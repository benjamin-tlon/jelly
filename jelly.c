#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "murmur3.h"
#include "xxh3.h"


// Internal murmur3 stuff //////////////////////////////////////////////////////

#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif

#define BIG_CONSTANT(x) (x##LLU)

/*
        Taken from murmur3.  Used to bit-mix direct atoms and
        shape-signatures.
*/
static FORCE_INLINE uint64_t fmix64 (uint64_t k) {
        k ^= k >> 33;
        k *= BIG_CONSTANT(0xff51afd7ed558ccd);
        k ^= k >> 33;
        k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
        k ^= k >> 33;

        return k;
}


/*
        Taken from stack overflow, derivative of boost::hash_combine,
        but for 64-bit hashes.
*/
size_t hash_combine(uint64_t x, uint64_t y) {
        x ^= y + BIG_CONSTANT(0x517cc1b727220a95) + (x << 6) + (x >> 2);
        return x;
}


// Misc Bin ////////////////////////////////////////////////////////////////////

/*
        What is the byte-width of a 64-bit word?

        For example: 0 takes 0 bytes, 127 takes 1 byte, 256 takes two
        bytes, etc.

        TODO: This should return 0 for 0, but it seems to return 1.
*/
uint32_t word64_bytes (uint64_t w) {
        uint32_t zeros = __builtin_clzll(w);
        uint32_t bits  = 64 - zeros;
        return (bits/8) + !!(bits%8);
}


// Types ///////////////////////////////////////////////////////////////////////

/*
    "Workspaces" are temporary `FanEntry`s used while building up a tree
    to be serialized.

    Pins : Map Hash256 Offset    <<  64 bits (word32 hash_fragment, word32 offset)
    Word : Map Word Offset       << 128 bits (word64 word, word64 offset)
    Nats : Map Natural Offset    << 192 bits (word64 leaves_hash, word64 length, word32 pointer)
    Bars : Map ByteString Offset << 192 bits (word64 leaves_hash, word64 length, word32 pointer)
*/

typedef struct hash256 {
    uint64_t a, b, c, d;
} hash256_t;

typedef struct { uint32_t ix; } pin_t;        //  Index into Jelly.pins
typedef struct { uint32_t ix; } bar_t;        //  Index into Jelly.bars
typedef struct { uint32_t ix; } nat_t;        //  Index into Jelly.nats
typedef struct { uint32_t ix; } fan_t;        //  Index into Jelly.fans
typedef struct { uint32_t ix; } treenode_t;   //  Index into Jelly.treenodes
typedef struct { uint32_t ix; } workspace_t;  //  Index into Jelly.workspaces

typedef struct fan_entry {
    uint32_t num_leaves;  // total number of interior (non-leaf) nodes in a tree.
    uint32_t num_bytes;   // total number of bytes in all leaves.
    uint64_t leaves_hash; // Summary hash of all leaf-data.
    uint64_t shape_hash;  // Summary hash of the tree shape.
    treenode_t pointer;   // If this is in a `workspaces` freelist,
                          // then this is instead an index into the
                          // workspaces array.
} FanEntry;

typedef struct treenode_value {
        uint64_t word;
} treenode_value;

typedef struct word_entry {
        uint64_t value;
        nat_t    index;
} word_entry_t;

/*
    10xx.. indicates a pin
    11xx.. indicates a bar
    0xxx.. indicates a nat
*/
typedef struct tagged_width {
        uint64_t word;
} tagged_width_t;

/*
        `offset` is an index into Jelly.pins, Jelly.bars, or Jelly.nats,
        depending on the tag in `tagged_width`

        If the byte-array width is less than nine, the data is directly
        inlined into `bytes`.
*/
typedef struct leaves_table_entry {
        uint64_t hash;
        tagged_width_t width;
        uint8_t *bytes;
        int32_t offset;
} leaves_table_entry_t;

/* A nat or a bar */
typedef struct leaf {
        int32_t width_bytes;
        uint8_t *bytes; // if width<9, then the bytes live in the pointer-space.
} leaf_t;

typedef struct jelly_ctx {
    // Ephemeral entries used during construction (recursive stack frames).
    FanEntry *workspaces;
    int32_t workspaces_width;
    int32_t workspaces_count;
    workspace_t workspaces_freelist;

    // One entry per tree-node
    treenode_value *treenodes;
    int32_t treenodes_width;
    int32_t treenodes_count;

    // Array of unique pin-leaves.
    hash256_t *pins;
    uint32_t pins_width;
    uint32_t pins_count;

    // Array of unique bar-leaves
    leaf_t   *bars;
    uint32_t bars_width;
    uint32_t bars_count;

    // Array of unique nat-leaves;
    leaf_t   *nats;
    uint32_t nats_width;
    uint32_t nats_count;

    // Leaf Deduplication table.
    leaves_table_entry_t *leaves_table;
    uint32_t leaves_table_width;
    uint32_t leaves_table_count;


    // Array of duplicate tree-nodes (and the top-level node).  Each one
    // is an index into `treenodes`.
    treenode_t *fans;
    uint32_t fans_width;
    uint32_t fans_count;
} Jelly;

#define TAG_PIN(x) (x | (4ULL << 61)) // 100...
#define TAG_FAN(x) (x | (7ULL << 61)) // 111...

FORCE_INLINE treenode_value TAG_BAR(bar_t bar) {
        uint64_t index = (uint64_t) bar.ix;
        uint64_t word  = (index | 5ULL << 61);
        return (treenode_value){ word };
}

FORCE_INLINE treenode_value TAG_NAT(nat_t nat) {
        uint64_t index = (uint64_t) nat.ix;
        uint64_t word  = (index | 6ULL << 61);
        return (treenode_value){ word };
}

// Pack two 32-bit indicies into one 64-bit words.  Since these indicies
// should both fit in 31 bits, we rely on the top-bit being set to zero
// and use that as a type-tag.
FORCE_INLINE treenode_value TAG_PAIR(treenode_t hed, treenode_t tel) {
        uint64_t hed_word  = (uint64_t) hed.ix;
        uint64_t tel_word  = (uint64_t) tel.ix;
        uint64_t result = ((hed_word << 32) | tel_word);
        return (treenode_value){ .word = result };
}

FORCE_INLINE uint64_t TREENODE_VAL_TAG(treenode_value v) {
        return (v.word >> 61);
}

FORCE_INLINE treenode_t TREENODE_VAL_HEAD(treenode_value v) {
        return (treenode_t){ .ix = (uint32_t) (v.word >> 32) };
}

// We rely on the cast to drop the high-bits.
FORCE_INLINE treenode_t TREENODE_VAL_TAIL(treenode_value v) {
        return (treenode_t){ .ix = (uint32_t) v.word };
}

// We rely on the cast to drop the high-bits.
FORCE_INLINE nat_t TREENODE_VAL_NAT(treenode_value v) {
    return (nat_t){ .ix = (uint32_t) v.word };
}

// We rely on the cast to drop the high-bits.
FORCE_INLINE bar_t TREENODE_VAL_BAR(treenode_value v) {
    return (bar_t){ .ix = (uint32_t) v.word };
}

// Memory Management ///////////////////////////////////////////////////////////

Jelly *new_jelly_ctx () {
        Jelly *res = malloc(sizeof(Jelly));
        res->workspaces          = malloc(32 * sizeof(FanEntry));
        res->workspaces_width    = 32;
        res->workspaces_count    = 0;
        res->workspaces_freelist = (workspace_t){-1};

        res->treenodes       = malloc(32 * sizeof(uint64_t));
        res->treenodes_width = 32;
        res->treenodes_count = 0;

        // Array of unique pin-leaves.
        res->pins       = malloc(8 * sizeof(hash256_t));
        res->pins_width = 8;
        res->pins_count = 0;

        // Array of unique bar-leaves
        res->bars       = malloc(8 * sizeof(leaf_t));
        res->bars_width = 8;
        res->bars_count = 0;

        // Array of unique nat-leaves
        res->nats       = malloc(8 * sizeof(leaf_t));
        res->nats_width = 8;
        res->nats_count = 0;

        // Deduplication table for leaves
        res->leaves_table = calloc(32, sizeof(leaves_table_entry_t));
        res->leaves_table_width = 32;
        res->leaves_table_count = 0;

        // Array of duplicate tree-nodes (and the top-level node).  Each one
        // is an index into `treenodes`.
        int32_t *fans      = malloc(64 * sizeof(uint32_t));
        int32_t fans_width = 64;
        int32_t fans_count = 0;

        return res;
}

void free_jelly_ctx (Jelly *ctx) {
        free(ctx->workspaces);
        free(ctx->treenodes);
        free(ctx->pins);
        free(ctx->bars);
        free(ctx->nats);
        free(ctx->leaves_table);
        free(ctx->fans);
        free(ctx);
}

FORCE_INLINE workspace_t alloc_workspace(Jelly *c) {
        /*
                (.pointer) is usually used to point into `treenodes`,
                but when it's in a freelist, we're using it to point to
                elements of the workspaces arena.

                That's why this block has all this weird casting.
        */
        if (c->workspaces_freelist.ix != -1) {
                int off = c->workspaces_freelist.ix;
                c->workspaces_freelist = (workspace_t){ c->workspaces[off].pointer.ix };
                return (workspace_t){ .ix = off };
        }

        uint32_t res = c->workspaces_count++;
        uint32_t wid = c->workspaces_width;

        if (res >= wid) {
                wid = wid + (wid / 2);
                c->workspaces = realloc(c->workspaces, wid * sizeof(FanEntry));
                c->workspaces_width = wid;
        }

        return (workspace_t){ .ix = res };
}

FORCE_INLINE void free_workspace(Jelly *c, workspace_t w) {
        c->workspaces[w.ix].pointer = (treenode_t){ c->workspaces_freelist.ix };
        c->workspaces_freelist = w;
        return;
}

FORCE_INLINE treenode_t alloc_treenode(Jelly *c, treenode_value v) {
        uint32_t res = c->treenodes_count++;
        uint32_t wid = c->treenodes_width;

        if (res >= wid) {
                wid = wid + (wid / 2);
                c->treenodes = reallocarray(c->treenodes, wid, sizeof(c->treenodes[0]));
                c->treenodes_width = wid;
        }

        c->treenodes[res] = v;

        return (treenode_t){ .ix = res };
}

FORCE_INLINE bar_t alloc_bar(Jelly *c) {
        uint32_t res = c->bars_count++;
        uint32_t wid = c->bars_width;

        if (res >= wid) {
                wid = wid + (wid / 2);
                c->bars = reallocarray(c->bars, wid, sizeof(c->bars[0]));
                c->bars_width = wid;
        }

        return (bar_t){ .ix = res };
}

FORCE_INLINE nat_t alloc_nat(Jelly *c) {
        uint32_t res = c->nats_count++;
        uint32_t wid = c->nats_width;

        if (res >= wid) {
                wid = wid + (wid / 2);
                c->nats = reallocarray(c->nats, wid, sizeof(c->nats[0]));
                c->nats_width = wid;
        }

        return (nat_t){ .ix = res };
}


// Tagged Widths ///////////////////////////////////////////////////////////////

/*
        0xxx.. indicates a nat
        10xx.. indicates a pin
        11xx.. indicates a bar
*/
FORCE_INLINE tagged_width_t NAT_TAGGED_WIDTH(uint32_t bytes) {
        uint64_t word = bytes;
        return (tagged_width_t){ .word = word };
}

FORCE_INLINE tagged_width_t PIN_TAGGED_WIDTH(uint32_t bytes) {
        uint64_t word = bytes;
        word |= 2ULL << 62;
        return (tagged_width_t){ .word = word };
}

FORCE_INLINE tagged_width_t BAR_TAGGED_WIDTH(uint32_t bytes) {
        uint64_t word = bytes;
        word |= 3ULL << 62;
        return (tagged_width_t){ .word = word };
}


// Inserting ///////////////////////////////////////////////////////////////////

void grow_leaves_if_full(Jelly *ctx) {
        uint32_t wid = ctx->leaves_table_width;
        uint32_t num = ctx->leaves_table_count;

        if ((num+1) >= wid) {
                wid *= 2;
                ctx->leaves_table_width = wid;
                ctx->leaves_table = reallocarray
                    ( ctx->leaves_table
                    , wid
                    , sizeof(leaves_table_entry_t)
                    );
                bzero( ctx->leaves_table + num
                     , (wid - num) * sizeof(leaves_table_entry_t)
                     );
        }

}

typedef struct leaf_insert_packed_t {
        tagged_width_t width;
        uint64_t hash;
        uint64_t word;
        uint32_t (*fill_entry)(Jelly*, leaf_t);
} LeafInsPacked;

uint32_t insert_packed(Jelly *ctx, LeafInsPacked ins) {
        uint64_t tag   = (ins.width.word >> 62);
        uint64_t bytes = (ins.width.word << 2) >> 2;

        printf(
            "\tinsert_packed(tag=%lu, wid=%lu, hash=%lu, %lu)\n",
            tag,
            bytes,
            ins.hash,
            ins.word
        );

        /*
                Do a linear-search over the leaves table, and return the
                index of the corresponding nat, if we find one.
        */

        grow_leaves_if_full(ctx); // Make sure there is space to insert,
                                  // if we do it later, we will invalidate
                                  // our pointer.

        leaves_table_entry_t *ent;

        for (ent = ctx->leaves_table; ent->hash; ent++) {
                bool match = ent->hash == ins.hash &&
                             ent->width.word == ins.width.word &&
                             ent->bytes == (uint8_t*) ins.word;

                if (match) {
                        printf("\t\tMATCH: n%d\n", ent->offset);
                        return ent->offset;
                }
        }

        ctx->leaves_table_count++;

        leaf_t leaf = (leaf_t){
                .width_bytes = bytes,
                .bytes = (uint8_t*) ins.word
        };

        uint32_t offset = ins.fill_entry(ctx, leaf);

        *ent = (leaves_table_entry_t){
            .hash  = ins.hash,
            .width = ins.width,
            .bytes = (uint8_t*) ins.word,
            .offset = offset,
        };
}

void print_nat_leaf(Jelly*, leaf_t);

int32_t bignat_insert(Jelly *ctx, bool isbar, uint64_t hash, leaf_t leaf) {
        printf("\tbignat_insert(wid=%u, ", leaf.width_bytes);
        print_nat_leaf(ctx, leaf);
        printf(")\n");

        /*
                Do a linear-search over the leaves table, and return the
                index of the corresponding nat, if we find one.
        */

        grow_leaves_if_full(ctx); // Make sure there is space to insert,
                                  // if we do it later, we will invalidate
                                  // our pointer.

        tagged_width_t tagged_width = isbar ? BAR_TAGGED_WIDTH(leaf.width_bytes)
                                            : NAT_TAGGED_WIDTH(leaf.width_bytes);
        leaves_table_entry_t *ent;

        for (ent = ctx->leaves_table; ent->hash; ent++) {
                bool hit = (ent->hash == hash) && (ent->width.word == tagged_width.word);
                if (!hit) continue;

                int o = memcmp(ent->bytes, leaf.bytes, leaf.width_bytes);
                if (o != 0) continue;

                printf("\t\tMATCH: n%d\n", ent->offset);
                return ent->offset;
        }

        ctx->leaves_table_count++;

        int32_t offset;
        char type;

        if (isbar) {
                bar_t bar = alloc_bar(ctx);
                ctx->bars[bar.ix] = leaf;
                offset = bar.ix;
                type = 'n';
        } else {
                nat_t nat = alloc_nat(ctx);
                ctx->nats[nat.ix] = leaf;
                offset = nat.ix;
                type = 'n';
        }

        printf("\t\tNEW_LEAF: %c%d\n", type, offset);

        *ent = (leaves_table_entry_t){
            .hash  = hash,
            .width = tagged_width,
            .bytes = leaf.bytes,
            .offset = offset,
        };

        return offset;
}

void push_fan_backref(Jelly *ctx, treenode_t tree) {
        uint32_t nex = ctx->fans_count++;
        uint32_t wid = ctx->fans_width;

        if (nex >= wid) {
                wid *= 2;
                ctx->fans = reallocarray(ctx->fans, wid, sizeof(ctx->fans[0]));
                ctx->fans_width = wid;
        }

        ctx->fans[nex] = tree;
}


// Building Up /////////////////////////////////////////////////////////////////

uint32_t fill_entry_word (Jelly *ctx, leaf_t leaf) {
        nat_t nat = alloc_nat(ctx);
        ctx->nats[nat.ix] = leaf;
        printf("\t\tNEW_LEAF: n%d\n", nat.ix);
        return nat.ix;
}

uint32_t fill_entry_bar (Jelly *ctx, leaf_t leaf) {
        bar_t bar = alloc_bar(ctx);
        ctx->bars[bar.ix] = leaf;
        printf("\t\tNEW_LEAF: b%d\n", bar.ix);
        return bar.ix;
}

/*
        TODO: Factor out the common code between packed nats and bars.

        I guess `fill_entry` could contain the alloc_treenode stuff as well?

        Have it return a treenode_value instead of an offset.

        And then the rest of the logic can live directly in
        `insert_packed`.
*/

workspace_t jelly_packed_nat(Jelly *ctx, uint32_t byte_width, uint64_t word) {
        printf("\tjelly_packed_nat(%lu, width=%u)\n", word, byte_width);

        uint64_t hash = fmix64(word);
        hash += !hash; // Zero hash means empty hashtable entry.

        tagged_width_t wid = NAT_TAGGED_WIDTH(byte_width);
        LeafInsPacked ins = { wid, hash, word, fill_entry_word };

        uint32_t offset = insert_packed(ctx, ins);
        treenode_value val = TAG_NAT((nat_t){ .ix = offset });

        treenode_t pointer = alloc_treenode(ctx, val);

        workspace_t res = alloc_workspace(ctx);
        FanEntry *ent = &(ctx->workspaces[res.ix]);

        ent->num_leaves  = 1;
        ent->num_bytes   = byte_width;
        ent->leaves_hash = hash;
        ent->shape_hash  = 0;
        ent->pointer     = pointer;

        return res;
}

workspace_t jelly_packed_bar(Jelly *ctx, uint32_t width_bytes, uint64_t word) {
        printf("\tjelly_packed_bar(%lu, width=%u)\n", word, width_bytes);

        uint64_t hash = fmix64(word);
        hash += !hash; // Zero hash means empty hashtable entry.

        tagged_width_t wid = BAR_TAGGED_WIDTH(width_bytes);
        LeafInsPacked ins = (LeafInsPacked){ wid, hash, word, fill_entry_bar };

        uint32_t offset = insert_packed(ctx, ins);

        treenode_value val = TAG_BAR((bar_t){offset});
        treenode_t pointer = alloc_treenode(ctx, val);

        workspace_t res = alloc_workspace(ctx);
        FanEntry *ent = &(ctx->workspaces[res.ix]);

        ent->num_leaves  = 1;
        ent->num_bytes   = width_bytes;
        ent->leaves_hash = hash;
        ent->shape_hash  = 0;
        ent->pointer     = pointer;

        return res;
}


/*
        TODO: Factor out the common bits, and split this into two
        functions, eliminating `isbar`.
*/
workspace_t jelly_nat(Jelly *ctx, bool isbar, leaf_t leaf) {
        printf("\tjelly_nat(width=%d, isbar=%u)\n", leaf.width_bytes, (unsigned)isbar);

        if (leaf.width_bytes < 9) {
                if (isbar) {
                        return jelly_packed_bar(ctx, leaf.width_bytes, (uint64_t)leaf.bytes);
                } else {
                        return jelly_packed_nat(ctx, leaf.width_bytes, (uint64_t)leaf.bytes);
                }
        }

        uint64_t hash = XXH3_64bits(leaf.bytes, leaf.width_bytes);
        hash += !hash;

        int32_t offset  = bignat_insert(ctx, isbar, hash, leaf);

        treenode_value val = isbar ? TAG_BAR((bar_t){ .ix = offset })
                                   : TAG_NAT((nat_t){ .ix = offset });

        treenode_t node = alloc_treenode(ctx, val);
        workspace_t res = alloc_workspace(ctx);
        FanEntry *ent   = &(ctx->workspaces[res.ix]);

        ent->num_leaves  = 1;
        ent->num_bytes   = leaf.width_bytes;
        ent->leaves_hash = hash;
        ent->shape_hash  = 0;
        ent->pointer     = node;

        return res;
}

workspace_t jelly_cons(Jelly *ctx, workspace_t hed, workspace_t tel) {
        printf("\tjelly_cons(%i, %i)\n", hed.ix, tel.ix);

        FanEntry h = ctx->workspaces[hed.ix];
        FanEntry t = ctx->workspaces[tel.ix];

        uint32_t num_leaves  = h.num_leaves + t.num_leaves;
        uint32_t num_bytes   = h.num_bytes  + t.num_bytes;
        uint64_t leaves_hash = hash_combine(h.leaves_hash, t.leaves_hash);
        uint64_t shape_width = (num_leaves * 2) - 1;

        /*
           The shape-hash is the outline of the tree.  Cells are 1-bits and
           leaves are 0-bits.

               The shape of 9:             0b0
               The shape of (9 8):         0b100
               The shape of ((9 8) 7):     0b11000
               The shape of ((9 8) (7 6)): 0b1100100

           These shapes quickly exceed our 64-bit word sizes (as soon as we
           are working with sub-trees that contain 33 leaves or more).
           At that point, we shuffle the bits using Murmur3's fmix64 and
           switch to just using a hash-combination routine.
        */
        uint32_t head_bits = (h.num_leaves * 2) - 1;
        uint32_t tail_bits = (t.num_leaves * 2) - 1;
        uint64_t shape_hash;
        if (shape_width > 64) {
                if (head_bits <= 64) h.shape_hash = fmix64(h.shape_hash);
                if (tail_bits <= 64) t.shape_hash = fmix64(t.shape_hash);
                shape_hash = hash_combine(h.shape_hash, t.shape_hash);
        } else {
                shape_hash = (1ULL << (shape_width-1))
                           | (h.shape_hash << tail_bits)
                           | t.shape_hash
                           ;
        }

        treenode_t pointer = alloc_treenode(ctx, TAG_PAIR(h.pointer, t.pointer));

        // Since we "consume" both of the workspace entries we've been given,
        // we can directly re-use the memory from the head.  The tail is
        // released.

        FanEntry *target = ctx->workspaces + hed.ix;
        target->num_leaves  = num_leaves;
        target->num_bytes   = num_bytes;
        target->leaves_hash = leaves_hash;
        target->shape_hash  = shape_hash;
        target->pointer     = pointer;

        free_workspace(ctx, tel);

        return hed;
}


// Testing /////////////////////////////////////////////////////////////////////

#define die(...) (fprintf(stderr, __VA_ARGS__),exit(1))

uint8_t hex_value(char a, char b) {
    if (b == EOF) die("partial hex literal");

    uint8_t top_byte = (isdigit(a) ? (a - '0') : (tolower(a) - ('a' - 10)));
    uint8_t low_byte = (isdigit(b) ? (b - '0') : (tolower(b) - ('a' - 10)));
    uint8_t result   = (top_byte << 4) | low_byte;

    return result;
}

uint64_t pack_bytes_msb(size_t width, char *bytes) {
        uint64_t res = 0;
        for (int i=0; i<width; i++) {
                uint64_t tmp = (uint8_t) bytes[i];
                res = (res<<8) | tmp;
        }
        return res;
}

uint64_t pack_bytes_lsb(size_t width, char *bytes) {
        printf("\t\tpack_bytes_lsb(width=%lu, ", width);
        uint64_t res = 0;
        for (int i=0; i<width; i++) {
                uint64_t tmp = bytes[i];
                printf("%02lx(%d)", tmp, i*8);
                res = res | (tmp << (i*8));
        }
        printf(")\n");
        return res;
}

workspace_t read_one(Jelly*);
workspace_t read_many(Jelly*, workspace_t);
workspace_t read_some(Jelly*);

leaf_t read_hex() {
        size_t count=0, wid=256;
        char tmp[256], *buf=tmp;

    loop:

        char c = getchar();

        if (isalnum(c)) {
                char d = getchar();

                uint8_t byte = hex_value(c, d);

                if (count >= wid) {
                        wid *= 4;
                        buf = ((tmp==buf) ? malloc(wid) : realloc(buf,wid));
                }

                buf[count++] = byte;

                goto loop;
        }

        ungetc(c, stdin);

        leaf_t result = { .width_bytes = count, .bytes = NULL };

        if (count < 8) {
                uint8_t **hack_slot = &(result.bytes);
                *((uint64_t*)hack_slot) = pack_bytes_msb(count, buf);
        } else {
                result.bytes = malloc(count);
                for (int i=0, j=(count-1); i<count; i++, j--) {
                    result.bytes[i] = buf[j];
                }
        }

        if (buf != tmp) free(buf);

        return result;
}

uint64_t read_word(uint64_t acc) {
        printf("\tread_word(%lu)\n", acc);
        int c;

        while (isdigit(c = getchar())) {
                printf("\t    acc=%lu\t('%c' is a digit)\n", acc, c);
                acc *= 10;
                acc += (uint64_t)(c - '0');
        }

        printf("\t    acc=%lu\t('%c' is not a digit)\n", acc, c);

        if (c != EOF) ungetc(c, stdin);
        return acc;
}

leaf_t read_string() {
        char *buf = malloc(1024);
        size_t wid = 1024;
        int ix = 0;

    loop:

        char c = getchar();
        switch (c) {
            case '"':
            case EOF:
                goto end;
            case '\\': {
                char d = getchar();
                if (d == EOF) goto end;
                c = d;
                break;
            }
        }

        buf[ix++] = c;
        if (ix >= wid) { wid *= 2; buf = realloc(buf, wid); }
        goto loop;

    end:

        size_t count = ix;
        leaf_t leaf;
        leaf.width_bytes = count;

        printf("\tread_string() -> \"%s\" (%lu)\n", buf, count);


        if (count < 9) {
                leaf.bytes = (uint8_t*) pack_bytes_lsb(count, buf);
        } else {
                leaf.bytes = calloc(count+1, 1);
                memcpy(leaf.bytes, buf, count);
        }

        free(buf);
        return leaf;

}

workspace_t read_some(Jelly *ctx) {
        return read_many(ctx, read_one(ctx));
}

workspace_t read_many(Jelly *ctx, workspace_t acc) {
    loop:
        int c = getchar();

        switch (c) {
            case '"': {
                workspace_t elmt = jelly_nat(ctx, true, read_string());
                acc = jelly_cons(ctx, acc, elmt);
                goto loop;
            }
            case '(':
                acc = jelly_cons(ctx, acc, read_many(ctx, read_one(ctx)));
                goto loop;
            case EOF:
            case ')':
                return acc;
            case ' ':
            case '\n':
            case '\r':
            case '\t':
                goto loop;
            case '0': {
                // TODO Break this out so that it can be implemented in
                // `read_one` also.
                c = getchar();
                if (c == 'x') {
                        workspace_t elmt = jelly_nat(ctx, false, read_hex());
                        acc = jelly_cons(ctx, acc, elmt);
                        goto loop;
                } else if (isdigit(c)) {
                        die("Non-zero numbers must not be zero-prefixed");
                } else {
                        ungetc(c, stdin);
                        uint32_t width = word64_bytes(0);
                        workspace_t elmt = jelly_packed_nat(ctx, width, 0);
                        acc = jelly_cons(ctx, acc, elmt);
                        goto loop;
                }
            }
            default:
                if isdigit(c) {
                        uint64_t word = read_word((uint64_t)(c - '0'));
                        uint32_t width = word64_bytes(word);
                        workspace_t elmt = jelly_packed_nat(ctx, width, word);
                        acc = jelly_cons(ctx, acc, elmt);
                        goto loop;
                } else {
                        die("Unexpected character: %c (read_many)\n", c);
                }
        }
}


workspace_t read_one(Jelly *ctx) {
        int c = getchar();
    loop:
        if (isdigit(c)) {
                uint64_t word = read_word((uint64_t)(c - '0'));
                uint32_t width = word64_bytes(word);
                return jelly_packed_nat(ctx, width, word);
        }

        switch (c) {
            case '"':
                return jelly_nat(ctx, true, read_string());
            case '(':
                return read_many(ctx, read_one(ctx));
            case ' ':
            case '\n':
            case '\r':
            case '\t':
                goto loop;
            case EOF:
                die("Unexpected EOF (read_one)\n");
            default:
                die("Unexpected character: %c (read_one)\n", c);
        }
}

int jelly_buffer_size(Jelly *ctx) {
        return 1 + strlen("\n\nGet fucked, lol\n\n");
}

void jelly_dump(Jelly *ctx, char *buf) {
        strcpy(buf, "\n\nGet fucked, lol\n\n");
}

void jelly_push_final(Jelly *ctx, workspace_t val) {
        treenode_t tree = ctx->workspaces[val.ix].pointer;
        push_fan_backref(ctx, tree);
}


void print_nat_leaf(Jelly *ctx, leaf_t l) {
        if (l.width_bytes > 8) {
                printf("0x");
                for (int i = l.width_bytes - 1; i>=0; i--) {
                        uint8_t byte = l.bytes[i];
                        printf("%02x", (unsigned) byte);
                }
        } else {
                uint64_t value = (uint64_t) l.bytes;
                printf("%lu", value);
        }
}

void print_nat(Jelly *ctx, nat_t nat) {
        leaf_t l = ctx->nats[nat.ix];
        if (l.width_bytes > 8) {
                printf("0x");
                for (int i = l.width_bytes - 1; i>=0; i--) {
                        uint8_t byte = l.bytes[i];
                        printf("%02x", (unsigned) byte);
                }
        } else {
                uint64_t value = (uint64_t) l.bytes;
                printf("%lu", value);
        }
}

void print_bar(Jelly *ctx, bar_t bar) {
        leaf_t l = ctx->bars[bar.ix];
        char *buf = 0;
        if (l.width_bytes > 8) {
                buf = l.bytes;
        } else {
                buf = (char*) &(l.bytes);
        }

        printf("\"");
        for (int i = 0; i < l.width_bytes; i++) {
                uint8_t byte = buf[i];
                printf("%c", byte);
        }
        printf("\"");
}


void print_tree(Jelly*, treenode_t);

void print_tree_list(Jelly *ctx, treenode_t tree) {
        treenode_value val = ctx->treenodes[tree.ix];

        switch (TREENODE_VAL_TAG(val)) {
            case 4: printf("pin"); break;
            case 5: print_bar(ctx, TREENODE_VAL_BAR(val)); break;
            case 6: print_nat(ctx, TREENODE_VAL_NAT(val)); break;
            case 7: printf("fan"); break;
            default: {
                treenode_t hed = TREENODE_VAL_HEAD(val);
                treenode_t tel = TREENODE_VAL_TAIL(val);
                print_tree_list(ctx, hed);
                putchar(' ');
                print_tree(ctx, tel);
            }
        }
}

void print_tree(Jelly *ctx, treenode_t tree) {
        treenode_value val = ctx->treenodes[tree.ix];

        switch (TREENODE_VAL_TAG(val)) {
            case 4: printf("pin"); break;
            case 5: print_bar(ctx, TREENODE_VAL_BAR(val)); break;
            case 6: print_nat(ctx, TREENODE_VAL_NAT(val)); break;
            case 7: printf("fan"); break;
            default: {
                treenode_t hed = TREENODE_VAL_HEAD(val);
                treenode_t tel = TREENODE_VAL_TAIL(val);
                putchar('(');
                print_tree_list(ctx, hed);
                putchar(' ');
                print_tree(ctx, tel);
                putchar(')');
            }
        }
}

int treenode_depth(Jelly *ctx, treenode_t node) {
        treenode_value val = ctx->treenodes[node.ix];
        if (val.word >> 63) return 1;
        int hed_depth = treenode_depth(ctx, TREENODE_VAL_HEAD(val));
        int tel_depth = treenode_depth(ctx, TREENODE_VAL_TAIL(val));
        return (1 + ((hed_depth > tel_depth) ? hed_depth : tel_depth));
}

void jelly_debug(Jelly *ctx) {
        printf("\njelly_debug():\n");

        {
                int free_workspaces = 0;
                int32_t tmp = ctx->workspaces_freelist.ix;
                while (tmp != -1) {
                        tmp = ctx->workspaces[tmp].pointer.ix;
                        free_workspaces++;
                }
                printf("\n\tworkspaces: (width=%u count=%u free=%d)\n",
                       ctx->workspaces_width,
                       ctx->workspaces_count,
                       free_workspaces);
        }

        {
                printf("\n\tleaves: (width=%u, count=%u)\n\n",
                       ctx->leaves_table_width,
                       ctx->leaves_table_count);

                int num = ctx->leaves_table_count;

                for (int i=0; i<num; i++) {
                        leaves_table_entry_t ent = ctx->leaves_table[i];

                        uint64_t width = ((ent.width.word << 2) >> 2);

                        char type;
                        switch (ent.width.word >> 62) {
                            case 2:  type = 'p'; break;
                            case 3:  type = 'b'; break;
                            default: type = 'n'; break;
                        }

                        printf("\t\t[%04d] %c%d = ", i, type, ent.offset);
                        if (type == 'n') {
                                print_nat(ctx, (nat_t){ .ix = ent.offset });
                        } else if (type == 'b') {
                                print_bar(ctx, (bar_t){ .ix = ent.offset });
                        } else {
                                die("todo: pins\n");
                        }
                        if (width < 9) {
                                uint64_t word = (uint64_t) ctx->nats[i].bytes;
                                if (word < 999) printf("\t");
                        }
                        printf("\t(hash=%lu, width=%lu)\n", ent.hash, width);
                }
        }

        uint32_t count = ctx->fans_count;

        {
                printf("\n\tbars: (width=%u, count=%u)\n\n",
                       ctx->bars_width,
                       ctx->bars_count);
                int num = ctx->bars_count;
                for (int i=0; i<num; i++) {
                        printf("\t\tn%d = ", i);
                        print_bar(ctx, (bar_t){ .ix = i });
                        printf("\n");
                }
        }


        {
                printf("\n\tnats: (width=%u, count=%u)\n\n",
                       ctx->nats_width,
                       ctx->nats_count);
                int num = ctx->nats_count;
                for (int i=0; i<num; i++) {
                        printf("\t\tn%d = ", i);
                        print_nat(ctx, (nat_t){ .ix = i });
                        printf("\n");
                }
        }

        {
                printf("\n\ttreenodes: (width=%u, count=%u)\n\n",
                       ctx->treenodes_width,
                       ctx->treenodes_count);
                int num = ctx->treenodes_count;
                for (int i=0; i<num; i++) {
                        treenode_value val = ctx->treenodes[i];
                        switch (val.word >> 61) {
                            case 4:
                                printf("\t\t[%04d]: p%u\n", i, (uint32_t) val.word);
                                continue;
                            case 5:
                                printf("\t\t[%04d]: b%u\n", i, (uint32_t) val.word);
                                continue;
                            case 6:
                                printf("\t\t[%04d]: n%u\n", i, (uint32_t) val.word);
                                continue;
                            case 7:
                                printf("\t\t[%04d]: f%u\n", i, (uint32_t) val.word);
                                continue;
                            default: {
                                treenode_t hed = TREENODE_VAL_HEAD(val);
                                treenode_t tel = TREENODE_VAL_TAIL(val);
                                printf("\t\t[%04d]: (%u, %u)\n", i, hed.ix, tel.ix);
                                continue;
                            }
                        }
                }
        }

        if (count == 0) die("jelly_debug: no value to print\n");
        if (count > 1) die("TODO: Print backreferences");

        treenode_t last = ctx->fans[count - 1];
        int depth = treenode_depth(ctx, last);
        printf("\n\tfinal value (treenode=%d depth=%d):\n\n\t\t", last.ix, depth);
        print_tree(ctx, last);
}


int main () {
        Jelly *ctx = new_jelly_ctx();

        workspace_t top = read_some(ctx);
        jelly_push_final(ctx, top);
        free_workspace(ctx, top);

        int wid = jelly_buffer_size(ctx);
        uint8_t *buf = calloc(wid, sizeof(uint8_t));

        jelly_debug(ctx);

        jelly_dump(ctx, buf);
        fwrite(buf, 1, wid, stdout);
        free_jelly_ctx(ctx);
        free(buf);

        return 0;
}
