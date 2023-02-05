#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "murmur3.h"
#include "xxh3.h"
#include "libbase58.h"


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

        k += !k; // Never return 0

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

#define die(...) (fprintf(stderr, __VA_ARGS__),exit(1))

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
    "deks" are temporary `FanEntry`s used while building up a tree to
    be serialized.

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

typedef struct treenode_value {
        uint64_t word;
} treenode_value;

typedef struct {
    treenode_value val;   // (Word32, Word32)
    uint32_t hash;        // Bit mixed + truncated version of `val.word`
    treenode_t pointer;   // If this is in a `desks` freelist,
                          // then this is instead an index into the
                          // desks array.
} FanEntry;

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
        If the byte-array width is less than nine, the data is directly
        inlined into `bytes`.
*/
typedef struct leaves_table_entry {
        uint64_t hash;
        tagged_width_t width;
        uint8_t *bytes;
        treenode_t pointer;
} leaves_table_entry_t;

/* A nat or a bar */
typedef struct leaf {
        int32_t width_bytes;
        uint8_t *bytes; // if width<9, then the bytes live in the pointer-space.
} leaf_t;

typedef struct jelly_ctx {
    // One entry per tree-node
    treenode_value *treenodes;
    int32_t treenodes_width;
    int32_t treenodes_count;

    // Array of unique pin-leaves.
    hash256_t **pins;
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
    uint64_t leaves_table_mask;

    // Interior-Node Deduplication table.
    FanEntry *nodes_table;
    uint32_t nodes_table_width;
    uint32_t nodes_table_count;
    uint64_t nodes_table_mask;

    // Array of duplicate tree-nodes (and the top-level node).  Each one
    // is an index into `treenodes`.
    treenode_t *fans;
    uint32_t fans_width;
    uint32_t fans_count;
} Jelly;

void print_tree(Jelly*, treenode_t);


FORCE_INLINE treenode_value TAG_FAN(fan_t fan) {
        uint64_t index = (uint64_t) fan.ix;
        uint64_t word  = (index | 7ULL << 61);
        return (treenode_value){ word };
}

FORCE_INLINE treenode_value TAG_PIN(pin_t pin) {
        uint64_t index = (uint64_t) pin.ix;
        uint64_t word  = (index | 4ULL << 61);
        return (treenode_value){ word };
}

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

FORCE_INLINE uint64_t NODEVAL_TAG(treenode_value v) {
        return (v.word >> 61);
}

FORCE_INLINE treenode_t NODEVAL_HEAD(treenode_value v) {
        return (treenode_t){ .ix = (uint32_t) (v.word >> 32) };
}

// We rely on the cast to drop the high-bits.
FORCE_INLINE treenode_t NODEVAL_TAIL(treenode_value v) {
        return (treenode_t){ .ix = (uint32_t) v.word };
}

// We rely on the cast to drop the high-bits.
FORCE_INLINE nat_t NODEVAL_NAT(treenode_value v) {
        return (nat_t){ .ix = (uint32_t) v.word };
}

// We rely on the cast to drop the high-bits.
FORCE_INLINE bar_t NODEVAL_BAR(treenode_value v) {
        return (bar_t){ .ix = (uint32_t) v.word };
}

// We rely on the cast to drop the high-bits.
FORCE_INLINE pin_t NODEVAL_PIN(treenode_value v) {
        return (pin_t){ .ix = (uint32_t) v.word };
}

// Memory Management ///////////////////////////////////////////////////////////

Jelly *new_jelly_ctx () {
        Jelly *res = malloc(sizeof(Jelly));

        // One entry per unique node (both leaves and interior nodes)
        res->treenodes       = calloc(32, sizeof(res->treenodes[0]));
        res->treenodes_width = 32;
        res->treenodes_count = 0;

        // Array of unique pin-leaves.
        res->pins       = calloc(8, sizeof(res->pins[0]));
        res->pins_width = 8;
        res->pins_count = 0;

        // Array of unique bar-leaves
        res->bars       = calloc(8, sizeof(res->bars[0]));
        res->bars_width = 8;
        res->bars_count = 0;

        // Array of unique nat-leaves
        res->nats       = calloc(8, sizeof(res->nats[0]));
        res->nats_width = 8;
        res->nats_count = 0;

        // Deduplication table for leaves
        res->leaves_table       = calloc((1<<2), sizeof(res->leaves_table[0]));
        res->leaves_table_width = (1<<2);
        res->leaves_table_mask  = ((1<<2) - 1);
        res->leaves_table_count = 0;

        // Deduplication table for interior nodes
        size_t nodes_table_byte_width = (1<<2) * sizeof(res->nodes_table[0]);
        res->nodes_table       = malloc(nodes_table_byte_width);
        res->nodes_table_width = (1<<2);
        res->nodes_table_mask  = ((1<<2) - 1);
        res->nodes_table_count = 0;

        memset(res->nodes_table, 255, nodes_table_byte_width);

        // Array of duplicate tree-nodes (and the top-level node).  Each one
        // is an index into `treenodes`.
        res->fans       = calloc(16, sizeof(res->fans[0]));
        res->fans_width = 16;
        res->fans_count = 0;

        return res;
}

void free_jelly_ctx (Jelly *ctx) {
        free(ctx->treenodes);
        free(ctx->pins);
        free(ctx->bars);
        free(ctx->nats);
        free(ctx->leaves_table);
        free(ctx->nodes_table);
        free(ctx->fans);
        free(ctx);
}


FORCE_INLINE treenode_t alloc_treenode(Jelly *c, treenode_value v) {
        uint32_t res = c->treenodes_count++;
        uint32_t wid = c->treenodes_width;

        if (res >= wid) {
                wid *= 2;
                c->treenodes = reallocarray(c->treenodes, wid, sizeof(c->treenodes[0]));
                c->treenodes_width = wid;
        }

        c->treenodes[res] = v;

        return (treenode_t){ .ix = res };
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

typedef struct {
        tagged_width_t width;
        uint64_t hash;
        uint64_t word;
        treenode_t (*new_leaf)(Jelly*, leaf_t);
} PackedInsertRequest;

typedef struct {
        tagged_width_t width;
        uint64_t hash;
        leaf_t leaf;
        treenode_t (*new_leaf)(Jelly*, leaf_t);
} IndirectInsertRequest;

void jelly_debug_leaves(Jelly*, bool);
void jelly_debug_interior_nodes(Jelly*);

void rehash_nodes_if_full(Jelly *ctx) {
        uint32_t oldwid = ctx->nodes_table_width;
        uint32_t num    = ctx->nodes_table_count;

        // If capacity is >= 50%, resize.
        if (num*2 < oldwid) return;

        uint32_t newwid = oldwid*2;

        printf("rehash (old=%u new=%u)\n", oldwid, newwid);

        // printf("OLD TABLE:\n");
        // jelly_debug_interior_nodes(ctx);
        // printf("\n");

        FanEntry *oldtab = ctx->nodes_table;
        FanEntry *newtab = malloc(newwid * sizeof(*newtab));

        memset(newtab, 255, newwid * sizeof(*newtab));

        // e.g. 0b111 -> 0b1111
        uint64_t newmask = (ctx->nodes_table_mask << 1) | 1;

        /*
                Loop over the whole of oldtab, and re-insert every
                non-empty slot.  Inserts are guarenteed to be unique,
                so we just need to, starting at the correct bucket,
                scan for an empty slot and write the value there.
        */
        for (uint64_t i=0; i<oldwid; i++) {
                FanEntry ent = oldtab[i];

                uint64_t j = ent.hash;

                if (ent.val.word == UINT64_MAX) continue; // empty slot

                for (;; j++) {
                        j &= newmask;
                        FanEntry *tar = newtab + j;
                        if (tar->val.word == UINT64_MAX) {
                                *tar = ent;
                                // printf("\t\t%lu -> %lu\n", i, j);
                                break;
                        } // else {
                                // printf("\t\t\t(%lu -> %lu) is taken\n", i, j);
                        // }
                }
        }

        free(ctx->nodes_table);

        ctx->nodes_table_width = newwid;
        ctx->nodes_table_mask  = newmask;
        ctx->nodes_table       = newtab;

        // printf("NEW TABLE:\n");
        // jelly_debug_interior_nodes(ctx);
        // printf("\n");
}

void rehash_leaves_if_full(Jelly *ctx) {
        uint32_t oldwid = ctx->leaves_table_width;
        uint32_t num    = ctx->leaves_table_count;

        // If capacity is >= 50%, resize.
        if (num*2 < oldwid) return;

        uint32_t newwid = oldwid*2;

        printf("rehash (old=%u new=%u)\n", oldwid, newwid);

        // printf("OLD TABLE:\n");
        // jelly_debug_leaves(ctx, false);
        // printf("\n");

        leaves_table_entry_t *oldtab = ctx->leaves_table;
        leaves_table_entry_t *newtab = calloc(newwid, sizeof(*newtab));

        // printf("\tcalloc(%u * %lu = %lu)\n", newwid, sizeof(*newtab), newwid*sizeof(*newtab));

        // e.g. 0b111 -> 0b1111
        uint64_t newmask = (ctx->leaves_table_mask << 1) | 1;

        /*
                Loop over the whole of oldtab, and re-insert every
                non-empty slot.  Inserts are guarenteed to be unique,
                so we just need to, starting at the correct bucket,
                scan for an empty slot and write the value there.
        */
        for (uint64_t i=0; i<oldwid; i++) {
                leaves_table_entry_t ent = oldtab[i];

                uint64_t j = ent.hash;

                if (j == 0) continue; // empty slot

                for (;; j++) {
                        j &= newmask;
                        leaves_table_entry_t *tar = newtab + j;
                        if (tar->hash == 0) {
                                *tar = ent;
                                // printf("\t\t%lu -> %lu\n", i, j);
                                break;
                        } // else {
                                // printf("\t\t\t(%lu -> %lu) is taken\n", i, j);
                        // }
                }
        }

        free(ctx->leaves_table);

        ctx->leaves_table_width = newwid;
        ctx->leaves_table_mask  = newmask;
        ctx->leaves_table       = newtab;

        // printf("NEW TABLE:\n");
        // jelly_debug_leaves(ctx, false);
        // printf("\n");
}


treenode_t insert_packed_leaf(Jelly *ctx, PackedInsertRequest req) {
        uint64_t tag   = (req.width.word >> 62);
        uint64_t bytes = (req.width.word << 2) >> 2;

        printf(
            "\tpacked_insert(tag=%lu, wid=%lu, hash=%lu, %lu)\n",
            tag,
            bytes,
            req.hash,
            req.word
        );

        if (req.hash == 0) {
                die("Broken invariant: hashes must be non-zero\n");
        }

        /*
                Do a linear-search over the leaves table, and use the
                index of the corresponding nat, if we find one.

                We make sure ahead-of-time that there is enough space
                for us to insert.  This way we don't have to move things
                around while were are inserting.
        */

        rehash_leaves_if_full(ctx);

        leaves_table_entry_t *ent = NULL;

        uint64_t mask = ctx->leaves_table_mask;
        uint64_t ix = req.hash;

        for (;; ix++) {
                ix &= mask;

                ent = &(ctx->leaves_table[ix]);

                if (ent->hash == 0) break;

                bool match = ent->hash == req.hash &&
                             ent->width.word == req.width.word &&
                             ent->bytes == (uint8_t*) req.word;

                if (!match) continue;

                printf("\t\tPACKED_LEAF_MATCH: ");
                print_tree(ctx, ent->pointer);
                printf("\n");

                return ent->pointer;
        }

        /*
                We didn't find any matching entries, but now `ent`
                is pointing into an empty slot, so we fill that and
                return it.
        */

        ctx->leaves_table_count++;

        leaf_t leaf = (leaf_t){
                .width_bytes = bytes,
                .bytes = (uint8_t*) req.word
        };

        treenode_t pointer = req.new_leaf(ctx, leaf);

        printf("\t\tNEW_PACKED_LEAF: ");
        print_tree(ctx, pointer);
        printf("\n");

        *ent = (leaves_table_entry_t){
            .hash  = req.hash,
            .width = req.width,
            .bytes = (uint8_t*) req.word,
            .pointer = pointer,
        };

        return pointer;
}

void print_nat_leaf(Jelly*, leaf_t);

treenode_t insert_indirect_leaf(Jelly *ctx, IndirectInsertRequest req) {
        printf("\tinsert_indirect_leaf(wid=%u, ", req.leaf.width_bytes);
        print_nat_leaf(ctx, req.leaf);
        printf(")\n");

        if (req.hash == 0) {
                die("Broken invariant: hashes must be non-zero\n");
        }

        /*
                Do a linear-search over the leaves table, and return the
                index of the corresponding nat, if we find one.
        */

        rehash_leaves_if_full(ctx); // Make sure there is space to insert,
                                    // if we do it later, we will invalidate
                                    // our pointer.

        leaves_table_entry_t *ent;

        uint64_t mask = ctx->leaves_table_mask;
        uint64_t ix = req.hash;

        for (;; ix++) {
                ix &= mask;

                ent = &(ctx->leaves_table[ix]);

                if (ent->hash == 0) break;

                bool hit = (ent->hash == req.hash) && (ent->width.word == req.width.word);
                if (!hit) continue;

                int o = memcmp(ent->bytes, req.leaf.bytes, req.leaf.width_bytes);
                if (o != 0) continue;

                printf("\t\tINDIRECT_LEAF_MATCH: ");
                print_tree(ctx, ent->pointer);
                printf("\n");

                return ent->pointer;
        }

        ctx->leaves_table_count++;

        treenode_t pointer = req.new_leaf(ctx, req.leaf);

        printf("\t\tNEW_INDIRECT_LEAF: ");
        print_tree(ctx, pointer);
        printf("\n");

        *ent = (leaves_table_entry_t){
            .hash  = req.hash,
            .width = req.width,
            .bytes = req.leaf.bytes,
            .pointer = pointer,
        };

        return pointer;
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


// Inserting Nats //////////////////////////////////////////////////////////////

FORCE_INLINE nat_t alloc_nat(Jelly *c) {
        uint32_t res = c->nats_count++;
        uint32_t wid = c->nats_width;

        if (res >= wid) {
                wid *= 2;
                c->nats = reallocarray(c->nats, wid, sizeof(c->nats[0]));
                c->nats_width = wid;
        }

        return (nat_t){ .ix = res };
}

static treenode_t
new_nat (Jelly *ctx, leaf_t leaf) {
        nat_t nat = alloc_nat(ctx);
        ctx->nats[nat.ix] = leaf;
        return alloc_treenode(ctx, TAG_NAT(nat));
}

FORCE_INLINE static treenode_t
jelly_packed_nat(Jelly *ctx, uint32_t byte_width, uint64_t word) {
        printf("\tjelly_packed_nat(%lu, width=%u)\n", word, byte_width);

        return insert_packed_leaf(ctx,
            (PackedInsertRequest){
                .width = NAT_TAGGED_WIDTH(byte_width),
                .hash = fmix64(word),
                .word = word,
                .new_leaf = new_nat,
            }
        );
}

treenode_t jelly_nat(Jelly *ctx, leaf_t leaf) {
        printf("\tjelly_nat(width=%d)\n", leaf.width_bytes);

        if (leaf.width_bytes < 9) {
                return jelly_packed_nat(ctx, leaf.width_bytes, (uint64_t)leaf.bytes);
        } else {
                uint64_t hash = XXH3_64bits(leaf.bytes, leaf.width_bytes);
                hash += !hash;

                return insert_indirect_leaf(ctx,
                    (IndirectInsertRequest){
                        .width = NAT_TAGGED_WIDTH(leaf.width_bytes),
                        .hash = hash,
                        .leaf = leaf,
                        .new_leaf = new_nat,
                    }
                );
        }
}




// Inserting Bars //////////////////////////////////////////////////////////////

FORCE_INLINE bar_t alloc_bar(Jelly *c) {
        uint32_t res = c->bars_count++;
        uint32_t wid = c->bars_width;

        if (res >= wid) {
                wid *= wid;
                c->bars = reallocarray(c->bars, wid, sizeof(c->bars[0]));
                c->bars_width = wid;
        }

        return (bar_t){ .ix = res };
}

static treenode_t
new_bar (Jelly *ctx, leaf_t leaf) {
        bar_t bar = alloc_bar(ctx);
        ctx->bars[bar.ix] = leaf;
        return alloc_treenode(ctx, TAG_BAR(bar));
}

treenode_t jelly_bar(Jelly *ctx, leaf_t leaf) {
        printf("\tjelly_bar(width=%d)\n", leaf.width_bytes);

        if (leaf.width_bytes < 9) {
                uint64_t word = (uint64_t) leaf.bytes;

                printf("\tjelly_packed_bar(%lu, width=%u)\n", word, leaf.width_bytes);

                return insert_packed_leaf(ctx,
                    (PackedInsertRequest){
                        .width = BAR_TAGGED_WIDTH(leaf.width_bytes),
                        .hash = fmix64(word),
                        .word = word,
                        .new_leaf = new_bar,
                    }
                );
        } else {
                uint64_t hash = XXH3_64bits(leaf.bytes, leaf.width_bytes);
                hash += !hash;

                return insert_indirect_leaf(ctx,
                    (IndirectInsertRequest){
                        .width = BAR_TAGGED_WIDTH(leaf.width_bytes),
                        .hash = hash,
                        .leaf = leaf,
                        .new_leaf = new_bar,
                    }
                );
        }
}


// Inserting Pins //////////////////////////////////////////////////////////////

FORCE_INLINE pin_t alloc_pin(Jelly *c) {
        uint32_t res = c->pins_count++;
        uint32_t wid = c->pins_width;

        if (res >= wid) {
                wid *= 2;
                c->pins = reallocarray(c->pins, wid, sizeof(c->pins[0]));
                c->pins_width = wid;
        }

        return (pin_t){ .ix = res };
}

static treenode_t new_pin (Jelly *ctx, leaf_t leaf) {
        hash256_t *hash = (hash256_t*) leaf.bytes;
        pin_t pin = alloc_pin(ctx);
        ctx->pins[pin.ix] = hash;
        return alloc_treenode(ctx, TAG_PIN(pin));
}

treenode_t jelly_pin(Jelly *ctx, hash256_t *pin) {
        printf("\tjelly_pin(hash=%lx)\n", pin->a);

        uint64_t hash = pin->a;
        hash += !hash;

        return insert_indirect_leaf(ctx,
            (IndirectInsertRequest){
                .width = PIN_TAGGED_WIDTH(32),
                .hash = hash,
                .leaf = (leaf_t) { .width_bytes = 32, .bytes = (uint8_t*) pin },
                .new_leaf = new_pin,
            }
        );
}


// Inserting Nats //////////////////////////////////////////////////////////////

/*
        Either returns a pointer to a matching treenode, or returns a
        pointer into an empty slot that may be written to.
*/
treenode_t jelly_cons(Jelly *ctx, treenode_t hed, treenode_t tel) {
        rehash_nodes_if_full(ctx);

        treenode_value val = TAG_PAIR(hed, tel);

        uint64_t target = val.word;
        uint64_t mask   = ctx->nodes_table_mask;

        uint32_t hash = (uint32_t) fmix64(val.word);
        uint64_t i = hash;

        for (;; i++) {
                i &= mask;

                FanEntry *cur = (ctx->nodes_table + i);

                uint64_t word = cur->val.word;

                if (word == target) {
                        printf("\t\tTREE MATCH @[%4d] = ", cur->pointer.ix);
                        print_tree(ctx, cur->pointer);
                        printf("\n");

                        return cur->pointer;
                }

                if (word == UINT64_MAX) {
                        treenode_t pointer = alloc_treenode(ctx, val);

                        cur->val     = val;
                        cur->hash    = (uint32_t) hash;
                        cur->pointer = pointer;

                        ctx->nodes_table_count++;

                        return pointer;
                }
        }
}


// Testing /////////////////////////////////////////////////////////////////////

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

treenode_t read_one(Jelly*);
treenode_t read_many(Jelly*, treenode_t);
treenode_t read_some(Jelly*);

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

treenode_t read_some(Jelly *ctx) {
        return read_many(ctx, read_one(ctx));
}

treenode_t read_many(Jelly *ctx, treenode_t acc) {
    loop:
        int c = getchar();

        switch (c) {
            case '"': {
                treenode_t elmt = jelly_bar(ctx, read_string());
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
                        treenode_t elmt = jelly_nat(ctx, read_hex());
                        acc = jelly_cons(ctx, acc, elmt);
                        goto loop;
                } else if (isdigit(c)) {
                        die("Non-zero numbers must not be zero-prefixed");
                } else {
                        ungetc(c, stdin);
                        uint32_t width = word64_bytes(0);
                        treenode_t elmt = jelly_packed_nat(ctx, width, 0);
                        acc = jelly_cons(ctx, acc, elmt);
                        goto loop;
                }
            }
            default:
                if isdigit(c) {
                        uint64_t word = read_word((uint64_t)(c - '0'));
                        uint32_t width = word64_bytes(word);
                        treenode_t elmt = jelly_packed_nat(ctx, width, word);
                        acc = jelly_cons(ctx, acc, elmt);
                        goto loop;
                } else {
                        die("Unexpected character: %c (read_many)\n", c);
                }
        }
}


treenode_t read_one(Jelly *ctx) {
        int c = getchar();
    loop:
        if (isdigit(c)) {
                uint64_t word = read_word((uint64_t)(c - '0'));
                uint32_t width = word64_bytes(word);
                return jelly_packed_nat(ctx, width, word);
        }

        switch (c) {
            case '"':
                return jelly_bar(ctx, read_string());
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

void jelly_dump(Jelly *ctx, uint8_t *buf) {
        strcpy((char*)buf, "\n\nGet fucked, lol\n\n");
}

void jelly_push_final(Jelly *ctx, treenode_t val) {
        push_fan_backref(ctx, val);
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
        uint8_t *buf = 0;
        if (l.width_bytes > 8) {
                buf = l.bytes;
        } else {
                buf = (uint8_t*) &(l.bytes);
        }

        printf("\"");
        for (int i = 0; i < l.width_bytes; i++) {
                uint8_t byte = buf[i];
                printf("%c", byte);
        }
        printf("\"");
}

size_t print_pin(Jelly *ctx, pin_t pin) {
        hash256_t h = *(ctx->pins[pin.ix]);
        char tmp[256];
        size_t output_size;
        b58enc(tmp, &output_size, &h, 32);
        printf("%s", tmp);
        return output_size;
}

void print_tree_outline(Jelly*, treenode_t);

void print_tree_outline_list(Jelly *ctx, treenode_t tree) {
        treenode_value val = ctx->treenodes[tree.ix];

        uint32_t offset = (uint32_t) val.word;

        switch (NODEVAL_TAG(val)) {
            case 4: printf("p%u", offset); break;
            case 5: printf("b%u", offset); break;
            case 6: printf("n%u", offset); break;
            case 7: printf("f%u", offset); break;
            default: {
                treenode_t hed = NODEVAL_HEAD(val);
                treenode_t tel = NODEVAL_TAIL(val);
                print_tree_outline_list(ctx, hed);
                putchar(' ');
                print_tree_outline(ctx, tel);
            }
        }
}

void print_tree_outline(Jelly *ctx, treenode_t tree) {
        treenode_value val = ctx->treenodes[tree.ix];

        uint32_t offset = (uint32_t) val.word;

        switch (NODEVAL_TAG(val)) {
            case 4: printf("p%u", offset); break;
            case 5: printf("b%u", offset); break;
            case 6: printf("n%u", offset); break;
            case 7: printf("f%u", offset); break;
            default: {
                treenode_t hed = NODEVAL_HEAD(val);
                treenode_t tel = NODEVAL_TAIL(val);
                putchar('(');
                print_tree_outline_list(ctx, hed);
                putchar(' ');
                print_tree_outline(ctx, tel);
                putchar(')');
            }
        }
}

void print_tree(Jelly*, treenode_t);

void print_tree_list(Jelly *ctx, treenode_t tree) {
        treenode_value val = ctx->treenodes[tree.ix];

        switch (NODEVAL_TAG(val)) {
            case 4: print_pin(ctx, NODEVAL_PIN(val)); break;
            case 5: print_bar(ctx, NODEVAL_BAR(val)); break;
            case 6: print_nat(ctx, NODEVAL_NAT(val)); break;
            case 7: printf("fan"); break;
            default: {
                treenode_t hed = NODEVAL_HEAD(val);
                treenode_t tel = NODEVAL_TAIL(val);
                print_tree_list(ctx, hed);
                putchar(' ');
                print_tree(ctx, tel);
            }
        }
}

void print_tree(Jelly *ctx, treenode_t tree) {
        treenode_value val = ctx->treenodes[tree.ix];

        switch (NODEVAL_TAG(val)) {
            case 4: print_pin(ctx, NODEVAL_PIN(val)); break;
            case 5: print_bar(ctx, NODEVAL_BAR(val)); break;
            case 6: print_nat(ctx, NODEVAL_NAT(val)); break;
            case 7: printf("fan"); break;
            default: {
                treenode_t hed = NODEVAL_HEAD(val);
                treenode_t tel = NODEVAL_TAIL(val);
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
        int hed_depth = treenode_depth(ctx, NODEVAL_HEAD(val));
        int tel_depth = treenode_depth(ctx, NODEVAL_TAIL(val));
        return (1 + ((hed_depth > tel_depth) ? hed_depth : tel_depth));
}

void jelly_debug_leaves(Jelly *ctx, bool details) {
        printf("\n\tleaves: (width=%u, count=%u, mask=%lu)\n\n",
               ctx->leaves_table_width,
               ctx->leaves_table_count,
               ctx->leaves_table_mask);

        int wid = ctx->leaves_table_width;

        uint64_t mask = ctx->leaves_table_mask;

        for (int i=0; i<wid; i++) {
                leaves_table_entry_t ent = ctx->leaves_table[i];

                // Empty slot
                if (ent.hash == 0) continue;

                printf("\t\t%4d = ", i);
                print_tree_outline(ctx, ent.pointer);
                printf("\t(val=");
                print_tree(ctx, ent.pointer);
                printf(")\n");

                if (details) {
                        uint64_t width = ((ent.width.word << 2) >> 2);
                        uint64_t bin = ent.hash & mask;
                        uint64_t distance = (((uint64_t)i) - bin);
                        printf("\t\t    bytes: %-4lu\n", width);
                        printf("\t\t    bin: %lu [dist=%lu]\n", bin, distance);
                        printf("\t\t    hash: 0x%016lx\n", ent.hash);
                }
        }
}


void jelly_debug_interior_nodes(Jelly *ctx) {
        printf("\n\tinterior_nodes_table: (width=%u, count=%u, mask=%lu)\n\n",
               ctx->nodes_table_width,
               ctx->nodes_table_count,
               ctx->nodes_table_mask);

        int wid = ctx->nodes_table_width;

        uint64_t mask = ctx->nodes_table_mask;

        for (int i=0; i<wid; i++) {
                FanEntry ent = ctx->nodes_table[i];

                // Empty slot
                if (ent.val.word == UINT64_MAX) continue;

                uint64_t bin = ent.hash & mask;
                uint64_t distance = (((uint64_t)i) - bin);

                printf("\t\t%4d = @%u\tbin=%lu\tdist=%lu\thash=%08x\n", i, ent.pointer.ix, bin, distance, ent.hash);

                // printf("\n\n\t\t  ");
                // print_tree(ctx, ent.pointer);
                // print_tree_outline(ctx, ent.pointer);
                // printf("\n\n");
        }
}


void jelly_debug(Jelly *ctx) {
        printf("\njelly_debug():\n");

        jelly_debug_leaves(ctx, false);

        jelly_debug_interior_nodes(ctx);

        uint32_t count = ctx->fans_count;

        {
                printf("\n\tpins: (width=%u, count=%u)\n\n",
                       ctx->pins_width,
                       ctx->pins_count);
                int num = ctx->pins_count;
                for (int i=0; i<num; i++) {
                        printf("\t\tp%d = ", i);
                        print_pin(ctx, (pin_t){ .ix = i });
                        printf("\n");
                }
        }


        {
                printf("\n\tbars: (width=%u, count=%u)\n\n",
                       ctx->bars_width,
                       ctx->bars_count);
                int num = ctx->bars_count;
                for (int i=0; i<num; i++) {
                        printf("\t\tb%d = ", i);
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
                                printf("\t\tt%d = p%u\n", i, (uint32_t) val.word);
                                continue;
                            case 5:
                                printf("\t\tt%d = b%u\n", i, (uint32_t) val.word);
                                continue;
                            case 6:
                                printf("\t\tt%d = n%u\n", i, (uint32_t) val.word);
                                continue;
                            case 7:
                                printf("\t\tt%d = f%u\n", i, (uint32_t) val.word);
                                continue;
                            default: {
                                treenode_t hed = NODEVAL_HEAD(val);
                                treenode_t tel = NODEVAL_TAIL(val);
                                printf("\t\tt%d = (%u, %u)\n", i, hed.ix, tel.ix);
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

        hash256_t dumb_hash_1 = { fmix64(111111), fmix64(65535), fmix64(9),  65536 };
        hash256_t dumb_hash_2 = { fmix64(222222), fmix64(33333), (0ULL - 1), (65535ULL << 12) };
        treenode_t top = read_some(ctx);
        treenode_t tmp = jelly_pin(ctx, &dumb_hash_1);
        top = jelly_cons(ctx, tmp, top);
        tmp = jelly_pin(ctx, &dumb_hash_1);
        top = jelly_cons(ctx, tmp, top);
        tmp = jelly_pin(ctx, &dumb_hash_2);
        top = jelly_cons(ctx, tmp, top);

        jelly_push_final(ctx, top);

        int wid = jelly_buffer_size(ctx);
        uint8_t *buf = calloc(wid, sizeof(uint8_t));

        jelly_debug(ctx);

        jelly_dump(ctx, buf);
        fwrite(buf, 1, wid, stdout);
        free_jelly_ctx(ctx);
        free(buf);

        return 0;
}
