#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "murmur3.h"


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

    // Deduplication table for direct atoms (atoms that fit in a 64 bit word)
    word_entry_t *words_table;
    uint32_t words_table_width;
    uint32_t words_table_count;

    // Array of duplicate tree-nodes (and the top-level node).  Each one
    // is an index into `treenodes`.
    treenode_t *fans;
    uint32_t fans_width;
    uint32_t fans_count;
} Jelly;

#define TAG_PIN(x) (x | (4ULL << 61)) // 100...
#define TAG_BAR(x) (x | (5ULL << 61)) // 101...
#define TAG_FAN(x) (x | (7ULL << 61)) // 111...

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

        // TODO: Deduplication table for pins

        // Array of unique bar-leaves
        res->bars       = malloc(8 * sizeof(leaf_t));
        res->bars_width = 8;
        res->bars_count = 0;

        // TODO: Deduplication table for bars

        // Array of unique nat-leaves
        res->nats       = malloc(8 * sizeof(leaf_t));
        res->nats_width = 8;
        res->nats_count = 0;

        // Deduplication table for direct-atoms.
        res->words_table = malloc(8 * sizeof(word_entry_t));
        res->words_table_width = 8;
        res->words_table_count = 0;

        // TODO: Deduplication for indirect-atoms.
        res->words_table = malloc(8 * sizeof(word_entry_t));
        res->words_table_width = 8;
        res->words_table_count = 0;

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
        free(ctx->words_table);
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


// Inserting ///////////////////////////////////////////////////////////////////

nat_t word_insert(Jelly *ctx, uint32_t bytes, uint64_t word) {
        printf("\tword_insert(wid=%u, %lu)\n", bytes, word);
        uint32_t end = ctx->words_table_count;

        /*
                Do a linear-search over the words table, and return the
                index of the corresponding nat, if we find one.
        */

        for (uint32_t i=0; i<end; i++) {
                if (ctx->words_table[i].value == word) {
                        return ctx->words_table[i].index;
                }
        }

        /*
                The word is unique, append the word to the end of the
                nats array.

                Because the input is a word, we can always store it in
                the `leaf_t` as a direct atom.
        */

        nat_t nat = alloc_nat(ctx);
        ctx->nats[nat.ix].width_bytes = bytes;
        leaf_t *l = (ctx->nats + nat.ix);
        uint8_t **hack_slot = &(l->bytes);
        *((uint64_t*)hack_slot) = word;

        /*
                If words_table isn't big enough to simply append the
                element, realloc it first.
        */

        uint32_t wid = ctx->words_table_width;
        if (end == wid) {
                wid = wid + (wid / 2);
                ctx->words_table = reallocarray(
                    ctx->words_table,
                    wid,
                    sizeof(ctx->words_table[0])
                );
                ctx->words_table_width = wid;
        }

        /*
                Append to the end of the words_table array.
        */

        ctx->words_table[end] = (word_entry_t){ .value = word, .index = nat };
        ctx->words_table_count = end+1;

        return nat;
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

workspace_t jelly_word(Jelly *ctx, uint64_t w) {
        printf("\tjelly_word(%lu)\n", w);
        uint64_t hash    = fmix64(w);
        uint32_t bytes   = word64_bytes(w);
        nat_t   nat      = word_insert(ctx, bytes, w);

        treenode_t pointer = alloc_treenode(ctx, TAG_NAT(nat));

        workspace_t res = alloc_workspace(ctx);
        FanEntry *ent = &(ctx->workspaces[res.ix]);

        ent->num_leaves  = 1;
        ent->num_bytes   = bytes;
        ent->leaves_hash = hash;
        ent->shape_hash  = 0;
        ent->pointer     = pointer;

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

workspace_t read_one(Jelly*);
workspace_t read_many(Jelly*, workspace_t);
workspace_t read_some(Jelly*);

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

workspace_t read_some(Jelly *ctx) {
        return read_many(ctx, read_one(ctx));
}

workspace_t read_many(Jelly *ctx, workspace_t acc) {
    loop:
        int c = getchar();

        if isdigit(c) {
                uint64_t word = read_word((uint64_t)(c - '0'));
                workspace_t elmt = jelly_word(ctx, word);
                acc = jelly_cons(ctx, acc, elmt);
                goto loop;
        }

        switch (c) {
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
            default:
                die("Unexpected character: %c (read_many)\n", c);
        }
}


workspace_t read_one(Jelly *ctx) {
        int c = getchar();
    loop:
        if (isdigit(c)) {
                uint64_t word = read_word((uint64_t)(c - '0'));
                return jelly_word(ctx, word);
        }

        switch (c) {
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

void print_nat(Jelly *ctx, nat_t nat) {
        leaf_t l = ctx->nats[nat.ix];
        if (l.width_bytes > 8) die("TODO");
        uint64_t value = (uint64_t) l.bytes;
        printf("%lu", value);
}

void print_tree(Jelly*, treenode_t);

void print_tree_list(Jelly *ctx, treenode_t tree) {
        treenode_value val = ctx->treenodes[tree.ix];

        switch (TREENODE_VAL_TAG(val)) {
            case 4: printf("pin"); break;
            case 5: printf("bar"); break;
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
            case 5: printf("bar"); break;
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

        uint32_t count = ctx->fans_count;

        {
                printf("\n\tnats: (width=%u, count=%u)\n\n",
                       ctx->nats_width,
                       ctx->nats_count);
                int num = ctx->nats_count;
                for (int i=0; i<num; i++) {
                        leaf_t l = ctx->nats[i];
                        if (l.width_bytes > 8) die("TODO");
                        uint64_t value = (uint64_t) l.bytes;
                        printf("\t\tn%d = %lu\n", i, value);
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
                            case 5:
                            case 6:
                            case 7:
                                printf("\t\t%d = f%u\n", i, (uint32_t) val.word);
                                continue;
                            default: {
                                treenode_t hed = TREENODE_VAL_HEAD(val);
                                treenode_t tel = TREENODE_VAL_TAIL(val);
                                printf("\t\t%d = (%u, %u)\n", i, hed.ix, tel.ix);
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
