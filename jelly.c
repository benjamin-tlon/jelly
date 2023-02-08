#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "murmur3.h"
#include "xxh3.h"
#include "libbase58.h"

#define DEBUG 0

#if DEBUG
#define debugf(...) printf(__VA_ARGS__)
#else
void pass() {}
#define debugf(...) pass(__VA_ARGS__)
#endif


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

#define die(...) (fprintf(stderr, __VA_ARGS__),exit(1))

// TODO: Testing.
FORCE_INLINE uint32_t word64_bits (uint64_t w) {
        return 64 - __builtin_clzll(w);
}

/*
        What is the byte-width of a 64-bit word?

        For example: 0 takes 0 bytes, 127 takes 1 byte, 256 takes two
        bytes, etc.

        TODO: This should return 0 for 0, but it seems to return 1.
*/
uint32_t word64_bytes (uint64_t w) {
        uint32_t bits = word64_bits(w);
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
typedef struct { uint32_t ix; } frag_t;       //  Index into Jelly.frags
typedef struct { uint32_t ix; } treenode_t;   //  Index into Jelly.treenodes

typedef struct {
        treenode_t head;
        treenode_t tail;
        uint32_t leaves;
} FragVal;

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
    uint32_t *refcounts;
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

    // Interior-Node Deduplication table.
    FanEntry *nodes_table;
    uint32_t nodes_table_width;
    uint32_t nodes_table_count;

    // Array of duplicate tree-nodes (and the top-level node).  Each one
    // is an index into `treenodes`.
    FragVal *frags;
    uint32_t frags_width;
    uint32_t frags_count;
} Jelly;

void print_tree_outline(Jelly*, treenode_t);
void print_tree(Jelly*, treenode_t);


FORCE_INLINE treenode_value TAG_FRAG(frag_t frag) {
        uint64_t index = (uint64_t) frag.ix;
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
FORCE_INLINE frag_t NODEVAL_FRAG(treenode_value v) {
        return (frag_t){ .ix = (uint32_t) v.word };
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
        res->refcounts       = calloc(32, sizeof(res->refcounts[0]));
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
        size_t leaves_table_bytes = (1<<2) * sizeof(res->leaves_table[0]);
        res->leaves_table       = malloc(leaves_table_bytes);
        res->leaves_table_width = (1<<2);
        res->leaves_table_count = 0;
        memset(res->leaves_table, 255, leaves_table_bytes);

        // Deduplication table for interior nodes
        size_t nodes_table_bytes = (1<<2) * sizeof(res->nodes_table[0]);
        res->nodes_table       = malloc(nodes_table_bytes);
        res->nodes_table_width = (1<<2);
        res->nodes_table_count = 0;
        memset(res->nodes_table, 255, nodes_table_bytes);

        // Array of duplicate tree-nodes (and the top-level node).  Each one
        // is an index into `treenodes`.
        res->frags       = calloc(16, sizeof(res->frags[0]));
        res->frags_width = 16;
        res->frags_count = 0;

        return res;
}

/*
        We don't free ny memory or shrink any tables, we just set the
        size-counts of everything to 0 and empty all the hashtable slots.
*/
void wipe_jelly_ctx (Jelly *ctx) {
        ctx->treenodes_count = 0;
        ctx->pins_count = 0;
        ctx->bars_count = 0;
        ctx->nats_count = 0;
        ctx->frags_count = 0;

        ctx->leaves_table_count = 0;
        ctx->nodes_table_count = 0;

        memset(ctx->leaves_table, 255,
               (sizeof(ctx->leaves_table[0]) * ctx->leaves_table_width));

        memset(ctx->nodes_table, 255,
               (sizeof(ctx->nodes_table[0]) * ctx->nodes_table_width));
}

void free_jelly_ctx (Jelly *ctx) {
        free(ctx->treenodes);
        free(ctx->refcounts);
        free(ctx->pins);
        free(ctx->bars);
        free(ctx->nats);
        free(ctx->leaves_table);
        free(ctx->nodes_table);
        free(ctx->frags);
        free(ctx);
}


FORCE_INLINE treenode_t alloc_treenode(Jelly *c, treenode_value v) {
        uint32_t res = c->treenodes_count++;
        uint32_t wid = c->treenodes_width;

        if (res >= wid) {
                wid *= 2;
                c->treenodes = reallocarray(c->treenodes, wid, sizeof(c->treenodes[0]));
                c->refcounts = reallocarray(c->refcounts, wid, sizeof(c->refcounts[0]));
                c->treenodes_width = wid;
        }

        c->treenodes[res] = v;
        c->refcounts[res] = 0;

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

        debugf("\t\tREHASH_NODES (old=%u new=%u)\n", oldwid, newwid);

        // debugf("OLD TABLE:\n");
        // jelly_debug_interior_nodes(ctx);
        // debugf("\n");

        FanEntry *oldtab = ctx->nodes_table;
        FanEntry *newtab = malloc(newwid * sizeof(*newtab));

        memset(newtab, 255, newwid * sizeof(*newtab));

        uint64_t newmask = newwid - 1;

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
                                // debugf("\t\t%lu -> %lu\n", i, j);
                                break;
                        } // else {
                                // debugf("\t\t\t(%lu -> %lu) is taken\n", i, j);
                        // }
                }
        }

        free(ctx->nodes_table);

        ctx->nodes_table_width = newwid;
        ctx->nodes_table       = newtab;

        // debugf("NEW TABLE:\n");
        // jelly_debug_interior_nodes(ctx);
        // debugf("\n");
}

void rehash_leaves_if_full(Jelly *ctx) {
        uint32_t oldwid = ctx->leaves_table_width;
        uint32_t num    = ctx->leaves_table_count;

        // If capacity is >= 50%, resize.
        if (num*2 < oldwid) return;

        uint32_t newwid = oldwid*2;

        debugf("\t\tREHASH_LEAVES (old=%u new=%u)\n", oldwid, newwid);

        // debugf("OLD TABLE:\n");
        // jelly_debug_leaves(ctx, false);
        // debugf("\n");

        leaves_table_entry_t *oldtab = ctx->leaves_table;
        leaves_table_entry_t *newtab = malloc(newwid * sizeof(*newtab));

        memset(newtab, 255, newwid * sizeof(*newtab));

        uint64_t newmask = newwid - 1;

        /*
                Loop over the whole of oldtab, and re-insert every
                non-empty slot.  Inserts are guarenteed to be unique,
                so we just need to, starting at the correct bucket,
                scan for an empty slot and write the value there.
        */
        for (uint64_t i=0; i<oldwid; i++) {
                leaves_table_entry_t ent = oldtab[i];

                // empty slot
                if (ent.pointer.ix == UINT32_MAX) continue;

                uint64_t j = ent.hash;
                for (;; j++) {
                        j &= newmask;
                        leaves_table_entry_t *tar = newtab + j;
                        if (tar->pointer.ix == UINT32_MAX) {
                                *tar = ent;
                                debugf("\t\t%lu -> %lu\n", i, j);
                                break;
                        } else {
                                debugf("\t\t\t(%lu -> %lu) is taken\n", i, j);
                        }
                }
        }

        free(ctx->leaves_table);

        ctx->leaves_table_width = newwid;
        ctx->leaves_table       = newtab;

        // debugf("NEW TABLE:\n");
        // jelly_debug_leaves(ctx, false);
        // debugf("\n");
}


treenode_t insert_packed_leaf(Jelly *ctx, PackedInsertRequest req) {
        uint64_t tag   = (req.width.word >> 62);
        uint64_t bytes = (req.width.word << 2) >> 2;

        debugf(
            "\tpacked_insert(tag=%lu, wid=%lu, hash=%lu, %lu)\n",
            tag,
            bytes,
            req.hash,
            req.word
        );

        /*
                Do a linear-search over the leaves table, and use the
                index of the corresponding nat, if we find one.

                We make sure ahead-of-time that there is enough space
                for us to insert.  This way we don't have to move things
                around while were are inserting.
        */

        rehash_leaves_if_full(ctx);

        leaves_table_entry_t *ent = NULL;

        uint64_t mask = ctx->leaves_table_width - 1;
        uint64_t ix = req.hash;

        for (;; ix++) {
                ix &= mask;

                ent = &(ctx->leaves_table[ix]);

                if (ent->pointer.ix == UINT32_MAX) break;

                bool match = ent->hash == req.hash &&
                             ent->width.word == req.width.word &&
                             ent->bytes == (uint8_t*) req.word;

                if (!match) continue;

                debugf("\t\tPACKED_LEAF_MATCH:\n\t\t\t");
                ctx->refcounts[ent->pointer.ix]++;
                if (DEBUG) print_tree_outline(ctx, ent->pointer);
                debugf(" = ");
                if (DEBUG) print_tree(ctx, ent->pointer);
                debugf("\n");

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

        debugf("\t\tNEW_PACKED_LEAF:\n\t\t\t");
        if (DEBUG) print_tree_outline(ctx, pointer);
        debugf(" = ");
        if (DEBUG) print_tree(ctx, pointer);
        debugf("\n");

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
        debugf("\tinsert_indirect_leaf(wid=%u)\n", req.leaf.width_bytes);

        /*
                Do a linear-search over the leaves table, and return the
                index of the corresponding nat, if we find one.
        */

        rehash_leaves_if_full(ctx); // Make sure there is space to insert,
                                    // if we do it later, we will invalidate
                                    // our pointer.

        leaves_table_entry_t *ent;

        uint64_t mask = ctx->leaves_table_width - 1;
        uint64_t ix = req.hash;

        int32_t byt_wid = req.leaf.width_bytes;
        uint8_t *byt    = req.leaf.bytes;

        for (;; ix++) {
                ix &= mask;

                ent = &(ctx->leaves_table[ix]);

                // An empty slot ends the search.
                if (ent->pointer.ix == UINT32_MAX) break;

                if (ent->hash != req.hash) continue;

                if (ent->width.word != req.width.word) continue;

                if (0 != memcmp(ent->bytes, byt, byt_wid)) continue;

                debugf("\t\tINDIRECT_LEAF_MATCH:\n\t\t\t");
                ctx->refcounts[ent->pointer.ix]++;
                if (DEBUG) print_tree_outline(ctx, ent->pointer);
                debugf(" = ");
                if (DEBUG) print_tree(ctx, ent->pointer);
                debugf("\n");

                return ent->pointer;
        }

        ctx->leaves_table_count++;

        treenode_t pointer = req.new_leaf(ctx, req.leaf);

        debugf("\t\tNEW_INDIRECT_LEAF:\n\t\t\t");
        if (DEBUG) print_tree_outline(ctx, pointer);
        debugf(" = ");
        if (DEBUG) print_tree(ctx, pointer);
        debugf("\n");

        *ent = (leaves_table_entry_t){
            .hash  = req.hash,
            .width = req.width,
            .bytes = req.leaf.bytes,
            .pointer = pointer,
        };

        return pointer;
}

frag_t alloc_frag(Jelly *ctx, FragVal frag) {
        uint32_t nex = ctx->frags_count++;
        uint32_t wid = ctx->frags_width;

        if (nex >= wid) {
                wid *= 2;
                ctx->frags = reallocarray(ctx->frags, wid, sizeof(ctx->frags[0]));
                ctx->frags_width = wid;
        }

        ctx->frags[nex] = frag;

        return (frag_t){ .ix = nex };
}


/*
        When we encounter our second reference to an interned interior
        node, we shatter it.

        We mutate it's value to be a reference into the fragments table,
        and then we insert
*/

void frag(Jelly *ctx, treenode_t tree, uint32_t leaf_count) {
        treenode_value val = ctx->treenodes[tree.ix];

        // Ignore leaves and already-fragmented values.
        switch (NODEVAL_TAG(val)) {
            case 4:
            case 5:
            case 6:
            case 7:
                return;
        }

        treenode_t hed = NODEVAL_HEAD(val);
        treenode_t tel = NODEVAL_TAIL(val);

        frag_t frag = alloc_frag(ctx, (FragVal){ .head=hed,
                                                 .tail=tel,
                                                 .leaves=leaf_count,
                                               });

        debugf("%u", leaf_count);

        ctx->treenodes[tree.ix] = TAG_FRAG(frag);
}

struct shatter {
        uint32_t refs;
        uint32_t leaves;
};

struct shatter shatter(Jelly *ctx, treenode_t tree, bool top) {
        treenode_value val = ctx->treenodes[tree.ix];

        uint32_t refs = ctx->refcounts[tree.ix];

        // Leaf or existing fragment.
        switch (NODEVAL_TAG(val)) {
            case 4:
            case 5:
            case 6:
            case 7:
                return (struct shatter){ .refs = refs, .leaves = 1 };
        }

        treenode_t head = NODEVAL_HEAD(val);
        treenode_t tail = NODEVAL_TAIL(val);

        struct shatter hed = shatter(ctx, head, false);
        struct shatter tel = shatter(ctx, tail, false);

        if (hed.refs > refs) { frag(ctx, head, hed.leaves); hed.leaves = 1; }
        if (tel.refs > refs) { frag(ctx, tail, tel.leaves); tel.leaves = 1; }

        uint32_t leaves = hed.leaves + tel.leaves;

        if (top) frag(ctx, tree, leaves);

        return (struct shatter){ .refs = refs, .leaves = leaves };
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
        debugf("\tjelly_packed_nat(%lu, width=%u)\n", word, byte_width);

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
        debugf("\tjelly_nat(width=%d)\n", leaf.width_bytes);

        if (leaf.width_bytes < 9) {
                return jelly_packed_nat(ctx, leaf.width_bytes, (uint64_t)leaf.bytes);
        } else {
                uint64_t hash = XXH3_64bits(leaf.bytes, leaf.width_bytes);

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
        debugf("\tjelly_bar(width=%d)\n", leaf.width_bytes);

        if (leaf.width_bytes < 9) {
                uint64_t word = (uint64_t) leaf.bytes;

                debugf("\tjelly_packed_bar(%lu, width=%u)\n", word, leaf.width_bytes);

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
        debugf("\tjelly_pin(hash=%lx)\n", pin->a);

        uint64_t hash = pin->a;

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
        uint64_t mask   = ctx->nodes_table_width - 1;

        uint32_t hash = (uint32_t) fmix64(val.word);
        uint64_t i = hash;

        for (;; i++) {
                i &= mask;

                FanEntry *cur = (ctx->nodes_table + i);

                uint64_t word = cur->val.word;

                if (word == target) {
                        debugf("\t\tTREE MATCH\n\t\t\tt%d = ", cur->pointer.ix);
                        if (DEBUG) print_tree(ctx, cur->pointer);
                        debugf("\n");

                        ctx->refcounts[cur->pointer.ix]++;

                        return cur->pointer;
                }

                if (word == UINT64_MAX) {
                        treenode_t pointer = alloc_treenode(ctx, val);

                        cur->val     = val;
                        cur->hash    = hash;
                        cur->pointer = pointer;

                        ctx->nodes_table_count++;

                        return pointer;
                }
        }
}


// Serializing /////////////////////////////////////////////////////////////////

// 127        -> 0b01111111
// 128        -> 0b10000001_10000000
// 255        -> 0b10000001_11111111
// 256        -> 0b10000001_00000000_00000001
// UINT64_MAX -> 0b10001000_11111111_(x8)
uint64_t word_dumpsize(uint64_t word) {
        if (word < 128) return 1;
        return 1 + word64_bytes(word);
}

// TODO: Test test test
void word_dump(uint8_t **ptr, uint64_t word) {
        uint8_t *buf = *ptr;
        if (word < 128) {
                *buf++ = (uint8_t) word;
        } else {
                uint8_t width = (uint8_t) word64_bytes(word);

                uint8_t tagged_width = width | (2<<6); // 0b10xxxxxx

                debugf("\tword_dump: width = %u\n", width);
                debugf("\tword_dump: tagged_width = %u\n", tagged_width);

                *buf++ = tagged_width;

                // TODO Can I use memcpy here?
                for (int i=0; i<width; i++) {
                        uint8_t byte = word;
                        debugf("\tword_dump: byte[%d] = %u (rest=%lu)\n", i, byte, word);
                        *buf++ = byte;
                        word = word >> 8;
                }
        }

        *ptr = buf;
}

// width=63  -> 0b10111111_xxxxxxxx(x63)
// width=64  -> 0b11000001_01000000_xxxxxxxx*63
// width=64  -> 0b11000001_01000000_xxxxxxxx*64
// width=255 -> 0b11000001_11111111_xxxxxxxx*255
// width=256 -> 0b11000010_00000000_000000001_xxxxxxxx(*256)
uint64_t leaf_dumpsize(leaf_t leaf) {
        if (leaf.width_bytes < 9) {
                // debugf("\tDIRECT (word=%lu)\n", (uint64_t) leaf.bytes);
                return word_dumpsize((uint64_t) leaf.bytes);
        }

        // single-byte length field
        if (leaf.width_bytes < 64) {
                // debugf("\tSHORT\n");
                return 1 + leaf.width_bytes;
        }

        // variable-byte length field
        // debugf("\tLONG %d %u\n", leaf.width_bytes, word64_bytes(leaf.width_bytes));
        return 1 + word64_bytes(leaf.width_bytes) + leaf.width_bytes;
}

// TODO: Test test test
uint8_t *dump_leaf(uint8_t *buf, leaf_t leaf) {
        uint64_t wid = leaf.width_bytes;

        if (wid < 9) {
                debugf("\tdump direct\n");
                word_dump(&buf, (uint64_t) leaf.bytes);
                return buf;
        }

        if (wid < 64) {
                debugf("\tdump small [wid=%lu]\n", wid);
                uint8_t wid_byte = wid;
                uint8_t wid_tagged = (wid_byte | 128); // 0b10xxxxxx
                debugf("\t\twidth_tagged = %u\n", wid_tagged);
                *buf++ = wid_tagged;
                memcpy(buf, leaf.bytes, wid);
                buf += wid;
                return buf;
        }

        {
                debugf("\tdump big\n");
                uint8_t widwid = word64_bytes(wid);

                *buf++ = (widwid | 192); // 0b11xxxxxx

                // TODO Can I use memcpy here too?
                uint64_t tmp = wid;
                for (int i=0; i<widwid; i++) {
                        *buf++ = (uint8_t) tmp;
                        tmp = tmp >> 8;
                }

                memcpy(buf, leaf.bytes, wid);
                buf += wid;
                return buf;
        }
}

void print_bar(Jelly*, bar_t);
void print_nat(Jelly*, nat_t);
void print_fragment_outline(Jelly*, frag_t);

struct ser {
        char *buf;
        size_t wid;
};

void showbits(char *key, int wid, int num, uint64_t bits) {
        if (!DEBUG) return;
        if (key) debugf(" %s:", key);
        else debugf(" (");

        int extra = wid - num;
        for (int i=0; i<extra; i++) putchar('_');

        for (num--; num >= 0; num--) {
                putchar(((bits>>num)&1) ? '1' : '0');
        }
        if (!key) putchar(')');
}

struct frag_state {
        uint8_t acc;
        uint8_t fil;
        treenode_t *stack;
        uint8_t offs[4];
        uint64_t refbits;
        uint8_t *out;
};

void serialize_frag(Jelly *ctx, struct frag_state *st, FragVal frag) {
        treenode_t *stack = st->stack;
        uint8_t fil       = st->fil;
        uint8_t acc       = st->acc;
        uint8_t *offs     = st->offs;
        uint64_t refbits  = st->refbits;
        uint8_t *out      = st->out;

        int sp = 0;

        stack[sp++] = frag.tail;
        stack[sp]   = frag.head;

        int parens = 1;

        while (sp >= 0) {

                treenode_t t = stack[sp];

                treenode_value val = ctx->treenodes[t.ix];

                if (val.word >> 63) {
                        debugf("\n");
                        for (int x=0; x<sp; x++) debugf(" ");
                        if (DEBUG) print_tree_outline(ctx, t);
                        //debugf("\n");

                        debugf("\n\t\t\t\t");
                        showbits(NULL, 8, fil, acc);

                        // Output a zero bit
                        fil = (fil + 1) % 8;
                        if (!fil) { showbits("flush", 8, 8, acc); *out++ = acc; acc = 0; }

                        showbits("tag", 8, fil, acc);

                        uint8_t ptr = (uint8_t) val.word;
                        uint64_t tag = ((val.word << 1) >> 62);
                        uint32_t bits = ptr + offs[tag];

                        uint8_t new_bits = (bits << fil);
                        uint8_t overflow = (bits >> (8 - fil));

                        // debugf(" [extra=%u]", extra);
                        showbits("v", refbits, refbits, bits);

                        // int new_count = ((fil+refbits)>8) ? 8 : (fil + refbits);
                        // int ovo_count = ((fil+refbits)>8) ? (8 - (fil + refbits)) : 0;
                        // showbits("mask", 8, new_count, new_bits);

                        acc |= new_bits;
                        fil += refbits;

                        if (fil >= 8) {
                                showbits("overflow", (fil-8), (fil-8), overflow);
                                showbits("flush", 8, 8, acc);
                                *out++ = acc;
                                acc = overflow;
                                fil -= 8;
                        }

                        showbits(NULL, 8, fil, acc);

                        debugf("\n");
                        sp--;
                } else {
                        debugf("\n");
                        for (int x=0; x<sp; x++) debugf(" ");
                        debugf("-");
                        parens++;

                        debugf("\n\t\t\t\t");
                        showbits(NULL, 8, fil, acc);

                        // Output a 1 bit.
                        acc |= (1 << fil);
                        showbits("tag", 8, (fil+1), acc);
                        fil = (fil + 1) % 8;
                        if (!fil) { showbits("flush", 8, 8, acc); *out++ = acc; acc = 0; }

                        showbits(NULL, 8, fil, acc);

                        // Replaces the current stack pointer.
                        stack[sp++] = NODEVAL_TAIL(val);
                        stack[sp]   = NODEVAL_HEAD(val);
                }
        }

        st->fil = fil;
        st->acc = acc;
        st->out = out;
}

struct ser serialize(Jelly *ctx) {
        uint64_t width = 0;
        uint64_t refrs = 0;

        uint32_t numpins  = ctx->pins_count;
        uint32_t numbars  = ctx->bars_count;
        uint32_t numnats  = ctx->nats_count;
        uint32_t numfrags = ctx->frags_count;

        refrs += numpins;
        width += word_dumpsize(numpins);
        debugf("buffer_bytes (after pin_width) = %lu\n", width);
        width += numpins * sizeof(hash256_t);
        debugf("buffer_bytes (after pins) = %lu\n", width);

        refrs += numbars;
        width += word_dumpsize(numbars);
        debugf("buffer_bytes (bar_width) = %lu\n", width);
        for (int i=0; i<numbars; i++) {
                width += leaf_dumpsize(ctx->bars[i]);
                debugf("buffer_bytes (after b%d) = %lu\n", i, width);
                debugf("\n\tb%d = ", i);
                if (DEBUG) print_bar(ctx, (bar_t){ .ix = i });
                debugf("\n\n");
        }

        refrs += numnats;
        width += word_dumpsize(numnats);
        debugf("buffer_bytes (after nats count) = %lu\n", width);
        for (int i=0; i<numnats; i++) {
                width += leaf_dumpsize(ctx->nats[i]);
                debugf("buffer_bytes (after n%d) = %lu\n", i, width);
                debugf("\n\tn%d = ", i);
                if (DEBUG) print_nat(ctx, (nat_t){ .ix = i });
                debugf("\n\n");
        }

        width += word_dumpsize(ctx->frags_count);
        debugf("buffer_bytes (after frags count) = %lu\n", width);

        uint64_t treebits = 0;

        for (int i=0; i<numfrags;  i++) {
                uint32_t maxref     = refrs - 1;
                uint32_t leaf_width = word64_bits(maxref);
                uint32_t leaves     = ctx->frags[i].leaves;
                uint32_t frag_bits  = (leaves*leaf_width) + (leaves * 2) - 2;
                    // (... - 2) because we can omit the outermost
                    // one bit.  Every fragment is a pair, so the
                    // parser just reads two forms per fragment:
                    // (head, tail).

                treebits += frag_bits;
                refrs++;

                debugf("tree_bits (frags) = %lu\n", treebits);
                debugf("\n\t[maxref=%u leafwid=%u leaves=%u bits=%u]\n", maxref, leaf_width, leaves, frag_bits);
                debugf("\n\t\tf%d = ", i);
                if (DEBUG) print_fragment_outline(ctx, (frag_t){i});
                debugf("\n\n");
        }

        // Tree-bits is padded to be a multiple of 8 (treat as a byte-array);
        uint64_t hanging_bits = treebits % 8;
        if (hanging_bits) treebits += (8 - hanging_bits);

        // Add in the tree bits;
        width += (treebits/8);

        debugf("total_byte_width = %lu\n", width);

        size_t total_byte_width = width;

        // Result is always a multiple of 8 (so we can treat it as an
        // array of 64-bit words);
        uint64_t hanging_bytes = width % 8;
        if (hanging_bytes) width += (8 - hanging_bytes);

        debugf("padded byte width = %lu\n", width);

        struct ser result = (struct ser) { .buf = calloc(1, width), .wid = total_byte_width };

        // Dumping

        uint8_t *out = (uint8_t*) result.buf;

        // Dumping Leaves

        word_dump(&out, numpins);
        for (int i=0; i<numpins; i++) {
                debugf("p%d\n", i);
                hash256_t *pin = ctx->pins[i];
                memcpy(out, pin, 32);
                out += 32;
        }

        word_dump(&out, numbars);
        for (int i=0; i<numbars; i++) {
                debugf("b%d\n", i);
                out = dump_leaf(out, ctx->bars[i]);
        }

        word_dump(&out, numnats);
        for (int i=0; i<numnats; i++) {
                debugf("n%d\n", i);
                out = dump_leaf(out, ctx->nats[i]);
        }

        word_dump(&out, numfrags);

        uint64_t maxdepth = 128; // TODO

        uint8_t acc = 0;
        uint8_t fil = 0;
        treenode_t stack[maxdepth];

        struct frag_state st;
        st.acc = 0;
        st.fil = 0;
        st.stack = stack;
        st.out = out;

        for (int i=0; i<numfrags; i++) {
                uint8_t pinoff   = 0;
                uint8_t baroff   = pinoff + numpins;
                uint8_t natoff   = baroff + numbars;
                uint8_t frgoff   = natoff + numnats;
                uint64_t maxref  = frgoff + i - 1;
                uint64_t refbits = word64_bits(maxref);

                debugf("\nFRAG(%d) [maxref=%lu refbits=%lu]\n\n", i, maxref, refbits);

                st.offs[0] = pinoff;
                st.offs[1] = baroff;
                st.offs[2] = natoff;
                st.offs[3] = frgoff;
                st.refbits = refbits;

                serialize_frag(ctx, &st, ctx->frags[i]);
        }

        out = st.out;
        acc = st.acc;
        fil = st.fil;

        if (fil > 0) {
                debugf("\n");
                showbits("final_flush", 8, 8, acc);
                debugf("\n\n");
                *out++ = acc;
                acc = 0;
                fil = 0;
        }

        // Dumping Tree Fragments (TODO)

        return result;
}


FORCE_INLINE uint8_t load_byte(struct ser *st) {
    debugf("\tload_byte() [remain=%lu]\n", st->wid);

    // debugf("\t    [");
    // for (int i=0; i<10; i++) debugf("%u%s", (uint8_t) st->buf[i], (i==9 ? "" : " "));
    // debugf("]\n");

    if (st->wid == 0) die("load_byte(): EOF");
    uint8_t byte = *(st->buf++);
    st->wid--;
    debugf("\t  -> %u [remain=%lu]\n", byte, st->wid);
    return byte;
}

uint64_t load_word(struct ser *st) {
    uint8_t byt = load_byte(st);

    if (byt < 128) {
            return (uint64_t) byt;
    }

    uint8_t flag = byt & (1<<6);
    uint8_t len = byt & ((1<<6)-1);

    if (flag || len > 8) {
            die("load_word(): length-byte is too big");
    }

    if (st->wid < len) {
            die("load_word(): EOF after length byte");
    }

    uint64_t result = 0;

    for (int i=0; i<len; i++) {
            uint64_t new = *(st->buf++);
            result |= (new << (8*i));
    }

    st->wid -= len;

    return result;
}

FORCE_INLINE leaf_t load_leaf(struct ser *st) {
        debugf("load_leaf:\n");

        leaf_t result;

        uint8_t byt = load_byte(st);

        if (byt < 128) {
                debugf("\tSingle-byte leaf: %u\n", byt);
                uint64_t word = byt;
                return (leaf_t) {
                        .width_bytes = (word == 0 ? 0 : 1),
                        .bytes       = (uint8_t*) word,
                };
        }

        uint8_t flag = byt & (1<<6);
        uint8_t len = byt & ((1<<6)-1);

        debugf("\tflag=%u len=%u\n", flag, len);

        if (st->wid < len) {
                die("load_word(): EOF after length byte");
        }

        // Packed Leaf
        if (!flag && len < 9) {

                debugf("\tpacked\n");

                uint64_t word = 0;

                for (int i=0; i<len; i++) {
                        uint8_t new_byte = *(st->buf++);
                        uint64_t new = new_byte;
                        debugf("\t\tword = %lu, new= %lu\n", word, new);
                        word |= (new << (8*i));
                }

                debugf("\t\tword = %lu\n", word);

                st->wid -= len;

                result.bytes = (uint8_t*) word;
                result.width_bytes = len;
                return result;
        }

        // Short indirect leaf, byte indicates length.
        if (!flag) {
                debugf("\tshort indirect\n");
                result.bytes = (uint8_t*) st->buf;
                result.width_bytes = len;

                st->wid -= len;
                st->buf += len;

                return result;
        }

        debugf("\tlong indirect\n");

        if (len > 8) die("impossibly big length-of-length\n");

        uint64_t actual_length = 0;

        for (int i=0; i<len; i++) {
                uint64_t new = *(st->buf++);
                actual_length |= (new << (8*i));
        }
        st->wid -= len;

        if (st->wid < actual_length) {
                die("load_word(): EOF after variable-width length\n");
        }

        result.bytes = (uint8_t*) st->buf;
        result.width_bytes = actual_length;

        st->wid -= actual_length;
        st->buf += actual_length;

        return result;
}

/*
        Alright, what's the algorithm here?

        To read a bit:

            mask the against `1<<red`.
            bit = (acc & (1<<red));
            red = (red+1) & 8;
            if (!red) { WIDTH_CHECK; acc = *buf++; }

        To read a tree:

            -   Read a bit.

            -   If the bit is zero, this is a leaf
                -   Read n bits.
                -   n is the bit-width of the maximum backref at this point.
                -   Use the resulting word, to construct a TreeVal.

                    Can we do better than the following?

                        CHECK_IF_WITHIN_MAXREF()
                        if (bits < num_pins) return TAG_PIN(bits);
                        bits -= num_pins
                        if (bits < num_bars) return TAG_BAR(bits);
                        bits -= num_bars;
                        if (bits < num_nats) return TAG_NAT(bits);
                        bits -= num_nats;
                        return TAG_FRAG(bits);

            -   Otherwise, this is a cell.

                Here, just use naive recursion, this is not
                performance-sensitive code.

                Just call load_frag() directly, since that gets a pair.

                    (Uninline it and make it into a function againn);

            -   That's it!  This is not terribly complicated.
*/

struct frag_loader_state {
        Jelly *ctx;
        uint8_t *buf; // Pointer into the remaining bytes.
        uint64_t wid; // Remaining bytes in the buffer.
        uint8_t acc;  // The last byte read from the buffer.
        uint64_t red; // The number of bits of `acc` that have been consumed.

        uint64_t ref_bits; // The number of bits per leaf
        uint64_t max_ref;  // The maximum ref that we can accept;
        uint64_t num_pins;
        uint64_t num_bars;
        uint64_t num_nats;
};

struct load_fragtree_result {
        treenode_t tree;
        uint32_t leaves;
};

struct load_fragtree_result load_fragtree(struct frag_loader_state *s);

FragVal load_fragment(struct frag_loader_state *s) {
        struct load_fragtree_result hed = load_fragtree(s);
        struct load_fragtree_result tel = load_fragtree(s);

        uint32_t leaves = hed.leaves + tel.leaves;

        FragVal fv = { .head=hed.tree, .tail=tel.tree, .leaves=leaves };
        return fv;
}

treenode_value decode_leaf(struct frag_loader_state *s, uint32_t leaf) {
        if (leaf < s->num_pins) return TAG_PIN((pin_t){leaf});

        leaf -= s->num_pins;
        if (leaf < s->num_bars) {
                return TAG_BAR((bar_t){leaf});
        }

        leaf -= s->num_bars;
        if (leaf < s->num_nats) {
                return TAG_NAT((nat_t){leaf});
        }

        leaf -= s->num_nats;
        return TAG_FRAG((frag_t){leaf});
}

struct load_fragtree_result
load_fragtree(struct frag_loader_state *s) {
        int refbits = s->ref_bits;

        uint8_t bit = (s->acc & (1 << s->red));

        s->red = (s->red + 1) % 8;
        if (!s->red) {
                // WIDTH_CHECK
                s->acc = *(s->buf)++;
        }

        if (bit) {
                debugf("cell\n");
                FragVal res = load_fragment(s);
                treenode_t tr = alloc_treenode(s->ctx, TAG_PAIR(res.head, res.tail));
                return (struct load_fragtree_result){ .tree=tr, .leaves=res.leaves };
        }

        debugf("leaf\n");
        debugf("refbits=%u\n", refbits);

        //      -   Read n bits.
        //      -   n is the bit-width of the maximum backref at this point.
        uint32_t leaf_mask = (1 << refbits) - 1;

        uint32_t leaf = (s->acc >> s->red) & leaf_mask;

        debugf("acc=%u red=%lu | mask=%u leaf=%u\n", s->acc, s->red, leaf_mask, leaf);

        int oldred = s->red;
        debugf("[[refbits=%d oldred=%d]]\n", refbits, oldred);

        s->red += refbits;

        if (s->red >= 8) {
                int extra   = (oldred + refbits) - 8;
                int remain  = 8-extra;
                int already = refbits - extra;

                uint8_t nex = s->buf[0];

                uint8_t why = nex & ((1<<extra) - 1);
                uint8_t more = why << already;

                debugf("[[nex=%u extra=%u remain=%u already=%u more=%u why=%u]]\n", nex, extra, remain, already, more, why);

                leaf |= more;

                s->red -= 8;
                s->acc = nex;
                s->buf++;
        }

        if (leaf > s->max_ref) {
                die("leaf val is out-of-bounds\n");
        }

        treenode_value v = decode_leaf(s, leaf);

        treenode_t t = alloc_treenode(s->ctx, v);

        return (struct load_fragtree_result){ .tree = t, .leaves = 0 };

/*
                -   Use the resulting word, to construct a TreeVal.

                    Can we do better than the following?

                        if (bits < num_pins) return TAG_PIN(bits);
                        bits -= num_pins
                        if (bits < num_bars) return TAG_BAR(bits);
                        bits -= num_bars;
                        if (bits < num_nats) return TAG_NAT(bits);
                        bits -= num_nats;
                        return TAG_FRAG(bits);

*/
}

/*
        Note that for pins, indirect atoms, and indirect bars we do not
        copy, we just slice the input buffer.
*/
void deserialize(Jelly *ctx, struct ser st) {
        debugf("\n");

        uint64_t num_pins = load_word(&st);

        debugf("num_pins = %lu\n", num_pins);

        for (int i=0; i<num_pins; i++) {
                debugf("\tload_pin() [remain=%lu]\n", st.wid);
                if (st.wid < 32) die("load_pin(): EOF");
                ctx->pins[alloc_pin(ctx).ix] = (hash256_t*) st.buf;
                st.buf += 32;
                st.wid -= 32;
        }

        uint64_t num_bars = load_word(&st);
        debugf("num_bars = %lu\n", num_bars);
        for (int i=0; i<num_bars; i++) {
                ctx->bars[alloc_bar(ctx).ix] = load_leaf(&st);
        }

        uint64_t num_nats = load_word(&st);

        debugf("num_nats = %lu\n", num_nats);
        for (int i=0; i<num_nats; i++) {
                ctx->nats[alloc_nat(ctx).ix] = load_leaf(&st);
        }

        uint64_t num_frags = load_word(&st);
        debugf("num_nats = %lu\n", num_nats);

        if (num_frags) {
                uint8_t acc = *st.buf++;

                struct frag_loader_state s = {
                        .ctx = ctx,
                        .buf = (uint8_t*) st.buf,
                        .acc = acc,
                        .red = 0,
                        .ref_bits = 0,
                        .max_ref = 0,
                        .num_pins = num_pins,
                        .num_bars = num_bars,
                        .num_nats = num_nats,
                };

                for (int i=0; i<num_frags; i++) {
                        int num_refs = num_pins + num_bars + num_nats + i;
                        s.max_ref  = num_refs - 1;
                        s.ref_bits = word64_bits(s.max_ref);

                        FragVal fv = load_fragment(&s);
                        frag_t frag = alloc_frag(ctx, fv);
                        debugf("loaded frag %u\n", frag.ix);
                }
        } else {
                uint32_t num_leaves = num_pins + num_bars + num_nats;

                if (num_leaves == 0) {
                        die("No pins and no leaves.\n");
                }

                if (num_leaves > 1) {
                        die("No frags, but multiple leaves.\n");
                }

                treenode_value v;
                if (num_pins) { v = TAG_PIN((pin_t){0}); }
                else if (num_bars) { v = TAG_BAR((bar_t){0}); }
                else { v = TAG_NAT((nat_t){0}); }
                alloc_treenode(ctx, v);
        }

        if (st.wid != 0) {
                debugf("EXTRA STUFF %u bytes unread!\n", st.wid);
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
        debugf("\t\tpack_bytes_lsb(width=%lu, ", width);
        uint64_t res = 0;
        for (int i=0; i<width; i++) {
                uint64_t tmp = bytes[i];
                debugf("%02lx(%d)", tmp, i*8);
                res = res | (tmp << (i*8));
        }
        debugf(")\n");
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
        // debugf("read_word(%lu)\n", acc);
        debugf("read_word()\n");
        int c;

        while (isdigit(c = getchar())) {
                // debugf("\t    acc=%lu\t('%c' is a digit)\n", acc, c);
                acc *= 10;
                acc += (uint64_t)(c - '0');
        }

        // debugf("\t    acc=%lu\t('%c' is not a digit)\n", acc, c);

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

        debugf("\tread_string() -> \"%s\" (%lu)\n", buf, count);
        // debugf("read_string()\n");


        if (count < 9) {
                leaf.bytes = (uint8_t*) pack_bytes_lsb(count, buf);
                debugf("\t\t(direct)\n");
        } else {
                leaf.bytes = calloc(count+1, 1);
                memcpy(leaf.bytes, buf, count);
                debugf("\t\t(indirect)\n");
        }

        debugf("\t\t(width=%d)\n", leaf.width_bytes);
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
void print_fragment_outline(Jelly*, frag_t);
void print_fragment(Jelly*, frag_t);
void print_tree(Jelly*, treenode_t);

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

void print_tree_list(Jelly *ctx, treenode_t tree) {
        treenode_value val = ctx->treenodes[tree.ix];

        switch (NODEVAL_TAG(val)) {
            case 4: print_pin(ctx, NODEVAL_PIN(val)); break;
            case 5: print_bar(ctx, NODEVAL_BAR(val)); break;
            case 6: print_nat(ctx, NODEVAL_NAT(val)); break;
            case 7: {
                FragVal frag = ctx->frags[NODEVAL_FRAG(val).ix];
                print_tree_list(ctx, frag.head);
                putchar(' ');
                print_tree(ctx, frag.tail);
                break;
            }
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
            case 7: print_fragment(ctx, NODEVAL_FRAG(val)); break;
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

void print_fragment_outline(Jelly *ctx, frag_t ref) {
        FragVal frag = ctx->frags[ref.ix];
        putchar('(');
        print_tree_outline_list(ctx, frag.head);
        putchar(' ');
        print_tree_outline(ctx, frag.tail);
        putchar(')');
}

void print_fragment(Jelly *ctx, frag_t ref) {
        FragVal frag = ctx->frags[ref.ix];
        putchar('(');
        print_tree_list(ctx, frag.head);
        putchar(' ');
        print_tree(ctx, frag.tail);
        putchar(')');
}

int treenode_depth(Jelly *ctx, treenode_t node) {
        treenode_value val = ctx->treenodes[node.ix];
        if (val.word >> 63) return 1;
        int hed_depth = treenode_depth(ctx, NODEVAL_HEAD(val));
        int tel_depth = treenode_depth(ctx, NODEVAL_TAIL(val));
        return (1 + ((hed_depth > tel_depth) ? hed_depth : tel_depth));
}

int fragment_depth(Jelly *ctx, FragVal frag) {
        int hed_depth = treenode_depth(ctx, frag.head);
        int tel_depth = treenode_depth(ctx, frag.tail);
        return (1 + ((hed_depth > tel_depth) ? hed_depth : tel_depth));
}

void jelly_debug_leaves(Jelly *ctx, bool details) {
        debugf("\n\tleaves: (width=%u, count=%u)\n\n",
               ctx->leaves_table_width,
               ctx->leaves_table_count);

        int wid = ctx->leaves_table_width;

        uint64_t mask = ctx->leaves_table_width - 1;

        for (int i=0; i<wid; i++) {
                leaves_table_entry_t ent = ctx->leaves_table[i];

                // Empty slot
                if (ent.pointer.ix == UINT32_MAX) continue;

                debugf("\t\t%4d = ", i);
                print_tree_outline(ctx, ent.pointer);
                debugf("\t(val=");
                print_tree(ctx, ent.pointer);
                debugf(")\n");

                if (details) {
                        uint64_t width = ((ent.width.word << 2) >> 2);
                        uint64_t bin = ent.hash & mask;
                        uint64_t distance = (((uint64_t)i) - bin);
                        debugf("\t\t    bytes: %-4lu\n", width);
                        debugf("\t\t    bin: %lu [dist=%lu]\n", bin, distance);
                        debugf("\t\t    hash: 0x%016lx\n", ent.hash);
                }
        }
}


void jelly_debug_interior_nodes(Jelly *ctx) {
        debugf("\n\tinterior_nodes_table: (width=%u, count=%u)\n\n",
               ctx->nodes_table_width,
               ctx->nodes_table_count);

        int wid = ctx->nodes_table_width;

        uint64_t mask = ctx->nodes_table_width - 1;

        for (int i=0; i<wid; i++) {
                FanEntry ent = ctx->nodes_table[i];

                // Empty slot
                if (ent.val.word == UINT64_MAX) continue;

                uint64_t bin = ent.hash & mask;
                uint64_t distance = (((uint64_t)i) - bin);

                debugf("\t\t%4d = t%u\tbin=%lu\tdist=%lu\thash=%08x\n", i, ent.pointer.ix, bin, distance, ent.hash);

                // debugf("\n\n\t\t  ");
                // print_tree(ctx, ent.pointer);
                // print_tree_outline(ctx, ent.pointer);
                // debugf("\n\n");
        }
}


void jelly_debug(Jelly *ctx) {
        if (!DEBUG) return;

        debugf("\njelly_debug():\n");

        jelly_debug_leaves(ctx, false);

        jelly_debug_interior_nodes(ctx);

        {
                debugf("\n\tpins: (width=%u, count=%u)\n\n",
                       ctx->pins_width,
                       ctx->pins_count);
                int num = ctx->pins_count;
                for (int i=0; i<num; i++) {
                        debugf("\t\tp%d = ", i);
                        print_pin(ctx, (pin_t){ .ix = i });
                        debugf("\n");
                }
        }


        {
                debugf("\n\tbars: (width=%u, count=%u)\n\n",
                       ctx->bars_width,
                       ctx->bars_count);
                int num = ctx->bars_count;
                for (int i=0; i<num; i++) {
                        debugf("\t\tb%d = ", i);
                        print_bar(ctx, (bar_t){ .ix = i });
                        debugf("\n");
                }
        }


        {
                debugf("\n\tnats: (width=%u, count=%u)\n\n",
                       ctx->nats_width,
                       ctx->nats_count);
                int num = ctx->nats_count;
                for (int i=0; i<num; i++) {
                        debugf("\t\tn%d = ", i);
                        print_nat(ctx, (nat_t){ .ix = i });
                        debugf("\n");
                }
        }

        {
                debugf("\n\ttreenodes: (width=%u, count=%u)\n\n",
                       ctx->treenodes_width,
                       ctx->treenodes_count);
                int num = ctx->treenodes_count;
                for (int i=0; i<num; i++) {
                        treenode_value val = ctx->treenodes[i];
                        switch (val.word >> 61) {
                            case 4:
                                debugf("\t\tt%d = p%u\n", i, (uint32_t) val.word);
                                continue;
                            case 5:
                                debugf("\t\tt%d = b%u\n", i, (uint32_t) val.word);
                                continue;
                            case 6:
                                debugf("\t\tt%d = n%u\n", i, (uint32_t) val.word);
                                continue;
                            case 7:
                                debugf("\t\tt%d = f%u\n", i, (uint32_t) val.word);
                                continue;
                            default: {
                                treenode_t hed = NODEVAL_HEAD(val);
                                treenode_t tel = NODEVAL_TAIL(val);
                                debugf("\t\tt%d = (%u, %u)\n", i, hed.ix, tel.ix);
                                continue;
                            }
                        }
                }
        }

        debugf("\nJelly Fragments:\n\n");

        uint32_t count = ctx->frags_count;
        for (int i=0; i<count; i++) {
                FragVal cur = ctx->frags[i];
                int depth = fragment_depth(ctx, cur);
                debugf("\tFragment[%d]: (depth=%d)\n\n\t\t(t%d, t%d) = ", i, depth, cur.head.ix, cur.tail.ix);
                print_fragment_outline(ctx, (frag_t){i});

                debugf("\n\n\t\t    ");
                print_fragment(ctx, (frag_t){i});
                debugf("\n\n");
        }
}

void jelly_print(Jelly *ctx) {
        uint32_t frags = ctx->frags_count;

        if (frags) {
                print_fragment(ctx, (frag_t){ .ix = (frags-1) });
        } else {

                if (!ctx->treenodes_count) {
                        die("No frags and no leaves\n");
                }

                if (ctx->treenodes_count > 1) {
                        die("No frags but %d leaves.  Nonsense!\n", ctx->treenodes_count);
                }

                print_tree(ctx, (treenode_t){0});
        }

        putchar('\n');
}


int main () {
        Jelly *ctx = new_jelly_ctx();

        // hash256_t dumb_hash_1 = { fmix64(111111), fmix64(65535), fmix64(9),  65536 };
        // hash256_t dumb_hash_2 = { fmix64(222222), fmix64(33333), (0ULL - 1), (65535ULL << 12) };
        treenode_t top = read_some(ctx);
        // treenode_t tmp = jelly_pin(ctx, &dumb_hash_1);
        // top = jelly_cons(ctx, tmp, top);
        // tmp = jelly_pin(ctx, &dumb_hash_1);
        // top = jelly_cons(ctx, tmp, top);
        // tmp = jelly_pin(ctx, &dumb_hash_2);
        // top = jelly_cons(ctx, tmp, top);

        debugf("# Shatter\n");

        shatter(ctx, top, true);

        jelly_debug(ctx);

        debugf("# Serialize\n");

        struct ser ser = serialize(ctx);

        debugf("# Clean Slate\n");
        wipe_jelly_ctx(ctx);
        // jelly_debug(ctx);

        debugf("# De-serialize\n");
        deserialize(ctx, ser);

        jelly_debug(ctx);

        jelly_print(ctx);

        free(ser.buf);
        free_jelly_ctx(ctx);

        return 0;
}
