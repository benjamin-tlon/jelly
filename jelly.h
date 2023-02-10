#define DEBUG 0

#if DEBUG
#define debugf(...) printf(__VA_ARGS__)
#else
static inline void pass() {}
#define debugf(...) pass(__VA_ARGS__)
#endif

#define die(...) (fprintf(stderr, __VA_ARGS__),exit(1))

#define INLINE __attribute__((always_inline)) inline

#define BIG_CONSTANT(x) (x##ULL)


// Inlines /////////////////////////////////////////////////////////////////////

/*
        Taken from murmur3.  Used to bit-mix direct atoms and
        shape-signatures.
*/
static INLINE uint64_t fmix64 (uint64_t k) {
        k ^= k >> 33;
        k *= BIG_CONSTANT(0xff51afd7ed558ccd);
        k ^= k >> 33;
        k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
        k ^= k >> 33;

        return k;
}

static inline int
word64_bits (uint64_t w) {
        if (!w) { return 0; }
        return 64 - __builtin_clzll(w);
}

static inline uint32_t
word64_bytes (uint64_t w) {
        uint32_t bits = word64_bits(w);
        uint32_t aligned = bits/8;
        if (bits%8) return aligned + 1;
        return aligned;
}


// Types ///////////////////////////////////////////////////////////////////////


typedef struct hash256 { uint64_t a, b, c, d; } hash256_t;

typedef struct { uint32_t ix; } pin_t;
typedef struct { uint32_t ix; } bar_t;
typedef struct { uint32_t ix; } nat_t;
typedef struct { uint32_t ix; } frag_t;
typedef struct { uint32_t ix; } treenode_t;

typedef struct jelly_ctx *Jelly;

// If width<9, then the bytes are packed directly into the `bytes` slot.
typedef struct leaf {
        int32_t width_bytes;
        uint8_t *bytes;
} leaf_t;

struct ser {
        uint8_t *buf;
        size_t wid;
};


// Functions ///////////////////////////////////////////////////////////////////

Jelly jelly_make();

void jelly_wipe(Jelly);
void jelly_free(Jelly);
void jelly_done(Jelly);
void jelly_dbug(Jelly);
void jelly_show(Jelly);

// TODO: Tell us if it was unique!
treenode_t jelly_pin(Jelly, bool*, uint8_t*);
treenode_t jelly_bar(Jelly, size_t, uint8_t*);
treenode_t jelly_nat(Jelly, size_t, uint8_t*);
treenode_t jelly_cons(Jelly, treenode_t, treenode_t);
treenode_t jelly_word(Jelly, uint64_t);

size_t jelly_size(Jelly);

void jelly_dump(Jelly, size_t, uint8_t*);
void jelly_load(Jelly, size_t, uint8_t*);
