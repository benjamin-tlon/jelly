#define DEBUG 0

#if DEBUG
#define debugf(...) printf(__VA_ARGS__)
#else
static inline void pass() {}
#define debugf(...) pass(__VA_ARGS__)
#endif

#define die(...) (fprintf(stderr, __VA_ARGS__),exit(1))


// Inlines /////////////////////////////////////////////////////////////////////

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

Jelly jelly_new_ctx();
void jelly_wipe_ctx(Jelly);
void jelly_free_ctx(Jelly);

void jelly_finalize(Jelly);
void jelly_debug(Jelly);
void jelly_print(Jelly);

treenode_t jelly_bar(Jelly, leaf_t);
treenode_t jelly_nat(Jelly, leaf_t);
treenode_t jelly_cons(Jelly, treenode_t, treenode_t);
treenode_t jelly_packed_nat(Jelly, uint32_t, uint64_t);

struct ser jelly_serialize(Jelly);

void jelly_deserialize(Jelly, struct ser);

