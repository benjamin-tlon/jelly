-   TODO Track the recursion depth.

    Just store a ctx->maximum_depth field.

    This will be an over-estimate in some cases (if there are
    back-references, then it's possible that the deepest part will live
    in a fragment).  However, it is still a valid upper-bound and it is
    still a reasonable estimate.

    -   TODO Use the maximum_depth field to pre-allocate a stack during
        bit-serialization, don't use recursion.

-   DONE Internior nodes table.
-   TODO treenode_equals
-   TODO Interior node deduplication (and backrefs)
-   TODO Printer shows backreferenes table.

-   TODO Calculate resulting buffer size.
-   TODO Actually serialize.
-   TODO Deserialize to a fresh jelly context.
-   TODO Final test-harness print should be the deserialized context.


## Misc

-   DONE Q: Can we reclaim space when we find duplicates?

    When we find a duplicate tree-node and turn that into a
    back-reference, all of the nodes in the duplicate tree "leak".
    They take up space, but will never be referenced.

    Can we immediatly release them, and re-use that space?

    Intuitively, these nodes should be the *last n* nodes for some value
    of `n`, in which case, we could simply decrease the "used" count.

    Is that true?  And if so, what is the that of n?

    is is `F(num_leaves)` for some function F?

    Let's look at an example:

        #((0 0) (0 0))    ->
        (#(0 0) (0 0))    ->
        ((#0 0) (0 0))    ->    []                  n0=0
        ((0 #0) (0 0))    ->    []                  n0=0 n0=0
        ((0 0) #(0 0))    ->    []                  n0=0 n0=0 (0,1)
        ((0 0) (#0 0))    ->    []                  n0=0 n0=0 (0,1)
        ((0 0) (0 #0))    ->    []                  n0=0 n0=0 (0,1) n0=0 n0=0
        ((0 0) ^(0 0))    ->    []                  n0=0 n0=0 (0,1) n0=0 n0=0 (0,1)

                          ->    [f0=(0,1)]          n0=0 n0=0 {f0} n0=0 n0=0 f0
                          ->    [f0=(0,1)]          n0=0 n0=0 f0 {n0=0 n0=0 f0}
                          ->    [f0=(0,1)]          n0=0 n0=0 f0
                 [here, 3 nodes can be reclaimed (2*2 - 1)]
                 [And they are the last 3 nodes]
                 [So, indeed, pop((leaves*2)-1) seems to work]

        ^((0 0) (0 0))    ->    [f0=(0,1)]          n0=0 n0=0 f0
        ^((0 0) (0 0))    ->    [f0=(0,1) f1=(2,2)] n0=0 n0=0 f0 (2,2)

    So, first of all: yes that seems to work great.

    But also, it seems like we're missing an opportunity to apply the
    same logic to matching laves.  Why duplicate those nodes?

        #((0 0) (0 0))    ->
        (#(0 0) (0 0))    ->
        ((^0 0) (0 0))    -> t:{0:n0}
        ((t0 ^0) (0 0))   -> t:{0:n0}
        (^(t0 t0) (0 0))  -> t:{0:n0, 1:(t0 t0)}
        (t1 #(0 0))       -> t:{0:n0, 1:(t0 t0)}
        (t1 (^0 0))       -> t:{0:n0, 1:(t0 t0)}
        (t1 (t0 ^0))      -> t:{0:n0, 1:(t0 t0)}
        ^(t1 t1)          -> t:{0:n0, 1:(t0 t0)}
        t2                -> t:{0:n0, 1:(t0 t0), 2:(t1 t1)}

    Alright, sick.  That seems quite possible.

    However, the number of nodes to release upon (interior-node-match)
    is now more complicated, since we didn't allocate every leaf.

    So, what we really want, I think, is to track the number of treenodes
    that have been allocated by this particular recursion-point.

    That goes up when we allocate a node (find a unique leaf or unique
    pair).  Whenever we find a matching tree, we pop(n) for the number
    of allocations that we own.

    Simple!

    This optimization will not only avoid wasting memory, but it will
    also improve the memory locality of node-equality tests.

    Furthermore, I believe it can maybe let us short-circuit tree-equality
    checks, since every equal sub-node will always be pointer-equals
    (except for the top-most newly-allocated entry).

        For example, if we are comparing:

            ((0 1) (2 3))

        with

            ((0 1) (2 3))

        Then we will have already compared (0 1) and (2 3) for equality.

        So, we are actually just comparing:

            (f0 f1) against (f0 f1)

        If we ever ask if two sub-cells are equal, we automatically
        know they are NOT EQUAL, since they would already be (equal)
        backreferences if they were equal trees.

        This allows us to completely short-circuit the equality check.

        Yes, I think even the treenode_value_t (packed word) will be
        directly word-equals if the two trees are equal.

        Let's go through the example again:

            #((4 8) (4 8))    -> 
            (#(4 8) (4 8))    -> 
            ((@0 ^8) (4 8))   ->  n:{4}   t:{n0}
            (^(@0 @1) (4 8))  ->  n:{4 8} t:{n0 n1}
            (@2 #(4 8))       ->  n:{4 8} t:{n0 n1 {0,1}}
            (@2 (^4 8))       ->  n:{4 8} t:{n0 n1 {0,1}}
            (@2 (@0 ^8))      ->  n:{4 8} t:{n0 n1 {0,1}}
            (@2 ^(@0 @1))     ->  n:{4 8} t:{n0 n1 {0,1}}
                                                           Here we see
                                                           that we are
                                                           comparing {0,1}
                                                           with {0,1},
                                                           which is a
                                                           word-equal
                                                           comparison.
            (@2 @2)           ->  n:{4 8} b:{(0,1)} t:{n0 n1 b0}

    Conclusion:

        The internal-nodes deduplication table need only be:

            HashMap (Word32, Word32) Int

        Even the mutation that happens when a backreference was found
        does not invalidate this key, because the value that it matches
        against is still an (x,y) pair.

        We can implement this in C as:

            //
            // Use pair = UINT64_MAX to indicate an empty
            // marker.  Use `memset(buf, 255, size)` to initialize.
            //
            // The `hash` field is only used for re-hashing (growing
            // the table).  When we scan the table we only ask:
            //
            //     NodeEntry *cur = ctx->node_table + i;
            //     uint64_t pair = cur->pair
            //     if (pair == target || pair == UINT64_MAX) {
            //             return cur;
            //     }
            //
            typedef struct {
                uint64_t pair;
                uint32_t hash;    // fmix64 truncated to 32 bits.
                uint32_t offset;
            } NodeEntry;

-   TODO Parser support for pins.

-   TODO Make sure no entity hash is ever zero.

    -   TODO Is there some other field we could look at?

    -   TODO This is an annoying invariant to enforce, error-prone.

-   TODO The front-end should own all the binary data for bars, bignats,
    and pins.  It shouldn't be freed with the Jelly* context, and instead
    should be freed directly from main().


## Encoding and Decoding

-   TODO Write code to dump a larger byte-string.

    -   number<128 -> one byte
    -   size<128   -> one byte length + bytes
    -   otherwise  -> one byte length-of-length + length + bytes

    TODO: Fast way to byte-dump a word?  Just mask and shift (n<8) times?

-   TODO Write code to decode a larger byte-string.

    - high bit is 0xxxxxxx? -> read one byte
    - high two bits are 10xxxxx? -> read 0bxxxxxx byte value
    - high two bits are 11xxxxx? -> read 0bxxxxxx byte length (and then length-byte value)

-   TODO Write code to encode the leaf-set:

    -   length_prefixed(pins)
    -   length_prefixed(varchar(bar) * bars)
    -   length_prefixed(varchar(nat) * nats)

-   TODO Write code to decode the leaf-set:

-   TODO Make sure that malicious inputs cannot cause us to read past
    the end of a buffer.

    That's easy enough, I think.  Just store the end position, and make
    sure we don't go past that.

    I guess we need to keep track of that at every step, which is
    annoyingly expensive.

-   DONE Determine the maximum stack-size required for serializaion.

    It's just the tree depth.

-   TODO Write code to bit-encode a fan with fixed-size leaves.

    -   Do this byte-wise at first.  Leave word-wise encoding as an optimization.
    -   As a short-term hack, only support working with leaf sizes <9 bits.

    -   Algorithm:

        -   One function call, use goto for flow control.

        -   Precalculate maximum stack size (workspace entry should record this).

        -   Maintain State

            -   a stack (array of int)
            -   current stack pointer (index into array).
            -   A pointer into the output buffer.
            -   An accumulator byte.
            -   An overflow byte.
            -   The number of bits that have been written to the accumulator.

        -   To write n bits:

                acc          |= ((value << bits_written) >> bits_written)
                overflow      = (value >> (8 - bits_written))
                bits_written += value_width

                if (bits_written >= 8) {
                    *buf++ = acc;
                    acc = overflow
                    bits_written -= 8.
                }

        -   To advance the position (by one bit).

                bits_written = (bits_written + 1) % 8;
                if (!bits_written) {
                    *buf++ = acc;
                    acc = 0;
                }

        -   If we encounter a tree-node, set the current bit.

                acc |= (1 << bits_written);
                ADVANCE;
                push(tail);
                push(head);
                goto loop;

        -   If we encounter a leaf-node.

                offsets = [0, pins_count, pins_count+bars_count, pins_count+bars_count+nats_coun]
                ADVANCE;
                val = node_value + offsets[node_type]
                WRITE_BITS;

        -   Finally

                if (bits_written) {
                    *buf = acc
                }

    -   Supporting word buffering:

        -   Make sure that the buffer-size that we allocate has a size that's a multiple of 8.

        -   Calculate the next word that contains the starting byte.

        -   Calculate the bit-index into that word where the writing should start.

        -   Initialize the accumulator by reading the current state of that work.

        -   Replace 8 with 64 as the bit-offset in which we flush to the buffer.
