-   TODO Make sure no entity hash is ever zero.

-   TODO The front-end should own all the binary data for bars, bignats,
    and pins.  It shouldn't be freed with the Jelly* context, and instead
    should be freed directly from main().

-   TODO The `leaf_table_insert` routines (direct and indirect) should
    be passed a continuation.

    This way the type-specific logic can be factored out.  The overhead
    of an indirect call should be minor, since we are only doing it once
    per (unique) leaf.

-   DONE: Support bars.
-   TODO: Support pins.

-   TODO: Make the leaves table into an actual hash table.

    Hash table capacity will be a power of two.

    Hash table capacity will double whenever it becomes half full.

    We calculate the offset by using bit-wise and against the hash.
    We store the mask in the environment to avoid recalculating it.

    We use a zero-hash to indicate an empty field.  Our hash function
    never produces zero, because we catch zero hashes and replace them
    with one.

    So, to insert:

            uint64_t index = hash
        loop:
            index &= mask
            if (table[index]->hash == 0) {
                INSERT;
                load++;
                if ((load * 2) > capacity) GROW;
            } else if (MATCH(value, table[index])) {
                FOUND;
            } else {
                index++;
                goto loop;
            }

    The specifics of FOUND and INSERT are specialized per entry-type,
    but the rest of the logic stays the same.

-   TODO: Make the leaves table into an actual hash table.

## Deduplication

-   TODO Deduplication support for pins.
-   TODO Parser/printer support for pins.

-   TODO Deduplication support for fan backreferences.
-   TODO Parser/printer support for fan backreferences.

-   TODO Calculate resulting buffer size.


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
