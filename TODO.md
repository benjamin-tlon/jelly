-   TODO The front-end should own all the binary data for bars, bignats,
    and pins.

    That data should not be freed with the context.

-   TODO Maybe we should have a single `leaves` array, and have each
    be tagged with it's type.

    This would simplify the code for deduplication of naturals, at the
    expense of slightly worse bucketing.

    Instead of having a table of words, a table of pins, a table of bars,
    and a table of nats, we would have a table of tagged leaves:

        struct leaves_table_entry {
            uint64_t hash;
            tagged_width_t width;
            uint8_t *bytes;
            int64_t offset; // Index into pins/bars/nats (depending on the tag).
        }

    So, you would hash, tag the width, calculate the starting offset,
    and then just scan the table.

    The duplicates case would be slightly slower for direct words, but
    not much.  The `bytes` field would still be inline, and we would have
    a custom code-path for searching for data the fits in a word.

    In the search-for-word path, we would ignore the hash and just
    directly compare the `bytes` and `width` fields.  Most failing cases
    are still only one compare, successful cases is two comparisions.

    It also doubles the width of the words table, but it saves us from
    needing to be constantly searching accross four tables, which has
    worse memory locality.

    This is also slightly slower for pins, but barely.  The hash (just
    the first word of the pin) will almost never match for equal values,
    so we get single-compare checks in most cases.  A successful compare
    is (hash + tagged_width + bytes[1] + bytes[2] + bytes[3]).

    Yeah, that's interesting, we can still have specialized insert logic
    per-type, but still have a single data structure, a single routine
    to grow hash tables, etc.

            uint64_t hash;
            tagged_width_t width;
            uint8_t *bytes;
            int64_t offset; // Index into pins/bars/nats (depending on the tag).

-   TODO: Decide about details of hash table implementation.

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

## Deduplication

-   TODO Deduplication support for large nats.
-   TODO Parser/printer support for large nats.

-   TODO Deduplication support for bars.
-   TODO Parser/printer support for bars.

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
