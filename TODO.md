-   TODO Cannibalize the refcounts table.

    The refcounts taable is (32bits * num_nodes), which is always bigger
    than (32bits * max_depth).

    We only need it for shatter(), so we can just canibalize it, and
    use that memory for the fragment-serialization stack.

    -   TODO Remove the `stack[128] /* TODO */` nonsense, and just reuse
        memory like the above suggests.

-   TODO Deserialize to a fresh jelly context.

    -   Just read each leaf, and append directly to the respective tables

    -   Do not populate the dedupe tables.

    -   Then, read each fan.

    -   All references get looked up in the prior tables.

    -   Reading the fans will definitely be the trickiest part here,
        since we will need to maintain bit-offset counts, like we do
        on serialization.

-   TODO Final test-harness print should be the deserialized context.


## Misc

-   TODO Parser support for pins.

-   TODO The front-end should own all the binary data for bars, bignats,
    and pins.  It shouldn't be freed with the Jelly* context, and instead
    should be freed directly from main().


## Encoding and Decoding

-   TODO Write code to decode a larger byte-string.

    - high bit is 0xxxxxxx? -> read one byte
    - high two bits are 10xxxxx? -> read 0bxxxxxx byte value
    - high two bits are 11xxxxx? -> read 0bxxxxxx byte length (and then length-byte value)

-   TODO Write code to decode the leaf-set:

-   TODO Make sure that malicious inputs cannot cause us to read past
    the end of a buffer.

    That's easy enough, I think.  Just store the end position, and make
    sure we don't go past that.

    I guess we need to keep track of that at every step, which is
    annoyingly expensive.

-   TODO Fragment serializer uses words instead of bytes.

    -   DONE Make sure that the buffer-size that we allocate has a size that's a multiple of 8.

    -   TODO Calculate the next word that contains the starting byte.

    -   TODO Calculate the bit-index into that word where the writing should start.

    -   TODO Initialize the accumulator by reading the current state of that work.

    -   TODO Replace 8 with 64 as the bit-offset in which we flush to the buffer.
