#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "jelly.h"
#include "libbase58.h"


// Forward Declarations ////////////////////////////////////////////////////////

treenode_t read_one(Jelly);
treenode_t read_many(Jelly);


// Testing Harness /////////////////////////////////////////////////////////////

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

// TODO: Instead of reading pairs, use an accumulator byte like when
// we bit-deserialize.
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
        debugf("read_word()\n");
        int c;

        while (isdigit(c = getchar())) {
                acc *= 10;
                acc += (uint64_t)(c - '0');
        }

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

void eat_comment() {
    loop:
        int c = getchar();
        switch (c) {
            case EOF:
            case '\n':
                return;
            default:
                goto loop;
        }
}

void eat_space() {
    loop:
        int c = getchar();
        switch (c) {
            case '#':
                eat_comment();
                goto loop;
            case '\n':
            case '\t':
            case ' ':
                goto loop;
            default:
                ungetc(c, stdin);
                return;
        }
}

treenode_t read_leaf(Jelly ctx) {
        int c = getchar();

        switch (c) {
            case '"':
                return jelly_bar(ctx, read_string());
            case '[':
                die("TODO: Parse pins\n");
            case '0': {
                int d = getchar();
                if (d == 'x') {
                        return jelly_nat(ctx, read_hex());
                }
                ungetc(d, stdin);
            } // fallthrough
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                uint64_t word = read_word((uint64_t)(c - '0'));
                uint32_t width = word64_bytes(word);
                return jelly_packed_nat(ctx, width, word);
            default:
                die("Not a leaf: '%c'\n", c);
        }
}

treenode_t read_many(Jelly ctx) {
       treenode_t acc = read_one(ctx);
    loop:
        eat_space();

        int c = getchar();
        switch (c) {
            case '(':
                treenode_t list = read_many(ctx);
                acc = jelly_cons(ctx, acc, list);
                goto loop;
            case EOF:
            case ')':
                return acc;
            default:
                ungetc(c, stdin);
                treenode_t elmt = read_leaf(ctx);
                acc = jelly_cons(ctx, acc, elmt);
                goto loop;
        }
}

treenode_t read_one(Jelly ctx) {
        eat_space();
        int c = getchar();
        switch (c) {
            case '(':
                return read_many(ctx);
            case EOF:
                die("Unexpected EOF (read_one)\n");
            default:
                ungetc(c, stdin);
                return read_leaf(ctx);
        }
}



int main () {
        Jelly ctx = jelly_new_ctx();

        // hash256_t dumb_hash_1 = { fmix64(111111), fmix64(65535), fmix64(9),  65536 };
        // hash256_t dumb_hash_2 = { fmix64(222222), fmix64(33333), (0ULL - 1), (65535ULL << 12) };
        read_many(ctx);

        // treenode_t top = read_many(ctx);
        // treenode_t tmp = jelly_pin(ctx, &dumb_hash_1);
        // top = jelly_cons(ctx, tmp, top);
        // tmp = jelly_pin(ctx, &dumb_hash_1);
        // top = jelly_cons(ctx, tmp, top);
        // tmp = jelly_pin(ctx, &dumb_hash_2);
        // top = jelly_cons(ctx, tmp, top);

        debugf("# Shatter\n");

        jelly_finalize(ctx);

        jelly_debug(ctx);

        debugf("# Serialize\n");

        struct ser ser = jelly_serialize(ctx);

        debugf("# Clean Slate\n");
        jelly_wipe_ctx(ctx);
        // jelly_debug(ctx);

        debugf("# De-serialize\n");
        jelly_deserialize(ctx, ser);

        jelly_debug(ctx);

        jelly_print(ctx);

        free(ser.buf);
        jelly_free_ctx(ctx);

        return 0;
}
