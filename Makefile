a.out: jelly.o murmur3.o xxhash.o xxh_x86dispatch.o base58.o harness.o
	gcc -O3 -Werror -Wall $^

harness.o: harness.c jelly.h
	gcc -O3 -Werror -Wall -c $<

jelly.o: jelly.c murmur3.h libbase58.h xxhash.h
	gcc -O3 -Werror -Wall -c $<

murmur3.o: murmur3.c murmur3.h
	gcc -O3 -Werror -Wall -c $<

base58.o: base58.c libbase58.h
	gcc -O3 -Werror -Wall -c $<

xxhash.o: xxhash.c xxhash.h
	gcc -O3 -Werror -Wall -c $<

xxh_x86dispatch.o: xxh_x86dispatch.c xxhash.h
	gcc -O3 -Werror -Wall -c $<

clean:
	rm -f *.o a.out test.tmp test.tmp?
