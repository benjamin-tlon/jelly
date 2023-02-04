a.out: jelly.o murmur3.o xxhash.o xxh_x86dispatch.o base58.o
	gcc -Wall $^

jelly.o: jelly.c murmur3.h libbase58.h xxhash.h
	gcc -Wall -c $<

murmur3.o: murmur3.c murmur3.h
	gcc -Wall -c $<

base58.o: base58.c libbase58.h
	gcc -Wall -c $<

xxhash.o: xxhash.c xxhash.h
	gcc -Wall -c $<

xxh_x86dispatch.o: xxh_x86dispatch.c xxhash.h
	gcc -Wall -c $<

clean:
	rm -f *.o a.out
