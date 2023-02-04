a.out: jelly.o murmur3.o xxhash.o xxh_x86dispatch.o
	gcc -Wall $^

jelly.o: jelly.c
	gcc -Wall -c $<

murmur3.o: murmur3.c
	gcc -Wall -c $<

xxhash.o: xxhash.c
	gcc -Wall -c $<

xxh_x86dispatch.o: xxh_x86dispatch.c
	gcc -Wall -c $<

clean:
	rm -f *.o a.out
