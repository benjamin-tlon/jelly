a.out: jelly.o murmur3.o xxhash.o xxh_x86dispatch.o
	gcc $^

jelly.o: jelly.c
	gcc -c $<

murmur3.o: murmur3.c
	gcc -c $<

xxhash.o: xxhash.c
	gcc -c $<

xxh_x86dispatch.o: xxh_x86dispatch.c
	gcc -c $<

clean:
	rm -f *.o a.out
