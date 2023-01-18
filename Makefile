default: ray3.c
	gcc ray3.c -o ray3 -g -Wall -std=gnu99 -lm

clean:
	rm -f ray3 reference.png
