
CC=gcc
CFLAGS=
LIBS=-lm -lfftw3f

test_fftw: test_fftw.c
	$(CC) -std=c99 -o $@ test_fftw.c $(CFLAGS) $(LIBS)

.phony: clean

clean:
	rm -f test_fftw 

