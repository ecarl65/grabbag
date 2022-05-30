
CC=gcc
CFLAGS=-g
LIBS=-lm -lfftw3f -lfftw3

all: test_fftw osconv

test_fftw: test_fftw.c
	$(CC) -std=c99 -o $@ $< $(CFLAGS) $(LIBS)

osconv: osconv.c
	$(CC) -std=c99 -o $@ $< $(CFLAGS) $(LIBS)

.phony: clean

clean:
	rm -f test_fftw osconv

