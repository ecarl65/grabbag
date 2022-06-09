
CC=gcc
# CFLAGS=-g -DDEBUG
CFLAGS=-g
LIBS=-lm -lfftw3f -lfftw3

all: test_fftw osconv dsconv poly udft ovudft

test_fftw: test_fftw.c
	$(CC) -std=c99 -o $@ $< $(CFLAGS) $(LIBS)

osconv: osconv.c
	$(CC) -std=c99 -o $@ $< $(CFLAGS) $(LIBS)

dsconv: dsconv.c
	$(CC) -std=c99 -o $@ $< $(CFLAGS) $(LIBS)

poly: poly.c
	$(CC) -std=c99 -o $@ $< $(CFLAGS) $(LIBS)

udft: udft.c
	$(CC) -std=c99 -o $@ $< $(CFLAGS) $(LIBS)

ovudft: ovudft.c
	$(CC) -std=c99 -o $@ $< $(CFLAGS) $(LIBS)

.phony: clean

clean:
	rm -f *.bin test_fftw osconv dsconv poly udft ovudft

