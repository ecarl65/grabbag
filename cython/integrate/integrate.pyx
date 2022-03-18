#!/usr/bin/env python

import cython

@cython.cfunc
@cython.exceptval(-2, check=True)
@cython.boundscheck(False)
@cython.wraparound(False)
def f(x: cython.double) -> cython.double:
    return x ** 2 - x

@cython.boundscheck(False)
@cython.wraparound(False)
def integrate_f(a: cython.double, b: cython.double, N: cython.int) -> cython.double:
    i: cython.int
    s: cython.double
    dx: cython.double
    out: cython.double
    s = 0
    dx = (b - a) / N
    for i in range(N):
        s += f(a + i * dx)
    out = s * dx
    return out

if __name__ == '__main__':
    a = 0
    b = 1000
    N = 10000000
    res = integrate_f(a, b, N)

    print(f"Integral of x^2 - x from {a} to {b} with {N} points is {res}")
