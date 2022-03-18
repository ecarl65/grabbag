#!/usr/bin/env python

import integrate

if __name__ == '__main__':
    a = 0
    b = 1000
    N = 10000000
    res = integrate.integrate_f(a, b, N)

    print(f"Integral of x^2 - x from {a} to {b} with {N} points is {res}")

