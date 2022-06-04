#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def f(x):
    return 0.2 * x**2 + 0.5 * x - 5

x = np.linspace(-10, 10, 1000)
y = f(x)

r1 = optimize.newton(f, -5)
r2 = optimize.newton(f, 5)

print(f"Newton's Method Optimization:\nRoot 1: {r1}\nRoot 2: {r2}")

r = np.roots([0.2, 0.5, -5])
print(f"Root finder results: {r}")

plt.plot(x, y)

plt.show()
