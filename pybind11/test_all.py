#!/usr/bin/env python3

import numpy as np
import pet
import iteration_mod as i_mod
import example

a = np.random.normal(0, 1, (3, 4)).astype('f')
print("")
print("Testing 'pet':")
print("--------------")
p = pet.Pet("Cassie")
print(p.getName())
p.setName("Ginger")
print(p.getName())
print("Testing a * (2 + 1j)")
print(f"a = {a}")
print(p.nump(a))

print("")
print("Testing 'iteration_mod':")
print("------------------------")
print(a)
print(i_mod.nump(a))

b = np.random.randn(25)
c = np.random.randn(25)
d = np.zeros(25)
e = i_mod.add_arrays(b, c, d)
print("")
print("i_mod.add_arrays:")
#  print(f"in1: {b}")
#  print(f"in2: {c}")
#  print(f"in1 + in2: {b + c}")
#  print(f"i_mod.add_arrays(in1, in2): {d}")
print(f"sum(|delta|) = {np.sum(np.abs(d - b - c))}")
print(f"sum(|delta|) = {np.sum(np.abs(d - e))}")

print("")
print("Testing 'example':")
print("------------------")
print(example.add1(10, 20))
print(example.add1(j=20))
print(example.add1(i=20))
print(example.add2(10, 20))
print(example.add2(j=20))
print(example.add2(i=20))
print("")
