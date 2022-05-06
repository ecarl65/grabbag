#!/usr/bin/env python3.8

import pet
import iteration_mod as i_mod
import numpy as np
import example

print("")
print("Testing 'pet':")
print("--------------")
p = pet.Pet("Cassie")
print(p.getName())
p.setName("Ginger")
print(p.getName())

print("")
print("Testing 'iteration_mod':")
print("------------------------")
a = np.random.normal(0, 1, (3, 4))
print(a)
print(i_mod.nump(a))

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

