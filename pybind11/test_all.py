#!/usr/bin/env python3

import numpy as np
import pet
import iteration_mod as i_mod
import example

a = np.random.normal(0, 1, (3, 4)).astype('f')
p = pet.Pet("Cassie")
assert "Cassie" == p.getName(), "Pet name incorrect"
p.setName("Ginger")
assert "Ginger" == p.getName(), "Pet name incorrect"
a_t = a * (2 + 1j)
a_c = a + 0j
a_m = p.nump(a_c)
assert np.allclose(a_t, a_m), "Pet.nump incorrect"

a_t = a.astype('d') * 2
a_m = i_mod.nump(a.astype('d'))
assert np.allclose(a_t, a_m), "iteration_mod.nump incorrect"

b = np.random.randn(25)
c = np.random.randn(25)
d = np.zeros(25)
e = i_mod.add_arrays(b, c, d)
assert np.allclose(b + c, d), "iteration_mod.add_arrays correct filling in last argument"
assert np.allclose(b + c, e), "iteration_mod.add_arrays correct returning argument"
assert example.add1(10, 20) == 30, "example.add1(10, 20) producing right result"
assert example.add1(j=20) == 21, "example.add1(j=20) producing right result"
assert example.add1(i=20) == 22, "example.add1(i=20) producing right result"
assert example.add2(10, 20) == 30, "example.add2(10, 20) producing right result"
assert example.add2(j=20) == 21, "example.add2(j=20) producing right result"
assert example.add2(i=20) == 22, "example.add2(i=20) producing right result"

f = np.random.normal(0, 1, (3, 4)) + 1j * np.random.normal(0, 1, (3, 4))
f_t = f * (2 + 2j)
f_m = i_mod.cnump(f)
assert np.allclose(f_t, f_m), "iteration_mod.cnump produces correct output"

print("All tests pased")
