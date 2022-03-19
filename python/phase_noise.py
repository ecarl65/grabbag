### current bpython session - make changes and save to reevaluate session.
### lines beginning with ### will be ignored.
### To return to bpython without reevaluating make no changes to this file
### or save an empty file.
w 
# OUT: Traceback (most recent call last):
# OUT:   File "<input>", line 1, in <module>
# OUT:     w
# OUT: NameError: name 'w' is not defined
import numpy as np
w = np.dot(np.array([1, 1j]), np.random.randn(2, 1000))
plt.figure(); plt.plot(np.real(w[:100]), np.imag(w[:100]), 'b.'); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b92c2b370>]
plt.figure(); plt.plot(np.real(w), np.imag(w), 'b.'); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b90980610>]
np.std(w)
# OUT: 1.4603238234390632
w = np.dot(np.array([1, 1j]), np.random.randn(2, 1000)) / np.sqrt(2)
np.std(w)
# OUT: 0.9915283355851306
plt.figure(); plt.plot(np.real(w), np.imag(w), 'b.'); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b92c2a1a0>]
s = np.exp(1j * 2 * np.pi * 0.1 * np.arange(1000))
sw = s + w
plt.figure(); plt.plot(np.angle(np.unwrap(sw))); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: Traceback (most recent call last):
# OUT:   File "<input>", line 1, in <module>
# OUT:     plt.figure(); plt.plot(np.angle(np.unwrap(sw))); plt.show()
# OUT:   File "<__array_function__ internals>", line 5, in unwrap
# OUT:   File "/home/eric/pyvenv/lib64/python3.10/site-packages/numpy/lib/function_base.py", line 1585, in unwrap
# OUT:     ddmod = mod(dd - interval_low, period) + interval_low
# OUT: TypeError: ufunc 'remainder' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''
plt.figure(); plt.plot(np.unwrap(np.angle(sw))); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b92cad150>]
w = np.dot(np.array([1, 1j]), np.random.randn(2, 1000)) / np.sqrt(2) / 100
sw = s + w
plt.figure(); plt.plot(np.unwrap(np.angle(sw))); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b90982380>]
plt.figure(); plt.plot(np.unwrap(np.angle(sw))); plt.plot(2 * np.pi * 0.1 * np.arange(1000), 'r'); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b92cad780>]
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b92bfd5d0>]
phn = np.unwrap(np.angle(sw)) - (2 * np.pi * 0.1 * np.arange(1000))
plt.figure(); plt.plot(phn); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b8ff98970>]
np.std(w)
# OUT: 0.01005740750865678
plt.figure(); plt.plot(phn); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b8fdeb1f0>]
s = 30 * np.exp(1j * 2 * np.pi * 0.1 * np.arange(1000))
sw = s + w
phn = np.unwrap(np.angle(sw)) - (2 * np.pi * 0.1 * np.arange(1000))
plt.figure(); plt.plot(phn); plt.show()
# OUT: <Figure size 640x480 with 0 Axes>
# OUT: [<matplotlib.lines.Line2D object at 0x7f9b909b0ee0>]
np.std(w)
# OUT: 0.01005740750865678
np.std(w) / 30
# OUT: 0.00033524691695522603
np.std(phn)
# OUT: 0.0002427803094646874
np.std(w) / 30 / np.sqrt(2)
# OUT: 0.00023705536835092367
np.std(phn) / np.sqrt(30)
# OUT: 4.432541733729809e-05
np.std(phn) / 30 / np.sqrt(2)
# OUT: 5.722386772034967e-06
x = np.random.randn(100) * 17
np.std(x)
# OUT: 15.986614289291284
y = np.convolve(x, np.ones(10)/10)
np.std(y)
# OUT: 5.244826663507928
np.std(y[10:-10])
# OUT: 5.280433690468031
np.std(x) / np.sqrt(10)
# OUT: 5.055411322875441
### 
