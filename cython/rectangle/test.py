#!/usr/bin/env python3
import numpy as np
import rect

x0, y0, x1, y1 = 1, 2, 3, 4
rect_obj = rect.PyRectangle(x0, y0, x1, y1)

print(rect_obj.x0)
print(rect_obj.get_size())
print(rect_obj.get_area())
print(rect_obj.get_start())
print("Moving by (1, 2)")
rect_obj.move(1, 2)
print(rect_obj.get_start())

a = np.arange(4) + 20
rect_obj.set_integer_arr_ptr(a)

print(a)
