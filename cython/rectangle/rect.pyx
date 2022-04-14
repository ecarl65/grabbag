# distutils: language = c++

from Rectangle cimport Rectangle

# Create a Cython extension type which holds a C++ instance
# as an attribute and create a bunch of forwarding methods
# Python extension type.
cdef class PyRectangle:
    cdef Rectangle c_rect  # Hold a C++ instance which we're wrapping

    def __cinit__(self, int x0, int y0, int x1, int y1):
        self.c_rect = Rectangle(x0, y0, x1, y1)

    def get_area(self):
        return self.c_rect.getArea()

    def get_size(self):
        cdef int width, height
        self.c_rect.getSize(&width, &height)
        return width, height

    def move(self, dx, dy):
        self.c_rect.move(dx, dy)

    def get_start(self):
        cdef int x_0, y_0
        x_0 = self.c_rect.getX0()
        y_0 = self.c_rect.getY0()
        return x_0, y_0


def main():
   rec_ptr = new Rectangle(1, 2, 3, 4) #  Instantiate a Rectangle object on the heap
   try:
        rec_area = rec_ptr.getArea()
   finally:
        del rec_ptr #  delete heap allocated object

   cdef Rectangle rec_stack # Instantiate a Rectangle object on the stack
