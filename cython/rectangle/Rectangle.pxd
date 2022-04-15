cdef extern from "Rectangle.cpp":
    pass

# Decalre the class with cdef
cdef extern from "Rectangle.h" namespace "shapes":
    cdef cppclass Rectangle:
        Rectangle() except +
        Rectangle(int, int, int, int) except +
        int x0, y0, x1, y1
        int getArea()
        void getSize(int* width, int* height)
        int getX0()
        int getX1()
        int getY0()
        int getY1()
        void move(int, int)
        void set_integer_arr_ptr(int* a)
