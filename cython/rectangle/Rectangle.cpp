/*
 * =====================================================================================
 *
 *       Filename:  Rectangle.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/14/2022 02:41:45 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <vector>
#include <algorithm>

#include <iostream>
#include "Rectangle.h"

namespace shapes {

    // Default constructor
    Rectangle::Rectangle () {}

    // Overloaded constructor
    Rectangle::Rectangle (int x0, int y0, int x1, int y1) {
        this->x0 = x0;
        this->y0 = y0;
        this->x1 = x1;
        this->y1 = y1;
    }

    // Destructor
    Rectangle::~Rectangle () {}

    // Return the area of the rectangle
    int Rectangle::getArea () {
        return (this->x1 - this->x0) * (this->y1 - this->y0);
    }

    // Get the size of the rectangle.
    // Put the size in the pointer args
    void Rectangle::getSize (int *width, int *height) {
        (*width) = x1 - x0;
        (*height) = y1 - y0;
    }

    int Rectangle::getX0() {
      return this->x0;
    }
    int Rectangle::getX1() {
      return this->x1;
    }
    int Rectangle::getY0() {
      return this->y0;
    }
    int Rectangle::getY1() {
      return this->y1;
    }

    // Move the rectangle by dx dy
    void Rectangle::move (int dx, int dy) {
        this->x0 += dx;
        this->y0 += dy;
        this->x1 += dx;
        this->y1 += dy;
    }

    void Rectangle::set_integer_arr_ptr(int* a) {
      std::vector<int> d(4);
      d[0] = 1;
      d[1] = 3;
      d[2] = 4;
      d[3] = 7;

      std::copy(d.begin(), d.end(), a);
    }
}
