/*
 * =====================================================================================
 *
 *       Filename:  Rectangle.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/14/2022 02:41:06 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef RECTANGLE_H
#define RECTANGLE_H

namespace shapes {
  class Rectangle {
    public:
      int x0, y0, x1, y1;
      Rectangle();
      Rectangle(int x0, int y0, int x1, int y1);
      ~Rectangle();
      int getArea();
      void getSize(int* width, int* height);
      void move(int dx, int dy);
      int getX0();
      int getX1();
      int getY0();
      int getY1();
      void set_integer_arr_ptr(int* a);
  };
}

#endif
