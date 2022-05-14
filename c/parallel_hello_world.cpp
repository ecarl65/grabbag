/*
 * =====================================================================================
 *
 *       Filename:  openmp_example.c
 *
 *    Description: Openmp Example 
 *
 *        Version:  1.0
 *        Created:  05/11/2022 11:32:16 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv) {
  printf("Hello from process: %d\n", omp_get_thread_num());

  return 0;
}

