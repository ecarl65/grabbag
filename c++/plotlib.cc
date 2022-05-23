/*
 * =====================================================================================
 *
 *       Filename:  plotlib.cc
 *
 *    Description:  Testing the gdbplotlib
 *
 *        Version:  1.0
 *        Created:  05/17/2022 09:49:53 AM
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
#include <array>

int main()
{
    std::array<double, 6> x = {0.1, 0.9, 0.8, 0.7, 0.2, 0.1};

    int* y = new int[100];
    for (int i = 0; i < 100; ++i) {
        y[i] = 50 - i + int(5e-3 * i * i);
    }

    std::vector<std::array<int*, 10>> z(10);
    for (int i = 0; i < z.size(); ++i) {
        for (int j = 0; j < z[i].size(); ++j) {
            z[i][j] = new int[10];
            for (int k = 0; k < 10; ++k) {
                z[i][j][k] = i + 2*j + 3*k;
            }
        }
    }

    return 0;
}
