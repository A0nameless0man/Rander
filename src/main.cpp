#include <iostream>
#include "matrix.hpp"

int main(void)
{
    std::cout << "Hello" << std::endl;
    Matrix::NumArray<4, 4, int> a = {{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}}};
    a=a/2;
    a.print(std::cout);
    Matrix::Matrix<4, 4, int> m = {{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}}};
    m=m/2;
    m.print(std::cout);
    // BaseMatrix::NumView<4, 4, int> &v = m;
}