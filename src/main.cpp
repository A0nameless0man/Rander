#include <iostream>
#include "matrix.hpp"

int main(void)
{
    std::cout << "Hello" << std::endl;
    using Target = double;

    Matrix::NumArray<4, 4, Target> a = {{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}}};
    // a.ISwapRow(1, 3);
    a.print(std::cout);
    std::cout<<std::endl;
    a.IGauss(0.01);
    // a.I_FMA_VRR(-5,0,1,1);
    a.print(std::cout);
    std::cout<<std::endl;
    Matrix::Matrix<4, 4, Target> m = {{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}}};
    // m=m/2;
    m.print(std::cout);
    std::cout<<std::endl;
    Matrix::SubView<3, 3, Target> v(m, 3, 3);
    Matrix::SubView<3, 3, Target> v2(v, 2, 2);
    v.print(std::cout);
    std::cout<<std::endl;
    v2.print(std::cout);
    // BaseMatrix::NumView<4, 4, int> &v = m;
    Matrix::NumArray<3,3,Target> e = {{{2,1,-1,-3,-1,2,-2,1,2}}};
    // e.IGauss(0.01);
    std::cout<<e.GetDeterminant(0.01);
    std::cout<<std::endl;
    e.print(std::cout);
}