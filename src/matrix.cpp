#include <iostream>

#include "matrix.hpp"

namespace Matrix
{

    template <size_t H, size_t W, typename V>
    size_t NumView<H, W, V>::print(std::ostream &os) const
    {
        for (size_t i = 0; i < H; i++)
        {
            for (size_t j = 0; j < W; j++)
            {
                os << this->get(i, j) << " ";
            }
            os << "\n";
        }
        return H * W;
    }



    template <size_t W, typename V>
    Matrix<W, W, V> Matrix<W, W, V>::inv() const
    {
        return this;
    }

}