#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <iostream>
#include <array>
#include <type_traits>

namespace Matrix
{

    template <size_t H, size_t W, typename V = double>
    class NumView
    {
    public:
        virtual V get(size_t y, size_t x) const = 0;
        virtual NumView &set(size_t y, size_t x, V &v) = 0;
        virtual NumView &set(size_t y, size_t x, V &&v) = 0;
        template <size_t VH, size_t VW>
        class View
        {
            virtual const NumView<VH, VW, V> getSubView(size_t offsetY, size_t offsetX) const = 0;
        };
        virtual size_t print(std::ostream &os) const;
    };

    template <size_t H, size_t W, typename V = double>
    class NumArray : public NumView<H, W, V>
    {
    public:
        V val[H][W];
        NumArray() {}
        NumArray(const std::array<V, H * W> &);
        virtual V get(size_t y, size_t x) const;
        virtual NumArray &set(size_t y, size_t x, V &v);
        virtual NumArray &set(size_t y, size_t x, V &&v);
        NumArray operator/(const V &rop)
        {
            NumArray res;
#pragma omp parallel for
            for (size_t i = 0; i < H; i++)
            {
                for (size_t j = 0; j < W; j++)
                {
                    res.set(i, j, this->get(i, j) / rop);
                }
            }
            return res;
        }
        template <size_t VH, size_t VW>
        class ArrayView : public NumView<H, W, V>
        {
        public:
            size_t offsetW;
            size_t offsetH;
        };
    };

    template <size_t H, size_t W, typename V, typename SELF>
    class BaseMatrix : public NumArray<H, W, V>
    {
    public:
        BaseMatrix() {}
        BaseMatrix(const std::array<V, H * W> &arr) : NumArray<H, W, V>(arr) {}
        SELF operator/(const V &rop)
        {
            SELF res;
#pragma omp parallel for
            for (size_t i = 0; i < H; i++)
            {
                for (size_t j = 0; j < W; j++)
                {
                    res.set(i, j, this->get(i, j) / rop);
                }
            }
            return res;
        }
    };
    template <size_t H, size_t W, typename V = double>
    class Matrix : public BaseMatrix<H, W, V, Matrix<H, W, V>>
    {
    public:
        Matrix() {}
        Matrix(const std::array<V, H * W> &arr) : BaseMatrix<H, W, V, Matrix<H, W, V>>(arr) {}
    };

    template <size_t W, typename V>
    class Matrix<W, W, V> : public BaseMatrix<W, W, V, Matrix<W, W, V>>
    {
    public:
        Matrix() {}
        Matrix(const std::array<V, W * W> &arr) : BaseMatrix<W, W, V, Matrix<W, W, V>>(arr) {}
        Matrix<W, W, V> inv() const;
    };

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

    template <size_t H, size_t W, typename V>
    NumArray<H, W, V>::NumArray(const std::array<V, H * W> &arr)
    {
        auto it = arr.cbegin();
        for (size_t i = 0; i < H; i++)
        {
            for (size_t j = 0; j < W; j++)
            {
                val[i][j] = *it;
                ++it;
            }
        }
    }

    template <size_t H, size_t W, typename V>
    V NumArray<H, W, V>::get(size_t y, size_t x) const
    {
        return val[y][x];
    }

    template <size_t H, size_t W, typename V>
    NumArray<H, W, V> &NumArray<H, W, V>::set(size_t y, size_t x, V &v)
    {
        val[y][x] = v;
        return *this;
    }

    template <size_t H, size_t W, typename V>
    NumArray<H, W, V> &NumArray<H, W, V>::set(size_t y, size_t x, V &&v)
    {
        val[y][x] = v;
        return *this;
    }

    template <size_t W, typename V>
    Matrix<W, W, V> Matrix<W, W, V>::inv() const
    {
        return this;
    }
}

#endif