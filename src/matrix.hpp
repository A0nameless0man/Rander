#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <iostream>
#include <array>
#include <type_traits>
#include <cmath>

namespace Matrix
{
    //Method start with I means Inplace

    template <typename V = double>
    class UnSizedNumView
    {
    public:
        virtual V get(size_t y, size_t x) const = 0;
    };

    template <size_t H, size_t W, typename V = double>
    class NumView : public UnSizedNumView<V>
    {
    public:
        virtual size_t print(std::ostream &os) const
        {
            for (size_t i = 0; i < H; ++i)
            {
                for (size_t j = 0; j < W; ++j)
                {
                    os << this->get(i, j) << " ";
                }
                os << "\n";
            }
            return H * W;
        }
    };

    template <typename V = double>
    class UnSizedWritableNumView
    {
    public:
        virtual V &at(size_t y, size_t x) = 0;
        virtual UnSizedWritableNumView &set(size_t y, size_t x, V &v) = 0;
        virtual UnSizedWritableNumView &set(size_t y, size_t x, V &&v) = 0;
    };

    template <size_t H, size_t W, typename V = double>
    class NumWritableView : public UnSizedWritableNumView<V>
    {
    public:
    };
    template <size_t H, size_t W, typename V = double>
    class SubView : public NumView<H, W, V>
    {
    public:
        const UnSizedNumView<V> &ref;
        size_t offsetX;
        size_t offsetY;

        SubView(const UnSizedNumView<V> &ref, size_t offsetY, size_t offsetX) : ref(ref), offsetY(offsetY), offsetX(offsetX) {}
        virtual V get(size_t y, size_t x) const
        {
            return ref.get(y % H + offsetY, x % W + offsetX);
        }
    };

    template <size_t H, size_t W, typename V, typename SELF>
    class BasicNumArray : public NumView<H, W, V>, public NumWritableView<H, W, V>
    {
    public:
        V val[H][W];
        BasicNumArray() {}
        BasicNumArray(const std::array<V, H * W> &arr)
        {
            auto it = arr.cbegin();
            for (size_t i = 0; i < H; ++i)
            {
                for (size_t j = 0; j < W; ++j)
                {
                    val[i][j] = *it;
                    ++it;
                }
            }
        }
        virtual V get(size_t y, size_t x) const
        {
            return val[y % H][x % W];
        }
        virtual V &at(size_t y, size_t x)
        {
            return val[y % H][x % W];
        }
        virtual BasicNumArray &set(size_t y, size_t x, V &v)
        {
            val[y % H][x % W] = v;
            return *this;
        }
        virtual BasicNumArray &set(size_t y, size_t x, V &&v)
        {
            val[y % H][x % W] = v;
            return *this;
        }
        SELF operator/(const V &rop)
        {
            SELF res;
#pragma omp parallel for
            for (size_t i = 0; i < H; ++i)
            {
                for (size_t j = 0; j < W; ++j)
                {
                    res.set(i, j, this->get(i, j) / rop);
                }
            }
            return res;
        }
        BasicNumArray operator*(const V &rop)
        {
            SELF res;
#pragma omp parallel for
            for (size_t i = 0; i < H; ++i)
            {
                for (size_t j = 0; j < W; ++j)
                {
                    res.set(i, j, this->get(i, j) * rop);
                }
            }
            return res;
        }

        BasicNumArray &ISwapRow(size_t rop, size_t lop)
        {
            if (rop != lop)
            {
#pragma omp parallel for
                for (size_t i = 0; i < W; ++i)
                {
                    std::swap(at(rop, i), at(lop, i));
                }
            }
            return *this;
        }
        BasicNumArray &I_FMA_VRR(const V &x, size_t yr, size_t zr, size_t targetRow)
        {
#pragma omp parallel for
            for (size_t i = 0; i < W; ++i)
            {
                set(targetRow, i, std::fma(x, get(yr, i), get(zr, i)));
            }
            return *this;
        }

        size_t IGauss(const V &eps)
        {
            size_t swapCnt;
            this->print(std::cout);
            std::cout << std::endl;
            for (size_t i = 0; i < W; ++i)
            {
                size_t j = i;
                for (; j < W; ++j)
                {
                    if (std::abs(get(j, i)) > eps)
                    {
                        ISwapRow(i, j);
                        ++swapCnt;
                        break;
                    }
                }
                ++j;
                for (; j < W; ++j)
                {
                    if (std::abs(get(j, i)) > eps)
                    {
                        I_FMA_VRR(-(get(j, i) / get(i, i)), i, j, j);
                    }
                }
                this->print(std::cout);
                std::cout << std::endl;
            }
            return swapCnt;
        }
    };

    template <size_t H, size_t W, typename V>
    class NumArray : public BasicNumArray<H, W, V, NumArray<H, W, V>>
    {
        NumArray() {}
        NumArray(const std::array<V, W * W> &arr) : BasicNumArray<H, W, V, NumArray<H, W, V>>(arr) {}
    };

    template <size_t W, typename V>
    class NumArray<W, W, V> : public BasicNumArray<W, W, V, NumArray<W, W, V>>
    {
    public:
        NumArray() {}
        NumArray(const std::array<V, W * W> &arr) : BasicNumArray<W, W, V, NumArray<W, W, V>>(arr) {}
        V GetDeterminant(const V &eps) const
        {
            NumArray det(*this);
            V res = det.IGauss(eps) & 1 ? -1 : 1;
            for (size_t i = 0; i < W; ++i)
            {
                res *= det.at(i, i);
            }
            return res;
        }
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
            for (size_t i = 0; i < H; ++i)
            {
                for (size_t j = 0; j < W; ++j)
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
        Matrix<W, W, V> inv() const
        {
        }
    };

}

#endif