#include <cstddef>
#include <iostream>
#include <array>
namespace Matrix
{

    template <size_t H, size_t W, typename V = double>
    class NumView
    {
    public:
        virtual V get(size_t y, size_t x) const = 0;
        virtual V set(size_t y, size_t x, V &v) = 0;
        virtual V set(size_t y, size_t x, V &&v) = 0;
        template <size_t VH, size_t VW>
        class View
        {
            virtual const NumView<VH, VW, V> getSubView(size_t offsetY, size_t offsetX) const = 0;
        };
        virtual size_t print(std::ostream &os) const;
    };

    template <size_t H, size_t W, typename V = double>
    class NumArray : NumView<H, W, V>
    {
    protected:
        V val[H][W];

    public:
        NumArray(const std::array<V, H * W> &);
        template <size_t VH, size_t VW>
        class ArrayView : NumView<H, W, V>
        {
        protected:
            size_t offsetW;
            size_t offsetH;
        };
    };

    template <size_t H, size_t W, typename V = double>
    class Matrix : NumArray<H, W, V>
    {
    private:
    public:
        Matrix();
        ~Matrix();
    };

    template <size_t W, typename V>
    class Matrix<W, W, V>
    {
    public:
        Matrix<W, W, V> inv() const;
    };

}