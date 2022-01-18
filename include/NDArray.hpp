///
/// Created by F.Moitzi on 09.01.2022.
///
/// The memory layout of this class is designed to be the same in {\bf Fortran}.
///

#ifndef MULTIARRAY_HPP
#define MULTIARRAY_HPP

#include <cstddef>
#include <stdexcept>
#include <vector>
#include <memory>
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <array>

namespace detail {
  template<bool...>
  struct bool_pack;
  template<bool... bs>
  using all_true = std::is_same<bool_pack<bs..., true>, bool_pack<true, bs...>>;
}

template<typename... Ts>
using all_true = detail::all_true<Ts::value...>;

template<typename T, typename... Ts>
using all_same = all_true<std::is_same<T, Ts>...>;


namespace lsms {


  template<typename std::size_t Size>
  using array_t = std::array<std::size_t, Size>;

  template<typename T>
  using allocator_t = std::allocator<T>;

  using size_t = std::size_t;

  template<class T>
  class ArrayItem {

    using iterator = typename std::vector<ArrayItem<T>>::iterator;
    using const_iterator = typename std::vector<ArrayItem<T>>::const_iterator;

  private:

    bool is_array_;
    std::vector<ArrayItem<T>> v_;
    T val_;

  public:


    explicit ArrayItem(const std::vector<ArrayItem<T>> &a)
        : is_array_{true}, v_(a), val_{0} {
    }

    ArrayItem(std::initializer_list<ArrayItem<T>> list)
        : is_array_{true}, v_(list), val_{0} {
    }

    ArrayItem(T val) // Cannot be explicit
        : is_array_{false}, val_{val} {
    }

    ArrayItem() : is_array_{false} {};

    ArrayItem(const ArrayItem &) = default;

    ArrayItem(ArrayItem &&) noexcept = default;

    bool is_array() const { return is_array_; }

    size_t size() const { return v_.size(); }

    T value() const { return val_; }

    iterator begin() { return v_.begin(); }

    iterator end() { return v_.end(); }

    const_iterator begin() const { return v_.begin(); };

    const_iterator end() const { return v_.end(); };


  };


  class ColumnMajor {

  public:

    template<size_t N>
    static void calculate_strides(const array_t<N> &shape, array_t<N> &strides) {
      size_t size = 1;
      for (size_t i = 0; i < N; ++i) {
        strides[i] = size;
        size *= shape[i];
      }
    }

  };

  template<typename Pointer>
  inline typename std::pointer_traits<Pointer>::element_type *to_plain_pointer(Pointer ptr) {
    return (std::addressof(*ptr));
  }


  template<size_t N, typename... Indices>
  inline size_t get_offset(array_t<N> strides, Indices... indices) {

    const array_t<N> inds{indices...};

    size_t off{0};
    for (auto i = 0; i < N; i++) {
      off += strides[i] * inds[i];
    }

    return off;
  }


  template<typename T, size_t N>
  class NDArray {

  private:

    T *data_ = nullptr; // null pointer literal

    size_t size_{0}; // total number of elements

    array_t<N> shape_; // shape of the n-dims array (can be change)

    array_t<N> strides_; //

    allocator_t<T> allocator_; // allocator

    void dim_from_initializer_list(const ArrayItem<T> &init, size_t shape) {
      bool is_array = false;
      size_t size = 0;

      size_t i = 0;
      for (const auto &item : init) {

        if (i == 0) {
          is_array = item.is_array();
          size = item.size();
          if (shape < N) {
            shape_[shape++] = init.size();
            if (is_array) {
              dim_from_initializer_list(item, shape);
            }
          }
        } else {

          if (is_array) {

            if (!item.is_array() || item.size() != size) {

              throw std::invalid_argument("initializer list contains non-conforming shapes");

            }

          } else if (shape != N) {

            throw std::invalid_argument("initializer list incompatible with array dimensionality");

          }
        }

        ++i;
      }
    }


    void data_from_initializer_list(const ArrayItem<T> &init, array_t<N> &indices, size_t index) {
      size_t i = 0;
      for (const auto &item : init) {
        indices[index] = i;

        if (item.is_array()) {

          data_from_initializer_list(item, indices, index + 1);

        } else {

          size_t offset = get_offset<N>(strides_, indices);
          if (offset < size()) {
            data_[offset] = item.value();
          }

        }

        ++i;
      }
    }


  public:


    NDArray() {
      std::fill(std::begin(shape_), std::end(shape_), 0);
      std::fill(std::begin(strides_), std::end(strides_), 0);
    };

    template<typename... Args>
    explicit NDArray(size_t i, Args... args) : shape_{i, args...} {
      static_assert(sizeof...(args) != N, "Number of arguments are wrong!");
      size_ = std::accumulate(std::begin(shape_), std::end(shape_), 1, std::multiplies<decltype(N)>());
      ColumnMajor::calculate_strides<N>(shape_, strides_);
      data_ = allocator_.allocate(size_);
      T def{};
      std::fill_n(data_, size_, def);
    };

    ~NDArray() {
      allocator_.deallocate(data_, size_);
    }

    NDArray(const NDArray &other)
        : data_(nullptr),
          size_(other.size()),
          shape_{other.shape_},
          strides_{other.strides_} {
      data_ = allocator_.allocate(size_);
      std::copy(other.data_, other.data_ + other.size_, data_);
    }

    NDArray(NDArray &&other) noexcept
        : data_(std::exchange(other.data_, nullptr)),
          size_(std::move(other.size())),
          shape_(std::move(other.shape_)),
          strides_(std::move(other.strides_)) {
    };

    friend void assign(NDArray &lhs, const NDArray &rhs) {
      lhs.size_ = rhs.size_;
      lhs.shape_ = rhs.shape_;
      lhs.strides_ = rhs.strides_;
    }

    NDArray &operator=(const NDArray &other) {

      if (this != &other) {
        assign(*this, other);
        data_ = allocator_.allocate(other.size_);
        std::copy(other.data_, other.data_ + other.size_, data_);
      }
      return *this;
    }


    NDArray &operator=(NDArray &&other) noexcept {
      if (data_) allocator_.deallocate(data_, size_);
      data_ = std::exchange(other.data_, nullptr);
      size_ = std::move(other.size());
      shape_ = std::move(other.shape_);
      strides_ = std::move(other.strides_);
      return *this;
    };


    // get access operator
    template<typename... Indices>
    inline T &operator()(size_t i, Indices... indices) {
      static_assert(sizeof...(indices) != N, "Number of arguments are wrong!");
      return data_[get_offset<N>(strides_, i, indices...)];
    }

    template<typename... Indices>
    inline const T &operator()(size_t i, Indices... indices) const {
      static_assert(sizeof...(indices) != N, "Number of arguments are wrong!");
      return data_[get_offset<N>(strides_, i, indices...)];
    }

    inline NDArray &operator=(const T &val) {
      std::fill(data_, data_ + size_, val);
      return *this;
    }

    void scale(const T &val) {
      std::transform(data_, data_ + size_, data_, [val](T &c) { return c * val; });
    }

    inline NDArray &operator*=(const T &val) {
      std::transform(data_, data_ + size_, data_, [val](T &c) { return c * val; });
      return *this;
    }

    inline NDArray &operator/=(const T &val) {
      std::transform(data_, data_ + size_, data_, [val](T &c) { return c / val; });
      return *this;
    }


    template<typename Indices>
    inline T &operator[](Indices i) {
      return data_[i];
    }

    std::size_t size() const {
      return size_;
    }

    const array_t<N> &shape() const { return shape_; }

    const array_t<N> &strides() const { return strides_; }

    T *data() {
      return to_plain_pointer(data_);
    }

    const T *data() const {
      return to_plain_pointer(data_);
    }

    size_t shape(size_t i) const { return shape_[i]; }


    NDArray(std::initializer_list<ArrayItem<T>> list) {

      dim_from_initializer_list(list, 0);

      ColumnMajor::calculate_strides<N>(shape_, strides_);

      size_ = std::accumulate(std::begin(shape_), std::end(shape_), 1, std::multiplies<decltype(N)>());

      data_ = allocator_.allocate(size_);

      data_from_initializer_list(list, shape_, 0);

    }


  };


}

#endif //MULTIARRAY_HPP
