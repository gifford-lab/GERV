// Copyright 2014 Daniel Kang

#include <algorithm>

#include "../consts.h"

#ifndef CLASSES_ARRAYF_H_
#define CLASSES_ARRAYF_H_

template <class T, size_t N, size_t offset>
class Array {
  T *arr;
  const size_t SIZE;
  const size_t OFFSET;

  typedef T& reference;
  typedef const T& const_reference;

 public:
  explicit Array() :
      arr(((T *) calloc(N, sizeof(T))) + offset), SIZE(N), OFFSET(offset) {}

  reference operator[] (size_t n) { return *(arr + n); }

  const_reference operator[](size_t n) const { return *(arr + n); }

  void clear() { memset(arr - OFFSET, 0, SIZE * sizeof(T)); }

  void swap(Array& other) {
    assert(other.SIZE == SIZE);
    std::swap(arr, other.arr);
  }

  void copy(const Array& other) {
    assert(other.SIZE == SIZE);
    memcpy(arr - OFFSET, other.arr - OFFSET, SIZE * sizeof(T));
  }

  T* data() { return arr; }

  constexpr size_t size() { return SIZE; }

  ~Array() { free(arr - OFFSET); }

  Array(const Array&) = delete;
  Array& operator=(const Array&) = delete;
};

template <size_t N, size_t offset>
using ArrayF = Array<ftype, N, offset>;

#endif  // CLASSES_ARRAYF_H_
