// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

extern "C" {
#include <sys/mman.h>
#include <sys/types.h>
}

template <typename T>
class TlMmapAllocater {
 public:
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef T* pointer;
  typedef const T* const_pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T value_type;
  template <class U>
  struct rebind {
    typedef TlMmapAllocater<U> other;
  };

 public:
  TlMmapAllocater() throw() {}

  TlMmapAllocater(const TlMmapAllocater&) throw() {}

  template <class U>
  TlMmapAllocater(const TlMmapAllocater<U>&) throw() {}

  ~TlMmapAllocater() {}

 public:
  pointer address(reference x) const { return &x; }

  const_pointer address(const_reference x) const { return &x; }

  pointer allocate(size_type n, void* hint = 0) {
    // return static_cast<T*>(::operator new(n * sizeof(T)));
    pointer p = static_cast<pointer>(
        ::mmap(hint, n, PROT_READ | PROT_WRITE, MAP_ANONYMOUS, -1, 0));
    if ((p != NULL) && ((int)*p == -1)) {
      p = NULL;
    }
    return p;
  }

  void deallocate(pointer p, size_type n) {
    ::operator delete(static_cast<void*>(p));
  }

  size_type max_size() const throw() {
    return std::numeric_limits<size_type>::max() / sizeof(T);
  }

  void construct(pointer p, const T& val) {
    new (static_cast<void*>(p)) T(val);
  }

  void destroy(pointer p) { p->~T(); }
};

template <typename T>
bool operator==(const TlMmapAllocater<T>&, const TlMmapAllocater<T>&) {
  return true;
}

template <typename T>
bool operator!=(const TlMmapAllocater<T>&, const TlMmapAllocater<T>&) {
  return false;
}

template <>
class TlMmapAllocater<void> {
 public:
  typedef void* pointer;
  typedef const void* const_pointer;
  typedef void value_type;
  template <class U>
  struct rebind {
    typedef TlMmapAllocater<U> other;
  };
};
