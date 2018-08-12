#ifndef TL_VECTOR_OBJECT_H
#define TL_VECTOR_OBJECT_H

class TlVectorObject {
  // ---------------------------------------------------------------------------
 public:
  typedef signed int index_type;
  typedef signed int size_type;

  // ---------------------------------------------------------------------------
 public:
  /// vector要素を格納するための構造体
  ///
  /// TlCommunicate で通信可能
  struct VectorElement {
    typedef TlVectorObject::index_type index_type;

   public:
    VectorElement(index_type i = 0, double v = 0.0)
        : index(i), value(v) {}

   public:
    index_type index;
    double value;
  };
};

#endif // TL_VECTOR_OBJECT_H
