#ifndef NULLPTR_H

#if (__cplusplus < 201103L)

const  // const オブジェクトであって、
    class nullptr_t {
 public:
  template <class T>
  operator T*() const  // 任意の非メンバ型のヌルポインタや、
  {
    return 0;
  }

  template <class C, class T>
  operator T C::*() const  // 任意のメンバ型のヌルポインタに変換可能であって、
  {
    return 0;
  }

 private:
  void operator&() const;  // アドレスを取得することができない。

} nullptr = {};

#endif  // (__cplusplus < 201103L)

#endif  // NULLPTR_H
