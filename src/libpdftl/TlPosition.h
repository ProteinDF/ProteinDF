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

#ifndef TLPOSITION_H
#define TLPOSITION_H

#include <cassert>
#include <cmath>
#include "tl_dense_general_matrix_lapack.h"
#include "tl_dense_vector_lapack.h"

/// 3次元座標(カーテシアン)を扱うクラス
class TlPosition {
 public:
  /// コンストラクタ
  ///
  /// @param[in] x X座標
  /// @param[in] y Y座標
  /// @param[in] z Z座標
  explicit TlPosition(const double x = 0.0, const double y = 0.0,
                      const double z = 0.0);

  /// コピーコンストラクタ
  TlPosition(const TlPosition& p);

  /// デストラクタ
  ~TlPosition();

 public:
  /// X座標の値を返す
  double x() const;

  /// Y座標の値を返す
  double y() const;

  /// Z座標の値を返す
  double z() const;

  // 座標ごとに設定
  void x(double nx);
  void y(double ny);
  void z(double nz);

  ///
  const double& operator[](int index) const;

  double& operator[](int index);

  /// 原点に移動する
  void clear();

  /// 指定された座標に移動する
  ///
  /// @param[in] x X座標
  /// @param[in] y Y座標
  /// @param[in] z Z座標
  TlPosition& moveTo(const double x = 0.0, const double y = 0.0,
                     const double z = 0.0);

  /// 指定された座標に移動する
  ///
  /// @param[in] position TlPositionオブジェクト
  TlPosition& moveTo(const TlPosition& position);

  /// 指定された座標分、移動する
  ///
  /// @param[in] x X座標
  /// @param[in] y Y座標
  /// @param[in] z Z座標
  TlPosition& shiftBy(const double x = 0.0, const double y = 0.0,
                      const double z = 0.0);

  /// 指定された座標分、移動する
  ///
  /// @param[in] position TlPositionオブジェクト
  TlPosition& shiftBy(const TlPosition& position);

  /// 回転
  ///
  /// rot must be a 3x3 rotation matrix
  /// @param[in] rot 回転行列
  TlPosition& rotateBy(const TlDenseGeneralMatrix_Lapack& rot);

  // 原点からの距離を返す
  double distanceFrom() const;

  /// 指定された座標からの距離を返す
  ///
  /// @param[in] x X座標
  /// @param[in] y Y座標
  /// @param[in] z Z座標
  double distanceFrom(const double x, const double y, const double z) const;

  /// 指定された座標からの距離を返す
  ///
  /// @param[in] position TlPositionオブジェクト
  double distanceFrom(const TlPosition& position) const;

  /// 原点からの距離の二乗を返す
  double squareDistanceFrom() const;

  /// 指定された座標からの距離の二乗を返す
  ///
  /// @param[in] x X座標
  /// @param[in] y Y座標
  /// @param[in] z Z座標
  double squareDistanceFrom(const double x, const double y,
                            const double z) const;

  /// 指定された座標からの距離の二乗を返す
  ///
  /// @param[in] position TlPositionオブジェクト
  double squareDistanceFrom(const TlPosition& position) const;

  /// 単位ベクトル(原点からの距離を1)にする
  void unit();

  /// 指定された点から指定範囲内にあるかどうかを返す
  ///
  /// @param[in] pos 指定座標
  /// @param[in] range 距離
  /// @retval true 範囲内にある
  /// @retval false 範囲内に無い
  bool isWithinRange(const TlPosition& pos, double range) const;

  /// 代入演算子
  TlPosition& operator=(const TlPosition& rhs);

  /// 座標ベクトルの加算
  TlPosition& operator+=(const TlPosition& rhs);

  /// 座標ベクトルの減算
  TlPosition& operator-=(const TlPosition& rhs);

  /// 座標ベクトルの定数倍
  TlPosition& operator*=(const double a);

  /// 座標ベクトルの定数除算
  TlPosition& operator/=(const double a);

  bool operator==(const TlPosition& rhs) const;
  bool operator!=(const TlPosition& rhs) const {
    return !(this->operator==(rhs));
  }

 protected:
  static const double TOO_SMALL;
  static const double EQUIV_POSITION_RANGE;

  double v_[3];

  mutable bool m_bHasDistance;  /// m_distanceを計算しているかどうか
  mutable double m_distance;    /// 原点からの距離を保存

  friend TlPosition operator+(const TlPosition& rhs1, const TlPosition& rhs2);
  friend TlPosition operator-(const TlPosition& rhs1, const TlPosition& rhs2);
  friend TlPosition operator*(const TlPosition& rhs1, const double& rhs2);

  friend TlPosition operator*(const double& rhs1, const TlPosition& rhs2) {
    return (rhs2 * rhs1);
  };

  friend TlPosition operator/(const TlPosition& rhs1, const double& rhs2) {
    assert(std::fabs(rhs2) > TOO_SMALL);
    return (rhs1 * (1.0 / rhs2));
  };

  friend TlPosition operator*(const TlDenseGeneralMatrix_Lapack& rot,
                              const TlPosition& pos) {
    assert(rot.getNumOfRows() == 3);
    assert(rot.getNumOfCols() == 3);

    TlDenseVector_Lapack x(3);
    x.set(0, pos.v_[0]);
    x.set(1, pos.v_[1]);
    x.set(2, pos.v_[2]);
    x = rot * x;

    return TlPosition(x.get(0), x.get(1), x.get(2));
  };
};

// inline ======================================================================
inline const double& TlPosition::operator[](int index) const {
  assert((0 <= index) && (index < 3));
  return this->v_[index];
}

inline double& TlPosition::operator[](int index) {
  assert((0 <= index) && (index < 3));
  return this->v_[index];
}

inline double TlPosition::x() const { return this->v_[0]; }

inline double TlPosition::y() const { return this->v_[1]; }

inline double TlPosition::z() const { return this->v_[2]; }

inline void TlPosition::x(double nx) {
  this->v_[0] = nx;
  this->m_bHasDistance = false;
}

inline void TlPosition::y(double ny) {
  this->v_[1] = ny;
  this->m_bHasDistance = false;
}

inline void TlPosition::z(double nz) {
  this->v_[2] = nz;
  this->m_bHasDistance = false;
}

inline TlPosition& TlPosition::moveTo(const double x, const double y,
                                      const double z) {
  this->v_[0] = x;
  this->v_[1] = y;
  this->v_[2] = z;
  this->m_bHasDistance = false;

  return *this;
}

inline TlPosition& TlPosition::moveTo(const TlPosition& p) {
  return this->moveTo(p.x(), p.y(), p.z());
}

inline TlPosition& TlPosition::shiftBy(const double x, const double y,
                                       const double z) {
  this->v_[0] += x;
  this->v_[1] += y;
  this->v_[2] += z;
  this->m_bHasDistance = false;

  return *this;
}

inline TlPosition& TlPosition::shiftBy(const TlPosition& p) {
  return this->shiftBy(p.x(), p.y(), p.z());
}

inline double TlPosition::distanceFrom() const {
  return std::sqrt(this->squareDistanceFrom());
}

inline double TlPosition::distanceFrom(const TlPosition& p) const {
  return this->distanceFrom(p.x(), p.y(), p.z());
}

inline double TlPosition::squareDistanceFrom() const {
  if (this->m_bHasDistance == false) {
    const double dx = this->v_[0];
    const double dy = this->v_[1];
    const double dz = this->v_[2];
    this->m_distance = dx * dx + dy * dy + dz * dz;
    this->m_bHasDistance = true;
  }

  return this->m_distance;
}

inline double TlPosition::distanceFrom(const double x, const double y,
                                       const double z) const {
  return std::sqrt(this->squareDistanceFrom(x, y, z));
}

inline double TlPosition::squareDistanceFrom(const double x, const double y,
                                             const double z) const {
  const double dx = this->x() - x;
  const double dy = this->y() - y;
  const double dz = this->z() - z;
  return (dx * dx + dy * dy + dz * dz);
}

inline double TlPosition::squareDistanceFrom(const TlPosition& p) const {
  return this->squareDistanceFrom(p.x(), p.y(), p.z());
}

inline TlPosition operator+(const TlPosition& rhs1, const TlPosition& rhs2) {
  TlPosition p(rhs1);
  p.shiftBy(rhs2);

  return p;
}

inline TlPosition operator-(const TlPosition& rhs1, const TlPosition& rhs2) {
  TlPosition p(rhs1);
  p.shiftBy(-rhs2.x(), -rhs2.y(), -rhs2.z());

  return p;
}

inline TlPosition operator*(const TlPosition& rhs1, const double& rhs2) {
  TlPosition p(rhs1);
  p *= rhs2;

  return p;
}

#endif  // TLPOSITION_H
