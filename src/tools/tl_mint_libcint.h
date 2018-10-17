#ifndef TL_MINT_LIBCINT_H
#define TL_MINT_LIBCINT_H

#include <vector>
#include "tl_mint_object.h"

class TlMint_Libcint : public TlMintObject {
 public:
  TlMint_Libcint(const Fl_Geometry& geom, const TlOrbitalInfo& orbInfo);
  virtual ~TlMint_Libcint();

 public:
  virtual void calc_nuc(int shell_p, int shell_q);
  virtual void calc_kin(int shell_p, int shell_q);
  virtual void calc_ovp(int shell_p, int shell_q);
  virtual void calc_eri(int shell_p, int shell_q, int shell_r, int shell_s);

public:
  void showAtomTable() const;

 protected:
  void makeAtomTable(const Fl_Geometry& geom);
  void makeBasissetTable(const TlOrbitalInfo& orbInfo);

  int getValuesOffset() const;
  void addValue(const double value);

 private:
  static const double AU2ANG;
  static const int ATOM_BLOCK_SIZE;
  static const int BASISSET_BLOCK_SIZE;

  std::vector<int> atomTable_;
  std::vector<int> basissetTable_;

  int valuesOffset_;
  std::vector<double> values_;
};

#endif  // TL_MINT_LIBCINT_H
