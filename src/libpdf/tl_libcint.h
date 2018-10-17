#ifndef TLLIBCINT_H
#define TLLIBCINT_H

#include <vector>

class Fl_Geometry;
class TlOrbitalInfo;

class TlLibcint
{
  public:
    TlLibcint(const Fl_Geometry &geom, const TlOrbitalInfo &orbInfo);
    ~TlLibcint();

  public:
    void calc_ovp(const int shell_id_i, const int shell_id_j);
    void calc_nuc(const int shell_id_i, const int shell_id_j);
    void calc_kin(const int shell_id_i, const int shell_id_j);
    void calc_ERI(const int shell_id_i, const int shell_id_j,
                  const int shell_id_k, const int shell_id_l);

  protected:
    void makeAtomTable(const Fl_Geometry &geom);
    void makeBasissetTable(const TlOrbitalInfo &orbInfo);

    int getValuesOffset() const;
    void addValue(const double value);

  private:
    static const double AU2ANG;
    static const int ATOM_BLOCK_SIZE;
    static const int BASISSET_BLOCK_SIZE;

    const Fl_Geometry &geom_;
    const TlOrbitalInfo &orbInfo_;

    std::vector<int> atomTable_;
    std::vector<int> basissetTable_;

    int valuesOffset_;
    std::vector<double> values_;
};

#endif // TLLIBCCINT_H
