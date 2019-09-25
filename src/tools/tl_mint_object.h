#ifndef TL_MINT_OBJECT_H
#define TL_MINT_OBJECT_H

class Fl_Geometry;
class TlOrbitalInfo;

class TlMintObject {
   public:
    TlMintObject(const Fl_Geometry& geom, const TlOrbitalInfo& orbInfo);
    virtual ~TlMintObject();

   public:
    bool checkAO(int aoIndex) const;

    virtual void calc_ovp(int shell_p, int shell_q) = 0;
    virtual void calc_nuc(int shell_p, int shell_q) = 0;
    virtual void calc_kin(int shell_p, int shell_q) = 0;
    virtual void calc_eri(int shell_p, int shell_q, int shell_r,
                          int shell_s) = 0;

   protected:
    const Fl_Geometry& geom_;
    const TlOrbitalInfo& orbInfo_;
};

#endif  // TL_MINT_OBJECT_H
