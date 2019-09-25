#ifndef TL_MINT_PDF_H
#define TL_MINT_PDF_H

#include <vector>
#include "DfHpqEngine.h"
#include "tl_mint_object.h"

class TlMint_Pdf : public TlMintObject {
   public:
    TlMint_Pdf(const Fl_Geometry& geom, const TlOrbitalInfo& orbInfo);
    virtual ~TlMint_Pdf();

   public:
    virtual void calc_nuc(int shell_p, int shell_q);
    virtual void calc_kin(int shell_p, int shell_q);
    virtual void calc_ovp(int shell_p, int shell_q);
    virtual void calc_eri(int shell_p, int shell_q, int shell_r, int shell_s);

   protected:
    void makeAtomList();
    DfHpqEngine::PGTOs getPGTOs(const int shellIndex);

   private:
    std::vector<TlAtom> Cs_;
    std::vector<TlAtom> Xs_;
};

#endif  // TL_MINT_PDF_H
