#include <iostream>
#include <string>
#include <vector>

#include "Fl_Geometry.h"
#include "TlGetopt.h"
#include "TlMsgPack.h"
#include "TlOrbitalInfo.h"
#include "TlSerializeData.h"
#include "TlUtils.h"
#include "tl_mint_object.h"
#include "tl_mint_pdf.h"

#ifdef HAVE_LIBCINT
#include "tl_mint_libcint.h"
#endif // HAVE_LIBCINT

void showHelp() {
  std::cout << "Molecular INTegral program: " << std::endl;
  std::cout << "mint <mpac_file>" << std::endl;
}

// ----------------------------------------------------------------------------
// mint
// ----------------------------------------------------------------------------
void doMint(TlMintObject* pMint, const int shell_p, const int shell_q) {
  if (pMint->checkAO(shell_p) && pMint->checkAO(shell_q)) {
    std::cout << ">>>> ovp" << std::endl;
    pMint->calc_ovp(shell_p, shell_q);
    std::cout << ">>>> nuc" << std::endl;
    pMint->calc_nuc(shell_p, shell_q);
    std::cout << ">>>> kin" << std::endl;
    pMint->calc_kin(shell_p, shell_q);
  }
}

void doMint(TlMintObject* pMint, const int shell_p, const int shell_q, const int shell_r, const int shell_s) {
  if (pMint->checkAO(shell_p) && pMint->checkAO(shell_q) &&
      pMint->checkAO(shell_r) && pMint->checkAO(shell_s)) {
    std::cout << ">>>> eri" << std::endl;
    pMint->calc_eri(shell_p, shell_q, shell_r, shell_s);
  }
}
// ----------------------------------------------------------------------------
// main 
// ----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
  TlGetopt opt(argc, argv, "hv");

  if ((opt.getCount() < 2) || (opt["h"] == "defined")) {
    showHelp();
    return 1;
  }

  bool verbose = false;
  if (opt["v"] == "defined") {
    verbose = true;
  }

  const std::string mpacPath = opt[1];
  if (verbose) {
    std::cerr << "mpac file: " << mpacPath << std::endl;
  }

  int shell_p = 0;
  int shell_q = 0;
  int shell_r = 0;
  int shell_s = 0;
  bool enable_1e = true;
  bool enable_2e = false;
  if (opt.getCount() >= 3) {
    shell_p = std::atoi(opt[2].c_str());
  }
  if (opt.getCount() >= 4) {
    enable_1e = true;
    shell_q = std::atoi(opt[3].c_str());
  }
  if (opt.getCount() >= 5) {
    shell_r = std::atoi(opt[4].c_str());
  }
  if (opt.getCount() >= 6) {
    enable_2e = true;
    shell_s = std::atoi(opt[5].c_str());
  }
  if (verbose) {
    std::cerr << TlUtils::format("input: %d %d %d %d", shell_p, shell_q,
                                 shell_r, shell_s)
              << std::endl;
  }

  TlSerializeData data;
  {
    TlMsgPack mpac;
    mpac.load(mpacPath);
    data = mpac.getSerializeData();
  }

  Fl_Geometry geom(data["coordinates"]);
  // std::cerr << "atoms: " << geom.getNumOfAtoms() << std::endl;

  TlOrbitalInfo orbInfo(data["coordinates"], data["basis_set"]);
  // orbInfo.printCGTOs(std::cout);


{
  std::cout << "# pdf ====" << std::endl;
  TlMintObject* pMint = new TlMint_Pdf(geom, orbInfo);
  if (enable_1e) {
    doMint(pMint, shell_p, shell_q);
  }
  if (enable_2e) {
    doMint(pMint, shell_p, shell_q, shell_r, shell_s);
  }
  delete pMint;
  pMint = NULL;
}

#ifdef HAVE_LIBCINT
{
  std::cout << "# libcint ====" << std::endl;
  TlMint_Libcint* pMint = new TlMint_Libcint(geom, orbInfo);
  if (enable_1e) {
    doMint(pMint, shell_p, shell_q);
  }
  if (enable_2e) {
    doMint(pMint, shell_p, shell_q, shell_r, shell_s);
  }

  //pMint->showAtomTable();
  pMint->showBasissetTable();

  delete pMint;
  pMint = NULL;
}
#endif // HAVE_LIBCINT

  return 0;
}
