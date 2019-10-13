#include "tl_mint_pdf.h"
#include "DfEriEngine.h"
#include "DfOverlapEngine.h"
#include "Fl_Geometry.h"
#include "TlOrbitalInfo.h"
#include "TlPrdctbl.h"

TlMint_Pdf::TlMint_Pdf(const Fl_Geometry& geom, const TlOrbitalInfo& orbInfo)
    : TlMintObject(geom, orbInfo) {
    this->makeAtomList();
}
TlMint_Pdf::~TlMint_Pdf() {}

void TlMint_Pdf::calc_nuc(int shell_p, int shell_q) {
    const int shellTypeP = this->orbInfo_.getShellType(shell_p);
    const int shellTypeQ = this->orbInfo_.getShellType(shell_q);
    const int maxStepsP = 2 * shellTypeP + 1;
    const int maxStepsQ = 2 * shellTypeQ + 1;
    const TlPosition posP = this->orbInfo_.getPosition(shell_p);
    const TlPosition posQ = this->orbInfo_.getPosition(shell_q);
    const DfHpqEngine::PGTOs pgtosP = this->getPGTOs(shell_p);
    const DfHpqEngine::PGTOs pgtosQ = this->getPGTOs(shell_q);
    const DfHpqEngine::Query query(0, 0, shellTypeP, shellTypeQ);

    DfHpqEngine engine;
    engine.calc(query, posP, posQ, pgtosP, pgtosQ, this->Cs_, this->Xs_);
    int i = 0;
    for (int p = 0; p < maxStepsP; ++p) {
        for (int q = 0; q < maxStepsQ; ++q) {
            std::cout << TlUtils::format("(%d %d) = % f", shell_p + p,
                                         shell_q + q, engine.WORK_NUC[i])
                      << std::endl;
            ++i;
        }
    }
}

void TlMint_Pdf::calc_kin(int shell_p, int shell_q) {
    const int shellTypeP = this->orbInfo_.getShellType(shell_p);
    const int shellTypeQ = this->orbInfo_.getShellType(shell_q);
    const int maxStepsP = 2 * shellTypeP + 1;
    const int maxStepsQ = 2 * shellTypeQ + 1;
    const TlPosition posP = this->orbInfo_.getPosition(shell_p);
    const TlPosition posQ = this->orbInfo_.getPosition(shell_q);
    const DfHpqEngine::PGTOs pgtosP = this->getPGTOs(shell_p);
    const DfHpqEngine::PGTOs pgtosQ = this->getPGTOs(shell_q);
    const DfHpqEngine::Query query(0, 0, shellTypeP, shellTypeQ);

    DfHpqEngine engine;
    engine.calc(query, posP, posQ, pgtosP, pgtosQ, this->Cs_, this->Xs_);
    int i = 0;
    for (int p = 0; p < maxStepsP; ++p) {
        for (int q = 0; q < maxStepsQ; ++q) {
            std::cout << TlUtils::format("(%d %d) = % f", shell_p + p,
                                         shell_q + q, engine.WORK_KIN[i])
                      << std::endl;
            ++i;
        }
    }
}

void TlMint_Pdf::calc_ovp(int shell_p, int shell_q) {
    const int shellTypeP = this->orbInfo_.getShellType(shell_p);
    const int shellTypeQ = this->orbInfo_.getShellType(shell_q);
    const int maxStepsP = 2 * shellTypeP + 1;
    const int maxStepsQ = 2 * shellTypeQ + 1;

    DfOverlapEngine engine;
    engine.calc(0, this->orbInfo_, shell_p, 0, this->orbInfo_, shell_q, 0,
                this->orbInfo_, -1, 0, this->orbInfo_, -1);

    int i = 0;
    for (int p = 0; p < maxStepsP; ++p) {
        for (int q = 0; q < maxStepsQ; ++q) {
            std::cout << TlUtils::format("(%d %d) = % f", shell_p + p,
                                         shell_q + q, engine.WORK[i])
                      << std::endl;
            ++i;
        }
    }
}

void TlMint_Pdf::calc_eri(int shell_p, int shell_q, int shell_r, int shell_s) {
    const int shellTypeP = this->orbInfo_.getShellType(shell_p);
    const int shellTypeQ = this->orbInfo_.getShellType(shell_q);
    const int shellTypeR = this->orbInfo_.getShellType(shell_r);
    const int shellTypeS = this->orbInfo_.getShellType(shell_s);
    const int maxStepsP = 2 * shellTypeP + 1;
    const int maxStepsQ = 2 * shellTypeQ + 1;
    const int maxStepsR = 2 * shellTypeR + 1;
    const int maxStepsS = 2 * shellTypeS + 1;

    DfEriEngine engine;
    engine.calc(0, this->orbInfo_, shell_p, 0, this->orbInfo_, shell_q, 0,
                this->orbInfo_, shell_r, 0, this->orbInfo_, shell_s);

    int i = 0;
    for (int p = 0; p < maxStepsP; ++p) {
        for (int q = 0; q < maxStepsQ; ++q) {
            for (int r = 0; r < maxStepsR; ++r) {
                for (int s = 0; s < maxStepsS; ++s) {
                    std::cout
                        << TlUtils::format("(%d %d|%d %d) = % f", shell_p + p,
                                           shell_q + q, shell_r + r,
                                           shell_s + s, engine.WORK[i])
                        << std::endl;
                    ++i;
                }
            }
        }
    }
}

void TlMint_Pdf::makeAtomList() {
    this->Cs_.clear();
    this->Xs_.clear();
    const int numOfAtoms = this->geom_.getNumOfAtoms();
    for (int i = 0; i < numOfAtoms; ++i) {
        TlAtom atom = this->geom_.getAtom(i);
        if (atom.getSymbol() == "X") {
            this->Xs_.push_back(atom);
        } else {
            atom.setCharge(TlPrdctbl::getAtomicNumber(atom.getSymbol()));
            this->Cs_.push_back(atom);
        }
    }
}

DfHpqEngine::PGTOs TlMint_Pdf::getPGTOs(const int shellIndex) {
    const int numOfContractions = this->orbInfo_.getCgtoContraction(shellIndex);
    DfHpqEngine::PGTOs pgtos(numOfContractions);

    for (int i = 0; i < numOfContractions; ++i) {
        const DfHpqEngine::PGTO pgto(
            this->orbInfo_.getCoefficient(shellIndex, i),
            this->orbInfo_.getExponent(shellIndex, i));
        pgtos[i] = pgto;
    }

    return pgtos;
}
