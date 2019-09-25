#include "tl_mint_libcint.h"
#include "Fl_Geometry.h"
#include "TlOrbitalInfo.h"
#include "TlPrdctbl.h"

extern "C" {
#include "cint.h"
int cint1e_ovlp_cart(double *buf, int *shls, int *atm, int natm, int *bas,
                     int nbas, double *env);

int cint1e_ovlp_sph(double *buf, int *shls, int *atm, int natm, int *bas,
                    int nbas, double *env);
int cint1e_nuc_sph(double *buf, int *shls, int *atm, int natm, int *bas,
                   int nbas, double *env);
int cint1e_kin_sph(double *buf, int *shls, int *atm, int natm, int *bas,
                   int nbas, double *env);
};

const double TlMint_Libcint::AU2ANG = 0.52917721067;
const int TlMint_Libcint::ATOM_BLOCK_SIZE = 6;
const int TlMint_Libcint::BASISSET_BLOCK_SIZE = 8;

TlMint_Libcint::TlMint_Libcint(const Fl_Geometry &geom,
                               const TlOrbitalInfo &orbInfo)
    : TlMintObject(geom, orbInfo), valuesOffset_(0) {
    this->makeAtomTable(geom);
    this->makeBasissetTable(orbInfo);
}

TlMint_Libcint::~TlMint_Libcint() {}

void TlMint_Libcint::calc_ovp(int shell_p, int shell_q) {
    const int numOfAtoms = this->orbInfo_.getNumOfAtoms();
    const int numOfAOs = this->orbInfo_.getNumOfOrbitals();

    int shls[2];
    shls[0] = shell_p;
    shls[1] = shell_q;
    const int di = CINTcgto_cart(shell_p, &(this->basissetTable_[0]));
    const int dj = CINTcgto_cart(shell_q, &(this->basissetTable_[0]));

    std::vector<double> buf(di * dj);
    if (0 != cint1e_ovlp_sph(&(buf[0]), shls, &(this->atomTable_[0]),
                             numOfAtoms, &(this->basissetTable_[0]), numOfAOs,
                             &(this->values_[0]))) {
        for (int q = 0; q < dj; ++q) {
            for (int p = 0; p < di; ++p) {
                std::cout << TlUtils::format("(%d %d) = % f", shell_p + p,
                                             shell_q + q, buf[p + q * dj])
                          << std::endl;
            }
        }
    } else {
        std::cout << "This integral is 0." << std::endl;
    };
}

void TlMint_Libcint::calc_nuc(int shell_p, int shell_q) {
    const int numOfAtoms = this->orbInfo_.getNumOfAtoms();
    const int numOfAOs = this->orbInfo_.getNumOfOrbitals();

    int shls[2];
    shls[0] = shell_p;
    shls[1] = shell_q;
    const int di = CINTcgto_cart(shell_p, &(this->basissetTable_[0]));
    const int dj = CINTcgto_cart(shell_q, &(this->basissetTable_[0]));

    std::vector<double> buf(di * dj);
    if (0 != cint1e_nuc_sph(&(buf[0]), shls, &(this->atomTable_[0]), numOfAtoms,
                            &(this->basissetTable_[0]), numOfAOs,
                            &(this->values_[0]))) {
        for (int q = 0; q < dj; ++q) {
            for (int p = 0; p < di; ++p) {
                std::cout << TlUtils::format("(%d %d) = % f", shell_p + p,
                                             shell_q + q, buf[p + q * dj])
                          << std::endl;
            }
        }
    } else {
        std::cout << "This integral is 0." << std::endl;
    };
}

void TlMint_Libcint::calc_kin(int shell_p, int shell_q) {
    const int numOfAtoms = this->orbInfo_.getNumOfAtoms();
    const int numOfAOs = this->orbInfo_.getNumOfOrbitals();

    int shls[2];
    shls[0] = shell_p;
    shls[1] = shell_q;
    const int di = CINTcgto_cart(shell_p, &(this->basissetTable_[0]));
    const int dj = CINTcgto_cart(shell_q, &(this->basissetTable_[0]));

    std::vector<double> buf(di * dj);
    if (0 != cint1e_kin_sph(&(buf[0]), shls, &(this->atomTable_[0]), numOfAtoms,
                            &(this->basissetTable_[0]), numOfAOs,
                            &(this->values_[0]))) {
        for (int q = 0; q < dj; ++q) {
            for (int p = 0; p < di; ++p) {
                std::cout << TlUtils::format("(%d %d) = % f", shell_p + p,
                                             shell_q + q, buf[p + q * dj])
                          << std::endl;
            }
        }
    } else {
        std::cout << "This integral is 0." << std::endl;
    };
}

void TlMint_Libcint::calc_eri(int shell_p, int shell_q, int shell_r,
                              int shell_s) {
    const int numOfAtoms = this->orbInfo_.getNumOfAtoms();
    const int numOfAOs = this->orbInfo_.getNumOfOrbitals();

    // call two-electron cartesian integrals
    int shls[4];
    shls[0] = shell_p;
    shls[1] = shell_q;
    shls[2] = shell_r;
    shls[3] = shell_s;
    const int di = CINTcgto_cart(shell_p, &(this->basissetTable_[0]));
    const int dj = CINTcgto_cart(shell_q, &(this->basissetTable_[0]));
    const int dk = CINTcgto_cart(shell_r, &(this->basissetTable_[0]));
    const int dl = CINTcgto_cart(shell_s, &(this->basissetTable_[0]));

    std::vector<double> buf(di * dj * dk * dl);
    if (0 != cint2e_sph(&(buf[0]), shls, &(this->atomTable_[0]), numOfAtoms,
                        &(this->basissetTable_[0]), numOfAOs,
                        &(this->values_[0]), NULL)) {
        for (int s = 0; s < dl; ++s) {
            for (int r = 0; r < dk; ++r) {
                for (int q = 0; q < dj; ++q) {
                    for (int p = 0; p < di; ++p) {
                        // printf("(%d %d|%d %d) = % f\n", p, q, r, s,
                        // buf[((p*dj +q)*dk +r)*dl +s]);
                        std::cout << TlUtils::format(
                                         "(%d %d|%d %d) = % f", shell_p + p,
                                         shell_q + q, shell_r + r, shell_s + s,
                                         buf[p + (q + (r + s * dl) * dk) * dj])
                                  << std::endl;
                    }
                }
            }
        }
    } else {
        std::cout << "This integral is 0." << std::endl;
    }
}

void TlMint_Libcint::makeAtomTable(const Fl_Geometry &geom) {
    // ATM_SLOTS = 6;
    // atm[i*6+0]: 原子iの電荷
    // atm[i*6+1]: env配列のオフセット; (x, y, z) = (env[atm[i*6+1]+0],
    // env[atm[i*6+1]+1], env[atm[i*6+1]+2]) atm[i*6+2]: 原子iの原子核モデル; =2
    // gaussian核モデル atm[i*6+3]: env offset to save the nuclear charge
    // distribution parameter ζ atm[i*6+4]: 未使用 atm[i*6+5]: 未使用
    const int numOfAtoms = geom.getNumOfAtoms();

    this->atomTable_.clear();
    this->atomTable_.resize(numOfAtoms * TlMint_Libcint::ATOM_BLOCK_SIZE);

    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const int base = atomIndex * TlMint_Libcint::ATOM_BLOCK_SIZE;
        const TlAtom &atom = geom.getAtom(atomIndex);
        this->atomTable_[base + 0] =
            TlPrdctbl::getAtomicNumber(atom.getSymbol());
        this->atomTable_[base + 1] = this->getValuesOffset();
        this->addValue(atom.getPosition().x() * AU2ANG);
        this->addValue(atom.getPosition().y() * AU2ANG);
        this->addValue(atom.getPosition().z() * AU2ANG);
        this->atomTable_[base + 2] = 0;
        this->atomTable_[base + 3] = 0;
        this->atomTable_[base + 4] = 0;
        this->atomTable_[base + 5] = 0;
    }
}

void TlMint_Libcint::makeBasissetTable(const TlOrbitalInfo &orbInfo) {
    // BAS_SLOTS = 8;
    // bas[i*8+0]: 対応原子の添字(0-based)
    // bas[i*8+1]: angular momentum
    // bas[i*8+2]: i番目のbasisにおける原始-GTOの数
    // bas[i*8+3]: i番目のbasisにおける縮約-GTOの数
    // bas[i*8+4]: kappa for spinor GTO.
    //             < 0 the basis ~ j = l + 1/2.
    //             > 0 the basis ~ j = l - 1/2.
    //             = 0 the basis includes both j = l + 1/2 and j = l - 1/2
    // bas[i*8+5]: 原始GTO exp のオフセット
    //             (env[bas[i*8+5]+0], env[bas[i*8+5]+1],...)
    // bas[i*8+6]: 原始GTO coef のオフセット
    //             env offset to save column-major
    //             contraction coefficients.
    //             e.g. 10 primitive -> 5 contraction needs a 10 × 5 array
    //
    //             env[bas[i*8+6]  ] | env[bas[i*8+6]+10] |     |
    //             env[bas[i*8+6]+40] env[bas[i*8+6]+1] | env[bas[i*8+6]+11] | |
    //             env[bas[i*8+6]+41] .                 | .                  |
    //             ... | . .                 | .                  |     | .
    //             env[bas[i*8+6]+9] | env[bas[i*8+6]+19] |     |
    //             env[bas[i*8+6]+49]
    // bas[i*8+7]: 未使用
    const int numOfAOs = orbInfo.getNumOfOrbitals();

    this->basissetTable_.clear();
    this->basissetTable_.resize(numOfAOs * TlMint_Libcint::BASISSET_BLOCK_SIZE);
    for (int aoIndex = 0; aoIndex < numOfAOs; ++aoIndex) {
        const int base = aoIndex * TlMint_Libcint::BASISSET_BLOCK_SIZE;
        this->basissetTable_[base + 0] = orbInfo.getAtomIndex(aoIndex);
        this->basissetTable_[base + 1] = orbInfo.getShellType(aoIndex);

        const int numOfPGTOs = orbInfo.getCgtoContraction(aoIndex);
        this->basissetTable_[base + 2] = numOfPGTOs;
        this->basissetTable_[base + 3] = 1;
        this->basissetTable_[base + 4] = 0;
        this->basissetTable_[base + 5] = this->getValuesOffset();
        for (int pgtoIndex = 0; pgtoIndex < numOfPGTOs; ++pgtoIndex) {
            this->addValue(orbInfo.getExponent(aoIndex, pgtoIndex));
        }

        this->basissetTable_[base + 6] = this->getValuesOffset();
        for (int pgtoIndex = 0; pgtoIndex < numOfPGTOs; ++pgtoIndex) {
            double coef = orbInfo.getCoefficient(aoIndex, pgtoIndex);
            double coef2 =
                CINTgto_norm(orbInfo.getShellType(aoIndex),
                             orbInfo.getExponent(aoIndex, pgtoIndex));
            this->addValue(coef * coef2);
        }

        this->basissetTable_[base + 7] = 0;
    }
}

int TlMint_Libcint::getValuesOffset() const { return this->valuesOffset_; }

void TlMint_Libcint::addValue(const double value) {
    this->values_.push_back(value);
    this->valuesOffset_ += 1;
    assert(this->values_.size() == this->valuesOffset_);
}

void TlMint_Libcint::showAtomTable() const {
    const int atomTableSize = this->atomTable_.size();
    for (int i = 0; i < atomTableSize; i += 6) {
        std::cout << TlUtils::format("charge: %d", this->atomTable_[i + 0])
                  << std::endl;
        std::cout << TlUtils::format("offset: %d", this->atomTable_[i + 1])
                  << std::endl;
        const int atom_offset = this->atomTable_[i + 1];
        std::cout << TlUtils::format("(% f, % f, % f)",
                                     this->values_[atom_offset],
                                     this->values_[atom_offset + 1],
                                     this->values_[atom_offset + 2])
                  << std::endl;
        std::cout << TlUtils::format("model:  %d", this->atomTable_[i + 2])
                  << std::endl;
        std::cout << TlUtils::format("dist:   %d", this->atomTable_[i + 3])
                  << std::endl;
        std::cout << std::endl;
    }
}

void TlMint_Libcint::showBasissetTable() const {
    const int basissetTableSize = this->basissetTable_.size();
    for (int i = 0; i < basissetTableSize; i += 8) {
        std::cout << TlUtils::format("atom id: %2d",
                                     this->basissetTable_[i + 0])
                  << std::endl;
        std::cout << TlUtils::format("ang. momentum: %d",
                                     this->basissetTable_[i + 1])
                  << std::endl;
        std::cout << TlUtils::format("# pGTOs: %d", this->basissetTable_[i + 2])
                  << std::endl;
        std::cout << TlUtils::format("# CTOs:  %d", this->basissetTable_[i + 3])
                  << std::endl;
        std::cout << TlUtils::format("kappa: %d", this->basissetTable_[i + 4])
                  << std::endl;
        std::cout << TlUtils::format("exp offset: %d",
                                     this->basissetTable_[i + 5])
                  << std::endl;
        std::cout << TlUtils::format("coef offset: %d",
                                     this->basissetTable_[i + 6])
                  << std::endl;
        std::cout << std::endl;
    }
}
