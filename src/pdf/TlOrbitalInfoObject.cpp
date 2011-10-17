#include <iostream>
#include <cmath>
#include <limits>
#include "TlOrbitalInfoObject.h"
#include "Fl_Geometry.h"
#include "Fl_Gto.h"

const int TlOrbitalInfoObject::MAX_SHELL_TYPE = 3;
const double TlOrbitalInfoObject::INV_SQRT3 = 1.0 / std::sqrt(3.0);
const double TlOrbitalInfoObject::MAX_EXPONENT = std::log(std::numeric_limits<double>::max());
const char* TlOrbitalInfoObject::basisTypeNameTbl_[] = {
    "s",
    "px",
    "py",
    "pz",
    "dxy",
    "dxz",
    "dyz",
    "dx2-y2",
    "dz2",
    "---"
};


TlOrbitalInfoObject::TlOrbitalInfoObject()
{
}


TlOrbitalInfoObject::~TlOrbitalInfoObject()
{
}


void TlOrbitalInfoObject::setCGTO(const Fl_Gto& flGto)
{
    const int numOfCgto = flGto.getNumOfCGTOs();
    this->cgtos_.resize(numOfCgto);

    for (int cgto = 0; cgto < numOfCgto; ++cgto) {
        CGTO tmpCGTO;

        const int shell = flGto.getShell(cgto);
        int l = 0, m = 0, n = 0;
        switch (shell) {
        case 's':
            l = m = n = 0;
            break;
        case 'p':
            l = 1;
            m = n = 0;
            break;
        case 'd':
            l = m = 1;
            n = 0;
            break;
        default:
            std::cerr << "### Now we support 's' to 'd'. ###" << std::endl;
            exit(1);
        }

        const int cgtoContraction = flGto.getContraction(cgto);
        PGTOs tmpPGTOs(cgtoContraction);

        double factor = 0.0;
        const double normalizedFactor = flGto.getNormalizedfactor(cgto, l, m, n);
        for (int pgto = 0; pgto < cgtoContraction; ++pgto) {
            const double coef = 
                flGto.getCoefficient(cgto, pgto) *
                flGto.getNormalized(cgto, pgto, l, m, n) * normalizedFactor;
            tmpPGTOs[pgto].dCoefficient = coef;

            const double exponent = flGto.getExponent(cgto, pgto);
            tmpPGTOs[pgto].dExponent = exponent;

            if (exponent < TlOrbitalInfoObject::MAX_EXPONENT) {
                factor += coef * std::exp(exponent);
            }
        }

        // sort
        std::sort(tmpPGTOs.begin(), tmpPGTOs.end(), cgto_sort_functor());

        tmpCGTO.pgtos = tmpPGTOs;
        tmpCGTO.basisName = flGto.getBasisName(cgto);
        tmpCGTO.atomSymbol = flGto.getAtom(cgto);
        tmpCGTO.shell = shell;
        tmpCGTO.basisName = flGto.getLabel(cgto);
        tmpCGTO.factor = factor;
        
        this->cgtos_[cgto] = tmpCGTO;
    }
}


void TlOrbitalInfoObject::setCGTO_coulomb(const Fl_Gto& flGto)
{
    const int numOfCgto = flGto.getNumOfCGTOs();
    this->cgtos_.resize(numOfCgto);

    for (int cgto = 0; cgto < numOfCgto; ++cgto) {
        CGTO tmpCGTO;

        const int shell = flGto.getShell(cgto);
        int l = 0, m = 0, n = 0;
        switch (shell) {
        case 's':
            l = m = n = 0;
            break;
        case 'p':
            l = 1;
            m = n = 0;
            break;
        case 'd':
            l = m = 1;
            n = 0;
            break;
        default:
            std::cerr << "### Now we support 's' to 'd'. ###" << std::endl;
            exit(1);
        }

        const int cgtoContraction = 1;
        PGTOs tmpPGTOs(cgtoContraction);
        tmpPGTOs[0].dCoefficient = flGto.getCoulombnormalized(cgto, 0, l, m, n);
        tmpPGTOs[0].dExponent = flGto.getExponent(cgto, 0);

        tmpCGTO.pgtos = tmpPGTOs;
        tmpCGTO.basisName = flGto.getBasisName(cgto);
        tmpCGTO.atomSymbol = flGto.getAtom(cgto);
        tmpCGTO.shell = shell;
        tmpCGTO.basisName = flGto.getLabel(cgto);

        this->cgtos_[cgto] = tmpCGTO;
    }
}


void TlOrbitalInfoObject::setAtoms(const Fl_Geometry& flGeom)
{
    const int numOfAtoms = flGeom.getNumOfAtoms();

    this->atoms_.clear();
    this->atoms_.resize(numOfAtoms);

    //int basiscount = 0;
    for (int atom = 0; atom < numOfAtoms; ++atom) {
        const std::string atomName = flGeom.getAtom(atom);
        const std::string label = flGeom.getLabel(atom);

        this->atoms_[atom].setElement(flGeom.getAtom(atom));
        this->atoms_[atom].moveTo(flGeom.getCoordinate(atom));
        this->atoms_[atom].setCharge(flGeom.getCharge(atom));
        this->atoms_[atom].setName(label);
    }

    //assert(this->getNumOfOrbitals() == basiscount);
}


void TlOrbitalInfoObject::makeOrbitalTable()
{
    struct Orbital tmpOrbital;
    this->orbitals_.clear();
    index_type numOfOrbitals = 0;

    const int numOfCgtos = this->cgtos_.size();
    const int numOfAtoms = this->getNumOfAtoms();
    for (int atom = 0; atom < numOfAtoms; ++atom) {
        const std::string symbol = this->atoms_[atom].getSymbol();
        if (symbol == "X") {
            // dummy atom
            continue;
        }
        
        const std::string label = this->atoms_[atom].getName();
        tmpOrbital.atomIndex = atom;
        
        bool isFound = false;
        for (int cgto = 0; cgto < numOfCgtos; ++cgto) {
            if ((this->cgtos_[cgto].atomSymbol == symbol) && (this->cgtos_[cgto].label == label)) {
                tmpOrbital.cgtoIndex = cgto;
                
                unsigned int shellType = 0;
                const char shell = this->cgtos_[cgto].shell;
                switch (shell) {
                case 's':
                    shellType = 0;
                    break;
                case 'p':
                    shellType = 1;
                    break;
                case 'd':
                    shellType = 2;
                    break;
                default:
                    std::cerr << TlUtils::format("make table error. @ %s:%d\n",__FILE__, __LINE__) << std::endl;
                    break;
                }
                tmpOrbital.shellType = shellType;
                
                const int maxBasisType = shellType * 2  +1;
                const int shellType2 = shellType * shellType;
                for (int basisType = 0; basisType < maxBasisType; ++basisType) {
                    tmpOrbital.basisType = shellType2 + basisType; // s=0, px=1, py=2, ...
                    this->orbitals_.push_back(tmpOrbital);
                    ++numOfOrbitals;
                }

                isFound = true;
                //break;
            }
        }

        if (isFound == false) {
            std::cerr << TlUtils::format("basis not found: atom index = %d, symbol=\"%s\", label=\"%s\"",
                                         atom, symbol.c_str(), label.c_str())
                      << std::endl;
        }
    }
    assert((index_type)this->orbitals_.size() == numOfOrbitals);
}


std::vector<TlOrbitalInfoObject::index_type> TlOrbitalInfoObject::getStartIndexArrayOfShellGroup() const
{
    std::vector<index_type> answer;

    const index_type numOfOrbitals = this->getNumOfOrbitals();
    for (index_type i = 0; i < numOfOrbitals; ) {
        const int shellType = this->getShellType(i);
        answer.push_back(i);
        
        i += 2 * shellType + 1;
    }

    // swap techneque
    std::vector<index_type>(answer).swap(answer);

    return answer;
}




void TlOrbitalInfoObject::getPrefactor_S(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = 1.0;
}


void TlOrbitalInfoObject::getPrefactor_PX(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = pos.x();
}


void TlOrbitalInfoObject::getPrefactor_PY(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = pos.y();
}


void TlOrbitalInfoObject::getPrefactor_PZ(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = pos.z();
}


void TlOrbitalInfoObject::getPrefactor_DXY(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = pos.x() * pos.y();
}


void TlOrbitalInfoObject::getPrefactor_DXZ(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = pos.x() * pos.z();
}


void TlOrbitalInfoObject::getPrefactor_DYZ(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = pos.y() * pos.z();
}


void TlOrbitalInfoObject::getPrefactor_DXX_YY(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = 0.5 * (pos.x() * pos.x() - pos.y() * pos.y());
}


void TlOrbitalInfoObject::getPrefactor_DZZ(const TlPosition& pos, double* pValue)
{
    assert(pValue != NULL);
    *pValue = INV_SQRT3 * (pos.z()*pos.z() - 0.5*(pos.x()*pos.x() + pos.y()*pos.y()));
}


void TlOrbitalInfoObject::getGradPrefactor_S(const TlPosition& pos, const double alpha,
                                             double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    *pValueX = alpha2 * pos.x();
    *pValueY = alpha2 * pos.y();
    *pValueZ = alpha2 * pos.z();
}


void TlOrbitalInfoObject::getGradPrefactor_PX(const TlPosition& pos, const double alpha,
                                              double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    *pValueX = alpha2 * pos.x() * pos.x() -1.0;
    *pValueY = alpha2 * pos.x() * pos.y();
    *pValueZ = alpha2 * pos.x() * pos.z();
}


void TlOrbitalInfoObject::getGradPrefactor_PY(const TlPosition& pos, const double alpha,
                                              double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    *pValueX = alpha2 * pos.y() * pos.x();
    *pValueY = alpha2 * pos.y() * pos.y() -1.0;
    *pValueZ = alpha2 * pos.y() * pos.z();
}


void TlOrbitalInfoObject::getGradPrefactor_PZ(const TlPosition& pos, const double alpha,
                                              double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    *pValueX = alpha2 * pos.y() * pos.x();
    *pValueY = alpha2 * pos.y() * pos.y();
    *pValueZ = alpha2 * pos.y() * pos.z() -1.0;
}


void TlOrbitalInfoObject::getGradPrefactor_DXY(const TlPosition& pos, const double alpha,
                                               double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    const double xy = pos.x() * pos.y();
    *pValueX = alpha2 * xy * pos.x() - pos.y(); 
    *pValueY = alpha2 * xy * pos.y() - pos.x();
    *pValueZ = alpha2 * xy * pos.z();
}


void TlOrbitalInfoObject::getGradPrefactor_DXZ(const TlPosition& pos, const double alpha,
                                               double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    const double xz = pos.x() * pos.z();
    *pValueX = alpha2 * xz * pos.x() - pos.z(); 
    *pValueY = alpha2 * xz * pos.y();
    *pValueZ = alpha2 * xz * pos.z() - pos.x();
}


void TlOrbitalInfoObject::getGradPrefactor_DYZ(const TlPosition& pos, const double alpha,
                                               double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    const double yz = pos.y() * pos.z();
    *pValueX = alpha2 * yz * pos.x(); 
    *pValueY = alpha2 * yz * pos.y() - pos.z();
    *pValueZ = alpha2 * yz * pos.z() - pos.y();
}


void TlOrbitalInfoObject::getGradPrefactor_DXX_YY(const TlPosition& pos, const double alpha,
                                                  double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    const double xx = pos.x() * pos.x();
    const double xx_X = alpha2 * xx * pos.x() -2.0 * pos.x();
    const double xx_Y = alpha2 * xx * pos.y();
    const double xx_Z = alpha2 * xx * pos.z();
    
    const double yy = pos.y() * pos.y();
    const double yy_X = alpha2 * yy * pos.x();
    const double yy_Y = alpha2 * yy * pos.y() -2.0 * pos.y();
    const double yy_Z = alpha2 * yy * pos.z();
    
    *pValueX = 0.5 * (xx_X - yy_X);
    *pValueY = 0.5 * (xx_Y - yy_Y);
    *pValueZ = 0.5 * (xx_Z - yy_Z);
}


void TlOrbitalInfoObject::getGradPrefactor_DZZ(const TlPosition& pos, const double alpha,
                                               double* pValueX, double* pValueY, double* pValueZ)
{
    assert(pValueX != NULL);
    assert(pValueY != NULL);
    assert(pValueZ != NULL);
    const double alpha2 = 2.0 * alpha;
    const double xx = pos.x() * pos.x();
    const double xx_X = alpha2 * xx * pos.x() -2.0 * pos.x();
    const double xx_Y = alpha2 * xx * pos.y();
    const double xx_Z = alpha2 * xx * pos.z();
    
    const double yy = pos.y() * pos.y();
    const double yy_X = alpha2 * yy * pos.x();
    const double yy_Y = alpha2 * yy * pos.y() -2.0 * pos.y();
    const double yy_Z = alpha2 * yy * pos.z();
    
    const double zz = pos.z() * pos.z();
    const double zz_X = alpha2 * zz * pos.x();
    const double zz_Y = alpha2 * zz * pos.y();
    const double zz_Z = alpha2 * zz * pos.z() -2.0 * pos.z();
    
    *pValueX = INV_SQRT3 * (zz_X - 0.5 * (xx_X + yy_X));
    *pValueY = INV_SQRT3 * (zz_Y - 0.5 * (xx_Y + yy_Y));
    *pValueZ = INV_SQRT3 * (zz_Z - 0.5 * (xx_Z + yy_Z));
}


