#include <iostream>
#include "TlOrbitalInfo_Density.h"

TlOrbitalInfo_Density::TlOrbitalInfo_Density(const TlSerializeData& geomData,
                                             const TlSerializeData& basisData)
{
    this->setCGTO_coulomb(Fl_Gto(basisData));
    this->setAtoms(Fl_Geometry(geomData));
    this->makeOrbitalTable();
    
    // std::cerr << "TlOrbitalInfo_Density::TlOrbitalInfo_Density(): "
    //           << this->cgtos_.size() << "/"
    //           << this->getNumOfOrbitals()
    //           << std::endl;
}


TlOrbitalInfo_Density::~TlOrbitalInfo_Density()
{
}


// void TlOrbitalInfo_Density::setDataFromFlTbl()
// {
//     // set data from Fl_Gto
//     const int numOfCgto = this->flGto_.getNumOfCGTOs();
//     //std::cerr << "TlOrbitalInfo_Density::setDataFromFlTbl() numOfCgto = " << numOfCgto << std::endl;
//     this->cgto_.resize(numOfCgto);

//     for (int i = 0; i < numOfCgto; ++i) {
//         int l = 0, m = 0, n = 0;
//         switch (this->flGto_.getShell(i)) {
//         case 's':
//             l = m = n = 0;
//             break;
//         case 'p':
//             l = 1;
//             m = n = 0;
//             break;
//         case 'd':
//             l = m = 1;
//             n = 0;
//             break;
//         default:
//             std::cerr << "### Now we support 's' to 'd'. ###" << std::endl;
//             exit(1);
//         }

//         this->cgto_[i].resize(1);
//         this->cgto_[i][0].dCoefficient = this->flGto_.getCoulombnormalized(i, 0, l, m, n);
//         this->cgto_[i][0].dExponent = this->flGto_.getExponent(i, 0);
// //     const int cgtoContraction = 1; this->pFlGto_->getContraction(i);
// //     this->m_aCGTO[i].resize(cgtoContraction);
// //     const double normalizedFactor = 1.0; //this->pFlGto_->getNormalizedfactor(i, l, m, n);
// //     for (int k = 0; k < cgtoContraction; ++k) {
// //       this->m_aCGTO[i][k].dCoefficient =
// //  this->pFlGto_->getCoefficient(i, k) *
// //  this->pFlGto_->getNormalized(i, k, l, m, n) *
// //  normalizedFactor;
// //       this->m_aCGTO[i][k].dExponent = this->pFlGto_->getExponent(i, k);
// //     }
//     }

//     // sort by exponent
//     CGTO::const_iterator pCgtoEnd = this->cgto_.end();
//     for (CGTO::iterator pCGTO = this->cgto_.begin(); pCGTO != pCgtoEnd; ++pCGTO) {
//         std::sort(pCGTO->begin(), pCGTO->end(), cgto_sort_functor());
//     }

//     //
//     const std::size_t atomNum = this->flGeom_.getNumOfAtoms();
//     const int atomKindNum = this->flGeom_.getAtomKindNumber();

//     this->atoms_.resize(atomNum);
//     struct Orbital tmpOrbital;
//     int basiscount = 0;
//     for (std::size_t i = 0; i < atomNum; ++i) {
//         this->atoms_[i].setElement(this->flGeom_.getAtom(i));
//         this->atoms_[i].moveTo(this->flGeom_.getCoordinate(i));

//         tmpOrbital.nAtomIndex = i;
//         const std::string atomName = this->flGeom_.getAtom(i);
//         const std::string label = this->flGeom_.getLabel(i);

//         for (int j = 0; j < atomKindNum; ++j) {
//             const int m = this->flGto_.getStartposition(j);

//             if ((this->flGto_.getAtom(m) == atomName) &&
//                 (this->flGto_.getLabel(m) == label)) {
//                 const int start = m;
//                 const int to = m + this->flGto_.getTermnumber(j);
//                 for (int k = start; k < to; ++k) {
//                     const int cGtonum = k;
//                     tmpOrbital.nCgtoIndex = cGtonum;

//                     const char shell = this->flGto_.getShell(k);
//                     switch (shell) {
//                     case 's':
//                         tmpOrbital.nType = 0;
//                         break;
//                     case 'p':
//                         tmpOrbital.nType = 1;
//                         break;
//                     case 'd':
//                         tmpOrbital.nType = 2;
//                         break;
//                     default:
//                         std::cerr << TlUtils::format("make table error. @ %s:%d\n",__FILE__, __LINE__) << std::endl;
//                         break;
//                     }

//                     const int maxT = tmpOrbital.nType * 2  +1;
//                     for (int t = 0; t < maxT; ++t) {
//                         tmpOrbital.nBasisType = (tmpOrbital.nType * tmpOrbital.nType) + t; // s=0, px=1, py=2, ...
//                         this->orbitals_.push_back(tmpOrbital);
//                         ++basiscount;
//                     }
//                 }

//                 break;
//             }
//         }
//     }

//     assert(this->getNumOfOrbitals() == static_cast<std::size_t>(basiscount));
// }
