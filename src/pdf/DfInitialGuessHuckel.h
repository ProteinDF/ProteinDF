#ifndef DFINITIALGUESSHUCKEL_H
#define DFINITIALGUESSHUCKEL_H

#include <string>

#include "DfObject.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"

/** 拡張Huckel法による初期値を作成する
 */
class DfInitialGuessHuckel : public DfObject {
public:
    DfInitialGuessHuckel(TlSerializeData* pPdfParam);
    ~DfInitialGuessHuckel();

private:
    void initialize();

    void createGuess();
    TlSymmetricMatrix getHuckelMatrix();
    TlSymmetricMatrix getHpqMatrix();
    double getHii(const std::string& sAtomName, int nOrbitalType);

    TlVector generateOccupation(RUN_TYPE nRunType, int nNumOrbcut, int nNumOfElectrons);
    void generatePMatrix(const TlMatrix& C, const TlVector& occ,
                         TlSymmetricMatrix& P1, TlSymmetricMatrix& P2);

protected:
//     GUESS_TYPE m_nGuessType;

    static const double EPS;
    std::map<std::string, std::map<int, double> > m_Hii; // [atom_number][orbital_type] = value
};

#endif // DFINITIALGUESSHUCKEL_H
