#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "DfGridFreeXC.h"
#include "DfXCFunctional.h"
#include "DfFunctional_SVWN.h"
#include "DfFunctional_HFS.h"
#include "DfFunctional_Becke88.h"
#include "DfFunctional_B88LYP.h"
#include "DfFunctional_B3LYP.h"
#include "DfOverlapX.h"
#include "DfXMatrix.h"
#include "DfCD.h"
#include "TlTime.h"
#include "TlRowVectorMatrix2.h"
#include "TlSystem.h"
#include "TlUtils.h"

const int DfGridFreeXC::MAX_SHELL_TYPE = 2 + 1;
const double DfGridFreeXC::ONE_THIRD = 1.0 / 3.0;

DfGridFreeXC::DfGridFreeXC(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam), pOvpEngines_(NULL),
      orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_sets"]) {

    this->numOfPQs_ = this->m_nNumOfAOs * (this->m_nNumOfAOs + 1) / 2;

    this->tau_ = 1.0E-10;
    if ((*pPdfParam)["gridfree/CDAM_tau"].getStr().empty() != true) {
        this->tau_ = (*pPdfParam)["grid_free/CDAM_tau"].getDouble();
    }    

    this->epsilon_ = 1.0E-4;
    if ((*pPdfParam)["gridfree/CD_epsilon"].getStr().empty() != true) {
        this->epsilon_ = (*pPdfParam)["grid_free/CD_epsilon"].getDouble();
    }

    this->isCanonicalOrthogonalize_ = true;
    if ((*pPdfParam)["gridfree/orthogonalize_method"].getStr().empty() != true) {
        const std::string method = TlUtils::toUpper((*pPdfParam)["gridfree/orthogonalize_method"].getStr());
        if (method == "LOWDIN") {
            this->isCanonicalOrthogonalize_ = false;
        }
    }

    this->debugSaveM_ = (*pPdfParam)["debug/DfGridFreeXC/saveM"].getBoolean();
    if (this->debugSaveM_) {
        this->log_.info("using GAMESS formula");
    }
}


DfGridFreeXC::~DfGridFreeXC()
{
}

DfOverlapX* DfGridFreeXC::getDfOverlapObject()
{
    DfOverlapX* pDfOverlapX = new DfOverlapX(this->pPdfParam_);
    return pDfOverlapX;
}

DfXMatrix* DfGridFreeXC::getDfXMatrixObject()
{
    DfXMatrix* pDfXMatrix = new DfXMatrix(this->pPdfParam_);
    return pDfXMatrix;
}

// before SCF ==================================================================
void DfGridFreeXC::preprocessBeforeSCF()
{
    this->preprocessBeforeSCF_templ<DfOverlapX, DfXMatrix,
                                    TlSymmetricMatrix, TlMatrix>();
}

// in SCF ======================================================================
void DfGridFreeXC::buildFxc()
{
    const DfXCFunctional xcFunc(this->pPdfParam_);
    if (xcFunc.getXcType() == DfXCFunctional::HF) {
        // need not pure-DFT term
        return;
    }

    const DfXCFunctional::FUNCTIONAL_TYPE funcType = xcFunc.getFunctionalType();
    switch (funcType) {
    case DfXCFunctional::LDA:
        this->buildFxc_LDA();
        break;

    case DfXCFunctional::GGA:
        this->buildFxc_GGA();
        break;

    default:
        this->log_.critical("unknown XC functional type. stop.");
        CnErr.abort();
        break;
    }
}

void DfGridFreeXC::buildFxc_LDA()
{
    this->buildFxc_LDA_method<DfOverlapX, DfCD, TlSymmetricMatrix, TlMatrix>();
}

void DfGridFreeXC::createEngines()
{
    assert(this->pOvpEngines_ == NULL);
    
    this->log_.info(TlUtils::format("create ERI engine: %d",
                                    this->numOfThreads_));
    this->pOvpEngines_ = new DfOverlapEngine[this->numOfThreads_];
}


void DfGridFreeXC::destroyEngines()
{
    this->log_.info("delete OpenMP ERI engine");
    if (this->pOvpEngines_ != NULL) {
        delete[] this->pOvpEngines_;
    }
    this->pOvpEngines_ = NULL;
}


DfTaskCtrl* DfGridFreeXC::getDfTaskCtrlObject() const
{
    DfTaskCtrl* pDfTaskCtrl = new DfTaskCtrl(this->pPdfParam_);
    return pDfTaskCtrl;
}


void DfGridFreeXC::finalize(TlSymmetricMatrix* pMtx)
{
    // do nothing
}

void DfGridFreeXC::get_F_lamda(const TlVector lamda,
                               TlMatrixObject* pF_lamda,
                               TlMatrixObject* pE_lamda)
{
    const int dim = lamda.getSize();
    assert(pF_lamda->getNumOfRows() == dim);
    assert(pF_lamda->getNumOfCols() == dim);
    assert(pE_lamda->getNumOfRows() == dim);
    assert(pE_lamda->getNumOfCols() == dim);

    DfFunctional_LDA* pFunc = NULL;
    std::string checkXC = this->m_sXCFunctional;
    if (checkXC == "SVWN") {
        pFunc = new DfFunctional_SVWN();
    } else if (checkXC == "HFS") {
        pFunc = new DfFunctional_HFS();
    } else {
        this->log_.critical(TlUtils::format("not support functional: %s", checkXC.c_str()));
        abort();
    }

    double fv_a = 0.0;
    double fv_b = 0.0;
    for (int i = 0; i < dim; ++i) {
        const double v = lamda.get(i);
        if (v > 1.0E-16) {
            pFunc->getDerivativeFunctional(v, v, &fv_a, &fv_b);
            pF_lamda->set(i, i, fv_a);

            const double f = pFunc->getFunctional(v, v) / (2.0 * v);
            pE_lamda->set(i, i, f);
        }
    }

    delete pFunc;
    pFunc = NULL;
}


void DfGridFreeXC::getM_exact(const TlSymmetricMatrix& P, TlSymmetricMatrix* pM)
{
    assert(pM != NULL);
    TlMatrix M(this->m_nNumOfAOs, this->m_nNumOfAOs);
    pM->resize(this->m_nNumOfAOs);

    DfOverlapEngine engine;
    
    const TlOrbitalInfo orbitalInfo((*(this->pPdfParam_))["coordinates"],
                                    (*(this->pPdfParam_))["basis_sets"]);

    const ShellArrayTable shellArrayTable = this->makeShellArrayTable(orbitalInfo);
    // const ShellPairArrayTable shellPairArrayTable = this->getShellPairArrayTable(shellArrayTable);

    for (int shellTypeP = MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const int maxStepsP = 2 * shellTypeP + 1;
        const ShellArray shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeQ = MAX_SHELL_TYPE -1; shellTypeQ >= 0; --shellTypeQ) {
            const int maxStepsQ = 2 * shellTypeQ + 1;
            const ShellArray shellArrayQ = shellArrayTable[shellTypeQ];
            ShellArray::const_iterator qItEnd = shellArrayQ.end();
            
            for (int shellTypeR = MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
                const int maxStepsR = 2 * shellTypeR + 1;
                const ShellArray shellArrayR = shellArrayTable[shellTypeR];
                ShellArray::const_iterator rItEnd = shellArrayR.end();

                for (int shellTypeS = MAX_SHELL_TYPE -1; shellTypeS >= 0; --shellTypeS) {
                    const int maxStepsS = 2 * shellTypeS + 1;
                    const ShellArray shellArrayS = shellArrayTable[shellTypeS];
                    ShellArray::const_iterator sItEnd = shellArrayS.end();

                    const DfOverlapEngine::Query query(0, 0, 0, 0,
                                                       shellTypeP, shellTypeQ,
                                                       shellTypeR, shellTypeS);

                    for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                        const index_type shellIndexP = *pIt;
                        // const TlPosition posP = orbitalInfo.getPosition(shellIndexP);
                        // const DfOverlapEngine::PGTOs pgtosP = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexP);

                        for (ShellArray::const_iterator qIt = shellArrayQ.begin(); qIt != qItEnd; ++qIt) {
                            const index_type shellIndexQ = *qIt;
                            // const TlPosition posQ = orbitalInfo.getPosition(shellIndexQ);
                            // const DfOverlapEngine::PGTOs pgtosQ = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexQ);

                            for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                                const index_type shellIndexR = *rIt;
                                // const TlPosition posR = orbitalInfo.getPosition(shellIndexR);
                                // const DfOverlapEngine::PGTOs pgtosR = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexR);
                                
                                for (ShellArray::const_iterator sIt = shellArrayS.begin(); sIt != sItEnd; ++sIt) {
                                    const index_type shellIndexS = *sIt;
                                    // const TlPosition posS = orbitalInfo.getPosition(shellIndexS);
                                    // const DfOverlapEngine::PGTOs pgtosS = DfOverlapEngine::getPGTOs(orbitalInfo, shellIndexS);

                                    // engine.calc0(query,
                                    //              posP, posQ, posR, posS,
                                    //              pgtosP, pgtosQ, pgtosR, pgtosS);
                                    engine.calc(0, orbitalInfo, shellIndexP,
                                                0, orbitalInfo, shellIndexQ,
                                                0, orbitalInfo, shellIndexR,
                                                0, orbitalInfo, shellIndexS);

                                    int index = 0;
                                    for (int i = 0; i < maxStepsP; ++i) {
                                        const int indexP = shellIndexP + i;

                                        for (int j = 0; j < maxStepsQ; ++j) {
                                            const int indexQ = shellIndexQ + j;
                                            
                                            for (int k = 0; k < maxStepsR; ++k) {
                                                const int indexR = shellIndexR + k;

                                                for (int l = 0; l < maxStepsS; ++l) {
                                                    const int indexS = shellIndexS + l;
                                                    
                                                    const double P_rs = P.get(indexR, indexS);
                                                    const double value = engine.WORK[index];
                                                    M.add(indexP, indexQ, P_rs * value);

                                                    ++index;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    *pM = M;
}


DfGridFreeXC::ShellArrayTable DfGridFreeXC::makeShellArrayTable(const TlOrbitalInfoObject& orbitalInfo)
{
    ShellArrayTable shellArrayTable(MAX_SHELL_TYPE);
    const index_type maxShellIndex = orbitalInfo.getNumOfOrbitals();

    index_type shellIndex = 0;
    while (shellIndex < maxShellIndex) {
        // shellType: 0=s, 1=p, 2=d
        const int shellType = orbitalInfo.getShellType(shellIndex);
        const int steps = 2 * shellType +1;

        shellArrayTable[shellType].push_back(shellIndex);
        
        shellIndex += steps;
    }

    return shellArrayTable;
}


DfGridFreeXC::ShellPairArrayTable DfGridFreeXC::getShellPairArrayTable(const ShellArrayTable& shellArrayTable)
{
    ShellPairArrayTable shellPairArrayTable(MAX_SHELL_TYPE * MAX_SHELL_TYPE);

    for (int shellTypeP = MAX_SHELL_TYPE -1; shellTypeP >= 0; --shellTypeP) {
        const ShellArray& shellArrayP = shellArrayTable[shellTypeP];
        ShellArray::const_iterator pItEnd = shellArrayP.end();

        for (int shellTypeR = MAX_SHELL_TYPE -1; shellTypeR >= 0; --shellTypeR) {
            const ShellArray& shellArrayR = shellArrayTable[shellTypeR];
            ShellArray::const_iterator rItEnd = shellArrayR.end();

            const int shellPairType_PR = shellTypeP * MAX_SHELL_TYPE + shellTypeR;
            for (ShellArray::const_iterator pIt = shellArrayP.begin(); pIt != pItEnd; ++pIt) {
                const index_type indexP = *pIt;

                for (ShellArray::const_iterator rIt = shellArrayR.begin(); rIt != rItEnd; ++rIt) {
                    const index_type indexR = *rIt;

                    if (indexP >= indexR) {
                        ShellPair shellPair(indexP, indexR);
                        shellPairArrayTable[shellPairType_PR].push_back(shellPair);
                    }
                }
            }
        }
    }

    return shellPairArrayTable;
}

// TlSymmetricMatrix DfGridFreeXC::getPMatrix()
// {
//     TlSymmetricMatrix P = this->getPpqMatrix<TlSymmetricMatrix>(RUN_RKS, this->m_nIteration -1);
//     return P;
// }


TlMatrix DfGridFreeXC::getL()
{
    //TlMatrix L = DfObject::getLMatrix<TlMatrix>();
    TlMatrix L;
    L.load("GF_L.mat");

    return L;
}


DfGridFreeXC::PQ_PairArray DfGridFreeXC::getI2PQ()
{
    std::string filepath = this->getI2pqVtrPath();
    std::ifstream ifs;
    ifs.open(filepath.c_str(), std::ofstream::in | std::ofstream::binary);
    if (ifs.fail()) {
        abort();
    }

    std::size_t size = 0;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    PQ_PairArray answer(size);
    index_type shellIndex1 = 0;
    index_type shellIndex2 = 0;
    for (std::size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&shellIndex1), sizeof(index_type));
        ifs.read(reinterpret_cast<char*>(&shellIndex2), sizeof(index_type));
        answer[i] = IndexPair2(shellIndex1, shellIndex2);
    }

    ifs.close();
    return answer;
}


void DfGridFreeXC::divideCholeskyBasis(const index_type numOfCBs,
                                       index_type *pStart, index_type *pEnd)
{
    *pStart = 0;
    *pEnd = numOfCBs;
}


TlSymmetricMatrix DfGridFreeXC::getCholeskyVector(const TlVector& L_col,
                                                  const PQ_PairArray& I2PQ)
{
    const index_type numOfItilde = L_col.getSize();
    TlSymmetricMatrix answer(this->m_nNumOfAOs);
    for (index_type i = 0; i < numOfItilde; ++i) {
        answer.set(I2PQ[i].index1(),
                   I2PQ[i].index2(),
                   L_col[i]);
    }

    return answer;
}

// -----------------------------------------------------------------------------
void DfGridFreeXC::buildFxc_GGA()
{
    this->log_.info("build Fxc by grid-free method: functional type is GGA.");

    std::string basisset_param = "basis_sets";
    if (this->isDedicatedBasisForGridFree_) {
        basisset_param = "basis_sets_GF";
    }
    const TlOrbitalInfo orbitalInfo_GF((*(this->pPdfParam_))["coordinates"],
                                       (*(this->pPdfParam_))[basisset_param]);
    const index_type numOfAOs = this->m_nNumOfAOs;
    const index_type numOfGfOrbs = orbitalInfo_GF.getNumOfOrbitals();
    this->log_.info(TlUtils::format("AOs = %d", numOfAOs));
    this->log_.info(TlUtils::format("auxAOs for GF = %d", numOfGfOrbs));

    // RKS
    const TlSymmetricMatrix PA = 0.5 * DfObject::getPpqMatrix<TlSymmetricMatrix>(RUN_RKS,
                                                                                 this->m_nIteration -1);
    assert(PA.getNumOfRows() == numOfAOs);

    TlSymmetricMatrix M;
    if (this->XC_engine_ == XC_ENGINE_GRIDFREE_CD) {
        this->log_.info("begin to create M matrix based on CD.");
        {
            DfCD dfCD(this->pPdfParam_);
            dfCD.getM(PA, &M);
            // M *= 2.0;
        }
    } else {
        this->log_.info("begin to create M matrix using 4-center overlap.");
        DfOverlapX dfOvp(this->pPdfParam_);
        if (this->isDedicatedBasisForGridFree_) {
            dfOvp.getM_A(PA, &M);
        } else {
            dfOvp.getM(PA, &M);
        }
    }
    if (this->debugSaveM_) {
        M.save(TlUtils::format("fl_Work/debug_M.%d.mat", this->m_nIteration));
    }

    this->log_.info("begin to generate Fxc using grid-free method.");
    // M~(=V^t * M * V) および SVU(=SVU, but now only SV)の作成
    TlMatrix S;
    TlMatrix V;
    if (this->isDedicatedBasisForGridFree_) {
        S = DfObject::getGfStildeMatrix<TlMatrix>();
        V = DfObject::getGfVMatrix<TlMatrix>();
    } else {
        S = DfObject::getSpqMatrix<TlSymmetricMatrix>();
        V = DfObject::getXMatrix<TlMatrix>();
    }

    const index_type numOfGFOrthNormBasis = V.getNumOfCols();
    this->log_.info(TlUtils::format("orthonormal basis = %d", numOfGFOrthNormBasis));
    TlMatrix Vt = V;
    Vt.transpose();

    TlMatrix St = S;
    St.transpose();

    TlSymmetricMatrix Mtilde = Vt * M * V;
    //Mtilde.save("Mtilde.mat");

    TlVector lambda;
    TlMatrix U;
    Mtilde.diagonal(&lambda, &U);
    //lambda.save("lambda.vct");
    //U.save("U.mat");
    TlMatrix Ut = U;
    Ut.transpose();
    assert(lambda.getSize() == numOfGFOrthNormBasis);
    assert(U.getNumOfRows() == numOfGFOrthNormBasis);
    assert(U.getNumOfCols() == numOfGFOrthNormBasis);

    // M[rho^(-1/3)]
    TlSymmetricMatrix Mtilde_13;
    {
        TlSymmetricMatrix lambda_13(numOfGFOrthNormBasis);
        for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
            double v_13 = 0.0;
            const double v = lambda.get(i);
            if (v > 1.0E-16) {
                v_13 = std::pow(v, - ONE_THIRD);
            }
            lambda_13.set(i, i, v_13);
        }
        Mtilde_13 = U * lambda_13 * Ut;
    }
    //Mtilde_13.save("Mtilde_13.mat");

    // GGA用gradient
    TlMatrix Gx = this->getDipoleVelocityIntegralsXMatrix<TlMatrix>();
    TlMatrix Gy = this->getDipoleVelocityIntegralsYMatrix<TlMatrix>();
    TlMatrix Gz = this->getDipoleVelocityIntegralsZMatrix<TlMatrix>();
    Gx *= -1.0;
    Gy *= -1.0;
    Gz *= -1.0;
    const TlMatrix DX = Vt * Gx * V;
    const TlMatrix DY = Vt * Gy * V;
    const TlMatrix DZ = Vt * Gz * V;
    // DX.save("DX.mat");
    // DY.save("DY.mat");
    // DZ.save("DZ.mat");

    TlMatrix DXt = DX;
    DXt.transpose();
    TlMatrix DYt = DY;
    DYt.transpose();
    TlMatrix DZt = DZ;
    DZt.transpose();
    const TlMatrix RTX = 3.0*(DXt * Mtilde_13 + Mtilde_13 * DX);
    const TlMatrix RTY = 3.0*(DYt * Mtilde_13 + Mtilde_13 * DY);
    const TlMatrix RTZ = 3.0*(DZt * Mtilde_13 + Mtilde_13 * DZ);
    TlMatrix RTXt = RTX;
    RTXt.transpose();
    TlMatrix RTYt = RTY;
    RTYt.transpose();
    TlMatrix RTZt = RTZ;
    RTZt.transpose();

    // RX2 := M[{nabla rho / rho^(-4/3)}^2]
    const TlSymmetricMatrix RX2 = RTXt*RTX + RTYt*RTY + RTZt*RTZ;
    // RTX.save("RTX.mat");
    // RTY.save("RTY.mat");
    // RTZ.save("RTZ.mat");
    // RX2.save("RX2.mat");

    // RZ2 := [nabla rho / rho^(4/3)] * (DX + DY + DZ)
    const TlMatrix RZ2 = RTX*DX + RTY*DY + RTZ*DZ;
    //RZ2.save("RZ2.mat");
    TlMatrix RZ2t = RZ2;
    RZ2t.transpose();

    TlVector x2;
    TlMatrix Ux2;
    RX2.diagonal(&x2, &Ux2);
    //x2.save("x2.vct");
    //Ux2.save("Ux2.mat");
    TlMatrix Ux2t = Ux2;
    Ux2t.transpose();

    // ------------------
    assert(lambda.getSize() == numOfGFOrthNormBasis);
    // TlVector rhoAs(numOfGfOrbs);
    // TlVector xAs(numOfGfOrbs);
    TlVector rhoAs(numOfGFOrthNormBasis);
    TlVector xAs(numOfGFOrthNormBasis);
    //for (index_type i = 0; i < numOfGfOrbs; ++i) {
    for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
        const double rho_value = lambda[i];
        const double rho = (rho_value > 1.0E-16) ? rho_value : 0.0;
        rhoAs[i] = rho;
        
        const double x2_value = x2.get(i);
        const double x = (x2_value > 1.0E-16) ? std::sqrt(x2_value) : 0.0;
        xAs[i] = x;
    }
    const TlVector rhoBs = rhoAs;
    const TlVector xBs = xAs;
    
    DfFunctional_GGA* pFunc = this->getFunctionalGGA();
    // Fxc -------------------------------------------------------------
    TlSymmetricMatrix FxcA(numOfAOs); // alpha spin
    TlSymmetricMatrix FxcB(numOfAOs); // beta spin
    {
        DerivativeFunctionalSets dfs = pFunc->getDerivativeFunctional_GF(rhoAs, rhoBs, xAs, xBs);
        
        TlVector rhoAA43(numOfGFOrthNormBasis);
        TlVector rhoAB43(numOfGFOrthNormBasis);
        TlVector rhoBB43(numOfGFOrthNormBasis);
        for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
            const double rhoA = lambda[i];
            if (rhoA > 1.0E-16) {
                const double rhoA43 = std::pow(rhoA, 4.0/3.0);
                const double rhoB43 = rhoA43; // RKS
                rhoAA43[i] = rhoA43;
                rhoBB43[i] = rhoB43;
                rhoAB43[i] = std::sqrt(rhoA43 * rhoB43);
            }
        }

        TlMatrix FxcA_tilde(numOfGFOrthNormBasis, numOfGFOrthNormBasis);
        TlMatrix FxcB_tilde(numOfGFOrthNormBasis, numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_RAR(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_RAX(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_RBR(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_RBX(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_GAAR(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_GAAX(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_GABR(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_GABX(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_GBBR(numOfGFOrthNormBasis);
        TlSymmetricMatrix diag_GBBX(numOfGFOrthNormBasis);
        const int numOfTerms = pFunc->getNumOfDerivativeFunctionalTerms();
        for (int term = 0; term < numOfTerms; ++term) {
            for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
                diag_RAR(i, i) = dfs.rFrRhoA_R(term, i);
                diag_RAX(i, i) = dfs.rFrRhoA_X(term, i);
                diag_RBR(i, i) = dfs.rFrRhoB_R(term, i);
                diag_RBX(i, i) = dfs.rFrRhoB_X(term, i);

                diag_GAAR(i, i) = dfs.rFrGAA_R(term, i) * rhoAA43[i];
                diag_GAAX(i, i) = dfs.rFrGAA_X(term, i);
                diag_GABR(i, i) = dfs.rFrGAB_R(term, i) * rhoAB43[i];
                diag_GABX(i, i) = dfs.rFrGAB_X(term, i);
                diag_GBBR(i, i) = dfs.rFrGBB_R(term, i) * rhoBB43[i];
                diag_GBBX(i, i) = dfs.rFrGBB_X(term, i);
            }

            // alpha spin ------------
            {
                const TlSymmetricMatrix Fxc_RR = U * diag_RAR * Ut;
                const TlSymmetricMatrix Fxc_RX = Ux2 * diag_RAX * Ux2t;
                TlMatrix Fxc_tilde_term1 = 0.5 * (Fxc_RR * Fxc_RX + Fxc_RX * Fxc_RR);
                FxcA_tilde += Fxc_tilde_term1;
            }

            TlMatrix Fxc_GAA;
            {
                const TlMatrix Fxc_GAAR = U * diag_GAAR * Ut;
                const TlMatrix Fxc_GAAX = Ux2 * diag_GAAX * Ux2t;
                Fxc_GAA = 2.0 * 0.5 * (Fxc_GAAR * Fxc_GAAX + Fxc_GAAX * Fxc_GAAR);
            }

            TlMatrix Fxc_GAB;
            {
                const TlSymmetricMatrix Fxc_GABR = U * diag_GABR * Ut;
                const TlSymmetricMatrix Fxc_GABX = Ux2 * diag_GABX * Ux2t;
                Fxc_GAB = 0.5 * (Fxc_GABR * Fxc_GABX + Fxc_GABX * Fxc_GABR);
            }

            TlMatrix Fxc_GBB;
            {
                const TlSymmetricMatrix Fxc_GBBR = U * diag_GBBR * Ut;
                const TlSymmetricMatrix Fxc_GBBX = Ux2 * diag_GBBX * Ux2t;
                Fxc_GBB = 2.0 * 0.5 * (Fxc_GBBR * Fxc_GBBX + Fxc_GBBX * Fxc_GBBR);
            }

            {
                TlMatrix FxcA_term2 = Fxc_GAA + Fxc_GAB;
                TlMatrix FxcA_term2t = FxcA_term2;
                FxcA_term2t.transpose();
                TlMatrix FxcA_tilde_term2 = FxcA_term2t * RZ2 + RZ2t * FxcA_term2;
                FxcA_tilde += FxcA_tilde_term2;
            }
            {
                TlMatrix FxcB_term2 = Fxc_GBB + Fxc_GAB;
                TlMatrix FxcB_term2t = FxcB_term2;
                FxcB_term2t.transpose();
                TlMatrix FxcB_tilde_term2 = FxcB_term2t * RZ2 + RZ2t * FxcB_term2;
                FxcB_tilde += FxcB_tilde_term2;
            }
        }
        FxcA = S * V * FxcA_tilde * Vt * St;
        FxcB = S * V * FxcB_tilde * Vt * St;
    }
    DfObject::saveFxcMatrix(RUN_RKS, this->m_nIteration, TlSymmetricMatrix(FxcA));

    // Exc -------------------------------------------------------------
    TlSymmetricMatrix ExcA(numOfAOs);
    //TlSymmetricMatrix ExcB(numOfAOs);
    {
        const FunctionalSets fs = pFunc->getFunctional_GF(rhoAs, rhoBs, xAs, xBs);

        TlMatrix ExcA_tilde(numOfGFOrthNormBasis, numOfGFOrthNormBasis);
        //TlMatrix ExcB_tilde(numOfGFOrthNormBasis, numOfGFOrthNormBasis);
        {
            TlSymmetricMatrix diag_AR(numOfGFOrthNormBasis);
            TlSymmetricMatrix diag_AX(numOfGFOrthNormBasis);
            //TlSymmetricMatrix diag_BR(numOfGFOrthNormBasis);
            //TlSymmetricMatrix diag_BX(numOfGFOrthNormBasis);
            const int numOfTerms = pFunc->getNumOfFunctionalTerms();
            for (int term = 0; term < numOfTerms; ++term) {
                for (index_type i = 0; i < numOfGFOrthNormBasis; ++i) {
                    const double rho = rhoAs[i] + rhoBs[i];
                    if (rho > 1.0E-16) {
                        const double inv_rho = 1.0 / rho;
                        diag_AR(i, i) = fs.FA_termR(term, i) * inv_rho;
                        diag_AX(i, i) = fs.FA_termX(term, i);
                        // diag_BR(i, i) = fs.FB_termR(term, i) * inv_rho;
                        // diag_BX(i, i) = fs.FB_termX(term, i);
                    }
                }
                
                const TlSymmetricMatrix ExcA_R = U * diag_AR * Ut;
                const TlSymmetricMatrix ExcA_X = Ux2 * diag_AX * Ux2t;
                const TlMatrix ExcA_tilde_term = ExcA_R * ExcA_X + ExcA_X * ExcA_R;
                // const TlSymmetricMatrix ExcB_R = U * diag_BR * Ut;
                // const TlSymmetricMatrix ExcB_X = Ux2 * diag_BX * Ux2t;
                // const TlMatrix ExcB_tilde_term = ExcB_R * ExcB_X + ExcB_X * ExcB_R;
                
                ExcA_tilde += ExcA_tilde_term;
                // ExcB_tilde += ExcB_tilde_term;
            }
        }

        ExcA = S * V * ExcA_tilde * Vt * St;
        // ExcB = S * V * ExcB_tilde * Vt * St;
        // Exc *= 2.0; // means RKS
    }
    DfObject::saveExcMatrix(RUN_RKS, this->m_nIteration, TlSymmetricMatrix(ExcA));
    // DfObject::saveExcMatrix(RUN_UKS_ALPHA, this->m_nIteration, TlSymmetricMatrix(ExcA));
    // DfObject::saveExcMatrix(RUN_UKS_BETA,  this->m_nIteration, TlSymmetricMatrix(ExcB));
    
    delete pFunc;
    pFunc = NULL;
}

DfFunctional_GGA* DfGridFreeXC::getFunctionalGGA()
{
    DfFunctional_GGA* pFunc = NULL;
    
    DfXCFunctional xcFunc(this->pPdfParam_);
    const DfXCFunctional::XC_TYPE xcType = xcFunc.getXcType();
    switch (xcType) {
    case DfXCFunctional::HFB:
        pFunc = new DfFunctional_Becke88();
        break;

    case DfXCFunctional::BLYP:
        pFunc = new DfFunctional_B88LYP();
        break;

    case DfXCFunctional::B3LYP:
        pFunc = new DfFunctional_B3LYP();
        break;

    default:
        std::cerr << "DfGridFreeXC::getFunctionalGGA() " << xcType << std::endl;
        pFunc = NULL;
        abort();
        break;
    }
    
    return pFunc;
}

