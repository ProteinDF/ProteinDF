#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "DfHpq.h"
#include "DfHpq_driver.h"
#include "TlTime.h"
#include "TlSymmetricMatrix.h"

// 25 = 5 x 5 means D x D.
#define WORK_SIZE 25
#define DEFAULT_BLOCK_SIZE 1000

const double DfHpq::FPAI   = 2.0 * std::pow(M_PI, 2.5);
const double DfHpq::D33    = 1.0 / 0.03;
const double DfHpq::SQR3I  = 1.0 / std::sqrt(3.0);
const double DfHpq::SQR3I2 = 2.0 / std::sqrt(3.0);

const std::size_t DfHpq::TF_SIZE = 25;
const std::size_t DfHpq::RMI_SIZE = 24;
const std::size_t DfHpq::GA_SIZE = 6;
const std::size_t DfHpq::EDAT_SIZE = 1901;
const std::size_t DfHpq::SDAT_SIZE = 1901;
const std::size_t DfHpq::ADAT_SIZE = 60 * 1901;


DfHpq::DfHpq(TlSerializeData* pPdfParam) : DfObject(pPdfParam), pOrbitalInfo_(NULL)
{
    this->pOrbitalInfo_ = new TlOrbitalInfo((*pPdfParam)["coordinates"],
                                            (*pPdfParam)["basis_sets"]);

    const TlSerializeData& data = *(this->pPdfParam_);
    this->cutvalue = data["cut-value"].getDouble();

    this->m_nNumOfTotalAtoms = data["control-number-of-atoms"].getInt();

    this->blockSize_ = DEFAULT_BLOCK_SIZE;
    if (data["block-size"].getStr().empty() != true) {
        this->blockSize_ = data["block-size"].getInt();
    }

    this->getMemory();
    this->totalCounter_.resize(16); // 4 * shelltype(p) + shelltype(q): shelltypeの範囲は0-2(s, p, d)
    this->progressCounter_.resize(16);
    this->cutoffCounter_.resize(16);
}


DfHpq::~DfHpq()
{
    if (this->pOrbitalInfo_ != NULL) {
        delete this->pOrbitalInfo_;
    }

    delete[] ADAT;

    for (int i = 0; i < 7; ++i) {
        delete[] *(FDAT+i);
    }
    delete[] FDAT;
    delete[] SDAT;
    delete[] EDAT;
    delete[] GA;
    delete[] RMI;
    delete[] TF;
}


void DfHpq::getMemory()
{
    TF = new double[TF_SIZE];
    RMI = new double[RMI_SIZE];
    GA = new double[GA_SIZE];
    EDAT = new double[EDAT_SIZE];
    SDAT = new double[SDAT_SIZE];
    FDAT = new double*[7];
    for (int i = 0; i < 7; ++i) {
        *(FDAT+i) = new double[1901];
    }
    this->ADAT = new double[ADAT_SIZE];
}


void DfHpq::getHpq(TlSymmetricMatrix* pHpq, TlSymmetricMatrix* pHpq2)
{
    assert(pHpq != NULL);
    assert(pHpq2 != NULL);
    
    // prepare
    this->loggerTime(" start");
    this->makeTable();

    this->loggerTime(" set aux");
    this->auxSet();

    this->loggerTime(" integral");
    this->resetCounter();
    this->getHpq_core(pHpq, pHpq2);
    if (this->m_nChargeExtrapolateNumber != 0) {
        *pHpq2 /= static_cast<double>(this->m_nChargeExtrapolateNumber);
    }

    this->loggerTime(" finalize");
    this->loggerTime(" end");
}


void DfHpq::getHpq_core(TlMatrixObject* pHpq, TlMatrixObject* pHpq2)
{
    const std::size_t maxNumOfPairs = this->blockSize_;
    std::vector<IJShellPair> ij = this->getQueue(maxNumOfPairs, true);
    while (ij.empty() != true) {
        this->hpqcalc(ij, pHpq, pHpq2);
        ij = this->getQueue(maxNumOfPairs);
    }
}


std::vector<DfObject::IJShellPair>
DfHpq::getQueue(const int maxNumOfPairs, const bool initialize)
{
    static const int MAXTYPE = 3; // s, p, d
    static bool hasCalced = false;
    struct State {
        int ity;
        int jty;
        int i;
        int j;
    };
    static State state;
    if (initialize == true) {
        hasCalced = false;
        state.ity = MAXTYPE -1;
        state.jty = state.ity;
        state.i = -1;
        state.j = -1;
    }

    std::vector<IJShellPair> ijShellPairList;
    ijShellPairList.reserve(maxNumOfPairs);
    if (hasCalced == true) {
        return ijShellPairList;
    }

    int numOfPairs = 0;
    // I-shell Type ------------------------------------------------------
    for (int ity = state.ity; ity >= 0; --ity) {

        // J-shell Type ----------------------------------------------------
        for (int jty = std::min(ity, state.jty); jty >= 0; --jty) {

            for (std::size_t i = std::max(0, state.i); i < this->shellList_[ity].size(); ++i) {
                const std::size_t index_i = this->shellList_[ity][i];

                for (std::size_t j = std::max(0, state.j); j < this->shellList_[jty].size(); ++j) {

                    this->countProgress(ity, jty);
//        if (this->isCutoffUsingSchwartzInequality(i, j) == true) {
//      this->countupCutoff(ity, jty);
//      continue;
//    }

                    const std::size_t index_j = this->shellList_[jty][j];
                    if ((ity == jty) && (index_i < index_j)) {
                        continue;
                    }

                    ijShellPairList.push_back(IJShellPair(index_i,
                                                          index_j));
                    ++numOfPairs;

                    if (numOfPairs >= maxNumOfPairs) {
                        // 現状の保存
                        state.ity = ity;
                        state.jty = jty;
                        state.i = i;
                        state.j = j +1; // "+1" means to start next j-shell at next time.

                        return ijShellPairList;
                    }
                }
                state.j = -1;
            }
            state.i = -1;
        }
        state.jty = MAXTYPE;
    }

    // finalize
    hasCalced = true;
    return ijShellPairList;
}


////////////////////////////////////////////////////////////////////////////////
// counter
//
void DfHpq::resetCounter()
{
    static const int MAXTYPE = 3; // s, p, d
    std::fill(this->totalCounter_.begin(), this->totalCounter_.end(), 0);
    std::fill(this->cutoffCounter_.begin(), this->cutoffCounter_.end(), 0);

    // count total pairs
    for (int i = 0; i < MAXTYPE; ++i) {
        const std::size_t coef_A = this->shellList_[i].size();
        for (int j = 0; j <= i; ++j) {
            const std::size_t coef_B = this->shellList_[j].size();

            const int shellType = 4 * i + j;
            this->totalCounter_[shellType] = coef_A * coef_B;
        }
    }
}


void DfHpq::countProgress(const int p_type, const int q_type)
{
    const int shellType = p_type * 4 + q_type;
    ++(this->progressCounter_[shellType]);
}


void DfHpq::countCutoff(const int p_type, const int q_type)
{
    const int shellType = p_type * 4 + q_type;
    ++(this->cutoffCounter_[shellType]);
}


int DfHpq::getProgress() const
{
    std::size_t sum = 0;
    std::size_t calcd = 0;
    const int max_i = this->progressCounter_.size();
    assert(static_cast<size_t>(max_i) == this->totalCounter_.size());
    for (int i = 0; i < max_i; ++i) {
        sum += this->totalCounter_[i];
        calcd += this->progressCounter_[i];
    }

    if (sum == 0) {
        return -1;
    }

    return (double)calcd / (double)sum * 100.0;
}


void DfHpq::hpqcalc(const std::vector<IJShellPair>& IJShellList,
                    TlMatrixObject* Hpq, TlMatrixObject* Hpq2)
{
    const int workSize = 25; // 5(d) x 5(d)
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);
    const std::size_t numOfAtoms = flGeom.getNumOfAtoms();

    int numOfPairs = IJShellList.size();
//#pragma omp parallel for
    for (int pairIndex = 0; pairIndex < numOfPairs; ++pairIndex) {
        const int orbA = IJShellList[pairIndex].nIShell;
        const int orbB = IJShellList[pairIndex].nJShell;
        const int nqA = this->pOrbitalInfo_->getShellType(orbA);
        const int nqB = this->pOrbitalInfo_->getShellType(orbB);
        const int iAngA = nqA * 2  +1;
        const int iAngB = nqB * 2 + 1;

        const int npA = this->pOrbitalInfo_->getCgtoContraction(orbA);
        assert(npA > 0);
        const TlPosition posA = this->pOrbitalInfo_->getPosition(orbA);
        double A[3];
        A[0] = posA.x();
        A[1] = posA.y();
        A[2] = posA.z();
        double* Ca = new double[npA];
        double* Za = new double[npA];
        for (int i = 0; i < npA; ++i) {
            Ca[i] = this->pOrbitalInfo_->getCoefficient(orbA, i);
            Za[i] = this->pOrbitalInfo_->getExponent(orbA, i);
        }

        const int npB = this->pOrbitalInfo_->getCgtoContraction(orbB);
        assert(npB > 0);
        const TlPosition posB = this->pOrbitalInfo_->getPosition(orbB);
        double B[3];
        B[0] = posB.x();
        B[1] = posB.y();
        B[2] = posB.z();
        double* Cb = new double[npB];
        double* Zb = new double[npB];
        for (int i = 0; i < npB; ++i) {
            Cb[i] = this->pOrbitalInfo_->getCoefficient(orbB, i);
            Zb[i] = this->pOrbitalInfo_->getExponent(orbB, i);
        }

        double* E = new double[workSize];
        double* EW = new double[workSize];

        //-------  KINETIC ENERGY INTEGRALS
        {
            // initialize work region
            for (int i = 0; i < workSize; ++i) {
                E[i] = 0.0;
            }

            const int type = nqA * 4 + nqB;
            switch (type) {
            case 0: // ss
                this->kinSS(npA, npB, Za, Zb, Ca, Cb, A, B, E);
                break;

            case 4: // ps
                this->kinPS(npA, npB, Za, Zb, Ca, Cb, A, B, E);
                break;

            case 5: // pp
                this->kinPP(npA, npB, Za, Zb, Ca, Cb, A, B, E);
                break;

            case 8: // ds
                this->kinDS(npA, npB, Za, Zb, Ca, Cb, A, B, E);
                break;

            case 9: // dp
                this->kinDP(npA, npB, Za, Zb, Ca, Cb, A, B, E);
                break;

            case 10: // dd
                this->kinDD(npA, npB, Za, Zb, Ca, Cb, A, B, E);
                break;

            default:
                abort();
            }

#pragma omp critical (DfHpq_hpqcalc)
            {
                int iww = 0;
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    const std::size_t I = orbA + iAng;

                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        const std::size_t J = orbB + jAng;

                        if ((orbA != orbB) || (I >= J)) {
                            if (fabs(E[iww]) >= cutvalue) {
                                Hpq->add(I, J, E[iww]);
                            }
                        }
                        ++iww;
                    }
                }
            }
        }

        //-------  NUCLEAR ATTRACTION INTEGRALS
        {
            // initialize work region
            for (int i = 0; i < workSize; ++i) {
                E[i] = 0.0;
            }

            for (std::size_t nuC = 0; nuC < numOfAtoms; ++nuC) {
                if (flGeom.getAtom(nuC) == "X") {
                    continue;
                }

                const TlPosition posC = flGeom.getCoordinate(nuC);
                double C[3];
                C[0] = posC.x();
                C[1] = posC.y();
                C[2] = posC.z();

                for (int i = 0; i < workSize; ++i) {
                    EW[i] = 0.0;
                }

                DfHpqDrv_getNuclearAttractionIntegrals(nqA, nqB,
                                                       npA, npB,
                                                       Za, Ca, Zb, Cb,
                                                       A, B, C,
                                                       this->TF, this->ADAT,
                                                       this->RMI, this->GA, this->EDAT,
                                                       EW);
                const double charge = flGeom.getCharge(nuC);
                for (int i = 0; i < WORK_SIZE; ++i) {
                    E[i] -= charge * EW[i];
                }
            }

#pragma omp critical (DfHpq_hpqcalc)
            {
                int iww = 0;
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    const std::size_t I = orbA + iAng;

                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        const std::size_t J = orbB + jAng;

                        if ((orbA != orbB) || (I >= J)) {
                            if (std::fabs(E[iww]) >= this->cutvalue) {
                                Hpq->add(I, J, E[iww]);
                            }
                        }
                        ++iww;
                    }
                }
            }
        }

        // for Hpq2
        if (this->m_nNumOfTotalAtoms - this->m_nNumOfRealAtoms > 0) {
            assert(Hpq2 != NULL);
            for (int i = 0; i < workSize; ++i) {
                E[i] = 0.0;
            }

            const std::size_t numOfAtoms = flGeom.getNumOfAtoms();
            for (std::size_t nuC = 0; nuC < numOfAtoms; ++nuC) {
                if (flGeom.getAtom(nuC) != "X") {
                    continue;
                }

                const TlPosition posC = flGeom.getCoordinate(nuC);
                double C[3];
                C[0] = posC.x();
                C[1] = posC.y();
                C[2] = posC.z();

                for (int i = 0; i < workSize; ++i) {
                    EW[i] = 0.0;
                }

                DfHpqDrv_getNuclearAttractionIntegrals(nqA, nqB,
                                                       npA, npB,
                                                       Za, Ca, Zb, Cb,
                                                       A, B, C,
                                                       this->TF, this->ADAT,
                                                       this->RMI, this->GA, this->EDAT,
                                                       EW);
                const double charge = flGeom.getCharge(nuC);
                int iww = 0;
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        E[iww] -= charge * EW[iww];
                        ++iww;
                    }
                }
            }

#pragma omp critical (DfHpq_hpqcalc)
            {
                int iww = 0;
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    const std::size_t I = orbA + iAng;
                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        const std::size_t J = orbB + jAng;

                        if ((orbA != orbB) || (I >= J)) {
                            if (fabs(E[iww]) >= this->cutvalue) {
                                Hpq2->add(I, J, E[iww]);
                            }
                        }
                        ++iww;
                    }
                }
            }
        }

        delete[] EW;
        EW = NULL;
        delete[] E;
        E = NULL;

        delete[] Ca;
        Ca = NULL;
        delete[] Za;
        Za = NULL;
        delete[] Cb;
        Cb = NULL;
        delete[] Zb;
        Zb = NULL;
    }
}


/*****************************************************************************/
/*   make table for Hpq                                                      */
/*      this->m_nNumOfTotalAtoms  = total number of atoms                    */
/*      ncgto  =                 orbital CGTOs                               */
/*      ndgto  =                 density function GTOs                       */
/*      norbf  =                 orbital basis functions                     */
/*****************************************************************************/
void DfHpq::makeTable()
{
    const std::size_t numOfAOs = this->m_nNumOfAOs;
    this->shellList_.clear();
    this->shellList_.resize(3); // s, p, d support
    {
        std::size_t shell = 0;
        while (shell < numOfAOs) {
            const int shellType = this->pOrbitalInfo_->getShellType(shell);
            this->shellList_[shellType].push_back(shell);
            shell += shellType * 2 + 1;
        }
    }

    // Get Atom Coordinates
    Fl_Geometry fgeom((*this->pPdfParam_)["coordinates"]);
    this->m_nNumOfTotalAtoms = fgeom.getNumOfAtoms();
    this->m_nNumOfRealAtoms  = fgeom.getNumOfAtoms() - fgeom.getDummyatom();
}

/*****************************************************************************/
/*   aux set for ERI [pq|Alpha]                          */
/*  FPAI   = 2*PAI**2.5                          */
/*  D33    = 1/0.03                              */
/*  SQR3I  = 1/sqrt(3)                           */
/*  SQR3I2 = 2*SQR3I                             */
/*****************************************************************************/
int DfHpq::auxSet()
{
    // TF Data
    double TF0[25];
    TF0[ 0]=33.0;
    TF0[ 1]=37.0;
    TF0[ 2]=41.0;
    TF0[ 3]=43.0;
    TF0[ 4]=46.0;
    TF0[ 5]=49.0;
    TF0[ 6]=51.0;
    TF0[ 7]=54.0;
    TF0[ 8]=56.0;
    TF0[ 9]=58.0;
    TF0[10]=61.0;
    TF0[11]=63.0;
    TF0[12]=66.0;
    TF0[13]=68.0;
    TF0[14]=70.0;
    TF0[15]=72.0;
    TF0[16]=74.0;
    TF0[17]=80.0;
    TF0[18]=85.0;
    TF0[19]=90.0;
    TF0[20]=95.0;
    TF0[21]=95.0;
    TF0[22]=98.0;
    TF0[23]=99.0;
    TF0[24]=99.0;

    int nk[6];
    nk[0] = -1;
    nk[2] = 18;
    nk[4] = -48;
    nk[1] = nk[3] = nk[5] = 0;

    const double DD = 0.015;
    double D[6];
    D[1]=DD*DD;
    D[3]=D[1]*D[1];
    D[5]=D[1]*D[3];

    for (int k = 0; k < 25; ++k) {
        TF[k] = TF0[k];
    }

    for (int k = 1; k < 25; ++k) {
        RMI[k-1] = 1.0/(double)(2*k+1);
    }

    GA[0]=1.0;
    for (int k = 1; k < 6; ++k) {
        GA[k] = GA[k-1]/(double)k;
    }

    const double GA5 = GA[5]/192.0; // 192=6*32
    for (int k = 0; k < 6; ++k) {
        if (nk[k] != 0) {
            GA[k] -= nk[k] * GA5 * D[5-k];
        }
    }

    // The following loop may be meaningless, because the same loop exists
    // in member function, fmtSet.
    for (int k = 0; k < 1901; ++k) {
        EDAT[k] = exp(-0.03*k);
    }

    this->fmtSet(); // call fmtSet
    // copy *FDAT[6]  to *SDAT
    for (int k = 0; k < 1901; ++k) {
        SDAT[k] = FDAT[6][k];
    }
    this->fmtRecursive(14, 1901, 5); // call fmtRecursive

    //
    int mm,jm,mk,nkk;
    int jk,jj;
    double gm,dk,gdk;
    double DMI,t2;

    int mp = 2;
    int mnp = 1;

    for (int m = 8; m >= 0; --m) {
        mm=m*6;
        gm=1.0;

        if (m!=8) {
            mp--;
            if (mp ==0)mp =7;
            mnp--;
            if (mnp==0)mnp=7;
            // fmtRecursive(m+1, 1901, 0);
            DMI=1.0/(2.0*m+1.0);
            for (int k = 0; k < 1901; ++k) {
                t2            =0.06*(double) k;
                FDAT[mnp-1][k]=DMI*(t2*FDAT[mp-1][k]+EDAT[k]);
            }
            // fmtRecursive end
        }

        jm = mnp-1;
        if (jm == 0) {
            jm = 7;
        }
        for (int k = 0; k < 6; ++k) {
            mk = mm + k;
            if (k!=0)gm/=(double)k;
            nkk=nk[k];
            dk =D[5-k];
            jk =mnp+k;
            jj =jk/8;  // check here
            jk -= jj*7;// check here
            if (nkk==0) {
                for (int l = 0; l < 1901; ++l) {
                    //ADAT[mk][l] = gm*FDAT[jk-1][l];
                    ADAT[1901 * mk + l] = gm * FDAT[jk-1][l];
                }
                goto SKIP;
            }
            gdk=nkk*dk*GA5;

            for (int l = 0; l < 1901; ++l) {
                //ADAT[mk][l] = gm*FDAT[jk-1][l]-gdk*FDAT[jm-1][l];
                ADAT[1901 * mk +l] = gm * FDAT[jk-1][l] - gdk*FDAT[jm-1][l];
            }

SKIP:
            continue;

        }   // end-of k-for
    } // end-of m-for

    //  for(k=0;k<60;++k)cout<<"ADAT["<<k<<"][0] = "<<ADAT[k][0]<<endl;

    //  for(l=0;l<1901;++l)
    //    for(k=0;k<60;++k)
    //      ADAT[k/10][k%10][l]=ADAT[k][l];
    // check here !!

    return 0;

}

/*****************************************************************************/
/*  FMT calculation table set (m = 14)                   */
/*  Should calculate "int double" (real*16).                 */
/*****************************************************************************/
void DfHpq::fmtSet()
{
    int n1=0, n2=0, n3=0, ms;
    double DMI;
    double T[1901], T2[1901], FMT[1901], Phai[1901];

    for (int k=0; k<1901; ++k) {
        T [k]=0.03*(double)k;
        T2[k]=2*T[k];
        if (T[k]<= 2.0)n1=k;
        if (T[k]<=10.0)n2=k;
        if (T[k]<=20.0)n3=k;
    }

    for (int k=0; k<1901; ++k) {
        EDAT[k]=exp(-T[k]);
        FMT [k]=0.0;
        Phai[k]=0.0;
    }

    for (int m=50; m>=14; m--) {
        DMI=1.0/(2.0*m+1.0);
        for (int k=0; k<n1; ++k) {
            FMT[k]=DMI*(T2[k]*FMT[k]+EDAT[k]);
        }
    }

    for (int m=70; m>=14; m--) {
        DMI=1.0/(2.0*m+1.0);
        for (int k=n1; k<n2; ++k) {
            FMT[k]=DMI*(T2[k]*FMT[k]+EDAT[k]);
        }
    }

    for (int m=120; m>=14; m--) {
        DMI=1.0/(2.0*m+1.0);
        for (int k=n2; k<n3; ++k) {
            FMT[k]=DMI*(T2[k]*FMT[k]+EDAT[k]);
        }
    }

    for (int k=n3; k<1901; ++k) {
        FMT[k]=0.5*sqrt(M_PI/T[k]);
    }

    for (int m=1; m<=14; ++m) {
        for (int k=n3; k<1901; ++k) {
            FMT[k]*=(2.0*m-1.0)/T2[k];
        }
    }

    for (int k=n3; k<1901; ++k) {
        ms      =(int)(T[k]+0.5);
        ms     *=-1;
        Phai[k] =EDAT[k]/(T2[k]-2.0*ms);
        for (int m=ms+1; m<=14; ++m) {
            Phai[k]=((2.0*m-1.0)*Phai[k]+EDAT[k])/T2[k];
        }
    }

    for (int k=n3; k<1901; ++k) {
        FMT[k]-=Phai[k];
    }

    for (int k= 0; k<1901; ++k) {
        FDAT[6][k]=FMT[k]; // 0-5 are work area.
    }

    return;

}

void DfHpq::fmtRecursive(int m, int n, int mm)
{
    /*****************************************************************************/
    /*                                         */
    /*    FMT backward recursive relation                      */
    /*                                         */
    /*****************************************************************************/
    int k, kk, l;
    double DMI,t2;

    //  if(mm<0)
    if (mm<1)goto ERROR_Label;

    DMI=1.0/(2.0*m-1.0);
    for (k=0; k<n; ++k) {
        t2         =0.06*(double)k;
        FDAT[mm][k]=DMI*(t2*SDAT[k]+EDAT[k]);
    }

    //  if(mm==0)
    //  return;

    kk=m-1;
    for (k=1; k<=mm; ++k) {
        kk--;
        DMI=1.0/(2.0*kk+1.0);
        for (l=0; l<n; ++l) {
            t2           =0.06*(double) l;
            FDAT[mm-k][l]=DMI*(t2*FDAT[mm-k+1][l]+EDAT[l]);
        }
    }

    return;

ERROR_Label:
    std::cout<<"### Error in DfEri::fmtRecursive ###"<<std::endl;
    exit(1);

}


/******************************************************************************/
/*   Kinetic energy evaluation for type SS                                    */
/*      npA     ; number of first   primitive pair data                       */
/*      npB     ; number of seconde primitive pair data                       */
/*      Za      ; orbital exponents alpha                                     */
/*      Zb      ; orbital exponents beta                                      */
/*      Ca      ; contracton coefficient A                                    */
/*      Cb      ; contracton coefficient B                                    */
/*      Ax,Ay,Az; nuclear coordinates of atom A                               */
/*      Bx,By,Bz; nuclear coordinates of atom B                               */
/******************************************************************************/
void DfHpq::kinSS(const int npA, const int npB,
                  const double* Za, const double* Zb, const double* Ca, const double* Cb,
                  const double* A, const double* B, double* E)
{
    const double Ax = A[0];
    const double Ay = A[1];
    const double Az = A[2];
    const double Bx = B[0];
    const double By = B[1];
    const double Bz = B[2];

    const double ABx = Ax - Bx;
    const double ABy = Ay - By;
    const double ABz = Az - Bz;
    const double SqAB = ABx*ABx + ABy*ABy + ABz*ABz;

    const double CPAI = M_PI * std::sqrt(M_PI);
    double STS =0.0;
    for (int i = 0; i < npA; ++i) {
        for (int j = 0; j < npB; ++j) {
            const double Zp = Za[i] + Zb[j];
            const double ZpI = 1.0 / Zp;
            const double GGI = Za[i]*Zb[j] * ZpI;
            STS += GGI * (3.0 - 2.0 * GGI * SqAB)
                   * CPAI * ZpI * std::sqrt(ZpI) * exp(-GGI*SqAB)*Ca[i]*Cb[j];
        }
    }
    E[0] = STS;
}

/******************************************************************************/
/*   Kinetic energy evaluation for type PS                                    */
/*      npA     ; number of first   primitive pair data                       */
/*      npB     ; number of seconde primitive pair data                       */
/*      Za      ; orbital exponents alpha                                     */
/*      Zb      ; orbital exponents beta                                      */
/*      Ca      ; contracton coefficient A                                    */
/*      Cb      ; contracton coefficient B                                    */
/*      Ax,Ay,Az; nuclear coordinates of atom A                               */
/*      Bx,By,Bz; nuclear coordinates of atom B                               */
/******************************************************************************/
void DfHpq::kinPS(const int npA, const int npB,
                  const double* Za, const double* Zb, const double* Ca, const double* Cb,
                  const double* A, const double* B, double* E)
{
    const double Ax=A[0];
    const double Ay=A[1];
    const double Az=A[2];
    const double Bx=B[0];
    const double By=B[1];
    const double Bz=B[2];

    const double ABx=Ax-Bx;
    const double ABy=Ay-By;
    const double ABz=Az-Bz;
    const double SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

    const double CPAI=M_PI*sqrt(M_PI);
    double xTS =0.0;
    double yTS =0.0;
    double zTS =0.0;
    for (int i = 0; i < npA; ++i) {
        for (int j = 0; j < npB; ++j) {
            const double Zp  =Za[i]+Zb[j];
            const double ZpI = 1.0 / Zp;
            const double GGI =Za[i]*Zb[j] * ZpI;
            const double SS  =CPAI*ZpI*sqrt(ZpI)*exp(-GGI*SqAB)*Ca[i]*Cb[j];
            const double STS =GGI*(3.0-2.0*GGI*SqAB)*SS;
            const double Px  =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
            const double Py  =(Za[i]*Ay+Zb[j]*By)*ZpI;
            const double Pz  =(Za[i]*Az+Zb[j]*Bz)*ZpI;
            const double PAx =Px-Ax;
            const double PAy =Py-Ay;
            const double PAz =Pz-Az;
            const double GGI2=GGI+GGI;
            const double xS  =PAx*SS;
            const double yS  =PAy*SS;
            const double zS  =PAz*SS;
            xTS += PAx*STS+GGI2*xS;
            yTS += PAy*STS+GGI2*yS;
            zTS += PAz*STS+GGI2*zS;
        }
    }
    E[0]=xTS;
    E[1]=yTS;
    E[2]=zTS;
}

/******************************************************************************/
/*   Kinetic energy evaluation for type PP                                    */
/*      npA     ; number of first   primitive pair data                       */
/*      npB     ; number of seconde primitive pair data                       */
/*      Za      ; orbital exponents alpha                                     */
/*      Zb      ; orbital exponents beta                                      */
/*      Ca      ; contracton coefficient A                                    */
/*      Cb      ; contracton coefficient B                                    */
/*      Ax,Ay,Az; nuclear coordinates of atom A                               */
/*      Bx,By,Bz; nuclear coordinates of atom B                               */
/******************************************************************************/
void DfHpq::kinPP(const int npA, const int npB,
                  const double* Za, const double* Zb, const double* Ca, const double* Cb,
                  const double* A, const double* B, double* E)
{
    const double Ax=A[0];
    const double Ay=A[1];
    const double Az=A[2];
    const double Bx=B[0];
    const double By=B[1];
    const double Bz=B[2];

    const double ABx=Ax-Bx;
    const double ABy=Ay-By;
    const double ABz=Az-Bz;
    const double SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

    const double CPAI=M_PI*sqrt(M_PI);
    double xTx=0.0;
    double xTy=0.0;
    double xTz=0.0;
    double yTx=0.0;
    double yTy=0.0;
    double yTz=0.0;
    double zTx=0.0;
    double zTy=0.0;
    double zTz=0.0;
    for (int i = 0; i < npA; ++i) {
        for (int j = 0; j < npB; ++j) {
            const double Zp  =Za[i]+Zb[j];
            const double ZpI =1.0/Zp;
            const double GGI =Za[i]*Zb[j] * ZpI;
            const double SS  =CPAI*ZpI*sqrt(ZpI)*exp(-GGI*SqAB)*Ca[i]*Cb[j];
            const double STS =GGI*(3.0-2.0*GGI*SqAB)*SS;
            const double Px  =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
            const double Py  =(Za[i]*Ay+Zb[j]*By)*ZpI;
            const double Pz  =(Za[i]*Az+Zb[j]*Bz)*ZpI;
            const double PAx =Px-Ax;
            const double PAy =Py-Ay;
            const double PAz =Pz-Az;
            const double PBx =Px-Bx;
            const double PBy =Py-By;
            const double PBz =Pz-Bz;
            const double GGI2=GGI+GGI;
            const double ZpI5=0.5*ZpI;
            const double xS  =PAx*SS;
            const double yS  =PAy*SS;
            const double zS  =PAz*SS;
            const double xTS =PAx*STS + GGI2*xS;
            const double yTS =PAy*STS + GGI2*yS;
            const double zTS =PAz*STS + GGI2*zS;
            const double xx  =PBx*xS  + ZpI5*SS;
            const double xy  =PBy*xS;
            const double xz  =PBz*xS;
            const double yx  =PBx*yS;
            const double yy  =PBy*yS  + ZpI5*SS;
            const double yz  =PBz*yS;
            const double zx  =PBx*zS;
            const double zy  =PBy*zS;
            const double zz  =PBz*zS  + ZpI5*SS;
            xTx += PBx*xTS + ZpI5*STS + GGI2*xx;
            xTy += PBy*xTS            + GGI2*xy;
            xTz += PBz*xTS            + GGI2*xz;
            yTx += PBx*yTS            + GGI2*yx;
            yTy += PBy*yTS + ZpI5*STS + GGI2*yy;
            yTz += PBz*yTS            + GGI2*yz;
            zTx += PBx*zTS            + GGI2*zx;
            zTy += PBy*zTS            + GGI2*zy;
            zTz += PBz*zTS + ZpI5*STS + GGI2*zz;
        }
    }

    E[0]=xTx;
    E[1]=xTy;
    E[2]=xTz;
    E[3]=yTx;
    E[4]=yTy;
    E[5]=yTz;
    E[6]=zTx;
    E[7]=zTy;
    E[8]=zTz;
}

/******************************************************************************/
/*   Kinetic energy evaluation for type DS                                    */
/*      npA     ; number of first   primitive pair data                       */
/*      npB     ; number of seconde primitive pair data                       */
/*      Za      ; orbital exponents alpha                                     */
/*      Zb      ; orbital exponents beta                                      */
/*      Ca      ; contracton coefficient A                                    */
/*      Cb      ; contracton coefficient B                                    */
/*      Ax,Ay,Az; nuclear coordinates of atom A                               */
/*      Bx,By,Bz; nuclear coordinates of atom B                               */
/******************************************************************************/
void DfHpq::kinDS(const int npA, const int npB,
                  const double* Za, const double* Zb, const double* Ca, const double* Cb,
                  const double* A, const double* B, double* E)
{
    const double Ax=A[0];
    const double Ay=A[1];
    const double Az=A[2];
    const double Bx=B[0];
    const double By=B[1];
    const double Bz=B[2];

    const double ABx=Ax-Bx;
    const double ABy=Ay-By;
    const double ABz=Az-Bz;
    const double SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

    const double CPAI=M_PI*sqrt(M_PI);
    double ATS=0.0;
    double BTS=0.0;
    double CTS=0.0;
    double DTS=0.0;
    double ETS=0.0;
    for (int i = 0; i < npA; ++i) {
        for (int j = 0; j < npB; ++j) {
            const double Zp   =Za[i]+Zb[j];
            const double GGI  =Za[i]*Zb[j]/Zp;
            const double ZpI  =1.0/Zp;
            const double SS   =CPAI*ZpI*sqrt(ZpI)*exp(-GGI*SqAB)*Ca[i]*Cb[j];
            const double STS  =GGI*(3.0-2.0*GGI*SqAB)*SS;
            const double Px   =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
            const double Py   =(Za[i]*Ay+Zb[j]*By)*ZpI;
            const double Pz   =(Za[i]*Az+Zb[j]*Bz)*ZpI;
            const double PAx  =Px-Ax;
            const double PAy  =Py-Ay;
            const double PAz  =Pz-Az;
            const double GGI2 =GGI+GGI;
            // PS
            const double xS   =PAx*SS;
            const double yS   =PAy*SS;
            const double zS   =PAz*SS;
            const double xTS  =PAx*STS+GGI2*xS;
            const double yTS  =PAy*STS+GGI2*yS;
            const double zTS  =PAz*STS+GGI2*zS;
            // DS
            const double xxSW =PAx*xS;
            const double yySW =PAy*yS;
            const double zzSW =PAz*zS;
            const double AS0  =PAx*yS;
            const double BS0  =PAx*zS;
            const double CS0  =PAy*zS;
            const double xxTSW=PAx*xTS+GGI2*xxSW;
            const double yyTSW=PAy*yTS+GGI2*yySW;
            const double zzTSW=PAz*zTS+GGI2*zzSW;
            ATS += PAx*yTS+GGI2*AS0;
            BTS += PAx*zTS+GGI2*BS0;
            CTS += PAy*zTS+GGI2*CS0;
            DTS += 0.5*(xxTSW - yyTSW);
            ETS += SQR3I*(zzTSW - 0.5*(xxTSW + yyTSW));
        }
    }

    E[0]=ATS;
    E[1]=BTS;
    E[2]=CTS;
    E[3]=DTS;
    E[4]=ETS;
}

/******************************************************************************/
/*   Kinetic energy evaluation for type DP                                    */
/*      npA     ; number of first   primitive pair data                       */
/*      npB     ; number of seconde primitive pair data                       */
/*      Za      ; orbital exponents alpha                                     */
/*      Zb      ; orbital exponents beta                                      */
/*      Ca      ; contracton coefficient A                                    */
/*      Cb      ; contracton coefficient B                                    */
/*      Ax,Ay,Az; nuclear coordinates of atom A                               */
/*      Bx,By,Bz; nuclear coordinates of atom B                               */
/******************************************************************************/
void DfHpq::kinDP(const int npA, const int npB,
                  const double* Za, const double* Zb, const double* Ca, const double* Cb,
                  const double* A, const double* B, double* E)
{
    const double Ax=A[0];
    const double Ay=A[1];
    const double Az=A[2];
    const double Bx=B[0];
    const double By=B[1];
    const double Bz=B[2];

    const double ABx=Ax-Bx;
    const double ABy=Ay-By;
    const double ABz=Az-Bz;
    const double SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

    const double CPAI=M_PI*sqrt(M_PI);
    double ATx=0.0;
    double ATy=0.0;
    double ATz=0.0;
    double BTx=0.0;
    double BTy=0.0;
    double BTz=0.0;
    double CTx=0.0;
    double CTy=0.0;
    double CTz=0.0;
    double DTx=0.0;
    double DTy=0.0;
    double DTz=0.0;
    double ETx=0.0;
    double ETy=0.0;
    double ETz=0.0;
    for (int i = 0; i <npA; ++i) {
        for (int j = 0; j < npB; ++j) {
            const double Zp  =Za[i]+Zb[j];
            const double GGI =Za[i]*Zb[j]/Zp;
            const double ZpI =1.0/Zp;
            const double SS  =CPAI*ZpI*sqrt(ZpI)*exp(-GGI*SqAB)*Ca[i]*Cb[j];
            const double STS =GGI*(3.0-2.0*GGI*SqAB)*SS;
            const double Px  =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
            const double Py  =(Za[i]*Ay+Zb[j]*By)*ZpI;
            const double Pz  =(Za[i]*Az+Zb[j]*Bz)*ZpI;
            const double PAx =Px-Ax;
            const double PAy =Py-Ay;
            const double PAz =Pz-Az;
            const double PBx =Px-Bx;
            const double PBy =Py-By;
            const double PBz =Pz-Bz;
            const double GGI2=GGI+GGI;
            const double ZpI5=0.5*ZpI;
            // PS
            const double xS  =PAx*SS;
            const double yS  =PAy*SS;
            const double zS  =PAz*SS;
            const double xTS =PAx*STS+GGI2*xS;
            const double yTS =PAy*STS+GGI2*yS;
            const double zTS =PAz*STS+GGI2*zS;
            // DS
            const double xxSW =PAx*xS;
            const double yySW =PAy*yS;
            const double zzSW =PAz*zS;
            const double AS0  =PAx*yS;
            const double BS0  =PAx*zS;
            const double CS0  =PAy*zS;
            const double DS0  =0.5*(xxSW-yySW);
            const double ES0  =SQR3I*(zzSW-0.5*(xxSW+yySW));
            const double xxTSW=PAx*xTS+GGI2*xxSW;
            const double yyTSW=PAy*yTS+GGI2*yySW;
            const double zzTSW=PAz*zTS+GGI2*zzSW;
            const double ATS  =PAx*yTS+GGI2*AS0;
            const double BTS  =PAx*zTS+GGI2*BS0;
            const double CTS  =PAy*zTS+GGI2*CS0;
            const double DTS  =0.5*(xxTSW-yyTSW);
            const double ETS  =SQR3I*(zzTSW-0.5*(xxTSW+yyTSW));
            // DP
            const double ZpxS =ZpI5*xS;
            const double ZpyS =ZpI5*yS;
            const double ZpzS =ZpI5*zS;
            const double Ax0  =PBx*AS0+ZpyS;
            const double Ay0  =PBy*AS0+ZpxS;
            const double Az0  =PBz*AS0;
            const double Bx0  =PBx*BS0+ZpzS;
            const double By0  =PBy*BS0;
            const double Bz0  =PBz*BS0+ZpxS;
            const double Cx0  =PBx*CS0;
            const double Cy0  =PBy*CS0+ZpzS;
            const double Cz0  =PBz*CS0+ZpyS;
            const double Dx0  =PBx*DS0+ZpxS;
            const double Dy0  =PBy*DS0-ZpyS;
            const double Dz0  =PBz*DS0;
            const double Ex0  =PBx*ES0-ZpxS*SQR3I;
            const double Ey0  =PBy*ES0-ZpyS*SQR3I;
            const double Ez0  =PBz*ES0+ZpzS*SQR3I2;
            const double ZpxTS=ZpI5*xTS;
            const double ZpyTS=ZpI5*yTS;
            const double ZpzTS=ZpI5*zTS;
            ATx += PBx*ATS+ZpyTS       +GGI2*Ax0;
            ATy += PBy*ATS+ZpxTS       +GGI2*Ay0;
            ATz += PBz*ATS             +GGI2*Az0;
            BTx += PBx*BTS+ZpzTS       +GGI2*Bx0;
            BTy += PBy*BTS             +GGI2*By0;
            BTz += PBz*BTS+ZpxTS       +GGI2*Bz0;
            CTx += PBx*CTS             +GGI2*Cx0;
            CTy += PBy*CTS+ZpzTS       +GGI2*Cy0;
            CTz += PBz*CTS+ZpyTS       +GGI2*Cz0;
            DTx += PBx*DTS+ZpxTS       +GGI2*Dx0;
            DTy += PBy*DTS-ZpyTS       +GGI2*Dy0;
            DTz += PBz*DTS             +GGI2*Dz0;
            ETx += PBx*ETS-SQR3I *ZpxTS+GGI2*Ex0;
            ETy += PBy*ETS-SQR3I *ZpyTS+GGI2*Ey0;
            ETz += PBz*ETS+SQR3I2*ZpzTS+GGI2*Ez0;
        }
    }

    E[ 0]=ATx;
    E[ 1]=ATy;
    E[ 2]=ATz;
    E[ 3]=BTx;
    E[ 4]=BTy;
    E[ 5]=BTz;
    E[ 6]=CTx;
    E[ 7]=CTy;
    E[ 8]=CTz;
    E[ 9]=DTx;
    E[10]=DTy;
    E[11]=DTz;
    E[12]=ETx;
    E[13]=ETy;
    E[14]=ETz;
}

/******************************************************************************/
/*   Kinetic energy evaluation for type DD                                    */
/*      npA     ; number of first   primitive pair data                       */
/*      npB     ; number of seconde primitive pair data                       */
/*      Za      ; orbital exponents alpha                                     */
/*      Zb      ; orbital exponents beta                                      */
/*      Ca      ; contracton coefficient A                                    */
/*      Cb      ; contracton coefficient B                                    */
/*      Ax,Ay,Az; nuclear coordinates of atom A                               */
/*      Bx,By,Bz; nuclear coordinates of atom B                               */
/******************************************************************************/
void DfHpq::kinDD(const int npA, const int npB,
                  const double* Za, const double* Zb, const double* Ca, const double* Cb,
                  const double* A, const double* B, double* E)
{
    const double Ax=A[0];
    const double Ay=A[1];
    const double Az=A[2];
    const double Bx=B[0];
    const double By=B[1];
    const double Bz=B[2];

    const double ABx=Ax-Bx;
    const double ABy=Ay-By;
    const double ABz=Az-Bz;
    const double SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

    const double CPAI=M_PI*sqrt(M_PI);
    double ATA=0.0;
    double ATB=0.0;
    double ATC=0.0;
    double ATD=0.0;
    double ATE=0.0;
    double BTA=0.0;
    double BTB=0.0;
    double BTC=0.0;
    double BTD=0.0;
    double BTE=0.0;
    double CTA=0.0;
    double CTB=0.0;
    double CTC=0.0;
    double CTD=0.0;
    double CTE=0.0;
    double DTA=0.0;
    double DTB=0.0;
    double DTC=0.0;
    double DTD=0.0;
    double DTE=0.0;
    double ETA=0.0;
    double ETB=0.0;
    double ETC=0.0;
    double ETD=0.0;
    double ETE=0.0;
    for (int i = 0; i < npA; ++i) {
        for (int j = 0; j < npB; ++j) {
            const double Zp  =Za[i]+Zb[j];
            const double ZpI =1.0/Zp;
            const double GGI =Za[i]*Zb[j] * ZpI;
            const double SS  =CPAI*ZpI*sqrt(ZpI)*exp(-GGI*SqAB)*Ca[i]*Cb[j];
            const double STS =GGI*(3.0-2.0*GGI*SqAB)*SS;
            const double Px  =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
            const double Py  =(Za[i]*Ay+Zb[j]*By)*ZpI;
            const double Pz  =(Za[i]*Az+Zb[j]*Bz)*ZpI;
            const double PAx =Px-Ax;
            const double PAy =Py-Ay;
            const double PAz =Pz-Az;
            const double PBx =Px-Bx;
            const double PBy =Py-By;
            const double PBz =Pz-Bz;
            const double GGI2=GGI+GGI;
            const double ZpI5=0.5*ZpI;
            // PS
            const double xS =PAx*SS;
            const double yS =PAy*SS;
            const double zS =PAz*SS;
            const double xTS=PAx*STS+GGI2*xS;
            const double yTS=PAy*STS+GGI2*yS;
            const double zTS=PAz*STS+GGI2*zS;
            // PP
            const double ZpSS =ZpI5*SS;
            const double ZpSTS=ZpI5*STS;
            const double xx0  =PBx*xS+ZpSS;
            const double xy0  =PBy*xS;
            const double xz0  =PBz*xS;
            const double yx0  =PBx*yS;
            const double yy0  =PBy*yS+ZpSS;
            const double yz0  =PBz*yS;
            const double zx0  =PBx*zS;
            const double zy0  =PBy*zS;
            const double zz0  =PBz*zS+ZpSS;
            const double xTx  =PBx*xTS+ZpSTS+GGI2*xx0;
            const double xTy  =PBy*xTS      +GGI2*xy0;
            const double xTz  =PBz*xTS      +GGI2*xz0;
            const double yTx  =PBx*yTS      +GGI2*yx0;
            const double yTy  =PBy*yTS+ZpSTS+GGI2*yy0;
            const double yTz  =PBz*yTS      +GGI2*yz0;
            const double zTx  =PBx*zTS      +GGI2*zx0;
            const double zTy  =PBy*zTS      +GGI2*zy0;
            const double zTz  =PBz*zTS+ZpSTS+GGI2*zz0;
            // DS
            const double xxSW =PAx*xS;
            const double yySW =PAy*yS;
            const double zzSW =PAz*zS;
            const double AS0  =PAx*yS;
            const double BS0  =PAx*zS;
            const double CS0  =PAy*zS;
            const double DS0  =0.5*(xxSW-yySW);
            const double ES0  =SQR3I*(zzSW-0.5*(xxSW+yySW));
            const double xxTSW=PAx*xTS+GGI2*xxSW;
            const double yyTSW=PAy*yTS+GGI2*yySW;
            const double zzTSW=PAz*zTS+GGI2*zzSW;
            const double ATS  =PAx*yTS+GGI2*AS0;
            const double BTS  =PAx*zTS+GGI2*BS0;
            const double CTS  =PAy*zTS+GGI2*CS0;
            const double DTS  =0.5*(xxTSW-yyTSW);
            const double ETS  =SQR3I*(zzTSW-0.5*(xxTSW+yyTSW));
            // DP
            const double ZpxS =ZpI5*xS;
            const double ZpyS =ZpI5*yS;
            const double ZpzS =ZpI5*zS;
            const double Ax0  =PBx*AS0+ZpyS;
            const double Ay0  =PBy*AS0+ZpxS;
            const double Az0  =PBz*AS0;
            const double Bx0  =PBx*BS0+ZpzS;
            const double By0  =PBy*BS0;
            const double Bz0  =PBz*BS0+ZpxS;
            const double Cx0  =PBx*CS0;
            const double Cy0  =PBy*CS0+ZpzS;
            const double Cz0  =PBz*CS0+ZpyS;
            const double Dx0  =PBx*DS0+ZpxS;
            const double Dy0  =PBy*DS0-ZpyS;
            const double Dz0  =PBz*DS0;
            const double Ex0  =PBx*ES0-ZpxS*SQR3I;
            const double Ey0  =PBy*ES0-ZpyS*SQR3I;
            const double Ez0  =PBz*ES0+ZpzS*SQR3I2;
            const double ZpxTS=ZpI5*xTS;
            const double ZpyTS=ZpI5*yTS;
            const double ZpzTS=ZpI5*zTS;
            const double ATx  =PBx*ATS+ZpyTS        +GGI2*Ax0;
            const double ATy  =PBy*ATS+ZpxTS        +GGI2*Ay0;
            const double ATz  =PBz*ATS              +GGI2*Az0;
            const double BTx  =PBx*BTS+ZpzTS        +GGI2*Bx0;
            const double BTy  =PBy*BTS              +GGI2*By0;
            const double BTz  =PBz*BTS+ZpxTS        +GGI2*Bz0;
            const double CTx  =PBx*CTS              +GGI2*Cx0;
            const double CTy  =PBy*CTS+ZpzTS        +GGI2*Cy0;
            const double CTz  =PBz*CTS+ZpyTS        +GGI2*Cz0;
            const double DTx  =PBx*DTS+ZpxTS        +GGI2*Dx0;
            const double DTy  =PBy*DTS-ZpyTS        +GGI2*Dy0;
            const double DTz  =PBz*DTS              +GGI2*Dz0;
            const double ETx  =PBx*ETS-SQR3I *ZpxTS +GGI2*Ex0;
            const double ETy  =PBy*ETS-SQR3I *ZpyTS +GGI2*Ey0;
            const double ETz  =PBz*ETS+SQR3I2*ZpzTS +GGI2*Ez0;
            // DD-OV
            const double Zpxx0=ZpI5*xx0;
            const double Zpxy0=ZpI5*xy0;
            const double Zpxz0=ZpI5*xz0;
            const double Zpyx0=ZpI5*yx0;
            const double Zpyy0=ZpI5*yy0;
            const double Zpyz0=ZpI5*yz0;
            const double Zpzx0=ZpI5*zx0;
            const double Zpzy0=ZpI5*zy0;
            const double Zpzz0=ZpI5*zz0;
            const double AxxW =PBx*Ax0+Zpyx0;
            const double AyyW =PBy*Ay0+Zpxy0;
            const double AzzW =PBz*Az0;
            const double BxxW =PBx*Bx0+Zpzx0;
            const double ByyW =PBy*By0;
            const double BzzW =PBz*Bz0+Zpxz0;
            const double CxxW =PBx*Cx0;
            const double CyyW =PBy*Cy0+Zpzy0;
            const double CzzW =PBz*Cz0+Zpyz0;
            const double DxxW =PBx*Dx0+Zpxx0;
            const double DyyW =PBy*Dy0-Zpyy0;
            const double DzzW =PBz*Dz0;
            const double ExxW =PBx*Ex0-Zpxx0*SQR3I;
            const double EyyW =PBy*Ey0-Zpyy0*SQR3I;
            const double EzzW =PBz*Ez0+Zpzz0*SQR3I2;
            const double AA0  =PBx*Ay0+Zpyy0;
            const double AB0  =PBz*Ax0;
            const double AC0  =PBz*Ay0;
            const double BA0  =PBy*Bx0;
            const double BB0  =PBx*Bz0+Zpzz0;
            const double BC0  =PBy*Bz0;
            const double CA0  =PBx*Cy0;
            const double CB0  =PBx*Cz0;
            const double CC0  =PBy*Cz0+Zpzz0;
            const double DA0  =PBx*Dy0+Zpxy0;
            const double DB0  =PBx*Dz0+Zpxz0;
            const double DC0  =PBy*Dz0-Zpyz0;
            const double EA0  =PBx*Ey0-Zpxy0*SQR3I;
            const double EB0  =PBx*Ez0-Zpxz0*SQR3I;
            const double EC0  =PBy*Ez0-Zpyz0*SQR3I;
            // DD-KINETIC
            const double ZpxTx=ZpI5*xTx;
            const double ZpxTy=ZpI5*xTy;
            const double ZpxTz=ZpI5*xTz;
            const double ZpyTx=ZpI5*yTx;
            const double ZpyTy=ZpI5*yTy;
            const double ZpyTz=ZpI5*yTz;
            const double ZpzTx=ZpI5*zTx;
            const double ZpzTy=ZpI5*zTy;
            const double ZpzTz=ZpI5*zTz;
            const double ATxxW=PBx*ATx+GGI2*AxxW+ZpyTx;
            const double ATyyW=PBy*ATy+GGI2*AyyW+ZpxTy;
            const double ATzzW=PBz*ATz+GGI2*AzzW;
            const double BTxxW=PBx*BTx+GGI2*BxxW+ZpzTx;
            const double BTyyW=PBy*BTy+GGI2*ByyW;
            const double BTzzW=PBz*BTz+GGI2*BzzW+ZpxTz;
            const double CTxxW=PBx*CTx+GGI2*CxxW;
            const double CTyyW=PBy*CTy+GGI2*CyyW+ZpzTy;
            const double CTzzW=PBz*CTz+GGI2*CzzW+ZpyTz;
            const double DTxxW=PBx*DTx+GGI2*DxxW+ZpxTx;
            const double DTyyW=PBy*DTy+GGI2*DyyW-ZpyTy;
            const double DTzzW=PBz*DTz+GGI2*DzzW;
            const double ETxxW=PBx*ETx+GGI2*ExxW-ZpxTx*SQR3I;
            const double ETyyW=PBy*ETy+GGI2*EyyW-ZpyTy*SQR3I;
            const double ETzzW=PBz*ETz+GGI2*EzzW+ZpzTz*SQR3I2;
            ATA += PBx*ATy+GGI2*AA0+ZpyTy;
            ATB += PBz*ATx+GGI2*AB0;
            ATC += PBz*ATy+GGI2*AC0;
            ATD += 0.5*(ATxxW-ATyyW);
            ATE += SQR3I*(ATzzW-0.5*(ATxxW+ATyyW));
            BTA += PBy*BTx+GGI2*BA0;
            BTB += PBz*BTx+GGI2*BB0+ZpxTx;
            BTC += PBy*BTz+GGI2*BC0;
            BTD += 0.5*(BTxxW-BTyyW);
            BTE += SQR3I*(BTzzW-0.5*(BTxxW+BTyyW));
            CTA += PBx*CTy+GGI2*CA0;
            CTB += PBx*CTz+GGI2*CB0;
            CTC += PBy*CTz+GGI2*CC0+ZpzTz;
            CTD += 0.5*(CTxxW-CTyyW);
            CTE += SQR3I*(CTzzW-0.5*(CTxxW+CTyyW));
            DTA += PBx*DTy+GGI2*DA0+ZpxTy;
            DTB += PBz*DTx+GGI2*DB0;
            DTC += PBz*DTy+GGI2*DC0;
            DTD += 0.5*(DTxxW-DTyyW);
            DTE += SQR3I*(DTzzW-0.5*(DTxxW+DTyyW));
            ETA += PBx*ETy+GGI2*EA0-ZpxTy*SQR3I;
            ETB += PBx*ETz+GGI2*EB0-ZpxTz*SQR3I;
            ETC += PBy*ETz+GGI2*EC0-ZpyTz*SQR3I;
            ETD += 0.5*(ETxxW-ETyyW);
            ETE += SQR3I*(ETzzW-0.5*(ETxxW+ETyyW));
        }
    }

    E[ 0]=ATA;
    E[ 1]=ATB;
    E[ 2]=ATC;
    E[ 3]=ATD;
    E[ 4]=ATE;
    E[ 5]=BTA;
    E[ 6]=BTB;
    E[ 7]=BTC;
    E[ 8]=BTD;
    E[ 9]=BTE;
    E[10]=CTA;
    E[11]=CTB;
    E[12]=CTC;
    E[13]=CTD;
    E[14]=CTE;
    E[15]=DTA;
    E[16]=DTB;
    E[17]=DTC;
    E[18]=DTD;
    E[19]=DTE;
    E[20]=ETA;
    E[21]=ETB;
    E[22]=ETC;
    E[23]=ETD;
    E[24]=ETE;
}


std::vector<double> DfHpq::getESP(const TlMatrixObject* pPpq,
                                  const std::vector<TlPosition>& grids)
{
    // prepare
    this->makeTable();
    this->auxSet();

    std::vector<double> values(grids.size());
    const std::size_t maxNumOfPairs = this->blockSize_;
    std::vector<IJShellPair> ij = this->getQueue(maxNumOfPairs, true);
    while (ij.empty() != true) {
        std::vector<double> tmp_values = this->getESP_core(ij, pPpq, grids);
        std::transform(values.begin(), values.end(), tmp_values.begin(), values.begin(), std::plus<double>());
        ij = this->getQueue(maxNumOfPairs);
    }

    return values;
}


std::vector<double> DfHpq::getESP_core(const std::vector<IJShellPair>& IJShellList,
                                       const TlMatrixObject* pPpq,                                       
                                       const std::vector<TlPosition>& grids)
{
    const int numOfPairs = IJShellList.size();
    const std::size_t numOfGrids = grids.size();
    std::vector<double> values(numOfGrids, 0.0);
    
#pragma omp parallel for schedule(runtime)
    for (int pairIndex = 0; pairIndex < numOfPairs; ++pairIndex) {
        const std::size_t orbA = IJShellList[pairIndex].nIShell;
        const std::size_t orbB = IJShellList[pairIndex].nJShell;
        const int nqA = this->pOrbitalInfo_->getShellType(orbA);
        const int nqB = this->pOrbitalInfo_->getShellType(orbB);
        const int iAngA = nqA * 2 + 1;
        const int iAngB = nqB * 2 + 1;
        const int npA = this->pOrbitalInfo_->getCgtoContraction(orbA);
        const int npB = this->pOrbitalInfo_->getCgtoContraction(orbB);

        const TlPosition posA = this->pOrbitalInfo_->getPosition(orbA);
        double A[3];
        A[0] = posA.x();
        A[1] = posA.y();
        A[2] = posA.z();
        double* Ca = new double[npA];
        double* Za = new double[npA];
        for (int i = 0; i < npA; ++i) {
            Ca[i] = this->pOrbitalInfo_->getCoefficient(orbA, i);
            Za[i] = this->pOrbitalInfo_->getExponent(orbA, i);
        }

        const TlPosition posB = this->pOrbitalInfo_->getPosition(orbB);
        double B[3];
        B[0] = posB.x();
        B[1] = posB.y();
        B[2] = posB.z();
        double* Cb = new double[npB];
        double* Zb = new double[npB];
        for (int i = 0; i < npB; ++i) {
            Cb[i] = this->pOrbitalInfo_->getCoefficient(orbB, i);
            Zb[i] = this->pOrbitalInfo_->getExponent(orbB, i);
        }

        // cutoff
//     double cutvalue = 1.0E-10;
//     {
//       double value = 0.0;
//       const int max_i = this->pOrbitalInfo_->getShellType(orbA);
//       const int max_j = this->pOrbitalInfo_->getShellType(orbB);
//       for (int i = 0; i < max_i; ++i) {
//  const std::size_t I = orbA + i;
//  for (int j = 0; j < max_j; ++j) {
//    const std::size_t J = orbB + j;
//    value = std::max(value, std::fabs(pPpq->get(I, J)));
//  }
//       }
//       if (value < 1.0E-10) {
//  continue;
//       }
//     }
//     if (this->isCutoffUsingOverlap(orbA, orbB, cutvalue) == true) {
//       continue;
//     }

        //-------  NUCLEAR ATTRACTION INTEGRALS
        double EW[WORK_SIZE];

        double C[3];
        for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
            C[0] = grids[gridIndex].x();
            C[1] = grids[gridIndex].y();
            C[2] = grids[gridIndex].z();

            for (int i = 0; i < WORK_SIZE; ++i) {
                EW[i] = 0.0;
            }

            DfHpqDrv_getNuclearAttractionIntegrals(nqA, nqB,
                                                   npA, npB,
                                                   Za, Ca, Zb, Cb,
                                                   A, B, C,
                                                   this->TF, this->ADAT,
                                                   this->RMI, this->GA, this->EDAT,
                                                   EW);
            
            double value = 0.0;
            int iww = 0;
            if ((orbA != orbB)) {
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    const std::size_t I = orbA + iAng;
                    
                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        const std::size_t J = orbB + jAng;
                        
                        double tmp = pPpq->get(I, J) * EW[iww];
                        value -= 2.0 * tmp;
                        ++iww;
                    }
                }
            } else {
                for (int iAng = 0; iAng < iAngA; ++iAng) {
                    const std::size_t I = orbA + iAng;
                    
                    for (int jAng = 0; jAng < iAngB; ++jAng) {
                        const std::size_t J = orbB + jAng;
                        
                        double tmp = pPpq->get(I, J) * EW[iww];
                        if (I > J) {
                            value -= 2.0 * tmp;
                        } else if (I == J) {
                            value -= tmp;
                        } // nothing to do in the case of (I < J)

                        ++iww;
                    }
                }
            }

            
//#pragma omp critical (DfHpq__getESP_core)
#pragma omp atomic
            values[gridIndex] += value;
        }

        delete[] Ca;
        Ca = NULL;
        delete[] Za;
        Za = NULL;
        delete[] Cb;
        Cb = NULL;
        delete[] Zb;
        Zb = NULL;
    }

    return values;
}


bool DfHpq::isCutoffUsingOverlap(const std::size_t orbA, const std::size_t orbB,
                                 const double cutvalue) const
{
    bool bAnswer = false;

    const TlPosition A = this->pOrbitalInfo_->getPosition(orbA);
    const TlPosition B = this->pOrbitalInfo_->getPosition(orbB);
    const double absQ = A.squareDistanceFrom(B);
//   static const double CPAI = M_PI * sqrt(M_PI);

    const double zetaA = this->pOrbitalInfo_->getMinExponent(orbA);
    const double zetaB = this->pOrbitalInfo_->getMinExponent(orbB);
    const double zetaP = zetaA + zetaB;
    const double invZetaP = 1.0 / zetaP;
    //double STS = CPAI * invZetaP * std::sqrt(invZetaP) * std::exp(-absQ * zetaA * zetaB * invZetaP);
    double STS = -absQ * zetaA * zetaB * invZetaP;

    if (STS < -30) {
        bAnswer = true;
    }

    return bAnswer;
}


