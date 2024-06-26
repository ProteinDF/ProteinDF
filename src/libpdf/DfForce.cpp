// Copyright (C) 2002-2014 The ProteinDF project
// see also AUTHORS and README.
//
// This file is part of ProteinDF.
//
// ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

#include "DfForce.h"

#include "DfCalcGridX.h"
#include "DfEriX.h"
#include "DfFunctional_B3LYP.h"
#include "DfFunctional_B88LYP.h"
#include "DfFunctional_SVWN.h"
#include "DfGridFreeXC.h"
#include "DfHpqX.h"
#include "DfOverlapX.h"
#include "DfXCFunctional.h"
#include "Fl_Geometry.h"
#include "TlMsgPack.h"
#include "TlTime.h"

DfForce::DfForce(TlSerializeData* pPdfParam)
    : DfObject(pPdfParam),
      orbitalInfo_((*pPdfParam)["coordinates"], (*pPdfParam)["basis_set"]),
      orbitalInfoDens_((*pPdfParam)["coordinates"],
                       (*pPdfParam)["basis_set_j"]),
      pdfParamForForce_(*pPdfParam) {
    // initialize
    this->force_.resize(this->m_nNumOfAtoms, 3);

    if ((*pPdfParam)["force_cut_value"].getStr().empty() != true) {
        this->pdfParamForForce_["cut_value"] =
            (*pPdfParam)["force_cut_value"].getDouble();
    }

    // debug
    this->isDebugOutMatrix_ = false;
    if ((*pPdfParam)["debug/save_forces"].getStr().empty() != true) {
        this->isDebugOutMatrix_ =
            (*pPdfParam)["debug/save_forces"].getBoolean();
    }
}

DfForce::~DfForce() {}

void DfForce::calcForce() {
    RUN_TYPE runType = RUN_RKS;
    const int iteration = this->m_nIteration;

    this->calcForceFromNuclei();
    this->calcForceFromWS(runType);

    TlDenseSymmetricMatrix_Lapack P =
        this->getPInMatrix<TlDenseSymmetricMatrix_Lapack>(runType, iteration);
    this->calcForceFromHpq(P);

    this->calcForceFromCoulomb(runType);

    this->calcForceFromK(runType);

    if (this->isDFT_ == true) {
        switch (this->XC_engine_) {
            case XC_ENGINE_GRID:
                this->calcForceFromPureXC(runType);
                break;

            case XC_ENGINE_GRIDFREE:
            case XC_ENGINE_GRIDFREE_CD:
                this->calcForceFromPureXC_gridfree(runType);
                break;

            default:
                this->log_.critical("program error.");
                break;
        }
    }

    this->force_ *= -1.0;

    this->output();
    this->saveForce();
}

void DfForce::output() {
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);
    const int numOfAtoms = this->m_nNumOfAtoms;

    // MAX element
    const double max_val = this->force_.getMaxAbsoluteElement();

    // calc RMS
    double rms = 0;
    {
        TlDenseGeneralMatrix_Lapack force2 = this->force_;
        force2.dotInPlace(force2);
        rms = force2.sum();
        rms = std::sqrt(rms / (double(numOfAtoms) * 3.0));
    }

    // output for log
    {
        this->log_.info("=== FORCE ===");

        for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
            const TlAtom atom = flGeom.getAtom(atomIndex);
            this->log_.info(TlUtils::format(
                "%4d:[%-2s] % f, % f, % f", atomIndex, atom.getSymbol().c_str(),
                this->force_.get(atomIndex, 0), this->force_.get(atomIndex, 1),
                this->force_.get(atomIndex, 2)));
        }
        this->log_.info(TlUtils::format("MAX force: % f", max_val));
        this->log_.info(TlUtils::format("RMS force: % f", rms));
        this->log_.info("=============");
    }
    {
        this->log_.info("=== FORCE (without X) ===");
        TlDenseGeneralMatrix_Lapack force_woX = this->force_;
        force_woX -= this->force_Xonly_;

        double max_val_woX = force_woX.getMaxAbsoluteElement();
        TlDenseGeneralMatrix_Lapack force_woX2 = force_woX;
        force_woX2.dotInPlace(force_woX2);
        double rms_woX = force_woX2.sum();
        rms_woX = std::sqrt(rms_woX / (double(numOfAtoms) * 3.0));

        for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
            const TlAtom atom = flGeom.getAtom(atomIndex);
            this->log_.info(TlUtils::format(
                "%4d:[%-2s] % f, % f, % f", atomIndex, atom.getSymbol().c_str(),
                force_woX.get(atomIndex, 0), force_woX.get(atomIndex, 1),
                force_woX.get(atomIndex, 2)));
        }
        this->log_.info(TlUtils::format("MAX force: % f", max_val_woX));
        this->log_.info(TlUtils::format("RMS force: % f", rms_woX));
        this->log_.info("=============");
    }

    // output for pdfparam
    {
        TlSerializeData force_dat;
        force_dat.resize(numOfAtoms);
        for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
            TlSerializeData v;
            v.resize(3);
            v[0] = this->force_.get(atomIndex, 0);
            v[1] = this->force_.get(atomIndex, 1);
            v[2] = this->force_.get(atomIndex, 2);
            force_dat.setAt(atomIndex, v);
        }
        (*this->pPdfParam_)["force"] = force_dat;

        (*this->pPdfParam_)["force_max"] = max_val;
        (*this->pPdfParam_)["force_rms"] = rms;
    }
}

void DfForce::saveForce() {
    this->force_.save("force.mtx");
}

void DfForce::outputStartTitle(const std::string& stepName,
                               const char lineChar) {
    const std::string timeString = TlUtils::format(
        "[%s %s]", TlTime::getNowDate().c_str(), TlTime::getNowTime().c_str());

    std::string title = ">>>> " + stepName + " ";
    TlUtils::pad(title, (72 - timeString.length()), ' ');
    title += timeString;

    // 出力
    this->logger(title + "\n");
}

void DfForce::outputEndTitle(const std::string& stepName, const char lineChar) {
    const std::string timeString = TlUtils::format(
        "[%s %s]", TlTime::getNowDate().c_str(), TlTime::getNowTime().c_str());

    std::string title = "<<<< " + stepName + " ";
    TlUtils::pad(title, (72 - timeString.length()), ' ');
    title += timeString;

    // 出力
    this->logger(title + "\n");
}

void DfForce::calcForceFromNuclei() {
    this->logger("E grad in Nuclei");

    const int numOfAtoms = this->m_nNumOfAtoms;
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);

    TlDenseGeneralMatrix_Lapack F_nuc(this->m_nNumOfAtoms, 3);
    for (int i = 0; i < numOfAtoms; ++i) {
        const TlPosition posI = flGeom.getCoordinate(i);
        const double chargeI = flGeom.getCharge(i);

        for (int j = i + 1; j < numOfAtoms; ++j) {
            const TlPosition posJ = flGeom.getCoordinate(j);
            const double chargeJ = flGeom.getCharge(j);

            const TlPosition IJ = posJ - posI;
            const double distance = IJ.distanceFrom();
            const double invIJ3 = 1.0 / (distance * distance * distance);

            const double chargeIJ = chargeI * chargeJ;
            const double coef = -chargeIJ * invIJ3;
            const double fx = coef * IJ.x();
            const double fy = coef * IJ.y();
            const double fz = coef * IJ.z();

            F_nuc.add(i, X, fx);
            F_nuc.add(i, Y, fy);
            F_nuc.add(i, Z, fz);
            F_nuc.add(j, X, -fx);
            F_nuc.add(j, Y, -fy);
            F_nuc.add(j, Z, -fz);
        }
    }

    F_nuc *= -1.0;
    if (this->isDebugOutMatrix_ == true) {
        F_nuc.save("F_nuc.mtx");
    }

    this->force_ += F_nuc;
}

void DfForce::calcForceFromHpq(const TlDenseSymmetricMatrix_Lapack& P) {
    this->loggerTime("calc core-H");

    DfHpqX dfHpqX(&(this->pdfParamForForce_));
    TlDenseGeneralMatrix_Lapack force_Hpq(this->m_nNumOfAtoms, 3);
    TlDenseGeneralMatrix_Lapack force_Hpq_Xonly(this->m_nNumOfAtoms, 3);
    dfHpqX.getForce(P, &force_Hpq, &force_Hpq_Xonly);

    if (this->isDebugOutMatrix_ == true) {
        force_Hpq.save("F_h.mtx");
        force_Hpq_Xonly.save("F_h_Xonly.mtx");
    }

    this->force_ += force_Hpq;
    this->force_Xonly_ = force_Hpq_Xonly;
}

void DfForce::calcForceFromWS(RUN_TYPE runType) {
    this->logger("calc WS");

    DfOverlapX dfOvpX(&(this->pdfParamForForce_));

    const TlDenseSymmetricMatrix_Lapack W =
        this->getEnergyWeightedDensityMatrix(runType);
    if (this->isDebugOutMatrix_ == true) {
        W.save("W.mtx");
    }

    TlDenseGeneralMatrix_Lapack F_WS(this->m_nNumOfAtoms, 3);
    dfOvpX.getForce(W, &F_WS);

    F_WS *= -1.0;
    if (this->isDebugOutMatrix_ == true) {
        F_WS.save("F_WS.mtx");
    }

    if (runType == DfObject::RUN_RKS) {
        F_WS *= 2.0;
    }
    this->force_ += F_WS;
}

TlDenseGeneralMatrix_Lapack DfForce::getEnergyWeightedDensityMatrix(
    RUN_TYPE runType) {
    const int iteration = this->m_nIteration;
    // const int numOfAOs = this->m_nNumOfAOs;

    const TlDenseSymmetricMatrix_Lapack P =
        DfObject::getSpinDensityMatrix<TlDenseSymmetricMatrix_Lapack>(
            runType, iteration);
    const TlDenseSymmetricMatrix_Lapack F =
        DfObject::getFpqMatrix<TlDenseSymmetricMatrix_Lapack>(runType,
                                                              iteration);

    const TlDenseGeneralMatrix_Lapack W = P * F * P;

    return W;
}

void DfForce::calcForceFromCoulomb(RUN_TYPE runType) {
    if (this->J_engine_ == J_ENGINE_RI_J) {
        this->calcForceFromCoulomb_RIJ(runType);
    } else {
        this->calcForceFromCoulomb_exact(runType);
    }
}

void DfForce::calcForceFromCoulomb_exact(RUN_TYPE runType) {
    this->loggerTime("calc force from J");

    const int iteration = this->m_nIteration;
    const int numOfAtoms = this->m_nNumOfAtoms;

    const TlDenseSymmetricMatrix_Lapack P =
        this->getPInMatrix<TlDenseSymmetricMatrix_Lapack>(runType, iteration);

    DfEriX dfEri(&(this->pdfParamForForce_));

    // ((pq)'|(rs))
    TlDenseGeneralMatrix_Lapack F_J(numOfAtoms, 3);
    dfEri.getForceJ(P, &F_J);

    // F_J *= 0.5;
    if (this->isDebugOutMatrix_ == true) {
        F_J.save("F_J.mtx");
    }

    this->force_ += F_J;
}

void DfForce::calcForceFromCoulomb_RIJ(const RUN_TYPE runType) {
    this->loggerTime("calc force from RI_J");

    const int iteration = this->m_nIteration;
    const int numOfAtoms = this->m_nNumOfAtoms;

    const TlDenseVector_Lapack rho =
        this->getRho<TlDenseVector_Lapack>(runType, iteration);
    const TlDenseSymmetricMatrix_Lapack P =
        this->getPInMatrix<TlDenseSymmetricMatrix_Lapack>(runType, iteration);

    DfEriX dfEri(&(this->pdfParamForForce_));

    // ((pq)'|a)
    TlDenseGeneralMatrix_Lapack F_pqa(numOfAtoms, 3);
    dfEri.getForceJ(P, rho, &F_pqa);

    // (a'|b)
    TlDenseGeneralMatrix_Lapack F_ab(numOfAtoms, 3);
    dfEri.getForceJ(rho, &F_ab);

    if (this->isDebugOutMatrix_ == true) {
        F_pqa.save("F_pqa.mtx");
        F_ab.save("F_ab.mtx");
    }

    const TlDenseGeneralMatrix_Lapack F_J = (F_pqa - 0.5 * F_ab);
    this->force_ += F_J;
}

DfCalcGridX* DfForce::getCalcGridObj() {
    this->logger("create calc-grid engine");
    DfCalcGridX* pDfCalcGrid = new DfCalcGridX(&(this->pdfParamForForce_));

    return pDfCalcGrid;
}

void DfForce::calcForceFromPureXC(const RUN_TYPE runType) {
    this->loggerTime("calc XC(grid)");

    const int iteration = this->m_nIteration;
    const int numOfAtoms = this->m_nNumOfAtoms;

    TlDenseGeneralMatrix_Lapack Fxc(numOfAtoms, 3);

    // for RKS
    const TlDenseSymmetricMatrix_Lapack P = 0.5 * this->getPInMatrix<TlDenseSymmetricMatrix_Lapack>(runType, iteration);

    DfCalcGridX* pCalcGrid = this->getCalcGridObj();

    DfXCFunctional dfXCFunctional(&(this->pdfParamForForce_));
    DfXCFunctional::XC_TYPE xcType = dfXCFunctional.getXcType();
    switch (xcType) {
        case DfXCFunctional::SVWN: {
            DfFunctional_SVWN svwn;
            Fxc = pCalcGrid->energyGradient(P, &svwn);
        } break;

        case DfXCFunctional::BLYP: {
            DfFunctional_B88LYP blyp;
            Fxc = pCalcGrid->energyGradient(P, &blyp);
        } break;

        case DfXCFunctional::B3LYP: {
            DfFunctional_B3LYP b3lyp;
            Fxc = pCalcGrid->energyGradient(P, &b3lyp);
        } break;

        case DfXCFunctional::HF: {
            // do nothing
        } break;

        default:
            std::cerr << "unsupported functional. XC_TYPE=" << xcType
                      << std::endl;
            abort();
            break;
    }

    delete pCalcGrid;
    pCalcGrid = NULL;

    if (this->isDebugOutMatrix_ == true) {
        Fxc.save("F_xc.mtx");
    }
    this->force_ += Fxc;
}

void DfForce::calcForceFromPureXC_gridfree(RUN_TYPE runType) {
    this->loggerTime("calc XC(gridfree)");

    DfGridFreeXC dfGridFreeXC(&(this->pdfParamForForce_));
    TlDenseGeneralMatrix_Lapack force = dfGridFreeXC.getForce();

    const TlDenseGeneralMatrix_Lapack T = this->getTransformMatrix(force);
    force += T;

    if (this->isDebugOutMatrix_ == true) {
        force.save("F_xc_gf.mtx");
    }
    this->force_ += force;
}

void DfForce::calcForceFromK(RUN_TYPE runType) {
    this->loggerTime("calc K");

    const DfXCFunctional dfXCFunctional(&(this->pdfParamForForce_));

    if (dfXCFunctional.isHybridFunctional() == true) {
        const int iteration = this->m_nIteration;
        const int numOfAtoms = this->m_nNumOfAtoms;

        const TlDenseSymmetricMatrix_Lapack P =
            this->getPInMatrix<TlDenseSymmetricMatrix_Lapack>(runType, iteration);

        DfEriX dfEri(&(this->pdfParamForForce_));

        TlDenseGeneralMatrix_Lapack F_K(numOfAtoms, 3);
        // for RKS
        dfEri.getForceK(P, &F_K);
        if (runType == RUN_RKS) {
            F_K *= 0.5;
        }

        F_K *= -1.0;
        F_K *= dfXCFunctional.getFockExchangeCoefficient();  // for B3LYP

        if (this->isDebugOutMatrix_ == true) {
            F_K.save("F_K.mtx");
        }
        this->force_ += F_K;
    }
}

TlDenseGeneralMatrix_Lapack DfForce::getTransformMatrix(
    const TlDenseGeneralMatrix_Lapack& force) {
    const Fl_Geometry flGeom((*this->pPdfParam_)["coordinates"]);
    const int numOfAtoms = this->m_nNumOfAtoms;
    TlDenseGeneralMatrix_Lapack answer(numOfAtoms, 3);

    // 重心
    std::vector<TlPosition> X(numOfAtoms);
    {
        TlPosition center;
        for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
            center += flGeom.getCoordinate(atomIndex);
        }
        center /= numOfAtoms;
        this->log_.info(TlUtils::format("center: (% f, % f, % f)", center.x(),
                                        center.y(), center.z()));

        for (int i = 0; i < numOfAtoms; ++i) {
            X[i] = flGeom.getCoordinate(i) - center;
            this->log_.info(TlUtils::format("X[%d] (% e, % e, %e)", X[i].x(),
                                            X[i].y(), X[i].z()));
        }
    }

    TlDenseSymmetricMatrix_Lapack rot(3);
    for (int i = 0; i < numOfAtoms; ++i) {
        const TlPosition p = X[i];
        const double x = p.x();
        const double y = p.y();
        const double z = p.z();
        rot.add(0, 0, y * y + z * z);
        rot.add(0, 1, -x * y);
        rot.add(0, 2, -x * z);
        rot.add(1, 1, x * x + z * z);
        rot.add(1, 2, -y * z);
        rot.add(2, 2, x * x + y * y);
    }
    // rot.save("force_rot.mat");

    const double chk = rot.get(0, 0) * rot.get(1, 1) * rot.get(2, 2);
    this->log_.info(TlUtils::format("chk = % 8.3e", chk));

    //
    static const double TOO_SMALL = 1.0E-5;
    if (chk < TOO_SMALL) {
        if (rot.get(0, 0) > TOO_SMALL) {
            if (rot.get(1, 1) > TOO_SMALL) {
                // x, y != 0
                const double det = rot.get(0, 0) * rot.get(1, 1) -
                                   rot.get(0, 1) * rot.get(1, 0);
                const double inv_det = 1.0 / det;
                const double trp = rot.get(1, 1);
                rot.set(0, 0, rot.get(1, 1) * inv_det);
                rot.set(1, 1, trp * inv_det);
                rot.set(0, 1, -rot.get(0, 1) * inv_det);
            } else if (rot.get(2, 2) > TOO_SMALL) {
                // x, z != 0
                const double det = rot.get(0, 0) * rot.get(2, 2) -
                                   rot.get(0, 2) * rot.get(2, 0);
                const double inv_det = 1.0 / det;
                const double trp = rot.get(2, 2);
                rot.set(0, 0, rot.get(2, 2) * inv_det);
                rot.set(2, 2, trp * inv_det);
                rot.set(0, 2, -rot.get(0, 2) * inv_det);
            } else {
                // x != 0
                const double v = rot.get(0, 0);
                rot.set(0, 0, 1.0 / v);
            }
        } else if (rot.get(1, 1) > TOO_SMALL) {
            if (rot.get(2, 2) > TOO_SMALL) {
                // y, z != 0
                const double det = rot.get(2, 2) * rot.get(1, 1) -
                                   rot.get(2, 1) * rot.get(1, 2);
                const double inv_det = 1.0 / det;
                const double trp = rot.get(2, 2);
                rot.set(2, 2, rot.get(1, 1) * inv_det);
                rot.set(1, 1, trp * inv_det);
                rot.set(2, 1, -rot.get(2, 1) * inv_det);
            } else {
                // y != 0
                const double v = rot.get(1, 1);
                rot.set(1, 1, 1.0 / v);
            }
        } else if (rot.get(2, 2) > TOO_SMALL) {
            // z != 0
            const double v = rot.get(2, 2);
            rot.set(2, 2, 1.0 / v);
        } else {
            return answer;
        }
    } else {
        rot.inverse();
        // rot.save("force_rotinv.mat");
    }

    // <-- OK
    int TENS[3][3][3];
    TENS[0][0][0] = 0;
    TENS[1][0][0] = 0;
    TENS[2][0][0] = 0;
    TENS[0][1][0] = 0;
    TENS[1][1][0] = 0;
    TENS[2][1][0] = -1;
    TENS[0][2][0] = 0;
    TENS[1][2][0] = 1;
    TENS[2][2][0] = 0;

    TENS[0][0][1] = 0;
    TENS[1][0][1] = 0;
    TENS[2][0][1] = 1;
    TENS[0][1][1] = 0;
    TENS[1][1][1] = 0;
    TENS[2][1][1] = 0;
    TENS[0][2][1] = -1;
    TENS[1][2][1] = 0;
    TENS[2][2][1] = 0;

    TENS[0][0][2] = 0;
    TENS[1][0][2] = -1;
    TENS[2][0][2] = 0;
    TENS[0][1][2] = 1;
    TENS[1][1][2] = 0;
    TENS[2][1][2] = 0;
    TENS[0][2][2] = 0;
    TENS[1][2][2] = 0;
    TENS[2][2][2] = 0;

    const double invNumOfAtoms = 1.0 / double(numOfAtoms);
    for (int IP = 0; IP < numOfAtoms; ++IP) {
        const int KNDX = std::max(IP, 2 * IP - numOfAtoms);

        for (int JP = 0; JP <= IP; ++JP) {
            const int LNDX = std::max(JP, 2 * JP - numOfAtoms);

            for (int IC = 0; IC < 3; ++IC) {
                const int JEND = (JP == IP) ? (IC + 1) : 3;
                for (int JC = 0; JC < JEND; ++JC) {
                    double sum = 0.0;
                    for (int IA = 0; IA < 3; ++IA) {
                        for (int IB = 0; IB < 3; ++IB) {
                            if (TENS[IA][IB][IC] == 0) {
                                continue;
                            }
                            for (int JA = 0; JA < 3; ++JA) {
                                for (int JB = 0; JB < 3; ++JB) {
                                    if (TENS[JA][JB][JC] == 0) {
                                        continue;
                                    }
                                    sum += TENS[IA][IB][IC] * TENS[JA][JB][JC] *
                                           rot.get(IA, JA) * X[IP][IB] *
                                           X[JP][JB];
                                    // this->log_.info(TlUtils::format("rot(%d,%d)=%
                                    // e, X1(%d,%d)=% e, X2(%d,%d)=% e",
                                    //                                 IA,JA,rot.get(IA,
                                    //                                 JA),
                                    //                                 IP,IB,X[IP][IB],
                                    //                                 JP,JB,X[JP][JB]));
                                }
                            }
                        }
                    }
                    // this->log_.info(TlUtils::format("IP=%d, JP=%d, IC=%d,
                    // JC=%d, SUM=% e",
                    //                                 IP, JP, IC, JC, sum));

                    if (IC == JC) {
                        sum += invNumOfAtoms;
                    }

                    answer.add(LNDX, JC, -force.get(KNDX, IC) * sum);
                    if ((3 * LNDX + JC) != (3 * KNDX + IC)) {
                        answer.add(KNDX, IC, -force.get(LNDX, JC) * sum);
                    }
                }
            }
        }
    }

    return answer;
}
