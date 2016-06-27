#include <cassert>
#include <cstdlib>

#include "DfObject.h"
#include "DfHpqX.h"
#include "Fl_Geometry.h"
#include "TlSymmetricMatrix.h"
#include "TlGetopt.h"
#include "TlSerializeData.h"
#include "TlMsgPack.h"
#include "TlSymmetricMatrix.h"
#include "TlVector.h"
#include "TlUtils.h"

void help(const std::string& progName)
{
    std::cout << TlUtils::format("%s [options] args...", progName.c_str()) << std::endl;
    std::cout << "calc population of CEEM."
              << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << " -p PATH      set ProteinDF parameter file. default = pdfparam.mpac" << std::endl;
    std::cout << " -d PATH      set density matrix file. default is presumed by parameter file." << std::endl;
    std::cout << " -h           show help message (this)." << std::endl;
    std::cout << " -v           show message verbosely." << std::endl;
}


double getCharge(TlSerializeData param)
{
    double charge = 0.0;
    double chargeWithoutX = 0.0;

    Fl_Geometry flGeom(param["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();
    for (int i = 0; i < numOfAtoms; ++i) {
        const std::string symbol = flGeom.getAtomSymbol(i);
        const double currentCharge = flGeom.getCharge(i);

        charge += currentCharge;
        if (symbol != "X") {
            chargeWithoutX += currentCharge;
        }
    }

    return charge;
}


TlSymmetricMatrix makeA(TlSerializeData param)
{
    const Fl_Geometry flGeom(param["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();
    TlSymmetricMatrix A(numOfAtoms +1);
    for (int i = 0; i < numOfAtoms; ++i) {
        const TlPosition posA = flGeom.getCoordinate(i);
        for (int j = 0; j < i; ++j) {
            const TlPosition posB = flGeom.getCoordinate(j);
            const double d = posA.distanceFrom(posB);
            
            const double v = 1.0 / d;
            A.set(i, j, 0.5 * v);
            A.set(j, i, 0.5 * v);
        }

        // diagonal element
        A.set(i, i, 1.0);
    }

    
    // for lambda
    for (int i = 0; i < numOfAtoms; ++i) {
        A.set(numOfAtoms, i,  1.0);
        // A.set(i, numOfAtoms, -1.0); // not need because it's symmetric!
    }

    return A;
}


TlVector makeB(TlSerializeData param, const TlSymmetricMatrix& P)
{
    const Fl_Geometry flGeom(param["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();

    std::vector<TlPosition> grids(numOfAtoms);
    int atomCount = 0;
    for (int atomIndex = 0; atomIndex < numOfAtoms; ++atomIndex) {
        const std::string atomSymbol = flGeom.getAtomSymbol(atomIndex);
        if (atomSymbol == "X") {
            // not calculate in case of dummy charge
            continue;
        }

        const TlPosition p = flGeom.getCoordinate(atomIndex);
        std::cout << TlUtils::format("pos[%4d] (% 8.3f, % 8.3f % 8.3f)", atomIndex, p.x(), p.y(), p.z()) << std::endl;
        grids[atomIndex] = p;
        ++atomCount;
    }
    assert(atomCount == numOfAtoms);
    
    DfHpqX dfHpq(&param);
    std::vector<double> values = dfHpq.getESP(P, grids);

    TlVector B(numOfAtoms +1);
    for (int i = 0; i < numOfAtoms; ++i) {
        B.set(i, - values[i]);
    }

    // charge
    // const double q = getCharge(param);
    const double q = 0.0;
    B.set(numOfAtoms, - q);

    return B;
}


double estimateEnergy(TlSerializeData param, const TlVector& x)
{
    const Fl_Geometry flGeom(param["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();

    double e = 0.0;
    for (int i = 0; i < numOfAtoms; ++i) {
        const TlPosition posA = flGeom.getCoordinate(i);
        for (int j = 0; j < i; ++j) {
            const TlPosition posB = flGeom.getCoordinate(j);
            const double d = posA.distanceFrom(posB);
            e += x[j] * x[i] / d;
        }
    }

    return e;
}


TlVector getQ(TlSerializeData param, const TlVector& q)
{
    Fl_Geometry flGeom(param["coordinates"]);
    const int numOfAtoms = flGeom.getNumOfAtoms();

    TlVector Q(numOfAtoms);
    for (int i = 0; i < numOfAtoms; ++i) {
        const std::string symbol = flGeom.getAtomSymbol(i);
        const double currentCharge = flGeom.getCharge(i);
        
        // Q[i] = q[i] + currentCharge;
        Q[i] = q[i];
    }

    return Q;
}


int main(int argc, char* argv[])
{
    TlGetopt opt(argc, argv, "p:d:i:fcmo:hv");
    
    const bool verbose = (opt["v"] == "defined");
    if ((opt["h"] == "defined")) {
        help(opt[0]);
        return EXIT_SUCCESS;
    }

    std::string pdfParamPath = "pdfparam.mpac";
    if (opt["p"].empty() != true) {
        pdfParamPath = opt["p"];
    }
    std::string PMatrixFilePath = "";
    if (opt["d"].empty() != true) {
        PMatrixFilePath = opt["d"];
    }

    // パラメータファイルの読み込み
    TlSerializeData param;
    {
        TlMsgPack mpac;
        mpac.load(pdfParamPath);
        param = mpac.getSerializeData();
    }

    // 密度行列の読み込み
    TlSymmetricMatrix P;
    if (PMatrixFilePath == "") {
        const int iteration = param["model"]["iterations"].getInt();
        DfObject::RUN_TYPE runType = DfObject::RUN_RKS;
        DfObject dfObject(&param);
        PMatrixFilePath = dfObject.getPpqMatrixPath(runType, iteration);
    }
    P.load(PMatrixFilePath);

    // A
    TlSymmetricMatrix A = makeA(param);
    A.save("ceem_a.mat");

    // B; ESP計算
    TlVector B = makeB(param, P);
    B.save("ceem_b.vtr");
    
    // solve
    TlMatrix Ainv = A;
    Ainv.inverse();
    TlVector x = Ainv * B;
    x.save("ceem_x.vtr");

    TlVector Q = getQ(param, x);
    Q.print(std::cout);

    const double e = estimateEnergy(param, x);
    std::cout << TlUtils::format("estimate Ene: %16.10f", e) << std::endl;

    return EXIT_SUCCESS;
}


