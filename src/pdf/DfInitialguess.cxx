#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <cstring>

#include "DfInitialguess.h"

#include "Fl_Geometry.h"
#include "Fl_Gto_Density.h"
#include "Fl_Gto_Xcpot.h"
#include "Fl_Gto_Xcpot2.h"
#include "Fl_Gdb_Molecular.h"
#include "Fl_Gdb_Atom.h"
#include "Fl_Tbl_Density.h"
#include "Fl_Tbl_Xcpot.h"
#include "Fl_Tbl_Xcpot2.h"

#include "DfDensityFittingX.h"

#include "DfInvMatrix.h"

#include "TlSymmetricMatrix.h"
#include "TlPrdctbl.h"
#include "TlLogX.h"
#include "CnError.h"

int DfInitialguess::w1RouCount = 0;
int DfInitialguess::w1MyuCount = 0;
int DfInitialguess::w1NyuCount = 0;

DfInitialguess::DfInitialguess(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
    this->prepare();

    this->getMemory();
    this->setGuessDBfile();
    this->prepare2();
    this->setKeyword();
}


DfInitialguess::~DfInitialguess()
{
}

int DfInitialguess::dfGusMain()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    
    TlLogX& log = TlLogX::getInstance();
    this->tempflag = 0;

    const std::string sScfStartGuess = pdfParam["scf-start-guess"].getStr();
    if ((sScfStartGuess == "rho") || (sScfStartGuess == "atom_rho")) {
        this->readData();        // Read and Analyzy needData.
        this->calcMain();        // Mainfunction for Calculation.
    } else if (sScfStartGuess == "file_rho") {
        log << " ##### DfInitialguess, Treatment of user_rho is begun.####\n";
        this->setCoefRoudata();
    }

    this->VctNormfromIGuess2(); // Normalyzation rootine copied from DfIGuess2

    this->GusOutput();       // Output : InitialValue(Rou,Myu,Nyu)

    if (this->tempflag == 0) {
        if (pdfParam["myu-nyu-zero"].getStr() == "yes") {
            log << "\n@@@@ DfInitialguess FACE TO 1/1000(ZERO) vector for Myu, Nyu by EKA @@@@\n\n";

            const int niter = 1;
            // read myu vector
            TlVector myu;
            myu.load("fl_Work/fl_Vct_Myu" + TlUtils::xtos(niter));
            myu /= 1000.0;
            myu.save("fl_Work/fl_Vct_Myu" + TlUtils::xtos(niter));

            // read nyu vector
            TlVector nyu;
            nyu.load("fl_Work/fl_Vct_Nyu" + TlUtils::xtos(niter));
            nyu /= 1000.0;
            nyu.save("fl_Work/fl_Vct_Myu" + TlUtils::xtos(niter));
        }
    }

    return 0;
}


int DfInitialguess::readData()
{
    // DataSet : Coefficient of EXP's shoulder and BasisSetName.
    // DataSet : struct Auxtype's Number for each Atom.
    this->setStruct();

    // Make Table : Table is datafile for making Initialguess.
    this->makeStepTable();

    // This routine make TempGDBfile fromGuessDB.
    // if need, execute corresponding orbital method.
    this->setGDBdatafile();

    // DataSet : Input Atomdata's VectorCoef in struct Auxtype .
    this->setAtomCoef();

    return 0;
}

int DfInitialguess::prepare()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    
    Fl_Geometry FlGeom(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Density FlGtoDen;
    Fl_Gto_Xcpot FlGtoXcpot;

    //--get scftype data --------------------------------------------------|
    if ((pdfParam["method"].getStr() == "roks") ||
        (pdfParam["method"].getStr() == "sp")) {
        this->scftype = SP ;
        if (pdfParam["method"].getStr() == "roks") {
            // roks
            this->AlphaSpinNum = pdfParam["method/roks/electron-number"].getInt();
            this->BetaSpinNum  = pdfParam["method/roks/electron-number-alpha"].getInt();
            this->ElectronNum  = pdfParam["method/roks/electron-number-beta"].getInt();
        } else {
            // sp
            this->AlphaSpinNum = pdfParam["UKS/alphaElectrons"].getInt();
            this->BetaSpinNum  = pdfParam["UKS/betaElectrons"].getInt();
            this->ElectronNum  = AlphaSpinNum + BetaSpinNum ;
        }

        if (pdfParam["xc-potential"].getStr() == "xalpha") {
            this->nAlpha = pdfParam["xc-potential/xalpha/alpha-value"].getDouble();
        } else {
            this->nAlpha = pdfParam["xc-potential/gxalpha/alpha-value"].getDouble();
        }
    } else {
        // rks
        this->scftype = NSP ;
        this->ElectronNum  = pdfParam["RKS/electrons"].getInt();
        if (pdfParam["xc-potential"].getStr() == "xalpha") {
            this->nAlpha = pdfParam["xc-potential/xalpha/alpha-value"].getDouble();
        } else {
            this->nAlpha = pdfParam["xc-potential/gxalpha/alpha-value"].getDouble();
        }
    }

    this->AtomNum              = FlGeom.getNumOfAtoms();      // set Total Atom Number;
    this->AtomKindNumInclDummy = FlGeom.getAtomKindNumber();  // set Atom Kind Num including Dummy;
    this->DumyAtomNum          = FlGeom.getDummyatom();       // ser
    this->AtomKindNum          = this->AtomKindNumInclDummy;
    if (DumyAtomNum != 0) {
        AtomKindNum -=  1;
    }
    this->AtomNum -= this->DumyAtomNum;

    this->MaxTermRou = 0;
    for (int j=0; j<AtomNum; j++) {
        const std::string Atm = FlGeom.getAtom(j);
        const std::string Lb2 = FlGeom.getLabel(j);

        for (int i=0; i<AtomKindNum; i++) {
            const int cgtoPoint = FlGtoDen.getStartposition(i);
            if ((Atm == FlGtoDen.getAtom(cgtoPoint)) &&
                    (Lb2 == FlGtoDen.getLabel(cgtoPoint))) {
                for (int k = cgtoPoint; k<cgtoPoint+FlGtoDen.getTermnumber(i); k++) {
                    int cGto = 0;
                    switch (FlGtoDen.getShell(k)) {
                    case 's':
                        cGto=1;
                        break;
                    case 'p':
                        cGto=3;
                        break;
                    case 'd':
                        cGto=5;
                        break;
                    default:
                        this->GusError("Bad Shell Char !!");
                        break;
                    }
                    this->MaxTermRou += (FlGtoDen.getContraction(k)*cGto) ;
                }
                break;
            }
        }
    }

    this->MaxTermMyu = 0;
    for (int j=0; j<AtomNum; j++) {
        const std::string Atm = FlGeom.getAtom(j);
        const std::string Lb2 = FlGeom.getLabel(j);
        for (int i=0; i<AtomKindNum; i++) {
            const int cgtoPoint = FlGtoXcpot.getStartposition(i);
            if ((Atm == FlGtoXcpot.getAtom(cgtoPoint)) &&
                    (Lb2 == FlGtoXcpot.getLabel(cgtoPoint))) {
                for (int k = cgtoPoint; k < cgtoPoint+FlGtoXcpot.getTermnumber(i); k++) {
                    int cGto = 0;
                    switch (FlGtoXcpot.getShell(k)) {
                    case 's':
                        cGto=1;
                        break;
                    case 'p':
                        cGto=3;
                        break;
                    case 'd':
                        cGto=5;
                        break;
                    default:
                        GusError("Bad Shell Char !!");
                        break;
                    }
                    this->MaxTermMyu += (FlGtoXcpot.getContraction(k)*cGto);
                }

                break;
            }
        }
    }

    this->MaxTermNyu = this->MaxTermMyu;

    return 0;
}

int DfInitialguess::getMemory()
{
    //EachAD  = new eachAtomData[ AtomNum ];
    this->EachAD.resize(this->AtomNum);

    //Atomtype = new Auxtype[ AtomKindNum ];
    this->Atomtype.resize(this->AtomKindNum);

    //ρについての構造体中のベクトルのメモリを取るルーチン
    {
        Fl_Gto_Density FlGtoDen;
        int to = 0;
        int count = 0;
        for (int i = 0; i<AtomKindNum; i++) {
            int from = to;
            to  += FlGtoDen.getTermnumber(i);
            int term = 0;
            for (int j=from; j<to; j++) {
                int cGto = 0;
                switch (FlGtoDen.getShell(j)) {
                case 's':
                    cGto=1;
                    break;
                case 'p':
                    cGto=3;
                    break;
                case 'd':
                    cGto=5;
                    break;
                default:
                    GusError("Bad Shell Char !!");
                    break;
                }
                term += (FlGtoDen.getContraction(j) * cGto);
            }

            this->Atomtype[i].alpha.resize(term);
            this->Atomtype[i].rou.resize(term);

            term=0;
            count = to;
        }
    }

    //μについての構造体中のベクトルのメモリを取るルーチン
    {
        Fl_Gto_Xcpot FlGtoXcpot;
        int to=0;
        int count=0;
        for (int i=0; i<AtomKindNum; i++) {
            int from = to;
            to   = to + FlGtoXcpot.getTermnumber(i);
            int term=0;
            for (int j=from; j<to; j++) {
                int cGto = 0;
                switch (FlGtoXcpot.getShell(j)) {
                case 's':
                    cGto=1;
                    break;
                case 'p':
                    cGto=3;
                    break;
                case 'd':
                    cGto=5;
                    break;
                default:
                    GusError("Bad Shell Char !!");
                    break;
                }
                term += (FlGtoXcpot.getContraction(j)*cGto);
            }

            this->Atomtype[i].gamma.resize(term);
            this->Atomtype[i].myu.resize(term);

            term=0;
            count = to;
        }
    }

    //νについての構造体中のベクトルのメモリを取るルーチン
    {
        Fl_Gto_Xcpot2   FlGtoXcpot2;

        int to=0;
        int count=0;
        for (int i=0; i<AtomKindNum; i++) {
            int from = to;
            to   = to + FlGtoXcpot2.getTermnumber(i);
            int term=0;
            for (int j=from; j<to; j++) {
                int cGto = 0;
                switch (FlGtoXcpot2.getShell(j)) {
                case 's':
                    cGto=1;
                    break;
                case 'p':
                    cGto=3;
                    break;
                case 'd':
                    cGto=5;
                    break;
                default:
                    GusError("Bad Shell Char !!");
                    break;
                }
                term += (FlGtoXcpot2.getContraction(j)*cGto);
            }

//       Atomtype[i].gamma2 = new double[ term ];
//       Atomtype[i].nyu = new double[ term ];
            this->Atomtype[i].gamma2.resize(term);
            this->Atomtype[i].nyu.resize(term);

            term=0;
            count = to;
        }
    }

    // for Work Variable;
//   DbFile  = new char[MaxFname];     assert(DbFile  != 0);
//   DbFile2 = new char[MaxFname];     assert(DbFile2  != 0);
//   WorkDB  = new char[MaxFname];     assert(WorkDB  != 0);
//   serchFile = new char[MaxFname];   assert(serchFile != 0);

    // For  Database File
    this->rouNormOnOff.resize(AtomNum);
    this->myuNormOnOff.resize(AtomNum);
    this->nyuNormOnOff.resize(AtomNum);

    this->GuessStep.resize(MaxStepNum);

    // First : this get_memory is in EXE_MOLECULE
    DbAmino.Order.resize(Max1DbAtomNum);
    DbAmino.equalPaircount.resize(MaxEqualPairCount);
    DbAmino.equalPair.resize(MaxEqualPair);
    DbAmino.Atomdata.resize(Max1DbAtomNum);
    InpAmino.Order.resize(Max1DbAtomNum);
    InpAmino.Atomdata.resize(Max1DbAtomNum);

    for (int j=0; j<Max1DbAtomNum; j++) {
        DbAmino.Atomdata[j].afinmat.resize(4*4);    // afinMatrix = [4 x 4];
        DbAmino.Atomdata[j].CoefRou.resize(aboutAuxTerm);
        DbAmino.Atomdata[j].CoefMyu.resize(aboutAuxTerm);
        DbAmino.Atomdata[j].CoefNyu.resize(aboutAuxTerm);
        DbAmino.Atomdata[j].rouDbShellNum.resize(MaxSPDterm);
        DbAmino.Atomdata[j].myuDbShellNum.resize(MaxSPDterm);
        DbAmino.Atomdata[j].nyuDbShellNum.resize(MaxSPDterm);
        DbAmino.Atomdata[j].Conection.resize(MaxCntNum);
    }


    return 0;
}

void DfInitialguess::setGuessDBfile()
{
    char* PdfHome = getenv("PDF_HOME");
    std::string GdbDir;
    if (PdfHome) {
        GdbDir = TlUtils::format("%s/data", PdfHome);
    } else {
        GdbDir = ".";
    }

    this->GDBfile0 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Atom", GdbDir.c_str());
    this->GDBfile1 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Namino.header", GdbDir.c_str());
    this->GDBfile2 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Peptid.header", GdbDir.c_str());
    this->GDBfile3 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Camino.header", GdbDir.c_str());
    this->GDBfile4 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Molecular.header", GdbDir.c_str());
    this->GDBfile5 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_UserMolec.header", GdbDir.c_str());

    this->GDBfile1_2 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Namino", GdbDir.c_str());
    this->GDBfile2_2 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Peptid", GdbDir.c_str());
    this->GDBfile3_2 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Camino", GdbDir.c_str());
    this->GDBfile4_2 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_Molecular", GdbDir.c_str());
    this->GDBfile5_2 = TlUtils::format("%s/fl_GuessDB/fl_Gdb_UserMolec", GdbDir.c_str());
}

int DfInitialguess::prepare2()
{
    Fl_Geometry FlGeom(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Density FlGtoDen;
    Fl_Gto_Xcpot FlGtoXcpot;

    // Prepare dataset for Rou.
    {
        int Now_AuxLine=0;
        std::vector<int> AtomSearchVector(this->AtomNum);
        for (int j=0; j < this->AtomNum; j++) {
            const std::string Atm = FlGeom.getAtom(j);
            const std::string Lb2 = FlGeom.getLabel(j);

            for (int i=0; i<AtomKindNum; i++) {
                int cgtoPoint = FlGtoDen.getStartposition(i);
                if ((Atm == FlGtoDen.getAtom(cgtoPoint)) && (Lb2 == FlGtoDen.getLabel(cgtoPoint))) {
                    int S_Shell = FlGtoDen.getSnum(cgtoPoint);
                    int P_Shell = FlGtoDen.getPnum(cgtoPoint);
                    int D_Shell = FlGtoDen.getDnum(cgtoPoint);

                    AtomSearchVector[j] = Now_AuxLine ;
                    Now_AuxLine += (S_Shell + (P_Shell*3) + (D_Shell*5));
                    break;
                }
            }
        }

        // Copy : AtomSerchVector ----> struct eachAtomData's data_member.
        for (int j=0; j < this->AtomNum; j++) {
            this->EachAD[j].Rou_StartNum = AtomSearchVector[j];
        }
    }

    //    Prepare dataset for Myu.
    {
        int Now_AuxLine = 0;
        std::vector<int> AtomSearchVector(this->AtomNum);
        for (int j=0; j<AtomNum; j++) {
            const std::string Atm = FlGeom.getAtom(j);
            const std::string Lb2 = FlGeom.getLabel(j);
            for (int i=0; i<AtomKindNum; i++) {
                int cgtoPoint = FlGtoXcpot.getStartposition(i);
                if ((Atm == FlGtoXcpot.getAtom(cgtoPoint)) && (Lb2 == FlGtoXcpot.getLabel(cgtoPoint))) {
                    int S_Shell = FlGtoXcpot.getSnum(cgtoPoint);
                    int P_Shell = FlGtoXcpot.getPnum(cgtoPoint);
                    int D_Shell = FlGtoXcpot.getDnum(cgtoPoint);

                    AtomSearchVector[j] = Now_AuxLine ;
                    Now_AuxLine += (S_Shell + (P_Shell*3) + (D_Shell*5)) ;
                }
            }
        }

        // Copy : AtomSerchVector ----> struct eachAtomData's data_member.
        for (int j=0; j<AtomNum; j++) {
            this->EachAD[j].Myu_StartNum = AtomSearchVector[j];
            this->EachAD[j].Nyu_StartNum = AtomSearchVector[j];
            // Myu == Nyu for Auxiliary Number. But eachValue is not equal.
        }
    }

    return 0;
}

int DfInitialguess::setKeyword()
{
    const TlSerializeData& pdfParam = *(this->pPdfParam_);
    
    TlLogX& log = TlLogX::getInstance();
    int   FromAtom[MaxPartialNum];
    int   ToAtom[MaxPartialNum];
    double Ele[MaxPartialNum];

    std::vector<std::string> VctFname(MaxPartialNum, "");

    DfInitialguess::WhichRMN  Conbination[MaxPartialNum];

    // Matrix File Name in NSP.
    {
        const std::string tmp = pdfParam["guess/nsp-ppq"].getStr();
        this->PpqFname = tmp;
        if (tmp == "nil") {
            this->PpqFname = "";
        } else {
            if (this->scftype != NSP) {
                log<<"Bad appointment Matrix Input Filename\n";
                log<<"Calulation Method is SP , now .\n";
                log<<"Input 2 Filename in guess/sp-ppq\n";
                CnErr.abort();
            }
        }
    }

    //Matrix File Name in SP
    {
        std::string tmp = pdfParam["guess/sp-ppq"].getStr();
        if (tmp == "nil") {
            this->AlphaPpqFname = "";
            this->BetaPpqFname  = "";
        } else {
            if (this->scftype != SP) {
                log << "Bad appointment Matrix Input Filename\n";
                log << "Calulation Method is NSP , now .\n";
                log << "Input Filename in guess/nsp-ppq\n";
                CnErr.abort();
            }

            std::vector<std::string> sItems = TlUtils::split(tmp, " ");
            if (sItems.size() != 2) {
                log << "Bad appointment : MatrixInput for SP \n";
                log << "Input :  two filenames for alpha and beta\n";
                CnErr.abort();
            }
            this->AlphaPpqFname = sItems[0];
            this->BetaPpqFname  = sItems[1];
        }
    }

    //  Option : Threshold Angle Gosa for Trans Density
    {
        std::string dummy = pdfParam["guess/trans-angle-threshold"].getStr();
        std::vector<std::string> sItems = TlUtils::split(dummy, " ");

        this->FirstAveGosa = std::atof(sItems[0].c_str());
        if (this->FirstAveGosa < 1E-6) {
            log << " FirstAveGosa = ["<<FirstAveGosa<<"] is very small\n";
            log << " This suitable value is from 0.1 to 2.0 .\n";
            CnErr.abort();
        }

        this->FirstEachGosa = std::atof(sItems[1].c_str());
        if (this->FirstEachGosa < this->FirstAveGosa) {
            log << "You must appoint (FirstEachGosa > FirstAveGosa)\n";
            log << "But now , FirstEachGosa = " << this->FirstEachGosa << "\n";
            log << "          FirstAveGosa  = " << this->FirstAveGosa  << "\n";
            CnErr.abort();
        }

        this->LastAveGosa = std::atof(sItems[2].c_str());
        if (this->LastAveGosa < this->FirstAveGosa) {
            log <<"You must appoint (LastAveGosa > FirstAveGosa)\n";
            log <<"But now , FirstAveGosa = " << this->FirstAveGosa << "\n";
            log <<"          LastAveGosa  = " << this->LastAveGosa  << "\n";
            CnErr.abort();
        }

        this->LastEachGosa = std::atof(sItems[3].c_str());
        if (this->LastAveGosa > this->LastEachGosa) {
            log << "You must appoint (LastAveGosa > LastEachGosa)\n";
            log << "But now , LastAveGosa  = "  << this->LastAveGosa  << "\n";
            log << "          LastEachGosa  = " << this->LastEachGosa << "\n";
            CnErr.abort();
        }
    }

    // Method of making Myu and Nyu.
    {
        std::string dummy = pdfParam["guess/make-myu-nyu"].getStr();
        if (dummy == "meth0") {
            this->MethodMyuNyu = Meth0 ;
        } else if (dummy == "meth1") {
            this->MethodMyuNyu = Meth1 ;
        } else if (dummy == "meth2") {
            this->MethodMyuNyu = Meth2 ;
        } else if (dummy == "meth3") {
            this->MethodMyuNyu = Meth3 ;
        } else if (dummy == "meth4") {
            this->MethodMyuNyu = Meth4 ;
        } else {
            log << "Bad appointment guess/make-myu-nyu\n";
            CnErr.abort();
        }
    }

    // Appointment : Total Normalize for Vector(Rou,Myu,Nyu)
    {
        std::string dummy = pdfParam["guess/vct-normalize"].getStr();
        if (dummy.empty()) {
            log<<"Bad appointment vct-normalize "<<"\n";
            log<<"Input ON or OFF for Rou,Myu,Nyu."<<"\n";
            CnErr.abort();
        } else {
            std::vector<std::string> sItems = TlUtils::split(dummy, " ");

            if (TlUtils::toUpper(sItems[0]) == "ON") {
                this->VctRouNormalize = ON;
            } else if (TlUtils::toUpper(sItems[0]) == "OFF") {
                this->VctRouNormalize = OFF;
            } else {
                log<<"Bad appointment guess/vct-normalize (Rou)"<<"\n";
                CnErr.abort();
            }

            if (TlUtils::toUpper(sItems[1]) == "ON") {
                this->VctMyuNormalize = ON;
            } else if (TlUtils::toUpper(sItems[1]) == "OFF") {
                this->VctMyuNormalize = OFF;
            } else {
                log << "Bad appointment guess/vct-normalize (Myu)\n";
                CnErr.abort();
            }

            if (TlUtils::toUpper(sItems[2]) == "ON") {
                this->VctNyuNormalize = ON;
            } else if (TlUtils::toUpper(sItems[2]) == "OFF") {
                this->VctNyuNormalize = OFF;
            } else {
                log << "Bad appointment guess/vct-normalize (Nyu)\n";
                CnErr.abort();
            }

        }
    }

    // Appointment : Partial Normalize for Vector(Rou,Myu,Nyu)
    {
        std::string tmp = pdfParam["guess/part-normalize"].getStr();
        //std::cerr << "[TH] DfInitialguess::setKeyword() guess/part-normalize = " << tmp << "." << std::endl;
        if (tmp == "nil") {
            this->PartialRouNormNum = 0;
            this->PartialMyuNormNum = 0;
            this->PartialNyuNormNum = 0;
        } else {
            std::vector<std::string> sItems = TlUtils::split(tmp, " ");
            int blocknumber = std::atoi(sItems[0].c_str());
            if (blocknumber == 0) {
                log << "Cannot Transfer atoi(dumyCHAR2).\n";
                log << "Check : format for guess/part-normalize\n";
                CnErr.abort();
            }

            int rcount = 0;
            int mcount = 0;
            int ncount = 0;
            for (int i=0; i < blocknumber; i++) {
                std::string tmp2 = TlUtils::toUpper(sItems[i +1]);
                if (tmp2 == "ROU") {
                    rcount++;
                    Conbination[i] = Rou;
                } else if (tmp2 == "MYU") {
                    mcount++;
                    Conbination[i] = Myu;
                } else if (tmp2 == "NYU") {
                    ncount++;
                    Conbination[i] = Nyu ;
                } else {
                    log << "Bad appintment guess/part-normalize\n";
                    CnErr.abort();
                }

                FromAtom[i] = std::atoi(sItems[i +2].c_str());
                ToAtom[i]   = std::atoi(sItems[i +3].c_str());
                Ele[i]      = std::atof(sItems[i +4].c_str());
            }

            // ここから実際のデータメンバ（キーワード）に代入する。
            this->PartialRouNormNum = rcount ;
            this->PartialMyuNormNum = mcount ;
            this->PartialNyuNormNum = ncount ;

            if (this->PartialRouNormNum != 0) {
                this->PartialRouNormalize.resize(PartialRouNormNum);
            }
            if (this->PartialMyuNormNum != 0) {
                this->PartialMyuNormalize.resize(PartialMyuNormNum);
            }
            if (PartialNyuNormNum != 0) {
                this->PartialNyuNormalize.resize(PartialNyuNormNum);
            }

            rcount = mcount = ncount = 0;
            for (int i=0; i<blocknumber; i++) {
                if (Conbination[i] == Rou) {
                    this->PartialRouNormalize[rcount].From1   = FromAtom[i];
                    this->PartialRouNormalize[rcount].To1     = ToAtom[i];
                    this->PartialRouNormalize[rcount].ElecNum = Ele[i];
                    rcount++;
                } else if (Conbination[i] == Myu) {
                    this->PartialMyuNormalize[mcount].From1   = FromAtom[i];
                    this->PartialMyuNormalize[mcount].To1     = ToAtom[i];
                    this->PartialMyuNormalize[mcount].ElecNum = Ele[i];
                    mcount++;
                } else if (Conbination[i] == Nyu) {
                    this->PartialNyuNormalize[ncount].From1   = FromAtom[i];
                    this->PartialNyuNormalize[ncount].To1     = ToAtom[i];
                    this->PartialNyuNormalize[ncount].ElecNum = Ele[i];
                    ncount++;
                }
            }
        }
    }

    // Input Normalize flag(rouNormOnOff,..) for Rou,Myu,Nyu .
    for (int i=0; i<AtomNum; i++) {
        // Initialize;  ON:defined , OFF:non defined.
        this->rouNormOnOff[i] = OFF;
        this->myuNormOnOff[i] = OFF;
        this->nyuNormOnOff[i] = OFF;
    }

    for (int i=0; i < PartialRouNormNum; i++) {
        int startAtom = PartialRouNormalize[i].From1;
        int finalAtom = PartialRouNormalize[i].To1;
        for (int j = startAtom; j <= finalAtom; j++) {
            this->rouNormOnOff[j] = ON;
        }
    }

    for (int i=0; i < PartialMyuNormNum; i++) {
        int startAtom = PartialMyuNormalize[i].From1;
        int finalAtom = PartialMyuNormalize[i].To1;
        for (int j = startAtom; j <= finalAtom; j++) {
            this->myuNormOnOff[j] = ON;
        }
    }

    for (int i=0; i < PartialNyuNormNum; i++) {
        int startAtom = PartialNyuNormalize[i].From1;
        int finalAtom = PartialNyuNormalize[i].To1;
        for (int j=startAtom; j <= finalAtom; j++) {
            this->nyuNormOnOff[j] = ON;
        }
    }

    // Appointment : User difine Vector(Rou,Myu,Nyu)
    {
        std::string dummy = pdfParam["guess/user-vector"].getStr();
        //std::cerr << "[TH] DfInitialguess::setKeyword() guess/user-vector = " << dummy << "." << std::endl;
        if (dummy == "nil") {
            this->UserVctRouNum = 0;
            this->UserVctMyuNum = 0;
            this->UserVctNyuNum = 0;
        } else {
            std::vector<std::string> sItems = TlUtils::split(dummy, " ");

            int blocknumber = std::atoi(sItems[0].c_str()) ;
            if (blocknumber == 0) {
                log << "Cannot Transfer atoi(dumyCHAR2).\n";
                log << "Check : format for guess/user-vector\n";
                CnErr.abort();
            }

            int rcount = 0;
            int mcount = 0;
            int ncount = 0;

            for (int i=0; i<blocknumber; i++) {
                std::string tmp2 = TlUtils::toUpper(sItems[i +1]);
                if (tmp2 == "ROU") {
                    rcount++;
                    Conbination[i] = Rou;
                } else if (tmp2 == "MYU") {
                    mcount++;
                    Conbination[i] = Myu ;
                } else if (tmp2 == "NYU") {
                    ncount++;
                    Conbination[i] = Nyu ;
                } else {
                    log << "Bad appintment guess/part-normalize\n";
                    CnErr.abort();
                }

                FromAtom[i] = std::atoi(sItems[i +2].c_str());
                ToAtom[i]   = std::atoi(sItems[i +3].c_str());
                VctFname[i] = sItems[i +4];
            }

            // ここから実際のデータメンバ（キーワード）に代入する。
            this->UserVctRouNum = rcount ;
            this->UserVctMyuNum = mcount ;
            this->UserVctNyuNum = ncount ;

            if (this->UserVctRouNum != 0) {
                this->UserVctRou.resize(UserVctRouNum);
//  for (int i=0; i<UserVctRouNum; i++){
//    UserVctRou[i].Filename = new char[MaxFname];
//  }
            }
            if (this->UserVctMyuNum != 0) {
                this->UserVctMyu.resize(UserVctMyuNum);
//  for (int i=0; i<UserVctMyuNum; i++){
//    UserVctMyu[i].Filename = new char[MaxFname];
//  }
            }
            if (this->UserVctNyuNum != 0) {
                this->UserVctNyu.resize(UserVctNyuNum);
//  for (int i=0; i<UserVctNyuNum; i++){
//    UserVctNyu[i].Filename = new char[MaxFname];
//  }
            }

            rcount = mcount = ncount = 0;
            for (int i=0; i<blocknumber; i++) {
                if (Conbination[i] == Rou) {
                    UserVctRou[rcount].From1    = FromAtom[i];
                    UserVctRou[rcount].To1      = ToAtom[i];
                    UserVctRou[rcount].Filename = VctFname[i];
                    rcount++;
                } else if (Conbination[i] == Myu) {
                    UserVctMyu[mcount].From1    = FromAtom[i];
                    UserVctMyu[mcount].To1      = ToAtom[i];
                    UserVctMyu[mcount].Filename = VctFname[i];
                    mcount++;
                } else if (Conbination[i] == Nyu) {
                    UserVctNyu[ncount].From1    = FromAtom[i];
                    UserVctNyu[ncount].To1      = ToAtom[i];
                    UserVctNyu[ncount].Filename = VctFname[i];
                    ncount++;
                }
            }

        }
    }

    //ここで、各種指定ファイルが存在し、その中身の値がおかしくないかを
    //この時点で調べておく。
    // for Rou
    for (int i=0; i < this->UserVctRouNum; i++) {
        int startAtom = UserVctRou[i].From1 ;
        int finalAtom = UserVctRou[i].To1;
        int endaux;
        if (finalAtom != AtomNum-1) {
            endaux   = this->EachAD[finalAtom+1].Rou_StartNum ;
        } else {
            endaux = MaxTermRou ;
        }
        int startaux = this->EachAD[startAtom].Rou_StartNum;
        int auxlength = endaux - startaux ;
        {
            //Fl_GuessVector FGV(UserVctRou[i].Filename);
            TlVector FGV;
            FGV.load(UserVctRou[i].Filename);

            std::cout << "----following debug out by eka----" << std::endl;
            std::cout << "              DfInititalguess.cxx modified by eka" << std::endl;
            std::cout << AtomNum << std::endl;
            std::cout << startAtom << " "<< finalAtom << std::endl;
            std::cout << startaux << " " << endaux << " " << auxlength << std::endl;
            //cout << FGV.getAuxterm() << endl;
            std::cout << FGV.getSize() << std::endl;
            std::cout << "----------------------------------" << std::endl;
            //if( auxlength != FGV.getAuxterm() ){
            if (static_cast<TlVector::size_type>(auxlength) != FGV.getSize()) {
                log<<"Bad defined VectorFile [" <<UserVctRou[i].Filename<<"]\n";
                CnErr.abort();
            }
        }
    }

    // for Myu
    for (int i=0; i < this->UserVctMyuNum; i++) {
        int startAtom = UserVctMyu[i].From1 ;
        int finalAtom = UserVctMyu[i].To1;
        int endaux;
        if (finalAtom != AtomNum-1) {
            endaux   = this->EachAD[finalAtom+1].Myu_StartNum ;
        } else {
            endaux = MaxTermMyu ;
        }
        int startaux = this->EachAD[startAtom].Myu_StartNum;
        int auxlength = endaux - startaux ;
        {
            //Fl_GuessVector  FGV(UserVctMyu[i].Filename);
            TlVector FGV;
            FGV.load(UserVctMyu[i].Filename);

            //if (auxlength != FGV.getAuxterm() ){
            if (static_cast<TlVector::size_type>(auxlength) != FGV.getSize()) {
                log<<"Bad defined VectorFile [" <<UserVctMyu[i].Filename<<"]\n";
                CnErr.abort();
            }
        }
    }

    // for Nyu
    for (int i=0; i < this->UserVctNyuNum; i++) {
        int startAtom = UserVctNyu[i].From1 ;
        int finalAtom = UserVctNyu[i].To1;
        int endaux;
        if (finalAtom != AtomNum-1) {
            endaux   = this->EachAD[finalAtom+1].Nyu_StartNum ;
        } else {
            endaux = MaxTermMyu;
        }
        int startaux = this->EachAD[startAtom].Nyu_StartNum;
        int auxlength = endaux - startaux ;
        {
            //Fl_GuessVector  FGV(UserVctNyu[i].Filename);
            TlVector FGV;
            FGV.load(UserVctNyu[i].Filename);

            //if( auxlength != FGV.getAuxterm() ){
            if (static_cast<TlVector::size_type>(auxlength) != FGV.getSize()) {
                log<<"Bad defined VectorFile [" << UserVctNyu[i].Filename << "]\n";
                CnErr.abort();
            }
        }
    }

    return 0;
}

int DfInitialguess::setStruct()
{
    Fl_Geometry    FlGeom(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Density FlGtoDen;
    Fl_Gto_Xcpot   FlGtoXcpot;
    Fl_Gto_Xcpot2   FlGtoXcpot2;

    // DataSet : struct Auxtype for Rou.
    {
        int to =0;
        int count=0;
        for (int i=0; i<AtomKindNum; i++) {
            int hikisuu = FlGtoDen.getStartposition(i);
            this->Atomtype[i].AtomKind = FlGtoDen.getAtom(hikisuu);
            this->Atomtype[i].AuxSetName = FlGtoDen.getBasisName(hikisuu);

            this->Atomtype[i].Label2 = FlGtoDen.getLabel(hikisuu);

            for (int x=0; x<AtomNum; x++) {
                if ((Atomtype[i].AtomKind == FlGeom.getAtom(x)) &&
                        (Atomtype[i].Label2 == FlGeom.getLabel(x))) {
                    this->Atomtype[i].Label1 = FlGeom.getLabel(x);
                }
            }

            int from = to;
            to += FlGtoDen.getTermnumber(i);

            this->Atomtype[i].routerm.resize(4);
            for (int index = 0; index < 4; ++index) {
                this->Atomtype[i].routerm[index] = 0;
            }

            int cGto = 0;
            int point=0;
            for (int j = from; j < to; j++) {
                switch (FlGtoDen.getShell(j)) {
                case 's':
                    cGto=1;
                    this->Atomtype[i].routerm[1]++;
                    break;
                case 'p':
                    cGto=3;
                    this->Atomtype[i].routerm[2]++;
                    break;
                case 'd':
                    cGto=5;
                    this->Atomtype[i].routerm[3]++;
                    break;
                default:
                    GusError("Bad Shell Char !!");
                    break;
                }

                for (int m=0; m < cGto; m++) {
                    for (int k=0; k < FlGtoDen.getContraction(j); k++, point++)
                        this->Atomtype[i].alpha[point] = FlGtoDen.getExponent(j,k);
                }
            }

            count++;

            this->Atomtype[i].routerm[0] = this->Atomtype[i].routerm[1] * 1
                                           + this->Atomtype[i].routerm[2] * 3
                                           + this->Atomtype[i].routerm[3] * 5;
        }

    }

    // DataSet : struct Auxtype for Myu.
    {
        int to=0;
        //int count=0;
        for (int i=0; i < AtomKindNum; i++) {
            int from = to;
            to += FlGtoXcpot.getTermnumber(i);

            this->Atomtype[i].myuterm.resize(4);
            for (int index = 0; index < 4; ++index) {
                this->Atomtype[i].myuterm[index] = 0;
            }

            int cGto = 0;
            int point = 0;
            for (int j = from; j < to; j++) {
                switch (FlGtoXcpot.getShell(j)) {
                case 's':
                    cGto=1;
                    this->Atomtype[i].myuterm[1]++;
                    break;
                case 'p':
                    cGto=3;
                    this->Atomtype[i].myuterm[2]++;
                    break;
                case 'd':
                    cGto=5;
                    this->Atomtype[i].myuterm[3]++;
                    break;
                default:
                    GusError("Bad Shell Char !!");
                    break;
                }

                for (int m=0; m < cGto; m++) {
                    for (int k=0; k < FlGtoXcpot.getContraction(j); k++, point++) {
                        this->Atomtype[i].gamma[point] = FlGtoXcpot.getExponent(j,k);
                    }
                }
            }

            this->Atomtype[i].myuterm[0] =  this->Atomtype[i].myuterm[1] * 1
                                            + this->Atomtype[i].myuterm[2] * 3
                                            + this->Atomtype[i].myuterm[3] * 5;
        }
    }

    // DataSet : struct Auxtype for Nyu.
    {
        int to = 0;
        for (int i=0; i < AtomKindNum; i++) {
            int from = to;
            to += FlGtoXcpot2.getTermnumber(i);

            this->Atomtype[i].nyuterm.resize(4);
            for (int index = 0; index < 4; ++index) {
                this->Atomtype[i].nyuterm[index] = 0;
            }

            int cGto = 0;
            int point=0;
            for (int j = from; j < to; j++) {
                switch (FlGtoXcpot2.getShell(j)) {
                case 's':
                    cGto=1;
                    this->Atomtype[i].nyuterm[1]++;
                    break;
                case 'p':
                    cGto=3;
                    this->Atomtype[i].nyuterm[2]++;
                    break;
                case 'd':
                    cGto=5;
                    this->Atomtype[i].nyuterm[3]++;
                    break;
                default:
                    GusError("Bad Shell Char !!");
                    break;
                }

                for (int m=0; m < cGto; m++) {
                    for (int k=0; k < FlGtoXcpot2.getContraction(j); k++, point++) {
                        this->Atomtype[i].gamma2[point] = FlGtoXcpot2.getExponent(j,k);
                    }
                }
            }

            this->Atomtype[i].nyuterm[0] =  this->Atomtype[i].nyuterm[1] * 1
                                            + this->Atomtype[i].nyuterm[2] * 3
                                            + this->Atomtype[i].nyuterm[3] * 5;
        }
    }

    // ===  dataset : Number of struct_Auxtypefor each Atom.
    //int Protein_Nryoutai = FlGeom.getGroupnum();
    //  log<<" Protein_Nryoutai = [ "<<Protein_Nryoutai<<" ]"<<"\n";
    for (int j = 0; j < AtomNum; j++) {
        for (int i = 0; i < AtomKindNum; i++) {
            if ((FlGeom.getAtom(j) == this->Atomtype[i].AtomKind) &&
                    (FlGeom.getLabel(j) == this->Atomtype[i].Label2)) {

                this->EachAD[j].AuxTypeNumber = i ;

                if (strcmp(FlGeom.getLabel(j).c_str(),"") != 0) {

                    //char Label1[MaxLabelChar];
                    //strcpy(Label1, FlGeom.getLabel1(j).c_str());
                    std::string Label1 = FlGeom.getLabel(j);
                    if ((Label1[0]=='U')&&(Label1[1]=='D')&&(Label1[2]=='B')) {
                        this->EachAD[j].Unitfile = fl_User ;
                        break;
                    } else if ((Label1[0]=='D') && (Label1[1]=='B')) {
                        this->EachAD[j].Unitfile = fl_Molec ;
                        break;
                    } else {
                        this->EachAD[j].Unitfile = fl_Atom ;
                        break;
                    }
                } else {
                    this->EachAD[j].Unitfile = fl_Atom ;
                    break;
                }

                break;
            }
        }
    }

    return 0;
}


int DfInitialguess::makeStepTable()
{
    TlLogX& log = TlLogX::getInstance();

    Fl_Geometry FlGeom(Fl_Geometry::getDefaultFileName());

    int count = 0;
    for (int i=0; i<AtomNum;) {
        //     if( OutLevel < -4 ){
        //       log<<"AtomNumber = "<< i <<"     ";
        //     }
        int atmnum  = 0;
        //int dumynum = EachAD[i].AuxTypeNumber;
        {
            Fl_Gdb_Molecular FGM;

            switch (this->EachAD[i].Unitfile) {
            case fl_Atom:
                count++ ;
                atmnum  = 1;
                break;

            case fl_Namino:
                atmnum  = FGM.getHeaderAtomNum(GDBfile1, "", 0);
                count++;
                break;

            case fl_Pepamino:
                atmnum  = FGM.getHeaderAtomNum(GDBfile2, "", 0);
                count++;
                break;

            case fl_Camino:
                atmnum  = FGM.getHeaderAtomNum(GDBfile3, "", 0);
                count++;
                break;

            case fl_Molec:
                atmnum  = FGM.getHeaderAtomNum(GDBfile4, "", 0);
                count++;
                log<<"MOLEC"<<"\n";
                break;

            case fl_User:
                atmnum  = FGM.getHeaderAtomNum(GDBfile5, "", 0);
                count++;
                break;

            default:
                log<<"Bad information UnitFile in makeStep Table";
                CnErr.abort();
                break;
            }
        }
        i += atmnum ;
    }

    STEPnumber = count;
    if (STEPnumber > MaxStepNum) {
        log<<"Error in makeStepTable in DfInitialguess.\n";
        log<<"You must change MaxStepNum in DfInitialguess.h\n";
        log<<"Now , MaxStepNum = " << static_cast<int>(MaxStepNum) <<"\n";
        log<<"So , MaxStepNum > "  << STEPnumber <<"\n";
    }

    count=0;
    for (int i=0; i<AtomNum;) {
        Fl_Gdb_Molecular  FGM;
        int atmnum = 0;
        {
            switch (EachAD[i].Unitfile) {
            case fl_Atom: {
                this->GuessStep[count].Unitfile = fl_Atom ;
                int dumynum = EachAD[i].AuxTypeNumber;
                this->GuessStep[count].DbLabel = this->Atomtype[dumynum].AuxSetName;
                atmnum  = 1;
                this->GuessStep[count].AtomNum = atmnum;
                this->GuessStep[count].FirstAtomNum = i ;
                this->GuessStep[count].UserVctRouExist = serchUserVctExist(Rou,i,i);
                this->GuessStep[count].UserVctMyuExist = serchUserVctExist(Myu,i,i);
                this->GuessStep[count].UserVctNyuExist = serchUserVctExist(Nyu,i,i);
                count++ ;
            }
            break;
            case fl_Namino: {
                this->GuessStep[count].Unitfile = fl_Namino ;
                this->GuessStep[count].DbLabel = ""; //FlGeom.getResidue(i);
                atmnum  = FGM.getHeaderAtomNum(GDBfile1, "", 0);
                this->GuessStep[count].AtomNum = atmnum ;
                this->GuessStep[count].FirstAtomNum = i ;
                this->GuessStep[count].UserVctRouExist = serchUserVctExist(Rou,i,i+atmnum-1);
                this->GuessStep[count].UserVctMyuExist = serchUserVctExist(Myu,i,i+atmnum-1);
                this->GuessStep[count].UserVctNyuExist = serchUserVctExist(Nyu,i,i+atmnum-1);
                count++ ;
            }
            break;
            case fl_Pepamino: {
                this->GuessStep[count].Unitfile = fl_Pepamino ;
                this->GuessStep[count].DbLabel = ""; //FlGeom.getResidue(i);
                atmnum  = FGM.getHeaderAtomNum(GDBfile2, "", 0);
                this->GuessStep[count].AtomNum = atmnum ;
                this->GuessStep[count].FirstAtomNum = i ;
                this->GuessStep[count].UserVctRouExist = serchUserVctExist(Rou,i,i+atmnum-1);
                this->GuessStep[count].UserVctMyuExist = serchUserVctExist(Myu,i,i+atmnum-1);
                this->GuessStep[count].UserVctNyuExist = serchUserVctExist(Nyu,i,i+atmnum-1);
                count++ ;
            }
            break;
            case fl_Camino: {
                this->GuessStep[count].Unitfile = fl_Camino ;
                this->GuessStep[count].DbLabel = ""; //FlGeom.getResidue(i);
                atmnum  = FGM.getHeaderAtomNum(GDBfile3, "", 0);
                this->GuessStep[count].AtomNum = atmnum ;
                this->GuessStep[count].FirstAtomNum = i ;
                this->GuessStep[count].UserVctRouExist = serchUserVctExist(Rou,i,i+atmnum-1);
                this->GuessStep[count].UserVctMyuExist = serchUserVctExist(Myu,i,i+atmnum-1);
                this->GuessStep[count].UserVctNyuExist = serchUserVctExist(Nyu,i,i+atmnum-1);
                count++ ;
            }
            break;
            case fl_Molec: {
                this->GuessStep[count].Unitfile = fl_Molec ;
                this->GuessStep[count].DbLabel = FlGeom.getLabel(i);
                atmnum  = FGM.getHeaderAtomNum(GDBfile4,FlGeom.getLabel(i).c_str(),0);
                this->GuessStep[count].AtomNum = atmnum ;
                this->GuessStep[count].FirstAtomNum = i ;
                this->GuessStep[count].UserVctRouExist = serchUserVctExist(Rou,i,i+atmnum-1);
                this->GuessStep[count].UserVctMyuExist = serchUserVctExist(Myu,i,i+atmnum-1);
                this->GuessStep[count].UserVctNyuExist = serchUserVctExist(Nyu,i,i+atmnum-1);
                count++ ;
            }
            break;
            case fl_User: {
                this->GuessStep[count].Unitfile = fl_User ;
                this->GuessStep[count].DbLabel = FlGeom.getLabel(i);
                atmnum  = FGM.getHeaderAtomNum(GDBfile5,FlGeom.getLabel(i).c_str(),0);
                this->GuessStep[count].AtomNum = atmnum ;
                this->GuessStep[count].FirstAtomNum = i ;
                this->GuessStep[count].UserVctRouExist = serchUserVctExist(Rou,i,i+atmnum-1);
                this->GuessStep[count].UserVctMyuExist = serchUserVctExist(Myu,i,i+atmnum-1);
                this->GuessStep[count].UserVctNyuExist = serchUserVctExist(Nyu,i,i+atmnum-1);
                count++ ;
            }
            break;
            default: {
                log<<"Bad information UnitFile in makeStep Table";
                CnErr.abort();
            }
            break;
            }
        }
        i += atmnum;
    }

    if (count != STEPnumber) {
        log<<"Bad counter \n";
        CnErr.abort();
    }

    return 0;
}

DfInitialguess::KeyOnOff DfInitialguess::serchUserVctExist(WhichRMN RMN, int from, int to)
{
    TlLogX& log = TlLogX::getInstance();

    KeyOnOff  Hantei;
    int flag,i;

    flag=0;
    Hantei=OFF;

    if (RMN==Rou) {
        for (i=0; i<UserVctRouNum; i++) {
            if (((UserVctRou[i].From1<=from) && (from<UserVctRou[i].To1)) &&
                    ((UserVctRou[i].From1<to)   && (to  <=UserVctRou[i].To1))) {
                flag++;
            }
        }
        if (flag==0) {
            Hantei=OFF;
        } else if (flag==1) {
            Hantei=ON;
        } else if (flag>1) {
            log<<"Bad appointment UserVectorRou in serchUserVctExist"<<"\n";
            CnErr.abort();
        }
        return  Hantei;
    }

    if (RMN==Myu) {
        for (i=0; i<UserVctMyuNum; i++) {
            if (((UserVctMyu[i].From1<from) && (from<UserVctMyu[i].To1)) &&
                    ((UserVctMyu[i].From1<to)   && (to  <UserVctMyu[i].To1))) {
                flag++;
            }
        }
        if (flag==0) {
            Hantei=OFF;
        } else if (flag==1) {
            Hantei=ON;
        } else if (flag>1) {
            log<<"Bad appointment UserVectorMyu in serchUserVctExist"<<"\n";
            CnErr.abort();
        }
        return  Hantei;
    }

    if (RMN==Nyu) {
        for (i=0; i<UserVctNyuNum; i++) {
            if (((UserVctNyu[i].From1<from) && (from<UserVctNyu[i].To1)) &&
                    ((UserVctNyu[i].From1<to)   && (to  <UserVctNyu[i].To1))) {
                flag++;
            }
        }
        if (flag==0) {
            Hantei=OFF;
        } else if (flag==1) {
            Hantei=ON;
        } else if (flag>1) {
            log<<"Bad appointment UserVectorNyu in serchUserVctExist"<<"\n";
            //        log<<"\7"<<"\n";
            CnErr.abort();
        }
        return  Hantei;
    }
    return  Hantei;
}

int DfInitialguess::setGDBdatafile()
{

    // ここでは、"./fl_GuessDB"の下のファイル
    // fl_Gdb_Namino,fl_Gdb_Peptid,fl_Gdb_Camino,fl_Gdb_Molecular,
    // fl_Gdb_UserMolec のデータのうち、計算に必要なDBラベルのものだけを
    // 取り出す操作を行なう。
    // このとき、使用する補助関数セットが違う場合は、
    // コレスポンディング・オービタルの方法によって、展開係数を変換しておく。
    // 変換先のデータファイルは、それぞれ以下のファイル。
    // fl_GuessDB/tmp の下の、fl_Gdb_Nanimo2,fl_Gdb_Peptid2,fl_Gdb_Camino2,
    // fl_Gdb_Molecular2,fl_Gdb_UserMolec2 である。

    // 以上は、まだ計画の段階で、ver1.00 ではサポートしていない。

    // また、原子のデータベースfl_GuessDB/fl_Gdb_Atomについては、
    // このままのファイルを扱っている。
//   if (OutLevel < 0){
//     log << "This routine make Temp_GuessDB from fl_GuessDB/files .\n";
//     log << "Now , This member function is empty\n";
//   }

    CorrespOrbital(); // This is dumy.


    return 0;
}

int DfInitialguess::CorrespOrbital()
{
    return 0;
}

// fl_Gdb_Atomに入っているアトムのデータをメモリ上（構造体）に読み込む。
int DfInitialguess::setAtomCoef()
{
    TlLogX& log = TlLogX::getInstance();

    //char SetName[LineBuf];
    Fl_Gto_Density FlGtoDen;

    std::ifstream fi;
    fi.open(this->GDBfile0.c_str());
    if (!fi) {
        log <<"Cannot Open in setAtomcoef() "<< GDBfile0 << "\n";
        log.flush();
    }

    for (int i = 0; i < AtomKindNum; ++i) {
        //strcpy(SetName , this->Atomtype[i].AuxSetName.c_str());
        std::string SetName = this->Atomtype[i].AuxSetName;
        Fl_Gdb_Atom FGA;
        FGA.readData(GDBfile0.c_str() ,GDBfile0.c_str(), SetName, 0);
        log << TlUtils::format(" file: %s, name=%s\n",
                               GDBfile0.c_str(), SetName.c_str())
            << std::endl;
        
        // setdata for Rou.
        tempflagw1 = FGA.gettempflag();
        if (tempflagw1 == 1) {
            this->tempflag = 1;
            const int startnumber = FlGtoDen.getStartposition(i);
            const int atomcgtonumber = FlGtoDen.getTermnumber(i);

            log<< "\n";
            log<< "**** WARNING ****" <<"\n";
            log<< "ATOM = " <<Atomtype[i].AtomKind <<"\n";
            log<< "The basis set was not included in fl_Gdb_Atom"<<"\n";
            log<< "However, force to continue"<<"\n";
            log<< "**** INFORMATION FROM BASIS2 ****" <<"\n";
            log<< "Startposition = "<<startnumber<<"  Termnumber = "<<atomcgtonumber <<"\n";
            log<< "\n";

            std::cout<<std::endl;
            std::cout<< "**** WARNING ****" <<std::endl;
            std::cout<< "ATOM = "<<Atomtype[i].AtomKind <<std::endl;
            std::cout<< "The basis set was not included in fl_Gdb_Atom"<<std::endl;
            std::cout<< "However, force to continue"<<std::endl;
            std::cout<< "**** INFORMATION FROM BASIS2 ****" <<std::endl;
            std::cout<< "Startposition = "<<startnumber<<"  Termnumber = "<<atomcgtonumber <<std::endl;
            std::cout<<std::endl;

            // force to set the values of elements of rho
            log<<" only the coefficient of s-shell term is set to 1.0 "<<"\n";
            int itmp, rhotermnum; // ,inum

            int snumber = 0;
            int pnumber = 0;
            int dnumber = 0;
            int inum2 = 0;
            for (int inum = startnumber; inum < startnumber + atomcgtonumber; inum++) {
                const char shellinfo = FlGtoDen.getShell(inum);
                switch (shellinfo) {
                case 's':
                    Atomtype[i].rou[inum2]=1.0;
                    snumber=snumber+1;
                    inum2+=1;
                    log<<"s is detected"<<"\n";
                    break;
                case 'p':
                    for (itmp=inum2; itmp<inum2+3; itmp++) Atomtype[i].rou[itmp]=0.0;
                    pnumber=pnumber+1;
                    inum2+=3;
                    break;
                case 'd':
                    for (itmp=inum2; itmp<inum2+5; itmp++) Atomtype[i].rou[itmp]=0.0;
                    dnumber=dnumber+1;
                    inum2+=5;
                    break;
                default:
                    CnErr.abort("DfInitialguess","","","Something wrong in basis information");
                    break;
                }
            }

            rhotermnum=snumber+pnumber*3+dnumber*5;
            Atomtype[i].routerm[0] = rhotermnum;
            Atomtype[i].rDbShellNum[0]=snumber;
            Atomtype[i].rDbShellNum[1]=pnumber;
            Atomtype[i].rDbShellNum[2]=dnumber;
            Atomtype[i].rDbShellNum[3]=0;
        } else {
            if (Atomtype[i].routerm[0] != FGA.getRouterm(0)) {
                log << TlUtils::format("Bad term Number: Atomtype[%d].routerm != FGA.getRouterm(0)", i)
                    << std::endl;
                log << Atomtype[i].routerm[0] << ", " << FGA.getRouterm(0) <<"\n";
                log.flush();
                CnErr.abort();
            } else {
                FGA.getRouCoef(0 , Atomtype[i].rou);
                Atomtype[i].rDbShellNum = FGA.getRouShellterm(0);
            }

            // setdata for Myu.
            if (Atomtype[i].myuterm[0] != FGA.getMyuterm(0)) {
                log<<"Bad term Number "<<"\n";
                log<<"Atomtype["<<i<<"].myuterm != FGA.getMyuterm(0)"<<"\n";
                log << Atomtype[i].myuterm[0] << "\n";
                log << FGA.getMyuterm(0) << "\n";
                log.flush();

                CnErr.abort();
            } else {
                FGA.getMyuCoef(0 , Atomtype[i].myu);
                Atomtype[i].mDbShellNum = FGA.getMyuShellterm(0);
            }

            // setdata for Nyu.
//             if (Atomtype[i].nyuterm[0] != FGA.getNyuterm(0)) {
//                 log << TlUtils::format("Bad term Number: Atomtype[%d].nyuterm != FGA.getNyuterm(0)\n", i);
//                 log << "lhs = " << Atomtype[i].nyuterm[0] << "\n";
//                 log << "rhs = " << FGA.getNyuterm(0) << "\n";
//                 log.flush();
//             } else {
//                 FGA.getNyuCoef(0 , Atomtype[i].nyu);
//                 Atomtype[i].nDbShellNum = FGA.getNyuShellterm(0);
//             }
        }
    }

    fi.close();
    return 0;
}
//#############################################################################
int DfInitialguess::calcMain()
{
    TlLogX& log = TlLogX::getInstance();

    if (this->PpqFname != "") {
        if (scftype != NSP) {
            GusError("Bad appointment MatrixInput_PpqFname,now SP calulation");
        }
        readPpq();
        callDensityfitting();
        makeXcpotential();
    } else if ((this->AlphaPpqFname != "") && (this->BetaPpqFname != "")) {
        if (scftype != SP) {
            GusError("Bad appointment MatrixInput,now NSP calulation");
        }

        if ((this->AlphaPpqFname == "BINARY") &&
                (this->BetaPpqFname  == "BINARY")) {
            // バイナリファイル fl_Mtr_PpqDra,fl_Mtr_PpqDrb,を直接扱う。
            callDensityfitting();
            makeXcpotential();
            //goto end;
        } else {
            readPpq();
            callDensityfitting();
            makeXcpotential();
        }
    } else {
        for (int i=0; i<STEPnumber; i++) {
            if ((GuessStep[i].UserVctRouExist==OFF) ||
                    (GuessStep[i].UserVctMyuExist==OFF) ||
                    (GuessStep[i].UserVctNyuExist==OFF)) {

                if (GuessStep[i].Unitfile == fl_Atom) {
                    this->EXE_ATOM(i);
                } else {
                    this->EXE_MOLECULE(i);
                }
            }
        }// for(i)

        if (w1RouCount != MaxTermRou) {
            log<<" Bad counter and vector in EXE_MOLECULE_routine(). "<<"\n";
            log<<" w1RouCount = "<<w1RouCount<<"\n";
            log<<" MaxTermRou = "<<MaxTermRou<<"\n";
            //    log<<"\7"<<"\n";
            CnErr.abort();
        }
        if (w1MyuCount != MaxTermMyu) {
            log<<" Bad counter and vector in EXE_MOLECULE_routine(). "<<"\n";
            log<<" w1MyuCount = "<<w1MyuCount<<"\n";
            log<<" MaxTermMyu = "<<MaxTermMyu<<"\n";
            //    log<<"\7"<<"\n";
            CnErr.abort();
        }
//         if (w1NyuCount != MaxTermNyu) {
//             log<<" Bad counter and vector in EXE_MOLECULE_routine(). "<<"\n";
//             log<<" w1NyuCount = "<<w1NyuCount<<"\n";
//             log<<" MaxTermNyu = "<<MaxTermNyu<<"\n";
//         }

        //-----------------------------------------
        // Cover Vector which is defined by user.
        //-----------------------------------------
        CoverVector();
    }

    return 0;

}

void DfInitialguess::readPpq()
{
    const int niter = 0;  // output file number :  fl_Mtr_PpqDr0.
    // or fl_Mtr_PpqDra0 , fl_Mtr_PpqDrb0 .

    if (this->scftype == NSP) {
        //Fl_GuessMatrix FGMAT(PpqFname);
        TlSymmetricMatrix FGMAT;
        FGMAT.load(PpqFname);
        //TlSymmetricMatrix PpqD = FGMAT.getValue();
        //PpqD.save("fl_Work/fl_Mtr_PpqDr" + TlUtils::xtos(niter));
        FGMAT.save("fl_Work/fl_Mtr_PpqDr" + TlUtils::xtos(niter));
    } else {
        //============================================================
        // for scftype==SP ;

        //----------------------
        // for ALPHA DENSITY
        //----------------------
        {
            //Fl_GuessMatrix  FGMAT(AlphaPpqFname);
            TlSymmetricMatrix FGMAT;
            FGMAT.load(AlphaPpqFname);
            //TlSymmetricMatrix PpqDa = FGMAT.getValue();
            //PpqDa.save("fl_Work/fl_Mtr_PpqDra" + TlUtils::xtos(niter));
            FGMAT.save("fl_Work/fl_Mtr_PpqDra" + TlUtils::xtos(niter));
        }

        //----------------------
        // for BETA DENSITY
        //----------------------
        {
            //Fl_GuessMatrix  FGMAT(BetaPpqFname);
            TlSymmetricMatrix FGMAT;
            FGMAT.load(BetaPpqFname);
            //TlSymmetricMatrix PpqDb = FGMAT.getValue();
            //PpqDb.save("fl_Work/fl_Mtr_PpqDrb" + TlUtils::xtos(niter));
            FGMAT.save("fl_Work/fl_Mtr_PpqDrb" + TlUtils::xtos(niter));
        }
    }
}

void DfInitialguess::callDensityfitting()
{
    TlLogX& log = TlLogX::getInstance();
    log<<" ### START DensityFitting Called in DfInitialguess ###"<<"\n";

    //int Rnum;
    {
        this->densityFitOnly();
    }

    if (scftype==NSP) {
        TlVector rho1;
        rho1.load("fl_Work/fl_Vct_Rou1");
        rho1.save("fl_Work/fl_Vct_Rou0");
    } else if (scftype == SP) {
        {
            TlVector rho1;
            rho1.load("fl_Work/fl_Vct_Roua1");
            rho1.save("fl_Work/fl_Vct_Roua0");

            rho1.load("fl_Work/fl_Vct_Roub1");
            rho1.save("fl_Work/fl_Vct_Roub0");
        }
    }
}

void DfInitialguess::makeXcpotential()
{
    TlLogX& log = TlLogX::getInstance();

    if (tempflag==1) {
        log<<"### makeXcpotential in DfInitialguess is skipped ###"<<"\n";
        return;
    }

    if (MethodMyuNyu == Meth0) {
        for (int i=0; i < STEPnumber; i++) {
            this->EXE_ATOM(i);
        }
        if (scftype == NSP) {
            this->CoefMyu = this->w1Myu;
            this->CoefNyu = this->w1Nyu;
        } else if (scftype == SP) {
            this->CoefMyuAlpha = this->w1Myu;
            this->CoefMyuBeta  = this->w1Myu;
            this->CoefNyuAlpha = this->w1Nyu;
            this->CoefNyuBeta  = this->w1Nyu;
        }
    } else if (MethodMyuNyu == Meth1) {
        if (MaxTermRou != MaxTermMyu) {
            log<<"This routine cannot use this calculation"<<"\n";
            log<<"Because , Total Auxiliary Number for Rou is "
            <<"differrent from Total Auxiliary Number for Myu."<<"\n";
            log<<"So, you must appoint other method to make Xcpotential"
            <<"\n";
            CnErr.abort();
        }
        this->CoefRou.save("fl_Work/fl_Vct_Rou1");

        if (scftype==NSP) {
            this->CoefMyu = -1.0 * this->CoefRou;
            this->CoefNyu = this->CoefRou;
        } else if (scftype == SP) {
            this->CoefMyuAlpha = -1.0 * this->CoefRou;
            this->CoefMyuBeta  = -1.0 * this->CoefRou;
            this->CoefNyuAlpha = this->CoefRou;
            this->CoefNyuBeta  = this->CoefRou;
        }
    } else if (MethodMyuNyu == Meth2) {
        log<<"Apporoximation(Meth0) for All Space.(fitting);"<<"\n";
        log<<"This routine is dumy menber function."<<"\n";
        log<<"So, you must appoint other method to make Xcpotential"<<"\n";
        //  log<<"\7"<<"\n";
        CnErr.abort();
    } else if (MethodMyuNyu == Meth3) {
        if (MaxTermRou != MaxTermMyu) {
            log<<"This routine cannot use this calculation"<<"\n";
            log<<"Because , Total Auxiliary Number for Rou is "
            <<"differrent from Total Auxiliary Number for Myu."<<"\n";
            log<<"So, you must appoint other method to make Xcpotential"
            <<"\n";
            CnErr.abort();
        }
        {
            this->CoefMyu.load("fl_Work/fl_Vct_Myu999");
            this->CoefNyu.load("fl_Work/fl_Vct_Nyu999");

            if (scftype == SP) {
                this->CoefMyuAlpha = this->CoefMyu * 0.5;
                this->CoefMyuBeta  = this->CoefMyu * 0.5;
                this->CoefNyuAlpha = this->CoefNyu * 0.5;
                this->CoefNyuBeta  = this->CoefNyu * 0.5;
            }
        }
    } else if (MethodMyuNyu == Meth4) {
        log<<"This routine is dumy menber function."<<"\n";
        log<<"So, you must appoint other method to make Xcpotential"<<"\n";
        CnErr.abort();
    }
}

void DfInitialguess::EXE_ATOM(int num)
{
    const int AtmNum  = GuessStep[num].FirstAtomNum;
    const int dumynum = this->EachAD[AtmNum].AuxTypeNumber;

    const int rStartNum = this->EachAD[AtmNum].Rou_StartNum;
    const int Rtm = Atomtype[dumynum].routerm[0];

    this->w1Rou.resize(rStartNum + Rtm);
    for (int i=0; i < Rtm; i++) {
        this->w1Rou[ rStartNum + i ] = Atomtype[dumynum].rou[i] ;
    }
    w1RouCount += Rtm ;

    const int mStartNum = this->EachAD[AtmNum].Myu_StartNum;
    const int Mtm = Atomtype[dumynum].myuterm[0];

    this->w1Myu.resize(mStartNum + Mtm);
    for (int i=0; i < Mtm; i++) {
        this->w1Myu[ mStartNum + i ] = Atomtype[dumynum].myu[i] ;
    }
    w1MyuCount += Mtm ;

    const int nStartNum = this->EachAD[AtmNum].Nyu_StartNum;
    const int Ntm = Atomtype[dumynum].nyuterm[0];

    this->w1Nyu.resize(nStartNum + Ntm);
    for (int i=0; i < Ntm; i++) {
        this->w1Nyu[ nStartNum + i ] = Atomtype[dumynum].nyu[i] ;
    }
    w1NyuCount += Ntm ;

    Consider_Connect(); // This is dumy. Content is empty.
}

void DfInitialguess::EXE_MOLECULE(int num)
{
    TlLogX& log = TlLogX::getInstance();

    int  tAN;  // tAN:totalAtomNumber of One Molecule(Residue) ;
    int  count,showcountR,showcountM,showcountN;
    int j,k,CnNum; // int r,i,,dumynum2, ,dumynum1;
    int atmRouterm,atmMyuterm,atmNyuterm; //,stdnum;
    std::string MolecDBLabel;

    std::string filename;  // Honmono DB file.
    std::string filename2; // Header DB file.

    Fl_Geometry FlGeom(Fl_Geometry::getDefaultFileName());
    showcountR=0;
    showcountM=0;
    showcountN=0;

    Fl_Gdb_Molecular  FGM;  // OBJECT Construct.

    // ---------- appointment GuessDBfile.------------------
    if (GuessStep[num].Unitfile == fl_Namino) {
        filename  = this->GDBfile1_2;
        filename2 = this->GDBfile1;
    } else if (GuessStep[num].Unitfile == fl_Pepamino) {
        filename  = this->GDBfile2_2;
        filename2 = this->GDBfile2;
    } else if (GuessStep[num].Unitfile == fl_Camino) {
        filename  = this->GDBfile3_2;
        filename2 = this->GDBfile3;
    } else if (GuessStep[num].Unitfile == fl_Molec) {
        filename  = this->GDBfile4_2;
        filename2 = this->GDBfile4;
    } else if (GuessStep[num].Unitfile == fl_User) {
        filename  = this->GDBfile5_2;
        filename2 = this->GDBfile5;
    }

    MolecDBLabel = GuessStep[num].DbLabel;

    // Data Initialize for struct AtomData and struct AminoData.
    for (j=0; j<Max1DbAtomNum; j++) {
        DbAmino.Order[j]  = 0;
        InpAmino.Order[j] = 0;
        for (k=0; k<(4*4); k++) {
            DbAmino.Atomdata[j].afinmat[k]=0;
        }
        for (k=0; k<aboutAuxTerm; k++) {
            DbAmino.Atomdata[j].CoefRou[k]=0;
            DbAmino.Atomdata[j].CoefMyu[k]=0;
            DbAmino.Atomdata[j].CoefNyu[k]=0;
        }
        for (k=0; k<MaxSPDterm; k++) {
            DbAmino.Atomdata[j].rouDbShellNum[k]=0;
            DbAmino.Atomdata[j].myuDbShellNum[k]=0;
            DbAmino.Atomdata[j].nyuDbShellNum[k]=0;
        }
        for (k=0; k<MaxCntNum; k++) {
            DbAmino.Atomdata[j].Conection[k]=0;
        }
    }

    // CHANGE
    //    FGM.readData(filename,MolecDBLabel,1); // readData from OBJECT.
    FGM.readData(filename.c_str(), filename2.c_str(), MolecDBLabel.c_str(),1); // readData from OBJECT.

    tAN = InpAmino.Res1AtomNum = FGM.getTotalAtomNum();
    if (tAN > Max1DbAtomNum) {
        log<<" Bad : EXE_MOLECULE "<<"\n";
        log<<" You must change the number(Max1DbAtomNum)";
        log<<" in DfInitialguess.h ."<<"\n";
        log<<" You must chage Max1DbAtomNum > "<< tAN <<"\n";
        CnErr.abort();
    }
    DbAmino.Res1AtomNum = InpAmino.Res1AtomNum ;

    DbAmino.equalPairNum = FGM.getTotalEqualPairNum();
    if (DbAmino.equalPairNum > MaxEqualPair) {
        log<<" Bad : EXE_MOLECULE "<<"\n";
        log<<" You must change the number(MaxEqualPair)";
        log<<" in DfInitialguess.h ."<<"\n";
        log<<" You must chage MaxEqualPair > ";
        log<< DbAmino.equalPairNum <<"\n";
        CnErr.abort();
    }

    DbAmino.TotalRouterm = FGM.getTotalRouterm();
    DbAmino.TotalMyuterm = FGM.getTotalMyuterm();
    DbAmino.TotalNyuterm = FGM.getTotalNyuterm();
    for (j=0; j<tAN; j++) {
        DbAmino.Atomdata[j].Routerm = FGM.getRouterm(j);
        DbAmino.Atomdata[j].Myuterm = FGM.getMyuterm(j);
        DbAmino.Atomdata[j].Nyuterm = FGM.getNyuterm(j);
    }

    for (j=0; j<tAN; j++) {
        atmRouterm = DbAmino.Atomdata[j].Routerm;
        if (atmRouterm > aboutAuxTerm) {
            log<<" Bad : EXE_MOLECULE "<<"\n";
            log<<" You must change the number(aboutAuxTerm)";
            log<<" in DfInitialguess.h ."<<"\n";
            log<<" You must chage aboutAuxTerm > "<< atmRouterm <<"\n";
            CnErr.abort();
        }

        atmMyuterm = DbAmino.Atomdata[j].Myuterm;
        if (atmMyuterm > aboutAuxTerm) {
            log<<" Bad : EXE_MOLECULE "<<"\n";
            log<<" You must change the number(aboutAuxTerm)";
            log<<" in DfInitialguess.h ."<<"\n";
            log<<" You must chage aboutAuxTerm > "<< atmMyuterm <<"\n";
            CnErr.abort();
        }

        atmNyuterm = DbAmino.Atomdata[j].Nyuterm;
        if (atmNyuterm > aboutAuxTerm) {
            log<<" Bad : EXE_MOLECULE "<<"\n";
            log<<" You must change the number(aboutAuxTerm)";
            log<<" in DfInitialguess.h ."<<"\n";
            log<<" You must chage aboutAuxTerm > "<< atmNyuterm <<"\n";
            CnErr.abort();
        }

        CnNum = DbAmino.Atomdata[j].ConectNum = FGM.getConnectNum(j);
        if (CnNum > MaxCntNum) {
            log<<" Bad : EXE_MOLECULE "<<"\n";
            log<<" You must change the number(MaxCntNum)";
            log<<" in DfInitialguess.h ."<<"\n";
            log<<" You must chage MaxCntNum > "<< CnNum <<"\n";
            CnErr.abort();
        }

    }

    // Get Value from OBJECT.
    FGM.getOrder(DbAmino.Order);
    FGM.getEqualPairCount(DbAmino.equalPaircount);
    FGM.getEqualPair(DbAmino.equalPair);
    for (j=0; j<tAN; j++) {
        DbAmino.Atomdata[j].auxname = FGM.getAuxSetName(j);
        DbAmino.Atomdata[j].atomname = FGM.getAtom(j);
        DbAmino.Atomdata[j].x = FGM.getXCoord(j);
        DbAmino.Atomdata[j].y = FGM.getYCoord(j);
        DbAmino.Atomdata[j].z = FGM.getZCoord(j);
        FGM.getRouCoef(j, DbAmino.Atomdata[j].CoefRou);
        FGM.getMyuCoef(j, DbAmino.Atomdata[j].CoefMyu);
        FGM.getNyuCoef(j, DbAmino.Atomdata[j].CoefNyu);
        DbAmino.Atomdata[j].rouDbShellNum = FGM.getRouShellterm(j);
        DbAmino.Atomdata[j].myuDbShellNum = FGM.getMyuShellterm(j);
        DbAmino.Atomdata[j].nyuDbShellNum = FGM.getNyuShellterm(j);
        FGM.getConnection(j , DbAmino.Atomdata[j].Conection);
    }

    // Initialyze for AfinMatrix.
    for (j=0; j<tAN; j++) {
        for (k=0; k<16; k++) {
            DbAmino.Atomdata[j].afinmat[k]=0.0;
        }
        for (k=0; k<4; k++) {
            DbAmino.Atomdata[j].afinmat[k*4+k]=1.0;
        }
    }

    //------------ 入力データの読み込み。------
    //---注意:アミノ酸同士の間にダミー原子を置いてはいけない。
    //---    :タンパクを解く時は○量体の頭にのみつけてよいことにする
    count = GuessStep[num].FirstAtomNum;
    for (j=0; j<tAN; j++) {
        InpAmino.Atomdata[j].atomname = FlGeom.getAtom(j+count);
        InpAmino.Atomdata[j].x = FlGeom.getCoordinate(j+count).x();
        InpAmino.Atomdata[j].y = FlGeom.getCoordinate(j+count).y();
        InpAmino.Atomdata[j].z = FlGeom.getCoordinate(j+count).z();
        InpAmino.Order[j]      = 0 ;
        // Order :  1 --->   use afinmatrix at TransOrbital ;
        // Order : -1 ---> nouse afinmatrix at TransOrbital ;
        //                 Input Only S_orbital at TransOrbital ;
    }


    // 入力された分子の原子並びと、指定されたDB分子の原子並びが一致しているかを
    // チェックする。

    CheckDBFormat(num);

    //=============================================================
    //----------------------------------------------------------------
    //----- r量体i番目の残基の構造を一致させるアフィンマトリクスを作る。
    //----------------------------------------------------------------

    makeAfinmatrix(num);
    transOrbital(tAN);
}

void DfInitialguess::CheckDBFormat(int num)
{
    TlLogX& log = TlLogX::getInstance();

    for (int i=0; i<DbAmino.Res1AtomNum; i++) {
        if (DbAmino.Atomdata[i].atomname != InpAmino.Atomdata[i].atomname) {
            //        log<<" ### Bad Format ### " <<"\7"<<"\n";
            log<<" ### Bad Format ### " <<"  "<<"\n";
            log<<" STEPnumber = "<< num << "\n";
            int dumynum1 = GuessStep[num].FirstAtomNum;
            //int dumynum2 = this->EachAD[dumynum1].AuxTypeNumber;
            log<<" The DB_Molecule's fitst number = "<<dumynum1<<"\n";
            //CHANGE   log<<" wrong DB_Label = "<< Atomtype[dumynum2].Label1<<"\n";
            log<<" wrong AtomNAME = "<<InpAmino.Atomdata[i].atomname<<"\n";
            log<<" It's number is "<< i <<"\n";
            CnErr.abort();
        }
    }
}

int DfInitialguess::makeAfinmatrix(int stepnum)
{
/////////////////////////////////////////////////////////////////////////
//第１次バージョンとしてある原子を考慮するのに、その原子に結合している
//原子群の座標のずれを最も小さいようにアフィンマトリクスを決める。
//ただし、ペプチド結合部を考慮してないことになるので
//そこもなんとかしないといけない。
//等価な原子のペアは、考慮している原子１つに対して、
//１ペアということにしている
////////////////////////////////////////////////////////////////////////
    TlLogX& log = TlLogX::getInstance();

//   if( OutLevel < -3 ){
//     log<<"\n";
//     log<<"    -------------------------------------------"<<"\n";
//     log<<"    -----   makeAfinmatrix(int)     Start -----"<<"\n";
//     log<<"    -------------------------------------------"<<"\n";
//     log<<"\n";
//   }

    enum { AfinMatrixYouso = 4*4 }; // アフィンマトリクスの要素数。
    enum { MaxToukaAtomNum = 6 }; //１原子の等価原子の最大数。
    enum { MaxMawasi = 720 };     // (MaxToukaAtomNum)! となる数。
    enum { MaxUVN = 3 };          // 未知数の数。
    //現在は最適な回転角（θx、θy、θz）の３つ。

    //double CenterAtomX,CenterAtomY,CenterAtomZ;
    //double InpCtrAtomX,InpCtrAtomY,InpCtrAtomZ;
    double ConectAtomX[MaxCntNum];
    double ConectAtomY[MaxCntNum];
    double ConectAtomZ[MaxCntNum];
    double InpAtomX[MaxCntNum];
    double InpAtomY[MaxCntNum];
    double InpAtomZ[MaxCntNum];
    double dmX1[MaxCntNum];
    double dmY1[MaxCntNum];
    double dmZ1[MaxCntNum];
    double dmX2[MaxCntNum];
    double dmY2[MaxCntNum];
    double dmZ2[MaxCntNum];

    int ConectSerialAtomNum[MaxToukaAtomNum];
    // ConectSerialAtomNumには等価原子の並びが
    //ConectAtom{X,Y,Z}の何番目にあたるかを保存。

    // 等価な原子がある場合に、何通りの並び替えが
    // あるか。n原子のときには、mawasi = n! 通り。
    int mawasi =0;

    int PairBangou =0;   // 等価なペアのうち、何番めのペアか。
    int VectBangou;   // ベクトルDbAmino.equalPairの何番めからの
    // データになるか。
    int ToukaAtomNum = 0; // そのペアで、いくつの等価な原子があるか。

    int PairAtom[MaxToukaAtomNum];
    // その等価な原子の番号を入れるベクトル。

    double X,Y,Z;//,vecXd1,vecYd1,vecZd1,vecXi1,vecYi1,vecZi1;
    double x,y,z,Norm,hantei;
    double x0,y0,z0,x1,y1,z1;// ,L,M,N;
    //int MemgetFlag1,MemgetFlag2;
    // int i;
    //int w;//ss, j,k,,s ,t ,CntNum; // intCenterAtomNum;
    //int pre,post,Agflag;
    int CntAtomNum,flag,scount,thetaflg,row,col; // ,m
    double SinTheta,HousenX,HousenY,HousenZ; // vecNormD1,vecNormI1, CosTheta,
    double HousenNorm,SinTheta1,SinTheta2,CosTheta1,CosTheta2;
    double MaxAngleGosa;
    double TotalAngleGosa,AverageAngleGosa,kyori,CosAngle,SinAngle;

    double eachAngle[MaxCntNum];
    double henkanval;
    double radX,radY,radZ;
    double Vdx1,Vdy1,Vdz1,Vdx2,Vdy2,Vdz2,VdNorm1,VdNorm2;
    double Vix1,Viy1,Viz1,Vix2,Viy2,Viz2,ViNorm1,ViNorm2;
    double Hdx,Hdy,Hdz,HdNorm,Hix,Hiy,Hiz,HiNorm;
    double HousenX3,HousenY3,HousenZ3;
    double dbgX,dbgY,dbgZ,dbgNorm;

    double dumymat[AfinMatrixYouso];
    double TransMatrix[AfinMatrixYouso];
    double Matrix[AfinMatrixYouso];
    double HozonMatrix[MaxMawasi][AfinMatrixYouso];
    double EachStepX[MaxMawasi][MaxCntNum];
    double EachStepY[MaxMawasi][MaxCntNum];
    double EachStepZ[MaxMawasi][MaxCntNum];
    double HozonXcore[MaxCntNum];
    double HozonYcore[MaxCntNum];
    double HozonZcore[MaxCntNum];
    double fitness[MaxMawasi];
    double U0[MaxUVN];
    double OptAngle[MaxUVN*MaxMawasi];

    double  minFitness;
    int     minFitnessNum;

    //--------------- ################### --------------------//
    //--------------- # Threshold Value # --------------------//
    //--------------- ################### --------------------//
    //キーワードからローカル変数への代入。
    const double Threshold1     = FirstAveGosa ;
    const double EachThreshold1 = FirstEachGosa;
    const double Threshold2     = LastAveGosa;
    const double EachThreshold2 = LastEachGosa;

    // スレッシュオールドの表示
    //   if(stepnum==0){
    //     if( OutLevel < -3 ){
    //       //----------------------------------------------------------//
    //       log<<"\n";
    //       log<<"     ########## Threshold Table ########"<<"\n";
    //       log<<"        Threshold 1 (Average) = "<< Threshold1     <<"\n";
    //       log<<"        EachThreshold 1       = "<< EachThreshold1 <<"\n";
    //       log<<                                                      "\n";
    //       log<<"        Threshold 2 (Average) = "<< Threshold2     <<"\n";
    //       log<<"        EachThreshold 2       = "<< EachThreshold2 <<"\n";
    //       log<<"    ------------------------------------"<<"\n";
    //       //----------------------------------------------------------//
    //     }
    //   }

    // 未知変数の数をローカル変数へ代入。
    const int UVN = MaxUVN; // Unknown Variable Number;
    // Unknown Variable Number is three.(angleX,angleY,angleZ);


    Fl_Geometry FlGeom(Fl_Geometry::getDefaultFileName());
//   FlGeom.load();
//   FlGeom.open("fl_Geometry", "read");
//   FlGeom.read();
//   FlGeom.close();

    //   if( OutLevel < -4 ){
    //     log<<"\n";
    //     log<<"   +++++++++++++++++++++++++++++++++"<<"\n";
    //     log<<"     STEPnumber(now) ["<<stepnum<<"]"<< "\n";
    //     log<<"   +++++++++++++++++++++++++++++++++"<<"\n";
    //     log<<"\n";
    //   }

    const int tAN = InpAmino.Res1AtomNum;
    //   if( OutLevel < -4 ){
    //     log<<"tAN = "<<tAN<<"\n";
    //   }

    // 等価な原子ペアがいくつあるかを調べる。
    int ToukaPairNum = DbAmino.equalPairNum;
    //   if( OutLevel < -4 ){
    //     log<<" Total ToukaPairNum = "<< ToukaPairNum <<"\n";
    //   }

    //#######################################################################
    //#######################################################################
    // Istart
    for (int i=0; i < tAN; ++i) {
        int MemgetFlag1 = 0;
        int MemgetFlag2 = 0 ;

        int CenterAtomNum = DbAmino.Order[i] ;
        //     if( OutLevel < -5 ){
        //       log<<"\n";
        //       log<<"======================================="<<"\n";
        //       log<<"          AtomNumber = ["<<CenterAtomNum<<"]"<<"\n";
        //       log<<"======================================="<<"\n"<<"\n";
        //       log.flush();
        //     }

        int CntNum = DbAmino.Atomdata[CenterAtomNum].ConectNum;
        // CntNum : Connect Number; not Center Number;

        if (CntNum == 1) {
            CntAtomNum = DbAmino.Atomdata[CenterAtomNum].Conection[0] ;
            //CntAtomNum : Connect Atom Number;

            //       dainyu_afinmat(DbAmino.Atomdata[CenterAtomNum].afinmat,
            //             DbAmino.Atomdata[CntAtomNum].afinmat);
            DbAmino.Atomdata[CenterAtomNum].afinmat = DbAmino.Atomdata[CntAtomNum].afinmat;

            InpAmino.Order[CenterAtomNum] = InpAmino.Order[CntAtomNum] ;

            //       if( OutLevel < -6 ){
            //    log<<" ConectAtom ga 1 ko ;"<<"\n";
            //    log<<" sorede OK!! InpAmino.Order["<<CenterAtomNum
            //       <<"] = InpAmino.Order[ "
            //       <<CntAtomNum<<" ] = "<<InpAmino.Order[CenterAtomNum]<<"\n";
            //    log.flush();
            //       }

            //goto pre_iLoopEnd;
            continue;
        }

        MemgetFlag1 = 1 ;

        //データベースに関する入力。
        double CenterAtomX = DbAmino.Atomdata[CenterAtomNum].x;
        double CenterAtomY = DbAmino.Atomdata[CenterAtomNum].y;
        double CenterAtomZ = DbAmino.Atomdata[CenterAtomNum].z;
        //     double L = CenterAtomX;
        //     double M = CenterAtomY;
        //     double N = CenterAtomZ;

        for (int t = 0; t < CntNum; t++) {
            const int s = DbAmino.Atomdata[CenterAtomNum].Conection[t] ;
            ConectAtomX[t] = DbAmino.Atomdata[s].x ;
            ConectAtomY[t] = DbAmino.Atomdata[s].y ;
            ConectAtomZ[t] = DbAmino.Atomdata[s].z ;
        }

        //入力データに関する入力。
        double InpCtrAtomX = InpAmino.Atomdata[CenterAtomNum].x ;
        double InpCtrAtomY = InpAmino.Atomdata[CenterAtomNum].y ;
        double InpCtrAtomZ = InpAmino.Atomdata[CenterAtomNum].z ;
        // Ctr : CenterAtom no ryaku;
        for (int t=0; t<CntNum; t++) {
            const int s = DbAmino.Atomdata[CenterAtomNum].Conection[t] ;
            InpAtomX[t] = InpAmino.Atomdata[s].x ;
            InpAtomY[t] = InpAmino.Atomdata[s].y ;
            InpAtomZ[t] = InpAmino.Atomdata[s].z ;
        }

        //-------------------
        // DEBUG PRINT
        //-------------------
        //     if( OutLevel < -9 ){
        //       log<<"L = "<<L<<"\n";
        //       log<<"M = "<<M<<"\n";
        //       log<<"N = "<<N<<"\n";
        //       log<<"debug  CenterInpAtom XYZ = "<<InpCtrAtomX<<"  "
        //   <<InpCtrAtomY<<"  "<<InpCtrAtomZ<<"\n";
        //       log.flush();
        //       for(t=0;t<CntNum;t++){
        //  log<<"debug  InpAtomXYZ = "<<InpAtomX[t]<<"   "
        //     <<InpAtomY[t]<<"   "<<InpAtomZ[t]<<"\n";
        //  log.flush();
        //       }
        //       log<<"debug  CenterDB_Atom XYZ = "<<L<<"  "<<M<<"  "<<N<<"\n";
        //       for(t=0;t<CntNum;t++){
        //  log<<"debug  ConectAtomXYZ = "<<ConectAtomX[t]<<"   "
        //     <<ConectAtomY[t]<<"   "<<ConectAtomZ[t]<<"\n";
        //  log.flush();
        //       }
        //     }

        // ------------- 考慮するCenter原子を原点に平行移動させる
        // ------------- 同時に、結合している原子を平行移動させる
        for (int t = 0; t < CntNum; t++) {
            ConectAtomX[t] -= CenterAtomX;
            ConectAtomY[t] -= CenterAtomY;
            ConectAtomZ[t] -= CenterAtomZ;
            InpAtomX[t]    -= InpCtrAtomX;
            InpAtomY[t]    -= InpCtrAtomY;
            InpAtomZ[t]    -= InpCtrAtomZ;
        }

        CenterAtomX = CenterAtomY = CenterAtomZ = 0.00 ;
        InpCtrAtomX = InpCtrAtomY = InpCtrAtomZ = 0.00 ;

        for (int t = 0; t < CntNum; t++) {
            HozonXcore[t] = ConectAtomX[t];
            HozonYcore[t] = ConectAtomY[t];
            HozonZcore[t] = ConectAtomZ[t];
        }

        //-------------------
        // DEBUG PRINT
        //-------------------
        //     if( OutLevel < -9 ){
        //       log<<"After Heikou Idou"<<"\n";
        //       log<<"debug  CenterInpAtom XYZ = "<<InpCtrAtomX<<"  "
        //   <<InpCtrAtomY<<"  "<<InpCtrAtomZ<<"\n";
        //       log.flush();
        //       for(t=0;t<CntNum;t++){
        //  log<<"debug  InpAtomXYZ = "<<InpAtomX[t]<<"   "
        //     <<InpAtomY[t]<<"   "<<InpAtomZ[t]<<"\n";
        //  log.flush();
        //       }
        //       log<<"debug  CenterDB_Atom XYZ = "<<CenterAtomX<<"  "
        //   <<CenterAtomY<<"  "<<CenterAtomZ<<"\n";

        //       log<<"First ConnectAtom( after Heikou idou )"<<"\n";
        //       log.flush();
        //       for(t=0;t<CntNum;t++){
        //  log<<"debug  ConectAtomXYZ = "<<ConectAtomX[t]<<"   "
        //     <<ConectAtomY[t]<<"   "<<ConectAtomZ[t]<<"\n";
        //  log.flush();
        //       }
        //     }

        //----------------------------------------------------------------//
        // 考慮している中心原子にたいする結合原子の中に等価な原子が入って //
        // いないかどうか調べる。                                         //
        //----------------------------------------------------------------//
        if (ToukaPairNum == 0) {
            mawasi = 1 ;
            flag=0;
        } else {
            flag=0;
            // 現在考慮している中心原子に対して、等価ペアがあるか／ないか。
            // flag==0 : 等価原子がない。
            // flag==1 : 等価原子有り。
            scount=0;
            bool bBreak_loop =false;
            for (int j = 0; j < ToukaPairNum; j++) {
                //log<<"debug loop_k is ["<<DbAmino.equalPaircount[j]<<"]";
                //log<<"\n";
                for (int k = 0; k < DbAmino.equalPaircount[j]; k++) {
                    int ss = DbAmino.equalPair[scount + k];
                    //log<<"debug  ss = "<< ss <<"\n";

                    for (int m = 0; m < CntNum; m++) {
                        int s = DbAmino.Atomdata[CenterAtomNum].Conection[m] ;
                        //log<<"debug s= "<< s <<"\n";

                        if (ss == s) {
                            PairBangou = k;
                            //log<<"debug PairBangou = "<<PairBangou<<"\n";
                            VectBangou = scount;
                            //log<<"debug VectBangou = "<<VectBangou<<"\n";
                            flag=1;

                            //goto skip;
                            bBreak_loop = true;
                            break;
                        }
                    }

                    if (bBreak_loop == true) {
                        break;
                    }
                }

                if (bBreak_loop == true) {
                    break;
                }

                scount += DbAmino.equalPaircount[j];
            }
        }

        //skip:;

        //log<<"debug   skip no sita;       scount = "<<scount<<"\n";

        //------------------------------------------------------------------//
        //等価原子が存在することも考慮したときに、何回ループを回すかの算出。//
        //------------------------------------------------------------------//
        if (flag == 0) {
            //       if( OutLevel < -5 ){
            //    log<<"Touka Atom  Nasi."<<"\n";
            //       }
            //等価原子なし。
            mawasi = 1 ;
        } else if (flag==1) {
            //       if( OutLevel < -5 ){
            //    log<<"Touka Atom  Ari."<<"\n";
            //       }
            log.flush();
            //等価原子あり。
            mawasi = 1 ;

            ToukaAtomNum = DbAmino.equalPaircount[PairBangou];
            //       if( OutLevel < -5 ){
            //    log<<"sono atom no ToukaAtomNum = "<< ToukaAtomNum<<"\n";
            //    log<<"PairBangou "<<PairBangou<<"\n";
            //    log.flush();
            //       }
            // mawasi = (ToukaAtomNum)! no keisan;
            for (int j = 1; j <= ToukaAtomNum; j++) {
                mawasi *= j;
            }

            // PairAtomにはDB分子中の通し番号が入る。
            for (int j=0; j < ToukaAtomNum; j++) {
                PairAtom[j] = DbAmino.equalPair[scount + j];
                //  if( OutLevel < -5 ){
                //    log<<"  PairAtom["<<j<<"] = "<<PairAtom[j]<<"\n";
                //  }
                log.flush();
            }

            // ConectSerialAtomNumには等価原子の並びがConectAtom{X,Y,Z}の
            // 何番目にあたるかを保存。
            for (int j = 0; j < ToukaAtomNum; j++) {
                for (int t=0; t<CntNum; t++) {
                    int s = DbAmino.Atomdata[CenterAtomNum].Conection[t];
                    if (PairAtom[j] == s) {
                        ConectSerialAtomNum[j] = t ;
                    }
                }
            }

            // debug
            //       if( OutLevel < -5 ){
            //    for(j=0;j<ToukaAtomNum;j++){
            //      log<<" ConectSerialAtomNum["<<j<<"] = "
            //         <<ConectSerialAtomNum[j]<<"\n";
            //      log.flush();
            //    }
            //       }

            MemgetFlag2 = 1 ;
        }


        //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        //$
        // ############### Angle Optimyzation Start ########################
        //$
        //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        //##################################################################//
        // 中心原子に結合している原子の座標のずれを最小にする               //
        // 回転角(θx,θy,θz,)を捜し出す。                                 //
        //##################################################################//
        //Jstart
        //     if( OutLevel < -5 ){
        //       log<<"     MAWASI = "<<mawasi<<"\n";
        //       log<<"\n"<<"JSTART"<<"\n"<<"\n";
        //       log.flush();
        //     }

        for (int j=0; j < mawasi; j++) {
            this->dainyu_vector(EachStepX[j],ConectAtomX,CntNum);
            this->dainyu_vector(EachStepY[j],ConectAtomY,CntNum);
            this->dainyu_vector(EachStepZ[j],ConectAtomZ,CntNum);

            //       if( OutLevel < -7 ){
            //    log<<"J no 1 ban  sentou *****************************"<<"\n";
            //    log.flush();
            //    for(t=0;t<CntNum;t++){
            //      log<<"debug  ConectAtomXYZ = "<<ConectAtomX[t]<<"   "
            //         <<ConectAtomY[t]<<"   "<<ConectAtomZ[t]<<"\n";
            //      log.flush();
            //    }
            //       }

            //-------------------------------------------------------------
            // 結合している原子中の最初の２つの原子のうち、
            // １つめの原子と中心原子が作るベクトルを入力構造のベクトルと
            // 一致させ、さらに２つめの原子との３原子で作る平面を一致させる
            // ように全体の構造を回転させる。
            //-------------------------------------------------------------

            double vecXd1 = ConectAtomX[0];
            double vecXi1 = InpAtomX[0];
            double vecYd1 = ConectAtomY[0];
            double vecYi1 = InpAtomY[0];
            double vecZd1 = ConectAtomZ[0];
            double vecZi1 = InpAtomZ[0];

            double vecNormD1= sqrt(pow(vecXd1,2) + pow(vecYd1,2) + pow(vecZd1,2)) ;
            double vecNormI1= sqrt(pow(vecXi1,2) + pow(vecYi1,2) + pow(vecZi1,2)) ;
            if ((vecNormD1 < 1E-8) || (vecNormI1 < 1E-8)) {
                log<<"**** Error vecNorm is smaller than 1E-8"<<"\n";
                log<<"**** vecNnormD1 = "<<vecNormD1<<"\n";
                log<<"**** vecNnormI1 = "<<vecNormI1<<"\n";
                //                log<<"\7"<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   norm "<<"\n";
                log.flush();
                CnErr.abort();
            }

            vecXd1 /= vecNormD1;
            vecXi1 /= vecNormI1;
            vecYd1 /= vecNormD1;
            vecYi1 /= vecNormI1;
            vecZd1 /= vecNormD1;
            vecZi1 /= vecNormI1;

            double CosTheta = vecXd1*vecXi1 + vecYd1*vecYi1 +vecZd1*vecZi1;
            //       if( OutLevel < -5 ){
            //    log<<"CosTheta = "<<CosTheta<<"\n";
            //    log.flush();
            //       }

            //平行移動した時点で１つ目のベクトル方向が一致している場合。
            if (pow(fabs(CosTheta-1.0) , 2) < 1E-13) {
                log<<" Vecto awase : type CosTheta==1.0"<<"\n";
                log<<"CosTheta = "<<CosTheta<<"\n";
                initialyze_AfinMatrix(DbAmino.Atomdata[CenterAtomNum].afinmat);
                initialyze_AfinMatrix(TransMatrix);

                for (int t=0; t < CntNum; t++) {
                    dmX1[t] = ConectAtomX[t];
                    dmY1[t] = ConectAtomY[t];
                    dmZ1[t] = ConectAtomZ[t];
                }
                SinTheta=0.0;
                //  if( OutLevel < -5 ){
                //    log<<" type : CosTheta==1 "<<"\n";
                //    log.flush();
                //  }
                goto second_rotation;
            }

            //平行移動した時点で１つ目のベクトル方向が逆向きになっている場合。
            if (pow(fabs(CosTheta+1.0) , 2) < 1E-13) {
                log<<" Vecto awase : type CosTheta==-1.0"<<"\n";
                log<<"CosTheta = "<<CosTheta<<"\n";
                if (pow(fabs(vecXi1),2) >= pow(fabs(vecYi1),2)) {
                    if (pow(fabs(vecXi1),2) >= pow(fabs(vecZi1),2)) {
                        if (pow(fabs(vecYi1),2) >= pow(fabs(vecZi1),2)) {
                            // x >= y >= z no jyun.
                            HousenX3 = -sqrt(2.0);
                            HousenY3 = sqrt(3.0);
                            HousenZ3 = sqrt(5.0);
                        } else {
                            // x >= z > y no jyun.
                            HousenZ3 = -sqrt(2.0);
                            HousenX3 = sqrt(3.0);
                            HousenY3 = sqrt(5.0);
                        }
                    } else {
                        // z > x > y no jyun.
                        HousenZ3 = -sqrt(2.0);
                        HousenX3 = sqrt(3.0);
                        HousenY3 = sqrt(5.0);
                    }
                } else {
                    if (pow(fabs(vecXi1),2) < pow(fabs(vecZi1),2)) {
                        // z , y > x no jyun.
                        HousenZ3 = -sqrt(2.0);
                        HousenY3 = sqrt(3.0);
                        HousenX3 = sqrt(5.0);
                    } else {
                        //  y > x > z no jyun.
                        HousenY3 = -sqrt(2.0);
                        HousenX3 = sqrt(3.0);
                        HousenZ3 = sqrt(5.0);
                    }
                }

                //---- 仮の回転軸となる直線の法線ベクトルを求める。
                HousenX = vecYi1*HousenZ3 - vecZi1*HousenY3 ;
                HousenY = vecZi1*HousenZ3 - vecXi1*HousenZ3 ;
                HousenZ = vecXi1*HousenY3 - vecYi1*HousenZ3 ;

                SinTheta=0.0;
                //  if( OutLevel < -5 ){
                //    log<<" type : CosTheta==-1 "<<"\n";
                //    log.flush();
                //  }

                goto taiou;
            }// type CosTheta==-1;

            //平行移動した時点で１つ目のベクトル方向が一致していない場合。
            if ((CosTheta < -1) || (CosTheta > 1)) {
                //              log<<"Bad value :   CosTheta = ["<<CosTheta<<"]"<<"\7"<<"\n";
                log<<"Bad value :   CosTheta = ["<<CosTheta<<"]"<<"  "<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   costheta"<<"\n";
                log.flush();
                CnErr.abort();
            }
            SinTheta = sqrt(1 - pow(CosTheta,2));

            //---- 回転軸となる直線の法線ベクトルを求める。
            HousenX = vecYi1*vecZd1 - vecZi1*vecYd1 ;
            HousenY = vecZi1*vecXd1 - vecXi1*vecZd1 ;
            HousenZ = vecXi1*vecYd1 - vecYi1*vecXd1 ;

            ////////////
            //
taiou:;      //
            //
            ///////////

            HousenNorm = sqrt(pow(HousenX,2)+pow(HousenY,2)+pow(HousenZ,2));
            if (HousenNorm < 1E-7) {
                log<<"### Cannot define Housen_Vector###"<<"\n";
                log<<"HousenNorm = "<<HousenNorm<<"\n";
                //      log<<"\7"<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   Housen"<<"\n";
                log.flush();
                CnErr.abort();
            }
            HousenX /= HousenNorm;
            HousenY /= HousenNorm;
            HousenZ /= HousenNorm;

            //log <<"HousenX = "<<HousenX<<"\n";
            //log <<"HousenY = "<<HousenY<<"\n";
            //log <<"HousenZ = "<<HousenZ<<"\n";
            //log.flush();

            // 角の算出。
            SinTheta2 = sqrt(pow(HousenX,2) + pow(HousenZ,2)) ;
            SinTheta1 = HousenX / SinTheta2;
            CosTheta1 = HousenZ / SinTheta2;
            CosTheta2 = HousenY;
            thetaflg=0;

            //       if( OutLevel < -9 ){
            //    log <<"SinTheta1 = "<<SinTheta1<<"\n";
            //    log <<"SinTheta2 = "<<SinTheta2<<"\n";
            //    log <<"CosTheta1 = "<<CosTheta1<<"\n";
            //    log <<"CosTheta2 = "<<CosTheta2<<"\n";
            //    log.flush();
            //       }

            //###############//
            //
first_rotation1:;   //
            //
            //###############//

            if (thetaflg == 1) {
                SinTheta = -SinTheta;
            }

            initialyze_AfinMatrix(TransMatrix);
            //--------------------------------------------------------------
            //-- 回転の任意軸を、Y軸を回転軸としてYZ平面まで回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=0;
            col=0;
            dumymat[ row*4 + col ] =  CosTheta1;
            row=0;
            col=2;
            dumymat[ row*4 + col ] =  SinTheta1;
            row=2;
            col=0;
            dumymat[ row*4 + col ] = -SinTheta1;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta1;

            multiply_afinmat(TransMatrix,dumymat,Matrix);

            dainyu_afinmat(TransMatrix,Matrix);

            //--------------------------------------------------------------
            //-- 移動した任意軸を、X軸を回転軸としてY軸まで回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=1;
            col=1;
            dumymat[ row*4 + col ] =  CosTheta2;
            row=1;
            col=2;
            dumymat[ row*4 + col ] = -SinTheta2;
            row=2;
            col=1;
            dumymat[ row*4 + col ] =  SinTheta2;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta2;

            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);
            //--------------------------------------------------------------
            //-- 希望の回転角θだけ回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=0;
            col=0;
            dumymat[ row*4 + col ] =  CosTheta;
            row=0;
            col=2;
            dumymat[ row*4 + col ] =  SinTheta;
            row=2;
            col=0;
            dumymat[ row*4 + col ] = -SinTheta;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta;
            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);

            //--------------------------------------------------------------
            //-- 回転した任意軸を、X軸を回転軸として角度（-θ2）だけ回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=1;
            col=1;
            dumymat[ row*4 + col ] =  CosTheta2;
            row=1;
            col=2;
            dumymat[ row*4 + col ] =  SinTheta2;
            row=2;
            col=1;
            dumymat[ row*4 + col ] = -SinTheta2;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta2;

            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);
            //--------------------------------------------------------------
            //-- 回転した任意軸を、Y軸を回転軸として角度（-θ1）だけ回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=0;
            col=0;
            dumymat[ row*4 + col ] =  CosTheta1;
            row=0;
            col=2;
            dumymat[ row*4 + col ] = -SinTheta1;
            row=2;
            col=0;
            dumymat[ row*4 + col ] =  SinTheta1;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta1;

            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);
            //--------------------------------------------------------------

            //--------------------------------------------------------------
            // 一つめのベクトル(vecXd1,vecYd1,vecZd1)が一致したかどうかを
            // 確かめるために回転させてみる。
            // x=ConectAtomX[0]; y=ConectAtomY[0]; z=ConectAtomZ[0];
            // [X Y Z T] = [x y z t]*AFINMATRIX ;  // T,t is dumy.
            //--------------------------------------------------------------
            x=ConectAtomX[0];
            y=ConectAtomY[0];
            z=ConectAtomZ[0];
            X = x*TransMatrix[0] + y*TransMatrix[4] + z*TransMatrix[8];
            Y = x*TransMatrix[1] + y*TransMatrix[5] + z*TransMatrix[9];
            Z = x*TransMatrix[2] + y*TransMatrix[6] + z*TransMatrix[10];

            Norm = sqrt(pow(X,2) + pow(Y,2) + pow(Z,2));
            X /= Norm;
            Y /= Norm;
            Z /= Norm;
            hantei=sqrt(pow((vecXi1-X),2)+pow((vecYi1-Y),2)+pow((vecZi1-Z),2));
            if (hantei < 1E-10) {
                //---------------------------//
                //結合原子を全て回転させる   //
                //---------------------------//
                for (int t = 0; t < CntNum; t++) {
                    dmX1[t] =   ConectAtomX[t] * TransMatrix[0]
                                + ConectAtomY[t] * TransMatrix[4]
                                + ConectAtomZ[t] * TransMatrix[8] ;

                    dmY1[t] =   ConectAtomX[t] * TransMatrix[1]
                                + ConectAtomY[t] * TransMatrix[5]
                                + ConectAtomZ[t] * TransMatrix[9] ;

                    dmZ1[t] =   ConectAtomX[t] * TransMatrix[2]
                                + ConectAtomY[t] * TransMatrix[6]
                                + ConectAtomZ[t] * TransMatrix[10] ;
                }
                dainyu_afinmat(DbAmino.Atomdata[CenterAtomNum].afinmat ,
                               TransMatrix);
                goto second_rotation;
            } else {
                thetaflg++;
                if (thetaflg==2) {
                    log<<"#@#@#@ Not agree : Bad angular rotation ";
                    //            log<<"first_rotaion1"<< "\7" << "\n";
                    log<<"first_rotaion1"<< "  " << "\n";
                    log<<"in DfInitialguess::makeAfinmatrix  thetaflg"<<"\n";
                    log.flush();
                    CnErr.abort();
                }
                goto first_rotation1;
            }

            //-------------------------------------------------------
            // ここで、１つめのベクトルが一致したことになる。
            //-------------------------------------------------------

            //##############//
            //
second_rotation:;  //
            //
            //##############//

            //       if( OutLevel < -7 ){
            //    log<<" 1 tu me no vector ga itti sita toki no afinmat"<<"\n";
            //    display_afinmatrix(TransMatrix);
            //    log<<"\n";
            //    log<<"1tu me no Vector ga itti sita tyokugo. "<<"\n";
            //    for(t=0;t<CntNum;t++){
            //      log<<"debug   dmXYZ1[ "<<t<<" ] = "<<dmX1[t]<<"  "
            //         <<dmY1[t]<<"  "<<dmZ1[t]<<"\n";
            //    }
            //       }
            /*
            //-------------------------------------------------------
            // 次に平面を合わせるように２つ目のベクトルを合わせる。
            //-------------------------------------------------------
            //---------------------------------------------
            // 今度の回転軸は、一致したベクトル、つまり
            // (vecXi1,vecYi1,vecZi1)である。
            //平面を合わせるために、ぞれぞれの平面の法線ベクトルを求めて
            //そのベクトルの間の角度分だけ、上の軸を中心に回転させる。
            //----------------------------------------------------------
            //データベースの原子から作られる平面の法線ベクトルを求める。
            //----------------------------------------------------------
            */
            Vdx1 = dmX1[0];
            Vdx2 = dmX1[1];
            Vdy1 = dmY1[0];
            Vdy2 = dmY1[1];
            Vdz1 = dmZ1[0];
            Vdz2 = dmZ1[1];
            VdNorm1 = sqrt(pow(Vdx1,2.0) + pow(Vdy1,2.0) + pow(Vdz1,2.0));
            VdNorm2 = sqrt(pow(Vdx2,2.0) + pow(Vdy2,2.0) + pow(Vdz2,2.0));
            if ((VdNorm1 < 1E-8) || (VdNorm2 < 1E-8)) {
                log<<"**** Error vecNorm is smaller than 1E-8"<<"\n";
                log<<"**** VdNorm1 = "<<VdNorm1<<"\n";
                log<<"**** VdNorm2 = "<<VdNorm2<<"\n";
                //                log<<"\7"<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   norm "<<"\n";
                CnErr.abort();
            }

            Vdx1 /= VdNorm1 ;
            Vdx2 /= VdNorm2 ;
            Vdy1 /= VdNorm1 ;
            Vdy2 /= VdNorm2 ;
            Vdz1 /= VdNorm1 ;
            Vdz2 /= VdNorm2 ;

            Hdx = Vdy2*Vdz1 - Vdz2*Vdy1 ;
            Hdy = Vdz2*Vdx1 - Vdx2*Vdz1 ;
            Hdz = Vdx2*Vdy1 - Vdy2*Vdx1 ;

            HdNorm = sqrt(pow(Hdx,2.0) + pow(Hdy,2.0) + pow(Hdz,2.0));
            if (HdNorm < 1E-7) {
                log<<"### Cannot define Housen_Vector###"<<"\n";
                log<<"HdNorm = "<<HdNorm<<"\n";
                //      log<<"\7"<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   Housen"<<"\n";
                CnErr.abort();
            }
            Hdx /= HdNorm;
            Hdy /= HdNorm;
            Hdz /= HdNorm;

            //--------------------------------------------------------
            //入力分子中の原子から作られる平面の法線ベクトルを求める。
            //--------------------------------------------------------
            Vix1 = InpAtomX[0];
            Vix2 = InpAtomX[1];
            Viy1 = InpAtomY[0];
            Viy2 = InpAtomY[1];
            Viz1 = InpAtomZ[0];
            Viz2 = InpAtomZ[1];
            ViNorm1 = sqrt(pow(Vix1,2.0) + pow(Viy1,2.0) + pow(Viz1,2.0));
            ViNorm2 = sqrt(pow(Vix2,2.0) + pow(Viy2,2.0) + pow(Viz2,2.0));
            if ((ViNorm1 < 1E-8) || (ViNorm2 < 1E-8)) {
                log<<"**** Error vecNorm is smaller than 1E-8"<<"\n";
                log<<"**** ViNorm1 = "<<ViNorm1<<"\n";
                log<<"**** ViNorm2 = "<<ViNorm2<<"\n";
                //                log<<"\7"<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   norm "<<"\n";
                CnErr.abort();
            }
            Vix1 /= ViNorm1 ;
            Vix2 /= ViNorm2 ;
            Viy1 /= ViNorm1 ;
            Viy2 /= ViNorm2 ;
            Viz1 /= ViNorm1 ;
            Viz2 /= ViNorm2 ;

            Hix = Viy2*Viz1 - Viz2*Viy1 ;
            Hiy = Viz2*Vix1 - Vix2*Viz1 ;
            Hiz = Vix2*Viy1 - Viy2*Vix1 ;

            HiNorm = sqrt(pow(Hix,2.0) + pow(Hiy,2.0) + pow(Hiz,2.0));
            if (HiNorm < 1E-7) {
                log<<"### Cannot define Housen_Vector###"<<"\n";
                log<<"HiNorm = "<<HiNorm<<"\n";
                //      log<<"\7"<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   Housen"<<"\n";
                CnErr.abort();
            }
            Hix /= HiNorm;
            Hiy /= HiNorm;
            Hiz /= HiNorm;

            //---------------------------------------------------------
            //次に、上で求めたHd{x,y,z},Hi{x,y,z}の２つのベクトルの間の
            //角度を求める。
            //---------------------------------------------------------
            CosTheta = Hdx*Hix + Hdy*Hiy + Hdz*Hiz ;
            //       if( OutLevel < -5 ){
            //    log<<"CosTheta = "<<CosTheta<<"\n";
            //       }
            //１つ目のベクトル方向を合わせた時点で２つ目のベクトルが
            //入力原子の２つ目のベクトルと同一平面上（重なっている時）にあるとき。
            if (pow((CosTheta-1.0),2)<1E-13) {
                dainyu_afinmat(HozonMatrix[j],
                               DbAmino.Atomdata[CenterAtomNum].afinmat);
                for (int t=0; t < CntNum; t++) {
                    dmX2[t] = dmX1[t];
                    dmY2[t] = dmY1[t];
                    dmZ2[t] = dmZ1[t];
                }
                log<<" Heimen awase : type CosTheta==1.0"<<"\n";
                goto next_step1;
            }
            //１つ目のベクトル方向を合わせた時点で２つ目のベクトルが
            //入力原子の２つ目のベクトルと同一平面上（反対の平面）にあるとき。
            if (pow(fabs(CosTheta+1.0) , 2)<1E-13) {
                log<<" Heimen awase : type CosTheta==-1.0"<<"\n";
                //              dbgX = Vdx1;
                //              dbgY = Vdy1;
                //              dbgZ = Vdz1;
                dbgX = Vix1;
                dbgY = Viy1;
                dbgZ = Viz1;
                goto taiou2;
            }
            //１つ目のベクトル方向を合わせた時点で２つ目のベクトルが
            //入力原子の２つ目のベクトルと同一平面上にないとき。
            if ((CosTheta < -1) || (CosTheta > 1)) {
                //              log<<"Bad value :   CosTheta = ["<<CosTheta<<"]"<<"\7"<<"\n";
                log<<"Bad value :   CosTheta = ["<<CosTheta<<"]"<<"  "<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   costheta"<<"\n";
                CnErr.abort();
            }
            SinTheta = sqrt(1.0 - pow(CosTheta,2.0));
            if ((SinTheta < -1) || (SinTheta > 1)) {
                //              log<<"Bad value :   SinTheta = ["<<SinTheta<<"]"<<"\7"<<"\n";
                log<<"Bad value :   SinTheta = ["<<SinTheta<<"]"<<"  "<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   sintheta"<<"\n";
                CnErr.abort();
            }
            CosAngle = acos(CosTheta)*180.0/3.14159265358979;
            SinAngle = asin(SinTheta)*180.0/3.14159265358979;

            dbgX = Hdy*Hiz - Hdz*Hiy;
            dbgY = Hdz*Hix - Hdx*Hiz;
            dbgZ = Hdx*Hiy - Hdy*Hix;

            ////////////
taiou2:;     //
            ////////////

            dbgNorm = sqrt(pow(dbgX,2.0) + pow(dbgY,2.0) + pow(dbgZ,2.0));
            if (dbgNorm < 1E-8) {
                //              log<<"Bad value : dbgNorm = ["<<dbgNorm<<"]"<<"\7"<<"\n";
                log<<"Bad value : dbgNorm = ["<<dbgNorm<<"]"<<"  "<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix "<<"\n";
                log<<"You must Rotation Geometry for Calculation's Molecule."<<"\n";
                CnErr.abort();
            }

            dbgX /= dbgNorm ;
            dbgY /= dbgNorm ;
            dbgZ /= dbgNorm ;
            HousenX = dbgX;
            HousenY = dbgY;
            HousenZ = dbgZ;
            //-----------------------------------
            // 各種、三角関数（角度）の計算
            //-----------------------------------
            SinTheta2 = sqrt(pow(HousenX,2.0) + pow(HousenZ,2.0)) ;
            if (SinTheta2 < 1E-8) {
                //              Out<<"Bad value :   SinTheta2 = ["<<SinTheta2<<"]"<<"\7"<<"\n";
                log<<"Bad value :   SinTheta2 = ["<<SinTheta2<<"]"<<"  "<<"\n";
                log<<" in DfInitialguess::makeAfinmatrix   sintheta"<<"\n";
                log<<"You must Rotation Geometry for Calculation's Molecule."<<"\n";
                CnErr.abort();
            }
            SinTheta1 = HousenX / SinTheta2;
            CosTheta1 = HousenZ / SinTheta2;
            CosTheta2 = HousenY;

            //平面一致の判定の基準の一つを作る。
            kyori = sqrt(pow(dmX1[1]-InpAtomX[1],2.0) +
                         pow(dmY1[1]-InpAtomY[1],2.0) +
                         pow(dmZ1[1]-InpAtomZ[1],2.0)) ;

            //       if( OutLevel < -9 ){
            //    log<<"\n"<<"debug  @@@ KYORI(hajime[1]) = "<<kyori<<"\n";
            //       }

            thetaflg=0;
            SinTheta = -SinTheta ; //こっちを先にした方が早いみたい。

first_rotation2:;  //

            if (thetaflg==1) {
                SinTheta = -SinTheta;
            }

            //log<<"debug SinTheta = "<<SinTheta<<"\n";

            initialyze_AfinMatrix(TransMatrix);
            //--------------------------------------------------------------
            //-- 回転の任意軸を、Y軸を回転軸としてYZ平面まで回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=0;
            col=0;
            dumymat[ row*4 + col ] =  CosTheta1;
            row=0;
            col=2;
            dumymat[ row*4 + col ] =  SinTheta1;
            row=2;
            col=0;
            dumymat[ row*4 + col ] = -SinTheta1;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta1;

            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);
            //--------------------------------------------------------------
            //-- 移動した任意軸を、X軸を回転軸としてY軸まで回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=1;
            col=1;
            dumymat[ row*4 + col ] =  CosTheta2;
            row=1;
            col=2;
            dumymat[ row*4 + col ] = -SinTheta2;
            row=2;
            col=1;
            dumymat[ row*4 + col ] =  SinTheta2;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta2;

            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);
            //--------------------------------------------------------------
            //-- 希望の回転角θだけ回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=0;
            col=0;
            dumymat[ row*4 + col ] =  CosTheta;
            row=0;
            col=2;
            dumymat[ row*4 + col ] =  SinTheta;
            row=2;
            col=0;
            dumymat[ row*4 + col ] = -SinTheta;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta;

            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);

            //--------------------------------------------------------------
            //-- 回転した任意軸を、X軸を回転軸として角度（-θ2）だけ回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=1;
            col=1;
            dumymat[ row*4 + col ] =  CosTheta2;
            row=1;
            col=2;
            dumymat[ row*4 + col ] =  SinTheta2;
            row=2;
            col=1;
            dumymat[ row*4 + col ] = -SinTheta2;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta2;

            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);
            //--------------------------------------------------------------
            //-- 回転した任意軸を、Y軸を回転軸として角度（-θ1）だけ回転する。
            //--------------------------------------------------------------
            initialyze_AfinMatrix(dumymat);
            clear_AfinMatrix(Matrix);

            row=0;
            col=0;
            dumymat[ row*4 + col ] =  CosTheta1;
            row=0;
            col=2;
            dumymat[ row*4 + col ] = -SinTheta1;
            row=2;
            col=0;
            dumymat[ row*4 + col ] =  SinTheta1;
            row=2;
            col=2;
            dumymat[ row*4 + col ] =  CosTheta1;

            multiply_afinmat(TransMatrix,dumymat,Matrix);
            dainyu_afinmat(TransMatrix,Matrix);
            //--------------------------------------------------------------

            //--------------------------------------------------------------
            // 2つめのベクトル(vecXd1,vecYd1,vecZd1)が同一平面上にのっている
            // ことを確かめるために回転させてみる。
            // x=dmX1[1]; y=dmY1[1]; z=dmZ1[1];
            // [X Y Z T] = [x y z t]*AFINMATRIX ;  // T,t is dumy.
            //--------------------------------------------------------------
            x=dmX1[1];
            y=dmY1[1];
            z=dmZ1[1];
            X = x*TransMatrix[0] + y*TransMatrix[4] + z*TransMatrix[8];
            Y = x*TransMatrix[1] + y*TransMatrix[5] + z*TransMatrix[9];
            Z = x*TransMatrix[2] + y*TransMatrix[6] + z*TransMatrix[10];

            //--------------------------------------------------------------
            // 回転を実行させて、平面にのったかどうかの確認
            // 原点を通る平面は aX+bY+cZ=0 なので（a,b,cは法線ベクトル）、
            // 上で求めた(X,Y,Z)を入れてゼロになるか調べる。
            //--------------------------------------------------------------
            hantei = Hix*X + Hiy*Y + Hiz*Z;
            hantei = fabs(hantei);
            if ((hantei < 1E-6) &&
                    (kyori > sqrt(pow(X-InpAtomX[1],2) +
                                  pow(Y-InpAtomY[1],2) +
                                  pow(Z-InpAtomZ[1],2)))) {

                //  if( OutLevel < -9 ){
                //    log<<"        ### ### secondVector Heimen ni notta"<<"\n";
                //  }
                //---------------------------//
                //結合原子を全て回転させる   //
                //---------------------------//
                for (int t = 0; t < CntNum; t++) {
                    dmX2[t] =   dmX1[t] * TransMatrix[0]
                                + dmY1[t] * TransMatrix[4]
                                + dmZ1[t] * TransMatrix[8] ;

                    dmY2[t] =   dmX1[t] * TransMatrix[1]
                                + dmY1[t] * TransMatrix[5]
                                + dmZ1[t] * TransMatrix[9] ;

                    dmZ2[t] =   dmX1[t] * TransMatrix[2]
                                + dmY1[t] * TransMatrix[6]
                                + dmZ1[t] * TransMatrix[10] ;
                }

                //--------------------------//
                //debug
                //--------------------------//
                //  if( OutLevel < -7 ){
                //    log<<" After Trans Tamesi ni look"<<"\n";
                //    for(t=0;t<CntNum;t++){
                //      log<<"dmXYZ2 = "<<dmX2[t]<<"  "<<dmY2[t]<<"  "
                //         <<dmZ2[t]<<"\n";
                //    }
                //  }
                //---------------------------------------------
                //平面に乗った時点で、dm{X,Y,Z}2に座標値が入っている。
                //---------------------------------------------
                clear_AfinMatrix(dumymat);
                multiply_afinmat(DbAmino.Atomdata[CenterAtomNum].afinmat,
                                 TransMatrix,dumymat);
                //  if( OutLevel < -7 ){
                //    log<<"second Vector no tameno matrix"<<"\n";
                //    display_afinmatrix(DbAmino.Atomdata[CenterAtomNum].afinmat);
                //  }

                dainyu_afinmat(HozonMatrix[j],dumymat);

                //  if( OutLevel < -7 ){
                //    log<<"2Heimen Itti no tyokugo no Matrix "<<"\n";
                //    log<<"mawasi no j = [ "<<j<<" ]"<<"\n";
                //    display_afinmatrix(HozonMatrix[j]);
                //  }

                goto next_step1;

            } else {
                thetaflg++;

                if (thetaflg == 2) {
                    log<<"Warning : #@#@#@ Not agree : Bad angular rotation ";
                    //            log<<"Warning : first_rotaion2"<< "\7" << "\n";
                    log<<"Warning : first_rotaion2"<< "  " << "\n";
                    log<<"Warning : DfInitialguess::makeAfinmatrix rot2"<<"\n";
                    dainyu_vector(dmX2,dmX1,CntNum);
                    dainyu_vector(dmY2,dmY1,CntNum);
                    dainyu_vector(dmZ2,dmZ1,CntNum);
                    goto serch; // 95/03/22 Changed by Okazaki and Fukue.
                    //CnErr.abort();
                }
                goto first_rotation2;
            }
            //-------------------------------------------------------
            // ここで、最初の２原子と中心原子が平面にのったことになる。
            //-------------------------------------------------------

            //##########//
            //
next_step1:    //
            //
            //##########//

            //#######################################################
            //ここから、やっと原子座標のずれを最小にする角を探す。  #
            //#######################################################
            //-----------------------------------------
            //その前に、この時点で座標のずれを調べて、
            //threshold以下であれば探索は行なわない。
            //
            //それぞれのInpAminoの原子とDbAminoの結合原子が作る角度
            //を算出する。
            //
            //誤差角度の平均を算出する。
            //-----------------------------------------
            TotalAngleGosa=0.0;
            MaxAngleGosa=1E-10;
            for (int t = 0; t < CntNum; t++) {
                x0 = dmX2[t];
                y0 = dmY2[t];
                z0 = dmZ2[t];
                x1 = InpAtomX[t];
                y1 = InpAtomY[t];
                z1 = InpAtomZ[t];
                // angle is degree.
                eachAngle[t] = calc_angle(x0,y0,z0,x1,y1,z1);
//  if( OutLevel < -5 ){
//    log<<"eachAngle[ "<<t<<" ] = "<<eachAngle[t]<<"         ";
//    log<<"distance = "<< sqrt( pow(x0-x1,2) + pow(y0-y1,2)
//                   +pow(z0-z1,2) )<<"\n";
//  }

                if (eachAngle[t] > MaxAngleGosa) {
                    MaxAngleGosa = eachAngle[t] ;
                }

                TotalAngleGosa += eachAngle[t] ;
            }
            AverageAngleGosa = TotalAngleGosa / CntNum;

//       if( OutLevel < -5 ){
//  log<<" ------------- First judge ------------- "<<"\n";
//  log<<"  Total  AngleGosa  = "<<TotalAngleGosa<<"\n";
//  log<<"  AverageAngleGosa  = "<<AverageAngleGosa<<"\n";
//  log<<"  MaxAngleGosa      = "<<MaxAngleGosa <<"\n";
//  log.flush();
//       }

            //------------------------------------------------------------
            //条件を満たすかの判定。
            // 条件は、誤差の平均がThreshold1以下であること、かつ、
            // それぞれの誤差がEachThreshold1以下であること、の２つ。
            //------------------------------------------------------------
            {
                int Agflag = 0;
                if (AverageAngleGosa <= Threshold1) {
                    for (int t = 0; t < CntNum; t++) {
                        Agflag = 0;

                        if (eachAngle[t] <= EachThreshold1) {
                            Agflag=1;
                        } else {
                            //break;
                            goto serch;
                        }
                    }

                    if (Agflag == 1) {
                        InpAmino.Order[CenterAtomNum] = 1;
                        //    if( OutLevel < -4 ){
                        //      log<<"\n";
                        //      log<<"########################################"<<"\n";
                        //      log<<"#                                      #"<<"\n";
                        //      log<<"     OK!! InpAmino.Order[";
                        //      log<<      CenterAtomNum<<"] = 1 ;     "       <<"\n";
                        //      log<<"#                                      #"<<"\n";
                        //      log<<"########################################"<<"\n";
                        //      log.flush();
                        //    }
                        goto warp_end;
                    }
                }
            }
serch:
            ;

            //       if( OutLevel < -5 ){
            //    log<<"\n"<<"\n";
            //    log<<"                  ###########################"<<"\n";
            //    log<<"                  #  Powell  Method  Start  #"<<"\n";
            //    log<<"                  ###########################"<<"\n";
            //    log.flush();
            //       }
            // Initialuze parameter U0; U0 is three parameter.
            // It is TheataX,TheataY,TheataZ : Unit is  [degree].
            for (int t = 0; t < UVN; t++) {
                U0[t] = 0.0 ;
            }

            // ##### Call   Modify-Powell-Mehtod Routine  #####
            {
                powell(U0,dmX2,dmY2,dmZ2,InpAtomX,InpAtomY,InpAtomZ,CntNum);
            }
            fitness[j] = AfinFunction(U0,dmX2,dmY2,dmZ2,
                                      InpAtomX,InpAtomY,InpAtomZ,CntNum);

            //       if( OutLevel < -5 ){
            //    log<<"\n";
            //    log<<"after Modify Powell Method "<<"\n";
            //    log<<"Optimal Fitness = "<<fitness[j]<<"\n";
            //    log<<"Optimal Angle = "<<U0[0]<<" ,  "
            //       <<U0[1]<<" ,  "<<U0[2]<<"\n";
            //    log.flush();
            //       }

            for (int w=0; w < UVN; w++) {
                OptAngle[ j*UVN + w ] = U0[w] ;
            }

            // ############### for Touka Atom #####################
            // SWAP ;ここで入れ替えを行なう。
            if (j != (mawasi-1)) {
                //  if( OutLevel < -5 ){
                //    log<<"---- NARABIKAE   START ------"<<"\n";
                //  }

                if (ToukaAtomNum > 4) {
                    log<<" ToukaAtomNum is "<<ToukaAtomNum<<"\n";
                    log<<" Not support : ToukaAtomNum>4 "<<"\n";
                    log<<" You must write swap argorizm in this case."<<"\n";
                    CnErr.abort();
                }

                if (ToukaAtomNum==2) {
                    int pre  = ConectSerialAtomNum[0];
                    int post = ConectSerialAtomNum[1];
                    //      if( OutLevel < -5 ){
                    //        log<<"pre , post ==> "<< pre <<","<<post<<"\n";
                    //        log<<"Irekae Atom Number "<<pre<<" <-->  ";
                    //        log<<post<<"\n";
                    //      }
                    swap(pre,post,ConectAtomX,ConectAtomY,ConectAtomZ);
                } else if (ToukaAtomNum==3) {
                    //      if( OutLevel < -5 ){
                    //        log<<"debug ToukaAtomNum in 3 no type"<<"\n";
                    //      }

                    int pre = 0;
                    int post = 0;
                    switch (j) {
                    case 0:
                        pre  = ConectSerialAtomNum[0];
                        post = ConectSerialAtomNum[1];
                        break;
                    case 1:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 2:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 3:
                        pre  = ConectSerialAtomNum[0];
                        post = ConectSerialAtomNum[1];
                        break;
                    case 4:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 5:
                        break;
                    default:
                        log<<"Bad Number in ToukaAtomNum==3"<<"\n";
                        //                log<<"\7"<<"\n";
                        CnErr.abort();
                    }// switch

                    //      if( OutLevel < -5 ){
                    //        log<<"pre , post ==> "<< pre <<","<<post<<"\n";
                    //        log<<"Irekae Atom Number "<<pre<<" <-->  ";
                    //        log<<post<<"\n";
                    //      }

                    swap(pre,post,ConectAtomX,ConectAtomY,ConectAtomZ);

                } else if (ToukaAtomNum==4) {
                    //      if( OutLevel < -5 ){
                    //        log<<"debug ToukaAtomNum in 4 no type"<<"\n";
                    //      }
                    int pre = 0;
                    int post = 0;
                    switch (j) {
                    case 0:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 1:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 2:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 3:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 4:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 5:
                        pre  = ConectSerialAtomNum[0];
                        post = ConectSerialAtomNum[1];
                        break;
                    case 6:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 7:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 8:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 9:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 10:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 11:
                        pre  = ConectSerialAtomNum[0];
                        post = ConectSerialAtomNum[1];
                        break;
                    case 12:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 13:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 14:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 15:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 16:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 17:
                        pre  = ConectSerialAtomNum[0];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 18:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 19:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 20:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 21:
                        pre  = ConectSerialAtomNum[1];
                        post = ConectSerialAtomNum[2];
                        break;
                    case 22:
                        pre  = ConectSerialAtomNum[2];
                        post = ConectSerialAtomNum[3];
                        break;
                    case 23:
                        break;
                    default:
                        log<<"Bad Number in ToukaAtomNum==4"<<"\n";
                        //                 log<<"\7"<<"\n";
                        CnErr.abort();
                    }// switch

                    //      if( OutLevel < -5 ){
                    //        log<<"pre , post ==> "<< pre <<","<<post<<"\n";
                    //        log<<"Irekae Atom Number "<<pre<<" <-->  ";
                    //        log<<post<<"\n";
                    //      }
                    swap(pre,post,ConectAtomX,ConectAtomY,ConectAtomZ);

                }// if(ToukaAtomNum==4)
            }// if(j != (mawasi-1))

            //       if( OutLevel < -8 ){
            //    log<<" After Swap conectAtom"<<"\n";
            //    for(t=0;t<CntNum;t++){
            //      log<<ConectAtomX[t]<<"  "
            //         <<ConectAtomY[t]<<"  "
            //         <<ConectAtomZ[t]<<"\n";
            //    }
            //       }
        }

        // serch Minimum fitness
        minFitnessNum=0;
        minFitness = 999999.999999;
        for (int j = 0; j < mawasi; j++) {
            if (minFitness > fitness[j]) {
                minFitness = fitness[j] ;
                minFitnessNum = j ;
            }
        }

        //     if( OutLevel < -5 ){
        //       log<<"   MinFitNum = "<<minFitnessNum<<"\n";
        //       log<<"   OptAngle = "<<OptAngle[ minFitnessNum*UVN + 0 ]<<"   "
        //   <<OptAngle[ minFitnessNum*UVN + 1 ]<<"   "
        //   <<OptAngle[ minFitnessNum*UVN + 2 ]<<"\n";
        //     }

        henkanval = 3.14159265358979/180.0 ;
        // Oprimal Rotation Angle is next.
        radX = OptAngle[ minFitnessNum*UVN + 0 ] * henkanval ;
        radY = OptAngle[ minFitnessNum*UVN + 1 ] * henkanval ;
        radZ = OptAngle[ minFitnessNum*UVN + 2 ] * henkanval ;
        //---------------------------------------------------------------
        //条件を満たしていれば、Y、X、Z軸回りの順でθy,θx,θzだけ
        //回転するようにアフィンマトリクスを追加してやる。
        //これで、ある原子（CenterAtomNumの番号の原子）のマトリクスが
        //求まったことになる。
        //さらに、InpAmino.Order[CenterAtomNum]に１を代入しておく。
        //最後に warp_end へ 飛び越し命令を出す。
        //---------------------------------------------------------------
        //---------------------------------------------------------------
        //Y軸回りに radY[radian] だけ回転させるマトリクスを作る。
        //---------------------------------------------------------------
        initialyze_AfinMatrix(TransMatrix);
        initialyze_AfinMatrix(dumymat);
        clear_AfinMatrix(Matrix);

        row=0;
        col=0;
        dumymat[ row*4 + col ] =  cos(radY);
        row=0;
        col=2;
        dumymat[ row*4 + col ] =  sin(radY);
        row=2;
        col=0;
        dumymat[ row*4 + col ] = -sin(radY);
        row=2;
        col=2;
        dumymat[ row*4 + col ] =  cos(radY);

        multiply_afinmat(TransMatrix,dumymat,Matrix);
        dainyu_afinmat(TransMatrix,Matrix);
        //---------------------------------------------------------------
        //X軸回りに radX[radian] だけ回転させるマトリクスを作る。
        //---------------------------------------------------------------
        initialyze_AfinMatrix(dumymat);
        clear_AfinMatrix(Matrix);

        row=1;
        col=1;
        dumymat[ row*4 + col ] =  cos(radX);
        row=1;
        col=2;
        dumymat[ row*4 + col ] = -sin(radX);
        row=2;
        col=1;
        dumymat[ row*4 + col ] =  sin(radX);
        row=2;
        col=2;
        dumymat[ row*4 + col ] =  cos(radX);

        multiply_afinmat(TransMatrix,dumymat,Matrix);
        dainyu_afinmat(TransMatrix,Matrix);
        //---------------------------------------------------------------
        //Z軸回りに radZ[radian] だけ回転させるマトリクスを作る。
        //---------------------------------------------------------------
        initialyze_AfinMatrix(dumymat);
        clear_AfinMatrix(Matrix);

        row=0;
        col=0;
        dumymat[ row*4 + col ] =  cos(radZ);
        row=0;
        col=1;
        dumymat[ row*4 + col ] = -sin(radZ);
        row=1;
        col=0;
        dumymat[ row*4 + col ] =  sin(radZ);
        row=1;
        col=1;
        dumymat[ row*4 + col ] =  cos(radZ);

        multiply_afinmat(TransMatrix,dumymat,Matrix);
        dainyu_afinmat(TransMatrix,Matrix);
        //---------------------------------------------------------------
        // 上で求めておいたDbAmino.Atomdata[CenterAtomNum].afinmat と
        // 上はだめ.
        // 真、HozonMatrix[minFitnessNum]と
        // TransMatrixの積を行ない、CenterAtomNumの原子の最終的な
        // アフィンマトリクスとする。
        //---------------------------------------------------------------
        clear_AfinMatrix(dumymat);

        //     if( OutLevel < -7 ){
        //       log<<"minFitnessNum no Hozonmatrix "<<"\n";
        //       display_afinmatrix(HozonMatrix[minFitnessNum]);

        //       log<<"kake zan no tame no TransMatrix"<<"\n";
        //       display_afinmatrix(TransMatrix);
        //     }

        multiply_afinmat(HozonMatrix[minFitnessNum],TransMatrix,dumymat);
        dainyu_afinmat(DbAmino.Atomdata[CenterAtomNum].afinmat,dumymat);

        //     if( OutLevel < -7 ){
        //       log<<" Kaketa ato no Matrix "<<"\n";
        //       log<<" This AfinMatrix is used at Translation Coefficient."<<"\n";
        //       display_afinmatrix(DbAmino.Atomdata[CenterAtomNum].afinmat);
        //     }

        //変換する前に一番fitnessが小さい時の等価原子の順番に並び変える。
        dainyu_vector(ConectAtomX,EachStepX[minFitnessNum],CntNum);
        dainyu_vector(ConectAtomY,EachStepY[minFitnessNum],CntNum);
        dainyu_vector(ConectAtomZ,EachStepZ[minFitnessNum],CntNum);


        // 計算前の出力
        //     if( OutLevel < -7 ){
        //       log<<" Henkan Matrix %%%%%%%%%"<<"\n";
        //       display_afinmatrix(dumymat);
        //       log<<"\n"<<" ConectAtom Coordinate  X   Y   Z"<<"\n";
        //       for(t=0;t<CntNum;t++){
        //  log<<ConectAtomX[t]<<" "<<ConectAtomY[t]
        //     <<" "<<ConectAtomY[t]<<"\n";
        //       }
        //       log<<"\n";
        //     }

        //実際に最適角度を用いて座標変換する。
        for (int t=0; t < CntNum; t++) {
            dmX1[t] = ConectAtomX[t]*dumymat[0]
                      + ConectAtomY[t]*dumymat[4]
                      + ConectAtomZ[t]*dumymat[8];

            dmY1[t] = ConectAtomX[t]*dumymat[1]
                      + ConectAtomY[t]*dumymat[5]
                      + ConectAtomZ[t]*dumymat[9];

            dmZ1[t] = ConectAtomX[t]*dumymat[2]
                      + ConectAtomY[t]*dumymat[6]
                      + ConectAtomZ[t]*dumymat[10];
        }

        //------------------------------------------------------------
        // 条件を満たすかどうかの判定。
        // 条件は、誤差の平均がThreshold2以下であること、かつ、
        // それぞれの誤差がEachThreshold2以下であること、の２つ。
        // まず、それぞれの判定角度を求める。
        //------------------------------------------------------------
        TotalAngleGosa=0.0;
        MaxAngleGosa=1E-10;
        for (int t = 0; t < CntNum; t++) {
            x0 = dmX1[t];
            y0 = dmY1[t];
            z0 = dmZ1[t];
            x1 = InpAtomX[t];
            y1 = InpAtomY[t];
            z1 = InpAtomZ[t];
            // angle is degree.

            eachAngle[t] = calc_angle(x0,y0,z0,x1,y1,z1);

            //       if( OutLevel < -5 ){
            //    log<<"eachAngle[ "<<t<<" ] = "<<eachAngle[t]<<"          ";
            //    log<<"distance = "<< sqrt( pow(x0-x1,2) + pow(y0-y1,2)
            //                   +pow(z0-z1,2) )<<"\n";
            //       }

            if (eachAngle[t] > MaxAngleGosa) {
                MaxAngleGosa = eachAngle[t] ;
            }

            TotalAngleGosa += eachAngle[t] ;
        }
        AverageAngleGosa = TotalAngleGosa / CntNum ;

        //     if( OutLevel < -5 ){
        //       log<<" ------------- Last judge ------------- "<<"\n";
        //       log<<"  Total  AngleGosa  = "<<TotalAngleGosa<<"\n";
        //       log<<"  AverageAngleGosa  = "<<AverageAngleGosa<<"\n";
        //       log<<"  MaxAngleGosa      = "<<MaxAngleGosa <<"\n";
        //     }

        //------------------------------------------------------------
        //条件を満たすかの判定。
        //------------------------------------------------------------
        if (AverageAngleGosa <= Threshold2) {
            int Agflag = 0;
            for (int t = 0; t < CntNum; t++) {
                Agflag=0;

                if (eachAngle[t] <= EachThreshold2) {
                    Agflag=1;
                } else {
                    InpAmino.Order[CenterAtomNum] = -1;
                    //      if( OutLevel < -4 ){
                    //        log<<"\n";
                    //        log<<"########################################"<<"\n";
                    //        log<<"#                                      #"<<"\n";
                    //        log<<"#     Not Accept,order !!              #"<<"\n";
                    //        log<<"#  Order [ "<<CenterAtomNum<<" ] = -1  #"<<"\n";
                    //        log<<"#                                      #"<<"\n";
                    //        log<<"########################################"<<"\n";
                    //      }
                    goto end_serch;
                }
            }

            if (Agflag==1) {
                InpAmino.Order[CenterAtomNum] = 1;
                //  if( OutLevel < -4 ){
                //    log<<"\n";
                //    log<<"##########################################"<<"\n";
                //    log<<"#                                        #"<<"\n";
                //    log<<"      OK!! InpAmino.Order[";
                //    log<<       CenterAtomNum<<"] = 1 ;            "<<"\n";
                //    log<<"#                                        #"<<"\n";
                //    log<<"##########################################"<<"\n";
                //  }
                goto end_serch;
            }
        } else {
            InpAmino.Order[CenterAtomNum] = -1;
            //       if( OutLevel < -4 ){
            //    log<<"\n";
            //    log<<"########################################"<<"\n";
            //    log<<"#                                      #"<<"\n";
            //    log<<"#     Not Accept,order !!              #"<<"\n";
            //    log<<"#  Order [ "<<CenterAtomNum<<" ] = -1  #"<<"\n";
            //    log<<"#                                      #"<<"\n";
            //    log<<"########################################"<<"\n";
            //       }
            goto end_serch;
        }

end_serch:;    //

warp_end:;     //

        //pre_iLoopEnd:; //

        //     if( OutLevel < -9 ){
        //       log<<"Loop i no Number is [ " <<i<<" ]"<<"\n";
        //     }

    }//for (i) ///////////////////////////////////////////////////////

    //   if( OutLevel < -3 ){
    //     log<<"\n"<<"\n";
    //     log<<"    -------------------------------------------"<<"\n";
    //     log<<"    ----- makeAfinmatrix(  int  ) Exit  -------"<<"\n";
    //     log<<"    -------------------------------------------"<<"\n"<<"\n";
    //   }

    return 0;
}
//#############################################################################
void DfInitialguess::powell(double* U0 ,
                            double* coreX , double* coreY , double* coreZ ,
                            double* inpX ,  double* inpY ,  double* inpZ ,
                            int CntNum)
{
    TlLogX& log = TlLogX::getInstance();

    int i,j,k; //int loop=1;
    int   iteration,MaxIteration;
    double alpha,F1,F2,F3,max,Inf,value; // double min_point
    int    UVN,M; // nknown Variable Number is Three(ThetaX,ThetaY,ThetaZ);
    enum { MaxUVN = 3 };
    M = UVN = MaxUVN ;  // Unkown Valraible number;
    Inf=1E+20 ;      // flag for variable infinity
    MaxIteration = 100;
    double Xi[MaxUVN][MaxUVN];
    double U[MaxUVN][MaxUVN];
    double U0_old[MaxUVN];
    double x[MaxUVN];
    double a[MaxUVN];
    double new_Xi[MaxUVN];
    double delta[MaxUVN];

    //int    kaiten;

//   if( OutLevel < -5 ){
//     log<<"First J(Kyori no gosa)in powell = "
//        <<AfinFunction(U0,coreX,coreY,coreZ,inpX,inpY,inpZ,CntNum)<<"\n";
//     log.flush();
//   }

    // Initialize data.
    for (i=0; i<M; i++) {
        for (j=0; j<M; j++) {
            Xi[i][j] = 0.0 ;
        }
        Xi[i][i] = 1.0 ;
        U0_old[i]=0.0;
    }

    iteration = 0;
    while (1) {
        iteration++;
        /*  ****** (1) ******   */
        for (i=0; i<M; i++) {
            for (j=0; j<M; j++) {
                x[j]=Xi[j][i];
            }
            a[i]=rsf(U0,0.0,x,coreX,coreY,coreZ,inpX,inpY,inpZ,CntNum);
            if (fabs(a[i])>=Inf) {
                log<<" go to infinity #1"<<"\n";
                log<<" one variable serch routine rsf(,,,,)"<<"\n";
                log<<" powell called rsf "<<"\n";
                log.flush();
                break;
            }
        }

        for (i=0; i<M; i++) {
            U[i][0]=U0[i]+a[0]*Xi[i][0];
            for (j=1; j<M; j++) {
                U[i][j]=U[i][j-1]+a[j]*Xi[i][j];
            }
        }

        /*   ****** (2) ******    */
        for (i=0; i<M; i++) {
            x[i]=U[i][0];
        }
        delta[0] = AfinFunction(U0,coreX,coreY,coreZ,
                                inpX,inpY,inpZ,CntNum) -
                   AfinFunction(x,coreX,coreY,coreZ,
                                inpX,inpY,inpZ,CntNum)   ;
        for (i=1; i<M; i++) {
            for (j=0; j<M; j++) {
                x[j]=U[j][i-1];
            }
            delta[i]=AfinFunction(x,coreX,coreY,coreZ,
                                  inpX,inpY,inpZ,CntNum);
            for (j=0; j<M; j++) {
                x[j]=U[j][i];
            }
            delta[i]=delta[i]-AfinFunction(x,coreX,coreY,coreZ,
                                           inpX,inpY,inpZ,CntNum);
        }
        max = delta[0];
        k=0;
        for (i=1; i<M; i++) {
            if (max<delta[i]) {
                k=i;
                max=delta[i];
            }
        }

        /*   ****** (3) ******   */
        F1 = AfinFunction(U0,coreX,coreY,coreZ,
                          inpX,inpY,inpZ,CntNum);
        for (i=0; i<M; i++) {
            x[i]=U[i][M-1];
        }
        F2 = AfinFunction(x,coreX,coreY,coreZ,
                          inpX,inpY,inpZ,CntNum);
        for (i=0; i<M; i++) {
            x[i] = 2*U[i][M-1]-U0[i];
        }
        F3 = AfinFunction(x,coreX,coreY,coreZ,
                          inpX,inpY,inpZ,CntNum);


        /*   ****** (4) ******   */
        if ((F1 <= F3) ||
                (1/2*max*(F1-F3)*(F1-F3))<=((F1-2*F2+F3)*(F1-F2-max)*(F1-F2-max))) {
            for (i=0; i<M; i++) {
                U0[i]=U[i][M-1];
            }
        } else {
            for (i=0; i<M; i++) {
                new_Xi[i]=U[i][M-1]-U0[i];
                x[i]=U[i][M-1];
            }
            alpha=rsf(x,0.0,new_Xi,coreX,coreY,coreZ,
                      inpX,inpY,inpZ,CntNum);
            if ((fabs(alpha)>=Inf)) {
                log<<" go to infinity #2"<<"\n";
                log<<" one variable serch routine rsf(,,,,)"<<"\n";
                log<<" powell called rsf "<<"\n";
                log.flush();
                break;
            }
            for (i=k; i<M-1; i++) {
                for (j=0; j<M; j++) {
                    Xi[j][i]=Xi[j][i+1];
                }
            }
            for (i=0; i<M; i++) {
                Xi[i][M-1]=new_Xi[i];
                U0[i]=U[i][M-1]+alpha*new_Xi[i];
            }
        }// else

        max = fabs(U0_old[0] - U0[0]);

        for (i=1; i<M; i++) {
            value = fabs(U0_old[i]-U0[i]) ;
            if (max < value) {
                max = value ;
            }
        }

        if (max >= Inf) {
            log<<" go to infinity #3"<<"\n";
            log<<" one variable serch routine rsf(,,,,)"<<"\n";
            log<<" powell called rsf "<<"\n";
            log.flush();
            break;
        } else if ((max <= 1e-4)) {
            log<<"\n";
            log<<" ######################"<<"\n";
            log<<" ## optimal solution ##"<<"\n";
            log<<" ######################"<<"\n";
            log<<"\n";
            log.flush();
            for (i=0; i<M; i++) {
                log<<"Uopt["<<i<<"] = "<<U0[i]<<"\n";
            }
            break;
        } else {
            AngleTrans(U0,M);
            for (i=0; i<M; i++) {
                U0_old[i] = U0[i];
            }
//       if( OutLevel < -7 ){
//  log<<"Now J(Kyori no gosa) = "
//     <<AfinFunction(U0_old,coreX,coreY,coreZ,
//            inpX,inpY,inpZ,CntNum)<<"\n";
//       }
        }//else

        if (iteration>MaxIteration) {
            log<<" Iteration is max."<<"\n";
            log<<" So , This ruotine  exit by force."<<"\n";
            break;
        }

        //  log<<"powell::loop no saigo"<<"\n";

    } //while(1);

//   if( OutLevel < -5 ){
//     log<<"Last Powell Mehtod Iteration is [ "<<iteration<<" ] "<<"\n";
//     log<<"Now J(Last no gosa) = "
//        <<AfinFunction(U0_old,coreX,coreY,coreZ,
//            inpX,inpY,inpZ,CntNum)<<"\n";
//   }

} // end.
//######################################################################
void DfInitialguess::AngleTrans(double* U0,int M)
{
    int i;
    double kaiten,henkan;

    henkan = 3.14159265358979/180.0 ;

    for (i=0; i<M; i++) {
        if (((U0[i]>=0.0) && (U0[i]<=180.0)) ||
                ((U0[i]>-180.0) && (U0[i]<0.0))) {
        } else {
            kaiten = floor(U0[i] / 180.0);
            if (sin(U0[i]*henkan) > 0.0) {
                if (U0[i]>0.0) {
                    U0[i] = U0[i] - kaiten*180.0 ;
                } else {
                    kaiten = fabs(kaiten);
                    U0[i] = U0[i] + kaiten*180.0 ;
                }
            } else {
                if (U0[i]>0.0) {
                    U0[i] = U0[i] - (kaiten+1)*180.0 ;
                } else {
                    kaiten = fabs(kaiten);
                    U0[i] = U0[i] + (kaiten-1)*180.0 ;
                }
            }

        }

    }
}
//#############################################################################
double DfInitialguess::AfinFunction(double* angle,
                                    double* X ,   double* Y ,   double* Z ,
                                    double* inpX ,double* inpY ,double* inpZ,
                                    int CntNum)
{
    // double* angle is degree unit.
    int i;
    double Tx,Ty,Tz;
    double  value,henkanval;
    double Xd[MaxCntNum];
    double Yd[MaxCntNum];
    double Zd[MaxCntNum];
    henkanval = 3.141592653538979/180.0 ;
    Tx = angle[0]*henkanval ;
    Ty = angle[1]*henkanval ;
    Tz = angle[2]*henkanval ;

    for (i=0; i<CntNum; i++) {
        Xd[i] =  X[i] * (cos(Ty)*cos(Tz) + sin(Tx)*sin(Ty)*sin(Tz))
                 + Y[i] * (cos(Tx)*sin(Tz))
                 + Z[i] * (-sin(Ty)*cos(Tz) + sin(Tx)*cos(Ty)*sin(Tz)) ;

        Yd[i] =  X[i] * (-cos(Ty)*sin(Tz) + sin(Tx)*sin(Ty)*cos(Tz))
                 + Y[i] * (cos(Tx)*cos(Tz))
                 + Z[i] * (sin(Ty)*sin(Tz) + sin(Tx)*cos(Ty)*cos(Tz)) ;

        Zd[i] =  X[i] * (+cos(Tx)*sin(Ty))
                 + Y[i] * (-sin(Tx))
                 + Z[i] * (cos(Tx)*cos(Ty)) ;
    }

    value = 0.0;
    for (i=0; i<CntNum; i++) {
        value += sqrt((Xd[i]-inpX[i])*(Xd[i]-inpX[i]) +
                      (Yd[i]-inpY[i])*(Yd[i]-inpY[i]) +
                      (Zd[i]-inpZ[i])*(Zd[i]-inpZ[i]));
    }

    return value;
}
//#############################################################################
double DfInitialguess::rsf(double* U , double a ,double* Xi,
                           double* coreX , double* coreY , double* coreZ ,
                           double* inpX  , double* inpY  , double* inpZ ,
                           int CntNum)
{
    int i;
    double  alpha,beta,delta,min,min_a,new_point,h,Inf;
    int    M ;
    enum { MaxUVN = 3 };
    M = MaxUVN; // Unknown Variable Number is Three(ThetaX,ThetaY,ThetaZ);
    Inf = 1E20;
    double x[MaxUVN];

    h = 1.0;         // first kizami haba.
    alpha = 5.0;     // serch +vector
    beta = 0.5;      // serch -vector
    delta = 1E-5; // Convergence criterion

    for (i=0; i<M; i++) {
        x[i]=U[i]+a*Xi[i];
    }
    min = AfinFunction(x,coreX,coreY,coreZ,inpX,inpY,inpZ,CntNum);
    for (i=0; i<M; i++) {
        x[i]=U[i]+(a+h)*Xi[i];
    }

    new_point = AfinFunction(x,coreX,coreY,coreZ,inpX,inpY,inpZ,CntNum);

    while (1) {
        if (new_point < min) {
            a = a+h;
            min = new_point;
            h = alpha*h;
        } else if (fabs(h)<delta) {
            min_a = a;
            break;
        } else {
            h = -beta*h;
        }

        for (i=0; i<M; i++) {
            x[i]=U[i]+(a+h)*Xi[i];
        }

        new_point=AfinFunction(x,coreX,coreY,coreZ,inpX,inpY,inpZ,CntNum);

        if (fabs(new_point) >= Inf) {
            //log<<new_point<<"\n";
            min_a = new_point;
            break;
        }

    } // while(1);

    return min_a;
}
//#############################################################################
void DfInitialguess::swap(int pre, int post, double* X, double* Y, double* Z)
{
//   if( OutLevel < -8 ){
//     log<<"in Swap Procedure"<<"\n";
//     log<<"Before Swap"<<"\n";
//     log<<X[pre] <<"     "<< Y[pre]<<"    "<<Z[pre]<<"\n";
//     log<<X[post]<<"     "<<Y[post]<<"    "<<Z[post]<<"\n";
//   }

//   double dumyval = X[pre];
//   X[pre]  = X[post];
//   X[post] = dumyval;
    std::swap(X[pre], X[post]);

//   dumyval = Y[pre];
//   Y[pre]  = Y[post];
//   Y[post] = dumyval;
    std::swap(Y[pre], Y[post]);

//   dumyval = Z[pre];
//   Z[pre]  = Z[post];
//   Z[post] = dumyval;
    std::swap(Z[pre], Z[post]);

//   if( OutLevel < -8 ){
//     log<<"After Swap"<<"\n";
//     log<<X[pre] <<"     "<< Y[pre]<<"    "<<Z[pre]<<"\n";
//     log<<X[post]<<"     "<<Y[post]<<"    "<<Z[post]<<"\n";
//   }

}
//#############################################################################
double DfInitialguess::calc_angle(double x0 , double y0 , double z0 ,
                                  double x1 , double y1 , double z1)
{
    TlLogX& log = TlLogX::getInstance();

    double  Angle_degree,Angle_radian;
    double  Naiseki,Norm0,Norm1,CosTheta;

    Norm0   = sqrt(pow(x0,2) + pow(y0,2) + pow(z0,2)) ;
    Norm1   = sqrt(pow(x1,2) + pow(y1,2) + pow(z1,2)) ;
    if (Norm0 < 1E-6) {
        log<<"in gusCalcAngle, Bad Norm0 = "<<Norm0<<"\n";
        log<<"vec = ("<<x0<<","<<y0<<","<<z0<<")   ,("<<x1<<","<<y1<<","<<z1
        <<")"<<"\n";
        //  log<<"\7"<<"\n";
        CnErr.abort();
    }
    if (Norm1 < 1E-6) {
        log<<"in gusCalcAngle, Bad Norm1 = "<<Norm1<<"\n";
        log<<"vec = ("<<x0<<","<<y0<<","<<z0<<")   ,("<<x1<<","<<y1<<","<<z1
        <<")"<<"\n";
        //  log<<"\7"<<"\n";
        CnErr.abort();
    }

    x0 /= Norm0;
    y0 /= Norm0;
    z0 /= Norm0;
    x1 /= Norm1;
    y1 /= Norm1;
    z1 /= Norm1;
    Naiseki = x0*x1 + y0*y1 + z0*z1 ;
    CosTheta = Naiseki ;

    if ((CosTheta < -1.0000001) || (CosTheta > 1.0000001)) {
        log<<"Error in DfInitialguess::calc_angle"<<"\n";
        //  log<<"Bad value CosTheta"<<"\7"<<"\n";
        log<<"Bad value CosTheta"<<"  "<<"\n";
        log<<"CosTheta = "<<CosTheta<<"\n";
        CnErr.abort();
    }
    Angle_radian = acos(CosTheta) ;
    Angle_degree = Angle_radian* 180.0 / 3.14159265358979;

    return (Angle_degree);

}
//#############################################################################
void DfInitialguess::multiply_afinmat(double* matA,double* matB,double* matC)
{
    int i,j,k;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            for (k=0; k<4; k++) {
                matC[i*4+j] += (matA[i*4+k]*matB[k*4+j]);
            }
        }
    }
}
void DfInitialguess::multiply_afinmat(std::vector<double>& matA,double* matB,double* matC)
{
    int i,j,k;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            for (k=0; k<4; k++) {
                matC[i*4+j] += (matA[i*4+k]*matB[k*4+j]);
            }
        }
    }
}
//#############################################################################
void DfInitialguess::dainyu_afinmat(double* matA,double* matB)
{
    int i,j;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            matA[i*4+j] = matB[i*4+j] ;
        }
    }
}
void DfInitialguess::dainyu_afinmat(std::vector<double>& matA, double* matB)
{
    int i,j;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            matA[i*4+j] = matB[i*4+j] ;
        }
    }
}
void DfInitialguess::dainyu_afinmat(double* matA, std::vector<double>& matB)
{
    int i,j;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            matA[i*4+j] = matB[i*4+j] ;
        }
    }
}
//#############################################################################
void DfInitialguess::initialyze_AfinMatrix(double* matA)
{
    int k;
    for (k=0; k<16; k++) {
        matA[k]=0.0;
    }
    for (k=0; k<4; k++) {
        matA[k*4+k]=1.0;
    }
}
void DfInitialguess::initialyze_AfinMatrix(std::vector<double>& matA)
{
    int k;
    for (k=0; k<16; k++) {
        matA[k]=0.0;
    }
    for (k=0; k<4; k++) {
        matA[k*4+k]=1.0;
    }
}

//#############################################################################
void DfInitialguess::clear_AfinMatrix(double* matA)
{
    int i;
    for (i=0; i<16; i++) {
        matA[i]=0.0;
    }
}
//#############################################################################
void DfInitialguess::display_afinmatrix(double* matA)
{
    TlLogX& log = TlLogX::getInstance();

    int i,j;
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            //log<<"   "<<setw(10)<< matA[i*4+j] <<"  ";
//       log.precision(10);
//       log.setf(ios::fixed,ios::floatfield);
            log<<"   "<< matA[i*4+j] ;
        }
        log<<"\n";
    }
    log<<"\n";
}
//#############################################################################
void DfInitialguess::dainyu_vector(double* vec1,double* vec2,int num)
{
    int i;
    for (i=0; i<num; i++) {
        vec1[i] = vec2[i] ;
    }
}
//#############################################################################
int DfInitialguess::transOrbital(int tAN)
{
    TlLogX& log = TlLogX::getInstance();

//   if( OutLevel < -3 ){
//     log<<"\n";
//     log<<"####################################################### "<<"\n";
//     log<<"##############   TransOrbital Start  ################## "<<"\n";
//     log<<"####################################################### "<<"\n";
//     log<<"\n";
//     log<<"         Total AtomNum of One Residue = "<<tAN << "\n"<< "\n";
//   }

//   if( OutLevel < -8 ){
//     log<<"Initial w1RouCount = "<<w1RouCount<<"\n";
//     log<<"Initial w1MyuCount = "<<w1MyuCount<<"\n";
//     log<<"Initial w1NyuCount = "<<w1NyuCount<<"\n";
//   }

    int   Sorbnum,Porbnum,Dorbnum,Si,Pi,Di,i,Snum,Pnum,Dnum,SPDnum;
    int   row,col,Rct,Mct,Nct;
    double Cxx,Cxy,Cxz,Cyx,Cyy,Cyz,Czx,Czy,Czz;
    int   NormNum;
    double NormFact,NormF_2,NormF_3;
    double Rs,InvRoot3;
    Rs = 1.0 / 6.0  ;          // for d_Orbtial;  ( 1/6 = 0.16666666..... )
    InvRoot3 = 1.0 / sqrt(3.0);  // for d_Orbital;  ( 1/sqrt(3) )

    Fl_Tbl_Density  FTD;
    Fl_Tbl_Xcpot    FTX;
    Fl_Tbl_Xcpot    FTX2;

    Fl_Gto_Density  FGD;
//   FGD.open("fl_Gto_Density","read");
//   FGD.read();
//   FGD.close();

    Fl_Gto_Xcpot    FGX;
//   FGX.open("fl_Gto_Xcpot","read");
//   FGX.read();
//   FGX.close();

    Fl_Gto_Xcpot2   FGX2;
    ///// modified by Fric(mizouchi) 2001/12/09 (s) /////
    //    FGX2.open("fl_Gto_Xcpot2","read");
//   FGX2.open("fl_Gto_Xcpot","read");
    ///// modified by Fric(mizouchi) 2001/12/09 (e) /////

//   FGX2.read();
//   FGX2.close();

    //----------------------------- Loop : i -----------------------
    // i : その残基の原子数で回る。
    //-------------------------------------------------------------

    for (i=0; i<tAN; i++) {
        Rct=0;
        Mct=0;
        Nct=0;
        //-----------------------------------------------------------------
        //
        // DbAmino.Atomdata[i].afinmatから回転マトリクスだけを取り出して
        // 次のように代入しておく。
        // +-                -+
        // |  Cxx   Cxy   Cxz |                           [0,0] [0,1] [0,2]
        // |  Cyx   Cyy   Cyz | <======================   [1,0] [1,1] [1,2]
        // |  Czx   Czy   Czz |                           [2,0] [2,1] [2,2]
        // +-                -+
        //
        //-----------------------------------------------------------------

        row=0;
        col=0;
        Cxx = DbAmino.Atomdata[i].afinmat[row*4+col];
        row=1;
        col=0;
        Cyx = DbAmino.Atomdata[i].afinmat[row*4+col];
        row=2;
        col=0;
        Czx = DbAmino.Atomdata[i].afinmat[row*4+col];

        row=0;
        col=1;
        Cxy = DbAmino.Atomdata[i].afinmat[row*4+col];
        row=1;
        col=1;
        Cyy = DbAmino.Atomdata[i].afinmat[row*4+col];
        row=2;
        col=1;
        Czy = DbAmino.Atomdata[i].afinmat[row*4+col];

        row=0;
        col=2;
        Cxz = DbAmino.Atomdata[i].afinmat[row*4+col];
        row=1;
        col=2;
        Cyz = DbAmino.Atomdata[i].afinmat[row*4+col];
        row=2;
        col=2;
        Czz = DbAmino.Atomdata[i].afinmat[row*4+col];

//     if( OutLevel < -7 ){
//       log<<"\n"<<"TransMatrix-----------------------------------"<<"\n";
//       log<<Cxx<<"    "<<Cxy<<"    "<<Cxz<<"\n";
//       log<<Cyx<<"    "<<Cyy<<"    "<<Cyz<<"\n";
//       log<<Czx<<"    "<<Czy<<"    "<<Czz<<"\n";
//     }

        if (InpAmino.Order[i] ==  1) {
//       if( OutLevel < -3 ){
//  log<<" Execute : Translation1  ---> ";
//  log<<" Use AfinMatrix (S,P,D)_Orbital Translation"<<"\n";
//       }
            goto Translation1;
        } else if (InpAmino.Order[i] == -1) {
//       if( OutLevel < -3 ){
//  log<<"Execute  : Translation2  ---> ";
//  log<<" ( S -> S   ,   P -> Zero  ,  D -> Zero )"<<"\n";
//       }
            goto Translation2;
        } else {
            log<<"Bad value : InpAmino.Order = "<<InpAmino.Order[i]<<"\n";
            //        log<<"\7"<<"\n";
            CnErr.abort();
        }

        //############//
        //
Translation1:;   //        InpAmino.Order is 1 .
        //
        //############//

        //***********************************************************
        // Rouについての変換。
        //***********************************************************
        //-------------------------------------------------------
        //-------------------------------------------------------
        //   ここから実際の軌道s,p,d（ρ、μ、νの展開係数）を変換する。
        //-------------------------------------------------------
        //-------------------------------------------------------
        Snum   = DbAmino.Atomdata[i].rouDbShellNum[0] ;
        Pnum   = DbAmino.Atomdata[i].rouDbShellNum[1] ;
        Dnum   = DbAmino.Atomdata[i].rouDbShellNum[2] ;
        SPDnum = DbAmino.Atomdata[i].rouDbShellNum[3] ;

#if VER1_0
        Sorbnum = Snum + SPDnum ;
#endif

#if VER1_1
        Sorbnum = Snum + SPDnum + SPDnum;
#endif

        Porbnum = Pnum + SPDnum ;
        Dorbnum = Dnum + SPDnum ;

        //-------------------------------------------------------
        //  S 軌道は変換なしに、そのままにする。   SSSSSSSSSSSSS
        //-------------------------------------------------------
        for (Si=0; Si<Sorbnum; Si++) {
            w1Rou[Si+w1RouCount] = DbAmino.Atomdata[i].CoefRou[Si] ;
        }
        w1RouCount += Sorbnum ;
        Rct += Sorbnum;
        //-------------------------------------------------------
        //  P 軌道の変換を行なう。                 PPPPPPPPPPPPP
        //                                          +-           -+
        //                                          | Cxx Cxy Cxz |
        //  [ Cpx' Cpy' Cpz' ] = [ Cpx Cpy Cpz ] x  | Cyx Cyy Cyz |
        //                                  | Czx Czy Czz |
        //                                  +-           -+
        //-------------------------------------------------------
        for (Pi=0; Pi<Porbnum; Pi++) {
            // ----- for Px -----;
            w1Rou[Pi*3+w1RouCount]  =
                DbAmino.Atomdata[i].CoefRou[Pi*3+Rct]  *Cxx
                + DbAmino.Atomdata[i].CoefRou[Pi*3+1+Rct]*Cyx
                + DbAmino.Atomdata[i].CoefRou[Pi*3+2+Rct]*Czx;

            // ----- for Py -----;
            w1Rou[Pi*3+1+w1RouCount]=
                DbAmino.Atomdata[i].CoefRou[Pi*3+Rct]  *Cxy
                + DbAmino.Atomdata[i].CoefRou[Pi*3+1+Rct]*Cyy
                + DbAmino.Atomdata[i].CoefRou[Pi*3+2+Rct]*Czy;

            // ----- for Pz -----;
            w1Rou[Pi*3+2+w1RouCount]=
                DbAmino.Atomdata[i].CoefRou[Pi*3+Rct]  *Cxz
                + DbAmino.Atomdata[i].CoefRou[Pi*3+1+Rct]*Cyz
                + DbAmino.Atomdata[i].CoefRou[Pi*3+2+Rct]*Czz;

        }
        w1RouCount += (Porbnum*3);
        Rct += (Porbnum*3) ;

        //---------------------------------------------------------------------
        //  D 軌道の変換を行なう。                 DDDDDDDDDDDDD
        //
        //  d(xy),d(xz),d(yz),d(x^2 - y^2),d(3z^2 - r^2)を変換することで
        //  それぞれ５つの成分に分解して、
        //  最後にそれぞれの成分を足したものを変換後の軌道とする。
        //
        //  ただし、規格化定数を考慮しなければならない。
        //  d(xy)を規格化する数をNormFactとしたとき、
        //  d(x2-y2)の規格化定数 NormF_2 は NormFact/2.0
        //  となる。さらにd(3z2-r2)の規格化定数 NormF_3 は
        //  NormFact/sqrt(3)となる。
        //
        //  ただし、ρはFl_Gtoのメンバ関数 getCoulonbnormalized(.....)を使い、
        //  μ、νはgetNormalized(.....)を用いて規格化する。
        //
        //---------------------------------------------------------------------
        //
        //  d(xy)' =     (CxxCyy + CxyCyx)            * d(xy)      * NormFact +
        //               (CxxCzy + CxyCzx)            * d(xz)      * NormFact +
        //               (CyxCzy + CyyCzx)            * d(yz)      * NormFact +
        //         1/2 * (CxxCxy - CyxCyy)            * (x^2 - y^2)* NormF_2  -
        //         1/6 * (CxxCxy + CyxCyy - 2CzxCzy)  * (3z^2 - r^2)* NormF_3 ;
        //
        //  求まったd(xy)'をさらに規格化定数で割ってやる。
        //
        //  d(xy)' = d(xy)'/NormFact ;
        //
        //---------------------------------------------------------------------
        //
        //  d(xz)' =     (CxxCyz + CxzCyx)            * d(xy)      * NormFact +
        //               (CxxCzz + CxzCzx)            * d(xz)      * NormFact +
        //               (CyxCzz + CyzCzx)            * d(yz)      * NormFact +
        //         1/2 * (CxxCxz - CyxCyz)            * (x^2 - y^2)* NormF_2  -
        //         1/6 * (CxxCxz + CyxCyz - 2CzxCzz)  * (3z^2 - r^2)* NormF_3 ;
        //
        //  求まったd(xz)'をさらに規格化定数で割ってやる。
        //
        //  d(xz)' = d(xz)'/NormFact ;
        //
        //---------------------------------------------------------------------
        //
        //  d(yz)' =     (CxyCyz + CxzCyy)            * d(xy)      * NormFact +
        //               (CxyCzz + CxzCzy)            * d(xz)      * NormFact +
        //               (CyyCzz + CyzCzy)            * d(yz)      * NormFact +
        //         1/2 * (CxyCxz - CyyCyz)            * (x^2 - y^2)* NormF_2  -
        //         1/6 * (CxyCxz + CyyCyz - 2CzyCzz)  * (3z^2 - r^2)* NormF_3 ;
        //
        //  求まったd(xz)'をさらに規格化定数で割ってやる。
        //
        //  d(yz)' = d(yz)'/NormFact ;
        //
        //---------------------------------------------------------------------
        //
        //  d(x^2 - y^2) =
        //               2 * (CxxCyx - CxyCyy)              * d(xy)* NormFact +
        //               2 * (CxxCzx - CxyCzy)              * d(xz)* NormFact +
        //               2 * (CyxCzx - CyyCzy)              * d(yz)* NormFact +
        //         1/2 * (CxxCxx - CyxCyx - CxyCxy +  CyyCyy)* (x^2 - y^2)
        //             * NormF_2  -
        //         1/6 * (CxxCxx + CyxCyx - 2CzxCzx -
        //                CxyCxy - CyyCyy + 2CzyCzy   )*(3z^2 - r^2) * NormF_3;
        //
        //  求まったd(x2-y2)'をさらに規格化定数で割ってやる。
        //
        //  d(x2-y2)' = d(x2-y2)'/NormF_2 ;
        //
        //
        //---------------------------------------------------------------------
        //
        //  d(3z^2 - r^2) = (4CxzCyz - 2CxxCyx - 2CxyCyy) * d(xy) * NormFact  +
        //                  (4CxzCzz - 2CxxCzx - 2CxyCzy) * d(xy) * NormFact  +
        //                  (4CyzCzz - 2CyxCzx - 2CyyCzy) * d(xy) * NormFact  +
        //           1/2 *  (2CxzCxz - 2CyzCyz - CxxCxx +
        //                    CyxCyx - CxyCxy  + CyyCyy ) * d(x2-y2)* NormF_2 -
        //           1/6 *  (2CxzCxz + 2CyzCyz - 4CzzCzz -
        //                   CxxCxx  - CyxCyx  + 2CzxCzx -
        //           CxyCxy  - CyyCyy  + 2CzyCzy) * d(3z2-r2)* NormF_3;
        //
        //  求まったd(3z2-r2)'をさらに規格化定数で割ってやる。
        //
        //  d(3z2-r2)' = d(3z2-r2)'/NormF_3 ;
        //
        //
        //---------------------------------------------------------------------
        for (Di=0; Di<Dorbnum; Di++) {
            NormNum  = FTD.getcgtonum(w1RouCount);
            NormFact = FGD.getCoulombnormalized(NormNum,0,1,1,0);
            NormF_2  = NormFact * 0.5 ;                    // NormFact/2 ;
            NormF_3  = NormFact * InvRoot3 * 0.5 ;         // NormFact/(2*sqrt(3));

            NormFact = 1.0 / NormFact ;
            NormF_2  = 1.0 / NormF_2  ;
            NormF_3  = 1.0 / NormF_3  ;
            //----- for d(xy) -----;
            w1Rou[Di*5+w1RouCount]   =
                (Cxx*Cyy + Cxy*Cyx)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+Rct])
                + (Cxx*Czy + Cxy*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+1+Rct])
                + (Cyx*Czy + Cyy*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+2+Rct])
                + (0.5)*NormF_2*(Cxx*Cxy - Cyx*Cyy)*
                (DbAmino.Atomdata[i].CoefRou[Di*5+3+Rct])
                -1*Rs*(Cxx*Cxy + Cyx*Cyy - 2*Czx*Czy)*NormF_3*
                (DbAmino.Atomdata[i].CoefRou[Di*5+4+Rct]) ;

            w1Rou[Di*5+w1RouCount] = w1Rou[Di*5+w1RouCount] / NormFact ;

            //----- for d(xz) -----;
            w1Rou[Di*5+1+w1RouCount] =
                (Cxx*Cyz + Cxz*Cyx)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+Rct])
                + (Cxx*Czz + Cxz*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+1+Rct])
                + (Cyx*Czz + Cyz*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+2+Rct])
                + (0.5)*NormF_2*(Cxx*Cxz - Cyx*Cyz)*
                (DbAmino.Atomdata[i].CoefRou[Di*5+3+Rct])
                -1*Rs*(Cxx*Cxz + Cyx*Cyz - 2*Czx*Czz)*NormF_3*
                (DbAmino.Atomdata[i].CoefRou[Di*5+4+Rct]);

            w1Rou[Di*5+1+w1RouCount] = w1Rou[Di*5+1+w1RouCount] / NormFact ;

            //----- for d(yz) -----;
            w1Rou[Di*5+2+w1RouCount] =
                (Cxy*Cyz + Cxz*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+Rct])
                + (Cxy*Czz + Cxz*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+1+Rct])
                + (Cyy*Czz + Cyz*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+2+Rct])
                + (0.5)*NormF_2*(Cxy*Cxz - Cyy*Cyz)*
                (DbAmino.Atomdata[i].CoefRou[Di*5+3+Rct])
                -1*Rs*(Cxy*Cxz + Cyy*Cyz - 2*Czy*Czz)*NormF_3*
                (DbAmino.Atomdata[i].CoefRou[Di*5+4+Rct]);

            w1Rou[Di*5+2+w1RouCount] = w1Rou[Di*5+2+w1RouCount] / NormFact ;

            //----- for d(x^2 - y^2) -----;
            w1Rou[Di*5+3+w1RouCount] =
                2*(Cxx*Cyx - Cxy*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+Rct])
                + 2*(Cxx*Czx - Cxy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+1+Rct])
                + 2*(Cyx*Czx - Cyy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+2+Rct])
                + (0.5)*NormF_2*(Cxx*Cxx - Cyx*Cyx - Cxy*Cxy + Cyy*Cyy)*
                (DbAmino.Atomdata[i].CoefRou[Di*5+3+Rct])
                -1*Rs*
                (Cxx*Cxx + Cyx*Cyx - 2*Czx*Czx - Cxy*Cxy - Cyy*Cyy + 2*Czy*Czy)*
                NormF_3*(DbAmino.Atomdata[i].CoefRou[Di*5+4+Rct]);

            w1Rou[Di*5+3+w1RouCount] = w1Rou[Di*5+3+w1RouCount] / NormF_2 ;

            //----- for d(3z^2 - r^2) -----;
            w1Rou[Di*5+4+w1RouCount] =
                (4*Cxz*Cyz - 2*Cxx*Cyx - 2*Cxy*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+Rct])
                + (4*Cxz*Czz - 2*Cxx*Czx - 2*Cxy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+1+Rct])
                + (4*Cyz*Czz - 2*Cyx*Czx - 2*Cyy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefRou[Di*5+2+Rct])
                + (0.5)*
                (2*Cxz*Cxz - 2*Cyz*Cyz - Cxx*Cxx + Cyx*Cyx - Cxy*Cxy + Cyy*Cyy)*
                NormF_2*(DbAmino.Atomdata[i].CoefRou[Di*5+3+Rct])
                -1*Rs*
                (2*Cxz*Cxz + 2*Cyz*Cyz - 4*Czz*Czz - Cxx*Cxx - Cyx*Cyx + 2*Czx*Czx
                 -Cxy*Cxy - Cyy*Cyy + 2*Czy*Czy)*NormF_3*
                (DbAmino.Atomdata[i].CoefRou[Di*5+4+Rct]);

            w1Rou[Di*5+4+w1RouCount] = w1Rou[Di*5+4+w1RouCount] /  NormF_3 ;

        }
        w1RouCount += (Dorbnum*5) ;
        Rct += (Dorbnum*5);

        //////////////////////////////////////////////////////////////////////////////
        //***********************************************************
        // Myuについての変換。
        //***********************************************************
        //-------------------------------------------------------
        //-------------------------------------------------------
        //   ここから実際の軌道s,p,d（ρ、μ、νの展開係数）を変換する。
        //-------------------------------------------------------
        //-------------------------------------------------------
        Snum   = DbAmino.Atomdata[i].myuDbShellNum[0] ;
        Pnum   = DbAmino.Atomdata[i].myuDbShellNum[1] ;
        Dnum   = DbAmino.Atomdata[i].myuDbShellNum[2] ;
        SPDnum = DbAmino.Atomdata[i].myuDbShellNum[3] ;

#if VER1_0
        Sorbnum = Snum + SPDnum ;
#endif

#if VER1_1
        Sorbnum = Snum + SPDnum + SPDnum;
#endif

        Porbnum = Pnum + SPDnum ;
        Dorbnum = Dnum + SPDnum ;
        //-------------------------------------------------------
        //  S 軌道は変換なしに、そのままにする。   SSSSSSSSSSSSS
        //-------------------------------------------------------
        for (Si=0; Si<Sorbnum; Si++) {
            w1Myu[Si+w1MyuCount] = DbAmino.Atomdata[i].CoefMyu[Si] ;
        }
        w1MyuCount += Sorbnum ;
        Mct += Sorbnum;
        //-------------------------------------------------------
        //  P 軌道の変換を行なう。                 PPPPPPPPPPPPP
        //                                          +-           -+
        //                                          | Cxx Cxy Cxz |
        //  [ Cpx' Cpy' Cpz' ] = [ Cpx Cpy Cpz ] x  | Cyx Cyy Cyz |
        //                                  | Czx Czy Czz |
        //                                  +-           -+
        //-------------------------------------------------------
        for (Pi=0; Pi<Porbnum; Pi++) {
            // ----- for Px -----;
            w1Myu[Pi*3+w1MyuCount]  =
                DbAmino.Atomdata[i].CoefMyu[Pi*3+Mct]  *Cxx
                + DbAmino.Atomdata[i].CoefMyu[Pi*3+1+Mct]*Cyx
                + DbAmino.Atomdata[i].CoefMyu[Pi*3+2+Mct]*Czx;

            // ----- for Py -----;
            w1Myu[Pi*3+1+w1MyuCount]=
                DbAmino.Atomdata[i].CoefMyu[Pi*3+Mct]  *Cxy
                + DbAmino.Atomdata[i].CoefMyu[Pi*3+1+Mct]*Cyy
                + DbAmino.Atomdata[i].CoefMyu[Pi*3+2+Mct]*Czy;

            // ----- for Pz -----;
            w1Myu[Pi*3+2+w1MyuCount]=
                DbAmino.Atomdata[i].CoefMyu[Pi*3+Mct]  *Cxz
                + DbAmino.Atomdata[i].CoefMyu[Pi*3+1+Mct]*Cyz
                + DbAmino.Atomdata[i].CoefMyu[Pi*3+2+Mct]*Czz;

        }
        w1MyuCount += (Porbnum*3);
        Mct += (Porbnum*3) ;
        //-------------------------------------------------------
        //  D 軌道の変換を行なう。                 DDDDDDDDDDDDD
        //
        //  ρと同様。ただし、規格化定数に注意。
        //
        //-------------------------------------------------------
        for (Di=0; Di<Dorbnum; Di++) {
            NormNum  = FTX.getcgtonum(w1MyuCount);
            NormFact = FGX.getNormalized(NormNum,0,1,1,0);
            NormF_2  = NormFact * 0.5 ;                    // NormFact/2 ;
            NormF_3  = NormFact * InvRoot3 * 0.5 ;         // NormFact/(2*sqrt(3));

            NormFact = 1.0 / NormFact ;
            NormF_2  = 1.0 / NormF_2  ;
            NormF_3  = 1.0 / NormF_3  ;

            //----- for d(xy) -----;
            w1Myu[Di*5+w1MyuCount]   =
                (Cxx*Cyy + Cxy*Cyx)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+Mct])
                + (Cxx*Czy + Cxy*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+1+Mct])
                + (Cyx*Czy + Cyy*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+2+Mct])
                + (0.5)*NormF_2*
                (Cxx*Cxy - Cyx*Cyy)*(DbAmino.Atomdata[i].CoefMyu[Di*5+3+Mct])
                -1*Rs*(Cxx*Cxy + Cyx*Cyy - 2*Czx*Czy)*NormF_3*
                (DbAmino.Atomdata[i].CoefMyu[Di*5+4+Mct]) ;

            w1Myu[Di*5+w1MyuCount]   = w1Myu[Di*5+w1MyuCount] / NormFact ;

            //----- for d(xz) -----;
            w1Myu[Di*5+1+w1MyuCount] =
                (Cxx*Cyz + Cxz*Cyx)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+Mct])
                + (Cxx*Czz + Cxz*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+1+Mct])
                + (Cyx*Czz + Cyz*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+2+Mct])
                + (0.5)*NormF_2*
                (Cxx*Cxz - Cyx*Cyz)*(DbAmino.Atomdata[i].CoefMyu[Di*5+3+Mct])
                -1*Rs*(Cxx*Cxz + Cyx*Cyz - 2*Czx*Czz)*NormF_3*
                (DbAmino.Atomdata[i].CoefMyu[Di*5+4+Mct]);

            w1Myu[Di*5+1+w1MyuCount] = w1Myu[Di*5+1+w1MyuCount] / NormFact ;


            //----- for d(yz) -----;
            w1Myu[Di*5+2+w1MyuCount] =
                (Cxy*Cyz + Cxz*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+Mct])
                + (Cxy*Czz + Cxz*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+1+Mct])
                + (Cyy*Czz + Cyz*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+2+Mct])
                + (0.5)*NormF_2*
                (Cxy*Cxz - Cyy*Cyz)*(DbAmino.Atomdata[i].CoefMyu[Di*5+3+Mct])
                -1*Rs*(Cxy*Cxz + Cyy*Cyz - 2*Czy*Czz)*NormF_3*
                (DbAmino.Atomdata[i].CoefMyu[Di*5+4+Mct]);

            w1Myu[Di*5+2+w1MyuCount] = w1Myu[Di*5+2+w1MyuCount] / NormFact ;


            //----- for d(x^2 - y^2) -----;
            w1Myu[Di*5+3+w1MyuCount] =
                2*(Cxx*Cyx - Cxy*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+Mct])
                + 2*(Cxx*Czx - Cxy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+1+Mct])
                + 2*(Cyx*Czx - Cyy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+2+Mct])
                + (0.5)*NormF_2*(Cxx*Cxx - Cyx*Cyx - Cxy*Cxy + Cyy*Cyy)*
                (DbAmino.Atomdata[i].CoefMyu[Di*5+3+Mct])
                -1*Rs*
                (Cxx*Cxx + Cyx*Cyx - 2*Czx*Czx - Cxy*Cxy - Cyy*Cyy + 2*Czy*Czy)*
                NormF_3*(DbAmino.Atomdata[i].CoefMyu[Di*5+4+Mct]);

            w1Myu[Di*5+3+w1MyuCount] =  w1Myu[Di*5+3+w1MyuCount] / NormF_2 ;


            //----- for d(3z^2 - r^2) -----;
            w1Myu[Di*5+4+w1MyuCount] =
                (4*Cxz*Cyz - 2*Cxx*Cyx - 2*Cxy*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+Mct])
                + (4*Cxz*Czz - 2*Cxx*Czx - 2*Cxy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+1+Mct])
                + (4*Cyz*Czz - 2*Cyx*Czx - 2*Cyy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefMyu[Di*5+2+Mct])
                + (0.5)*
                (2*Cxz*Cxz - 2*Cyz*Cyz - Cxx*Cxx + Cyx*Cyx - Cxy*Cxy + Cyy*Cyy)*
                NormF_2*(DbAmino.Atomdata[i].CoefMyu[Di*5+3+Mct])
                -1*Rs*
                (2*Cxz*Cxz + 2*Cyz*Cyz - 4*Czz*Czz - Cxx*Cxx - Cyx*Cyx + 2*Czx*Czx
                 -Cxy*Cxy - Cyy*Cyy + 2*Czy*Czy)*NormF_3*
                (DbAmino.Atomdata[i].CoefMyu[Di*5+4+Mct]);

            w1Myu[Di*5+4+w1MyuCount] = w1Myu[Di*5+4+w1MyuCount] / NormF_3 ;
        }
        w1MyuCount += (Dorbnum*5) ;
        Mct += (Dorbnum*5);
        //////////////////////////////////////////////////////////////////////////////
        //***********************************************************
        // Nyuについての変換。
        //***********************************************************
        //-------------------------------------------------------
        //-------------------------------------------------------
        //   ここから実際の軌道s,p,d（ρ、μ、νの展開係数）を変換する。
        //-------------------------------------------------------
        //-------------------------------------------------------
        Snum   = DbAmino.Atomdata[i].nyuDbShellNum[0] ;
        Pnum   = DbAmino.Atomdata[i].nyuDbShellNum[1] ;
        Dnum   = DbAmino.Atomdata[i].nyuDbShellNum[2] ;
        SPDnum = DbAmino.Atomdata[i].nyuDbShellNum[3] ;

#if VER1_0
        Sorbnum = Snum + SPDnum ;
#endif

#if VER1_1
        Sorbnum = Snum + SPDnum + SPDnum;
#endif

        Porbnum = Pnum + SPDnum ;
        Dorbnum = Dnum + SPDnum ;
        //-------------------------------------------------------
        //  S 軌道は変換なしに、そのままにする。   SSSSSSSSSSSSS
        //-------------------------------------------------------
        for (Si=0; Si<Sorbnum; Si++) {
            w1Nyu[Si+w1NyuCount] = DbAmino.Atomdata[i].CoefNyu[Si] ;
        }
        w1NyuCount += Sorbnum ;
        Nct += Sorbnum;
        //-------------------------------------------------------
        //  P 軌道の変換を行なう。                 PPPPPPPPPPPPP
        //                                          +-           -+
        //                                          | Cxx Cxy Cxz |
        //  [ Cpx' Cpy' Cpz' ] = [ Cpx Cpy Cpz ] x  | Cyx Cyy Cyz |
        //                                  | Czx Czy Czz |
        //                                  +-           -+
        //-------------------------------------------------------
        for (Pi=0; Pi<Porbnum; Pi++) {
            // ----- for Px -----;
            w1Nyu[Pi*3+w1NyuCount]  =
                DbAmino.Atomdata[i].CoefNyu[Pi*3+Nct]  *Cxx
                + DbAmino.Atomdata[i].CoefNyu[Pi*3+1+Nct]*Cyx
                + DbAmino.Atomdata[i].CoefNyu[Pi*3+2+Nct]*Czx;

            // ----- for Py -----;
            w1Nyu[Pi*3+1+w1NyuCount]=
                DbAmino.Atomdata[i].CoefNyu[Pi*3+Nct]  *Cxy
                + DbAmino.Atomdata[i].CoefNyu[Pi*3+1+Nct]*Cyy
                + DbAmino.Atomdata[i].CoefNyu[Pi*3+2+Nct]*Czy;

            // ----- for Pz -----;
            w1Nyu[Pi*3+2+w1NyuCount]=
                DbAmino.Atomdata[i].CoefNyu[Pi*3+Nct]  *Cxz
                + DbAmino.Atomdata[i].CoefNyu[Pi*3+1+Nct]*Cyz
                + DbAmino.Atomdata[i].CoefNyu[Pi*3+2+Nct]*Czz;

        }
        w1NyuCount += (Porbnum*3);
        Nct += (Porbnum*3) ;
        //-------------------------------------------------------
        //  D 軌道の変換を行なう。                 DDDDDDDDDDDDD
        //
        //  μと同様。規格化定数はμのものと同じ。
        //
        //-------------------------------------------------------
        for (Di=0; Di<Dorbnum; Di++) {
            NormNum  = FTX2.getcgtonum(w1NyuCount);
            NormFact = FGX2.getNormalized(NormNum,0,1,1,0);
            NormF_2  = NormFact * 0.5 ;                    // NormFact/2 ;
            NormF_3  = NormFact * InvRoot3 * 0.5 ;         // NormFact/(2*sqrt(3));


            NormFact = 1.0 / NormFact ;
            NormF_2  = 1.0 / NormF_2  ;
            NormF_3  = 1.0 / NormF_3  ;

            //----- for d(xy) -----;
            w1Nyu[Di*5+w1NyuCount]   =
                (Cxx*Cyy + Cxy*Cyx)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+Nct])
                + (Cxx*Czy + Cxy*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+1+Nct])
                + (Cyx*Czy + Cyy*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+2+Nct])
                + (0.5)*NormF_2*
                (Cxx*Cxy - Cyx*Cyy)*(DbAmino.Atomdata[i].CoefNyu[Di*5+3+Nct])
                -1*Rs*(Cxx*Cxy + Cyx*Cyy - 2*Czx*Czy)*NormF_3*
                (DbAmino.Atomdata[i].CoefNyu[Di*5+4+Nct]) ;

            w1Nyu[Di*5+w1NyuCount]  = w1Nyu[Di*5+w1NyuCount] / NormFact ;

            //----- for d(xz) -----;
            w1Nyu[Di*5+1+w1NyuCount] =
                (Cxx*Cyz + Cxz*Cyx)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+Nct])
                + (Cxx*Czz + Cxz*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+1+Nct])
                + (Cyx*Czz + Cyz*Czx)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+2+Nct])
                + (0.5)*NormF_2*
                (Cxx*Cxz - Cyx*Cyz)*(DbAmino.Atomdata[i].CoefNyu[Di*5+3+Nct])
                -1*Rs*(Cxx*Cxz + Cyx*Cyz - 2*Czx*Czz)*NormF_3*
                (DbAmino.Atomdata[i].CoefNyu[Di*5+4+Nct]);

            w1Nyu[Di*5+1+w1NyuCount] = w1Nyu[Di*5+1+w1NyuCount] / NormFact ;

            //----- for d(yz) -----;
            w1Nyu[Di*5+2+w1NyuCount] =
                (Cxy*Cyz + Cxz*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+Nct])
                + (Cxy*Czz + Cxz*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+1+Nct])
                + (Cyy*Czz + Cyz*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+2+Nct])
                + (0.5)*NormF_2*
                (Cxy*Cxz - Cyy*Cyz)*(DbAmino.Atomdata[i].CoefNyu[Di*5+3+Nct])
                -1*Rs*(Cxy*Cxz + Cyy*Cyz - 2*Czy*Czz)*NormF_3*
                (DbAmino.Atomdata[i].CoefNyu[Di*5+4+Nct]);

            w1Nyu[Di*5+2+w1NyuCount] = w1Nyu[Di*5+2+w1NyuCount] / NormFact ;


            //----- for d(x^2 - y^2) -----;
            w1Nyu[Di*5+3+w1NyuCount] =
                2*(Cxx*Cyx - Cxy*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+Nct])
                + 2*(Cxx*Czx - Cxy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+1+Nct])
                + 2*(Cyx*Czx - Cyy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+2+Nct])
                + (0.5)*NormF_2*(Cxx*Cxx - Cyx*Cyx - Cxy*Cxy + Cyy*Cyy)*
                (DbAmino.Atomdata[i].CoefNyu[Di*5+3+Nct])
                -1*Rs*
                (Cxx*Cxx + Cyx*Cyx - 2*Czx*Czx - Cxy*Cxy - Cyy*Cyy + 2*Czy*Czy)*
                NormF_3*(DbAmino.Atomdata[i].CoefNyu[Di*5+4+Nct]);

            w1Nyu[Di*5+3+w1NyuCount] = w1Nyu[Di*5+3+w1NyuCount] / NormF_2 ;


            //----- for d(3z^2 - r^2) -----;
            w1Nyu[Di*5+4+w1NyuCount] =
                (4*Cxz*Cyz - 2*Cxx*Cyx - 2*Cxy*Cyy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+Nct])
                + (4*Cxz*Czz - 2*Cxx*Czx - 2*Cxy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+1+Nct])
                + (4*Cyz*Czz - 2*Cyx*Czx - 2*Cyy*Czy)*
                (NormFact*DbAmino.Atomdata[i].CoefNyu[Di*5+2+Nct])
                + (0.5)*
                (2*Cxz*Cxz - 2*Cyz*Cyz - Cxx*Cxx + Cyx*Cyx - Cxy*Cxy + Cyy*Cyy)*
                NormF_2*(DbAmino.Atomdata[i].CoefNyu[Di*5+3+Nct])
                -1*Rs*
                (2*Cxz*Cxz + 2*Cyz*Cyz - 4*Czz*Czz - Cxx*Cxx - Cyx*Cyx + 2*Czx*Czx
                 -Cxy*Cxy - Cyy*Cyy + 2*Czy*Czy)*NormF_3*
                (DbAmino.Atomdata[i].CoefNyu[Di*5+4+Nct]);

            w1Nyu[Di*5+4+w1NyuCount] = w1Nyu[Di*5+4+w1NyuCount] / NormF_3 ;
        }
        w1NyuCount += (Dorbnum*5) ;
        Nct += (Dorbnum*5);

        // InpAmino.Order is -1 no type no Translation kokode owari.

        goto  before_iLoopEnd;

        //############//
        //
Translation2:;   //      InpAmino.Order is -1 .
        //
        //############//

        //###############################################################
        //
        //  Rou Rou Rou Rou Rou Rou Rou Rou Rou Rou Rou Rou Rou Rou Rou
        //
        //###############################################################
        Snum   = DbAmino.Atomdata[i].rouDbShellNum[0] ;
        Pnum   = DbAmino.Atomdata[i].rouDbShellNum[1] ;
        Dnum   = DbAmino.Atomdata[i].rouDbShellNum[2] ;
        SPDnum = DbAmino.Atomdata[i].rouDbShellNum[3] ;

#if VER1_0
        Sorbnum = Snum + SPDnum ;
#endif

#if VER1_1
        Sorbnum = Snum + SPDnum + SPDnum;
#endif

        Porbnum = Pnum + SPDnum ;
        Dorbnum = Dnum + SPDnum ;
        //-------------------------------------------------------
        //  for S orbital
        //-------------------------------------------------------
        for (Si=0; Si<Sorbnum; Si++) {
            w1Rou[Si+w1RouCount] = DbAmino.Atomdata[i].CoefRou[Si] ;
        }
        w1RouCount += Sorbnum ;
        Rct += Sorbnum;
        //-------------------------------------------------------
        //  for P orbital
        //-------------------------------------------------------
        for (Pi=0; Pi<Porbnum; Pi++) {
            // ----- for Px -----;
            w1Rou[Pi*3+w1RouCount]  = 0.0 ;
            // ----- for Py -----;
            w1Rou[Pi*3+1+w1RouCount]= 0.0 ;
            // ----- for Pz -----;
            w1Rou[Pi*3+2+w1RouCount]= 0.0 ;
        }
        w1RouCount += (Porbnum*3);
        Rct += (Porbnum*3) ;
        //-------------------------------------------------------
        //  for D orbital
        //-------------------------------------------------------
        for (Di=0; Di<Dorbnum; Di++) {
            //----- for d(xy) -----;
            w1Rou[Di*5+w1RouCount]   = 0.0 ;
            //----- for d(xz) -----;
            w1Rou[Di*5+1+w1RouCount] = 0.0 ;
            //----- for d(yz) -----;
            w1Rou[Di*5+2+w1RouCount] = 0.0 ;
            //----- for d(x^2 - y^2) -----;
            w1Rou[Di*5+3+w1RouCount] = 0.0 ;
            //----- for d(3z^2 - r^2) -----;
            w1Rou[Di*5+4+w1RouCount] = 0.0 ;
        }
        w1RouCount += (Dorbnum*5) ;
        Rct += (Dorbnum*5);

        //###############################################################
        //
        //  Myu Myu Myu Myu Myu Myu Myu Myu Myu Myu Myu Myu Myu Myu Myu
        //
        //###############################################################
        Snum   = DbAmino.Atomdata[i].myuDbShellNum[0] ;
        Pnum   = DbAmino.Atomdata[i].myuDbShellNum[1] ;
        Dnum   = DbAmino.Atomdata[i].myuDbShellNum[2] ;
        SPDnum = DbAmino.Atomdata[i].myuDbShellNum[3] ;

#if VER1_0
        Sorbnum = Snum + SPDnum ;
#endif

#if VER1_1
        Sorbnum = Snum + SPDnum + SPDnum;
#endif

        Porbnum = Pnum + SPDnum ;
        Dorbnum = Dnum + SPDnum ;
        //-------------------------------------------------------
        //  for S orbital
        //-------------------------------------------------------
        for (Si=0; Si<Sorbnum; Si++) {
            w1Myu[Si+w1MyuCount] = DbAmino.Atomdata[i].CoefMyu[Si] ;
        }
        w1MyuCount += Sorbnum ;
        Mct += Sorbnum;
        //-------------------------------------------------------
        //  for P orbital
        //-------------------------------------------------------
        for (Pi=0; Pi<Porbnum; Pi++) {
            // ----- for Px -----;
            w1Myu[Pi*3+w1MyuCount]  = 0.0 ;
            // ----- for Py -----;
            w1Myu[Pi*3+1+w1MyuCount]= 0.0 ;
            // ----- for Pz -----;
            w1Myu[Pi*3+2+w1MyuCount]= 0.0 ;
        }
        w1MyuCount += (Porbnum*3);
        Mct += (Porbnum*3) ;
        //-------------------------------------------------------
        //  for D orbital
        //-------------------------------------------------------
        for (Di=0; Di<Dorbnum; Di++) {
            //----- for d(xy) -----;
            w1Myu[Di*5+w1MyuCount]   = 0.0 ;
            //----- for d(xz) -----;
            w1Myu[Di*5+1+w1MyuCount] = 0.0 ;
            //----- for d(yz) -----;
            w1Myu[Di*5+2+w1MyuCount] = 0.0 ;
            //----- for d(x^2 - y^2) -----;
            w1Myu[Di*5+3+w1MyuCount] = 0.0 ;
            //----- for d(3z^2 - r^2) -----;
            w1Myu[Di*5+4+w1MyuCount] = 0.0 ;
        }
        w1MyuCount += (Dorbnum*5) ;
        Mct += (Dorbnum*5);

        //###############################################################
        //
        //  Nyu Nyu Nyu Nyu Nyu Nyu Nyu Nyu Nyu Nyu Nyu Nyu Nyu Nyu Nyu
        //
        //###############################################################
        Snum   = DbAmino.Atomdata[i].nyuDbShellNum[0] ;
        Pnum   = DbAmino.Atomdata[i].nyuDbShellNum[1] ;
        Dnum   = DbAmino.Atomdata[i].nyuDbShellNum[2] ;
        SPDnum = DbAmino.Atomdata[i].nyuDbShellNum[3] ;

#if VER1_0
        Sorbnum = Snum + SPDnum ;
#endif

#if VER1_1
        Sorbnum = Snum + SPDnum + SPDnum;
#endif

        Porbnum = Pnum + SPDnum ;
        Dorbnum = Dnum + SPDnum ;
        //-------------------------------------------------------
        //  for S orbital
        //-------------------------------------------------------
        for (Si=0; Si<Sorbnum; Si++) {
            w1Nyu[Si+w1NyuCount] = DbAmino.Atomdata[i].CoefNyu[Si] ;
        }
        w1NyuCount += Sorbnum ;
        Nct += Sorbnum;
        //-------------------------------------------------------
        //  for P orbital
        //-------------------------------------------------------

        for (Pi=0; Pi<Porbnum; Pi++) {
            // ----- for Px -----;
            w1Nyu[Pi*3+w1NyuCount]  = 0.0 ;
            // ----- for Py -----;
            w1Nyu[Pi*3+1+w1NyuCount]= 0.0 ;
            // ----- for Pz -----;
            w1Nyu[Pi*3+2+w1NyuCount]= 0.0 ;
        }
        w1NyuCount += (Porbnum*3);
        Nct += (Porbnum*3) ;
        //-------------------------------------------------------
        //  for D orbital
        //-------------------------------------------------------
        for (Di=0; Di<Dorbnum; Di++) {
            //----- for d(xy) -----;
            w1Nyu[Di*5+w1NyuCount]   = 0.0 ;
            //----- for d(xz) -----;
            w1Nyu[Di*5+1+w1NyuCount] = 0.0 ;
            //----- for d(yz) -----;
            w1Nyu[Di*5+2+w1NyuCount] = 0.0 ;
            //----- for d(x^2 - y^2) -----;
            w1Nyu[Di*5+3+w1NyuCount] = 0.0 ;
            //----- for d(3z^2 - r^2) -----;
            w1Nyu[Di*5+4+w1NyuCount] = 0.0 ;
        }
        w1NyuCount += (Dorbnum*5) ;
        Nct += (Dorbnum*5);

        // InpAmino.Order is -1 no type no Translation kokode owari.

        goto  before_iLoopEnd;

        //#############################################################

before_iLoopEnd:
        ;

    }// loop_i ;

//   if( OutLevel < -3 ){
//     log<<"\n";
//     log<<"######################################################## "<<"\n";
//     log<<"###############   TransOrbital  Exit  ################## "<<"\n";
//     log<<"######################################################## "<<"\n";
//     log<<"\n";
//   }
    return 0;

}
//#############################################################################
void DfInitialguess::CoverVector()
{
    TlLogX& log = TlLogX::getInstance();

//   if( OutLevel < -3 ){
//     log<<"\n";
//     log<<"### Cover Vector(User defined)  START### "<<"\n";
//     log<<"### and Copy Vector(CoefRou,CoefMyu,CoefNyu) ###"<<"\n";
//   }

    //int startAtom,finalAtom,endaux,startaux,auxlength;
    //int i,j;

    // for Rou
    for (int i=0; i < UserVctRouNum; i++) {
        const int startAtom = UserVctRou[i].From1 ;
        const int finalAtom = UserVctRou[i].To1;
        int endaux;
        if (finalAtom != AtomNum) {
            endaux  =this->EachAD[finalAtom+1].Rou_StartNum ;
        } else {
            endaux = MaxTermRou;
        }
        const int startaux = this->EachAD[startAtom].Rou_StartNum;
        const int auxlength = endaux - startaux ;
//     if (OutLevel < -5){
//       log<<"DEBUG auxlength(rou) = "<<auxlength<<"\n";
//     }

        {
            //Fl_GuessVector  FGV(UserVctRou[i].Filename);
            TlVector FGV;
            FGV.load(UserVctRou[i].Filename);

            //if (auxlength != FGV.getAuxterm()){
            if (static_cast<TlVector::size_type>(auxlength) != FGV.getSize()) {
                log << "Bad defined VectorFile ["
                << UserVctRou[i].Filename << "]" << "\n";
                CnErr.abort();
            }
            for (int j=0; j < auxlength; j++) {
                //w1Rou[j+startaux] = FGV.getCoef(j) ;
                w1Rou[j+startaux] = FGV[j];
            }
        }
    }

    if (tempflag==0) {
        // for Myu
        for (int i=0; i < UserVctMyuNum; i++) {
            const int startAtom = UserVctMyu[i].From1 ;
            const int finalAtom = UserVctMyu[i].To1;
            int endaux;
            if (finalAtom != AtomNum) {
                endaux   = this->EachAD[finalAtom+1].Myu_StartNum ;
            } else {
                endaux = MaxTermMyu;
            }
            const int startaux = this->EachAD[startAtom].Myu_StartNum;
            const int auxlength = endaux - startaux ;
//       if (OutLevel < -5){
//  log << "DEBUG auxlength(myu) = " << auxlength << "\n";
//       }

            {
                //Fl_GuessVector FGV(UserVctMyu[i].Filename);
                TlVector FGV;
                FGV.load(UserVctMyu[i].Filename);
                //if (auxlength != FGV.getAuxterm()){
                if (static_cast<TlVector::size_type>(auxlength) != FGV.getSize()) {
                    log << "Bad defined VectorFile ["
                    << UserVctMyu[i].Filename << "]" << "\n";
                    CnErr.abort();
                }
                for (int j=0; j < auxlength; j++) {
                    //w1Myu[j+startaux] = FGV.getCoef(j) ;
                    w1Myu[j+startaux] = FGV[j];
                }
            }
        }

        // for Nyu
        for (int i=0; i < UserVctNyuNum; i++) {
            const int startAtom = UserVctNyu[i].From1 ;
            const int finalAtom = UserVctNyu[i].To1;
            int endaux;
            if (finalAtom != AtomNum) {
                endaux   = this->EachAD[finalAtom+1].Nyu_StartNum ;
            } else {
                endaux = MaxTermNyu ;
            }
            const int startaux = this->EachAD[startAtom].Nyu_StartNum;
            const int auxlength = endaux - startaux ;
//       if( OutLevel < -5 ){
//  log<<"DEBUG auxlength(nyu) = "<<auxlength<<"\n";
//       }

            {
                //Fl_GuessVector  FGV(UserVctNyu[i].Filename);
                TlVector FGV;
                FGV.load(UserVctNyu[i].Filename);
                //if (auxlength != FGV.getAuxterm()){
                if (static_cast<TlVector::size_type>(auxlength) != FGV.getSize()) {
                    log << "Bad defined VectorFile ["
                    << UserVctNyu[i].Filename << "]" << "\n";
                    CnErr.abort();
                }
                for (int j=0; j < auxlength; j++) {
                    //w1Nyu[j+startaux] = FGV.getCoef(j) ;
                    w1Nyu[j+startaux] = FGV[j];
                }
            }
        }
    }

    if (scftype == NSP) {
        // Copy vector : w1Rou --> CoefRou
//     for(i=0;i<MaxTermRou;i++){
//       CoefRou[i] = w1Rou[i];
//     }
        this->CoefRou = this->w1Rou;

        if (tempflag==0) {
            // Copy vector : w1Myu --> CoefMyu
//       for(i=0;i<MaxTermMyu;i++){
//  CoefMyu[i] = w1Myu[i];
//       }
            this->CoefMyu = this->w1Myu;

            // Copy vector : w1Nyu --> CoefNyu
//       for(i=0;i<MaxTermNyu;i++){
//  CoefNyu[i] = w1Nyu[i];
//       }
            this->CoefNyu = this->w1Nyu;
        }

        if (tempflag==0) {
            if (MethodMyuNyu == Meth1) {
//  if( OutLevel < -5 ){
//    log<<"This routine make Xcpotential."<<"\n";
//    log<<"The Method is  (Myu) = Pow( Rou , 1/3 ) ."<<"\n";
//    log<<"attention :  MaxTermRou and MaxTermMyu is equal."<<"\n";
//  }
                if (MaxTermRou != MaxTermMyu) {
                    log<<"This routine cannot use this calculation"<<"\n";
                    log<<"Because , Total Auxiliary Number for Rou is "
                    <<"differrent from Total Auxiliary Number for Myu."<<"\n";
                    log<<"So, you must appoint other method to make Xcpotential"
                    <<"\n";
                    CnErr.abort();
                }

//  for(i=0;i<MaxTermRou;i++){
//    CoefMyu[i] = -1* CoefRou[i];
//    CoefNyu[i] = CoefRou[i];
//  }
                this->CoefMyu = -1.0 * this->CoefRou;
                this->CoefNyu = this->CoefRou;
            }//meth1

            if (MethodMyuNyu == Meth3) {
                log<<"Bad appointment"<<"\n";
                log<<"MethodMyuNyu meth3 is SP calculation option."<<"\n";
                log<<"You must change option's keyword make-myu-nyu "<<"\n";
                log<<"Keyword is [ meth0 or meth1 ]."<<"\n";
                CnErr.abort();
            }
        }
    } else if (scftype == SP) {
        // Copy vector : w1Rou --> CoefRou{Alpha,Beta}
//     for(i=0;i<MaxTermRou;i++){
//       CoefRouAlpha[i] = w1Rou[i];
//       CoefRouBeta[i]  = w1Rou[i];
//     }
        this->CoefRouAlpha = this->w1Rou;
        this->CoefRouBeta  = this->w1Rou;

        if (tempflag==0) {
            // Copy vector : w1Myu --> CoefMyu{Alpha,Beta}
//       for(i=0;i<MaxTermMyu;i++){
//  CoefMyuAlpha[i] = w1Myu[i];
//  CoefMyuBeta[i]  = w1Myu[i];
//       }
            this->CoefMyuAlpha = this->w1Myu;

            // Copy vector : w1Nyu --> CoefNyu{Alpha,Beta}
//       for(i=0;i<MaxTermNyu;i++){
//  CoefNyuAlpha[i] = w1Nyu[i];
//  CoefNyuBeta[i]  = w1Nyu[i];
//       }
            this->CoefNyuAlpha = this->w1Nyu;
            this->CoefNyuBeta  = this->w1Nyu;
        }

        // option for make XCpotential.

        if (tempflag==0) {
            if (MethodMyuNyu == Meth1) {
//  if( OutLevel < -5 ){
//    log<<"This routine make Xcpotential."<<"\n";
//    log<<"The Method is  (Myu) = Pow( Rou , 1/3 ) ."<<"\n";
//    log<<"attention :  MaxTermRou and MaxTermMyu is equal."<<"\n";
//  }
                if (MaxTermRou != MaxTermMyu) {
                    log<<"This routine cannot use this calculation"<<"\n";
                    log<<"Because , Total Auxiliary Number for Rou is "
                    <<"differrent from Total Auxiliary Number for Myu."<<"\n";
                    log<<"So, you must appoint other method to make Xcpotential"
                    <<"\n";
                    //        log<<"\7"<<"\n";
                    CnErr.abort();
                }

//  for(i=0;i<MaxTermRou;i++){
//    CoefMyuAlpha[i] = -1* CoefRouAlpha[i];
//    CoefMyuBeta[i]  = -1* CoefRouBeta[i];
//    CoefNyuAlpha[i] = CoefRouAlpha[i];
//    CoefNyuBeta[i]  = CoefRouBeta[i];
//  }
                this->CoefMyuAlpha = -1.0 * this->CoefRouAlpha;
                this->CoefMyuBeta  = -1.0 * this->CoefRouBeta;
                this->CoefNyuAlpha = this->CoefRouAlpha;
                this->CoefNyuBeta  = this->CoefRouBeta;
            }

            if (MethodMyuNyu == Meth3) {
                if (MaxTermRou != MaxTermMyu) {
                    log<<"This routine cannot use this calculation"<<"\n";
                    log<<"Because , Total Auxiliary Number for Rou is "
                    <<"differrent from Total Auxiliary Number for Myu."<<"\n";
                    log<<"So, you must appoint other method to make Xcpotential"
                    <<"\n";
                    CnErr.abort();
                }

                {
                    //      int Mnum,Nnum;

//    Fl_Vct_Myu  fvm(999);
//    Fl_Vct_Nyu  fvn(999);

//    fvm.open("fl_Work",fvm.getfilename(),"read");
//    fvm.getelemnum(&Mnum);
//    fvm.read(Mnum,CoefMyu);
//    fvm.close();
                    this->CoefMyu.load("fl_Work/fl_Vct_Myu999");

//    fvn.open("fl_Work",fvn.getfilename(),"read");
//    fvn.getelemnum(&Nnum);
//    fvn.read(Nnum,CoefNyu);
//    fvn.close();
                    this->CoefNyu.load("fl_Work/fl_Vct_Nyu999");

                    /* // for debug
                       log<<"in gusCoverVector +++ Meth3 debug"<<"\n";
                       for(i=0;i<Mnum;i++){
                       log<<"CoefMyu["<<i<<"] = "<<CoefMyu[i]<<"\n";
                       }
                       for(i=0;i<Nnum;i++){
                       log<<"CoefNyu["<<i<<"] = "<<CoefNyu[i]<<"\n";
                       }
                    */

                    assert(this->CoefMyu.getSize() == this->CoefNyu.getSize());
                    assert(static_cast<TlVector::size_type>(MaxTermRou) == this->CoefMyu.getSize());

//    if( (Mnum!=Nnum) || (MaxTermRou!=Mnum) ){
//      log<<"This routine cannot use this calculation"<<"\n";
//      log<<"Because , Total Auxiliary Number for Rou is "
//         <<"differrent from Total Auxiliary Number for Myu."<<"\n";
//      log<<"You must make fl_Vct_Myu999,fl_Vct_Nyu999"<<"\n";
//      CnErr.abort();
//    }

//    for(i=0;i<Mnum;i++){
//      CoefMyuAlpha[i] = CoefMyu[i]*0.5;
//      CoefMyuBeta[i]  = CoefMyu[i]*0.5;
//      CoefNyuAlpha[i] = CoefNyu[i]*0.5;
//      CoefNyuBeta[i]  = CoefNyu[i]*0.5;
//    }
                    this->CoefMyuAlpha = this->CoefMyu * 0.5;
                    this->CoefMyuBeta  = this->CoefMyu * 0.5;
                    this->CoefNyuAlpha = this->CoefNyu * 0.5;
                    this->CoefNyuBeta  = this->CoefNyu * 0.5;
                }
            }
        }
    }

//   if( OutLevel < -8 ){

//     // ------ DEBUG PRINT -----------
//     if(scftype==NSP){
//       for (int i=0;i<MaxTermRou;i++){
//  log<<"CoefRou["<<i<<"] = "<<CoefRou[i]<<"\n";
//       }

//       if (tempflag==0){
//  for (int i=0;i<MaxTermMyu;i++){
//    log<<"CoefMyu["<<i<<"] = "<<CoefMyu[i]<<"\n";
//  }

//  for (int i=0;i<MaxTermNyu;i++){
//    log<<"CoefNyu["<<i<<"] = "<<CoefNyu[i]<<"\n";
//  }
//       }
//     }

//     if (scftype == SP){
//       for (int i=0;i<MaxTermRou;i++){
//  log<<"CoefRouAlpha["<<i<<"] = "<<CoefRouAlpha[i]<<"\n";
//       }
//       for (int i=0;i<MaxTermRou;i++){
//  log<<"CoefRouBeta["<<i<<"] = "<<CoefRouBeta[i]<<"\n";
//       }

//       if (tempflag==0){
//  for (int i=0;i<MaxTermMyu;i++){
//    log<<"CoefMyuAlpha["<<i<<"] = "<<CoefMyuAlpha[i]<<"\n";
//  }
//  for (int i=0;i<MaxTermMyu;i++){
//    log<<"CoefMyuBeta["<<i<<"] = "<<CoefMyuBeta[i]<<"\n";
//  }
//  for (int i=0;i<MaxTermNyu;i++){
//    log<<"CoefNyuAlpha["<<i<<"] = "<<CoefNyuAlpha[i]<<"\n";
//  }
//  for (int i=0;i<MaxTermNyu;i++){
//    log<<"CoefNyuBeta["<<i<<"] = "<<CoefNyuBeta[i]<<"\n";
//  }
//       }
//     }
//   }

//   if( OutLevel < -3 ){
//     log<<"### Cover Vector(User defined)  end.### "<<"\n"<<"\n";
//   }
}
//#############################################################################
int DfInitialguess::VectorNormalyze()
{
    TlLogX& log = TlLogX::getInstance();

    if (VctRouNormalize==OFF) {
        if (PartialRouNormNum!=0) {
            log<<"Warning : VctRouNormalize is OFF(=Not Normalize)"<<"\n";
            log<<"        : But , PartialRouNormNum isnot Zero "<<"\n";
            log<<"        : PartialRouNormNum = "<<PartialRouNormNum<<"\n";
        }
    }
    if (VctMyuNormalize==OFF) {
        if (PartialMyuNormNum!=0) {
            log<<"Warning : VctMyuNormalize is OFF(=Not Normalize)"<<"\n";
            log<<"        : But , PartialMyuNormNum isnot Zero "<<"\n";
            log<<"        : PartialMyuNormNum = "<<PartialMyuNormNum<<"\n";
        }
    }
    if (VctNyuNormalize==OFF) {
        if (PartialNyuNormNum!=0) {
            log<<"Warning : VctNyuNormalize is OFF(=Not Normalize)"<<"\n";
            log<<"        : But , PartialNyuNormNum isnot Zero "<<"\n";
            log<<"        : PartialNyuNormNum = "<<PartialNyuNormNum<<"\n";
        }
    }

    if ((VctRouNormalize==OFF) && (VctMyuNormalize==OFF) &&
            (VctNyuNormalize==OFF)) {
        goto vctNormalizeEND;
    }

    //---------------------------------------------------
//   if( OutLevel < 0 ){
//     log<<"  ----- Normalyze rootine start ----- "<<"\n";
//   }

    if (VctRouNormalize==ON) {
        VctRouNorm();
    }
    if (VctMyuNormalize==ON) {
        VctMyuNorm();
    }
    if (VctNyuNormalize==ON) {
        VctNyuNorm();
    }

//   if( OutLevel < 0 ){
//     log<<"\n";
//     log<<"  ----- Normalyze rootine exit ----- "<<"\n";
//   }

vctNormalizeEND:
    ;

    return 0 ;

}
//#############################################################################
void DfInitialguess::VctRouNorm()
{
    TlLogX& log = TlLogX::getInstance();

    int    i,j,k,mcnt,hikisuu,from,to;
    double  myuCC,Pai,N_one3rd,N_two3rd,kata_alpha,Sum,PreSum,InvmyuCC;
    double  FlVctNalpha,invPai;
    //char Atom[20],Lab2[20];
    std::string Atom, Lab2;
    double defElec[MaxPartialNum];
    double  defElecTotalNum; // defined atom total electron num.
    int x,f2,t2;
    double Dens;

    Fl_Geometry FG(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Density FGD;

    Pai = 3.14159265358979;
    invPai = 1.0/Pai ;
    myuCC  = -(1.5)*nAlpha*pow((3.0/Pai) , (1.0/3.0));
    InvmyuCC = 1.0 / myuCC ;
    //myuCC = (-3/2)*nAlpha(3/Pai)^(1/3)
    //myu ---> CoefMyu[i] * pow( (2*Pai/kata_gamma) , (1.0/3.0) ) / myuCC ;
    N_one3rd = pow(ElectronNum , (1.0/3.0));
    N_two3rd = pow(ElectronNum , (2.0/3.0));

//   FG.open("fl_Geometry","read");
//   FG.read();
//   FG.close();
//   FG.load();

//   FGD.open("fl_Gto_Density","read");
//   FGD.read();
//   FGD.close();



//   if( OutLevel < 0 ){
//     log<<"  ----- Rou Normalyze rootine start ----- "<<"\n";
//   }

    //-----------------------------------------------------
    // Normalyze for CoefRou[i]  <----- Coefficient Rou
    //-----------------------------------------------------

    defElecTotalNum=0.0;
    for (i=0; i<PartialRouNormNum; i++) {
        defElecTotalNum += PartialRouNormalize[i].ElecNum;
    }

    //-----------------------------------------------------
    // Count ElectronNumber for Non-defined-Atom.
    //-----------------------------------------------------
    Sum = PreSum = 0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Rou_StartNum;
        //strcpy(Atom, FG.getAtom(i).c_str());
        //strcpy(Lab2, FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGD.getStartposition(j);
            // if ((strcmp(Atom, FGD.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2, FGD.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGD.getAtom(hikisuu)) &&
                (Lab2 == FGD.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGD.getTermnumber(j);
                //log<<"debug From = "<<from<<" to = "<<to<<"\n";
                for (k=from; k<to; k++) {
                    if (FGD.getShell(k)=='s') {
                        if (rouNormOnOff[i]==OFF) {
                            kata_alpha = FGD.getExponent(k,0);
                            FlVctNalpha = pow((Pai/(2*kata_alpha)) , 0.25) ;
                            Sum += (FlVctNalpha * CoefRou[mcnt]);
                        }
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;

    //-----------------------------------------------
    // Normalize for Non-defined-Atom.
    //-----------------------------------------------
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Rou_StartNum;
        //strcpy(Atom,FG.getAtom(i).c_str());
        //strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGD.getStartposition(j);
            // if ((strcmp(Atom,FGD.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGD.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGD.getAtom(hikisuu)) &&
                (Lab2 == FGD.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGD.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGD.getShell(k)=='s') {
                        if (rouNormOnOff[i]==OFF) {
                            CoefRou[mcnt]*=((ElectronNum-defElecTotalNum)/Sum);
                        }
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;


    //-----------------------------------------------------
    // Count ElectronNumber for defined-Atom.
    //-----------------------------------------------------
    for (x=0; x<PartialRouNormNum; x++) {
        defElec[x] = 0.0 ;
        f2 = PartialRouNormalize[x].From1;
        t2 = PartialRouNormalize[x].To1;
        for (i=f2; i<=t2; i++) {
            mcnt = this->EachAD[i].Rou_StartNum;
            //strcpy(Atom,FG.getAtom(i).c_str());
            //strcpy(Lab2,FG.getLabel(i).c_str());
            Atom = FG.getAtom(i);
            Lab2 = FG.getLabel(i);
            to = 0;
            for (j=0; j<AtomKindNum; j++) {
                hikisuu = FGD.getStartposition(j);
                // if ((strcmp(Atom,FGD.getAtom(hikisuu).c_str())==0) &&
                //         (strcmp(Lab2,FGD.getLabel(hikisuu).c_str())==0)) {
                if ((Atom == FGD.getAtom(hikisuu)) &&
                    (Lab2 == FGD.getLabel(hikisuu))) {
                    from = hikisuu ;
                    to   = from + FGD.getTermnumber(j);
                    for (k=from; k<to; k++) {
                        if (FGD.getShell(k)=='s') {
                            if (rouNormOnOff[i]==ON) {
                                kata_alpha = FGD.getExponent(k,0);
                                FlVctNalpha = pow((Pai/(2*kata_alpha)),0.25);
                                defElec[x] += (FlVctNalpha * CoefRou[mcnt]);
                            }
                            mcnt++;
                        }
                    } //k;
                } //if;
            } //j;
        } //i;
    }    // x;

    //-----------------------------------------------
    // Normalize for defined-Atom.
    //-----------------------------------------------
    for (x=0; x<PartialRouNormNum; x++) {
        f2 = PartialRouNormalize[x].From1;
        t2 = PartialRouNormalize[x].To1;
        for (i=f2; i<=t2; i++) {
            mcnt = this->EachAD[i].Rou_StartNum;
            //strcpy(Atom,FG.getAtom(i).c_str());
            //strcpy(Lab2,FG.getLabel(i).c_str());
            Atom = FG.getAtom(i);
            Lab2 = FG.getLabel(i);
            to = 0;
            for (j=0; j<AtomKindNum; j++) {
                hikisuu = FGD.getStartposition(j);
                // if ((strcmp(Atom,FGD.getAtom(hikisuu).c_str())==0) &&
                //         (strcmp(Lab2,FGD.getLabel(hikisuu).c_str())==0)) {
                if ((Atom == FGD.getAtom(hikisuu)) &&
                    (Lab2 == FGD.getLabel(hikisuu))) {
                    from = hikisuu ;
                    to   = from + FGD.getTermnumber(j);
                    for (k=from; k<to; k++) {
                        if (FGD.getShell(k)=='s') {
                            if (rouNormOnOff[i]==ON) {
                                Dens = PartialRouNormalize[x].ElecNum;
                                CoefRou[mcnt] *= (Dens / defElec[x]);
                            }
                            mcnt++;
                        }
                    } //k;
                } //if;
            } //j;
        } //i;
    }    // x;


    ///////////////////////////////////////////
    //-----------------------------------------
    // check electron number
    //-----------------------------------------
    ///////////////////////////////////////////

    Sum=0.0;
    PreSum=0.0;

    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Rou_StartNum;
        //strcpy(Atom,FG.getAtom(i).c_str());
        //strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGD.getStartposition(j);
            // if ((strcmp(Atom,FGD.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGD.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGD.getAtom(hikisuu)) &&
                (Lab2 == FGD.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGD.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGD.getShell(k)=='s') {
                        kata_alpha = FGD.getExponent(k,0);
                        Sum += pow((Pai/(2*kata_alpha)),0.25)*CoefRou[mcnt];
                        mcnt++;
                    } else if (FGD.getShell(k)=='p') {
                        mcnt=mcnt+3;
                    } else if (FGD.getShell(k)=='d') {
                        mcnt=mcnt+5;
                    }
                } //k;
            } //if;
        } //j;
//     if( OutLevel < -3 ){
//       log<<"Electron  atom["<<i<<"] =  "
//   << Atom<<" =  "<< Sum - PreSum <<"\n";
//     }
        PreSum = Sum ;
    } //i;

    // modified and commented out by AS(koike.s) 2003/05/23
    //    if( OutLevel < -2 ){
    log<<"Total  CoefRou Electron Num = [" << Sum <<" ]"<<"\n";
    log<<"Real ElectronNum = [ "<<ElectronNum<<" ]"<<"\n";

    std::cout <<"Total CoefRou Electron Num = [" << Sum <<" ]"<<"\n";
    std::cout <<"Real ElectronNum = [ "<<ElectronNum<<" ]"<<"\n" ;
    //    }

    if (pow((Sum-ElectronNum) , 2) > 1E-6) {
        log<<"BAD Normalize (Rou) "<<"\n";
        log<<"ElectronNumber is different from ElecNumber after";
        log<<" Normalized Coefficient"<<"\n";
        //  log<<"\7"<<"\n";

        log << "#### CoefRou ####\n";
        for (i=0; i<MaxTermRou; i++) {
            log<<"CoefRou [ "<<i<<" ]   = "<<CoefRou[i]<<"\n";
        }

        CnErr.abort();
    }

    /*
    // for debug
    for(i=0;i<MaxTermRou;i++){
    log<<"CoefRou [ "<<i<<" ]   = "<<CoefRou[i]<<"\n";
    }
    */

//   if( OutLevel < 0 ){
//     log<<"\n";
//     log<<"  ----- Rou Normalyze rootine exit ----- "<<"\n";
//   }

}
//#############################################################################
void DfInitialguess::VctMyuNorm()
{
    TlLogX& log = TlLogX::getInstance();

    int    i,j,k,mcnt,hikisuu,from,to;
    double  myuCC,Pai,N_one3rd,N_two3rd,kata_alpha,Sum,PreSum,InvmyuCC;
    double  invPai; // double FlVctNalpha;
    //char Atom[20],Lab2[20];
    std::string Atom, Lab2;
    double defElec[MaxPartialNum];
    double defElecTotalNum;
    int x,f2,t2;
    double Dens;
    //int Rnum;

    Fl_Geometry FG(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Xcpot FGX;

    Pai = 3.14159265358979;
    invPai = 1.0/Pai ;
    myuCC  = -(1.5)*nAlpha*pow((3.0/Pai) , (1.0/3.0));
    InvmyuCC = 1.0 / myuCC ;
    //myuCC = (-3/2)*nAlpha(3/Pai)^(1/3)
    //myu ---> CoefMyu[i] * pow( (2*Pai/kata_gamma) , (1.0/3.0) ) / myuCC ;
    N_one3rd = pow(ElectronNum , (1.0/3.0));
    N_two3rd = pow(ElectronNum , (2.0/3.0));

//   FG.load();
//   FG.open("fl_Geometry","read");
//   FG.read();
//   FG.close();

//   FGX.open("fl_Gto_Xcpot","read");
//   FGX.read();
//   FGX.close();

//   if( OutLevel < 0 ){
//     log<<"  ----- Myu Normalyze rootine start ----- "<<"\n";
//   }
    //-----------------------------------------------------
    // Normalyze for CoefMyu[i]  <----- Coefficient Myu
    //-----------------------------------------------------
    defElecTotalNum = 0.0;
    for (i=0; i<PartialMyuNormNum; i++) {
        defElecTotalNum += PartialMyuNormalize[i].ElecNum;
    }

    //-----------------------------------------------------
    // Count Electron Number for Non-defined-Atom.
    //-----------------------------------------------------
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Myu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX.getStartposition(j);
            // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX.getAtom(hikisuu)) &&
                (Lab2 == FGX.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX.getShell(k)=='s') {
                        if (myuNormOnOff[i]==OFF) {
                            kata_alpha = FGX.getExponent(k,0);
                            Sum = Sum + CoefMyu[mcnt] *
                                  pow((2*Pai/kata_alpha) , 0.75) * InvmyuCC ;
                        }
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
        //log<<"Sum = "<<Sum<<"\n";
    } //i;
    //log<<"Total Sum  = "<<Sum<<"\n";

    //-----------------------------------------------------
    // Normalize for Non-defined-Atom.
    //-----------------------------------------------------
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Myu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX.getStartposition(j);
            // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX.getAtom(hikisuu)) &&
                (Lab2 == FGX.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX.getShell(k)=='s') {
                        if (myuNormOnOff[i]==OFF) {
                            CoefMyu[mcnt] *= ((N_one3rd-defElecTotalNum)/Sum);
                        }
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;


    //-----------------------------------------------------
    // Count Electron Number for defined-Atom.
    //-----------------------------------------------------
    for (x=0; x<PartialMyuNormNum; x++) {
        defElec[x]=0.0;
        f2 = PartialMyuNormalize[x].From1;
        t2 = PartialMyuNormalize[x].To1;
        for (i=f2; i<=t2; i++) {
            mcnt = this->EachAD[i].Myu_StartNum;
            // strcpy(Atom,FG.getAtom(i).c_str());
            // strcpy(Lab2,FG.getLabel(i).c_str());
            Atom = FG.getAtom(i);
            Lab2 = FG.getLabel(i);
            to = 0;
            for (j=0; j<AtomKindNum; j++) {
                hikisuu = FGX.getStartposition(j);
                // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
                //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
                if ((Atom == FGX.getAtom(hikisuu)) &&
                    (Lab2 == FGX.getLabel(hikisuu))) {
                    from = hikisuu ;
                    to   = from + FGX.getTermnumber(j);
                    for (k=from; k<to; k++) {
                        if (FGX.getShell(k)=='s') {
                            if (myuNormOnOff[i]==ON) {
                                kata_alpha = FGX.getExponent(k,0);
                                defElec[x] = defElec[x] + CoefMyu[mcnt] *
                                             pow((2*Pai/kata_alpha) , 0.75) * InvmyuCC;
                            }
                            mcnt++;
                        }
                    } //k;
                } //if;
            } //j;
        } //i;
    } //x;

    //-----------------------------------------------------
    // Normalize for Non-defined-Atom.
    //-----------------------------------------------------
    for (x=0; x<PartialMyuNormNum; x++) {
        f2 = PartialMyuNormalize[x].From1;
        t2 = PartialMyuNormalize[x].To1;
        for (i=f2; i<=t2; i++) {
            mcnt = this->EachAD[i].Myu_StartNum;
            // strcpy(Atom,FG.getAtom(i).c_str());
            // strcpy(Lab2,FG.getLabel(i).c_str());
            Atom = FG.getAtom(i);
            Lab2 = FG.getLabel(i);
            to = 0;
            for (j=0; j<AtomKindNum; j++) {
                hikisuu = FGX.getStartposition(j);
                // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
                //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
                if ((Atom == FGX.getAtom(hikisuu)) &&
                    (Lab2 == FGX.getLabel(hikisuu))) {
                    from = hikisuu ;
                    to   = from + FGX.getTermnumber(j);
                    for (k=from; k<to; k++) {
                        if (FGX.getShell(k)=='s') {
                            if (myuNormOnOff[i]==ON) {
                                Dens = PartialMyuNormalize[x].ElecNum;
                                CoefMyu[mcnt] *= (Dens / defElec[x]);
                            }
                            mcnt++;
                        }
                    } //k;
                } //if;
            } //j;
        } //i;
    } //x;


    ///////////////////////////////////////////
    //-----------------------------------------
    // check electron number
    //-----------------------------------------
    ///////////////////////////////////////////
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Myu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX.getStartposition(j);
            // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX.getAtom(hikisuu)) &&
                (Lab2 == FGX.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX.getShell(k)=='s') {
                        kata_alpha = FGX.getExponent(k,0);
                        Sum = Sum + CoefMyu[mcnt] *
                              pow((2*Pai/kata_alpha) , 0.75) * InvmyuCC ;
                        mcnt++;
                    } else if (FGX.getShell(k)=='p') {
                        mcnt += 3 ;
                    } else if (FGX.getShell(k)=='d') {
                        mcnt += 5 ;
                    }
                } //k;
            } //if;
        } //j;
//     if( OutLevel < -3 ){
//       log<<"Electron  atom["<<i<<"] =  "
//   << Atom<<" =  "<< Sum - PreSum <<"\n";
//     }
        PreSum = Sum ;
    } //i;

//   if( OutLevel < -2 ){
//     log<<"Total  CoefMyu Electron Num = [" << Sum <<" ]"<<"\n";
//     log<<"Rou no 1/3 jyou = [ "<< N_one3rd <<" ] "<<"\n";
//   }

    if (pow((Sum-N_one3rd) , 2) > 1E-6) {
        log<<"BAD Normalize (Myu) "<<"\n";
        log<<"ElectronNumber^(1/3) is different from ElecNumber^(1/3) after";
        log<<" Normalized Coefficient"<<"\n";
        //  log<<"\7"<<"\n";
        CnErr.abort();
    }


//   if( OutLevel < -9 ){
//     // for debug
//     for(i=0;i<MaxTermMyu;i++){
//       log<<"CoefMyu [ "<<i<<" ]   = "<<CoefMyu[i]<<"\n";
//     }
//   }

//   if( OutLevel < 0 ){
//     log<<"\n";
//     log<<"  ----- Myu Normalyze rootine exit ----- "<<"\n";
//   }

}
//#############################################################################
void DfInitialguess::VctNyuNorm()
{
    TlLogX& log = TlLogX::getInstance();

    int    i,j,k,mcnt,hikisuu,from,to;
    double  myuCC,Pai,N_one3rd,N_two3rd,kata_alpha,Sum,PreSum,InvmyuCC;
    double  invPai,Tb2,dumy; // FlVctNalpha
    //char Atom[20],Lab2[20];
    std::string Atom, Lab2;
    double defElec[MaxPartialNum];
    double defElecTotalNum;
    int x,f2,t2;
    double Dens;
    // int Rnum;

    Fl_Geometry FG(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Xcpot2 FGX2;

    Pai = 3.14159265358979;
    invPai = 1.0/Pai ;
    myuCC  = -(1.5)*nAlpha*pow((3.0/Pai) , (1.0/3.0));
    InvmyuCC = 1.0 / myuCC ;
    //myuCC = (-3/2)*nAlpha(3/Pai)^(1/3)
    //myu ---> CoefNyu[i] * pow( (2*Pai/kata_gamma) , (1.0/3.0) ) / myuCC ;
    N_one3rd = pow(ElectronNum , (1.0/3.0));
    Tb2 = 2.0/3.0;
    N_two3rd = pow(ElectronNum , (2.0/3.0));

//   FG.load();
//   FG.open("fl_Geometry","read");
//   FG.read();
//   FG.close();

    ///// modified by Fric(mizouchi) 2001/12/09 (s) /////
    //    FGX2.open("fl_Gto_Xcpot2","read");
//   FGX2.open("fl_Gto_Xcpot","read");
    ///// modified by Fric(mizouchi) 2001/12/09 (e) /////

//   FGX2.read();
//   FGX2.close();

//   if( OutLevel < 0 ){
//     log<<"  ----- Nyu Normalyze rootine start ----- "<<"\n";
//   }
    //-----------------------------------------------------
    // Normalyze for CoefNyu[i]  <----- Coefficient Nyu
    //-----------------------------------------------------
    defElecTotalNum = 0.0;
    for (i=0; i<PartialNyuNormNum; i++) {
        defElecTotalNum += PartialNyuNormalize[i].ElecNum;
    }

    //-----------------------------------------------------
    // Count Electron Number for Non-defined-Atom.
    //-----------------------------------------------------
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Nyu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX2.getStartposition(j);
            // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX2.getAtom(hikisuu)) &&
                (Lab2 == FGX2.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX2.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX2.getShell(k)=='s') {
                        if (nyuNormOnOff[i]==OFF) {
                            kata_alpha = FGX2.getExponent(k,0);
                            Sum = Sum + CoefNyu[mcnt] *
                                  pow((2*Pai/kata_alpha) , 0.75) ;
                        }
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;

//   if( OutLevel < -4 ){
//     log<<"DEBUG  Pre sum = "<<Sum<<"\n";
//   }

    //-----------------------------------------------------
    // Normalize for Non-defined-Atom.
    //-----------------------------------------------------
    dumy = pow(ElectronNum,(2./3.)) - defElecTotalNum ;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Nyu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX2.getStartposition(j);
            // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX2.getAtom(hikisuu)) &&
                (Lab2 == FGX2.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX2.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX2.getShell(k)=='s') {
                        if (nyuNormOnOff[i]==OFF) {
                            CoefNyu[mcnt] *= (dumy / Sum);
                        }
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;


    //-----------------------------------------------------
    // Count Electron Number for defined-Atom.
    //-----------------------------------------------------
    for (x=0; x<PartialNyuNormNum; x++) {
        defElec[x]=0.0;
        f2 = PartialNyuNormalize[x].From1;
        t2 = PartialNyuNormalize[x].To1;
        for (i=f2; i<=t2; i++) {
            mcnt = this->EachAD[i].Nyu_StartNum;
            // strcpy(Atom,FG.getAtom(i).c_str());
            // strcpy(Lab2,FG.getLabel(i).c_str());
            Atom = FG.getAtom(i);
            Lab2 = FG.getLabel(i);
            to = 0;
            for (j=0; j<AtomKindNum; j++) {
                hikisuu = FGX2.getStartposition(j);
                // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
                //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
                if ((Atom == FGX2.getAtom(hikisuu)) &&
                    (Lab2 == FGX2.getLabel(hikisuu))) {
                    from = hikisuu ;
                    to   = from + FGX2.getTermnumber(j);
                    for (k=from; k<to; k++) {
                        if (FGX2.getShell(k)=='s') {
                            if (nyuNormOnOff[i]==ON) {
                                kata_alpha = FGX2.getExponent(k,0);
                                defElec[x] = defElec[x] + CoefNyu[mcnt] *
                                             pow((2*Pai/kata_alpha) , 0.75) ;
                            }
                            mcnt++;
                        }
                    } //k;
                } //if;
            } //j;
        } //i;
    } //x;

    //-----------------------------------------------------
    // Normalize for defined-Atom.
    //-----------------------------------------------------
    for (x=0; x<PartialNyuNormNum; x++) {
        f2 = PartialNyuNormalize[x].From1;
        t2 = PartialNyuNormalize[x].To1;
        for (i=f2; i<=t2; i++) {
            mcnt = this->EachAD[i].Nyu_StartNum;
            // strcpy(Atom,FG.getAtom(i).c_str());
            // strcpy(Lab2,FG.getLabel(i).c_str());
            Atom = FG.getAtom(i);
            Lab2 = FG.getLabel(i);
            to = 0;
            for (j=0; j<AtomKindNum; j++) {
                hikisuu = FGX2.getStartposition(j);
                // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
                //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
                if ((Atom == FGX2.getAtom(hikisuu)) &&
                    (Lab2 == FGX2.getLabel(hikisuu))) {
                    from = hikisuu ;
                    to   = from + FGX2.getTermnumber(j);
                    for (k=from; k<to; k++) {
                        if (FGX2.getShell(k)=='s') {
                            if (nyuNormOnOff[i]==ON) {
                                Dens=PartialNyuNormalize[x].ElecNum;
                                CoefNyu[mcnt] *= (Dens / defElec[x]);
                            }
                            mcnt++;
                        }
                    } //k;
                } //if;
            } //j;
        } //i;
    } //x;

    ///////////////////////////////////////////
    //-----------------------------------------
    // check electron number
    //-----------------------------------------
    ///////////////////////////////////////////
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Nyu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX2.getStartposition(j);
            // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX2.getAtom(hikisuu)) &&
                (Lab2 == FGX2.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX2.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX2.getShell(k)=='s') {
                        kata_alpha = FGX2.getExponent(k,0);
                        Sum = Sum + CoefNyu[mcnt] *
                              pow((2*Pai/kata_alpha) , 0.75) ;
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
//     if( OutLevel < -3 ){
//       log<<"Electron  atom["<<i<<"] =  "
//   << Atom<<" =  "<< Sum - PreSum <<"\n";
//     }
        PreSum = Sum ;
    } //i;

//   if( OutLevel < -2 ){
//     log<<"Total  CoefNyu Electron Num = [" << Sum <<" ]"<<"\n";
//     log<<"Rou no 2/3 jyou = [ "<< N_two3rd <<" ]"<<"\n";
//   }

    if (pow((Sum-N_two3rd) , 2) > 1E-6) {
        log<<"BAD Normalize (Nyu) "<<"\n";
        log<<"ElectronNumber^(2/3) is different from ElecNumber^(2/3) after";
        log<<" Normalized Coefficient"<<"\n";
        //  log<<"\7"<<"\n";
        CnErr.abort();
    }

//   if( OutLevel < -9 ){
//     // for debug
//     for(i=0;i<MaxTermNyu;i++){
//       log<<"CoefNyu [ "<<i<<" ]   = "<<CoefNyu[i]<<"\n";
//     }
//   }

//   if( OutLevel < 0 ){
//     log<<"\n";
//     log<<"  ----- Nyu Normalyze rootine exit ----- "<<"\n";
//   }

}
//#############################################################################
int DfInitialguess::VectorNormalyze2(void)
{
//   if( OutLevel < 0 ){
//     log<<"  ----- Normalyze rootine start ----- "<<"\n";
//   }

    if (VctRouNormalize==ON) {
        VctRouNorm2();
    }
    if (VctMyuNormalize==ON) {
        VctMyuNorm2();
    }
    if (VctNyuNormalize==ON) {
        VctNyuNorm2();
    }

//   if( OutLevel < 0 ){
//     log<<"  ----- Normalyze rootine exit ----- "<<"\n";
//   }

    return 0 ;
}
//#############################################################################
void DfInitialguess::VctRouNorm2(void)
{

    int    i,j,k,mcnt,hikisuu,from,to;
    double  myuCC,Pai,N_one3rd,N_two3rd,kata_alpha,Sum,PreSum,InvmyuCC;
    double  invPai; // double FlVctNalpha
    //char Atom[20],Lab2[20];
    std::string Atom, Lab2;
    //int x,f2,t2;
    //double Dens;

    Fl_Geometry FG(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Density FGD;

    Pai = 3.14159265358979;
    invPai = 1.0/Pai ;
    myuCC  = -(1.5)*nAlpha*pow((3.0/Pai) , (1.0/3.0));
    InvmyuCC = 1.0 / myuCC ;
    //myuCC = (-3/2)*nAlpha(3/Pai)^(1/3)
    //myu ---> CoefMyu[i] * pow( (2*Pai/kata_gamma) , (1.0/3.0) ) / myuCC ;
    N_one3rd = pow(ElectronNum , (1.0/3.0));
    N_two3rd = pow(ElectronNum , (2.0/3.0));

//   FG.load();
//   FG.open("fl_Geometry","read");
//   FG.read();
//   FG.close();

//   FGD.open("fl_Gto_Density","read");
//   FGD.read();
//   FGD.close();


//   if( OutLevel < 0 ){
//     log<<"  ----- Rou Normalyze rootine start ----- "<<"\n";
//   }

    ////////////////////////////////////////////
    //  Count Number of alpha Density
    ////////////////////////////////////////////
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Rou_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGD.getStartposition(j);
            // if ((strcmp(Atom,FGD.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGD.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGD.getAtom(hikisuu)) &&
                (Lab2 == FGD.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGD.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGD.getShell(k)=='s') {
                        kata_alpha = FGD.getExponent(k,0);
                        Sum+=pow((Pai/(2*kata_alpha)),0.25)*CoefRouAlpha[mcnt];
                        mcnt++;
                    } else if (FGD.getShell(k)=='p') {
                        mcnt=mcnt+3;
                    } else if (FGD.getShell(k)=='d') {
                        mcnt=mcnt+5;
                    }
                } //k;
            } //if;
        } //j;
    } //i;

    ////////////////////////////////////////////
    //  Normalyze of alpha Density
    ////////////////////////////////////////////
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Rou_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGD.getStartposition(j);
            // if ((strcmp(Atom,FGD.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGD.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGD.getAtom(hikisuu)) &&
                (Lab2 == FGD.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGD.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGD.getShell(k)=='s') {
                        CoefRouAlpha[mcnt] *= (AlphaSpinNum / Sum);
                        mcnt++;
                    } else if (FGD.getShell(k)=='p') {
                        mcnt=mcnt+3;
                    } else if (FGD.getShell(k)=='d') {
                        mcnt=mcnt+5;
                    }
                } //k;
            } //if;
        } //j;
    } //i;


    //###################################################################

    ////////////////////////////////////////////
    //  Count Number of Beta Density
    ////////////////////////////////////////////
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Rou_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGD.getStartposition(j);
            // if ((strcmp(Atom,FGD.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGD.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGD.getAtom(hikisuu)) &&
                (Lab2 == FGD.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGD.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGD.getShell(k)=='s') {
                        kata_alpha = FGD.getExponent(k,0);
                        Sum+=pow((Pai/(2*kata_alpha)),0.25)*CoefRouBeta[mcnt];
                        mcnt++;
                    } else if (FGD.getShell(k)=='p') {
                        mcnt=mcnt+3;
                    } else if (FGD.getShell(k)=='d') {
                        mcnt=mcnt+5;
                    }
                } //k;
            } //if;
        } //j;
    } //i;

    ////////////////////////////////////////////
    //  Normalyze of Beta Density
    ////////////////////////////////////////////
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Rou_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGD.getStartposition(j);
            // if ((strcmp(Atom,FGD.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGD.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGD.getAtom(hikisuu)) &&
                (Lab2 == FGD.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGD.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGD.getShell(k)=='s') {
                        CoefRouBeta[mcnt] *= (BetaSpinNum / Sum);
                        mcnt++;
                    } else if (FGD.getShell(k)=='p') {
                        mcnt=mcnt+3;
                    } else if (FGD.getShell(k)=='d') {
                        mcnt=mcnt+5;
                    }
                } //k;
            } //if;
        } //j;
    } //i;

//   if( OutLevel < 0 ){
//     log<<"  ----- Rou Normalyze rootine end. ----- "<<"\n";
//   }
}
//#############################################################################
void DfInitialguess::VctMyuNorm2(void)
{

    int    i,j,k,mcnt,hikisuu,from,to;
    double  myuCC,Pai,N_one3rd,N_two3rd,kata_alpha,Sum,PreSum,InvmyuCC;
    double  aN_one3rd,bN_one3rd;
    double  invPai; // double FlVctNalpha;
    //char Atom[20],Lab2[20];
    std::string Atom, Lab2;
    // int x,f2,t2;
    //double Dens;

    Fl_Geometry FG(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Xcpot FGX;

    Pai = 3.14159265358979;
    invPai = 1.0/Pai ;
    myuCC  = -(1.5)*nAlpha*pow((3.0/Pai) , (1.0/3.0));
    InvmyuCC = 1.0 / myuCC ;
    //myuCC = (-3/2)*nAlpha(3/Pai)^(1/3)
    //myu ---> CoefMyu[i] * pow( (2*Pai/kata_gamma) , (1.0/3.0) ) / myuCC ;
    N_one3rd = pow(ElectronNum , (1.0/3.0));
    N_two3rd = pow(ElectronNum , (2.0/3.0));
    aN_one3rd = pow(AlphaSpinNum , (1.0/3.0));   // for alpha
    bN_one3rd = pow(BetaSpinNum  , (1.0/3.0));   // for beta

//   FG.load();
//   FG.open("fl_Geometry","read");
//   FG.read();
//   FG.close();

//   FGX.open("fl_Gto_Xcpot","read");
//   FGX.read();
//   FGX.close();

    //-----------------------------------------------------
    // Count Electron Number for Alpha Xcpotential
    //-----------------------------------------------------
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Myu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX.getStartposition(j);
            // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX.getAtom(hikisuu)) &&
                (Lab2 == FGX.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX.getShell(k)=='s') {
                        kata_alpha = FGX.getExponent(k,0);
                        Sum = Sum + CoefMyuAlpha[mcnt] *
                              pow((2*Pai/kata_alpha) , 0.75) * InvmyuCC ;
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;

    //-----------------------------------------------------
    // Normalize for Non-defined-Atom.
    //-----------------------------------------------------
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Myu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX.getStartposition(j);
            // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX.getAtom(hikisuu)) &&
                (Lab2 == FGX.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX.getShell(k)=='s') {
                        if (myuNormOnOff[i]==OFF) {
                            CoefMyuAlpha[mcnt] *= (aN_one3rd/Sum);
                        }
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;



    //==================================================================

    //-----------------------------------------------------
    // Count Electron Number for Beta Xcpotential
    //-----------------------------------------------------
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Myu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX.getStartposition(j);
            // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX.getAtom(hikisuu)) &&
                (Lab2 == FGX.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX.getShell(k)=='s') {
                        kata_alpha = FGX.getExponent(k,0);
                        Sum = Sum + CoefMyuBeta[mcnt] *
                              pow((2*Pai/kata_alpha) , 0.75) * InvmyuCC ;
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;

    //-----------------------------------------------------
    // Normalize for Non-defined-Atom.
    //-----------------------------------------------------
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Myu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX.getStartposition(j);
            // if ((strcmp(Atom,FGX.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX.getAtom(hikisuu)) &&
                (Lab2 == FGX.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX.getShell(k)=='s') {
                        if (myuNormOnOff[i]==OFF) {
                            CoefMyuBeta[mcnt] *= (bN_one3rd/Sum);
                        }
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;

//   if( OutLevel < 0 ){
//     log<<"\n";
//     log<<"  ----- Myu Normalyze rootine exit ----- "<<"\n";
//   }

}
//#############################################################################
void DfInitialguess::VctNyuNorm2(void)
{

    int    i,j,k,mcnt,hikisuu,from,to;
    double  myuCC,Pai,N_one3rd,N_two3rd,kata_alpha,Sum,PreSum,InvmyuCC;
    double  invPai,Tb2,dumy; // double FlVctNalpha;
    //char Atom[20],Lab2[20];
    std::string Atom, Lab2;
    //int x,f2,t2;
    //double Dens;

    Fl_Geometry FG(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Xcpot2 FGX2;

    Pai = 3.14159265358979;
    invPai = 1.0/Pai ;
    myuCC  = -(1.5)*nAlpha*pow((3.0/Pai) , (1.0/3.0));
    InvmyuCC = 1.0 / myuCC ;
    //myuCC = (-3/2)*nAlpha(3/Pai)^(1/3)
    //myu ---> CoefNyu[i] * pow( (2*Pai/kata_gamma) , (1.0/3.0) ) / myuCC ;
    N_one3rd = pow(ElectronNum , (1.0/3.0));
    Tb2 = 2.0/3.0;
    N_two3rd = pow(ElectronNum , (2.0/3.0));

//   FG.load();
//   FG.open("fl_Geometry","read");
//   FG.read();
//   FG.close();

    ///// modified by Fric(mizouchi) 2001/12/09 (s) /////
    //    FGX2.open("fl_Gto_Xcpot2","read");
//   FGX2.open("fl_Gto_Xcpot","read");
    ///// modified by Fric(mizouchi) 2001/12/09 (e) /////

//   FGX2.read();
//   FGX2.close();
    //-----------------------------------------------------
    // Count Electron Number for alpha
    //-----------------------------------------------------
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Nyu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX2.getStartposition(j);
            // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX2.getAtom(hikisuu)) &&
                (Lab2 == FGX2.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX2.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX2.getShell(k)=='s') {
                        kata_alpha = FGX2.getExponent(k,0);
                        Sum = Sum + CoefNyuAlpha[mcnt] *
                              pow((2*Pai/kata_alpha) , 0.75) ;
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;


    //-----------------------------------------------------
    // Normalize for alpha
    //-----------------------------------------------------
    dumy = pow(AlphaSpinNum,(2./3.));
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Nyu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX2.getStartposition(j);
            // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX2.getAtom(hikisuu)) &&
                (Lab2 == FGX2.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX2.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX2.getShell(k)=='s') {
                        CoefNyuAlpha[mcnt] *= (dumy / Sum);
                    }
                    mcnt++;
                }// k
            } // if
        } //j;
    } //i;

    //===================================================================

    //-----------------------------------------------------
    // Count Electron Number for beta
    //-----------------------------------------------------
    Sum=0.0;
    PreSum=0.0;
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Nyu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX2.getStartposition(j);
            // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX2.getAtom(hikisuu)) &&
                (Lab2 == FGX2.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX2.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX2.getShell(k)=='s') {
                        kata_alpha = FGX2.getExponent(k,0);
                        Sum = Sum + CoefNyuBeta[mcnt] *
                              pow((2*Pai/kata_alpha) , 0.75) ;
                        mcnt++;
                    }
                } //k;
            } //if;
        } //j;
    } //i;


    //-----------------------------------------------------
    // Normalize for alpha
    //-----------------------------------------------------
    dumy = pow(BetaSpinNum,(2./3.));
    for (i=0; i<AtomNum; i++) {
        mcnt = this->EachAD[i].Nyu_StartNum;
        // strcpy(Atom,FG.getAtom(i).c_str());
        // strcpy(Lab2,FG.getLabel(i).c_str());
        Atom = FG.getAtom(i);
        Lab2 = FG.getLabel(i);
        to = 0;
        for (j=0; j<AtomKindNum; j++) {
            hikisuu = FGX2.getStartposition(j);
            // if ((strcmp(Atom,FGX2.getAtom(hikisuu).c_str())==0) &&
            //         (strcmp(Lab2,FGX2.getLabel(hikisuu).c_str())==0)) {
            if ((Atom == FGX2.getAtom(hikisuu)) &&
                (Lab2 == FGX2.getLabel(hikisuu))) {
                from = hikisuu ;
                to   = from + FGX2.getTermnumber(j);
                for (k=from; k<to; k++) {
                    if (FGX2.getShell(k)=='s') {
                        CoefNyuBeta[mcnt] *= (dumy / Sum);
                    }
                    mcnt++;
                }// k
            } // if
        } //j;
    } //i;


//   if( OutLevel < 0 ){
//     log<<"\n";
//     log<<"  ----- Nyu Normalyze rootine exit ----- "<<"\n";
//   }

}

int DfInitialguess::GusOutput()
{
    if (this->scftype == NSP) {
        if (this->PpqFname.empty()) {
            //this->CoefRou.save("fl_Work/fl_Vct_Rou0");
            this->CoefRou.save("fl_Work/fl_Vct_Rou1");
        }

        if (tempflag==0) {
            //this->CoefMyu.save("fl_Work/fl_Vct_Myu0");
            this->CoefMyu.save("fl_Work/fl_Vct_Myu1");
            //this->CoefNyu.save("fl_Work/fl_Vct_Nyu0");
            this->CoefNyu.save("fl_Work/fl_Vct_Nyu1");
        }
    } else if (this->scftype == SP) {
        if ((this->AlphaPpqFname.empty()) || (this->BetaPpqFname.empty())) {
            //this->CoefRouAlpha.save("fl_Work/fl_Vct_Roua0");
            //this->CoefRouBeta.save("fl_Work/fl_Vct_Roub0");
            this->CoefRouAlpha.save("fl_Work/fl_Vct_Roua1");
            this->CoefRouBeta.save("fl_Work/fl_Vct_Roub1");
        }

        if (tempflag==0) {
            //this->CoefMyuAlpha.save("fl_Work/fl_Vct_Myua0");
            //this->CoefMyuBeta.save("fl_Work/fl_Vct_Myub0");
            this->CoefMyuAlpha.save("fl_Work/fl_Vct_Myua1");
            this->CoefMyuBeta.save("fl_Work/fl_Vct_Myub1");
            //this->CoefNyuAlpha.save("fl_Work/fl_Vct_Nyua0");
            //this->CoefNyuBeta.save("fl_Work/fl_Vct_Nyub0");
            this->CoefNyuAlpha.save("fl_Work/fl_Vct_Nyua1");
            this->CoefNyuBeta.save("fl_Work/fl_Vct_Nyub1");
        }
    }

    return 0;
}

void DfInitialguess::Print()
{
    TlLogX& log = TlLogX::getInstance();
    if (scftype == NSP) {
        TlVector rho;
        rho.load("fl_Work/fl_Vct_Rou1");

        log<<"\n";
        log<<"###########################################################"<<"\n";
        log<<"############ Initialguess Vector ( fl_Vct_Rou ) ###########"<<"\n";
        log<<"###########################################################"<<"\n";
        for (TlVector::size_type x = 0; x < rho.getSize(); x++) {
            log << "Coefficient_Rou[" << x << "]= " << rho[x] << "\n";
        }

        if (tempflag==0) {
            TlVector myu;
            myu.load("fl_Work/fl_Vct_Myu1");

            log<<"\n";
            log<<"###########################################################"<<"\n";
            log<<"############ Initialguess Vector ( fl_Vct_Myu ) ###########"<<"\n";
            log<<"###########################################################"<<"\n";
            for (TlVector::size_type x=0; x < myu.getSize(); x++) {
                log << TlUtils::format("Coefficient_Myu[%6ld] = %15.8lf\n", x, myu[x]);
            }

            TlVector nyu;
            nyu.load("fl_Work/fl_Vct_Nyu1");

            log<<"\n";
            log<<"###########################################################"<<"\n";
            log<<"############ Initialguess Vector ( fl_Vct_Nyu ) ###########"<<"\n";
            log<<"###########################################################"<<"\n";
            for (TlVector::size_type x = 0; x < nyu.getSize(); x++) {
                log << TlUtils::format("Coefficient_Nyu[%6ld] = %15.8lf\n",x, nyu[x]);
            }
        }
    } else if (scftype == SP) {
        //--- for UHF(SP) ------------
        {
            TlVector rhoa;
            rhoa.load("fl_Work/fl_Vct_Roua1");

            log<<"\n";
            log<<"###########################################################"<<"\n";
            log<<"####### Initialguess Vector ( fl_Vct_Rou_Alpha ) ##########"<<"\n";
            log<<"###########################################################"<<"\n";
            for (TlVector::size_type x = 0; x < rhoa.getSize(); x++) {
                log << TlUtils::format("Coefficient_Rou_Alpha[%6ld] = %15.8lf\n", x, rhoa[x]);
            }
        }

        {
            TlVector rhob;
            rhob.load("fl_Work/fl_Vct_Roub1");

            log<<"\n";
            log<<"###########################################################"<<"\n";
            log<<"####### Initialguess Vector ( fl_Vct_Rou_Beta ) ###########"<<"\n";
            log<<"###########################################################"<<"\n";
            for (TlVector::size_type x = 0; x < rhob.getSize(); x++) {
                log << TlUtils::format("Coefficient_Rou_Beta[%6ld] = %15.8lf\n", x, rhob[x]);
            }
        }

        if (tempflag==0) {
            {
                TlVector myua;
                myua.load("fl_Work/fl_Vct_Myua1");

                log<<"\n";
                log<<"###########################################################"<<"\n";
                log<<"####### Initialguess Vector ( fl_Vct_Myu_Alpha ) ##########"<<"\n";
                log<<"###########################################################"<<"\n";
                for (TlVector::size_type x = 0; x < myua.getSize(); x++) {
                    log << TlUtils::format("Coefficient_Myu_Alpha[%6ld] = %15.8lf\n", x, myua[x]);
                }
            }

            {
                TlVector myub;
                myub.load("fl_Work/fl_Vct_Myub1");

                log<<"\n";
                log<<"###########################################################"<<"\n";
                log<<"####### Initialguess Vector ( fl_Vct_Myu_Beta ) ###########"<<"\n";
                log<<"###########################################################"<<"\n";
                for (TlVector::size_type x=0; x < myub.getSize(); x++) {
                    log << TlUtils::format("Coefficient_Myu_Beta[%6ld] = %15.8lf\n" ,x, myub[x]);
                }
            }

            {
                TlVector nyua;
                nyua.load("fl_Work/fl_Vct_Nyua1");

                log<<"\n";
                log<<"###########################################################"<<"\n";
                log<<"####### Initialguess Vector ( fl_Vct_Nyu_Alpha ) ##########"<<"\n";
                log<<"###########################################################"<<"\n";
                for (TlVector::size_type x = 0; x < nyua.getSize(); x++) {
                    log << TlUtils::format("Coefficient_Nyu_Alpha[%6ld] = %15.8lf\n", x, nyua[x]);
                }
            }

            {
                TlVector nyub;
                nyub.load("fl_Work/fl_Vct_Nyub1");

                log<<"\n";
                log<<"###########################################################"<<"\n";
                log<<"####### Initialguess Vector ( fl_Vct_Nyu_Beta ) ###########"<<"\n";
                log<<"###########################################################"<<"\n";
                for (TlVector::size_type x=0; x < nyub.getSize(); x++) {
                    log << TlUtils::format("Coefficient_Nyu_Beta[%6ld] = %15.8lf\n", x, nyub[x]);
                }
            }
        }
        //-----
    }

}


// Essentially, this rootine is copied from DfIGuess2
// Variable name is partially modified in order to keep the
// consistency with DfInitiallguess class.
void DfInitialguess::setCoefRoudata()
{
    TlLogX& log = TlLogX::getInstance();

    //---- file open and read data ---------------------------------------------------
    tempflag = 1;

    std::ifstream  fi;
    fi.open("guess.rho", std::ios::in);
    if (!fi) {
        log << "Cannot open \"guess.rho\"" <<"\n";
        std::cout << "Cannot open \"guess.rho\"" << std::endl;
        CnErr.abort();
    }
    log << "    ./guess.rho is opened\n\n";
    log.flush();

    int term;
    fi >> term;

    this->CoefRou.resize(term);
    for (int k=0; k<term; k++) {
        fi >> this->CoefRou[k];
    }
    fi.close();

    if (scftype == SP) {
        this->CoefRouAlpha = 0.5 * this->CoefRou;
        this->CoefRouBeta  = 0.5 * this->CoefRou;
    }
}

//added by AS(koike.s) 2003/05/23
// This part is mainly moved from DfIguess2.
// The bugs in the case of SP in DfIguess2 has been modified.
void DfInitialguess::VctNormfromIGuess2()
{
    TlLogX& log = TlLogX::getInstance();

    tempflag = 1;

    TlVector coefA;
    coefA.load(this->getNalphaPath());

    if (coefA.getSize() != static_cast<TlVector::size_type>(MaxTermRou)) {
        log<<" Inconsistency detected between fl_Vct_Rou and fl_Vct_Nalpha "<<"\n";
        CnErr.abort();
    }

    if (scftype==NSP) {
        const double dumelenum = this->CoefRou * coefA;

        std::cout << "    number of electron = " << ElectronNum << "\n";
        std::cout << "    Rho * Nalpha       = " << dumelenum << "\n";
        std::cout << "    difference         = " << ElectronNum - dumelenum << "\n";
        std::cout << "    Rho is normalized by the number of electrons." << "\n";
        log  << "    number of electron = " << ElectronNum << "\n";
        log  << "    Rho * Nalpha       = " << dumelenum << "\n";
        log  << "    difference         = " << ElectronNum - dumelenum << "\n";
        log  << "    Rho is normalized by the number of electrons." << "\n";

        this->CoefRou *= (ElectronNum / dumelenum);

        // NSP case end
    } else if (scftype==SP) {
        const double dumelenumalpha = this->CoefRouAlpha * coefA;

        std::cout << "    number of alpha electron = " << AlphaSpinNum << "\n";
        std::cout << "    Rho * Nalpha       = " << dumelenumalpha << "\n";
        std::cout << "    difference         = " << AlphaSpinNum - dumelenumalpha << "\n";
        std::cout << "    Rho-alpha is normalized by the number of electrons." << "\n";
        std::cout <<"\n\n\n";
        log  << "    number of alpha electron = " << AlphaSpinNum << "\n";
        log  << "    Rho * Nalpha       = " << dumelenumalpha << "\n";
        log  << "    difference         = " << AlphaSpinNum - dumelenumalpha << "\n";
        log  << "    Rho-alpha is normalized by the number of electrons." << "\n";
        log <<"\n\n\n";

        this->CoefRouAlpha *= (AlphaSpinNum / dumelenumalpha);

        // for beta-spin
        const double dumelenumbeta = this->CoefRouBeta * coefA;

        std::cout << "    number of alpha electron = " << BetaSpinNum << "\n";
        std::cout << "    Rho * Nalpha       = " << dumelenumbeta << "\n";
        std::cout << "    difference         = " << BetaSpinNum - dumelenumbeta << "\n";
        std::cout << "    Rho-alpha is normalized by the number of electrons." << "\n";
        std::cout <<"\n\n\n";
        log  << "    number of alpha electron = " << BetaSpinNum << "\n";
        log  << "    Rho * Nalpha       = " << dumelenumbeta << "\n";
        log  << "    difference         = " << BetaSpinNum - dumelenumbeta << "\n";
        log  << "    Rho-alpha is normalized by the number of electrons." << "\n";
        log <<"\n\n\n";

        this->CoefRouBeta *= (BetaSpinNum / dumelenumbeta);
    } else {
        std::cout<< " scftype = "<<scftype<<std::endl;
        std::cout<< " DfInitialguess::VctNormfromIGuess2"<<std::endl;
        std::cout<< " scftype is somethig wrong "<<std::endl;
        CnErr.abort();
    }
}


int DfInitialguess::GusError(const std::string& str)
{
    TlLogX& log = TlLogX::getInstance();

    log << "**** Error in DfInitialguess ****\n";
    log << str << "\n";
    CnErr.abort();
    return 0;
}

void DfInitialguess::densityFitOnly()
{
    TlLogX& log = TlLogX::getInstance();

    log << "Densityfitonly\n";

    //diagonarize S-matrix (aux.basis) and obtain inverse matrix etc.
    {
        DfInvMatrix dfinvMatrix(this->pPdfParam_);    //evaluate Sinv(Sab and Sgd)
        dfinvMatrix.DfInvMain();
    }

    {
        // density fitting
        DfDensityFittingX dfDensityFitting(this->pPdfParam_);
        dfDensityFitting.exec();
    }

    log << "Densityfitonly end\n";
}

