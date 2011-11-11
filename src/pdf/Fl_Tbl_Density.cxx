#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cassert>

#include "Fl_Tbl_Density.h"
#include "FileX.h"
#include "TlLogging.h"
#include "CnError.h"

int  Fl_Tbl_Density::flag1=0;        // flag for Constructor;
// first : flag1=0; after : flag1=1;
int  Fl_Tbl_Density::flag2=0;        // flag for Destructor;
// notGetmemory : flag2=0;
// getMemory    : flag2=1;
// int Fl_Tbl_Density::cGtoTotalNum=0; // s*1,p*3,d*5;

Fl_Tbl_Density::Fl_Tbl_Density()
{
    flag2=0;
    if (flag1 == 0) {
        this->prepare();

        if (FileX::isExist("fl_Table/DensityTable") == false) {
            this->makeTable();
        } else {
            this->setData();
        }

        flag2=1;
        flag1=1;
    } else {
        setData();
    }
}


Fl_Tbl_Density::~Fl_Tbl_Density()
{
}


int Fl_Tbl_Density::prepare()
{
    return 0;
}


int Fl_Tbl_Density::makeTable()
{
    TlLogging& log = TlLogging::getInstance();

    Fl_Geometry      FlGeom(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Density   FlGtoDen;
    TlPrdctbl        TPobj;    //これは、使うか？;

    const int AtomNum = FlGeom.getNumOfAtoms();
    const int AtomKindNum = FlGeom.getAtomKindNumber();

    static const char* Sorb[20] = {"s()"};
    static const char* Porb[20] = {"p(x)","p(y)","p(z)"};
    static const char* Dorb[20] = {"d(xy)","d(xz)","d(yz)","d(x2-y2)","d(3z2-r2)"};

    std::ofstream fo;
    std::string TblFile = "fl_Table/DensityTable";
    fo.open(TblFile.c_str(), std::ios::out | std::ios::trunc);
    if (!fo) {
        log.error("Cannot Open " + TblFile);
        CnErr.abort();
    }

    int basiscount = 0;
    for (int i = 0; i < AtomNum; ++i) {
        //int flag = 0;
        bool isFound = false;
        std::string Atm = FlGeom.getAtom(i);
        std::string Lb2 = FlGeom.getLabel(i);

        if (Atm == "X") {
            continue;
        }

        for (int j = 0; j < AtomKindNum; ++j) {
            int m = FlGtoDen.getStartposition(j);
            if ((Atm == FlGtoDen.getAtom(m)) &&
                    (Lb2 == FlGtoDen.getLabel(m))) {
                int start = m;
                int to = m + FlGtoDen.getTermnumber(j);
                for (int k = start; k < to; ++k) {
                    int cGtonum = k;

                    int iter = 0;
                    const char shell = FlGtoDen.getShell(k);
                    switch (shell) {
                    case 's':
                        iter = 1;
                        break;

                    case 'p':
                        iter = 3;
                        break;

                    case 'd':
                        iter = 5;
                        break;

                    default:
                        log.error("Error (making Table for Yahito)!!");
                        CnErr.abort();
                    }

                    int rept = 0;
                    for (int t = 0; t < iter; ++t) {
                        fo << std::setw(8) <<  basiscount
                        << std::setw(6) <<  cGtonum
                        << std::setw(6) <<     i
                        << std::setw(6) <<    Atm      <<"    " ;
                        ++basiscount;

                        switch (shell) {
                        case 's':
                            fo <<" s    "
                            << std::setw(10) << Sorb[rept]
                            <<"   ["<< Lb2<<"]"<< "\n";
                            break;

                        case 'p':
                            fo <<" p    "
                            << std::setw(10) << Porb[rept]
                            <<"   ["<<Lb2<<"]"<< "\n";
                            break;

                        case 'd':
                            fo <<" d    "
                            << std::setw(10) << Dorb[rept]
                            <<"   ["<<Lb2<<"]"<< "\n";
                            break;

                        default:
                            log.error("error(making Table for Yahito(aux))!!");
                            CnErr.abort();
                            break;
                        }
                        rept++;
                    }
                }

                isFound = true;
                break;
            }
        }

        if (isFound == false) {
            std::string msg = TlUtils::format("Fl_Tbl_Density::load(): not found basis set. %s", Atm.c_str());
            if (Lb2.empty() == false) {
                msg += TlUtils::format("@%s", Lb2.c_str());
            }
            CnErr.abort(msg);
        }
    }

    fo.close();
    setData();

    return 0;
}


int Fl_Tbl_Density::setData()
{
    TlLogging& log = TlLogging::getInstance();

    std::ifstream fi;

    std::string TblFile = "fl_Table/DensityTable";
    fi.open(TblFile.c_str(), std::ios::in);
    if (!fi) {
        log.error("Cannot open " + TblFile);
        CnErr.abort();
    }

    this->BsTbl.clear();
    while (fi) {
        int inttmp = 0;
        BasisTable basisTable;
        fi >> inttmp
        >> basisTable.cGtonum
        >> basisTable.InpAtomnumber
        >> basisTable.Atom
        >> basisTable.shell
        >> basisTable.dumy1
        >> basisTable.dumy2;
        this->BsTbl.push_back(basisTable);
    }

    fi.close();

    return 0;
}










