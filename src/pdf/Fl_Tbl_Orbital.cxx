#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>

#include "Fl_Tbl_Orbital.h"
#include "FileX.h"
#include "TlUtils.h"
#include "CnError.h"

const std::string Fl_Tbl_Orbital::S_ORB_STR[] = {"s()"};
const std::string Fl_Tbl_Orbital::P_ORB_STR[] = {"p(x)", "p(y)", "p(z)"};
const std::string Fl_Tbl_Orbital::D_ORB_STR[] = {"d(xy)", "d(xz)", "d(yz)", "d(x2-y2)", "d(3z2-r2)"};

int Fl_Tbl_Orbital::flag1=0;        // flag for Constructor;
// first : flag1=0; after : flag1=1;
int Fl_Tbl_Orbital::cGtoTotalNum=0; // s*1,p*3,d*5;

Fl_Tbl_Orbital::Fl_Tbl_Orbital(const std::string& path) : m_sTableFilePath(path)
{
    if (flag1 == 0) {
        this->prepare();

        if (FileX::isExist(this->m_sTableFilePath) == false) {
            this->makeTable();
        } else {
            this->setData();
        }
    } else {
        this->setData();
    }
}

Fl_Tbl_Orbital::~Fl_Tbl_Orbital()
{
}

void Fl_Tbl_Orbital::prepare()
{
    Fl_Geometry FlGeom(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Orbital FlGtoOrb;

    AtomNum = FlGeom.getNumOfAtoms();
    AtomKindNum = FlGeom.getAtomKindNumber();

    cGtoTotalNum = 0;

    for (int i = 0; i < AtomNum; ++i) {
        std::string Atm = FlGeom.getAtom(i);
        std::string Lb2 = FlGeom.getLabel(i);

        for (int j = 0; j < AtomKindNum; ++j) {
            int m = FlGtoOrb.getStartposition(j);
            if (m < 0) {
                //std::cerr << "[TH] !(m >= 0) at Fl_Tbl_Orbital::prepare()." << std::endl;
                continue;
            }

            if ((Atm == FlGtoOrb.getAtom(m)) &&
                    (Lb2 == FlGtoOrb.getLabel(m))) {
                const int start = m;
                const int to = m + FlGtoOrb.getTermnumber(j);
                for (int k = start; k < to; ++k) {
                    const char shell = FlGtoOrb.getShell(k);
                    int iter = 0;
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
                        CnErr.abort(TlUtils::format("program error. %s: %d\n", __FILE__, __LINE__));
                    }

                    this->cGtoTotalNum += iter;
                }
                break;
            }
        }
    }
}

void Fl_Tbl_Orbital::makeTable()
{
    Fl_Geometry FlGeom(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Orbital FlGtoOrb;

    const int nAtomNum = FlGeom.getNumOfAtoms();
    const int nAtomKindNum = FlGeom.getAtomKindNumber();

    std::ofstream fo;
    fo.open(Fl_Tbl_Orbital::m_sTableFilePath.c_str(), std::ios::out | std::ios::trunc);
    if (!fo) {
        CnErr.abort(TlUtils::format("Cannot Open %s\n", Fl_Tbl_Orbital::m_sTableFilePath.c_str()));
    }

    int basiscount = 0;
    for (int i = 0; i < nAtomNum; ++i) {
        const std::string Atm = FlGeom.getAtom(i);
        const std::string Lb2 = FlGeom.getLabel(i);

        //std::cout << "Atm=\"" << Atm << "\", Lb2=\"" << Lb2 << "\"" << std::endl;

        for (int j = 0; j < nAtomKindNum; ++j) {
            int m = FlGtoOrb.getStartposition(j);
            if ((Atm == FlGtoOrb.getAtom(m)) &&
                    (Lb2 == FlGtoOrb.getLabel(m))) {
                const int start = m;
                const int to = m + FlGtoOrb.getTermnumber(j);
                for (int k = start; k < to; ++k) {
                    const int cGtonum = k;
                    const char shell = FlGtoOrb.getShell(k);
                    int iter = 0;
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
                        CnErr.abort(TlUtils::format("make table error. @ %s:%d\n",__FILE__, __LINE__));
                        break;
                    }

                    int rept = 0;
                    for (int t = 0; t < iter; ++t) {
                        fo << std::setw(8) <<  basiscount
                        << std::setw(6) <<  cGtonum
                        << std::setw(6) <<  i
                        << std::setw(6) <<  Atm
                        << "    ";
                        basiscount++;
                        switch (shell) {
                        case 's':
                            fo <<" s    "
                            << std::setw(10) << Fl_Tbl_Orbital::S_ORB_STR[rept]
                            <<"   ["<< Lb2 <<"]"<< "\n";
                            break;
                        case 'p':
                            fo <<" p    "
                            << std::setw(10) << Fl_Tbl_Orbital::P_ORB_STR[rept]
                            <<"   ["<< Lb2 <<"]"<< "\n";
                            break;
                        case 'd':
                            fo <<" d    "
                            << std::setw(10) << Fl_Tbl_Orbital::D_ORB_STR[rept]
                            <<"   ["<< Lb2 <<"]"<< "\n";
                            break;
                        default:
                            CnErr.abort(TlUtils::format("make table error. @ %s:%d\n",__FILE__, __LINE__));
                            break;
                        }
                        rept++;
                    }
                }

                break;
            }
        }
    }

    fo.close();

    //  this->setData();
}

void Fl_Tbl_Orbital::setData()
{
    std::ifstream ifs;
    ifs.open(this->m_sTableFilePath.c_str(), std::ios::in);
    if (!ifs) {
        CnErr.abort(TlUtils::format("Cannot open %s\n", Fl_Tbl_Orbital::m_sTableFilePath.c_str()));
    }

    int nBasisCount;
    //for (int intI = 0; intI < cGtoTotalNum; intI++){
    std::string sLine;
    while ((ifs) && (std::getline(ifs, sLine))) {
        std::istringstream iss(sLine);
        BasisTable basis;
        iss >> nBasisCount
        >> basis.cGtonum
        >> basis.InpAtomnumber
        >> basis.Atom
        >> basis.shell
        >> basis.dumy1
        >> basis.dumy2;

        this->m_BasisTable.push_back(basis);
    }

    this->cGtoTotalNum = this->m_BasisTable.size();

    ifs.close();
}







