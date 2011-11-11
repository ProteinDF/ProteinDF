#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <iomanip>
#include <cassert>

#include "Fl_Tbl_Xcpot2.h"
#include "FileX.h"
#include "TlLogging.h"
#include "CnError.h"

int  Fl_Tbl_Xcpot2::flag1=0;         // flag for Constructor;
// first : flag1=0; after : flag1=1;
int  Fl_Tbl_Xcpot2::flag2=0;         // flag for Destructor;
// notGetmemory : flag2=0;
// getMemory    : flag2=1;
int Fl_Tbl_Xcpot2::cGtoTotalNum=0; // s*1,p*3,d*5;

//#######################################################################
Fl_Tbl_Xcpot2::Fl_Tbl_Xcpot2()
{
    TlLogging& log = TlLogging::getInstance();
//   std::cerr << "Fl_Tbl_Xcpot2::Fl_Tbl_Xcpot2() constructer" << std::endl;

    int memory;

    flag2=0;
    if (flag1==0) {
        prepare();
        memory = Estmemval();
        if (memory > MaxFlBasisTblmemory) {
            log.error("Cannot get memory !!(in Fl_Tbl_Xcpot2.cxx no constractor)");
            log.error("Used Memory is " + TlUtils::xtos(8 * memory) + " [byte] ");
            CnErr.abort();
        } else {
//      Log <<"Used Memory is " << 8 * memory <<" [byte] " <<"\n";
            getMemory();

            if (FileX::isExist("fl_Table/XcpotTable2") == false) {
                makeTable();
            } else {
                this->setData();
            }

            flag2=1;
        }
        flag1=1;
    } else {
        getMemory();
        setData();
    }
}

Fl_Tbl_Xcpot2::~Fl_Tbl_Xcpot2()
{
    int intI;
    for (intI=0; intI<cGtoTotalNum; intI++) {
        delete[] BsTbl[intI].Atom ;
        delete[] BsTbl[intI].shell;
        delete[] BsTbl[intI].dumy1;
        delete[] BsTbl[intI].dumy2;
    }

    delete[] BsTbl;
}

//#######################################################################
int Fl_Tbl_Xcpot2::prepare()
{
    TlLogging& log = TlLogging::getInstance();

    int i,j,k,m,start,to,flag;
    char *Atm = new char[MaxAtomName];
    char *Lb2 = new char[MaxLabelname];
    char shell;
    Fl_Geometry FlGeom(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Xcpot2      FlGtoXcpot2;

//   FlGeom.load();
//   FlGeom.open("fl_Geometry","read");
//   FlGeom.read();
// //FlGeom.show();
//   FlGeom.close();

///// modified by Fric(mizouchi) 2001/12/09 (s) /////
//    FlGtoXcpot2.open("fl_Gto_Xcpot2","read");
//   FlGtoXcpot2.open("fl_Gto_Xcpot","read");
///// modified by Fric(mizouchi) 2001/12/09 (e) /////

//   FlGtoXcpot2.read();
//FlGtoXcpot2.show();
//   FlGtoXcpot2.close();

    AtomNum     = FlGeom.getNumOfAtoms();
    AtomKindNum = FlGeom.getAtomKindNumber();

    cGtoTotalNum=0;

    for (i=0; i<AtomNum; i++) {
        flag=0;
        std::strcpy(Atm,FlGeom.getAtom(i).c_str());
        std::strcpy(Lb2,FlGeom.getLabel(i).c_str());
        for (j=0; j<AtomKindNum; j++) {
            m = FlGtoXcpot2.getStartposition(j);
            if ((strcmp(Atm,FlGtoXcpot2.getAtom(m).c_str())==0)  &&
                    (strcmp(Lb2,FlGtoXcpot2.getLabel(m).c_str())==0)) {
                start = m;
                to    = m + FlGtoXcpot2.getTermnumber(j);
                for (k=start; k<to; k++) {
                    shell   = FlGtoXcpot2.getShell(k);
                    int iter = 0;
                    switch (shell) {
                    case 's' :
                        iter=1;
                        break;
                    case 'p' :
                        iter=3;
                        break;
                    case 'd' :
                        iter=5;
                        break;
                    default  :
                        log.error("Error (making Table for Yahito(Xcot)))");
                        CnErr.abort();
                    }

                    cGtoTotalNum += iter;
                }
                flag=1;
            }
            if (flag==1) break;
        }
    }

    delete[] Atm;
    delete[] Lb2;

    return 0;
}
//#######################################################################
int Fl_Tbl_Xcpot2::Estmemval(void)
{
    int EstMemory;

    EstMemory =  sizeof(BasisTable*) * cGtoTotalNum
                 + sizeof(char*) * MaxLabelname *4     ;

    return EstMemory;
}
//#######################################################################
int  Fl_Tbl_Xcpot2::getMemory(void)
{
    int intI;
    //TblFile = new char[Maxfilename];
    //assert(TblFile != 0);

    BsTbl = new BasisTable[cGtoTotalNum];
    assert(BsTbl != 0);

    for (intI=0; intI<cGtoTotalNum; intI++) {
        BsTbl[intI].Atom  = new char[MaxAtomName];
        assert(BsTbl[intI].Atom != 0);
        BsTbl[intI].shell = new char[MaxLabelname];
        assert(BsTbl[intI].shell != 0);
        BsTbl[intI].dumy1 = new char[MaxLabelname];
        assert(BsTbl[intI].dumy1 != 0);
        BsTbl[intI].dumy2 = new char[MaxLabelname];
        assert(BsTbl[intI].dumy2 != 0);
    }

    return 0;
}
//#######################################################################
int Fl_Tbl_Xcpot2::makeTable(void)
{
    TlLogging& log = TlLogging::getInstance();

    Fl_Geometry      FlGeom(Fl_Geometry::getDefaultFileName());
    Fl_Gto_Xcpot2   FlGtoXcpot2;
    TlPrdctbl        TPobj;    //これは、使うか？;

    int i,j,k,m,t,rept,start,to,cGtonum,flag;
    int basiscount;
//  char Atm[MaxAtomName],Lb2[MaxLabelname];
    char *Atm = new char[MaxAtomName];
    char *Lb2 = new char[MaxLabelname];
//  char *Atm,*Lb2;
    char shell;

    static const char* Sorb[20]={"s()"};
    static const char* Porb[20]={"p(x)","p(y)","p(z)"};
    static const char* Dorb[20]={"d(xy)","d(xz)","d(yz)","d(x2-y2)","d(3z2-r2)"};

//   FlGeom.load();
//   FlGeom.open("fl_Geometry","read");
//   FlGeom.read();
// //FlGeom.show();
//   FlGeom.close();

///// modified by Fric(mizouchi) 2001/12/09 (s) /////
//    FlGtoXcpot2.open("fl_Gto_Xcpot2","read");
//   FlGtoXcpot2.open("fl_Gto_Xcpot","read");
///// modified by Fric(mizouchi) 2001/12/09 (e) /////

//   FlGtoXcpot2.read();
//FlGtoXcpot2.show();
//   FlGtoXcpot2.close();


    std::ofstream fo;
///// modified by Fric(mizouchi) 2001/12/09 (s) /////
//  TblFile = "fl_Table/Xcpot2Table";
    TblFile = "fl_Table/XcpotTable2";
///// modified by Fric(mizouchi) 2001/12/09 (e) /////

    fo.open(TblFile.c_str(), std::ios::out | std::ios::trunc);
    if (!fo) {
        log.error("Cannot Open " + TblFile);
        CnErr.abort();
    }

    basiscount = 0;
    for (i=0; i<AtomNum; i++) {
        flag=0;
        strcpy(Atm,FlGeom.getAtom(i).c_str());
        strcpy(Lb2,FlGeom.getLabel(i).c_str());
        for (j=0; j<AtomKindNum; j++) {
            m = FlGtoXcpot2.getStartposition(j);
            if ((strcmp(Atm,FlGtoXcpot2.getAtom(m).c_str())==0)  &&
                    (strcmp(Lb2,FlGtoXcpot2.getLabel(m).c_str())==0)) {
                start = m;
                to    = m + FlGtoXcpot2.getTermnumber(j);
                for (k=start; k<to; k++) {
                    cGtonum = k;
                    shell   = FlGtoXcpot2.getShell(k);
                    int iter = 0;
                    switch (shell) {
                    case 's' :
                        iter=1;
                        break;
                    case 'p' :
                        iter=3;
                        break;
                    case 'd' :
                        iter=5;
                        break;
                    default  :
                        log.error("Error (making Table for Yahito(Xpot2))!!");
                        CnErr.abort();
                        break;
                    }
                    rept=0;
                    for (t=0; t<iter; t++) {
                        fo << std::setw(8) <<  basiscount
                        << std::setw(6) <<  cGtonum
                        << std::setw(6) <<     i
                        << std::setw(6) <<    Atm      <<"    " ;
                        basiscount++;
                        switch (shell) {
                        case 's' :
                            fo <<" s    "
                            << std::setw(10) << Sorb[rept]
                            <<"   ["<< Lb2<<"]"<< "\n";
                            break;
                        case 'p' :
                            fo <<" p    "
                            << std::setw(10) << Porb[rept]
                            <<"   ["<<Lb2<<"]"<< "\n";
                            break;
                        case 'd' :
                            fo <<" d    "
                            << std::setw(10) << Dorb[rept]
                            <<"   ["<<Lb2<<"]"<< "\n";
                            break;
                        default  :
                            log.error("error(making Table for Yahiro(Xcpot2))");
                            CnErr.abort();
                            break;
                        }       //switch(shell);
                        rept++;
                    }         //for( t );
                }           //for( k );
                flag=1;
            }     //if;
            if (flag==1) break;
        }   //for( j );
    }    //for( i );

    fo.close();
    setData();

    delete[] Atm;
    delete[] Lb2;

    std::cerr << "Fl_Tbl_Xcpot2::makeTable() exit" << std::endl;
    return 0;
}
//#########################################################################
int  Fl_Tbl_Xcpot2::setData(void)
{
    TlLogging& log = TlLogging::getInstance();
    int intI;
    std::ifstream fi;
    int inttmp;

    TblFile = "fl_Table/XcpotTable2";

    fi.open(TblFile.c_str(), std::ios::in);
    if (!fi) {
        log.error("Cannot open " + TblFile);
        CnErr.abort();
    }

    for (intI=0; intI<cGtoTotalNum; intI++) {
        fi >>inttmp >> BsTbl[intI].cGtonum >> BsTbl[intI].InpAtomnumber
           >> BsTbl[intI].Atom       >> BsTbl[intI].shell
           >> BsTbl[intI].dumy1      >> BsTbl[intI].dumy2     ;
    }

    fi.close();

    return 0;

}

//#########################################################################
