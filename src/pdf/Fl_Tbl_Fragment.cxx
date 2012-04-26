#include "Fl_Tbl_Fragment.h"
#include "TlParameter.h"
#include "FileX.h"
#include "TlLogging.h"

int  Fl_Tbl_Fragment :: flag = 0; // flag for Constructor
// first : flag = 0
// after : flag = 1


Fl_Tbl_Fragment::Fl_Tbl_Fragment(const Fl_Geometry& flGeom)
    : FlFrag(flGeom) {
//     FlFrag.open("fl_Input/fl_Fragment", "read");
//     FlFrag.read();

    number_qclo = 0;
    number_fragment = FlFrag.getNumOfFragments();
    number_fragmentqclo = new int[ number_fragment ];
    for (int i = 0; i < number_fragment; i++) {
        number_fragmentqclo[ i ] = FlFrag.getNumOfOrbitals(i);
        number_qclo += number_fragmentqclo[ i ];
    }

    FragTbl = new FragmentTable[ number_qclo ];
    for (int i = 0; i < number_qclo; i++) {
        FragTbl[ i ].qclo = i;
        FragTbl[ i ].fragment = 0;
        FragTbl[ i ].fragment_alpha = 0;
        FragTbl[ i ].fragment_beta = 0;
        FragTbl[ i ].fragmentqclo = i;
        FragTbl[ i ].fragmentqclo_alpha = i;
        FragTbl[ i ].fragmentqclo_beta = i;
    }

    if (flag == 0) {            // If the table was not written, write table.
        prepare();

        if (FileX::isExist("fl_Table/FragmentTable") == false) {
            this->makeTable();
        }

        flag = 1;               // The table was written.
    } else if (flag == 1) {      // If the table was already written.
        setData();
    }
//     FlFrag.close();
}


Fl_Tbl_Fragment::~Fl_Tbl_Fragment()
{
    if (FragTbl)             delete[] FragTbl;
    if (number_fragmentqclo) delete[] number_fragmentqclo;
}


void Fl_Tbl_Fragment::prepare()
{
    //Fl_Globalinput  GbiScf (">>>>SCF" );
    TlParameter flGbi;
    flGbi.load("fl_Input/fl_Globalinput");

    std::string method = flGbi["SCF"]["method"];

    int   qcloindex  = 0;
    int*  number_occ = new int[ number_fragment ];

    if (method == "nsp" || method == "roks") {
        for (int i = 0; i < number_fragment; i++) {
            number_occ[ i ] = FlFrag.getNumOfOccupiedOrbitals(i);
        }
        for (int i = 0; i < number_fragment; i++) {
            for (int j = 0; j < number_occ[ i ]; j++) {
                FragTbl[ qcloindex ].qclo = qcloindex;
                FragTbl[ qcloindex ].fragment = i;
                FragTbl[ qcloindex ].fragment_alpha = i;
                FragTbl[ qcloindex ].fragment_beta = i;
                FragTbl[ qcloindex ].fragmentqclo = j;
                FragTbl[ qcloindex ].fragmentqclo_alpha = j;
                FragTbl[ qcloindex ].fragmentqclo_beta = j;
                qcloindex++;
            }
        }
        for (int i = 0; i < number_fragment; i++) {
            for (int j = number_occ[ i ]; j < number_fragmentqclo[ i ]; j++) {
                FragTbl[ qcloindex ].qclo = qcloindex;
                FragTbl[ qcloindex ].fragment = i;
                FragTbl[ qcloindex ].fragment_alpha = i;
                FragTbl[ qcloindex ].fragment_beta = i;
                FragTbl[ qcloindex ].fragmentqclo = j;
                FragTbl[ qcloindex ].fragmentqclo_alpha = j;
                FragTbl[ qcloindex ].fragmentqclo_beta = j;
                qcloindex++;
            }
        }
    }

    else if (method == "sp") {
        for (int i = 0; i < number_fragment; i++) {
            number_occ[ i ] = FlFrag.getNumOfOccupiedAlphaOrbitals(i);
        }
        for (int i = 0; i < number_fragment; i++) {
            for (int j = 0; j < number_occ[ i ]; j++) {
                FragTbl[ qcloindex ].qclo = qcloindex;
                FragTbl[ qcloindex ].fragment_alpha = i;
                FragTbl[ qcloindex ].fragmentqclo_alpha = j;
                qcloindex++;
            }
        }
        for (int i = 0; i < number_fragment; i++) {
            for (int j = number_occ[ i ]; j < number_fragmentqclo[ i ]; j++) {
                FragTbl[ qcloindex ].qclo = qcloindex;
                FragTbl[ qcloindex ].fragment_alpha = i;
                FragTbl[ qcloindex ].fragmentqclo_alpha = j;
                qcloindex++;
            }
        }

        qcloindex = 0;
        for (int i = 0; i < number_fragment; i++) {
            number_occ[ i ] = FlFrag.getNumOfOccupiedBetaOrbitals(i);
        }
        for (int i = 0; i < number_fragment; i++) {
            for (int j = 0; j < number_occ[ i ]; j++) {
                FragTbl[ qcloindex ].fragment_beta = i;
                FragTbl[ qcloindex ].fragmentqclo_beta = j;
                qcloindex++;
            }
        }
        for (int i = 0; i < number_fragment; i++) {
            for (int j = number_occ[ i ]; j < number_fragmentqclo[ i ]; j++) {
                FragTbl[ qcloindex ].fragment_beta = i;
                FragTbl[ qcloindex ].fragmentqclo_beta = j;
                qcloindex++;
            }
        }
    }
    if (number_occ) delete[] number_occ;
}


void Fl_Tbl_Fragment::makeTable()
{
    TlLogging& log = TlLogging::getInstance();

    std::ofstream  fo;
    TblFile = "fl_Table/FragmentTable";
    fo.open(TblFile.c_str(), std::ios::out | std::ios::trunc);
    if (!fo) {
        log.error("Cannot open " + TblFile);
        CnErr.abort();
    }
    for (int i = 0; i < number_qclo; i++) {
        fo << std::setw(8) << FragTbl[ i ].qclo
        << std::setw(6) << FragTbl[ i ].fragment_alpha
        << std::setw(6) << FragTbl[ i ].fragmentqclo_alpha
        << std::setw(6) << FragTbl[ i ].fragment_beta
        << std::setw(6) << FragTbl[ i ].fragmentqclo_beta
        << std::endl;
    }
    fo.close();
}


void Fl_Tbl_Fragment::setData()
{
    TlLogging& log = TlLogging::getInstance();

    std::ifstream  fi;
    TblFile = "fl_Table/FragmentTable";
    fi.open(TblFile.c_str(), std::ios::in);
    if (!fi) {
        log.error("Cannot open " + TblFile);
        CnErr.abort("Fl_Tbl_Fragment", "", "setData", "Cannot open FragmentTable") ;
    }
    for (int i = 0; i < number_fragment; i++) {
        number_fragmentqclo[ i ] = 0;
    }
    for (int i = 0; i < number_qclo; i++) {
        fi >> FragTbl[ i ].qclo
        >> FragTbl[ i ].fragment_alpha
        >> FragTbl[ i ].fragmentqclo_alpha
        >> FragTbl[ i ].fragment_beta
        >> FragTbl[ i ].fragmentqclo_beta;
        FragTbl[ i ].fragment = FragTbl[ i ].fragment_alpha;
        FragTbl[ i ].fragmentqclo = FragTbl[ i ].fragmentqclo_alpha;
        number_fragmentqclo[ FragTbl[ i ].fragment ]++;
    }
    fi.close();
}


int  Fl_Tbl_Fragment :: getQclo(int fragindex, int fragqcloindex)
{
    for (int i = 0; i < number_qclo; i++) {
        if (FragTbl[ i ].fragment     == fragindex &&
                FragTbl[ i ].fragmentqclo == fragqcloindex)
            return i;
    }
    CnErr.abort("Fl_Tbl_Fragment", "", "getQclo" , "Cannot find such mo");
    return -1;
}


int  Fl_Tbl_Fragment :: getQcloAlpha(int fragindex, int fragqcloindex)
{
    for (int i = 0; i < number_qclo; i++) {
        if (FragTbl[ i ].fragment_alpha     == fragindex &&
                FragTbl[ i ].fragmentqclo_alpha == fragqcloindex)
            return i;
    }
    CnErr.abort("Fl_Tbl_Fragment", "", "getQcloAlpha" , "Cannot find such mo");
    return -1;
}


int  Fl_Tbl_Fragment :: getQcloBeta(int fragindex, int fragqcloindex)
{
    for (int i = 0; i < number_qclo; i++) {
        if (FragTbl[ i ].fragment_beta     == fragindex &&
                FragTbl[ i ].fragmentqclo_beta == fragqcloindex)
            return i;
    }
    CnErr.abort("Fl_Tbl_Fragment", "", "getQcloBeta" , "Cannot find such mo");
    return -1;
}


int  Fl_Tbl_Fragment :: getFragment(int qcloindex)
{
    return FragTbl[ qcloindex ].fragment;
}


int  Fl_Tbl_Fragment :: getFragmentAlpha(int qcloindex)
{
    return FragTbl[ qcloindex ].fragment_alpha;
}


int  Fl_Tbl_Fragment :: getFragmentBeta(int qcloindex)
{
    return FragTbl[ qcloindex ].fragment_beta;
}


int  Fl_Tbl_Fragment :: getFragmentqclo(int qcloindex)
{
    return FragTbl[ qcloindex ].fragmentqclo;
}


int  Fl_Tbl_Fragment :: getFragmentqcloAlpha(int qcloindex)
{
    return FragTbl[ qcloindex ].fragmentqclo_alpha;
}


int  Fl_Tbl_Fragment :: getFragmentqcloBeta(int qcloindex)
{
    return FragTbl[ qcloindex ].fragmentqclo_beta;
}


int  Fl_Tbl_Fragment :: getNumberFragmentqclo(int fragindex)
{
    return number_fragmentqclo[ fragindex ];
}


/* Fl_Tbl_Fragment.cxx */
/* end of file */
