#include <fstream>
#include <string>
#include <cmath>

#include "DfCqclomatrix.h"
#include "Fl_Tbl_AtomFragment.h"
#include "Fl_Tbl_Fragment.h"
#include "TlUtils.h"
#include "TlVector.h"
#include "TlMatrix.h"
#include "TlSymmetricMatrix.h"
#include "TlLogX.h"

DfCqclomatrix::DfCqclomatrix(TlSerializeData* pPdfParam) : DfObject(pPdfParam)
{
}


DfCqclomatrix::~DfCqclomatrix()
{
}


void DfCqclomatrix::main()
{
    //const TlSerializeData& pdfParam = *(this->pPdfParam_);
    TlLogX& Log = TlLogX::getInstance();

    //const int natom = atoi(this->m_flGbi["SCF"]["control-number-of-atoms"].c_str());
    const int natom = this->m_nNumOfAtoms;
    Fl_Fragment FlFrag;
    this->number_fragment = FlFrag.getNumOfFragments();

    int* atom_fragment        = new int[natom];
    int* basis_fragment       = new int[this->m_nNumOfAOs];
    //int* number_fragmentbasis = new int[number_fragment];
    std::vector<int> numOfFragmentBasis(this->m_nNumOfAOs, 0);
//   for (int i = 0; i < number_fragment; i++){
//     number_fragmentbasis[i] = 0;
//   }

    // read AtomFragmentTable
    {
        std::ifstream fi;
        int atom, fragment;

        fi.open("fl_Table/AtomFragmentTable");
        for (int i = 0; i < natom; i++) {
            fi >> atom >> fragment;
            atom_fragment[ atom ] = fragment;
        }
        fi.close();
    }

    // read OrbitalTable and make table which correlates basis and fragment
    {
        std::ifstream fi;
        int basis, atom, fragment, dummy;
        char linebuf[256];

        fi.open("fl_Table/OrbitalTable");
        for (int i = 0; i < m_nNumOfAOs; i++) {
            fi >> basis >> dummy >> atom;
            fragment = atom_fragment[ atom ];
            basis_fragment[ basis ] = fragment;
            numOfFragmentBasis[fragment]++;
            fi.getline(linebuf, 256);
        }
        fi.close();
    }

    // read S matrix
    TlSymmetricMatrix Spq;
    Spq.load("fl_Work/fl_Mtr_Spq.matrix");

    // fragment loop
    int count_basis = 0;
    double* eigval  = new double[m_nNumOfAOs];
    int*   eigfrag = new int[m_nNumOfAOs];

    for (int frag = 0; frag < number_fragment; frag++) {
        // If number of qclo in the fragment is not set in fl_Fragment,
        // the value is calculated and written to fl_Fragment.

        // create partial S matrix
        TlSymmetricMatrix fragmentS(numOfFragmentBasis[frag]);
        int counti = 0;
        for (int i = 0; i < m_nNumOfAOs; i++) {
            if (basis_fragment[i] == frag) {
                int countj = 0;
                for (int j = 0; j <= i; j++) {
                    if (basis_fragment[j] == frag) {
                        fragmentS(counti, countj) = Spq(i, j);
                        countj++;
                    }
                }
                counti++;
            }
        }
        Log << "partial S matrix " << frag << "\n";

        TlVector eigVal;
        TlMatrix eigVec;
        fragmentS.diagonal(&eigVal, &eigVec);

        Log << "eigenvalue of partial S matrix\n";
        eigVal.print(Log);

        for (int i = 0; i < numOfFragmentBasis[frag]; i++) {
            eigval[ count_basis+i ]  = eigVal[i];
            eigfrag[ count_basis+i ] = frag;
        }
        count_basis += numOfFragmentBasis[ frag ];
    }

    // read guess.lcao (temporary),
    // and write the number of independent basis
    {
        std::string type;
        switch (this->m_nMethodType) {
        case METHOD_RKS:
            type = "rks";        // RKS
            break;

        case METHOD_UKS:
            type = "uks-alpha";  // UKS alpha spin
            // type = "uks-beta";   // UKS beta spin
            break;

        case METHOD_ROKS:
            type = "roks";       // ROKS
            break;

        default:
            CnErr.abort();
            break;
        }

        std::ifstream fi;
        fi.open(std::string("./guess.lcao." + type).c_str(), std::ios::in);
        if (fi.rdstate()) {
            std::cerr << "cannot open file " << ("./guess.lcao." + type) << std::endl;
            abort();
        }

        std::string dummy_line;
        int row_dimension, col_dimension;

        fi >> dummy_line;
        fi >> row_dimension >> col_dimension;
        if (row_dimension != m_nNumOfAOs) {
            CnErr.abort("DfCqclomatrix", "", "main", "inputted guess lcao has illegal dimension");
        }
        this->m_nNumOfMOs = col_dimension;

        (*(this->pPdfParam_))["model"]["MOs"] = this->m_nNumOfMOs;
    }

    // count independant basis in each fragment
    // by the search of N smallest eigenvalues
    // ( N is the number of dependent basis )
    {
        for (int i = 0; i < this->m_nNumOfAOs - this->m_nNumOfMOs; i++) {
            for (int j = i+1; j < m_nNumOfAOs; j++) {
                if (eigval[j] < eigval[i]) {
                    double tmp1 = eigval[i];
                    eigval[i]   = eigval[j];
                    eigval[j]   = tmp1;
                    int tmp2   = eigfrag[i];
                    eigfrag[i]  = eigfrag[j];
                    eigfrag[j]  = tmp2;
                }
            }
            numOfFragmentBasis[ eigfrag[i] ]--;
        }
    }

    // put data and rewrite fl_Fragment
    {
        FlFrag.calcDefault(numOfFragmentBasis);
//     FlFrag.open( "fl_Input/fl_Fragment", "write" );
//     FlFrag.write();
//     FlFrag.close();
    }

    // write FragmentTable
    Fl_Tbl_Fragment  Tfrag;

    // write matrices
    switch (this->m_nMethodType) {
    case METHOD_RKS:
        this->main("rks");        // RKS
        break;

    case METHOD_UKS:
        this->main("uks-alpha");  // UKS alpha spin
        this->main("uks-beta");   // UKS beta spin
        break;

    case METHOD_ROKS:
        this->main("roks");         // ROKS
        break;

    default:
        CnErr.abort();
        break;
    }
    
    // release memory
    if (eigval)               delete [] eigval;
    if (eigfrag)              delete [] eigfrag;
    if (atom_fragment)        delete [] atom_fragment;
    if (basis_fragment)       delete [] basis_fragment;
    //if( number_fragmentbasis ) delete [] number_fragmentbasis;
}


void DfCqclomatrix::main(std::string type)
{
    TlLogX& Log = TlLogX::getInstance();

    TlMatrix guess_lcao;
    Fl_Tbl_Fragment  Tfrag;

    // read guess.lcao
    {
        std::ifstream fi;

        fi.open(std::string("./guess.lcao." + type).c_str(), std::ios::in);
        if (fi.rdstate()) {
            std::cerr << "cannot open file " << ("./guess.lcao." + type) << std::endl;
            abort();
        }

        std::string  dummy_line;
        int  row_dimension, col_dimension;

        fi >> dummy_line;
        fi >> row_dimension >> col_dimension;

        if (row_dimension != m_nNumOfAOs) {
            CnErr.abort("DfCqclomatrix", "", "main", "inputted guess lcao has illegal dimension");
        }

        this->m_nNumOfMOs = col_dimension;

        guess_lcao.resize(this->m_nNumOfAOs, this->m_nNumOfMOs);

        for (int i = 0; i < this->m_nNumOfAOs; i++) {
            for (int j = 0; j < this->m_nNumOfMOs; j++) {
                fi >> guess_lcao(i,j);
            }
        }

        Log << std::string("guess lcao " + type) << "\n";
        guess_lcao.print(Log);
    }

    // write Cqclo of each fragment
    for (int frag = 0; frag < number_fragment; frag++) {
        int  number_fragqclo = Tfrag.getNumberFragmentqclo(frag);
        TlMatrix Cqclo(m_nNumOfAOs, number_fragqclo);

        for (int fragqclo = 0; fragqclo < number_fragqclo; fragqclo++) {
            int qclo =0;
            if (type == "rks" || type == "roks") {
                qclo = Tfrag.getQclo(frag, fragqclo);
            } else if (type == "uks-alpha") {
                qclo = Tfrag.getQcloAlpha(frag, fragqclo);
            } else if (type == "uks-beta") {
                qclo = Tfrag.getQcloBeta(frag, fragqclo);
            }
            for (int basis = 0; basis < this->m_nNumOfAOs; basis++) {
                Cqclo(basis, fragqclo) = guess_lcao(basis, qclo);
            }
        }

        Cqclo.save("fl_Work/fl_Mtr_C.matrix.frag" + TlUtils::xtos(frag) + "." + type + "0");
    }

    // write Cprime.matrix.frag#.(type)0 for Level-Shift
    for (int frag = 0; frag < number_fragment; frag++) {
        int  number_fragqclo = Tfrag.getNumberFragmentqclo(frag);
        TlMatrix Cprime_frag(number_fragqclo, number_fragqclo);

        // Cprime_frag = Cqclo^-1 * Cqclo
        for (int i = 0; i < number_fragqclo; i++) {
            Cprime_frag(i, i) = 1.0;
        }

        std::string fname = "fl_Work/fl_Mtr_Cprime.matrix.frag" + TlUtils::xtos(frag) + "." + type + "0";
        Cprime_frag.save(fname);
    }
}

