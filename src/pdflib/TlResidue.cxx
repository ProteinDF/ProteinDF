#include <cassert>
#include "TlResidue.h"
#include "TlUtils.h"

const char* TlResidue::residueTable[] = {
    "ALA",   "ALAC-", "ALAN+", "ARG", "ASNN+",
    "ASN",   "ASP"  , "CYS",   "GLN", "GLU",
    "GLY",   "GLYC-", "HID",   "HIP", "HIS",
    "LEU",   "LEUN+", "LEUC-", "LYS", "MET",
    "PHE",   "PHEC-", "PRO",   "SER", "SERC-",
    "THR",   "THRN+", "THRC-", "TRP", "TYR",
    "TYRC-", "TYRN+", "VALC-", "VAL", "HEM"
};

int TlResidue::number()
{
    return (sizeof(TlResidue::residueTable) / sizeof(char*));
}

// [TH] ���ͤ��������������ϣ����ѹ����������褤��
std::string TlResidue::label(int p)
{
    assert((1 <= p)  && p <= this->number());

//   if( p>number() ){
//     CnErr.abort( "TlResidue", "", "label", TlUtils::format("illegal order of amino acid information table %d", p));
//   }
    //    printf("%s\n",residue[p-1]);
    return std::string(TlResidue::residueTable[p-1]);
}

TlResidue::TlResidue()
{
}

bool TlResidue::isResidue(const std::string& str)
{
    bool bAnswer = false;

    int nNumOfResidues = sizeof(TlResidue::residueTable) / sizeof(char*);

    for (int i = 0; i < nNumOfResidues; ++i) {
        if (str == std::string(TlResidue::residueTable[i])) {
            bAnswer = true;
            break;
        }
    }

    return bAnswer;
}
