#ifndef DFCQCLOMATRIX_H
#define DFCQCLOMATRIX_H

#include "DfObject.h"

//  function : enaluate X-matrix by diagonalizeing orbital based S-matrix
//
//  input : S (hermite,positive define matrix)
//
//  output : Cqclo (Unitary ,column/row oriented,unpacked)

class DfCqclomatrix : public DfObject {
public:
    // flGbi は書き換えられる
    DfCqclomatrix(TlSerializeData* pPdfParam);
    ~DfCqclomatrix();

    void main();

private:
    void main(std::string type);

private:
    int number_fragment;           // Number of fragments
};

#endif // DFCQCLOMATRIX_H

