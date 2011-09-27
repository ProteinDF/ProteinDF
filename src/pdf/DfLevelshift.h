#ifndef DFLEVELSHIFT_H
#define DFLEVELSHIFT_H

#include <string>
#include "DfObject.h"

/// Fock行列に対してレベルシフト法の処理を行うクラス
class DfLevelshift : public DfObject {
public:
    DfLevelshift(TlSerializeData* pPdfParam, int num_iter);
    virtual ~DfLevelshift();

public:
    void DfLshiftMain();

    /// フラグメントfragname のF 行列に対してレベルシフトの処置を施す
    void DfLshiftQclo(const std::string& fragname, int norbcut);

private:
    void main(const RUN_TYPE runType, int iteration, const std::string& fragname ="", bool bPdfQcloMode = false);
};

#endif // DFLEVELSHIFT_H
