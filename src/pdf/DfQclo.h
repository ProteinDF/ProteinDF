#ifndef DFQCLO_H
#define DFQCLO_H

#include <string>

#include "Fl_Fragment.h"
#include "DfObject.h"

/** QCLO based QCLO 法の実装クラス
 */
class DfQclo : public DfObject {
public:
    DfQclo(TlSerializeData* pPdfParam, int num_iter, bool bExecDiis);
    virtual ~DfQclo();

    void  DfQcloMain();

private:
    void combineCqclo(const std::string& runtype, int iteration);

protected:
    int number_fragment;      // Number of fragments
    bool m_bExecDiis;
    Fl_Fragment  FlFrag;
};

#endif // DFQCLO_H
