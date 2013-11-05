#include "DfFunctional.h"

////////////////////////////////////////////////////////////////////////
// DfFunctional

DfFunctional_LDA::DfFunctional_LDA()
{
}

DfFunctional_LDA::~DfFunctional_LDA()
{
}

////////////////////////////////////////////////////////////////////////
// DfFunctional_GGA

DfFunctional_GGA::DfFunctional_GGA()
{
    this->numOfFunctionalTerms_ = 1;
    this->numOfDerivativeFunctionalTerms_ = 1;
}

DfFunctional_GGA::~DfFunctional_GGA()
{
}

// for grid-free
FunctionalSets DfFunctional_GGA::getFunctional_GF(const TlVector& rhoAs, 
                                                  const TlVector& rhoBs,
                                                  const TlVector& xAs,
                                                  const TlVector& xBs)
{
    const index_type dim = rhoAs.getSize();
    assert(dim == rhoBs.getSize());
    assert(dim == xAs.getSize());
    assert(dim == xBs.getSize());

    const int numOfTerms = this->getNumOfFunctionalTerms();
    FunctionalSets answer(numOfTerms, dim);
    for (index_type i = 0; i < dim; ++i) {
        const TlMatrix fs = this->getFunctionalCore(rhoAs[i], rhoBs[i],
                                                    xAs[i], xBs[i]);
        
        for (int term = 0; term < numOfTerms; ++term) {
            answer.FA_termR.set(term, i, fs.get(FA_R, term));
            answer.FA_termX.set(term, i, fs.get(FA_X, term));
            answer.FB_termR.set(term, i, fs.get(FB_R, term));
            answer.FB_termX.set(term, i, fs.get(FB_X, term));
        }
    }

    return answer;
}

DerivativeFunctionalSets DfFunctional_GGA::getDerivativeFunctional_GF(const TlVector& rhoAs, 
                                                                      const TlVector& rhoBs,
                                                                      const TlVector& xAs,
                                                                      const TlVector& xBs)
{
    const index_type dim = rhoAs.getSize();
    assert(dim == rhoBs.getSize());
    assert(dim == xAs.getSize());
    assert(dim == xBs.getSize());

    const int numOfTerms = this->getNumOfDerivativeFunctionalTerms();
    DerivativeFunctionalSets answer(numOfTerms, dim);
    for (index_type i = 0; i < dim; ++i) {
        assert(std::fabs(rhoAs[i] - rhoBs[i]) < 1.0E-5);
        assert(std::fabs(xAs[i] - xBs[i]) < 1.0E-5);
        const TlMatrix dfs = this->getDerivativeFunctionalCore(rhoAs[i], rhoBs[i],
                                                               xAs[i], xBs[i]);
        
        for (int term = 0; term < numOfTerms; ++term) {
            answer.rFrRhoA_R.set(term, i, dfs.get(RA_R, term));
            answer.rFrRhoA_X.set(term, i, dfs.get(RA_X, term));
            answer.rFrRhoB_R.set(term, i, dfs.get(RB_R, term));
            answer.rFrRhoB_X.set(term, i, dfs.get(RB_X, term));
            answer.rFrGAA_R.set(term, i, dfs.get(GAA_R, term));
            answer.rFrGAA_X.set(term, i, dfs.get(GAA_X, term));
            answer.rFrGAB_R.set(term, i, dfs.get(GAB_R, term));
            answer.rFrGAB_X.set(term, i, dfs.get(GAB_X, term));
            answer.rFrGBB_R.set(term, i, dfs.get(GBB_R, term));
            answer.rFrGBB_X.set(term, i, dfs.get(GBB_X, term));
        }
    }

    return answer;
}

