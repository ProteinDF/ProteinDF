#ifndef TLESPPOP_H
#define TLESPPOP_H

#include <vector>
#include "Fl_Geometry.h"
#include "TlAtom.h"
#include "TlPosition.h"
#include "TlVector.h"

class TlEspPop {
 public:
  TlEspPop(const TlSerializeData& param);
  virtual ~TlEspPop();

  enum RESP_RESTRICTION { REST_NONE = 0, REST_QUADRIC, REST_HYPERBOLIC };

 public:
  void exec(std::string PMatrixFilePath = "");

  RESP_RESTRICTION getRespRestriction() const;
  void setRespRestriction(TlEspPop::RESP_RESTRICTION v);
  double getRestrictionParameterA() const;
  void setRestrictionParameterA(double a);
  double getRestrictionParameterB() const;
  void setRestrictionParameterB(double b);

  void verbose(bool yn);
  void saveMpacFilePath(const std::string& path);
  void saveDesignMatrixPath(const std::string& path);
  void savePredictedVectorPath(const std::string& path);
  void saveModelCoefVectorPath(const std::string& path);

 protected:
  std::vector<TlAtom> getRealAtoms();

  std::vector<TlPosition> getMerzKollmanGrids();
  std::vector<TlPosition> getMKGridsOnAtom(const TlPosition& center,
                                           const double radii);
  bool isInMolecule(const TlPosition& p, double coef);

  TlMatrix getInvDistanceMatrix();
  void makeDesignMatrix_MK(TlSymmetricMatrix* pDesignMat, TlVector* pPredicted);
  void makeDesignMatrix_quadric(const TlSymmetricMatrix& MK_designMat,
                                const TlVector& MK_predicted,
                                TlSymmetricMatrix* pDesignMat,
                                TlVector* pPredicted);
  void makeDesignMatrix_hyperbolic(const TlSymmetricMatrix& MK_designMat,
                                   const TlVector& MK_predicted,
                                   TlSymmetricMatrix* pDesignMat,
                                   TlVector* pPredicted);

  bool convCheck(const TlVector& modelCoef);
  void output(const TlVector& modelCoef);

 protected:
  static const double AU2ANG;
  static const double ANG2AU;

  TlSerializeData param_;
  Fl_Geometry flGeom_;

  std::vector<TlAtom> realAtoms_;
  std::vector<TlPosition> grids_;
  TlVector esp_;
  double sumOfCounterCharges_;
  double totalCharge_;

  // for RESP
  RESP_RESTRICTION resp_restriction_;
  double param_a_;
  double param_b_;
  TlVector expected_;

  // convcheck
  int itr_;
  int maxItr_;
  double maxErrorThreshold_;
  double rmsErrorThreshold_;
  TlVector prevModelCoef_;

  // debug
  bool verbose_;
  std::string saveMpacFilePath_;
  std::string saveDesignMatPath_;
  std::string savePredictedPath_;
  std::string saveModelCoefPath_;
};

#endif  // TLESPPOP_H
