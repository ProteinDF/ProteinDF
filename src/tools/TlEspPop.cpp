#include "TlEspPop.h"

#include "DfObject.h"
#include "Fl_Geometry.h"
#include "TlEspField.h"
#include "TlLebedevGrid.h"
#include "TlMath.h"
#include "TlMsgPack.h"
#include "TlPrdctbl.h"
#include "TlSymmetricMatrix.h"
#include "TlUtils.h"
#include "TlVector.h"

const double TlEspPop::AU2ANG = 0.5291772108;
const double TlEspPop::ANG2AU = 1.889762;

TlEspPop::TlEspPop(const TlSerializeData& param)
    : param_(param),
      flGeom_(param["coordinates"]),
      resp_restriction_(TlEspPop::REST_NONE),
      param_a_(0.0005),
      param_b_(0.1),
      itr_(1),
      maxItr_(100),
      maxErrorThreshold_(1.0E-3),
      rmsErrorThreshold_(1.0E-2),
      verbose_(false) {
  // initialize --------------------------------------------------------------
  // set real atoms
  this->realAtoms_ = this->getRealAtoms();

  this->expected_ = TlVector(this->realAtoms_.size());
}

TlEspPop::~TlEspPop() {}

TlEspPop::RESP_RESTRICTION TlEspPop::getRespRestriction() const {
  return this->resp_restriction_;
}

void TlEspPop::setRespRestriction(TlEspPop::RESP_RESTRICTION v) {
  this->resp_restriction_ = v;
}

double TlEspPop::getRestrictionParameterA() const { return this->param_a_; }

void TlEspPop::setRestrictionParameterA(const double a) { this->param_a_ = a; }

double TlEspPop::getRestrictionParameterB() const { return this->param_b_; }

void TlEspPop::setRestrictionParameterB(const double b) { this->param_b_ = b; }

void TlEspPop::verbose(bool yn) { this->verbose_ = yn; }

void TlEspPop::saveMpacFilePath(const std::string& path) {
  this->saveMpacFilePath_ = path;
}

void TlEspPop::saveDesignMatrixPath(const std::string& path) {
  this->saveDesignMatPath_ = path;
}

void TlEspPop::savePredictedVectorPath(const std::string& path) {
  this->savePredictedPath_ = path;
}

void TlEspPop::saveModelCoefVectorPath(const std::string& path) {
  this->saveModelCoefPath_ = path;
}

void TlEspPop::exec(std::string PMatrixFilePath) {
  // generate grids
  this->grids_ = this->getMerzKollmanGrids();
  const std::size_t numOfGrids = this->grids_.size();
  if (this->verbose_) {
    std::cerr << TlUtils::format("# grids: %ld", numOfGrids) << std::endl;
  }

  // 密度行列の読み込み
  TlSymmetricMatrix P;
  if (PMatrixFilePath.empty()) {
    const int iteration = this->param_["num_of_iterations"].getInt();
    DfObject::RUN_TYPE runType = DfObject::RUN_RKS;
    DfObject dfObject(&(this->param_));
    PMatrixFilePath = dfObject.getPpqMatrixPath(runType, iteration);
  }
  if (this->verbose_) {
    std::cerr << "loading: " << PMatrixFilePath << std::endl;
  }
  P.load(PMatrixFilePath);

  // ESP計算
  TlEspField espFld(this->param_);
  this->esp_ = espFld.makeEspFld(P, this->grids_);

  // save data by msgpack
  if (this->saveMpacFilePath_.empty() != true) {
    TlSerializeData output;
    output["version"] = "2014.0";
    output["num_of_grids"] = numOfGrids;
    for (std::size_t gridIndex = 0; gridIndex < numOfGrids; ++gridIndex) {
      TlSerializeData pos;
      pos.pushBack(this->grids_[gridIndex].x() * TlEspPop::AU2ANG);
      pos.pushBack(this->grids_[gridIndex].y() * TlEspPop::AU2ANG);
      pos.pushBack(this->grids_[gridIndex].z() * TlEspPop::AU2ANG);

      output["grids"].pushBack(pos);
    }
    output["grid_unit"] = "angstrom";

    for (std::size_t i = 0; i < numOfGrids; ++i) {
      output["ESP"].setAt(i, this->esp_[i]);
    }

    TlMsgPack mpac(output);
    if (this->saveMpacFilePath_.empty() != true) {
      if (this->verbose_) {
        std::cerr << "output mpac file: " << this->saveMpacFilePath_
                  << std::endl;
      }
      mpac.save(this->saveMpacFilePath_);
    }
  }

  // solve MK
  TlSymmetricMatrix designMat;
  TlVector predicted;
  TlVector modelCoef;
  // setup MK data
  this->makeDesignMatrix_MK(&designMat, &predicted);

  if (this->resp_restriction_ == REST_NONE) {
    // solve
    TlSymmetricMatrix invDesignMat = designMat;
    invDesignMat.inverse();

    modelCoef = invDesignMat * predicted;
  } else {
    TlSymmetricMatrix MK_designMat = designMat;
    TlVector MK_predicted = predicted;

    while (this->itr_ < this->maxItr_) {
      switch (this->resp_restriction_) {
        case REST_QUADRIC:
          this->makeDesignMatrix_quadric(MK_designMat, MK_predicted, &designMat,
                                         &predicted);
          break;

        case REST_HYPERBOLIC:
          this->makeDesignMatrix_hyperbolic(MK_designMat, MK_predicted,
                                            &designMat, &predicted);
          break;

        default:
          break;
      }

      // solve
      TlSymmetricMatrix invDesignMat = designMat;
      invDesignMat.inverse();

      modelCoef = invDesignMat * predicted;
      this->expected_ = modelCoef;

      if (this->convCheck(modelCoef) == true) {
        break;
      }

      ++(this->itr_);
    }
  }

  this->output(modelCoef);
  if (this->saveDesignMatPath_.empty() != true) {
    if (this->verbose_) {
      std::cerr << "design mat path: " << this->saveDesignMatPath_ << std::endl;
    }
    designMat.save(this->saveDesignMatPath_);
  }
  if (this->savePredictedPath_.empty() != true) {
    if (this->verbose_) {
      std::cerr << "predicetd vtr path: " << this->savePredictedPath_
                << std::endl;
    }
    predicted.save(this->savePredictedPath_);
  }
  if (this->saveModelCoefPath_.empty() != true) {
    if (this->verbose_) {
      std::cerr << "model coef path: " << this->saveModelCoefPath_ << std::endl;
    }
    modelCoef.save(this->saveModelCoefPath_);
  }
}

std::vector<TlAtom> TlEspPop::getRealAtoms() {
  // const Fl_Geometry flGeom(this->param_["coordinates"]); // 単位はa.u.
  const int numOfAllAtoms = this->flGeom_.getNumOfAtoms();
  const int numOfRealAtoms = numOfAllAtoms - this->flGeom_.getNumOfDummyAtoms();

  this->sumOfCounterCharges_ = 0.0;
  std::vector<TlAtom> realAtoms(numOfRealAtoms);
  {
    int realAtomIndex = 0;
    for (int i = 0; i < numOfAllAtoms; ++i) {
      const TlAtom atom = this->flGeom_.getAtom(i);
      if (atom.getSymbol() != "X") {
        realAtoms[realAtomIndex] = atom;
        ++realAtomIndex;
      } else {
        double charge = atom.getCharge();
        this->sumOfCounterCharges_ += charge;
      }
    }
    assert(realAtomIndex == numOfRealAtoms);
  }

  return realAtoms;
}

// radii: atomic unit
std::vector<TlPosition> TlEspPop::getMKGridsOnAtom(const TlPosition& center,
                                                   const double radii) {
  TlLebedevGrid lebGrid;
  const std::vector<int> gridList = lebGrid.getSupportedGridNumber();

  std::vector<TlPosition> grids;
  const double r = radii * AU2ANG;  // to angstroam unit

  const double area = 4.0 * TlMath::PI() * r * r;
  std::vector<int>::const_iterator it = std::upper_bound(
      gridList.begin(), gridList.end(), static_cast<int>(area));
  const int numOfGrids = *it;
  if (this->verbose_) {
    std::cerr << TlUtils::format("area=%8.3f ANG^2 grids=%d", area, numOfGrids)
              << std::endl;
  }

  std::vector<TlPosition> layerGrids;
  std::vector<double> layerWeights;
  lebGrid.getGrids(numOfGrids, &layerGrids, &layerWeights);
  for (int grid = 0; grid < numOfGrids; ++grid) {
    layerGrids[grid] *= radii;
    layerGrids[grid].shiftBy(center);
  }

  grids.insert(grids.end(), layerGrids.begin(), layerGrids.end());

  return grids;
}

bool TlEspPop::isInMolecule(const TlPosition& p, double coef) {
  bool answer = false;

  const int numOfAtoms = this->flGeom_.getNumOfAtoms();
  for (int i = 0; i < numOfAtoms; ++i) {
    const TlAtom& atom = this->flGeom_.getAtom(i);
    const std::string symbol = atom.getSymbol();
    if (symbol != "X") {
      const double distance = p.distanceFrom(atom.getPosition());
      const double vdwr =
          TlPrdctbl::getVdwRadii(TlPrdctbl::getAtomicNumber(symbol)) * ANG2AU;
      if (distance < vdwr * coef) {
        answer = true;
        break;
      }
    }
  }

  return answer;
}

std::vector<TlPosition> TlEspPop::getMerzKollmanGrids() {
  static const double layers[] = {1.4, 1.6, 1.8, 2.0};
  static const int numOfLayers = sizeof(layers) / sizeof(layers[0]);

  // Fl_Geometry flGeom(this->param_["coordinates"]); // 単位はa.u.

  std::vector<TlPosition> allGrids;
  std::size_t numOfGrids = 0;
  std::size_t screened = 0;
  const int numOfAtoms = this->flGeom_.getNumOfAtoms();
  for (int i = 0; i < numOfAtoms; ++i) {
    const TlAtom atom = this->flGeom_.getAtom(i);
    const std::string symbol = atom.getSymbol();
    if (symbol != "X") {
      const double vdwr =
          TlPrdctbl::getVdwRadii(TlPrdctbl::getAtomicNumber(symbol)) * ANG2AU;

      for (int layer_index = 0; layer_index < numOfLayers; ++layer_index) {
        const double coef = layers[layer_index];

        // define the number of grids
        std::vector<TlPosition> grids =
            this->getMKGridsOnAtom(atom.getPosition(), coef * vdwr);

        // check in molecule
        std::vector<TlPosition>::iterator itEnd = grids.end();
        for (std::vector<TlPosition>::iterator it = grids.begin(); it != itEnd;
             ++it) {
          if (this->isInMolecule(*it, coef) != true) {
            allGrids.push_back(*it);
          } else {
            ++screened;
          }
          ++numOfGrids;
        }
      }
    }
  }

  if (this->verbose_) {
    std::cerr << TlUtils::format("screened %ld/%ld", screened, numOfGrids)
              << std::endl;
  }
  return allGrids;
}

TlMatrix TlEspPop::getInvDistanceMatrix() {
  // make 1/r distance table
  const int numOfRealAtoms = this->realAtoms_.size();
  const int numOfGrids = this->grids_.size();
  TlMatrix d(numOfRealAtoms, numOfGrids);
  if (this->verbose_) {
    std::cerr << TlUtils::format("# of atoms: %d", numOfRealAtoms) << std::endl;
    std::cerr << TlUtils::format("# of grids: %d", numOfGrids) << std::endl;
  }

  for (int i = 0; i < numOfRealAtoms; ++i) {
    const TlPosition& posA = this->realAtoms_[i].getPosition();
    // std::cout << TlUtils::format("MK > % f, %f , %f", posA.x(), posA.y(),
    // posA.z()) << std::endl;
    for (int j = 0; j < numOfGrids; ++j) {
      const double r_ai = posA.distanceFrom(this->grids_[j]);

      d.set(i, j, 1.0 / r_ai);
    }
  }
  // d.save("d.mat");

  return d;
}

void TlEspPop::makeDesignMatrix_MK(TlSymmetricMatrix* pDesignMat,
                                   TlVector* pPredicted) {
  assert(pDesignMat != NULL);
  assert(pPredicted != NULL);
  assert(this->esp_.getSize() == static_cast<int>(this->grids_.size()));

  const TlMatrix d = this->getInvDistanceMatrix();

  const int numOfRealAtoms = this->realAtoms_.size();
  pDesignMat->resize(numOfRealAtoms + 1);
  pPredicted->resize(numOfRealAtoms + 1);
  for (int a = 0; a < numOfRealAtoms; ++a) {
    const TlVector r_a = d.getRowVector(a);
    assert(r_a.getSize() == static_cast<int>(this->grids_.size()));

    // a == b
    {
      TlVector r_a2 = r_a;
      r_a2.dot(r_a);
      pDesignMat->set(a, a, r_a2.sum());
    }

    // a != b
    for (int b = 0; b < a; ++b) {
      TlVector r_b = d.getRowVector(b);
      r_b.dot(r_a);
      pDesignMat->set(a, b, r_b.sum());
    }

    // for Lagurange
    pDesignMat->set(a, numOfRealAtoms, 1.0);
    pDesignMat->set(numOfRealAtoms, a, 1.0);

    //
    {
      TlVector r = r_a;
      (*pPredicted)[a] = r.dot(this->esp_).sum();
    }
  }

  if (this->verbose_) {
    std::cerr << TlUtils::format(
                     "total charge: % 8.2f = % 8.2f - % 8.2f",
                     this->totalCharge_ - this->sumOfCounterCharges_,
                     this->totalCharge_, this->sumOfCounterCharges_)
              << std::endl;
  }
  pPredicted->set(numOfRealAtoms,
                  this->totalCharge_ - this->sumOfCounterCharges_);
}

void TlEspPop::makeDesignMatrix_quadric(const TlSymmetricMatrix& MK_designMat,
                                        const TlVector& MK_predicted,
                                        TlSymmetricMatrix* pDesignMat,
                                        TlVector* pPredicted) {
  *pDesignMat = MK_designMat;
  *pPredicted = MK_predicted;

  const double param_a = this->param_a_;
  const TlVector& q = this->expected_;
  const double target = 0.0;  // TODO

  const int numOfRealAtoms = this->realAtoms_.size();
  for (int a = 0; a < numOfRealAtoms; ++a) {
    const double rest = -2.0 * param_a * (target - q[a]);

    pDesignMat->add(a, a, rest);
    (*pPredicted)[a] += target * rest;
  }
}

void TlEspPop::makeDesignMatrix_hyperbolic(
    const TlSymmetricMatrix& MK_designMat, const TlVector& MK_predicted,
    TlSymmetricMatrix* pDesignMat, TlVector* pPredicted) {
  *pDesignMat = MK_designMat;
  *pPredicted = MK_predicted;

  const double param_a = this->param_a_;
  const double param_b2 = this->param_b_ * this->param_b_;
  const TlVector& q = this->expected_;
  const double target = 0.0;  // TODO

  const int numOfRealAtoms = this->realAtoms_.size();
  for (int a = 0; a < numOfRealAtoms; ++a) {
    const double ya = q[a];
    const double rest = param_a * ya * (1.0 / std::sqrt(ya * ya + param_b2));

    pDesignMat->add(a, a, rest);
    (*pPredicted)[a] += target * rest;
  }
  pPredicted->set(numOfRealAtoms,
                  this->totalCharge_ - this->sumOfCounterCharges_);
}

bool TlEspPop::convCheck(const TlVector& modelCoef) {
  bool judge = false;
  if (this->itr_ > 1) {
    TlVector diff = modelCoef - this->prevModelCoef_;

    const double maxErr = diff.getMaxAbsoluteElement();
    const int numOfSize = modelCoef.getSize();
    const double rmsErr = diff.dot(diff).sum() / double(numOfSize);

    std::cerr << TlUtils::format("#%3d: MAX: % e RMS: % e", this->itr_, maxErr,
                                 rmsErr)
              << std::endl;
    if ((maxErr < this->maxErrorThreshold_) &&
        (rmsErr < this->rmsErrorThreshold_)) {
      judge = true;
    }
  }

  if (judge == false) {
    this->prevModelCoef_ = modelCoef;
  }

  return judge;
}

void TlEspPop::output(const TlVector& modelCoef) {
  const int numOfRealAtoms = this->realAtoms_.size();
  assert(modelCoef.getSize() == numOfRealAtoms + 1);

  switch (this->resp_restriction_) {
    case REST_QUADRIC:
      std::cout << TlUtils::format("RESP charge (rest: quadric; a=% f)",
                                   this->param_a_)
                << std::endl;
      break;

    case REST_HYPERBOLIC:
      std::cout << TlUtils::format(
                       "RESP charge (rest: hyperbolic; a=% f, b=% f)",
                       this->param_a_, this->param_b_)
                << std::endl;
      break;

    default:
      std::cout << ">>> ESP (MK) charge" << std::endl;
      break;
  }
  double totalCharge = 0.0;
  for (int i = 0; i < numOfRealAtoms; ++i) {
    const double charge = modelCoef[i];
    totalCharge += charge;
    std::cout << TlUtils::format("[%4d] % 8.3f", i, charge) << std::endl;
  }
  std::cout << TlUtils::format("total: % 8.3f", totalCharge) << std::endl;
}
