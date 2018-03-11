#ifndef TLLEBEDEVGRID_H
#define TLLEBEDEVGRID_H

#include <vector>
#include "TlPosition.h"

class TlLebedevGrid {
 public:
  TlLebedevGrid();

  void getGrids(const int numOfGrids, std::vector<TlPosition>* pGrids,
                std::vector<double>* pWeights);

  std::vector<int> getSupportedGridNumber() const;

 private:
  void getOh(const int type, double a, double b, double w,
             std::vector<TlPosition>* pGrids, std::vector<double>* pWeights);

  void getLeb0006(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0014(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0026(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0038(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0050(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0074(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0086(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0110(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0146(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0170(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0194(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0230(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0266(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0302(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0350(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0434(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0590(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0770(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb0974(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb1202(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb1454(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);
  void getLeb1730(std::vector<TlPosition>* pGrids,
                  std::vector<double>* pWeights);

 private:
  static const int supportedGrids[];
};

#endif  // TLLEBEDEVGRID_H
