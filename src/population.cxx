#include "ProteinDF.h"
#include "DfPopulation.h"

#include "TlGetopt.h"
#include "TlParameter.h"

int main(int argc, char *argv[]) {
  TlGetopt opt(argc, argv, "i:");

  TlParameter flGbi;
  flGbi.load("fl_Input/fl_Globalinput");

  const int nIteration = std::atoi(opt["i"].c_str());
  
  DfPopulation dfPop(flGbi, nIteration);
  dfPop.DfPopuMain();

  return 0;
}

