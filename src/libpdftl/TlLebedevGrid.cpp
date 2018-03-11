#include "TlLebedevGrid.h"
#include <iostream>

const int TlLebedevGrid::supportedGrids[] = {
    6,   14,  26,  38,  50,  74,  86,  110, 146,  170,  194,
    230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730};

TlLebedevGrid::TlLebedevGrid() {}

void TlLebedevGrid::getGrids(const int numOfGrids,
                             std::vector<TlPosition>* pGrids,
                             std::vector<double>* pWeights) {
  switch (numOfGrids) {
    case 6:
      this->getLeb0006(pGrids, pWeights);
      break;
    case 14:
      this->getLeb0014(pGrids, pWeights);
      break;
    case 26:
      this->getLeb0026(pGrids, pWeights);
      break;
    case 38:
      this->getLeb0038(pGrids, pWeights);
      break;
    case 50:
      this->getLeb0050(pGrids, pWeights);
      break;
    case 74:
      this->getLeb0074(pGrids, pWeights);
      break;
    case 86:
      this->getLeb0086(pGrids, pWeights);
      break;
    case 110:
      this->getLeb0110(pGrids, pWeights);
      break;
    case 146:
      this->getLeb0146(pGrids, pWeights);
      break;
    case 170:
      this->getLeb0170(pGrids, pWeights);
      break;
    case 194:
      this->getLeb0194(pGrids, pWeights);
      break;
    case 230:
      this->getLeb0230(pGrids, pWeights);
      break;
    case 266:
      this->getLeb0266(pGrids, pWeights);
      break;
    case 302:
      this->getLeb0302(pGrids, pWeights);
      break;
    case 350:
      this->getLeb0350(pGrids, pWeights);
      break;
    case 434:
      this->getLeb0434(pGrids, pWeights);
      break;
    case 590:
      this->getLeb0590(pGrids, pWeights);
      break;
    case 770:
      this->getLeb0770(pGrids, pWeights);
      break;
    case 974:
      this->getLeb0974(pGrids, pWeights);
      break;
    case 1202:
      this->getLeb1202(pGrids, pWeights);
      break;
    case 1454:
      this->getLeb1454(pGrids, pWeights);
      break;
    case 1730:
      this->getLeb1730(pGrids, pWeights);
      break;
    default:
      std::cerr << "unsupport grid number in Lebedev grid." << std::endl;
      break;
  }
}

std::vector<int> TlLebedevGrid::getSupportedGridNumber() const {
  const int entries = sizeof(TlLebedevGrid::supportedGrids) /
                      sizeof(TlLebedevGrid::supportedGrids[0]);
  return std::vector<int>(TlLebedevGrid::supportedGrids,
                          TlLebedevGrid::supportedGrids + entries);
}

void TlLebedevGrid::getOh(const int type, double a, double b, double w,
                          std::vector<TlPosition>* pGrids,
                          std::vector<double>* pWeights) {
  std::vector<TlPosition> A;

  switch (type) {
    case 1: {
      A.resize(6);
      a = 1.0;
      A[0] = TlPosition(a, 0.0, 0.0);
      A[1] = TlPosition(-a, 0.0, 0.0);
      A[2] = TlPosition(0.0, a, 0.0);
      A[3] = TlPosition(0.0, -a, 0.0);
      A[4] = TlPosition(0.0, 0.0, a);
      A[5] = TlPosition(0.0, 0.0, -a);
    } break;

    case 2: {
      A.resize(12);
      a = std::sqrt(0.5);
      A[0] = TlPosition(a, a, 0.0);
      A[1] = TlPosition(a, -a, 0.0);
      A[2] = TlPosition(-a, a, 0.0);
      A[3] = TlPosition(-a, -a, 0.0);
      A[4] = TlPosition(a, 0.0, a);
      A[5] = TlPosition(a, 0.0, -a);
      A[6] = TlPosition(-a, 0.0, a);
      A[7] = TlPosition(-a, 0.0, -a);
      A[8] = TlPosition(0.0, a, a);
      A[9] = TlPosition(0.0, a, -a);
      A[10] = TlPosition(0.0, -a, a);
      A[11] = TlPosition(0.0, -a, -a);
    } break;

    case 3: {
      A.resize(8);
      a = std::sqrt(1.0 / 3.0);
      A[0] = TlPosition(a, a, a);
      A[1] = TlPosition(a, a, -a);
      A[2] = TlPosition(a, -a, a);
      A[3] = TlPosition(a, -a, -a);
      A[4] = TlPosition(-a, a, a);
      A[5] = TlPosition(-a, a, -a);
      A[6] = TlPosition(-a, -a, a);
      A[7] = TlPosition(-a, -a, -a);
    } break;

    case 4: {
      A.resize(24);
      b = std::sqrt(1.0 - 2.0 * a * a);
      A[0] = TlPosition(a, a, b);
      A[1] = TlPosition(a, a, -b);
      A[2] = TlPosition(a, -a, b);
      A[3] = TlPosition(a, -a, -b);
      A[4] = TlPosition(-a, a, b);
      A[5] = TlPosition(-a, a, -b);
      A[6] = TlPosition(-a, -a, b);
      A[7] = TlPosition(-a, -a, -b);

      A[8] = TlPosition(a, b, a);
      A[9] = TlPosition(a, -b, a);
      A[10] = TlPosition(a, b, -a);
      A[11] = TlPosition(a, -b, -a);
      A[12] = TlPosition(-a, b, a);
      A[13] = TlPosition(-a, -b, a);
      A[14] = TlPosition(-a, b, -a);
      A[15] = TlPosition(-a, -b, -a);

      A[16] = TlPosition(b, a, a);
      A[17] = TlPosition(-b, a, a);
      A[18] = TlPosition(b, a, -a);
      A[19] = TlPosition(-b, a, -a);
      A[20] = TlPosition(b, -a, a);
      A[21] = TlPosition(-b, -a, a);
      A[22] = TlPosition(b, -a, -a);
      A[23] = TlPosition(-b, -a, -a);
    } break;

    case 5: {
      A.resize(24);
      b = std::sqrt(1.0 - a * a);
      A[0] = TlPosition(a, b, 0.0);
      A[1] = TlPosition(a, -b, 0.0);
      A[2] = TlPosition(-a, b, 0.0);
      A[3] = TlPosition(-a, -b, 0.0);
      A[4] = TlPosition(b, a, 0.0);
      A[5] = TlPosition(-b, a, 0.0);
      A[6] = TlPosition(b, -a, 0.0);
      A[7] = TlPosition(-b, -a, 0.0);

      A[8] = TlPosition(a, 0.0, b);
      A[9] = TlPosition(a, 0.0, -b);
      A[10] = TlPosition(-a, 0.0, b);
      A[11] = TlPosition(-a, 0.0, -b);
      A[12] = TlPosition(b, 0.0, a);
      A[13] = TlPosition(-b, 0.0, a);
      A[14] = TlPosition(b, 0.0, -a);
      A[15] = TlPosition(-b, 0.0, -a);

      A[16] = TlPosition(0.0, a, b);
      A[17] = TlPosition(0.0, a, -b);
      A[18] = TlPosition(0.0, -a, b);
      A[19] = TlPosition(0.0, -a, -b);
      A[20] = TlPosition(0.0, b, a);
      A[21] = TlPosition(0.0, -b, a);
      A[22] = TlPosition(0.0, b, -a);
      A[23] = TlPosition(0.0, -b, -a);
    } break;

    case 6: {
      A.resize(48);
      double c = std::sqrt(1.0 - a * a - b * b);
      A[0] = TlPosition(a, b, c);
      A[1] = TlPosition(a, -b, c);
      A[2] = TlPosition(-a, b, c);
      A[3] = TlPosition(-a, -b, c);
      A[4] = TlPosition(b, a, c);
      A[5] = TlPosition(-b, a, c);
      A[6] = TlPosition(b, -a, c);
      A[7] = TlPosition(-b, -a, c);

      A[8] = TlPosition(a, c, b);
      A[9] = TlPosition(a, c, -b);
      A[10] = TlPosition(-a, c, b);
      A[11] = TlPosition(-a, c, -b);
      A[12] = TlPosition(b, c, a);
      A[13] = TlPosition(-b, c, a);
      A[14] = TlPosition(b, c, -a);
      A[15] = TlPosition(-b, c, -a);

      A[16] = TlPosition(c, a, b);
      A[17] = TlPosition(c, a, -b);
      A[18] = TlPosition(c, -a, b);
      A[19] = TlPosition(c, -a, -b);
      A[20] = TlPosition(c, b, a);
      A[21] = TlPosition(c, -b, a);
      A[22] = TlPosition(c, b, -a);
      A[23] = TlPosition(c, -b, -a);

      A[24] = TlPosition(a, b, -c);
      A[25] = TlPosition(a, -b, -c);
      A[26] = TlPosition(-a, b, -c);
      A[27] = TlPosition(-a, -b, -c);
      A[28] = TlPosition(b, a, -c);
      A[29] = TlPosition(-b, a, -c);
      A[30] = TlPosition(b, -a, -c);
      A[31] = TlPosition(-b, -a, -c);

      A[32] = TlPosition(a, -c, b);
      A[33] = TlPosition(a, -c, -b);
      A[34] = TlPosition(-a, -c, b);
      A[35] = TlPosition(-a, -c, -b);
      A[36] = TlPosition(b, -c, a);
      A[37] = TlPosition(-b, -c, a);
      A[38] = TlPosition(b, -c, -a);
      A[39] = TlPosition(-b, -c, -a);

      A[40] = TlPosition(-c, a, b);
      A[41] = TlPosition(-c, a, -b);
      A[42] = TlPosition(-c, -a, b);
      A[43] = TlPosition(-c, -a, -b);
      A[44] = TlPosition(-c, b, a);
      A[45] = TlPosition(-c, -b, a);
      A[46] = TlPosition(-c, b, -a);
      A[47] = TlPosition(-c, -b, -a);
    } break;
  }

  const int numOfGrids = A.size();
  std::vector<double> W(numOfGrids);
  for (int i = 0; i < numOfGrids; ++i) {
    W[i] = w;
  }

  pGrids->insert(pGrids->end(), A.begin(), A.end());
  pWeights->insert(pWeights->end(), W.begin(), W.end());
}

void TlLebedevGrid::getLeb0006(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.1666666666666667;
  this->getOh(1, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 6);
  assert(pWeights->size() == 6);
}

void TlLebedevGrid::getLeb0014(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.6666666666666667E-1;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.7500000000000000E-1;
  this->getOh(3, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 14);
  assert(pWeights->size() == 14);
}

void TlLebedevGrid::getLeb0026(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.4761904761904762E-1;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.3809523809523810E-1;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.3214285714285714E-1;
  this->getOh(3, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 26);
  assert(pWeights->size() == 26);
}

void TlLebedevGrid::getLeb0038(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.9523809523809524E-2;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.3214285714285714E-1;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.4597008433809831E+0;
  w = 0.2857142857142857E-1;
  this->getOh(5, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 38);
  assert(pWeights->size() == 38);
}

void TlLebedevGrid::getLeb0050(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.1269841269841270E-1;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.2257495590828924E-1;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.2109375000000000E-1;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.3015113445777636E+0;
  w = 0.2017333553791887E-1;
  this->getOh(4, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 50);
  assert(pWeights->size() == 50);
}

void TlLebedevGrid::getLeb0074(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.5130671797338464E-3;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.1660406956574204E-1;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = -0.2958603896103896E-1;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.4803844614152614E+0;
  w = 0.2657620708215946E-1;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3207726489807764E+0;
  w = 0.1652217099371571E-1;
  this->getOh(5, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 74);
  assert(pWeights->size() == 74);
}

void TlLebedevGrid::getLeb0086(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.1154401154401154E-1;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.1194390908585628E-1;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.3696028464541502E+0;
  w = 0.1111055571060340E-1;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6943540066026664E+0;
  w = 0.1187650129453714E-1;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3742430390903412E+0;
  w = 0.1181230374690448E-1;
  this->getOh(5, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 86);
  assert(pWeights->size() == 86);
}

void TlLebedevGrid::getLeb0110(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.3828270494937162E-2;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.9793737512487512E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.1851156353447362E+0;
  w = 0.8211737283191111E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6904210483822922E+0;
  w = 0.9942814891178103E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3956894730559419E+0;
  w = 0.9595471336070963E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4783690288121502E+0;
  w = 0.9694996361663028E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 110);
  assert(pWeights->size() == 110);
}

void TlLebedevGrid::getLeb0146(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.5996313688621381E-3;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.7372999718620756E-2;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.7210515360144488E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.6764410400114264E+0;
  w = 0.7116355493117555E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4174961227965453E+0;
  w = 0.6753829486314477E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1574676672039082E+0;
  w = 0.7574394159054034E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1403553811713183E+0;
  b = 0.4493328323269557E+0;
  w = 0.6991087353303262E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 146);
  assert(pWeights->size() == 146);
}

void TlLebedevGrid::getLeb0170(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.5544842902037365E-2;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.6071332770670752E-2;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.6383674773515093E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.2551252621114134E+0;
  w = 0.5183387587747790E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6743601460362766E+0;
  w = 0.6317929009813725E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4318910696719410E+0;
  w = 0.6201670006589077E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2613931360335988E+0;
  w = 0.5477143385137348E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.4990453161796037E+0;
  b = 0.1446630744325115E+0;
  w = 0.5968383987681156E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 170);
  assert(pWeights->size() == 170);
}

void TlLebedevGrid::getLeb0194(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.1782340447244611E-2;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.5716905949977102E-2;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.5573383178848738E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.6712973442695226E+0;
  w = 0.5608704082587997E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2892465627575439E+0;
  w = 0.5158237711805383E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4446933178717437E+0;
  w = 0.5518771467273614E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1299335447650067E+0;
  w = 0.4106777028169394E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3457702197611283E+0;
  w = 0.5051846064614808E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.1590417105383530E+0;
  b = 0.8360360154824589E+0;
  w = 0.5530248916233094E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 194);
  assert(pWeights->size() == 194);
}

void TlLebedevGrid::getLeb0230(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = -0.5522639919727325E-1;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.4450274607445226E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.4492044687397611E+0;
  w = 0.4496841067921404E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2520419490210201E+0;
  w = 0.5049153450478750E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6981906658447242E+0;
  w = 0.3976408018051883E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6587405243460960E+0;
  w = 0.4401400650381014E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4038544050097660E-1;
  w = 0.1724544350544401E-1;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.5823842309715585E+0;
  w = 0.4231083095357343E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.3545877390518688E+0;
  w = 0.5198069864064399E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.2272181808998187E+0;
  b = 0.4864661535886647E+0;
  w = 0.4695720972568883E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 230);
  assert(pWeights->size() == 230);
}

void TlLebedevGrid::getLeb0266(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = -0.1313769127326952E-2;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = -0.2522728704859336E-2;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.4186853881700583E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.7039373391585475E+0;
  w = 0.5315167977810885E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1012526248572414E+0;
  w = 0.4047142377086219E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4647448726420539E+0;
  w = 0.4112482394406990E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3277420654971629E+0;
  w = 0.3595584899758782E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6620338663699974E+0;
  w = 0.4256131351428158E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.8506508083520399E+0;
  w = 0.4229582700647240E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.3233484542692899E+0;
  b = 0.1153112011009701E+0;
  w = 0.4080914225780505E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2314790158712601E+0;
  b = 0.5244939240922365E+0;
  w = 0.4071467593830964E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 266);
  assert(pWeights->size() == 266);
}

void TlLebedevGrid::getLeb0302(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.8545911725128148E-3;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.3599119285025571E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.3515640345570105E+0;
  w = 0.3449788424305883E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6566329410219612E+0;
  w = 0.3604822601419882E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4729054132581005E+0;
  w = 0.3576729661743367E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.9618308522614784E-1;
  w = 0.2352101413689164E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2219645236294178E+0;
  w = 0.3108953122413675E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.7011766416089545E+0;
  w = 0.3650045807677255E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2644152887060663E+0;
  w = 0.2982344963171804E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.5718955891878961E+0;
  w = 0.3600820932216460E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.2510034751770465E+0;
  b = 0.8000727494073952E+0;
  w = 0.3571540554273387E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1233548532583327E+0;
  b = 0.4127724083168531E+0;
  w = 0.3392312205006170E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 302);
  assert(pWeights->size() == 302);
}

void TlLebedevGrid::getLeb0350(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.3006796749453936E-2;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.3050627745650771E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.7068965463912316E+0;
  w = 0.1621104600288991E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4794682625712025E+0;
  w = 0.3005701484901752E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1927533154878019E+0;
  w = 0.2990992529653774E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6930357961327123E+0;
  w = 0.2982170644107595E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3608302115520091E+0;
  w = 0.2721564237310992E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6498486161496169E+0;
  w = 0.3033513795811141E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1932945013230339E+0;
  w = 0.3007949555218533E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.3800494919899303E+0;
  w = 0.2881964603055307E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.2899558825499574E+0;
  b = 0.7934537856582316E+0;
  w = 0.2958357626535696E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.9684121455103957E-1;
  b = 0.8280801506686862E+0;
  w = 0.3036020026407088E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1833434647041659E+0;
  b = 0.9074658265305127E+0;
  w = 0.2832187403926303E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 350);
  assert(pWeights->size() == 350);
}

void TlLebedevGrid::getLeb0434(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.5265897968224436E-3;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.2548219972002607E-2;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.2512317418927307E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.6909346307509111E+0;
  w = 0.2530403801186355E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1774836054609158E+0;
  w = 0.2014279020918528E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4914342637784746E+0;
  w = 0.2501725168402936E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6456664707424256E+0;
  w = 0.2513267174597564E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2861289010307638E+0;
  w = 0.2302694782227416E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.7568084367178018E-1;
  w = 0.1462495621594614E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3927259763368002E+0;
  w = 0.2445373437312980E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.8818132877794288E+0;
  w = 0.2417442375638981E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.9776428111182649E+0;
  w = 0.1910951282179532E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.2054823696403044E+0;
  b = 0.8689460322872412E+0;
  w = 0.2416930044324775E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5905157048925271E+0;
  b = 0.7999278543857286E+0;
  w = 0.2512236854563495E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5550152361076807E+0;
  b = 0.7717462626915901E+0;
  w = 0.2496644054553086E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.9371809858553722E+0;
  b = 0.3344363145343455E+0;
  w = 0.2236607760437849E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 434);
  assert(pWeights->size() == 434);
}

void TlLebedevGrid::getLeb0590(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.3095121295306187E-3;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.1852379698597489E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.7040954938227469E+0;
  w = 0.1871790639277744E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6807744066455243E+0;
  w = 0.1858812585438317E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6372546939258752E+0;
  w = 0.1852028828296213E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.5044419707800358E+0;
  w = 0.1846715956151242E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4215761784010967E+0;
  w = 0.1818471778162769E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3317920736472123E+0;
  w = 0.1749564657281154E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2384736701421887E+0;
  w = 0.1617210647254411E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1459036449157763E+0;
  w = 0.1384737234851692E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6095034115507196E-1;
  w = 0.9764331165051050E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6116843442009876E+0;
  w = 0.1857161196774078E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.3964755348199858E+0;
  w = 0.1705153996395864E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.1724782009907724E+0;
  w = 0.1300321685886048E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.5610263808622060E+0;
  b = 0.3518280927733519E+0;
  w = 0.1842866472905286E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4742392842551980E+0;
  b = 0.2634716655937950E+0;
  w = 0.1802658934377451E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5984126497885380E+0;
  b = 0.1816640840360209E+0;
  w = 0.1849830560443660E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3791035407695563E+0;
  b = 0.1720795225656878E+0;
  w = 0.1713904507106709E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2778673190586244E+0;
  b = 0.8213021581932511E-1;
  w = 0.1555213603396808E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5033564271075117E+0;
  b = 0.8999205842074875E-1;
  w = 0.1802239128008525E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 590);
  assert(pWeights->size() == 590);
}

void TlLebedevGrid::getLeb0770(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.2192942088181184E-3;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.1436433617319080E-2;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.1421940344335877E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.5087204410502360E-1;
  w = 0.6798123511050502E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1228198790178831E+0;
  w = 0.9913184235294912E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2026890814408786E+0;
  w = 0.1180207833238949E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2847745156464294E+0;
  w = 0.1296599602080921E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3656719078978026E+0;
  w = 0.1365871427428316E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4428264886713469E+0;
  w = 0.1402988604775325E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.5140619627249735E+0;
  w = 0.1418645563595609E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6306401219166803E+0;
  w = 0.1421376741851662E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6716883332022612E+0;
  w = 0.1423996475490962E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6979792685336881E+0;
  w = 0.1431554042178567E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1446865674195309E+0;
  w = 0.9254401499865368E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.3390263475411216E+0;
  w = 0.1250239995053509E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.5335804651263506E+0;
  w = 0.1394365843329230E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.6944024393349413E-1;
  b = 0.2355187894242326E+0;
  w = 0.1127089094671749E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2269004109529460E+0;
  b = 0.4102182474045730E+0;
  w = 0.1345753760910670E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.8025574607775339E-1;
  b = 0.6214302417481605E+0;
  w = 0.1424957283316783E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1467999527896572E+0;
  b = 0.3245284345717394E+0;
  w = 0.1261523341237750E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1571507769824727E+0;
  b = 0.5224482189696630E+0;
  w = 0.1392547106052696E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2365702993157246E+0;
  b = 0.6017546634089558E+0;
  w = 0.1418761677877656E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.7714815866765732E-1;
  b = 0.4346575516141163E+0;
  w = 0.1338366684479554E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3062936666210730E+0;
  b = 0.4908826589037616E+0;
  w = 0.1393700862676131E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3822477379524787E+0;
  b = 0.5648768149099500E+0;
  w = 0.1415914757466932E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 770);
  assert(pWeights->size() == 770);
}

void TlLebedevGrid::getLeb0974(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.1438294190527431E-3;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.1125772288287004E-2;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.4292963545341347E-1;
  w = 0.4948029341949241E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1051426854086404E+0;
  w = 0.7357990109125470E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1750024867623087E+0;
  w = 0.8889132771304384E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2477653379650257E+0;
  w = 0.9888347838921435E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3206567123955957E+0;
  w = 0.1053299681709471E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3916520749849983E+0;
  w = 0.1092778807014578E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4590825874187624E+0;
  w = 0.1114389394063227E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.5214563888415861E+0;
  w = 0.1123724788051555E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6253170244654199E+0;
  w = 0.1125239325243814E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6637926744523170E+0;
  w = 0.1126153271815905E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6910410398498301E+0;
  w = 0.1130286931123841E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.7052907007457760E+0;
  w = 0.1134986534363955E-2;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1236686762657990E+0;
  w = 0.6823367927109931E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.2940777114468387E+0;
  w = 0.9454158160447096E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.4697753849207649E+0;
  w = 0.1074429975385679E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.6334563241139567E+0;
  w = 0.1129300086569132E-2;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.5974048614181342E-1;
  b = 0.2029128752777523E+0;
  w = 0.8436884500901954E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1375760408473636E+0;
  b = 0.4602621942484054E+0;
  w = 0.1075255720448885E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3391016526336286E+0;
  b = 0.5030673999662036E+0;
  w = 0.1108577236864462E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1271675191439820E+0;
  b = 0.2817606422442134E+0;
  w = 0.9566475323783357E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2693120740413512E+0;
  b = 0.4331561291720157E+0;
  w = 0.1080663250717391E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1419786452601918E+0;
  b = 0.6256167358580814E+0;
  w = 0.1126797131196295E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.6709284600738255E-1;
  b = 0.3798395216859157E+0;
  w = 0.1022568715358061E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.7057738183256172E-1;
  b = 0.5517505421423520E+0;
  w = 0.1108960267713108E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2783888477882155E+0;
  b = 0.6029619156159187E+0;
  w = 0.1122790653435766E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1979578938917407E+0;
  b = 0.3589606329589096E+0;
  w = 0.1032401847117460E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2087307061103274E+0;
  b = 0.5348666438135476E+0;
  w = 0.1107249382283854E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4055122137872836E+0;
  b = 0.5674997546074373E+0;
  w = 0.1121780048519972E-2;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 974);
  assert(pWeights->size() == 974);
}

void TlLebedevGrid::getLeb1202(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.1105189233267572E-3;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.9205232738090741E-3;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.9133159786443561E-3;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.3712636449657089E-1;
  w = 0.3690421898017899E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.9140060412262223E-1;
  w = 0.5603990928680660E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1531077852469906E+0;
  w = 0.6865297629282609E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2180928891660612E+0;
  w = 0.7720338551145630E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2839874532200175E+0;
  w = 0.8301545958894795E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3491177600963764E+0;
  w = 0.8686692550179628E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4121431461444309E+0;
  w = 0.8927076285846890E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4718993627149127E+0;
  w = 0.9060820238568219E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.5273145452842337E+0;
  w = 0.9119777254940867E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6209475332444019E+0;
  w = 0.9128720138604181E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6569722711857291E+0;
  w = 0.9130714935691735E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6841788309070143E+0;
  w = 0.9152873784554116E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.7012604330123631E+0;
  w = 0.9187436274321654E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1072382215478166E+0;
  w = 0.5176977312965694E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.2582068959496968E+0;
  w = 0.7331143682101417E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.4172752955306717E+0;
  w = 0.8463232836379928E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.5700366911792503E+0;
  w = 0.9031122694253992E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.9827986018263947E+0;
  b = 0.1771774022615325E+0;
  w = 0.6485778453163257E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.9624249230326228E+0;
  b = 0.2475716463426288E+0;
  w = 0.7435030910982369E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.9402007994128811E+0;
  b = 0.3354616289066489E+0;
  w = 0.7998527891839054E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.9320822040143202E+0;
  b = 0.3173615246611977E+0;
  w = 0.8101731497468018E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.9043674199393299E+0;
  b = 0.4090268427085357E+0;
  w = 0.8483389574594331E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.8912407560074747E+0;
  b = 0.3854291150669224E+0;
  w = 0.8556299257311812E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.8676435628462708E+0;
  b = 0.4932221184851285E+0;
  w = 0.8803208679738260E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.8581979986041619E+0;
  b = 0.4785320675922435E+0;
  w = 0.8811048182425720E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.8396753624049856E+0;
  b = 0.4507422593157064E+0;
  w = 0.8850282341265444E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.8165288564022188E+0;
  b = 0.5632123020762100E+0;
  w = 0.9021342299040653E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.8015469370783529E+0;
  b = 0.5434303569693900E+0;
  w = 0.9010091677105086E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.7773563069070351E+0;
  b = 0.5123518486419871E+0;
  w = 0.9022692938426915E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.7661621213900394E+0;
  b = 0.6394279634749102E+0;
  w = 0.9158016174693465E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.7553584143533510E+0;
  b = 0.6269805509024392E+0;
  w = 0.9131578003189435E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.7344305757559503E+0;
  b = 0.6031161693096310E+0;
  w = 0.9107813579482705E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.7043837184021765E+0;
  b = 0.5693702498468441E+0;
  w = 0.9105760258970126E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 1202);
  assert(pWeights->size() == 1202);
}

void TlLebedevGrid::getLeb1454(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.7777160743261247E-4;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.7557646413004701E-3;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.3229290663413854E-1;
  w = 0.2841633806090617E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.8036733271462222E-1;
  w = 0.4374419127053555E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1354289960531653E+0;
  w = 0.5417174740872172E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1938963861114426E+0;
  w = 0.6148000891358593E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2537343715011275E+0;
  w = 0.6664394485800705E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3135251434752570E+0;
  w = 0.7025039356923220E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3721558339375338E+0;
  w = 0.7268511789249627E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4286809575195696E+0;
  w = 0.7422637534208629E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4822510128282994E+0;
  w = 0.7509545035841214E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.5320679333566263E+0;
  w = 0.7548535057718401E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6172998195394274E+0;
  w = 0.7554088969774001E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6510679849127481E+0;
  w = 0.7553147174442808E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6777315251687360E+0;
  w = 0.7564767653292297E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6963109410648741E+0;
  w = 0.7587991808518730E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.7058935009831749E+0;
  w = 0.7608261832033027E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.9955546194091857E+0;
  w = 0.4021680447874916E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.9734115901794209E+0;
  w = 0.5804871793945964E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.9275693732388626E+0;
  w = 0.6792151955945159E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.8568022422795103E+0;
  w = 0.7336741211286294E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.7623495553719372E+0;
  w = 0.7581866300989608E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.5707522908892223E+0;
  b = 0.4387028039889501E+0;
  w = 0.7538257859800743E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5196463388403083E+0;
  b = 0.3858908414762617E+0;
  w = 0.7483517247053123E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4646337531215351E+0;
  b = 0.3301937372343854E+0;
  w = 0.7371763661112059E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4063901697557691E+0;
  b = 0.2725423573563777E+0;
  w = 0.7183448895756934E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3456329466643087E+0;
  b = 0.2139510237495250E+0;
  w = 0.6895815529822191E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2831395121050332E+0;
  b = 0.1555922309786647E+0;
  w = 0.6480105801792886E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2197682022925330E+0;
  b = 0.9892878979686097E-1;
  w = 0.5897558896594636E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1564696098650355E+0;
  b = 0.4598642910675510E-1;
  w = 0.5095708849247346E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.6027356673721295E+0;
  b = 0.3376625140173426E+0;
  w = 0.7536906428909755E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5496032320255096E+0;
  b = 0.2822301309727988E+0;
  w = 0.7472505965575118E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4921707755234567E+0;
  b = 0.2248632342592540E+0;
  w = 0.7343017132279698E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4309422998598483E+0;
  b = 0.1666224723456479E+0;
  w = 0.7130871582177445E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3664108182313672E+0;
  b = 0.1086964901822169E+0;
  w = 0.6817022032112776E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2990189057758436E+0;
  b = 0.5251989784120085E-1;
  w = 0.6380941145604121E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.6268724013144998E+0;
  b = 0.2297523657550023E+0;
  w = 0.7550381377920310E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5707324144834607E+0;
  b = 0.1723080607093800E+0;
  w = 0.7478646640144802E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5096360901960365E+0;
  b = 0.1140238465390513E+0;
  w = 0.7335918720601220E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4438729938312456E+0;
  b = 0.5611522095882537E-1;
  w = 0.7110120527658118E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.6419978471082389E+0;
  b = 0.1164174423140873E+0;
  w = 0.7571363978689501E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5817218061802611E+0;
  b = 0.5797589531445219E-1;
  w = 0.7489908329079234E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 1454);
  assert(pWeights->size() == 1454);
}

void TlLebedevGrid::getLeb1730(std::vector<TlPosition>* pGrids,
                               std::vector<double>* pWeights) {
  pGrids->clear();
  pWeights->clear();

  double a = 0.0;
  double b = 0.0;
  double w = 0.0;

  w = 0.6309049437420976E-4;
  this->getOh(1, a, b, w, pGrids, pWeights);
  w = 0.6398287705571748E-3;
  this->getOh(2, a, b, w, pGrids, pWeights);
  w = 0.6357185073530720E-3;
  this->getOh(3, a, b, w, pGrids, pWeights);
  a = 0.2860923126194662E-1;
  w = 0.2221207162188168E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.7142556767711522E-1;
  w = 0.3475784022286848E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1209199540995559E+0;
  w = 0.4350742443589804E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.1738673106594379E+0;
  w = 0.4978569136522127E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2284645438467734E+0;
  w = 0.5435036221998053E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.2834807671701512E+0;
  w = 0.5765913388219542E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3379680145467339E+0;
  w = 0.6001200359226003E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.3911355454819537E+0;
  w = 0.6162178172717512E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4422860353001403E+0;
  w = 0.6265218152438485E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.4907781568726057E+0;
  w = 0.6323987160974212E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.5360006153211468E+0;
  w = 0.6350767851540569E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6142105973596603E+0;
  w = 0.6354362775297107E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6459300387977504E+0;
  w = 0.6352302462706235E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6718056125089225E+0;
  w = 0.6358117881417972E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.6910888533186254E+0;
  w = 0.6373101590310117E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.7030467416823252E+0;
  w = 0.6390428961368665E-3;
  this->getOh(4, a, b, w, pGrids, pWeights);
  a = 0.8354951166354646E-1;
  w = 0.3186913449946576E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.2050143009099486E+0;
  w = 0.4678028558591711E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.3370208290706637E+0;
  w = 0.5538829697598626E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.4689051484233963E+0;
  w = 0.6044475907190476E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.5939400424557334E+0;
  w = 0.6313575103509012E-3;
  this->getOh(5, a, b, w, pGrids, pWeights);
  a = 0.1394983311832261E+0;
  b = 0.4097581162050343E-1;
  w = 0.4078626431855630E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.1967999180485014E+0;
  b = 0.8851987391293348E-1;
  w = 0.4759933057812725E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2546183732548967E+0;
  b = 0.1397680182969819E+0;
  w = 0.5268151186413440E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3121281074713875E+0;
  b = 0.1929452542226526E+0;
  w = 0.5643048560507316E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3685981078502492E+0;
  b = 0.2467898337061562E+0;
  w = 0.5914501076613073E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4233760321547856E+0;
  b = 0.3003104124785409E+0;
  w = 0.6104561257874195E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4758671236059246E+0;
  b = 0.3526684328175033E+0;
  w = 0.6230252860707806E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5255178579796463E+0;
  b = 0.4031134861145713E+0;
  w = 0.6305618761760796E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5718025633734589E+0;
  b = 0.4509426448342351E+0;
  w = 0.6343092767597889E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.2686927772723415E+0;
  b = 0.4711322502423248E-1;
  w = 0.5176268945737826E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3306006819904809E+0;
  b = 0.9784487303942695E-1;
  w = 0.5564840313313692E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.3904906850594983E+0;
  b = 0.1505395810025273E+0;
  w = 0.5856426671038980E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4479957951904390E+0;
  b = 0.2039728156296050E+0;
  w = 0.6066386925777091E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5027076848919780E+0;
  b = 0.2571529941121107E+0;
  w = 0.6208824962234458E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5542087392260217E+0;
  b = 0.3092191375815670E+0;
  w = 0.6296314297822907E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.6020850887375187E+0;
  b = 0.3593807506130276E+0;
  w = 0.6340423756791859E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4019851409179594E+0;
  b = 0.5063389934378671E-1;
  w = 0.5829627677107342E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.4635614567449800E+0;
  b = 0.1032422269160612E+0;
  w = 0.6048693376081110E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5215860931591575E+0;
  b = 0.1566322094006254E+0;
  w = 0.6202362317732461E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5758202499099271E+0;
  b = 0.2098082827491099E+0;
  w = 0.6299005328403779E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.6259893683876795E+0;
  b = 0.2618824114553391E+0;
  w = 0.6347722390609353E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5313795124811891E+0;
  b = 0.5263245019338556E-1;
  w = 0.6203778981238834E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.5893317955931995E+0;
  b = 0.1061059730982005E+0;
  w = 0.6308414671239979E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.6426246321215801E+0;
  b = 0.1594171564034221E+0;
  w = 0.6362706466959498E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);
  a = 0.6511904367376113E+0;
  b = 0.5354789536565540E-1;
  w = 0.6375414170333233E-3;
  this->getOh(6, a, b, w, pGrids, pWeights);

  assert(pGrids->size() == 1730);
  assert(pWeights->size() == 1730);
}
