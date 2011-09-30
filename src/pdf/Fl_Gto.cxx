#include <cstdlib>
#include <fstream>
#include <string>
#include <cassert>

#include "CnError.h"
#include "FileX.h"
#include "Fl_Gto.h"
#include "TlMath.h"
#include "TlUtils.h"
#include "TlLogX.h"
#include "TlStringTokenizer.h"
#include "TlFile.h"
#include "PdfUtils.h"
#include "TlPseudoYaml.h"

Fl_Gto::Fl_Gto(const std::string& str) : m_sFileName(str), m_isUpdate(false)
{
    if ((this->m_sFileName != "") &&
        (TlFile::isExist(this->m_sFileName) == true)) {
        this->load();
    }
}


Fl_Gto::Fl_Gto(const TlSerializeData& gtoData) : m_isUpdate(false)
{
    this->setup(gtoData);
}


Fl_Gto::Fl_Gto(const Fl_Gto& rhs)
    : cgto(rhs.cgto), m_sFileName(rhs.m_sFileName), m_isUpdate(false) {
}

Fl_Gto::~Fl_Gto()
{
    if ((this->m_sFileName != "") &&
        (this->m_isUpdate == true)) {
        this->save();
    }
}

Fl_Gto Fl_Gto::operator=(const Fl_Gto& rhs)
{
    if (this != &rhs) {
        this->cgto = rhs.cgto;
        this->m_sFileName = rhs.m_sFileName;
        this->m_isUpdate = false;
    }

    return *this;
}

int Fl_Gto::getSnum(int i) const
{
    return this->cgto[i].Snum;
}

int Fl_Gto::getPnum(int i) const
{
    return this->cgto[i].Pnum;
}

int Fl_Gto::getDnum(int i) const
{
    return this->cgto[i].Dnum;
}

int Fl_Gto::getNumOfCGTOs() const
{
    return this->cgto.size();
}

std::string Fl_Gto::getAtom(int s) const
{
    if (s < 0) {
// #ifndef NDEBUG
//     std::cerr << "Assertion! (s>0) at Fl_Gto::getAtom(). s = " << s << std::endl;
// #endif // NDEBUG
        return "";
    }

    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));
    return this->cgto[s].atom;
}

std::string Fl_Gto::getLabel(int s) const
{
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));

    return this->cgto[s].label;
}

std::string Fl_Gto::getShellname(int s) const
{
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));

    return this->cgto[s].shellname;
}

char Fl_Gto::getShell(int s) const
{
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));

    return this->cgto[s].shell;
}

double Fl_Gto::getScalefactor(int s) const
{
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));

    return this->cgto[s].scalefactor;
}

int Fl_Gto::getContraction(int s) const
{
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));

    return this->cgto[s].pgto.size();
}

double Fl_Gto::getExponent(int s, int p) const
{
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));
    assert((0 <= p) && (p < static_cast<int>(this->cgto[s].contraction())));

    return this->cgto[s].pgto[p].exponent;
}

double Fl_Gto::getCoefficient(int s, int p) const
{
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));
    assert((0 <= p) && (p < static_cast<int>(this->cgto[s].contraction())));

    return this->cgto[s].pgto[p].coefficient;
}

// y番目が属している基底関数の俗称を与える。
std::string Fl_Gto::getBasisName(int y) const
{
    assert((0 <= y) && (y < static_cast<int>(this->cgto.size())));

    return this->cgto[y].basisName;
}

void Fl_Gto::set(int i, const Cgto& data)
{
    assert(0 <= i);

    if (static_cast<int>(this->cgto.size()) <= i) {
        this->cgto.resize(i +1);
    }

    this->cgto[i] = data;
    this->m_isUpdate = true;
}

void Fl_Gto::push_back(const Cgto& data)
{
    this->cgto.push_back(data);
    this->m_isUpdate = true;
}


// for serialize =============================================================
//
void Fl_Gto::save(const std::string& newSavePath)
{
    if (newSavePath != "") {
        this->m_sFileName = newSavePath;
    }
    const std::string savePath = this->m_sFileName;
    
    std::ofstream fs(savePath.c_str(), std::ios::out | std::ios::trunc);
    if (fs.is_open() != true) {
        CnErr.abort(TlUtils::format("could not open file \"%s\" at Fl_Gto::write()", this->m_sFileName.c_str()));
    }

    fs.precision(16);
    fs.setf(std::ios::scientific, std::ios::floatfield);

    fs << " ------------------------------------------------------------------------\n";
    fs << " numbercgto     : " << this->getNumOfCGTOs()     << "\n";
    fs << " ------------------------------------------------------------------------\n";
    for (int i = 0; i < this->getNumOfCGTOs(); ++i) {
        fs << " --- " << i << " th CGTO --------------------\n";
        fs << "Snum [" << this->cgto[i].Snum << "]\n";
        fs << "Pnum [" << this->cgto[i].Pnum << "]\n";
        fs << "Dnum [" << this->cgto[i].Dnum << "]\n";
        fs << "basis            [" << cgto[i].basisName        << "]\n";
        fs << "atom             [" << cgto[i].atom             << "]\n";
        fs << "label            [" << cgto[i].label            << "]\n";
        fs << "shellname        [" << cgto[i].shellname        << "]\n";
        fs << "shell             " << cgto[i].shell            << " \n";
        fs << "scalefactor       " << cgto[i].scalefactor      << " \n";
        fs << "contraction       " << cgto[i].contraction()      << " \n";
        fs << " // primitive GTO //\n";

        for (int j = 0; j < cgto[i].contraction(); ++j) {
            fs << this->cgto[i].pgto[j].exponent   << "  "
            << this->cgto[i].pgto[j].coefficient << "\n";
        }
    }
    fs << " ------------------------------------------------------------------------\n";

    fs.close();
}

void Fl_Gto::showAMOSS(const std::string& title) const
{
    TlLogX& Log = TlLogX::getInstance();
    Log << "\n\n";
    Log << "  " << title << "\n";
    Log << "\n";

    std::string previous_atom  = "-1";
    std::string previous_label = "-1";

    for (int i=0; i < this->getNumOfCGTOs(); i++) {
        if (previous_atom != cgto[i].atom || previous_label != cgto[i].label) {
            Log << "\n";
            Log << TlUtils::format("%5s %5s = NULL ;; %s\n",
                                   cgto[i].atom.c_str(), cgto[i].label.c_str(), cgto[i].basisName.c_str());

            previous_atom  = cgto[i].atom;
            previous_label = cgto[i].label;
        }

        Log << "      &  " << cgto[i].shell << TlUtils::format("  ( %3.1lf )\n", cgto[i].scalefactor);

        for (int j = 0; j < cgto[i].contraction(); j++) {
            Log << TlUtils::format("      %18.8lE %18.8lE\n",
                                   this->cgto[i].pgto[j].exponent, this->cgto[i].pgto[j].coefficient);
        }
    }
}

void Fl_Gto::showGAMESS() const
{
    TlLogX& Log = TlLogX::getInstance();
    std::string prevAtom  = "-1";
    std::string prevBasis = "-1";

    for (int i = 0; i < this->getNumOfCGTOs(); ++i) {
        if (prevAtom   != this->cgto[i].atom ||
                prevBasis  != this->cgto[i].basisName) {
            Log << "\n";
            Log << TlUtils::format("%2s %s\n", cgto[i].atom.c_str(), cgto[i].basisName.c_str());
            prevAtom  = cgto[i].atom;
            prevBasis = cgto[i].basisName;
        }

        const char sShell       = cgto[i].shell;
        const int nContraction = cgto[i].contraction();
        Log << TlUtils::format("   %c %d\n", sShell, nContraction);

        for (int j = 0; j < cgto[i].contraction(); j++) {
            Log << TlUtils::format(" %2d %18.8lE %18.8lE\n",
                                   j+1, this->cgto[i].pgto[j].exponent, this->cgto[i].pgto[j].coefficient);
        }
    }
}

void Fl_Gto::show(const std::string& title) const
{
    TlLogX& Log = TlLogX::getInstance();
    Log << "\n\n";
    Log << "  " << title << "\n";
    Log << "\n";
    this->show();
}


void Fl_Gto::show() const
{
    TlLogX& Log = TlLogX::getInstance();
    Log << "\n\n";
    int angular;

    Log << " ------------------------------------------------------------------------\n";
//   Log << " MaxContraction : " << int(MaxContraction) << "\n";
//   Log << " MaxCgto        : " << int(MaxCgto)        << "\n";
    Log << " numbercgto     : " << this->getNumOfCGTOs()<< "\n";
    Log << " ------------------------------------------------------------------------\n";
    for (int i=0; i < this->getNumOfCGTOs(); i++) {
        Log << " --- " << i << " th CGTO --------------------\n";
        Log << "Snum [" << this->cgto[i].Snum << "]\n";
        Log << "Pnum [" << this->cgto[i].Pnum << "]\n";
        Log << "Dnum [" << this->cgto[i].Dnum << "]\n";
        Log << "basis             [" << cgto[i].basisName        << "]\n";
        Log << "atom              [" << cgto[i].atom             << "]\n";
        Log << "label             [" << cgto[i].label            << "]\n";
        Log << "shellname         [" << cgto[i].shellname        << "]\n";
        Log << "shell              " << cgto[i].shell            << " \n";

        Log << TlUtils::format("  scale factor     = %18.8lf\n", cgto[i].scalefactor);
        Log << TlUtils::format("  contraction      = %18ld\n",   cgto[i].contraction());
        Log << "  <normalizedfactor> ";

        if (getShell(i) == 's') {
            angular=0;
        } else if (getShell(i) == 'p') {
            angular=1;
        } else if (getShell(i) == 'd') {
            angular=2;
        } else if (getShell(i) == 'f') {
            angular=3;
        } else {
            Log << " *** Fl_Gto::show() not supported such large shell.\n";
            Log << "     +++ calculate in case of l=m=n=0, and continue.\n";
            angular=0;
        }
        for (int l=0; l<=angular; l++) {
            for (int m=0; m<=angular; m++) {
                for (int n=0; n<=angular; n++) {
                    if (l+m+n == angular) {
                        Log << TlUtils::format("  (%1ld,%1ld,%1ld) = %18.8lE ", l, m, n, getNormalizedfactor(i, l,m,n));
                    }
                }
            }
        }
        Log << "\n";
        Log << "    // primitive GTO //\n";

        for (int j=0; j<cgto[i].contraction(); j++) {
            Log << TlUtils::format("      %18.8lE %18.8lE ",
                                   this->cgto[i].pgto[j].exponent, this->cgto[i].pgto[j].coefficient);
            Log << "  <normalized> ";

            if (getShell(i)=='s') angular=0;
            else if (getShell(i)=='p') angular=1;
            else if (getShell(i)=='d') angular=2;
            else if (getShell(i)=='f') angular=3;
            else {
                Log << " *** Fl_Gto::show() not supported such large shell.\n";
                Log << "     +++ calculate in case of l=m=n=0, and continue.\n";
                angular=0;
            }
            for (int l=0; l <= angular; l++) {
                for (int m=0; m <= angular; m++) {
                    for (int n=0; n <= angular; n++) {
                        if (l +m +n == angular) {
                            Log << TlUtils::format("  (%1ld,%1ld,%1ld) = %18.8lE %18.8lE",
                                                   l, m, n,
                                                   getNormalized(i,j, l,m,n),
                                                   getCoulombnormalized(i,j, l,m,n)
                                                  );
                        }
                    }
                }
            }
            Log << "\n";
        }
    }
    Log << " ------------------------------------------------------------------------\n";
}

void Fl_Gto::load()
{
    if (FileX::isExist(this->m_sFileName) == false) {
        std::cerr << this->m_sFileName << " is not found." << std::endl;
        //return false;
        return;
    }

    std::ifstream fs(this->m_sFileName.c_str(), std::ios::in);
    if (fs.fail()) {
        fs.clear(); // if stream at eof, must clear error.
    }
    fs.seekg(0, std::ios::beg);

    // read header
    int numOfCGTOs = -1;
    {
        std::string line;
        while (std::getline(fs, line)) {
            if (PdfUtils::isComment(line) == true) {
                continue;
            }

            TlUtils::trim_ws(line);

            TlStringTokenizer token(line);
            if (token.countTokens() == 3) {
                std::string term1 = token.nextToken();
                std::string term2 = token.nextToken();
                std::string term3 = token.nextToken();
                if (term1 == "numbercgto") {
                    // ex)  numbercgto     : 17
                    numOfCGTOs = std::atoi(term3.c_str());
                    break;
                }
            }
        }
    }
    if (numOfCGTOs == -1) {
        CnErr.abort(TlUtils::format("not found numbercgto line at %s", this->m_sFileName.c_str()));
    }

    // read body
    this->cgto.resize(numOfCGTOs);
    for (int c = 0; c < numOfCGTOs; ++c) {
        Cgto cgto;

        bool isReadPGTO = false;
        int contraction = 0;
        int currentContraction = 0;
        std::string line;
        while (std::getline(fs, line)) {
            if (PdfUtils::isComment(line) == true) {
                continue;
            }
            TlUtils::trim_ws(line);

            std::string term1 = TlUtils::getPdfParam(line);
            std::string term2 = TlUtils::getPdfParam(line);

            if (isReadPGTO != true) {
                if (term1 == "Snum") {
                    cgto.Snum = std::atoi(TlUtils::getPdfParam(term2).c_str());
                } else if (term1 == "Pnum") {
                    cgto.Pnum = std::atoi(TlUtils::getPdfParam(term2).c_str());
                } else if (term1 == "Dnum") {
                    cgto.Dnum = std::atoi(TlUtils::getPdfParam(term2).c_str());
                } else if (term1 == "basis") {
                    cgto.basisName = term2;
                } else if (term1 == "atom") {
                    cgto.atom = TlUtils::getPdfParam(term2);
                } else if (term1 == "label") {
                    cgto.label = TlUtils::getPdfParam(term2);
                } else if (term1 == "shellname") {
                    cgto.shellname = TlUtils::getPdfParam(term2);
                } else if (term1 == "shell") {
                    cgto.shell = TlUtils::getPdfParam(term2)[0];
                } else if (term1 == "scalefactor") {
                    cgto.scalefactor = std::atof(TlUtils::getPdfParam(term2).c_str());
                } else if (term1 == "contraction") {
                    contraction = std::atoi(TlUtils::getPdfParam(term2).c_str());
                    cgto.pgto.resize(contraction);
                    currentContraction = 0;
                    isReadPGTO = true;
                } else {
                    // error!
                }
            } else {
                cgto.pgto[currentContraction].exponent = std::atof(term1.c_str());
                cgto.pgto[currentContraction].coefficient = std::atof(term2.c_str());

                ++currentContraction;
                if (currentContraction == contraction) {
                    break;
                }
            }
        }

        this->cgto[c] = cgto;
    }
}


void Fl_Gto::setup(const TlSerializeData& basisSet)
{
//     TlPseudoYaml yaml(basisSet);
//     std::cerr << yaml.str() << std::endl;
    this->cgto.clear();
    
    //const TlSerializeData& basisSet = data["basis_set"];
    TlSerializeData::MapConstIterator pEnd = basisSet.endMap();
    for (TlSerializeData::MapConstIterator p = basisSet.beginMap(); p != pEnd; ++p) {
        const TlSerializeData contents = p->second;

        std::string atomLabel = p->first.getStr();
//      std::cerr << "atom label = " << atomLabel << std::endl;
        const std::string::size_type pos = atomLabel.find('@', 0);
        const std::string atom = atomLabel.substr(0, pos);
        std::string label = "";
        if ((pos != std::string::npos) &&
            (pos +1) <= atomLabel.size()) {
            label = atomLabel.substr(pos +1);
        }

        const std::string basisName = contents["name"].getStr();
        
        TlSerializeData::ArrayConstIterator qEnd = contents["cGTOs"].endArray();
        for (TlSerializeData::ArrayConstIterator q = contents["cGTOs"].beginArray(); q != qEnd; ++q) {
            Cgto cgto;
            cgto.basisName = basisName;
            cgto.atom = atom;
            cgto.label = label;
            cgto.scalefactor = (*q)["scale_factor"].getDouble();
            std::string shellStr = (*q)["shell_type"].getStr();
            if (shellStr.size() > 0) {
                cgto.shell = shellStr[0];
            }

//             std::cerr << TlUtils::format("scale_factor=%e, shell=%c", cgto.scalefactor, cgto.shell)
//                       << std::endl;
            
            const TlSerializeData& pGTOs = (*q)["pGTOs"];
            TlSerializeData::ArrayConstIterator rEnd = pGTOs.endArray();
            for (TlSerializeData::ArrayConstIterator r = pGTOs.beginArray(); r != rEnd; ++r) {
                double coef = (*r)["coefficient"].getDouble();
                double exp = (*r)["exponent"].getDouble();
                Cgto::Pgto pgto;
                pgto.exponent = exp;
                pgto.coefficient = coef;
                cgto.pgto.push_back(pgto);

//                 std::cerr << TlUtils::format("coef=%e, exp=%e", coef, exp)
//                           << std::endl;
            }

            this->push_back(cgto);
        }
    }
}


// int Fl_Gto::getNumOfCGTOs(const std::string& atomSymbol,
//                           const std::string& label) const
// {
//     const int numOfCGTOs = this->getNumOfCGTOs();
//     int count = 0;
//     for (int i = 0; i < numOfCGTOs; ++i) {
//         if ((this->cgto[i].atom == atomSymbol) &&
//             (this->cgto[i].label = label)) {
//             ++count;
//         }
//     }

//     return count;
// }


// std::string Fl_Gto::getBasisName(const std::string& atomSymbol,
//                                  const std::string& label,
//                                  int cGTO_index) const
// {
// }


// normalized factor =========================================================
//
double Fl_Gto::getNormalizedfactor(int s, int l, int m, int n) const
{
    TlLogX& Log = TlLogX::getInstance();

    // index check for s.
    if (s > this->getNumOfCGTOs() -1) {
        Log << " *** Fl_Gto::getNormalizedfactor()\n";
        Log << "     force to load(get) basis data over avairable maximun(getNumbercgto()), Log of Range.\n";
        Log << "     this is program error.\n";
        CnErr.abort();
    }

    // index chech l, m, n.
    if (l<0 || m<0 || n<0) {
        Log << " *** Fl_Gto::getNormalizedfactor() illegal angular momentum.\n";
        CnErr.abort();
    }

    // chech consistent with shell.
    int maxangular = 0;
    switch (getShell(s)) {
    case 's':
        maxangular=0;
        break;
    case 'p':
        maxangular=1;
        break;
    case 'd':
        maxangular=2;
        break;
    case 'f':
        maxangular=3;
        break;
    default:
        Log << " *** Fl_Gto::getNormalizedfactor() not support over f angular momentum, sorry.\n";
        CnErr.abort();
        break;
    }

    if (l+m+n>maxangular) {
        Log << " *** Fl_Gto::getNormalizedfactor() illegal angular momentum, inconsistent with shell.\n";
        CnErr.abort();
    }

    // calculate normalized factor of CGTO, of which angular part is X^l Y^m Z^n type.
    TlMath Math;
    double pwr = (double)(l+m+n) + 3.0 / 2.0;
    double ans = 0.0;

    const int contractions = getContraction(s);
    for (int a = 0; a < contractions; a++) {
        const double coef_a = this->getCoefficient(s, a);
        const double norm_a = this->getNormalized(s, a, l, m, n);
        const double exp_a = this->getExponent(s, a);

        for (int b = 0; b < contractions; b++) {
            const double coef_b = this->getCoefficient(s, b);
            const double norm_b = this->getNormalized(s, b, l, m, n);
            const double exp_b = this->getExponent(s, b);

            double trm  = coef_a * coef_b;
            trm *= norm_a * norm_b;
            trm *= std::pow(exp_a + exp_b, -pwr);
            ans += trm;
        }
    }
    ans *= Math.dbfact(2*l-1) * Math.dbfact(2*m-1) * Math.dbfact(2*n-1);
    ans *= std::pow(2.0, -(double)(l+m+n));
    ans *= std::pow(Math.PI(), 3.0 / 2.0);
    ans =  Math.sqrt(1.0 / ans);

    return ans;
}


/// PGTOのnormalize
double Fl_Gto::getNormalized(const int s, const int p, const int l, const int m, const int n) const
{
    //   std::cout << TlUtils::format("Fl_Gto::getNormalized s=%d, p=%d, l=%d, m=%d, n=%d) ",
    //                   s, p, l, m, n)
    //        << TlUtils::format("cgto_size = %d, contraction = %d",
    //                   this->cgto.size(), static_cast<int>(this->cgto[s].contraction()))
    //        << std::endl;
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));
    assert((0 <= p) && (p < static_cast<int>(this->cgto[s].contraction())));
    assert(l >= 0);
    assert(m >= 0);
    assert(n >= 0);

    TlLogX& Log = TlLogX::getInstance();

    // chech consistent with shell.
    int maxangular = 0;
    switch (getShell(s)) {
    case 's':
        maxangular=0;
        break;
    case 'p':
        maxangular=1;
        break;
    case 'd':
        maxangular=2;
        break;
    case 'f':
        maxangular=3;
        break;
    default:
        Log << " *** Fl_Gto::getNormalized() not support over f angular momentum, sorry.\n";
        CnErr.abort();
        break;
    }

    if (l+m+n > maxangular) {
        Log << " *** Fl_Gto::getNormalized() illegal angular momentum, inconsistent with shell.\n";
        CnErr.abort();
    }

    // calculate normalized factor of pGTO, of which angular part is X^l Y^m Z^n type.
    TlMath Math;
    double pwr  = (double)(l+m+n);
    double ans  = std::pow(2.0, pwr);
    ans *= std::pow(Math.dbfact(2*l-1)*Math.dbfact(2*m-1)*Math.dbfact(2*n-1), -1.0/2.0);
    ans *= std::pow(2.0 / Math.PI(), 3.0 / 4.0);
    ans *= std::pow(this->getExponent(s,p), (pwr+3.0/2.0)/2.0);

    return ans;
}

//福江君に依頼された関数その１
int Fl_Gto::getBasiskindnumber() const
{
    std::string atom ="";
    std::string label ="";
    std::string residue ="";

    int count = 0;
    for (int p=0; p < this->getNumOfCGTOs(); p++) {
        if ((atom == cgto[p].atom) && (label == cgto[p].label)) {
            continue;
        }

        count++;
        atom  = cgto[p].atom;
        label = cgto[p].label;
    }

    return count;
}

//福江君に依頼された関数その２ (引数は0〜である)
int Fl_Gto::getTermnumber(int k) const
{
    std::string atom = "";
    std::string label = "";
    std::string residue = "";

    int position = -1;
    int p = 0;
    for (; p < this->getNumOfCGTOs(); ++p) {
        if ((atom == cgto[p].atom) && (label == cgto[p].label)) {
            continue;
        }
        ++position;

        if (position == k) {
            break;
        }
        atom = cgto[p].atom;
        label = cgto[p].label;
    }

    if (p == this->getNumOfCGTOs()) {
        return -1; // 引数が間違っている
    }

    // 引数で示されるk番目の原子種に属するCGTOはcgtoのp番目の要素番号から始まっている
    atom = cgto[p].atom;
    label = cgto[p].label;
    int count = 0;
    for (int r = p; r < this->getNumOfCGTOs(); ++r) {
        if ((atom == cgto[r].atom) && (label == cgto[r].label)) {
            ++count;
        }
    }

    return count;
}


int Fl_Gto::getStartposition(int k) const
{
    std::string atom = "";
    std::string label = "";

    int position = -1;
    int p = 0;
    for (; p < this->getNumOfCGTOs(); ++p) {
        if ((atom == cgto[p].atom) && (label == cgto[p].label)) {
            continue;
        }
        ++position;
        if (position == k) {
            break;
        }
        atom = cgto[p].atom;
        label = cgto[p].label;
    }

    if (p == this->getNumOfCGTOs()) {
        return -1; // 引数が間違っている
    }

    return p;
}

double Fl_Gto::getCoulombnormalized(const int s, const int p, const int l, const int m, const int n) const
{
    assert((0 <= s) && (s < static_cast<int>(this->cgto.size())));
    assert((0 <= p) && (p < static_cast<int>(this->cgto[s].contraction())));
    assert(l >= 0);
    assert(m >= 0);
    assert(n >= 0);
    //assert((l +m +n) < MaxAngular);

    // check consistent with shell.
    int maxangular = 0;
    switch (this->getShell(s)) {
    case 's':
        maxangular=0;
        break;
    case 'p':
        maxangular=1;
        break;
    case 'd':
        maxangular=2;
        break;
    case 'f':
        maxangular=3;
        break;
    default:
        CnErr.abort("Fl_Gto", "", "double getNormalized(int,int,int,int,int)",
                    "not support over f angular momentum, sorry.");
        break;
    }

    if (l+m+n > maxangular) {
        CnErr.abort("Fl_Gto", "", "double getNormalized(int,int,int,int,int)",
                    "illegal angular momentum, inconsistent with shell.");
    }

    // calculate normalized factor of pGTO, of which angular part is X^l Y^m Z^n type.
    TlMath Math;

    const int angular = l +m +n;
    double ans = 0.0;
    switch (angular) {
    case  0:
        /* s type pGTO */
        ans  = std::pow(2.0, -1.0/4.0);
        ans *= std::pow(Math.PI()/getExponent(s,p), -5.0/4.0);
        break;

    case 1:
        /* px, py, pz type pGTO */
        ans  = std::pow(2.0, -1.0/4.0);
        ans *= std::pow(Math.PI()/getExponent(s,p), -5.0/4.0);
        ans *= std::pow(12.0*getExponent(s,p), 1.0/2.0);
        break;

    case 2:
        if (l==2 || m==2 || n==2) {
            /* dxx, dyy, dzz type pGTO */
            ans  = std::pow(2.0, -1.0/4.0);
            ans *= std::pow(Math.PI()/getExponent(s,p), -5.0/4.0);
            ans *= std::pow(60.0*4.0*getExponent(s,p)*getExponent(s,p)/49.0, 1.0/2.0);
        } else {
            /* dxy, dyz, dzx type pGTO */
            ans  = std::pow(2.0, -1.0/4.0);
            ans *= std::pow(Math.PI()/getExponent(s,p), -5.0/4.0);
            ans *= std::pow(4.0*20.0*getExponent(s,p)*getExponent(s,p), 1.0/2.0);
        }
        break;
    default:
        CnErr.abort("Fl_Gto", "", "double getNormalized(int,int,int,int,int)",
                    "not support over d angular momentum, sorry.(H6/05/27)");
        break;
    }

    return ans;
}

