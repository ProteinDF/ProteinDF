#include <iostream>
#include <fstream>

#include "Fl_Integrals.h"
#include "Common.h"

Fl_Integrals::Fl_Integrals(const std::string& sFileName) : ofile(sFileName)
{
}

Fl_Integrals::~Fl_Integrals()
{
}

// read fixed length
bool Fl_Integrals::read(int* icont, int* nABC,   int* nC,     int* indexC,
                        int* nAB,   int* indexA, int* indexB, double* ABC)
{
    std::ifstream ifs(this->ofile.c_str(), std::ios::in);

    ifs.read((char*) icont, (int)sizeof(int) * 1);
    ifs.read((char*) nABC, (int)sizeof(int) * 1);
    ifs.read((char*) nC, (int)sizeof(int) * 1);
    ifs.read((char*) indexC, (int)sizeof(int) * MAXNA);
    ifs.read((char*) nAB, (int)sizeof(int) * MAXNA);
    ifs.read((char*) indexA, (int)sizeof(int) * MAXNPQA);
    ifs.read((char*) indexB, (int)sizeof(int) * MAXNPQA);
    ifs.read((char*) ABC, (int)sizeof(double) * MAXNPQA);

    ifs.close();

    return true;
}

// write function
bool Fl_Integrals::write(int* icont, int* nABC,   int* nC,     int* indexC,
                         int* nAB,   int* indexA, int* indexB, double* ABC)
{
    std::ofstream ofs(this->ofile.c_str(), std::ios::out | std::ios::trunc);

    ofs.write((char*) icont, (int)sizeof(int) * 1);
    ofs.write((char*) nABC, (int)sizeof(int) * 1);
    ofs.write((char*) nC, (int)sizeof(int) * 1);
    ofs.write((char*) indexC, (int)sizeof(int) * MAXNA);
    ofs.write((char*) nAB, (int)sizeof(int) * MAXNA);
    ofs.write((char*) indexA, (int)sizeof(int) * MAXNPQA);
    ofs.write((char*) indexB, (int)sizeof(int) * MAXNPQA);
    ofs.write((char*) ABC, (int)sizeof(double) * MAXNPQA);

    return true;
}


