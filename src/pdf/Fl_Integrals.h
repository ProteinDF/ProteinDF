#ifndef FL_INTEGRALS_H
#define FL_INTEGRALS_H

#include <string>

class Fl_Integrals {
public:
    Fl_Integrals(const std::string& sFileName = "");
    virtual ~Fl_Integrals();

public:
    // read  fixed length
    bool read(int*, int*, int*, int*, int*, int*, int*, double*);

    // write fixed length
    bool write(int*, int*, int*, int*, int*, int*, int*, double*);

private:
    std::string ofile;
};

#endif // FL_INTEGRALS_H

