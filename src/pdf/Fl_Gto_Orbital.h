#ifndef FL_GTO_ORBITAL_H
#define FL_GTO_ORBITAL_H

#include <string>
#include "Fl_Gto.h"

class Fl_Gto_Orbital : public Fl_Gto {
public:
    explicit Fl_Gto_Orbital(const std::string& path = Fl_Gto_Orbital::getDefaultFileName()) : Fl_Gto(path) {
    }

    virtual ~Fl_Gto_Orbital() {
    }

public:
    static std::string getDefaultFileName() {
        return "fl_Input/fl_Gto_Orbital";
    }
};

#endif // FL_GTO_ORBITAL_H
