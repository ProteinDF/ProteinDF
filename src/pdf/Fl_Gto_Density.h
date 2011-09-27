#ifndef CLASS_FL_GTO_DENSITY
#define CLASS_FL_GTO_DENSITY

#include <string>
#include "Fl_Gto.h"

class Fl_Gto_Density : public Fl_Gto {
public:
    explicit Fl_Gto_Density(const std::string& path = Fl_Gto_Density::getDefaultFileName()) : Fl_Gto(path) {
    }

    virtual ~Fl_Gto_Density() {
    }

public:
    static std::string getDefaultFileName() {
        return "fl_Input/fl_Gto_Density";
    }
};

#endif // CLASS_FL_GTO_DENSITY

