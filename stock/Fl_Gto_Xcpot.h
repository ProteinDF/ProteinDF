#ifndef FL_GTO_XCPOT_H
#define FL_GTO_XCPOT_H

#include <string>
#include  "Fl_Gto.h"

class Fl_Gto_Xcpot : public Fl_Gto {
public:
    explicit Fl_Gto_Xcpot(const std::string& path = Fl_Gto_Xcpot::getDefaultFileName()) : Fl_Gto(path) {
    }

    virtual ~Fl_Gto_Xcpot() {
    }

public:
    static std::string getDefaultFileName() {
        return "fl_Input/Fl_Gto_Xcpot";
    }
};

#endif // FL_GTO_XCPOT_H

