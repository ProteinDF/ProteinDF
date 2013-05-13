#ifndef FL_GTO_XCPOT2_H
#define FL_GTO_XCPOT2_H

#include  "Fl_Gto.h"

class Fl_Gto_Xcpot2 : public Fl_Gto {
public:
    Fl_Gto_Xcpot2(const std::string& path = Fl_Gto_Xcpot2::getDefaultFileName()) : Fl_Gto(path) {
    }

    virtual ~Fl_Gto_Xcpot2() {
    }

public:
    static std::string getDefaultFileName() {
        return "fl_Input/Fl_Gto_Xcpot2";
    }
    
};

#endif // FL_GTO_XCPOT2_H

