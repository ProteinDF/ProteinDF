#ifndef FL_OUT_H
#define FL_OUT_H

#include "FileX.h"

class Fl_Out : public FileX {
public:
    Fl_Out(const std::string& rFilePath) : FileX(rFilePath, std::ios_base::out | std::ios_base::app) {
        this->open();
    };

    virtual ~Fl_Out() {
        this->flush();
        this->close();
    };
};


class Fl_Out_Sys : public Fl_Out {
public:
    Fl_Out_Sys() : Fl_Out("fl_Out_Sys") {
        // TH NOTE for debug
        //std::cerr << "[TH] object Fl_Out_Sys is created" << std::endl;
    }
};


class Fl_Out_Std : public Fl_Out {
public:
    Fl_Out_Std() : Fl_Out("fl_Out_Std") {
        // TH NOTE for debug
        //std::cerr << "[TH] object Fl_Out_Std is created." << std::endl;
    }
    virtual ~Fl_Out_Std() {
        //std::cerr << "[TH] object Fl_Out_Std is destructed." << std::endl;
    }
};


class Fl_Out_Arc : public Fl_Out {
public:
    Fl_Out_Arc() : Fl_Out("fl_Out_Arc") {
        // TH NOTE for debug
        //std::cerr << "[TH] object Fl_Out_Arc is created." << std::endl;
    }
};

#endif // FL_OUT_H
