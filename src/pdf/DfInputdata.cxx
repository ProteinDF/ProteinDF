#include "DfInputdata.h"

#include "CnError.h"
#include "TlPrdctbl.h"
#include "TlResidue.h"
#include "Fl_Geometry.h"
#include "Fl_Gto.h"
#include "Fl_Gto_Orbital.h"
#include "Fl_Gto_Density.h"
#include "Fl_Gto_Xcpot.h"
#include "Fl_Gto_Xcpot2.h"
#include "Fl_Db_Basis.h"
#include "Fl_Tbl_Orbital.h"
#include "Fl_Tbl_Density.h"
#include "Fl_Tbl_Xcpot.h"
#include "TlUtils.h"
#include "TlLogX.h"
#include "TlFile.h"

#include "PdfKeyword.h"
#include "PdfUserInput.h"
#include "TlMsgPack.h"

DfInputdata::DfInputdata()
{
}


DfInputdata::~DfInputdata()
{
}


TlSerializeData DfInputdata::main()
{
    // load default parameters
    PdfKeyword pdfKwd;
    TlSerializeData param = pdfKwd.getDefault();
    
    // include user input parameters
    const std::string mpacFilePath = "pdfparam.mpac";
    if (TlFile::isExist(mpacFilePath) == true) {
        TlMsgPack msgPack;
        msgPack.load(mpacFilePath);
        TlSerializeData tmpParam = msgPack.getSerializeData();
        param.merge(tmpParam);
    }

    // fl_Userinputを読み取る
    PdfUserInput pdfUserInput;
    pdfUserInput.load();
    TlSerializeData inputParam = pdfUserInput.getSerializeData();
    pdfKwd.convertAlias(&(inputParam["model"]));
    
    // default値にユーザー入力値を上書きする
    param.merge(inputParam);

    // keyword check
    pdfKwd.convertAlias(&(param["model"]));
    pdfKwd.checkInputParam(param["model"]);
    
    // テーブルの作成
    {
        Fl_Tbl_Orbital Tbl;
        param["model"]["AOs"] = Tbl.getcGtoTotalNum();
    }

    {
        Fl_Tbl_Density Tbl;
        param["model"]["AuxCDs"] = Tbl.getcGtoTotalNum();
    }

    {
        Fl_Tbl_Xcpot Tbl;
        param["model"]["AuxXCs"] = Tbl.getcGtoTotalNum();
    }

    {
        Fl_Geometry  Geom(Fl_Geometry::getDefaultFileName());
        param["model"]["atoms"] = Geom.getNumOfAtoms();
        param["model"]["dummyAtoms"] = Geom.getDummyatom();
    }

    // 保存
    this->data_ = param;
    TlMsgPack msgPack(this->data_);
    msgPack.save(mpacFilePath);

    // 表示
    this->show(this->data_);

    return this->data_;
}


void DfInputdata::show(const TlSerializeData& data) const
{
    TlLogX& log = TlLogX::getInstance();

    PdfKeyword kwd;

    // print comment
    log << "  >>>> Comment <<<<\n";
    log << data["model"]["comment"].getStr() << "\n";
    log << "\n";
    
    // print input keyword list
    if (data["model"]["show_keyword"].getBoolean() == true) {
        log << " >>>> Available Keywords <<<<\n";
        kwd.printDefault(log);
    } else {
        log << " printing available keywords is rejected.\n";
        log << " to print the keywords, please input 'show_keyword = yes'.\n";
    }
    log << "\n";

    // print global input
    if (data["model"]["show_input"].getBoolean() == true) {
        log << " >>>> Input Parameters <<<<\n";
        kwd.print(log, data["model"]);
    } else {
        log << " printing available keywords is rejected.\n";
        log << " to print the keywords, please input 'show_input = yes'.\n";
    }
    log << "\n";

    // print coordinates
    if (data["model"]["show_coordinates"].getBoolean() == true) {
        log << " >>>> Coordinates <<<<\n";
        log << " symbol charge (coord) [label]\n";

        Fl_Geometry geom(data["model"]["coordinates"]);
        const int numOfAtoms = geom.getNumOfAtoms();
        for (int i = 0; i < numOfAtoms; ++i) {
            const std::string symbol = geom.getAtom(i);
            const double charge = geom.getCharge(i);
            const std::string label = geom.getLabel(i);
            const TlPosition pos = geom.getCoordinate(i);
            
            log << TlUtils::format(" %-2s %+8.3e (%+16.10e %+16.10e %+16.10e) [%s] \n",
                                   symbol.c_str(), charge,
                                   pos.x(), pos.y(), pos.z(), label.c_str());
        }
    } else {
        log << " printing coordinates is rejected.\n";
        log << " to print the keywords, please input 'show_coordinates = yes'.\n";
    }
    log << "\n";
    
    // print basis set
    {
        const std::string showOrbitalBasis = TlUtils::toUpper(data["model"]["show_orbital_basis"].getStr());

        const Fl_Gto_Orbital flGtoOrbital;
        if (showOrbitalBasis == "GAMESS") {
            log << " >>>> Inputted Orbital Basis Set (GAMESS format) <<<<\n";
            flGtoOrbital.showGAMESS();
            log << "\n";
        } else if (showOrbitalBasis == "AMOSS") {
            log << " >>>> Inputted Orbital Basis Set (AMOSS format) <<<<\n";
            flGtoOrbital.showAMOSS();
            log << "\n";
        }

        const Fl_Gto_Density flGtoDensity;
        if (showOrbitalBasis == "GAMESS") {
            log << " >>>> Inputted Density Orbital Basis Set (GAMESS format) <<<<\n";
            flGtoDensity.showGAMESS();
            log << "\n";
        } else if (showOrbitalBasis == "AMOSS") {
            log << " >>>> Inputted Density Orbital Basis Set (AMOSS format) <<<<\n";
            flGtoDensity.showAMOSS();
            log << "\n";
        }
    }
}

