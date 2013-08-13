#include <sstream>
#include "DfInputdata.h"

#include "CnError.h"
#include "TlPrdctbl.h"
#include "TlResidue.h"
#include "Fl_Geometry.h"
#include "Fl_Gto.h"

#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"
#include "TlOrbitalInfo_XC.h"

#include "TlUtils.h"
#include "TlLogging.h"
#include "TlFile.h"

#include "PdfKeyword.h"
#include "PdfUserInput.h"
#include "TlOrbitalInfo.h"
#include "TlOrbitalInfo_Density.h"
#include "TlOrbitalInfo_XC.h"
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
    pdfKwd.convertAlias(&(inputParam));
    
    // default値にユーザー入力値を上書きする
    param.merge(inputParam);

    // keyword check
    pdfKwd.convertAlias(&(param));
    pdfKwd.checkInputParam(param);

    // 
    {
        const Fl_Geometry flGeom(param["coordinates"]);
        const TlOrbitalInfo orbInfo(param["coordinates"], param["basis_sets"]);
        const TlOrbitalInfo_Density orbInfo_J(param["coordinates"], param["basis_sets_j"]);
        const TlOrbitalInfo_XC orbInfo_XC(param["coordinates"], param["basis_sets_k"]);

        param["num_of_atoms"] = flGeom.getNumOfAtoms();
        param["num_of_dummy_atoms"] = flGeom.getNumOfDummyAtoms();
        param["num_of_AOs"] = orbInfo.getNumOfOrbitals();
        param["num_of_auxCDs"] = orbInfo_J.getNumOfOrbitals();
        param["num_of_auxXCs"] = orbInfo_XC.getNumOfOrbitals();
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
    TlLogging& log = TlLogging::getInstance();

    PdfKeyword kwd;

    // print comment
    log.info(TlUtils::format(" Comment: %s", data["comment"].getStr().c_str()));
    
    // print input keyword list
    if (data["show_keyword"].getBoolean() == true) {
        log.info(" >>>> Available Keywords <<<<");

        std::stringstream ss;
        kwd.printDefault(ss);
        log.info(ss.str());
    } else {
        log.info(" printing available keywords is rejected.");
        log.info(" to print the keywords, please input 'show_keyword = yes'.");
    }
    log.info("\n");

    // print global input
    if (data["show_input"].getBoolean() == true) {
        log.info(" >>>> Input Parameters <<<<");

        std::stringstream ss;
        kwd.print(ss, data);
        log.info(ss.str());
    } else {
        log.info(" printing available keywords is rejected.");
        log.info(" to print the keywords, please input 'show_input = yes'.");
    }
    log.info("\n");

    // print coordinates
    if (data["show_coordinates"].getBoolean() == true) {
        log.info(" >>>> Coordinates <<<<");
        log.info(" symbol charge (coord) [label]");

        Fl_Geometry geom(data["coordinates"]);
        const int numOfAtoms = geom.getNumOfAtoms();
        for (int i = 0; i < numOfAtoms; ++i) {
            const std::string symbol = geom.getAtom(i);
            const double charge = geom.getCharge(i);
            const std::string label = geom.getLabel(i);
            const TlPosition pos = geom.getCoordinate(i);
            
            log.info(TlUtils::format(" %-2s %+8.3e (%+16.10e %+16.10e %+16.10e) [%s] \n",
                                     symbol.c_str(), charge,
                                     pos.x(), pos.y(), pos.z(), label.c_str()));
        }
    } else {
        log.info(" printing coordinates is rejected.");
        log.info(" to print the keywords, please input 'show_coordinates = yes'.");
    }
    log.info("\n");
    
    // print basis set
    {
        const std::string showOrbitalBasis = TlUtils::toUpper(data["show_orbital_basis"].getStr());

        const Fl_Gto orb(data["basis_sets"]);
        const Fl_Gto orb_j(data["basis_sets_j"]);
        const Fl_Gto orb_k(data["basis_sets_k"]);
        const Fl_Gto orb_GF(data["basis_sets_GF"]);
        const bool isUseGF = data["grid_free"].getBoolean();
        if (showOrbitalBasis == "GAMESS") {
            log.info(" >>>> Inputted Orbital Basis Set (GAMESS format) <<<<");
            log.info(orb.getStr_GAMESS());
            log.info("\n");
            log.info(" >>>> Inputted Density Orbital Basis Set (GAMESS format) <<<<");
            log.info(orb_j.getStr_GAMESS());
            log.info("\n");
            log.info(" >>>> Inputted XC Orbital Basis Set (GAMESS format) <<<<");
            log.info(orb_k.getStr_GAMESS());
            log.info("\n");
            if (isUseGF) {
                log.info(" >>>> Inputted GridFree Orbital Basis Set (GAMESS format) <<<<");
                log.info(orb_GF.getStr_GAMESS());
                log.info("\n");
            }
        } else if (showOrbitalBasis == "AMOSS") {
            log.info(" >>>> Inputted Orbital Basis Set (AMOSS format) <<<<");
            log.info(orb.getStr_AMOSS());
            log.info("\n");
            log.info(" >>>> Inputted Density Orbital Basis Set (GAMESS format) <<<<");
            log.info(orb_j.getStr_AMOSS());
            log.info("\n");
            log.info(" >>>> Inputted XC Orbital Basis Set (GAMESS format) <<<<");
            log.info(orb_k.getStr_AMOSS());
            log.info("\n");
            if (isUseGF) {
                log.info(" >>>> Inputted XC Orbital Basis Set (GAMESS format) <<<<");
                log.info(orb_GF.getStr_AMOSS());
                log.info("\n");
            }
        }
    }
}

