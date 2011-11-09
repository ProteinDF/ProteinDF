#include <iostream>

#include "DfIGuess2.h"
#include "TlVector.h"
#include "TlParameter.h"
#include "CnError.h"

void DfIGuess2::main()
{
    //---- file open and read data ---------------------------------------------------
    std::ifstream  fi;
    fi.open("guess2.txt", std::ios::in);
    if (!fi) {
        this->log_.error("Cannot open \"guess2.txt\"");
        CnErr.abort();
    }
    this->log_.info("./guess2.txt is opened");

    int term;
    fi >> term;

    TlVector coefR(term);
    TlVector coefM(term);
    TlVector coefN(term);
    TlVector coefA(term);
    TlVector coefRA(term);
    TlVector coefRB(term);

    for (int k=0; k<term; k++) {
        double tmp;
        fi >> tmp;
        coefR[k] = tmp;
        fi >> tmp;
        coefM[k] = tmp;
        fi >> tmp;
        coefN[k] = tmp;
    }
    fi.close();

    TlVector nalpha;
    nalpha.load("fl_Work/fl_Vct_Nalpha");

    this->log_.info("./guess2.txt is closed");

    //---- Normalized by the number of electrons -------------------------------------
    int elenum =0;
    int aspnum =0, bspnum =0;

    TlParameter flGbi;
    flGbi.load("fl_Input/fl_Globalinput");
    //Fl_Globalinput   Gbi(">>>>SCF");
    std::string method = flGbi["SCF"]["method"];

    if ("roks" == method) {
        elenum = atoi(flGbi["SCF"]["method/roks/electron-number-beta"].c_str());
    } else if ("sp" == method) {
        aspnum = atoi(flGbi["SCF"]["UKS/alpha_electrons"].c_str());
        bspnum = atoi(flGbi["SCF"]["UKS/beta_electrons"].c_str());
        elenum = aspnum + bspnum;
    } else {
        elenum = atoi(flGbi["SCF"]["RKS/num_electrons"].c_str());
    }

    double dumelenum = 0.0;
    for (int k=0; k<term; k++) {
        dumelenum += coefR[k]*coefA[k];
    }

    std::cout << "    number of electron = " << elenum << "\n";
    std::cout << "    Rho * Nalpha       = " << dumelenum << "\n";
    std::cout << "    difference         = " << elenum - dumelenum << "\n";
    std::cout << "    Rho is normalized by the number of electrons." << "\n";
    Log  << "    number of electron = " << elenum << "\n";
    Log  << "    Rho * Nalpha       = " << dumelenum << "\n";
    Log  << "    difference         = " << elenum - dumelenum << "\n";
    Log  << "    Rho is normalized by the number of electrons." << "\n";

    for (int k=0; k<term; k++) {
        coefR[k] *= elenum / dumelenum;
    }

    Log << "    SCF method is " << std::string(method) << "\n\n";

    if ("sp" == method) {
        for (int k=0; k<term; k++) {
            coefRA[k] *= aspnum/elenum;
            coefRB[k] *= bspnum/elenum;
        }
    }

    //---- output fl_Vct_Rou, Myu, Nyu -----------------------------------------------
    if ("sp" == method) {
        Log << " output fl_Vct_Roua0, Roub0, Myua0, Myub0, Nyua0, Nyub0\n\n";
        {
//       Fl_Vct_Rou fl_Vct_Rou_Alpha(0,"a");
//       fl_Vct_Rou_Alpha.open("fl_Work", fl_Vct_Rou_Alpha.getfilename(), "write");
//       fl_Vct_Rou_Alpha.putelemnum(&term);
//       fl_Vct_Rou_Alpha.write(term,coefRA);
//       fl_Vct_Rou_Alpha.close();
            coefRA.save("fl_Work/fl_Vct_Roua0");

//       Fl_Vct_Myu fl_Vct_Myu_Alpha(0,"a");
//       fl_Vct_Myu_Alpha.open("fl_Work", fl_Vct_Myu_Alpha.getfilename(), "write");
//       fl_Vct_Myu_Alpha.putelemnum(&term);
//       fl_Vct_Myu_Alpha.write(term,coefM);
//       fl_Vct_Myu_Alpha.close();
            coefM.save("fl_Work/fl_Vct_Myua0");

//       Fl_Vct_Nyu fl_Vct_Nyu_Alpha(0,"a");
//       fl_Vct_Nyu_Alpha.open("fl_Work", fl_Vct_Nyu_Alpha.getfilename(), "write");
//       fl_Vct_Nyu_Alpha.putelemnum(&term);
//       fl_Vct_Nyu_Alpha.write(term,coefN);
//       fl_Vct_Nyu_Alpha.close();
            coefN.save("fl_Work/fl_Vct_Nyua0");

//       Fl_Vct_Rou fl_Vct_Rou_Beta(0,"b");
//       fl_Vct_Rou_Beta.open("fl_Work", fl_Vct_Rou_Beta.getfilename(), "write");
//       fl_Vct_Rou_Beta.putelemnum(&term);
//       fl_Vct_Rou_Beta.write(term,coefRB);
//       fl_Vct_Rou_Beta.close();
            coefRB.save("fl_Work/fl_Vct_Roub0");

//       Fl_Vct_Myu fl_Vct_Myu_Beta(0,"b");
//       fl_Vct_Myu_Beta.open("fl_Work", fl_Vct_Myu_Beta.getfilename(), "write");
//       fl_Vct_Myu_Beta.putelemnum(&term);
//       fl_Vct_Myu_Beta.write(term,coefM);
//       fl_Vct_Myu_Beta.close();
            coefM.save("fl_Work/fl_Vct_Myub0");

//       Fl_Vct_Nyu fl_Vct_Nyu_Beta(0,"b");
//       fl_Vct_Nyu_Beta.open("fl_Work", fl_Vct_Nyu_Beta.getfilename(), "write");
//       fl_Vct_Nyu_Beta.putelemnum(&term);
//       fl_Vct_Nyu_Beta.write(term,coefN);
//       fl_Vct_Nyu_Beta.close();
            coefN.save("fl_Work/fl_Vct_Nyub0");
        }

        Log << " output fl_Vct_Roua1, Roub1, Myua1, Myub1, Nyua1, Nyub1\n\n";
        {
//       Fl_Vct_Rou       fl_Vct_Rou_Alpha(1,"a");
//       fl_Vct_Rou_Alpha.open("fl_Work", fl_Vct_Rou_Alpha.getfilename(), "write");
//       fl_Vct_Rou_Alpha.putelemnum(&term);
//       fl_Vct_Rou_Alpha.write(term,coefRA);
//       fl_Vct_Rou_Alpha.close();
            coefRA.save("fl_Work/fl_Vct_Roua1");

//       Fl_Vct_Myu       fl_Vct_Myu_Alpha(1,"a");
//       fl_Vct_Myu_Alpha.open("fl_Work", fl_Vct_Myu_Alpha.getfilename(), "write");
//       fl_Vct_Myu_Alpha.putelemnum(&term);
//       fl_Vct_Myu_Alpha.write(term,coefM);
//       fl_Vct_Myu_Alpha.close();
            coefM.save("fl_Work/fl_Vct_Myua1");

//       Fl_Vct_Nyu       fl_Vct_Nyu_Alpha(1,"a");
//       fl_Vct_Nyu_Alpha.open("fl_Work", fl_Vct_Nyu_Alpha.getfilename(), "write");
//       fl_Vct_Nyu_Alpha.putelemnum(&term);
//       fl_Vct_Nyu_Alpha.write(term,coefN);
//       fl_Vct_Nyu_Alpha.close();
            coefN.save("fl_Work/fl_Vct_Nyua1");

//       Fl_Vct_Rou      fl_Vct_Rou_Beta(1,"b");
//       fl_Vct_Rou_Beta.open("fl_Work", fl_Vct_Rou_Beta.getfilename(), "write");
//       fl_Vct_Rou_Beta.putelemnum(&term);
//       fl_Vct_Rou_Beta.write(term,coefRB);
//       fl_Vct_Rou_Beta.close();
            coefRB.save("fl_Work/fl_Vct_Roub1");

//       Fl_Vct_Myu      fl_Vct_Myu_Beta(1,"b");
//       fl_Vct_Myu_Beta.open("fl_Work", fl_Vct_Myu_Beta.getfilename(), "write");
//       fl_Vct_Myu_Beta.putelemnum(&term);
//       fl_Vct_Myu_Beta.write(term,coefM);
//       fl_Vct_Myu_Beta.close();
            coefM.save("fl_Work/fl_Vct_Myub1");

//       Fl_Vct_Nyu      fl_Vct_Nyu_Beta(1,"b");
//       fl_Vct_Nyu_Beta.open("fl_Work", fl_Vct_Nyu_Beta.getfilename(), "write");
//       fl_Vct_Nyu_Beta.putelemnum(&term);
//       fl_Vct_Nyu_Beta.write(term,coefN);
//       fl_Vct_Nyu_Beta.close();
            coefN.save("fl_Work/fl_Vct_Nyub1");
        }
    } else if ("nsp" == method) {
        Log << " output fl_Vct_Rou0, Myu0, Nyu0\n\n";
        {
//       Fl_Vct_Rou  fl_Vct_Rou(0);
//       fl_Vct_Rou.open("fl_Work", fl_Vct_Rou.getfilename(), "write");
//       fl_Vct_Rou.putelemnum(&term);
//       fl_Vct_Rou.write(term,coefR);
//       fl_Vct_Rou.close();
            coefR.save("fl_Work/fl_Vct_Rou0");

//       Fl_Vct_Myu  fl_Vct_Myu(0);
//       fl_Vct_Myu.open("fl_Work", fl_Vct_Myu.getfilename(), "write");
//       fl_Vct_Myu.putelemnum(&term);
//       fl_Vct_Myu.write(term,coefM);
//       fl_Vct_Myu.close();
            coefM.save("fl_Work/fl_Vct_Myu0");

//       Fl_Vct_Nyu  fl_Vct_Nyu(0);
//       fl_Vct_Nyu.open("fl_Work", fl_Vct_Nyu.getfilename(), "write");
//       fl_Vct_Nyu.putelemnum(&term);
//       fl_Vct_Nyu.write(term,coefN);
//       fl_Vct_Nyu.close();
            coefN.save("fl_work/fl_Vct_Nyu0");
        }
        Log << " output fl_Vct_Rou1, Myu1, Nyu1\n\n";
        {
//       Fl_Vct_Rou  fl_Vct_Rou(1);
//       fl_Vct_Rou.open("fl_Work", fl_Vct_Rou.getfilename(), "write");
//       fl_Vct_Rou.putelemnum(&term);
//       fl_Vct_Rou.write(term,coefR);
//       fl_Vct_Rou.close();
            coefR.save("fl_Work/fl_Vct_Rou1");

//       Fl_Vct_Myu  fl_Vct_Myu(1);
//       fl_Vct_Myu.open("fl_Work", fl_Vct_Myu.getfilename(), "write");
//       fl_Vct_Myu.putelemnum(&term);
//       fl_Vct_Myu.write(term,coefM);
//       fl_Vct_Myu.close();
            coefM.save("fl_Work/fl_Vct_Myu1");

//       Fl_Vct_Nyu  fl_Vct_Nyu(1);
//       fl_Vct_Nyu.open("fl_Work", fl_Vct_Nyu.getfilename(), "write");
//       fl_Vct_Nyu.putelemnum(&term);
//       fl_Vct_Nyu.write(term,coefN);
//       fl_Vct_Nyu.close();
            coefN.save("fl_Work/fl_Vct_Nyu1");
        }
    } else {
        std::cout << "illegal method is specified." << std::endl;
        CnErr.abort();
    }

//   delete [] coefRB;
//   delete [] coefRA;
//   delete [] coefR;
//   delete [] coefM;
//   delete [] coefN;
//   delete [] coefA;

}
