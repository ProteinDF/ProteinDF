#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cassert>

#include "Fl_Database.h"
#include "Fl_Db_Basis.h"
#include "TlUtils.h"
#include "TlLogX.h"

Fl_Db_Basis::Fl_Db_Basis(const std::string& sBasisName) : m_sBasisName(sBasisName)
{
    this->setData();
}


Fl_Db_Basis::~Fl_Db_Basis()
{
}


int Fl_Db_Basis::setData()
{
    //TlLogX& Log = TlLogX::getInstance();

    std::string DbFile = "basis2";
    const char* PdfHome = std::getenv("PDF_HOME");
    if (PdfHome != NULL) {
        DbFile = TlUtils::format("%s/data/basis2", PdfHome);
    }

    std::ifstream fi;
    fi.open(DbFile.c_str(), std::ios::in);
    if (!fi) {
        CnErr.abort(TlUtils::format("Cannot open %s", DbFile.c_str()));
    }

    //int flag = 0;
    if (this->m_sBasisName.compare(0, 2, "O-") == 0) {
        bool isFound = false;
        while (fi) {
            std::string line = "";
            std::getline(fi, line);

            if (this->m_sBasisName.compare(line) == 0) {
                int ScGto = 0;
                int PcGto = 0;
                int DcGto = 0;
                fi >> ScGto >> PcGto >> DcGto ;

                int Snum = 0;
                for (int i = 0; i < ScGto; ++i) {
                    Cgto tmp_cgto;
                    tmp_cgto.shell = 's';
                    tmp_cgto.scalefactor = 1.0;
                    tmp_cgto.normalizedfactor = 1.0;

                    int GtoNum = 0;
                    fi >> GtoNum;
                    tmp_cgto.contraction.resize(GtoNum);
                    for (int j = 0; j < GtoNum; ++j) {
                        fi >> tmp_cgto.contraction[j].exponent;
                        fi >> tmp_cgto.contraction[j].coefficient;
                    }
                    this->cgto.push_back(tmp_cgto);
                    Snum += GtoNum;
                }

                int Pnum = 0;
                for (int i = 0; i < PcGto; ++i) {
                    Cgto tmp_cgto;
                    tmp_cgto.shell = 'p';
                    tmp_cgto.scalefactor = 1.0;
                    tmp_cgto.normalizedfactor = 1.0;

                    int GtoNum = 0;
                    fi >> GtoNum;
                    tmp_cgto.contraction.resize(GtoNum);
                    for (int j = 0; j < GtoNum; ++j) {
                        fi >> tmp_cgto.contraction[j].exponent;
                        fi >> tmp_cgto.contraction[j].coefficient;
                    }
                    this->cgto.push_back(tmp_cgto);
                    Pnum += GtoNum;
                }

                int Dnum = 0;
                for (int i = 0; i < DcGto; i++) {
                    Cgto tmp_cgto;
                    tmp_cgto.shell = 'd';
                    tmp_cgto.scalefactor = 1.0;
                    tmp_cgto.normalizedfactor = 1.0;

                    int GtoNum = 0;
                    fi >> GtoNum;
                    tmp_cgto.contraction.resize(GtoNum);
                    for (int j = 0; j < GtoNum; j++) {
                        fi >> tmp_cgto.contraction[j].exponent;
                        fi >> tmp_cgto.contraction[j].coefficient;
                    }
                    this->cgto.push_back(tmp_cgto);
                    Dnum+=GtoNum;
                }

                isFound = true;
                break;
            }
        }

        if (isFound == false) {
            CnErr.abort(TlUtils::format("The orbital basis set \'%s\' is not found. Please check Basis set name.",
                                        this->m_sBasisName.c_str()));
        }

    } else if (this->m_sBasisName.compare(0, 2, "A-") == 0) {
        bool isFound = false;
        while (fi) {
            std::string line = "";
            std::getline(fi, line);

            if (this->m_sBasisName.compare(line) == 0) {
                int ScGto = 0;
                int PcGto = 0;
                int DcGto = 0;
                int SPDcGto = 0;
                fi >> ScGto >> PcGto >> DcGto >>SPDcGto;

                if (SPDcGto == 0) {
                    int Snum = 0;
                    if (ScGto != 0) {
                        fi >>  Snum;
                        for (int i = 0; i < ScGto; i++) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell = 's';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            fi >> tmp_cgto.contraction[0].exponent;
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->rhoCgto.push_back(tmp_cgto);
                        }
                    }

                    int Pnum = 0;
                    if (PcGto != 0) {
                        fi >> Pnum;
                        for (int i = 0; i < PcGto; i++) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell = 'p';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            fi >> tmp_cgto.contraction[0].exponent;
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->rhoCgto.push_back(tmp_cgto);
                        }
                    }

                    int Dnum = 0;
                    if (DcGto != 0) {
                        fi >> Dnum;

                        for (int i = 0; i < DcGto; i++) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='d';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            fi >> tmp_cgto.contraction[0].exponent;
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->rhoCgto.push_back(tmp_cgto);
                        }
                    }

                } else {
                    if ((PcGto == 0) && (DcGto == 0)) {
                        int Snum = 0;
                        if (ScGto != 0) {
                            fi >>  Snum;
                            for (int i = 0; i < ScGto; i++) {
                                Cgto tmp_cgto;
                                tmp_cgto.shell = 's';
                                tmp_cgto.scalefactor = 1.0;
                                tmp_cgto.normalizedfactor = 1.0;

                                tmp_cgto.contraction.resize(1);
                                fi >> tmp_cgto.contraction[0].exponent;
                                tmp_cgto.contraction[0].coefficient = 1.0;

                                this->rhoCgto.push_back(tmp_cgto);
                            }
                        }

                        int SPDnum = 0;
                        fi >> SPDnum;
                        std::vector<double> spd_exponent(SPDnum, 0.0);
                        for (int i = 0; i < SPDnum; ++i) {
                            fi >> spd_exponent[i];
                        }

                        for (int i = 0; i < SPDnum; ++i) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='s';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            tmp_cgto.contraction[0].exponent = spd_exponent[i];
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->rhoCgto.push_back(tmp_cgto);
                        }
                        for (int i = 0; i < SPDnum; ++i) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='p';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            tmp_cgto.contraction[0].exponent = spd_exponent[i];
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->rhoCgto.push_back(tmp_cgto);
                        }
                        for (int i = 0; i < SPDnum; i++) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='d';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            tmp_cgto.contraction[0].exponent = spd_exponent[i];
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->rhoCgto.push_back(tmp_cgto);
                        }
                    }    /* if(PcGto==0 && DcGto==0) */
                }     /* else */

                //----------------------Myu Myu Myu
                fi >> ScGto >> PcGto >> DcGto >>SPDcGto;
                if (SPDcGto == 0) {
                    int Snum = 0;
                    if (ScGto != 0) {
                        fi >>  Snum;
                        for (int i = 0; i < ScGto; ++i) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='s';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            fi >> tmp_cgto.contraction[0].exponent;
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->myuCgto.push_back(tmp_cgto);
                        }
                    }

                    int Pnum = 0;
                    if (PcGto != 0) {
                        fi >> Pnum;
                        for (int i = 0; i < PcGto; i++) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='p';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            fi >> tmp_cgto.contraction[0].exponent;
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->myuCgto.push_back(tmp_cgto);
                        }
                    }

                    int Dnum = 0;
                    if (DcGto != 0) {
                        fi >> Dnum;
                        for (int i = 0; i < DcGto; i++) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='d';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            fi >> tmp_cgto.contraction[0].exponent;
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->myuCgto.push_back(tmp_cgto);
                        }
                    }

                } else {
                    if ((PcGto == 0) && (DcGto == 0)) {
                        int Snum=0;
                        if (ScGto != 0) {
                            fi >>  Snum;
                            for (int i = 0; i < ScGto; ++i) {
                                Cgto tmp_cgto;
                                tmp_cgto.shell='s';
                                tmp_cgto.scalefactor = 1.0;
                                tmp_cgto.normalizedfactor = 1.0;

                                tmp_cgto.contraction.resize(1);
                                fi >> tmp_cgto.contraction[0].exponent;
                                tmp_cgto.contraction[0].coefficient = 1.0;

                                this->myuCgto.push_back(tmp_cgto);
                            }
                        }

                        int SPDnum = 0;
                        fi >> SPDnum;
                        std::vector<double> spd_exponent(SPDnum, 0.0);
                        for (int i = 0; i < SPDnum; ++i) {
                            fi >> spd_exponent[i];
                        }

                        for (int i = 0; i < SPDnum; ++i) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='s';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            tmp_cgto.contraction[0].exponent = spd_exponent[i];
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->myuCgto.push_back(tmp_cgto);
                        }
                        for (int i = 0; i < SPDnum; ++i) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='p';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            tmp_cgto.contraction[0].exponent = spd_exponent[i];
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->myuCgto.push_back(tmp_cgto);
                        }
                        for (int i = 0; i < SPDnum; ++i) {
                            Cgto tmp_cgto;
                            tmp_cgto.shell='d';
                            tmp_cgto.scalefactor = 1.0;
                            tmp_cgto.normalizedfactor = 1.0;

                            tmp_cgto.contraction.resize(1);
                            tmp_cgto.contraction[0].exponent = spd_exponent[i];
                            tmp_cgto.contraction[0].coefficient = 1.0;

                            this->myuCgto.push_back(tmp_cgto);
                        }


                    }    /* if(PcGto==0 && DcGto==0) */

                }     /* else */

                isFound = true;
                break;
            }
        }

        if (isFound == false) {
            CnErr.abort(TlUtils::format("The aux orbital basis set \'%s\' is not found. Please check Basis set name.",
                                        this->m_sBasisName.c_str()));
        }

    } else if (this->m_sBasisName.compare(0, 2, "P-") == 0) {
        //今は空。
    }

    fi.close();

    return 0;
}

