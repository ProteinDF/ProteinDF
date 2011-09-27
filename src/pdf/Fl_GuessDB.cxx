#include <cassert>
#include "Fl_GuessDB.h"
#include "TlLogX.h"

Fl_GuessDB::Fl_GuessDB()
        : DbFile(""), HdFile(""), DbLabel("")
{

    DBM.TotalRouterm = DBM.TotalMyuterm = DBM.TotalNyuterm = 0;
    Atm_or_Mol = 0;
    DBM.TotalAtomNum = 0;
}


Fl_GuessDB::~Fl_GuessDB()
{
}


int Fl_GuessDB::serchLine(const std::string& filename, const std::string& dblabel) const
{
    return this->getHeaderAtomNum(filename, dblabel, 1);
}


int Fl_GuessDB::getHeaderAtomNum(const std::string& filename, const std::string& dblabel, int hantei) const
{
    TlLogX& Log = TlLogX::getInstance();

    std::ifstream fi;
    fi.open(filename.c_str(), std::ios::in);
    if (!fi) {
        Log <<"Cannot Open "<< filename << "\n";
        Log <<"\7"<<"\n";
        CnErr.abort();
    }

    if (dblabel.empty()) {
        Log<<"Bad appointment DBLABEL\n";
        Log<<"in Fl_GuessDB::getHeaderAtomNum\n";
        CnErr.abort();
    }

    int  AtomNum,aboutLineNumber;
    for (;;) {
        //fi >> DLB;
        std::string DLB;
        std::getline(fi, DLB);
        if (DLB == dblabel) {
            fi >> AtomNum >> aboutLineNumber ;
            break;
        }
    }
    fi.close();

    if (hantei == 0) {
        return  AtomNum;
    } else if (hantei == 1) {
        return  aboutLineNumber;
    }

    return 0;
}


int Fl_GuessDB::readData(const std::string& FileName, const std::string& HeaderFile,
                         const std::string& str, int ATOMorMolecular)
{
    TlLogX& Log = TlLogX::getInstance();
    this->Atm_or_Mol = ATOMorMolecular ;

    this->DbFile = FileName;
    this->HdFile = HeaderFile;
    this->DbLabel = str;

    if (this->Atm_or_Mol == 0) {
        //  Log <<"AuxName  = [ "<< DbLabel << " ]"<< "\n";
    } else if (this->Atm_or_Mol == 1) {
        //  Log <<"DbLabel  = [ "<< DbLabel << " ]"<< "\n";
        //int passline = serchLine(HdFile,DbLabel);
    } else {
        Log<<"Bad appointment parameter : Atm_or_Mol"<<"\n";
        CnErr.abort();
    }

    std::ifstream fi;
    fi.open(this->DbFile.c_str(), std::ios::in);
    if (!fi) {
        Log <<"Cannot Open "<< DbFile << "\n";
        CnErr.abort();
    }

    this->tempflagw = 0;

    for (;;) {
        //       fi.getline(line, LineBuf);
        std::string line;
        std::getline(fi, line);

        if (fi.eof()) {
            this->tempflagw = 1;
            goto skip;
        }

        if (this->DbLabel == line) {
            if (this->Atm_or_Mol == 0) {
                this->DBM.TotalAtomNum = 1;
            } else if (this->Atm_or_Mol == 1) {
                fi >> this->DBM.TotalAtomNum;
            }
            //Log<<"DBM.TotalAtomNum = "<<DBM.TotalAtomNum<<"\n";

            //this->getMemory();
            this->DBM.AtomData.resize(this->DBM.TotalAtomNum);

            if (this->Atm_or_Mol == 1) {
                // read ( Num  Atom  X  Y  Z )
                for (int i=0; i < this->DBM.TotalAtomNum; i++) {
                    int tmp_int;
                    fi >> tmp_int
                    >> this->DBM.AtomData[i].Atom
                    >> this->DBM.AtomData[i].X
                    >> this->DBM.AtomData[i].Y
                    >> this->DBM.AtomData[i].Z;
                }

                // read Connect information.
                // read ( AtomNum  ConnectNum  Connect_AtomNumber* )
                for (int i = 0; i < this->DBM.TotalAtomNum; i++) {
                    int dumynum = 0;
                    fi >> dumynum;  // atom-number
                    int dConnectNum = 0;
                    fi >> dConnectNum; //this->DBM.AtomData[i].ConectNum;   // Conect no kazu
                    this->DBM.AtomData[i].Conection.resize(dConnectNum);
                    for (int j = 0; j < dConnectNum; j++) {
                        fi >> DBM.AtomData[i].Conection[j];
                    }
                }

                // read Touka_Atom information.
                // read ( TotalEqualPariNum ).
                // read ( PairNum   Pair_AtomNumber* ).
                //int sum = 0;
                int dEqualPairNum = 0;
                fi >> dEqualPairNum; //this->DBM.equalPairNum;
                this->DBM.equalPaircount.resize(dEqualPairNum);
                for (int i = 0; i < dEqualPairNum; i++) {
                    int dEqualPairCount = 0;
                    fi >> dEqualPairCount; //
                    DBM.equalPaircount[i] = dEqualPairCount;

                    for (int j = 0; j < dEqualPairCount; j++) {
                        int nEqualPair = 0;
                        this->DBM.equalPair.push_back(nEqualPair);
                        //fi >> DBM.equalPair[sum + j];
                    }
                    //sum += DBM.equalPaircount[i] ;
                }

                DBM.Order.resize(this->DBM.TotalAtomNum);
                for (int i = 0; i < this->DBM.TotalAtomNum; i++) {
                    fi >> DBM.Order[i];
                }
            }//if(Atm....)

            for (int i = 0; i < this->DBM.TotalAtomNum; i++) {
                if (this->Atm_or_Mol == 1) {
                    int tmp_int = 0;
                    fi >> tmp_int ;
                    //fi.getline(line,LineBuf);
                    std::getline(fi, line);
                    //strcpy(dumyString,line);
                    std::string dumyString = line;
                    //length = strlen(dumyString);
                    int length = dumyString.length();
                    int spaceCharNum = this->Count_HeadSpaceNum(dumyString);
                    this->charShift(dumyString, spaceCharNum);
                    this->add_strend(dumyString,length-spaceCharNum);
                    this->DBM.AtomData[i].AuxSetName = dumyString;
                    //Log<<"AUX_SET = ["<<DBM.AtomData[i].AuxSetName<<;
                    //Log<<"]"<<"\n";
                    if (tmp_int != i) {
                        Log<<"Bad GuessDatabase Format";
                        Log<<"(each Atom's data)"<<"\n";
                        Log<<"Database_File = ["<<DbFile<<"]"<<"\n";
                        Log<<"DB_Label = ["<<DbLabel<<"]"<<"\n";
                        Log<<"inDEBUG tmp_int ="<<tmp_int<<"\n";
                        Log<<"inDEBUG i = "<<i<<"\n";
                        Log<<"\7"<<"\n";
                        CnErr.abort();
                    }
                }//if(Atm....)

                //------ read  Rou_Data -----------------------------
                {
                    this->DBM.AtomData[i].rouDbShellNum.resize(4);
                    fi >> this->DBM.AtomData[i].rouDbShellNum[0]
                    >> this->DBM.AtomData[i].rouDbShellNum[1]
                    >> this->DBM.AtomData[i].rouDbShellNum[2]
                    >> this->DBM.AtomData[i].rouDbShellNum[3] ;

                    const int Snum   = DBM.AtomData[i].rouDbShellNum[0];
                    const int Pnum   = DBM.AtomData[i].rouDbShellNum[1];
                    const int Dnum   = DBM.AtomData[i].rouDbShellNum[2];
                    const int SPDnum = DBM.AtomData[i].rouDbShellNum[3];

#if VER1_0
                    this->DBM.AtomData[i].Routerm =
                        Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 ;
#endif

#if VER1_1
                    this->DBM.AtomData[i].Routerm =
                        Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 + SPDnum ;
#endif

                    this->DBM.AtomData[i].CoefRou.resize(this->DBM.AtomData[i].Routerm);
                    for (int j = 0; j < this->DBM.AtomData[i].Routerm; j++) {
                        fi >> DBM.AtomData[i].CoefRou[j] ;
                    }
                    this->DBM.TotalRouterm += this->DBM.AtomData[i].Routerm;
                }
                //------ read  Myu_Data -----------------------------
                {
                    this->DBM.AtomData[i].myuDbShellNum.resize(4);
                    fi >> this->DBM.AtomData[i].myuDbShellNum[0]
                    >> this->DBM.AtomData[i].myuDbShellNum[1]
                    >> this->DBM.AtomData[i].myuDbShellNum[2]
                    >> this->DBM.AtomData[i].myuDbShellNum[3] ;

                    const int Snum   = DBM.AtomData[i].myuDbShellNum[0];
                    const int Pnum   = DBM.AtomData[i].myuDbShellNum[1];
                    const int Dnum   = DBM.AtomData[i].myuDbShellNum[2];
                    const int SPDnum = DBM.AtomData[i].myuDbShellNum[3];

#if VER1_0
                    this->DBM.AtomData[i].Myuterm =
                        Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 ;
#endif

#if VER1_1
                    this->DBM.AtomData[i].Myuterm =
                        Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 + SPDnum ;
#endif

                    this->DBM.AtomData[i].CoefMyu.resize(this->DBM.AtomData[i].Myuterm);
                    for (int j = 0; j < this->DBM.AtomData[i].Myuterm; j++) {
                        fi >> DBM.AtomData[i].CoefMyu[j];
                    }
                    this->DBM.TotalMyuterm += this->DBM.AtomData[i].Myuterm;
                }
                //------ read  Nyu_Data -----------------------------
                {
                    this->DBM.AtomData[i].nyuDbShellNum.resize(4);
                    fi >> this->DBM.AtomData[i].nyuDbShellNum[0]
                    >> this->DBM.AtomData[i].nyuDbShellNum[1]
                    >> this->DBM.AtomData[i].nyuDbShellNum[2]
                    >> this->DBM.AtomData[i].nyuDbShellNum[3] ;

                    const int Snum   = DBM.AtomData[i].nyuDbShellNum[0];
                    const int Pnum   = DBM.AtomData[i].nyuDbShellNum[1];
                    const int Dnum   = DBM.AtomData[i].nyuDbShellNum[2];
                    const int SPDnum = DBM.AtomData[i].nyuDbShellNum[3];

#if VER1_0
                    this->DBM.AtomData[i].Nyuterm =
                        Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 ;
#endif

#if VER1_1
                    this->DBM.AtomData[i].Nyuterm =
                        Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 + SPDnum ;
#endif

                    this->DBM.AtomData[i].CoefNyu.resize(this->DBM.AtomData[i].Nyuterm);
                    for (int j = 0; j < this->DBM.AtomData[i].Nyuterm; j++) {
                        fi >> this->DBM.AtomData[i].CoefNyu[j];
                    }
                    this->DBM.TotalNyuterm += this->DBM.AtomData[i].Nyuterm;
                }
            } // for i_loop;

            fi.close();
            goto skip;
        } //  strcmp(DbLabel,line)==0 no line;
    } // for(;;) no line;


skip:
    ;
    return 0;
}


void Fl_GuessDB::Print() const
{
    TlLogX& Log = TlLogX::getInstance();

    Log<<"\n"<<"\n"<<"Print Start for DEBUG"<<"\n"<<"\n"<<"\n";
//     int Snum,Pnum,Dnum,SPDnum,i,k;
    int DRouterm,DMyuterm,DNyuterm;;

    int TotalAtomNum = DBM.TotalAtomNum ;
    Log <<"TotalAtomNum = "<< TotalAtomNum << "\n";

    if (this->Atm_or_Mol) {
        Log<<"\n"<<" ### ATOM DATA ###"<<"\n";
        for (int j = 0; j < TotalAtomNum; j++) {
            Log << j                     << "  "
            << DBM.AtomData[j].Atom << "   "
            << DBM.AtomData[j].X << "   "
            << DBM.AtomData[j].Y << "   "
            << DBM.AtomData[j].Z << "   " << "\n" ;
        }

        Log<<"  ### Connection Table ### "<<"\n"<<"\n";
        for (int i=0; i<TotalAtomNum; i++) {
            Log<< i <<"   "<<DBM.AtomData[i].ConectNum<<"      ";
            for (int j=0; j<DBM.AtomData[i].ConectNum; j++) {
                Log<< DBM.AtomData[i].Conection[j] <<"   ";
            }
            Log<<"\n";
        }
        Log<<"\n";

        int sum=0;
        Log<<" ### EqualAtomData ### "<<"\n"<<"\n";
        Log<<" Equal atom pair Num = "<< DBM.equalPairNum << "\n";
        for (int i=0; i<DBM.equalPairNum; i++) {
            Log<< DBM.equalPaircount[i] <<"    ";
            for (int j=0; j<DBM.equalPaircount[i]; j++) {
                Log<< DBM.equalPair[sum+j] <<"  ";
            }
            sum +=  DBM.equalPaircount[i] ;
            Log<<"\n";
        }
        Log<<"\n";

        for (int j=0; j<TotalAtomNum; j++) {
            Log << DBM.Order[j] << "  " ;
        }
        Log << "\n";
    }
    //-----------------------------------------------------
    Log <<"\n"<<"Rou(s,p,d,spd)"<<"\n";
    for (int j=0; j<TotalAtomNum; j++) {
        int Snum   = DBM.AtomData[j].rouDbShellNum[0] ;
        int Pnum   = DBM.AtomData[j].rouDbShellNum[1] ;
        int Dnum   = DBM.AtomData[j].rouDbShellNum[2] ;
        int SPDnum = DBM.AtomData[j].rouDbShellNum[3] ;

        if (Atm_or_Mol==1) {
            Log<<"AuxSetName = ["<<DBM.AtomData[j].AuxSetName <<"]"<<"\n";
        }
        Log <<"ATOM_No[ "<< j <<" ]" << " "
        << Snum   << "  "
        << Pnum   << "  "
        << Dnum   << "  "
        << SPDnum << "  "
        << "\n" ;

#if VER1_0
        DRouterm = Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 ;
#endif

#if VER1_1
        DRouterm = Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 + SPDnum ;
#endif

        for (int k=0; k<DRouterm; k++) {
            Log <<"[ "<<k<<" ]"<<DBM.AtomData[j].CoefRou[k] << "\n";
        }

    }


    //-----------------------------------------------------
    Log <<"\n"<<"Myu(s,p,d,spd)"<<"\n";
    for (int j=0; j<TotalAtomNum; j++) {
        int Snum   = DBM.AtomData[j].myuDbShellNum[0] ;
        int Pnum   = DBM.AtomData[j].myuDbShellNum[1] ;
        int Dnum   = DBM.AtomData[j].myuDbShellNum[2] ;
        int SPDnum = DBM.AtomData[j].myuDbShellNum[3] ;

        Log <<"ATOM_No[ "<< j <<" ]" << " "
        << Snum   << "  "
        << Pnum   << "  "
        << Dnum   << "  "
        << SPDnum << "  "
        << "\n" ;

#if VER1_0
        DMyuterm = Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 ;
#endif

#if VER1_1
        DMyuterm = Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 + SPDnum ;
#endif

        for (int k=0; k<DMyuterm; k++) {
            Log <<"[ "<<k<<" ]"<<DBM.AtomData[j].CoefMyu[k] << "\n";
        }

    }


    //-----------------------------------------------------
    Log <<"\n"<<"Nyu(s,p,d,spd)"<<"\n";
    for (int j=0; j<TotalAtomNum; j++) {
        int Snum   = DBM.AtomData[j].nyuDbShellNum[0] ;
        int Pnum   = DBM.AtomData[j].nyuDbShellNum[1] ;
        int Dnum   = DBM.AtomData[j].nyuDbShellNum[2] ;
        int SPDnum = DBM.AtomData[j].nyuDbShellNum[3] ;

        Log <<"ATOM_No[ "<< j <<" ]" << " "
        << Snum   << "  "
        << Pnum   << "  "
        << Dnum   << "  "
        << SPDnum << "  "
        << "\n" ;

#if VER1_0
        DNyuterm = Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 ;
#endif

#if VER1_1
        DNyuterm = Snum*1 + Pnum*3 + Dnum*5 + SPDnum*9 + SPDnum ;
#endif

        for (int k=0; k<DNyuterm; k++) {
            Log <<"[ "<<k<<" ]"<<DBM.AtomData[j].CoefNyu[k] << "\n";
        }

    }
}


// Set Member Function  for  AtomDatabase
void Fl_GuessDB::openAtomDBfile(const std::string& sFile)
{
    TlLogX& Log = TlLogX::getInstance();
    this->DbFile = sFile;

    SETATOM.open(DbFile.c_str(), std::ios::app);
    if (!SETATOM) {
        Log<<"Cannot open file ["<<DbFile<<"]."<<"\n";
        Log<<"Cannot append AtomDatabase"<<"\n";
        Log<<"\7"<<"\n";
        CnErr.abort();
    }
}


void Fl_GuessDB::closeAtomDBfile()
{
    SETATOM<<"\n";
    SETATOM.close();
}

// get member function
int Fl_GuessDB::getTotalAtomNum() const
{
    return  DBM.TotalAtomNum;
}


std::string Fl_GuessDB::getAuxSetName(int num) const
{
    return DBM.AtomData[num].AuxSetName;
}


std::string Fl_GuessDB::getAtom(int num) const
{
    return DBM.AtomData[num].Atom;
}


double Fl_GuessDB::getXCoord(int num) const
{
    return  DBM.AtomData[num].X;
}


double Fl_GuessDB::getYCoord(int num)const
{
    return  DBM.AtomData[num].Y;
}


double Fl_GuessDB::getZCoord(int num)const
{
    return  DBM.AtomData[num].Z;
}


int Fl_GuessDB::getConnectNum(int  num)const
{
    return  DBM.AtomData[num].ConectNum;
}


void Fl_GuessDB::getConnection(int num, std::vector<int>& AtmNum) const
{
    assert(AtmNum.size() >= static_cast<unsigned int>(DBM.AtomData[num].ConectNum));

    for (int i=0; i < DBM.AtomData[num].ConectNum; i++) {
        AtmNum[i] = DBM.AtomData[num].Conection[i];
    }
}


int Fl_GuessDB::getTotalEqualPairNum() const
{
    return DBM.equalPairNum;
}


void Fl_GuessDB::getEqualPairCount(std::vector<int>& val) const
{
    assert(val.size() >= static_cast<unsigned int>(DBM.equalPairNum));

    for (int i=0; i < DBM.equalPairNum; i++) {
        val[i] = DBM.equalPaircount[i];
    }
}


void Fl_GuessDB::getEqualPair(std::vector<int>& val) const
{

    int sum=0;
    for (int i=0; i < DBM.equalPairNum; i++) {
        for (int j=0; j < DBM.equalPaircount[i]; j++) {
            val[ i*DBM.equalPaircount[i]+j ] = DBM.equalPair[sum+j];
        }
        sum += DBM.equalPaircount[i];
    }
}


void Fl_GuessDB::getOrder(std::vector<int>& val) const
{
    assert(val.size() >= static_cast<unsigned int>(DBM.TotalAtomNum));

    for (int i=0; i < DBM.TotalAtomNum; i++) {
        val[i] = DBM.Order[i];
    }
}


std::vector<int> Fl_GuessDB::getRouShellterm(int num) const
{
    std::vector<int> term(MaxSPDterm);

    for (int i=0; i < MaxSPDterm; i++) {
        term[i] = DBM.AtomData[num].rouDbShellNum[i];
    }

    return term;
}


std::vector<int> Fl_GuessDB::getMyuShellterm(int num) const
{
    std::vector<int> term(MaxSPDterm);

    for (int i=0; i<MaxSPDterm; i++) {
        term[i] = DBM.AtomData[num].myuDbShellNum[i];
    }

    return term;
}


std::vector<int> Fl_GuessDB::getNyuShellterm(int num) const
{
    std::vector<int> term(MaxSPDterm);

    for (int i=0; i<MaxSPDterm; i++) {
        term[i] = DBM.AtomData[num].nyuDbShellNum[i];
    }

    return term;
}


void Fl_GuessDB::getRouCoef(int num, std::vector<double>& val) const
{
    assert(val.size() >= static_cast<unsigned int>(DBM.AtomData[num].Routerm));

    for (int i=0; i < DBM.AtomData[num].Routerm; i++) {
        val[i] = DBM.AtomData[num].CoefRou[i] ;
    }
}


void Fl_GuessDB::getMyuCoef(int num, std::vector<double>& val) const
{
    assert(val.size() >= static_cast<unsigned int>(DBM.AtomData[num].Myuterm));

    for (int i=0; i < DBM.AtomData[num].Myuterm; i++) {
        val[i] = DBM.AtomData[num].CoefMyu[i];
    }
}


void Fl_GuessDB::getNyuCoef(int num, std::vector<double>& val) const
{
    assert(val.size() >= static_cast<unsigned int>(DBM.AtomData[num].Nyuterm));

    for (int i=0; i < DBM.AtomData[num].Nyuterm; i++) {
        val[i] = DBM.AtomData[num].CoefNyu[i];
    }
}


int Fl_GuessDB::getTotalRouterm() const
{
    return DBM.TotalRouterm;
}


int Fl_GuessDB::getTotalMyuterm() const
{
    return DBM.TotalMyuterm;
}


int Fl_GuessDB::getTotalNyuterm() const
{
    return DBM.TotalNyuterm;
}


int Fl_GuessDB::getRouterm(int num) const
{
    return  DBM.AtomData[num].Routerm;
}


int Fl_GuessDB::getMyuterm(int num) const
{
    return  DBM.AtomData[num].Myuterm;
}


int Fl_GuessDB::getNyuterm(int num) const
{
    return  DBM.AtomData[num].Nyuterm;
}


// Translation  for String
int Fl_GuessDB::Count_HeadSpaceNum(const std::string& str) const
{
    int ptr = 0;

    for (;;) {
        if (str[ptr] == ' ') {
            ptr++;
        } else {
            break;
        }
    }

    return ptr;
}


void Fl_GuessDB::charShift(std::string& str, const int SpaceNum) const
{
    int length = str.length() - SpaceNum;
    for (int i=0; i < length ; i++) {
        str[i] = str[i+SpaceNum];
    }
}


void Fl_GuessDB::add_strend(std::string& str, const int length) const
{
    str[length] = '\0';
}

