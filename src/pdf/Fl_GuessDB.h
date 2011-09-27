#ifndef FL_GUESSDB_H
#define FL_GUESSDB_H

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>

#include "CnError.h"

#define VER1_0  1  // for non add 1s
#define VER1_1  0  // for add 1s( by d_orbital ).

class Fl_GuessDB {
public:
    Fl_GuessDB();
    virtual ~Fl_GuessDB();

public:
    int gettempflag() const {
        return tempflagw;
    };

public:
    void Print() const;

    // Common memberfunction for Atom and Molecular .
    int readData(const std::string&, const std::string&, const std::string&, int);

    void openAtomDBfile(const std::string&);
    void closeAtomDBfile();

public:
    int  getHeaderAtomNum(const std::string&, const std::string&, int) const;

    int   getTotalAtomNum() const ;
    std::string getAuxSetName(int) const;
    std::string getAtom(int) const;
    double  getXCoord(int) const;
    double  getYCoord(int) const;
    double  getZCoord(int) const;

    int getConnectNum(int) const;
    void getConnection(int, std::vector<int>&) const;
    int getTotalEqualPairNum() const;
    void getEqualPairCount(std::vector<int>&) const;
    void getEqualPair(std::vector<int>&) const;
    void getOrder(std::vector<int>&) const;

    std::vector<int> getRouShellterm(int) const;
    std::vector<int> getMyuShellterm(int) const;
    std::vector<int> getNyuShellterm(int) const;

    void  getRouCoef(int, std::vector<double>&) const;
    void  getMyuCoef(int, std::vector<double>&) const;
    void  getNyuCoef(int, std::vector<double>&) const;

    int   getRouterm(int) const;
    int   getMyuterm(int) const;
    int   getNyuterm(int) const;
    int   getTotalRouterm() const;
    int   getTotalMyuterm() const;
    int   getTotalNyuterm() const;

private:
    int serchLine(const std::string&, const std::string&) const;

    int  Count_HeadSpaceNum(const std::string&) const;
    void charShift(std::string&, int) const;
    void add_strend(std::string&, int) const;

private:
    enum { MaxSPDterm = 4};

    int  Atm_or_Mol;
    std::string  DbFile;
    std::string  HdFile;
    std::string  DbLabel;

    int tempflagw;

    std::ofstream SETATOM; // AtomDatabase wo Tukuru tame no file_pointer.

    struct AtomSet {
        std::string AuxSetName;
        std::string Atom;
        double X;
        double Y;
        double Z;

        std::vector<int> rouDbShellNum;
        std::vector<int> myuDbShellNum;
        std::vector<int> nyuDbShellNum;

        int Routerm;
        int Myuterm;
        int Nyuterm;

        std::vector<double> CoefRou;
        std::vector<double> CoefMyu;
        std::vector<double> CoefNyu;

        int ConectNum;
        std::vector<int> Conection;    // conect table
    };

    struct MolecularData {
        int TotalAtomNum;
        std::vector<AtomSet> AtomData;
        int equalPairNum;
        std::vector<int> equalPaircount;
        std::vector<int> equalPair;
        std::vector<int> Order;        // itti saseru jyunban
        int TotalRouterm;
        int TotalMyuterm;
        int TotalNyuterm;
    };
    MolecularData DBM;
};

#endif // FL_GUESSDB_H

