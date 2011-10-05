#ifndef FL_DB_BASIS_H
#define FL_DB_BASIS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Fl_Database.h"
#include "CnError.h"

// Basisset名に対応した値を、データベースから
// 取ってくるクラス。データベースの名前は今の所basis2である
class Fl_Db_Basis : public Fl_Database {
private:
    struct Cgto {
        char shell;              // Shelltype { s , p , d };
        double scalefactor;        // Scalefactor;
        double normalizedfactor;   // Normalazedfactor;
        struct Pgto {
            double exponent;          // Expornent;
            double coefficient;        // Coefficient;
        };
        std::vector<Pgto> contraction;
    };

public:
    Fl_Db_Basis(const std::string&);
    virtual ~Fl_Db_Basis();

public:
    int getTotalcgto() const {
        return this->cgto.size();
    }

    int getrouSnum() const {
        int nAnswer = 0;
        for (unsigned int i = 0; i < this->rhoCgto.size(); ++i) {
            if (this->rhoCgto[i].shell == 's') {
                ++nAnswer;
            }
        }
        return static_cast<int>(nAnswer);
    }

    int getrouPnum() const {
        int nAnswer = 0;
        for (unsigned int i = 0; i < this->rhoCgto.size(); ++i) {
            if (this->rhoCgto[i].shell == 'p') {
                ++nAnswer;
            }
        }
        return static_cast<int>(nAnswer);
    }

    int getrouDnum() const {
        int nAnswer = 0;
        for (unsigned int i = 0; i < this->rhoCgto.size(); ++i) {
            if (this->rhoCgto[i].shell == 'd') {
                ++nAnswer;
            }
        }

        return static_cast<int>(nAnswer);
    }

    int getrouTotalnum() const {
        return this->rhoCgto.size();
    }

    int getmyuSnum() const {
        int nAnswer = 0;
        for (unsigned int i = 0; i < this->myuCgto.size(); ++i) {
            if (this->myuCgto[i].shell == 's') {
                ++nAnswer;
            }
        }
        return static_cast<int>(nAnswer);
    }

    int getmyuPnum() const {
        int nAnswer = 0;
        for (unsigned int i = 0; i < this->myuCgto.size(); ++i) {
            if (this->myuCgto[i].shell == 'p') {
                ++nAnswer;
            }
        }
        return static_cast<int>(nAnswer);
    }

    int getmyuDnum() const {
        int nAnswer = 0;
        for (unsigned int i = 0; i < this->myuCgto.size(); ++i) {
            if (this->myuCgto[i].shell == 'd') {
                ++nAnswer;
            }
        }
        return static_cast<int>(nAnswer);
    }

    int getmyuTotalnum() const {
        return this->myuCgto.size();
    }

    char getShell(int s) const {
        return this->cgto[s].shell;
    }

    double getScalefactor(int s) {
        return   cgto[s].scalefactor;
    }

    double getNormalizedfactor(int s) {
        return   cgto[s].normalizedfactor;
    }

    int getContraction(int s) {
        return this->cgto[s].contraction.size();
    }

    double getExpornent(int s,int p) {
        return this->cgto[s].contraction[p].exponent;
    }

    double getCoefficient(int s,int p) {
        return this->cgto[s].contraction[p].coefficient;
    }

    char getrouShell(int s) {
        return this->rhoCgto[s].shell;
    }

    double getrouScalefactor(int s) {
        return this->rhoCgto[s].scalefactor;
    }

    double getrouNormalizedfactor(int s) {
        return this->rhoCgto[s].normalizedfactor;
    }

    int getrouContraction(int s) {
        return this->rhoCgto[s].contraction.size();
    }

    double getrouExpornent(int s,int p) {
        return this->rhoCgto[s].contraction[p].exponent;
    }

    double getrouCoefficient(int s,int p) {
        return this->rhoCgto[s].contraction[p].coefficient;
    }

    char getmyuShell(int s) {
        return this->myuCgto[s].shell;
    }

    double getmyuScalefactor(int s) {
        return this->myuCgto[s].scalefactor;
    }

    double getmyuNormalizedfactor(int s) {
        return this->myuCgto[s].normalizedfactor;
    }

    int getmyuContraction(int s) {
        return this->myuCgto[s].contraction.size();
    }

    double getmyuExpornent(int s,int p) {
        return this->myuCgto[s].contraction[p].exponent;
    }

    double getmyuCoefficient(int s,int p) {
        return this->myuCgto[s].contraction[p].coefficient;
    }

private:
    void setData();
    void read_cGTOs(const int numOfCGTOs,
                    const std::string& shellType,
                    std::ifstream& fi,
                    std::vector<Cgto>* pCGTOs);

    std::vector<std::string> getWords(const std::string& str);

private:
    std::string m_sBasisName;

    std::vector<Cgto> cgto;
    std::vector<Cgto> rhoCgto;
    std::vector<Cgto> myuCgto;
};

#endif // FL_DB_BASIS_H



