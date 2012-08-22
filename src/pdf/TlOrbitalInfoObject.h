#ifndef TLORBITALINFOOBJECT_H
#define TLORBITALINFOOBJECT_H

#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <iostream>

#include "TlPosition.h"
#include "TlAtom.h"
#include "TlUtils.h"
#include "Fl_Gto.h"
#include "Fl_Geometry.h"

/// 軌道の情報を管理するクラス
///
/// PGTOは指数の小さい順にソートされる
class TlOrbitalInfoObject {
public:
    typedef int index_type;
    typedef long size_type;

public:
    typedef std::vector<TlAtom> Atoms;
    
    /// primitive-GaussTypeOrbitalの情報を保持する構造体
    struct PGTO {
        double dCoefficient;
        double dExponent;
    };
    typedef std::vector<PGTO> PGTOs;

    struct CGTO {
        std::string basisName;
        std::string atomSymbol; // H, He, ...
        std::string label;
        char shell;
        PGTOs pgtos;

        double factor; // sum of exp(-a)
    };

    typedef std::vector<CGTO> CGTOs;

    
    /// CGTO内のPGTOをソートする関数オブジェクト
    struct cgto_sort_functor {
        bool operator()(const PGTO& a, const PGTO& b) const {
            return (a.dExponent < b.dExponent);
        }
    };

protected:
    /// 軌道の情報を保持する構造体
    struct Orbital {
        int shellType;      // 0: s, 1: p, 2:d
        int basisType; // s=0, px=1, py=2, ...
        int atomIndex; // 対応するatoms_のインデックス
        int cgtoIndex; // 対応するcgtos_のインデックス
        
        double prefactor;
        double prefactor_gradx;
        double prefactor_grady;
        double prefactor_gradz;
    };

public:
    /// 原子数を返す
    int getNumOfAtoms() const {
        return this->atoms_.size();
    }
    
    /// 系の全軌道数を返す
    index_type getNumOfOrbitals() const {
        return this->orbitals_.size();
    }


    /// 与えられた軌道のシェルの型を返す.
    ///
    /// @param[in] index 軌道番号
    /// @retval 0 s型
    /// @retval 1 p型
    /// @retval 2 d型
    int getShellType(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->orbitals_[AO].shellType;
    }


    static int getMaxShellType() {
        return TlOrbitalInfoObject::MAX_SHELL_TYPE;
    }

    /// シェル毎にグループ化された軌道の最初の軌道番号の配列を返す
    ///
    /// 例えば 0:s, 1:s, 2:px, 3:py, 4:pz, 5:s, ...ならば
    /// 0, 1, 2, 5, ...を返す
    std::vector<index_type> getStartIndexArrayOfShellGroup() const;
    
    /// 与えられた軌道のシェル(角運動量)の型を返す.
    ///
    /// @param[in] index 軌道番号
    /// @retval 0 s
    /// @retval 1 px
    /// @retval 2 py
    /// @retval 3 pz
    /// @retval 4 dxy
    /// @retval 5 dxz
    /// @retval 6 dyz
    /// @retval 7 dxx-yy
    /// @retval 8 dzz
    int getBasisType(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->orbitals_[AO].basisType;
    }

    /// 与えられたインデックスのシェルインデックスを返す
    index_type getShellIndex(const index_type AO) const {
        static const int basisTypeBase[] = {0, 1, 4}; // s, px, dxy
        const int shellType = this->getShellType(AO);
        const int basisType = this->getBasisType(AO);
        return AO - (basisType - basisTypeBase[shellType]);
    }

    std::string getBasisTypeName(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return std::string(TlOrbitalInfoObject::basisTypeNameTbl_[this->getBasisType(AO)]);
    }
    
    /// 与えられた軌道が属する原子インデックス(原子番号ではない)を返す.
    ///
    /// @param[in] index 軌道番号
    /// @return 原子のインデックス
    int getAtomIndex(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->orbitals_[AO].atomIndex;
    }


    /// 与えられた軌道が属する原子名を返す.
    ///
    /// @param[in] index 軌道番号
    /// @return 原子名
    std::string getAtomName(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->atoms_[this->getAtomIndex(AO)].getSymbol();
    }


    std::string getAtomLabel(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->atoms_[this->getAtomIndex(AO)].getName();
    }
   

    /// 与えられた軌道が属する原子の位置を返す.
    ///
    /// @param[in] index 軌道番号
    /// @return TlPositionオブジェクト
    TlPosition getPosition(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->atoms_[this->getAtomIndex(AO)].getPosition();
    }

    /// 与えられた軌道が属する原子の電荷を返す.
    double getAtomCharge(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->atoms_[this->getAtomIndex(AO)].getCharge();
    }

    /// 短縮数を返す
    int getCgtoContraction(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->cgtos_[this->getCgtoIndex(AO)].pgtos.size();
    }


    double getCoefficient(const index_type AO, const int pGTO) const {
        assert(AO < this->getNumOfOrbitals());
        assert(pGTO < this->getCgtoContraction(AO));
        return this->cgtos_[this->getCgtoIndex(AO)].pgtos[pGTO].dCoefficient;
    }


    double getExponent(const index_type AO, const int pGTO) const {
        assert(AO < this->getNumOfOrbitals());
        assert(pGTO < this->getCgtoContraction(AO));
        return this->cgtos_[this->getCgtoIndex(AO)].pgtos[pGTO].dExponent;
    }


    /// 最小の pGTO 指数を返す
    double getMinExponent(const index_type AO) const {
        // PGTOは小さい順にソートされているので、最初が最小
        return this->getExponent(AO, 0);
    }

    
    double getSumOfCoef(const index_type AO) const {
        return this->cgtos_[this->getCgtoIndex(AO)].factor;
    }
    
    

public:
    // for debug
    template<typename T>
    void printCGTOs(T& out) const;

    template<typename T>
    void printPGTOs(const PGTOs rPGTOs, T& out) const;

    
protected:
    /// コンストラクタ
    ///
    TlOrbitalInfoObject();

    /// デストラクタ
    virtual ~TlOrbitalInfoObject();

    
protected:
    int getCgtoIndex(const index_type AO) const {
        assert(AO < this->getNumOfOrbitals());
        return this->orbitals_[AO].cgtoIndex;
    }

    void setCGTO(const Fl_Gto& flGto);
    void setCGTO_coulomb(const Fl_Gto& flGto);
    void setAtoms(const Fl_Geometry& flGeom);
    void makeOrbitalTable();

protected:
    static void getPrefactor_S(const TlPosition& pos, double* pValue);
    static void getPrefactor_PX(const TlPosition& pos, double* pValue);
    static void getPrefactor_PY(const TlPosition& pos, double* pValue);
    static void getPrefactor_PZ(const TlPosition& pos, double* pValue);
    static void getPrefactor_DXY(const TlPosition& pos, double* pValue);
    static void getPrefactor_DXZ(const TlPosition& pos, double* pValue);
    static void getPrefactor_DYZ(const TlPosition& pos, double* pValue);
    static void getPrefactor_DXX_YY(const TlPosition& pos, double* pValue);
    static void getPrefactor_DZZ(const TlPosition& pos, double* pValue);
    
    static void getGradPrefactor_S(const TlPosition& pos, const double alpha,
                                   double* pDX, double* pDY, double* pDZ);
    static void getGradPrefactor_PX(const TlPosition& pos, const double alpha,
                                    double* pDX, double* pDY, double* pDZ);
    static void getGradPrefactor_PY(const TlPosition& pos, const double alpha,
                                    double* pDX, double* pDY, double* pDZ);
    static void getGradPrefactor_PZ(const TlPosition& pos, const double alpha,
                                    double* pDX, double* pDY, double* pDZ);
    static void getGradPrefactor_DXY(const TlPosition& pos, const double alpha,
                                     double* pDX, double* pDY, double* pDZ);
    static void getGradPrefactor_DXZ(const TlPosition& pos, const double alpha,
                                     double* pDX, double* pDY, double* pDZ);
    static void getGradPrefactor_DYZ(const TlPosition& pos, const double alpha,
                                     double* pDX, double* pDY, double* pDZ);
    static void getGradPrefactor_DXX_YY(const TlPosition& pos, const double alpha,
                                        double* pDX, double* pDY, double* pDZ);
    static void getGradPrefactor_DZZ(const TlPosition& pos, const double alpha,
                                     double* pDX, double* pDY, double* pDZ);
    
public:
    /// 扱うことのできるシェルの型の総数
    /// @sa getShellType()
    static const int MAX_SHELL_TYPE;

    static const double INV_SQRT3;
    static const double MAX_EXPONENT;

    static const char* basisTypeNameTbl_[];
    
protected:
    Atoms atoms_;
    CGTOs cgtos_;
    std::vector<Orbital> orbitals_;
};


template<typename T>
void TlOrbitalInfoObject::printCGTOs(T& out) const
{
    const int nMaxCgtoSize = this->cgtos_.size();
    out << "number of CGTO = " << nMaxCgtoSize << std::endl;

    for (int i = 0; i < nMaxCgtoSize; ++i) {
        out << TlUtils::format("--- CGTO ID = %d\n", i);
        this->printPGTOs(this->cgtos_[i].pgtos, out);
        out << std::endl;
    }
}


template<typename T>
void TlOrbitalInfoObject::printPGTOs(const PGTOs rPGTOs, T& out) const
{
    const int nMaxPgtoSize = rPGTOs.size();

    for (int i = 0; i < nMaxPgtoSize; ++i) {
        out << TlUtils::format("exp= % .9e, coef= % .9e\n", rPGTOs[i].dExponent, rPGTOs[i].dCoefficient);
    }

}

#endif // TLORBITALINFOOBJECT_H
