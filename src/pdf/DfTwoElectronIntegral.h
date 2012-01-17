#ifndef DFTWOELECTRONINTEGRAL_H
#define DFTWOELECTRONINTEGRAL_H

#ifdef HAVE_CONFIG_H
#include "config.h"    // this file created by autotools
#endif // HAVE_CONFIG_H

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <cmath>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <functional>

#include "DfObject.h"
#include "DfTEI.h"

#include "TlSimpleVector.h"
#include "TlPosition.h"
#include "TlSymmetricMatrix.h"
#include "TlSparseSymmetricMatrix.h"
#include "TlOrbitalInfo.h"

/**
 *  4中心積分クラス
 *
 *  制限：
 *  o intは32bit以上であること
 *  o 基底関数の数は4,294,967,296の半分、すなわち2,147,483,648以下であること
 */
class DfTwoElectronIntegral : public DfObject {
protected:
    class ShellPairIndex {
    public:
        ShellPairIndex(int i, int j) : m_IShell(i), m_JShell(j) {
            if (m_IShell > m_JShell) {
                std::swap(m_IShell, m_JShell);
            }
        }

    public:
        friend bool operator<(const ShellPairIndex& a, const ShellPairIndex& b) {
            return ((a.m_IShell < b.m_IShell) ||
                    ((a.m_IShell == b.m_IShell) && (a.m_JShell < b.m_JShell)));
        }

    private:
        int m_IShell;
        int m_JShell;
    };

    // for I-K shell pair
    struct IKShellPair {
        explicit IKShellPair(int i =0, int k =0): nIShell(i), nKShell(k) {
        }

        int nIShell;
        int nKShell;
    };


public:
    DfTwoElectronIntegral(TlSerializeData* pPdfParam);
    virtual ~DfTwoElectronIntegral();

public:
    virtual void getContractKMatrixByIntegralDriven(const TlSymmetricMatrix& P,
                                                    TlSymmetricMatrix* pK);
    virtual void getContractKMatrixByRTmethod(const TlSymmetricMatrix& P,
                                              TlSymmetricMatrix* pK);

protected:
    // this->m_IKShellPairListから計算すべきIK shellの組を返す
    // 並列化用(シリアルでは全部返す)
    virtual std::vector<IKShellPair> getLocalIKShellPairList(const int nIKShellPairIndex);

    virtual void finalize(TlSymmetricMatrix& K);

protected:
    /// RT法によりFockの交換項を求める
    ///
    /// @param[in] P 密度行列
    void getContractKMatrixByRTmethodCore(const TlSymmetricMatrix& P,
                                          TlSymmetricMatrix* pK);

    /// s, p, d各軌道の先頭のインデックスを検索し、m_ShellList[nShellType]に格納する
    ///
    /// m_ShellList[nShellType] には(nShellType)に対応する軌道indexが順に入っている
    void makeShellList();

    /// (pq|**)において、計算すべきpqの組を計算する
    ///
    /// m_DistributionShellList[Shell][ShellType]には
    /// Shell番目の軌道に対してShellType型の軌道が積の大きい順に入っている
    void screenDistributions();

    /// 計算すべきi-shell, k-shellの組み合わせを計算する
    ///
    /// m_IKShellPairList[nIKShellPairIndex]を作成する.
    /// リストを求める際にカットオフは行わない
    void getShellList_density();

    /// 計算すべきi-shell, k-shellの組み合わせを計算する
    ///
    /// m_IKShellPairList[nIKShellPairIndex]を作成する.
    /// リストを求める際にカットオフを行う
    void getShellList_density_nocut();

    std::vector<int> getShellList(const std::vector<int>& rTargetShellList,
                                  const int nMaxShell);

    /// integral-driven
    void getContractKMatrixByIntegralDrivenCore(const TlSymmetricMatrix& P,
                                                TlSymmetricMatrix* pK);

    /// integral-drive
    void IntegralDrivenCore(const unsigned int eriType,
                            const int iShell, const int kShell,
                            const TlSymmetricMatrix& P,
                            TlSymmetricMatrix* pK);


    /// RT計算法のコア部分
    ///
    /// @param P_i[in] P(i, *)の要素を格納した密度行列
    /// @param P_k[in] P(k, *)の要素を格納した密度行列
    template<class DensitySymmetricMatrixType, class ExchangeSymmetricMatrixType>
    void RT_Core(const unsigned int eriType,
                 const int iShell, const int kShell,
                 const DensitySymmetricMatrixType& P_i,
                 const DensitySymmetricMatrixType& P_k,
                 ExchangeSymmetricMatrixType* pK);

    /// 計算すべきi-shell, j-shellの組み合わせ(<ij|**>)を計算する
    std::vector<int> getShellList_distribute(const std::vector<int>& rShellList,
                                             int nTargetShell, int nTargetShellType);
    
    std::size_t getShellList_distribute(const std::vector<int>& rShellList,
                                        const int nTargetShell, const int nTargetShellType,
                                        index_type* pOutShellList);

    double getMaxDensityMatrixElement(int nShell1, int nStep1,
                                      int nShell2, int nStep2,
                                      const TlSymmetricMatrix& P);

public:
    DfTEI::ShellPair getShellPair(int nIShell, int nJShell);
    void prepare_ERI();

protected:
    /// 積分値と密度行列要素の積をK行列に足しこむ
    void storeKmatrixByIntegralDriven(const int ishell, const int istep,
                                      const int jshell, const int jstep,
                                      const int kshell, const int kstep,
                                      const int lshell, const int lstep,
                                      const TlMatrix& P, const double* pERI,
                                      TlMatrix* pKmat);

    /// 積分値と密度行列要素の積をK行列に足しこむ
    template <class DensityMatrixType, class KMatrixType>
    void storeKmatrixByRTmethod(const int ishell, const int istep,
                                const int jshell, const int jstep,
                                const int kshell, const int kstep,
                                const int lshell, const int lstep,
                                const DensityMatrixType& P_i,
                                const DensityMatrixType& P_k,
                                const double* pERI,
                                KMatrixType* pKmat);
protected:
    void cutoffReport(const std::string& shell, const long cutoffCount, const long totalCount);


protected:
    static const double CK;
//   static const double INV_SQRT3;
    static const int MAX_SHELL_TYPE; // =3 (s, p, d)
    static const double CONTRIBUTE_COEF;

    /// 軌道の情報を保持したオブジェクト
    TlOrbitalInfo m_TlOrbitalInfo;

protected:
    /// ユーザー入力のカットオフ値
    double m_dInputCutoffThreshold;

    /// 動的変化によって得られたカットオフ値
    double m_dCutoffThreshold;

    // for primitive shell integral
    double m_dPrimitiveIntegralsCutoffThreshold;

    ///
    bool isDensityCutoff_;

    std::size_t JShellListSize_; /// pJShellList_の有効サイズ
    index_type* pJShellList_; /// Iインデックスに対する計算すべきJのリスト
    std::size_t LShellListSize_; /// pJShellList_の有効サイズ
    index_type* pLShellList_; /// Iインデックスに対する計算すべきJのリスト
    
    /// length scale paramter
    double lengthScaleParameter_;

    std::map<int, std::vector<int> > m_ShellList; // m_ShellList[nShellType]
    std::map<int, std::map<int, std::vector<int> > > m_DistributionShellList;

    std::map<ShellPairIndex, DfTEI::ShellPair> m_storedShellPair;

protected:
    // for primitive shell

    /// PGTOの計算が十分に小さいか判断し, その結果を返す
    //bool isPrimitiveShellsCutoff() const;

protected:
    // for parallel
    std::vector<std::vector<IKShellPair> > m_IKShellPairList; // indexは(4*iShellType +kShellType)

protected:
    // for debug
//   template <typename T>
//   void printIntegrals(int nIShell, int nJShell, int nKShell, int nLShell,
//            T& out) const;
    bool isOutputTEI_;
    std::string outputTEI_path_; // ２電子積分の積分値を出力するファイル名
    //bool isOutputHGP_;

    //bool isUseNewEngine_;
};


// ---------------------------------------------------------------------
// templates
// ---------------------------------------------------------------------

template<class DensitySymmetricMatrixType, class ExchangeSymmetricMatrixType>
void DfTwoElectronIntegral::RT_Core(const unsigned int eriType,
                                    const int iShell, const int kShell,
                                    const DensitySymmetricMatrixType& P_i,
                                    const DensitySymmetricMatrixType& P_k,
                                    ExchangeSymmetricMatrixType* pK)
{
    assert(pK != NULL);

    // shell type
    const int iShellType = int((eriType >> 6) & 3);
    const int jShellType = int((eriType >> 4) & 3);
    const int kShellType = int((eriType >> 2) & 3);
    const int lShellType = int(eriType        & 3);

    // steps
    const int iSteps = 2*iShellType +1;
    const int jSteps = 2*jShellType +1;
    const int kSteps = 2*kShellType +1;
    const int lSteps = 2*lShellType +1;

    // j-Shell (from 0 to i-shell)
//     const std::vector<int> JShellList
//     = this->getShellList_distribute(this->getShellList(this->m_ShellList[jShellType], iShell +1),
//                                     iShell, jShellType);
    // 行列の次元数はINT_MAXが最大なので、
    // 受け取る要素数(indexJ_max)はINT_MAX以下になる。
//     const int indexJ_max = JShellList.size();
    const index_type indexJ_max = this->getShellList_distribute(this->getShellList(this->m_ShellList[jShellType], iShell +1),
                                                                iShell, jShellType,
                                                                this->pJShellList_);

#pragma omp parallel
    {
        DfTEI dfTEI(this->m_dPrimitiveIntegralsCutoffThreshold);

#pragma omp for schedule(runtime)
        for (index_type indexJ = 0; indexJ < indexJ_max; ++indexJ) {
            //const int jShell = JShellList[indexJ];
            const index_type jShell = this->pJShellList_[indexJ];
            const DfTEI::ShellPair IJ = this->getShellPair(iShell, jShell);
            const double dIJIJ = IJ.dSchwartzValue;
            //const double dP_jk = this->getMaxDensityMatrixElement(jShell, jSteps, kShell, kSteps, P);

            // l-shell loop from 0 to i-shell
//             const std::vector<int> LShellList =
//                 this->getShellList_distribute(this->getShellList(this->m_ShellList[lShellType], iShell +1),
//                                               kShell, lShellType);
//             const std::size_t indexL_max = LShellList.size();
            const index_type indexL_max = this->getShellList_distribute(this->getShellList(this->m_ShellList[lShellType], iShell +1),
                                                                        kShell, lShellType,
                                                                        this->pLShellList_);
            for (index_type indexL = 0; indexL < indexL_max; ++indexL) {
                //const int lShell = LShellList[indexL];
                const index_type lShell = this->pLShellList_[indexL];
                const DfTEI::ShellPair KL = this->getShellPair(kShell, lShell);
                const double dKLKL = KL.dSchwartzValue;
                //const double dP_il = this->getMaxDensityMatrixElement(iShell, iSteps, lShell, lSteps, P);

                // shwartz cutoff
                //const double dMaxP = std::max(dP_jk, dP_il);
                //if ((dMaxP * dIJIJ * dKLKL) < this->m_dCutoffThreshold) {
                if ((dIJIJ * dKLKL) < this->m_dCutoffThreshold) {
                    continue;
                }

                // calc
                dfTEI.calc(eriType, IJ, KL);

                // store K matrix
#pragma omp critical (storeK_RT)
                {
                    this->storeKmatrixByRTmethod(iShell, iSteps, jShell, jSteps,
                                                 kShell, kSteps, lShell, lSteps,
                                                 P_i, P_k,
                                                 dfTEI.ERI, pK);
                }
            }

        }
    }
}


// 密度行列 P と、積分値 xxxx から得られた値を、K 行列に足しこむ
// RT method 用
template <class DensityMatrixType, class KMatrixType>
void DfTwoElectronIntegral::storeKmatrixByRTmethod(const int ishell, const int istep, const int jshell, const int jstep,
                                                   const int kshell, const int kstep, const int lshell, const int lstep,
                                                   const DensityMatrixType& P_i, const DensityMatrixType& P_k,
                                                   const double* pERI, KMatrixType* pKmat)
{
    assert(pKmat != NULL);

    int index = 0;
    for (int i = 0; i < istep; ++i) {
        const int I = ishell + i;

        for (int j = 0; j < jstep; ++j) {
            const int J = jshell + j;
            if (J > I) {
                // bypass
                index += (kstep * lstep);
                continue;
            }

            for (int k = 0; k < kstep; ++k) {
                const int K = kshell + k;
                if (K > I) {
                    // bypass
                    index += lstep;
                    continue;
                }

                for (int l = 0; l < lstep; ++l) {
                    const int L = lshell + l;
                    if (L > I) {
                        // bypass
                        ++index;
                        continue;
                    }

                    const double value = -1.0 * pERI[index];

                    // 代入
                    //(*pKmat)(I, L) += P(K, J) * value;
                    pKmat->add(I, L, P_k.get(K, J) * value);
//                     if (P_k.hasKey(K, J) == false) {
//                         std::cerr << TlUtils::format("P_k(%3d, %3d) [NG] value = %f, kshell=%d, kstep=%d.",
//                                                      K, J, P_k.get(K, J), kshell, kstep)
//                                   << std::endl;
//                     }

                    // (ij|kl)の対称性 1
                    if ((K <= J) && (I > J)) {
                        //(*pKmat)(J, K) += P(I, L) * value;
                        pKmat->add(J, K, P_i.get(I, L) * value);
//                         if (P_i.hasKey(I, L) == false) {
//                             std::cerr << TlUtils::format("P_i(%3d, %3d) not found. value = %f, ishell=%d, istep=%d.",
//                                                          I, L, P_i.get(I, L), ishell, istep)
//                                       << std::endl;
//                         }
                    }

                    // (ij|kl)の対称性 2
                    if ((K >= J) && (I > K) && (I > L)) {
                        //pKmat->add(K, J, P(I, L) * value);
                        pKmat->add(J, K, P_i.get(I, L) * value);
//                         if (P_i.hasKey(I, L) == false) {
//                             std::cerr << TlUtils::format("P_i(%3d, %3d) not found. value = %f, ishell=%d, istep=%d.",
//                                                          I, L, P_i.get(I, L), ishell, istep)
//                                       << std::endl;
//                         }
                    }

                    ++index;
                }
            }
        }
    }
}


#endif // DFTWOELECTRONINTEGRAL_H

