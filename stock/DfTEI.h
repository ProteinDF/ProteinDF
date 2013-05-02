#ifndef DFTEI_H
#define DFTEI_H

#ifdef USE_OLD_TEI_ENGINE

#include <vector>
#include <string>
#include <iostream>

#include "TlPosition.h"

class DfTEI {
public:
    struct PrimitiveShellPair {
public:
        PrimitiveShellPair(double z = 0.0, double k = 0.0,
                           const TlPosition& pos = TlPosition(0.0, 0.0, 0.0))
                : zeta(z), dKab(k), P(pos) {
        }

public:
        double zeta;
        double dKab;
        TlPosition P;
    };


    typedef std::vector<PrimitiveShellPair> PrimitiveShellPairList;


    /// PrimitiveShellPair ソート用関数オブジェクト
    ///
    /// PrimitiveShellPair.dKab の絶対値の大きい順にソート
    struct sort_primitiveShellPair_cmp {
        bool operator()(const PrimitiveShellPair& a, const PrimitiveShellPair& b) {
            return (std::fabs(a.dKab) > std::fabs(b.dKab));
        }
    };

    struct ShellPair {
public:
        TlPosition A;
        TlPosition AB;
        PrimitiveShellPairList pList;
        double dSchwartzValue; // max(sqrt(abs(ij|ij)))の値

public:
        void swap() {
            A = A - AB;
            AB *= -1.0;
        }
    };


    // ERI の型
    // 関数テーブルなどに利用する
    //
    // s: 0, p: 1, d:2, f:3
    // p_shell = 4^3(=64), q_shell = 4^2(=16), r_shell = 4^1(=4), s_shell = 4^0
    enum {
        ERI_SSSS =   0, ERI_SSSP, ERI_SSSD, ERI_SSSF, // 0000
        ERI_SSPS =   4, ERI_SSPP, ERI_SSPD, ERI_SSPF, // 0010
        ERI_SSDS =   8, ERI_SSDP, ERI_SSDD, ERI_SSDF, // 0020
        ERI_SSFS =  12, ERI_SSFP, ERI_SSFD, ERI_SSFF, // 0030
        ERI_SPSS =  16, ERI_SPSP, ERI_SPSD, ERI_SPSF, // 0100
        ERI_SPPS =  20, ERI_SPPP, ERI_SPPD, ERI_SPPF, // 0110
        ERI_SPDS =  24, ERI_SPDP, ERI_SPDD, ERI_SPDF,
        ERI_SPFS =  28, ERI_SPFP, ERI_SPFD, ERI_SPFF,
        ERI_SDSS =  32, ERI_SDSP, ERI_SDSD, ERI_SDSF, // 0200
        ERI_SDPS =  36, ERI_SDPP, ERI_SDPD, ERI_SDPF,
        ERI_SDDS =  40, ERI_SDDP, ERI_SDDD, ERI_SDDF,
        ERI_SDFS =  44, ERI_SDFP, ERI_SDFD, ERI_SDFF,
        ERI_SFSS =  48, ERI_SFSP, ERI_SFSD, ERI_SFSF, // 0300
        ERI_SFPS =  52, ERI_SFPP, ERI_SFPD, ERI_SFPF,
        ERI_SFDS =  56, ERI_SFDP, ERI_SFDD, ERI_SFDF,
        ERI_SFFS =  60, ERI_SFFP, ERI_SFFD, ERI_SFFF,
        //
        ERI_PSSS =  64, ERI_PSSP, ERI_PSSD, ERI_PSSF, // 1000
        ERI_PSPS =  68, ERI_PSPP, ERI_PSPD, ERI_PSPF, // 1010
        ERI_PSDS =  72, ERI_PSDP, ERI_PSDD, ERI_PSDF, // 1020
        ERI_PSFS =  76, ERI_PSFP, ERI_PSFD, ERI_PSFF, // 1030
        ERI_PPSS =  80, ERI_PPSP, ERI_PPSD, ERI_PPSF, // 1100
        ERI_PPPS =  84, ERI_PPPP, ERI_PPPD, ERI_PPPF, // 1110
        ERI_PPDS =  88, ERI_PPDP, ERI_PPDD, ERI_PPDF,
        ERI_PPFS =  92, ERI_PPFP, ERI_PPFD, ERI_PPFF,
        ERI_PDSS =  96, ERI_PDSP, ERI_PDSD, ERI_PDSF, // 1200
        ERI_PDPS = 100, ERI_PDPP, ERI_PDPD, ERI_PDPF,
        ERI_PDDS = 104, ERI_PDDP, ERI_PDDD, ERI_PDDF,
        ERI_PDFS = 108, ERI_PDFP, ERI_PDFD, ERI_PDFF,
        ERI_PFSS = 112, ERI_PFSP, ERI_PFSD, ERI_PFSF, // 1300
        ERI_PFPS = 116, ERI_PFPP, ERI_PFPD, ERI_PFPF,
        ERI_PFDS = 120, ERI_PFDP, ERI_PFDD, ERI_PFDF,
        ERI_PFFS = 124, ERI_PFFP, ERI_PFFD, ERI_PFFF,
        //
        ERI_DSSS = 128, ERI_DSSP, ERI_DSSD, ERI_DSSF, // 2000
        ERI_DSPS = 132, ERI_DSPP, ERI_DSPD, ERI_DSPF, // 2010
        ERI_DSDS = 136, ERI_DSDP, ERI_DSDD, ERI_DSDF, // 2020
        ERI_DSFS = 140, ERI_DSFP, ERI_DSFD, ERI_DSFF, // 2030
        ERI_DPSS = 144, ERI_DPSP, ERI_DPSD, ERI_DPSF, // 2100
        ERI_DPPS = 148, ERI_DPPP, ERI_DPPD, ERI_DPPF, // 2110
        ERI_DPDS = 152, ERI_DPDP, ERI_DPDD, ERI_DPDF,
        ERI_DPFS = 156, ERI_DPFP, ERI_DPFD, ERI_DPFF,
        ERI_DDSS = 160, ERI_DDSP, ERI_DDSD, ERI_DDSF, // 2200
        ERI_DDPS = 164, ERI_DDPP, ERI_DDPD, ERI_DDPF,
        ERI_DDDS = 168, ERI_DDDP, ERI_DDDD, ERI_DDDF,
        ERI_DDFS = 172, ERI_DDFP, ERI_DDFD, ERI_DDFF,
        ERI_DFSS = 176, ERI_DFSP, ERI_DFSD, ERI_DFSF, // 2300
        ERI_DFPS = 180, ERI_DFPP, ERI_DFPD, ERI_DFPF,
        ERI_DFDS = 184, ERI_DFDP, ERI_DFDD, ERI_DFDF,
        ERI_DFFS = 188, ERI_DFFP, ERI_DFFD, ERI_DFFF,
        //
        //
        ERI_FSSS = 192, ERI_FSSP, ERI_FSSD, ERI_FSSF, // 3000
        ERI_FSPS = 196, ERI_FSPP, ERI_FSPD, ERI_FSPF, // 3010
        ERI_FSDS = 200, ERI_FSDP, ERI_FSDD, ERI_FSDF, // 3020
        ERI_FSFS = 204, ERI_FSFP, ERI_FSFD, ERI_FSFF, // 3030
        ERI_FPSS = 208, ERI_FPSP, ERI_FPSD, ERI_FPSF, // 3100
        ERI_FPPS = 212, ERI_FPPP, ERI_FPPD, ERI_FPPF, // 3110
        ERI_FPDS = 216, ERI_FPDP, ERI_FPDD, ERI_FPDF,
        ERI_FPFS = 220, ERI_FPFP, ERI_FPFD, ERI_FPFF,
        ERI_FDSS = 224, ERI_FDSP, ERI_FDSD, ERI_FDSF, // 3200
        ERI_FDPS = 228, ERI_FDPP, ERI_FDPD, ERI_FDPF,
        ERI_FDDS = 232, ERI_FDDP, ERI_FDDD, ERI_FDDF,
        ERI_FDFS = 236, ERI_FDFP, ERI_FDFD, ERI_FDFF,
        ERI_FFSS = 240, ERI_FFSP, ERI_FFSD, ERI_FFSF, // 3300
        ERI_FFPS = 244, ERI_FFPP, ERI_FFPD, ERI_FFPF,
        ERI_FFDS = 248, ERI_FFDP, ERI_FFDD, ERI_FFDF,
        ERI_FFFS = 252, ERI_FFFP, ERI_FFFD, ERI_FFFF,
        //
        ERI_MAXSHELLS
    };


public:
    DfTEI(double value = 1.0E-20);
    ~DfTEI();

    void setCutoffValue(const double value) {
        this->primitiveIntegralsCutoffThreshold_ = value;
    }
    void calc(unsigned int eriType, const ShellPair& IJ, const ShellPair& KL);

    void debug() {
        std::cerr << "DfTEI: " << primitiveIntegralsCutoffThreshold_ << std::endl;
    }

    // type(0:s, 1:p, 2:d)
    static unsigned int getEriType(const int itype, const int jtype,
                                   const int ktype, const int ltype) {
        return ((((itype * 4) + jtype) * 4 + ktype) * 4 + ltype);
    }

private:
    void contractSSSS(const ShellPair& IJ, const ShellPair& KL);
    void contractSSSP(const ShellPair& IJ, const ShellPair& KL);
    void contractSSSD(const ShellPair& IJ, const ShellPair& KL);

    void contractSSPS(const ShellPair& IJ, const ShellPair& KL);
    void contractSSPP(const ShellPair& IJ, const ShellPair& KL);
    void contractSSPD(const ShellPair& IJ, const ShellPair& KL);

    void contractSSDS(const ShellPair& IJ, const ShellPair& KL);
    void contractSSDP(const ShellPair& IJ, const ShellPair& KL);
    void contractSSDD(const ShellPair& IJ, const ShellPair& KL);

    void contractSPSS(const ShellPair& IJ, const ShellPair& KL);
    void contractSPSP(const ShellPair& IJ, const ShellPair& KL);
    void contractSPSD(const ShellPair& IJ, const ShellPair& KL);

    void contractSPPS(const ShellPair& IJ, const ShellPair& KL);
    void contractSPPP(const ShellPair& IJ, const ShellPair& KL);
    void contractSPPD(const ShellPair& IJ, const ShellPair& KL);

    void contractSPDS(const ShellPair& IJ, const ShellPair& KL);
    void contractSPDP(const ShellPair& IJ, const ShellPair& KL);
    void contractSPDD(const ShellPair& IJ, const ShellPair& KL);

    void contractSDSS(const ShellPair& IJ, const ShellPair& KL);
    void contractSDSP(const ShellPair& IJ, const ShellPair& KL);
    void contractSDSD(const ShellPair& IJ, const ShellPair& KL);

    void contractSDPS(const ShellPair& IJ, const ShellPair& KL);
    void contractSDPP(const ShellPair& IJ, const ShellPair& KL);
    void contractSDPD(const ShellPair& IJ, const ShellPair& KL);

    void contractSDDS(const ShellPair& IJ, const ShellPair& KL);
    void contractSDDP(const ShellPair& IJ, const ShellPair& KL);
    void contractSDDD(const ShellPair& IJ, const ShellPair& KL);

    void contractPSSS(const ShellPair& IJ, const ShellPair& KL);
    void contractPSSP(const ShellPair& IJ, const ShellPair& KL);
    void contractPSSD(const ShellPair& IJ, const ShellPair& KL);

    void contractPSPS(const ShellPair& IJ, const ShellPair& KL);
    void contractPSPP(const ShellPair& IJ, const ShellPair& KL);
    void contractPSPD(const ShellPair& IJ, const ShellPair& KL);

    void contractPSDS(const ShellPair& IJ, const ShellPair& KL);
    void contractPSDP(const ShellPair& IJ, const ShellPair& KL);
    void contractPSDD(const ShellPair& IJ, const ShellPair& KL);

    void contractPPSS(const ShellPair& IJ, const ShellPair& KL);
    void contractPPSP(const ShellPair& IJ, const ShellPair& KL);
    void contractPPSD(const ShellPair& IJ, const ShellPair& KL);

    void contractPPPS(const ShellPair& IJ, const ShellPair& KL);
    void contractPPPP(const ShellPair& IJ, const ShellPair& KL);
    void contractPPPD(const ShellPair& IJ, const ShellPair& KL);

    void contractPPDS(const ShellPair& IJ, const ShellPair& KL);
    void contractPPDP(const ShellPair& IJ, const ShellPair& KL);
    void contractPPDD(const ShellPair& IJ, const ShellPair& KL);

    void contractPDSS(const ShellPair& IJ, const ShellPair& KL);
    void contractPDSP(const ShellPair& IJ, const ShellPair& KL);
    void contractPDSD(const ShellPair& IJ, const ShellPair& KL);

    void contractPDPS(const ShellPair& IJ, const ShellPair& KL);
    void contractPDPP(const ShellPair& IJ, const ShellPair& KL);
    void contractPDPD(const ShellPair& IJ, const ShellPair& KL);

    void contractPDDS(const ShellPair& IJ, const ShellPair& KL);
    void contractPDDP(const ShellPair& IJ, const ShellPair& KL);
    void contractPDDD(const ShellPair& IJ, const ShellPair& KL);

    void contractDSSS(const ShellPair& IJ, const ShellPair& KL);
    void contractDSSP(const ShellPair& IJ, const ShellPair& KL);
    void contractDSSD(const ShellPair& IJ, const ShellPair& KL);

    void contractDSPS(const ShellPair& IJ, const ShellPair& KL);
    void contractDSPP(const ShellPair& IJ, const ShellPair& KL);
    void contractDSPD(const ShellPair& IJ, const ShellPair& KL);

    void contractDSDS(const ShellPair& IJ, const ShellPair& KL);
    void contractDSDP(const ShellPair& IJ, const ShellPair& KL);
    void contractDSDD(const ShellPair& IJ, const ShellPair& KL);

    void contractDPSS(const ShellPair& IJ, const ShellPair& KL);
    void contractDPSP(const ShellPair& IJ, const ShellPair& KL);
    void contractDPSD(const ShellPair& IJ, const ShellPair& KL);

    void contractDPPS(const ShellPair& IJ, const ShellPair& KL);
    void contractDPPP(const ShellPair& IJ, const ShellPair& KL);
    void contractDPPD(const ShellPair& IJ, const ShellPair& KL);

    void contractDPDS(const ShellPair& IJ, const ShellPair& KL);
    void contractDPDP(const ShellPair& IJ, const ShellPair& KL);
    void contractDPDD(const ShellPair& IJ, const ShellPair& KL);

    void contractDDSS(const ShellPair& IJ, const ShellPair& KL);
    void contractDDSP(const ShellPair& IJ, const ShellPair& KL);
    void contractDDSD(const ShellPair& IJ, const ShellPair& KL);

    void contractDDPS(const ShellPair& IJ, const ShellPair& KL);
    void contractDDPP(const ShellPair& IJ, const ShellPair& KL);
    void contractDDPD(const ShellPair& IJ, const ShellPair& KL);

    void contractDDDS(const ShellPair& IJ, const ShellPair& KL);
    void contractDDDP(const ShellPair& IJ, const ShellPair& KL);
    void contractDDDD(const ShellPair& IJ, const ShellPair& KL);

private:
    /// m x n の要素を転置する
    void transform(const int nBeforeNumOfColumns,
                   const int nBeforeNumOfRows,
                   double* pMatrix);

    ShellPair swapShellPair(const ShellPair& sp);

    bool isPrimitiveShellsCutoff() const;

private:
    void HRR_DPXX(const double ABx, const double ABy, const double ABz,
                  const double* FSXX, const double* DSXX,
                  const int nMaxK, const int nMaxL, double* DPXX);
    void HRR_FPXX(const double ABx, const double ABy, const double ABz,
                  const double* GSXX, const double* FSXX,
                  const int nMaxK, const int nMaxL, double* FPXX);
    void HRR_DDXX(const double ABx, const double ABy, const double ABz,
                  const double* FPXX, const double* DPXX,
                  const int nMaxK, const int nMaxL, double* DDXX);
    void HRR_XXPP(const double CDx, const double CDy, const double CDz,
                  const double* XXDS, const double* XXPS,
                  const int nMaxI, const int nMaxJ, double* XXPP);
    void HRR_XXDP(const double CDx, const double CDy, const double CDz,
                  const double* XXFS, const double* XXDS,
                  const int nMaxI, const int nMaxJ, double* XXDP);
    void HRR_XXFP(const double CDx, const double CDy, const double CDz,
                  const double* XXGS, const double* XXFS,
                  const int nMaxI, const int nMaxJ, double* XXFP);
    void HRR_XXDD(const double CDx, const double CDy, const double CDz,
                  const double* XXFP, const double* XXDP,
                  const int nMaxI, const int nMaxJ, double* XXDD);

private:
    /// convert integrals based on 6D to 5D at I-shell
    void convert6Dto5D_I(const int maxJ, const int maxK, const int maxL,
                         const double* pIn, double* pOut);

    /// convert integrals based on 6D to 5D at J-shell
    void convert6Dto5D_J(const int maxI, const int maxK, const int maxL,
                         const double* pIn, double* pOut);

    /// convert integrals based on 6D to 5D at K-shell
    void convert6Dto5D_K(const int maxI, const int maxJ, const int maxL,
                         const double* pIn, double* pOut);

    /// convert integrals based on 6D to 5D at L-shell
    void convert6Dto5D_L(const int maxI, const int maxJ, const int maxK,
                         const double* pIn, double* pOut);

private:
    void primitiveSSSS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveSSPS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveSSDS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitivePSSS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitivePSPS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitivePSDS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitivePSFS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveDSSS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveDSPS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveDSDS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveDSFS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveDSGS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveFSSS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveGSSS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveFSPS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveGSPS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveFSDS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveGSDS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveFSFS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveFSGS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveGSFS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);
    void primitiveGSGS(const PrimitiveShellPair& ij, const PrimitiveShellPair& kl,
                       const int nEndM, const int nStartM =0);

public:
    // contract ERIを格納する配列
    double ERI[625]; // 5*5*5*5

private:
    static const double INV_SQRT3;

    /// 関数テーブルに使う、型の宣言
    typedef void (DfTEI::*ContractXXXX)(const ShellPair&, const ShellPair&);

    // contractGTO積分用関数テーブル
    //
    // インデックスには ERI_TYPE を利用する。
    // 即ち、(*m_contractXXXX[ERI_TYPE])で、ERI_TYPEに該当するメンバ関数が呼ばれる。
    ContractXXXX contractXXXX_[256];

    /// 積分種類の文字列(eg. ERI_SSSS)を格納したテーブル
    ///
    /// 配列の添字にはgetEriで得られた数を用いる
    static const std::string STR_XXXX[];

private:
    double primitiveIntegralsCutoffThreshold_;

    double m_dRho;
    double m_divZetaEta;
    TlPosition posP;
    TlPosition WP;
    TlPosition PA;
    TlPosition posQ;
    TlPosition WQ;
    TlPosition QC;
    TlPosition posW;

    double pFmTBuf[9]; // FmT計算用バッファ

    double pSSSS[ 9];      //
    double pSSPS[ 4][  3]; //  3
    double pSSDS[ 3][  6]; //  3

    double pPSSS[ 8][  3]; //  3
    double pPSPS[ 4][  9]; //  3* 3
    double pPSDS[ 3][ 18]; //  3* 6
    double pPSFS[ 2][ 30]; //  3* 10

    double pDSSS[ 7][  6]; //  6
    double pDSPS[ 4][ 18]; //  6* 3
    double pDSDS[ 3][ 36]; //  6* 6
    double pDSFS[ 2][ 60]; //  6* 10
    double pDSGS[ 1][ 90]; //  6* 15

    double pFSSS[ 6][ 10]; // 10
    double pFSPS[ 4][ 30]; // 10* 3
    double pFSDS[ 3][ 60]; // 10* 6
    double pFSFS[ 2][100]; // 10*10
    double pFSGS[ 1][150]; // 10*15

    double pGSSS[ 5][ 15]; // 15
    double pGSPS[ 4][ 45]; // 15* 3

    double pGSDS[ 3][ 90]; // 15* 6
    double pGSFS[ 2][150]; // 15*10
    double pGSGS[ 1][225]; // 15*15
};

#endif // USE_OLD_TEI_ENGINE

#endif // DFTEI_H
