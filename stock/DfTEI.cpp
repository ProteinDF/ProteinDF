#ifdef USE_OLD_TEI_ENGINE

#include "DfTEI.h"
#include "TlFmt.h"

const double DfTEI::INV_SQRT3 = 1.0 / sqrt(3.0);
// const double DfTwoElectronIntegral::CK          = sqrt(2.0 * sqrt(M_PI)) * M_PI;
// const int DfTwoElectronIntegral::MAX_SHELL_TYPE = 3;
// const double DfTwoElectronIntegral::CONTRIBUTE_COEF = 2.0 * std::pow(M_PI, 2.5);

const std::string DfTEI::STR_XXXX[] = {
    "ERI_SSSS", "ERI_SSSP", "ERI_SSSD", "ERI_SSSF", // 0000
    "ERI_SSPS", "ERI_SSPP", "ERI_SSPD", "ERI_SSPF", // 0010
    "ERI_SSDS", "ERI_SSDP", "ERI_SSDD", "ERI_SSDF", // 0020
    "ERI_SSFS", "ERI_SSFP", "ERI_SSFD", "ERI_SSFF", // 0030
    "ERI_SPSS", "ERI_SPSP", "ERI_SPSD", "ERI_SPSF", // 0100
    "ERI_SPPS", "ERI_SPPP", "ERI_SPPD", "ERI_SPPF", // 0110
    "ERI_SPDS", "ERI_SPDP", "ERI_SPDD", "ERI_SPDF",
    "ERI_SPFS", "ERI_SPFP", "ERI_SPFD", "ERI_SPFF",
    "ERI_SDSS", "ERI_SDSP", "ERI_SDSD", "ERI_SDSF", // 0200
    "ERI_SDPS", "ERI_SDPP", "ERI_SDPD", "ERI_SDPF",
    "ERI_SDDS", "ERI_SDDP", "ERI_SDDD", "ERI_SDDF",
    "ERI_SDFS", "ERI_SDFP", "ERI_SDFD", "ERI_SDFF",
    "ERI_SFSS", "ERI_SFSP", "ERI_SFSD", "ERI_SFSF", // 0300
    "ERI_SFPS", "ERI_SFPP", "ERI_SFPD", "ERI_SFPF",
    "ERI_SFDS", "ERI_SFDP", "ERI_SFDD", "ERI_SFDF",
    "ERI_SFFS", "ERI_SFFP", "ERI_SFFD", "ERI_SFFF",
    //
    "ERI_PSSS", "ERI_PSSP", "ERI_PSSD", "ERI_PSSF", // 1000
    "ERI_PSPS", "ERI_PSPP", "ERI_PSPD", "ERI_PSPF", // 1010
    "ERI_PSDS", "ERI_PSDP", "ERI_PSDD", "ERI_PSDF", // 1020
    "ERI_PSFS", "ERI_PSFP", "ERI_PSFD", "ERI_PSFF", // 1030
    "ERI_PPSS", "ERI_PPSP", "ERI_PPSD", "ERI_PPSF", // 1100
    "ERI_PPPS", "ERI_PPPP", "ERI_PPPD", "ERI_PPPF", // 1110
    "ERI_PPDS", "ERI_PPDP", "ERI_PPDD", "ERI_PPDF",
    "ERI_PPFS", "ERI_PPFP", "ERI_PPFD", "ERI_PPFF",
    "ERI_PDSS", "ERI_PDSP", "ERI_PDSD", "ERI_PDSF", // 1200
    "ERI_PDPS", "ERI_PDPP", "ERI_PDPD", "ERI_PDPF",
    "ERI_PDDS", "ERI_PDDP", "ERI_PDDD", "ERI_PDDF",
    "ERI_PDFS", "ERI_PDFP", "ERI_PDFD", "ERI_PDFF",
    "ERI_PFSS", "ERI_PFSP", "ERI_PFSD", "ERI_PFSF", // 1300
    "ERI_PFPS", "ERI_PFPP", "ERI_PFPD", "ERI_PFPF",
    "ERI_PFDS", "ERI_PFDP", "ERI_PFDD", "ERI_PFDF",
    "ERI_PFFS", "ERI_PFFP", "ERI_PFFD", "ERI_PFFF",
    //
    "ERI_DSSS", "ERI_DSSP", "ERI_DSSD", "ERI_DSSF", // 2000
    "ERI_DSPS", "ERI_DSPP", "ERI_DSPD", "ERI_DSPF", // 2010
    "ERI_DSDS", "ERI_DSDP", "ERI_DSDD", "ERI_DSDF", // 2020
    "ERI_DSFS", "ERI_DSFP", "ERI_DSFD", "ERI_DSFF", // 2030
    "ERI_DPSS", "ERI_DPSP", "ERI_DPSD", "ERI_DPSF", // 2100
    "ERI_DPPS", "ERI_DPPP", "ERI_DPPD", "ERI_DPPF", // 2110
    "ERI_DPDS", "ERI_DPDP", "ERI_DPDD", "ERI_DPDF",
    "ERI_DPFS", "ERI_DPFP", "ERI_DPFD", "ERI_DPFF",
    "ERI_DDSS", "ERI_DDSP", "ERI_DDSD", "ERI_DDSF", // 2200
    "ERI_DDPS", "ERI_DDPP", "ERI_DDPD", "ERI_DDPF",
    "ERI_DDDS", "ERI_DDDP", "ERI_DDDD", "ERI_DDDF",
    "ERI_DDFS", "ERI_DDFP", "ERI_DDFD", "ERI_DDFF",
    "ERI_DFSS", "ERI_DFSP", "ERI_DFSD", "ERI_DFSF", // 2300
    "ERI_DFPS", "ERI_DFPP", "ERI_DFPD", "ERI_DFPF",
    "ERI_DFDS", "ERI_DFDP", "ERI_DFDD", "ERI_DFDF",
    "ERI_DFFS", "ERI_DFFP", "ERI_DFFD", "ERI_DFFF",
    //
    "ERI_FSSS", "ERI_FSSP", "ERI_FSSD", "ERI_FSSF", // 3000
    "ERI_FSPS", "ERI_FSPP", "ERI_FSPD", "ERI_FSPF", // 3010
    "ERI_FSDS", "ERI_FSDP", "ERI_FSDD", "ERI_FSDF", // 3020
    "ERI_FSFS", "ERI_FSFP", "ERI_FSFD", "ERI_FSFF", // 3030
    "ERI_FPSS", "ERI_FPSP", "ERI_FPSD", "ERI_FPSF", // 3100
    "ERI_FPPS", "ERI_FPPP", "ERI_FPPD", "ERI_FPPF", // 3110
    "ERI_FPDS", "ERI_FPDP", "ERI_FPDD", "ERI_FPDF",
    "ERI_FPFS", "ERI_FPFP", "ERI_FPFD", "ERI_FPFF",
    "ERI_FDSS", "ERI_FDSP", "ERI_FDSD", "ERI_FDSF", // 3200
    "ERI_FDPS", "ERI_FDPP", "ERI_FDPD", "ERI_FDPF",
    "ERI_FDDS", "ERI_FDDP", "ERI_FDDD", "ERI_FDDF",
    "ERI_FDFS", "ERI_FDFP", "ERI_FDFD", "ERI_FDFF",
    "ERI_FFSS", "ERI_FFSP", "ERI_FFSD", "ERI_FFSF", // 3300
    "ERI_FFPS", "ERI_FFPP", "ERI_FFPD", "ERI_FFPF",
    "ERI_FFDS", "ERI_FFDP", "ERI_FFDD", "ERI_FFDF",
    "ERI_FFFS", "ERI_FFFP", "ERI_FFFD", "ERI_FFFF"
};


DfTEI::DfTEI(double value) : primitiveIntegralsCutoffThreshold_(value)
{
    // 関数テーブルの作成
    this->contractXXXX_[ERI_SSSS] = &DfTEI::contractSSSS;
    this->contractXXXX_[ERI_SSSP] = &DfTEI::contractSSSP;
    this->contractXXXX_[ERI_SSSD] = &DfTEI::contractSSSD;
    this->contractXXXX_[ERI_SSPS] = &DfTEI::contractSSPS;
    this->contractXXXX_[ERI_SSPP] = &DfTEI::contractSSPP;
    this->contractXXXX_[ERI_SSPD] = &DfTEI::contractSSPD;
    this->contractXXXX_[ERI_SSDS] = &DfTEI::contractSSDS;
    this->contractXXXX_[ERI_SSDP] = &DfTEI::contractSSDP;
    this->contractXXXX_[ERI_SSDD] = &DfTEI::contractSSDD;

    this->contractXXXX_[ERI_SPSS] = &DfTEI::contractSPSS;
    this->contractXXXX_[ERI_SPSP] = &DfTEI::contractSPSP;
    this->contractXXXX_[ERI_SPSD] = &DfTEI::contractSPSD;
    this->contractXXXX_[ERI_SPPS] = &DfTEI::contractSPPS;
    this->contractXXXX_[ERI_SPPP] = &DfTEI::contractSPPP;
    this->contractXXXX_[ERI_SPPD] = &DfTEI::contractSPPD;
    this->contractXXXX_[ERI_SPDS] = &DfTEI::contractSPDS;
    this->contractXXXX_[ERI_SPDP] = &DfTEI::contractSPDP;
    this->contractXXXX_[ERI_SPDD] = &DfTEI::contractSPDD;

    this->contractXXXX_[ERI_SDSS] = &DfTEI::contractSDSS;
    this->contractXXXX_[ERI_SDSP] = &DfTEI::contractSDSP;
    this->contractXXXX_[ERI_SDSD] = &DfTEI::contractSDSD;
    this->contractXXXX_[ERI_SDPS] = &DfTEI::contractSDPS;
    this->contractXXXX_[ERI_SDPP] = &DfTEI::contractSDPP;
    this->contractXXXX_[ERI_SDPD] = &DfTEI::contractSDPD;
    this->contractXXXX_[ERI_SDDS] = &DfTEI::contractSDDS;
    this->contractXXXX_[ERI_SDDP] = &DfTEI::contractSDDP;
    this->contractXXXX_[ERI_SDDD] = &DfTEI::contractSDDD;

    this->contractXXXX_[ERI_PSSS] = &DfTEI::contractPSSS;
    this->contractXXXX_[ERI_PSSP] = &DfTEI::contractPSSP;
    this->contractXXXX_[ERI_PSSD] = &DfTEI::contractPSSD;
    this->contractXXXX_[ERI_PSPS] = &DfTEI::contractPSPS;
    this->contractXXXX_[ERI_PSPP] = &DfTEI::contractPSPP;
    this->contractXXXX_[ERI_PSPD] = &DfTEI::contractPSPD;
    this->contractXXXX_[ERI_PSDS] = &DfTEI::contractPSDS;
    this->contractXXXX_[ERI_PSDP] = &DfTEI::contractPSDP;
    this->contractXXXX_[ERI_PSDD] = &DfTEI::contractPSDD;

    this->contractXXXX_[ERI_PPSS] = &DfTEI::contractPPSS;
    this->contractXXXX_[ERI_PPSP] = &DfTEI::contractPPSP;
    this->contractXXXX_[ERI_PPSD] = &DfTEI::contractPPSD;
    this->contractXXXX_[ERI_PPPS] = &DfTEI::contractPPPS;
    this->contractXXXX_[ERI_PPPP] = &DfTEI::contractPPPP;
    this->contractXXXX_[ERI_PPPD] = &DfTEI::contractPPPD;
    this->contractXXXX_[ERI_PPDS] = &DfTEI::contractPPDS;
    this->contractXXXX_[ERI_PPDP] = &DfTEI::contractPPDP;
    this->contractXXXX_[ERI_PPDD] = &DfTEI::contractPPDD;

    this->contractXXXX_[ERI_PDSS] = &DfTEI::contractPDSS;
    this->contractXXXX_[ERI_PDSP] = &DfTEI::contractPDSP;
    this->contractXXXX_[ERI_PDSD] = &DfTEI::contractPDSD;
    this->contractXXXX_[ERI_PDPS] = &DfTEI::contractPDPS;
    this->contractXXXX_[ERI_PDPP] = &DfTEI::contractPDPP;
    this->contractXXXX_[ERI_PDPD] = &DfTEI::contractPDPD;
    this->contractXXXX_[ERI_PDDS] = &DfTEI::contractPDDS;
    this->contractXXXX_[ERI_PDDP] = &DfTEI::contractPDDP;
    this->contractXXXX_[ERI_PDDD] = &DfTEI::contractPDDD;

    this->contractXXXX_[ERI_DSSS] = &DfTEI::contractDSSS;
    this->contractXXXX_[ERI_DSSP] = &DfTEI::contractDSSP;
    this->contractXXXX_[ERI_DSSD] = &DfTEI::contractDSSD;
    this->contractXXXX_[ERI_DSPS] = &DfTEI::contractDSPS;
    this->contractXXXX_[ERI_DSPP] = &DfTEI::contractDSPP;
    this->contractXXXX_[ERI_DSPD] = &DfTEI::contractDSPD;
    this->contractXXXX_[ERI_DSDS] = &DfTEI::contractDSDS;
    this->contractXXXX_[ERI_DSDP] = &DfTEI::contractDSDP;
    this->contractXXXX_[ERI_DSDD] = &DfTEI::contractDSDD;

    this->contractXXXX_[ERI_DPSS] = &DfTEI::contractDPSS;
    this->contractXXXX_[ERI_DPSP] = &DfTEI::contractDPSP;
    this->contractXXXX_[ERI_DPSD] = &DfTEI::contractDPSD;
    this->contractXXXX_[ERI_DPPS] = &DfTEI::contractDPPS;
    this->contractXXXX_[ERI_DPPP] = &DfTEI::contractDPPP;
    this->contractXXXX_[ERI_DPPD] = &DfTEI::contractDPPD;
    this->contractXXXX_[ERI_DPDS] = &DfTEI::contractDPDS;
    this->contractXXXX_[ERI_DPDP] = &DfTEI::contractDPDP;
    this->contractXXXX_[ERI_DPDD] = &DfTEI::contractDPDD;

    this->contractXXXX_[ERI_DDSS] = &DfTEI::contractDDSS;
    this->contractXXXX_[ERI_DDSP] = &DfTEI::contractDDSP;
    this->contractXXXX_[ERI_DDSD] = &DfTEI::contractDDSD;
    this->contractXXXX_[ERI_DDPS] = &DfTEI::contractDDPS;
    this->contractXXXX_[ERI_DDPP] = &DfTEI::contractDDPP;
    this->contractXXXX_[ERI_DDPD] = &DfTEI::contractDDPD;
    this->contractXXXX_[ERI_DDDS] = &DfTEI::contractDDDS;
    this->contractXXXX_[ERI_DDDP] = &DfTEI::contractDDDP;
    this->contractXXXX_[ERI_DDDD] = &DfTEI::contractDDDD;
}


DfTEI::~DfTEI()
{
}


void DfTEI::calc(unsigned int eriType, const ShellPair& IJ, const ShellPair& KL)
{
    (this->*contractXXXX_[eriType])(IJ, KL);
}

void DfTEI::transform(const int nBeforeNumOfColumns,
                      const int nBeforeNumOfRows, double* pMatrix)
{
    const int nNumOfElements = nBeforeNumOfColumns * nBeforeNumOfRows;
    double* pBuf = new double[nNumOfElements];

    for (int i = 0; i < nNumOfElements; ++i) {
        pBuf[i] = pMatrix[i];
    }

    for (int m = 0; m < nBeforeNumOfColumns; ++m) {
        for (int n = 0; n < nBeforeNumOfRows; ++n) {
            pMatrix[n*nBeforeNumOfColumns +m] = pBuf[m*nBeforeNumOfRows +n];
        }
    }

    delete[] pBuf;
    pBuf = NULL;
}


DfTEI::ShellPair DfTEI::swapShellPair(const ShellPair& sp)
{
    ShellPair ans = sp;
    ans.swap();

    return ans;
}


// HRR DPXX
void DfTEI::HRR_DPXX(const double ABx, const double ABy, const double ABz,
                     const double* FSXX, const double* DSXX,
                     const int nMaxK, const int nMaxL, double* DPXX)
{
    // 0:xxx 1:xxy 2:xxz 3:xyy 4:xyz 5:xzz
    // 6:yyy 7:yyz 8:yzz 9:zzz
    const int nJ = nMaxK * nMaxL;
    for (int k = 0; k < nMaxK; ++k) {
        const int nK = nMaxL*k;
        for (int l = 0; l < nMaxL; ++l) {
            DPXX[3*nJ*0 +nJ*0 +nK +l] = FSXX[nJ* 0 +nK +l] +ABx*DSXX[nJ* 0 +nK +l]; // i=x a=xx
            DPXX[3*nJ*0 +nJ*1 +nK +l] = FSXX[nJ* 1 +nK +l] +ABy*DSXX[nJ* 0 +nK +l]; // i=y
            DPXX[3*nJ*0 +nJ*2 +nK +l] = FSXX[nJ* 2 +nK +l] +ABz*DSXX[nJ* 0 +nK +l]; // i=z

            DPXX[3*nJ*1 +nJ*0 +nK +l] = FSXX[nJ* 1 +nK +l] +ABx*DSXX[nJ* 1 +nK +l]; // i=x a=xy
            DPXX[3*nJ*1 +nJ*1 +nK +l] = FSXX[nJ* 3 +nK +l] +ABy*DSXX[nJ* 1 +nK +l]; // i=y
            DPXX[3*nJ*1 +nJ*2 +nK +l] = FSXX[nJ* 4 +nK +l] +ABz*DSXX[nJ* 1 +nK +l]; // i=z

            DPXX[3*nJ*2 +nJ*0 +nK +l] = FSXX[nJ* 2 +nK +l] +ABx*DSXX[nJ* 2 +nK +l]; // i=x a=xz
            DPXX[3*nJ*2 +nJ*1 +nK +l] = FSXX[nJ* 4 +nK +l] +ABy*DSXX[nJ* 2 +nK +l]; // i=y
            DPXX[3*nJ*2 +nJ*2 +nK +l] = FSXX[nJ* 5 +nK +l] +ABz*DSXX[nJ* 2 +nK +l]; // i=z

            DPXX[3*nJ*3 +nJ*0 +nK +l] = FSXX[nJ* 3 +nK +l] +ABx*DSXX[nJ* 3 +nK +l]; // i=x a=yy
            DPXX[3*nJ*3 +nJ*1 +nK +l] = FSXX[nJ* 6 +nK +l] +ABy*DSXX[nJ* 3 +nK +l]; // i=y
            DPXX[3*nJ*3 +nJ*2 +nK +l] = FSXX[nJ* 7 +nK +l] +ABz*DSXX[nJ* 3 +nK +l]; // i=z

            DPXX[3*nJ*4 +nJ*0 +nK +l] = FSXX[nJ* 4 +nK +l] +ABx*DSXX[nJ* 4 +nK +l]; // i=x a=yz
            DPXX[3*nJ*4 +nJ*1 +nK +l] = FSXX[nJ* 7 +nK +l] +ABy*DSXX[nJ* 4 +nK +l]; // i=y
            DPXX[3*nJ*4 +nJ*2 +nK +l] = FSXX[nJ* 8 +nK +l] +ABz*DSXX[nJ* 4 +nK +l]; // i=z

            DPXX[3*nJ*5 +nJ*0 +nK +l] = FSXX[nJ* 5 +nK +l] +ABx*DSXX[nJ* 5 +nK +l]; // i=x a=zz
            DPXX[3*nJ*5 +nJ*1 +nK +l] = FSXX[nJ* 8 +nK +l] +ABy*DSXX[nJ* 5 +nK +l]; // i=y
            DPXX[3*nJ*5 +nJ*2 +nK +l] = FSXX[nJ* 9 +nK +l] +ABz*DSXX[nJ* 5 +nK +l]; // i=z
        }
    }
}


// HRR FPXX
void DfTEI::HRR_FPXX(const double ABx, const double ABy, const double ABz,
                     const double* GSXX, const double* FSXX,
                     const int nMaxK, const int nMaxL, double* FPXX)
{
    // 0:xxxx  1:xxxy  2:xxxz  3:xxyy  4:xxyz
    // 5:xxzz  6:xyyy  7:xyyz  8:xyzz  9:xzzz
    //10:yyyy 11:yyyz 12:yyzz 13:yzzz 14:zzzz
    const int nJ = nMaxK * nMaxL;
    for (int k = 0; k < nMaxK; ++k) {
        const int nK = nMaxL*k;
        for (int l = 0; l < nMaxL; ++l) {
            FPXX[3*nJ*0 +nJ*0 +nK +l] = GSXX[nJ* 0 +nK +l] + ABx*FSXX[nJ*0 +nK +l]; // i=x a=xxx
            FPXX[3*nJ*0 +nJ*1 +nK +l] = GSXX[nJ* 1 +nK +l] + ABy*FSXX[nJ*0 +nK +l]; // i=y
            FPXX[3*nJ*0 +nJ*2 +nK +l] = GSXX[nJ* 2 +nK +l] + ABz*FSXX[nJ*0 +nK +l]; // i=z

            FPXX[3*nJ*1 +nJ*0 +nK +l] = GSXX[nJ* 1 +nK +l] + ABx*FSXX[nJ*1 +nK +l]; // i=x a=xxy
            FPXX[3*nJ*1 +nJ*1 +nK +l] = GSXX[nJ* 3 +nK +l] + ABy*FSXX[nJ*1 +nK +l]; // i=y
            FPXX[3*nJ*1 +nJ*2 +nK +l] = GSXX[nJ* 4 +nK +l] + ABz*FSXX[nJ*1 +nK +l]; // i=z

            FPXX[3*nJ*2 +nJ*0 +nK +l] = GSXX[nJ* 2 +nK +l] + ABx*FSXX[nJ*2 +nK +l]; // i=x a=xxz
            FPXX[3*nJ*2 +nJ*1 +nK +l] = GSXX[nJ* 4 +nK +l] + ABy*FSXX[nJ*2 +nK +l]; // i=y
            FPXX[3*nJ*2 +nJ*2 +nK +l] = GSXX[nJ* 5 +nK +l] + ABz*FSXX[nJ*2 +nK +l]; // i=z

            FPXX[3*nJ*3 +nJ*0 +nK +l] = GSXX[nJ* 3 +nK +l] + ABx*FSXX[nJ*3 +nK +l]; // i=x a=xyy
            FPXX[3*nJ*3 +nJ*1 +nK +l] = GSXX[nJ* 6 +nK +l] + ABy*FSXX[nJ*3 +nK +l]; // i=y
            FPXX[3*nJ*3 +nJ*2 +nK +l] = GSXX[nJ* 7 +nK +l] + ABz*FSXX[nJ*3 +nK +l]; // i=z

            FPXX[3*nJ*4 +nJ*0 +nK +l] = GSXX[nJ* 4 +nK +l] + ABx*FSXX[nJ*4 +nK +l]; // i=x a=xyz
            FPXX[3*nJ*4 +nJ*1 +nK +l] = GSXX[nJ* 7 +nK +l] + ABy*FSXX[nJ*4 +nK +l]; // i=y
            FPXX[3*nJ*4 +nJ*2 +nK +l] = GSXX[nJ* 8 +nK +l] + ABz*FSXX[nJ*4 +nK +l]; // i=z

            FPXX[3*nJ*5 +nJ*0 +nK +l] = GSXX[nJ* 5 +nK +l] + ABx*FSXX[nJ*5 +nK +l]; // i=x a=xzz
            FPXX[3*nJ*5 +nJ*1 +nK +l] = GSXX[nJ* 8 +nK +l] + ABy*FSXX[nJ*5 +nK +l]; // i=y
            FPXX[3*nJ*5 +nJ*2 +nK +l] = GSXX[nJ* 9 +nK +l] + ABz*FSXX[nJ*5 +nK +l]; // i=z

            FPXX[3*nJ*6 +nJ*0 +nK +l] = GSXX[nJ* 6 +nK +l] + ABx*FSXX[nJ*6 +nK +l]; // i=x a=yyy
            FPXX[3*nJ*6 +nJ*1 +nK +l] = GSXX[nJ*10 +nK +l] + ABy*FSXX[nJ*6 +nK +l]; // i=y
            FPXX[3*nJ*6 +nJ*2 +nK +l] = GSXX[nJ*11 +nK +l] + ABz*FSXX[nJ*6 +nK +l]; // i=z

            FPXX[3*nJ*7 +nJ*0 +nK +l] = GSXX[nJ* 7 +nK +l] + ABx*FSXX[nJ*7 +nK +l]; // i=x a=yyz
            FPXX[3*nJ*7 +nJ*1 +nK +l] = GSXX[nJ*11 +nK +l] + ABy*FSXX[nJ*7 +nK +l]; // i=y
            FPXX[3*nJ*7 +nJ*2 +nK +l] = GSXX[nJ*12 +nK +l] + ABz*FSXX[nJ*7 +nK +l]; // i=z

            FPXX[3*nJ*8 +nJ*0 +nK +l] = GSXX[nJ* 8 +nK +l] + ABx*FSXX[nJ*8 +nK +l]; // i=x a=yzz
            FPXX[3*nJ*8 +nJ*1 +nK +l] = GSXX[nJ*12 +nK +l] + ABy*FSXX[nJ*8 +nK +l]; // i=y
            FPXX[3*nJ*8 +nJ*2 +nK +l] = GSXX[nJ*13 +nK +l] + ABz*FSXX[nJ*8 +nK +l]; // i=z

            FPXX[3*nJ*9 +nJ*0 +nK +l] = GSXX[nJ* 9 +nK +l] + ABx*FSXX[nJ*9 +nK +l]; // i=x a=zzz
            FPXX[3*nJ*9 +nJ*1 +nK +l] = GSXX[nJ*13 +nK +l] + ABy*FSXX[nJ*9 +nK +l]; // i=y
            FPXX[3*nJ*9 +nJ*2 +nK +l] = GSXX[nJ*14 +nK +l] + ABz*FSXX[nJ*9 +nK +l]; // i=z
        }
    }
}


// HRR DDXX
void DfTEI::HRR_DDXX(const double ABx, const double ABy, const double ABz,
                     const double* FPXX, const double* DPXX,
                     const int nMaxK, const int nMaxL, double* DDXX)
{
    // 0xxx 1xxy 2xxz 3xyy 4xyz 5xzz
    // 6yyy 7yyz 8yzz 9zzz
    const int nJ = nMaxK * nMaxL;
    for (int k = 0; k < nMaxK; ++k) {
        const int nK = nMaxL*k;
        for (int l = 0; l < nMaxL; ++l) {
            DDXX[6*nJ*0 +nJ*0 +nK +l] = FPXX[3*nJ*0 +nJ*0 +nK +l] +ABx*DPXX[3*nJ*0 +nJ*0 +nK +l]; // i=x a=xx b=x
            DDXX[6*nJ*0 +nJ*1 +nK +l] = FPXX[3*nJ*1 +nJ*0 +nK +l] +ABy*DPXX[3*nJ*0 +nJ*0 +nK +l]; // i=y a=xx b=x
            DDXX[6*nJ*0 +nJ*2 +nK +l] = FPXX[3*nJ*2 +nJ*0 +nK +l] +ABz*DPXX[3*nJ*0 +nJ*0 +nK +l]; // i=z a=xx b=x
            DDXX[6*nJ*0 +nJ*3 +nK +l] = FPXX[3*nJ*1 +nJ*1 +nK +l] +ABy*DPXX[3*nJ*0 +nJ*1 +nK +l]; // i=y a=xx b=y
            DDXX[6*nJ*0 +nJ*4 +nK +l] = FPXX[3*nJ*2 +nJ*1 +nK +l] +ABz*DPXX[3*nJ*0 +nJ*1 +nK +l]; // i=z a=xx b=y
            DDXX[6*nJ*0 +nJ*5 +nK +l] = FPXX[3*nJ*2 +nJ*2 +nK +l] +ABz*DPXX[3*nJ*0 +nJ*2 +nK +l]; // i=z a=xx b=z

            DDXX[6*nJ*1 +nJ*0 +nK +l] = FPXX[3*nJ*1 +nJ*0 +nK +l] +ABx*DPXX[3*nJ*1 +nJ*0 +nK +l]; // i=x a=xy b=x
            DDXX[6*nJ*1 +nJ*1 +nK +l] = FPXX[3*nJ*3 +nJ*0 +nK +l] +ABy*DPXX[3*nJ*1 +nJ*0 +nK +l]; // i=y a=xy b=x
            DDXX[6*nJ*1 +nJ*2 +nK +l] = FPXX[3*nJ*4 +nJ*0 +nK +l] +ABz*DPXX[3*nJ*1 +nJ*0 +nK +l]; // i=z a=xy b=x
            DDXX[6*nJ*1 +nJ*3 +nK +l] = FPXX[3*nJ*3 +nJ*1 +nK +l] +ABy*DPXX[3*nJ*1 +nJ*1 +nK +l]; // i=y a=xy b=y
            DDXX[6*nJ*1 +nJ*4 +nK +l] = FPXX[3*nJ*4 +nJ*1 +nK +l] +ABz*DPXX[3*nJ*1 +nJ*1 +nK +l]; // i=z a=xy b=y
            DDXX[6*nJ*1 +nJ*5 +nK +l] = FPXX[3*nJ*4 +nJ*2 +nK +l] +ABz*DPXX[3*nJ*1 +nJ*2 +nK +l]; // i=z a=xy b=z

            DDXX[6*nJ*2 +nJ*0 +nK +l] = FPXX[3*nJ*2 +nJ*0 +nK +l] +ABx*DPXX[3*nJ*2 +nJ*0 +nK +l]; // i=x a=xz b=x
            DDXX[6*nJ*2 +nJ*1 +nK +l] = FPXX[3*nJ*4 +nJ*0 +nK +l] +ABy*DPXX[3*nJ*2 +nJ*0 +nK +l]; // i=y a=xz b=x
            DDXX[6*nJ*2 +nJ*2 +nK +l] = FPXX[3*nJ*5 +nJ*0 +nK +l] +ABz*DPXX[3*nJ*2 +nJ*0 +nK +l]; // i=z a=xz b=x
            DDXX[6*nJ*2 +nJ*3 +nK +l] = FPXX[3*nJ*4 +nJ*1 +nK +l] +ABy*DPXX[3*nJ*2 +nJ*1 +nK +l]; // i=y a=xz b=y
            DDXX[6*nJ*2 +nJ*4 +nK +l] = FPXX[3*nJ*5 +nJ*1 +nK +l] +ABz*DPXX[3*nJ*2 +nJ*1 +nK +l]; // i=z a=xz b=y
            DDXX[6*nJ*2 +nJ*5 +nK +l] = FPXX[3*nJ*5 +nJ*2 +nK +l] +ABz*DPXX[3*nJ*2 +nJ*2 +nK +l]; // i=z a=xz b=z

            DDXX[6*nJ*3 +nJ*0 +nK +l] = FPXX[3*nJ*3 +nJ*0 +nK +l] +ABx*DPXX[3*nJ*3 +nJ*0 +nK +l]; // i=x a=yy b=x
            DDXX[6*nJ*3 +nJ*1 +nK +l] = FPXX[3*nJ*6 +nJ*0 +nK +l] +ABy*DPXX[3*nJ*3 +nJ*0 +nK +l]; // i=y a=yy b=x
            DDXX[6*nJ*3 +nJ*2 +nK +l] = FPXX[3*nJ*7 +nJ*0 +nK +l] +ABz*DPXX[3*nJ*3 +nJ*0 +nK +l]; // i=z a=yy b=x
            DDXX[6*nJ*3 +nJ*3 +nK +l] = FPXX[3*nJ*6 +nJ*1 +nK +l] +ABy*DPXX[3*nJ*3 +nJ*1 +nK +l]; // i=y a=yy b=y
            DDXX[6*nJ*3 +nJ*4 +nK +l] = FPXX[3*nJ*7 +nJ*1 +nK +l] +ABz*DPXX[3*nJ*3 +nJ*1 +nK +l]; // i=z a=yy b=y
            DDXX[6*nJ*3 +nJ*5 +nK +l] = FPXX[3*nJ*7 +nJ*2 +nK +l] +ABz*DPXX[3*nJ*3 +nJ*2 +nK +l]; // i=z a=yy b=z

            DDXX[6*nJ*4 +nJ*0 +nK +l] = FPXX[3*nJ*4 +nJ*0 +nK +l] +ABx*DPXX[3*nJ*4 +nJ*0 +nK +l]; // i=x a=yz b=x
            DDXX[6*nJ*4 +nJ*1 +nK +l] = FPXX[3*nJ*7 +nJ*0 +nK +l] +ABy*DPXX[3*nJ*4 +nJ*0 +nK +l]; // i=y a=yz b=x
            DDXX[6*nJ*4 +nJ*2 +nK +l] = FPXX[3*nJ*8 +nJ*0 +nK +l] +ABz*DPXX[3*nJ*4 +nJ*0 +nK +l]; // i=z a=yz b=x
            DDXX[6*nJ*4 +nJ*3 +nK +l] = FPXX[3*nJ*7 +nJ*1 +nK +l] +ABy*DPXX[3*nJ*4 +nJ*1 +nK +l]; // i=y a=yz b=y
            DDXX[6*nJ*4 +nJ*4 +nK +l] = FPXX[3*nJ*8 +nJ*1 +nK +l] +ABz*DPXX[3*nJ*4 +nJ*1 +nK +l]; // i=z a=yz b=y
            DDXX[6*nJ*4 +nJ*5 +nK +l] = FPXX[3*nJ*8 +nJ*2 +nK +l] +ABz*DPXX[3*nJ*4 +nJ*2 +nK +l]; // i=z a=yz b=z

            DDXX[6*nJ*5 +nJ*0 +nK +l] = FPXX[3*nJ*5 +nJ*0 +nK +l] +ABx*DPXX[3*nJ*5 +nJ*0 +nK +l]; // i=x a=zz b=x
            DDXX[6*nJ*5 +nJ*1 +nK +l] = FPXX[3*nJ*8 +nJ*0 +nK +l] +ABy*DPXX[3*nJ*5 +nJ*0 +nK +l]; // i=y a=zz b=x
            DDXX[6*nJ*5 +nJ*2 +nK +l] = FPXX[3*nJ*9 +nJ*0 +nK +l] +ABz*DPXX[3*nJ*5 +nJ*0 +nK +l]; // i=z a=zz b=x
            DDXX[6*nJ*5 +nJ*3 +nK +l] = FPXX[3*nJ*8 +nJ*1 +nK +l] +ABy*DPXX[3*nJ*5 +nJ*1 +nK +l]; // i=y a=zz b=y
            DDXX[6*nJ*5 +nJ*4 +nK +l] = FPXX[3*nJ*9 +nJ*1 +nK +l] +ABz*DPXX[3*nJ*5 +nJ*1 +nK +l]; // i=z a=zz b=y
            DDXX[6*nJ*5 +nJ*5 +nK +l] = FPXX[3*nJ*9 +nJ*2 +nK +l] +ABz*DPXX[3*nJ*5 +nJ*2 +nK +l]; // i=z a=zz b=z
        }
    }
}


void DfTEI::HRR_XXPP(const double CDx, const double CDy, const double CDz,
                     const double* XXDS, const double* XXPS,
                     const int nMaxI, const int nMaxJ, double* XXPP)
{
    const int nI = nMaxJ *3 *3;
    const int nID = nMaxJ *6;
    const int nIP = nMaxJ *3;
    for (int i = 0; i < nMaxI; ++i) {
        for (int j = 0; j < nMaxJ; ++j) {
            XXPP[nI*i +9*j +3*0 +0] = XXDS[nID*i +6*j +0] + CDx*XXPS[nIP*i +3*j +0]; // -- x x
            XXPP[nI*i +9*j +3*0 +1] = XXDS[nID*i +6*j +1] + CDy*XXPS[nIP*i +3*j +0]; // -- x y
            XXPP[nI*i +9*j +3*0 +2] = XXDS[nID*i +6*j +2] + CDz*XXPS[nIP*i +3*j +0]; // -- x z
            XXPP[nI*i +9*j +3*1 +0] = XXDS[nID*i +6*j +1] + CDx*XXPS[nIP*i +3*j +1]; // -- y x
            XXPP[nI*i +9*j +3*1 +1] = XXDS[nID*i +6*j +3] + CDy*XXPS[nIP*i +3*j +1]; // -- y y
            XXPP[nI*i +9*j +3*1 +2] = XXDS[nID*i +6*j +4] + CDz*XXPS[nIP*i +3*j +1]; // -- y z
            XXPP[nI*i +9*j +3*2 +0] = XXDS[nID*i +6*j +2] + CDx*XXPS[nIP*i +3*j +2]; // -- z x
            XXPP[nI*i +9*j +3*2 +1] = XXDS[nID*i +6*j +4] + CDy*XXPS[nIP*i +3*j +2]; // -- z y
            XXPP[nI*i +9*j +3*2 +2] = XXDS[nID*i +6*j +5] + CDz*XXPS[nIP*i +3*j +2]; // -- z z
        }
    }
}


void DfTEI::HRR_XXDP(const double CDx, const double CDy, const double CDz,
                     const double* XXFS, const double* XXDS,
                     const int nMaxI, const int nMaxJ, double* XXDP)
{
    // 0xxx 1xxy 2xxz 3xyy 4xyz 5xzz
    // 6yyy 7yyz 8yzz 9zzz
    const int nI = nMaxJ *6 *3;
    const int nIF = nMaxJ *10;
    const int nID = nMaxJ *6;
    for (int i = 0; i < nMaxI; ++i) {
        for (int j = 0; j < nMaxJ; ++j) {
            XXDP[nI*i +18*j +3*0 +0] = XXFS[nIF*i +10*j +0] + CDx*XXDS[nID*i +6*j +0]; // i=x c=xx
            XXDP[nI*i +18*j +3*0 +1] = XXFS[nIF*i +10*j +1] + CDy*XXDS[nID*i +6*j +0]; // i=y
            XXDP[nI*i +18*j +3*0 +2] = XXFS[nIF*i +10*j +2] + CDz*XXDS[nID*i +6*j +0]; // i=z

            XXDP[nI*i +18*j +3*1 +0] = XXFS[nIF*i +10*j +1] + CDx*XXDS[nID*i +6*j +1]; // i=x c=xy
            XXDP[nI*i +18*j +3*1 +1] = XXFS[nIF*i +10*j +3] + CDy*XXDS[nID*i +6*j +1]; // i=y
            XXDP[nI*i +18*j +3*1 +2] = XXFS[nIF*i +10*j +4] + CDz*XXDS[nID*i +6*j +1]; // i=z

            XXDP[nI*i +18*j +3*2 +0] = XXFS[nIF*i +10*j +2] + CDx*XXDS[nID*i +6*j +2]; // i=x c=xz
            XXDP[nI*i +18*j +3*2 +1] = XXFS[nIF*i +10*j +4] + CDy*XXDS[nID*i +6*j +2]; // i=y
            XXDP[nI*i +18*j +3*2 +2] = XXFS[nIF*i +10*j +5] + CDz*XXDS[nID*i +6*j +2]; // i=z

            XXDP[nI*i +18*j +3*3 +0] = XXFS[nIF*i +10*j +3] + CDx*XXDS[nID*i +6*j +3]; // i=x c=yy
            XXDP[nI*i +18*j +3*3 +1] = XXFS[nIF*i +10*j +6] + CDy*XXDS[nID*i +6*j +3]; // i=y
            XXDP[nI*i +18*j +3*3 +2] = XXFS[nIF*i +10*j +7] + CDz*XXDS[nID*i +6*j +3]; // i=z

            XXDP[nI*i +18*j +3*4 +0] = XXFS[nIF*i +10*j +4] + CDx*XXDS[nID*i +6*j +4]; // i=x c=yz
            XXDP[nI*i +18*j +3*4 +1] = XXFS[nIF*i +10*j +7] + CDy*XXDS[nID*i +6*j +4]; // i=y
            XXDP[nI*i +18*j +3*4 +2] = XXFS[nIF*i +10*j +8] + CDz*XXDS[nID*i +6*j +4]; // i=z

            XXDP[nI*i +18*j +3*5 +0] = XXFS[nIF*i +10*j +5] + CDx*XXDS[nID*i +6*j +5]; // i=x c=zz
            XXDP[nI*i +18*j +3*5 +1] = XXFS[nIF*i +10*j +8] + CDy*XXDS[nID*i +6*j +5]; // i=y
            XXDP[nI*i +18*j +3*5 +2] = XXFS[nIF*i +10*j +9] + CDz*XXDS[nID*i +6*j +5]; // i=z
        }
    }
}


void DfTEI::HRR_XXFP(const double CDx, const double CDy, const double CDz,
                     const double* XXGS, const double* XXFS,
                     const int nMaxI, const int nMaxJ, double* XXFP)
{
    const int nI = nMaxJ *10 *3;
    const int nIF = nMaxJ *15;
    const int nID = nMaxJ *10;

    //  0:xxxx  1:xxxy  2:xxxz  3:xxyy  4:xxyz
    //  5:xxzz  6:xyyy  7:xyyz  8:xyzz  9:xzzz
    // 10:yyyy 11:yyyz 12:yyzz 13:yzzz 14:zzzz

    for (int i = 0; i < nMaxI; ++i) {
        for (int j = 0; j < nMaxJ; ++j) {
            XXFP[nI*i +30*j +3*0 +0] = XXGS[nIF*i +15*j + 0] + CDx*XXFS[nID*i +10*j +0]; // i=x c=xxx
            XXFP[nI*i +30*j +3*0 +1] = XXGS[nIF*i +15*j + 1] + CDy*XXFS[nID*i +10*j +0]; // i=y
            XXFP[nI*i +30*j +3*0 +2] = XXGS[nIF*i +15*j + 2] + CDz*XXFS[nID*i +10*j +0]; // i=z

            XXFP[nI*i +30*j +3*1 +0] = XXGS[nIF*i +15*j + 1] + CDx*XXFS[nID*i +10*j +1]; // i=x c=xxy
            XXFP[nI*i +30*j +3*1 +1] = XXGS[nIF*i +15*j + 3] + CDy*XXFS[nID*i +10*j +1]; // i=y
            XXFP[nI*i +30*j +3*1 +2] = XXGS[nIF*i +15*j + 4] + CDz*XXFS[nID*i +10*j +1]; // i=z

            XXFP[nI*i +30*j +3*2 +0] = XXGS[nIF*i +15*j + 2] + CDx*XXFS[nID*i +10*j +2]; // i=x c=xxz
            XXFP[nI*i +30*j +3*2 +1] = XXGS[nIF*i +15*j + 4] + CDy*XXFS[nID*i +10*j +2]; // i=y
            XXFP[nI*i +30*j +3*2 +2] = XXGS[nIF*i +15*j + 5] + CDz*XXFS[nID*i +10*j +2]; // i=z

            XXFP[nI*i +30*j +3*3 +0] = XXGS[nIF*i +15*j + 3] + CDx*XXFS[nID*i +10*j +3]; // i=x c=xyy
            XXFP[nI*i +30*j +3*3 +1] = XXGS[nIF*i +15*j + 6] + CDy*XXFS[nID*i +10*j +3]; // i=y
            XXFP[nI*i +30*j +3*3 +2] = XXGS[nIF*i +15*j + 7] + CDz*XXFS[nID*i +10*j +3]; // i=z

            XXFP[nI*i +30*j +3*4 +0] = XXGS[nIF*i +15*j + 4] + CDx*XXFS[nID*i +10*j +4]; // i=x c=xyz
            XXFP[nI*i +30*j +3*4 +1] = XXGS[nIF*i +15*j + 7] + CDy*XXFS[nID*i +10*j +4]; // i=y
            XXFP[nI*i +30*j +3*4 +2] = XXGS[nIF*i +15*j + 8] + CDz*XXFS[nID*i +10*j +4]; // i=z

            XXFP[nI*i +30*j +3*5 +0] = XXGS[nIF*i +15*j + 5] + CDx*XXFS[nID*i +10*j +5]; // i=x c=xzz
            XXFP[nI*i +30*j +3*5 +1] = XXGS[nIF*i +15*j + 8] + CDy*XXFS[nID*i +10*j +5]; // i=y
            XXFP[nI*i +30*j +3*5 +2] = XXGS[nIF*i +15*j + 9] + CDz*XXFS[nID*i +10*j +5]; // i=z

            XXFP[nI*i +30*j +3*6 +0] = XXGS[nIF*i +15*j + 6] + CDx*XXFS[nID*i +10*j +6]; // i=x c=yyy
            XXFP[nI*i +30*j +3*6 +1] = XXGS[nIF*i +15*j +10] + CDy*XXFS[nID*i +10*j +6]; // i=y
            XXFP[nI*i +30*j +3*6 +2] = XXGS[nIF*i +15*j +11] + CDz*XXFS[nID*i +10*j +6]; // i=z

            XXFP[nI*i +30*j +3*7 +0] = XXGS[nIF*i +15*j + 7] + CDx*XXFS[nID*i +10*j +7]; // i=x c=yyz
            XXFP[nI*i +30*j +3*7 +1] = XXGS[nIF*i +15*j +11] + CDy*XXFS[nID*i +10*j +7]; // i=y
            XXFP[nI*i +30*j +3*7 +2] = XXGS[nIF*i +15*j +12] + CDz*XXFS[nID*i +10*j +7]; // i=z

            XXFP[nI*i +30*j +3*8 +0] = XXGS[nIF*i +15*j + 8] + CDx*XXFS[nID*i +10*j +8]; // i=x c=yzz
            XXFP[nI*i +30*j +3*8 +1] = XXGS[nIF*i +15*j +12] + CDy*XXFS[nID*i +10*j +8]; // i=y
            XXFP[nI*i +30*j +3*8 +2] = XXGS[nIF*i +15*j +13] + CDz*XXFS[nID*i +10*j +8]; // i=z

            XXFP[nI*i +30*j +3*9 +0] = XXGS[nIF*i +15*j + 9] + CDx*XXFS[nID*i +10*j +9]; // i=x c=zzz
            XXFP[nI*i +30*j +3*9 +1] = XXGS[nIF*i +15*j +13] + CDy*XXFS[nID*i +10*j +9]; // i=y
            XXFP[nI*i +30*j +3*9 +2] = XXGS[nIF*i +15*j +14] + CDz*XXFS[nID*i +10*j +9]; // i=z
        }
    }
}


void DfTEI::HRR_XXDD(const double CDx, const double CDy, const double CDz,
                     const double* XXFP, const double* XXDP,
                     const int nMaxI, const int nMaxJ, double* XXDD)
{
    const int nI = nMaxJ *6 *6;
    const int nIF = nMaxJ *3 *10;
    const int nID = nMaxJ *3 *6;
    for (int i = 0; i < nMaxI; ++i) {
        for (int j = 0; j < nMaxJ; ++j) {
            XXDD[nI*i +36*j +6*0 +0] = XXFP[nIF*i +30*j +3*0 +0] + CDx*XXDP[nID*i +18*j +3*0 +0]; // i=x, xx xx
            XXDD[nI*i +36*j +6*0 +1] = XXFP[nIF*i +30*j +3*1 +0] + CDy*XXDP[nID*i +18*j +3*0 +0]; // i=y, xx xy
            XXDD[nI*i +36*j +6*0 +2] = XXFP[nIF*i +30*j +3*2 +0] + CDz*XXDP[nID*i +18*j +3*0 +0]; // i=z, xx xz
            XXDD[nI*i +36*j +6*0 +3] = XXFP[nIF*i +30*j +3*1 +1] + CDy*XXDP[nID*i +18*j +3*0 +1]; // i=y, xx yy
            XXDD[nI*i +36*j +6*0 +4] = XXFP[nIF*i +30*j +3*2 +1] + CDz*XXDP[nID*i +18*j +3*0 +1]; // i=y, xx yz
            XXDD[nI*i +36*j +6*0 +5] = XXFP[nIF*i +30*j +3*2 +2] + CDz*XXDP[nID*i +18*j +3*0 +2]; // i=z, xx zz

            XXDD[nI*i +36*j +6*1 +0] = XXFP[nIF*i +30*j +3*1 +0] + CDx*XXDP[nID*i +18*j +3*1 +0]; // i=x, xy xx
            XXDD[nI*i +36*j +6*1 +1] = XXFP[nIF*i +30*j +3*3 +0] + CDy*XXDP[nID*i +18*j +3*1 +0]; // i=y, xy xy
            XXDD[nI*i +36*j +6*1 +2] = XXFP[nIF*i +30*j +3*4 +0] + CDz*XXDP[nID*i +18*j +3*1 +0]; // i=z, xy xz
            XXDD[nI*i +36*j +6*1 +3] = XXFP[nIF*i +30*j +3*3 +1] + CDy*XXDP[nID*i +18*j +3*1 +1]; // i=y, xy yy
            XXDD[nI*i +36*j +6*1 +4] = XXFP[nIF*i +30*j +3*4 +1] + CDz*XXDP[nID*i +18*j +3*1 +1]; // i=y, xy yz
            XXDD[nI*i +36*j +6*1 +5] = XXFP[nIF*i +30*j +3*4 +2] + CDz*XXDP[nID*i +18*j +3*1 +2]; // i=z, xy zz

            XXDD[nI*i +36*j +6*2 +0] = XXFP[nIF*i +30*j +3*2 +0] + CDx*XXDP[nID*i +18*j +3*2 +0]; // i=x, xz xx
            XXDD[nI*i +36*j +6*2 +1] = XXFP[nIF*i +30*j +3*4 +0] + CDy*XXDP[nID*i +18*j +3*2 +0]; // i=y, xz xy
            XXDD[nI*i +36*j +6*2 +2] = XXFP[nIF*i +30*j +3*5 +0] + CDz*XXDP[nID*i +18*j +3*2 +0]; // i=z, xz xz
            XXDD[nI*i +36*j +6*2 +3] = XXFP[nIF*i +30*j +3*4 +1] + CDy*XXDP[nID*i +18*j +3*2 +1]; // i=y, xz yy
            XXDD[nI*i +36*j +6*2 +4] = XXFP[nIF*i +30*j +3*5 +1] + CDz*XXDP[nID*i +18*j +3*2 +1]; // i=y, xz yz
            XXDD[nI*i +36*j +6*2 +5] = XXFP[nIF*i +30*j +3*5 +2] + CDz*XXDP[nID*i +18*j +3*2 +2]; // i=z, xz zz

            XXDD[nI*i +36*j +6*3 +0] = XXFP[nIF*i +30*j +3*3 +0] + CDx*XXDP[nID*i +18*j +3*3 +0]; // i=x, yy xx
            XXDD[nI*i +36*j +6*3 +1] = XXFP[nIF*i +30*j +3*6 +0] + CDy*XXDP[nID*i +18*j +3*3 +0]; // i=y, yy xy
            XXDD[nI*i +36*j +6*3 +2] = XXFP[nIF*i +30*j +3*7 +0] + CDz*XXDP[nID*i +18*j +3*3 +0]; // i=z, yy xz
            XXDD[nI*i +36*j +6*3 +3] = XXFP[nIF*i +30*j +3*6 +1] + CDy*XXDP[nID*i +18*j +3*3 +1]; // i=y, yy yy
            XXDD[nI*i +36*j +6*3 +4] = XXFP[nIF*i +30*j +3*7 +1] + CDz*XXDP[nID*i +18*j +3*3 +1]; // i=y, yy yz
            XXDD[nI*i +36*j +6*3 +5] = XXFP[nIF*i +30*j +3*7 +2] + CDz*XXDP[nID*i +18*j +3*3 +2]; // i=z, yy zz

            XXDD[nI*i +36*j +6*4 +0] = XXFP[nIF*i +30*j +3*4 +0] + CDx*XXDP[nID*i +18*j +3*4 +0]; // i=x, yz xx
            XXDD[nI*i +36*j +6*4 +1] = XXFP[nIF*i +30*j +3*7 +0] + CDy*XXDP[nID*i +18*j +3*4 +0]; // i=y, yz xy
            XXDD[nI*i +36*j +6*4 +2] = XXFP[nIF*i +30*j +3*8 +0] + CDz*XXDP[nID*i +18*j +3*4 +0]; // i=z, yz xz
            XXDD[nI*i +36*j +6*4 +3] = XXFP[nIF*i +30*j +3*7 +1] + CDy*XXDP[nID*i +18*j +3*4 +1]; // i=y, yz yy
            XXDD[nI*i +36*j +6*4 +4] = XXFP[nIF*i +30*j +3*8 +1] + CDz*XXDP[nID*i +18*j +3*4 +1]; // i=y, yz yz
            XXDD[nI*i +36*j +6*4 +5] = XXFP[nIF*i +30*j +3*8 +2] + CDz*XXDP[nID*i +18*j +3*4 +2]; // i=z, yz zz

            XXDD[nI*i +36*j +6*5 +0] = XXFP[nIF*i +30*j +3*5 +0] + CDx*XXDP[nID*i +18*j +3*5 +0]; // i=x, zz xx
            XXDD[nI*i +36*j +6*5 +1] = XXFP[nIF*i +30*j +3*8 +0] + CDy*XXDP[nID*i +18*j +3*5 +0]; // i=y, zz xy
            XXDD[nI*i +36*j +6*5 +2] = XXFP[nIF*i +30*j +3*9 +0] + CDz*XXDP[nID*i +18*j +3*5 +0]; // i=z, zz xz
            XXDD[nI*i +36*j +6*5 +3] = XXFP[nIF*i +30*j +3*8 +1] + CDy*XXDP[nID*i +18*j +3*5 +1]; // i=y, zz yy
            XXDD[nI*i +36*j +6*5 +4] = XXFP[nIF*i +30*j +3*9 +1] + CDz*XXDP[nID*i +18*j +3*5 +1]; // i=y, zz yz
            XXDD[nI*i +36*j +6*5 +5] = XXFP[nIF*i +30*j +3*9 +2] + CDz*XXDP[nID*i +18*j +3*5 +2]; // i=z, zz zz
        }
    }
}


void DfTEI::convert6Dto5D_I(const int nMaxJ, const int nMaxK, const int nMaxL,
                            const double* pIn, double* pOut)
{

    const int nMaxKL = nMaxK * nMaxL;
    const int nMaxJKL = nMaxJ * nMaxKL;
    for (int j = 0; j < nMaxJ; ++j) {
        for (int k = 0; k < nMaxK; ++k) {
            for (int l = 0; l < nMaxL; ++l) {
                pOut[nMaxJKL*0 +nMaxKL*j +nMaxL*k +l] = pIn[nMaxJKL*1 +nMaxKL*j +nMaxL*k +l]; // xy
                pOut[nMaxJKL*1 +nMaxKL*j +nMaxL*k +l] = pIn[nMaxJKL*2 +nMaxKL*j +nMaxL*k +l]; // xz
                pOut[nMaxJKL*2 +nMaxKL*j +nMaxL*k +l] = pIn[nMaxJKL*4 +nMaxKL*j +nMaxL*k +l]; // yz
                pOut[nMaxJKL*3 +nMaxKL*j +nMaxL*k +l] = 0.5*(+pIn[nMaxJKL*0 +nMaxKL*j +nMaxL*k +l]
                                                             -pIn[nMaxJKL*3 +nMaxKL*j +nMaxL*k +l]); // xx-yy
                pOut[nMaxJKL*4 +nMaxKL*j +nMaxL*k +l] = INV_SQRT3*(+pIn[nMaxJKL*5 +nMaxKL*j +nMaxL*k +l]
                                                                   -0.5*(+pIn[nMaxJKL*0 +nMaxKL*j +nMaxL*k +l]
                                                                         +pIn[nMaxJKL*3 +nMaxKL*j +nMaxL*k +l])); // 2zz-(xx+yy)
            }
        }
    }
}


void DfTEI::convert6Dto5D_J(const int nMaxI, const int nMaxK, const int nMaxL,
                            const double* pIn, double* pOut)
{

    const int nMaxKL = nMaxK * nMaxL;
    const int nMaxJKL_Out = 5 * nMaxKL;
    const int nMaxJKL_In = 6 * nMaxKL;
    for (int i = 0; i < nMaxI; ++i) {
        for (int k = 0; k < nMaxK; ++k) {
            for (int l = 0; l < nMaxL; ++l) {
                pOut[nMaxJKL_Out*i +nMaxKL*0 +nMaxL*k +l] = pIn[nMaxJKL_In*i +nMaxKL*1 +nMaxL*k +l]; // xy
                pOut[nMaxJKL_Out*i +nMaxKL*1 +nMaxL*k +l] = pIn[nMaxJKL_In*i +nMaxKL*2 +nMaxL*k +l]; // xz
                pOut[nMaxJKL_Out*i +nMaxKL*2 +nMaxL*k +l] = pIn[nMaxJKL_In*i +nMaxKL*4 +nMaxL*k +l]; // yz
                pOut[nMaxJKL_Out*i +nMaxKL*3 +nMaxL*k +l] = 0.5*(+pIn[nMaxJKL_In*i +nMaxKL*0 +nMaxL*k +l]
                                                                 -pIn[nMaxJKL_In*i +nMaxKL*3 +nMaxL*k +l]); // xx-yy
                pOut[nMaxJKL_Out*i +nMaxKL*4 +nMaxL*k +l] = INV_SQRT3*(+pIn[nMaxJKL_In*i +nMaxKL*5 +nMaxL*k +l]
                                                                       -0.5*(+pIn[nMaxJKL_In*i +nMaxKL*0 +nMaxL*k +l]
                                                                             +pIn[nMaxJKL_In*i +nMaxKL*3 +nMaxL*k +l])); // 2zz-(xx+yy)
            }
        }
    }
}


void DfTEI::convert6Dto5D_K(const int nMaxI, const int nMaxJ, const int nMaxL,
                            const double* pIn, double* pOut)
{

    const int nMaxKL_Out = 5 * nMaxL;
    const int nMaxKL_In = 6 * nMaxL;
    const int nMaxJKL_Out = nMaxJ * nMaxKL_Out;
    const int nMaxJKL_In = nMaxJ * nMaxKL_In;
    for (int i = 0; i < nMaxI; ++i) {
        for (int j = 0; j < nMaxJ; ++j) {
            for (int l = 0; l < nMaxL; ++l) {
                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +nMaxL*0 +l] = pIn[nMaxJKL_In*i +nMaxKL_In*j +nMaxL*1 +l]; // xy
                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +nMaxL*1 +l] = pIn[nMaxJKL_In*i +nMaxKL_In*j +nMaxL*2 +l]; // xz
                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +nMaxL*2 +l] = pIn[nMaxJKL_In*i +nMaxKL_In*j +nMaxL*4 +l]; // yz
                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +nMaxL*3 +l] = 0.5*(+pIn[nMaxJKL_In*i +nMaxKL_In*j +nMaxL*0 +l]
                                                                     -pIn[nMaxJKL_In*i +nMaxKL_In*j +nMaxL*3 +l]); // xx-yy
                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +nMaxL*4 +l] = INV_SQRT3*(+pIn[nMaxJKL_In*i +nMaxKL_In*j +nMaxL*5 +l]
                                                                           -0.5*(+pIn[nMaxJKL_In*i +nMaxKL_In*j +nMaxL*0 +l]
                                                                                 +pIn[nMaxJKL_In*i +nMaxKL_In*j +nMaxL*3 +l])); // 2zz-(xx+yy)
            }
        }
    }
}


void DfTEI::convert6Dto5D_L(const int nMaxI, const int nMaxJ, const int nMaxK,
                            const double* pIn, double* pOut)
{
    const int nMaxKL_Out = nMaxK * 5;
    const int nMaxKL_In  = nMaxK * 6;
    const int nMaxJKL_Out = nMaxJ * nMaxKL_Out;
    const int nMaxJKL_In  = nMaxJ * nMaxKL_In;

    for (int i = 0; i < nMaxI; ++i) {
        for (int j = 0; j < nMaxJ; ++j) {
            for (int k = 0; k < nMaxK; ++k) {
                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +5*k +0] = pIn[nMaxJKL_In*i +nMaxKL_In*j +6*k +1]; // xy
                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +5*k +1] = pIn[nMaxJKL_In*i +nMaxKL_In*j +6*k +2]; // xz
                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +5*k +2] = pIn[nMaxJKL_In*i +nMaxKL_In*j +6*k +4]; // yz

                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +5*k +3] = 0.5*(+pIn[nMaxJKL_In*i +nMaxKL_In*j +6*k +0]
                                                                 -pIn[nMaxJKL_In*i +nMaxKL_In*j +6*k +3]); // xx-yy

                pOut[nMaxJKL_Out*i +nMaxKL_Out*j +5*k +4] = INV_SQRT3*(+pIn[nMaxJKL_In*i +nMaxKL_In*j +6*k +5]
                                                                       -0.5*(+pIn[nMaxJKL_In*i +nMaxKL_In*j +6*k +0]
                                                                             +pIn[nMaxJKL_In*i +nMaxKL_In*j +6*k +3])); // 2zz-(xx+yy)
            }
        }
    }
}


bool DfTEI::isPrimitiveShellsCutoff() const
{
    bool bAnswer = false;
    if (std::fabs(this->pSSSS[0]) < this->primitiveIntegralsCutoffThreshold_) {
        bAnswer = true;
    }

    return bAnswer;
}


/////////////////////////////////////////////////////////////////////////
// contract SSSS
void DfTEI::contractSSSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double SSSS = 0.0;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {
        this->posP = ij->P;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 1);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            SSSS += this->pSSSS[0];
        }
    }

    this->ERI[0] = SSSS;
}


/////////////////////////////////////////////////////////////////////////
// contract SSSP
void DfTEI::contractSSSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSSPS(IJ, this->swapShellPair(KL));
}


/////////////////////////////////////////////////////////////////////////
// constract SSSD
void DfTEI::contractSSSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSSDS(IJ, this->swapShellPair(KL));
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// constract SSPS
void DfTEI::contractSSPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSSS(KL, IJ); // exchange
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// constract SSPP
void DfTEI::contractSSPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPPSS(KL, IJ); // exchange
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// constract SSPD
void DfTEI::contractSSPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSSDP(IJ, this->swapShellPair(KL));
    this->transform(5, 3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// constract SSDS
void DfTEI::contractSSDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSSS(KL, IJ); // << exchange
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// constract SSDP
void DfTEI::contractSSDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPSS(KL, IJ); // << exchange
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// constract SSDD
void DfTEI::contractSSDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDSS(KL, IJ); // << exchange
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SPSS
void DfTEI::contractSPSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSSS(this->swapShellPair(IJ), KL);
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SPSP
void DfTEI::contractSPSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSPPS(IJ, this->swapShellPair(KL));
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// // contract SPSD
void DfTEI::contractSPSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSPDS(IJ, this->swapShellPair(KL));
    // no transform
}

/////////////////////////////////////////////////////////////////////////
// contract SPPS
void DfTEI::contractSPPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSPS(this->swapShellPair(IJ), KL);
    // no transform
}

/////////////////////////////////////////////////////////////////////////
// contract SPPP
void DfTEI::contractSPPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSPP(this->swapShellPair(IJ), KL);
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SPPD
void DfTEI::contractSPPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSPDP(IJ, this->swapShellPair(KL));
    for (int i = 0; i < 3; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
}


/////////////////////////////////////////////////////////////////////////
// contract SPDS
void DfTEI::contractSPDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSDS(this->swapShellPair(IJ), KL);
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SPDP
void DfTEI::contractSPDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSDP(this->swapShellPair(IJ), KL);
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SPDD
void DfTEI::contractSPDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSDD(this->swapShellPair(IJ), KL);
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// constract SDSS
void DfTEI::contractSDSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSSS(this->swapShellPair(IJ), KL);
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SDSP
void DfTEI::contractSDSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSDPS(IJ, this->swapShellPair(KL));
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SDSD
void DfTEI::contractSDSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSDDS(IJ, this->swapShellPair(KL));
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SDPS
void DfTEI::contractSDPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSPS(this->swapShellPair(IJ), KL);
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SDPP
void DfTEI::contractSDPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSPP(this->swapShellPair(IJ), KL);
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract SDPD
void DfTEI::contractSDPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractSDDP(IJ, this->swapShellPair(KL));
    for (int i = 0; i < 5; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
}


/////////////////////////////////////////////////////////////////////////
// contract SDDS
void DfTEI::contractSDDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSDS(this->swapShellPair(IJ), KL);
}


/////////////////////////////////////////////////////////////////////////
// contract SDDP
void DfTEI::contractSDDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSDP(this->swapShellPair(IJ), KL);
}


/////////////////////////////////////////////////////////////////////////
// contract SDDD
void DfTEI::contractSDDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSDD(this->swapShellPair(IJ), KL);
}


/////////////////////////////////////////////////////////////////////////
// contract PSSS
void DfTEI::contractPSSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double PSSS[3];
    for (int i = 0; i < 3; ++i) {
        PSSS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {
        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 2);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->WP = this->posW - this->posP;

            this->primitivePSSS(*ij, *kl, 1);

            for (int i = 0; i < 3; ++i) {
                PSSS[i] += this->pPSSS[0][i];
            }
        }
    }

    for (int i = 0; i < 3; ++i) {
        this->ERI[i] = PSSS[i];
    }
}


/////////////////////////////////////////////////////////////////////////
// contract PSSP
void DfTEI::contractPSSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSPS(IJ, this->swapShellPair(KL));
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract PSSD
void DfTEI::contractPSSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSDS(IJ, this->swapShellPair(KL));
}


/////////////////////////////////////////////////////////////////////////
// contract PSPS
void DfTEI::contractPSPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double PSPS[9];
    for (int i = 0; i < 9; ++i) {
        PSPS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {
        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 3);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->WP = this->posW - this->posP;

            this->primitivePSSS(*ij, *kl, 2);

            this->QC = this->posQ - posC;
            this->WQ = this->posW - this->posQ;

            this->primitivePSPS(*ij, *kl, 1);

            for (int i = 0; i < 9; ++i) {
                PSPS[i] += this->pPSPS[0][i];
            }
        }
    }

    for (int i = 0; i < 9; ++i) {
        this->ERI[i] = PSPS[i];
    }
}


/////////////////////////////////////////////////////////////////////////
// contract PSPP
void DfTEI::contractPSPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPPPS(KL, IJ); // exchange
    this->transform(3*3, 3, this->ERI); // transform ERI
}


/////////////////////////////////////////////////////////////////////////
// contract PSPD
void DfTEI::contractPSPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPSDP(IJ, this->swapShellPair(KL));
    for (int i = 0; i < 3; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
}


/////////////////////////////////////////////////////////////////////////
// contract PSDS
void DfTEI::contractPSDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSPS(KL, IJ);
    this->transform(5, 3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// constract PSDP
void DfTEI::contractPSDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPPS(KL, IJ);
    this->transform(5*3, 3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// constract PSDD
void DfTEI::contractPSDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDPS(KL, IJ);
    this->transform(5*5, 3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// constract PPSS
void DfTEI::contractPPSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double PSSS[3];
    for (int i = 0; i < 3; ++i) {
        PSSS[i] = 0.0;
    }
    double DSSS[6];
    for (int i = 0; i < 6; ++i) {
        DSSS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    //const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {
        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 3);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->WP = this->posW - this->posP;
            //this->WQ = this->posW - this->posQ;
            //this->QC = this->posQ - this->posC;

            this->primitivePSSS(*ij, *kl, 2);
            this->primitiveDSSS(*ij, *kl, 1);

            for (int i = 0; i < 3; ++i) {
                PSSS[i] += this->pPSSS[0][i];
            }
            for (int i = 0; i < 6; ++i) {
                DSSS[i] += this->pDSSS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();

    // HRR(PPSS)
    double PPSS[9]; // 3*3
    PPSS[0] = DSSS[0] + ABx * PSSS[0]; // x x
    PPSS[1] = DSSS[1] + ABy * PSSS[0]; // x y
    PPSS[2] = DSSS[2] + ABz * PSSS[0]; // x z
    PPSS[3] = DSSS[1] + ABx * PSSS[1]; // y x
    PPSS[4] = DSSS[3] + ABy * PSSS[1]; // y y
    PPSS[5] = DSSS[4] + ABz * PSSS[1]; // y z
    PPSS[6] = DSSS[2] + ABx * PSSS[2]; // z x
    PPSS[7] = DSSS[4] + ABy * PSSS[2]; // z y
    PPSS[8] = DSSS[5] + ABz * PSSS[2]; // z z

    // ERI
    for (int i = 0; i < 9; ++i) {
        this->ERI[i] = PPSS[i];
    }
}


/////////////////////////////////////////////////////////////////////////
// contract PPSP
void DfTEI::contractPPSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPPPS(IJ, this->swapShellPair(KL));
    // no transform
}


/////////////////////////////////////////////////////////////////////////
// contract PPSD
void DfTEI::contractPPSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPPDS(IJ, this->swapShellPair(KL));
}


/////////////////////////////////////////////////////////////////////////
// contract PPPS
void DfTEI::contractPPPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double PSPS[9];
    for (int i = 0; i < 9; ++i) {
        PSPS[i] = 0.0;
    }
    double DSPS[18];
    for (int i = 0; i < 18; ++i) {
        DSPS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {
        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 4);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->WP = this->posW - this->posP;

            this->primitivePSSS(*ij, *kl, 3);
            this->primitiveDSSS(*ij, *kl, 2);

            this->QC = this->posQ - posC;
            this->WQ = this->posW - this->posQ;

            this->primitivePSPS(*ij, *kl, 1);
            this->primitiveDSPS(*ij, *kl, 1);

            for (int i = 0; i < 9; ++i) {
                PSPS[i] += this->pPSPS[0][i];
            }
            for (int i = 0; i < 18; ++i) {
                DSPS[i] += this->pDSPS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();

    // HRR
    double PPPS[27]; // 3 * 3 * 3
    PPPS[9*0 +3*0 +0] = DSPS[3*0 +0] + ABx * PSPS[3*0 +0]; // xxx
    PPPS[9*0 +3*0 +1] = DSPS[3*0 +1] + ABx * PSPS[3*0 +1]; // xxy
    PPPS[9*0 +3*0 +2] = DSPS[3*0 +2] + ABx * PSPS[3*0 +2]; // xxz
    PPPS[9*0 +3*1 +0] = DSPS[3*1 +0] + ABy * PSPS[3*0 +0]; // xyx
    PPPS[9*0 +3*1 +1] = DSPS[3*1 +1] + ABy * PSPS[3*0 +1]; // xyy
    PPPS[9*0 +3*1 +2] = DSPS[3*1 +2] + ABy * PSPS[3*0 +2]; // xyz
    PPPS[9*0 +3*2 +0] = DSPS[3*2 +0] + ABz * PSPS[3*0 +0]; // xzx
    PPPS[9*0 +3*2 +1] = DSPS[3*2 +1] + ABz * PSPS[3*0 +1]; // xzy
    PPPS[9*0 +3*2 +2] = DSPS[3*2 +2] + ABz * PSPS[3*0 +2]; // xzz

    PPPS[9*1 +3*0 +0] = DSPS[3*1 +0] + ABx * PSPS[3*1 +0]; // yxx
    PPPS[9*1 +3*0 +1] = DSPS[3*1 +1] + ABx * PSPS[3*1 +1]; // yxy
    PPPS[9*1 +3*0 +2] = DSPS[3*1 +2] + ABx * PSPS[3*1 +2]; // yxz
    PPPS[9*1 +3*1 +0] = DSPS[3*3 +0] + ABy * PSPS[3*1 +0]; // yyx
    PPPS[9*1 +3*1 +1] = DSPS[3*3 +1] + ABy * PSPS[3*1 +1]; // yyy
    PPPS[9*1 +3*1 +2] = DSPS[3*3 +2] + ABy * PSPS[3*1 +2]; // yyz
    PPPS[9*1 +3*2 +0] = DSPS[3*4 +0] + ABz * PSPS[3*1 +0]; // yzx
    PPPS[9*1 +3*2 +1] = DSPS[3*4 +1] + ABz * PSPS[3*1 +1]; // yzy
    PPPS[9*1 +3*2 +2] = DSPS[3*4 +2] + ABz * PSPS[3*1 +2]; // yzz

    PPPS[9*2 +3*0 +0] = DSPS[3*2 +0] + ABx * PSPS[3*2 +0]; // zxx
    PPPS[9*2 +3*0 +1] = DSPS[3*2 +1] + ABx * PSPS[3*2 +1]; // zxy
    PPPS[9*2 +3*0 +2] = DSPS[3*2 +2] + ABx * PSPS[3*2 +2]; // zxz
    PPPS[9*2 +3*1 +0] = DSPS[3*4 +0] + ABy * PSPS[3*2 +0]; // zyx
    PPPS[9*2 +3*1 +1] = DSPS[3*4 +1] + ABy * PSPS[3*2 +1]; // zyy
    PPPS[9*2 +3*1 +2] = DSPS[3*4 +2] + ABy * PSPS[3*2 +2]; // zyz
    PPPS[9*2 +3*2 +0] = DSPS[3*5 +0] + ABz * PSPS[3*2 +0]; // zzx
    PPPS[9*2 +3*2 +1] = DSPS[3*5 +1] + ABz * PSPS[3*2 +1]; // zzy
    PPPS[9*2 +3*2 +2] = DSPS[3*5 +2] + ABz * PSPS[3*2 +2]; // zzz

    for (int i = 0; i < 27; ++i) {
        this->ERI[i] = PPPS[i];
    }
}


/////////////////////////////////////////////////////////////////////////
// contract PPPP
void DfTEI::contractPPPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double PSPS[9];
    for (int i = 0; i < 9; ++i) {
        PSPS[i] = 0.0;
    }
    double PSDS[18];
    double DSPS[18];
    for (int i = 0; i < 18; ++i) {
        PSDS[i] = 0.0;
        DSPS[i] = 0.0;
    }
    double DSDS[36];
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {
        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 5);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 4);
            this->primitiveDSSS(*ij, *kl, 3);
            this->primitivePSPS(*ij, *kl, 2); // for DSDS & PSDS
            this->primitiveDSPS(*ij, *kl, 2);
            this->primitiveDSDS(*ij, *kl, 1);

            //this->primitiveSSSS(*ij, *kl, 4);
            //this->primitivePSSS(*ij, *kl, 3);
            //this->primitivePSPS(*ij, *kl, 2); // for DSDS & PSDS
            this->primitiveSSPS(*ij, *kl, 2); // for PSDS
            this->primitivePSDS(*ij, *kl, 1);

            for (int i = 0; i < 9; ++i) {
                PSPS[i] += this->pPSPS[0][i];
            }
            for (int i = 0; i < 18; ++i) {
                PSDS[i] += this->pPSDS[0][i];
                DSPS[i] += this->pDSPS[0][i];
            }
            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
    const double CDx = KL.AB.x();
    const double CDy = KL.AB.y();
    const double CDz = KL.AB.z();

    // HRR(PSPP)
    double PSPP[27]; // 3*3*3
    for (int i = 0; i < 3; ++i) {
        PSPP[9*i +3*0 +0] = PSDS[6*i +0] + CDx*PSPS[3*i +0]; // x xx
        PSPP[9*i +3*0 +1] = PSDS[6*i +1] + CDy*PSPS[3*i +0]; // x xy
        PSPP[9*i +3*0 +2] = PSDS[6*i +2] + CDz*PSPS[3*i +0]; // x xz
        PSPP[9*i +3*1 +0] = PSDS[6*i +1] + CDx*PSPS[3*i +1]; // x yx
        PSPP[9*i +3*1 +1] = PSDS[6*i +3] + CDy*PSPS[3*i +1]; // x yy
        PSPP[9*i +3*1 +2] = PSDS[6*i +4] + CDz*PSPS[3*i +1]; // x yz
        PSPP[9*i +3*2 +0] = PSDS[6*i +2] + CDx*PSPS[3*i +2]; // x zx
        PSPP[9*i +3*2 +1] = PSDS[6*i +4] + CDy*PSPS[3*i +2]; // x zy
        PSPP[9*i +3*2 +2] = PSDS[6*i +5] + CDz*PSPS[3*i +2]; // x zz
    }
    // std::cout << "PSPP=" << PSPP[0] << std::endl;

    // HRR(DSPP)
    double DSPP[54]; //6*3*3
    for (int i = 0; i < 6; ++i) {
        DSPP[9*i +3*0 +0] = DSDS[6*i +0] + CDx*DSPS[3*i +0]; // -- x x
        DSPP[9*i +3*0 +1] = DSDS[6*i +1] + CDy*DSPS[3*i +0]; // -- x y
        DSPP[9*i +3*0 +2] = DSDS[6*i +2] + CDz*DSPS[3*i +0]; // -- x z
        DSPP[9*i +3*1 +0] = DSDS[6*i +1] + CDx*DSPS[3*i +1]; // -- y x
        DSPP[9*i +3*1 +1] = DSDS[6*i +3] + CDy*DSPS[3*i +1]; // -- y y
        DSPP[9*i +3*1 +2] = DSDS[6*i +4] + CDz*DSPS[3*i +1]; // -- y z
        DSPP[9*i +3*2 +0] = DSDS[6*i +2] + CDx*DSPS[3*i +2]; // -- z x
        DSPP[9*i +3*2 +1] = DSDS[6*i +4] + CDy*DSPS[3*i +2]; // -- z y
        DSPP[9*i +3*2 +2] = DSDS[6*i +5] + CDz*DSPS[3*i +2]; // -- z z
    }

    // HRR(PPPP)
    double PPPP[81]; // 3*3*3*3
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            PPPP[27*0 +9*0 +3*i +j] = DSPP[9*0 +3*i +j] + ABx*PSPP[9*0 +3*i +j]; // xx --
            PPPP[27*0 +9*1 +3*i +j] = DSPP[9*1 +3*i +j] + ABy*PSPP[9*0 +3*i +j]; // xy --
            PPPP[27*0 +9*2 +3*i +j] = DSPP[9*2 +3*i +j] + ABz*PSPP[9*0 +3*i +j]; // xz --
            PPPP[27*1 +9*0 +3*i +j] = DSPP[9*1 +3*i +j] + ABx*PSPP[9*1 +3*i +j]; // yx --
            PPPP[27*1 +9*1 +3*i +j] = DSPP[9*3 +3*i +j] + ABy*PSPP[9*1 +3*i +j]; // yy --
            PPPP[27*1 +9*2 +3*i +j] = DSPP[9*4 +3*i +j] + ABz*PSPP[9*1 +3*i +j]; // yz --
            PPPP[27*2 +9*0 +3*i +j] = DSPP[9*2 +3*i +j] + ABx*PSPP[9*2 +3*i +j]; // xx --
            PPPP[27*2 +9*1 +3*i +j] = DSPP[9*4 +3*i +j] + ABy*PSPP[9*2 +3*i +j]; // xy --
            PPPP[27*2 +9*2 +3*i +j] = DSPP[9*5 +3*i +j] + ABz*PSPP[9*2 +3*i +j]; // xz --
        }
    }

    for (int i = 0; i < 81; ++i) {
        this->ERI[i] = PPPP[i];
    }
}


/////////////////////////////////////////////////////////////////////////
// contract PPPD
void DfTEI::contractPPPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPPDP(IJ, this->swapShellPair(KL));
    for (int i = 0; i < 3*3; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
}


/////////////////////////////////////////////////////////////////////////
// contract PPDS
void DfTEI::contractPPDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSPP(KL, IJ);
    this->transform(5, 3*3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// contract PPDP
void DfTEI::contractPPDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPPP(KL, IJ);
    this->transform(5*3, 3*3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// contract PPDD
void DfTEI::contractPPDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDPP(KL, IJ);
    this->transform(5*5, 3*3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// contract PDSS
void DfTEI::contractPDSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPSS(this->swapShellPair(IJ), KL);
    this->transform(5, 3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// contract PDSP
void DfTEI::contractPDSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPDPS(IJ, this->swapShellPair(KL));
}


/////////////////////////////////////////////////////////////////////////
// contract PDSD
void DfTEI::contractPDSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPDDS(IJ, this->swapShellPair(KL));
}


/////////////////////////////////////////////////////////////////////////
// contract PDPS
void DfTEI::contractPDPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPPS(this->swapShellPair(IJ), KL);
    // transform
    this->transform(5*3, 3, this->ERI);
    for (int i = 0; i < 3; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
    this->transform(3, 5*3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// contract PDPP
void DfTEI::contractPDPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPPP(this->swapShellPair(IJ), KL);
    // transform
    this->transform(5*3, 3*3, this->ERI);
    for (int i = 0; i < 3*3; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
    this->transform(3*3, 5*3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// contract PDPD
void DfTEI::contractPDPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractPDDP(IJ, this->swapShellPair(KL));
    // transform
    for (int i = 0; i < 5*3; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
}


/////////////////////////////////////////////////////////////////////////
// contract PDDS
void DfTEI::contractPDDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPDS(this->swapShellPair(IJ), KL);
    // transform
    this->transform(5*3, 5, this->ERI);
    for (int i = 0; i < 5; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
    this->transform(5, 5*3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// contract PDDP
void DfTEI::contractPDDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPDP(this->swapShellPair(IJ), KL);
    // transform
    this->transform(5*3, 5*3, this->ERI);
    for (int i = 0; i < 5*3; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
    this->transform(5*3, 5*3, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// contract PDDD
void DfTEI::contractPDDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDPD(KL, IJ);
    this->transform(5*5, 3*5, this->ERI);
}


/////////////////////////////////////////////////////////////////////////
// constract DSSS
void DfTEI::contractDSSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSSS[6];
    for (int i = 0; i < 6; ++i) {
        DSSS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    //const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 3);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            //this->QC = this->posQ - this->posC;
            this->WP = this->posW - this->posP;
            //this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 2);
            this->primitiveDSSS(*ij, *kl, 1);

            for (int i = 0; i < 6; ++i) {
                DSSS[i] += this->pDSSS[0][i];
            }
        }
    }

    // convert 6D to 5D
    this->ERI[ 0] = DSSS[ 1];                                         // xy
    this->ERI[ 1] = DSSS[ 2];                                         // xz
    this->ERI[ 2] = DSSS[ 4];                                         // yz
    this->ERI[ 3] = 0.5*(DSSS[ 0] - DSSS[ 3]);                        // xx - yy
    this->ERI[ 4] = INV_SQRT3*(DSSS[ 5] - 0.5*(DSSS[ 0] + DSSS[ 3])); // 2zz - (xx + yy)
}


/////////////////////////////////////////////////////////////////////////
// contract DSSP
void DfTEI::contractDSSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSPS(IJ, this->swapShellPair(KL));
}


/////////////////////////////////////////////////////////////////////////
// contract DSSD
void DfTEI::contractDSSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSDS(IJ, this->swapShellPair(KL));
}


/////////////////////////////////////////////////////////////////////////
// contract DSPS
void DfTEI::contractDSPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSPS[18];
    for (int i = 0; i < 18; ++i) {
        DSPS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 4);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 3);
            this->primitiveDSSS(*ij, *kl, 2);
            this->primitiveDSPS(*ij, *kl, 1);

            for (int i = 0; i < 18; ++i) {
                DSPS[i] += this->pDSPS[0][i];
            }
        }
    }

    // convert 6D to 5D
    for (int r = 0; r < 3; ++r) {
        this->ERI[3*0 +r] = DSPS[3*1 +r];                                                 // xy
        this->ERI[3*1 +r] = DSPS[3*2 +r];                                                 // xz
        this->ERI[3*2 +r] = DSPS[3*4 +r];                                                 // yz
        this->ERI[3*3 +r] = 0.5*(DSPS[3*0 +r] - DSPS[3*3 +r]);                            // xx - yy
        this->ERI[3*4 +r] = INV_SQRT3*(DSPS[3*5 +r] - 0.5*(DSPS[3*0 +r] + DSPS[3*3 +r])); // 2zz - (xx + yy)
    }
}


// contract DSPP
void DfTEI::contractDSPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSPS[18];
    for (int i = 0; i < 18; ++i) {
        DSPS[i] = 0.0;
    }
    double DSDS[36];
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 5);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 4);
            this->primitiveDSSS(*ij, *kl, 3);

            this->primitivePSPS(*ij, *kl, 2);
            this->primitiveDSPS(*ij, *kl, 2);
            this->primitiveDSDS(*ij, *kl, 1);

            for (int i = 0; i < 18; ++i) {
                DSPS[i] += this->pDSPS[0][i];
            }
            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
        }
    }

    const double CDx = KL.AB.x();
    const double CDy = KL.AB.y();
    const double CDz = KL.AB.z();

    // HRR (DSPP)
    double DSPP[54]; // 6*3*3
    for (int p = 0; p < 6; ++p) {
        DSPP[9*p +3*0 + 0] = DSDS[6*p +0] + CDx*DSPS[3*p +0]; // xx i=x
        DSPP[9*p +3*0 + 1] = DSDS[6*p +1] + CDy*DSPS[3*p +0]; // xy i=y
        DSPP[9*p +3*0 + 2] = DSDS[6*p +2] + CDz*DSPS[3*p +0]; // xz i=z
        DSPP[9*p +3*1 + 0] = DSDS[6*p +1] + CDx*DSPS[3*p +1]; // yx i=x
        DSPP[9*p +3*1 + 1] = DSDS[6*p +3] + CDy*DSPS[3*p +1]; // yy i=y
        DSPP[9*p +3*1 + 2] = DSDS[6*p +4] + CDz*DSPS[3*p +1]; // yz i=z
        DSPP[9*p +3*2 + 0] = DSDS[6*p +2] + CDx*DSPS[3*p +2]; // zx i=x
        DSPP[9*p +3*2 + 1] = DSDS[6*p +4] + CDy*DSPS[3*p +2]; // zy i=y
        DSPP[9*p +3*2 + 2] = DSDS[6*p +5] + CDz*DSPS[3*p +2]; // zz i=z
    }

    // convert 6D to 5D
    for (int r = 0; r < 3; ++r) {
        for (int s = 0; s < 3; ++s) {
            this->ERI[9*0 +3*r +s] = DSPP[9*1 +3*r +s];                                                           // xy
            this->ERI[9*1 +3*r +s] = DSPP[9*2 +3*r +s];                                                           // xz
            this->ERI[9*2 +3*r +s] = DSPP[9*4 +3*r +s];                                                           // yz
            this->ERI[9*3 +3*r +s] = 0.5*(DSPP[9*0 +3*r +s] - DSPP[9*3 +3*r +s]);                                 // xx - yy
            this->ERI[9*4 +3*r +s] = INV_SQRT3*(DSPP[9*5 +3*r +s] - 0.5*(DSPP[9*0 +3*r +s] + DSPP[9*3 +3*r +s])); // 2zz - (xx + yy)
        }
    }
}


////////////////////////////////////////////////////////////////////////
// contract DSPD
void DfTEI::contractDSPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDSDP(IJ, this->swapShellPair(KL));
    for (int i = 0; i < 5; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
}


////////////////////////////////////////////////////////////////////////
// contract DSDS
void DfTEI::contractDSDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSDS[36];
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 5);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 4);
            this->primitiveDSSS(*ij, *kl, 3);

            this->primitivePSPS(*ij, *kl, 2);
            this->primitiveDSPS(*ij, *kl, 2);
            this->primitiveDSDS(*ij, *kl, 1);

            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
        }
    }

    // convert 6D to 5D
    double Dsds[30]; // 6*5
    for (int i = 0; i < 6; ++i) {
        Dsds[5*i +0] = DSDS[6*i +1]; // xy
        Dsds[5*i +1] = DSDS[6*i +2]; // xz
        Dsds[5*i +2] = DSDS[6*i +4]; // yz
        Dsds[5*i +3] = 0.5*(DSDS[6*i +0] -DSDS[6*i +3]); // xx-yy
        Dsds[5*i +4] = INV_SQRT3*(DSDS[6*i +5] -0.5*(DSDS[6*i +0] +DSDS[6*i +3])); // 2zz-(xx+yy)
    }

    for (int j = 0; j < 5; ++j) {
        this->ERI[5*0 +j] = Dsds[5*1 +j]; // xy
        this->ERI[5*1 +j] = Dsds[5*2 +j]; // xz
        this->ERI[5*2 +j] = Dsds[5*4 +j]; // yz
        this->ERI[5*3 +j] = 0.5*(Dsds[5*0 +j] -Dsds[5*3 +j]); // xx-yy
        this->ERI[5*4 +j] = INV_SQRT3*(Dsds[5*5 +j] -0.5*(Dsds[5*0 +j] +Dsds[5*3 +j])); // 2zz-(xx+yy)
    }

    // for debug
//   {
//     std::ofstream ofs("DSDS.txt", std::ios::out | std::ios::app);
//     int index = 0;
//     ofs << "6D" << std::endl;
//     for (int i = 0; i < 6; ++i){
//       for (int j = 0; j < 1; ++j){
//  for (int k = 0; k < 6; ++k){
//    for (int l = 0; l < 1; ++l){
//      ofs << TlUtils::format("(%2d %2d|%2d %2d) = % .8e", i, j, k, l, DSDS[index])
//      << std::endl;
//      ++index;
//    }
//  }
//       }
//     }
//
//     ofs << "5D" << std::endl;
//     for (int i = 0; i < 25; ++i){
//       ofs << TlUtils::format("%2d: % .8e", i, this->ERI[i]) << std::endl;
//     }
//   }


}


////////////////////////////////////////////////////////////////////////
// contract DSDP
void DfTEI::contractDSDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPDS(KL, IJ);
    this->transform(5*3, 5, this->ERI);
}


////////////////////////////////////////////////////////////////////////
// contract DSDD
void DfTEI::contractDSDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDDS(KL, IJ);
    this->transform(5*5, 5, this->ERI);
}


////////////////////////////////////////////////////////////////////////
// contract DPSS
void DfTEI::contractDPSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSSS[6];
    for (int i = 0; i < 6; ++i) {
        DSSS[i] = 0.0;
    }
    double FSSS[10];
    for (int i = 0; i < 10; ++i) {
        FSSS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    //const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {
        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 4);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->WP = this->posW - this->posP;
            //this->WQ = this->posW - this->posQ;
            //this->QC = this->posQ - this->posC;

            this->primitivePSSS(*ij, *kl, 3);
            this->primitiveDSSS(*ij, *kl, 2);
            this->primitiveFSSS(*ij, *kl, 1);

            for (int i = 0; i < 6; ++i) {
                DSSS[i] += this->pDSSS[0][i];
            }
            for (int i = 0; i < 10; ++i) {
                FSSS[i] += this->pFSSS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();

    // HRR (DPSS)
    // [a(b+1i)|cd](m) = [(a+1i)b|cd](m) +(Ai-Bi)[ab|cd]
    double DPSS[18]; // 6*3
    DPSS[3*0 +0] = FSSS[ 0] +ABx*DSSS[ 0]; // i=x, a=xx b=0
    DPSS[3*0 +1] = FSSS[ 1] +ABy*DSSS[ 0]; // i=y, a=xx b=0
    DPSS[3*0 +2] = FSSS[ 2] +ABz*DSSS[ 0]; // i=z, a=xx b=0
    DPSS[3*1 +0] = FSSS[ 1] +ABx*DSSS[ 1]; // i=x, a=xy b=0
    DPSS[3*1 +1] = FSSS[ 3] +ABy*DSSS[ 1]; // i=y,   xy
    DPSS[3*1 +2] = FSSS[ 4] +ABz*DSSS[ 1]; // i=z,   xy
    DPSS[3*2 +0] = FSSS[ 2] +ABx*DSSS[ 2]; // i=x, a=xz
    DPSS[3*2 +1] = FSSS[ 4] +ABy*DSSS[ 2]; // i=y,   xz
    DPSS[3*2 +2] = FSSS[ 5] +ABz*DSSS[ 2]; // i=z,   xz
    DPSS[3*3 +0] = FSSS[ 3] +ABx*DSSS[ 3]; // i=x,   yy
    DPSS[3*3 +1] = FSSS[ 6] +ABy*DSSS[ 3]; // i=y,   yy
    DPSS[3*3 +2] = FSSS[ 7] +ABz*DSSS[ 3]; // i=z,   yy
    DPSS[3*4 +0] = FSSS[ 4] +ABx*DSSS[ 4]; // i=x,   yz
    DPSS[3*4 +1] = FSSS[ 7] +ABy*DSSS[ 4]; // i=y,   yz
    DPSS[3*4 +2] = FSSS[ 8] +ABz*DSSS[ 4]; // i=z,   yz
    DPSS[3*5 +0] = FSSS[ 5] +ABx*DSSS[ 5]; // i=x,   zz
    DPSS[3*5 +1] = FSSS[ 8] +ABy*DSSS[ 5]; // i=y,   zz
    DPSS[3*5 +2] = FSSS[ 9] +ABz*DSSS[ 5]; // i=z,   zz

    // convert 6D to 5D ==========
    for (int j = 0; j < 3; ++j) {
        this->ERI[3*0 +j] = DPSS[3*1 +j];                                                 // xy
        this->ERI[3*1 +j] = DPSS[3*2 +j];                                                 // xz
        this->ERI[3*2 +j] = DPSS[3*4 +j];                                                 // yz
        this->ERI[3*3 +j] = 0.5*(DPSS[3*0 +j] - DPSS[3*3 +j]);                            // xx - yy
        this->ERI[3*4 +j] = INV_SQRT3*(DPSS[3*5 +j] - 0.5*(DPSS[3*0 +j] + DPSS[3*3 +j])); // 2zz - (xx + yy)
    }
}


////////////////////////////////////////////////////////////////////////
// constract DPSP
void DfTEI::contractDPSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPPS(IJ, this->swapShellPair(KL));
    // no transform
}


////////////////////////////////////////////////////////////////////////
// constract DPSD
void DfTEI::contractDPSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPDS(IJ, this->swapShellPair(KL));
    // no transform
}


////////////////////////////////////////////////////////////////////////
// constract DPPS
void DfTEI::contractDPPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSPS[18]; // 6*3
    for (int i = 0; i < 18; ++i) {
        DSPS[i] = 0.0;
    }
    double FSPS[30]; // 10*3
    for (int i = 0; i < 30; ++i) {
        FSPS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 5);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 4);
            this->primitiveDSSS(*ij, *kl, 3);
            this->primitiveFSSS(*ij, *kl, 2);
            this->primitiveFSPS(*ij, *kl, 1);

            this->primitiveDSPS(*ij, *kl, 1);

            for (int i = 0; i < 18; ++i) {
                DSPS[i] += this->pDSPS[0][i];
            }
            for (int i = 0; i < 30; ++i) {
                FSPS[i] += this->pFSPS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();

    // HRR (DPPS)
    double DPPS[54]; // 6*3*3
    for (int k = 0; k < 3; ++k) {
        DPPS[9*0 +3*0 +k] = FSPS[3*0 +k] +ABx*DSPS[3*0 +k]; // i=x a=xx xxx=0
        DPPS[9*0 +3*1 +k] = FSPS[3*1 +k] +ABy*DSPS[3*0 +k]; // i=y a=xx xxy=1
        DPPS[9*0 +3*2 +k] = FSPS[3*2 +k] +ABz*DSPS[3*0 +k]; // i=z a=xx xxz=2

        DPPS[9*1 +3*0 +k] = FSPS[3*1 +k] +ABx*DSPS[3*1 +k]; // i=x a=xy
        DPPS[9*1 +3*1 +k] = FSPS[3*3 +k] +ABy*DSPS[3*1 +k]; // i=y a=xy xyy=3
        DPPS[9*1 +3*2 +k] = FSPS[3*4 +k] +ABz*DSPS[3*1 +k]; // i=z a=xy xyz=4

        DPPS[9*2 +3*0 +k] = FSPS[3*2 +k] +ABx*DSPS[3*2 +k]; // i=x a=xz
        DPPS[9*2 +3*1 +k] = FSPS[3*4 +k] +ABy*DSPS[3*2 +k]; // i=y a=xz
        DPPS[9*2 +3*2 +k] = FSPS[3*5 +k] +ABz*DSPS[3*2 +k]; // i=z a=xz xzz=5

        DPPS[9*3 +3*0 +k] = FSPS[3*3 +k] +ABx*DSPS[3*3 +k]; // i=x a=yy
        DPPS[9*3 +3*1 +k] = FSPS[3*6 +k] +ABy*DSPS[3*3 +k]; // i=y a=yy yyy=6
        DPPS[9*3 +3*2 +k] = FSPS[3*7 +k] +ABz*DSPS[3*3 +k]; // i=z a=yy yyz=7

        DPPS[9*4 +3*0 +k] = FSPS[3*4 +k] +ABx*DSPS[3*4 +k]; // i=x a=yz
        DPPS[9*4 +3*1 +k] = FSPS[3*7 +k] +ABy*DSPS[3*4 +k]; // i=y a=yz
        DPPS[9*4 +3*2 +k] = FSPS[3*8 +k] +ABz*DSPS[3*4 +k]; // i=z a=yz yzz=8

        DPPS[9*5 +3*0 +k] = FSPS[3*5 +k] +ABx*DSPS[3*5 +k]; // i=x a=zz
        DPPS[9*5 +3*1 +k] = FSPS[3*8 +k] +ABy*DSPS[3*5 +k]; // i=y a=zz
        DPPS[9*5 +3*2 +k] = FSPS[3*9 +k] +ABz*DSPS[3*5 +k]; // i=z a=zz zzz=9
    }

    // convert 6D to 5D ==========
    for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
            this->ERI[9*0 +3*j +k] = DPPS[9*1 +3*j +k];                                                         // xy
            this->ERI[9*1 +3*j +k] = DPPS[9*2 +3*j +k];                                                         // xz
            this->ERI[9*2 +3*j +k] = DPPS[9*4 +3*j +k];                                                         // yz
            this->ERI[9*3 +3*j +k] = 0.5*(DPPS[9*0 +3*j +k] -DPPS[9*3 +3*j +k]);                                // xx - yy
            this->ERI[9*4 +3*j +k] = INV_SQRT3*(DPPS[9*5 +3*j +k] -0.5*(DPPS[9*0 +3*j +k] +DPPS[9*3 +3*j +k])); // 2zz - (xx + yy)
        }
    }
}


////////////////////////////////////////////////////////////////////////
// contract DPPP
void DfTEI::contractDPPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSPS[18]; // 6*3
    for (int i = 0; i < 18; ++i) {
        DSPS[i] = 0.0;
    }
    double DSDS[36]; // 6*6
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }
    double FSPS[30]; // 10*3
    for (int i = 0; i < 30; ++i) {
        FSPS[i] = 0.0;
    }
    double FSDS[60]; // 10*6
    for (int i = 0; i < 60; ++i) {
        FSDS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 6);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 5);
            this->primitiveDSSS(*ij, *kl, 4);
            this->primitiveFSSS(*ij, *kl, 3);

            this->primitivePSPS(*ij, *kl, 2);

            this->primitiveDSPS(*ij, *kl, 2);
            this->primitiveDSDS(*ij, *kl, 1); // needs DSPS, PSPS

            this->primitiveFSPS(*ij, *kl, 2);
            this->primitiveFSDS(*ij, *kl, 1); // needs FSPS, DSPS

            for (int i = 0; i < 18; ++i) {
                DSPS[i] += this->pDSPS[0][i];
            }
            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
            for (int i = 0; i < 30; ++i) {
                FSPS[i] += this->pFSPS[0][i];
            }
            for (int i = 0; i < 60; ++i) {
                FSDS[i] += this->pFSDS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
    const double CDx = KL.AB.x();
    const double CDy = KL.AB.y();
    const double CDz = KL.AB.z();

    // HRR(DSPP)
    double DSPP[54]; //6*3*3
    for (int i = 0; i < 6; ++i) {
        DSPP[9*i +3*0 +0] = DSDS[6*i +0] + CDx*DSPS[3*i +0]; // -- x x
        DSPP[9*i +3*0 +1] = DSDS[6*i +1] + CDy*DSPS[3*i +0]; // -- x y
        DSPP[9*i +3*0 +2] = DSDS[6*i +2] + CDz*DSPS[3*i +0]; // -- x z
        DSPP[9*i +3*1 +0] = DSDS[6*i +1] + CDx*DSPS[3*i +1]; // -- y x
        DSPP[9*i +3*1 +1] = DSDS[6*i +3] + CDy*DSPS[3*i +1]; // -- y y
        DSPP[9*i +3*1 +2] = DSDS[6*i +4] + CDz*DSPS[3*i +1]; // -- y z
        DSPP[9*i +3*2 +0] = DSDS[6*i +2] + CDx*DSPS[3*i +2]; // -- z x
        DSPP[9*i +3*2 +1] = DSDS[6*i +4] + CDy*DSPS[3*i +2]; // -- z y
        DSPP[9*i +3*2 +2] = DSDS[6*i +5] + CDz*DSPS[3*i +2]; // -- z z
    }

    // HRR (FSPP)
    double FSPP[90]; // 10*3*3
    for (int i = 0; i < 10; ++i) {
        FSPP[9*i +3*0 +0] = FSDS[6*i +0] + CDx*FSPS[3*i +0]; // -- x x
        FSPP[9*i +3*0 +1] = FSDS[6*i +1] + CDy*FSPS[3*i +0]; // -- x y
        FSPP[9*i +3*0 +2] = FSDS[6*i +2] + CDz*FSPS[3*i +0]; // -- x z
        FSPP[9*i +3*1 +0] = FSDS[6*i +1] + CDx*FSPS[3*i +1]; // -- y x
        FSPP[9*i +3*1 +1] = FSDS[6*i +3] + CDy*FSPS[3*i +1]; // -- y y
        FSPP[9*i +3*1 +2] = FSDS[6*i +4] + CDz*FSPS[3*i +1]; // -- y z
        FSPP[9*i +3*2 +0] = FSDS[6*i +2] + CDx*FSPS[3*i +2]; // -- z x
        FSPP[9*i +3*2 +1] = FSDS[6*i +4] + CDy*FSPS[3*i +2]; // -- z y
        FSPP[9*i +3*2 +2] = FSDS[6*i +5] + CDz*FSPS[3*i +2]; // -- z z
    }

    // HRR (DPPP)
    double DPPP[162]; // 6*3*3*3
    for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
            DPPP[27*0 +9*0 +3*k +l] = FSPP[9*0 +3*k +l] +ABx*DSPP[9*0 +3*k +l]; // i=x, xx x
            DPPP[27*0 +9*1 +3*k +l] = FSPP[9*1 +3*k +l] +ABy*DSPP[9*0 +3*k +l]; // i=y, xx y
            DPPP[27*0 +9*2 +3*k +l] = FSPP[9*2 +3*k +l] +ABz*DSPP[9*0 +3*k +l]; // i=z, xx z
            DPPP[27*1 +9*0 +3*k +l] = FSPP[9*1 +3*k +l] +ABx*DSPP[9*1 +3*k +l]; // i=x, xy x
            DPPP[27*1 +9*1 +3*k +l] = FSPP[9*3 +3*k +l] +ABy*DSPP[9*1 +3*k +l]; // i=y, xy y
            DPPP[27*1 +9*2 +3*k +l] = FSPP[9*4 +3*k +l] +ABz*DSPP[9*1 +3*k +l]; // i=z, xy z
            DPPP[27*2 +9*0 +3*k +l] = FSPP[9*2 +3*k +l] +ABx*DSPP[9*2 +3*k +l]; // i=x, xz x
            DPPP[27*2 +9*1 +3*k +l] = FSPP[9*4 +3*k +l] +ABy*DSPP[9*2 +3*k +l]; // i=y, xz y
            DPPP[27*2 +9*2 +3*k +l] = FSPP[9*5 +3*k +l] +ABz*DSPP[9*2 +3*k +l]; // i=z, xz z
            DPPP[27*3 +9*0 +3*k +l] = FSPP[9*3 +3*k +l] +ABx*DSPP[9*3 +3*k +l]; // i=x, yy x
            DPPP[27*3 +9*1 +3*k +l] = FSPP[9*6 +3*k +l] +ABy*DSPP[9*3 +3*k +l]; // i=y, yy y
            DPPP[27*3 +9*2 +3*k +l] = FSPP[9*7 +3*k +l] +ABz*DSPP[9*3 +3*k +l]; // i=z, yy z
            DPPP[27*4 +9*0 +3*k +l] = FSPP[9*4 +3*k +l] +ABx*DSPP[9*4 +3*k +l]; // i=x, yz x
            DPPP[27*4 +9*1 +3*k +l] = FSPP[9*7 +3*k +l] +ABy*DSPP[9*4 +3*k +l]; // i=y, yz y
            DPPP[27*4 +9*2 +3*k +l] = FSPP[9*8 +3*k +l] +ABz*DSPP[9*4 +3*k +l]; // i=z, yz z
            DPPP[27*5 +9*0 +3*k +l] = FSPP[9*5 +3*k +l] +ABx*DSPP[9*5 +3*k +l]; // i=x, zz x
            DPPP[27*5 +9*1 +3*k +l] = FSPP[9*8 +3*k +l] +ABy*DSPP[9*5 +3*k +l]; // i=y, zz y
            DPPP[27*5 +9*2 +3*k +l] = FSPP[9*9 +3*k +l] +ABz*DSPP[9*5 +3*k +l]; // i=z, zz z
        }
    }

//   // HRR (DPPS)
//   double DPPS[54]; // 6*3*3
//   for (int k = 0; k < 3; ++k){
//     DPPS[9*0 +3*0 +k] = FSPS[3*0 +k] +ABx*DSPS[3*0 +k]; // i=x, a=xx
//     DPPS[9*0 +3*1 +k] = FSPS[3*1 +k] +ABy*DSPS[3*0 +k]; // i=y, a=xx
//     DPPS[9*0 +3*2 +k] = FSPS[3*2 +k] +ABz*DSPS[3*0 +k]; // i=z, a=xx
//     DPPS[9*1 +3*0 +k] = FSPS[3*1 +k] +ABx*DSPS[3*1 +k]; // i=x, a=xy
//     DPPS[9*1 +3*1 +k] = FSPS[3*3 +k] +ABy*DSPS[3*1 +k]; // i=y, a=xy
//     DPPS[9*1 +3*2 +k] = FSPS[3*4 +k] +ABz*DSPS[3*1 +k]; // i=z, a=xy
//     DPPS[9*2 +3*0 +k] = FSPS[3*2 +k] +ABx*DSPS[3*2 +k]; // i=x, a=xz
//     DPPS[9*2 +3*1 +k] = FSPS[3*4 +k] +ABy*DSPS[3*2 +k]; // i=y, a=xz
//     DPPS[9*2 +3*2 +k] = FSPS[3*5 +k] +ABz*DSPS[3*2 +k]; // i=z, a=xz
//     DPPS[9*3 +3*0 +k] = FSPS[3*3 +k] +ABx*DSPS[3*3 +k]; // i=x, a=yy
//     DPPS[9*3 +3*1 +k] = FSPS[3*6 +k] +ABy*DSPS[3*3 +k]; // i=y, a=yy
//     DPPS[9*3 +3*2 +k] = FSPS[3*7 +k] +ABz*DSPS[3*3 +k]; // i=z, a=yy
//     DPPS[9*4 +3*0 +k] = FSPS[3*4 +k] +ABx*DSPS[3*4 +k]; // i=x, a=yz
//     DPPS[9*4 +3*1 +k] = FSPS[3*7 +k] +ABy*DSPS[3*4 +k]; // i=y, a=yz
//     DPPS[9*4 +3*2 +k] = FSPS[3*8 +k] +ABz*DSPS[3*4 +k]; // i=z, a=yz
//     DPPS[9*5 +3*0 +k] = FSPS[3*5 +k] +ABx*DSPS[3*5 +k]; // i=x, a=zz
//     DPPS[9*5 +3*1 +k] = FSPS[3*8 +k] +ABy*DSPS[3*5 +k]; // i=y, a=zz
//     DPPS[9*5 +3*2 +k] = FSPS[3*9 +k] +ABz*DSPS[3*5 +k]; // i=z, a=zz
//   }

//   // HRR (DPDS)
//   double DPDS[108]; // 6*3*6
//   for (int k = 0; k < 6; ++k){
//     DPDS[18*0 +6*0 +k] = FSDS[6*0 +k] +ABx*DSDS[6*0 +k]; // i=x, a=xx
//     DPDS[18*0 +6*1 +k] = FSDS[6*1 +k] +ABy*DSDS[6*0 +k]; // i=y, a=xx
//     DPDS[18*0 +6*2 +k] = FSDS[6*2 +k] +ABz*DSDS[6*0 +k]; // i=z, a=xx
//     DPDS[18*1 +6*0 +k] = FSDS[6*1 +k] +ABx*DSDS[6*1 +k]; // i=x, a=xy
//     DPDS[18*1 +6*1 +k] = FSDS[6*3 +k] +ABy*DSDS[6*1 +k]; // i=y, a=xy
//     DPDS[18*1 +6*2 +k] = FSDS[6*4 +k] +ABz*DSDS[6*1 +k]; // i=z, a=xy
//     DPDS[18*2 +6*0 +k] = FSDS[6*2 +k] +ABx*DSDS[6*2 +k]; // i=x, xz x
//     DPDS[18*2 +6*1 +k] = FSDS[6*4 +k] +ABy*DSDS[6*2 +k]; // i=y, xz y
//     DPDS[18*2 +6*2 +k] = FSDS[6*5 +k] +ABz*DSDS[6*2 +k]; // i=z, xz z
//     DPDS[18*3 +6*0 +k] = FSDS[6*3 +k] +ABx*DSDS[6*3 +k]; // i=x, yy x
//     DPDS[18*3 +6*1 +k] = FSDS[6*6 +k] +ABy*DSDS[6*3 +k]; // i=y, yy y
//     DPDS[18*3 +6*2 +k] = FSDS[6*7 +k] +ABz*DSDS[6*3 +k]; // i=z, yy z
//     DPDS[18*4 +6*0 +k] = FSDS[6*4 +k] +ABx*DSDS[6*4 +k]; // i=x, yz x
//     DPDS[18*4 +6*1 +k] = FSDS[6*7 +k] +ABy*DSDS[6*4 +k]; // i=y, yz y
//     DPDS[18*4 +6*2 +k] = FSDS[6*8 +k] +ABz*DSDS[6*4 +k]; // i=z, yz z
//     DPDS[18*5 +6*0 +k] = FSDS[6*5 +k] +ABx*DSDS[6*5 +k]; // i=x, zz x
//     DPDS[18*5 +6*1 +k] = FSDS[6*8 +k] +ABy*DSDS[6*5 +k]; // i=y, zz y
//     DPDS[18*5 +6*2 +k] = FSDS[6*9 +k] +ABz*DSDS[6*5 +k]; // i=z, zz z
//   }

//   // HRR (DPPP)
//   double DPPP[162]; // 6*3*3*3
//   for (int i = 0; i < 6; ++i){
//     for (int j = 0; j < 3; ++j){
//       DPPP[27*i +9*j +3*0 +0] = DPDS[18*i +6*j +0] +CDx*DPPS[9*i +3*j +0]; //i=x, c=x
//       DPPP[27*i +9*j +3*0 +1] = DPDS[18*i +6*j +1] +CDy*DPPS[9*i +3*j +0]; //i=y, c=x
//       DPPP[27*i +9*j +3*0 +2] = DPDS[18*i +6*j +2] +CDz*DPPS[9*i +3*j +0]; //i=z, c=x
//       DPPP[27*i +9*j +3*1 +0] = DPDS[18*i +6*j +1] +CDx*DPPS[9*i +3*j +1]; //i=x, c=y
//       DPPP[27*i +9*j +3*1 +1] = DPDS[18*i +6*j +3] +CDy*DPPS[9*i +3*j +1]; //i=y, c=y
//       DPPP[27*i +9*j +3*1 +2] = DPDS[18*i +6*j +4] +CDz*DPPS[9*i +3*j +1]; //i=z, c=y
//       DPPP[27*i +9*j +3*2 +0] = DPDS[18*i +6*j +2] +CDx*DPPS[9*i +3*j +2]; //i=x, c=z
//       DPPP[27*i +9*j +3*2 +1] = DPDS[18*i +6*j +4] +CDy*DPPS[9*i +3*j +2]; //i=y, c=z
//       DPPP[27*i +9*j +3*2 +2] = DPDS[18*i +6*j +5] +CDz*DPPS[9*i +3*j +2]; //i=z, c=z
//     }
//   }


    // convert 6D to 5D ==========
    for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l) {
                this->ERI[27*0 +9*j +3*k +l] = DPPP[27*1 +9*j +3*k +l];                                          // xy
                this->ERI[27*1 +9*j +3*k +l] = DPPP[27*2 +9*j +3*k +l];                                          // xz
                this->ERI[27*2 +9*j +3*k +l] = DPPP[27*4 +9*j +3*k +l];                                          // yz
                this->ERI[27*3 +9*j +3*k +l] = 0.5*(DPPP[27*0 +9*j +3*k +l] - DPPP[27*3 +9*j +3*k +l]);          // xx - yy
                this->ERI[27*4 +9*j +3*k +l] = INV_SQRT3*(+DPPP[27*5 +9*j +3*k +l]
                                                          -0.5*(DPPP[27*0 +9*j +3*k +l] +DPPP[27*3 +9*j +3*k +l])); // 2zz - (xx + yy)
            }
        }
    }
}


////////////////////////////////////////////////////////////////////////
// constract DPPD
void DfTEI::contractDPPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDPDP(IJ, this->swapShellPair(KL));
    for (int i = 0; i < 5*3; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
}


////////////////////////////////////////////////////////////////////////
// constract DPDS
void DfTEI::contractDPDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSDS[36]; // 6*6
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }
    double FSDS[60]; // 6*10
    for (int i = 0; i < 60; ++i) {
        FSDS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 6);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 5);
            this->primitiveDSSS(*ij, *kl, 4);
            this->primitiveFSSS(*ij, *kl, 3);

            this->primitivePSPS(*ij, *kl, 2);

            this->primitiveDSPS(*ij, *kl, 2);
            this->primitiveDSDS(*ij, *kl, 1); // needs DSPS, PSPS

            this->primitiveFSPS(*ij, *kl, 2);
            this->primitiveFSDS(*ij, *kl, 1); // needs FSPS, DSPS

            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
            for (int i = 0; i < 60; ++i) {
                FSDS[i] += this->pFSDS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
//   const double CDx = KL.AB.x();
//   const double CDy = KL.AB.y();
//   const double CDz = KL.AB.z();

    // xxx xxy xxz xyy xyz xzz
    // yyy yyz yzz zzz
    // HRR (DPDS)
    double DPDS[108]; // 6*3*6
    for (int k = 0; k < 6; ++k) {
        DPDS[18*0 +6*0 +k] = FSDS[6*0 +k] +ABx*DSDS[6*0 +k]; // i=x, a=xx
        DPDS[18*0 +6*1 +k] = FSDS[6*1 +k] +ABy*DSDS[6*0 +k]; // i=y, a=xx
        DPDS[18*0 +6*2 +k] = FSDS[6*2 +k] +ABz*DSDS[6*0 +k]; // i=z, a=xx

        DPDS[18*1 +6*0 +k] = FSDS[6*1 +k] +ABx*DSDS[6*1 +k]; // i=x, a=xy
        DPDS[18*1 +6*1 +k] = FSDS[6*3 +k] +ABy*DSDS[6*1 +k]; // i=y, a=xy
        DPDS[18*1 +6*2 +k] = FSDS[6*4 +k] +ABz*DSDS[6*1 +k]; // i=z, a=xy

        DPDS[18*2 +6*0 +k] = FSDS[6*2 +k] +ABx*DSDS[6*2 +k]; // i=x, a=xz
        DPDS[18*2 +6*1 +k] = FSDS[6*4 +k] +ABy*DSDS[6*2 +k]; // i=y, a=xz
        DPDS[18*2 +6*2 +k] = FSDS[6*5 +k] +ABz*DSDS[6*2 +k]; // i=z, a=xz

        DPDS[18*3 +6*0 +k] = FSDS[6*3 +k] +ABx*DSDS[6*3 +k]; // i=x, yy x
        DPDS[18*3 +6*1 +k] = FSDS[6*6 +k] +ABy*DSDS[6*3 +k]; // i=y, yy y
        DPDS[18*3 +6*2 +k] = FSDS[6*7 +k] +ABz*DSDS[6*3 +k]; // i=z, yy z

        DPDS[18*4 +6*0 +k] = FSDS[6*4 +k] +ABx*DSDS[6*4 +k]; // i=x, yz x
        DPDS[18*4 +6*1 +k] = FSDS[6*7 +k] +ABy*DSDS[6*4 +k]; // i=y, yz y
        DPDS[18*4 +6*2 +k] = FSDS[6*8 +k] +ABz*DSDS[6*4 +k]; // i=z, yz z

        DPDS[18*5 +6*0 +k] = FSDS[6*5 +k] +ABx*DSDS[6*5 +k]; // i=x, zz x
        DPDS[18*5 +6*1 +k] = FSDS[6*8 +k] +ABy*DSDS[6*5 +k]; // i=y, zz y
        DPDS[18*5 +6*2 +k] = FSDS[6*9 +k] +ABz*DSDS[6*5 +k]; // i=z, zz z
    }

    // convert 6D to 5D
    double Dpds[90]; // 6*3*5
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 3; ++j) {
            Dpds[15*i +5*j +0] = DPDS[18*i +6*j +1]; // xy
            Dpds[15*i +5*j +1] = DPDS[18*i +6*j +2]; // xz
            Dpds[15*i +5*j +2] = DPDS[18*i +6*j +4]; // yz
            Dpds[15*i +5*j +3] = 0.5*(DPDS[18*i +6*j +0] -DPDS[18*i +6*j +3]); // xx-yy
            Dpds[15*i +5*j +4] = INV_SQRT3*(DPDS[18*i +6*j +5] -0.5*(DPDS[18*i +6*j +0] +DPDS[18*i +6*j +3])); // 2zz-(xx+yy)
        }
    }

    for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 5; ++k) {
            this->ERI[15*0 +5*j +k] = Dpds[15*1 +5*j +k]; // xy
            this->ERI[15*1 +5*j +k] = Dpds[15*2 +5*j +k]; // xz
            this->ERI[15*2 +5*j +k] = Dpds[15*4 +5*j +k]; // yz
            this->ERI[15*3 +5*j +k] = 0.5*(Dpds[15*0 +5*j +k] -Dpds[15*3 +5*j +k]); // xx-yy
            this->ERI[15*4 +5*j +k] = INV_SQRT3*(Dpds[15*5 +5*j +k] -0.5*(Dpds[15*0 +5*j +k] +Dpds[15*3 +5*j +k])); // 2zz-(xx+yy)
        }
    }
}


////////////////////////////////////////////////////////////////////////
// constract DPDP
void DfTEI::contractDPDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSDS[36]; // 6*6
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }
    double DSFS[60]; // 6*10
    for (int i = 0; i < 60; ++i) {
        DSFS[i] = 0.0;
    }
    double FSDS[60]; // 6*10
    for (int i = 0; i < 60; ++i) {
        FSDS[i] = 0.0;
    }
    double FSFS[100]; // 10*10
    for (int i = 0; i < 100; ++i) {
        FSFS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 7);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 6);
            this->primitiveDSSS(*ij, *kl, 5);
            this->primitiveFSSS(*ij, *kl, 4); //

            this->primitiveSSPS(*ij, *kl, 3);

            this->primitivePSPS(*ij, *kl, 3);
            this->primitivePSDS(*ij, *kl, 2); // needs SSPS

            this->primitiveDSPS(*ij, *kl, 3);
            this->primitiveDSDS(*ij, *kl, 2); // needs DSPS, PSPS
            this->primitiveDSFS(*ij, *kl, 1); // needs DSDS, DSPS, PSDS

            this->primitiveFSPS(*ij, *kl, 3); //
            this->primitiveFSDS(*ij, *kl, 2); // needs FSPS, DSPS
            this->primitiveFSFS(*ij, *kl, 1); // needs FSDS, FSPS, DSDS

            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
            for (int i = 0; i < 60; ++i) {
                DSFS[i] += this->pDSFS[0][i];
            }
            for (int i = 0; i < 60; ++i) {
                FSDS[i] += this->pFSDS[0][i];
            }
            for (int i = 0; i < 100; ++i) {
                FSFS[i] += this->pFSFS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
    const double CDx = KL.AB.x();
    const double CDy = KL.AB.y();
    const double CDz = KL.AB.z();

    // HRR (DPFS)
    double DPFS[180]; //6*3*10
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSFS, DSFS, 10, 1, DPFS);

    // HRR (DPDS)
    double DPDS[108]; // 6*3*6
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSDS, DSDS, 6, 1, DPDS);

    // HRR (DPDP)
    double DPDP[324]; // 6*3*6*3
    this->HRR_XXDP(CDx, CDy, CDz,
                   DPFS, DPDS, 6, 3, DPDP);





    // convert 6D to 5D
    double Dpdp[270]; // 6*3*5*3
    this->convert6Dto5D_K(6, 3, 3, DPDP, Dpdp);

    this->convert6Dto5D_I(3, 5, 3, Dpdp, this->ERI);
}


////////////////////////////////////////////////////////////////////////
// constract DPDD
void DfTEI::contractDPDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDDP(KL, IJ);
    this->transform(5*5, 5*3, this->ERI);
}


////////////////////////////////////////////////////////////////////////
// constract DDSS
void DfTEI::contractDDSS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSSS[6]; // 6
    for (int i = 0; i < 6; ++i) {
        DSSS[i] = 0.0;
    }
    double FSSS[10]; //10
    for (int i = 0; i < 10; ++i) {
        FSSS[i] = 0.0;
    }
    double GSSS[15]; //15
    for (int i = 0; i < 15; ++i) {
        GSSS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    //const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 5);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            //this->QC = this->posQ - this->posC;
            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->WP = this->posW - this->posP;
            //this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 4);
            this->primitiveDSSS(*ij, *kl, 3);
            this->primitiveFSSS(*ij, *kl, 2);
            this->primitiveGSSS(*ij, *kl, 1);

            for (int i = 0; i < 6; ++i) {
                DSSS[i] += this->pDSSS[0][i];
            }
            for (int i = 0; i < 10; ++i) {
                FSSS[i] += this->pFSSS[0][i];
            }
            for (int i = 0; i < 15; ++i) {
                GSSS[i] += this->pGSSS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
//   const double CDx = KL.AB.x();
//   const double CDy = KL.AB.y();
//   const double CDz = KL.AB.z();

    // HRR (DPSS)
    double DPSS[18]; // 6*3
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSSS, DSSS, 1, 1, DPSS);

    // HRR (FPSS)
    double FPSS[30]; // 10*3
    this->HRR_FPXX(ABx, ABy, ABz,
                   GSSS, FSSS, 1, 1, FPSS);

    // HRR (DDSS)
    double DDSS[36]; // 6*6
    this->HRR_DDXX(ABx, ABy, ABz,
                   FPSS, DPSS, 1, 1, DDSS);

    // convert 6D to 5D
    double Ddss[30]; // 6*5
    for (int i = 0; i < 6; ++i) {
        Ddss[5*i +0] = DDSS[6*i +1]; // xy
        Ddss[5*i +1] = DDSS[6*i +2]; // xz
        Ddss[5*i +2] = DDSS[6*i +4]; // yz
        Ddss[5*i +3] = 0.5*(DDSS[6*i +0] -DDSS[6*i +3]); // xx-yy
        Ddss[5*i +4] = INV_SQRT3*(DDSS[6*i +5] -0.5*(DDSS[6*i +0] +DDSS[6*i +3])); // 2zz-(xx+yy)
    }

    for (int j = 0; j < 5; ++j) {
        this->ERI[5*0 +j] = Ddss[5*1 +j]; // xy
        this->ERI[5*1 +j] = Ddss[5*2 +j]; // xz
        this->ERI[5*2 +j] = Ddss[5*4 +j]; // yz
        this->ERI[5*3 +j] = 0.5*(Ddss[5*0 +j] -Ddss[5*3 +j]); // xx-yy
        this->ERI[5*4 +j] = INV_SQRT3*(Ddss[5*5 +j] -0.5*(Ddss[5*0 +j] +Ddss[5*3 +j])); // 2zz-(xx+yy)
    }
}


////////////////////////////////////////////////////////////////////////
// constract DDSP
void DfTEI::contractDDSP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDPS(IJ, this->swapShellPair(KL));
    // no transform
}


////////////////////////////////////////////////////////////////////////
// constract DDSD
void DfTEI::contractDDSD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDDS(IJ, this->swapShellPair(KL));
    // no transform
}


////////////////////////////////////////////////////////////////////////
// constract DDPS
void DfTEI::contractDDPS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSPS[18]; // 6*3
    for (int i = 0; i < 18; ++i) {
        DSPS[i] = 0.0;
    }
    double FSPS[30]; //10*3
    for (int i = 0; i < 30; ++i) {
        FSPS[i] = 0.0;
    }
    double GSPS[45]; //15*3
    for (int i = 0; i < 45; ++i) {
        GSPS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 6);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 5);
            this->primitiveDSSS(*ij, *kl, 4);
            this->primitiveFSSS(*ij, *kl, 3);
            this->primitiveGSSS(*ij, *kl, 2);

            this->primitiveDSPS(*ij, *kl, 1);
            this->primitiveFSPS(*ij, *kl, 1);
            this->primitiveGSPS(*ij, *kl, 1);

            for (int i = 0; i < 18; ++i) {
                DSPS[i] += this->pDSPS[0][i];
            }
            for (int i = 0; i < 30; ++i) {
                FSPS[i] += this->pFSPS[0][i];
            }
            for (int i = 0; i < 45; ++i) {
                GSPS[i] += this->pGSPS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
//   const double CDx = KL.AB.x();
//   const double CDy = KL.AB.y();
//   const double CDz = KL.AB.z();

    // HRR (DPPS)
    double DPPS[54]; // 6*3*3
    for (int k = 0; k < 3; ++k) {
        DPPS[9*0 +3*0 +k] = FSPS[3*0 +k] +ABx*DSPS[3*0 +k]; // i=x, xx x
        DPPS[9*0 +3*1 +k] = FSPS[3*1 +k] +ABy*DSPS[3*0 +k]; // i=y, xx y
        DPPS[9*0 +3*2 +k] = FSPS[3*2 +k] +ABz*DSPS[3*0 +k]; // i=z, xx z
        DPPS[9*1 +3*0 +k] = FSPS[3*1 +k] +ABx*DSPS[3*1 +k]; // i=x, xy x
        DPPS[9*1 +3*1 +k] = FSPS[3*3 +k] +ABy*DSPS[3*1 +k]; // i=y, xy y
        DPPS[9*1 +3*2 +k] = FSPS[3*4 +k] +ABz*DSPS[3*1 +k]; // i=z, xy z
        DPPS[9*2 +3*0 +k] = FSPS[3*2 +k] +ABx*DSPS[3*2 +k]; // i=x, xz x
        DPPS[9*2 +3*1 +k] = FSPS[3*4 +k] +ABy*DSPS[3*2 +k]; // i=y, xz y
        DPPS[9*2 +3*2 +k] = FSPS[3*5 +k] +ABz*DSPS[3*2 +k]; // i=z, xz z
        DPPS[9*3 +3*0 +k] = FSPS[3*3 +k] +ABx*DSPS[3*3 +k]; // i=x, yy x
        DPPS[9*3 +3*1 +k] = FSPS[3*6 +k] +ABy*DSPS[3*3 +k]; // i=y, yy y
        DPPS[9*3 +3*2 +k] = FSPS[3*7 +k] +ABz*DSPS[3*3 +k]; // i=z, yy z
        DPPS[9*4 +3*0 +k] = FSPS[3*4 +k] +ABx*DSPS[3*4 +k]; // i=x, yz x
        DPPS[9*4 +3*1 +k] = FSPS[3*7 +k] +ABy*DSPS[3*4 +k]; // i=y, yz y
        DPPS[9*4 +3*2 +k] = FSPS[3*8 +k] +ABz*DSPS[3*4 +k]; // i=z, yz z
        DPPS[9*5 +3*0 +k] = FSPS[3*5 +k] +ABx*DSPS[3*5 +k]; // i=x, zz x
        DPPS[9*5 +3*1 +k] = FSPS[3*8 +k] +ABy*DSPS[3*5 +k]; // i=y, zz y
        DPPS[9*5 +3*2 +k] = FSPS[3*9 +k] +ABz*DSPS[3*5 +k]; // i=z, zz z
    }

    // HRR (FPPS)
    double FPPS[90]; // 10*3*3
    for (int k = 0; k < 3; ++k) {
        FPPS[9*0 +3*0 +k] = GSPS[3* 0 +k] + ABx*FSPS[3*0 +k]; // i=x, xxx x
        FPPS[9*0 +3*1 +k] = GSPS[3* 1 +k] + ABy*FSPS[3*0 +k]; // i=y, xxx y
        FPPS[9*0 +3*2 +k] = GSPS[3* 2 +k] + ABz*FSPS[3*0 +k]; // i=z, xxx z
        FPPS[9*1 +3*0 +k] = GSPS[3* 1 +k] + ABx*FSPS[3*1 +k]; // i=x, xxy x
        FPPS[9*1 +3*1 +k] = GSPS[3* 3 +k] + ABy*FSPS[3*1 +k]; // i=y, xxy y
        FPPS[9*1 +3*2 +k] = GSPS[3* 4 +k] + ABz*FSPS[3*1 +k]; // i=z, xxy z
        FPPS[9*2 +3*0 +k] = GSPS[3* 2 +k] + ABx*FSPS[3*2 +k]; // i=x, xxz x
        FPPS[9*2 +3*1 +k] = GSPS[3* 4 +k] + ABy*FSPS[3*2 +k]; // i=y, xxz y
        FPPS[9*2 +3*2 +k] = GSPS[3* 5 +k] + ABz*FSPS[3*2 +k]; // i=z, xxz z
        FPPS[9*3 +3*0 +k] = GSPS[3* 3 +k] + ABx*FSPS[3*3 +k]; // i=x, xyy x
        FPPS[9*3 +3*1 +k] = GSPS[3* 6 +k] + ABy*FSPS[3*3 +k]; // i=y, xyy y
        FPPS[9*3 +3*2 +k] = GSPS[3* 7 +k] + ABz*FSPS[3*3 +k]; // i=z, xyy z
        FPPS[9*4 +3*0 +k] = GSPS[3* 4 +k] + ABx*FSPS[3*4 +k]; // i=x, xyz x
        FPPS[9*4 +3*1 +k] = GSPS[3* 7 +k] + ABy*FSPS[3*4 +k]; // i=y, xyz y
        FPPS[9*4 +3*2 +k] = GSPS[3* 8 +k] + ABz*FSPS[3*4 +k]; // i=z, xyz z
        FPPS[9*5 +3*0 +k] = GSPS[3* 5 +k] + ABx*FSPS[3*5 +k]; // i=x, xzz x
        FPPS[9*5 +3*1 +k] = GSPS[3* 8 +k] + ABy*FSPS[3*5 +k]; // i=y, xzz y
        FPPS[9*5 +3*2 +k] = GSPS[3* 9 +k] + ABz*FSPS[3*5 +k]; // i=z, xzz z
        FPPS[9*6 +3*0 +k] = GSPS[3* 6 +k] + ABx*FSPS[3*6 +k]; // i=x, yyy x
        FPPS[9*6 +3*1 +k] = GSPS[3*10 +k] + ABy*FSPS[3*6 +k]; // i=y, yyy y
        FPPS[9*6 +3*2 +k] = GSPS[3*11 +k] + ABz*FSPS[3*6 +k]; // i=z, yyy z
        FPPS[9*7 +3*0 +k] = GSPS[3* 7 +k] + ABx*FSPS[3*7 +k]; // i=x, yyz x
        FPPS[9*7 +3*1 +k] = GSPS[3*11 +k] + ABy*FSPS[3*7 +k]; // i=y, yyz y
        FPPS[9*7 +3*2 +k] = GSPS[3*12 +k] + ABz*FSPS[3*7 +k]; // i=z, yyz z
        FPPS[9*8 +3*0 +k] = GSPS[3* 8 +k] + ABx*FSPS[3*8 +k]; // i=x, yzz x
        FPPS[9*8 +3*1 +k] = GSPS[3*12 +k] + ABy*FSPS[3*8 +k]; // i=y, yzz y
        FPPS[9*8 +3*2 +k] = GSPS[3*13 +k] + ABz*FSPS[3*8 +k]; // i=z, yzz z
        FPPS[9*9 +3*0 +k] = GSPS[3* 9 +k] + ABx*FSPS[3*9 +k]; // i=x, zzz x
        FPPS[9*9 +3*1 +k] = GSPS[3*13 +k] + ABy*FSPS[3*9 +k]; // i=y, zzz y
        FPPS[9*9 +3*2 +k] = GSPS[3*14 +k] + ABz*FSPS[3*9 +k]; // i=z, zzz z
    }

    // HRR (DDPS)
    double DDPS[108]; // 6*6*3
    this->HRR_DDXX(ABx, ABy, ABz,
                   FPPS, DPPS, 3, 1, DDPS);

    // convert 6D to 5D
    double Ddps[90]; // 6*5*3
    for (int i = 0; i < 6; ++i) {
        for (int l = 0; l < 3; ++l) {
            Ddps[15*i +3*0 +l] = DDPS[18*i +3*1 +l]; // xy
            Ddps[15*i +3*1 +l] = DDPS[18*i +3*2 +l]; // xz
            Ddps[15*i +3*2 +l] = DDPS[18*i +3*4 +l]; // yz
            Ddps[15*i +3*3 +l] = 0.5*(DDPS[18*i +3*0 +l] -DDPS[18*i +3*3 +l]); // xx-yy
            Ddps[15*i +3*4 +l] = INV_SQRT3*(DDPS[18*i +3*5 +l] -0.5*(DDPS[18*i +3*0 +l] +DDPS[18*i +3*3 +l])); // 2zz-(xx+yy)
        }
    }

    for (int j = 0; j < 5; ++j) {
        for (int l = 0; l < 3; ++l) {
            this->ERI[15*0 +3*j +l] = Ddps[15*1 +3*j +l]; // xy
            this->ERI[15*1 +3*j +l] = Ddps[15*2 +3*j +l]; // xz
            this->ERI[15*2 +3*j +l] = Ddps[15*4 +3*j +l]; // yz
            this->ERI[15*3 +3*j +l] = 0.5*(Ddps[15*0 +3*j +l] -Ddps[15*3 +3*j +l]); // xx-yy
            this->ERI[15*4 +3*j +l] = INV_SQRT3*(Ddps[15*5 +3*j +l] -0.5*(Ddps[15*0 +3*j +l] +Ddps[15*3 +3*j +l])); // 2zz-(xx+yy)
        }
    }
}


////////////////////////////////////////////////////////////////////////
// constract DDPP
void DfTEI::contractDDPP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSPS[18]; // 6*3
    for (int i = 0; i < 18; ++i) {
        DSPS[i] = 0.0;
    }
    double DSDS[36]; // 6*6
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }
    double FSPS[30]; //10*3
    for (int i = 0; i < 30; ++i) {
        FSPS[i] = 0.0;
    }
    double FSDS[60]; //10*6
    for (int i = 0; i < 60; ++i) {
        FSDS[i] = 0.0;
    }
    double GSPS[45]; //15*3
    for (int i = 0; i < 45; ++i) {
        GSPS[i] = 0.0;
    }
    double GSDS[90]; //15*6
    for (int i = 0; i < 90; ++i) {
        GSDS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 7);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 6);
            this->primitiveDSSS(*ij, *kl, 5);
            this->primitiveFSSS(*ij, *kl, 4);
            this->primitiveGSSS(*ij, *kl, 3);

            this->primitivePSPS(*ij, *kl, 2);

            this->primitiveDSPS(*ij, *kl, 2);
            this->primitiveDSDS(*ij, *kl, 1); // needs DSPS, PSPS

            this->primitiveFSPS(*ij, *kl, 2);
            this->primitiveFSDS(*ij, *kl, 1); // needs FSPS, DSPS

            this->primitiveGSPS(*ij, *kl, 2);
            this->primitiveGSDS(*ij, *kl, 1); // needs GSPS, FSPS

            for (int i = 0; i < 18; ++i) {
                DSPS[i] += this->pDSPS[0][i];
            }
            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
            for (int i = 0; i < 30; ++i) {
                FSPS[i] += this->pFSPS[0][i];
            }
            for (int i = 0; i < 60; ++i) {
                FSDS[i] += this->pFSDS[0][i];
            }
            for (int i = 0; i < 45; ++i) {
                GSPS[i] += this->pGSPS[0][i];
            }
            for (int i = 0; i < 90; ++i) {
                GSDS[i] += this->pGSDS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
    const double CDx = KL.AB.x();
    const double CDy = KL.AB.y();
    const double CDz = KL.AB.z();

    // HRR(DSPP)
    double DSPP[54]; //6*3*3
    this->HRR_XXPP(CDx, CDy, CDz,
                   DSDS, DSPS, 6, 1, DSPP);
//   for (int i = 0; i < 6; ++i){
//     DSPP[9*i +3*0 +0] = DSDS[6*i +0] + CDx*DSPS[3*i +0]; // -- x x
//     DSPP[9*i +3*0 +1] = DSDS[6*i +1] + CDy*DSPS[3*i +0]; // -- x y
//     DSPP[9*i +3*0 +2] = DSDS[6*i +2] + CDz*DSPS[3*i +0]; // -- x z
//     DSPP[9*i +3*1 +0] = DSDS[6*i +1] + CDx*DSPS[3*i +1]; // -- y x
//     DSPP[9*i +3*1 +1] = DSDS[6*i +3] + CDy*DSPS[3*i +1]; // -- y y
//     DSPP[9*i +3*1 +2] = DSDS[6*i +4] + CDz*DSPS[3*i +1]; // -- y z
//     DSPP[9*i +3*2 +0] = DSDS[6*i +2] + CDx*DSPS[3*i +2]; // -- z x
//     DSPP[9*i +3*2 +1] = DSDS[6*i +4] + CDy*DSPS[3*i +2]; // -- z y
//     DSPP[9*i +3*2 +2] = DSDS[6*i +5] + CDz*DSPS[3*i +2]; // -- z z
//   }

    // HRR (FSPP)
    double FSPP[90]; // 10*3*3
    this->HRR_XXPP(CDx, CDy, CDz,
                   FSDS, FSPS, 10, 1, FSPP);
//   for (int i = 0; i < 10; ++i){
//     FSPP[9*i +3*0 +0] = FSDS[6*i +0] + CDx*FSPS[3*i +0]; // -- x x
//     FSPP[9*i +3*0 +1] = FSDS[6*i +1] + CDy*FSPS[3*i +0]; // -- x y
//     FSPP[9*i +3*0 +2] = FSDS[6*i +2] + CDz*FSPS[3*i +0]; // -- x z
//     FSPP[9*i +3*1 +0] = FSDS[6*i +1] + CDx*FSPS[3*i +1]; // -- y x
//     FSPP[9*i +3*1 +1] = FSDS[6*i +3] + CDy*FSPS[3*i +1]; // -- y y
//     FSPP[9*i +3*1 +2] = FSDS[6*i +4] + CDz*FSPS[3*i +1]; // -- y z
//     FSPP[9*i +3*2 +0] = FSDS[6*i +2] + CDx*FSPS[3*i +2]; // -- z x
//     FSPP[9*i +3*2 +1] = FSDS[6*i +4] + CDy*FSPS[3*i +2]; // -- z y
//     FSPP[9*i +3*2 +2] = FSDS[6*i +5] + CDz*FSPS[3*i +2]; // -- z z
//   }

    // HRR (GSPP)
    double GSPP[135]; // 15*3*3
    this->HRR_XXPP(CDx, CDy, CDz,
                   GSDS, GSPS, 15, 1, GSPP);
//   for (int i = 0; i < 15; ++i){
//     GSPP[9*i +3*0 +0] = GSDS[6*i +0] + CDx*GSPS[3*i +0]; // -- x x
//     GSPP[9*i +3*0 +1] = GSDS[6*i +1] + CDy*GSPS[3*i +0]; // -- x y
//     GSPP[9*i +3*0 +2] = GSDS[6*i +2] + CDz*GSPS[3*i +0]; // -- x z
//     GSPP[9*i +3*1 +0] = GSDS[6*i +1] + CDx*GSPS[3*i +1]; // -- y x
//     GSPP[9*i +3*1 +1] = GSDS[6*i +3] + CDy*GSPS[3*i +1]; // -- y y
//     GSPP[9*i +3*1 +2] = GSDS[6*i +4] + CDz*GSPS[3*i +1]; // -- y z
//     GSPP[9*i +3*2 +0] = GSDS[6*i +2] + CDx*GSPS[3*i +2]; // -- z x
//     GSPP[9*i +3*2 +1] = GSDS[6*i +4] + CDy*GSPS[3*i +2]; // -- z y
//     GSPP[9*i +3*2 +2] = GSDS[6*i +5] + CDz*GSPS[3*i +2]; // -- z z
//   }

    // HRR (DPPP)
    double DPPP[162]; // 6*3*3*3
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSPP, DSPP, 3, 3, DPPP);
//   for (int k = 0; k < 3; ++k){
//     for (int l = 0; l < 3; ++l){
//       DPPP[27*0 +9*0 +3*k +l] = FSPP[9*0 +3*k +l] +ABx*DSPP[9*0 +3*k +l]; // i=x, xx x
//       DPPP[27*0 +9*1 +3*k +l] = FSPP[9*1 +3*k +l] +ABy*DSPP[9*0 +3*k +l]; // i=y, xx y
//       DPPP[27*0 +9*2 +3*k +l] = FSPP[9*2 +3*k +l] +ABz*DSPP[9*0 +3*k +l]; // i=z, xx z
//       DPPP[27*1 +9*0 +3*k +l] = FSPP[9*1 +3*k +l] +ABx*DSPP[9*1 +3*k +l]; // i=x, xy x
//       DPPP[27*1 +9*1 +3*k +l] = FSPP[9*3 +3*k +l] +ABy*DSPP[9*1 +3*k +l]; // i=y, xy y
//       DPPP[27*1 +9*2 +3*k +l] = FSPP[9*4 +3*k +l] +ABz*DSPP[9*1 +3*k +l]; // i=z, xy z
//       DPPP[27*2 +9*0 +3*k +l] = FSPP[9*2 +3*k +l] +ABx*DSPP[9*2 +3*k +l]; // i=x, xz x
//       DPPP[27*2 +9*1 +3*k +l] = FSPP[9*4 +3*k +l] +ABy*DSPP[9*2 +3*k +l]; // i=y, xz y
//       DPPP[27*2 +9*2 +3*k +l] = FSPP[9*5 +3*k +l] +ABz*DSPP[9*2 +3*k +l]; // i=z, xz z
//       DPPP[27*3 +9*0 +3*k +l] = FSPP[9*3 +3*k +l] +ABx*DSPP[9*3 +3*k +l]; // i=x, yy x
//       DPPP[27*3 +9*1 +3*k +l] = FSPP[9*6 +3*k +l] +ABy*DSPP[9*3 +3*k +l]; // i=y, yy y
//       DPPP[27*3 +9*2 +3*k +l] = FSPP[9*7 +3*k +l] +ABz*DSPP[9*3 +3*k +l]; // i=z, yy z
//       DPPP[27*4 +9*0 +3*k +l] = FSPP[9*4 +3*k +l] +ABx*DSPP[9*4 +3*k +l]; // i=x, yz x
//       DPPP[27*4 +9*1 +3*k +l] = FSPP[9*7 +3*k +l] +ABy*DSPP[9*4 +3*k +l]; // i=y, yz y
//       DPPP[27*4 +9*2 +3*k +l] = FSPP[9*8 +3*k +l] +ABz*DSPP[9*4 +3*k +l]; // i=z, yz z
//       DPPP[27*5 +9*0 +3*k +l] = FSPP[9*5 +3*k +l] +ABx*DSPP[9*5 +3*k +l]; // i=x, zz x
//       DPPP[27*5 +9*1 +3*k +l] = FSPP[9*8 +3*k +l] +ABy*DSPP[9*5 +3*k +l]; // i=y, zz y
//       DPPP[27*5 +9*2 +3*k +l] = FSPP[9*9 +3*k +l] +ABz*DSPP[9*5 +3*k +l]; // i=z, zz z
//     }
//   }

    // HRR (FPPP)
    double FPPP[270]; // 10*3*3*3
    this->HRR_FPXX(ABx, ABy, ABz,
                   GSPP, FSPP, 3, 3, FPPP);
//   for (int k = 0; k < 3; ++k){
//     for (int l = 0; l < 3; ++l){
//       FPPP[27*0 +9*0 +3*k +l] = GSPP[9* 0 +3*k +l] + ABx*FSPP[9*0 +3*k +l]; // i=x, xxx x
//       FPPP[27*0 +9*1 +3*k +l] = GSPP[9* 1 +3*k +l] + ABy*FSPP[9*0 +3*k +l]; // i=y, xxx y
//       FPPP[27*0 +9*2 +3*k +l] = GSPP[9* 2 +3*k +l] + ABz*FSPP[9*0 +3*k +l]; // i=z, xxx z
//       FPPP[27*1 +9*0 +3*k +l] = GSPP[9* 1 +3*k +l] + ABx*FSPP[9*1 +3*k +l]; // i=x, xxy x
//       FPPP[27*1 +9*1 +3*k +l] = GSPP[9* 3 +3*k +l] + ABy*FSPP[9*1 +3*k +l]; // i=y, xxy y
//       FPPP[27*1 +9*2 +3*k +l] = GSPP[9* 4 +3*k +l] + ABz*FSPP[9*1 +3*k +l]; // i=z, xxy z
//       FPPP[27*2 +9*0 +3*k +l] = GSPP[9* 2 +3*k +l] + ABx*FSPP[9*2 +3*k +l]; // i=x, xxz x
//       FPPP[27*2 +9*1 +3*k +l] = GSPP[9* 4 +3*k +l] + ABy*FSPP[9*2 +3*k +l]; // i=y, xxz y
//       FPPP[27*2 +9*2 +3*k +l] = GSPP[9* 5 +3*k +l] + ABz*FSPP[9*2 +3*k +l]; // i=z, xxz z
//       FPPP[27*3 +9*0 +3*k +l] = GSPP[9* 3 +3*k +l] + ABx*FSPP[9*3 +3*k +l]; // i=x, xyy x
//       FPPP[27*3 +9*1 +3*k +l] = GSPP[9* 6 +3*k +l] + ABy*FSPP[9*3 +3*k +l]; // i=y, xyy y
//       FPPP[27*3 +9*2 +3*k +l] = GSPP[9* 7 +3*k +l] + ABz*FSPP[9*3 +3*k +l]; // i=z, xyy z
//       FPPP[27*4 +9*0 +3*k +l] = GSPP[9* 4 +3*k +l] + ABx*FSPP[9*4 +3*k +l]; // i=x, xyz x
//       FPPP[27*4 +9*1 +3*k +l] = GSPP[9* 7 +3*k +l] + ABy*FSPP[9*4 +3*k +l]; // i=y, xyz y
//       FPPP[27*4 +9*2 +3*k +l] = GSPP[9* 8 +3*k +l] + ABz*FSPP[9*4 +3*k +l]; // i=z, xyz z
//       FPPP[27*5 +9*0 +3*k +l] = GSPP[9* 5 +3*k +l] + ABx*FSPP[9*5 +3*k +l]; // i=x, xzz x
//       FPPP[27*5 +9*1 +3*k +l] = GSPP[9* 8 +3*k +l] + ABy*FSPP[9*5 +3*k +l]; // i=y, xzz y
//       FPPP[27*5 +9*2 +3*k +l] = GSPP[9* 9 +3*k +l] + ABz*FSPP[9*5 +3*k +l]; // i=z, xzz z
//       FPPP[27*6 +9*0 +3*k +l] = GSPP[9* 6 +3*k +l] + ABx*FSPP[9*6 +3*k +l]; // i=x, yyy x
//       FPPP[27*6 +9*1 +3*k +l] = GSPP[9*10 +3*k +l] + ABy*FSPP[9*6 +3*k +l]; // i=y, yyy y
//       FPPP[27*6 +9*2 +3*k +l] = GSPP[9*11 +3*k +l] + ABz*FSPP[9*6 +3*k +l]; // i=z, yyy z
//       FPPP[27*7 +9*0 +3*k +l] = GSPP[9* 7 +3*k +l] + ABx*FSPP[9*7 +3*k +l]; // i=x, yyz x
//       FPPP[27*7 +9*1 +3*k +l] = GSPP[9*11 +3*k +l] + ABy*FSPP[9*7 +3*k +l]; // i=y, yyz y
//       FPPP[27*7 +9*2 +3*k +l] = GSPP[9*12 +3*k +l] + ABz*FSPP[9*7 +3*k +l]; // i=z, yyz z
//       FPPP[27*8 +9*0 +3*k +l] = GSPP[9* 8 +3*k +l] + ABx*FSPP[9*8 +3*k +l]; // i=x, yzz x
//       FPPP[27*8 +9*1 +3*k +l] = GSPP[9*12 +3*k +l] + ABy*FSPP[9*8 +3*k +l]; // i=y, yzz y
//       FPPP[27*8 +9*2 +3*k +l] = GSPP[9*13 +3*k +l] + ABz*FSPP[9*8 +3*k +l]; // i=z, yzz z
//       FPPP[27*9 +9*0 +3*k +l] = GSPP[9* 9 +3*k +l] + ABx*FSPP[9*9 +3*k +l]; // i=x, zzz x
//       FPPP[27*9 +9*1 +3*k +l] = GSPP[9*13 +3*k +l] + ABy*FSPP[9*9 +3*k +l]; // i=y, zzz y
//       FPPP[27*9 +9*2 +3*k +l] = GSPP[9*14 +3*k +l] + ABz*FSPP[9*9 +3*k +l]; // i=z, zzz z
//     }
//   }

    // HRR (DDPP)
    double DDPP[324]; // 6*6*3*3
    this->HRR_DDXX(ABx, ABy, ABz,
                   FPPP, DPPP, 3, 3, DDPP);
//   for (int k = 0; k < 3; ++k){
//     for (int l = 0; l < 3; ++l){
//       DDPP[54*0 +9*0 +3*k +l] = FPPP[27*0 +9*0 +3*k +l] + ABx*DPPP[27*0 +9*0 +3*k +l]; // i=x, xx xx
//       DDPP[54*0 +9*1 +3*k +l] = FPPP[27*0 +9*1 +3*k +l] + ABy*DPPP[27*0 +9*0 +3*k +l]; // i=y, xx xy
//       DDPP[54*0 +9*2 +3*k +l] = FPPP[27*0 +9*2 +3*k +l] + ABz*DPPP[27*0 +9*0 +3*k +l]; // i=z, xx xz
//       DDPP[54*0 +9*3 +3*k +l] = FPPP[27*1 +9*1 +3*k +l] + ABy*DPPP[27*0 +9*1 +3*k +l]; // i=y, xx yy
//       DDPP[54*0 +9*4 +3*k +l] = FPPP[27*1 +9*2 +3*k +l] + ABz*DPPP[27*0 +9*1 +3*k +l]; // i=z, xx yz
//       DDPP[54*0 +9*5 +3*k +l] = FPPP[27*2 +9*2 +3*k +l] + ABz*DPPP[27*0 +9*2 +3*k +l]; // i=z, xx zz

//       DDPP[54*1 +9*0 +3*k +l] = FPPP[27*1 +9*0 +3*k +l] + ABx*DPPP[27*1 +9*0 +3*k +l]; // i=x, xy xx
//       DDPP[54*1 +9*1 +3*k +l] = FPPP[27*1 +9*1 +3*k +l] + ABy*DPPP[27*1 +9*0 +3*k +l]; // i=y, xy xy
//       DDPP[54*1 +9*2 +3*k +l] = FPPP[27*1 +9*2 +3*k +l] + ABz*DPPP[27*1 +9*0 +3*k +l]; // i=z, xy xz
//       DDPP[54*1 +9*3 +3*k +l] = FPPP[27*3 +9*1 +3*k +l] + ABy*DPPP[27*1 +9*1 +3*k +l]; // i=y, xy yy
//       DDPP[54*1 +9*4 +3*k +l] = FPPP[27*3 +9*2 +3*k +l] + ABz*DPPP[27*1 +9*1 +3*k +l]; // i=z, xy yz
//       DDPP[54*1 +9*5 +3*k +l] = FPPP[27*4 +9*2 +3*k +l] + ABz*DPPP[27*1 +9*2 +3*k +l]; // i=z, xy zz

//       DDPP[54*2 +9*0 +3*k +l] = FPPP[27*2 +9*0 +3*k +l] + ABx*DPPP[27*2 +9*0 +3*k +l]; // i=x, xz xx
//       DDPP[54*2 +9*1 +3*k +l] = FPPP[27*2 +9*1 +3*k +l] + ABy*DPPP[27*2 +9*0 +3*k +l]; // i=y, xz xy
//       DDPP[54*2 +9*2 +3*k +l] = FPPP[27*2 +9*2 +3*k +l] + ABz*DPPP[27*2 +9*0 +3*k +l]; // i=z, xz xz
//       DDPP[54*2 +9*3 +3*k +l] = FPPP[27*4 +9*1 +3*k +l] + ABy*DPPP[27*2 +9*1 +3*k +l]; // i=y, xz yy
//       DDPP[54*2 +9*4 +3*k +l] = FPPP[27*4 +9*2 +3*k +l] + ABz*DPPP[27*2 +9*1 +3*k +l]; // i=z, xz yz
//       DDPP[54*2 +9*5 +3*k +l] = FPPP[27*5 +9*2 +3*k +l] + ABz*DPPP[27*2 +9*2 +3*k +l]; // i=z, xz zz

//       DDPP[54*3 +9*0 +3*k +l] = FPPP[27*3 +9*0 +3*k +l] + ABx*DPPP[27*3 +9*0 +3*k +l]; // i=x, yy xx
//       DDPP[54*3 +9*1 +3*k +l] = FPPP[27*3 +9*1 +3*k +l] + ABy*DPPP[27*3 +9*0 +3*k +l]; // i=y, yy xy
//       DDPP[54*3 +9*2 +3*k +l] = FPPP[27*3 +9*2 +3*k +l] + ABz*DPPP[27*3 +9*0 +3*k +l]; // i=z, yy xz
//       DDPP[54*3 +9*3 +3*k +l] = FPPP[27*6 +9*1 +3*k +l] + ABy*DPPP[27*3 +9*1 +3*k +l]; // i=y, yy yy
//       DDPP[54*3 +9*4 +3*k +l] = FPPP[27*6 +9*2 +3*k +l] + ABz*DPPP[27*3 +9*1 +3*k +l]; // i=z, yy yz
//       DDPP[54*3 +9*5 +3*k +l] = FPPP[27*7 +9*2 +3*k +l] + ABz*DPPP[27*3 +9*2 +3*k +l]; // i=z, yy zz

//       DDPP[54*4 +9*0 +3*k +l] = FPPP[27*4 +9*0 +3*k +l] + ABx*DPPP[27*4 +9*0 +3*k +l]; // i=x, yz xx
//       DDPP[54*4 +9*1 +3*k +l] = FPPP[27*4 +9*1 +3*k +l] + ABy*DPPP[27*4 +9*0 +3*k +l]; // i=y, yz xy
//       DDPP[54*4 +9*2 +3*k +l] = FPPP[27*4 +9*2 +3*k +l] + ABz*DPPP[27*4 +9*0 +3*k +l]; // i=z, yz xz
//       DDPP[54*4 +9*3 +3*k +l] = FPPP[27*7 +9*1 +3*k +l] + ABy*DPPP[27*4 +9*1 +3*k +l]; // i=y, yz yy
//       DDPP[54*4 +9*4 +3*k +l] = FPPP[27*7 +9*2 +3*k +l] + ABz*DPPP[27*4 +9*1 +3*k +l]; // i=z, yz yz
//       DDPP[54*4 +9*5 +3*k +l] = FPPP[27*8 +9*2 +3*k +l] + ABz*DPPP[27*4 +9*2 +3*k +l]; // i=z, yz zz

//       DDPP[54*5 +9*0 +3*k +l] = FPPP[27*5 +9*0 +3*k +l] + ABx*DPPP[27*5 +9*0 +3*k +l]; // i=x, zz xx
//       DDPP[54*5 +9*1 +3*k +l] = FPPP[27*5 +9*1 +3*k +l] + ABy*DPPP[27*5 +9*0 +3*k +l]; // i=y, zz xy
//       DDPP[54*5 +9*2 +3*k +l] = FPPP[27*5 +9*2 +3*k +l] + ABz*DPPP[27*5 +9*0 +3*k +l]; // i=z, zz xz
//       DDPP[54*5 +9*3 +3*k +l] = FPPP[27*8 +9*1 +3*k +l] + ABy*DPPP[27*5 +9*1 +3*k +l]; // i=y, zz yy
//       DDPP[54*5 +9*4 +3*k +l] = FPPP[27*8 +9*2 +3*k +l] + ABz*DPPP[27*5 +9*1 +3*k +l]; // i=z, zz yz
//       DDPP[54*5 +9*5 +3*k +l] = FPPP[27*9 +9*2 +3*k +l] + ABz*DPPP[27*5 +9*2 +3*k +l]; // i=z, zz zz
//     }
//   }

    // convert 6D to 5D
    double Ddpp[270]; // 6*5*3*3
    this->convert6Dto5D_J(6, 3, 3, DDPP, Ddpp);
//   for (int i = 0; i < 6; ++i){
//     for (int k = 0; k < 3; ++k){
//       for (int l = 0; l < 3; ++l){
//  Ddpp[45*i +9*0 +3*k +l] = DDPP[54*i +9*1 +3*k +l]; // xy
//  Ddpp[45*i +9*1 +3*k +l] = DDPP[54*i +9*2 +3*k +l]; // xz
//  Ddpp[45*i +9*2 +3*k +l] = DDPP[54*i +9*4 +3*k +l]; // yz
//  Ddpp[45*i +9*3 +3*k +l] = 0.5*(DDPP[54*i +9*0 +3*k +l] -DDPP[54*i +9*3 +3*k +l]); // xx-yy
//  Ddpp[45*i +9*4 +3*k +l] = INV_SQRT3*(DDPP[54*i +9*5 +3*k +l]
//                       -0.5*(DDPP[54*i +9*0 +3*k +l] +DDPP[54*i +9*3 +3*k +l])); // 2zz-(xx+yy)
//       }
//     }
//   }

    this->convert6Dto5D_I(5, 3, 3, Ddpp, this->ERI);
//   for (int j = 0; j < 5; ++j){
//     for (int k = 0; k < 3; ++k){
//       for (int l = 0; l < 3; ++l){
//  this->ERI[45*0 +9*j +3*k +l] = Ddpp[45*1 +9*j +3*k +l]; // xy
//  this->ERI[45*1 +9*j +3*k +l] = Ddpp[45*2 +9*j +3*k +l]; // xz
//  this->ERI[45*2 +9*j +3*k +l] = Ddpp[45*4 +9*j +3*k +l]; // yz
//  this->ERI[45*3 +9*j +3*k +l] = 0.5*(Ddpp[45*0 +9*j +3*k +l] -Ddpp[45*3 +9*j +3*k +l]); // xx-yy
//  this->ERI[45*4 +9*j +3*k +l] = INV_SQRT3*(Ddpp[45*5 +9*j +3*k +l]
//                        -0.5*(Ddpp[45*0 +9*j +3*k +l] +Ddpp[45*3 +9*j +3*k +l])); // 2zz-(xx+yy)
//       }
//     }
//}
}


////////////////////////////////////////////////////////////////////////
// constract DDPD
void DfTEI::contractDDPD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    this->contractDDDP(IJ, this->swapShellPair(KL));
    for (int i = 0; i < 5*5; ++i) {
        this->transform(5, 3, &(this->ERI[15*i]));
    }
}


////////////////////////////////////////////////////////////////////////
// constract DDDS
void DfTEI::contractDDDS(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSDS[36]; // 6*6
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }

    double FSDS[60]; //10*6
    for (int i = 0; i < 60; ++i) {
        FSDS[i] = 0.0;
    }

    double GSDS[90]; //15*6
    for (int i = 0; i < 90; ++i) {
        GSDS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 7);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 6);
            this->primitiveDSSS(*ij, *kl, 5);
            this->primitiveFSSS(*ij, *kl, 4);
            this->primitiveGSSS(*ij, *kl, 3);

            this->primitivePSPS(*ij, *kl, 2);

            this->primitiveDSPS(*ij, *kl, 2);
            this->primitiveDSDS(*ij, *kl, 1); // needs DSPS, PSPS

            this->primitiveFSPS(*ij, *kl, 2);
            this->primitiveFSDS(*ij, *kl, 1); // needs FSPS, DSPS

            this->primitiveGSPS(*ij, *kl, 2);
            this->primitiveGSDS(*ij, *kl, 1); // needs GSPS, FSPS

            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
            for (int i = 0; i < 60; ++i) {
                FSDS[i] += this->pFSDS[0][i];
            }
            for (int i = 0; i < 90; ++i) {
                GSDS[i] += this->pGSDS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
//   const double CDx = KL.AB.x();
//   const double CDy = KL.AB.y();
//   const double CDz = KL.AB.z();

    // HRR (DPDS)
    double DPDS[108]; // 6*3*6
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSDS, DSDS, 6, 1, DPDS);

    // HRR (FPDS)
    double FPDS[180]; // 10*3*6
    this->HRR_FPXX(ABx, ABy, ABz,
                   GSDS, FSDS, 6, 1, FPDS);

    // HRR (DDDS)
    double DDDS[216]; // 6*6*6
    this->HRR_DDXX(ABx, ABy, ABz,
                   FPDS, DPDS, 6, 1, DDDS);

    // convert 6D to 5D
    double DDds[180]; // 6*6*5
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            DDds[30*i +5*j +0] = DDDS[36*i +6*j +1]; // xy
            DDds[30*i +5*j +1] = DDDS[36*i +6*j +2]; // yz
            DDds[30*i +5*j +2] = DDDS[36*i +6*j +4]; // xz
            DDds[30*i +5*j +3] = 0.5*(DDDS[36*i +6*j +0] -DDDS[36*i +6*j +3]); // xx-yy
            DDds[30*i +5*j +4] = INV_SQRT3*(DDDS[36*i +6*j +5]
                                            -0.5*(DDDS[36*i +6*j +0] +DDDS[36*i +6*j +3])); // 2zz-(xx+yy)
        }
    }

    double Ddds[150]; // 6*5*5
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 5; ++k) {
            Ddds[25*i +5*0 +k] = DDds[30*i +5*1 +k]; // xy
            Ddds[25*i +5*1 +k] = DDds[30*i +5*2 +k]; // yz
            Ddds[25*i +5*2 +k] = DDds[30*i +5*4 +k]; // xz
            Ddds[25*i +5*3 +k] = 0.5*(DDds[30*i +5*0 +k] -DDds[30*i +5*3 +k]); // xx-yy
            Ddds[25*i +5*4 +k] = INV_SQRT3*(DDds[30*i +5*5 +k]
                                            -0.5*(DDds[30*i +5*0 +k] +DDds[30*i +5*3 +k])); // 2zz-(xx+yy)
        }
    }

    for (int j = 0; j < 5; ++j) {
        for (int k = 0; k < 5; ++k) {
            this->ERI[25*0 +5*j +k] = Ddds[25*1 +5*j +k]; // xy
            this->ERI[25*1 +5*j +k] = Ddds[25*2 +5*j +k]; // yz
            this->ERI[25*2 +5*j +k] = Ddds[25*4 +5*j +k]; // xz
            this->ERI[25*3 +5*j +k] = 0.5*(Ddds[25*0 +5*j +k] -Ddds[25*3 +5*j +k]); // xx-yy
            this->ERI[25*4 +5*j +k] = INV_SQRT3*(Ddds[25*5 +5*j +k]
                                                 -0.5*(Ddds[25*0 +5*j +k] +Ddds[25*3 +5*j +k])); // 2zz-(xx+yy)
        }
    }
}


////////////////////////////////////////////////////////////////////////
// constract DDDP
void DfTEI::contractDDDP(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSDS[36]; // 6*6
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }
    double DSFS[60]; // 6*10
    for (int i = 0; i < 60; ++i) {
        DSFS[i] = 0.0;
    }
    double FSDS[60]; //10*6
    for (int i = 0; i < 60; ++i) {
        FSDS[i] = 0.0;
    }
    double FSFS[100]; //10*10
    for (int i = 0; i < 100; ++i) {
        FSFS[i] = 0.0;
    }
    double GSDS[90]; //15*6
    for (int i = 0; i < 90; ++i) {
        GSDS[i] = 0.0;
    }
    double GSFS[150]; //15*10
    for (int i = 0; i < 150; ++i) {
        GSFS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 8);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 7);
            this->primitiveDSSS(*ij, *kl, 6);
            this->primitiveFSSS(*ij, *kl, 5);
            this->primitiveGSSS(*ij, *kl, 4);

            this->primitiveSSPS(*ij, *kl, 3, 2);

            this->primitivePSPS(*ij, *kl, 3);
            this->primitivePSDS(*ij, *kl, 2); // needs ssps

            this->primitiveDSPS(*ij, *kl, 3);
            this->primitiveDSDS(*ij, *kl, 2);
            this->primitiveDSFS(*ij, *kl, 1); // needs dsds,psds

            this->primitiveFSPS(*ij, *kl, 3);
            this->primitiveFSDS(*ij, *kl, 2); // needs fsps,dsps
            this->primitiveFSFS(*ij, *kl, 1); // needs fsds,dsds

            this->primitiveGSPS(*ij, *kl, 3);
            this->primitiveGSDS(*ij, *kl, 2); // needs gsps,fsps
            this->primitiveGSFS(*ij, *kl, 1); // needs gsds,fsds

            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i]; //
            }
            for (int i = 0; i < 60; ++i) {
                DSFS[i] += this->pDSFS[0][i]; //
            }
            for (int i = 0; i < 60; ++i) {
                FSDS[i] += this->pFSDS[0][i]; //
            }
            for (int i = 0; i < 100; ++i) {
                FSFS[i] += this->pFSFS[0][i]; //
            }
            for (int i = 0; i < 90; ++i) {
                GSDS[i] += this->pGSDS[0][i]; //
            }
            for (int i = 0; i < 150; ++i) {
                GSFS[i] += this->pGSFS[0][i]; //
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
    const double CDx = KL.AB.x();
    const double CDy = KL.AB.y();
    const double CDz = KL.AB.z();

    // HRR (DSDP)
    double DSDP[108]; // 6*6*3
    this->HRR_XXDP(CDx, CDy, CDz,
                   DSFS, DSDS, 6, 1, DSDP);

    // HRR (FSDP)
    double FSDP[180]; // 10*6*3
    this->HRR_XXDP(CDx, CDy, CDz,
                   FSFS, FSDS, 10, 1, FSDP);

    // HRR (GSDP)
    double GSDP[270]; // 15*6*3
    this->HRR_XXDP(CDx, CDy, CDz,
                   GSFS, GSDS, 15, 1, GSDP);

    // HRR (DPDP)
    double DPDP[324]; // 6*3*6*3
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSDP, DSDP, 6, 3, DPDP);

    // HRR (FPDP)
    double FPDP[540]; // 10*3*6*3
    this->HRR_FPXX(ABx, ABy, ABz,
                   GSDP, FSDP, 6, 3, FPDP);

    // HRR (DDDP)
    double DDDP[648]; // 6*6*6*3
    this->HRR_DDXX(ABx, ABy, ABz,
                   FPDP, DPDP, 6, 3, DDDP);

    // convert 6D to 5D
    double DDdp[540]; // 6*6*5*3
    this->convert6Dto5D_K(6, 6, 3, DDDP, DDdp);

    double Dddp[450]; // 6*5*5*3
    this->convert6Dto5D_J(6, 5, 3, DDdp, Dddp);

    this->convert6Dto5D_I(5, 5, 3, Dddp, this->ERI);
}


////////////////////////////////////////////////////////////////////////
// constract DDDD
void DfTEI::contractDDDD(const ShellPair& IJ,
                         const ShellPair& KL)
{
    // initialize
    double DSDS[36]; // 6*6
    for (int i = 0; i < 36; ++i) {
        DSDS[i] = 0.0;
    }
    double DSFS[60]; // 6*10
    for (int i = 0; i < 60; ++i) {
        DSFS[i] = 0.0;
    }
    double DSGS[90]; // 6*15
    for (int i = 0; i < 90; ++i) {
        DSGS[i] = 0.0;
    }

    double FSDS[60]; //10*6
    for (int i = 0; i < 60; ++i) {
        FSDS[i] = 0.0;
    }
    double FSFS[100]; //10*10
    for (int i = 0; i < 100; ++i) {
        FSFS[i] = 0.0;
    }
    double FSGS[150]; //10*15
    for (int i = 0; i < 150; ++i) {
        FSGS[i] = 0.0;
    }

    double GSDS[90]; //15*6
    for (int i = 0; i < 90; ++i) {
        GSDS[i] = 0.0;
    }
    double GSFS[150]; //15*10
    for (int i = 0; i < 150; ++i) {
        GSFS[i] = 0.0;
    }
    double GSGS[225]; //15*15
    for (int i = 0; i < 225; ++i) {
        GSGS[i] = 0.0;
    }

    const TlPosition posA = IJ.A;
    const TlPosition posC = KL.A;

    // start ij-shell loop
    PrimitiveShellPairList::const_iterator ijEnd = IJ.pList.end();
    for (PrimitiveShellPairList::const_iterator ij = IJ.pList.begin(); ij != ijEnd; ++ij) {

        this->posP = ij->P;
        this->PA = this->posP - posA;

        // start kl-shell loop
        PrimitiveShellPairList::const_iterator klEnd = KL.pList.end();
        for (PrimitiveShellPairList::const_iterator kl = KL.pList.begin(); kl != klEnd; ++kl) {
            this->posQ = kl->P;
            this->primitiveSSSS(*ij, *kl, 9);

            if (this->isPrimitiveShellsCutoff() == true) {
                break;
            }

            this->posW = (ij->zeta*this->posP + kl->zeta*this->posQ) * this->m_divZetaEta;
            this->QC = this->posQ - posC;
            this->WP = this->posW - this->posP;
            this->WQ = this->posW - this->posQ;

            this->primitivePSSS(*ij, *kl, 8);
            this->primitiveDSSS(*ij, *kl, 7);
            this->primitiveFSSS(*ij, *kl, 6);
            this->primitiveGSSS(*ij, *kl, 5);

            this->primitiveSSPS(*ij, *kl, 4);
            this->primitiveSSDS(*ij, *kl, 3);

            this->primitivePSPS(*ij, *kl, 4);
            this->primitivePSDS(*ij, *kl, 3);
            this->primitivePSFS(*ij, *kl, 2); // needs SSDS

            this->primitiveDSPS(*ij, *kl, 4);
            this->primitiveDSDS(*ij, *kl, 3); // needs DSPS
            this->primitiveDSFS(*ij, *kl, 2); // needs DSDS, PSDS
            this->primitiveDSGS(*ij, *kl, 1); // needs DSFS, PSFS

            this->primitiveFSPS(*ij, *kl, 4);
            this->primitiveFSDS(*ij, *kl, 3); // needs FSPS, DSPS
            this->primitiveFSFS(*ij, *kl, 2); // needs FSDS, DSDS
            this->primitiveFSGS(*ij, *kl, 1); // needs FSFS, DSFS

            this->primitiveGSPS(*ij, *kl, 4);
            this->primitiveGSDS(*ij, *kl, 3); // needs GSPS, FSPS
            this->primitiveGSFS(*ij, *kl, 2); // needs GSDS, FSDS
            this->primitiveGSGS(*ij, *kl, 1); // needs GSFS, FSFS

            for (int i = 0; i < 36; ++i) {
                DSDS[i] += this->pDSDS[0][i];
            }
            for (int i = 0; i < 60; ++i) {
                DSFS[i] += this->pDSFS[0][i];
            }
            for (int i = 0; i < 90; ++i) {
                DSGS[i] += this->pDSGS[0][i];
            }

            for (int i = 0; i < 60; ++i) {
                FSDS[i] += this->pFSDS[0][i];
            }
            for (int i = 0; i < 100; ++i) {
                FSFS[i] += this->pFSFS[0][i];
            }
            for (int i = 0; i < 150; ++i) {
                FSGS[i] += this->pFSGS[0][i];
            }

            for (int i = 0; i < 90; ++i) {
                GSDS[i] += this->pGSDS[0][i];
            }
            for (int i = 0; i < 150; ++i) {
                GSFS[i] += this->pGSFS[0][i];
            }
            for (int i = 0; i < 225; ++i) {
                GSGS[i] += this->pGSGS[0][i];
            }
        }
    }

    const double ABx = IJ.AB.x();
    const double ABy = IJ.AB.y();
    const double ABz = IJ.AB.z();
    const double CDx = KL.AB.x();
    const double CDy = KL.AB.y();
    const double CDz = KL.AB.z();

    double DPDS[108];
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSDS, DSDS,  6, 1, DPDS);
    double DPFS[180];
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSFS, DSFS, 10, 1, DPFS);
    double DPGS[270];
    this->HRR_DPXX(ABx, ABy, ABz,
                   FSGS, DSGS, 15, 1, DPGS);

    double FPDS[180];
    this->HRR_FPXX(ABx, ABy, ABz,
                   GSDS, FSDS, 6, 1, FPDS);
    double FPFS[300];
    this->HRR_FPXX(ABx, ABy, ABz,
                   GSFS, FSFS, 10, 1, FPFS);
    double FPGS[450];
    this->HRR_FPXX(ABx, ABy, ABz,
                   GSGS, FSGS, 15, 1, FPGS);

    double DDDS[216];
    this->HRR_DDXX(ABx, ABy, ABz,
                   FPDS, DPDS, 6, 1, DDDS);
    double DDFS[360];
    this->HRR_DDXX(ABx, ABy, ABz,
                   FPFS, DPFS, 10, 1, DDFS);
    double DDGS[540];
    this->HRR_DDXX(ABx, ABy, ABz,
                   FPGS, DPGS, 15, 1, DDGS);

    double DDDP[648]; //
    this->HRR_XXDP(CDx, CDy, CDz,
                   DDFS, DDDS, 6, 6, DDDP);
    double DDFP[1080]; //
    this->HRR_XXFP(CDx, CDy, CDz,
                   DDGS, DDFS, 6, 6, DDFP);

    double DDDD[1296]; //
    this->HRR_XXDD(CDx, CDy, CDz,
                   DDFP, DDDP, 6, 6, DDDD);

//   {
//     std::ofstream ofs("DDDD_DDDP.txt", std::ios::out | std::ios::app);

//     int index = 0;
//     ofs << "6D" << std::endl;
//     for (int i = 0; i < 6; ++i){
//       for (int j = 0; j < 6; ++j){
//  for (int k = 0; k < 6; ++k){
//    for (int l = 0; l < 3; ++l){
//      ofs << TlUtils::format("(%2d %2d|%2d %2d) = % .8e", i, j, k, l, DDDP[index])
//      << std::endl;
//      ++index;
//    }
//  }
//       }
//     }

//     this->contractDDDP2(IJ, KL);
//   }

    // convert 6D to 5D
    double DDDd[1080]; // 6*6*6*5
    this->convert6Dto5D_L(6, 6, 6, DDDD, DDDd);

    double DDdd[900]; // 6*6*5*5
    this->convert6Dto5D_K(6, 6, 5, DDDd, DDdd);

    double Dddd[750]; // 6*5*5*5
    this->convert6Dto5D_J(6, 5, 5, DDdd, Dddd);

    this->convert6Dto5D_I(5, 5, 5, Dddd, this->ERI);
}

// =============================================================================
// primitive
// =============================================================================

void DfTEI::primitiveSSSS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double zeta = ij.zeta;
    const double eta = kl.zeta;

    const double dZetaEta = zeta + eta;
    this->m_divZetaEta = 1.0 / dZetaEta;
    this->m_dRho = zeta * eta * this->m_divZetaEta; // zeta*eta/(zeta+eta)
    const double t = this->m_dRho * (this->posP.squareDistanceFrom(this->posQ));

    //const std::vector<double> f = this->m_FmT.getFmT((nEndM -1), t);
    TlFmt& FmT = TlFmt::getInstance();
    FmT.getFmT((nEndM -1), t, this->pFmTBuf);

    const double dCoef = (1.0 / sqrt(dZetaEta)) * ij.dKab * kl.dKab;

    for (int m = nStartM; m < nEndM; ++m) {
        //this->pSSSS[m] = dCoef * f[m];
        this->pSSSS[m] = dCoef * this->pFmTBuf[m];
    }
}


void DfTEI::primitiveSSPS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;
        this->pSSPS[m][0] = QCx * this->pSSSS[m] + WQx * this->pSSSS[m1];
        this->pSSPS[m][1] = QCy * this->pSSSS[m] + WQy * this->pSSSS[m1];
        this->pSSPS[m][2] = QCz * this->pSSSS[m] + WQz * this->pSSSS[m1];
    }
}


void DfTEI::primitiveSSDS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;
        this->pSSDS[m][0] = QCx * this->pSSPS[m][0] + WQx * this->pSSPS[m1][0];
        this->pSSDS[m][1] = QCx * this->pSSPS[m][1] + WQx * this->pSSPS[m1][1];
        this->pSSDS[m][2] = QCx * this->pSSPS[m][2] + WQx * this->pSSPS[m1][2];
        this->pSSDS[m][3] = QCy * this->pSSPS[m][1] + WQy * this->pSSPS[m1][1];
        this->pSSDS[m][4] = QCy * this->pSSPS[m][2] + WQy * this->pSSPS[m1][2];
        this->pSSDS[m][5] = QCz * this->pSSPS[m][2] + WQz * this->pSSPS[m1][2];

        const double zssss = eta2 * (this->pSSSS[m] - re * this->pSSSS[m1]);
        this->pSSDS[m][0] += zssss;
        this->pSSDS[m][3] += zssss;
        this->pSSDS[m][5] += zssss;
    }
}


void DfTEI::primitivePSSS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double WPx = this->WP.x();
    const double WPy = this->WP.y();
    const double WPz = this->WP.z();
    const double PAx = this->PA.x();
    const double PAy = this->PA.y();
    const double PAz = this->PA.z();

    for (int m = nStartM; m < nEndM; ++m) {
        const double SSSS0 = this->pSSSS[m];
        const double SSSS1 = this->pSSSS[m +1];

        this->pPSSS[m][0] = PAx * SSSS0 + WPx * SSSS1;
        this->pPSSS[m][1] = PAy * SSSS0 + WPy * SSSS1;
        this->pPSSS[m][2] = PAz * SSSS0 + WPz * SSSS1;
    }
}


void DfTEI::primitiveDSSS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double WPx = this->WP.x();
    const double WPy = this->WP.y();
    const double WPz = this->WP.z();
    const double PAx = this->PA.x();
    const double PAy = this->PA.y();
    const double PAz = this->PA.z();
    const double zeta2 = 1.0 / (2.0 * ij.zeta);
    const double rz = this->m_dRho / ij.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;
        this->pDSSS[m][0] = PAx * this->pPSSS[m][0] + WPx * this->pPSSS[m1][0];
        this->pDSSS[m][1] = PAx * this->pPSSS[m][1] + WPx * this->pPSSS[m1][1];
        this->pDSSS[m][2] = PAx * this->pPSSS[m][2] + WPx * this->pPSSS[m1][2];
        this->pDSSS[m][3] = PAy * this->pPSSS[m][1] + WPy * this->pPSSS[m1][1];
        this->pDSSS[m][4] = PAy * this->pPSSS[m][2] + WPy * this->pPSSS[m1][2];
        this->pDSSS[m][5] = PAz * this->pPSSS[m][2] + WPz * this->pPSSS[m1][2];

        const double zssss = zeta2 * (this->pSSSS[m] - rz * this->pSSSS[m1]);
        this->pDSSS[m][0] += zssss;
        this->pDSSS[m][3] += zssss;
        this->pDSSS[m][5] += zssss;
    }
}


void DfTEI::primitiveFSSS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double WPx = this->WP.x();
    const double WPy = this->WP.y();
    const double WPz = this->WP.z();
    const double PAx = this->PA.x();
    const double PAy = this->PA.y();
    const double PAz = this->PA.z();
    const double zeta2 = 1.0 / (2.0 * ij.zeta);
    const double rz = this->m_dRho / ij.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        this->pFSSS[m][0] = PAx*this->pDSSS[m][0] +WPx*this->pDSSS[m1][0]; // i=x, a=xx
        this->pFSSS[m][1] = PAx*this->pDSSS[m][1] +WPx*this->pDSSS[m1][1]; // i=x, a=xy
        this->pFSSS[m][2] = PAx*this->pDSSS[m][2] +WPx*this->pDSSS[m1][2]; // i=x, a=xz
        this->pFSSS[m][3] = PAx*this->pDSSS[m][3] +WPx*this->pDSSS[m1][3]; // i=x, a=yy
        this->pFSSS[m][4] = PAx*this->pDSSS[m][4] +WPx*this->pDSSS[m1][4]; // xyz
        this->pFSSS[m][5] = PAx*this->pDSSS[m][5] +WPx*this->pDSSS[m1][5]; // xzz
        this->pFSSS[m][6] = PAy*this->pDSSS[m][3] +WPy*this->pDSSS[m1][3]; // yyy
        this->pFSSS[m][7] = PAy*this->pDSSS[m][4] +WPy*this->pDSSS[m1][4]; // yyz
        this->pFSSS[m][8] = PAy*this->pDSSS[m][5] +WPy*this->pDSSS[m1][5]; // yzz
        this->pFSSS[m][9] = PAz*this->pDSSS[m][5] +WPz*this->pDSSS[m1][5]; // zzz

        const double zpsss0 = zeta2*(this->pPSSS[m][0] - rz*this->pPSSS[m1][0]);
        const double zpsss1 = zeta2*(this->pPSSS[m][1] - rz*this->pPSSS[m1][1]);
        const double zpsss2 = zeta2*(this->pPSSS[m][2] - rz*this->pPSSS[m1][2]);
        this->pFSSS[m][0] += 2.0 * zpsss0; // i=x, a=xx
        this->pFSSS[m][1] +=       zpsss1; // i=x, a=xy
        this->pFSSS[m][2] +=       zpsss2; // i=x, a=xz
        //this->pFSSS[m][3] += 0.0;        // i=x, a=yy
        //this->pFSSS[m][4] += 0.0;        // i=x, a=yz
        //this->pFSSS[m][5] += 0.0;        // i=x, a=zz
        this->pFSSS[m][6] += 2.0 * zpsss1; // i=y, a=yy
        this->pFSSS[m][7] +=       zpsss2; // i=y, a=yz
        //this->pFSSS[m][8] += 0.0;        // i=y, a=zz
        this->pFSSS[m][9] += 2.0 * zpsss2; // i=z, a=zz
    }
}


void DfTEI::primitiveGSSS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double WPx = this->WP.x();
    const double WPy = this->WP.y();
    const double WPz = this->WP.z();
    const double PAx = this->PA.x();
    const double PAy = this->PA.y();
    const double PAz = this->PA.z();
    const double zeta2 = 1.0 / (2.0 * ij.zeta);
    const double rz = this->m_dRho / ij.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        this->pGSSS[m][ 0] = PAx* this->pFSSS[m][0] +WPx*this->pFSSS[m1][0]; // xxxx *
        this->pGSSS[m][ 1] = PAx* this->pFSSS[m][1] +WPx*this->pFSSS[m1][1]; // xxxy
        this->pGSSS[m][ 2] = PAx* this->pFSSS[m][2] +WPx*this->pFSSS[m1][2]; // xxxz
        this->pGSSS[m][ 3] = PAx* this->pFSSS[m][3] +WPx*this->pFSSS[m1][3]; // xxyy *
        this->pGSSS[m][ 4] = PAx* this->pFSSS[m][4] +WPx*this->pFSSS[m1][4]; // xxyz
        this->pGSSS[m][ 5] = PAx* this->pFSSS[m][5] +WPx*this->pFSSS[m1][5]; // xxzz *
        this->pGSSS[m][ 6] = PAx* this->pFSSS[m][6] +WPx*this->pFSSS[m1][6]; // xyyy
        this->pGSSS[m][ 7] = PAx* this->pFSSS[m][7] +WPx*this->pFSSS[m1][7]; // xyyz
        this->pGSSS[m][ 8] = PAx* this->pFSSS[m][8] +WPx*this->pFSSS[m1][8]; // xyzz
        this->pGSSS[m][ 9] = PAx* this->pFSSS[m][9] +WPx*this->pFSSS[m1][9]; // xzzz
        this->pGSSS[m][10] = PAy* this->pFSSS[m][6] +WPy*this->pFSSS[m1][6]; // yyyy *
        this->pGSSS[m][11] = PAy* this->pFSSS[m][7] +WPy*this->pFSSS[m1][7]; // yyyz
        this->pGSSS[m][12] = PAy* this->pFSSS[m][8] +WPy*this->pFSSS[m1][8]; // yyzz *
        this->pGSSS[m][13] = PAy* this->pFSSS[m][9] +WPy*this->pFSSS[m1][9]; // yzzz
        this->pGSSS[m][14] = PAz* this->pFSSS[m][9] +WPz*this->pFSSS[m1][9]; // zzzz *

        const double zdsss0 = zeta2 * (this->pDSSS[m][0] -rz*this->pDSSS[m1][0]); // xx
        const double zdsss1 = zeta2 * (this->pDSSS[m][1] -rz*this->pDSSS[m1][1]); // xy
        const double zdsss2 = zeta2 * (this->pDSSS[m][2] -rz*this->pDSSS[m1][2]); // xz
        const double zdsss3 = zeta2 * (this->pDSSS[m][3] -rz*this->pDSSS[m1][3]); // yy
        const double zdsss4 = zeta2 * (this->pDSSS[m][4] -rz*this->pDSSS[m1][4]); // yz
        const double zdsss5 = zeta2 * (this->pDSSS[m][5] -rz*this->pDSSS[m1][5]); // zz

        this->pGSSS[m][ 0] += 3.0 * zdsss0; // xxxx (i=x, a=xxx) *
        this->pGSSS[m][ 1] += 2.0 * zdsss1; // xxxy (i=x, a=xxy)
        this->pGSSS[m][ 2] += 2.0 * zdsss2; // xxxz (i=x, a=xxz)
        this->pGSSS[m][ 3] +=       zdsss3; // xxyy (i=x, a=xyy) *
        this->pGSSS[m][ 4] +=       zdsss4; // xxyz (i=x, a=xyz)
        this->pGSSS[m][ 5] +=       zdsss5; // xxzz (i=x, a=xzz) *
        //this->pGSSS[m][ 6] += 0.0;        // xyyy (i=x, a=yyy)
        //this->pGSSS[m][ 7] += 0.0;        // xyyz (i=x, a=yyz)
        //this->pGSSS[m][ 8] += 0.0;        // xyzz (i=x, a=yzz)
        //this->pGSSS[m][ 9] += 0.0;        // xzzz (i=x, a=zzz)
        this->pGSSS[m][10] += 3.0 * zdsss3; // yyyy (i=y, a=yyy) *
        this->pGSSS[m][11] += 2.0 * zdsss4; // yyyz (i=y, a=yyz)
        this->pGSSS[m][12] +=       zdsss5; // yyzz (i=y, a=yzz) *
        //this->pGSSS[m][13] += 0.0;        // yzzz (i=y, a=zzz)
        this->pGSSS[m][14] += 3.0 * zdsss5; // zzzz (i=z, a=zzz) *
    }
}


void DfTEI::primitivePSPS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 3; ++i) {
            this->pPSPS[m][3*i +0] = QCx*this->pPSSS[m][i] +WQx*this->pPSSS[m1][i];
            this->pPSPS[m][3*i +1] = QCy*this->pPSSS[m][i] +WQy*this->pPSSS[m1][i];
            this->pPSPS[m][3*i +2] = QCz*this->pPSSS[m][i] +WQz*this->pPSSS[m1][i];
        }

        const double zssss = ze2 * this->pSSSS[m1];
        this->pPSPS[m][3*0 +0] += zssss;
        this->pPSPS[m][3*1 +1] += zssss;
        this->pPSPS[m][3*2 +2] += zssss;
    }
}


void DfTEI::primitivePSDS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = QC.x();
    const double QCy = QC.y();
    const double QCz = QC.z();
    const double WQx = WQ.x();
    const double WQy = WQ.y();
    const double WQz = WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 3; ++i) {
            this->pPSDS[m][6*i +0] = QCx * this->pPSPS[m][3*i +0] + WQx * this->pPSPS[m1][3*i +0];
            this->pPSDS[m][6*i +1] = QCx * this->pPSPS[m][3*i +1] + WQx * this->pPSPS[m1][3*i +1];
            this->pPSDS[m][6*i +2] = QCx * this->pPSPS[m][3*i +2] + WQx * this->pPSPS[m1][3*i +2];
            this->pPSDS[m][6*i +3] = QCy * this->pPSPS[m][3*i +1] + WQy * this->pPSPS[m1][3*i +1];
            this->pPSDS[m][6*i +4] = QCy * this->pPSPS[m][3*i +2] + WQy * this->pPSPS[m1][3*i +2];
            this->pPSDS[m][6*i +5] = QCz * this->pPSPS[m][3*i +2] + WQz * this->pPSPS[m1][3*i +2];
        }

        for (int i = 0; i < 3; ++i) {
            const double zpsss = eta2 * (this->pPSSS[m][i] - re * this->pPSSS[m1][i]);
            this->pPSDS[m][6*i +0] += zpsss;
            this->pPSDS[m][6*i +3] += zpsss;
            this->pPSDS[m][6*i +5] += zpsss;
        }

        const double zssps_0 = ze2 * this->pSSPS[m1][0];
        const double zssps_1 = ze2 * this->pSSPS[m1][1];
        const double zssps_2 = ze2 * this->pSSPS[m1][2];
        this->pPSDS[m][6*0 +0] += zssps_0;
        this->pPSDS[m][6*0 +1] += zssps_1;
        this->pPSDS[m][6*0 +2] += zssps_2;
        this->pPSDS[m][6*1 +3] += zssps_1;
        this->pPSDS[m][6*1 +4] += zssps_2;
        this->pPSDS[m][6*2 +5] += zssps_2;
    }
}


void DfTEI::primitivePSFS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = QC.x();
    const double QCy = QC.y();
    const double QCz = QC.z();
    const double WQx = WQ.x();
    const double WQy = WQ.y();
    const double WQz = WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 3; ++i) {
            this->pPSFS[m][10*i +0] = QCx * this->pPSDS[m][6*i +0] + WQx * this->pPSDS[m1][6*i +0];
            this->pPSFS[m][10*i +1] = QCx * this->pPSDS[m][6*i +1] + WQx * this->pPSDS[m1][6*i +1];
            this->pPSFS[m][10*i +2] = QCx * this->pPSDS[m][6*i +2] + WQx * this->pPSDS[m1][6*i +2];
            this->pPSFS[m][10*i +3] = QCx * this->pPSDS[m][6*i +3] + WQx * this->pPSDS[m1][6*i +3];
            this->pPSFS[m][10*i +4] = QCx * this->pPSDS[m][6*i +4] + WQx * this->pPSDS[m1][6*i +4];
            this->pPSFS[m][10*i +5] = QCx * this->pPSDS[m][6*i +5] + WQx * this->pPSDS[m1][6*i +5];
            this->pPSFS[m][10*i +6] = QCy * this->pPSDS[m][6*i +3] + WQy * this->pPSDS[m1][6*i +3];
            this->pPSFS[m][10*i +7] = QCy * this->pPSDS[m][6*i +4] + WQy * this->pPSDS[m1][6*i +4];
            this->pPSFS[m][10*i +8] = QCy * this->pPSDS[m][6*i +5] + WQy * this->pPSDS[m1][6*i +5];
            this->pPSFS[m][10*i +9] = QCz * this->pPSDS[m][6*i +5] + WQz * this->pPSDS[m1][6*i +5];
        }

        for (int i = 0; i < 3; ++i) {
            const double zpsps_0 = eta2 * (this->pPSPS[m][3*i +0] - re * this->pPSPS[m1][3*i +0]);
            const double zpsps_1 = eta2 * (this->pPSPS[m][3*i +1] - re * this->pPSPS[m1][3*i +1]);
            const double zpsps_2 = eta2 * (this->pPSPS[m][3*i +2] - re * this->pPSPS[m1][3*i +2]);

            this->pPSFS[m][10*i +0] += 2.0 * zpsps_0; // c= xx
            this->pPSFS[m][10*i +1] +=       zpsps_1; // xy
            this->pPSFS[m][10*i +2] +=       zpsps_2; // xz
            //this->pPSFS[m][10*i +3] += 0.0;         // yy
            //this->pPSFS[m][10*i +4] += 0.0;         // yz
            //this->pPSFS[m][10*i +5] += 0.0;         // zz
            this->pPSFS[m][10*i +6] += 2.0 * zpsps_1; // yy
            this->pPSFS[m][10*i +7] +=       zpsps_2; // yz
            //this->pPSFS[m][10*i +8] += 0.0;         // zz
            this->pPSFS[m][10*i +9] += 2.0 * zpsps_2; // zz
        }

        const double zssds_0 = ze2 * this->pSSDS[m1][0];
        const double zssds_1 = ze2 * this->pSSDS[m1][1];
        const double zssds_2 = ze2 * this->pSSDS[m1][2];
        const double zssds_3 = ze2 * this->pSSDS[m1][3];
        const double zssds_4 = ze2 * this->pSSDS[m1][4];
        const double zssds_5 = ze2 * this->pSSDS[m1][5];

        this->pPSFS[m][10*0 +0] +=       zssds_0; // i=x a=x
        this->pPSFS[m][10*0 +1] +=       zssds_1;
        this->pPSFS[m][10*0 +2] +=       zssds_2;
        this->pPSFS[m][10*0 +3] +=       zssds_3;
        this->pPSFS[m][10*0 +4] +=       zssds_4;
        this->pPSFS[m][10*0 +5] +=       zssds_5;
        //this->pPSFS[m][10*0 +6] += 0.0;
        //this->pPSFS[m][10*0 +7] += 0.0;
        //this->pPSFS[m][10*0 +8] += 0.0;
        //this->pPSFS[m][10*0 +9] += 0.0;

        //this->pPSFS[m][10*1 +0] += 0.0;
        //this->pPSFS[m][10*1 +1] += 0.0;
        //this->pPSFS[m][10*1 +2] += 0.0;
        //this->pPSFS[m][10*1 +3] += 0.0;
        //this->pPSFS[m][10*1 +4] += 0.0;
        //this->pPSFS[m][10*1 +5] += 0.0;
        this->pPSFS[m][10*1 +6] +=       zssds_3; // i=y a=y
        this->pPSFS[m][10*1 +7] +=       zssds_4;
        this->pPSFS[m][10*1 +8] +=       zssds_5;
        //this->pPSFS[m][10*1 +9] += 0.0;

        //this->pPSFS[m][10*0 +0] += 0.0;
        //this->pPSFS[m][10*2 +1] += 0.0;
        //this->pPSFS[m][10*2 +2] += 0.0;
        //this->pPSFS[m][10*2 +3] += 0.0;
        //this->pPSFS[m][10*2 +4] += 0.0;
        //this->pPSFS[m][10*2 +5] += 0.0;
        //this->pPSFS[m][10*2 +6] += 0.0;
        //this->pPSFS[m][10*2 +7] += 0.0;
        //this->pPSFS[m][10*2 +8] += 0.0;
        this->pPSFS[m][10*2 +9] +=       zssds_5; // i=z a=z
    }
}


void DfTEI::primitiveDSPS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    //const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    //const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;
        for (int i = 0; i < 6; ++i) {
            this->pDSPS[m][3*i +0] = QCx*this->pDSSS[m][i] +WQx*this->pDSSS[m1][i];
            this->pDSPS[m][3*i +1] = QCy*this->pDSSS[m][i] +WQy*this->pDSSS[m1][i];
            this->pDSPS[m][3*i +2] = QCz*this->pDSSS[m][i] +WQz*this->pDSSS[m1][i];
        }

        const double zpsss_0 = ze2 * this->pPSSS[m1][0];
        const double zpsss_1 = ze2 * this->pPSSS[m1][1];
        const double zpsss_2 = ze2 * this->pPSSS[m1][2];
        this->pDSPS[m][3*0 +0] += 2.0 * zpsss_0; // i=x a=xx
        //this->pDSPS[m][3*0 +1] += 0.0;         // i=y a=xx
        //this->pDSPS[m][3*0 +2] += 0.0;         // i=z a=xx

        this->pDSPS[m][3*1 +0] +=       zpsss_1; // i=x a=xy
        this->pDSPS[m][3*1 +1] +=       zpsss_0; // i=y a=xy
        //this->pDSPS[m][3*1 +2] += 0.0;         // i=z a=xy

        this->pDSPS[m][3*2 +0] +=       zpsss_2; // i=x a=xz
        //this->pDSPS[m][3*2 +1] += 0.0;         // i=y a=xz
        this->pDSPS[m][3*2 +2] +=       zpsss_0; // i=z a=xz

        //this->pDSPS[m][3*3 +0] += 0.0;         // i=x a=yy
        this->pDSPS[m][3*3 +1] += 2.0 * zpsss_1; // i=y a=yy
        //this->pDSPS[m][3*3 +2] += 0.0;         // i=z a=yy

        //this->pDSPS[m][3*4 +0] += 0.0;         // i=x a=yz
        this->pDSPS[m][3*4 +1] +=       zpsss_2; // i=y a=yz
        this->pDSPS[m][3*4 +2] +=       zpsss_1; // i=z a=yz

        //this->pDSPS[m][3*5 +0] += 0.0;         // i=x a=zz
        //this->pDSPS[m][3*5 +1] += 0.0;         // i=y a=zz
        this->pDSPS[m][3*5 +2] += 2.0 * zpsss_2; // i=z a=zz
    }
}


void DfTEI::primitiveDSDS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;
        for (int i = 0; i < 6; ++i) {
            this->pDSDS[m][6*i +0] = QCx * this->pDSPS[m][3*i +0] + WQx * this->pDSPS[m1][3*i +0];
            this->pDSDS[m][6*i +1] = QCx * this->pDSPS[m][3*i +1] + WQx * this->pDSPS[m1][3*i +1];
            this->pDSDS[m][6*i +2] = QCx * this->pDSPS[m][3*i +2] + WQx * this->pDSPS[m1][3*i +2];
            this->pDSDS[m][6*i +3] = QCy * this->pDSPS[m][3*i +1] + WQy * this->pDSPS[m1][3*i +1];
            this->pDSDS[m][6*i +4] = QCy * this->pDSPS[m][3*i +2] + WQy * this->pDSPS[m1][3*i +2];
            this->pDSDS[m][6*i +5] = QCz * this->pDSPS[m][3*i +2] + WQz * this->pDSPS[m1][3*i +2];
        }

        for (int i = 0; i < 6; ++i) {
            const double zdsss = eta2 * (this->pDSSS[m][i] -re * this->pDSSS[m1][i]);
            this->pDSDS[m][6*i +0] += zdsss;
            this->pDSDS[m][6*i +3] += zdsss;
            this->pDSDS[m][6*i +5] += zdsss;
        }

        const double zpsps_00 = ze2 * this->pPSPS[m1][3*0 +0];
        const double zpsps_01 = ze2 * this->pPSPS[m1][3*0 +1];
        const double zpsps_02 = ze2 * this->pPSPS[m1][3*0 +2];
        const double zpsps_10 = ze2 * this->pPSPS[m1][3*1 +0];
        const double zpsps_11 = ze2 * this->pPSPS[m1][3*1 +1];
        const double zpsps_12 = ze2 * this->pPSPS[m1][3*1 +2];
        const double zpsps_20 = ze2 * this->pPSPS[m1][3*2 +0];
        const double zpsps_21 = ze2 * this->pPSPS[m1][3*2 +1];
        const double zpsps_22 = ze2 * this->pPSPS[m1][3*2 +2];

//     this->pDSDS[m][6*0 +0] += 2.0 * zpsps_00; // i=x a=xx
//     this->pDSDS[m][6*1 +0] +=       zpsps_10; // i=x a=xy
//     this->pDSDS[m][6*2 +0] +=       zpsps_20; // i=x a=xz
//     this->pDSDS[m][6*0 +1] += 2.0 * zpsps_01; // i=y a=xx
//     this->pDSDS[m][6*1 +1] +=       zpsps_11;
//     this->pDSDS[m][6*2 +1] +=       zpsps_21;
//     this->pDSDS[m][6*0 +2] += 2.0 * zpsps_02;
//     this->pDSDS[m][6*1 +2] +=       zpsps_12;
//     this->pDSDS[m][6*2 +2] +=       zpsps_22;
//     this->pDSDS[m][6*1 +3] +=       zpsps_01;
//     this->pDSDS[m][6*3 +3] += 2.0 * zpsps_11;
//     this->pDSDS[m][6*4 +3] +=       zpsps_21;
//     this->pDSDS[m][6*1 +4] +=       zpsps_02;
//     this->pDSDS[m][6*3 +4] += 2.0 * zpsps_12;
//     this->pDSDS[m][6*4 +4] +=       zpsps_22;
//     this->pDSDS[m][6*2 +5] +=       zpsps_02;
//     this->pDSDS[m][6*4 +5] +=       zpsps_12;
//     this->pDSDS[m][6*5 +5] += 2.0 * zpsps_22;

        this->pDSDS[m][6*0 +0] += 2.0 * zpsps_00; // i=x a=xx
        this->pDSDS[m][6*0 +1] += 2.0 * zpsps_01; // i=y a=xx
        this->pDSDS[m][6*0 +2] += 2.0 * zpsps_02;

        this->pDSDS[m][6*1 +0] +=       zpsps_10; // i=x a=xy
        this->pDSDS[m][6*1 +1] +=       zpsps_11;
        this->pDSDS[m][6*1 +2] +=       zpsps_12;
        this->pDSDS[m][6*1 +3] +=       zpsps_01; // i=y
        this->pDSDS[m][6*1 +4] +=       zpsps_02;

        this->pDSDS[m][6*2 +0] +=       zpsps_20; // i=x a=xz
        this->pDSDS[m][6*2 +1] +=       zpsps_21;
        this->pDSDS[m][6*2 +2] +=       zpsps_22;
        this->pDSDS[m][6*2 +5] +=       zpsps_02; // i=x

        this->pDSDS[m][6*3 +3] += 2.0 * zpsps_11; // i=y a=yy
        this->pDSDS[m][6*3 +4] += 2.0 * zpsps_12;

        this->pDSDS[m][6*4 +3] +=       zpsps_21; // i=y a=yz
        this->pDSDS[m][6*4 +4] +=       zpsps_22;
        this->pDSDS[m][6*4 +5] +=       zpsps_12; // i=z

        this->pDSDS[m][6*5 +5] += 2.0 * zpsps_22; // i=z a=zz
    }
}


void DfTEI::primitiveDSFS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;
        for (int i = 0; i < 6; ++i) {
            this->pDSFS[m][10*i +0] = QCx * this->pDSDS[m][6*i +0] +WQx*this->pDSDS[m1][6*i +0]; // xxx
            this->pDSFS[m][10*i +1] = QCx * this->pDSDS[m][6*i +1] +WQx*this->pDSDS[m1][6*i +1]; // xxy
            this->pDSFS[m][10*i +2] = QCx * this->pDSDS[m][6*i +2] +WQx*this->pDSDS[m1][6*i +2]; // xxz
            this->pDSFS[m][10*i +3] = QCx * this->pDSDS[m][6*i +3] +WQx*this->pDSDS[m1][6*i +3]; // xyy
            this->pDSFS[m][10*i +4] = QCx * this->pDSDS[m][6*i +4] +WQx*this->pDSDS[m1][6*i +4]; // xyz
            this->pDSFS[m][10*i +5] = QCx * this->pDSDS[m][6*i +5] +WQx*this->pDSDS[m1][6*i +5]; // xzz
            this->pDSFS[m][10*i +6] = QCy * this->pDSDS[m][6*i +3] +WQy*this->pDSDS[m1][6*i +3]; // yyy
            this->pDSFS[m][10*i +7] = QCy * this->pDSDS[m][6*i +4] +WQy*this->pDSDS[m1][6*i +4]; // yyz
            this->pDSFS[m][10*i +8] = QCy * this->pDSDS[m][6*i +5] +WQy*this->pDSDS[m1][6*i +5]; // yzz
            this->pDSFS[m][10*i +9] = QCz * this->pDSDS[m][6*i +5] +WQz*this->pDSDS[m1][6*i +5]; // zzz
        }

        for (int i = 0; i < 6; ++i) {
            const double zdsps_0 = eta2 * (this->pDSPS[m][3*i +0] -re * this->pDSPS[m1][3*i +0]);
            const double zdsps_1 = eta2 * (this->pDSPS[m][3*i +1] -re * this->pDSPS[m1][3*i +1]);
            const double zdsps_2 = eta2 * (this->pDSPS[m][3*i +2] -re * this->pDSPS[m1][3*i +2]);

            this->pDSFS[m][10*i +0] += 2.0 * zdsps_0; // i=x c=xx
            this->pDSFS[m][10*i +1] +=       zdsps_1; // i=x c=xy
            this->pDSFS[m][10*i +2] +=       zdsps_2; // i=x c=xz
            //this->pDSFS[m][10*i +3] += 0.0;         // i=x c=yy
            //this->pDSFS[m][10*i +4] += 0.0;         // i=x c=yz
            //this->pDSFS[m][10*i +5] += 0.0;         // i=x c=zz
            this->pDSFS[m][10*i +6] += 2.0 * zdsps_1; // i=y c=yy
            this->pDSFS[m][10*i +7] +=       zdsps_2; // i=y c=yz
            //this->pDSFS[m][10*i +8] += 0.0;         // i=y c=zz
            this->pDSFS[m][10*i +9] += 2.0 * zdsps_2; // i=z c=zz
        }

        const double zpsds_00 = ze2 * this->pPSDS[m1][6*0 +0];
        const double zpsds_01 = ze2 * this->pPSDS[m1][6*0 +1];
        const double zpsds_02 = ze2 * this->pPSDS[m1][6*0 +2];
        const double zpsds_03 = ze2 * this->pPSDS[m1][6*0 +3];
        const double zpsds_04 = ze2 * this->pPSDS[m1][6*0 +4];
        const double zpsds_05 = ze2 * this->pPSDS[m1][6*0 +5];
        const double zpsds_10 = ze2 * this->pPSDS[m1][6*1 +0];
        const double zpsds_11 = ze2 * this->pPSDS[m1][6*1 +1];
        const double zpsds_12 = ze2 * this->pPSDS[m1][6*1 +2];
        const double zpsds_13 = ze2 * this->pPSDS[m1][6*1 +3];
        const double zpsds_14 = ze2 * this->pPSDS[m1][6*1 +4];
        const double zpsds_15 = ze2 * this->pPSDS[m1][6*1 +5];
        const double zpsds_20 = ze2 * this->pPSDS[m1][6*2 +0];
        const double zpsds_21 = ze2 * this->pPSDS[m1][6*2 +1];
        const double zpsds_22 = ze2 * this->pPSDS[m1][6*2 +2];
        const double zpsds_23 = ze2 * this->pPSDS[m1][6*2 +3];
        const double zpsds_24 = ze2 * this->pPSDS[m1][6*2 +4];
        const double zpsds_25 = ze2 * this->pPSDS[m1][6*2 +5];

        this->pDSFS[m][10*0 +0] += 2.0 * zpsds_00; // i=x a=xx
        this->pDSFS[m][10*0 +1] += 2.0 * zpsds_01; // i=x
        this->pDSFS[m][10*0 +2] += 2.0 * zpsds_02; // i=x
        this->pDSFS[m][10*0 +3] += 2.0 * zpsds_03; // i=x
        this->pDSFS[m][10*0 +4] += 2.0 * zpsds_04; // i=x
        this->pDSFS[m][10*0 +5] += 2.0 * zpsds_05; // i=x
        //this->pDSFS[m][10*0 +6] += 0.0;          // i=y
        //this->pDSFS[m][10*0 +7] += 0.0;          // i=y
        //this->pDSFS[m][10*0 +8] += 0.0;          // i=y
        //this->pDSFS[m][10*0 +9] += 0.0;          // i=z

        this->pDSFS[m][10*1 +0] +=       zpsds_10; // i=x a=xy
        this->pDSFS[m][10*1 +1] +=       zpsds_11; // i=x
        this->pDSFS[m][10*1 +2] +=       zpsds_12; // i=x
        this->pDSFS[m][10*1 +3] +=       zpsds_13; // i=x
        this->pDSFS[m][10*1 +4] +=       zpsds_14; // i=x
        this->pDSFS[m][10*1 +5] +=       zpsds_15; // i=x
        this->pDSFS[m][10*1 +6] +=       zpsds_03; // i=y
        this->pDSFS[m][10*1 +7] +=       zpsds_04; // i=y
        this->pDSFS[m][10*1 +8] +=       zpsds_05; // i=y
        //this->pDSFS[m][10*1 +9] += 0.0;          // i=z

        this->pDSFS[m][10*2 +0] +=       zpsds_20; // i=x a=xz
        this->pDSFS[m][10*2 +1] +=       zpsds_21; // i=x
        this->pDSFS[m][10*2 +2] +=       zpsds_22; // i=x
        this->pDSFS[m][10*2 +3] +=       zpsds_23; // i=x
        this->pDSFS[m][10*2 +4] +=       zpsds_24; // i=x
        this->pDSFS[m][10*2 +5] +=       zpsds_25; // i=x
        //this->pDSFS[m][10*2 +6] += 0.0;          // i=y
        //this->pDSFS[m][10*2 +7] += 0.0;          // i=y
        //this->pDSFS[m][10*2 +8] += 0.0;          // i=y
        this->pDSFS[m][10*2 +9] +=       zpsds_05; // i=z // ??

        //this->pDSFS[m][10*3 +0] += 0.0;          // i=x a=yy
        //this->pDSFS[m][10*3 +1] += 0.0;          // i=x
        //this->pDSFS[m][10*3 +2] += 0.0;          // i=x
        //this->pDSFS[m][10*3 +3] += 0.0;          // i=x
        //this->pDSFS[m][10*3 +4] += 0.0;          // i=x
        //this->pDSFS[m][10*3 +5] += 0.0;          // i=x
        this->pDSFS[m][10*3 +6] += 2.0 * zpsds_13; // i=y
        this->pDSFS[m][10*3 +7] += 2.0 * zpsds_14; // i=y
        this->pDSFS[m][10*3 +8] += 2.0 * zpsds_15; // i=y
        //this->pDSFS[m][10*3 +9] += 0.0;          // i=z

        //this->pDSFS[m][10*4 +0] += 0.0;          // i=x a=yz
        //this->pDSFS[m][10*4 +1] += 0.0;          // i=x
        //this->pDSFS[m][10*4 +2] += 0.0;          // i=x
        //this->pDSFS[m][10*4 +3] += 0.0;          // i=x
        //this->pDSFS[m][10*4 +4] += 0.0;          // i=x
        //this->pDSFS[m][10*4 +5] += 0.0;          // i=x
        this->pDSFS[m][10*4 +6] +=       zpsds_23; // i=y
        this->pDSFS[m][10*4 +7] +=       zpsds_24; // i=y
        this->pDSFS[m][10*4 +8] +=       zpsds_25; // i=y
        this->pDSFS[m][10*4 +9] +=       zpsds_15; // i=z

        //this->pDSFS[m][10*5 +0] += 0.0;          // i=x a=zz
        //this->pDSFS[m][10*5 +1] += 0.0;          // i=x
        //this->pDSFS[m][10*5 +2] += 0.0;          // i=x
        //this->pDSFS[m][10*5 +3] += 0.0;          // i=x
        //this->pDSFS[m][10*5 +4] += 0.0;          // i=x
        //this->pDSFS[m][10*5 +5] += 0.0;          // i=x
        //this->pDSFS[m][10*5 +6] += 0.0;          // i=y
        //this->pDSFS[m][10*5 +7] += 0.0;          // i=y
        //this->pDSFS[m][10*5 +8] += 0.0;          // i=y
        this->pDSFS[m][10*5 +9] += 2.0 * zpsds_25; // i=z
    }
}


void DfTEI::primitiveDSGS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;
        for (int i = 0; i < 6; ++i) {
            this->pDSGS[m][15*i + 0] = QCx * this->pDSFS[m][10*i +0] + WQx * this->pDSFS[m1][10*i +0]; // xxxx
            this->pDSGS[m][15*i + 1] = QCx * this->pDSFS[m][10*i +1] + WQx * this->pDSFS[m1][10*i +1]; // xxxy
            this->pDSGS[m][15*i + 2] = QCx * this->pDSFS[m][10*i +2] + WQx * this->pDSFS[m1][10*i +2]; // xxxz
            this->pDSGS[m][15*i + 3] = QCx * this->pDSFS[m][10*i +3] + WQx * this->pDSFS[m1][10*i +3]; // xxyy
            this->pDSGS[m][15*i + 4] = QCx * this->pDSFS[m][10*i +4] + WQx * this->pDSFS[m1][10*i +4]; // xxyz
            this->pDSGS[m][15*i + 5] = QCx * this->pDSFS[m][10*i +5] + WQx * this->pDSFS[m1][10*i +5]; // xxzz
            this->pDSGS[m][15*i + 6] = QCx * this->pDSFS[m][10*i +6] + WQx * this->pDSFS[m1][10*i +6]; // xyyy
            this->pDSGS[m][15*i + 7] = QCx * this->pDSFS[m][10*i +7] + WQx * this->pDSFS[m1][10*i +7]; // xyyz
            this->pDSGS[m][15*i + 8] = QCx * this->pDSFS[m][10*i +8] + WQx * this->pDSFS[m1][10*i +8]; // xyzz
            this->pDSGS[m][15*i + 9] = QCx * this->pDSFS[m][10*i +9] + WQx * this->pDSFS[m1][10*i +9]; // xzzz
            this->pDSGS[m][15*i +10] = QCy * this->pDSFS[m][10*i +6] + WQy * this->pDSFS[m1][10*i +6]; // yyyy
            this->pDSGS[m][15*i +11] = QCy * this->pDSFS[m][10*i +7] + WQy * this->pDSFS[m1][10*i +7]; // yyyz
            this->pDSGS[m][15*i +12] = QCy * this->pDSFS[m][10*i +8] + WQy * this->pDSFS[m1][10*i +8]; // yyzz
            this->pDSGS[m][15*i +13] = QCy * this->pDSFS[m][10*i +9] + WQy * this->pDSFS[m1][10*i +9]; // yzzz
            this->pDSGS[m][15*i +14] = QCz * this->pDSFS[m][10*i +9] + WQz * this->pDSFS[m1][10*i +9]; // zzzz
        }

        for (int i = 0; i < 6; ++i) {
            const double zdsds_0 = eta2 * (this->pDSDS[m][6*i +0] -re * this->pDSDS[m1][6*i +0]);
            const double zdsds_1 = eta2 * (this->pDSDS[m][6*i +1] -re * this->pDSDS[m1][6*i +1]);
            const double zdsds_2 = eta2 * (this->pDSDS[m][6*i +2] -re * this->pDSDS[m1][6*i +2]);
            const double zdsds_3 = eta2 * (this->pDSDS[m][6*i +3] -re * this->pDSDS[m1][6*i +3]);
            const double zdsds_4 = eta2 * (this->pDSDS[m][6*i +4] -re * this->pDSDS[m1][6*i +4]);
            const double zdsds_5 = eta2 * (this->pDSDS[m][6*i +5] -re * this->pDSDS[m1][6*i +5]);

            this->pDSGS[m][15*i + 0] += 3.0 * zdsds_0; // i=x c=xxx
            this->pDSGS[m][15*i + 1] += 2.0 * zdsds_1; // i=x xxy
            this->pDSGS[m][15*i + 2] += 2.0 * zdsds_2; // i=x xxz
            this->pDSGS[m][15*i + 3] +=       zdsds_3; // i=x xyy
            this->pDSGS[m][15*i + 4] +=       zdsds_4; // i=x xyz
            this->pDSGS[m][15*i + 5] +=       zdsds_5; // i=x xzz
            //this->pDSGS[m][15*i + 6] += 0.0;         // i=x yyy
            //this->pDSGS[m][15*i + 7] += 0.0;         // i=x yyz
            //this->pDSGS[m][15*i + 8] += 0.0;         // i=x yzz
            //this->pDSGS[m][15*i + 9] += 0.0;         // i=x zzz
            this->pDSGS[m][15*i +10] += 3.0 * zdsds_3; // i=y yyy
            this->pDSGS[m][15*i +11] += 2.0 * zdsds_4; // i=y yyz
            this->pDSGS[m][15*i +12] +=       zdsds_5; // i=y yzz
            //this->pDSGS[m][15*i +13] += 0.0;         // i=y zzz
            this->pDSGS[m][15*i +14] += 3.0 * zdsds_5; // i=z zzz
        }

        const double zpsfs_00 = ze2 * this->pPSFS[m1][10*0 +0];
        const double zpsfs_01 = ze2 * this->pPSFS[m1][10*0 +1];
        const double zpsfs_02 = ze2 * this->pPSFS[m1][10*0 +2];
        const double zpsfs_03 = ze2 * this->pPSFS[m1][10*0 +3];
        const double zpsfs_04 = ze2 * this->pPSFS[m1][10*0 +4];
        const double zpsfs_05 = ze2 * this->pPSFS[m1][10*0 +5];
        const double zpsfs_06 = ze2 * this->pPSFS[m1][10*0 +6];
        const double zpsfs_07 = ze2 * this->pPSFS[m1][10*0 +7];
        const double zpsfs_08 = ze2 * this->pPSFS[m1][10*0 +8];
        const double zpsfs_09 = ze2 * this->pPSFS[m1][10*0 +9];

        const double zpsfs_10 = ze2 * this->pPSFS[m1][10*1 +0];
        const double zpsfs_11 = ze2 * this->pPSFS[m1][10*1 +1];
        const double zpsfs_12 = ze2 * this->pPSFS[m1][10*1 +2];
        const double zpsfs_13 = ze2 * this->pPSFS[m1][10*1 +3];
        const double zpsfs_14 = ze2 * this->pPSFS[m1][10*1 +4];
        const double zpsfs_15 = ze2 * this->pPSFS[m1][10*1 +5];
        const double zpsfs_16 = ze2 * this->pPSFS[m1][10*1 +6];
        const double zpsfs_17 = ze2 * this->pPSFS[m1][10*1 +7];
        const double zpsfs_18 = ze2 * this->pPSFS[m1][10*1 +8];
        const double zpsfs_19 = ze2 * this->pPSFS[m1][10*1 +9];

        const double zpsfs_20 = ze2 * this->pPSFS[m1][10*2 +0];
        const double zpsfs_21 = ze2 * this->pPSFS[m1][10*2 +1];
        const double zpsfs_22 = ze2 * this->pPSFS[m1][10*2 +2];
        const double zpsfs_23 = ze2 * this->pPSFS[m1][10*2 +3];
        const double zpsfs_24 = ze2 * this->pPSFS[m1][10*2 +4];
        const double zpsfs_25 = ze2 * this->pPSFS[m1][10*2 +5];
        const double zpsfs_26 = ze2 * this->pPSFS[m1][10*2 +6];
        const double zpsfs_27 = ze2 * this->pPSFS[m1][10*2 +7];
        const double zpsfs_28 = ze2 * this->pPSFS[m1][10*2 +8];
        const double zpsfs_29 = ze2 * this->pPSFS[m1][10*2 +9];

        this->pDSGS[m][15*0 + 0] += 2.0 * zpsfs_00; // i=x a=xx
        this->pDSGS[m][15*0 + 1] += 2.0 * zpsfs_01; // i=x
        this->pDSGS[m][15*0 + 2] += 2.0 * zpsfs_02; // i=x
        this->pDSGS[m][15*0 + 3] += 2.0 * zpsfs_03; // i=x
        this->pDSGS[m][15*0 + 4] += 2.0 * zpsfs_04; // i=x
        this->pDSGS[m][15*0 + 5] += 2.0 * zpsfs_05; // i=x
        this->pDSGS[m][15*0 + 6] += 2.0 * zpsfs_06; // i=x
        this->pDSGS[m][15*0 + 7] += 2.0 * zpsfs_07; // i=x
        this->pDSGS[m][15*0 + 8] += 2.0 * zpsfs_08; // i=x
        this->pDSGS[m][15*0 + 9] += 2.0 * zpsfs_09; // i=x
        //this->pDSGS[m][15*0 +10] += 0.0;          // i=y
        //this->pDSGS[m][15*0 +11] += 0.0;          // i=y
        //this->pDSGS[m][15*0 +12] += 0.0;          // i=y
        //this->pDSGS[m][15*0 +13] += 0.0;          // i=y
        //this->pDSGS[m][15*0 +14] += 0.0;          // i=z

        this->pDSGS[m][15*1 + 0] +=       zpsfs_10; // i=x a=xy
        this->pDSGS[m][15*1 + 1] +=       zpsfs_11; // i=x
        this->pDSGS[m][15*1 + 2] +=       zpsfs_12; // i=x
        this->pDSGS[m][15*1 + 3] +=       zpsfs_13; // i=x
        this->pDSGS[m][15*1 + 4] +=       zpsfs_14; // i=x
        this->pDSGS[m][15*1 + 5] +=       zpsfs_15; // i=x
        this->pDSGS[m][15*1 + 6] +=       zpsfs_16; // i=x
        this->pDSGS[m][15*1 + 7] +=       zpsfs_17; // i=x
        this->pDSGS[m][15*1 + 8] +=       zpsfs_18; // i=x
        this->pDSGS[m][15*1 + 9] +=       zpsfs_19; // i=x
        this->pDSGS[m][15*1 +10] +=       zpsfs_06; // i=y
        this->pDSGS[m][15*1 +11] +=       zpsfs_07; // i=y
        this->pDSGS[m][15*1 +12] +=       zpsfs_08; // i=y
        this->pDSGS[m][15*1 +13] +=       zpsfs_09; // i=y
        //this->pDSGS[m][15*1 +14] += 0.0;          // i=z

        this->pDSGS[m][15*2 + 0] +=       zpsfs_20; // i=x a=xz
        this->pDSGS[m][15*2 + 1] +=       zpsfs_21; // i=x
        this->pDSGS[m][15*2 + 2] +=       zpsfs_22; // i=x
        this->pDSGS[m][15*2 + 3] +=       zpsfs_23; // i=x
        this->pDSGS[m][15*2 + 4] +=       zpsfs_24; // i=x
        this->pDSGS[m][15*2 + 5] +=       zpsfs_25; // i=x
        this->pDSGS[m][15*2 + 6] +=       zpsfs_26; // i=x
        this->pDSGS[m][15*2 + 7] +=       zpsfs_27; // i=x
        this->pDSGS[m][15*2 + 8] +=       zpsfs_28; // i=x
        this->pDSGS[m][15*2 + 9] +=       zpsfs_29; // i=x
        //this->pDSGS[m][15*2 +10] += 0.0;          // i=y
        //this->pDSGS[m][15*2 +11] += 0.0;          // i=y
        //this->pDSGS[m][15*2 +12] += 0.0;          // i=y
        //this->pDSGS[m][15*2 +13] += 0.0;          // i=y
        this->pDSGS[m][15*2 +14] +=       zpsfs_09; // i=z

        //this->pDSGS[m][15*3 + 0] += 0.0;          // i=x a=yy
        //this->pDSGS[m][15*3 + 1] += 0.0;          // i=x
        //this->pDSGS[m][15*3 + 2] += 0.0;          // i=x
        //this->pDSGS[m][15*3 + 3] += 0.0;          // i=x
        //this->pDSGS[m][15*3 + 4] += 0.0;          // i=x
        //this->pDSGS[m][15*3 + 5] += 0.0;          // i=x
        //this->pDSGS[m][15*3 + 6] += 0.0;          // i=x
        //this->pDSGS[m][15*3 + 7] += 0.0;          // i=x
        //this->pDSGS[m][15*3 + 8] += 0.0;          // i=x
        //this->pDSGS[m][15*3 + 9] += 0.0;          // i=x
        this->pDSGS[m][15*3 +10] += 2.0 * zpsfs_16; // i=y
        this->pDSGS[m][15*3 +11] += 2.0 * zpsfs_17; // i=y
        this->pDSGS[m][15*3 +12] += 2.0 * zpsfs_18; // i=y
        this->pDSGS[m][15*3 +13] += 2.0 * zpsfs_19; // i=y
        //this->pDSGS[m][15*3 +14] += 0.0;          // i=z

        //this->pDSGS[m][15*4 + 0] += 0.0;          // i=x a=yz
        //this->pDSGS[m][15*4 + 1] += 0.0;          // i=x
        //this->pDSGS[m][15*4 + 2] += 0.0;          // i=x
        //this->pDSGS[m][15*4 + 3] += 0.0;          // i=x
        //this->pDSGS[m][15*4 + 4] += 0.0;          // i=x
        //this->pDSGS[m][15*4 + 5] += 0.0;          // i=x
        //this->pDSGS[m][15*4 + 6] += 0.0;          // i=x
        //this->pDSGS[m][15*4 + 7] += 0.0;          // i=x
        //this->pDSGS[m][15*4 + 8] += 0.0;          // i=x
        //this->pDSGS[m][15*4 + 9] += 0.0;          // i=x
        this->pDSGS[m][15*4 +10] +=       zpsfs_26; // i=y
        this->pDSGS[m][15*4 +11] +=       zpsfs_27; // i=y
        this->pDSGS[m][15*4 +12] +=       zpsfs_28; // i=y
        this->pDSGS[m][15*4 +13] +=       zpsfs_29; // i=y
        this->pDSGS[m][15*4 +14] +=       zpsfs_19; // i=z

        //this->pDSGS[m][15*5 + 0] += 0.0;          // i=x a=zz
        //this->pDSGS[m][15*5 + 1] += 0.0;          // i=x
        //this->pDSGS[m][15*5 + 2] += 0.0;          // i=x
        //this->pDSGS[m][15*5 + 3] += 0.0;          // i=x
        //this->pDSGS[m][15*5 + 4] += 0.0;          // i=x
        //this->pDSGS[m][15*5 + 5] += 0.0;          // i=x
        //this->pDSGS[m][15*5 + 6] += 0.0;          // i=x
        //this->pDSGS[m][15*5 + 7] += 0.0;          // i=x
        //this->pDSGS[m][15*5 + 8] += 0.0;          // i=x
        //this->pDSGS[m][15*5 + 9] += 0.0;          // i=x
        //this->pDSGS[m][15*5 +10] += 0.0;          // i=y
        //this->pDSGS[m][15*5 +11] += 0.0;          // i=y
        //this->pDSGS[m][15*5 +12] += 0.0;          // i=y
        //this->pDSGS[m][15*5 +13] += 0.0;          // i=y
        this->pDSGS[m][15*5 +14] += 2.0 * zpsfs_29; // i=z
    }
}


void DfTEI::primitiveFSPS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    //const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    //const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 10; ++i) {
            this->pFSPS[m][3*i +0] = QCx * this->pFSSS[m][i] +WQx * this->pFSSS[m1][i]; // --- x
            this->pFSPS[m][3*i +1] = QCy * this->pFSSS[m][i] +WQy * this->pFSSS[m1][i]; // --- y
            this->pFSPS[m][3*i +2] = QCz * this->pFSSS[m][i] +WQz * this->pFSSS[m1][i]; // --- x
        }

        // [ab|(c+1i)d](m) = QC[ab|cd](m) +WQ[ab|cd](m+1)
        //                  +1/2(zeta+ita){ai[(a-1i)b|cd](m+1)}
        const double zdsss_0 = ze2 * this->pDSSS[m1][0]; // xx
        const double zdsss_1 = ze2 * this->pDSSS[m1][1]; // xy
        const double zdsss_2 = ze2 * this->pDSSS[m1][2]; // xz
        const double zdsss_3 = ze2 * this->pDSSS[m1][3]; // yy
        const double zdsss_4 = ze2 * this->pDSSS[m1][4]; // yz
        const double zdsss_5 = ze2 * this->pDSSS[m1][5]; // zz

        this->pFSPS[m][3*0 +0] += 3.0 * zdsss_0; // i=x, a=xxx
        //this->pFSPS[m][3*0 +1] += 0.0;         // i=y, a=xxx
        //this->pFSPS[m][3*0 +2] += 0.0;         // i=z, a=xxx

        this->pFSPS[m][3*1 +0] += 2.0 * zdsss_1; // i=x, a=xxy
        this->pFSPS[m][3*1 +1] +=       zdsss_0; // i=y, a=xxy
        //this->pFSPS[m][3*1 +2] += 0.0;         // i=z, a=xxy

        this->pFSPS[m][3*2 +0] += 2.0 * zdsss_2; // i=x, a=xxz
        //this->pFSPS[m][3*2 +1] += 0.0;         // i=y, a=xxz
        this->pFSPS[m][3*2 +2] +=       zdsss_0; // i=z, a=xxz

        this->pFSPS[m][3*3 +0] +=       zdsss_3; // i=x, a=xyy
        this->pFSPS[m][3*3 +1] += 2.0 * zdsss_1; // i=y, a=xyy
        //this->pFSPS[m][3*3 +2] += 0.0;         // i=z, a=xyy

        this->pFSPS[m][3*4 +0] +=       zdsss_4; // xyz x (i=x, a=xyz)
        this->pFSPS[m][3*4 +1] +=       zdsss_2; // xyz y (i=y, a=xyz)
        this->pFSPS[m][3*4 +2] +=       zdsss_1; // xyz z (i=z, a=xyz)

        this->pFSPS[m][3*5 +0] +=       zdsss_5; // xzz x (i=x, a=xzz)
        //this->pFSPS[m][3*5 +1] += 0.0;         // xzz y (i=y, a=xzz)
        this->pFSPS[m][3*5 +2] += 2.0 * zdsss_2; // xzz z (i=z, a=xzz)

        //this->pFSPS[m][3*6 +0] += 0.0;         // yyy x (i=x, a=yyy)
        this->pFSPS[m][3*6 +1] += 3.0 * zdsss_3; // yyy y (i=y, a=yyy)
        //this->pFSPS[m][3*6 +2] += 0.0;         // yyy z (i=z, a=yyy)

        //this->pFSPS[m][3*7 +0] += 0.0;         // yyz x (i=x, a=yyz)
        this->pFSPS[m][3*7 +1] += 2.0 * zdsss_4; // yyz y (i=y, a=yyz)
        this->pFSPS[m][3*7 +2] +=       zdsss_3; // yyz z (i=z, a=yyz)

        //this->pFSPS[m][3*8 +0] += 0.0;         // yzz x (i=x, a=yzz)
        this->pFSPS[m][3*8 +1] +=       zdsss_5; // yzz y (i=y, a=yzz)
        this->pFSPS[m][3*8 +2] += 2.0 * zdsss_4; // yzz z (i=z, a=yzz)

        //this->pFSPS[m][3*9 +0] += 0.0;         // zzz x (i=x, a=zzz)
        //this->pFSPS[m][3*9 +1] += 0.0;         // zzz y (i=y, a=zzz)
        this->pFSPS[m][3*9 +2] += 3.0 * zdsss_5; // zzz z (i=z, a=zzz)
    }
}


void DfTEI::primitiveFSDS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 10; ++i) {
            this->pFSDS[m][6*i +0] = QCx*this->pFSPS[m][3*i +0] +WQx*this->pFSPS[m1][3*i +0];
            this->pFSDS[m][6*i +1] = QCx*this->pFSPS[m][3*i +1] +WQx*this->pFSPS[m1][3*i +1];
            this->pFSDS[m][6*i +2] = QCx*this->pFSPS[m][3*i +2] +WQx*this->pFSPS[m1][3*i +2];
            this->pFSDS[m][6*i +3] = QCy*this->pFSPS[m][3*i +1] +WQy*this->pFSPS[m1][3*i +1];
            this->pFSDS[m][6*i +4] = QCy*this->pFSPS[m][3*i +2] +WQy*this->pFSPS[m1][3*i +2];
            this->pFSDS[m][6*i +5] = QCz*this->pFSPS[m][3*i +2] +WQz*this->pFSPS[m1][3*i +2];
        }

        for (int i = 0; i < 10; ++i) {
            const double zfsss = eta2 * (this->pFSSS[m][i] -re * this->pFSSS[m1][i]);
            this->pFSDS[m][6*i +0] += zfsss; //
            this->pFSDS[m][6*i +3] += zfsss; //
            this->pFSDS[m][6*i +5] += zfsss; //
        }

        // ai*(a-1i)
        const double zdsps_00 = ze2 * this->pDSPS[m1][3*0 +0]; // xx x
        const double zdsps_01 = ze2 * this->pDSPS[m1][3*0 +1]; // xx y
        const double zdsps_02 = ze2 * this->pDSPS[m1][3*0 +2]; // xx z
        const double zdsps_10 = ze2 * this->pDSPS[m1][3*1 +0]; // xy x
        const double zdsps_11 = ze2 * this->pDSPS[m1][3*1 +1]; // xy y
        const double zdsps_12 = ze2 * this->pDSPS[m1][3*1 +2]; // xy z
        const double zdsps_20 = ze2 * this->pDSPS[m1][3*2 +0]; // xz x
        const double zdsps_21 = ze2 * this->pDSPS[m1][3*2 +1]; // xz y
        const double zdsps_22 = ze2 * this->pDSPS[m1][3*2 +2]; // xz z
        const double zdsps_30 = ze2 * this->pDSPS[m1][3*3 +0]; //
        const double zdsps_31 = ze2 * this->pDSPS[m1][3*3 +1];
        const double zdsps_32 = ze2 * this->pDSPS[m1][3*3 +2];
        const double zdsps_40 = ze2 * this->pDSPS[m1][3*4 +0];
        const double zdsps_41 = ze2 * this->pDSPS[m1][3*4 +1];
        const double zdsps_42 = ze2 * this->pDSPS[m1][3*4 +2];
        const double zdsps_50 = ze2 * this->pDSPS[m1][3*5 +0];
        const double zdsps_51 = ze2 * this->pDSPS[m1][3*5 +1];
        const double zdsps_52 = ze2 * this->pDSPS[m1][3*5 +2];

        this->pFSDS[m][6*0 +0] += 3.0 * zdsps_00; // xxx xx (i=x, a=xxx)
        this->pFSDS[m][6*0 +1] += 3.0 * zdsps_01; // xxx xy (i=x, a=xxx)
        this->pFSDS[m][6*0 +2] += 3.0 * zdsps_02; // xxx xz (i=x, a=xxx)
        //this->pFSDS[m][6*0 +3] += 0.0; // xxx yy (i=y, a=xxx)
        //this->pFSDS[m][6*0 +4] += 0.0; // xxx yz (i=y, a=xxx)
        //this->pFSDS[m][6*0 +5] += 0.0; // xxx zz (i=z, a=xxx)

        this->pFSDS[m][6*1 +0] += 2.0 * zdsps_10; // xxy xx (i=x, a=xxy)
        this->pFSDS[m][6*1 +1] += 2.0 * zdsps_11; // xxy xy (i=x, a=xxy)
        this->pFSDS[m][6*1 +2] += 2.0 * zdsps_12; // xxy xz (i=x, a=xxy)
        this->pFSDS[m][6*1 +3] +=       zdsps_01; // xxy yy (i=y, a=xxy)
        this->pFSDS[m][6*1 +4] +=       zdsps_02; // xxy yz (i=y, a=xxy)
        //this->pFSDS[m][6*1 +5] += 0.0; // xxy zz (i=z, a=xxy)

        this->pFSDS[m][6*2 +0] += 2.0 * zdsps_20; // xxz xx (i=x, a=xxz)
        this->pFSDS[m][6*2 +1] += 2.0 * zdsps_21; // xxz xy (i=x, a=xxz)
        this->pFSDS[m][6*2 +2] += 2.0 * zdsps_22; // xxz xz (i=x, a=xxz)
        //this->pFSDS[m][6*2 +3] += 0.0;          // xxz yy (i=y, a=xxz)
        //this->pFSDS[m][6*2 +4] += 0.0;          // xxz yz (i=y, a=xxz)
        this->pFSDS[m][6*2 +5] +=       zdsps_02; // xxz zz (i=z, a=xxz)

        this->pFSDS[m][6*3 +0] +=       zdsps_30; // xyy xx (i=x, a=xyy)
        this->pFSDS[m][6*3 +1] +=       zdsps_31; // xyy xy (i=x, a=xyy)
        this->pFSDS[m][6*3 +2] +=       zdsps_32; // xyy xz (i=x, a=xyy)
        this->pFSDS[m][6*3 +3] += 2.0 * zdsps_11; // xyy yy (i=y, a=xyy)
        this->pFSDS[m][6*3 +4] += 2.0 * zdsps_12; // xyy yz (i=y, a=xyy)
        //this->pFSDS[m][6*3 +5] += 0.0;          // xyy zz (i=z, a=xyy)

        this->pFSDS[m][6*4 +0] +=       zdsps_40; // xyz xx (i=x, a=xyz)
        this->pFSDS[m][6*4 +1] +=       zdsps_41; // xyz xy (i=x, a=xyz)
        this->pFSDS[m][6*4 +2] +=       zdsps_42; // xyz xz (i=x, a=xyz)
        this->pFSDS[m][6*4 +3] +=       zdsps_21; // xyz yy (i=y, a=xyz)
        this->pFSDS[m][6*4 +4] +=       zdsps_22; // xyz yz (i=y, a=xyz)
        this->pFSDS[m][6*4 +5] +=       zdsps_12; // xyz zz (i=z, a=xyz)

        this->pFSDS[m][6*5 +0] +=       zdsps_50; // xzz xx (i=x, a=xzz)
        this->pFSDS[m][6*5 +1] +=       zdsps_51; // xzz xy (i=x, a=xzz)
        this->pFSDS[m][6*5 +2] +=       zdsps_52; // xzz xz (i=x, a=xzz)
        //this->pFSDS[m][6*5 +3] += 0.0;          // xzz yy (i=y, a=xzz)
        //this->pFSDS[m][6*5 +4] += 0.0;          // xzz yz (i=y, a=xzz)
        this->pFSDS[m][6*5 +5] += 2.0 * zdsps_22; // xzz zz (i=z, a=xzz)

        //this->pFSDS[m][6*6 +0] += 0.0;          // yyy xx (i=x, a=yyy)
        //this->pFSDS[m][6*6 +1] += 0.0;          // yyy xy (i=x, a=yyy)
        //this->pFSDS[m][6*6 +2] += 0.0;          // yyy xz (i=x, a=yyy)
        this->pFSDS[m][6*6 +3] += 3.0 * zdsps_31; // yyy yy (i=y, a=yyy)
        this->pFSDS[m][6*6 +4] += 3.0 * zdsps_32; // yyy yz (i=y, a=yyy)
        //this->pFSDS[m][6*6 +5] += 0.0;          // yyy zz (i=z, a=yyy)

        //this->pFSDS[m][6*7 +0] += 0.0;          // yyz xx (i=x, a=yyz)
        //this->pFSDS[m][6*7 +1] += 0.0;          // yyz xy (i=x, a=yyz)
        //this->pFSDS[m][6*7 +2] += 0.0;          // yyz xz (i=x, a=yyz)
        this->pFSDS[m][6*7 +3] += 2.0 * zdsps_41; // yyz yy (i=y, a=yyz)
        this->pFSDS[m][6*7 +4] += 2.0 * zdsps_42; // yyz yz (i=y, a=yyz)
        this->pFSDS[m][6*7 +5] +=       zdsps_32; // yyz zz (i=z, a=yyz)

        //this->pFSDS[m][6*8 +0] += 0.0;          // yzz xx (i=x, a=yzz)
        //this->pFSDS[m][6*8 +1] += 0.0;          // yzz xy (i=x, a=yzz)
        //this->pFSDS[m][6*8 +2] += 0.0;          // yzz xz (i=x, a=yzz)
        this->pFSDS[m][6*8 +3] +=       zdsps_51; // yzz yy (i=y, a=yzz)
        this->pFSDS[m][6*8 +4] +=       zdsps_52; // yzz yz (i=y, a=yzz)
        this->pFSDS[m][6*8 +5] += 2.0 * zdsps_42; // yzz zz (i=z, a=yzz)

        //this->pFSDS[m][6*9 +0] += 0.0;          // zzz xx (i=x, a=zzz)
        //this->pFSDS[m][6*9 +1] += 0.0;          // zzz xy (i=x, a=zzz)
        //this->pFSDS[m][6*9 +2] += 0.0;          // zzz xz (i=x, a=zzz)
        //this->pFSDS[m][6*9 +3] += 0.0;          // zzz yy (i=y, a=zzz)
        //this->pFSDS[m][6*9 +4] += 0.0;          // zzz yz (i=y, a=zzz)
        this->pFSDS[m][6*9 +5] += 3.0 * zdsps_52; // zzz zz (i=z, a=zzz)
    }
}


void DfTEI::primitiveFSFS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 10; ++i) {
            this->pFSFS[m][10*i +0] = QCx * this->pFSDS[m][6*i +0] +WQx * this->pFSDS[m1][6*i +0]; // xxx
            this->pFSFS[m][10*i +1] = QCx * this->pFSDS[m][6*i +1] +WQx * this->pFSDS[m1][6*i +1]; // xxy
            this->pFSFS[m][10*i +2] = QCx * this->pFSDS[m][6*i +2] +WQx * this->pFSDS[m1][6*i +2]; // xxz
            this->pFSFS[m][10*i +3] = QCx * this->pFSDS[m][6*i +3] +WQx * this->pFSDS[m1][6*i +3]; // xyy
            this->pFSFS[m][10*i +4] = QCx * this->pFSDS[m][6*i +4] +WQx * this->pFSDS[m1][6*i +4]; // xyz
            this->pFSFS[m][10*i +5] = QCx * this->pFSDS[m][6*i +5] +WQx * this->pFSDS[m1][6*i +5]; // xzz
            this->pFSFS[m][10*i +6] = QCy * this->pFSDS[m][6*i +3] +WQy * this->pFSDS[m1][6*i +3]; // yyy
            this->pFSFS[m][10*i +7] = QCy * this->pFSDS[m][6*i +4] +WQy * this->pFSDS[m1][6*i +4]; // yyz
            this->pFSFS[m][10*i +8] = QCy * this->pFSDS[m][6*i +5] +WQy * this->pFSDS[m1][6*i +5]; // yzz
            this->pFSFS[m][10*i +9] = QCz * this->pFSDS[m][6*i +5] +WQz * this->pFSDS[m1][6*i +5]; // zzz
        }

        for (int i = 0; i < 10; ++i) {
            const double zfsps_0 = eta2 * (this->pFSPS[m][3*i +0] -re * this->pFSPS[m1][3*i +0]); // x
            const double zfsps_1 = eta2 * (this->pFSPS[m][3*i +1] -re * this->pFSPS[m1][3*i +1]); // y
            const double zfsps_2 = eta2 * (this->pFSPS[m][3*i +2] -re * this->pFSPS[m1][3*i +2]); // z

            this->pFSFS[m][10*i +0] += 2.0 * zfsps_0; // c=xx, i=x
            this->pFSFS[m][10*i +1] +=       zfsps_1; // c=xy, i=x
            this->pFSFS[m][10*i +2] +=       zfsps_2; // c=xz, i=x
            //this->pFSFS[m][10*i +3] += 0.0;         // c=yy, i=x
            //this->pFSFS[m][10*i +4] += 0.0;         // c=yz, i=x
            //this->pFSFS[m][10*i +5] += 0.0;         // c=zz, i=x
            this->pFSFS[m][10*i +6] += 2.0 * zfsps_1; // c=yy, i=y
            this->pFSFS[m][10*i +7] +=       zfsps_2; // c=yz, i=y
            //this->pFSFS[m][10*i +8] += 0.0;         // c=zz, i=y
            this->pFSFS[m][10*i +9] += 2.0 * zfsps_2; // c=zz, i=z
        }

        const double zdsds_00 = ze2 * this->pDSDS[m1][6*0 +0]; // xx xx
        const double zdsds_01 = ze2 * this->pDSDS[m1][6*0 +1]; // xx xy
        const double zdsds_02 = ze2 * this->pDSDS[m1][6*0 +2]; // xx xz
        const double zdsds_03 = ze2 * this->pDSDS[m1][6*0 +3]; // xx yy
        const double zdsds_04 = ze2 * this->pDSDS[m1][6*0 +4]; // xx yz
        const double zdsds_05 = ze2 * this->pDSDS[m1][6*0 +5]; // xx zz
        const double zdsds_10 = ze2 * this->pDSDS[m1][6*1 +0]; // xy
        const double zdsds_11 = ze2 * this->pDSDS[m1][6*1 +1]; // xy
        const double zdsds_12 = ze2 * this->pDSDS[m1][6*1 +2]; // xy
        const double zdsds_13 = ze2 * this->pDSDS[m1][6*1 +3]; // xy
        const double zdsds_14 = ze2 * this->pDSDS[m1][6*1 +4]; // xy
        const double zdsds_15 = ze2 * this->pDSDS[m1][6*1 +5]; // xy
        const double zdsds_20 = ze2 * this->pDSDS[m1][6*2 +0]; // xz
        const double zdsds_21 = ze2 * this->pDSDS[m1][6*2 +1]; // xz
        const double zdsds_22 = ze2 * this->pDSDS[m1][6*2 +2]; // xz
        const double zdsds_23 = ze2 * this->pDSDS[m1][6*2 +3]; // xz
        const double zdsds_24 = ze2 * this->pDSDS[m1][6*2 +4]; // xz
        const double zdsds_25 = ze2 * this->pDSDS[m1][6*2 +5]; // xz
        const double zdsds_30 = ze2 * this->pDSDS[m1][6*3 +0]; // yy
        const double zdsds_31 = ze2 * this->pDSDS[m1][6*3 +1]; // yy
        const double zdsds_32 = ze2 * this->pDSDS[m1][6*3 +2]; // yy
        const double zdsds_33 = ze2 * this->pDSDS[m1][6*3 +3]; // yy
        const double zdsds_34 = ze2 * this->pDSDS[m1][6*3 +4]; // yy
        const double zdsds_35 = ze2 * this->pDSDS[m1][6*3 +5]; // yy
        const double zdsds_40 = ze2 * this->pDSDS[m1][6*4 +0]; // yz
        const double zdsds_41 = ze2 * this->pDSDS[m1][6*4 +1]; // yz
        const double zdsds_42 = ze2 * this->pDSDS[m1][6*4 +2]; // yz
        const double zdsds_43 = ze2 * this->pDSDS[m1][6*4 +3]; // yz
        const double zdsds_44 = ze2 * this->pDSDS[m1][6*4 +4]; // yz
        const double zdsds_45 = ze2 * this->pDSDS[m1][6*4 +5]; // yz
        const double zdsds_50 = ze2 * this->pDSDS[m1][6*5 +0]; // zz
        const double zdsds_51 = ze2 * this->pDSDS[m1][6*5 +1]; // zz
        const double zdsds_52 = ze2 * this->pDSDS[m1][6*5 +2]; // zz
        const double zdsds_53 = ze2 * this->pDSDS[m1][6*5 +3]; // zz
        const double zdsds_54 = ze2 * this->pDSDS[m1][6*5 +4]; // zz
        const double zdsds_55 = ze2 * this->pDSDS[m1][6*5 +5]; // zz

        this->pFSFS[m][10*0 +0] += 3.0 * zdsds_00; // i=x, a=xxx
        this->pFSFS[m][10*0 +1] += 3.0 * zdsds_01; // i=x
        this->pFSFS[m][10*0 +2] += 3.0 * zdsds_02; // xxx xxz (i=x, a=xxx)
        this->pFSFS[m][10*0 +3] += 3.0 * zdsds_03; // xxx xyy (i=x, a=xxx)
        this->pFSFS[m][10*0 +4] += 3.0 * zdsds_04; // xxx xyz (i=x, a=xxx)
        this->pFSFS[m][10*0 +5] += 3.0 * zdsds_05; // xxx xzz (i=x, a=xxx)
        //this->pFSFS[m][10*0 +6] += 0.0;          // xxx yyy (i=y, a=xxx)
        //this->pFSFS[m][10*0 +7] += 0.0;          // xxx yyz (i=y, a=xxx)
        //this->pFSFS[m][10*0 +8] += 0.0;          // xxx yzz (i=y, a=xxx)
        //this->pFSFS[m][10*0 +9] += 0.0;          // xxx zzz (i=z, a=xxx)

        this->pFSFS[m][10*1 +0] += 2.0 * zdsds_10; // xxy xxx (i=x, a=xxy)
        this->pFSFS[m][10*1 +1] += 2.0 * zdsds_11; // xxy xxy (i=x, a=xxy)
        this->pFSFS[m][10*1 +2] += 2.0 * zdsds_12; // xxy xxz (i=x, a=xxy)
        this->pFSFS[m][10*1 +3] += 2.0 * zdsds_13; // xxy xyy (i=x, a=xxy)
        this->pFSFS[m][10*1 +4] += 2.0 * zdsds_14; // xxy xyz (i=x, a=xxy)
        this->pFSFS[m][10*1 +5] += 2.0 * zdsds_15; // xxy xzz (i=x, a=xxy)
        this->pFSFS[m][10*1 +6] +=       zdsds_03; // xxy yyy (i=y, a=xxy)
        this->pFSFS[m][10*1 +7] +=       zdsds_04; // xxy yyz (i=y, a=xxy)
        this->pFSFS[m][10*1 +8] +=       zdsds_05; // xxy yzz (i=y, a=xxy)
        //this->pFSFS[m][10*1 +9] += 0.0;          // xxy zzz (i=z, a=xxy)

        this->pFSFS[m][10*2 +0] += 2.0 * zdsds_20; // xxz xxx (i=x, a=xxz)
        this->pFSFS[m][10*2 +1] += 2.0 * zdsds_21; // xxz xxy (i=x, a=xxz)
        this->pFSFS[m][10*2 +2] += 2.0 * zdsds_22; // xxz xxz (i=x, a=xxz)
        this->pFSFS[m][10*2 +3] += 2.0 * zdsds_23; // xxz xyy (i=x, a=xxz)
        this->pFSFS[m][10*2 +4] += 2.0 * zdsds_24; // xxz xyz (i=x, a=xxz)
        this->pFSFS[m][10*2 +5] += 2.0 * zdsds_25; // xxz xzz (i=x, a=xxz)
        //this->pFSFS[m][10*2 +6] += 0.0;          // xxz yyy (i=y, a=xxz)
        //this->pFSFS[m][10*2 +7] += 0.0;          // xxz yyz (i=y, a=xxz)
        //this->pFSFS[m][10*2 +8] += 0.0;          // xxz yzz (i=y, a=xxz)
        this->pFSFS[m][10*2 +9] +=       zdsds_05; // xxz zzz (i=z, a=xxz)

        this->pFSFS[m][10*3 +0] +=       zdsds_30; // xyy xxx (i=x, a=xyy)
        this->pFSFS[m][10*3 +1] +=       zdsds_31; // xyy xxy (i=x, a=xyy)
        this->pFSFS[m][10*3 +2] +=       zdsds_32; // xyy xxz (i=x, a=xyy)
        this->pFSFS[m][10*3 +3] +=       zdsds_33; // xyy xyy (i=x, a=xyy)
        this->pFSFS[m][10*3 +4] +=       zdsds_34; // xyy xyz (i=x, a=xyy)
        this->pFSFS[m][10*3 +5] +=       zdsds_35; // xyy xzz (i=x, a=xyy)
        this->pFSFS[m][10*3 +6] += 2.0 * zdsds_13; // xyy yyy (i=y, a=xyy)
        this->pFSFS[m][10*3 +7] += 2.0 * zdsds_14; // xyy yyz (i=y, a=xyy)
        this->pFSFS[m][10*3 +8] += 2.0 * zdsds_15; // xyy yzz (i=y, a=xyy)
        //this->pFSFS[m][10*3 +9] += 0.0;          // xyy zzz (i=z, a=xyy)

        this->pFSFS[m][10*4 +0] +=       zdsds_40; // xyz xxx (i=x, a=xyz)
        this->pFSFS[m][10*4 +1] +=       zdsds_41; // xyz xxy (i=x, a=xyz)
        this->pFSFS[m][10*4 +2] +=       zdsds_42; // xyz xxz (i=x, a=xyz)
        this->pFSFS[m][10*4 +3] +=       zdsds_43; // xyz xyy (i=x, a=xyz)
        this->pFSFS[m][10*4 +4] +=       zdsds_44; // xyz xyz (i=x, a=xyz)
        this->pFSFS[m][10*4 +5] +=       zdsds_45; // xyz xzz (i=x, a=xyz)
        this->pFSFS[m][10*4 +6] +=       zdsds_23; // xyz yyy (i=y, a=xyz)
        this->pFSFS[m][10*4 +7] +=       zdsds_24; // xyz yyz (i=y, a=xyz)
        this->pFSFS[m][10*4 +8] +=       zdsds_25; // xyz yzz (i=y, a=xyz)
        this->pFSFS[m][10*4 +9] +=       zdsds_15; // xyz zzz (i=z, a=xyz)

        this->pFSFS[m][10*5 +0] +=       zdsds_50; // xzz xxx (i=x, a=xzz)
        this->pFSFS[m][10*5 +1] +=       zdsds_51; // xzz xxy (i=x, a=xzz)
        this->pFSFS[m][10*5 +2] +=       zdsds_52; // xzz xxz (i=x, a=xzz)
        this->pFSFS[m][10*5 +3] +=       zdsds_53; // xzz xyy (i=x, a=xzz)
        this->pFSFS[m][10*5 +4] +=       zdsds_54; // xzz xyz (i=x, a=xzz)
        this->pFSFS[m][10*5 +5] +=       zdsds_55; // xzz xzz (i=x, a=xzz)
        //this->pFSFS[m][10*5 +6] += 0.0;          // xzz yyy (i=y, a=xzz)
        //this->pFSFS[m][10*5 +7] += 0.0;          // xzz yyz (i=y, a=xzz)
        //this->pFSFS[m][10*5 +8] += 0.0;          // xzz yzz (i=y, a=xzz)
        this->pFSFS[m][10*5 +9] += 2.0 * zdsds_25; // xzz zzz (i=z, a=xzz)

        //this->pFSFS[m][10*6 +0] += 0.0;          // yyy xxx (i=x, a=yyy)
        //this->pFSFS[m][10*6 +1] += 0.0;          // yyy xxy (i=x, a=yyy)
        //this->pFSFS[m][10*6 +2] += 0.0;          // yyy xxz (i=x, a=yyy)
        //this->pFSFS[m][10*6 +3] += 0.0;          // yyy xyy (i=x, a=yyy)
        //this->pFSFS[m][10*6 +4] += 0.0;          // yyy xyz (i=x, a=yyy)
        //this->pFSFS[m][10*6 +5] += 0.0;          // yyy xzz (i=x, a=yyy)
        this->pFSFS[m][10*6 +6] += 3.0 * zdsds_33; // yyy yyy (i=y, a=yyy)
        this->pFSFS[m][10*6 +7] += 3.0 * zdsds_34; // yyy yyz (i=y, a=yyy)
        this->pFSFS[m][10*6 +8] += 3.0 * zdsds_35; // yyy yzz (i=y, a=yyy)
        //this->pFSFS[m][10*6 +9] += 0.0;          // yyy zzz (i=z, a=yyy)

        //this->pFSFS[m][10*7 +0] += 0.0;          // yyz xxx (i=x, a=yyz)
        //this->pFSFS[m][10*7 +1] += 0.0;          // yyz xxy (i=x, a=yyz)
        //this->pFSFS[m][10*7 +2] += 0.0;          // yyz xxz (i=x, a=yyz)
        //this->pFSFS[m][10*7 +3] += 0.0;          // yyz xyy (i=x, a=yyz)
        //this->pFSFS[m][10*7 +4] += 0.0;          // yyz xyz (i=x, a=yyz)
        //this->pFSFS[m][10*7 +5] += 0.0;          // yyz xzz (i=x, a=yyz)
        this->pFSFS[m][10*7 +6] += 2.0 * zdsds_43; // yyz yyy (i=y, a=yyz)
        this->pFSFS[m][10*7 +7] += 2.0 * zdsds_44; // yyz yyz (i=y, a=yyz)
        this->pFSFS[m][10*7 +8] += 2.0 * zdsds_45; // yyz yzz (i=y, a=yyz)
        this->pFSFS[m][10*7 +9] +=       zdsds_35; // yyz zzz (i=z, a=yyz)

        //this->pFSFS[m][10*8 +0] += 0.0;          // yzz xxx (i=x, a=yzz)
        //this->pFSFS[m][10*8 +1] += 0.0;          // yzz xxy (i=x, a=yzz)
        //this->pFSFS[m][10*8 +2] += 0.0;          // yzz xxz (i=x, a=yzz)
        //this->pFSFS[m][10*8 +3] += 0.0;          // yzz xyy (i=x, a=yzz)
        //this->pFSFS[m][10*8 +4] += 0.0;          // yzz xyz (i=x, a=yzz)
        //this->pFSFS[m][10*8 +5] += 0.0;          // yzz xzz (i=x, a=yzz)
        this->pFSFS[m][10*8 +6] +=       zdsds_53; // yzz yyy (i=y, a=yzz)
        this->pFSFS[m][10*8 +7] +=       zdsds_54; // yzz yyz (i=y, a=yzz)
        this->pFSFS[m][10*8 +8] +=       zdsds_55; // yzz yzz (i=y, a=yzz)
        this->pFSFS[m][10*8 +9] += 2.0 * zdsds_45; // yzz zzz (i=z, a=yzz)

        //this->pFSFS[m][10*9 +0] += 0.0;          // zzz xxx (i=x, a=zzz)
        //this->pFSFS[m][10*9 +1] += 0.0;          // zzz xxy (i=x, a=zzz)
        //this->pFSFS[m][10*9 +2] += 0.0;          // zzz xxz (i=x, a=zzz)
        //this->pFSFS[m][10*9 +3] += 0.0;          // zzz xyy (i=x, a=zzz)
        //this->pFSFS[m][10*9 +4] += 0.0;          // zzz xyz (i=x, a=zzz)
        //this->pFSFS[m][10*9 +5] += 0.0;          // zzz xzz (i=x, a=zzz)
        //this->pFSFS[m][10*9 +6] += 0.0;          // zzz yyy (i=y, a=zzz)
        //this->pFSFS[m][10*9 +7] += 0.0;          // zzz yyz (i=y, a=zzz)
        //this->pFSFS[m][10*9 +8] += 0.0;          // zzz yzz (i=y, a=zzz)
        this->pFSFS[m][10*9 +9] += 3.0 * zdsds_55; // zzz zzz (i=z, a=zzz)
    }
}


void DfTEI::primitiveFSGS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 10; ++i) {
            this->pFSGS[m][15*i + 0] = QCx * this->pFSFS[m][10*i +0] +WQx * this->pFSFS[m1][10*i +0]; // xxx
            this->pFSGS[m][15*i + 1] = QCx * this->pFSFS[m][10*i +1] +WQx * this->pFSFS[m1][10*i +1]; // xxy
            this->pFSGS[m][15*i + 2] = QCx * this->pFSFS[m][10*i +2] +WQx * this->pFSFS[m1][10*i +2]; // xxz
            this->pFSGS[m][15*i + 3] = QCx * this->pFSFS[m][10*i +3] +WQx * this->pFSFS[m1][10*i +3]; // xyy
            this->pFSGS[m][15*i + 4] = QCx * this->pFSFS[m][10*i +4] +WQx * this->pFSFS[m1][10*i +4]; // xyz
            this->pFSGS[m][15*i + 5] = QCx * this->pFSFS[m][10*i +5] +WQx * this->pFSFS[m1][10*i +5]; // xzz
            this->pFSGS[m][15*i + 6] = QCx * this->pFSFS[m][10*i +6] +WQx * this->pFSFS[m1][10*i +6]; // yyy
            this->pFSGS[m][15*i + 7] = QCx * this->pFSFS[m][10*i +7] +WQx * this->pFSFS[m1][10*i +7]; // yyz
            this->pFSGS[m][15*i + 8] = QCx * this->pFSFS[m][10*i +8] +WQx * this->pFSFS[m1][10*i +8]; // yzz
            this->pFSGS[m][15*i + 9] = QCx * this->pFSFS[m][10*i +9] +WQx * this->pFSFS[m1][10*i +9]; // zzz
            this->pFSGS[m][15*i +10] = QCy * this->pFSFS[m][10*i +6] +WQy * this->pFSFS[m1][10*i +6]; // xzz
            this->pFSGS[m][15*i +11] = QCy * this->pFSFS[m][10*i +7] +WQy * this->pFSFS[m1][10*i +7]; // yyy
            this->pFSGS[m][15*i +12] = QCy * this->pFSFS[m][10*i +8] +WQy * this->pFSFS[m1][10*i +8]; // yyz
            this->pFSGS[m][15*i +13] = QCy * this->pFSFS[m][10*i +9] +WQy * this->pFSFS[m1][10*i +9]; // yzz
            this->pFSGS[m][15*i +14] = QCz * this->pFSFS[m][10*i +9] +WQz * this->pFSFS[m1][10*i +9]; // zzz
        }

        for (int i = 0; i < 10; ++i) {
            const double zfsds_0 = eta2 * (this->pFSDS[m][6*i +0] -re * this->pFSDS[m1][6*i +0]); // xx
            const double zfsds_1 = eta2 * (this->pFSDS[m][6*i +1] -re * this->pFSDS[m1][6*i +1]); // xy
            const double zfsds_2 = eta2 * (this->pFSDS[m][6*i +2] -re * this->pFSDS[m1][6*i +2]); // xz
            const double zfsds_3 = eta2 * (this->pFSDS[m][6*i +3] -re * this->pFSDS[m1][6*i +3]); // yy
            const double zfsds_4 = eta2 * (this->pFSDS[m][6*i +4] -re * this->pFSDS[m1][6*i +4]); // yz
            const double zfsds_5 = eta2 * (this->pFSDS[m][6*i +5] -re * this->pFSDS[m1][6*i +5]); // zz

            this->pFSGS[m][15*i + 0] += 3.0 * zfsds_0;  // c=xxx, i=x
            this->pFSGS[m][15*i + 1] += 2.0 * zfsds_1;  // c=xxy, i=x
            this->pFSGS[m][15*i + 2] += 2.0 * zfsds_2;  // c=xxz, i=x
            this->pFSGS[m][15*i + 3] +=       zfsds_3;  // c=xyy, i=x
            this->pFSGS[m][15*i + 4] +=       zfsds_4;  // c=xyz, i=x
            this->pFSGS[m][15*i + 5] +=       zfsds_5;  // c=xzz, i=x
            //this->pFSGS[m][15*i + 6] += 0.0;         // c=yyy, i=x
            //this->pFSGS[m][15*i + 7] += 0.0;         // c=yyz, i=x
            //this->pFSGS[m][15*i + 8] += 0.0;         // c=yzz, i=x
            //this->pFSGS[m][15*i + 9] += 0.0;         // c=zzz, i=x
            this->pFSGS[m][15*i +10] += 3.0 * zfsds_3; // c=yyy, i=y
            this->pFSGS[m][15*i +11] += 2.0 * zfsds_4; // c=yyz, i=y
            this->pFSGS[m][15*i +12] +=       zfsds_5; // c=yzz, i=y
            //this->pFSGS[m][15*i +13] += 0.0;         // c=zzz, i=y
            this->pFSGS[m][15*i +14] += 3.0 * zfsds_5; // c=zzz, i=z
        }

        const double zdsfs_00 = ze2 * this->pDSFS[m1][10*0 +0]; // xx xxx
        const double zdsfs_01 = ze2 * this->pDSFS[m1][10*0 +1]; // xx xxy
        const double zdsfs_02 = ze2 * this->pDSFS[m1][10*0 +2]; // xx xxz
        const double zdsfs_03 = ze2 * this->pDSFS[m1][10*0 +3]; // xx xyy
        const double zdsfs_04 = ze2 * this->pDSFS[m1][10*0 +4]; // xx xyz
        const double zdsfs_05 = ze2 * this->pDSFS[m1][10*0 +5]; // xx xzz
        const double zdsfs_06 = ze2 * this->pDSFS[m1][10*0 +6]; // xx yyy
        const double zdsfs_07 = ze2 * this->pDSFS[m1][10*0 +7]; // xx yyz
        const double zdsfs_08 = ze2 * this->pDSFS[m1][10*0 +8]; // xx yzz
        const double zdsfs_09 = ze2 * this->pDSFS[m1][10*0 +9]; // xx zzz

        const double zdsfs_10 = ze2 * this->pDSFS[m1][10*1 +0]; // xy xxx
        const double zdsfs_11 = ze2 * this->pDSFS[m1][10*1 +1]; // xy xxy
        const double zdsfs_12 = ze2 * this->pDSFS[m1][10*1 +2]; // xy xxz
        const double zdsfs_13 = ze2 * this->pDSFS[m1][10*1 +3]; // xy xyy
        const double zdsfs_14 = ze2 * this->pDSFS[m1][10*1 +4]; // xy xyz
        const double zdsfs_15 = ze2 * this->pDSFS[m1][10*1 +5]; // xy xzz
        const double zdsfs_16 = ze2 * this->pDSFS[m1][10*1 +6]; // xy yyy
        const double zdsfs_17 = ze2 * this->pDSFS[m1][10*1 +7]; // xy yyz
        const double zdsfs_18 = ze2 * this->pDSFS[m1][10*1 +8]; // xy yzz
        const double zdsfs_19 = ze2 * this->pDSFS[m1][10*1 +9]; // xy zzz

        const double zdsfs_20 = ze2 * this->pDSFS[m1][10*2 +0]; // xz xxx
        const double zdsfs_21 = ze2 * this->pDSFS[m1][10*2 +1]; // xz xxy
        const double zdsfs_22 = ze2 * this->pDSFS[m1][10*2 +2]; // xz xxz
        const double zdsfs_23 = ze2 * this->pDSFS[m1][10*2 +3]; // xz xyy
        const double zdsfs_24 = ze2 * this->pDSFS[m1][10*2 +4]; // xz xyz
        const double zdsfs_25 = ze2 * this->pDSFS[m1][10*2 +5]; // xz xzz
        const double zdsfs_26 = ze2 * this->pDSFS[m1][10*2 +6]; // xz yyy
        const double zdsfs_27 = ze2 * this->pDSFS[m1][10*2 +7]; // xz yyz
        const double zdsfs_28 = ze2 * this->pDSFS[m1][10*2 +8]; // xz yzz
        const double zdsfs_29 = ze2 * this->pDSFS[m1][10*2 +9]; // xz zzz

        const double zdsfs_30 = ze2 * this->pDSFS[m1][10*3 +0]; // yy xxx
        const double zdsfs_31 = ze2 * this->pDSFS[m1][10*3 +1]; // yy xxy
        const double zdsfs_32 = ze2 * this->pDSFS[m1][10*3 +2]; // yy xxz
        const double zdsfs_33 = ze2 * this->pDSFS[m1][10*3 +3]; // yy xyy
        const double zdsfs_34 = ze2 * this->pDSFS[m1][10*3 +4]; // yy xyz
        const double zdsfs_35 = ze2 * this->pDSFS[m1][10*3 +5]; // yy xzz
        const double zdsfs_36 = ze2 * this->pDSFS[m1][10*3 +6]; // yy yyy
        const double zdsfs_37 = ze2 * this->pDSFS[m1][10*3 +7]; // yy yyz
        const double zdsfs_38 = ze2 * this->pDSFS[m1][10*3 +8]; // yy yzz
        const double zdsfs_39 = ze2 * this->pDSFS[m1][10*3 +9]; // yy zzz

        const double zdsfs_40 = ze2 * this->pDSFS[m1][10*4 +0]; // yz xxx
        const double zdsfs_41 = ze2 * this->pDSFS[m1][10*4 +1]; // yz xxy
        const double zdsfs_42 = ze2 * this->pDSFS[m1][10*4 +2]; // yz xxz
        const double zdsfs_43 = ze2 * this->pDSFS[m1][10*4 +3]; // yz xyy
        const double zdsfs_44 = ze2 * this->pDSFS[m1][10*4 +4]; // yz xyz
        const double zdsfs_45 = ze2 * this->pDSFS[m1][10*4 +5]; // yz xzz
        const double zdsfs_46 = ze2 * this->pDSFS[m1][10*4 +6]; // yz yyy
        const double zdsfs_47 = ze2 * this->pDSFS[m1][10*4 +7]; // yz yyz
        const double zdsfs_48 = ze2 * this->pDSFS[m1][10*4 +8]; // yz yzz
        const double zdsfs_49 = ze2 * this->pDSFS[m1][10*4 +9]; // yz zzz

        const double zdsfs_50 = ze2 * this->pDSFS[m1][10*5 +0]; // yz xxx
        const double zdsfs_51 = ze2 * this->pDSFS[m1][10*5 +1]; // yz xxy
        const double zdsfs_52 = ze2 * this->pDSFS[m1][10*5 +2]; // yz xxz
        const double zdsfs_53 = ze2 * this->pDSFS[m1][10*5 +3]; // yz xyy
        const double zdsfs_54 = ze2 * this->pDSFS[m1][10*5 +4]; // yz xyz
        const double zdsfs_55 = ze2 * this->pDSFS[m1][10*5 +5]; // yz xzz
        const double zdsfs_56 = ze2 * this->pDSFS[m1][10*5 +6]; // yz yyy
        const double zdsfs_57 = ze2 * this->pDSFS[m1][10*5 +7]; // yz yyz
        const double zdsfs_58 = ze2 * this->pDSFS[m1][10*5 +8]; // yz yzz
        const double zdsfs_59 = ze2 * this->pDSFS[m1][10*5 +9]; // yz zzz

        this->pFSGS[m][15*0 + 0] += 3.0 * zdsfs_00; // i=x a=xxx
        this->pFSGS[m][15*0 + 1] += 3.0 * zdsfs_01; // xxx xxxy (i=x)
        this->pFSGS[m][15*0 + 2] += 3.0 * zdsfs_02; // xxx xxxz (i=x)
        this->pFSGS[m][15*0 + 3] += 3.0 * zdsfs_03; // xxx xxyy (i=x)
        this->pFSGS[m][15*0 + 4] += 3.0 * zdsfs_04; // xxx xxyz (i=x)
        this->pFSGS[m][15*0 + 5] += 3.0 * zdsfs_05; // xxx xxzz (i=x)
        this->pFSGS[m][15*0 + 6] += 3.0 * zdsfs_06; // xxx xyyy (i=x)
        this->pFSGS[m][15*0 + 7] += 3.0 * zdsfs_07; // xxx xyyz (i=x)
        this->pFSGS[m][15*0 + 8] += 3.0 * zdsfs_08; // xxx xyzz (i=x)
        this->pFSGS[m][15*0 + 9] += 3.0 * zdsfs_09; // xxx xzzz (i=x)
        //this->pFSGS[m][15*0 +10] += 0.0;          // xxx yyyy (i=y)
        //this->pFSGS[m][15*0 +11] += 0.0;          // xxx yyyz (i=y)
        //this->pFSGS[m][15*0 +12] += 0.0;          // xxx yyzz (i=y)
        //this->pFSGS[m][15*0 +13] += 0.0;          // xxx yzzz (i=y)
        //this->pFSGS[m][15*0 +14] += 0.0;          // xxx zzzz (i=z)

        this->pFSGS[m][15*1 + 0] += 2.0 * zdsfs_10; // i=x a=xxy
        this->pFSGS[m][15*1 + 1] += 2.0 * zdsfs_11; // xxy xxxy (i=x)
        this->pFSGS[m][15*1 + 2] += 2.0 * zdsfs_12; // xxy xxxz (i=x)
        this->pFSGS[m][15*1 + 3] += 2.0 * zdsfs_13; // xxy xxyy (i=x)
        this->pFSGS[m][15*1 + 4] += 2.0 * zdsfs_14; // xxy xxyz (i=x)
        this->pFSGS[m][15*1 + 5] += 2.0 * zdsfs_15; // xxy xxzz (i=x)
        this->pFSGS[m][15*1 + 6] += 2.0 * zdsfs_16; // xxy xyyy (i=x)
        this->pFSGS[m][15*1 + 7] += 2.0 * zdsfs_17; // xxy xyyz (i=x)
        this->pFSGS[m][15*1 + 8] += 2.0 * zdsfs_18; // xxy xyzz (i=x)
        this->pFSGS[m][15*1 + 9] += 2.0 * zdsfs_19; // xxy xzzz (i=x)
        this->pFSGS[m][15*1 +10] +=       zdsfs_06; // xxy yyyy (i=y)
        this->pFSGS[m][15*1 +11] +=       zdsfs_07; // xxy yyyz (i=y)
        this->pFSGS[m][15*1 +12] +=       zdsfs_08; // xxy yyzz (i=y)
        this->pFSGS[m][15*1 +13] +=       zdsfs_09; // xxy yzzz (i=y)
        //this->pFSGS[m][15*1 +14] += 0.0;          // xxy zzzz (i=z)

        this->pFSGS[m][15*2 + 0] += 2.0 * zdsfs_20; // i=x a=xxz
        this->pFSGS[m][15*2 + 1] += 2.0 * zdsfs_21; // xxz xxxy (i=x)
        this->pFSGS[m][15*2 + 2] += 2.0 * zdsfs_22; // xxz xxxz (i=x)
        this->pFSGS[m][15*2 + 3] += 2.0 * zdsfs_23; // xxz xxyy (i=x)
        this->pFSGS[m][15*2 + 4] += 2.0 * zdsfs_24; // xxz xxyz (i=x)
        this->pFSGS[m][15*2 + 5] += 2.0 * zdsfs_25; // xxz xxzz (i=x)
        this->pFSGS[m][15*2 + 6] += 2.0 * zdsfs_26; // xxz xyyy (i=x)
        this->pFSGS[m][15*2 + 7] += 2.0 * zdsfs_27; // xxz xyyz (i=x)
        this->pFSGS[m][15*2 + 8] += 2.0 * zdsfs_28; // xxz xyzz (i=x)
        this->pFSGS[m][15*2 + 9] += 2.0 * zdsfs_29; // xxz xzzz (i=x)
        //this->pFSGS[m][15*2 +10] += 0.0;          // xxz yyyy (i=y)
        //this->pFSGS[m][15*2 +11] += 0.0;          // xxz yyyz (i=y)
        //this->pFSGS[m][15*2 +12] += 0.0;          // xxz yyzz (i=y)
        //this->pFSGS[m][15*2 +13] += 0.0;          // xxz yzzz (i=y)
        this->pFSGS[m][15*2 +14] +=       zdsfs_09; // xxz zzzz (i=z)

        this->pFSGS[m][15*3 + 0] +=       zdsfs_30; // i=x a=xyy
        this->pFSGS[m][15*3 + 1] +=       zdsfs_31; // xyy xxxy (i=x)
        this->pFSGS[m][15*3 + 2] +=       zdsfs_32; // xyy xxxz (i=x)
        this->pFSGS[m][15*3 + 3] +=       zdsfs_33; // xyy xxyy (i=x)
        this->pFSGS[m][15*3 + 4] +=       zdsfs_34; // xyy xxyz (i=x)
        this->pFSGS[m][15*3 + 5] +=       zdsfs_35; // xyy xxzz (i=x)
        this->pFSGS[m][15*3 + 6] +=       zdsfs_36; // xyy xyyy (i=x)
        this->pFSGS[m][15*3 + 7] +=       zdsfs_37; // xyy xyyz (i=x)
        this->pFSGS[m][15*3 + 8] +=       zdsfs_38; // xyy xyzz (i=x)
        this->pFSGS[m][15*3 + 9] +=       zdsfs_39; // xyy xzzz (i=x)
        this->pFSGS[m][15*3 +10] += 2.0 * zdsfs_16; // xyy yyyy (i=y)
        this->pFSGS[m][15*3 +11] += 2.0 * zdsfs_17; // xyy yyyz (i=y)
        this->pFSGS[m][15*3 +12] += 2.0 * zdsfs_18; // xyy yyzz (i=y)
        this->pFSGS[m][15*3 +13] += 2.0 * zdsfs_19; // xyy yzzz (i=y)
        //this->pFSGS[m][15*3 +14] += 0.0;          // xyy zzzz (i=z)

        this->pFSGS[m][15*4 + 0] +=       zdsfs_40; // i=x a=xyz
        this->pFSGS[m][15*4 + 1] +=       zdsfs_41; // xyz xxxy (i=x)
        this->pFSGS[m][15*4 + 2] +=       zdsfs_42; // xyz xxxz (i=x)
        this->pFSGS[m][15*4 + 3] +=       zdsfs_43; // xyz xxyy (i=x)
        this->pFSGS[m][15*4 + 4] +=       zdsfs_44; // xyz xxyz (i=x)
        this->pFSGS[m][15*4 + 5] +=       zdsfs_45; // xyz xxzz (i=x)
        this->pFSGS[m][15*4 + 6] +=       zdsfs_46; // xyz xyyy (i=x)
        this->pFSGS[m][15*4 + 7] +=       zdsfs_47; // xyz xyyz (i=x)
        this->pFSGS[m][15*4 + 8] +=       zdsfs_48; // xyz xyzz (i=x)
        this->pFSGS[m][15*4 + 9] +=       zdsfs_49; // xyz xzzz (i=x)
        this->pFSGS[m][15*4 +10] +=       zdsfs_26; // xyz yyyy (i=y)
        this->pFSGS[m][15*4 +11] +=       zdsfs_27; // xyz yyyz (i=y)
        this->pFSGS[m][15*4 +12] +=       zdsfs_28; // xyz yyzz (i=y)
        this->pFSGS[m][15*4 +13] +=       zdsfs_29; // xyz yzzz (i=y)
        this->pFSGS[m][15*4 +14] +=       zdsfs_19; // xyz zzzz (i=z)

        this->pFSGS[m][15*5 + 0] +=       zdsfs_50; // i=x a=xzz
        this->pFSGS[m][15*5 + 1] +=       zdsfs_51; // xzz xxxy (i=x)
        this->pFSGS[m][15*5 + 2] +=       zdsfs_52; // xzz xxxz (i=x)
        this->pFSGS[m][15*5 + 3] +=       zdsfs_53; // xzz xxyy (i=x)
        this->pFSGS[m][15*5 + 4] +=       zdsfs_54; // xzz xxyz (i=x)
        this->pFSGS[m][15*5 + 5] +=       zdsfs_55; // xzz xxzz (i=x)
        this->pFSGS[m][15*5 + 6] +=       zdsfs_56; // xzz xyyy (i=x)
        this->pFSGS[m][15*5 + 7] +=       zdsfs_57; // xzz xyyz (i=x)
        this->pFSGS[m][15*5 + 8] +=       zdsfs_58; // xzz xyzz (i=x)
        this->pFSGS[m][15*5 + 9] +=       zdsfs_59; // xzz xzzz (i=x)
        //this->pFSGS[m][15*5 +10] += 0.0;          // xzz yyyy (i=y)
        //this->pFSGS[m][15*5 +11] += 0.0;          // xzz yyyz (i=y)
        //this->pFSGS[m][15*5 +12] += 0.0;          // xzz yyzz (i=y)
        //this->pFSGS[m][15*5 +13] += 0.0;          // xzz yzzz (i=y)
        this->pFSGS[m][15*5 +14] += 2.0 * zdsfs_29; // xzz zzzz (i=z)

        //this->pFSGS[m][15*6 + 0] += 0.0;          // i=x a=yyy
        //this->pFSGS[m][15*6 + 1] += 0.0;          // yyy xxxy (i=x)
        //this->pFSGS[m][15*6 + 2] += 0.0;          // yyy xxxz (i=x)
        //this->pFSGS[m][15*6 + 3] += 0.0;          // yyy xxyy (i=x)
        //this->pFSGS[m][15*6 + 4] += 0.0;          // yyy xxyz (i=x)
        //this->pFSGS[m][15*6 + 5] += 0.0;          // yyy xxzz (i=x)
        //this->pFSGS[m][15*6 + 6] += 0.0;          // yyy xyyy (i=x)
        //this->pFSGS[m][15*6 + 7] += 0.0;          // yyy xyyz (i=x)
        //this->pFSGS[m][15*6 + 8] += 0.0;          // yyy xyzz (i=x)
        //this->pFSGS[m][15*6 + 9] += 0.0;          // yyy xzzz (i=x)
        this->pFSGS[m][15*6 +10] += 3.0 * zdsfs_36; // yyy yyyy (i=y)
        this->pFSGS[m][15*6 +11] += 3.0 * zdsfs_37; // yyy yyyz (i=y)
        this->pFSGS[m][15*6 +12] += 3.0 * zdsfs_38; // yyy yyzz (i=y)
        this->pFSGS[m][15*6 +13] += 3.0 * zdsfs_39; // yyy yzzz (i=y)
        //this->pFSGS[m][15*6 +14] += 0.0;          // yyy zzzz (i=z)

        //this->pFSGS[m][15*7 + 0] += 0.0;          // i=x a=yyz
        //this->pFSGS[m][15*7 + 1] += 0.0;          // yyz xxxy (i=x)
        //this->pFSGS[m][15*7 + 2] += 0.0;          // yyz xxxz (i=x)
        //this->pFSGS[m][15*7 + 3] += 0.0;          // yyz xxyy (i=x)
        //this->pFSGS[m][15*7 + 4] += 0.0;          // yyz xxyz (i=x)
        //this->pFSGS[m][15*7 + 5] += 0.0;          // yyz xxzz (i=x)
        //this->pFSGS[m][15*7 + 6] += 0.0;          // yyz xyyy (i=x)
        //this->pFSGS[m][15*7 + 7] += 0.0;          // yyz xyyz (i=x)
        //this->pFSGS[m][15*7 + 8] += 0.0;          // yyz xyzz (i=x)
        //this->pFSGS[m][15*7 + 9] += 0.0;          // yyz xzzz (i=x)
        this->pFSGS[m][15*7 +10] += 2.0 * zdsfs_46; // yyz yyyy (i=y)
        this->pFSGS[m][15*7 +11] += 2.0 * zdsfs_47; // yyz yyyz (i=y)
        this->pFSGS[m][15*7 +12] += 2.0 * zdsfs_48; // yyz yyzz (i=y)
        this->pFSGS[m][15*7 +13] += 2.0 * zdsfs_49; // yyz yzzz (i=y)
        this->pFSGS[m][15*7 +14] +=       zdsfs_39; // yyz zzzz (i=z)

        //this->pFSGS[m][15*8 + 0] += 0.0;          // i=x a=yzz
        //this->pFSGS[m][15*8 + 1] += 0.0;          // yzz xxxy (i=x)
        //this->pFSGS[m][15*8 + 2] += 0.0;          // yzz xxxz (i=x)
        //this->pFSGS[m][15*8 + 3] += 0.0;          // yzz xxyy (i=x)
        //this->pFSGS[m][15*8 + 4] += 0.0;          // yzz xxyz (i=x)
        //this->pFSGS[m][15*8 + 5] += 0.0;          // yzz xxzz (i=x)
        //this->pFSGS[m][15*8 + 6] += 0.0;          // yzz xyyy (i=x)
        //this->pFSGS[m][15*8 + 7] += 0.0;          // yzz xyyz (i=x)
        //this->pFSGS[m][15*8 + 8] += 0.0;          // yzz xyzz (i=x)
        //this->pFSGS[m][15*8 + 9] += 0.0;          // yzz xzzz (i=x)
        this->pFSGS[m][15*8 +10] +=       zdsfs_56; // yzz yyyy (i=y)
        this->pFSGS[m][15*8 +11] +=       zdsfs_57; // yzz yyyz (i=y)
        this->pFSGS[m][15*8 +12] +=       zdsfs_58; // yzz yyzz (i=y)
        this->pFSGS[m][15*8 +13] +=       zdsfs_59; // yzz yzzz (i=y)
        this->pFSGS[m][15*8 +14] += 2.0 * zdsfs_49; // yzz zzzz (i=z)

        //this->pFSGS[m][15*9 + 0] += 0.0;          // i=x a=zzz
        //this->pFSGS[m][15*9 + 1] += 0.0;          // zzz xxxy (i=x)
        //this->pFSGS[m][15*9 + 2] += 0.0;          // zzz xxxz (i=x)
        //this->pFSGS[m][15*9 + 3] += 0.0;          // zzz xxyy (i=x)
        //this->pFSGS[m][15*9 + 4] += 0.0;          // zzz xxyz (i=x)
        //this->pFSGS[m][15*9 + 5] += 0.0;          // zzz xxzz (i=x)
        //this->pFSGS[m][15*9 + 6] += 0.0;          // zzz xyyy (i=x)
        //this->pFSGS[m][15*9 + 7] += 0.0;          // zzz xyyz (i=x)
        //this->pFSGS[m][15*9 + 8] += 0.0;          // zzz xyzz (i=x)
        //this->pFSGS[m][15*9 + 9] += 0.0;          // zzz xzzz (i=x)
        //this->pFSGS[m][15*9 +10] += 0.0;          // zzz yyyy (i=y)
        //this->pFSGS[m][15*9 +11] += 0.0;          // zzz yyyz (i=y)
        //this->pFSGS[m][15*9 +12] += 0.0;          // zzz yyzz (i=y)
        //this->pFSGS[m][15*9 +13] += 0.0;          // zzz yzzz (i=y)
        this->pFSGS[m][15*9 +14] += 3.0 * zdsfs_59; // zzz zzzz (i=z)
    }
}


void DfTEI::primitiveGSPS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    //const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    //const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 15; ++i) {
            this->pGSPS[m][3*i +0] = QCx * this->pGSSS[m][i] +WQx * this->pGSSS[m1][i];
            this->pGSPS[m][3*i +1] = QCy * this->pGSSS[m][i] +WQy * this->pGSSS[m1][i];
            this->pGSPS[m][3*i +2] = QCz * this->pGSSS[m][i] +WQz * this->pGSSS[m1][i];
        }

        const double zfsss_0 = ze2 * this->pFSSS[m1][0]; // xxx
        const double zfsss_1 = ze2 * this->pFSSS[m1][1]; // xxy
        const double zfsss_2 = ze2 * this->pFSSS[m1][2]; // xxz
        const double zfsss_3 = ze2 * this->pFSSS[m1][3]; // xyy
        const double zfsss_4 = ze2 * this->pFSSS[m1][4]; // xyz
        const double zfsss_5 = ze2 * this->pFSSS[m1][5]; // xzz
        const double zfsss_6 = ze2 * this->pFSSS[m1][6]; // yyy
        const double zfsss_7 = ze2 * this->pFSSS[m1][7]; // yyz
        const double zfsss_8 = ze2 * this->pFSSS[m1][8]; // yzz
        const double zfsss_9 = ze2 * this->pFSSS[m1][9]; // zzz

        this->pGSPS[m][3* 0 +0] += 4.0 * zfsss_0; // i=x, a=xxxx
        //this->pGSPS[m][3* 0 +1] += 0.0;         // i=y, a=xxxx
        //this->pGSPS[m][3* 0 +2] += 0.0;         // i=z, a=xxxx

        this->pGSPS[m][3* 1 +0] += 3.0 * zfsss_1; // i=x, a=xxxy
        this->pGSPS[m][3* 1 +1] +=       zfsss_0; // i=y, a=xxxy
        //this->pGSPS[m][3* 1 +2] += 0.0;         // i=z, a=xxxy

        this->pGSPS[m][3* 2 +0] += 3.0 * zfsss_2; // i=x, a=xxxz
        //this->pGSPS[m][3* 2 +1] += 0.0;         // i=y, a=xxxz
        this->pGSPS[m][3* 2 +2] +=       zfsss_0; // i=z, a=xxxz

        this->pGSPS[m][3* 3 +0] += 2.0 * zfsss_3; // i=x, a=xxyy
        this->pGSPS[m][3* 3 +1] += 2.0 * zfsss_1; // i=y, a=xxyy
        //this->pGSPS[m][3* 3 +2] += 0.0;         // i=z, a=xxyy

        this->pGSPS[m][3* 4 +0] += 2.0 * zfsss_4; // i=x, a=xxyz
        this->pGSPS[m][3* 4 +1] +=       zfsss_2; // i=y, a=xxyz
        this->pGSPS[m][3* 4 +2] +=       zfsss_1; // i=z, a=xxyz

        this->pGSPS[m][3* 5 +0] += 2.0 * zfsss_5; // i=x, a=xxzz
        //this->pGSPS[m][3* 5 +1] += 0.0;         // i=y, a=xxzz
        this->pGSPS[m][3* 5 +2] += 2.0 * zfsss_2; // i=z, a=xxzz

        this->pGSPS[m][3* 6 +0] +=       zfsss_6; // i=x, a=xyyy
        this->pGSPS[m][3* 6 +1] += 3.0 * zfsss_3; // i=y, a=xyyy
        //this->pGSPS[m][3* 6 +2] += 0.0;         // i=z, a=xyyy

        this->pGSPS[m][3* 7 +0] +=       zfsss_7; // i=x, a=xyyz
        this->pGSPS[m][3* 7 +1] += 2.0 * zfsss_4; // i=y, a=xyyz
        this->pGSPS[m][3* 7 +2] +=       zfsss_3; // i=z, a=xyyz

        this->pGSPS[m][3* 8 +0] +=       zfsss_8; // i=x, a=xyzz
        this->pGSPS[m][3* 8 +1] +=       zfsss_5; // i=y, a=xyzz
        this->pGSPS[m][3* 8 +2] += 2.0 * zfsss_4; // i=z, a=xyzz

        this->pGSPS[m][3* 9 +0] +=       zfsss_9; // i=x, a=xzzz
        //this->pGSPS[m][3* 9 +1] += 0.0;         // i=y, a=xzzz
        this->pGSPS[m][3* 9 +2] += 3.0 * zfsss_5; // i=z, a=xzzz

        //this->pGSPS[m][3*10 +0] += 0.0;         // i=x, a=yyyy
        this->pGSPS[m][3*10 +1] += 4.0 * zfsss_6; // i=y, a=yyyy
        //this->pGSPS[m][3*10 +2] += 0.0;         // i=z, a=yyyy

        //this->pGSPS[m][3*11 +0] += 0.0;         // i=x, a=yyyz
        this->pGSPS[m][3*11 +1] += 3.0 * zfsss_7; // i=y, a=yyyz
        this->pGSPS[m][3*11 +2] +=       zfsss_6; // i=z, a=yyyz

        //this->pGSPS[m][3*12 +0] += 0.0;         // i=x, a=yyzz
        this->pGSPS[m][3*12 +1] += 2.0 * zfsss_8; // i=y, a=yyzz
        this->pGSPS[m][3*12 +2] += 2.0 * zfsss_7; // i=z, a=yyzz

        //this->pGSPS[m][3*13 +0] += 0.0;         // i=x, a=yzzz
        this->pGSPS[m][3*13 +1] +=       zfsss_9; // i=y, a=yzzz
        this->pGSPS[m][3*13 +2] += 3.0 * zfsss_8; // i=z, a=yzzz

        //this->pGSPS[m][3*14 +0] += 0.0;         // i=x, a=zzzz
        //this->pGSPS[m][3*14 +1] += 0.0;         // i=y, a=zzzz
        this->pGSPS[m][3*14 +2] += 4.0 * zfsss_9; // i=z, a=zzzz
    }
}


void DfTEI::primitiveGSDS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 15; ++i) {
            this->pGSDS[m][6*i +0] = QCx * this->pGSPS[m][3*i +0] +WQx * this->pGSPS[m1][3*i +0];
            this->pGSDS[m][6*i +1] = QCx * this->pGSPS[m][3*i +1] +WQx * this->pGSPS[m1][3*i +1];
            this->pGSDS[m][6*i +2] = QCx * this->pGSPS[m][3*i +2] +WQx * this->pGSPS[m1][3*i +2];
            this->pGSDS[m][6*i +3] = QCy * this->pGSPS[m][3*i +1] +WQy * this->pGSPS[m1][3*i +1];
            this->pGSDS[m][6*i +4] = QCy * this->pGSPS[m][3*i +2] +WQy * this->pGSPS[m1][3*i +2];
            this->pGSDS[m][6*i +5] = QCz * this->pGSPS[m][3*i +2] +WQz * this->pGSPS[m1][3*i +2];
        }

        for (int i = 0; i < 15; ++i) {
            const double zgsss = eta2 * (this->pGSSS[m][i] -re * this->pGSSS[m1][i]);
            this->pGSDS[m][6*i +0] += zgsss; //
            this->pGSDS[m][6*i +3] += zgsss; //
            this->pGSDS[m][6*i +5] += zgsss; //
        }

        const double zfsps_00 = ze2 * this->pFSPS[m1][3*0 +0]; // xxx x
        const double zfsps_01 = ze2 * this->pFSPS[m1][3*0 +1]; // xxx y
        const double zfsps_02 = ze2 * this->pFSPS[m1][3*0 +2]; // xxx z
        const double zfsps_10 = ze2 * this->pFSPS[m1][3*1 +0]; // xxy x
        const double zfsps_11 = ze2 * this->pFSPS[m1][3*1 +1]; // xxy y
        const double zfsps_12 = ze2 * this->pFSPS[m1][3*1 +2]; // xxy z
        const double zfsps_20 = ze2 * this->pFSPS[m1][3*2 +0]; // xxz x
        const double zfsps_21 = ze2 * this->pFSPS[m1][3*2 +1]; // xxz y
        const double zfsps_22 = ze2 * this->pFSPS[m1][3*2 +2]; // xxz z
        const double zfsps_30 = ze2 * this->pFSPS[m1][3*3 +0]; // xyy x
        const double zfsps_31 = ze2 * this->pFSPS[m1][3*3 +1]; // xyy y
        const double zfsps_32 = ze2 * this->pFSPS[m1][3*3 +2]; // xyy z
        const double zfsps_40 = ze2 * this->pFSPS[m1][3*4 +0]; // xyz x
        const double zfsps_41 = ze2 * this->pFSPS[m1][3*4 +1]; // xyz y
        const double zfsps_42 = ze2 * this->pFSPS[m1][3*4 +2]; // xyz z
        const double zfsps_50 = ze2 * this->pFSPS[m1][3*5 +0]; // xzz x
        const double zfsps_51 = ze2 * this->pFSPS[m1][3*5 +1]; // xzz y
        const double zfsps_52 = ze2 * this->pFSPS[m1][3*5 +2]; // xzz z
        const double zfsps_60 = ze2 * this->pFSPS[m1][3*6 +0]; // yyy x
        const double zfsps_61 = ze2 * this->pFSPS[m1][3*6 +1]; // yyy y
        const double zfsps_62 = ze2 * this->pFSPS[m1][3*6 +2]; // yyy z
        const double zfsps_70 = ze2 * this->pFSPS[m1][3*7 +0]; // yyz x
        const double zfsps_71 = ze2 * this->pFSPS[m1][3*7 +1]; // yyz y
        const double zfsps_72 = ze2 * this->pFSPS[m1][3*7 +2]; // yyz z
        const double zfsps_80 = ze2 * this->pFSPS[m1][3*8 +0]; // yzz x
        const double zfsps_81 = ze2 * this->pFSPS[m1][3*8 +1]; // yzz y
        const double zfsps_82 = ze2 * this->pFSPS[m1][3*8 +2]; // yzz z
        const double zfsps_90 = ze2 * this->pFSPS[m1][3*9 +0]; // zzz x
        const double zfsps_91 = ze2 * this->pFSPS[m1][3*9 +1]; // zzz y
        const double zfsps_92 = ze2 * this->pFSPS[m1][3*9 +2]; // zzz z

        this->pGSDS[m][6* 0 +0] += 4.0 * zfsps_00; // i=x, a=xxxx
        this->pGSDS[m][6* 0 +1] += 4.0 * zfsps_01; // i=x, a=xxxx
        this->pGSDS[m][6* 0 +2] += 4.0 * zfsps_02; // i=x, a=xxxx
        //this->pGSDS[m][6* 0 +3] += 0.0;          // i=y, a=xxxx
        //this->pGSDS[m][6* 0 +4] += 0.0;          // i=y, a=xxxx
        //this->pGSDS[m][6* 0 +5] += 0.0;          // i=z, a=xxxx

        this->pGSDS[m][6* 1 +0] += 3.0 * zfsps_10; // i=x, a=xxxy
        this->pGSDS[m][6* 1 +1] += 3.0 * zfsps_11; // i=x, a=xxxy
        this->pGSDS[m][6* 1 +2] += 3.0 * zfsps_12; // i=x, a=xxxy
        this->pGSDS[m][6* 1 +3] +=       zfsps_01; // i=y, a=xxxy
        this->pGSDS[m][6* 1 +4] +=       zfsps_02; // i=y, a=xxxy
        //this->pGSDS[m][6* 1 +5] += 0.0;          // i=z, a=xxxy

        this->pGSDS[m][6* 2 +0] += 3.0 * zfsps_20; // i=x, a=xxxz
        this->pGSDS[m][6* 2 +1] += 3.0 * zfsps_21; // i=x, a=xxxz
        this->pGSDS[m][6* 2 +2] += 3.0 * zfsps_22; // i=x, a=xxxz
        //this->pGSDS[m][6* 2 +3] += 0.0;          // i=y, a=xxxz
        //this->pGSDS[m][6* 2 +4] += 0.0;          // i=y, a=xxxz
        this->pGSDS[m][6* 2 +5] +=       zfsps_02; // i=z, a=xxxz

        this->pGSDS[m][6* 3 +0] += 2.0 * zfsps_30; // i=x, a=xxyy
        this->pGSDS[m][6* 3 +1] += 2.0 * zfsps_31; // i=x, a=xxyy
        this->pGSDS[m][6* 3 +2] += 2.0 * zfsps_32; // i=x, a=xxyy
        this->pGSDS[m][6* 3 +3] += 2.0 * zfsps_11; // i=y, a=xxyy
        this->pGSDS[m][6* 3 +4] += 2.0 * zfsps_12; // i=y, a=xxyy
        //this->pGSDS[m][6* 3 +5] += 0.0;          // i=z, a=xxyy

        this->pGSDS[m][6* 4 +0] += 2.0 * zfsps_40; // i=x, a=xxyz
        this->pGSDS[m][6* 4 +1] += 2.0 * zfsps_41; // i=x, a=xxyz
        this->pGSDS[m][6* 4 +2] += 2.0 * zfsps_42; // i=x, a=xxyz
        this->pGSDS[m][6* 4 +3] +=       zfsps_21; // i=y, a=xxyz
        this->pGSDS[m][6* 4 +4] +=       zfsps_22; // i=y, a=xxyz
        this->pGSDS[m][6* 4 +5] +=       zfsps_12; // i=z, a=xxyz

        this->pGSDS[m][6* 5 +0] += 2.0 * zfsps_50; // i=x, a=xxzz
        this->pGSDS[m][6* 5 +1] += 2.0 * zfsps_51; // i=x, a=xxzz
        this->pGSDS[m][6* 5 +2] += 2.0 * zfsps_52; // i=x, a=xxzz
        //this->pGSDS[m][6* 5 +3] += 0.0;          // i=y, a=xxzz
        //this->pGSDS[m][6* 5 +4] += 0.0;          // i=y, a=xxzz
        this->pGSDS[m][6* 5 +5] += 2.0 * zfsps_22; // i=z, a=xxzz

        this->pGSDS[m][6* 6 +0] +=       zfsps_60; // i=x, a=xyyy
        this->pGSDS[m][6* 6 +1] +=       zfsps_61; // i=x, a=xyyy
        this->pGSDS[m][6* 6 +2] +=       zfsps_62; // i=x, a=xyyy
        this->pGSDS[m][6* 6 +3] += 3.0 * zfsps_31; // i=y, a=xyyy
        this->pGSDS[m][6* 6 +4] += 3.0 * zfsps_32; // i=y, a=xyyy
        //this->pGSDS[m][6* 6 +5] += 0.0;          // i=z, a=xyyy

        this->pGSDS[m][6* 7 +0] +=       zfsps_70; // i=x, a=xyyz
        this->pGSDS[m][6* 7 +1] +=       zfsps_71; // i=x, a=xyyz
        this->pGSDS[m][6* 7 +2] +=       zfsps_72; // i=x, a=xyyz
        this->pGSDS[m][6* 7 +3] += 2.0 * zfsps_41; // i=y, a=xyyz
        this->pGSDS[m][6* 7 +4] += 2.0 * zfsps_42; // i=y, a=xyyz
        this->pGSDS[m][6* 7 +5] +=       zfsps_32; // i=z, a=xyyz

        this->pGSDS[m][6* 8 +0] +=       zfsps_80; // i=x, a=xyzz
        this->pGSDS[m][6* 8 +1] +=       zfsps_81; // i=x, a=xyzz
        this->pGSDS[m][6* 8 +2] +=       zfsps_82; // i=x, a=xyzz
        this->pGSDS[m][6* 8 +3] +=       zfsps_51; // i=y, a=xyzz
        this->pGSDS[m][6* 8 +4] +=       zfsps_52; // i=y, a=xyzz
        this->pGSDS[m][6* 8 +5] += 2.0 * zfsps_42; // i=z, a=xyzz

        this->pGSDS[m][6* 9 +0] +=       zfsps_90; // i=x, a=xzzz
        this->pGSDS[m][6* 9 +1] +=       zfsps_91; // i=x, a=xzzz
        this->pGSDS[m][6* 9 +2] +=       zfsps_92; // i=x, a=xzzz
        //this->pGSDS[m][6* 9 +3] += 0.0;          // i=y, a=xzzz
        //this->pGSDS[m][6* 9 +4] += 0.0;          // i=y, a=xzzz
        this->pGSDS[m][6* 9 +5] += 3.0 * zfsps_52; // i=z, a=xzzz

        //this->pGSDS[m][6*10 +0] += 0.0;          // i=x, a=yyyy
        //this->pGSDS[m][6*10 +1] += 0.0;          // i=x, a=yyyy
        //this->pGSDS[m][6*10 +2] += 0.0;          // i=x, a=yyyy
        this->pGSDS[m][6*10 +3] += 4.0 * zfsps_61; // i=y, a=yyyy
        this->pGSDS[m][6*10 +4] += 4.0 * zfsps_62; // i=y, a=yyyy
        //this->pGSDS[m][6*10 +5] += 0.0;          // i=z, a=yyyy

        //this->pGSDS[m][6*11 +0] += 0.0;          // i=x, a=yyyz
        //this->pGSDS[m][6*11 +1] += 0.0;          // i=x, a=yyyz
        //this->pGSDS[m][6*11 +2] += 0.0;          // i=x, a=yyyz
        this->pGSDS[m][6*11 +3] += 3.0 * zfsps_71; // i=y, a=yyyz
        this->pGSDS[m][6*11 +4] += 3.0 * zfsps_72; // i=y, a=yyyz
        this->pGSDS[m][6*11 +5] +=       zfsps_62; // i=z, a=yyyz

        //this->pGSDS[m][6*12 +0] += 0.0;          // i=x, a=yyzz
        //this->pGSDS[m][6*12 +1] += 0.0;          // i=x, a=yyzz
        //this->pGSDS[m][6*12 +2] += 0.0;          // i=x, a=yyzz
        this->pGSDS[m][6*12 +3] += 2.0 * zfsps_81; // i=y, a=yyzz
        this->pGSDS[m][6*12 +4] += 2.0 * zfsps_82; // i=y, a=yyzz
        this->pGSDS[m][6*12 +5] += 2.0 * zfsps_72; // i=z, a=yyzz

        //this->pGSDS[m][6*13 +0] += 0.0;          // i=x, a=yzzz
        //this->pGSDS[m][6*13 +1] += 0.0;          // i=x, a=yzzz
        //this->pGSDS[m][6*13 +2] += 0.0;          // i=x, a=yzzz
        this->pGSDS[m][6*13 +3] +=       zfsps_91; // i=y, a=yzzz
        this->pGSDS[m][6*13 +4] +=       zfsps_92; // i=y, a=yzzz
        this->pGSDS[m][6*13 +5] += 3.0 * zfsps_82; // i=z, a=yzzz

        //this->pGSDS[m][6*14 +0] += 0.0;          // i=x, a=zzzz
        //this->pGSDS[m][6*14 +1] += 0.0;          // i=x, a=zzzz
        //this->pGSDS[m][6*14 +2] += 0.0;          // i=x, a=zzzz
        //this->pGSDS[m][6*14 +3] += 0.0;          // i=y, a=zzzz
        //this->pGSDS[m][6*14 +4] += 0.0;          // i=y, a=zzzz
        this->pGSDS[m][6*14 +5] += 4.0 * zfsps_92; // i=z, a=zzzz
    }
}


void DfTEI::primitiveGSFS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 15; ++i) {
            this->pGSFS[m][10*i +0] = QCx * this->pGSDS[m][6*i +0] +WQx * this->pGSDS[m1][6*i +0]; // xxx
            this->pGSFS[m][10*i +1] = QCx * this->pGSDS[m][6*i +1] +WQx * this->pGSDS[m1][6*i +1]; // xxy
            this->pGSFS[m][10*i +2] = QCx * this->pGSDS[m][6*i +2] +WQx * this->pGSDS[m1][6*i +2]; // xxz
            this->pGSFS[m][10*i +3] = QCx * this->pGSDS[m][6*i +3] +WQx * this->pGSDS[m1][6*i +3]; // xyy
            this->pGSFS[m][10*i +4] = QCx * this->pGSDS[m][6*i +4] +WQx * this->pGSDS[m1][6*i +4]; // xyz
            this->pGSFS[m][10*i +5] = QCx * this->pGSDS[m][6*i +5] +WQx * this->pGSDS[m1][6*i +5]; // xzz
            this->pGSFS[m][10*i +6] = QCy * this->pGSDS[m][6*i +3] +WQy * this->pGSDS[m1][6*i +3]; // yyy
            this->pGSFS[m][10*i +7] = QCy * this->pGSDS[m][6*i +4] +WQy * this->pGSDS[m1][6*i +4]; // yyz
            this->pGSFS[m][10*i +8] = QCy * this->pGSDS[m][6*i +5] +WQy * this->pGSDS[m1][6*i +5]; // yzz
            this->pGSFS[m][10*i +9] = QCz * this->pGSDS[m][6*i +5] +WQz * this->pGSDS[m1][6*i +5]; // zzz
        }

        for (int i = 0; i < 15; ++i) {
            const double zgsps_0 = eta2 * (this->pGSPS[m][3*i +0] -re * this->pGSPS[m1][3*i +0]);
            const double zgsps_1 = eta2 * (this->pGSPS[m][3*i +1] -re * this->pGSPS[m1][3*i +1]);
            const double zgsps_2 = eta2 * (this->pGSPS[m][3*i +2] -re * this->pGSPS[m1][3*i +2]);
            this->pGSFS[m][10*i +0] += 2.0 * zgsps_0; // i=x, c=xx
            this->pGSFS[m][10*i +1] +=       zgsps_1; // i=x, c=xy
            this->pGSFS[m][10*i +2] +=       zgsps_2; // i=x, c=xz
            //this->pGSFS[m][10*i +3] += zgsps; // i=x, c=yy
            //this->pGSFS[m][10*i +4] += zgsps; // i=x, c=yz
            //this->pGSFS[m][10*i +5] += zgsps; // i=x, c=zz
            this->pGSFS[m][10*i +6] += 2.0 * zgsps_1; // i=y, c=yy
            this->pGSFS[m][10*i +7] +=       zgsps_2; // i=y, c=yz
            //this->pGSFS[m][10*i +8] += zgsps; // i=y, c=zz
            this->pGSFS[m][10*i +9] += 2.0 * zgsps_2; // i=z, c=zz
        }

        //--------------------
        const double zfsds_00 = ze2 * this->pFSDS[m1][6*0 +0]; // xxx xx
        const double zfsds_01 = ze2 * this->pFSDS[m1][6*0 +1]; // xxx xy
        const double zfsds_02 = ze2 * this->pFSDS[m1][6*0 +2]; // xxx xz
        const double zfsds_03 = ze2 * this->pFSDS[m1][6*0 +3]; // xxx yy
        const double zfsds_04 = ze2 * this->pFSDS[m1][6*0 +4]; // xxx yz
        const double zfsds_05 = ze2 * this->pFSDS[m1][6*0 +5]; // xxx zz
        const double zfsds_10 = ze2 * this->pFSDS[m1][6*1 +0]; // xxy xx
        const double zfsds_11 = ze2 * this->pFSDS[m1][6*1 +1]; // xxy xy
        const double zfsds_12 = ze2 * this->pFSDS[m1][6*1 +2]; // xxy xz
        const double zfsds_13 = ze2 * this->pFSDS[m1][6*1 +3]; // xxy yy
        const double zfsds_14 = ze2 * this->pFSDS[m1][6*1 +4]; // xxy yz
        const double zfsds_15 = ze2 * this->pFSDS[m1][6*1 +5]; // xxy zz
        const double zfsds_20 = ze2 * this->pFSDS[m1][6*2 +0]; // xxz xx
        const double zfsds_21 = ze2 * this->pFSDS[m1][6*2 +1]; // xxz xy
        const double zfsds_22 = ze2 * this->pFSDS[m1][6*2 +2]; // xxz xz
        const double zfsds_23 = ze2 * this->pFSDS[m1][6*2 +3]; // xxz yy
        const double zfsds_24 = ze2 * this->pFSDS[m1][6*2 +4]; // xxz yz
        const double zfsds_25 = ze2 * this->pFSDS[m1][6*2 +5]; // xxz zz
        const double zfsds_30 = ze2 * this->pFSDS[m1][6*3 +0]; // xyy xx
        const double zfsds_31 = ze2 * this->pFSDS[m1][6*3 +1]; // xyy xy
        const double zfsds_32 = ze2 * this->pFSDS[m1][6*3 +2]; // xyy xz
        const double zfsds_33 = ze2 * this->pFSDS[m1][6*3 +3]; // xyy yy
        const double zfsds_34 = ze2 * this->pFSDS[m1][6*3 +4]; // xyy yz
        const double zfsds_35 = ze2 * this->pFSDS[m1][6*3 +5]; // xyy zz
        const double zfsds_40 = ze2 * this->pFSDS[m1][6*4 +0]; // xyz xx
        const double zfsds_41 = ze2 * this->pFSDS[m1][6*4 +1]; // xyz xy
        const double zfsds_42 = ze2 * this->pFSDS[m1][6*4 +2]; // xyz xz
        const double zfsds_43 = ze2 * this->pFSDS[m1][6*4 +3]; // xyz yy
        const double zfsds_44 = ze2 * this->pFSDS[m1][6*4 +4]; // xyz yz
        const double zfsds_45 = ze2 * this->pFSDS[m1][6*4 +5]; // xyz zz
        const double zfsds_50 = ze2 * this->pFSDS[m1][6*5 +0]; // xzz xx
        const double zfsds_51 = ze2 * this->pFSDS[m1][6*5 +1]; // xzz xy
        const double zfsds_52 = ze2 * this->pFSDS[m1][6*5 +2]; // xzz xz
        const double zfsds_53 = ze2 * this->pFSDS[m1][6*5 +3]; // xzz yy
        const double zfsds_54 = ze2 * this->pFSDS[m1][6*5 +4]; // xzz yz
        const double zfsds_55 = ze2 * this->pFSDS[m1][6*5 +5]; // xzz zz
        const double zfsds_60 = ze2 * this->pFSDS[m1][6*6 +0]; // yyy xx
        const double zfsds_61 = ze2 * this->pFSDS[m1][6*6 +1]; // yyy xy
        const double zfsds_62 = ze2 * this->pFSDS[m1][6*6 +2]; // yyy xz
        const double zfsds_63 = ze2 * this->pFSDS[m1][6*6 +3]; // yyy yy
        const double zfsds_64 = ze2 * this->pFSDS[m1][6*6 +4]; // yyy yz
        const double zfsds_65 = ze2 * this->pFSDS[m1][6*6 +5]; // yyy zz
        const double zfsds_70 = ze2 * this->pFSDS[m1][6*7 +0]; // yyz xx
        const double zfsds_71 = ze2 * this->pFSDS[m1][6*7 +1]; // yyz xy
        const double zfsds_72 = ze2 * this->pFSDS[m1][6*7 +2]; // yyz xz
        const double zfsds_73 = ze2 * this->pFSDS[m1][6*7 +3]; // yyz yy
        const double zfsds_74 = ze2 * this->pFSDS[m1][6*7 +4]; // yyz yz
        const double zfsds_75 = ze2 * this->pFSDS[m1][6*7 +5]; // yyz zz
        const double zfsds_80 = ze2 * this->pFSDS[m1][6*8 +0]; // yzz xx
        const double zfsds_81 = ze2 * this->pFSDS[m1][6*8 +1]; // yzz xy
        const double zfsds_82 = ze2 * this->pFSDS[m1][6*8 +2]; // yzz xz
        const double zfsds_83 = ze2 * this->pFSDS[m1][6*8 +3]; // yzz yy
        const double zfsds_84 = ze2 * this->pFSDS[m1][6*8 +4]; // yzz yz
        const double zfsds_85 = ze2 * this->pFSDS[m1][6*8 +5]; // yzz zz
        const double zfsds_90 = ze2 * this->pFSDS[m1][6*9 +0]; // zzz xx
        const double zfsds_91 = ze2 * this->pFSDS[m1][6*9 +1]; // zzz xy
        const double zfsds_92 = ze2 * this->pFSDS[m1][6*9 +2]; // zzz xz
        const double zfsds_93 = ze2 * this->pFSDS[m1][6*9 +3]; // zzz yy
        const double zfsds_94 = ze2 * this->pFSDS[m1][6*9 +4]; // zzz yz
        const double zfsds_95 = ze2 * this->pFSDS[m1][6*9 +5]; // zzz zz

        this->pGSFS[m][10* 0 +0] += 4.0 * zfsds_00; // i=x, a=xxxx
        this->pGSFS[m][10* 0 +1] += 4.0 * zfsds_01; // i=x, a=xxxx
        this->pGSFS[m][10* 0 +2] += 4.0 * zfsds_02; // i=x, a=xxxx
        this->pGSFS[m][10* 0 +3] += 4.0 * zfsds_03; // i=x, a=xxxx
        this->pGSFS[m][10* 0 +4] += 4.0 * zfsds_04; // i=x, a=xxxx
        this->pGSFS[m][10* 0 +5] += 4.0 * zfsds_05; // i=x, a=xxxx
        //this->pGSFS[m][10* 0 +6] += 0.0;          // i=y, a=xxxx
        //this->pGSFS[m][10* 0 +7] += 0.0;          // i=y, a=xxxx
        //this->pGSFS[m][10* 0 +8] += 0.0;          // i=y, a=xxxx
        //this->pGSFS[m][10* 0 +9] += 0.0;          // i=z, a=xxxx

        this->pGSFS[m][10* 1 +0] += 3.0 * zfsds_10; // i=x, a=xxxy
        this->pGSFS[m][10* 1 +1] += 3.0 * zfsds_11; // i=x, a=xxxy
        this->pGSFS[m][10* 1 +2] += 3.0 * zfsds_12; // i=x, a=xxxy
        this->pGSFS[m][10* 1 +3] += 3.0 * zfsds_13; // i=x, a=xxxy
        this->pGSFS[m][10* 1 +4] += 3.0 * zfsds_14; // i=x, a=xxxy
        this->pGSFS[m][10* 1 +5] += 3.0 * zfsds_15; // i=x, a=xxxy
        this->pGSFS[m][10* 1 +6] +=       zfsds_03; // i=y, a=xxxy
        this->pGSFS[m][10* 1 +7] +=       zfsds_04; // i=y, a=xxxy
        this->pGSFS[m][10* 1 +8] +=       zfsds_05; // i=y, a=xxxy
        //this->pGSFS[m][10* 1 +9] += 0.0;          // i=z, a=xxxy

        this->pGSFS[m][10* 2 +0] += 3.0 * zfsds_20; // i=x, a=xxxz
        this->pGSFS[m][10* 2 +1] += 3.0 * zfsds_21; // i=x, a=xxxz
        this->pGSFS[m][10* 2 +2] += 3.0 * zfsds_22; // i=x, a=xxxz
        this->pGSFS[m][10* 2 +3] += 3.0 * zfsds_23; // i=x, a=xxxz
        this->pGSFS[m][10* 2 +4] += 3.0 * zfsds_24; // i=x, a=xxxz
        this->pGSFS[m][10* 2 +5] += 3.0 * zfsds_25; // i=x, a=xxxz
        //this->pGSFS[m][10* 2 +6] += 0.0;          // i=y, a=xxxz
        //this->pGSFS[m][10* 2 +7] += 0.0;          // i=y, a=xxxz
        //this->pGSFS[m][10* 2 +8] += 0.0;          // i=y, a=xxxz
        this->pGSFS[m][10* 2 +9] +=       zfsds_05; // i=z, a=xxxz

        this->pGSFS[m][10* 3 +0] += 2.0 * zfsds_30; // i=x, a=xxyy
        this->pGSFS[m][10* 3 +1] += 2.0 * zfsds_31; // i=x, a=xxyy
        this->pGSFS[m][10* 3 +2] += 2.0 * zfsds_32; // i=x, a=xxyy
        this->pGSFS[m][10* 3 +3] += 2.0 * zfsds_33; // i=x, a=xxyy
        this->pGSFS[m][10* 3 +4] += 2.0 * zfsds_34; // i=x, a=xxyy
        this->pGSFS[m][10* 3 +5] += 2.0 * zfsds_35; // i=x, a=xxyy
        this->pGSFS[m][10* 3 +6] += 2.0 * zfsds_13; // i=y, a=xxyy
        this->pGSFS[m][10* 3 +7] += 2.0 * zfsds_14; // i=y, a=xxyy
        this->pGSFS[m][10* 3 +8] += 2.0 * zfsds_15; // i=y, a=xxyy
        //this->pGSFS[m][10* 3 +9] += 0.0;          // i=z, a=xxyy

        this->pGSFS[m][10* 4 +0] += 2.0 * zfsds_40; // i=x, a=xxyz
        this->pGSFS[m][10* 4 +1] += 2.0 * zfsds_41; // i=x, a=xxyz
        this->pGSFS[m][10* 4 +2] += 2.0 * zfsds_42; // i=x, a=xxyz
        this->pGSFS[m][10* 4 +3] += 2.0 * zfsds_43; // i=x, a=xxyz
        this->pGSFS[m][10* 4 +4] += 2.0 * zfsds_44; // i=x, a=xxyz
        this->pGSFS[m][10* 4 +5] += 2.0 * zfsds_45; // i=x, a=xxyz
        this->pGSFS[m][10* 4 +6] +=       zfsds_23; // i=y, a=xxyz
        this->pGSFS[m][10* 4 +7] +=       zfsds_24; // i=y, a=xxyz
        this->pGSFS[m][10* 4 +8] +=       zfsds_25; // i=y, a=xxyz
        this->pGSFS[m][10* 4 +9] +=       zfsds_15; // i=z, a=xxyz

        this->pGSFS[m][10* 5 +0] += 2.0 * zfsds_50; // i=x, a=xxzz
        this->pGSFS[m][10* 5 +1] += 2.0 * zfsds_51; // i=x, a=xxzz
        this->pGSFS[m][10* 5 +2] += 2.0 * zfsds_52; // i=x, a=xxzz
        this->pGSFS[m][10* 5 +3] += 2.0 * zfsds_53; // i=x, a=xxzz
        this->pGSFS[m][10* 5 +4] += 2.0 * zfsds_54; // i=x, a=xxzz
        this->pGSFS[m][10* 5 +5] += 2.0 * zfsds_55; // i=x, a=xxzz
        //this->pGSFS[m][10* 5 +6] += 0.0;          // i=y, a=xxzz
        //this->pGSFS[m][10* 5 +7] += 0.0;          // i=y, a=xxzz
        //this->pGSFS[m][10* 5 +8] += 0.0;          // i=y, a=xxzz
        this->pGSFS[m][10* 5 +9] += 2.0 * zfsds_25; // i=z, a=xxzz

        this->pGSFS[m][10* 6 +0] +=       zfsds_60; // i=x, a=xyyy
        this->pGSFS[m][10* 6 +1] +=       zfsds_61; // i=x, a=xyyy
        this->pGSFS[m][10* 6 +2] +=       zfsds_62; // i=x, a=xyyy
        this->pGSFS[m][10* 6 +3] +=       zfsds_63; // i=x, a=xyyy
        this->pGSFS[m][10* 6 +4] +=       zfsds_64; // i=x, a=xyyy
        this->pGSFS[m][10* 6 +5] +=       zfsds_65; // i=x, a=xyyy
        this->pGSFS[m][10* 6 +6] += 3.0 * zfsds_33; // i=y, a=xyyy
        this->pGSFS[m][10* 6 +7] += 3.0 * zfsds_34; // i=y, a=xyyy
        this->pGSFS[m][10* 6 +8] += 3.0 * zfsds_35; // i=y, a=xyyy
        //this->pGSFS[m][10* 6 +9] += 0.0           // i=z, a=xyyy

        this->pGSFS[m][10* 7 +0] +=       zfsds_70; // i=x, a=xyyz
        this->pGSFS[m][10* 7 +1] +=       zfsds_71; // i=x, a=xyyz
        this->pGSFS[m][10* 7 +2] +=       zfsds_72; // i=x, a=xyyz
        this->pGSFS[m][10* 7 +3] +=       zfsds_73; // i=x, a=xyyz
        this->pGSFS[m][10* 7 +4] +=       zfsds_74; // i=x, a=xyyz
        this->pGSFS[m][10* 7 +5] +=       zfsds_75; // i=x, a=xyyz
        this->pGSFS[m][10* 7 +6] += 2.0 * zfsds_43; // i=y, a=xyyz
        this->pGSFS[m][10* 7 +7] += 2.0 * zfsds_44; // i=y, a=xyyz
        this->pGSFS[m][10* 7 +8] += 2.0 * zfsds_45; // i=y, a=xyyz
        this->pGSFS[m][10* 7 +9] +=       zfsds_35; // i=z, a=xyyz

        this->pGSFS[m][10* 8 +0] +=       zfsds_80; // i=x, a=xyzz
        this->pGSFS[m][10* 8 +1] +=       zfsds_81; // i=x, a=xyzz
        this->pGSFS[m][10* 8 +2] +=       zfsds_82; // i=x, a=xyzz
        this->pGSFS[m][10* 8 +3] +=       zfsds_83; // i=x, a=xyzz
        this->pGSFS[m][10* 8 +4] +=       zfsds_84; // i=x, a=xyzz
        this->pGSFS[m][10* 8 +5] +=       zfsds_85; // i=x, a=xyzz
        this->pGSFS[m][10* 8 +6] +=       zfsds_53; // i=y, a=xyzz
        this->pGSFS[m][10* 8 +7] +=       zfsds_54; // i=y, a=xyzz
        this->pGSFS[m][10* 8 +8] +=       zfsds_55; // i=y, a=xyzz
        this->pGSFS[m][10* 8 +9] += 2.0 * zfsds_45; // i=z, a=xyzz

        this->pGSFS[m][10* 9 +0] +=       zfsds_90; // i=x, a=xzzz
        this->pGSFS[m][10* 9 +1] +=       zfsds_91; // i=x, a=xzzz
        this->pGSFS[m][10* 9 +2] +=       zfsds_92; // i=x, a=xzzz
        this->pGSFS[m][10* 9 +3] +=       zfsds_93; // i=x, a=xzzz
        this->pGSFS[m][10* 9 +4] +=       zfsds_94; // i=x, a=xzzz
        this->pGSFS[m][10* 9 +5] +=       zfsds_95; // i=x, a=xzzz
        //this->pGSFS[m][10* 9 +6] += 0.0;          // i=y, a=xzzz
        //this->pGSFS[m][10* 9 +7] += 0.0;          // i=y, a=xzzz
        //this->pGSFS[m][10* 9 +8] += 0.0;          // i=y, a=xzzz
        this->pGSFS[m][10* 9 +9] += 3.0 * zfsds_55; // i=z, a=xzzz

        //this->pGSFS[m][10*10 +0] += 0.0;          // i=x, a=yyyy
        //this->pGSFS[m][10*10 +1] += 0.0;          // i=x, a=yyyy
        //this->pGSFS[m][10*10 +2] += 0.0;          // i=x, a=yyyy
        //this->pGSFS[m][10*10 +3] += 0.0;          // i=x, a=yyyy
        //this->pGSFS[m][10*10 +4] += 0.0;          // i=x, a=yyyy
        //this->pGSFS[m][10*10 +5] += 0.0;          // i=x, a=yyyy
        this->pGSFS[m][10*10 +6] += 4.0 * zfsds_63; // i=y, a=yyyy
        this->pGSFS[m][10*10 +7] += 4.0 * zfsds_64; // i=y, a=yyyy
        this->pGSFS[m][10*10 +8] += 4.0 * zfsds_65; // i=y, a=yyyy
        //this->pGSFS[m][10*10 +9] += 0.0;          // i=z, a=yyyy

        //this->pGSFS[m][10*11 +0] += 0.0;          // i=x, a=yyyz
        //this->pGSFS[m][10*11 +1] += 0.0;          // i=x, a=yyyz
        //this->pGSFS[m][10*11 +2] += 0.0;          // i=x, a=yyyz
        //this->pGSFS[m][10*11 +3] += 0.0;          // i=x, a=yyyz
        //this->pGSFS[m][10*11 +4] += 0.0;          // i=x, a=yyyz
        //this->pGSFS[m][10*11 +5] += 0.0;          // i=x, a=yyyz
        this->pGSFS[m][10*11 +6] += 3.0 * zfsds_73; // i=y, a=yyyz
        this->pGSFS[m][10*11 +7] += 3.0 * zfsds_74; // i=y, a=yyyz
        this->pGSFS[m][10*11 +8] += 3.0 * zfsds_75; // i=y, a=yyyz
        this->pGSFS[m][10*11 +9] +=       zfsds_65; // i=z, a=yyyz

        //this->pGSFS[m][10*12 +0] += 0.0;          // i=x, a=yyzz
        //this->pGSFS[m][10*12 +1] += 0.0;          // i=x, a=yyzz
        //this->pGSFS[m][10*12 +2] += 0.0;          // i=x, a=yyzz
        //this->pGSFS[m][10*12 +3] += 0.0;          // i=x, a=yyzz
        //this->pGSFS[m][10*12 +4] += 0.0;          // i=x, a=yyzz
        //this->pGSFS[m][10*12 +5] += 0.0;          // i=x, a=yyzz
        this->pGSFS[m][10*12 +6] += 2.0 * zfsds_83; // i=y, a=yyzz
        this->pGSFS[m][10*12 +7] += 2.0 * zfsds_84; // i=y, a=yyzz
        this->pGSFS[m][10*12 +8] += 2.0 * zfsds_85; // i=y, a=yyzz
        this->pGSFS[m][10*12 +9] += 2.0 * zfsds_75; // i=z, a=yyzz

        //this->pGSFS[m][10*13 +0] += 0.0;          // i=x, a=yzzz
        //this->pGSFS[m][10*13 +1] += 0.0;          // i=x, a=yzzz
        //this->pGSFS[m][10*13 +2] += 0.0;          // i=x, a=yzzz
        //this->pGSFS[m][10*13 +3] += 0.0;          // i=x, a=yzzz
        //this->pGSFS[m][10*13 +4] += 0.0;          // i=x, a=yzzz
        //this->pGSFS[m][10*13 +5] += 0.0;          // i=x, a=yzzz
        this->pGSFS[m][10*13 +6] +=       zfsds_93; // i=y, a=yzzz
        this->pGSFS[m][10*13 +7] +=       zfsds_94; // i=y, a=yzzz
        this->pGSFS[m][10*13 +8] +=       zfsds_95; // i=y, a=yzzz
        this->pGSFS[m][10*13 +9] += 3.0 * zfsds_85; // i=z, a=yzzz

        //this->pGSFS[m][10*14 +0] += 0.0;          // i=x, a=zzzz
        //this->pGSFS[m][10*14 +1] += 0.0;          // i=x, a=zzzz
        //this->pGSFS[m][10*14 +2] += 0.0;          // i=x, a=zzzz
        //this->pGSFS[m][10*14 +3] += 0.0;          // i=x, a=zzzz
        //this->pGSFS[m][10*14 +4] += 0.0;          // i=x, a=zzzz
        //this->pGSFS[m][10*14 +5] += 0.0;          // i=x, a=zzzz
        //this->pGSFS[m][10*14 +6] += 0.0;          // i=y, a=zzzz
        //this->pGSFS[m][10*14 +7] += 0.0;          // i=y, a=zzzz
        //this->pGSFS[m][10*14 +8] += 0.0;          // i=y, a=zzzz
        this->pGSFS[m][10*14 +9] += 4.0 * zfsds_95; // i=z, a=zzzz
    }
}


void DfTEI::primitiveGSGS(const PrimitiveShellPair& ij,
                          const PrimitiveShellPair& kl,
                          const int nEndM, const int nStartM)
{
    const double QCx = this->QC.x();
    const double QCy = this->QC.y();
    const double QCz = this->QC.z();
    const double WQx = this->WQ.x();
    const double WQy = this->WQ.y();
    const double WQz = this->WQ.z();
    const double ze2 = 1.0 / (2.0 * (ij.zeta + kl.zeta)); // 0.5/(zeta+eta)
    const double eta2 = 1.0 / (2.0 * kl.zeta); //0.5 * inv_eta;
    const double re = this->m_dRho / kl.zeta;

    for (int m = nStartM; m < nEndM; ++m) {
        const int m1 = m +1;

        for (int i = 0; i < 15; ++i) {
            this->pGSGS[m][15*i + 0] = QCx * this->pGSFS[m][10*i +0] +WQx * this->pGSFS[m1][10*i +0]; // xxxx
            this->pGSGS[m][15*i + 1] = QCx * this->pGSFS[m][10*i +1] +WQx * this->pGSFS[m1][10*i +1]; // xxxy
            this->pGSGS[m][15*i + 2] = QCx * this->pGSFS[m][10*i +2] +WQx * this->pGSFS[m1][10*i +2]; // xxxz
            this->pGSGS[m][15*i + 3] = QCx * this->pGSFS[m][10*i +3] +WQx * this->pGSFS[m1][10*i +3]; // xxyy
            this->pGSGS[m][15*i + 4] = QCx * this->pGSFS[m][10*i +4] +WQx * this->pGSFS[m1][10*i +4]; // xxyz
            this->pGSGS[m][15*i + 5] = QCx * this->pGSFS[m][10*i +5] +WQx * this->pGSFS[m1][10*i +5]; // xxzz
            this->pGSGS[m][15*i + 6] = QCx * this->pGSFS[m][10*i +6] +WQx * this->pGSFS[m1][10*i +6]; // xyyy
            this->pGSGS[m][15*i + 7] = QCx * this->pGSFS[m][10*i +7] +WQx * this->pGSFS[m1][10*i +7]; // xyyz
            this->pGSGS[m][15*i + 8] = QCx * this->pGSFS[m][10*i +8] +WQx * this->pGSFS[m1][10*i +8]; // xyzz
            this->pGSGS[m][15*i + 9] = QCx * this->pGSFS[m][10*i +9] +WQx * this->pGSFS[m1][10*i +9]; // xzzz
            this->pGSGS[m][15*i +10] = QCy * this->pGSFS[m][10*i +6] +WQy * this->pGSFS[m1][10*i +6]; // yyyy
            this->pGSGS[m][15*i +11] = QCy * this->pGSFS[m][10*i +7] +WQy * this->pGSFS[m1][10*i +7]; // yyyz
            this->pGSGS[m][15*i +12] = QCy * this->pGSFS[m][10*i +8] +WQy * this->pGSFS[m1][10*i +8]; // yyzz
            this->pGSGS[m][15*i +13] = QCy * this->pGSFS[m][10*i +9] +WQy * this->pGSFS[m1][10*i +9]; // yzzz
            this->pGSGS[m][15*i +14] = QCz * this->pGSFS[m][10*i +9] +WQz * this->pGSFS[m1][10*i +9]; // zzzz
        }

        for (int i = 0; i < 15; ++i) {
            const double zgsds_0 = eta2 * (this->pGSDS[m][6*i +0] -re * this->pGSDS[m1][6*i +0]);
            const double zgsds_1 = eta2 * (this->pGSDS[m][6*i +1] -re * this->pGSDS[m1][6*i +1]);
            const double zgsds_2 = eta2 * (this->pGSDS[m][6*i +2] -re * this->pGSDS[m1][6*i +2]);
            const double zgsds_3 = eta2 * (this->pGSDS[m][6*i +3] -re * this->pGSDS[m1][6*i +3]);
            const double zgsds_4 = eta2 * (this->pGSDS[m][6*i +4] -re * this->pGSDS[m1][6*i +4]);
            const double zgsds_5 = eta2 * (this->pGSDS[m][6*i +5] -re * this->pGSDS[m1][6*i +5]);

            this->pGSGS[m][15*i + 0] += 3.0 * zgsds_0; // i=x, c=xxx
            this->pGSGS[m][15*i + 1] += 2.0 * zgsds_1; // i=x, c=xxy
            this->pGSGS[m][15*i + 2] += 2.0 * zgsds_2; // i=x, c=xxz
            this->pGSGS[m][15*i + 3] +=       zgsds_3; // i=x, c=xyy
            this->pGSGS[m][15*i + 4] +=       zgsds_4; // i=x, c=xyz
            this->pGSGS[m][15*i + 5] +=       zgsds_5; // i=x, c=xzz
            //this->pGSGS[m][15*i + 6] += 0.0;         // i=x, c=yyy
            //this->pGSGS[m][15*i + 7] += 0.0;         // i=x, c=yyz
            //this->pGSGS[m][15*i + 8] += 0.0;         // i=x, c=yzz
            //this->pGSGS[m][15*i + 9] += 0.0;         // i=x, c=zzz
            this->pGSGS[m][15*i +10] += 3.0 * zgsds_3; // i=y, c=yyy
            this->pGSGS[m][15*i +11] += 2.0 * zgsds_4; // i=y, c=yyz
            this->pGSGS[m][15*i +12] +=       zgsds_5; // i=y, c=yzz
            //this->pGSGS[m][15*i +13] += 0.0;         // i=y, c=yzz
            this->pGSGS[m][15*i +14] += 3.0 * zgsds_5; // i=z, c=zzz
        }

        //--------------------
        const double zfsfs_00 = ze2 * this->pFSFS[m1][10*0 +0]; // xxx xxx
        const double zfsfs_01 = ze2 * this->pFSFS[m1][10*0 +1]; // xxx xxy
        const double zfsfs_02 = ze2 * this->pFSFS[m1][10*0 +2]; // xxx xxz
        const double zfsfs_03 = ze2 * this->pFSFS[m1][10*0 +3]; // xxx xyy
        const double zfsfs_04 = ze2 * this->pFSFS[m1][10*0 +4]; // xxx xyz
        const double zfsfs_05 = ze2 * this->pFSFS[m1][10*0 +5]; // xxx xzz
        const double zfsfs_06 = ze2 * this->pFSFS[m1][10*0 +6]; // xxx yyy
        const double zfsfs_07 = ze2 * this->pFSFS[m1][10*0 +7]; // xxx yyz
        const double zfsfs_08 = ze2 * this->pFSFS[m1][10*0 +8]; // xxx yzz
        const double zfsfs_09 = ze2 * this->pFSFS[m1][10*0 +9]; // xxx zzz
        const double zfsfs_10 = ze2 * this->pFSFS[m1][10*1 +0]; // xxy xxx
        const double zfsfs_11 = ze2 * this->pFSFS[m1][10*1 +1]; // xxy xxy
        const double zfsfs_12 = ze2 * this->pFSFS[m1][10*1 +2]; // xxy xxz
        const double zfsfs_13 = ze2 * this->pFSFS[m1][10*1 +3]; // xxy xyy
        const double zfsfs_14 = ze2 * this->pFSFS[m1][10*1 +4]; // xxy xyz
        const double zfsfs_15 = ze2 * this->pFSFS[m1][10*1 +5]; // xxy xzz
        const double zfsfs_16 = ze2 * this->pFSFS[m1][10*1 +6]; // xxy yyy
        const double zfsfs_17 = ze2 * this->pFSFS[m1][10*1 +7]; // xxy yyz
        const double zfsfs_18 = ze2 * this->pFSFS[m1][10*1 +8]; // xxy yzz
        const double zfsfs_19 = ze2 * this->pFSFS[m1][10*1 +9]; // xxy zzz
        const double zfsfs_20 = ze2 * this->pFSFS[m1][10*2 +0]; // xxz xxx
        const double zfsfs_21 = ze2 * this->pFSFS[m1][10*2 +1]; // xxz xxy
        const double zfsfs_22 = ze2 * this->pFSFS[m1][10*2 +2]; // xxz xxz
        const double zfsfs_23 = ze2 * this->pFSFS[m1][10*2 +3]; // xxz xyy
        const double zfsfs_24 = ze2 * this->pFSFS[m1][10*2 +4]; // xxz xyz
        const double zfsfs_25 = ze2 * this->pFSFS[m1][10*2 +5]; // xxz xzz
        const double zfsfs_26 = ze2 * this->pFSFS[m1][10*2 +6]; // xxz yyy
        const double zfsfs_27 = ze2 * this->pFSFS[m1][10*2 +7]; // xxz yyz
        const double zfsfs_28 = ze2 * this->pFSFS[m1][10*2 +8]; // xxz yzz
        const double zfsfs_29 = ze2 * this->pFSFS[m1][10*2 +9]; // xxz zzz
        const double zfsfs_30 = ze2 * this->pFSFS[m1][10*3 +0]; // xyy xxx
        const double zfsfs_31 = ze2 * this->pFSFS[m1][10*3 +1]; // xyy xxy
        const double zfsfs_32 = ze2 * this->pFSFS[m1][10*3 +2]; // xyy xxz
        const double zfsfs_33 = ze2 * this->pFSFS[m1][10*3 +3]; // xyy xyy
        const double zfsfs_34 = ze2 * this->pFSFS[m1][10*3 +4]; // xyy xyz
        const double zfsfs_35 = ze2 * this->pFSFS[m1][10*3 +5]; // xyy xzz
        const double zfsfs_36 = ze2 * this->pFSFS[m1][10*3 +6]; // xyy yyy
        const double zfsfs_37 = ze2 * this->pFSFS[m1][10*3 +7]; // xyy yyz
        const double zfsfs_38 = ze2 * this->pFSFS[m1][10*3 +8]; // xyy yzz
        const double zfsfs_39 = ze2 * this->pFSFS[m1][10*3 +9]; // xyy zzz
        const double zfsfs_40 = ze2 * this->pFSFS[m1][10*4 +0]; // xyz xxx
        const double zfsfs_41 = ze2 * this->pFSFS[m1][10*4 +1]; // xyz xxy
        const double zfsfs_42 = ze2 * this->pFSFS[m1][10*4 +2]; // xyz xxz
        const double zfsfs_43 = ze2 * this->pFSFS[m1][10*4 +3]; // xyz xyy
        const double zfsfs_44 = ze2 * this->pFSFS[m1][10*4 +4]; // xyz xyz
        const double zfsfs_45 = ze2 * this->pFSFS[m1][10*4 +5]; // xyz xzz
        const double zfsfs_46 = ze2 * this->pFSFS[m1][10*4 +6]; // xyz yyy
        const double zfsfs_47 = ze2 * this->pFSFS[m1][10*4 +7]; // xyz yyz
        const double zfsfs_48 = ze2 * this->pFSFS[m1][10*4 +8]; // xyz yzz
        const double zfsfs_49 = ze2 * this->pFSFS[m1][10*4 +9]; // xyz zzz
        const double zfsfs_50 = ze2 * this->pFSFS[m1][10*5 +0]; // xzz xxx
        const double zfsfs_51 = ze2 * this->pFSFS[m1][10*5 +1]; // xzz xxy
        const double zfsfs_52 = ze2 * this->pFSFS[m1][10*5 +2]; // xzz xxz
        const double zfsfs_53 = ze2 * this->pFSFS[m1][10*5 +3]; // xzz xyy
        const double zfsfs_54 = ze2 * this->pFSFS[m1][10*5 +4]; // xzz xyz
        const double zfsfs_55 = ze2 * this->pFSFS[m1][10*5 +5]; // xzz xzz
        const double zfsfs_56 = ze2 * this->pFSFS[m1][10*5 +6]; // xzz yyy
        const double zfsfs_57 = ze2 * this->pFSFS[m1][10*5 +7]; // xzz yyz
        const double zfsfs_58 = ze2 * this->pFSFS[m1][10*5 +8]; // xzz yzz
        const double zfsfs_59 = ze2 * this->pFSFS[m1][10*5 +9]; // xzz zzz
        const double zfsfs_60 = ze2 * this->pFSFS[m1][10*6 +0]; // yyy xxx
        const double zfsfs_61 = ze2 * this->pFSFS[m1][10*6 +1]; // yyy xxy
        const double zfsfs_62 = ze2 * this->pFSFS[m1][10*6 +2]; // yyy xxz
        const double zfsfs_63 = ze2 * this->pFSFS[m1][10*6 +3]; // yyy xyy
        const double zfsfs_64 = ze2 * this->pFSFS[m1][10*6 +4]; // yyy xyz
        const double zfsfs_65 = ze2 * this->pFSFS[m1][10*6 +5]; // yyy xzz
        const double zfsfs_66 = ze2 * this->pFSFS[m1][10*6 +6]; // yyy yyy
        const double zfsfs_67 = ze2 * this->pFSFS[m1][10*6 +7]; // yyy yyz
        const double zfsfs_68 = ze2 * this->pFSFS[m1][10*6 +8]; // yyy yzz
        const double zfsfs_69 = ze2 * this->pFSFS[m1][10*6 +9]; // yyy zzz
        const double zfsfs_70 = ze2 * this->pFSFS[m1][10*7 +0]; // yyz xxx
        const double zfsfs_71 = ze2 * this->pFSFS[m1][10*7 +1]; // yyz xxy
        const double zfsfs_72 = ze2 * this->pFSFS[m1][10*7 +2]; // yyz xxz
        const double zfsfs_73 = ze2 * this->pFSFS[m1][10*7 +3]; // yyz xyy
        const double zfsfs_74 = ze2 * this->pFSFS[m1][10*7 +4]; // yyz xyz
        const double zfsfs_75 = ze2 * this->pFSFS[m1][10*7 +5]; // yyz xzz
        const double zfsfs_76 = ze2 * this->pFSFS[m1][10*7 +6]; // yyz yyy
        const double zfsfs_77 = ze2 * this->pFSFS[m1][10*7 +7]; // yyz yyz
        const double zfsfs_78 = ze2 * this->pFSFS[m1][10*7 +8]; // yyz yzz
        const double zfsfs_79 = ze2 * this->pFSFS[m1][10*7 +9]; // yyz zzz
        const double zfsfs_80 = ze2 * this->pFSFS[m1][10*8 +0]; // yzz xxx
        const double zfsfs_81 = ze2 * this->pFSFS[m1][10*8 +1]; // yzz xxy
        const double zfsfs_82 = ze2 * this->pFSFS[m1][10*8 +2]; // yzz xxz
        const double zfsfs_83 = ze2 * this->pFSFS[m1][10*8 +3]; // yzz xyy
        const double zfsfs_84 = ze2 * this->pFSFS[m1][10*8 +4]; // yzz xyz
        const double zfsfs_85 = ze2 * this->pFSFS[m1][10*8 +5]; // yzz xzz
        const double zfsfs_86 = ze2 * this->pFSFS[m1][10*8 +6]; // yzz yyy
        const double zfsfs_87 = ze2 * this->pFSFS[m1][10*8 +7]; // yzz yyz
        const double zfsfs_88 = ze2 * this->pFSFS[m1][10*8 +8]; // yzz yzz
        const double zfsfs_89 = ze2 * this->pFSFS[m1][10*8 +9]; // yzz zzz
        const double zfsfs_90 = ze2 * this->pFSFS[m1][10*9 +0]; // zzz xxx
        const double zfsfs_91 = ze2 * this->pFSFS[m1][10*9 +1]; // zzz xxy
        const double zfsfs_92 = ze2 * this->pFSFS[m1][10*9 +2]; // zzz xxz
        const double zfsfs_93 = ze2 * this->pFSFS[m1][10*9 +3]; // zzz xyy
        const double zfsfs_94 = ze2 * this->pFSFS[m1][10*9 +4]; // zzz xyz
        const double zfsfs_95 = ze2 * this->pFSFS[m1][10*9 +5]; // zzz xzz
        const double zfsfs_96 = ze2 * this->pFSFS[m1][10*9 +6]; // zzz yyy
        const double zfsfs_97 = ze2 * this->pFSFS[m1][10*9 +7]; // zzz yyz
        const double zfsfs_98 = ze2 * this->pFSFS[m1][10*9 +8]; // zzz yzz
        const double zfsfs_99 = ze2 * this->pFSFS[m1][10*9 +9]; // zzz zzz

        this->pGSGS[m][15* 0 + 0] += 4.0 * zfsfs_00; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 1] += 4.0 * zfsfs_01; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 2] += 4.0 * zfsfs_02; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 3] += 4.0 * zfsfs_03; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 4] += 4.0 * zfsfs_04; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 5] += 4.0 * zfsfs_05; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 6] += 4.0 * zfsfs_06; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 7] += 4.0 * zfsfs_07; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 8] += 4.0 * zfsfs_08; // i=x, a=xxxx
        this->pGSGS[m][15* 0 + 9] += 4.0 * zfsfs_09; // i=x, a=xxxx
        //this->pGSGS[m][15* 0 +10] += 0.0;          // i=y, a=xxxx
        //this->pGSGS[m][15* 0 +11] += 0.0;          // i=y, a=xxxx
        //this->pGSGS[m][15* 0 +12] += 0.0;          // i=y, a=xxxx
        //this->pGSGS[m][15* 0 +13] += 0.0;          // i=y, a=xxxx
        //this->pGSGS[m][15* 0 +14] += 0.0;          // i=z, a=xxxx

        this->pGSGS[m][15* 1 + 0] += 3.0 * zfsfs_10; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 1] += 3.0 * zfsfs_11; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 2] += 3.0 * zfsfs_12; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 3] += 3.0 * zfsfs_13; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 4] += 3.0 * zfsfs_14; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 5] += 3.0 * zfsfs_15; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 6] += 3.0 * zfsfs_16; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 7] += 3.0 * zfsfs_17; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 8] += 3.0 * zfsfs_18; // i=x, a=xxxy
        this->pGSGS[m][15* 1 + 9] += 3.0 * zfsfs_19; // i=x, a=xxxy
        this->pGSGS[m][15* 1 +10] +=       zfsfs_06; // i=y, a=xxxy
        this->pGSGS[m][15* 1 +11] +=       zfsfs_07; // i=y, a=xxxy
        this->pGSGS[m][15* 1 +12] +=       zfsfs_08; // i=y, a=xxxy
        this->pGSGS[m][15* 1 +13] +=       zfsfs_09; // i=y, a=xxxy
        //this->pGSGS[m][15* 1 +14] += 0.0;          // i=z, a=xxxy

        this->pGSGS[m][15* 2 + 0] += 3.0 * zfsfs_20; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 1] += 3.0 * zfsfs_21; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 2] += 3.0 * zfsfs_22; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 3] += 3.0 * zfsfs_23; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 4] += 3.0 * zfsfs_24; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 5] += 3.0 * zfsfs_25; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 6] += 3.0 * zfsfs_26; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 7] += 3.0 * zfsfs_27; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 8] += 3.0 * zfsfs_28; // i=x, a=xxxz
        this->pGSGS[m][15* 2 + 9] += 3.0 * zfsfs_29; // i=x, a=xxxz
        //this->pGSGS[m][15* 2 +10] += 0.0;          // i=y, a=xxxz
        //this->pGSGS[m][15* 2 +11] += 0.0;          // i=y, a=xxxz
        //this->pGSGS[m][15* 2 +12] += 0.0;          // i=y, a=xxxz
        //this->pGSGS[m][15* 2 +13] += 0.0;          // i=y, a=xxxz
        this->pGSGS[m][15* 2 +14] +=       zfsfs_09; // i=z, a=xxxz

        this->pGSGS[m][15* 3 + 0] += 2.0 * zfsfs_30; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 1] += 2.0 * zfsfs_31; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 2] += 2.0 * zfsfs_32; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 3] += 2.0 * zfsfs_33; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 4] += 2.0 * zfsfs_34; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 5] += 2.0 * zfsfs_35; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 6] += 2.0 * zfsfs_36; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 7] += 2.0 * zfsfs_37; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 8] += 2.0 * zfsfs_38; // i=x, a=xxyy
        this->pGSGS[m][15* 3 + 9] += 2.0 * zfsfs_39; // i=x, a=xxyy
        this->pGSGS[m][15* 3 +10] += 2.0 * zfsfs_16; // i=y, a=xxyy
        this->pGSGS[m][15* 3 +11] += 2.0 * zfsfs_17; // i=y, a=xxyy
        this->pGSGS[m][15* 3 +12] += 2.0 * zfsfs_18; // i=y, a=xxyy
        this->pGSGS[m][15* 3 +13] += 2.0 * zfsfs_19; // i=y, a=xxyy
        //this->pGSGS[m][15* 3 +14] += 0.0;          // i=z, a=xxyy

        this->pGSGS[m][15* 4 + 0] += 2.0 * zfsfs_40; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 1] += 2.0 * zfsfs_41; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 2] += 2.0 * zfsfs_42; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 3] += 2.0 * zfsfs_43; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 4] += 2.0 * zfsfs_44; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 5] += 2.0 * zfsfs_45; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 6] += 2.0 * zfsfs_46; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 7] += 2.0 * zfsfs_47; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 8] += 2.0 * zfsfs_48; // i=x, a=xxyz
        this->pGSGS[m][15* 4 + 9] += 2.0 * zfsfs_49; // i=x, a=xxyz
        this->pGSGS[m][15* 4 +10] +=       zfsfs_26; // i=y, a=xxyz
        this->pGSGS[m][15* 4 +11] +=       zfsfs_27; // i=y, a=xxyz
        this->pGSGS[m][15* 4 +12] +=       zfsfs_28; // i=y, a=xxyz
        this->pGSGS[m][15* 4 +13] +=       zfsfs_29; // i=y, a=xxyz
        this->pGSGS[m][15* 4 +14] +=       zfsfs_19; // i=z, a=xxyz

        this->pGSGS[m][15* 5 + 0] += 2.0 * zfsfs_50; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 1] += 2.0 * zfsfs_51; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 2] += 2.0 * zfsfs_52; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 3] += 2.0 * zfsfs_53; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 4] += 2.0 * zfsfs_54; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 5] += 2.0 * zfsfs_55; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 6] += 2.0 * zfsfs_56; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 7] += 2.0 * zfsfs_57; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 8] += 2.0 * zfsfs_58; // i=x, a=xxzz
        this->pGSGS[m][15* 5 + 9] += 2.0 * zfsfs_59; // i=x, a=xxzz
        //this->pGSGS[m][15* 5 +10] += 0.0;          // i=y, a=xxzz
        //this->pGSGS[m][15* 5 +11] += 0.0;          // i=y, a=xxzz
        //this->pGSGS[m][15* 5 +12] += 0.0;          // i=y, a=xxzz
        //this->pGSGS[m][15* 5 +13] += 0.0;          // i=y, a=xxzz
        this->pGSGS[m][15* 5 +14] += 2.0 * zfsfs_29; // i=z, a=xxzz

        this->pGSGS[m][15* 6 + 0] +=       zfsfs_60; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 1] +=       zfsfs_61; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 2] +=       zfsfs_62; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 3] +=       zfsfs_63; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 4] +=       zfsfs_64; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 5] +=       zfsfs_65; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 6] +=       zfsfs_66; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 7] +=       zfsfs_67; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 8] +=       zfsfs_68; // i=x, a=xyyy
        this->pGSGS[m][15* 6 + 9] +=       zfsfs_69; // i=x, a=xyyy
        this->pGSGS[m][15* 6 +10] += 3.0 * zfsfs_36; // i=y, a=xyyy
        this->pGSGS[m][15* 6 +11] += 3.0 * zfsfs_37; // i=y, a=xyyy
        this->pGSGS[m][15* 6 +12] += 3.0 * zfsfs_38; // i=y, a=xyyy
        this->pGSGS[m][15* 6 +13] += 3.0 * zfsfs_39; // i=y, a=xyyy
        //this->pGSGS[m][15* 6 +14] += 0.0;          // i=z, a=xyyy

        this->pGSGS[m][15* 7 + 0] +=       zfsfs_70; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 1] +=       zfsfs_71; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 2] +=       zfsfs_72; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 3] +=       zfsfs_73; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 4] +=       zfsfs_74; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 5] +=       zfsfs_75; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 6] +=       zfsfs_76; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 7] +=       zfsfs_77; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 8] +=       zfsfs_78; // i=x, a=xyyz
        this->pGSGS[m][15* 7 + 9] +=       zfsfs_79; // i=x, a=xyyz
        this->pGSGS[m][15* 7 +10] += 2.0 * zfsfs_46; // i=y, a=xyyz
        this->pGSGS[m][15* 7 +11] += 2.0 * zfsfs_47; // i=y, a=xyyz
        this->pGSGS[m][15* 7 +12] += 2.0 * zfsfs_48; // i=y, a=xyyz
        this->pGSGS[m][15* 7 +13] += 2.0 * zfsfs_49; // i=y, a=xyyz
        this->pGSGS[m][15* 7 +14] +=       zfsfs_39; // i=z, a=xyyz

        this->pGSGS[m][15* 8 + 0] +=       zfsfs_80; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 1] +=       zfsfs_81; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 2] +=       zfsfs_82; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 3] +=       zfsfs_83; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 4] +=       zfsfs_84; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 5] +=       zfsfs_85; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 6] +=       zfsfs_86; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 7] +=       zfsfs_87; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 8] +=       zfsfs_88; // i=x, a=xyzz
        this->pGSGS[m][15* 8 + 9] +=       zfsfs_89; // i=x, a=xyzz
        this->pGSGS[m][15* 8 +10] +=       zfsfs_56; // i=y, a=xyzz
        this->pGSGS[m][15* 8 +11] +=       zfsfs_57; // i=y, a=xyzz
        this->pGSGS[m][15* 8 +12] +=       zfsfs_58; // i=y, a=xyzz
        this->pGSGS[m][15* 8 +13] +=       zfsfs_59; // i=y, a=xyzz
        this->pGSGS[m][15* 8 +14] += 2.0 * zfsfs_49; // i=z, a=xyzz

        this->pGSGS[m][15* 9 + 0] +=       zfsfs_90; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 1] +=       zfsfs_91; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 2] +=       zfsfs_92; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 3] +=       zfsfs_93; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 4] +=       zfsfs_94; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 5] +=       zfsfs_95; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 6] +=       zfsfs_96; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 7] +=       zfsfs_97; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 8] +=       zfsfs_98; // i=x, a=xzzz
        this->pGSGS[m][15* 9 + 9] +=       zfsfs_99; // i=x, a=xzzz
        //this->pGSGS[m][15* 9 +10] += 0.0;          // i=y, a=xzzz
        //this->pGSGS[m][15* 9 +11] += 0.0;          // i=y, a=xzzz
        //this->pGSGS[m][15* 9 +12] += 0.0;          // i=y, a=xzzz
        //this->pGSGS[m][15* 9 +13] += 0.0;          // i=y, a=xzzz
        this->pGSGS[m][15* 9 +14] += 3.0 * zfsfs_59; // i=z, a=xzzz

        //this->pGSGS[m][15*10 + 0] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 1] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 2] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 3] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 4] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 5] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 6] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 7] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 8] += 0.0;          // i=x, a=yyyy
        //this->pGSGS[m][15*10 + 9] += 0.0;          // i=x, a=yyyy
        this->pGSGS[m][15*10 +10] += 4.0 * zfsfs_66; // i=y, a=yyyy
        this->pGSGS[m][15*10 +11] += 4.0 * zfsfs_67; // i=y, a=yyyy
        this->pGSGS[m][15*10 +12] += 4.0 * zfsfs_68; // i=y, a=yyyy
        this->pGSGS[m][15*10 +13] += 4.0 * zfsfs_69; // i=y, a=yyyy
        //this->pGSGS[m][15*10 +14] += 0.0;          // i=z, a=yyyy

        //this->pGSGS[m][15*11 + 0] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 1] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 2] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 3] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 4] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 5] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 6] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 7] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 8] += 0.0;          // i=x, a=yyyz
        //this->pGSGS[m][15*11 + 9] += 0.0;          // i=x, a=yyyz
        this->pGSGS[m][15*11 +10] += 3.0 * zfsfs_76; // i=y, a=yyyz
        this->pGSGS[m][15*11 +11] += 3.0 * zfsfs_77; // i=y, a=yyyz
        this->pGSGS[m][15*11 +12] += 3.0 * zfsfs_78; // i=y, a=yyyz
        this->pGSGS[m][15*11 +13] += 3.0 * zfsfs_79; // i=y, a=yyyz
        this->pGSGS[m][15*11 +14] +=       zfsfs_69; // i=z, a=yyyz

        //this->pGSGS[m][15*12 + 0] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 1] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 2] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 3] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 4] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 5] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 6] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 7] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 8] += 0.0;          // i=x, a=yyzz
        //this->pGSGS[m][15*12 + 9] += 0.0;          // i=x, a=yyzz
        this->pGSGS[m][15*12 +10] += 2.0 * zfsfs_86; // i=y, a=yyzz
        this->pGSGS[m][15*12 +11] += 2.0 * zfsfs_87; // i=y, a=yyzz
        this->pGSGS[m][15*12 +12] += 2.0 * zfsfs_88; // i=y, a=yyzz
        this->pGSGS[m][15*12 +13] += 2.0 * zfsfs_89; // i=y, a=yyzz
        this->pGSGS[m][15*12 +14] += 2.0 * zfsfs_79; // i=z, a=yyzz

        //this->pGSGS[m][15*13 + 0] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 1] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 2] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 3] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 4] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 5] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 6] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 7] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 8] += 0.0;          // i=x, a=yzzz
        //this->pGSGS[m][15*13 + 9] += 0.0;          // i=x, a=yzzz
        this->pGSGS[m][15*13 +10] +=       zfsfs_96; // i=y, a=yzzz
        this->pGSGS[m][15*13 +11] +=       zfsfs_97; // i=y, a=yzzz
        this->pGSGS[m][15*13 +12] +=       zfsfs_98; // i=y, a=yzzz
        this->pGSGS[m][15*13 +13] +=       zfsfs_99; // i=y, a=yzzz
        this->pGSGS[m][15*13 +14] += 3.0 * zfsfs_89; // i=z, a=yzzz

        //this->pGSGS[m][15*14 + 0] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 1] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 2] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 3] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 4] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 5] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 6] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 7] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 8] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 + 9] += 0.0;          // i=x, a=zzzz
        //this->pGSGS[m][15*14 +10] += 0.0;          // i=y, a=zzzz
        //this->pGSGS[m][15*14 +11] += 0.0;          // i=y, a=zzzz
        //this->pGSGS[m][15*14 +12] += 0.0;          // i=y, a=zzzz
        //this->pGSGS[m][15*14 +13] += 0.0;          // i=y, a=zzzz
        this->pGSGS[m][15*14 +14] += 4.0 * zfsfs_99; // i=z, a=zzzz
    }
}

#endif // USE_OLD_TEI_ENGINE
