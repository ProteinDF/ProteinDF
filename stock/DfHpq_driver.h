#ifndef DFHPQ_DRIVER_H
#define DFHPQ_DRIVER_H

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#ifdef HAVE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#define CUDA_GLOBAL __global__
#define CUDA_DEVICE __device__
#define DOUBLE float
#else
#define CUDA_GLOBAL
#define CUDA_DEVICE
#define DOUBLE double
#endif // HAVE_CUDA

    CUDA_GLOBAL void DfHpqDrv_getNuclearAttractionIntegrals(const int nqA, const int nqB,
    const int npA, const int npB,
    const DOUBLE* Za, const DOUBLE* Ca,
    const DOUBLE* Zb, const DOUBLE* Cb,
    const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
    const DOUBLE* TF, const DOUBLE* ADAT,
    const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
    DOUBLE* pE);

    /******************************************************************************/
    /*   Nuclear attraction integral evaluation for type SS                       */
    /*      npA     ; number of first   primitive pair data                       */
    /*      npB     ; number of seconde primitive pair data                       */
    /*      Za      ; orbital exponents alpha                                     */
    /*      Zb      ; orbital exponents beta                                      */
    /*      Ca      ; contracton coefficient A                                    */
    /*      Cb      ; contracton coefficient B                                    */
    /*      Ax,Ay,Az; nuclear coordinates of atom A                               */
    /*      Bx,By,Bz; nuclear coordinates of atom B                               */
    /*      Cx,Cy,Cz; nuclear coordinates of atom C                               */
    /******************************************************************************/
    CUDA_DEVICE void DfHpq_nucSS(const int npA, const int npB,
                                 const DOUBLE* Za, const DOUBLE* Zb,
                                 const DOUBLE* Ca, const DOUBLE* Cb,
                                 const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
                                 const DOUBLE* TF, const DOUBLE* ADAT,
                                 const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
                                 DOUBLE* pE);

    /******************************************************************************/
    /*   Nuclear attraction integral evaluation for type PS                       */
    /*      npA     ; number of first   primitive pair data                       */
    /*      npB     ; number of seconde primitive pair data                       */
    /*      Za      ; orbital exponents alpha                                     */
    /*      Zb      ; orbital exponents beta                                      */
    /*      Ca      ; contracton coefficient A                                    */
    /*      Cb      ; contracton coefficient B                                    */
    /*      Ax,Ay,Az; nuclear coordinates of atom A                               */
    /*      Bx,By,Bz; nuclear coordinates of atom B                               */
    /*      Cx,Cy,Cz; nuclear coordinates of atom C                               */
    /******************************************************************************/
    CUDA_DEVICE void DfHpq_nucPS(const int npA, const int npB,
                                 const DOUBLE* Za, const DOUBLE* Zb,
                                 const DOUBLE* Ca, const DOUBLE* Cb,
                                 const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
                                 const DOUBLE* TF, const DOUBLE* ADAT,
                                 const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
                                 DOUBLE* pE);

    /******************************************************************************/
    /*   Nuclear attraction integral evaluation for type PP                       */
    /*      npA     ; number of first   primitive pair data                       */
    /*      npB     ; number of seconde primitive pair data                       */
    /*      Za      ; orbital exponents alpha                                     */
    /*      Zb      ; orbital exponents beta                                      */
    /*      Ca      ; contracton coefficient A                                    */
    /*      Cb      ; contracton coefficient B                                    */
    /*      Ax,Ay,Az; nuclear coordinates of atom A                               */
    /*      Bx,By,Bz; nuclear coordinates of atom B                               */
    /*      Cx,Cy,Cz; nuclear coordinates of atom C                               */
    /******************************************************************************/
    CUDA_DEVICE void DfHpq_nucPP(const int npA, const int npB,
                                 const DOUBLE* Za, const DOUBLE* Zb,
                                 const DOUBLE* Ca, const DOUBLE* Cb,
                                 const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
                                 const DOUBLE* TF, const DOUBLE* ADAT,
                                 const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
                                 DOUBLE* pE);

    /******************************************************************************/
    /*   Nuclear attraction integral evaluation for type  DS                      */
    /*      npA     ; number of first   primitive pair data                       */
    /*      npB     ; number of seconde primitive pair data                       */
    /*      Za      ; orbital exponents alpha                                     */
    /*      Zb      ; orbital exponents beta                                      */
    /*      Ca      ; contracton coefficient A                                    */
    /*      Cb      ; contracton coefficient B                                    */
    /*      Ax,Ay,Az; nuclear coordinates of atom A                               */
    /*      Bx,By,Bz; nuclear coordinates of atom B                               */
    /*      Cx,Cy,Cz; nuclear coordinates of atom C                               */
    /******************************************************************************/
    CUDA_DEVICE void DfHpq_nucDS(const int npA, int npB,
                                 const DOUBLE* Za, const DOUBLE* Zb,
                                 const DOUBLE* Ca, const DOUBLE* Cb,
                                 const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
                                 const DOUBLE* TF, const DOUBLE* ADAT,
                                 const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
                                 DOUBLE* pE);

    /******************************************************************************/
    /*   Nuclear attraction integral evaluation for type  DP                      */
    /*      npA     ; number of first   primitive pair data                       */
    /*      npB     ; number of seconde primitive pair data                       */
    /*      Za      ; orbital exponents alpha                                     */
    /*      Zb      ; orbital exponents beta                                      */
    /*      Ca      ; contracton coefficient A                                    */
    /*      Cb      ; contracton coefficient B                                    */
    /*      Ax,Ay,Az; nuclear coordinates of atom A                               */
    /*      Bx,By,Bz; nuclear coordinates of atom B                               */
    /*      Cx,Cy,Cz; nuclear coordinates of atom C                               */
    /******************************************************************************/
    CUDA_DEVICE void DfHpq_nucDP(const int npA, const int npB,
                                 const DOUBLE* Za, const DOUBLE* Zb,
                                 const DOUBLE* Ca, const DOUBLE* Cb,
                                 const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
                                 const DOUBLE* TF, const DOUBLE* ADAT,
                                 const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
                                 DOUBLE* pE);


    /******************************************************************************/
    /*   Nuclear attraction integral evaluation for type  DD                      */
    /*      npA     ; number of first   primitive pair data                       */
    /*      npB     ; number of seconde primitive pair data                       */
    /*      Za      ; orbital exponents alpha                                     */
    /*      Zb      ; orbital exponents beta                                      */
    /*      Ca      ; contracton coefficient A                                    */
    /*      Cb      ; contracton coefficient B                                    */
    /*      Ax,Ay,Az; nuclear coordinates of atom A                               */
    /*      Bx,By,Bz; nuclear coordinates of atom B                               */
    /*      Cx,Cy,Cz; nuclear coordinates of atom C                               */
    /******************************************************************************/
    CUDA_DEVICE void DfHpq_nucDD(const int npA, const int npB,
                                 const DOUBLE* Za, const DOUBLE* Zb,
                                 const DOUBLE* Ca, const DOUBLE* Cb,
                                 const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
                                 const DOUBLE* TF, const DOUBLE* ADAT,
                                 const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
                                 DOUBLE* pE);


#ifdef __cplusplus
}
#endif // __cplusplus

#endif // DFHPQ_DRIVER_H
