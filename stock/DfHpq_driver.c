#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "DfHpq_driver.h"

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

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
  
static const DOUBLE D33 = 1.0 / 0.03;
static const DOUBLE SQR3I = 0.57735026918962576450914878050196; //1.0 / sqrt(3.0);
static const DOUBLE SQR3I2 = 1.1547005383792515290182975610039; //2.0 / sqrt(3.0);
  
CUDA_GLOBAL void DfHpqDrv_getNuclearAttractionIntegrals(const int nqA, const int nqB,
							const int npA, const int npB,
							const DOUBLE* Za, const DOUBLE* Ca,
							const DOUBLE* Zb, const DOUBLE* Cb,
							const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
							const DOUBLE* TF, const DOUBLE* ADAT,
							const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
							DOUBLE* pE) {
#ifdef HAVE_CUDA
  DOUBLE newC[3];
  newC[0] = C[0] + blockDim.x * blockIdx.x + threadIdx.x;
  newC[1] = C[1] + blockDim.y * blockIdx.y + threadIdx.y;
  newC[2] = C[2] + threadIdx.z;
#else
  DOUBLE newC[3];
  newC[0] = C[0];
  newC[1] = C[1];
  newC[2] = C[2];
#endif

  // s, p, dをサポート
  const int type = 4 * nqA + nqB;

  switch (type) {
  case 0: // SS
    DfHpq_nucSS(npA, npB, Za, Zb, Ca, Cb, A, B, newC,
		TF, ADAT,
		RMI, GA, EDAT,
		pE);
    break;

  case 4: // PS
    DfHpq_nucPS(npA, npB, Za, Zb, Ca, Cb, A, B, newC,
		TF, ADAT,
		RMI, GA, EDAT,
		pE);
    break;
    
  case 5: // PP
    DfHpq_nucPP(npA, npB, Za, Zb, Ca, Cb, A, B, newC,
		TF, ADAT,
		RMI, GA, EDAT,
		pE);
    break;

  case 8: // DS
    DfHpq_nucDS(npA, npB, Za, Zb, Ca, Cb, A, B, newC,
		TF, ADAT,
		RMI, GA, EDAT,
		pE);
    break;

  case 9: // DP
    DfHpq_nucDP(npA, npB, Za, Zb, Ca, Cb, A, B, newC,
		TF, ADAT,
		RMI, GA, EDAT,
		pE);
    break;

  case 10: // DD
    DfHpq_nucDD(npA, npB, Za, Zb, Ca, Cb, A, B, newC,
		TF, ADAT,
		RMI, GA, EDAT,
		pE);
    break;

  default:
#ifndef HAVE_CUDA
    printf("program error: file %s, line %d\n", __FILE__, __LINE__);
    abort();
#endif // HAVE_CUDA
  }
}


CUDA_DEVICE void DfHpq_nucSS(const int npA, const int npB,
			     const DOUBLE* Za, const DOUBLE* Zb,
			     const DOUBLE* Ca, const DOUBLE* Cb,
			     const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
			     const DOUBLE* TF, const DOUBLE* ADAT,
			     const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
			     DOUBLE* pE) {
  const DOUBLE Ax = A[0];
  const DOUBLE Ay = A[1];
  const DOUBLE Az = A[2];
  const DOUBLE Bx = B[0];
  const DOUBLE By = B[1];
  const DOUBLE Bz = B[2];
  const DOUBLE Cx = C[0];
  const DOUBLE Cy = C[1];
  const DOUBLE Cz = C[2];

  const DOUBLE ABx = Ax - Bx;
  const DOUBLE ABy = Ay - By;
  const DOUBLE ABz = Az - Bz;
  const DOUBLE SqAB = ABx*ABx + ABy*ABy + ABz*ABz;

  const DOUBLE TF0 = TF[0];
  const DOUBLE CPAI = M_PI + M_PI;
  DOUBLE SNS = 0.0;
  int i, j;
  for (i = 0; i < npA; ++i) {
    for (j = 0; j < npB; ++j) {
      const DOUBLE Zp = Za[i] + Zb[j];
      const DOUBLE ZpI = 1.0 / Zp;
      const DOUBLE GGI = Za[i]*Zb[j] * ZpI;
      const DOUBLE HP = CPAI*ZpI*exp(-GGI*SqAB)*Ca[i]*Cb[j];
      const DOUBLE Px = (Za[i]*Ax+Zb[j]*Bx)*ZpI;
      const DOUBLE Py = (Za[i]*Ay+Zb[j]*By)*ZpI;
      const DOUBLE Pz = (Za[i]*Az+Zb[j]*Bz)*ZpI;
      const DOUBLE PCx = Px-Cx;
      const DOUBLE PCy = Py-Cy;
      const DOUBLE PCz = Pz-Cz;
      const DOUBLE T = Zp*(PCx*PCx+PCy*PCy+PCz*PCz);
      DOUBLE F0;
      if (T <= TF0) {
	const int IT= (int)((T+0.015)*D33);
	const DOUBLE DT = 0.03 * IT - T;
/* 	F0 = ((((ADAT[ 5][IT] *DT+ADAT[ 4][IT])*DT */
/* 		+ADAT[ 3][IT])*DT+ADAT[ 2][IT])*DT */
/* 	      +  ADAT[ 1][IT])*DT+ADAT[ 0][IT]; */
	F0 = ((((ADAT[ 5*1901+IT] *DT+ADAT[ 4*1901+IT])*DT
		+ADAT[ 3*1901+IT])*DT+ADAT[ 2*1901+IT])*DT
	      +  ADAT[ 1*1901+IT])*DT+ADAT[ 0*1901+IT];
      } else {
	const DOUBLE tinv = 1.0 / T; 
	F0 = sqrt(M_PI * tinv) * 0.5;
      }
      SNS += HP * F0;
    }
  }

  pE[0] = SNS;
}


CUDA_DEVICE void DfHpq_nucPS(const int npA, const int npB,
			     const DOUBLE* Za, const DOUBLE* Zb,
			     const DOUBLE* Ca, const DOUBLE* Cb,
			     const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
			     const DOUBLE* TF, const DOUBLE* ADAT,
			     const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
			     DOUBLE* pE) {
  const DOUBLE Ax=A[0];
  const DOUBLE Ay=A[1];
  const DOUBLE Az=A[2];
  const DOUBLE Bx=B[0];
  const DOUBLE By=B[1];
  const DOUBLE Bz=B[2];
  const DOUBLE Cx = C[0];
  const DOUBLE Cy = C[1];
  const DOUBLE Cz = C[2];

  const DOUBLE ABx=Ax-Bx;
  const DOUBLE ABy=Ay-By;
  const DOUBLE ABz=Az-Bz;
  const DOUBLE SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

  const DOUBLE TFW = TF[1];
  const DOUBLE CPAI = M_PI + M_PI;
  DOUBLE xNS=0.0;
  DOUBLE yNS=0.0;
  DOUBLE zNS=0.0;
  int i, j;
  for (i = 0; i < npA; ++i) {
    for (j = 0; j < npB; ++j) {
      const DOUBLE Zp = Za[i]+Zb[j];
      const DOUBLE ZpI = 1.0 / Zp;
      const DOUBLE GGI = Za[i] * Zb[j] * ZpI;
      const DOUBLE HP =CPAI*ZpI*exp(-GGI*SqAB)*Ca[i]*Cb[j];
      const DOUBLE Px =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
      const DOUBLE Py =(Za[i]*Ay+Zb[j]*By)*ZpI;
      const DOUBLE Pz =(Za[i]*Az+Zb[j]*Bz)*ZpI;
      const DOUBLE PAx=Px-Ax;
      const DOUBLE PAy=Py-Ay;
      const DOUBLE PAz=Pz-Az;
      const DOUBLE PCx=Px-Cx;
      const DOUBLE PCy=Py-Cy;
      const DOUBLE PCz=Pz-Cz;
      const DOUBLE T = Zp*(PCx*PCx+PCy*PCy+PCz*PCz);
      DOUBLE F0, F1;
      if (T <= TFW) {
	const int IT = (int)((T+0.015)*D33);
	const DOUBLE DT = 0.03 * IT - T;
/* 	F1 = ((((ADAT[11][IT] *DT+ADAT[10][IT])*DT */
/* 		+ADAT[ 9][IT])*DT+ADAT[ 8][IT])*DT */
/* 	      +  ADAT[ 7][IT])*DT+ADAT[ 6][IT];  */
	F1 = ((((ADAT[11*1901+IT] *DT+ADAT[10*1901+IT])*DT
		+ADAT[ 9*1901+IT])*DT+ADAT[ 8*1901+IT])*DT
	      +  ADAT[ 7*1901+IT])*DT+ADAT[ 6*1901+IT]; 
	const DOUBLE EED = ((((GA[5]*DT+GA[4])*DT+GA[3])*DT+GA[2])*DT+GA[1])*DT+GA[0];
	const DOUBLE EE = EDAT[IT]*EED;
	F0 = 2.0 * T * F1 + EE;
      } else { 
	const DOUBLE TINV=1.0 / T;
	F0 =sqrt(M_PI * TINV) * 0.5;
	F1 =0.5*F0*TINV;
      }
      const DOUBLE SS0 = HP * F0;
      const DOUBLE SS1 = HP * F1;

      xNS += PAx*SS0 - PCx*SS1;
      yNS += PAy*SS0 - PCy*SS1;
      zNS += PAz*SS0 - PCz*SS1;
    }
  }

  pE[0] = xNS;
  pE[1] = yNS;
  pE[2] = zNS;
}


CUDA_DEVICE void DfHpq_nucPP(const int npA, const int npB,
			     const DOUBLE* Za, const DOUBLE* Zb,
			     const DOUBLE* Ca, const DOUBLE* Cb,
			     const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
			     const DOUBLE* TF, const DOUBLE* ADAT,
			     const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
			     DOUBLE* pE) {
  const DOUBLE Ax=A[0];
  const DOUBLE Ay=A[1];
  const DOUBLE Az=A[2];
  const DOUBLE Bx=B[0];
  const DOUBLE By=B[1];
  const DOUBLE Bz=B[2];
  const DOUBLE Cx = C[0];
  const DOUBLE Cy = C[1];
  const DOUBLE Cz = C[2];

  const DOUBLE ABx=Ax-Bx;
  const DOUBLE ABy=Ay-By;
  const DOUBLE ABz=Az-Bz;
  const DOUBLE SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

  const DOUBLE TFW =TF[2];
  const DOUBLE RMI1=RMI[0];
  const DOUBLE CPAI=M_PI+M_PI;
  DOUBLE xNx =0.0;
  DOUBLE xNy =0.0;
  DOUBLE xNz =0.0;
  DOUBLE yNx =0.0;
  DOUBLE yNy =0.0;
  DOUBLE yNz =0.0;
  DOUBLE zNx =0.0;
  DOUBLE zNy =0.0;
  DOUBLE zNz =0.0;
  int i, j;
  for (i = 0; i < npA; ++i) {
    for (j = 0; j < npB; ++j) {
      const DOUBLE Zp =Za[i]+Zb[j];
      const DOUBLE ZpI=1.0/Zp;
      const DOUBLE GGI=Za[i]*Zb[j] * ZpI;
      const DOUBLE HP =CPAI*ZpI*exp(-GGI*SqAB)*Ca[i]*Cb[j];
      const DOUBLE Px =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
      const DOUBLE Py =(Za[i]*Ay+Zb[j]*By)*ZpI;
      const DOUBLE Pz =(Za[i]*Az+Zb[j]*Bz)*ZpI;
      const DOUBLE PAx=Px-Ax;
      const DOUBLE PAy=Py-Ay;
      const DOUBLE PAz=Pz-Az;
      const DOUBLE PBx=Px-Bx;
      const DOUBLE PBy=Py-By;
      const DOUBLE PBz=Pz-Bz;
      const DOUBLE PCx=Px-Cx;
      const DOUBLE PCy=Py-Cy;
      const DOUBLE PCz=Pz-Cz;
      const DOUBLE T =Zp*(PCx*PCx+PCy*PCy+PCz*PCz);
      DOUBLE F0, F1, F2;
      if (T <= TFW) {
	const int IT =(int)((T+0.015)*D33);
	const DOUBLE DT =0.03*IT-T;
/* 	F2 =((((ADAT[17][IT] *DT + ADAT[16][IT])*DT */
/* 	       +ADAT[15][IT])*DT + ADAT[14][IT])*DT */
/* 	     +  ADAT[13][IT])*DT + ADAT[12][IT]; */
	F2 =((((ADAT[17*1901+IT] *DT + ADAT[16*1901+IT])*DT
	       +ADAT[15*1901+IT])*DT + ADAT[14*1901+IT])*DT
	     +  ADAT[13*1901+IT])*DT + ADAT[12*1901+IT];
	const DOUBLE EED=((((GA[5]*DT+GA[4])*DT+GA[3])*DT+GA[2])*DT+GA[1])*DT+GA[0];
	const DOUBLE EE =EDAT[IT]*EED;
	const DOUBLE T2 =T * 2.0;
	F1 =(T2*F2 + EE)*RMI1;
	F0 =T2*F1 + EE;
      }else{ 
	const DOUBLE TINV=1.0/T;
	F0  =sqrt(M_PI*TINV)*0.5;
	F1  =0.5*F0*TINV;
	F2  =1.5*F1*TINV;
      }
      const DOUBLE SS0=HP*F0;
      const DOUBLE SS1=HP*F1;
      const DOUBLE SS2=HP*F2;
      // PS
      const DOUBLE xS0=PAx*SS0-PCx*SS1;
      const DOUBLE yS0=PAy*SS0-PCy*SS1;
      const DOUBLE zS0=PAz*SS0-PCz*SS1;
      const DOUBLE xS1=PAx*SS1-PCx*SS2;
      const DOUBLE yS1=PAy*SS1-PCy*SS2;
      const DOUBLE zS1=PAz*SS1-PCz*SS2;
      // PP
      const DOUBLE ZpI5 =0.5*ZpI;
      const DOUBLE ZpS01=ZpI5*(SS0 - SS1);
      xNx += PBx*xS0-PCx*xS1+ZpS01;
      xNy += PBy*xS0-PCy*xS1;
      xNz += PBz*xS0-PCz*xS1;
      yNx += PBx*yS0-PCx*yS1;
      yNy += PBy*yS0-PCy*yS1+ZpS01;
      yNz += PBz*yS0-PCz*yS1;
      zNx += PBx*zS0-PCx*zS1;
      zNy += PBy*zS0-PCy*zS1;
      zNz += PBz*zS0-PCz*zS1+ZpS01;
    }
  }     

  pE[0]=xNx;
  pE[1]=xNy;
  pE[2]=xNz;
  pE[3]=yNx;
  pE[4]=yNy;
  pE[5]=yNz;
  pE[6]=zNx;
  pE[7]=zNy;
  pE[8]=zNz;
}


CUDA_DEVICE void DfHpq_nucDS(const int npA, int npB,
			     const DOUBLE* Za, const DOUBLE* Zb,
			     const DOUBLE* Ca, const DOUBLE* Cb,
			     const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
			     const DOUBLE* TF, const DOUBLE* ADAT,
			     const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
			     DOUBLE* pE) {
  const DOUBLE Ax=A[0];
  const DOUBLE Ay=A[1];
  const DOUBLE Az=A[2];
  const DOUBLE Bx=B[0];
  const DOUBLE By=B[1];
  const DOUBLE Bz=B[2];
  const DOUBLE Cx = C[0];
  const DOUBLE Cy = C[1];
  const DOUBLE Cz = C[2];

  const DOUBLE ABx=Ax-Bx;
  const DOUBLE ABy=Ay-By;
  const DOUBLE ABz=Az-Bz;
  const DOUBLE SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

  const DOUBLE TFW =TF[2];
  const DOUBLE RMI1=RMI[0];
  const DOUBLE CPAI=M_PI+M_PI;
  DOUBLE AS0 =0.0;
  DOUBLE BS0 =0.0;
  DOUBLE CS0 =0.0;
  DOUBLE DS0 =0.0;
  DOUBLE ES0 =0.0;
  int i, j;
  for (i = 0; i < npA; ++i) {
    for (j = 0; j < npB; ++j) {
      const DOUBLE Zp =Za[i]+Zb[j];
      const DOUBLE ZpI=1.0/Zp;
      const DOUBLE GGI=Za[i]*Zb[j]*ZpI;
      const DOUBLE HP =CPAI*ZpI*exp(-GGI*SqAB)*Ca[i]*Cb[j];
      const DOUBLE Px =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
      const DOUBLE Py =(Za[i]*Ay+Zb[j]*By)*ZpI;
      const DOUBLE Pz =(Za[i]*Az+Zb[j]*Bz)*ZpI;
      const DOUBLE PAx=Px-Ax;
      const DOUBLE PAy=Py-Ay;
      const DOUBLE PAz=Pz-Az;
      const DOUBLE PCx=Px-Cx;
      const DOUBLE PCy=Py-Cy;
      const DOUBLE PCz=Pz-Cz;
      const DOUBLE T  =Zp*(PCx*PCx+PCy*PCy+PCz*PCz);
      DOUBLE F0, F1, F2;
      if (T <= TFW) {
	int IT =(int)((T+0.015)*D33);
	const DOUBLE DT =0.03*IT-T;
/* 	F2 =((((ADAT[17][IT] *DT+ADAT[16][IT])*DT */
/* 	       +ADAT[15][IT])*DT+ADAT[14][IT])*DT */
/* 	     +  ADAT[13][IT])*DT+ADAT[12][IT]; */
	F2 =((((ADAT[17*1901+IT] *DT+ADAT[16*1901+IT])*DT
	       +ADAT[15*1901+IT])*DT+ADAT[14*1901+IT])*DT
	     +  ADAT[13*1901+IT])*DT+ADAT[12*1901+IT];
	const DOUBLE EED=((((GA[5]*DT+GA[4])*DT+GA[3])*DT+GA[2])*DT+GA[1])*DT+GA[0];
	const DOUBLE EE =EDAT[IT]*EED;
	const DOUBLE T2 =T*2.0;
	F1 =(T2*F2+EE)*RMI1;
	F0 =T2*F1+EE;
      } else {
	const DOUBLE TINV=1.0/T;
	F0  =sqrt(M_PI*TINV)*0.5;
	F1  =0.5*F0*TINV;
	F2  =1.5*F1*TINV;
      }
      const DOUBLE SS0 =HP*F0;
      const DOUBLE SS1 =HP*F1;
      const DOUBLE SS2 =HP*F2;
      // PS
      const DOUBLE xS0=PAx*SS0-PCx*SS1;
      const DOUBLE yS0=PAy*SS0-PCy*SS1;
      const DOUBLE zS0=PAz*SS0-PCz*SS1;
      const DOUBLE xS1=PAx*SS1-PCx*SS2;
      const DOUBLE yS1=PAy*SS1-PCy*SS2;
      const DOUBLE zS1=PAz*SS1-PCz*SS2;
      // DS     A=XY, B=XZ, C=YZ, D=XX-YY, E=3ZZ-RR
      const DOUBLE xxSW=PAx*xS0-PCx*xS1;
      const DOUBLE yySW=PAy*yS0-PCy*yS1;
      const DOUBLE zzSW=PAz*zS0-PCz*zS1;
      AS0 += PAx*yS0-PCx*yS1;
      BS0 += PAx*zS0-PCx*zS1;
      CS0 += PAy*zS0-PCy*zS1;
      DS0 += 0.5*(xxSW-yySW);
      ES0 += SQR3I*(zzSW-0.5*(xxSW+yySW));
    }
  }

  pE[0]=AS0;
  pE[1]=BS0;
  pE[2]=CS0;
  pE[3]=DS0;
  pE[4]=ES0;
}


CUDA_DEVICE void DfHpq_nucDP(const int npA, const int npB,
			     const DOUBLE* Za, const DOUBLE* Zb,
			     const DOUBLE* Ca, const DOUBLE* Cb,
			     const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
			     const DOUBLE* TF, const DOUBLE* ADAT,
			     const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
			     DOUBLE* pE) {
  const DOUBLE Ax=A[0];
  const DOUBLE Ay=A[1];
  const DOUBLE Az=A[2];
  const DOUBLE Bx=B[0];
  const DOUBLE By=B[1];
  const DOUBLE Bz=B[2];
  const DOUBLE Cx = C[0];
  const DOUBLE Cy = C[1];
  const DOUBLE Cz = C[2];

  const DOUBLE ABx=Ax-Bx;
  const DOUBLE ABy=Ay-By;
  const DOUBLE ABz=Az-Bz;
  const DOUBLE SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

  const DOUBLE TFW =TF[3];
  const DOUBLE RMI1=RMI[0];
  const DOUBLE RMI2=RMI[1];
  const DOUBLE CPAI=M_PI+M_PI;
  DOUBLE Ax0 =0.0;
  DOUBLE Ay0 =0.0;
  DOUBLE Az0 =0.0;
  DOUBLE Bx0 =0.0;
  DOUBLE By0 =0.0;
  DOUBLE Bz0 =0.0;
  DOUBLE Cx0 =0.0;
  DOUBLE Cy0 =0.0;
  DOUBLE Cz0 =0.0;
  DOUBLE Dx0 =0.0;
  DOUBLE Dy0 =0.0;
  DOUBLE Dz0 =0.0;
  DOUBLE Ex0 =0.0;
  DOUBLE Ey0 =0.0;
  DOUBLE Ez0 =0.0;
  int i, j;
  for (i = 0; i < npA; ++i) {
    for (j = 0; j < npB; ++j) {
      const DOUBLE Zp =Za[i]+Zb[j];
      const DOUBLE ZpI=1.0/Zp;
      const DOUBLE GGI=Za[i]*Zb[j] * ZpI;
      const DOUBLE HP =CPAI*ZpI*exp(-GGI*SqAB)*Ca[i]*Cb[j];
      const DOUBLE Px =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
      const DOUBLE Py =(Za[i]*Ay+Zb[j]*By)*ZpI;
      const DOUBLE Pz =(Za[i]*Az+Zb[j]*Bz)*ZpI;
      const DOUBLE PAx=Px-Ax;
      const DOUBLE PAy=Py-Ay;
      const DOUBLE PAz=Pz-Az;
      const DOUBLE PBx=Px-Bx;
      const DOUBLE PBy=Py-By;
      const DOUBLE PBz=Pz-Bz;
      const DOUBLE PCx=Px-Cx;
      const DOUBLE PCy=Py-Cy;
      const DOUBLE PCz=Pz-Cz;
      const DOUBLE T  =Zp*(PCx*PCx+PCy*PCy+PCz*PCz);
      DOUBLE F0, F1, F2, F3;
      if (T <= TFW) {
	const int IT =(int)((T+0.015)*D33);
	const DOUBLE DT =0.03*IT-T;
/* 	F3 =((((ADAT[23][IT] *DT+ADAT[22][IT])*DT */
/* 	       +ADAT[21][IT])*DT+ADAT[20][IT])*DT */
/* 	     +  ADAT[19][IT])*DT+ADAT[18][IT]; */
	F3 =((((ADAT[23*1901+IT] *DT+ADAT[22*1901+IT])*DT
	       +ADAT[21*1901+IT])*DT+ADAT[20*1901+IT])*DT
	     +  ADAT[19*1901+IT])*DT+ADAT[18*1901+IT];
	const DOUBLE EED=((((GA[5]*DT+GA[4])*DT+GA[3])*DT+GA[2])*DT+GA[1])*DT+GA[0];
	const DOUBLE EE =EDAT[IT]*EED;
	const DOUBLE T2 =T*2.0;
	F2 =(T2*F3+EE)*RMI2;
	F1 =(T2*F2+EE)*RMI1;
	F0 =T2*F1+EE;
      } else {
	const DOUBLE TINV=1.0/T;
	F0  =sqrt(M_PI*TINV)*0.5;
	F1  =0.5*F0*TINV;
	F2  =1.5*F1*TINV;
	F3  =2.5*F2*TINV;
      }
      const DOUBLE SS0 =HP*F0;
      const DOUBLE SS1 =HP*F1;
      const DOUBLE SS2 =HP*F2;
      const DOUBLE SS3 =HP*F3;
      const DOUBLE ZpI5=ZpI*0.5;
      // PS
      const DOUBLE xS0=PAx*SS0-PCx*SS1;
      const DOUBLE yS0=PAy*SS0-PCy*SS1;
      const DOUBLE zS0=PAz*SS0-PCz*SS1;
      const DOUBLE xS1=PAx*SS1-PCx*SS2;
      const DOUBLE yS1=PAy*SS1-PCy*SS2;
      const DOUBLE zS1=PAz*SS1-PCz*SS2;
      const DOUBLE xS2=PAx*SS2-PCx*SS3;
      const DOUBLE yS2=PAy*SS2-PCy*SS3;
      const DOUBLE zS2=PAz*SS2-PCz*SS3;
      // DS     A=XY, B=XZ, C=YZ, D=XX-YY, E=3ZZ-RR
      const DOUBLE xxSW1=PAx*xS0-PCx*xS1;
      const DOUBLE yySW1=PAy*yS0-PCy*yS1;
      const DOUBLE zzSW1=PAz*zS0-PCz*zS1;
      const DOUBLE AS0 =PAx*yS0-PCx*yS1;
      const DOUBLE BS0 =PAx*zS0-PCx*zS1;
      const DOUBLE CS0 =PAy*zS0-PCy*zS1;
      const DOUBLE DS0 =0.5*(xxSW1-yySW1);
      const DOUBLE ES0 =SQR3I*(zzSW1-0.5*(xxSW1+yySW1));
      const DOUBLE xxSW2=PAx*xS1-PCx*xS2;
      const DOUBLE yySW2=PAy*yS1-PCy*yS2;
      const DOUBLE zzSW2=PAz*zS1-PCz*zS2;
      const DOUBLE AS1 =PAx*yS1-PCx*yS2;
      const DOUBLE BS1 =PAx*zS1-PCx*zS2;
      const DOUBLE CS1 =PAy*zS1-PCy*zS2;
      const DOUBLE DS1 =0.5*(xxSW2-yySW2);
      const DOUBLE ES1 =SQR3I*(zzSW2-0.5*(xxSW2+yySW2));
      // DP
      const DOUBLE RxS01=ZpI5*(xS0-xS1);
      const DOUBLE RyS01=ZpI5*(yS0-yS1);
      const DOUBLE RzS01=ZpI5*(zS0-zS1);
      Ax0 += PBx*AS0-PCx*AS1+RyS01;
      Ay0 += PBy*AS0-PCy*AS1+RxS01;
      Az0 += PBz*AS0-PCz*AS1;
      Bx0 += PBx*BS0-PCx*BS1+RzS01;
      By0 += PBy*BS0-PCy*BS1;
      Bz0 += PBz*BS0-PCz*BS1+RxS01;
      Cx0 += PBx*CS0-PCx*CS1;
      Cy0 += PBy*CS0-PCy*CS1+RzS01;
      Cz0 += PBz*CS0-PCz*CS1+RyS01;
      Dx0 += PBx*DS0-PCx*DS1+RxS01;
      Dy0 += PBy*DS0-PCy*DS1-RyS01;
      Dz0 += PBz*DS0-PCz*DS1;
      Ex0 += PBx*ES0-PCx*ES1-RxS01*SQR3I;
      Ey0 += PBy*ES0-PCy*ES1-RyS01*SQR3I;
      Ez0 += PBz*ES0-PCz*ES1+RzS01*SQR3I2;
    }
  }

  pE[ 0]=Ax0;
  pE[ 1]=Ay0;
  pE[ 2]=Az0;
  pE[ 3]=Bx0;
  pE[ 4]=By0;
  pE[ 5]=Bz0;
  pE[ 6]=Cx0;
  pE[ 7]=Cy0;
  pE[ 8]=Cz0;
  pE[ 9]=Dx0;
  pE[10]=Dy0;
  pE[11]=Dz0;
  pE[12]=Ex0;
  pE[13]=Ey0;
  pE[14]=Ez0;
}


CUDA_DEVICE void DfHpq_nucDD(const int npA, const int npB,
			     const DOUBLE* Za, const DOUBLE* Zb,
			     const DOUBLE* Ca, const DOUBLE* Cb,
			     const DOUBLE* A, const DOUBLE* B, const DOUBLE* C,
			     const DOUBLE* TF, const DOUBLE* ADAT,
			     const DOUBLE* RMI, const DOUBLE* GA, const DOUBLE* EDAT,
			     DOUBLE* pE) {
  const DOUBLE Ax=A[0];
  const DOUBLE Ay=A[1];
  const DOUBLE Az=A[2];
  const DOUBLE Bx=B[0];
  const DOUBLE By=B[1];
  const DOUBLE Bz=B[2];
  const DOUBLE Cx = C[0];
  const DOUBLE Cy = C[1];
  const DOUBLE Cz = C[2];

  const DOUBLE ABx=Ax-Bx;
  const DOUBLE ABy=Ay-By;
  const DOUBLE ABz=Az-Bz;
  const DOUBLE SqAB=ABx*ABx+ABy*ABy+ABz*ABz;

  const DOUBLE TFW =TF[4];
  const DOUBLE RMI1=RMI[0];
  const DOUBLE RMI2=RMI[1];
  const DOUBLE RMI3=RMI[2];
  const DOUBLE CPAI=M_PI+M_PI;
  DOUBLE AA0 =0.0;
  DOUBLE AB0 =0.0;
  DOUBLE AC0 =0.0;
  DOUBLE AD0 =0.0;
  DOUBLE AE0 =0.0;
  DOUBLE BA0 =0.0;
  DOUBLE BB0 =0.0;
  DOUBLE BC0 =0.0;
  DOUBLE BD0 =0.0;
  DOUBLE BE0 =0.0;
  DOUBLE CA0 =0.0;
  DOUBLE CB0 =0.0;
  DOUBLE CC0 =0.0;
  DOUBLE CD0 =0.0;
  DOUBLE CE0 =0.0;
  DOUBLE DA0 =0.0;
  DOUBLE DB0 =0.0;
  DOUBLE DC0 =0.0;
  DOUBLE DD0 =0.0;
  DOUBLE DE0 =0.0;
  DOUBLE EA0 =0.0;
  DOUBLE EB0 =0.0;
  DOUBLE EC0 =0.0;
  DOUBLE ED0 =0.0;
  DOUBLE EE0 =0.0;
  int i, j;
  for (i = 0; i < npA; ++i) {
    for (j = 0; j < npB; ++j) {
      const DOUBLE Zp =Za[i]+Zb[j];
      const DOUBLE ZpI=1.0/Zp;
      const DOUBLE GGI=Za[i]*Zb[j] * ZpI;
      const DOUBLE HP =CPAI*ZpI*exp(-GGI*SqAB)*Ca[i]*Cb[j];
      const DOUBLE Px =(Za[i]*Ax+Zb[j]*Bx)*ZpI;
      const DOUBLE Py =(Za[i]*Ay+Zb[j]*By)*ZpI;
      const DOUBLE Pz =(Za[i]*Az+Zb[j]*Bz)*ZpI;
      const DOUBLE PAx=Px-Ax;
      const DOUBLE PAy=Py-Ay;
      const DOUBLE PAz=Pz-Az;
      const DOUBLE PBx=Px-Bx;
      const DOUBLE PBy=Py-By;
      const DOUBLE PBz=Pz-Bz;
      const DOUBLE PCx=Px-Cx;
      const DOUBLE PCy=Py-Cy;
      const DOUBLE PCz=Pz-Cz;
      const DOUBLE T =Zp*(PCx*PCx+PCy*PCy+PCz*PCz);
      DOUBLE F0, F1, F2, F3, F4;
      if (T <= TFW) { 
	const int IT =(int)((T+0.015)*D33);
	const DOUBLE DT =0.03*IT-T;
/* 	F4 =((((ADAT[29][IT] *DT  + ADAT[28][IT])*DT */
/* 	       +ADAT[27][IT])*DT  + ADAT[26][IT])*DT */
/* 	     +  ADAT[25][IT])*DT  + ADAT[24][IT]; */
	F4 =((((ADAT[29*1901+IT] *DT  + ADAT[28*1901+IT])*DT
	       +ADAT[27*1901+IT])*DT  + ADAT[26*1901+IT])*DT
	     +  ADAT[25*1901+IT])*DT  + ADAT[24*1901+IT];
	const DOUBLE EED=((((GA[5]*DT+GA[4])*DT+GA[3])*DT+GA[2])*DT+GA[1])*DT+GA[0];
	const DOUBLE EE =EDAT[IT]*EED;
	const DOUBLE T2 =T*2.0;
	F3 =(T2*F4+EE)*RMI3;
	F2 =(T2*F3+EE)*RMI2;
	F1 =(T2*F2+EE)*RMI1;
	F0 = T2*F1+EE;
      } else {
	const DOUBLE TINV = 1.0 / T;
	F0 = sqrt(M_PI*TINV)*0.5;
	F1 = 0.5*F0*TINV;
	F2 = 1.5*F1*TINV;
	F3 = 2.5*F2*TINV;
	F4 = 3.5*F3*TINV;
      }
      const DOUBLE SS0 = HP*F0;
      const DOUBLE SS1 = HP*F1;
      const DOUBLE SS2 = HP*F2;
      const DOUBLE SS3 = HP*F3;
      const DOUBLE SS4 = HP*F4;
      const DOUBLE ZpI5= ZpI*0.5;
      // PS
      const DOUBLE xS0=PAx*SS0-PCx*SS1;
      const DOUBLE yS0=PAy*SS0-PCy*SS1;
      const DOUBLE zS0=PAz*SS0-PCz*SS1;
      const DOUBLE xS1=PAx*SS1-PCx*SS2;
      const DOUBLE yS1=PAy*SS1-PCy*SS2;
      const DOUBLE zS1=PAz*SS1-PCz*SS2;
      const DOUBLE xS2=PAx*SS2-PCx*SS3;
      const DOUBLE yS2=PAy*SS2-PCy*SS3;
      const DOUBLE zS2=PAz*SS2-PCz*SS3;
      const DOUBLE xS3=PAx*SS3-PCx*SS4;
      const DOUBLE yS3=PAy*SS3-PCy*SS4;
      const DOUBLE zS3=PAz*SS3-PCz*SS4;
      // PP
      const DOUBLE ZpS01=ZpI5*(SS0-SS1);
      const DOUBLE xx0  =PBx*xS0-PCx*xS1+ZpS01;
      const DOUBLE xy0  =PBy*xS0-PCy*xS1;
      const DOUBLE xz0  =PBz*xS0-PCz*xS1;
      const DOUBLE yx0  =PBx*yS0-PCx*yS1;
      const DOUBLE yy0  =PBy*yS0-PCy*yS1+ZpS01;
      const DOUBLE yz0  =PBz*yS0-PCz*yS1;
      const DOUBLE zx0  =PBx*zS0-PCx*zS1;
      const DOUBLE zy0  =PBy*zS0-PCy*zS1;
      const DOUBLE zz0  =PBz*zS0-PCz*zS1+ZpS01;
      const DOUBLE ZpS12=ZpI5*(SS1-SS2);
      const DOUBLE xx1  =PBx*xS1-PCx*xS2+ZpS12;
      const DOUBLE xy1  =PBy*xS1-PCy*xS2;
      const DOUBLE xz1  =PBz*xS1-PCz*xS2;
      const DOUBLE yx1  =PBx*yS1-PCx*yS2;
      const DOUBLE yy1  =PBy*yS1-PCy*yS2+ZpS12;
      const DOUBLE yz1  =PBz*yS1-PCz*yS2;
      const DOUBLE zx1  =PBx*zS1-PCx*zS2;
      const DOUBLE zy1  =PBy*zS1-PCy*zS2;
      const DOUBLE zz1  =PBz*zS1-PCz*zS2+ZpS12;
      // DS     A=XY, B=XZ, C=YZ, D=XX-YY, E=3ZZ-RR
      const DOUBLE xxSW1=PAx*xS0-PCx*xS1;
      const DOUBLE yySW1=PAy*yS0-PCy*yS1;
      const DOUBLE zzSW1=PAz*zS0-PCz*zS1;
      const DOUBLE AS0 =PAx*yS0-PCx*yS1;
      const DOUBLE BS0 =PAx*zS0-PCx*zS1;
      const DOUBLE CS0 =PAy*zS0-PCy*zS1;
      const DOUBLE DS0 =0.5*(xxSW1-yySW1);
      const DOUBLE ES0 =SQR3I*(zzSW1-0.5*(xxSW1+yySW1));
      const DOUBLE xxSW2=PAx*xS1-PCx*xS2;
      const DOUBLE yySW2=PAy*yS1-PCy*yS2;
      const DOUBLE zzSW2=PAz*zS1-PCz*zS2;
      const DOUBLE AS1 =PAx*yS1-PCx*yS2;
      const DOUBLE BS1 =PAx*zS1-PCx*zS2;
      const DOUBLE CS1 =PAy*zS1-PCy*zS2;
      const DOUBLE DS1 =0.5*(xxSW2-yySW2);
      const DOUBLE ES1 =SQR3I*(zzSW2-0.5*(xxSW2+yySW2));
      const DOUBLE xxSW3=PAx*xS2-PCx*xS3;
      const DOUBLE yySW3=PAy*yS2-PCy*yS3;
      const DOUBLE zzSW3=PAz*zS2-PCz*zS3;
      const DOUBLE AS2 =PAx*yS2-PCx*yS3;
      const DOUBLE BS2 =PAx*zS2-PCx*zS3;
      const DOUBLE CS2 =PAy*zS2-PCy*zS3;
      const DOUBLE DS2 =0.5*(xxSW3-yySW3);
      const DOUBLE ES2 =SQR3I*(zzSW3-0.5*(xxSW3+yySW3));
      // DP
      const DOUBLE RxS01=ZpI5*(xS0-xS1);
      const DOUBLE RyS01=ZpI5*(yS0-yS1);
      const DOUBLE RzS01=ZpI5*(zS0-zS1);
      const DOUBLE Ax0  =PBx*AS0-PCx*AS1+RyS01;
      const DOUBLE Ay0  =PBy*AS0-PCy*AS1+RxS01;
      const DOUBLE Az0  =PBz*AS0-PCz*AS1;
      const DOUBLE Bx0  =PBx*BS0-PCx*BS1+RzS01;
      const DOUBLE By0  =PBy*BS0-PCy*BS1;
      const DOUBLE Bz0  =PBz*BS0-PCz*BS1+RxS01;
      const DOUBLE Cx0  =PBx*CS0-PCx*CS1;
      const DOUBLE Cy0  =PBy*CS0-PCy*CS1+RzS01;
      const DOUBLE Cz0  =PBz*CS0-PCz*CS1+RyS01;
      const DOUBLE Dx0  =PBx*DS0-PCx*DS1+RxS01;
      const DOUBLE Dy0  =PBy*DS0-PCy*DS1-RyS01;
      const DOUBLE Dz0  =PBz*DS0-PCz*DS1;
      const DOUBLE Ex0  =PBx*ES0-PCx*ES1-RxS01*SQR3I;
      const DOUBLE Ey0  =PBy*ES0-PCy*ES1-RyS01*SQR3I;
      const DOUBLE Ez0  =PBz*ES0-PCz*ES1+RzS01*SQR3I2;
      const DOUBLE RxS12=ZpI5*(xS1-xS2);
      const DOUBLE RyS12=ZpI5*(yS1-yS2);
      const DOUBLE RzS12=ZpI5*(zS1-zS2);
      const DOUBLE Ax1  =PBx*AS1-PCx*AS2+RyS12;
      const DOUBLE Ay1  =PBy*AS1-PCy*AS2+RxS12;
      const DOUBLE Az1  =PBz*AS1-PCz*AS2;
      const DOUBLE Bx1  =PBx*BS1-PCx*BS2+RzS12;
      const DOUBLE By1  =PBy*BS1-PCy*BS2;
      const DOUBLE Bz1  =PBz*BS1-PCz*BS2+RxS12;
      const DOUBLE Cx1  =PBx*CS1-PCx*CS2;
      const DOUBLE Cy1  =PBy*CS1-PCy*CS2+RzS12;
      const DOUBLE Cz1  =PBz*CS1-PCz*CS2+RyS12;
      const DOUBLE Dx1  =PBx*DS1-PCx*DS2+RxS12;
      const DOUBLE Dy1  =PBy*DS1-PCy*DS2-RyS12;
      const DOUBLE Dz1  =PBz*DS1-PCz*DS2;
      const DOUBLE Ex1  =PBx*ES1-PCx*ES2-RxS12*SQR3I;
      const DOUBLE Ey1  =PBy*ES1-PCy*ES2-RyS12*SQR3I;
      const DOUBLE Ez1  =PBz*ES1-PCz*ES2+RzS12*SQR3I2;
      // DD
      const DOUBLE Rxx01=ZpI5*(xx0-xx1);
      const DOUBLE Rxy01=ZpI5*(xy0-xy1);
      const DOUBLE Rxz01=ZpI5*(xz0-xz1);
      const DOUBLE Ryx01=ZpI5*(yx0-yx1);
      const DOUBLE Ryy01=ZpI5*(yy0-yy1);
      const DOUBLE Ryz01=ZpI5*(yz0-yz1);
      const DOUBLE Rzx01=ZpI5*(zx0-zx1);
      const DOUBLE Rzy01=ZpI5*(zy0-zy1);
      const DOUBLE Rzz01=ZpI5*(zz0-zz1);
      const DOUBLE AxxW =PBx*Ax0-PCx*Ax1+Ryx01;
      const DOUBLE AyyW =PBy*Ay0-PCy*Ay1+Rxy01;
      const DOUBLE AzzW =PBz*Az0-PCz*Az1;
      const DOUBLE BxxW =PBx*Bx0-PCx*Bx1+Rzx01;
      const DOUBLE ByyW =PBy*By0-PCy*By1;
      const DOUBLE BzzW =PBz*Bz0-PCz*Bz1+Rxz01;
      const DOUBLE CxxW =PBx*Cx0-PCx*Cx1;
      const DOUBLE CyyW =PBy*Cy0-PCy*Cy1+Rzy01;
      const DOUBLE CzzW =PBz*Cz0-PCz*Cz1+Ryz01;
      const DOUBLE DxxW =PBx*Dx0-PCx*Dx1+Rxx01;
      const DOUBLE DyyW =PBy*Dy0-PCy*Dy1-Ryy01;
      const DOUBLE DzzW =PBz*Dz0-PCz*Dz1;
      const DOUBLE ExxW =PBx*Ex0-PCx*Ex1-Rxx01*SQR3I;
      const DOUBLE EyyW =PBy*Ey0-PCy*Ey1-Ryy01*SQR3I;
      const DOUBLE EzzW =PBz*Ez0-PCz*Ez1+Rzz01*SQR3I2;

      AA0 += PBx*Ay0-PCx*Ay1+Ryy01;
      AB0 += PBx*Az0-PCx*Az1+Ryz01;
      AC0 += PBz*Ay0-PCz*Ay1;
      AD0 += 0.5*(AxxW-AyyW);
      AE0 += SQR3I*(AzzW-0.5*(AxxW+AyyW));
      BA0 += PBx*By0-PCx*By1+Rzy01;
      BB0 += PBx*Bz0-PCx*Bz1+Rzz01;
      BC0 += PBz*By0-PCz*By1+Rxy01;
      BD0 += 0.5*(BxxW-ByyW);
      BE0 += SQR3I*(BzzW-0.5*(BxxW+ByyW));
      CA0 += PBx*Cy0-PCx*Cy1;
      CB0 += PBx*Cz0-PCx*Cz1;
      CC0 += PBz*Cy0-PCz*Cy1+Ryy01;
      CD0 += 0.5*(CxxW-CyyW);
      CE0 += SQR3I*(CzzW-0.5*(CxxW+CyyW));
      DA0 += PBx*Dy0-PCx*Dy1+Rxy01;
      DB0 += PBx*Dz0-PCx*Dz1+Rxz01;
      DC0 += PBz*Dy0-PCz*Dy1;
      DD0 += 0.5*(DxxW-DyyW);
      DE0 += SQR3I*(DzzW-0.5*(DxxW+DyyW));
      EA0 += PBx*Ey0-PCx*Ey1-Rxy01*SQR3I;
      EB0 += PBx*Ez0-PCx*Ez1-Rxz01*SQR3I;
      EC0 += PBz*Ey0-PCz*Ey1+Rzy01*SQR3I2;
      ED0 += 0.5*(ExxW-EyyW);
      EE0 += SQR3I*(EzzW-0.5*(ExxW+EyyW));
    }
  } 

  pE[ 0]=AA0;
  pE[ 1]=AB0;
  pE[ 2]=AC0;
  pE[ 3]=AD0;
  pE[ 4]=AE0;
  pE[ 5]=BA0;
  pE[ 6]=BB0;
  pE[ 7]=BC0;
  pE[ 8]=BD0;
  pE[ 9]=BE0;
  pE[10]=CA0;
  pE[11]=CB0;
  pE[12]=CC0;
  pE[13]=CD0;
  pE[14]=CE0;
  pE[15]=DA0;
  pE[16]=DB0;
  pE[17]=DC0;
  pE[18]=DD0;
  pE[19]=DE0;
  pE[20]=EA0;
  pE[21]=EB0;
  pE[22]=EC0;
  pE[23]=ED0;
  pE[24]=EE0;
}

#ifdef __cplusplus
}
#endif //__cplusplus
