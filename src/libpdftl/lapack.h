#ifndef LAPACK_H
#define LAPACK_H

extern "C" {
void daxpy_(const int* n, const double* alpha, const double* x, const int* incx,
            double* y, const int* incy);
void dscal_(const int* n, const double* a, double* x, const int* incx);

void dgemv_(const char* TRANS, const int* M, const int* N, const double* ALPHA,
            const double* A, const int* LDA, const double* X, const int* INCX,
            const double* BETA, double* Y, const int* INCY);

void dgemm_(const char* transa, const char* transb, const int* m, const int* n,
            const int* k, const double* alpha, const double* A, const int* ldA,
            const double* B, const int* ldB, const double* beta, double* C,
            const int* ldC);

void dgetrf_(const int* M, const int* N, double* A, const int* LDA, int* IPIV,
             int* INFO);
void dgetri_(const int* N, double* A, const int* LDA, int* IPIV, double* WORK,
             int* LWORK, int* INFO);

void dgelss_(const int* M, const int* N, const int* NRHS, double* A,
             const int* LDA, double* B, const int* LDB, double* S,
             const double* RCOND, int* RANK, double* WORK, const int* LWORK,
             int* INFO);

void dtpttr_(const char* UPLO, const int* N, const double* AP, double* A,
             const int* LDA, int* INFO);
void dtrttp_(const char* UPLO, const int* N, const double* A, const int* LDA,
             double* AP, int* INFO);

void dsymm_(const char* SIDE, const char* UPLO, const int* M, const int* N,
            const double* ALPHA, const double* A, const int* LDA,
            const double* B, const int* LDB, const double* BETA, double* C,
            const int* LDC);

void dspmv_(const char* UPLO, const int* N, const double* ALPHA,
            const double* AP, const double* X, const int* INCX,
            const double* BETA, double* Y, const int* INCY);

void dsptrf_(const char* UPLO, const int* N, double* AP, int* IPIV, int* INFO);
void dsptri_(const char* UPLO, const int* N, double* AP, const int* IPIV,
             double* WORK, int* INFO);

void dspev_(const char* JOBZ, const char* UPLO, const int* N, double* AP,
            double* W, double* Z, const int* LDZ, double* WORK, int* INFO);
}

#endif  // LAPACK_H
