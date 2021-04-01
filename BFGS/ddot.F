#include "molpro_blas.fh"
#ifdef MOLPRO_BLAS_DDOT
#include "ddot.fh"
#endif

      DOUBLE PRECISION FUNCTION DDOT_X(N_,DX,INCX_,DY,INCY_)
      use molpro_int_mod
      IMPLICIT NONE
      INTEGER                        :: INCX_,INCY_,N_
      DOUBLE PRECISION, DIMENSION(*) :: DX,DY
      DOUBLE PRECISION, EXTERNAL     :: DDOT
#if defined(MOLPRO_BLAS_DDOT) || MOLPRO_INT == BLAS_INT
      DDOT_X=DDOT(N_,DX,INCX_,DY,INCY_)
#else
      INTEGER*BLAS_INT               :: INCX ,INCY ,N
      call mint(INCX_,INCX)
      call mint(INCY_,INCY)
      call mint(N_,N)
      DDOT_X=DDOT(N,DX,INCX,DY,INCY)
#endif
      RETURN
      END
