cstart none
c;c--------------------------------------------------------------------
c;      subroutine mxma(a,mcola,mrowa,b,mcolb,mrowb,
c;     *r,mcolr,mrowr,ncol,nlink,nrow)
c;c.....simple inner prodoct version
c;c--------------------------------------------------------------------
c;      implicit double precision (a-h,o-z)
c;      dimension r(1),a(1),b(1)
c;c... r(ncol,nrow)=a(ncol,nlink)*b(nlink,nrow) matrix mult
c;      iaa=1
c;      irr=1
c;      do 300 i=1,ncol
c;      ibb=1
c;      ir=irr
c;      do 200 j=1,nrow
c;      ib=ibb
c;      ia=iaa
c;      s=0
c;      do 100 k=1,nlink
c;      s=a(ia)*b(ib)+s
c;      ib=ib+mcolb
c;100   ia=ia+mrowa
c;      r(ir)=s
c;      ir=ir+mrowr
c;200   ibb=ibb+mrowb
c;      irr=irr+mcolr
c;300   iaa=iaa+mcola
c;      return
c;      end
c;c--------------------------------------------------------------------
c;      subroutine mxmb(a,mcola,mrowa,b,mcolb,mrowb,
c;     *r,mcolr,mrowr,ncol,nlink,nrow)
c;c.....simple inner prodoct version
c;c--------------------------------------------------------------------
c;      implicit double precision (a-h,o-z)
c;      dimension r(1),a(1),b(1)
c;c... r(ncol,nrow)=r(ncol,nrow)+a(ncol,nlink)*b(nlink,nrow) matrix mult
c;      iaa=1
c;      irr=1
c;      do 300 i=1,ncol
c;      ibb=1
c;      ir=irr
c;      do 200 j=1,nrow
c;      ib=ibb
c;      ia=iaa
c;      s=0
c;      do 100 k=1,nlink
c;      s=a(ia)*b(ib)+s
c;      ib=ib+mcolb
c;100   ia=ia+mrowa
c;      r(ir)=r(ir)+s
c;      ir=ir+mrowr
c;200   ibb=ibb+mrowb
c;      irr=irr+mcolr
c;300   iaa=iaa+mcola
c;      return
c;      end
c;c--------------------------------------------------------------------
c;      subroutine mxman(a,mcola,mrowa,b,mcolb,mrowb,
c;     *r,mcolr,mrowr,ncol,nlink,nrow)
c;c.....simple inner prodoct version
c;c--------------------------------------------------------------------
c;      implicit double precision (a-h,o-z)
c;      dimension r(1),a(1),b(1)
c;c... r(ncol,nrow)=-a(ncol,nlink)*b(nlink,nrow) matrix mult
c;      iaa=1
c;      irr=1
c;      do 300 i=1,ncol
c;      ibb=1
c;      ir=irr
c;      do 200 j=1,nrow
c;      ib=ibb
c;      ia=iaa
c;      s=0
c;      do 100 k=1,nlink
c;      s=a(ia)*b(ib)+s
c;      ib=ib+mcolb
c;100   ia=ia+mrowa
c;      r(ir)=-s
c;      ir=ir+mrowr
c;200   ibb=ibb+mrowb
c;      irr=irr+mcolr
c;300   iaa=iaa+mcola
c;      return
c;      end
c;c--------------------------------------------------------------------
c;      subroutine mxmbn(a,mcola,mrowa,b,mcolb,mrowb,
c;     *r,mcolr,mrowr,ncol,nlink,nrow)
c;c.....simple inner prodoct version
c;c--------------------------------------------------------------------
c;      implicit double precision (a-h,o-z)
c;      dimension r(1),a(1),b(1)
c;c... r(ncol,nrow)=r(ncol,nrow)-a(ncol,nlink)*b(nlink,nrow) matrix mult
c;      iaa=1
c;      irr=1
c;      do 300 i=1,ncol
c;      ibb=1
c;      ir=irr
c;      do 200 j=1,nrow
c;      ib=ibb
c;      ia=iaa
c;      s=0
c;      do 100 k=1,nlink
c;      s=a(ia)*b(ib)+s
c;      ib=ib+mcolb
c;100   ia=ia+mrowa
c;      r(ir)=r(ir)-s
c;      ir=ir+mrowr
c;200   ibb=ibb+mrowb
c;      irr=irr+mcolr
c;300   iaa=iaa+mcola
c;      return
c;      end
cend

#ifdef MOLPRO_sxf90
c--------------------------------------------------------------------
      subroutine mxma(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
c... r(ncol,nrow)=a(ncol,nlink)*b(nlink,nrow) matrix mult
c     calling NEC internal vector routine
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
      if(ncol.eq.0.or.nrow.eq.0) return
      if(nlink.eq.0) then
        do j=1,nrow
          do i=1,ncol
            r((i-1)*mcolr+(j-1)*mrowr+1) = 0.0d0
          enddo
        enddo
      else
        call VDMXMA(a,mcola,mrowa,b,mcolb,mrowb,
     &              r,mcolr,mrowr,ncol,nlink,nrow)
      end if
      return
      end
c--------------------------------------------------------------------
      subroutine mxmb(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
c... r(ncol,nrow)=r(ncol,nrow)+a(ncol,nlink)*b(nlink,nrow) matrix mult
c     calling NEC internal vector routine
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
      if( ncol.eq.0.or.nrow.eq.0.or.nlink.eq.0) return
        call VDMXQA(a,mcola,mrowa,b,mcolb,mrowb,
     &              r,mcolr,mrowr,ncol,nlink,nrow)
      return
      end
c--------------------------------------------------------------------
      subroutine mxman(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
c... r(ncol,nrow)=-a(ncol,nlink)*b(nlink,nrow) matrix mult
c     primitive coding, let the compiler insert vector routines
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
      if(ncol.eq.0.or.nrow.eq.0) return
      do j=1,nrow
        do i=1,ncol
          r((i-1)*mcolr+(j-1)*mrowr+1) = 0.0d0
        enddo
      enddo
      if(nlink.eq.0) return
      do j=1,nrow
        do l=1,nlink
          do i=1,ncol
            r((i-1)*mcolr+(j-1)*mrowr+1)=
     1      r((i-1)*mcolr+(j-1)*mrowr+1) -
     2      a((i-1)*mcola+(l-1)*mrowa+1) *
     3      b((l-1)*mcolb+(j-1)*mrowb+1)
          enddo
        enddo
      enddo
      return
      end
c--------------------------------------------------------------------
      subroutine mxmbn(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
c... r(ncol,nrow)=r(ncol,nrow)-a(ncol,nlink)*b(nlink,nrow) matrix mult
c     primitive coding, let the compiler insert vector routines
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
      if( ncol.eq.0.or.nrow.eq.0.or.nlink.eq.0) return
      do j=1,nrow
        do l=1,nlink
          do i=1,ncol
            r((i-1)*mcolr+(j-1)*mrowr+1)=
     1      r((i-1)*mcolr+(j-1)*mrowr+1) -
     2      a((i-1)*mcola+(l-1)*mrowa+1) *
     3      b((l-1)*mcolb+(j-1)*mrowb+1)
          enddo
        enddo
      enddo
      return
      end
#else

!---------------------------------------------
!     Commented by me down to here
!---------------------------------------------




!> matrix-matrix product
!>
!> calculates the \f${\tt ncol}\times{\tt nrow}\f$ matrix product
!> \f$ \bf R = A B\f$.
!>
!> \param a real array containing the matrix \f$A\f$; note that the
!>     matrix may be supplied transposed, with non-unit stride, and with
!>     enlarged leading dimension, under the control of the parameters
!>     \c mcola_, \c mrowa described below.
!> \param mcola_ memory increment in \c a between adjacent rows
!>     of \f$A\f$.
!> \param mrowa_ memory increment in \c a between adjacent columns of
!>     \f$A\f$.
!> \param b real array containing the matrix \f$B\f$
!> \param mcolb_ memory increment in \c b between adjacent rows
!>     of \f$B\f$.
!> \param mrowb_ memory increment in \c b between adjacent columns of
!>     \f$B\f$.
!> \param r real array containing the matrix \f$R\f$; \c r is cleared
!>     to zero at the start of \c mxma, and so any previous contents
!>     are lost.
!> \param mcolr_ memory increment in \c r between adjacent rows
!>     of \f$R\f$.
!> \param mrowr_ memory increment in \c r between adjacent columns of
!>     \f$R\f$.
!> \param ncol_ number of rows in \f$R\f$ and \f$A\f$.
!> \param nlink_ dimension which is summed over: number of rows in
!>     \f$B\f$ and number of columns in \f$A\f$.
!> \param nrow_ number of columns in \f$R\f$ and \f$B\f$.
!>
!> Thus standard use of \c mxma would take the form \f$ A*B=R: A(ncol,nlink)*B(nlink,nrow)=R(ncol,nrow)\f$:
!> \verbatim call mxma(a,1,ncol,b,1,nlink,r,1,ncol,ncol,nlink,nrow) \endverbatim
!> whilst, for example \f$ A*B(T)=R: A(ncol,nlink)*B(nrow,nlink)=R(ncol,nrow)\f$:
!> \verbatim call mxma(a,1,ncol,b,nrow,1,r,1,ncol,ncol,nlink,nrow) \endverbatim
!> and \f$ A(T)*B=R:  A(nlink,ncol)*B(nlink,nrow)=R(ncol,nrow)\f$:
!> \verbatim call mxma(a,nlink,1, b,1, nlink, r,1,ncol, ncol,nlink,nrow) \endverbatim
c--------------------------------------------------------------------
      subroutine mxma(a,mcola_,mrowa_,b,mcolb_,mrowb_,
     1                r,mcolr_,mrowr_,ncol_,nlink_,nrow_)
c--------------------------------------------------------------------
c.....unrolled version for ibm6000 and other risc machines
      implicit double precision (a-h,o-z)
      integer mcola_,mrowa_,mcolb_,mrowb_,mcolr_,mrowr_,
     >        ncol_,nlink_,nrow_
    !  include "common/clseg"
      parameter (maxb=128)
      dimension buf(maxb*maxb)
      double precision, dimension(*) :: r, a, b
      if(ncol_.eq.0.or.nrow_.eq.0) return
      nrow=nrow_
      ncol=ncol_
      nlink=nlink_
      mcolr=mcolr_
      mrowr=mrowr_
      if(nlink.eq.0) then
        ijj=1
        do 6201 j=1,nrow
        ij=ijj
        do 6101 i=1,ncol
        r(ij)=0.0d0
6101    ij=ij+mcolr
6201    ijj=ijj+mrowr
        return
      end if
      mcola=mcola_
      mrowa=mrowa_
      mcolb=mcolb_
      mrowb=mrowb_
      if(ncol.gt.5.or.nrow.gt.5) goto 6000
      goto (1000,2000,3000,4000,5000),nrow
c
c....nrow=1
1000  goto (1001,1002,1003,1004,1005),ncol
1001    ia1=1
        ib1=1
        s1=0.0d0
        do 1010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
1010    ib1=ib1+mcolb
        r(1)=s1
        return
c
1002    ia1=1
        ib1=1
        ia2=1+mcola
        s1=0.0d0
        s2=0.0d0
        do 1020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
1020    ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        return
c
1003    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        do 1030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
1030    ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        return
c
1004    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        do 1040 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
1040    ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        r(1+3*mcolr)=s4
        return
c
1005    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        do 1050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
        s5=a(ia5)*b(ib1)+s5
        ia5=ia5+mrowa
1050    ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        r(1+3*mcolr)=s4
        r(1+4*mcolr)=s5
        return
c
c....nrow=2
2000  goto (2001,2002,2003,2004,2005),ncol
2001    ia1=1
        ib1=1
        ib2=1+mrowb
        s1=0.0d0
        s2=0.0d0
        do 2010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
2010    ia1=ia1+mrowa
        r(1)=s1
        r(1+mrowr)=s2
        return
c
2002    ia1=1
        ib1=1
        ia2=1+mcola
        ib2=1+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        do 2020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia2)*b(ib1)+s2
        ib1=ib1+mcolb
        s3=a(ia1)*b(ib2)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib2)+s4
        ib2=ib2+mcolb
2020    ia2=ia2+mrowa
        r(1)=s1
        r(1+mcolr)=s2
        r(1+mrowr)=s3
        r(1+mcolr+mrowr)=s4
        return
c
2003    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        do 2030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s4=a(ia1)*b(ib2)+s4
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        s5=a(ia2)*b(ib2)+s5
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ib1=ib1+mcolb
        s6=a(ia3)*b(ib2)+s6
        ia3=ia3+mrowa
2030    ib2=ib2+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        r(1+mrowr)=s4
        r(1+mrowr+mcolr)=s5
        r(1+mrowr+2*mcolr)=s6
        return
c
2004    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        do 2040 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s42=a(ia4)*b(ib2)+s42
        ia4=ia4+mrowa
2040    ib2=ib2+mcolb
        r(1)=s11
        r(1+mcolr)=s21
        r(1+2*mcolr)=s31
        r(1+3*mcolr)=s41
        r(1+mrowr)=s12
        r(1+mrowr+mcolr)=s22
        r(1+mrowr+2*mcolr)=s32
        r(1+mrowr+3*mcolr)=s42
        return
c
2005    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s51=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s52=0.0d0
        do 2050 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        s42=a(ia4)*b(ib2)+s42
        ia4=ia4+mrowa
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s52=a(ia5)*b(ib2)+s52
        ia5=ia5+mrowa
2050    ib2=ib2+mcolb
        r(1)=s11
        r(1+mcolr)=s21
        r(1+2*mcolr)=s31
        r(1+3*mcolr)=s41
        r(1+4*mcolr)=s51
        r(1+mrowr)=s12
        r(1+mrowr+mcolr)=s22
        r(1+mrowr+2*mcolr)=s32
        r(1+mrowr+3*mcolr)=s42
        r(1+mrowr+4*mcolr)=s52
        return
c
c....nrow=3
3000  goto (3001,3002,3003,3004,3005),ncol
3001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        do 3010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
3010    ib3=ib3+mcolb
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        return
c
3002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        do 3020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s4=a(ia2)*b(ib1)+s4
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        s5=a(ia2)*b(ib2)+s5
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib3)+s6
        ib3=ib3+mcolb
3020    ia2=ia2+mrowa
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+mcolr)=s4
        r(1+mcolr+mrowr)=s5
        r(1+mcolr+2*mrowr)=s6
        return
c
3003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        s7=0.0d0
        s8=0.0d0
        s9=0.0d0
        do 3030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib1)+s4
        s5=a(ia2)*b(ib2)+s5
        s6=a(ia2)*b(ib3)+s6
        ia2=ia2+mrowa
        s7=a(ia3)*b(ib1)+s7
        ib1=ib1+mcolb
        s8=a(ia3)*b(ib2)+s8
        ib2=ib2+mcolb
        s9=a(ia3)*b(ib3)+s9
        ia3=ia3+mrowa
3030    ib3=ib3+mcolb
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+mcolr)=s4
        r(1+mcolr+mrowr)=s5
        r(1+mcolr+2*mrowr)=s6
        r(1+2*mcolr)=s7
        r(1+2*mcolr+mrowr)=s8
        r(1+2*mcolr+2*mrowr)=s9
        return
c
3004    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        s7=0.0d0
        s8=0.0d0
        s9=0.0d0
        s10=0.0d0
        s11=0.0d0
        s12=0.0d0
        do 3040 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib1)+s4
        s5=a(ia2)*b(ib2)+s5
        s6=a(ia2)*b(ib3)+s6
        ia2=ia2+mrowa
        s7=a(ia3)*b(ib1)+s7
        s8=a(ia3)*b(ib2)+s8
        s9=a(ia3)*b(ib3)+s9
        ia3=ia3+mrowa
        s10=a(ia4)*b(ib1)+s10
        ib1=ib1+mcolb
        s11=a(ia4)*b(ib2)+s11
        ib2=ib2+mcolb
        s12=a(ia4)*b(ib3)+s12
        ia4=ia4+mrowa
3040    ib3=ib3+mcolb
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+mcolr)=s4
        r(1+mcolr+mrowr)=s5
        r(1+mcolr+2*mrowr)=s6
        r(1+2*mcolr)=s7
        r(1+2*mcolr+mrowr)=s8
        r(1+2*mcolr+2*mrowr)=s9
        r(1+3*mcolr)=s10
        r(1+3*mcolr+mrowr)=s11
        r(1+3*mcolr+2*mrowr)=s12
        return
c
3005    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s51=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s52=0.0d0
        s13=0.0d0
        s23=0.0d0
        s33=0.0d0
        s43=0.0d0
        s53=0.0d0
        do 3050 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        s13=a(ia1)*b(ib3)+s13
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        s23=a(ia2)*b(ib3)+s23
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        s33=a(ia3)*b(ib3)+s33
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        s42=a(ia4)*b(ib2)+s42
        s43=a(ia4)*b(ib3)+s43
        ia4=ia4+mrowa
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        s53=a(ia5)*b(ib3)+s53
        ia5=ia5+mrowa
3050    ib3=ib3+mcolb
        r(1)=s11
        r(1+mcolr)=s21
        r(1+2*mcolr)=s31
        r(1+3*mcolr)=s41
        r(1+4*mcolr)=s51
        r(1+mrowr)=s12
        r(1+mrowr+mcolr)=s22
        r(1+mrowr+2*mcolr)=s32
        r(1+mrowr+3*mcolr)=s42
        r(1+mrowr+4*mcolr)=s52
        r(1+2*mrowr)=s13
        r(1+2*mrowr+mcolr)=s23
        r(1+2*mrowr+2*mcolr)=s33
        r(1+2*mrowr+3*mcolr)=s43
        r(1+2*mrowr+4*mcolr)=s53
        return
c
c.....nrow=4
4000  goto(4001,4002,4003,4004,4005),ncol
4001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        do 4010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
        s4=a(ia1)*b(ib4)+s4
        ib4=ib4+mcolb
4010    ia1=ia1+mrowa
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+3*mrowr)=s4
        return
c
4002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        do 4020 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ib4=ib4+mcolb
4020    ia2=ia2+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        return
c
4003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s31=0.0d0
        s32=0.0d0
        s33=0.0d0
        s34=0.0d0
        do 4030 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ia2=ia2+mrowa
        s34=a(ia3)*b(ib4)+s34
        ib4=ib4+mcolb
4030    ia3=ia3+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        r(1+2*mcolr)=s31
        r(1+2*mcolr+mrowr)=s32
        r(1+2*mcolr+2*mrowr)=s33
        r(1+2*mcolr+3*mrowr)=s34
        return
c
4004    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s31=0.0d0
        s32=0.0d0
        s33=0.0d0
        s34=0.0d0
        s41=0.0d0
        s42=0.0d0
        s43=0.0d0
        s44=0.0d0
        do 4040 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ia2=ia2+mrowa
        s34=a(ia3)*b(ib4)+s34
        ia3=ia3+mrowa
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
4040    ia4=ia4+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        r(1+2*mcolr)=s31
        r(1+2*mcolr+mrowr)=s32
        r(1+2*mcolr+2*mrowr)=s33
        r(1+2*mcolr+3*mrowr)=s34
        r(1+3*mcolr)=s41
        r(1+3*mcolr+mrowr)=s42
        r(1+3*mcolr+2*mrowr)=s43
        r(1+3*mcolr+3*mrowr)=s44
        return
c
4005    ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s6=0.0d0
        s7=0.0d0
        s8=0.0d0
        s9=0.0d0
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s16=0.0d0
        s17=0.0d0
        s18=0.0d0
        s19=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        do 4050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        s4=a(ia1)*b(ib4)+s4
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib1)+s6
        s7=a(ia2)*b(ib2)+s7
        s8=a(ia2)*b(ib3)+s8
        s9=a(ia2)*b(ib4)+s9
        ia2=ia2+mrowa
        s11=a(ia3)*b(ib1)+s11
        s12=a(ia3)*b(ib2)+s12
        s13=a(ia3)*b(ib3)+s13
        s14=a(ia3)*b(ib4)+s14
        ia3=ia3+mrowa
        s16=a(ia4)*b(ib1)+s16
        s17=a(ia4)*b(ib2)+s17
        s18=a(ia4)*b(ib3)+s18
        s19=a(ia4)*b(ib4)+s19
        ia4=ia4+mrowa
        s21=a(ia5)*b(ib1)+s21
        s22=a(ia5)*b(ib2)+s22
        s23=a(ia5)*b(ib3)+s23
        s24=a(ia5)*b(ib4)+s24
        ia5=ia5+mrowa
        ib1=ib1+mcolb
        ib2=ib2+mcolb
        ib3=ib3+mcolb
        ib4=ib4+mcolb
4050    continue
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+3*mrowr)=s4
        r(1+mcolr)=s6
        r(1+mcolr+mrowr)=s7
        r(1+mcolr+2*mrowr)=s8
        r(1+mcolr+3*mrowr)=s9
        r(1+2*mcolr)=s11
        r(1+2*mcolr+mrowr)=s12
        r(1+2*mcolr+2*mrowr)=s13
        r(1+2*mcolr+3*mrowr)=s14
        r(1+3*mcolr)=s16
        r(1+3*mcolr+mrowr)=s17
        r(1+3*mcolr+2*mrowr)=s18
        r(1+3*mcolr+3*mrowr)=s19
        r(1+4*mcolr)=s21
        r(1+4*mcolr+mrowr)=s22
        r(1+4*mcolr+2*mrowr)=s23
        r(1+4*mcolr+3*mrowr)=s24
        return
c
c.....nrow=5
5000  goto(5001,5002,5003,5004,5004),ncol
5001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        do 5010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
        s4=a(ia1)*b(ib4)+s4
        ib4=ib4+mcolb
        s5=a(ia1)*b(ib5)+s5
        ib5=ib5+mcolb
5010    ia1=ia1+mrowa
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+3*mrowr)=s4
        r(1+4*mrowr)=s5
        return
c
5002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=1
        ia2=ia1+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s15=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s25=0.0d0
        do 5020 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        ia1=ia1+mrowa
        s25=a(ia2)*b(ib5)+s25
        ib5=ib5+mcolb
5020    ia2=ia2+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+4*mrowr)=s15
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        r(1+mcolr+4*mrowr)=s25
        return
c
5003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s15=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s25=0.0d0
        s31=0.0d0
        s32=0.0d0
        s33=0.0d0
        s34=0.0d0
        s35=0.0d0
        do 5030 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        ia1=ia1+mrowa
        s25=a(ia2)*b(ib5)+s25
        ia2=ia2+mrowa
        s35=a(ia3)*b(ib5)+s35
        ib5=ib5+mcolb
5030    ia3=ia3+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+4*mrowr)=s15
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        r(1+mcolr+4*mrowr)=s25
        r(1+2*mcolr)=s31
        r(1+2*mcolr+mrowr)=s32
        r(1+2*mcolr+2*mrowr)=s33
        r(1+2*mcolr+3*mrowr)=s34
        r(1+2*mcolr+4*mrowr)=s35
        return
c
5004    ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        s7=0.0d0
        s8=0.0d0
        s9=0.0d0
        s10=0.0d0
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s15=0.0d0
        s16=0.0d0
        s17=0.0d0
        s18=0.0d0
        s19=0.0d0
        s20=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s25=0.0d0
        do 5050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        s4=a(ia1)*b(ib4)+s4
        s5=a(ia1)*b(ib5)+s5
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib1)+s6
        s7=a(ia2)*b(ib2)+s7
        s8=a(ia2)*b(ib3)+s8
        s9=a(ia2)*b(ib4)+s9
        s10=a(ia2)*b(ib5)+s10
        ia2=ia2+mrowa
        s11=a(ia3)*b(ib1)+s11
        s12=a(ia3)*b(ib2)+s12
        s13=a(ia3)*b(ib3)+s13
        s14=a(ia3)*b(ib4)+s14
        s15=a(ia3)*b(ib5)+s15
        ia3=ia3+mrowa
        s16=a(ia4)*b(ib1)+s16
        s17=a(ia4)*b(ib2)+s17
        s18=a(ia4)*b(ib3)+s18
        s19=a(ia4)*b(ib4)+s19
        s20=a(ia4)*b(ib5)+s20
        ia4=ia4+mrowa
        if(ncol.eq.4) goto 5040
        s21=a(ia5)*b(ib1)+s21
        s22=a(ia5)*b(ib2)+s22
        s23=a(ia5)*b(ib3)+s23
        s24=a(ia5)*b(ib4)+s24
        s25=a(ia5)*b(ib5)+s25
        ia5=ia5+mrowa
5040    ib1=ib1+mcolb
        ib2=ib2+mcolb
        ib3=ib3+mcolb
        ib4=ib4+mcolb
5050    ib5=ib5+mcolb
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+3*mrowr)=s4
        r(1+4*mrowr)=s5
        r(1+mcolr)=s6
        r(1+mcolr+mrowr)=s7
        r(1+mcolr+2*mrowr)=s8
        r(1+mcolr+3*mrowr)=s9
        r(1+mcolr+4*mrowr)=s10
        r(1+2*mcolr)=s11
        r(1+2*mcolr+mrowr)=s12
        r(1+2*mcolr+2*mrowr)=s13
        r(1+2*mcolr+3*mrowr)=s14
        r(1+2*mcolr+4*mrowr)=s15
        r(1+3*mcolr)=s16
        r(1+3*mcolr+mrowr)=s17
        r(1+3*mcolr+2*mrowr)=s18
        r(1+3*mcolr+3*mrowr)=s19
        r(1+3*mcolr+4*mrowr)=s20
        if(ncol.eq.4) return
        r(1+4*mcolr)=s21
        r(1+4*mcolr+mrowr)=s22
        r(1+4*mcolr+2*mrowr)=s23
        r(1+4*mcolr+3*mrowr)=s24
        r(1+4*mcolr+4*mrowr)=s25
        return
c
6000  continue
      if(nlink.le.6) then
        goto (1,2,3,4,5,6),nlink
1         call mxml1_(a,mcola,mrowa,b,mcolb,mrowb,
     >               r,mcolr,mrowr,ncol,nrow)
          return
2         call mxml2_(a,mcola,mrowa,b,mcolb,mrowb,
     >               r,mcolr,mrowr,ncol,nrow)
          return
3         call mxml3_(a,mcola,mrowa,b,mcolb,mrowb,
     >               r,mcolr,mrowr,ncol,nrow)
          return
4         call mxml4_(a,mcola,mrowa,b,mcolb,mrowb,
     >               r,mcolr,mrowr,ncol,nrow)
          return
5         call mxml5_(a,mcola,mrowa,b,mcolb,mrowb,
     >               r,mcolr,mrowr,ncol,nrow)
          return
6         call mxml6_(a,mcola,mrowa,b,mcolb,mrowb,
     >               r,mcolr,mrowr,ncol,nrow)
          return
      end if
#ifdef MOLPRO_BLAS
      if(noblas.ne.0) goto 6010
      if(mcolb.eq.1.and.mrowb.eq.1) goto 6010
      if(mcolr.eq.1) then
        if(nlink.ge.mindgm.and.nrow.ge.mindgm.and.ncol.ge.mindgm)
     1     goto 6001
        if(nlink.ge.mindgl.and.ncol.ge.mindgm2.and.nrow.ge.mindgm2)
     1     goto 6001
        if(nrow.ge.mindgr.and.ncol.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6001
        if(ncol.ge.mindgc.and.nrow.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6001
        goto 6010
6001    mr=mrowr
        if(nrow.eq.1) mr=max(ncol,mr)
        if(mr.lt.ncol) call error('mr.lt.ncol','mxma')
        call zeromat(r,mr,ncol,nrow)
        if(mcola.eq.1.and.mcolb.eq.1.and.mrowa.ge.ncol) then
          call dgemm_x('N','N',ncol,nrow,nlink,1.0d0,a,mrowa,
     1                    b,max(nlink,mrowb),0.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mcolb.eq.1) then
          call dgemm_x('T','N',ncol,nrow,nlink,1.0d0,a,mcola,
     1                    b,max(nlink,mrowb),0.0d0,r,mr)
          return
        else if(mcola.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','T',ncol,nrow,nlink,1.0d0,a,mrowa,
     1                    b,mcolb,0.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mrowb.eq.1) then
          call dgemm_x('T','T',ncol,nrow,nlink,1.0d0,a,mcola,
     1                    b,mcolb,0.0d0,r,mr)
          return
        else
          call mxma_dgm(a,mcola,mrowa,b,mcolb,mrowb,
     1                  r,mcolr,mrowr,ncol,nlink,nrow)
          return
        end if
      else if(mrowr.eq.1) then
        if(nlink.ge.mindgm.and.nrow.ge.mindgm.and.ncol.ge.mindgm)
     1     goto 6002
        if(nlink.ge.mindgl.and.ncol.ge.mindgm2.and.nrow.ge.mindgm2)
     1     goto 6002
        if(nrow.ge.mindgc.and.ncol.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6002
        if(ncol.ge.mindgr.and.nrow.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6002
        goto 6010
6002    mr=mcolr
        if(ncol.eq.1) mr=max(nrow,mr)
        if(mr.lt.nrow) call error('mr.lt.nrow','mxma')
        call zeromat(r,mr,nrow,ncol)
        if(mcola.eq.1.and.mcolb.eq.1) then
          call dgemm_x('T','T',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb),
     1                     a,mrowa,0.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mcolb.eq.1) then
          call dgemm_x('T','N',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb),
     1                     a,mcola,0.0d0,r,mr)
          return
        else if(mcola.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','T',nrow,ncol,nlink,1.0d0,b,mcolb,
     1                     a,mrowa,0.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','N',nrow,ncol,nlink,1.0d0,b,mcolb,
     1                     a,mcola,0.0d0,r,mr)
          return
        end if
      end if
      if(nlink.lt.mindgm.or.nrow.lt.mindgm.or.ncol.lt.mindgm) goto 6010
      call mxma_dgm(a,mcola,mrowa,b,mcolb,mrowb,
     1              r,mcolr,mrowr,ncol,nlink,nrow)
      return
6010  continue
#endif
#ifdef MOLPRO_SunOS
#else
      i0=0
      if(nroll.le.1) then
           call mxma2_(a,mcola,mrowa,b,mcolb,mrowb,
     1                 r,mcolr,mrowr,ncol,nlink,nrow,i0)
           return
      else if(nroll.eq.2) then
        if(mcola.eq.1) then
          if(mcolb.eq.1) then
             call mxma2_nn(a,mrowa,b,mrowb,
     1                     r,mcolr,mrowr,ncol,nlink,nrow,i0)
             return
          else if(mrowb.eq.1) then
             call mxma2_nt(a,mrowa,b,mcolb,
     1                     r,mcolr,mrowr,ncol,nlink,nrow,i0)
             return
          else
             call mxma2_(a,mcola,mrowa,b,mcolb,mrowb,
     1                   r,mcolr,mrowr,ncol,nlink,nrow,i0)
             return
          end if
        else if(mrowa.eq.1.and.mcolb.eq.1) then
             call mxma2_tn(a,mcola,b,mrowb,
     1                     r,mcolr,mrowr,ncol,nlink,nrow,i0)
             return
        else
             call mxma2_(a,mcola,mrowa,b,mcolb,mrowb,
     1                   r,mcolr,mrowr,ncol,nlink,nrow,i0)
             return
        end if
      end if
      if(nlink*(ncol+nrow)+nrow*ncol.le.ncache) then
        if(nroll.eq.3) then
          call mxma3_(a,mcola,mrowa,b,mcolb,mrowb,
     1                r,mcolr,mrowr,ncol,nlink,nrow)
          return
        else
          call mxma4_(a,mcola,mrowa,b,mcolb,mrowb,
     1                r,mcolr,mrowr,ncol,nlink,nrow)
          return
        end if
      end if
#endif
      i1=1
      if(mrowa.ne.1) goto 6200
      mxa=mxmblk
      mxb=min(maxb,mxmblk)
      nkb=min(maxb,mxmbln)
      ke=0
      do 60 k=1,nlink,nkb
      nr=nlink-ke
      if(nr.eq.0) goto 60
      nk=min(nkb,nr)
      if(nr-nk.gt.0.and.nr-nk.lt.nkb) nk=min(maxb,(nr-1)/2+1)
      je=0
      do 70 j=1,nrow,mxb
      nr=nrow-je
      if(nr.eq.0) goto 70
      nj=min(mxb,nr)
      if(nr-nj.gt.0.and.nr-nj.lt.mxb) nj=min(maxb,(nr-1)/2+1)
      ie=0
      ib1=je*mrowb+ke*mcolb+1
      call copmx(b(ib1),mcolb,mrowb,buf,i1,nk, nk,nj)
      do 80 i=1,ncol,mxa
      nr=ncol-ie
      if(nr.eq.0) goto 80
      ni=min(mxa,nr)
      if(nr-ni.gt.0.and.nr-ni.lt.mxa) ni=(nr-1)/2+1
      ia1=ke*mrowa+ie*mcola+1
      ir1=je*mrowr+ie*mcolr+1
      if(nroll.le.3) then
        if(k.eq.1) then
           call mxma3_(a(ia1),mcola,mrowa,buf,i1,nk,r(ir1),
     1                 mcolr, mrowr,ni,nk,nj)
         else
           call mxmb3_(a(ia1),mcola,mrowa,buf,i1,nk,r(ir1),
     1                 mcolr,mrowr, ni,nk,nj)
         end if
      else
        if(k.eq.1) then
          call mxma4_(a(ia1),mcola,mrowa,buf,i1,nk,r(ir1),
     1                mcolr, mrowr,ni,nk,nj)
        else
          call mxmb4_(a(ia1),mcola,mrowa,buf,i1,nk,r(ir1),
     1                mcolr,mrowr, ni,nk,nj)
        end if
      end if
80    ie=ie+ni
70    je=je+nj
60    ke=ke+nk
      return
6200  ke=0
      mxa=min(maxb,mxmblk)
      nkb=min(maxb,mxmbln)
      mxb=mxmblk
      do 61 k=1,nlink,nkb
      nr=nlink-ke
      if(nr.eq.0) goto 61
      nk=min(nkb,nr)
      if(nr-nk.gt.0.and.nr-nk.lt.nkb) nk=min(maxb,(nr-1)/2+1)
      ie=0
      do 81 i=1,ncol,mxa
      nr=ncol-ie
      if(nr.eq.0) goto 81
      ni=min(mxa,nr)
      if(nr-ni.gt.0.and.nr-ni.lt.mxa) ni=min(maxb,(nr-1)/2+1)
      ia1=ke*mrowa+ie*mcola+1
      call copmx(a(ia1),mcola,mrowa, buf,nk,i1, ni,nk)
      je=0
      do 71 j=1,nrow,mxb
      nr=nrow-je
      if(nr.eq.0) goto 71
      nj=min(mxb,nr)
      if(nr-nj.gt.0.and.nr-nj.lt.mxb) nj=(nr-1)/2+1
      ib1=je*mrowb+ke*mcolb+1
      ir1=je*mrowr+ie*mcolr+1
      if(nroll.le.3) then
        if(k.eq.1) then
          call mxma3_(buf,nk,i1,b(ib1),mcolb,mrowb,r(ir1),
     1                mcolr, mrowr,ni,nk,nj)
        else
          call mxmb3_(buf,nk,i1,b(ib1),mcolb,mrowb,r(ir1),
     1                mcolr,mrowr, ni,nk,nj)
        end if
      else
        if(k.eq.1) then
          call mxma4_(buf,nk,i1,b(ib1),mcolb,mrowb,r(ir1),
     1                mcolr, mrowr,ni,nk,nj)
        else
          call mxmb4_(buf,nk,i1,b(ib1),mcolb,mrowb,r(ir1),
     1                mcolr,mrowr, ni,nk,nj)
        end if
      end if
71    je=je+nj
81    ie=ie+ni
61    ke=ke+nk
      return
      end

!> matrix-matrix product
!>
!> calculates the \f${\tt ncol}\times{\tt nrow}\f$ matrix product
!> \f$\bf R = R + A B\f$.
!> Parameters are exactly as for \c mxma; the only difference is that
!>       \c r is not cleared to zero at the start, but the matrix
!>       product is added to the previous contents of \c r.
c--------------------------------------------------------------------
      subroutine mxmb(a,mcola_,mrowa_,b,mcolb_,mrowb_,
     1                r,mcolr_,mrowr_,ncol_,nlink_,nrow_)
c--------------------------------------------------------------------
c.....unrolled version for ibm6000 and other risc machines
      implicit double precision (a-h,o-z)
      integer mcola_,mrowa_,mcolb_,mrowb_,mcolr_,mrowr_,
     >        ncol_,nlink_,nrow_
    !  include "common/clseg"
      parameter (maxb=128)
      dimension buf(maxb*maxb)
      dimension r(*),a(*),b(*)
      if(ncol_.eq.0.or.nrow_.eq.0.or.nlink_.eq.0) return
cmgs start
      nrow=nrow_
      ncol=ncol_
      nlink=nlink_
      mcola=mcola_
      mrowa=mrowa_
      mcolb=mcolb_
      mrowb=mrowb_
      mcolr=mcolr_
      mrowr=mrowr_
      if(nlink.eq.1) then
        jj=1
        ijj=1
        do 6210 j=1,nrow
        ii=1
        ij=ijj
        if(b(jj).ne.0.0d0) then
        do 6110 i=1,ncol
        r(ij)=r(ij)+a(ii)*b(jj)
        ii=ii+mcola
 6110   ij=ij+mcolr
        end if
        ijj=ijj+mrowr
 6210   jj=jj+mrowb
        return
      end if
cmgs end
      if(ncol.gt.5.or.nrow.gt.5) goto 6000
      goto (1000,2000,3000,4000,5000),nrow
c
c....nrow=1
1000  goto (1001,1002,1003,1004,1005),ncol
1001    ia1=1
        ib1=1
        s1=0.0d0
        do 1010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
1010    ib1=ib1+mcolb
        r(1)=r(1)+s1
        return
c
1002    ia1=1
        ib1=1
        ia2=1+mcola
        s1=0.0d0
        s2=0.0d0
        do 1020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
1020    ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        return
c
1003    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        do 1030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
1030    ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        return
c
1004    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        do 1040 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
1040    ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        r(1+3*mcolr)=r(1+3*mcolr)+s4
        return
c
1005    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        do 1050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
        s5=a(ia5)*b(ib1)+s5
        ia5=ia5+mrowa
1050    ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        r(1+3*mcolr)=r(1+3*mcolr)+s4
        r(1+4*mcolr)=r(1+4*mcolr)+s5
        return
c
c....nrow=2
2000  goto (2001,2002,2003,2004,2005),ncol
2001    ia1=1
        ib1=1
        ib2=1+mrowb
        s1=0.0d0
        s2=0.0d0
        do 2010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
2010    ia1=ia1+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        return
c
2002    ia1=1
        ib1=1
        ia2=1+mcola
        ib2=1+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        do 2020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia2)*b(ib1)+s2
        ib1=ib1+mcolb
        s3=a(ia1)*b(ib2)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib2)+s4
        ib2=ib2+mcolb
2020    ia2=ia2+mrowa
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+mrowr)=r(1+mrowr)+s3
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s4
        return
c
2003    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        do 2030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s4=a(ia1)*b(ib2)+s4
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        s5=a(ia2)*b(ib2)+s5
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ib1=ib1+mcolb
        s6=a(ia3)*b(ib2)+s6
        ia3=ia3+mrowa
2030    ib2=ib2+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        r(1+mrowr)=r(1+mrowr)+s4
        r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s5
        r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s6
        return
c
2004    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        do 2040 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s42=a(ia4)*b(ib2)+s42
        ia4=ia4+mrowa
2040    ib2=ib2+mcolb
        r(1)=r(1)+s11
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+3*mcolr)=r(1+3*mcolr)+s41
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s22
        r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s32
        r(1+mrowr+3*mcolr)=r(1+mrowr+3*mcolr)+s42
        return
c
2005    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s51=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s52=0.0d0
        do 2050 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        s42=a(ia4)*b(ib2)+s42
        ia4=ia4+mrowa
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s52=a(ia5)*b(ib2)+s52
        ia5=ia5+mrowa
2050    ib2=ib2+mcolb
        r(1)=r(1)+s11
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+3*mcolr)=r(1+3*mcolr)+s41
        r(1+4*mcolr)=r(1+4*mcolr)+s51
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s22
        r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s32
        r(1+mrowr+3*mcolr)=r(1+mrowr+3*mcolr)+s42
        r(1+mrowr+4*mcolr)=r(1+mrowr+4*mcolr)+s52
        return
c
c....nrow=3
3000  goto (3001,3002,3003,3004,3005),ncol
3001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        do 3010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
3010    ia1=ia1+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        return
c
3002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        do 3020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s4=a(ia2)*b(ib1)+s4
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        s5=a(ia2)*b(ib2)+s5
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib3)+s6
        ib3=ib3+mcolb
3020    ia2=ia2+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+mcolr)=r(1+mcolr)+s4
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s5
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s6
        return
c
3003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        s7=0.0d0
        s8=0.0d0
        s9=0.0d0
        do 3030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib1)+s4
        s5=a(ia2)*b(ib2)+s5
        s6=a(ia2)*b(ib3)+s6
        ia2=ia2+mrowa
        s7=a(ia3)*b(ib1)+s7
        ib1=ib1+mcolb
        s8=a(ia3)*b(ib2)+s8
        ib2=ib2+mcolb
        s9=a(ia3)*b(ib3)+s9
        ia3=ia3+mrowa
3030    ib3=ib3+mcolb
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+mcolr)=r(1+mcolr)+s4
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s5
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s6
        r(1+2*mcolr)=r(1+2*mcolr)+s7
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s8
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s9
        return
c
3004    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        s7=0.0d0
        s8=0.0d0
        s9=0.0d0
        s10=0.0d0
        s11=0.0d0
        s12=0.0d0
        do 3040 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib1)+s4
        s5=a(ia2)*b(ib2)+s5
        s6=a(ia2)*b(ib3)+s6
        ia2=ia2+mrowa
        s7=a(ia3)*b(ib1)+s7
        s8=a(ia3)*b(ib2)+s8
        s9=a(ia3)*b(ib3)+s9
        ia3=ia3+mrowa
        s10=a(ia4)*b(ib1)+s10
        ib1=ib1+mcolb
        s11=a(ia4)*b(ib2)+s11
        ib2=ib2+mcolb
        s12=a(ia4)*b(ib3)+s12
        ia4=ia4+mrowa
3040    ib3=ib3+mcolb
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+mcolr)=r(1+mcolr)+s4
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s5
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s6
        r(1+2*mcolr)=r(1+2*mcolr)+s7
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s8
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s9
        r(1+3*mcolr)=r(1+3*mcolr)+s10
        r(1+3*mcolr+mrowr)=r(1+3*mcolr+mrowr)+s11
        r(1+3*mcolr+2*mrowr)=r(1+3*mcolr+2*mrowr)+s12
        return
c
3005    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s51=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s52=0.0d0
        s13=0.0d0
        s23=0.0d0
        s33=0.0d0
        s43=0.0d0
        s53=0.0d0
        do 3050 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        s13=a(ia1)*b(ib3)+s13
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        s23=a(ia2)*b(ib3)+s23
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        s33=a(ia3)*b(ib3)+s33
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        s42=a(ia4)*b(ib2)+s42
        s43=a(ia4)*b(ib3)+s43
        ia4=ia4+mrowa
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        s53=a(ia5)*b(ib3)+s53
        ia5=ia5+mrowa
3050    ib3=ib3+mcolb
        r(1)=r(1)+s11
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+3*mcolr)=r(1+3*mcolr)+s41
        r(1+4*mcolr)=r(1+4*mcolr)+s51
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s22
        r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s32
        r(1+mrowr+3*mcolr)=r(1+mrowr+3*mcolr)+s42
        r(1+mrowr+4*mcolr)=r(1+mrowr+4*mcolr)+s52
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+2*mrowr+mcolr)=r(1+2*mrowr+mcolr)+s23
        r(1+2*mrowr+2*mcolr)=r(1+2*mrowr+2*mcolr)+s33
        r(1+2*mrowr+3*mcolr)=r(1+2*mrowr+3*mcolr)+s43
        r(1+2*mrowr+4*mcolr)=r(1+2*mrowr+4*mcolr)+s53
        return
c
c.....nrow=4
4000  goto(4001,4002,4003,4004,4005),ncol
4001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        do 4010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
        s4=a(ia1)*b(ib4)+s4
        ib4=ib4+mcolb
4010    ia1=ia1+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+3*mrowr)=r(1+3*mrowr)+s4
        return
c
4002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        do 4020 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ib4=ib4+mcolb
4020    ia2=ia2+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        return
c
4003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s31=0.0d0
        s32=0.0d0
        s33=0.0d0
        s34=0.0d0
        do 4030 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ia2=ia2+mrowa
        s34=a(ia3)*b(ib4)+s34
        ib4=ib4+mcolb
4030    ia3=ia3+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s32
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s33
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s34
        return
c
4004    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s31=0.0d0
        s32=0.0d0
        s33=0.0d0
        s34=0.0d0
        s41=0.0d0
        s42=0.0d0
        s43=0.0d0
        s44=0.0d0
        do 4040 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ia2=ia2+mrowa
        s34=a(ia3)*b(ib4)+s34
        ia3=ia3+mrowa
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
4040    ia4=ia4+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s32
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s33
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s34
        r(1+3*mcolr)=r(1+3*mcolr)+s41
        r(1+3*mcolr+mrowr)=r(1+3*mcolr+mrowr)+s42
        r(1+3*mcolr+2*mrowr)=r(1+3*mcolr+2*mrowr)+s43
        r(1+3*mcolr+3*mrowr)=r(1+3*mcolr+3*mrowr)+s44
        return
c
4005    ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s6=0.0d0
        s7=0.0d0
        s8=0.0d0
        s9=0.0d0
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s16=0.0d0
        s17=0.0d0
        s18=0.0d0
        s19=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        do 4050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        s4=a(ia1)*b(ib4)+s4
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib1)+s6
        s7=a(ia2)*b(ib2)+s7
        s8=a(ia2)*b(ib3)+s8
        s9=a(ia2)*b(ib4)+s9
        ia2=ia2+mrowa
        s11=a(ia3)*b(ib1)+s11
        s12=a(ia3)*b(ib2)+s12
        s13=a(ia3)*b(ib3)+s13
        s14=a(ia3)*b(ib4)+s14
        ia3=ia3+mrowa
        s16=a(ia4)*b(ib1)+s16
        s17=a(ia4)*b(ib2)+s17
        s18=a(ia4)*b(ib3)+s18
        s19=a(ia4)*b(ib4)+s19
        ia4=ia4+mrowa
        s21=a(ia5)*b(ib1)+s21
        s22=a(ia5)*b(ib2)+s22
        s23=a(ia5)*b(ib3)+s23
        s24=a(ia5)*b(ib4)+s24
        ia5=ia5+mrowa
        ib1=ib1+mcolb
        ib2=ib2+mcolb
        ib3=ib3+mcolb
        ib4=ib4+mcolb
4050    continue
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+3*mrowr)=r(1+3*mrowr)+s4
        r(1+mcolr)=r(1+mcolr)+s6
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s7
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s8
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s9
        r(1+2*mcolr)=r(1+2*mcolr)+s11
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s12
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s13
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s14
        r(1+3*mcolr)=r(1+3*mcolr)+s16
        r(1+3*mcolr+mrowr)=r(1+3*mcolr+mrowr)+s17
        r(1+3*mcolr+2*mrowr)=r(1+3*mcolr+2*mrowr)+s18
        r(1+3*mcolr+3*mrowr)=r(1+3*mcolr+3*mrowr)+s19
        r(1+4*mcolr)=r(1+4*mcolr)+s21
        r(1+4*mcolr+mrowr)=r(1+4*mcolr+mrowr)+s22
        r(1+4*mcolr+2*mrowr)=r(1+4*mcolr+2*mrowr)+s23
        r(1+4*mcolr+3*mrowr)=r(1+4*mcolr+3*mrowr)+s24
        return
c
c.....nrow=5
5000  goto(5001,5002,5003,5004,5004),ncol
5001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        do 5010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
        s4=a(ia1)*b(ib4)+s4
        ib4=ib4+mcolb
        s5=a(ia1)*b(ib5)+s5
        ib5=ib5+mcolb
5010    ia1=ia1+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+3*mrowr)=r(1+3*mrowr)+s4
        r(1+4*mrowr)=r(1+4*mrowr)+s5
        return
c
5002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=1
        ia2=ia1+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s15=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s25=0.0d0
        do 5020 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        ia1=ia1+mrowa
        s25=a(ia2)*b(ib5)+s25
        ib5=ib5+mcolb
5020    ia2=ia2+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+4*mrowr)=r(1+4*mrowr)+s15
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        r(1+mcolr+4*mrowr)=r(1+mcolr+4*mrowr)+s25
        return
c
5003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s15=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s25=0.0d0
        s31=0.0d0
        s32=0.0d0
        s33=0.0d0
        s34=0.0d0
        s35=0.0d0
        do 5030 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        ia1=ia1+mrowa
        s25=a(ia2)*b(ib5)+s25
        ia2=ia2+mrowa
        s35=a(ia3)*b(ib5)+s35
        ib5=ib5+mcolb
5030    ia3=ia3+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+4*mrowr)=r(1+4*mrowr)+s15
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        r(1+mcolr+4*mrowr)=r(1+mcolr+4*mrowr)+s25
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s32
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s33
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s34
        r(1+2*mcolr+4*mrowr)=r(1+2*mcolr+4*mrowr)+s35
        return
c
5004    ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        s4=0.0d0
        s5=0.0d0
        s6=0.0d0
        s7=0.0d0
        s8=0.0d0
        s9=0.0d0
        s10=0.0d0
        s11=0.0d0
        s12=0.0d0
        s13=0.0d0
        s14=0.0d0
        s15=0.0d0
        s16=0.0d0
        s17=0.0d0
        s18=0.0d0
        s19=0.0d0
        s20=0.0d0
        s21=0.0d0
        s22=0.0d0
        s23=0.0d0
        s24=0.0d0
        s25=0.0d0
        do 5050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        s4=a(ia1)*b(ib4)+s4
        s5=a(ia1)*b(ib5)+s5
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib1)+s6
        s7=a(ia2)*b(ib2)+s7
        s8=a(ia2)*b(ib3)+s8
        s9=a(ia2)*b(ib4)+s9
        s10=a(ia2)*b(ib5)+s10
        ia2=ia2+mrowa
        s11=a(ia3)*b(ib1)+s11
        s12=a(ia3)*b(ib2)+s12
        s13=a(ia3)*b(ib3)+s13
        s14=a(ia3)*b(ib4)+s14
        s15=a(ia3)*b(ib5)+s15
        ia3=ia3+mrowa
        s16=a(ia4)*b(ib1)+s16
        s17=a(ia4)*b(ib2)+s17
        s18=a(ia4)*b(ib3)+s18
        s19=a(ia4)*b(ib4)+s19
        s20=a(ia4)*b(ib5)+s20
        ia4=ia4+mrowa
        if(ncol.eq.4) goto 5040
        s21=a(ia5)*b(ib1)+s21
        s22=a(ia5)*b(ib2)+s22
        s23=a(ia5)*b(ib3)+s23
        s24=a(ia5)*b(ib4)+s24
        s25=a(ia5)*b(ib5)+s25
        ia5=ia5+mrowa
5040    ib1=ib1+mcolb
        ib2=ib2+mcolb
        ib3=ib3+mcolb
        ib4=ib4+mcolb
5050    ib5=ib5+mcolb
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+3*mrowr)=r(1+3*mrowr)+s4
        r(1+4*mrowr)=r(1+4*mrowr)+s5
        r(1+mcolr)=r(1+mcolr)+s6
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s7
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s8
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s9
        r(1+mcolr+4*mrowr)=r(1+mcolr+4*mrowr)+s10
        r(1+2*mcolr)=r(1+2*mcolr)+s11
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s12
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s13
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s14
        r(1+2*mcolr+4*mrowr)=r(1+2*mcolr+4*mrowr)+s15
        r(1+3*mcolr)=r(1+3*mcolr)+s16
        r(1+3*mcolr+mrowr)=r(1+3*mcolr+mrowr)+s17
        r(1+3*mcolr+2*mrowr)=r(1+3*mcolr+2*mrowr)+s18
        r(1+3*mcolr+3*mrowr)=r(1+3*mcolr+3*mrowr)+s19
        r(1+3*mcolr+4*mrowr)=r(1+3*mcolr+4*mrowr)+s20
        if(ncol.eq.4) return
        r(1+4*mcolr)=r(1+4*mcolr)+s21
        r(1+4*mcolr+mrowr)=r(1+4*mcolr+mrowr)+s22
        r(1+4*mcolr+2*mrowr)=r(1+4*mcolr+2*mrowr)+s23
        r(1+4*mcolr+3*mrowr)=r(1+4*mcolr+3*mrowr)+s24
        r(1+4*mcolr+4*mrowr)=r(1+4*mcolr+4*mrowr)+s25
        return
c
6000  continue
      if(nlink.le.6) then
        goto (1,2,3,4,5,6),nlink
1         call mxml1b_(a,mcola,mrowa,b,mcolb,mrowb,
     >                r,mcolr,mrowr,ncol,nrow)
          return
2         call mxml2b_(a,mcola,mrowa,b,mcolb,mrowb,
     >                r,mcolr,mrowr,ncol,nrow)
          return
3         call mxml3b_(a,mcola,mrowa,b,mcolb,mrowb,
     >                r,mcolr,mrowr,ncol,nrow)
          return
4         call mxml4b_(a,mcola,mrowa,b,mcolb,mrowb,
     >                r,mcolr,mrowr,ncol,nrow)
          return
5         call mxml5b_(a,mcola,mrowa,b,mcolb,mrowb,
     >                r,mcolr,mrowr,ncol,nrow)
          return
6         call mxml6b_(a,mcola,mrowa,b,mcolb,mrowb,
     >                r,mcolr,mrowr,ncol,nrow)
          return
      end if
#ifdef MOLPRO_BLAS
      if(noblas.ne.0) goto 6010
      if(mcolb.eq.1.and.mrowb.eq.1) goto 6010
      if(mcolr.eq.1) then
        if(nlink.ge.mindgm.and.nrow.ge.mindgm.and.ncol.ge.mindgm)
     1     goto 6001
        if(nlink.ge.mindgl.and.ncol.ge.mindgm2.and.nrow.ge.mindgm2)
     1     goto 6001
        if(nrow.ge.mindgr.and.ncol.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6001
        if(ncol.ge.mindgc.and.nrow.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6001
        goto 6010
6001    mr=mrowr
        if(nrow.eq.1) mr=max(ncol,mr)
        if(mr.lt.ncol) call error('mr.lt.ncol','mxmb')
        if(mcola.eq.1.and.mcolb.eq.1.and.mrowa.ge.ncol) then
          call dgemm_x('N','N',ncol,nrow,nlink,1.0d0,a,mrowa,
     1                     b,max(nlink,mrowb),1.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mcolb.eq.1) then
          call dgemm_x('T','N',ncol,nrow,nlink,1.0d0,a,mcola,
     1                     b,max(nlink,mrowb),1.0d0,r,mr)
          return
        else if(mcola.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','T',ncol,nrow,nlink,1.0d0,a,mrowa,
     1                     b,mcolb,1.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mrowb.eq.1) then
          call dgemm_x('T','T',ncol,nrow,nlink,1.0d0,a,mcola,
     1                     b,mcolb,1.0d0,r,mr)
          return
        else
          call mxmb_dgm(a,mcola,mrowa,b,mcolb,mrowb,
     1                  r,mcolr,mrowr,ncol,nlink,nrow)
          return
        end if
      else if(mrowr.eq.1) then
        if(nlink.ge.mindgm.and.nrow.ge.mindgm.and.ncol.ge.mindgm)
     1     goto 6002
        if(nlink.ge.mindgl.and.ncol.ge.mindgm2.and.nrow.ge.mindgm2)
     1     goto 6002
        if(nrow.ge.mindgc.and.ncol.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6002
        if(ncol.ge.mindgr.and.nrow.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6002
        goto 6010
6002    mr=mcolr
        if(ncol.eq.1) mr=max(nrow,mr)
        if(mr.lt.nrow) call error('mr.lt.nrow','mxmb')
        if(mcola.eq.1.and.mcolb.eq.1) then
          call dgemm_x('T','T',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb),
     1                     a,mrowa,1.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mcolb.eq.1) then
          call dgemm_x('T','N',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb),
     1                     a,mcola,1.0d0,r,mr)
          return
        else if(mcola.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','T',nrow,ncol,nlink,1.0d0,b,mcolb,
     1                     a,mrowa,1.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','N',nrow,ncol,nlink,1.0d0,b,mcolb,
     1                     a,mcola,1.0d0,r,mr)
          return
        end if
      end if
      if(nlink.lt.mindgm.or.nrow.lt.mindgm.or.ncol.lt.mindgm) goto 6010
      call mxmb_dgm(a,mcola,mrowa,b,mcolb,mrowb,
     1              r,mcolr,mrowr,ncol,nlink,nrow)
      return
6010  continue
#endif
      i1=1
      if(nroll.le.1) then
         call mxma2_(a,mcola,mrowa,b,mcolb,mrowb,
     1               r,mcolr,mrowr,ncol,nlink,nrow,i1)
         return
      else if(nroll.eq.2) then
        if(mcola.eq.1) then
          if(mcolb.eq.1) then
             call mxma2_nn(a,mrowa,b,mrowb,
     1                     r,mcolr,mrowr,ncol,nlink,nrow,i1)
             return
          else if(mrowb.eq.1) then
             call mxma2_nt(a,mrowa,b,mcolb,
     1                     r,mcolr,mrowr,ncol,nlink,nrow,i1)
             return
          else
             call mxma2_(a,mcola,mrowa,b,mcolb,mrowb,
     1                   r,mcolr,mrowr,ncol,nlink,nrow,i1)
             return
          end if
        else if(mrowa.eq.1.and.mcolb.eq.1) then
             call mxma2_tn(a,mcola,b,mrowb,
     1                   r,mcolr,mrowr,ncol,nlink,nrow,i1)
             return
        else
             call mxma2_(a,mcola,mrowa,b,mcolb,mrowb,
     1                   r,mcolr,mrowr,ncol,nlink,nrow,i1)
             return
        end if
      end if
#ifdef MOLPRO_SunOS
#else
      if(nlink*(ncol+nrow)+nrow*ncol.le.ncache) then
        if(nroll.eq.3) then
          call mxmb3_(a,mcola,mrowa,b,mcolb,mrowb,
     1                r,mcolr,mrowr,ncol,nlink,nrow)
          return
        else
          call mxmb4_(a,mcola,mrowa,b,mcolb,mrowb,
     1                 r,mcolr,mrowr,ncol,nlink,nrow)
          return
        end if
      end if
#endif
      mxa=mxmblk
      mxb=mxmblk
      nkb=mxmbln
#ifdef MOLPRO_SunOS
      if(mrowa.ne.1.and.mcolb.ne.1) then
        mxa=mxmblk/2
        mxb=mxmblk/2
        nkb=mxmbln/2
      else if(mrowa.ne.1) then
        mxa=mxmblk/2
      else if(mcolb.ne.1) then
        mxb=mxmblk/2
      else
        mxa=mxmblk/2
        mxb=mxmblk/2
        nkb=min(512,nlink)
      end if
#else
      if(mrowa.eq.1.and.mcolb.eq.1) nkb=nkb*2
      if(mrowa.ne.1.and.mcolb.ne.1) nkb=nkb/2
#endif
      ke=0
      i1=1
      if(mrowa.ne.1) goto 6200
      mxa=mxmblk
      mxb=min(maxb,mxmblk)
      nkb=min(maxb,mxmbln)
      do 60 k=1,nlink,nkb
      nr=nlink-ke
      if(nr.eq.0) goto 60
      nk=min(nkb,nr)
      if(nr-nk.gt.0.and.nr-nk.lt.nkb) nk=min(maxb,(nr-1)/2+1)
      je=0
      do 70 j=1,nrow,mxb
      nr=nrow-je
      if(nr.eq.0) goto 70
      nj=min(mxb,nr)
      if(nr-nj.gt.0.and.nr-nj.lt.mxb) nj=min(maxb,(nr-1)/2+1)
      ib1=je*mrowb+ke*mcolb+1
      call copmx(b(ib1),mcolb,mrowb,buf,i1,nk,nk,nj)
      ie=0
      do 80 i=1,ncol,mxa
      nr=ncol-ie
      if(nr.eq.0) goto 80
      ni=min(mxa,nr)
      if(nr-ni.gt.0.and.nr-ni.lt.mxa) ni=(nr-1)/2+1
      ia1=ke*mrowa+ie*mcola+1
      ir1=je*mrowr+ie*mcolr+1
      if(nroll.le.3) then
        call mxmb3_(a(ia1),mcola,mrowa,buf,i1,nk,r(ir1),
     1              mcolr,mrowr, ni,nk,nj)
      else
        call mxmb4_(a(ia1),mcola,mrowa,buf,i1,nk,r(ir1),
     1              mcolr,mrowr, ni,nk,nj)
      end if
80    ie=ie+ni
70    je=je+nj
60    ke=ke+nk
      return
6200  ke=0
      mxa=min(maxb,mxmblk)
      nkb=min(maxb,mxmbln)
      mxb=mxmblk
      do 61 k=1,nlink,nkb
      nr=nlink-ke
      if(nr.eq.0) goto 61
      nk=min(nkb,nr)
      if(nr-nk.gt.0.and.nr-nk.lt.nkb) nk=min(maxb,(nr-1)/2+1)
      ie=0
      do 81 i=1,ncol,mxa
      nr=ncol-ie
      if(nr.eq.0) goto 81
      ni=min(mxa,nr)
      if(nr-ni.gt.0.and.nr-ni.lt.mxa) ni=min(maxb,(nr-1)/2+1)
      ia1=ke*mrowa+ie*mcola+1
      call copmx(a(ia1),mcola,mrowa, buf,nk,i1, ni,nk)
      je=0
      do 71 j=1,nrow,mxb
      nr=nrow-je
      if(nr.eq.0) goto 71
      nj=min(mxb,nr)
      if(nr-nj.gt.0.and.nr-nj.lt.mxb) nj=(nr-1)/2+1
      ib1=je*mrowb+ke*mcolb+1
      ir1=je*mrowr+ie*mcolr+1
      if(nroll.le.3) then
        call mxmb3_(buf,nk,i1,b(ib1),mcolb,mrowb,r(ir1),
     1              mcolr,mrowr, ni,nk,nj)
      else
        call mxmb4_(buf,nk,i1,b(ib1),mcolb,mrowb,r(ir1),
     1              mcolr,mrowr, ni,nk,nj)
      end if
71    je=je+nj
81    ie=ie+ni
61    ke=ke+nk
      return
      end
c--------------------------------------------------------------------
      subroutine mxman(a,mcola_,mrowa_,b,mcolb_,mrowb_,
     1                r,mcolr_,mrowr_,ncol_,nlink_,nrow_)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
    !  include "common/clseg"
      integer mcola_,mrowa_,mcolb_,mrowb_,mcolr_,mrowr_,
     >        ncol_,nlink_,nrow_
      parameter (i2=2,i4=4)
      dimension r(1),a(1),b(1)
c... r(ncol,nrow)=-a(ncol,nlink)*b(nlink,nrow) matrix mult
      nrow=nrow_
      ncol=ncol_
      nlink=nlink_
      mcola=mcola_
      mrowa=mrowa_
      mcolb=mcolb_
      mrowb=mrowb_
      mcolr=mcolr_
      mrowr=mrowr_
#ifdef MOLPRO_BLAS
      if(noblas.ne.0) goto 6010
      if(mcolb.eq.1.and.mrowb.eq.1) goto 6010
      if(mcolr.eq.1) then
        if(nlink.ge.mindgm.and.nrow.ge.mindgm.and.ncol.ge.mindgm)
     1     goto 6001
        if(nlink.ge.mindgl.and.ncol.ge.mindgm2.and.nrow.ge.mindgm2)
     1     goto 6001
        if(nrow.ge.mindgr.and.ncol.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6001
        if(ncol.ge.mindgc.and.nrow.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6001
        goto 6010
6001    mr=mrowr
        if(nrow.eq.1) mr=max(ncol,mr)
        if(mr.lt.ncol) call error('mr.lt.ncol','mxma')
        call zeromat(r,mr,ncol,nrow)
        if(mcola.eq.1.and.mcolb.eq.1.and.mrowa.ge.ncol) then
          call dgemm_x('N','N',ncol,nrow,nlink,-1.0d0,a,mrowa,
     1                    b,max(nlink,mrowb),0.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mcolb.eq.1) then
          call dgemm_x('T','N',ncol,nrow,nlink,-1.0d0,a,mcola,
     1                    b,max(nlink,mrowb),0.0d0,r,mr)
          return
        else if(mcola.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','T',ncol,nrow,nlink,-1.0d0,a,mrowa,
     1                    b,mcolb,0.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mrowb.eq.1) then
          call dgemm_x('T','T',ncol,nrow,nlink,-1.0d0,a,mcola,
     1                    b,mcolb,0.0d0,r,mr)
          return
        else
          call mxman_dgm(a,mcola,mrowa,b,mcolb,mrowb,
     1                   r,mcolr,mrowr,ncol,nlink,nrow)
          return
        end if
      else if(mrowr.eq.1) then
        if(nlink.ge.mindgm.and.nrow.ge.mindgm.and.ncol.ge.mindgm)
     1     goto 6002
        if(nlink.ge.mindgl.and.ncol.ge.mindgm2.and.nrow.ge.mindgm2)
     1     goto 6002
        if(nrow.ge.mindgc.and.ncol.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6002
        if(ncol.ge.mindgr.and.nrow.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6002
        goto 6010
6002    mr=mcolr
        if(ncol.eq.1) mr=max(nrow,mr)
        if(mr.lt.nrow) call error('mr.lt.nrow','mxma')
        call zeromat(r,mr,nrow,ncol)
        if(mcola.eq.1.and.mcolb.eq.1) then
         call dgemm_x('T','T',nrow,ncol,nlink,-1.0d0,b,max(nlink,mrowb),
     1                     a,mrowa,0.0d0,r,mr)
         return
        else if(mrowa.eq.1.and.mcolb.eq.1) then
         call dgemm_x('T','N',nrow,ncol,nlink,-1.0d0,b,max(nlink,mrowb),
     1                     a,mcola,0.0d0,r,mr)
         return
        else if(mcola.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','T',nrow,ncol,nlink,-1.0d0,b,mcolb,
     1                     a,mrowa,0.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','N',nrow,ncol,nlink,-1.0d0,b,mcolb,
     1                     a,mcola,0.0d0,r,mr)
          return
        end if
      end if
      if(nlink.lt.mindgm.or.nrow.lt.mindgm.or.ncol.lt.mindgm) goto 6010
      call mxman_dgm(a,mcola,mrowa,b,mcolb,mrowb,
     1               r,mcolr,mrowr,ncol,nlink,nrow)
      return
6010  continue
#endif
      iaa=1
      irr=1
      nrest=mod(nrow,i4)
      if(nrest.eq.1.and.nrow.gt.1) nrest=5
      nre=nrow-nrest
      ncest=mod(ncol,i4)
      if(ncest.eq.1.and.ncol.gt.1) ncest=5
      nce=ncol-ncest
      do 300 i=1,nce,4
      ibb=1
      ir1=irr
      ir2=ir1+mcolr
      ir3=ir2+mcolr
      ir4=ir3+mcolr
      irr=ir4+mcolr
      do 200 j=1,nre,4
      ib1=ibb
      ib2=ib1+mrowb
      ib3=ib2+mrowb
      ib4=ib3+mrowb
      ibb=ib4+mrowb
      ia1=iaa
      ia2=ia1+mcola
      ia3=ia2+mcola
      ia4=ia3+mcola
      s11=0.0d0
      s21=0.0d0
      s31=0.0d0
      s41=0.0d0
      s12=0.0d0
      s22=0.0d0
      s32=0.0d0
      s42=0.0d0
      s13=0.0d0
      s23=0.0d0
      s33=0.0d0
      s43=0.0d0
      s14=0.0d0
      s24=0.0d0
      s34=0.0d0
      s44=0.0d0
      do 100 k=1,nlink
      s11=a(ia1)*b(ib1)+s11
      s21=a(ia2)*b(ib1)+s21
      s31=a(ia3)*b(ib1)+s31
      s41=a(ia4)*b(ib1)+s41
      ib1=ib1+mcolb
      s12=a(ia1)*b(ib2)+s12
      s22=a(ia2)*b(ib2)+s22
      s32=a(ia3)*b(ib2)+s32
      s42=a(ia4)*b(ib2)+s42
      ib2=ib2+mcolb
      s13=a(ia1)*b(ib3)+s13
      s23=a(ia2)*b(ib3)+s23
      s33=a(ia3)*b(ib3)+s33
      s43=a(ia4)*b(ib3)+s43
      ib3=ib3+mcolb
      s14=a(ia1)*b(ib4)+s14
      s24=a(ia2)*b(ib4)+s24
      s34=a(ia3)*b(ib4)+s34
      s44=a(ia4)*b(ib4)+s44
      ib4=ib4+mcolb
      ia1=ia1+mrowa
      ia2=ia2+mrowa
      ia3=ia3+mrowa
100   ia4=ia4+mrowa
      r(ir1)=-s11
      r(ir2)=-s21
      r(ir3)=-s31
      r(ir4)=-s41
      r(ir1+mrowr)=-s12
      r(ir2+mrowr)=-s22
      r(ir3+mrowr)=-s32
      r(ir4+mrowr)=-s42
      r(ir1+2*mrowr)=-s13
      r(ir2+2*mrowr)=-s23
      r(ir3+2*mrowr)=-s33
      r(ir4+2*mrowr)=-s43
      r(ir1+3*mrowr)=-s14
      r(ir2+3*mrowr)=-s24
      r(ir3+3*mrowr)=-s34
      r(ir4+3*mrowr)=-s44
      ir1=ir1+4*mrowr
      ir2=ir2+4*mrowr
      ir3=ir3+4*mrowr
      ir4=ir4+4*mrowr
200   continue
      if(nrest.eq.0) goto 300
      goto (201,202,203,300,205),nrest
201     ib1=ibb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        do 101 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
101     ia4=ia4+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        goto 300
c
202     ib1=ibb
        ib2=ib1+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        do 102 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
102     ia4=ia4+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        r(ir4+mrowr)=-s42
        goto 300
c
203     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s13=0.0d0
        s23=0.0d0
        s33=0.0d0
        s43=0.0d0
        do 103 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
103     ia4=ia4+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        r(ir4+mrowr)=-s42
        r(ir1+2*mrowr)=-s13
        r(ir2+2*mrowr)=-s23
        r(ir3+2*mrowr)=-s33
        r(ir4+2*mrowr)=-s43
        goto 300
c
205     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s13=0.0d0
        s23=0.0d0
        s33=0.0d0
        s43=0.0d0
        s14=0.0d0
        s24=0.0d0
        s34=0.0d0
        s44=0.0d0
        s15=0.0d0
        s25=0.0d0
        s35=0.0d0
        s45=0.0d0
        do 105 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        s25=a(ia2)*b(ib5)+s25
        s35=a(ia3)*b(ib5)+s35
        s45=a(ia4)*b(ib5)+s45
        ib5=ib5+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
105     ia4=ia4+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        r(ir4+mrowr)=-s42
        r(ir1+2*mrowr)=-s13
        r(ir2+2*mrowr)=-s23
        r(ir3+2*mrowr)=-s33
        r(ir4+2*mrowr)=-s43
        r(ir1+3*mrowr)=-s14
        r(ir2+3*mrowr)=-s24
        r(ir3+3*mrowr)=-s34
        r(ir4+3*mrowr)=-s44
        r(ir1+4*mrowr)=-s15
        r(ir2+4*mrowr)=-s25
        r(ir3+4*mrowr)=-s35
        r(ir4+4*mrowr)=-s45
300   iaa=iaa+4*mcola
      if(ncest.eq.0) return
      goto (301,302,303,304,305),ncest
301     ibb=1
        ir1=irr
        do 2001 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        s11=0.0d0
        s12=0.0d0
        do 1001 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        ib2=ib2+mcolb
1001    ia1=ia1+mrowa
        r(ir1)=-s11
        r(ir1+mrowr)=-s12
        ir1=ir1+2*mrowr
2001    continue
        if(mod(nrow,i2).eq.0) return
        s11=0.0d0
        do 1011 k=1,nlink
        s11=a(iaa)*b(ibb)+s11
        ibb=ibb+mcolb
1011    iaa=iaa+mrowa
        r(ir1)=-s11
        return
c
302     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        do 2002 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        s11=0.0d0
        s21=0.0d0
        s12=0.0d0
        s22=0.0d0
        do 1002 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        ia1=ia1+mrowa
1002    ia2=ia2+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
2002    continue
        if(mod(nrow,i2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        s11=0.0d0
        s21=0.0d0
        do 1012 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        ia1=ia1+mrowa
        s21=a(ia2)*b(ibb)+s21
        ia2=ia2+mrowa
1012    ibb=ibb+mcolb
        r(ir1)=-s11
        r(ir2)=-s21
        return
c
303     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        do 2003 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        do 1003 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1003    ia3=ia3+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
2003    continue
        if(mod(nrow,i2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        do 1013 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1013    ia3=ia3+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
304     return
305     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        ir4=ir3+mcolr
        ir5=ir4+mcolr
        do 2005 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s51=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s52=0.0d0
        do 1005 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
1005    ia5=ia5+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir5)=-s51
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        r(ir4+mrowr)=-s42
        r(ir5+mrowr)=-s52
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
        ir4=ir4+2*mrowr
        ir5=ir5+2*mrowr
2005    continue
        if(mod(nrow,i2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s51=0.0d0
        do 1015 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        s41=a(ia4)*b(ibb)+s41
        s51=a(ia5)*b(ibb)+s51
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
        ia5=ia5+mrowa
1015    continue
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir5)=-s51
        return
      end
c--------------------------------------------------------------------
      subroutine mxmbn(a,mcola_,mrowa_,b,mcolb_,mrowb_,
     1                r,mcolr_,mrowr_,ncol_,nlink_,nrow_)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer mcola_,mrowa_,mcolb_,mrowb_,mcolr_,mrowr_,
     >        ncol_,nlink_,nrow_
    !  include "common/clseg"
      parameter (i2=2,i4=4)
      dimension r(1),a(1),b(1)
      if(ncol_.eq.0.or.nrow_.eq.0.or.nlink_.eq.0) return
      nrow=nrow_
      ncol=ncol_
      nlink=nlink_
      mcola=mcola_
      mrowa=mrowa_
      mcolb=mcolb_
      mrowb=mrowb_
      mcolr=mcolr_
      mrowr=mrowr_
#ifdef MOLPRO_BLAS
      if(noblas.ne.0) goto 6010
      if(mcolb.eq.1.and.mrowb.eq.1) goto 6010
      if(mcolr.eq.1) then
        if(nlink.ge.mindgm.and.nrow.ge.mindgm.and.ncol.ge.mindgm)
     1     goto 6001
        if(nlink.ge.mindgl.and.ncol.ge.mindgm2.and.nrow.ge.mindgm2)
     1     goto 6001
        if(nrow.ge.mindgr.and.ncol.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6001
        if(ncol.ge.mindgc.and.nrow.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6001
        goto 6010
6001    mr=mrowr
        if(nrow.eq.1) mr=max(ncol,mr)
        if(mr.lt.ncol) call error('mr.lt.ncol','mxma')
        if(mcola.eq.1.and.mcolb.eq.1.and.mrowa.ge.ncol) then
          call dgemm_x('N','N',ncol,nrow,nlink,-1.0d0,a,mrowa,
     1                    b,max(nlink,mrowb),1.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mcolb.eq.1) then
          call dgemm_x('T','N',ncol,nrow,nlink,-1.0d0,a,mcola,
     1                    b,max(nlink,mrowb),1.0d0,r,mr)
          return
        else if(mcola.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','T',ncol,nrow,nlink,-1.0d0,a,mrowa,
     1                    b,mcolb,1.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mrowb.eq.1) then
          call dgemm_x('T','T',ncol,nrow,nlink,-1.0d0,a,mcola,
     1                    b,mcolb,1.0d0,r,mr)
          return
        else
          call mxmbn_dgm(a,mcola,mrowa,b,mcolb,mrowb,
     1                  r,mcolr,mrowr,ncol,nlink,nrow)
          return
        end if
      else if(mrowr.eq.1) then
        if(nlink.ge.mindgm.and.nrow.ge.mindgm.and.ncol.ge.mindgm)
     1     goto 6002
        if(nlink.ge.mindgl.and.ncol.ge.mindgm2.and.nrow.ge.mindgm2)
     1     goto 6002
        if(nrow.ge.mindgc.and.ncol.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6002
        if(ncol.ge.mindgr.and.nrow.ge.mindgm2.and.nlink.ge.mindgm2)
     1     goto 6002
        goto 6010
6002    mr=mcolr
        if(ncol.eq.1) mr=max(nrow,mr)
        if(mr.lt.nrow) call error('mr.lt.nrow','mxma')
        if(mcola.eq.1.and.mcolb.eq.1) then
         call dgemm_x('T','T',nrow,ncol,nlink,-1.0d0,b,max(nlink,mrowb),
     1                     a,mrowa,1.0d0,r,mr)
         return
        else if(mrowa.eq.1.and.mcolb.eq.1) then
         call dgemm_x('T','N',nrow,ncol,nlink,-1.0d0,b,max(nlink,mrowb),
     1                     a,mcola,1.0d0,r,mr)
         return
        else if(mcola.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','T',nrow,ncol,nlink,-1.0d0,b,mcolb,
     1                     a,mrowa,1.0d0,r,mr)
          return
        else if(mrowa.eq.1.and.mrowb.eq.1) then
          call dgemm_x('N','N',nrow,ncol,nlink,-1.0d0,b,mcolb,
     1                     a,mcola,1.0d0,r,mr)
          return
        end if
      end if
      if(nlink.lt.mindgm.or.nrow.lt.mindgm.or.ncol.lt.mindgm) goto 6010
      call mxmbn_dgm(a,mcola,mrowa,b,mcolb,mrowb,
     1               r,mcolr,mrowr,ncol,nlink,nrow)
      return
6010  continue
#endif
c
      if(nlink.eq.1) then
        jj=1
        ijj=1
        do 6210 j=1,nrow
        ii=1
        ij=ijj
        if(b(jj).ne.0.0d0) then
        do 6110 i=1,ncol
        r(ij)=r(ij)-a(ii)*b(jj)
        ii=ii+mcola
 6110   ij=ij+mcolr
        end if
        ijj=ijj+mrowr
 6210   jj=jj+mrowb
        return
      end if
c... r(ncol,nrow)=r(ncol,nrow)-a(ncol,nlink)*b(nlink,nrow) matrix mult
      iaa=1
      irr=1
      nrest=mod(nrow,i4)
      if(nrest.eq.1.and.nrow.gt.1) nrest=5
      nre=nrow-nrest
      ncest=mod(ncol,i4)
      if(ncest.eq.1.and.ncol.gt.1) ncest=5
      nce=ncol-ncest
      do 300 i=1,nce,4
      ibb=1
      ir1=irr
      ir2=ir1+mcolr
      ir3=ir2+mcolr
      ir4=ir3+mcolr
      irr=ir4+mcolr
      do 200 j=1,nre,4
      ib1=ibb
      ib2=ib1+mrowb
      ib3=ib2+mrowb
      ib4=ib3+mrowb
      ibb=ib4+mrowb
      ia1=iaa
      ia2=ia1+mcola
      ia3=ia2+mcola
      ia4=ia3+mcola
      s11=0.0d0
      s21=0.0d0
      s31=0.0d0
      s41=0.0d0
      s12=0.0d0
      s22=0.0d0
      s32=0.0d0
      s42=0.0d0
      s13=0.0d0
      s23=0.0d0
      s33=0.0d0
      s43=0.0d0
      s14=0.0d0
      s24=0.0d0
      s34=0.0d0
      s44=0.0d0
      do 100 k=1,nlink
      s11=a(ia1)*b(ib1)+s11
      s21=a(ia2)*b(ib1)+s21
      s31=a(ia3)*b(ib1)+s31
      s41=a(ia4)*b(ib1)+s41
      ib1=ib1+mcolb
      s12=a(ia1)*b(ib2)+s12
      s22=a(ia2)*b(ib2)+s22
      s32=a(ia3)*b(ib2)+s32
      s42=a(ia4)*b(ib2)+s42
      ib2=ib2+mcolb
      s13=a(ia1)*b(ib3)+s13
      s23=a(ia2)*b(ib3)+s23
      s33=a(ia3)*b(ib3)+s33
      s43=a(ia4)*b(ib3)+s43
      ib3=ib3+mcolb
      s14=a(ia1)*b(ib4)+s14
      s24=a(ia2)*b(ib4)+s24
      s34=a(ia3)*b(ib4)+s34
      s44=a(ia4)*b(ib4)+s44
      ib4=ib4+mcolb
      ia1=ia1+mrowa
      ia2=ia2+mrowa
      ia3=ia3+mrowa
100   ia4=ia4+mrowa
      r(ir1)=r(ir1)-s11
      r(ir2)=r(ir2)-s21
      r(ir3)=r(ir3)-s31
      r(ir4)=r(ir4)-s41
      r(ir1+mrowr)=r(ir1+mrowr)-s12
      r(ir2+mrowr)=r(ir2+mrowr)-s22
      r(ir3+mrowr)=r(ir3+mrowr)-s32
      r(ir4+mrowr)=r(ir4+mrowr)-s42
      r(ir1+2*mrowr)=r(ir1+2*mrowr)-s13
      r(ir2+2*mrowr)=r(ir2+2*mrowr)-s23
      r(ir3+2*mrowr)=r(ir3+2*mrowr)-s33
      r(ir4+2*mrowr)=r(ir4+2*mrowr)-s43
      r(ir1+3*mrowr)=r(ir1+3*mrowr)-s14
      r(ir2+3*mrowr)=r(ir2+3*mrowr)-s24
      r(ir3+3*mrowr)=r(ir3+3*mrowr)-s34
      r(ir4+3*mrowr)=r(ir4+3*mrowr)-s44
      ir1=ir1+4*mrowr
      ir2=ir2+4*mrowr
      ir3=ir3+4*mrowr
      ir4=ir4+4*mrowr
200   continue
      if(nrest.eq.0) goto 300
      goto (201,202,203,300,205),nrest
201     ib1=ibb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        do 101 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
101     ia4=ia4+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        goto 300
c
202     ib1=ibb
        ib2=ib1+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        do 102 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
102     ia4=ia4+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        r(ir4+mrowr)=r(ir4+mrowr)-s42
        goto 300
c
203     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s13=0.0d0
        s23=0.0d0
        s33=0.0d0
        s43=0.0d0
        do 103 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
103     ia4=ia4+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        r(ir4+mrowr)=r(ir4+mrowr)-s42
        r(ir1+2*mrowr)=r(ir1+2*mrowr)-s13
        r(ir2+2*mrowr)=r(ir2+2*mrowr)-s23
        r(ir3+2*mrowr)=r(ir3+2*mrowr)-s33
        r(ir4+2*mrowr)=r(ir4+2*mrowr)-s43
        goto 300
c
205     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s13=0.0d0
        s23=0.0d0
        s33=0.0d0
        s43=0.0d0
        s14=0.0d0
        s24=0.0d0
        s34=0.0d0
        s44=0.0d0
        s15=0.0d0
        s25=0.0d0
        s35=0.0d0
        s45=0.0d0
        do 105 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        s25=a(ia2)*b(ib5)+s25
        s35=a(ia3)*b(ib5)+s35
        s45=a(ia4)*b(ib5)+s45
        ib5=ib5+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
105     ia4=ia4+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        r(ir4+mrowr)=r(ir4+mrowr)-s42
        r(ir1+2*mrowr)=r(ir1+2*mrowr)-s13
        r(ir2+2*mrowr)=r(ir2+2*mrowr)-s23
        r(ir3+2*mrowr)=r(ir3+2*mrowr)-s33
        r(ir4+2*mrowr)=r(ir4+2*mrowr)-s43
        r(ir1+3*mrowr)=r(ir1+3*mrowr)-s14
        r(ir2+3*mrowr)=r(ir2+3*mrowr)-s24
        r(ir3+3*mrowr)=r(ir3+3*mrowr)-s34
        r(ir4+3*mrowr)=r(ir4+3*mrowr)-s44
        r(ir1+4*mrowr)=r(ir1+4*mrowr)-s15
        r(ir2+4*mrowr)=r(ir2+4*mrowr)-s25
        r(ir3+4*mrowr)=r(ir3+4*mrowr)-s35
        r(ir4+4*mrowr)=r(ir4+4*mrowr)-s45
300   iaa=iaa+4*mcola
      if(ncest.eq.0) return
      goto (301,302,303,304,305),ncest
301     ibb=1
        ir1=irr
        do 2001 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        s11=0.0d0
        s12=0.0d0
        do 1001 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        ib2=ib2+mcolb
1001    ia1=ia1+mrowa
        r(ir1)=r(ir1)-s11
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        ir1=ir1+2*mrowr
2001    continue
        if(mod(nrow,i2).eq.0) return
        s11=0.0d0
        do 1011 k=1,nlink
        s11=a(iaa)*b(ibb)+s11
        ibb=ibb+mcolb
1011    iaa=iaa+mrowa
        r(ir1)=r(ir1)-s11
        return
c
302     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        do 2002 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        s11=0.0d0
        s21=0.0d0
        s12=0.0d0
        s22=0.0d0
        do 1002 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        ia1=ia1+mrowa
1002    ia2=ia2+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
2002    continue
        if(mod(nrow,i2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        s11=0.0d0
        s21=0.0d0
        do 1012 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        ia1=ia1+mrowa
        s21=a(ia2)*b(ibb)+s21
        ia2=ia2+mrowa
1012    ibb=ibb+mcolb
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        return
c
303     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        do 2003 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        do 1003 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1003    ia3=ia3+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
2003    continue
        if(mod(nrow,i2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        do 1013 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1013    ia3=ia3+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
304     return
305     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        ir4=ir3+mcolr
        ir5=ir4+mcolr
        do 2005 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s51=0.0d0
        s12=0.0d0
        s22=0.0d0
        s32=0.0d0
        s42=0.0d0
        s52=0.0d0
        do 1005 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
1005    ia5=ia5+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir5)=r(ir5)-s51
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        r(ir4+mrowr)=r(ir4+mrowr)-s42
        r(ir5+mrowr)=r(ir5+mrowr)-s52
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
        ir4=ir4+2*mrowr
        ir5=ir5+2*mrowr
2005    continue
        if(mod(nrow,i2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0.0d0
        s21=0.0d0
        s31=0.0d0
        s41=0.0d0
        s51=0.0d0
        do 1015 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        s41=a(ia4)*b(ibb)+s41
        s51=a(ia5)*b(ibb)+s51
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
        ia5=ia5+mrowa
1015    continue
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir5)=r(ir5)-s51
        return
      end

#endif


!> vector-vector outer product
!>
!> calculates the \f${\tt ncol}\times{\tt nrow}\f$ matrix product
!> \f$ \bf R = A B^\dagger\f$.
!> \param a real array containing the vector \f$\bf A\f$; note that the
!>       vector may be supplied with non-unit stride, and with
!>       enlarged leading dimension, under the control of the parameter
!>       \c mrowa described below.
!> \param mrowa memory increment in \c a between adjacent elements
!>       of \f$\bf A\f$.
!> \param b real array containing the vector \f$\bf B\f$
!> \param mrowb memory increment in \c b between adjacent elements
!>       of \f$\bf B\f$.
!> \param r real array containing the matrix \f$\bf R\f$; \c r is cleared
!>       to zero at the start of \c vxva, and so any previous contents
!>       are lost.
!> \param mcolr memory increment in \c r between adjacent rows
!>       of \f$\bf R\f$.
!> \param mrowr memory increment in \c r between adjacent columns of
!>       \f$\bf R\f$.
!> \param ncol number of rows in \f$\bf R\f$ and \f$\bf A\f$.
!> \param nrow number of columns in \f$\bf R\f$ and \f$\bf B^\dagger\f$.
!>
!> Thus standard use of \c vxva would take the form
!> \verbatim call vxva (a, 1, b, 1, r, 1, ncol, ncol, nrow). \endverbatim
c--------------------------------------------------------------------
      subroutine vxva (a, mrowa, b, mrowb, r, mcolr, mrowr,
     +     ncol, nrow)
c--------------------------------------------------------------------
      implicit double precision (a-h, o-z)
      integer mrowa, mrowb, mcolr, mrowr, ncol, nrow
      double precision, dimension(*) :: r, a, b

      irowb = 1

      do i = 1, nrow
       irowa = 1

       belt = b(irowb)          ! Element of B for this iteration

       do j = 1, ncol
        ioffset = 1 + ((i - 1) * mrowr) + ((j - 1) * mcolr)
        r(ioffset) = belt * a(irowa)
        irowa = irowa + mrowa
       end do

       irowb = irowb + mrowb
      end do

      end

!> vector-vector outer product
!>
!> calculates the \f${\tt ncol}\times{\tt nrow}\f$ matrix product
!> \f$\bf R = R + A B^\dagger\f$.
!> Parameters are exactly as for \c vxva; the only difference is that
!>       \c r is not cleared to zero at the start, but the matrix
!>       product is added to the previous contents of \c r.
c--------------------------------------------------------------------
      subroutine vxvb (a, mrowa, b, mrowb, r, mcolr, mrowr,
     +     ncol, nrow)
c--------------------------------------------------------------------
      implicit double precision (a-h, o-z)
      dimension r(1), a(1), b(1)

      irowb = 1

      do i = 1, nrow
       irowa = 1

       belt = b(irowb)          ! Element of B for this iteration

       do j = 1, ncol
        ioffset = 1 + ((i - 1) * mrowr) + ((j - 1) * mcolr)
        r(ioffset) = r(ioffset) + (belt * a(irowa))
        irowa = irowa + mrowa
       end do

       irowb = irowb + mrowb
      end do

      end
