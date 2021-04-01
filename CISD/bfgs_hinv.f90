  subroutine bfgs_hinv(hinveq,ninter,direction,dgra,log,hu_column)

!****************************************************************************
!*     Broyden-Fletcher-Goldfarb-Shanno update of inverse hessian           *
!*                                                                          *
!*     HInv(k+1) = HInv(k) + ( 1. + <dg|HInv(k)|dg> / <dx|dg> ) *           *
!*                                                                          *
!*     |dx><dx| / <dx|dg> - ( |dx><dg|HInv(k) + HInv(k)|dg><dx| ) / <dx|dg> *
!****************************************************************************
  implicit double precision (a-h,o-z)
  double precision :: hinveq(ninter,ninter), dgra(ninter)
  double precision :: dgra_column(ninter,1), dgra_row(1,ninter)
  double precision :: direction(ninter),direction_column(ninter,1),direction_row(1,ninter)
  double precision :: hinveq_temp(ninter,ninter),hu(ninter),hu_column(ninter,1),hu_row(1,ninter)
  double precision :: uhu(1,1), alpha(1,1), a, out_product_a(ninter,ninter), out_product_b(ninter,ninter)
  double precision :: temp1(ninter,ninter), temp2(ninter,ninter), temp3(ninter,ninter)
  Integer :: log,ninter

  log = 0

  dgra_row = 0.0d0
  dgra_column = 0.0d0
  direction_row = 0.0d0
  direction_column = 0.0d0

  do i=1,ninter
    dgra_column(i,1) = dgra(i)
    dgra_row(1,i) = dgra(i)
    direction_column(i,1) = direction(i)
    direction_row(1,i) = direction(i)
  end do

!-----------<dx|dg>---------------------------

  !alpha = 0.0d0
  !call DGEMM('N','N',1,1,ninter,1,direction_row,1,dgra_column,ninter,1,alpha,1)

  alpha = matmul(direction_row,dgra_column)

  a = 1.d0/alpha(1,1)

  if(dabs(alpha(1,1)).lt.1.d-14) then
    write(*,*) 'BFGS <dx|dg> too small'
    return
  end if
!---------------------------------------------

!*     form HInvEq(k)|dg>

  !hu_column = 0.0d0
  !call DGEMM('N','N',ninter,1,ninter,1,hinveq,ninter,dgra_column,ninter,1,hu_column,ninter)

  !hu = 0.d0
  !call DGEMV('N',ninter,ninter,1,hinveq,ninter,dgra,1,1,hu,1)

  hu_column = matmul(hinveq,dgra_column)


!*     form <dg|HInvEq(k)|dg>

!  uhu = 0.0d0
!  call DGEMM('N','N',1,1,ninter,1,dgra_row,1,hu_column,ninter,1,uhu,1)

  uhu = matmul(dgra_row,hu_column)

  if (uhu(1,1) .gt. 1.d2) then
    log = 1
  end if

!-----------|dx><dg|---------------------------

  !out_product_a = 0.0d0
  !call DGEMM('N','N',ninter,ninter,1,1,direction_column,ninter,dgra_row,1,1,out_product_a,ninter)

  out_product_a = matmul(direction_column,dgra_row)

!-----------|dg><dx|---------------------------

  !out_product_b = 0.0d0
  !call DGEMM('N','N',ninter,ninter,1,1,dgra_column,ninter,direction_row,1,1,out_product_b,ninter)
  out_product_b = matmul(dgra_column,direction_row)


!  temp1 = 0.0d0
!  call DGEMM('N','N',ninter,ninter,1,1,direction_column,ninter,direction_row,1,1,temp1,ninter)
  temp1 = matmul(direction_column,direction_row)

!  temp2 = 0.0d0
!  call DGEMM('N','N',ninter,ninter,ninter,1,out_product_a,ninter,hinveq,ninter,1,temp2,ninter)
  temp2 = matmul(out_product_a,hinveq)

!  temp3 = 0.0d0
!  call DGEMM('N','N',ninter,ninter,ninter,1,hinveq,ninter,out_product_b,ninter,1,temp3,ninter)
  temp3 = matmul(hinveq,out_product_b)

!-------------Updated Hessian------------------

  hinveq = hinveq + (a + uhu(1,1)*a**2)*temp1 - a*(temp2 + temp3)



return
end
