      subroutine bfgs_oep(hinveq,ninter,af,di_1,dgra,hu)

!************************************************************************
!*     Broyden-Fletcher-Goldfarb-Shanno update of inverse hessian       *
!*                                                                      *
!*     HInv(k+1) = HInv(k) + ( 1. + <dg|HInv(k)|dg> / <dx|dg> ) *       *
!*                                                                      *
!*     |dx><dx| / <dx|dg> - ( |dx><dg|HInv(k) + HInv(k)|dg><dx| ) / <dx|*
!************************************************************************
  implicit double precision (a-h,o-z)
  double precision hinveq(ninter*ninter), hu(ninter), af(ninter**2),di_1(ninter), dgra(ninter)
!*
!      call qenter('BFGS')
!      Call recPrt(' BFGS: HInvEq before update',' ',
!     &            HInvEq,nInter,nInter)
!*
!*     form <dx|dg>
!*

      if(abs(ddot(ninter,di_1,1,dgra,1)).lt.1.d-14) then
!         write(6,*) 'BFGS <dx|dg> too small',
!     >        abs(ddot_X(ninter,di_1,1,dgra,1))
         return
      end if

      alpha = 1.d0/ ddot(ninter,di_1,1,dgra,1)

!*
!*     form HInvEq(k)|dg>
!*

 !call DGEMM('N','N',ninter,ninter,ninter,1,hinveq,ninter,dgra,ninter,1,hu,ninter)
      call mxva (hinveq, 1, ninter, dgra, 1, hu, 1, ninter, ninter)
!      Call RecPrt(' BFGS: Hu',' ',Hu,nInter,1)
!*
!*     form |dx><dg|HInvEq(k) and correct hessian
!*

     call vxva (di_1, 1, hu, 1, af, 1, ninter, ninter, ninter)
      call daxpy(ninter*ninter,-alpha,af,1,hinveq,1)
!*
!*     form HInvEq(k)|dg><dx| and correct hessian
!*
      call vxva (hu, 1, di_1, 1, af, 1, ninter, ninter, ninter)
      call daxpy(ninter*ninter,-alpha,af,1,hinveq,1)
!*
!*     form <dg|HInvEq(k)|dg>
!*
      uhu = ddot(ninter,dgra,1,hu,1)
!      Write (6,*) ' uHu=',uHu
!*
!*     form |dx><dx| and correct hessian
!*
      call vxva (di_1, 1, di_1, 1, af, 1, ninter, ninter, ninter)
      call daxpy(ninter*ninter,(1.d0+uhu*alpha)*alpha,af,1,hinveq,1)
!      Write (6,*) '(One+uHu*Alpha)*Alpha=',(One+uHu*Alpha)*Alpha

!      Call recPrt('In BFGS: HInvEq in new basis',' ',
!     &            HInvEq,nInter,nInter)
!      call qexit('BFGS')
      return
      end
