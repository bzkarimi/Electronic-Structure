  Subroutine bfgs_oep
    USE array
!====================================================
!                   OEP
!====================================================

  write(*,*) "OEP Calculation Started!"

!------------Transform the Kinetic Matrix to Orthognalized AO-------

  KE_t_temp = 0.0d0
  call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,KE_t,norbit_t,beta1,KE_t_temp,norbit_t)

  KE_t_prime = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,KE_t_temp,norbit_t,X,norbit_t,beta1,KE_t_prime,norbit_t)

  step = 0

!------Initial Guess for Potential-----------------

  PE_t_new = F_Prime - KE_t_prime

! KE_t_temp = 0.0d0
! call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,PE_t_new,norbit_t,beta1,KE_t_temp,norbit_t)

! PE_t_new = 0.0d0
! call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,KE_t_temp,norbit_t,X,norbit_t,beta1,PE_t_new,norbit_t)

!  write(*,*) "Initial F Matrix in AO orthogonalized basis"
!  do i=1,norbit_t
!    do j=1,norbit_t
!      write(*, "(f16.8)", advance = "No") F_Prime(i,j)
!    end do
!    write(*,"(f16.8)")
!  end do

!  write(*,*) "Initial K Matrix in AO orthogonalized basis"
!  do i=1,norbit_t
!    do j=1,norbit_t
!      write(*, "(f16.8)", advance = "No") KE_t_prime(i,j)
!    end do
!    write(*,"(f16.8)")
!  end do

!  write(*,*) "Initial P Matrix in AO orthogonalized basis"
!  do i=1,norbit_t
!    do j=1,norbit_t
!      write(*, "(f16.8)", advance = "No") PE_t_new(i,j)
!    end do
!    write(*,"(f16.8)")
!  end do


! PE_t_new = 0.0d0
  PE_t_init = PE_t_new
  Ws_old = 0.0d0
  Ws_new = 0.0d0
  Gradient_new = 0.0d0
  direction = 0.0d0

!------Checking the basis orders------------------

  P_HF = 0.0d0
  do k=1,nelec_A/2
    do i=1,norbit_t
      do j=1,norbit_t
        P_HF(i,j,k) = P_HF(i,j,k) + 2.0d0*Eigen_C(i,k)*Eigen_C(j,k)
      end do
    end do
  end do

!	P_fci_A = P_new

!	Open(Unit = 9, File = 'densityin.txt')
!	read(9,*)
!	do i=1,norbit_t
!		read(9,*) (P_new(i,j), j=1,norbit_t)
!	end do

!-------------Set Hessian as I-------------------------

  Hessian_inv = 0.0d0
  do i=1,ninter
    Hessian_inv(i,i) = 1.0d0
  end do

!========================================================================

  write(2,*) "OEP Calculation:"
  Open(Unit = 8, File = 'Ws.txt')
  write(8,*) 'Step','                   ','delta_OEP','         ','Ws'
  write(8,*) '----------','          ','----------','          ','----------'


  KE_f = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_fci_A,norbit_t,KE_t_prime,norbit_t,beta1,KE_f,norbit_t)

  !Hint !
  Tr_K = 0.0d0
  do i=1,norbit
    Tr_K = Tr_K + KE_f(i,i)
  end do

  Tr_K_fci = Tr_K
  write(*,*) "Ws Max:", Tr_K
  write(*,*)

!=============================================================
!                       OEP zeroth step
!=============================================================

  write(2,*) "OEP Step:",step

!-------------------------Vext * (P_new-P_fci_A)-----------------

  del_P = P_orth-P_fci_A
  VP = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,del_P,norbit_t,PE_t_new,norbit_t,beta1,VP,norbit_t)

  KE_f = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_orth,norbit_t,KE_t_prime,norbit_t,beta1,KE_f,norbit_t)

  !Hint !

  Tr_K = 0.0d0
  do i=1,norbit
    Tr_K = Tr_K + KE_f(i,i)
  end do

  !Hint !
  Tr_P = 0.0d0
  do i=1,norbit
    Tr_P = Tr_P + VP(i,i)
!Tr_P = PE_t_new(1,1) * (P_orth(1,1)-P_fci_A(1,1))

  end do

  Ws_old = Ws_new
  Ws_new = Tr_P + Tr_K

! Note !
  delta_OEP = 0.0d0
  do i=1,norbit_t
    if(i<norbit+1) then
      do j=1,norbit
        delta_OEP = delta_OEP+(P_orth(i,j)-P_fci_A(i,j))**2.0d0
      end do
  !  else
  !    do j=3,4
  !      delta_OEP = delta_OEP+(P_orth(i,j)-P_fci_A(i,j))**2.0d0
  !    end do
    end if
  end do

  delta_OEP = dsqrt(delta_OEP/(norbit_t*norbit_t))

  write(8,*) Step,delta_OEP,Ws_new

  write(2,*) 'Ws = ', Ws_new
  write(2,*) 'integral(V*(p-pin)) = ', Tr_P
  write(2,*) 'Kinetic Energy = ', Tr_K

  write(2,*) 'delta_OEP = ', delta_OEP

  write(2,*) "Initial Density Matrix AO orthogonalized"
  do i=1,norbit_t
    do j=1,norbit_t
      write(2, "(f16.8)", advance = "No") P_orth(i,j)
    end do
    write(2,"(f16.8)")
  end do

!	write(2,*) "Initial Hessian_inv"
!	do i=1,ninter
!		do j=1,ninter
!			write(2, "(f16.8)", advance = "No") Hessian_inv(i,j)
!		end do
!		write(2,"(f16.8)")
!	end do

  write(2,*) '==============================='
  write(2,*) '==============================='



!=============================================================
!                         OEP Cycle
!=============================================================

  50 continue

  step = step +1

  write(2,*) "OEP Step:",step

!---------------------Calculating Gradient----------------

  Gradient_old = Gradient_new

!----------------------------------------------------------------------
!		Using Upper triangle to make sure the potential remains symmetric
!----------------------------------------------------------------------

! Note !
  counti = 1
  do i=1,norbit_t
    if (i<norbit+1) then
      do j=1,norbit
        Gradient_new(counti) = P_orth(i,j) - P_fci_A(i,j)
        counti = counti + 1
      end do
  !  else
  !    do j=3,4
  !      Gradient_new(counti) = P_orth(i,j) - P_fci_A(i,j)
  !      counti = counti + 1
  !    end do
    end if
  end do


!-------------delta gradient-------------------------------

  dgrad = Gradient_new - Gradient_old

!--------------- increasing the rank of grad matrix--------

  do i=1,ninter
    Grad_new_column(i,1) = Gradient_new(i)
  end do

	write(2,*)  "Gradient:"
	do i=1,ninter
		write(2,*) Gradient_new(i)
	end do

!---------------direction = H^(-1)*Grad--------

  direction_column = 0.0d0
  call DGEMM('N','N',ninter,1,ninter,alpha1,Hessian_inv,ninter,Grad_new_column,ninter,beta1,direction_column,ninter)

!---------------decreasing the rank of direction matrix--------

  counti = 1
  do i=1,ninter
    direction_1(counti) = direction_column(counti,1)
    counti = counti + 1
  end do

	write(2,*)  "direction_column:"
	do i=1,ninter
		write(2,*) direction_column(i,1)
	end do

! Note !

  counti = 1
  do i=1,norbit_t
    if(i<norbit+1) then
      do j=1,norbit
        direction(i,j) = direction_column(counti,1)
        counti = counti + 1
      end do
    !else
    !  do j=3,4
    !    direction(i,j) = direction_column(counti,1)
    !    counti = counti + 1
    !  end do
    end if
  end do

	write(2,*) "Upper-Triangular direction matrix:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") direction(i,j)
		end do
		write(2,"(f16.8)")
	end do

  PE_t_old = PE_t_new


!--------------Updating potential------------------------------
  stp = 1.0d0

  80 continue

  PE_t_new = PE_t_old

  write(2,*) "Old Potential Matrix AO ortho:"
  do i=1,norbit_t
    do j=1,norbit_t
      write(2, "(f16.8)", advance = "No") PE_t_new(i,j)
    end do
    write(2,"(f16.8)")
  end do

! Note !

  do i=1,norbit_t
    if (i<norbit+1) then
      do j=1,norbit
        PE_t_new(i,j) = PE_t_new(i,j) + stp*direction(i,j)
      end do
  !  else
  !    do j=3,4
  !      PE_t_new(i,j) = PE_t_new(i,j) + stp*direction(i,j)
  !    end do
    end if
  end do


  write(2,*) "Symmetric Potential Matrix:"
  do i=1,norbit_t
    do j=1,norbit_t
      write(2, "(f16.8)", advance = "No") PE_t_new(i,j)
    end do
    write(2,"(f16.8)")
  end do

  ! Note !

!  do j=1,norbit_t
!  	do i=j+1,norbit_t
!  		PE_t_new(i,j) = PE_t_new(j,i)
!  	end do
!  end do


!  write(2,*) "Symmetric Potential Matrix:"
!  do i=1,norbit_t
!    do j=1,norbit_t
!      write(2, "(f16.8)", advance = "No") PE_t_new(i,j)
!    end do
!    write(2,"(f16.8)")
!  end do

  write(2,*) "Hessian_inv*Grad:"
  do i=1,norbit_t
    do j=1,norbit_t
      write(2, "(f16.8)", advance = "No") direction(i,j)
    end do
    write(2,"(f16.8)")
  end do
  write(2,"(f16.8)")

!------Transform back to Non-orthogonalized AO--------

!  KE_t_temp = 0.0d0
!  call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X_inv,norbit_t,PE_t_new,norbit_t,beta1,KE_t_temp,norbit_t)

!  PE_t_new = 0.0d0
!  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,KE_t_temp,norbit_t,X_inv,norbit_t,beta1,PE_t_new,norbit_t)

!  write(2,*) "Potential Matrix Non-orthogonalized AO"
!  do i=1,norbit_t
!    do j=1,norbit_t
!      write(2, "(f16.8)", advance = "No") PE_t_new(i,j)
!    end do
!    write(2,"(f16.8)")
!  end do

  F_Prime = KE_t_prime + PE_t_new


!  F = KE_t + PE_t_new

	write(2,*) "F' Matrix"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") F_Prime(i,j)
		end do
		write(2,"(f16.8)")
	end do
!----------------------Transform and Diagonalize---------

!----------------------G = F*X----------------------------

!  G = 0.0d0
!  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,F,norbit_t,X,norbit_t,beta1,G,norbit_t)

!----------------------F_Prime = XT*G = XT*F*X--------------------

!  F_Prime = 0.0d0
!  call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,G,norbit_t,beta1,F_Prime,norbit_t)

  do i=1,norbit_t
    do j=1,norbit_t
      Eigen_C_prime(i,j) = F_Prime(i,j)
    end do
  end do

  write(2,*) "Transformed Fock Matrix"
  do i=1,norbit_t
    do j=1,norbit_t
      write(2, "(f16.8)", advance = "No") F_Prime(i,j)
    end do
    write(2,"(f16.8)")
  end do

!----------------------Eigenvalues and Eigenvectors of F_Prime---------

  Eigen_value = 0.0d0
  call DSYEV('V', 'U', norbit_t, Eigen_C_prime, norbit_t, Eigen_value, WORK, LWORK, INFO)

!---------------------- X*C_Prime = C---------------------------------

  Eigen_C = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,Eigen_C_prime,norbit_t,beta1,Eigen_C,norbit_t)


  write(2,*) "Eigenvectors of F' Matrix"
  do i=1,norbit_t
    do j=1,norbit_t
      write(2, "(f16.8)", advance = "No") Eigen_C_prime(i,j)
    end do
    write(2,"(f16.8)")
  end do


!====================================================
!                Check the order of MO
!====================================================

  delta_den = 0.0d0
  minimum = 1.0d4

  do k=1,norbit_t

  	P_new = 0.0d0
  	do i=1,norbit_t
  		do j=1,norbit_t
  			P_new(i,j) = P_new(i,j) + 2.0d0*Eigen_C(i,k)*Eigen_C(j,k)

        delta_den(k) = delta_den(k)+dabs(P_new(i,j)-P_HF(i,j,1))
			end do
		end do

		if(delta_den(k) .lt. minimum) then
			minimum = delta_den(k)
			ktemp = k
		end if

	end do

	write(2,*) "********************"

	write(2,*) "delta density"
	do i=1,norbit_t
		write(2,*) delta_den(i)
	end do
	write(2,*) "********************"

!	write(*,*) 'step=',step-1,'fragment orbital(should be 1)=',ktemp

	write(2,*) "Old C' Matrix"
	do i=1,norbit_t
		do j=1,norbit_t
		  write(2, "(f16.10)", advance = "No") Eigen_C_prime(i,j)
	  end do
		write(2,*)
	end do

!------------------Swap columns if there is any change in order-----

	if (ktemp /= 1) then

		70 continue

		do i=1,norbit_t
			Call Swap(Eigen_C(i,ktemp),Eigen_C(i,1))
      Call Swap(Eigen_C_prime(i,ktemp),Eigen_C_prime(i,1))
	  end do

		ktemp = ktemp-1
		if (ktemp /= 1) then
			goto 70
		end if
	end if

	write(2,*) "New C' Matrix"
	do i=1,norbit_t
  	do j=1,norbit_t
  		write(2, "(f16.10)", advance = "No") Eigen_C_prime(i,j)
  	end do
  write(2,*)
  end do


  P_new = 0.0d0
  do i=1,norbit_t
    do j=1,norbit_t
      do k=1,nelec/2
        P_new(i,j) = P_new(i,j) + 2.0d0*Eigen_C(i,k)*Eigen_C(j,k)
      end do
    end do
  end do

  P_orth = 0.0d0
  do i=1,norbit_t
    do j=1,norbit_t
      do k=1,nelec/2
        P_orth(i,j) = P_orth(i,j) + 2.0d0*Eigen_C_prime(i,k)*Eigen_C_prime(j,k)
      end do
    end do
  end do

  PS_check = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_new,norbit_t,S_M,norbit_t,beta1,PS_check,norbit_t)

  Tr_check = 0.0d0
  do i=1,norbit_t
    Tr_check = Tr_check + PS_check(i,i)
  end do
  write(*,*) "Trace(P*S):"
  write(*,"(f16.8)") Tr_check

!  do j=1,int(norbit_t/2)
!    do i=1,norbit_t
!      call swap(P_orth(i,j),P_orth(i,norbit_t-j+1))
!    end do
!  end do
!  do i=1,int(norbit_t/2)
!    do j=1,norbit_t
!      call swap(P_orth(i,j),P_orth(norbit_t-i+1,j))
!    end do
!  end do

  write(2,*) "Orthogonalized New Density Matrix"
  do i=1,norbit_t
    do j=1,norbit_t
      write(2, "(f16.8)", advance = "No") P_orth(i,j)
    end do
    write(2,"(f16.8)")
  end do


!  KE_t_temp = 0.0d0
!  call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,PE_t_new,norbit_t,beta1,KE_t_temp,norbit_t)

!  PE_t_new = 0.0d0
!  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,KE_t_temp,norbit_t,X,norbit_t,beta1,PE_t_new,norbit_t)

  write(2,*) "Potential AO ortho"
  do i=1,norbit_t
    do j=1,norbit_t
      write(2, "(f16.8)", advance = "No") PE_t_new(i,j)
    end do
    write(2,"(f16.8)")
  end do

!-------------------------Vext * (P_new-P_fci_A)-----------------

  del_P = P_orth-P_fci_A
  VP = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,del_P,norbit_t,PE_t_new,norbit_t,beta1,VP,norbit_t)

  KE_f = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_orth,norbit_t,KE_t_prime,norbit_t,beta1,KE_f,norbit_t)

  ! Hint !

  Tr_K = 0.0d0
  do i=1,norbit
    Tr_K = Tr_K + KE_f(i,i)
  end do

  ! Hint !
  Tr_P = 0.0d0
  do i=1,norbit
    Tr_P = Tr_P + VP(i,i)
    !Tr_P = PE_t_new(1,1) * (P_orth(1,1)-P_fci_A(1,1))
  end do

  Ws_old = Ws_new
  Ws_new = Tr_P + Tr_K

  if (Ws_new > Tr_K_fci) then
    PE_t_new = PE_t_old
    goto 100
  end if
  !Ws_new = Tr_K_fci - (Ws_new-Tr_K_fci)
  !	stp = 0.5d0*stp

  !	if (stp >= 1.d-14) then
  !		goto 80
  !	else
  !		write(*,*) 'OEP does not converge!', Ws_new
  !		goto 100
  !	end if

  if ( (Ws_new - Ws_old) < 0.0d0) then

    stp = 0.5d0*stp
    Ws_new = Ws_old

    if (stp >= 1.d-14) then
      goto 80
    else
      write(*,*) 'You are at the maximum!', Ws_new
      goto 100
    end if

  end if

! Note!
  delta_OEP = 0.0d0
  do i=1,norbit_t
    if (i<norbit+1) then
      do j=1,norbit
        delta_OEP = delta_OEP+(P_orth(i,j)-P_fci_A(i,j))**2.0d0
      end do
  !  else
  !    do j=3,4
  !      delta_OEP = delta_OEP+(P_orth(i,j)-P_fci_A(i,j))**2.0d0
  !    end do
    end if
  end do

  delta_OEP = dsqrt(delta_OEP/(norbit_t*norbit_t))

  write(8,*) Step,delta_OEP,Ws_new

  write(2,*) 'Ws = ', Ws_new
  write(2,*) 'integral(V*(p-pin)) = ', Tr_P
  write(2,*) 'delta_OEP = ', delta_OEP
  write(2,*) 'Kinetic Energy = ', Tr_K

!****************************************************************************
!*     Broyden-Fletcher-Goldfarb-Shanno update of inverse hessian           *
!*                                                                          *
!*     HInv(k+1) = HInv(k) + ( 1. + <dg|HInv(k)|dg> / <dx|dg> ) *           *
!*                                                                          *
!*     |dx><dx| / <dx|dg> - ( |dx><dg|HInv(k) + HInv(k)|dg><dx| ) / <dx|dg> *
!****************************************************************************

  do i=1,ninter
    dgra_column(i,1) = dgrad(i)
    dgra_row(1,i) = dgrad(i)
    direction_row(1,i) = direction_column(i,1)
  end do

!-----------<dx|dg>---------------------------

  alpha_bracket = 0.0d0
  call DGEMM('N','N',1,1,ninter,alpha1,direction_row,1,dgra_column,ninter,beta1,alpha_bracket,1)


  write(2,*) "dgra_column:"
  do i=1,ninter
    write(2,*) dgra_column(i,1)
  end do
  write(2,*) '******************'
  write(2,*)

  write(2,*) "direction_row:"
  do i=1,ninter
  	write(2,*) direction_row(1,i)
  end do
  write(2,*) '******************'

!-----------1/<dx|dg>---------------------------

  alpha_inv = 1.0d0/alpha_bracket(1,1)

  write(2,*) 'step=',step,'<dx|dg>=',alpha_bracket(1,1)

  if(dabs(alpha_bracket(1,1)).lt.1.d-14) then
  		write(*,*) 'BFGS <dx|dg> too small',alpha_bracket(1,1)
  end if

!-------------------HInvEq(k)|dg>--------------------------

  hu_column = 0.0d0
  call DGEMM('N','N',ninter,1,ninter,alpha1,Hessian_inv,ninter,dgra_column,ninter,beta1,hu_column,ninter)

!-------------------<dg|HInvEq(k)|dg>-----------------------

  uhu = 0.0d0
  call DGEMM('N','N',1,1,ninter,alpha1,dgra_row,1,hu_column,ninter,beta1,uhu,1)

  write(2,*) '<dg|HInvEq(k)|dg> =', uhu(1,1)


!-----------|dx><dg|---------------------------

  out_product_a = 0.0d0
  call DGEMM('N','N',ninter,ninter,1,alpha1,direction_column,ninter,dgra_row,1,beta1,out_product_a,ninter)

!-----------|dg><dx|---------------------------

  out_product_b = 0.0d0
  call DGEMM('N','N',ninter,ninter,1,alpha1,dgra_column,ninter,direction_row,1,beta1,out_product_b,ninter)

!-----------|dx><dx|---------------------------

  temp1 = 0.0d0
  call DGEMM('N','N',ninter,ninter,1,alpha1,direction_column,ninter,direction_row,1,beta1,temp1,ninter)

!-----------|dx><dg|HInvEq---------------------------

  temp2 = 0.0d0
  call DGEMM('N','N',ninter,ninter,ninter,alpha1,out_product_a,ninter,Hessian_inv,ninter,beta1,temp2,ninter)

!-----------HInvEq|dg><dx|----------------------------

  temp3 = 0.0d0
  call DGEMM('N','N',ninter,ninter,ninter,alpha1,Hessian_inv,ninter,out_product_b,ninter,beta1,temp3,ninter)


!  write(2,*) "|dx><dx|"
!  do i=1,ninter
!    do j=1,ninter
!  		write(2, "(f16.8)", advance = "No") temp1(i,j)
!    end do
!    write(2,"(f16.8)")
!  end do


!  write(2,*) "|dx><dg|HInvEq"
!  do i=1,ninter
!    do j=1,ninter
!      write(2, "(f16.8)", advance = "No") temp2(i,j)
!    end do
!    write(2,"(f16.8)")
!  end do

!  write(2,*) "HInvEq|dg><dx|"
!  do i=1,ninter
!    do j=1,ninter
!      write(2, "(f16.8)", advance = "No") temp3(i,j)
!    end do
!    write(2,"(f16.8)")
!  end do

!  write(2,*) "old Hessian_inv"
!  do i=1,ninter
!    do j=1,ninter
!    	write(2, "(f16.8)", advance = "No") Hessian_inv(i,j)
!    end do
!    write(2,"(f16.8)")
!  end do

!-------------Updated Hessian------------------

do i=1,ninter
  do j=1,ninter
    Hessian_inv(i,j) = Hessian_inv(i,j) + (alpha_inv*(1.0d0 + uhu(1,1)*alpha_inv))*temp1(i,j) &
    - alpha_inv*(temp2(i,j) + temp3(i,j))

  !  write(*,*) 'update1', (alpha_inv*(1.0d0 + uhu(1,1)*alpha_inv))*temp1(i,j)
  !  write(*,*)  'update2',alpha_inv*(temp2(i,j) + temp3(i,j))
  end do
  !write(*,*) '------------------------------------'

end do

write(*,*) '========================================'



!--------------Updating Hessian--------------------------

  !	write(2,*) "Delta Gradient:"
  !	do i=1,ninter
  !		write(2,*) dgrad(i)
  !	end do

  !	counti = 1
  !	do i=1,ninter
  !		dgra_column(counti,1) = dgrad(counti)
  !		counti = counti + 1
  !	end do



  !	write(2,*) "Delta Gradient:"
  !	do i=1,ninter
  !		write(2,*) dgra_column(i,1)
  !	end do


!		write(2,*) "hu_column:"
!		do i=1,ninter
!			write(2,*) hu_column(i,1)
!		end do

!		write(2,"(f16.8)")


!	write(2,*) "Updated Hessian_inv"
!	do i=1,ninter
!		do j=1,ninter
!			write(2, "(f16.8)", advance = "No") Hessian_inv(i,j)
!		end do
!		write(2,"(f16.8)")
!	end do
  write(2,*) "****************************************"
  write(2,*) "  ************************************  "

  if ( delta_OEP .gt. threshold_OEP .and. step .le. max_it) then

    goto 50

  else if (step .gt. max_it) then
    write(*,*) "It does not converge in ",max_it," steps!"
  else
    write(*,*) "OEP Calculation Converged!"
    write(2,*) "step=", step

  end if

!==============================================================

  100 continue
  return
end
