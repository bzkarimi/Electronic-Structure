	Subroutine rhf
		USE array


!=====================================================
!             RHF PROGRAM (CLOSED SHELL)
!=====================================================

!=====================================================
! RHF calculation is done based on the book,
! "modern quantum chemistry" by Szabu and Ostland.
! to better undersand the notation please refer to chapters
! 2 & 3 of this book.
!=====================================================


!=====================================================
!        ORTHOGONALIZATION OF THE BASIS(CANONICAL)
!        X = Us^-0.5, where s = diagonalized S
!-----------------------------------------------------

!------------- Diagonalizing S ----------------

	U = 0.0d0
	do i=1,norbit_t
	  do j=1,norbit_t
	    U(i,j) = S_M(i,j)
	  end do
	end do

	Eigen_value = 0.0d0
	call DSYEV('V', 'U', norbit_t, U, norbit_t, Eigen_value, WORK, LWORK, INFO)

!------------- Matrix Transformation X ----------------

	do i=1,norbit_t
	  Eigen_value_half(i) = (Eigen_value(i))**(0.5d0)
	end do

	do i=1,norbit_t
	  Eigen_value(i) = (Eigen_value(i))**(-0.5d0)
	end do

	X = 0.0d0
	do i=1,norbit_t
	  do j=1,norbit_t
	    X(i,j) = Eigen_value(j)*U(i,j)
	  end do
	end do


!-------------- Making S^-0.5  --------------------

	S_half_diag = 0.0d0

	do i=1,norbit_t
	  S_half_diag(i,i) = Eigen_value(i)
	end do

	S_half_temp = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,U,norbit_t,S_half_diag,norbit_t,beta1,S_half_temp,norbit_t)

	S_half_inv = 0.0d0
	call DGEMM('N','T',norbit_t,norbit_t,norbit_t,alpha1,S_half_temp,norbit_t,U,norbit_t,beta1,S_half_inv,norbit_t)

!--------------Making S^0.5 --------------------

	S_half_diag = 0.0d0

	do i=1,norbit_t
    S_half_diag(i,i) = Eigen_value_half(i)
	end do

	S_half_temp = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,U,norbit_t,S_half_diag,norbit_t,beta1,S_half_temp,norbit_t)

	S_half = 0.0d0
	call DGEMM('N','T',norbit_t,norbit_t,norbit_t,alpha1,S_half_temp,norbit_t,U,norbit_t,beta1,S_half,norbit_t)

	X = S_half_inv
	X_inv = S_half

!	call DGETRF(norbit_t, norbit_t, X_inv, norbit_t, IPIV2, INFO2)
!	call DGETRI(norbit_t, X_inv, norbit_t, IPIV2, WORK2, LWORK2, INFO2)

!--------------Making S^(-1) --------------------

	S_inv = 0.0d0
	call DGEMM('N','T',norbit_t,norbit_t,norbit_t,alpha1,S_half_inv,norbit_t,S_half_inv,norbit_t,beta1,S_inv,norbit_t)
!-------------Writing the Transformation Matrix ----

!	write(2,*) "X Matrix"
!	do i=1,norbit_t
!	  do j=1,norbit_t
!	    write(2, "(f16.8)", advance = "No") X(i,j)
!	  end do
!	write(2,*)
!	end do

!====================================================
!                 SCF ITERATION
!====================================================

	write(*,*) "HF SCF iteration Started"

	write(2,*) "-------------HF SCF iteration--------------"

	step = 0

	do i=1,norbit_t
	  do j=1,norbit_t
	    P_old(i,j) = 0.0d0
	    P_new(i,j) = 0.0d0
	  end do
	end do

	write(2,*) "**********************************************"

	EN_old = 1.0d0

	20 continue

	write(2,*) "Step", step
	step = step+1

!	stp_hf = 1.0d0

!	10 continue

!---------------------Form the 2e part of Fock Matrix-----------
!	write(2,*) "G Matrix"
	do i=1,norbit_t
	  do j=1,norbit_t
			G(i,j) = 0.0d0
	    	do k=1,norbit_t
	      	do l=1,norbit_t
	        	G(i,j) = G(i,j) + P_new(k,l)*(M_2e(i,j,k,l)-0.5d0*M_2e(i,l,k,j))
	      	end do
	    	end do
!	    write(2, "(f16.8)", advance = "No") G(i,j)
	  end do
!	write(2,*)
	end do

!------------------Constructing Fock Matrix--------------


!	write(2,*) "Fock Matrix"
	F = 0.0d0
	do i=1,norbit_t
	  do j=1,norbit_t
	    F(i,j) = H_core(i,j)+G(i,j)
!	    write(2, "(f16.8)", advance = "No") F(i,j)
	  end do
!	write(2,*)
	end do
!	F(1,1) = 0.0d0
	write(2,*)

!-----------------------Electronic Energy----------------

	EN = 0.0d0
	do i=1,norbit_t
	  do j=1,norbit_t
	    EN = EN + 0.5d0*P_new(i,j)*(H_core(i,j)+F(i,j))
	  end do
	end do
	write(2,*) "Electronic Energy = ",EN
	write(2,*)
!	EN_old = EN

	!EN_new = EN
	!if (EN_new-EN_old > 0.d0) then
	!	P_new = P_old*0.3d0 + P_new*0.7d0
		!stp_hf = stp_hf*0.88d0
	!	goto 10
	!end if

!----------------------Transform and Diagonalize---------


!----------------------G = F*X----------------------------

	G = 0.0d0
  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,F,norbit_t,X,norbit_t,beta1,G,norbit_t)

!----------------------F_Prime = XT*G = XT*F*X--------------------

	F_Prime = 0.0d0
  call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,G,norbit_t,beta1,F_Prime,norbit_t)

	do i=1,norbit_t
	  do j=1,norbit_t
	    Eigen_C_prime(i,j) = F_Prime(i,j)
	  end do
	end do

!----------------------Eigenvalues and Eigenvectors of F_Prime---------

	Eigen_value = 0.0d0
	call DSYEV('V', 'U', norbit_t, Eigen_C_prime, norbit_t, Eigen_value, WORK, LWORK, INFO)

!---------------------- X*C_Prime = C---------------------------------

	Eigen_C = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,Eigen_C_prime,norbit_t,beta1,Eigen_C,norbit_t)

!------------------Creating New Density Matrix----------

	P_old = P_new

	do i=1,norbit_t
	  do j=1,norbit_t
	    P_new(i,j) = 0.0d0
	    do k=1,nelec/2
	      	P_new(i,j) = P_new(i,j) + 2.0d0*Eigen_C(i,k)*Eigen_C(j,k)
	    end do
	  end do
	end do

!---------------------Writing F', C', Energy matrix, C, P--------------------

!	write(2,*) "Transformed Fock Matrix"
!	do i=1,norbit_t
!	  do j=1,norbit_t
!	    write(2, "(f16.8)", advance = "No") F_Prime(i,j)
!	  end do
!	write(2,*)
!	end do
!
!	write(2,*) "Transformed C Matrix"
!	do i=1,norbit_t
!	  do j=1,norbit_t
!	    write(2, "(f16.8)", advance = "No") Eigen_C_Prime(i,j)
!	  end do
!	write(2,*)
!	end do
!
!	write(2,*) "Diagonalized F Matrix"
!	do i=1,norbit_t
!	  do j=1,norbit_t
!	    if (i.eq.j) then
!	      write(2, "(f16.8)", advance = "No") Eigen_value(i)
!	    else
!	      write(2, "(f16.8)", advance = "No") 0.0d0
!	    end if
!	  end do
!	write(2,*)
!	end do
!
!	write(2,*) "C Matrix"
!	do i=1,norbit_t
!	  do j=1,norbit_t
!	    write(2, "(f16.8)", advance = "No") Eigen_C(i,j)
!	  end do
!	write(2,*)
!	end do
!

!	write(*,*) "HF Density Matrix:"
!	do i=1,norbit_t
!	  do j=1,norbit_t
!	    write(*, "(f16.8)", advance = "No") P_new(i,j)
!	  end do
!	write(*,*)
!	end do


!---------------------Convergence of Density Matrix(Using Standard Deviation)---------

	delta_HF = 0.0d0
	do i=1,norbit_t
	  do j=1,norbit_t
	    delta_HF = delta_HF +(P_new(i,j)-P_old(i,j))**2.0d0
	  end do
	end do

	delta_HF = dsqrt(delta_HF/(norbit_t*norbit_t))
	write(2,*) "delta_density = ",delta_HF


	delta_E = dabs(EN-EN_old)
	write(2,*) "delta_E = ",delta_E


	write(2,*) "*********************************************"

!-----------------Check for Convergence-------------

	if (delta_HF < threshold_HF .and. delta_E < E_thr) then

		EN_t = EN
		do i=1,natom-1
			do j=i+1,natom
				call Distance(R(i,1),R(i,2),R(i,3),R(j,1),R(j,2),R(j,3),d)
				EN_t = EN_t + Z(i)*Z(j)/d
			end do
		end do

		EN_t = EN_t/natom
		write(2,*) "Total Energy per Atom = ",EN_t

	  open(Unit = 7,File = 'Energy')
		write(7,*) a, EN_t


	  write(2,*) "HF Electronic Energy = ",EN
	  write(2,*) "HF Total Energy per Atom = ",EN_t
		write(2,*) "========================================="
		write(*,*) "Congratulations! HF Calculation Converged"
		write(*,*) "DMET Calculation Started"
	else if(step .lt. max_it) then
		EN_old = EN
		goto 20
	else
		write(2,*) "HF calculation does not converge in ",max_it," steps!"
		bool = 0
	end if

!==============================================================



	return
	end
