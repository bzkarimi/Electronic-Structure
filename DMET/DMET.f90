	Program DMET

	USE array

	IMPLICIT NONE

	INTERFACE
		subroutine alloc_array
		end subroutine alloc_array
	end INTERFACE

!====================================================
!                  MAIN PROGRAM
!====================================================

!====================================================
!          INPUT VALUES + ALLOCATE ARRAYS
!====================================================

	call alloc_array


	call input_geo

!----------------Reading coordinates from file----------------

	open(Unit = 4,File ="Input.xyz")
	do i=1,natom
	  read(4,*) atomname,R(i,1),R(i,2),R(i,3)
	end do

!----------------Reading input density from file----------------

!	open(Unit = 5,File ="densityin.txt")
!	read(5,*)
!	do i=1,norbit_t
!		read(5,*) (P_in(i,j), j=1, norbit_t)
!	end do

!---------------------Calculate all integrals------------------

	call all_integrals

!---------------------RHF Calculation--------------------------

	call rhf

!---------------------Checking the RHF convergence-------------

	if (bool == 0) then
		write(*,*) "RHF does not converge."
		goto 100
	end if


!=======================================================
!								          DMET
!=======================================================

	write(2,*) "-----------DMET Calculations---------------"
	write(2,*)
	write(2,*) "-----------Orbital Localization------------"

!--------------Orbital Localization on Site A-----------
!---------atom_basis = num of basis on site A-----------

!---------------------------------------
!				Loop over fragments
!---------------------------------------


	do i=1,norbit
		do j=1,norbit
			S_A(i,j) = S_M(i,j)
		end do
	end do

	do i=1,norbit
		do j=1,norbit_t
			S_prime(i,j) = S_M(i,j)
		end do
	end do

!		write(2,*) 'S_prime:'

!		do i=1,norbit
!			do j=1,norbit_t
!				write(2,"(f16.6)", advance = "No") S_prime(i,j)
!			end do
!			write(2,*)
!		end do


!--------------P = Projection Operator-------------------

!--------------P = S'(T) * S_A(-1) * S'------------------
	call DGETRF(norbit, norbit, S_A, norbit, IPIV, INFO)
	call DGETRI(norbit, S_A, norbit, IPIV, WORK, LWORK, INFO)

	!	S_prime_T = transpose(S_prime)

	Project_temp = 0.0d0
call DGEMM('T','N',norbit_t,norbit,norbit,alpha1,S_prime,norbit,S_A,norbit,beta1,Project_temp,norbit_t)

!Project_temp = matmul(S_prime_T,S_A)

	!	write(2,*) 'Project_temp:'

	!	do i=1,norbit_t
	!		do j=1,norbit
	!			write(2,"(f16.6)", advance = "No") Project_temp(i,j)
	!		end do
	!		write(2,*)
	!	end do

	Project_new = 0.0d0
call DGEMM('N','N',norbit_t,norbit_t,norbit,alpha1,Project_temp,norbit_t,S_prime,norbit,beta1,Project_new,norbit_t)

!Project_new = matmul(Project_temp,S_prime)

!	write(2,*) 'Projection Operator:'

!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(2,"(f16.6)", advance = "No") Project_new(i,j)
!		end do
!		write(2,*)
!	end do

!---------------------------------------------------------------------------
! Change the basis of projection matrix to MO basis	Using occupied orbitals

!-------------------- Gama = C(T) * P * C--------------------


	do i=1,norbit_t
		do j=1,nelec/2
			eigen_c_occ(i,j) = Eigen_C(i,j)
		end do
	end do

	!	eigen_c_occ_T = Transpose(eigen_c_occ)
	Gama_temp = 0.0d0
call DGEMM('T','N',nelec/2,norbit_t,norbit_t,alpha1,eigen_c_occ,norbit_t,Project_new,norbit_t,beta1,Gama_temp,nelec/2)

!Gama_temp = matmul(eigen_c_occ_T,Project_new)

	Gama = 0.0d0
call DGEMM('N','N',nelec/2,nelec/2,norbit_t,alpha1,Gama_temp,nelec/2,eigen_c_occ,norbit_t,beta1,Gama,nelec/2)

!		Gama = matmul(Gama_temp,eigen_c_occ)

	E_value = 0.0d0
	Eigen_Gama = Gama
	call DSYEV('V', 'U', nelec/2, Eigen_Gama, nelec/2, E_value, WORK, LWORK, INFO)

	write(*,*) 'Gama(Projection Operator in Occupied MO basis):'
	do i=1,nelec/2
		do j=1,nelec/2
				write(*,"(E16.8)", advance = "No") Gama(i,j)
		end do
	write(*,*)
	end do
	write(*,*)

!	write(2,*) 'Diagonalized Gama:'
!	do i=1,nelec/2
!		do j=1,nelec/2
!			if (i .eq. j) then
!				Gama_diag(i,j) = E_value(i)
!				write(2,"(E16.8)", advance = "No") E_value(i)
!			else
!				Gama_diag(i,j) = 0.0d0
!				write(2,"(E16.8)", advance = "No") 0.0d0
!			end if
!		end do
!		write(2,*)
!	end do

!	write(2,*)

	write(*,*) 'Gama Eigenvectors:'
	do i=1,nelec/2
		do j=1,nelec/2
			write(*,"(E16.8)", advance = "No") Eigen_Gama(i,j)
		end do
		write(*,*)
	end do
	write(*,*)


!------------------------------------------------------------
!					Localization of the occupied orbitals
!------------------------------------------------------------


	eigen_c_occ_local = 0.0d0
call DGEMM('N','N',norbit_t,nelec/2,nelec/2,alpha1,eigen_c_occ,norbit_t,Eigen_Gama,nelec/2,beta1,eigen_c_occ_local,norbit_t)

!		eigen_c_occ_local = matmul(eigen_c_occ,Eigen_Gama)

	write(*,*) "Localized Orbitals:"
	do i=1,norbit_t
		do j=1,nelec/2
			write(*, "(f16.8)", advance = "No") eigen_c_occ_local(i,j)
		end do
	write(*,"(f16.8)")
	end do
	write(*,*) "Occupied Orbitals:"
	do i=1,norbit_t
		do j=1,nelec/2
			write(*, "(f16.8)", advance = "No") eigen_c_occ(i,j)
		end do
	write(*,"(f16.8)")
	end do
!============================================================
!								   CHECK THE PROGRAM I
!============================================================

! This checks whether the orbitals are correctly localized on
! the fragment or not.
! The Trace(P*S) = num of e on fragment, where P is the density
! made from localized orbitals.


	P_check = 0.0d0

!----------------- Constructing Localized density just for site 1---------------
!---------------- Note: the last column of eigen-vectors = for the first atom---
!---------------- Since they are in the same order as their eigen-values--------

	do i=1,norbit_t
		do j=1,norbit_t
			do k=1,nelec/2
				if (k .eq. 1) then
					P_check(i,j) = P_check(i,j) + 2*eigen_c_occ_local(i,k)*eigen_c_occ_local(j,k)
				end if
			end do
		end do
	end do

	PS_check = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_check,norbit_t,S_M,norbit_t,beta1,PS_check,norbit_t)

!!!	PS_check = matmul(P_check,S_M)

	Tr_check = 0.0d0
	do i=1,norbit_t
		Tr_check = Tr_check + PS_check(i,i)
	end do
	write(*,*) "local Trace(P*S):"
	write(*,"(f16.8)") Tr_check

!	write(2,*) "Localized Density Matrix"
!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(2, "(f16.8)", advance = "No") P_check(i,j)
!		end do
!		write(2,"(f16.8)")
!	end do


!============================================================
!								   Constructing the new Fock (OLD)
!============================================================

!	Gama_B = 0.0d0
!	do i=1,norbit_t
!		do j=1,norbit_t
!			do k=1,nelec/2-1
!				Gama_B(i,j) = Gama_B(i,j) + 2.0d0*eigen_c_occ_local(i,k)*eigen_c_occ_local(j,k)
!			end do
!		end do
!	end do

!-------------------- P = S * Gama_B * S--------------------------


!	Project_new_temp = 0.0d0
!call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_M,norbit_t,Gama_B,norbit_t,beta1,Project_new_temp,norbit_t)

!	Project_new = 0.0d0
!call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,Project_new_temp,norbit_t,S_M,norbit_t,beta1,Project_new,norbit_t)
!!!	Project_new_temp = matmul(S_M_T,Gama_B)
!!!	Project_new = matmul(Project_new_temp,S_M)

!	write(2,*) 'Projection matrix for localized orbitals except the ones on fragment:'
!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(2,"(f16.6)", advance = "No") Project_new(i,j)
!		end do
!		write(2,*)
!	end do

!-------------------- New Fock--------------------------

!	F_A = F + mu*Project_new

!----------------------Transform and Diagonalize---------

!----------------------G = F_A*X----------------------------

!	G = 0.0d0
!  call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,F_A,norbit_t,X,norbit_t,beta1,G,norbit_t)

!----------------------F_Prime = XT*G = XT*F_A*X--------------------

!	F_Prime = 0.0d0
!  call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,G,norbit_t,beta1,F_Prime,norbit_t)

!	Eigen_F_A_prime = F_Prime

!----------------------Eigenvalues and Eigenvectors of F_Prime---------

!	E_value_F_A = 0.0d0
!	call DSYEV('V', 'U', norbit_t, Eigen_F_A_prime, norbit_t, E_value_F_A, WORK, LWORK, INFO)

!---------------------- X*C_Prime = C---------------------------------
!	Eigen_F_A = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,Eigen_F_A_prime,norbit_t,beta1,Eigen_F_A,norbit_t)

!	write(2,*) "Diagonalized Fock(F_A) Matrix:"
!	do i=1,norbit_t
!		do j=1,norbit_t
!			if (i.eq.j) then
!				write(2,"(f16.6)", advance = "No") E_value_F_A(i)
!			else
!				write(2,"(f16.6)", advance = "No") 0.0d0
!			end if
!		end do
!		write(2,*)
!	end do

!	write(2,*) "Fock(F_A) Eigen vectors:"
!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(2,"(f16.6)", advance = "No") Eigen_F_A(i,j)
!		end do
!		write(2,*)
!	end do

!============================================================
!								   CHECK THE PROGRAM II
!============================================================


!	P_check = 0.0d0

!----------------- Tr(PS) = num of e on site A----------------------
! P = density matrix constructed from eigenvectors of F_A corresponding to negative e-value
!-------------------------------------------------------------------

!	do i=1,norbit_t
!		do j=1,norbit_t
!			do k=1,nelec/2
!				if (k .eq. 1) then
!					P_check(i,j) = P_check(i,j) + 2*Eigen_F_A(i,k)*Eigen_F_A(j,k)
!				end if
!			end do
!		end do
!	end do

!	PS_check = 0.0d0
!call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_check,norbit_t,S_M,norbit_t,beta1,PS_check,norbit_t)

!!!		PS_check = matmul(P_check,S_M)

!	Tr_check = 0.0d0
!	do i=1,norbit_t
!		Tr_check = Tr_check + PS_check(i,i)
!	end do
!	write(2,*) "Trace(P*S):"
!	write(2,"(f16.8)") Tr_check

!	write(2,*) "F_A Density Matrix corresponding to the negative eigenvalue"
!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(2, "(f16.8)", advance = "No") P_check(i,j)
!		end do
!		write(2,"(f16.8)")
!	end do


!==========================================================

!==========================================================
!							Garnet Bath Orbital Construction
!==========================================================

	P_garnet = 0.0d0

	do i=1,norbit_t
		do j=1,norbit_t
			do k=1,nelec/2
				P_garnet(i,j) = P_garnet(i,j) + 2.0d0*eigen_c_occ_local(i,k)*eigen_c_occ_local(j,k)
			end do
		end do
	end do

!--------------------------------------------
!			Change the basis to natural orbitals
!							S1/2 * D * S1/2
!--------------------------------------------

	P_garnet_temp = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_garnet,norbit_t,S_half,norbit_t,beta1,P_garnet_temp,norbit_t)

	P_garnet_nat = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_half,norbit_t,P_garnet_temp,norbit_t,beta1,P_garnet_nat,norbit_t)

	write(*,*) "Density calculated from localized orbitals in AO basis:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(*, "(f16.8)", advance = "No") P_garnet(i,j)
		end do
		write(*,"(f16.8)")
	end do

	!write(2,*) "S^1/2*D*S^1/2 (natural orbital basis):"
	!do i=1,norbit_t
	!	do j=1,norbit_t
	!	write(2, "(f16.8)", advance = "No") P_garnet_nat(i,j)
	!end do
	!write(2,"(f16.8)")
	!end do

!--------------------------------------------
!					Check the basis change
!--------------------------------------------

!	do i=1,norbit_t
!		do j=1,norbit_t
!			eigen_C_nat(i,j) = P_garnet_nat(i,j)
!		end do
!	end do

!	P_garnet_evalue = 0.0d0
!	call DSYEV('V', 'U', norbit_t, eigen_C_nat, norbit_t, P_garnet_evalue, WORK, LWORK, INFO)

!	write(*,*) " (natural orbital basis):"
!	do i=1,norbit_t
!		do j=1,norbit_t
!		write(*, "(f16.8)", advance = "No") eigen_C_nat(i,j)
!	end do
!	write(*,"(f16.8)")
!	end do

!--------------SWAP the eigenvectors to be in the desired order-----------------

!	do j=1,int(norbit_t/2)
!		do i=1,norbit_t
!			call swap(eigen_C_nat(i,j),eigen_C_nat(i,norbit_t-j+1))
!		end do
!	end do

!	write(*,*) " (natural orbital basis):"
!	do i=1,norbit_t
!		do j=1,norbit_t
!		write(*, "(f16.8)", advance = "No") eigen_C_nat(i,j)
!	end do
!	write(*,"(f16.8)")
!	end do

!------------Multiply by S^-0.5 to get in AO basis----------------

!	eigen_c_nat_ao = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_half_inv,norbit_t,eigen_C_nat,norbit_t,beta1,eigen_c_nat_ao,norbit_t)

!	PS_check = 0.0d0
!	do i=1,norbit_t
!		do j=1,norbit_t
!			do k=1,nelec/2
!				PS_check(i,j) = PS_check(i,j) + 2.0d0*eigen_c_nat_ao(i,k)*eigen_c_nat_ao(j,k)
!			end do
!		end do
!	end do


!----------Zero out row/column------------------------------

	eigen_C_nat = 0.0d0

	do i=1,norbit_t
		do j=1,norbit_t
			if (i>norbit .and. j>norbit) then
 				eigen_C_nat(i,j) = P_garnet_nat(i,j)
			end if
		end do
	end do

	write(2,*) "P_garnet_natural(n-norbitxn-norbit):"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") eigen_C_nat(i,j)
		end do
		write(2,"(f16.8)")
	end do

	P_garnet_evalue = 0.0d0
	call DSYEV('V', 'U', norbit_t, eigen_C_nat, norbit_t, P_garnet_evalue, WORK, LWORK, INFO)

	countenv = 0.0d0
	write(2,*) "Diagonalized Garnet's Density matrix:"
	do i=1,norbit_t
		do j=1,norbit_t
			if (i==j) then
				write(2, "(f16.8)", advance = "No") P_garnet_evalue(i)

				if(dabs(P_garnet_evalue(i)-2.0d0)<=1.0d-8 ) then
						!.or. dabs(P_garnet_evalue(i)-0.0d0)<=1.0d-8
						countenv = countenv + 1
				else if (dabs(P_garnet_evalue(i)-0.0d0)>=1.0d-8) then
					bath = bath + 1
				end if

			else
				write(2,"(f16.8)",advance = "No") 0.0d0
			end if
		end do
		write(2,"(f16.8)")
	end do

!--------------SWAP the eigenvectors to be in the desired order-----------------

	do j=1,int(norbit_t/2)
		do i=1,norbit_t
			call swap(eigen_C_nat(i,j),eigen_C_nat(i,norbit_t-j+1))
		end do
	end do

	write(2,*) "Natural orbital eigenvectors:"
	do i=1,norbit_t
		do j=1,norbit_t
		write(2, "(f16.8)", advance = "No") eigen_C_nat(i,j)
	end do
	write(2,"(f16.8)")
	end do

!------------Multiply by S^-0.5 to get in AO basis----------------

	eigen_c_nat_ao = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_half_inv,norbit_t,eigen_C_nat,norbit_t,beta1,eigen_c_nat_ao,norbit_t)

!	write(2,*) "Natural orbital eigenvectors in AO basis:"
!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(2, "(f16.8)", advance = "No") eigen_C_nat_ao(i,j)
!		end do
!		write(2,"(f16.8)")
!	end do

!----------------Computing Gamma_B-------------------------------

	Gamma_garnet = 0.0d0
	do i=1,norbit_t
		do j=1,norbit_t
			do k=1,countenv
				!if(dabs(P_garnet_evalue(k)-2.0d0)<=3.0d-1 ) then
					!.or. dabs(P_garnet_evalue(k)-1.0d0)<=1.0d-8
					Gamma_garnet(i,j) = Gamma_garnet(i,j) + 2.0d0*eigen_C_nat_ao(i,k)*eigen_C_nat_ao(j,k)
				!end if
			end do
		end do
	end do
	write(*,*) "unentangled environment=",countenv
	write(*,*) "bath orbitals=",bath

	do i=1,norbit_t
		write(*, "(f16.8)") P_garnet_evalue(i)
	end do

	write(2,*) "Gamma_B:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") Gamma_garnet(i,j)
		end do
		write(2,"(f16.8)")
	end do

	PS_check = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,Gamma_garnet,norbit_t,S_M,norbit_t,beta1,PS_check,norbit_t)

	Tr_check = 0.0d0
	do i=1,norbit_t
		Tr_check = Tr_check + PS_check(i,i)
	end do
	write(2,*) 'Trace(Gamma_B*S) = ', Tr_check
!-------------------- Project_garnet = S * Gama_garnet * S--------------------------

	Project_garnet_temp = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,Gamma_garnet,norbit_t,S_M,norbit_t,beta1,Project_garnet_temp,norbit_t)

	Project_garnet = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_M,norbit_t,Project_garnet_temp,norbit_t,beta1,Project_garnet,norbit_t)


	!write(*,*) "Garnet's Projection matrix:"
	!do i=1,norbit_t
	!	do j=1,norbit_t
	!		write(*, "(f16.8)", advance = "No") Project_garnet(i,j)
	!	end do
	!	write(*,"(f16.8)")
	!end do

!==========================================================
!								New Fock Matrix Using Project_garnet
!==========================================================

!-------------------- New Fock--------------------------

	F_A = F + mu*Project_garnet

!----------------------Transform and Diagonalize---------

!----------------------G = F_A*X----------------------------

	G = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,F_A,norbit_t,X,norbit_t,beta1,G,norbit_t)

!----------------------F_Prime = XT*G = XT*F_A*X--------------------

	F_A_Prime = 0.0d0
	call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,G,norbit_t,beta1,F_A_Prime,norbit_t)

	Eigen_F_A_prime = F_A_Prime

!----------------------Eigenvalues and Eigenvectors of F_Prime---------

	E_value_F_A = 0.0d0
	call DSYEV('V', 'U', norbit_t, Eigen_F_A_prime, norbit_t, E_value_F_A, WORK, LWORK, INFO)

!---------------------- X*C_Prime = C---------------------------------
	Eigen_F_A = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,Eigen_F_A_prime,norbit_t,beta1,Eigen_F_A,norbit_t)

	write(2,*) "Diagonalized Fock(F_A) Matrix:"
	do i=1,norbit_t
		do j=1,norbit_t
			if (i.eq.j) then
				write(2,"(f16.6)", advance = "No") E_value_F_A(i)
			else
				write(2,"(f16.6)", advance = "No") 0.0d0
			end if
		end do
		write(2,*)
	end do

	write(2,*) "Fock(F_A) Eigen vectors:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2,"(f16.6)", advance = "No") Eigen_F_A(i,j)
		end do
		write(2,*)
	end do


	!Eigen_F_A = matmul(S_half,Eigen_F_A)
	!Eigen_F_A = matmul(Eigen_F_A,S_half)

!-------------------- Check--------------------------

	PS_check = 0.0d0
	do i=1,norbit_t
		do j=1,norbit_t
			do k=1,nelec_A/2
				if (E_value_F_A(k)<-1.0d-8) then
					PS_check(i,j) = PS_check(i,j) + 2.0d0*Eigen_F_A(i,k)*Eigen_F_A(j,k)
				end if
			end do
		end do
	end do

	write(2,*) "P_A:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") PS_check(i,j)
		end do
		write(2,"(f16.8)")
	end do
	write(2,"(f16.8)")

	write(2,*) "Gamma_B:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") Gamma_garnet(i,j)
		end do
		write(2,"(f16.8)")
	end do
	write(2,"(f16.8)")

	PS_check = PS_check + Gamma_garnet

	write(2,*) "P_A+Gamma_B:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") PS_check(i,j)
		end do
		write(2,"(f16.8)")
	end do
	write(2,"(f16.8)")

	write(2,*) "P_HF:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") P_new(i,j)
		end do
		write(2,"(f16.8)")
	end do


!---------------------high-level calculation (fci)-----------------

	!norbit_fci = norbit
	call fci_2_in_n

!---------------------Mulliken Population Analysis-----------------


					!	  write(7,*) EN_t

					!	  if (a .lt. max_bondlength) then

					!	    close(4,Status = 'delete')
					!	    a = a + bond_incr
					!	    EN = 0.0d0
					!	    EN_t = 0.0d0
					!	    write(2,*) "****************************************"
					!	    write(2,*) "  ************************************  "
					!	    go to 30

					!	  end if




!====================================================
! 			Matching the elements simultaneously
!----------------------------------------------------

			!do i=3,norbit_t
			!	do j=3,norbit_t
 			!		P_fci_A(i,j) = P_fci_A(i-norbit,j-norbit)
			!	end do
			!end do

			write(2,*) "High-level Density matrix"
			do i=1,norbit_t
				do j=1,norbit_t
					write(2, "(f16.8)", advance = "No") P_fci_A(i,j)
				end do
				write(2,*)
			end do


!====================================================
!          Intrinsic Atomic Orbitals (IAO)
!====================================================


!	S1_inv = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_half_inv,norbit_t,S_half_inv,norbit_t,beta1,S1_inv,norbit_t)

!	S2_inv = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_half_inv,norbit_t,S_half_inv,norbit_t,beta1,S2_inv,norbit_t)

!	S21 = Transpose(S12)

!	P12 = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S1_inv,norbit_t,S12,norbit_t,beta1,P12,norbit_t)

!	S1_inv_temp = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S1_inv,norbit_t,S12,norbit_t,beta1,S1_inv_temp,norbit_t)

!	S2_inv_temp = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S2_inv,norbit_t,S21,norbit_t,beta1,S2_inv_temp,norbit_t)

!	S_temp = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S1_inv_temp,norbit_t,S2_inv_temp,norbit_t,beta1,S_temp,norbit_t)

!	Eigen_C_temp = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_temp,norbit_t,Eigen_C,norbit_t,beta1,Eigen_C_temp,norbit_t)

!	Eigen_C_temp1 = 0.0d0
!	call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C_temp,norbit_t,S1,norbit_t,beta1,Eigen_C_temp1,norbit_t)

!	Eigen_C_coeff = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C_temp1,norbit_t,Eigen_C_temp,norbit_t,beta1,Eigen_C_coeff,norbit_t)

!	Eigen_C_tilde = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C,norbit_t,Eigen_C_coeff,norbit_t,beta1,Eigen_C_tilde,norbit_t)

!===============================
!	A1_temp = 0.0d0
!	call DGEMM('N','T',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C,norbit_t,Eigen_C,norbit_t,beta1,A1_temp,norbit_t)

!	A1 = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,A1_temp,norbit_t,S1,norbit_t,beta1,A1,norbit_t)

!	A2_temp = 0.0d0
!	call DGEMM('N','T',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C_tilde,norbit_t,Eigen_C_tilde,norbit_t,beta1,A2_temp,norbit_t)

!	A2 = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,A2_temp,norbit_t,S1,norbit_t,beta1,A2,norbit_t)

!	A3_temp = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,A1,norbit_t,A2,norbit_t,beta1,A3_temp,norbit_t)

!	A3 = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,A3_temp,norbit_t,P12,norbit_t,beta1,A3,norbit_t)

!	A1 = 1.0d0 - A1
!	A2 = 1.0d0 - A2

!	A1_temp = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,A1,norbit_t,A2,norbit_t,beta1,A1_temp,norbit_t)

!	A4 = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,A1_temp,norbit_t,P12,norbit_t,beta1,A4,norbit_t)

!	rho = A3 + A4

!	rho_temp = 0.0d0
!	call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,rho,norbit_t,S1,norbit_t,beta1,rho_temp,norbit_t)

!	rho_coeff = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,rho_temp,norbit_t,rho,norbit_t,beta1,rho_coeff,norbit_t)

!	rho_orth = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,rho,norbit_t,rho_temp1,norbit_t,beta1,rho_orth,norbit_t)


!===========================================

!	write(*,*) "Orthognalized S:"
!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(*, "(f16.8)", advance = "No") S_orth(i,j)
!		end do
!		write(*,"(f16.8)")
!	end do


!	S_half_temp = 0.0d0
!	call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C,norbit_t,S_half_inv,norbit_t,beta1,S_half_temp,norbit_t)
!	S_half_orth = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_half_temp,norbit_t,Eigen_C,norbit_t,beta1,S_half_orth,norbit_t)

!	write(*,*) "S^-1/2:"
!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(*, "(f16.8)", advance = "No") S_half_inv(i,j)
!		end do
!		write(*,"(f16.8)")
!	end do

!	Eigen_C_orth = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C,norbit_t,S_half_inv,norbit_t,beta1,Eigen_C_orth,norbit_t)

!	write(*,*) "Orthognalized Eigenvectors:"
!	do i=1,norbit_t
!		do j=1,norbit_t
!			write(*, "(f16.8)", advance = "No") Eigen_C_prime(i,j)
!		end do
!		write(*,"(f16.8)")
!	end do

!-------------XT * S * D * S * X---------------------------------------

	P_temp = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_M,norbit_t,P_new,norbit_t,beta1,P_temp,norbit_t)

	P_orth = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_temp,norbit_t,S_M,norbit_t,beta1,P_orth,norbit_t)

	P_temp = 0.0d0
	call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,P_orth,norbit_t,beta1,P_temp,norbit_t)

	P_orth = 0.0d0
	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_temp,norbit_t,X,norbit_t,beta1,P_orth,norbit_t)

!	do j=1,int(norbit_t/2)
!		do i=1,norbit_t
!			call swap(P_orth(i,j),P_orth(i,norbit_t-j+1))
!		end do
!	end do
!	do i=1,int(norbit_t/2)
!    do j=1,norbit_t
!      call swap(P_orth(i,j),P_orth(norbit_t-i+1,j))
!    end do
!  end do
!-----------------------------------------------------------------------


!	P_temp = 0.0d0
!	call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C,norbit_t,P_orth,norbit_t,beta1,P_temp,norbit_t)

!	P_orth = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_temp,norbit_t,Eigen_C,norbit_t,beta1,P_orth,norbit_t)


!Eigen_C_prime = 0.0d0
!call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_half_inv,norbit_t,Eigen_C,norbit_t,beta1,Eigen_C_prime,norbit_t)

	P_new_orth = 0.0d0
	do i=1,norbit_t
		do j=1,norbit_t
			do k=1,nelec/2
				P_new_orth(i,j) = P_new_orth(i,j) + 2.0d0*Eigen_C_prime(i,k)*Eigen_C_prime(j,k)
			end do
		end do
	end do

!	P_temp = 0.0d0
!	call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,Eigen_C,norbit_t,P_new_orth,norbit_t,beta1,P_temp,norbit_t)

!	P_new_orth = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_temp,norbit_t,Eigen_C,norbit_t,beta1,P_new_orth,norbit_t)


	write(2,*) "Orthognalized HF density USING XT*S*P*S*X:"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") P_orth(i,j)
		end do
		write(2,"(f16.8)")
	end do

	write(2,*) "Orthognalized HF density USING C':"
	do i=1,norbit_t
		do j=1,norbit_t
			write(2, "(f16.8)", advance = "No") P_new_orth(i,j)
		end do
		write(2,"(f16.8)")
	end do

!	P_temp = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_M,norbit_t,P_fci_A,norbit_t,beta1,P_temp,norbit_t)

!	P_orth = 0.0d0
!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_temp,norbit_t,S_M,norbit_t,beta1,P_orth,norbit_t)

!!	P_temp = 0.0d0
!!	call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,P_fci_A_mo,norbit_t,beta1,P_temp,norbit_t)

!!	P_orth = 0.0d0
!!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_temp,norbit_t,X,norbit_t,beta1,P_orth,norbit_t)

!!	write(*,*) "FCI density:"
!!	do i=1,norbit_t
!!		do j=1,norbit_t
!!			write(*, "(f16.8)", advance = "No") P_orth(i,j)
!!		end do
!!		write(*,"(f16.8)")
!!	end do


	call bfgs_oep

		100 continue

!==============================================================
!						Calculating the Final Energy
!==============================================================

		PE_t_final = PE_t_new

		do i=1,norbit
			do j=1,norbit
				PE_t_final(i+norbit,j+norbit)= PE_t_new(i,j)
			end do
		end do

	!	PE_t_final = PE_t_new
		F_Prime = KE_t_prime + PE_t_final


		write(2,*) "Initial Potential Matrix:"
		do i=1,norbit_t
			do j=1,norbit_t
				write(2, "(f16.8)", advance = "No") PE_t_init(i,j)
			end do
			write(2,"(f16.8)")
		end do

		write(2,*) "Final Potential Matrix after OEP:"
		do i=1,norbit_t
			do j=1,norbit_t
				write(2, "(f16.8)", advance = "No") PE_t_new(i,j)
			end do
			write(2,"(f16.8)")
		end do

		write(2,*) "Final Potential Matrix:"
		do i=1,norbit_t
			do j=1,norbit_t
				write(2, "(f16.8)", advance = "No") PE_t_final(i,j)
			end do
			write(2,"(f16.8)")
		end do
!----------------------Transform and Diagonalize---------

!----------------------G = F*X----------------------------

!		G = 0.0d0
!		call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X_inv,norbit_t,F_Prime,norbit_t,beta1,G,norbit_t)

!----------------------F_Prime = XT*G = XT*F*X--------------------

!		F = 0.0d0
!		call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,G,norbit_t,X_inv,norbit_t,beta1,F,norbit_t)

		do i=1,norbit_t
			do j=1,norbit_t
				Eigen_C_prime(i,j) = F_Prime(i,j)
			end do
		end do

	!	write(2,*) "Transformed Fock Matrix"
	!	do i=1,norbit_t
	!		do j=1,norbit_t
	!			write(2, "(f16.8)", advance = "No") F_Prime(i,j)
	!		end do
	!		write(2,"(f16.8)")
	!	end do

!----------------------Eigenvalues and Eigenvectors of F_Prime---------

		Eigen_value = 0.0d0
		call DSYEV('V', 'U', norbit_t, Eigen_C_prime, norbit_t, Eigen_value, WORK, LWORK, INFO)

!---------------------- X*C_Prime = C---------------------------------

	!	Eigen_C = 0.0d0
	!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,Eigen_C_prime,norbit_t,beta1,Eigen_C,norbit_t)

		P_orth = 0.0d0
		do i=1,norbit_t
			do j=1,norbit_t
				do k=1,nelec/2
					P_orth(i,j) = P_orth(i,j) + 2.0d0*Eigen_C_prime(i,k)*Eigen_C_prime(j,k)
				end do
			end do
		end do

!		do i=3,norbit_t
!			do j=3,norbit_t
!				P_orth(i,j) = P_orth(i-norbit,j-norbit)
!			end do
!		end do

		write(2,*) "Final Density Matrix:"
		do i=1,norbit_t
			do j=1,norbit_t
				write(2, "(f16.8)", advance = "No") P_orth(i,j)
			end do
			write(2,"(f16.8)")
		end do

		PS_check = 0.0d0
		call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_new,norbit_t,S_M,norbit_t,beta1,PS_check,norbit_t)



		Tr_check = 0.0d0
		do i=1,norbit_t
			Tr_check = Tr_check + P_orth(i,i)
		end do
		write(*,*) "FINAL Trace(P*S):"
		write(*,"(f16.8)") Tr_check
!-----------------------Electronic Energy----------------


		G = 0.0d0
		call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,H_core,norbit_t,X,norbit_t,beta1,G,norbit_t)

!----------------------H_core_orth = XT*H_core*X--------------------

		H_core_orth = 0.0d0
		call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,G,norbit_t,beta1,H_core_orth,norbit_t)

		EN = 0.0d0
		do i=1,norbit_t
		  do j=1,norbit_t
		    EN = EN + 0.5d0*P_orth(i,j)*(H_core_orth(i,j)+F_Prime(i,j))
		  end do
		end do


		EN_t = EN
		do i=1,natom-1
			do j=i+1,natom
				call Distance(R(i,1),R(i,2),R(i,3),R(j,1),R(j,2),R(j,3),d)
				EN_t = EN_t + Z(i)*Z(j)/d
			end do
		end do

		open(Unit = 7,File = 'Energy')
		write(7,*) "DMET energy:",a, EN_t


		write(2,*) "DMET Electronic Energy = ",EN
		write(*,*) "DMET Total Energy per Atom = ",EN_t/natom
		write(*,*) "========================================="
		!write(*,*) "Congratulations! DMET Calculation Converged"

		write(*,*) "Orthogonalized AO HF Density Matrix:"
		do i=1,norbit_t
		  do j=1,norbit_t
		    write(*, "(f16.8)", advance = "No") P_orth(i,j)
		  end do
		write(*,*)
		end do


		write(*,*) "Orthogonalized AO FCI Density matrix:"
		do i=1,norbit_t
			do j=1,norbit_t
				write(*, "(f16.8)", advance = "No") P_fci_A(i,j)
			end do
			write(*,*)
		end do
		write(*,*) "****************************************"
		write(*,*) "  ************************************  "


  End Program

!====================================================
!                 END OF THE MAIN PROGRAM
!====================================================
