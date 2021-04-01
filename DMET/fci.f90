  Subroutine fci_2_in_n
    USE array

!====================================================
!                 CI CALCULATION
!====================================================

!=====================================================
! CISD calculation is done based on the book,
! "modern quantum chemistry" by Szabu and Ostland.
! to better undersand the notation please refer to chapters
! 3 & 4 of this book.
!=====================================================

		write(*,*) "CI Calculation Started"
		write(2,*) "======================================"

		write(2,*) "------------CI Calculations-------------"

!====================================================
!              H_core in MO basis(F_A: fragment)
!----------------------------------------------------

	H_core_new = 0.0d0
	do i=1,norbit_t
	  do j=1,norbit_t
			do k=1,norbit_t
				do l=1,norbit_t
		      H_core_new(i,j) = H_core_new(i,j) + H_core(k,l)*Eigen_F_A(k,i)*Eigen_F_A(l,j)
				end do
			end do
		end do
	end do

 write(*,*) "Hello1"
!====================================================
!           2e-INTEGRALS in MO basis(F_A: fragment)
!----------------------------------------------------

M_2e_new = 0.0d0
	do i=1,norbit_t
		do j=1,norbit_t
			do k=1,norbit_t
		    do l=1,norbit_t
!-----------------------------------------------------
					do m=1,norbit_t
					  do n=1,norbit_t
							do o=1,norbit_t
						  	do p=1,norbit_t
							  	M_2e_new(i,j,k,l) = M_2e_new(i,j,k,l) + M_2e(m,n,o,p)*Eigen_F_A(m,i)*Eigen_F_A(n,j)*Eigen_F_A(o,k)*Eigen_F_A(p,l)
					      end do
							end do
		  			end do
					end do
        end do
      end do
    end do
  end do

write(*,*) "Hello2"
!====================================================
!     CALCULATING Full CI: Hamiltonian
!             2e in n orbitals
!----------------------------------------------------

!----------------------------------------------------
!								S1 = Singlet Single excitation
!								D1 =  (a,a----r,r) Singlet
!								D2 = (a,a----r,s)  Singlet
!								D3 = (a,a----r,s)  Triplet
!----------------------------------------------------
	H_t = 0.0d0
	counti = 1
	countj = 1

!----------------------------------------------------
!										Ground * Ground
!----------------------------------------------------


	H_t(counti,countj) = 2.0d0*H_core_new(1,1) + M_2e_new(1,1,1,1)

!----------------------------------------------------

	counti = 2
	countj = 2
!------------Ground * Single Excitations-------------
!											ZERO
!----------------------------------------------------


!----------------------------------------------------
!								     S1 * S1
!----------------------------------------------------
	do i=2,norbit_fci
		do j=2,norbit_fci
			if (i==j) then
				H_t(counti,countj) = H_core_new(1,1) + H_core_new(i,i) + M_2e_new(1,1,i,i) + M_2e_new(1,i,i,1)
			else

				H_t(counti,countj) = H_core_new(i,j) + M_2e_new(i,j,1,1) + M_2e_new(1,j,i,1)
			end if
		countj = countj+1
		end do
		countj = 2

		counti = counti+1

	end do

	!------------------------------------------------------
	!											S1 * D1
	!------------------------------------------------------


		counti = 2
		countj = norbit_fci+1

		do i=2,norbit_fci
				do j=2,norbit_fci
					if (i==j) then
						H_t(counti,countj) =  dsqrt(2.0d0)*(H_core_new(1,i)+M_2e_new(1,i,i,i))
					else
						H_t(counti,countj) =  dsqrt(2.0d0)*M_2e_new(i,j,1,j)
					end if

						H_t(countj,counti) = H_t(counti,countj)
						countj = countj + 1

					end do
	!--------------------------------------------------------
			counti = counti + 1
			countj = norbit_fci+1

		end do

!------------------------------------------------------
!											S1 * D2
!------------------------------------------------------


	counti = 2
	countj = 2*norbit_fci

	do i=2,norbit_fci
			do j=2,norbit_fci-1
				do k=j+1,norbit_fci

!------------ (a----r)*(a,a----r',s') ---------------

					if (i/=j .and. i/=k) then
						H_t(counti,countj) =  M_2e_new(1,j,i,k) + M_2e_new(1,k,i,j)

!------------ (a----r)*(a,a----r,s) ---------------

				else if(i==j) then
						H_t(counti,countj) =  H_core_new(1,k) +M_2e_new(1,k,i,i)+M_2e_new(i,k,1,i)

					else if(i==k) then
						H_t(counti,countj) =  H_core_new(1,j) +M_2e_new(1,j,i,i)+M_2e_new(i,j,1,i)

					end if

					H_t(countj,counti) = H_t(counti,countj)
					countj = countj + 1

				end do
!--------------------------------------------------------
		end do
		counti = counti + 1
		countj = 2*norbit_fci

	end do


!------------Ground * Double Excitations-------------
!----------------------------------------------------
!								     Ground * D1
!----------------------------------------------------

!------------ (a,a----r,r) --------------------------
counti = 1
countj = norbit_fci+1

		do i=2,norbit_fci
			H_t(1,countj) = M_2e_new(i,1,1,i)
			H_t(countj,1) = H_t(1,countj)

			countj = countj + 1
		end do
!----------------------------------------------------
!											Ground * D2
!----------------------------------------------------
!------------ (a,a----r,s) --------------------------

		do i=2,norbit_fci-1
			do j=i+1,norbit_fci

				H_t(1,countj) = dsqrt(2.0d0)*M_2e_new(1,i,1,j)

				H_t(countj,1) = H_t(1,countj)
				countj = countj + 1

			end do
		end do

!------------------------------------------------------
!												D1 * D1
!------------------------------------------------------

	counti = norbit_fci+1
	countj = norbit_fci+1

	do i=2,norbit_fci
		do j=2,norbit_fci

!------------ Diagonal Elements ---------------------

			if (i==j) then

			H_t(counti,countj) = 2.0d0*H_core_new(i,i) + M_2e_new(i,i,i,i)

!------------ (a,a----r,r)*(a,a----r',r') -------------

			else if (i/=j) then
				H_t(counti,countj) =  M_2e_new(i,j,i,j)
			end if

			countj = countj + 1
		end do
		counti = counti + 1
		countj = norbit_fci+1
	end do


!----------------------------------------------------
!											D2 * D2
!----------------------------------------------------


	counti = 2*norbit_fci
	countj = 2*norbit_fci

	do i=2,norbit_fci-1
		do j=i+1,norbit_fci

			do k=2,norbit_fci-1
				do l=k+1,norbit_fci

!------------ Diagonal Elements ---------------------

					if(i==k .and. j==l) then

H_t(counti,countj) = H_core_new(i,i)+H_core_new(j,j) +M_2e_new(i,i,j,j)+M_2e_new(i,j,j,i)


!------------ (a,a----r,s)*(a,a----r,s') ------------

					else	if (i==k .and. j/=l) then

						H_t(counti,countj) = H_core_new(j,l)+M_2e_new(j,l,i,i)+M_2e_new(i,l,j,i)
					else	if (i==l .and. j/=k) then

						H_t(counti,countj) = H_core_new(j,k)+M_2e_new(j,k,i,i)+M_2e_new(i,k,j,i)
					else	if (i/=k .and. j==l) then

						H_t(counti,countj) = H_core_new(i,k)+M_2e_new(i,k,j,j)+M_2e_new(j,k,i,j)
					else	if (i/=l .and. j==k) then

						H_t(counti,countj) = H_core_new(i,l)+M_2e_new(i,l,j,j)+M_2e_new(j,l,i,j)
!------------ (a,a----r,s)*(a,a----r',s') -----------

					else if (i/=k .and. j/=l .and. i/=l .and. j/=k) then
						H_t(counti,countj) = M_2e_new(i,k,j,l)+M_2e_new(i,l,j,k)
					end if
!----------------------------------------------------

					countj = countj + 1

				end do
			end do
			countj = 2*norbit_fci
			counti = counti + 1

		end do
	end do



	!------------------------------------------------------
	!											D1 * D2
	!------------------------------------------------------


		counti = norbit_fci+1
		countj = 2*norbit_fci

		do i=2,norbit_fci
				do j=2,norbit_fci-1
					do k=j+1,norbit_fci

	!------------ (a,a----r,r)*(a,a----r',s') ---------------

						if (i/=j .and. i/=k) then
							H_t(counti,countj) =  dsqrt(2.0d0)*M_2e_new(i,j,i,k)

	!------------ (a,a----r,r)*(a,a----r,s) ---------------

						else if(i==j) then
							H_t(counti,countj) =  dsqrt(2.0d0)*(H_core_new(i,k) +M_2e_new(i,k,i,i))

						else if(i==k) then
							H_t(counti,countj) =  dsqrt(2.0d0)*(H_core_new(i,j) +M_2e_new(i,j,i,i))

						end if

						H_t(countj,counti) = H_t(counti,countj)
						countj = countj + 1

					end do
	!--------------------------------------------------------
			end do
			counti = counti + 1
			countj = 2*norbit_fci

		end do

	!----------------------------------------------------
	!											D3 * D3
	!----------------------------------------------------


	!	counti = norbit_fci + (norbit_fci-1)*(norbit_fci-2)/2 +1
	!	countj = norbit_fci + (norbit_fci-1)*(norbit_fci-2)/2 +1

	!	do i=2,norbit_fci-1
	!		do j=i+1,norbit_fci

	!			do k=2,norbit_fci-1
	!				do l=k+1,norbit_fci

	!------------ Diagonal Elements ---------------------

	!					if(i==k .and. j==l) then

!	H_t(counti,countj) = (Eigen_value(i)+Eigen_value(j)) - 2.0d0*Eigen_value(1) + M_2e_new(1,1,1,1) + M_2e_new(i,i,j,j) &
!				+ M_2e_new(i,j,j,i) - 2.0d0*M_2e_new(j,j,1,1)- 2.0d0*M_2e_new(i,i,1,1) &
!				+ M_2e_new(i,1,1,i) + M_2e_new(j,1,1,j)
!H_t(counti,countj) = H_core_new(i,i)+H_core_new(j,j) + (M_2e_new(i,i,j,j)-M_2e_new(i,j,j,i))

!write(2,*) 'Remaining Diagonal:',counti,countj,H_t(counti,countj)



	!------------ (a,a----r,s)*(a,a----r,s') ------------

!						else	if (i==k .and. j/=l) then

!							H_t(counti,countj) = H_core_new(j,l)+M_2e_new(j,l,i,i)-M_2e_new(j,i,i,l)
!						else	if (i==l .and. j/=k) then

!						H_t(counti,countj) = H_core_new(j,k)+M_2e_new(j,k,i,i)-M_2e_new(j,i,i,k)
!						else	if (i/=k .and. j==l) then

!							H_t(counti,countj) = H_core_new(i,k)+M_2e_new(i,k,j,j)-M_2e_new(i,j,j,k)
!						else	if (i/=l .and. j==k) then

!							H_t(counti,countj) = H_core_new(i,l)+M_2e_new(i,l,j,j)-M_2e_new(i,j,j,l)
	!------------ (a,a----r,s)*(a,a----r',s') -----------

!						else if (i/=k .and. j/=l) then
!							H_t(counti,countj) = M_2e_new(i,k,j,l)-M_2e_new(i,l,j,k)
!						end if
	!----------------------------------------------------

!						countj = countj + 1

!					end do
!				end do
!				countj = norbit_fci + (norbit_fci-1)*(norbit_fci-2)/2 +1
!				counti = counti + 1

!			end do
!		end do




	!------------------------------------------------------
	!											D1 * D3
	!------------------------------------------------------


!		counti = 2
!		countj = 2*norbit_fci

!		do i=2,norbit_fci
!				do j=2,norbit_fci-1
!					do k=j+1,norbit_fci

!						if (i/=j .and. i/=k) then
!							H_t(counti,countj) = 0.0d0

	!------------ (a,a----r,r)*(a,a----r,s) ---------------

!						else if(i==k) then
!							H_t(counti,countj) =  dsqrt(2.0d0)*H_core_new(i,j)
!						else if(i==j) then
!							H_t(counti,countj) =  dsqrt(2.0d0)*H_core_new(i,k)
!						end if

!						H_t(countj,counti) = H_t(counti,countj)
!						countj = countj + 1

!					end do
	!--------------------------------------------------------
!			end do
!			counti = counti + 1
!			countj = 2*norbit_fci

!		end do

!----------------------------------------------------
!											D2 * D3
!----------------------------------------------------


!	counti = norbit_fci+1
!	countj = 2*norbit_fci

!	do i=2,norbit_fci-1
!		do j=i+1,norbit_fci
!
!			do k=2,norbit_fci-1
!				do l=k+1,norbit_fci

		!------------ (a,a----r,s)*(a,a----r,s) ------------

!					if (i==k .and. j==l) then

!						H_t(counti,countj) = H_core_new(i,i)+H_core_new(j,j)

!					else if(i==k .and. j/=l) then
!						H_t(counti,countj) = H_core_new(j,l)

!					else if(i==l .and. j/=k) then
!						H_t(counti,countj) = H_core_new(j,k)

!					else if(i/=k .and. j==l) then
!						H_t(counti,countj) = H_core_new(i,k)

!					else if(i/=l .and. j==k) then
!						H_t(counti,countj) = H_core_new(i,l)

!					end if
		!----------------------------------------------------
!		H_t(countj,counti) = H_t(counti,countj)
!		countj = countj + 1

!				end do
!			end do
!			countj = 2*norbit_fci
!			counti = counti + 1

!		end do
!	end do


!---------------------Diagonalizing H-----------------------

	E_CI = 0.0d0
	Eigen_C_H_t = H_t
	call DSYEV('V', 'U', mat_dim, Eigen_C_H_t, mat_dim, E_CI, WORK_FCI, LWORK_FCI, INFO)
!------------------------------------------------------------

!---------------------Total Energy---------------------------

	write(2,*) "Eigenvalues of FCI Hamiltonian: "
	do i=1,mat_dim
		write(2,*)	E_CI(i)
	end do
	write(2,*)

	E_corr = E_CI(1)-H_t(1,1)
  write(*,*) "H_t(1,1)=", H_t(1,1)
  write(*,*) "E_CI(1)=", E_CI(1)
	write(*,*) "total correlation energy:", (nelec/2)*E_corr
	write(*,*) "****************************************"
	write(*,*) "  ************************************  "

!---------------------FCI Eigenvectors------------------
	write(2,*) "Eigenvectors of FCI Hamiltonian:"
	do i=1,mat_dim
		do j=1,mat_dim
			write(2, "(f16.8)", advance = "No") Eigen_C_H_t(i,j)
		end do
		write(2,*)
	end do
	write(2,*) "****************************************"
	write(2,*) "  ************************************  "

	write(*,*) "RHF energy per atom:", EN_t
	write(*,*) "FCI energy per atom:", ((EN_t*natom)+(nelec/2)*E_corr)/natom

!---------------------fci eigenvectors Matrix------------------

	write(2,*) "Eigenvectors of FCI Hamiltonian:"
	do i=1,mat_dim
		do j=1,mat_dim
			write(2, "(f16.8)", advance = "No") Eigen_C_H_t(i,j)
		end do
		write(2,*)
	end do
	write(*,*) "CI Calculation is done on the fragment!"
	write(2,*) "----------------------------------------"

!-------------Density matrix using CI Vector---------
!     			P_fci_mo(p,q)	= <Di|ap aq + ap aq |Dj>
!---------------------------------------------------------

		P_fci_A = 0.0d0
		P_fci_A_mo = 0.0d0


!-----------------------------------------------------------------------
!										Generalized Construction
!-----------------------------------------------------------------------

		k = 1
		Slater_coeff =0.0d0
		do i=1,norbit_fci
			do j=i,norbit_fci
				if(i==j .and. i==2) then
					Slater_coeff(i,j) = 0.0d0
					k = k+norbit_fci-1
				else if(i/=j) then
					Slater_coeff(i,j) = Eigen_C_H_t(k,1)
					k = k+1
					Slater_coeff(j,i)=Slater_coeff(i,j)
				else if(i==j .and. i==1) then
					Slater_coeff(i,j) = Eigen_C_H_t(k,1)
					k = k+1
					Slater_coeff(j,i)=Slater_coeff(i,j)
				end if
			end do
		end do

		k =norbit_fci+1
		do i=2,norbit_fci
			Slater_coeff(i,i) = Eigen_C_H_t(k,1)
			k = k+1
		end do

		write(2,*) "Slater coefficients:"
		do i=1,norbit_fci
			do j=1,norbit_fci
				write(2, "(f16.8)", advance = "No") Slater_coeff(i,j)
			end do
			write(2,*)
		end do
		write(2,*) "****************************************"
		write(2,*) "  ************************************  "

		P_fci_A_mo =0.0d0

		do i=1,norbit_fci
			do j=i,norbit_fci

	!--------------------------------------------------------------------
			 	if (i==j) then
					do k=1,norbit_fci
						if (i/=k) then
							P_fci_A_mo(i,j) =P_fci_A_mo(i,j)+Slater_coeff(i,k)**2.0d0
						end if
						if (i==k) then
							P_fci_A_mo(i,j) =P_fci_A_mo(i,j)+ 2.0d0*Slater_coeff(i,k)**2.0d0
						end if
					end do
		!--------------------------------------------------------------------
				else
					P_fci_A_mo(i,j) = dsqrt(2.0d0)*(Slater_coeff(i,j)*Slater_coeff(i,i)+Slater_coeff(j,j)*Slater_coeff(i,j))

					do k=1,norbit_fci
						if (i/=k .and. j/=k) then
							P_fci_A_mo(i,j) = P_fci_A_mo(i,j) + Slater_coeff(k,j)*Slater_coeff(i,k)
						end if
					end do

					P_fci_A_mo(j,i)= P_fci_A_mo(i,j)
				end if
			end do
		end do

    !---------------------CI Vector------------------
    		write(2,*) "CI Vector:"
    		do i=1,mat_dim
    				write(2, "(f16.8)") Eigen_C_H_t(i,1)
    			end do
    		write(2,*) "****************************************"
    		write(2,*) "  ************************************  "



    Eigen_F_A_inv = Eigen_F_A
  !  Eigen_F_A_inv (1,1) = 1.0d0
  !  Eigen_F_A_inv (2,2) = 1.0d0
  !  Eigen_F_A_inv (3,3) = 1.0d0
  !  Eigen_F_A_inv (4,4) = 1.0d0

  call DGETRF(norbit_t, norbit_t, Eigen_F_A_inv, norbit_t, IPIV1, INFO1)
  call DGETRI(norbit_t, Eigen_F_A_inv, norbit_t, IPIV1, WORK1, LWORK1, INFO1)

!-------------Converting from MO to AO basis--------------

		P_fci_A_temp = 0.0d0
		call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,Eigen_F_A,norbit_t,P_fci_A_mo,norbit_t,beta1,P_fci_A_temp,norbit_t)

		P_fci_A = 0.0d0
		call DGEMM('N','T',norbit_t,norbit_t,norbit_t,alpha1,P_fci_A_temp,norbit_t,Eigen_F_A,norbit_t,beta1,P_fci_A,norbit_t)

  !  P_fci_A_temp = 0.0d0
	!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_inv,norbit_t,P_fci_A,norbit_t,beta1,P_fci_A_temp,norbit_t)

	!	P_fci_A = 0.0d0
	!	call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_fci_A_temp,norbit_t,S_inv,norbit_t,beta1,P_fci_A,norbit_t)

    PS_check = 0.0d0
		call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,X_inv,norbit_t,X,norbit_t,beta1,PS_check,norbit_t)

!---------------------Check the inversion----------
    write(*,*) "X*Xinv:"
    do i=1,norbit_t
      do j=1,norbit_t
        write(*, "(f16.8)", advance = "No") PS_check(i,j)
      end do
      write(*,*)
    end do
    write(*,*) "****************************************"
    write(*,*) "  ************************************  "
  !  write(*,*) "info=",INFO1

!---------------------fci Density Matrix AO basis----------
    write(*,*) "Non-orthogonalized AO FCI Density matrix:"
    do i=1,norbit_t
    	do j=1,norbit_t
    		write(*, "(f16.8)", advance = "No") P_fci_A(i,j)
    	end do
    	write(*,*)
    end do
    write(*,*) "****************************************"
  	write(*,*) "  ************************************  "

!-------------XT * S * D * S * X---------------------------------------

    P_fci_A_temp = 0.0d0
		call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_fci_A,norbit_t,S_M,norbit_t,beta1,P_fci_A_temp,norbit_t)

    P_fci_A = 0.0d0
    call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,S_M,norbit_t,P_fci_A_temp,norbit_t,beta1,P_fci_A,norbit_t)

    P_fci_A_temp = 0.0d0
    call DGEMM('T','N',norbit_t,norbit_t,norbit_t,alpha1,X,norbit_t,P_fci_A,norbit_t,beta1,P_fci_A_temp,norbit_t)

    P_fci_A = 0.0d0
    call DGEMM('N','N',norbit_t,norbit_t,norbit_t,alpha1,P_fci_A_temp,norbit_t,X,norbit_t,beta1,P_fci_A,norbit_t)

  !  do j=1,int(norbit_t/2)
  !  	do i=1,norbit_t
  !  		call swap(P_fci_A(i,j),P_fci_A(i,norbit_t-j+1))
  !  	end do
  !  end do
  !  do i=1,int(norbit_t/2)
  !    do j=1,norbit_t
  !      call swap(P_fci_A(i,j),P_fci_A(norbit_t-i+1,j))
  !    end do
  !  end do
!---------------------fci Density Matrix MO basis----------
		write(2,*) "FCI Density matrix MO basis on fragment A:"
		do i=1,norbit_t
			do j=1,norbit_t
				write(2, "(f16.8)", advance = "No") P_fci_A_mo(i,j)
			end do
			write(2,*)
		end do
		write(2,*) "****************************************"
		write(2,*) "  ************************************  "

!---------------------fci Density Matrix AO basis----------
		write(2,*) "Orthogonalized AO FCI Density matrix:"
		do i=1,norbit_t
			do j=1,norbit_t
				write(2, "(f16.8)", advance = "No") P_fci_A(i,j)
			end do
			write(2,*)
		end do
		write(2,*) "****************************************"
		write(2,*) "  ************************************  "

!---------------------HF Density Matrix-----------------

		write(2,*) "Non-orthogonalized total HF Density Matrix:"
		do i=1,norbit_t
			do j=1,norbit_t
				write(2, "(f16.8)", advance = "No") P_new(i,j)
			end do
			write(2,*)
		end do
		write(2,*) "****************************************"
		write(2,*) "  ************************************  "

!---------------------Overlap Matrix------------------
			!	write(2,*) "Overlap Matrix S"
			!	do i=1,norbit_t
			!	  do j=1,norbit_t
			!	    write(2, "(f16.8)", advance = "No") S_M(i,j)
			!	  end do
			!		write(2,*)
			!	end do


!---------------------Core Hamiltonian-----------------

					!	write(2,*) "H_core Matrix"
					!	do i=1,norbit_t
					!	  do j=1,norbit_t
					!	    write(2, "(f16.10)", advance = "No") H_core(i,j)
					!	  end do
					!	write(2,*)
					!	end do
!---------------------H_t Matrix------------------
					write(2,*) "H_t Matrix"
					do i=1,mat_dim
					  do j=1,mat_dim
			    		write(2, "(f16.8)", advance = "No") H_t(i,j)
					  end do
						write(2,*)
					end do
!---------------------FCI Orbital Energies------------------

				write(2,*) "Eigenvalues of FCI Hamiltonian on fragment A: "
				do i=1,mat_dim
					write(2,"(f16.8)")  E_CI(i)
				end do
				write(2,*)
!---------------------HF Orbital Energies------------------

					write(2,*) "HF Orbital Energies"
					do i=1,norbit_t
				    write(2, "(f16.8)") Eigen_value(i)
					end do
					write(2,*) "========================================="


    return
  end
