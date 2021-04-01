	Program FullCI

!====================================================
!                  SOME GENERAL COMMENTS
!----------------------------------------------------
! 1. All calculations were done in atomic units.
! 2. 1s primitive Gaussian functions were used.
! 3. The basis set is in a file named "Basis.txt".
! 4. The coordinates of atoms are in a file named "Input.xyz".
! 5. The final results for each bond length will be saved in a file named "Output.txt".
! 6. Total Energy per atom as a function of distance will be saved in a file named "Energy"

!====================================================

!====================================================
!                  VARIABLES(ALPHABETICAL ORDER)
!----------------------------------------------------
!
! a = side length of the polygon(bond length)
! alpha(i,j) = Guassian exponent for atom "i" and primitive "j"
! alpha1 = Scaling factor for DGEMM subroutine
! atomname = name of the atom
! beta1 = Scaling factor for DGEMM subroutine
! bond_incr = bond increment (for plotting total energy vs bond length)
! coeff(i,j) = Contraction coefficient for atom "i" and primitive "j"
! Converge = whenever (delta < Converge) this means that our result is in the desired interval
! d,d1,d2,d3 = distances between centers
! delta = standard deviation for change in density matrix
! Eigen_C/Eigen_C_prime = Eigen vectors of Fock/transformed Fock Matrix
! Eigen_value() = Eigenvalues of matrix
! EN = Electronic energy
! EN_t = total energy
! F(),F_Prime() = Fock Matrix before/after transformation
! G() = 2e part of Fock Matrix
! geo = if geo = 1 the geometry would be polygon, if geo = 2 the geometry will be grid
! H_core() = core H matrix
! Integral = output of 2e-integral
! itemp = temporary counter for total basis
! jtemp = temporary counter for total basis
! KE,KE_t = Kinetic energy(primitive)/total
! ktemp = temporary counter for total basis
! ltemp = temporary counter for total basis
! M_2e() = Total matrix regarding to 2e-integrals
! max_it = maximum iteration
! max_bondlength = maximum bond length for which we are going to calculate the total energy
! natom = number of atom
! nbas = number of basis on each atom
! nbas_t = total number of basis
! N_elec = number of electrons calculated by Tr(PS)
! nelec = number of electrons
! nside = number of polygon sides
! P_new() = New Density Matrix
! P_old() = Old Density Matrix
! PE,PE_t = Potential energy(primitive)/total
! R(i,j) = Distance between the nuclei "i" and "j"
! S,S_M = Overlap Integral : each primitive/Matrix
! step = the number of iterations
! title = title of the basis
! U() = Unitary Matrix
! X() = Matrix Transformation
! xp,yp,zp = Coordinates of P
! xq,yq,zq = Coordinates of Q
! Z(i) = Nuclear charge for atom i
! Zeta(i) = Slater exponent for atom i
!
!====================================================

!====================================================
!             SUBROUTINES(ALPHABETICAL ORDER)
!----------------------------------------------------
!

! 1. DGEMM('N','N',M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC) = Lapack Subroutine
!    Calculates C = alpha*A*B + beta*C
!
! 2. Distance(x1,y1,z1,x2,y2,z2,d) = Calculates the distance between atom 1 and 2
!
! 3. DSYEV('V', 'U', natom, U, natom, Eigen_value, WORK, LWORK, INFO) = Lapack Subroutine
!    Calculates the eigenvalues and eigenvectors of U
!
! 4. F0(t,F_0) = A function related to error function to calculate integrals
!
! 5. Grid(m,n,r) = Produces a grid with m rows and n columns with bond length of r
!
! 6. Input(natom,nbas) = Reads number of atoms and basis
!
! 7. Integral_2e(a,b,c,d,r1,r2,r3,Integral) = Calculates 2e integrals using 4 exponents: a,b,c,d and 3 distances:
!    r1 = R(A,B), r2 = R(C,D), r3 = R(P,Q) see fig below
!
! 8. Kinetic(a,b,d,KE) = Calculates Kinetic Energy integral using 2 exponents: a,b, and 1 distance:
!    r1 = R(A,B) see fig below
!
! 9. Overlap(a,b,d,S) = Calculates the overlap integral using exponents a,b and the distance d
!
! 10. Polygon(n,a) = It generates an equilateral polygon with n sides and length a with the origin
!    at the center of circumscrimbed circle, and saves the results in a file named "Input.xyz"
!
! 11. Potential(a,b,Z,d1,d2,PE) = Calculates Potential Energy integral using 2 exponents: a,b, and 2 distances:
!     d1 = R(A,B), d2 = R(P,C) see fig below


!====================================================
!                  FIGURE OF THE PROBLEM
!----------------------------------------------------
!        A------------P----B
!        C----------Q------D
!
!     A,B,C,D are centers for integration
!====================================================
  IMPLICIT NONE

  Integer ::  column, geo, i, itemp, j, jtemp, k, ktemp, l, ltemp, mat_dim, iatom,jatom, log
	Integer ::  m, n, nbas_t, nelec, o, p, nside, max_it, natom, nbas,norbit_fci,katom,latom, ninter
	Integer ::  row, step, lwork,LWORK_FCI, info, CI, primitive, norbit, counti,countj,countdim
	Integer,Allocatable :: 	Sub_basis(:)

  Double Precision :: a, alpha1, beta1, bond_incr, Converge, del,delta, d, d1, d2, d3, minimum
  Double Precision :: EN, EN_t, Integral, KE, max_bondlength, PE, pi, E1,E2,E_corr, num, stp
  Double Precision :: S, xp, yp, zp, xq, yq, zq, N_elec, sum, temp, H_atom, delta_OEP, incr, Tr_K_fci
	Double Precision :: Tr_P, Tr_K, Ws_old, Ws_new, E_thr, EN_old, delta_E, OEP_thr, result,alpha_inv

  Double Precision,Allocatable :: alpha(:,:), coeff(:,:), KE_t(:,:), R(:,:), KE_f(:,:), PE_f(:,:)
  Double Precision,Allocatable :: S_M(:,:), Z(:), Zeta(:),X(:,:), P_fci(:,:), del_P(:,:), H_core_new(:,:)
  Double Precision,Allocatable :: H_core(:,:), G(:,:), F(:,:), Work(:), E_CI(:), VP(:,:), den(:,:,:)
  Double Precision,Allocatable :: F_Prime(:,:), P_new(:,:), P_old(:,:), Eigen_C_H_t(:,:),WORK_FCI(:),E_CI_hamilton(:)
  Double Precision,Allocatable :: M_2e(:,:,:,:), PE_t(:,:),Eigen_value(:), H_t(:,:),Hamilton(:,:),Eigen_Hamilton(:,:)
  Double Precision,Allocatable :: Eigen_C(:,:), Eigen_C_Prime(:,:),U(:,:), PE_t_new(:,:)
	Double Precision,Allocatable :: P_fci_mo(:,:), P_fci_temp(:,:),KE_t_new(:,:),M_2e_new(:,:,:,:)
	Double Precision,Allocatable :: Slater_coeff(:,:), Gradient_new(:),Gradient_old(:),Hessian_inv(:,:)
	Double Precision,Allocatable :: direction(:,:),direction_1(:), hu(:), af(:), hu_column(:,:)
	Double Precision,Allocatable :: Grad_new_column(:,:), direction_column(:,:),dgrad(:),dgra_column(:,:)

	Double Precision,Allocatable :: dgra(:), P_HF(:,:,:),delta_den(:), S_orth(:,:)
	Double Precision,Allocatable :: dgra_row(:,:), PE_t_old(:,:), S_half_inv(:,:)
	Double Precision,Allocatable :: direction_row(:,:), Eigen_value_half(:), S_half(:,:)
	Double Precision,Allocatable :: hu_row(:,:),alpha_bracket(:,:), S_half_diag(:,:)
	Double Precision,Allocatable :: uhu(:,:), out_product_a(:,:), out_product_b(:,:)
	Double Precision,Allocatable :: temp1(:,:), temp2(:,:), temp3(:,:), S_half_temp(:,:)
	Double Precision,Allocatable :: Eigen_C_orth(:,:), S_half_orth(:,:),P_new_orth(:,:)

  Character(70) :: title,atomname

  pi = 4.0d0*datan(1.0d0)

  Converge = 1.0d-15
	OEP_thr = 1.0d-8
	E_thr = 1.0d-15
  max_it = 100
  max_bondlength = 2.1d0
  bond_incr = 0.1d0
  alpha1 = 1.0d0
  beta1 = 1.0d0
	incr = 0.5d0


!====================================================
!                  MAIN PROGRAM
!====================================================

	call Input(natom,nelec,nbas,CI)

	norbit = nbas*natom
	norbit_fci = norbit
	mat_dim = 2*norbit_fci-1 + (norbit_fci-1)*(norbit_fci-2)/2


	Allocate(Sub_basis(norbit))
	Allocate(Z(natom))
	Allocate(alpha(norbit,30))
	Allocate(coeff(norbit,30))

!====================================================
!                  INPUT BASIS
!====================================================


!----------------Reading basis from file----------------



	open(Unit = 3,File ="Basis.txt")
	read(3,*)
	read(3,"(A)") title
	read(3,*)
	do i=1,natom
		read(3,*) atomname, Z(i)
	end do
	read(3,*)


	nbas_t = 0
	do i=1,norbit

		read(3,*) atomname,Sub_basis(i)

		do j=1,Sub_basis(i)
		  read (3,*) alpha(i,j), coeff(i,j)
		end do

		nbas_t = nbas_t + Sub_basis(i)
	end do

	close(3)

	lwork = 3*nbas_t
	LWORK_FCI = 3*mat_dim
	ninter = norbit*(norbit+1)/2

	Allocate(S_M(norbit,norbit))
	Allocate(U(norbit,norbit))
	Allocate(KE_t(norbit,norbit))
	Allocate(KE_f(norbit,norbit))

	Allocate(PE_t(norbit,norbit))
	Allocate(PE_t_old(norbit,norbit))
	Allocate(PE_f(norbit,norbit))

	Allocate(KE_t_new(norbit,norbit))
	Allocate(PE_t_new(norbit,norbit))
	Allocate(H_core(norbit,norbit))
	Allocate(H_core_new(norbit,norbit))
	Allocate(Hamilton(16,16))
	Allocate(Eigen_Hamilton(16,16))

	Allocate(Eigen_C_H_t(mat_dim,mat_dim))
	Allocate(H_t(mat_dim,mat_dim))

	Allocate(P_new(norbit,norbit))
	Allocate(delta_den(norbit))
	Allocate(P_HF(norbit,norbit,nelec/2))

	Allocate(Slater_coeff(norbit,norbit))

	Allocate(S_half_inv(norbit,norbit))
	Allocate(Eigen_value_half(norbit))
	Allocate(S_half(norbit,norbit))
	Allocate(S_half_diag(norbit,norbit))
	Allocate(S_half_temp(norbit,norbit))
	Allocate(S_orth(norbit,norbit))
	Allocate(S_half_orth(norbit,norbit))
	Allocate(P_new_orth(norbit,norbit))


	Allocate(den(norbit_fci,norbit,norbit))

	Allocate(P_old(norbit,norbit))
	Allocate(P_fci(norbit,norbit))
	Allocate(P_fci_mo(norbit,norbit))
	Allocate(P_fci_temp(norbit,norbit))

	Allocate(del_P(norbit,norbit))
	Allocate(VP(norbit,norbit))

	Allocate(direction(norbit,norbit))
	Allocate(direction_1(ninter))
	Allocate(direction_column(ninter,1))

	Allocate(Hessian_inv(ninter,ninter))
 	Allocate(Gradient_new(ninter))
	Allocate(Grad_new_column(ninter,1))
	Allocate(Gradient_old(ninter))
	Allocate(dgrad(ninter))
	Allocate(dgra_column(ninter,1))
	Allocate(hu(ninter))
	Allocate(hu_column(ninter,1))
	Allocate(af(ninter))

	Allocate(dgra(ninter))
	Allocate(dgra_row(1,ninter))
	Allocate(direction_row(1,ninter))
	Allocate(hu_row(1,ninter))

	Allocate(uhu(1,1))
	Allocate(alpha_bracket(1,1))
	Allocate(out_product_a(ninter,ninter))
	Allocate(out_product_b(ninter,ninter))
	Allocate(temp1(ninter,ninter))
	Allocate(temp2(ninter,ninter))
	Allocate(temp3(ninter,ninter))

	Allocate(F_Prime(norbit,norbit))
	Allocate(F(norbit,norbit))
	Allocate(G(norbit,norbit))
	Allocate(X(norbit,norbit))
	Allocate(Eigen_C_Prime(norbit,norbit))
	Allocate(Eigen_C(norbit,norbit))
	Allocate(Eigen_C_orth(norbit,norbit))
	Allocate(M_2e(norbit,norbit,norbit,norbit))
	Allocate(M_2e_new(norbit,norbit,norbit,norbit))

	Allocate(Eigen_value(norbit))
	Allocate(Work(lwork))
	Allocate(WORK_FCI(LWORK_FCI))
	Allocate(Zeta(natom))
	Allocate(R(natom,3))
	Allocate(E_CI(mat_dim))
	Allocate(E_CI_hamilton(16))


!====================================================
!         SCALING of Coefficients & Exponents
!----------------------------------------------------

	do i=1,norbit
	  do j=1,Sub_basis(i)
!	    alpha(i,j) = alpha(i,j)*(zeta(i)**2.0d0)
	    coeff(i,j) = coeff(i,j)*((2.0d0*alpha(i,j)/pi)**0.75d0)
	  end do
	end do

!====================================================
!                  INPUT GEOMETRY
!====================================================
!	40 continue
!	write(*,*) 'Which Geometry are you interested in:'
!	write(*,*) '1. Polygon'
!	write(*,*) '2. Grid'
!	read(*,*) geo

!	if (geo.eq.1) then
!		write(*,*) 'Please insert the number of Polygon side:'
!		read(*,*) nside
!		write(*,*) 'Please insert the Polygon side length(bond length):'
!		read(*,*) a
!	else if (geo.eq.2) then
!		write(*,*) 'Please insert the number of rows:'
!		read(*,*) row
!		write(*,*) 'Please insert the number of columns:'
!		read(*,*) column
!		write(*,*) 'Please insert the bond length:'
!		read(*,*) a
!	else
!		write(*,*) 'Input data is invalid!'
!		goto 40
!	end if

!	30 continue

!	if (geo.eq.1) then
!  	call Polygon(nside,a)
!  else if (geo.eq.2) then
!		call Grid(row,column,a)
!	end if



!----------------Reading coordinates from file----------------

	open(Unit = 4,File ="Input.xyz")
	do i=1,natom
	  read(4,*) atomname,R(i,1),R(i,2),R(i,3)
	end do


!====================================================
!                  OVERLAP INTEGRAL
!----------------------------------------------------

	open(Unit = 2,File ="Output.txt" )
	write(2,*) title
!-------------loop over atoms------------------------
	do i=1,norbit
		iatom = int((i-1)/nbas)+1
	  do j=1,norbit
			S_M(i,j) = 0.0d0
			jatom = int((j-1)/nbas)+1
!-------------loop over primitives-------------------
	    do k=1,Sub_basis(i)
	      do l=1,Sub_basis(j)
					call Distance(R(iatom,1),R(iatom,2),R(iatom,3),R(jatom,1),R(jatom,2),R(jatom,3),d)
	        call Overlap(alpha(i,k),alpha(j,l),d,S)
		      S_M(i,j) = S_M(i,j) + S*coeff(i,k)*coeff(j,l)
	      end do
	    end do
!----------------------------------------------------

		end do

	end do

!====================================================
!                  KINETIC ENERGY
!----------------------------------------------------

!-------------loop over atoms------------------------
	H_core = 0.0d0

	do i=1,norbit
		iatom = int((i-1)/nbas)+1
	  do j=1,norbit
			KE_t(i,j) = 0.0d0
			jatom = int((j-1)/nbas)+1

!-------------loop over primitives-------------------
	    do k=1,Sub_basis(i)
	      do l=1,Sub_basis(j)
					call Distance(R(iatom,1),R(iatom,2),R(iatom,3),R(jatom,1),R(jatom,2),R(jatom,3),d)
	        call Kinetic(alpha(i,k),alpha(j,l),d,KE)

		      KE_t(i,j) = KE_t(i,j) + KE*coeff(i,k)*coeff(j,l)

	      end do
	    end do
!----------------------------------------------------
			H_core(i,j) = H_core(i,j) + KE_t(i,j)

		end do
	end do

!====================================================
!                  POTENTIAL ENERGY
!----------------------------------------------------
!----------------------------------------------------
!                (A|sum(Zc/(r-R))|B)
!----------------------------------------------------


!-------------loop over atoms------------------------

	do i=1,norbit
		iatom = int((i-1)/nbas)+1

	  do j=1,norbit
			PE_t(i,j) = 0.0d0
			jatom = int((j-1)/nbas)+1

!-------------loop over primitives-------------------
	    do l=1,Sub_basis(i)
	      do m=1,Sub_basis(j)

!-----------------Calculating the P Coordinates------

					xp = (alpha(i,l)*R(iatom,1)+alpha(j,m)*R(jatom,1) )/(alpha(i,l)+alpha(j,m))
					yp = (alpha(i,l)*R(iatom,2)+alpha(j,m)*R(jatom,2) )/(alpha(i,l)+alpha(j,m))
					zp = (alpha(i,l)*R(iatom,3)+alpha(j,m)*R(jatom,3) )/(alpha(i,l)+alpha(j,m))

!----------------loop over all nuclei(for Zc)--------

					do k=1,natom

!----------------------------------------------------
!    d1 is the distance between P and C
!    d is the distance between A and B
!----------------------------------------------------

					call Distance(xp,yp,zp,R(k,1),R(k,2),R(k,3),d1)
					call Distance(R(iatom,1),R(iatom,2),R(iatom,3),R(jatom,1),R(jatom,2),R(jatom,3),d)
					call Potential(alpha(i,l),alpha(j,m),Z(k),d,d1,PE)

		  			PE_t(i,j) = PE_t(i,j) + PE*coeff(i,l)*coeff(j,m)

	        end do


	      end do
	    end do
			H_core(i,j) = H_core(i,j) + PE_t(i,j)

!----------------------------------------------------
	  end do
	end do


!====================================================
!                 2e-Integrals
!----------------------------------------------------

  M_2e = 0.0d0

  do i=1,norbit
		iatom = int((i-1)/nbas)+1

    do j=1,norbit
			jatom = int((j-1)/nbas)+1

      do k=1,norbit
				katom = int((k-1)/nbas)+1

        do l=1,norbit
					latom = int((l-1)/nbas)+1

!-------------loop over primitives-------------------
          do m=1,Sub_basis(i)
            do n=1,Sub_basis(j)
	      			do o=1,Sub_basis(k)
	        			do p=1,Sub_basis(l)

!-----------------Calculating the P Coordinates------

		  						xp = (alpha(i,m)*R(iatom,1)+alpha(j,n)*R(jatom,1) )/(alpha(i,m)+alpha(j,n))
		  						yp = (alpha(i,m)*R(iatom,2)+alpha(j,n)*R(jatom,2) )/(alpha(i,m)+alpha(j,n))
		  						zp = (alpha(i,m)*R(iatom,3)+alpha(j,n)*R(jatom,3) )/(alpha(i,m)+alpha(j,n))

!-----------------Calculating the Q Coordinates------

		  						xq = (alpha(k,o)*R(katom,1)+alpha(l,p)*R(latom,1) )/(alpha(k,o)+alpha(l,p))
		  						yq = (alpha(k,o)*R(katom,2)+alpha(l,p)*R(latom,2) )/(alpha(k,o)+alpha(l,p))
		  						zq = (alpha(k,o)*R(katom,3)+alpha(l,p)*R(latom,3) )/(alpha(k,o)+alpha(l,p))

!--------------------------------------------------------------------------------------
!      d1 = distance(A,B) , d2 = distance(C,D) , d3 = distance(P,Q)
!--------------------------------------------------------------------------------------

		  						call Distance(R(iatom,1),R(iatom,2),R(iatom,3),R(jatom,1),R(jatom,2),R(jatom,3),d1)
		  						call Distance(R(katom,1),R(katom,2),R(katom,3),R(latom,1),R(latom,2),R(latom,3),d2)
		  						call Distance(xp,yp,zp,xq,yq,zq,d3)

		  						call Integral_2e(alpha(i,m),alpha(j,n),alpha(k,o),alpha(l,p),d1,d2,d3,Integral)

		  						M_2e(i,j,k,l) = M_2e(i,j,k,l) + Integral*coeff(i,m)*coeff(j,n)*coeff(k,o)*coeff(l,p)
                end do
	      			end do
	    			end do
	  			end do
        end do
      end do
    end do
  end do

!===========================================================================
!                   WRITING THE WHOLE 2e-MATRIX
!---------------------------------------------------------------------------

!	write(2,*) "All 2e Matrix"
!	do i=1,norbit
!	  do j=1,norbit
!	    do k=1,norbit
!	      do l=1,norbit
!	        write(2, "(A1,I3,I3,A1,I3,I3,A1,f16.10)") '(',i,j,'|',k,l,')',M_2e(i,j,k,l)
!	      end do
!	    end do
!	  end do
!	end do

!	write(2,*) "**********************************************"





!=====================================================
!             RHF PROGRAM (CLOSED SHELL)
!=====================================================

!=====================================================
!        ORTHOGONALIZATION OF THE BASIS(CANONICAL)
!        X = Us^-0.5, where s = diagonalized S
!-----------------------------------------------------

!-------------Diagonalizing S and producing s^-0.5----

	U = 0.0d0
	do i=1,norbit
	  do j=1,norbit
	    U(i,j) = S_M(i,j)
	  end do
	end do

!-------------Calculating eigenvalues & eigenvectors ----

	Eigen_value = 0.0d0
	call DSYEV('V', 'U', norbit, U, norbit, Eigen_value, WORK, LWORK, INFO)

	do i=1,norbit
	  Eigen_value_half(i) = (Eigen_value(i))**(0.5d0)
	end do

	do i=1,norbit
	  Eigen_value(i) = (Eigen_value(i))**(-0.5d0)
	end do

!-----------------Transformation----------------------

	X = 0.0d0
	do i=1,norbit
	  do j=1,norbit
	    X(i,j) = Eigen_value(j)*U(i,j)
	  end do
	end do

!	write(2,*) "X Matrix"
!	do i=1,norbit
!	  do j=1,norbit
!	    write(2, "(f16.10)", advance = "No") X(i,j)
!	  end do
!	write(2,*)
!	end do

!--------------S^-0.5  for future use!!!--------------------

	S_half_diag = 0.0d0

	do i=1,norbit
	  S_half_diag(i,i) = Eigen_value(i)
	end do

	S_half_temp = 0.0d0
	call DGEMM('N','N',norbit,norbit,norbit,alpha1,U,norbit,S_half_diag,norbit,beta1,S_half_temp,norbit)

	S_half_inv = 0.0d0
	call DGEMM('N','T',norbit,norbit,norbit,alpha1,S_half_temp,norbit,U,norbit,beta1,S_half_inv,norbit)

!--------------S^0.5  for future use!!!--------------------

	S_half_diag = 0.0d0

	do i=1,norbit
    S_half_diag(i,i) = Eigen_value_half(i)
	end do

	S_half_temp = 0.0d0
	call DGEMM('N','N',norbit,norbit,norbit,alpha1,U,norbit,S_half_diag,norbit,beta1,S_half_temp,norbit)

	S_half = 0.0d0
	call DGEMM('N','T',norbit,norbit,norbit,alpha1,S_half_temp,norbit,U,norbit,beta1,S_half,norbit)



!---------------------Overlap Matrix------------------
		  	write(2,*) "Overlap matrix"
				do i=1,norbit
				  do j=1,norbit
				    write(2, "(f16.8)", advance = "No") S_M(i,j)
				  end do
					write(2,*)
				end do
				write(2,*) "****************************************"
				write(2,*) "  ************************************  "
!---------------------Kinetic Matrix------------------
		  	write(2,*) "Kinetic matrix"
				do i=1,norbit
				  do j=1,norbit
				    write(2, "(f16.8)", advance = "No") KE_t(i,j)
				  end do
					write(2,*)
				end do
				write(2,*) "****************************************"
				write(2,*) "  ************************************  "

!---------------------Potential Matrix------------------
				write(2,*) "Potential matrix"
				do i=1,norbit
				  do j=1,norbit
			    write(2, "(f16.8)", advance = "No") PE_t(i,j)
					  end do
					write(2,*)
					end do
			write(2,*) "****************************************"
			write(2,*) "  ************************************  "
!====================================================
!                 SCF CALCULATION
!====================================================

	step = 0

	write(2,*) "Initial Density Matrix"
	do i=1,norbit
	  do j=1,norbit
	    P_old(i,j) = 0.0d0
	    P_new(i,j) = 0.0d0
	    write(2, "(f16.10)", advance = "No") P_old(i,j)
	  end do
	write(2,*)
	end do

	write(2,*) "**********************************************"

	20 continue


	write(2,*) "Step", step
	step = step+1

!---------------------Form the 2e part of Fock Matrix-----------

!	write(2,*) "G Matrix"
	do i=1,norbit
	  do j=1,norbit
			G(i,j) = 0.0d0
	    do k=1,norbit
	      do l=1,norbit
	        G(i,j) = G(i,j) + P_new(k,l)*(M_2e(i,j,k,l)-0.5d0*M_2e(i,l,k,j))
	      end do
	    end do
!	    write(2, "(f16.10)", advance = "No") G(i,j)
	  end do
!	write(2,*)
	end do

!------------------Constructing Fock Matrix--------------

!	write(2,*) "Fock Matrix"
	F = 0.0d0
	do i=1,norbit
	  do j=1,norbit
	    F(i,j) = H_core(i,j)+G(i,j)
!	    write(2, "(f16.10)", advance = "No") F(i,j)
	  end do
!	write(2,*)
	end do
!	F(1,1) = 0.0d0
!	write(2,*)

!-----------------------Electronic Energy----------------

	EN = 0.0d0
	do i=1,norbit
	  do j=1,norbit
	    EN = EN + 0.5d0*P_new(i,j)*(H_core(i,j)+F(i,j))
	  end do
	end do
	write(2,*) "Electronic Energy = ",EN
	write(2,*)
	EN_old = EN

!----------------------Transform and Diagonalize---------

!----------------------G = F*X----------------------------

	G = 0.0d0
  call DGEMM('N','N',norbit,norbit,norbit,alpha1,F,norbit,X,norbit,beta1,G,norbit)

!----------------------F_Prime = XT*G = XT*F*X--------------------

	F_Prime = 0.0d0
	call DGEMM('T','N',norbit,norbit,norbit,alpha1,X,norbit,G,norbit,beta1,F_Prime,norbit)

	do i=1,norbit
	  do j=1,norbit
	    Eigen_C_prime(i,j) = F_Prime(i,j)
	  end do
	end do

!----------------------Eigenvalues and Eigenvectors of F_Prime---------

	Eigen_value = 0.0d0
	call DSYEV('V', 'U', norbit, Eigen_C_prime, norbit, Eigen_value, WORK, LWORK, INFO)

!---------------------- X*C_Prime = C---------------------------------

	Eigen_C = 0.0d0
	call DGEMM('N','N',norbit,norbit,norbit,alpha1,X,norbit,Eigen_C_prime,norbit,beta1,Eigen_C,norbit)

!------------------Creating New Density Matrix----------

	P_old = P_new

	do i=1,norbit
	  do j=1,norbit
	    P_new(i,j) = 0.0d0
	    do k=1,nelec/2
	      	P_new(i,j) = P_new(i,j) + 2*Eigen_C(i,k)*Eigen_C(j,k)
	    end do
	  end do
	end do

!---------------------Writing F', C', Energy matrix, C, P--------------------

!	write(2,*) "Transformed Fock Matrix"
!	do i=1,norbit
!	  do j=1,norbit
!	    write(2, "(f16.10)", advance = "No") F_Prime(i,j)
!	  end do
!	write(2,*)
!	end do
!
!	write(2,*) "Transformed C Matrix"
!	do i=1,norbit
!	  do j=1,norbit
!	    write(2, "(f16.10)", advance = "No") Eigen_C_Prime(i,j)
!	  end do
!	write(2,*)
!	end do
!
!	write(2,*) "E Matrix"
!	do i=1,norbit
!	  do j=1,norbit
!	    if (i.eq.j) then
!	      write(2, "(f16.10)", advance = "No") Eigen_value(i)
!	    else
!	      write(2, "(f16.10)", advance = "No") 0.0d0
!	    end if
!	  end do
!	write(2,*)
!	end do
!
!	write(2,*) "C Matrix"
!	do i=1,norbit
!	  do j=1,norbit
!	    write(2, "(f16.10)", advance = "No") Eigen_C(i,j)
!	  end do
!	write(2,*)
!	end do
!
!	write(2,*) "Density Matrix"
!	do i=1,norbit
!	  do j=1,norbit
!	    write(2, "(f16.10)", advance = "No") P_new(i,j)
!	  end do
!	write(2,*)
!	end do


!---------------------Convergence of Density Matrix(Using Standard Deviation)---------

	delta = 0.0d0
	do i=1,norbit
	  do j=1,norbit
	    delta = delta+(P_new(i,j)-P_old(i,j))**2
	  end do
	end do

	delta = dsqrt(delta/(ninter))
	write(2,*) "Delta in Density = ",delta

	EN = 0.0d0
	do i=1,norbit
	  do j=1,norbit
	    EN = EN + 0.5d0*P_new(i,j)*(H_core(i,j)+F(i,j))
	  end do
	end do
	write(2,*) "Electronic Energy = ",EN
	write(2,*)

	delta_E = dabs(EN-EN_old)
	write(2,*) "Delta in Energy= ",delta_E


!		EN_t = EN_t/natom
!		write(2,*) "Total Energy = ",EN_t


	write(2,*) "*********************************************"

!-----------------Check for Convergence-------------

	if (delta .lt. Converge .and. delta_E .lt. E_thr) then
	  EN_t = EN
	  do i=1,natom-1
	    do j=i+1,natom
	      call Distance(R(i,1),R(i,2),R(i,3),R(j,1),R(j,2),R(j,3),d)
	      EN_t = EN_t + Z(i)*Z(j)/d
	    end do
	  end do

!	  EN_t = EN_t/natom
		write(*,*) "HF Calculation Converged!"
		write(*,*) "CI Calculation Started!"

	  open(Unit = 7,File = 'Energy')
		write(2,*) "Hartree-Fock Calculation: "
	  write(2,*) "HF Electronic Energy = ",EN
	  write(2,*) "Total Energy = ",EN_t
		write(2,*) "*********************************************"

	else if(step .lt. max_it) then
		EN_old = EN
	  goto 20
	else
	  write(2,*) "It does not converge in ",max_it," steps!"
!		goto 100
	end if

!====================================================
!                 CI CALCULATION
!====================================================


!====================================================
!              H_core in MO basis
!----------------------------------------------------

H_core_new = 0.0d0

	do i=1,norbit
	  do j=1,norbit
			do k=1,norbit
				do l=1,norbit
		      H_core_new(i,j) = H_core_new(i,j) + H_core(k,l)*Eigen_C(k,i)*Eigen_C(l,j)
				end do
			end do
		end do
	end do


!====================================================
!              2e-INTEGRALS OF FINAL HF in MO basis
!----------------------------------------------------

  M_2e_new = 0.0d0
  do i=1,norbit
    do j=1,norbit
      do k=1,norbit
        do l=1,norbit
!-----------------------------------------------------
          do m=1,norbit
            do n=1,norbit
	      			do o=1,norbit
	        			do p=1,norbit
		  						M_2e_new(i,j,k,l) = M_2e_new(i,j,k,l) + M_2e(m,n,o,p)*Eigen_C(m,i)*Eigen_C(n,j)*Eigen_C(o,k)*Eigen_C(p,l)
                end do
	      			end do
	    			end do
	  			end do
        end do
      end do
    end do
  end do



!====================================================
!     CALCULATING CID: Hamiltonian
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
		! Double Check this:
				H_t(counti,countj) = H_core_new(1,1) + H_core_new(i,i) + M_2e_new(1,1,i,i) + M_2e_new(1,i,i,1)
			else
				! Double Check this:

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

	!		H_t(counti,countj) = 2.0d0*(Eigen_value(i)-Eigen_value(1)) + M_2e_new(1,1,1,1) + M_2e_new(i,i,i,i) &
	!		- 4.0d0*M_2e_new(1,1,i,i) + 2.0d0*M_2e_new(i,1,1,i)
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

!H_t(counti,countj) = (Eigen_value(i)+Eigen_value(j)) - 2.0d0*Eigen_value(1) + M_2e_new(1,1,1,1) + M_2e_new(i,i,j,j) &
!			+ M_2e_new(i,j,j,i) - 2.0d0*M_2e_new(j,j,1,1)- 2.0d0*M_2e_new(i,i,1,1) &
!			+ M_2e_new(i,1,1,i) + M_2e_new(j,1,1,j)
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

!---------------------Reading Hamiltonian from molpro-----------------------


!open(Unit = 10,File ="Hamilton.txt")
!read(10,*)
!do i=1,16
!	read(10,*)
!	read(10,*)
!	do j=1,16
!		read(10,*) Hamilton(j,i)
!	end do
!end do



!---------------------Diagonalizing H Molpro-----------------------

	!	E_CI_hamilton = 0.0d0
	!	Eigen_Hamilton = Hamilton
	!	call DSYEV('V', 'U', 16, Eigen_Hamilton, 16, E_CI_hamilton, WORK_FCI, LWORK_FCI, INFO)
!------------------------------------------------------------

	!---------------------Total Energy---------------------------

!	write(2,*) "Eigenvalues of FCI Hamiltonian: "
!	do i=1,mat_dim
!		write(2,*)	E_CI(i)
!	end do
!	write(2,*)



	write(*,*) "correlation energy:", E_CI(1)-H_t(1,1)
	write(2,*) "****************************************"
	write(2,*) "  ************************************  "

!---------------------FCI Eigenvectors------------------
!	write(2,*) "Eigenvectors of FCI Hamiltonian:"
!	do i=1,mat_dim
!		do j=1,mat_dim
!			write(2, "(f16.8)", advance = "No") Eigen_C_H_t(i,j)
!		end do
!		write(2,*)
!	end do
!	write(2,*) "****************************************"
!	write(2,*) "  ************************************  "

!-------------Density matrix using CI Vector---------
!     			P_fci_mo(p,q)	= <Di|ap aq + ap aq |Dj>
!---------------------------------------------------------

		P_fci = 0.0d0
		P_fci_mo = 0.0d0


!-----------------------------------------------------------------------
!										Manual Construction for 6-31g
!-----------------------------------------------------------------------

!	P_fci_mo(1,1) = 2.0d0*Eigen_C_H_t(1,1)**2.0d0+Eigen_C_H_t(2,1)**2.0d0 &
!	+Eigen_C_H_t(3,1)**2.0d0+Eigen_C_H_t(4,1)**2.0d0

!	P_fci_mo(1,2) = (Eigen_C_H_t(2,1)*Eigen_C_H_t(5,1)+Eigen_C_H_t(1,1)*Eigen_C_H_t(2,1))*dsqrt(2.0d0) &
!	+ Eigen_C_H_t(3,1)*Eigen_C_H_t(8,1)+ Eigen_C_H_t(4,1)*Eigen_C_H_t(9,1)


!	P_fci_mo(1,3) = (Eigen_C_H_t(3,1)*Eigen_C_H_t(6,1)+Eigen_C_H_t(1,1)*Eigen_C_H_t(3,1))*dsqrt(2.0d0) &
!	+ Eigen_C_H_t(2,1)*Eigen_C_H_t(8,1)+ Eigen_C_H_t(4,1)*Eigen_C_H_t(10,1)

!	P_fci_mo(1,4) = (Eigen_C_H_t(4,1)*Eigen_C_H_t(7,1)+Eigen_C_H_t(1,1)*Eigen_C_H_t(4,1))*dsqrt(2.0d0) &
!	+ Eigen_C_H_t(2,1)*Eigen_C_H_t(9,1)+ Eigen_C_H_t(3,1)*Eigen_C_H_t(10,1)

!	P_fci_mo(2,1) = P_fci_mo(1,2)


!	P_fci_mo(2,2) = 2.0d0*Eigen_C_H_t(5,1)**2.0d0+Eigen_C_H_t(2,1)**2.0d0 &
!	+Eigen_C_H_t(8,1)**2.0d0+Eigen_C_H_t(9,1)**2.0d0

!	P_fci_mo(2,3) = (Eigen_C_H_t(8,1)*Eigen_C_H_t(6,1)+Eigen_C_H_t(5,1)*Eigen_C_H_t(8,1))*dsqrt(2.0d0) &
!	+ Eigen_C_H_t(2,1)*Eigen_C_H_t(3,1)+ Eigen_C_H_t(9,1)*Eigen_C_H_t(10,1)


!	P_fci_mo(2,4) = (Eigen_C_H_t(7,1)*Eigen_C_H_t(9,1)+Eigen_C_H_t(5,1)*Eigen_C_H_t(9,1))*dsqrt(2.0d0) &
!	+ Eigen_C_H_t(2,1)*Eigen_C_H_t(4,1)+ Eigen_C_H_t(8,1)*Eigen_C_H_t(9,1)

!	P_fci_mo(3,1) = P_fci_mo(1,3)

!	P_fci_mo(3,2) = P_fci_mo(2,3)

!	P_fci_mo(3,3) = 2.0d0*Eigen_C_H_t(6,1)**2.0d0+Eigen_C_H_t(3,1)**2.0d0 &
!	+Eigen_C_H_t(8,1)**2.0d0+Eigen_C_H_t(10,1)**2.0d0

!	P_fci_mo(3,4) = (Eigen_C_H_t(7,1)*Eigen_C_H_t(10,1)+Eigen_C_H_t(6,1)*Eigen_C_H_t(10,1))*dsqrt(2.0d0) &
!	+ Eigen_C_H_t(3,1)*Eigen_C_H_t(4,1)+ Eigen_C_H_t(8,1)*Eigen_C_H_t(9,1)

!	P_fci_mo(4,1) = P_fci_mo(1,4)

!	P_fci_mo(4,2) = P_fci_mo(2,4)

!	P_fci_mo(4,3) = P_fci_mo(3,4)

!	P_fci_mo(4,4) = 2.0d0*Eigen_C_H_t(7,1)**2.0d0+Eigen_C_H_t(4,1)**2.0d0 &
!	+Eigen_C_H_t(9,1)**2.0d0+Eigen_C_H_t(10,1)**2.0d0
!------------------------------------------------------------------------------


!---------------------fci Density Matrix MO basis----------
!		write(2,*) "Manual FCI Density matrix MO basis 6-31g:"
!		do i=1,norbit
!			do j=1,norbit
!				write(2, "(f16.8)", advance = "No") P_fci_mo(i,j)
!			end do
!			write(2,*)
!		end do
!		write(2,*) "****************************************"
!		write(2,*) "  ************************************  "


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
	do i=1,norbit
		do j=1,norbit
			write(2, "(f16.8)", advance = "No") Slater_coeff(i,j)
		end do
		write(2,*)
	end do
	write(2,*) "****************************************"
	write(2,*) "  ************************************  "

	P_fci_mo =0.0d0

	do i=1,norbit_fci
		do j=i,norbit_fci

!--------------------------------------------------------------------
		 	if (i==j) then
				do k=1,norbit_fci
					if (i/=k) then
						P_fci_mo(i,j) =P_fci_mo(i,j)+Slater_coeff(i,k)**2.0d0
					end if
					if (i==k) then
						P_fci_mo(i,j) =P_fci_mo(i,j)+ 2.0d0*Slater_coeff(i,k)**2.0d0
					end if
				end do
!--------------------------------------------------------------------
			else
				P_fci_mo(i,j) = dsqrt(2.0d0)*(Slater_coeff(i,j)*Slater_coeff(i,i)+Slater_coeff(j,j)*Slater_coeff(i,j))

				do k=1,norbit_fci
					if (i/=k .and. j/=k) then
						P_fci_mo(i,j) = P_fci_mo(i,j) + Slater_coeff(k,j)*Slater_coeff(i,k)
					end if
				end do

				P_fci_mo(j,i)= P_fci_mo(i,j)
			end if
		end do
	end do


	P_fci_temp = 0.0d0
call DGEMM('N','N',norbit,norbit,norbit,alpha1,Eigen_C,norbit,P_fci_mo,norbit,beta1,P_fci_temp,norbit)

	P_fci = 0.0d0
call DGEMM('N','T',norbit,norbit,norbit,alpha1,P_fci_temp,norbit,Eigen_C,norbit,beta1,P_fci,norbit)


!---------------------CI Vector------------------
	write(2,*) "CI Vector:"
	do i=1,mat_dim
			write(2, "(f16.8)") Eigen_C_H_t(i,1)
		end do
	write(2,*) "****************************************"
	write(2,*) "  ************************************  "

	!---------------------fci Density Matrix MO basis----------
			write(2,*) "FCI Density matrix MO basis:"
			do i=1,norbit
				do j=1,norbit
					write(2, "(f16.8)", advance = "No") P_fci_mo(i,j)
				end do
				write(2,*)
			end do
			write(2,*) "****************************************"
			write(2,*) "  ************************************  "

!---------------------fci Density Matrix AO basis----------
		write(2,*) "FCI Density matrix AO basis:"
		do i=1,norbit
			do j=1,norbit
				write(2, "(f16.8)", advance = "No") P_fci(i,j)
			end do
			write(2,*)
		end do
		write(2,*) "****************************************"
		write(2,*) "  ************************************  "

!---------------------HF Density Matrix-----------------

!				write(2,*) "HF Density Matrix"
!					do i=1,norbit
!						do j=1,norbit
!							write(2, "(f20.16)", advance = "No") P_new(i,j)
!						end do
!						write(2,*)
!					end do
!					write(2,*) "****************************************"
!					write(2,*) "  ************************************  "


				!---------------------H_t Matrix------------------
		!			write(2,*) "H_t Matrix"
		!			do i=1,mat_dim
		!			  do j=1,mat_dim
		!			    write(2, "(f16.8)", advance = "No") H_t(i,j)
		!			  end do
		!				write(2,*)
		!			end do
		!			write(2,*) "****************************************"
		!			write(2,*) "  ************************************  "
					!---------------------HF Orbital Energies------------------
		!				write(2,*) "HF Orbital Energies"
		!				do i=1,norbit
		!				    write(2, "(f16.8)") Eigen_value(i)
		!				end do
		!				write(2,*) "****************************************"
		!				write(2,*) "  ************************************  "

		!				write(*,*) "CI Calculation is Done!"

!						goto 100
		!---------------------Mulliken Population Analysis-----------------

!		N_elec = 0.0d0
!		do i=1,norbit
!			do j=1,norbit
!				N_elec = N_elec + P_new(i,j)*S_M(j,i)
!			end do
!		end do
!		write(2,*) "Num of electrons via Trace(PS): ", N_elec


	  write(7,*) EN_t

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
!          Intrinsic Atomic Orbitals (IAO)
!====================================================

!----------------------------------------------------
!          Determining the fragment core orbitals
!----------------------------------------------------

	S_half_temp = 0.0d0
	call DGEMM('T','N',norbit,norbit,norbit,alpha1,Eigen_C,norbit,S_M,norbit,beta1,S_half_temp,norbit)
	S_orth = 0.0d0
	call DGEMM('N','N',norbit,norbit,norbit,alpha1,S_half_temp,norbit,Eigen_C,norbit,beta1,S_orth,norbit)


	write(*,*) "Orthognalized S:"
	do i=1,norbit
		do j=1,norbit
			write(*, "(f16.8)", advance = "No") S_orth(i,j)
		end do
		write(*,"(f16.8)")
	end do


!	S_half_temp = 0.0d0
!	call DGEMM('T','N',norbit,norbit,norbit,alpha1,Eigen_C,norbit,S_half_inv,norbit,beta1,S_half_temp,norbit)
!	S_half_orth = 0.0d0
!	call DGEMM('N','N',norbit,norbit,norbit,alpha1,S_half_temp,norbit,Eigen_C,norbit,beta1,S_half_orth,norbit)

	write(*,*) "S^-1/2:"
	do i=1,norbit
		do j=1,norbit
			write(*, "(f16.8)", advance = "No") S_half_inv(i,j)
		end do
		write(*,"(f16.8)")
	end do

	Eigen_C_orth = 0.0d0
	call DGEMM('N','N',norbit,norbit,norbit,alpha1,Eigen_C,norbit,S_half_inv,norbit,beta1,Eigen_C_orth,norbit)

	write(*,*) "Orthognalized Eigenvectors:"
	do i=1,norbit
		do j=1,norbit
			write(*, "(f16.8)", advance = "No") Eigen_C_prime(i,j)
		end do
		write(*,"(f16.8)")
	end do

	!S_orth = 0.0d0
	!call DGEMM('N','N',norbit,norbit,norbit,alpha1,S_half_temp,norbit,Eigen_C,norbit,beta1,S_orth,norbit)


	do i=1,norbit
		do j=1,norbit
			P_new_orth(i,j) = P_new_orth(i,j) + 2.0d0*Eigen_C_prime(i,k)*Eigen_C_prime(j,k)
		end do
	end do

	write(*,*) "Orthognalized HF density:"
	do i=1,norbit
		do j=1,norbit
			write(*, "(f16.8)", advance = "No") P_new_orth(i,j)
		end do
		write(*,"(f16.8)")
	end do
	write(*,*) "HF density:"
	do i=1,norbit
		do j=1,norbit
			write(*, "(f16.8)", advance = "No") P_new(i,j)
		end do
		write(*,"(f16.8)")
	end do
!====================================================
!                   OEP
!====================================================

	write(*,*) "OEP Calculation Started!"

	step = 0
	PE_t_new = F - KE_t
	Ws_old = 0.0d0
 	Ws_new = 0.0d0
	Gradient_new = 0.0d0
	Gradient_old = 0.0d0

	P_HF = 0.0d0
	do k=1,nelec/2
		do i=1,norbit
			do j=1,norbit
				P_HF(i,j,k) = P_HF(i,j,k) + 2.0d0*Eigen_C(i,k)*Eigen_C(j,k)
			end do
		end do
	end do

!	P_fci = P_new

!	Open(Unit = 9, File = 'densityin.txt')
!	read(9,*)
!	do i=1,norbit
!		read(9,*) (P_new(i,j), j=1,norbit)
!	end do

!-------------Set Hessian as I-------------------------

	Hessian_inv = 0.0d0
	do i=1,ninter
		Hessian_inv(i,i) = 1.0d0
	end do
!------------------------------------------------------

	write(2,*) "OEP Calculation:"
	Open(Unit = 8, File = 'Ws.txt')
	write(8,*) 'Step','                   ','delta_OEP','         ','Ws'
	write(8,*) '----------','          ','----------','          ','----------'

	KE_f = 0.0d0
	call DGEMM('N','N',norbit,norbit,norbit,alpha1,P_fci,norbit,KE_t,norbit,beta1,KE_f,norbit)

	Tr_K = 0.0d0
	do i=1,norbit
		Tr_K = Tr_K + KE_f(i,i)
	end do

	Tr_K_fci = Tr_K
	write(*,*) "Ws max", Tr_K
	write(*,*)

!	VP = 0.0d0
!	call DGEMM('N','N',norbit,norbit,norbit,alpha1,P_fci,norbit,PE_t_new,norbit,beta1,VP,norbit)

!	Tr_P = 0.0d0
!	do i=1,norbit
!		Tr_P = Tr_P + VP(i,i)
!	end do

!	write(*,*) "FCI Potential energy:", Tr_P
!	write(*,*)

!-------------------------OEP zeroth step---------------------

	write(2,*) "OEP Step:",step


!-------------------------Vext * (P_new-P_fci)-----------------
	del_P = P_new-P_fci
	VP = 0.0d0
	call DGEMM('N','N',norbit,norbit,norbit,alpha1,del_P,norbit,PE_t_new,norbit,beta1,VP,norbit)


	KE_f = 0.0d0
	call DGEMM('N','N',norbit,norbit,norbit,alpha1,P_new,norbit,KE_t,norbit,beta1,KE_f,norbit)

	Tr_K = 0.0d0
	do i=1,norbit
		Tr_K = Tr_K + KE_f(i,i)
	end do

	Tr_P = 0.0d0
	do i=1,norbit
		Tr_P = Tr_P + VP(i,i)
	end do



	Ws_old = Ws_new
	Ws_new = Tr_P + Tr_K

	delta_OEP = 0.0d0
	do i=1,1
		do j=1,1
			delta_OEP = delta_OEP+(P_new(i,j)-P_fci(i,j))**2.0d0
		end do
	end do

	delta_OEP = dsqrt(delta_OEP/(1*1))

	write(8,*) Step,delta_OEP,Ws_new

	write(2,*) 'Ws = ', Ws_new
	write(2,*) 'integral(V*(p-pin)) = ', Tr_P
	write(2,*) 'Kinetic Energy = ', Tr_K

	write(2,*) 'delta_OEP = ', delta_OEP

	write(2,*) "Initial Density Matrix"
	do i=1,norbit
		do j=1,norbit
			write(2, "(f16.8)", advance = "No") P_new(i,j)
		end do
		write(2,"(f16.8)")
	end do

	write(2,*) "Initial Hessian_inv"
	do i=1,ninter
		do j=1,ninter
			write(2, "(f16.8)", advance = "No") Hessian_inv(i,j)
		end do
		write(2,"(f16.8)")
	end do

	write(2,*) '-------------------------------'
	write(2,*) '-------------------------------'



!---------------------OEP Cycle---------------------------


	50 continue

	step = step +1

	!PE_t_new = F - KE_t

	write(2,*) "OEP Step:",step


!---------------------Calculating Gradient----------------

	Gradient_old = Gradient_new

!----------------------------------------------------------------------
!		Using Upper triangle to make sure the potential remains symmetric
!----------------------------------------------------------------------

	counti = 1
	do i=1,1
		do j=i,1
			Gradient_new(counti) = P_new(i,j) - P_fci(i,j)
			counti = counti + 1
		end do
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
!	direction_column = 0.0d0
!	call DGEMM('N','N',ninter,1,ninter,1,Hessian_inv,ninter,Grad_new_column,ninter,1,direction_column,ninter)


	direction_column = matmul(Hessian_inv,Grad_new_column)

!--------------- decreasing the rank of direction matrix--------

	counti = 1
	do i=1,ninter
		direction_1(counti) = direction_column(counti,1)
		counti = counti + 1
	end do

	write(2,*)  "direction_column:"
	do i=1,ninter
		write(2,*) direction_column(i,1)
	end do

	direction = 0.0d0
	counti = 1
	do i=1,1
		do j=i,1
			direction(i,j) = direction_column(counti,1)
			counti = counti + 1
		end do
	end do

	write(2,*) "Upper-Triangular direction matrix:"
	do i=1,norbit
		do j=1,norbit
			write(2, "(f16.8)", advance = "No") direction(i,j)
		end do
		write(2,"(f16.8)")
	end do

	PE_t_old = PE_t_new


!--------------Updating potential------------------------------
	stp = 1.0d0

	80 continue

	PE_t_new = PE_t_old

	do i=1,1
		do j=i,1
			PE_t_new(i,j) = PE_t_new(i,j) + stp*direction(i,j)
		end do
	end do

	do j=1,2
		do i=j+1,2
			PE_t_new(i,j) = PE_t_new(j,i)
		end do
	end do


		write(2,*) "Symmetric Potential Matrix:"
		do i=1,norbit
			do j=1,norbit
				write(2, "(f16.8)", advance = "No") PE_t_new(i,j)
			end do
			write(2,"(f16.8)")
		end do

			write(2,*) "Hessian_inv*Grad:"
			do i=1,norbit
				do j=1,norbit
				write(2, "(f16.8)", advance = "No") direction(i,j)
			end do
			write(2,"(f16.8)")
		end do
		write(2,"(f16.8)")


!-------------------------Vext * (P_new-P_fci)-----------------
!			del_P = P_new-P_fci
!			VP = 0.0d0
!			call DGEMM('N','N',norbit,norbit,norbit,alpha1,del_P,norbit,PE_t_new,norbit,beta1,VP,norbit)


!			KE_f = 0.0d0
!			call DGEMM('N','N',norbit,norbit,norbit,alpha1,P_new,norbit,KE_t,norbit,beta1,KE_f,norbit)

!			Tr_K = 0.0d0
!			do i=1,norbit
!				Tr_K = Tr_K + KE_f(i,i)
!			end do

!			Tr_P = 0.0d0
!			do i=1,norbit
!				Tr_P = Tr_P + VP(i,i)
!			end do



!			Ws_old = Ws_new
!			Ws_new = Tr_P + Tr_K


		F = KE_t + PE_t_new

		write(2,*) "Fock Matrix"
		do i=1,norbit
			do j=1,norbit
				write(2, "(f16.8)", advance = "No") F(i,j)
			end do
			write(2,"(f16.8)")
		end do
!----------------------Transform and Diagonalize---------

!----------------------G = F*X----------------------------

		G = 0.0d0
		call DGEMM('N','N',norbit,norbit,norbit,alpha1,F,norbit,X,norbit,beta1,G,norbit)

!----------------------F_Prime = XT*G = XT*F*X--------------------

		F_Prime = 0.0d0
		call DGEMM('T','N',norbit,norbit,norbit,alpha1,X,norbit,G,norbit,beta1,F_Prime,norbit)

		do i=1,norbit
			do j=1,norbit
				Eigen_C_prime(i,j) = F_Prime(i,j)
			end do
		end do

		write(2,*) "Transformed Fock Matrix"
		do i=1,norbit
			do j=1,norbit
				write(2, "(f16.8)", advance = "No") F_Prime(i,j)
			end do
			write(2,"(f16.8)")
		end do

!----------------------Eigenvalues and Eigenvectors of F_Prime---------

		Eigen_value = 0.0d0
		call DSYEV('V', 'U', norbit, Eigen_C_prime, norbit, Eigen_value, WORK, LWORK, INFO)

!---------------------- X*C_Prime = C---------------------------------

		Eigen_C = 0.0d0
		call DGEMM('N','N',norbit,norbit,norbit,alpha1,X,norbit,Eigen_C_prime,norbit,beta1,Eigen_C,norbit)


		!write(2,*) "Eigenvectors of F Matrix"
		!do i=1,norbit
		!	do j=1,norbit
		!		write(2, "(f16.8)", advance = "No") Eigen_C(i,j)
		!	end do
		!	write(2,"(f16.8)")
		!end do


!====================================================
!                Check the order of MO
!====================================================

		delta_den = 0.0d0
		minimum = 1.0d4

		do k=1,norbit

			P_new = 0.0d0
			do i=1,norbit
				do j=1,norbit
					P_new(i,j) = P_new(i,j) + 2*Eigen_C(i,k)*Eigen_C(j,k)

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
		do i=1,norbit
			write(2,*) delta_den(i)
		end do
		write(2,*) "********************"

	!	write(*,*) 'step=',step-1,'min-orbital=',ktemp

		write(2,*) "Old C Matrix"
		do i=1,norbit
		  do j=1,norbit
		    write(2, "(f16.10)", advance = "No") Eigen_C(i,j)
		  end do
		write(2,*)
		end do

!------------------Swap columns if there is any change in order-----

		if (ktemp /= 1) then

			70 continue

			do i=1,norbit
				Call Swap(Eigen_C(i,ktemp),Eigen_C(i,1))
			end do

			ktemp = ktemp-1
			if (ktemp /= 1) then
				goto 70
			end if
		end if

		write(2,*) "New C Matrix"
		do i=1,norbit
			do j=1,norbit
				write(2, "(f16.10)", advance = "No") Eigen_C(i,j)
			end do
		write(2,*)
		end do


		P_new = 0.0d0
		do i=1,norbit
			do j=1,norbit
				do k=1,nelec/2
					P_new(i,j) = P_new(i,j) + 2*Eigen_C(i,k)*Eigen_C(j,k)
				end do
			end do
		end do

		write(2,*) "New Density Matrix"
		do i=1,norbit
			do j=1,norbit
				write(2, "(f20.16)", advance = "No") P_new(i,j)
			end do
			write(2,"(f16.8)")
		end do

	!-------------------------Vext * (P_new-P_fci)-----------------
		del_P = P_new-P_fci
		VP = 0.0d0
		call DGEMM('N','N',norbit,norbit,norbit,alpha1,del_P,norbit,PE_t_new,norbit,beta1,VP,norbit)


		KE_f = 0.0d0
		call DGEMM('N','N',norbit,norbit,norbit,alpha1,P_new,norbit,KE_t,norbit,beta1,KE_f,norbit)

		Tr_K = 0.0d0
			do i=1,norbit
				Tr_K = Tr_K + KE_f(i,i)
			end do

		Tr_P = 0.0d0
		do i=1,norbit
			Tr_P = Tr_P + VP(i,i)
		end do

		Ws_old = Ws_new
		Ws_new = Tr_P + Tr_K

		if ( (Ws_new - Ws_old) .lt. 0.d0) then

			stp = 0.5d0*stp
			Ws_new = Ws_old

			if (stp >= 1.d-14) then
				goto 80
			else
				write(*,*) 'You are at the maximum!', Ws_new
				goto 100
			end if

		end if

		delta_OEP = 0.0d0
		do i=1,1
			do j=1,1
				delta_OEP = delta_OEP+(P_new(i,j)-P_fci(i,j))**2.0d0
			end do
		end do

		delta_OEP = dsqrt(delta_OEP/(1*1))

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


		log = 0

		do i=1,ninter
			dgra_column(i,1) = dgrad(i)
		  dgra_row(1,i) = dgrad(i)
	    direction_row(1,i) = direction_column(i,1)
	  end do

!-----------<dx|dg>---------------------------

		!alpha_bracket = 0.0d0
		!call DGEMM('N','N',1,1,ninter,1,direction_row,1,dgra_column,ninter,1,alpha_bracket,1)

		alpha_bracket = matmul(direction_row,dgra_column)
		!do i=1,ninter
		!	write(*,*) dgra_column(i,1)
		!end do
		!write(*,*) '******************'
		!do i=1,ninter
		!	write(*,*) direction_row(1,i)
		!end do
		!write(*,*) '******************'
		!alpha_bracket = direction_row*dgra_column

		alpha_inv = 1.d0/alpha_bracket(1,1)

		!write(*,*) step,'<dx|dg>=',alpha_bracket(1,1)

	!	if(dabs(alpha_bracket(1,1)).lt.1.d-14) then
	!		write(*,*) 'BFGS <dx|dg> too small',alpha_bracket(1,1)
	!	end if
!---------------------------------------------

		!*     form HInvEq(k)|dg>

		!hu_column = 0.0d0
		!call DGEMM('N','N',ninter,1,ninter,1,hinveq,ninter,dgra_column,ninter,1,hu_column,ninter)

		!hu = 0.d0
		!call DGEMV('N',ninter,ninter,1,hinveq,ninter,dgra,1,1,hu,1)

		hu_column = matmul(Hessian_inv,dgra_column)


		!*     form <dg|HInvEq(k)|dg>

		!  uhu = 0.0d0
		!  call DGEMM('N','N',1,1,ninter,1,dgra_row,1,hu_column,ninter,1,uhu,1)

		uhu = matmul(dgra_row,hu_column)

		write(*,*) 'uhu =', uhu(1,1)

	!	if (uhu(1,1) .gt. 1.d2) then
	!		log = 1
	!	end if

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
		temp2 = matmul(out_product_a,Hessian_inv)

	!  temp3 = 0.0d0
	!  call DGEMM('N','N',ninter,ninter,ninter,1,hinveq,ninter,out_product_b,ninter,1,temp3,ninter)
		temp3 = matmul(Hessian_inv,out_product_b)



		!	write(2,*) "old Hessian_inv"
		!	do i=1,ninter
		!		do j=1,ninter
		!			write(2, "(f16.8)", advance = "No") Hessian_inv(i,j)
		!		end do
		!		write(2,"(f16.8)")
		!	end do

!-------------Updated Hessian------------------

	Hessian_inv = Hessian_inv + (alpha_inv*(1.0d0 + uhu(1,1)*alpha_inv))*temp1 - alpha_inv*(temp2 + temp3)


!--------------Updating Hessian--------------------------


	!	call bfgs_hinv(Hessian_inv,ninter,direction_1,dgrad,log,hu_column)

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
!
!hu_column = 0.0d0
!call DGEMM('N','N',ninter,1,ninter,1,Hessian_inv,ninter,dgra_column,ninter,1,hu_column,ninter)

!	hu_column = matmul(Hessian_inv,dgra_column)

!		write(2,*) "hu_column:"
!		do i=1,ninter
!			write(2,*) hu_column(i,1)
!		end do

!		write(2,"(f16.8)")



	!	if (log == 1) then
	!		write(*,*) step,'uhu'
	!	end if


		!call bfgs_oep(Hessian_inv,ninter,af,direction_1,dgrad,hu)


		write(2,*) "Updated Hessian_inv"
		do i=1,ninter
			do j=1,ninter
				write(2, "(f16.8)", advance = "No") Hessian_inv(i,j)
			end do
			write(2,"(f16.8)")
		end do
		write(2,*) "****************************************"
		write(2,*) "  ************************************  "

		if ( delta_OEP .gt. OEP_thr .and. step .le. max_it) then

				goto 50

			else if (step .gt. max_it) then
				write(*,*) "It does not converge in ",max_it," steps!"
			else
				write(*,*) "OEP Calculation Converged!"
				write(2,*) "step=", step

			end if

	!==============================================================

	100 continue

  End Program

!====================================================
!                 END OF THE MAIN PROGRAM
!====================================================

!====================================================
!                   SUBROUTINES
!----------------------------------------------------

!====================================================
!                      INPUT
!====================================================

  Subroutine Input(natom,nelec,nbas,CI)
	Integer :: nbas,nelec,natom,CI


	write(*,*)
	write(*,*) "|--------------Full CI Calculations-------------|"
	write(*,*)
	write(*,*) "Please specify the number of atoms:"
	read(*,*) natom
	write(*,*) "Please specify the number of electrons:"
	read(*,*) nelec
	write(*,*) "Please specify the number of orbital for each atom:"
	read(*,*) nbas
	write(*,*) "Please specify the number of Excitations you are going to use:"
	read(*,*) CI

  return
  end



!====================================================
!                  DISTANCE
!====================================================

  Subroutine Distance(x1,y1,z1,x2,y2,z2,d)
	Double Precision :: d,x1,y1,z1,x2,y2,z2



!-------------Calculating distance between points "1" and "2"------------

	d = dsqrt( (x1-x2) **2.0d0+(y1-y2) **2.0d0+(z1-z2) **2.0d0)

  return
  end

!====================================================
!                  OVERAP INTEGRAL (UNNORMALIZED)
!====================================================

  Subroutine Overlap(a,b,d,S)
	Double Precision :: a,b,d,pi,S
	pi = 4.0d0*datan(1.0d0)

	S = ( ( pi/(a+b) )**1.5d0)*exp((-a*b/(a+b))*d**2.0d0)
  return
  end
!====================================================
!                  KINETIC PART (UNNORMALIZED)
!====================================================

  Subroutine Kinetic(a,b,d,KE)
	Double Precision :: a,b,d,KE,pi
	pi = 4.0d0*datan(1.0d0)

	KE = a*b/(a+b)*(3.0d0-2.0d0*a*b/(a+b)*d**2)*((pi/(a+b))**1.5d0)*exp(-a*b/(a+b)*d**2.0d0)
  return
  end

!====================================================
!                  F0 FUNCTION
!====================================================

  Subroutine F0(t,F_0)
  	Double Precision :: pi,t,F_0
	pi = 4.0d0*datan(1.0d0)

!----------prevent from dividing by zero-----------

	if (t .gt. 1.0d-6) then
	  F_0 = 0.5d0*(pi/t)**0.5d0 *erf(t**0.5d0)
	else
	  F_0 = 1.0d0-t/3.0d0
	end if
  return
  end

!====================================================
!                  POTENTIAL PART (UNNORMALIZED)
!====================================================

  Subroutine Potential(a,b,Z,d1,d2,PE)
	Double Precision :: a,b,d1,d2,F_0,PE,pi,Z
	pi = 4.0d0*datan(1.0d0)

	call F0((a+b)*d2**2.0d0,F_0)
	PE = -2.0d0*pi/(a+b)*Z*exp(-a*b/(a+b)*d1**2.0d0)*F_0
  return
  end
!====================================================
!              TWO ELECTRON INTEGRAL (UNNORMALIZED)
!====================================================

  Subroutine Integral_2e(a,b,c,d,r1,r2,r3,Integral)
	Double Precision :: a,b,c,d,F_0,Integral,pi,r1,r2,r3
	pi = 4.0d0*datan(1.0d0)

	call F0( (a+b)*(c+d)/(a+b+c+d)*r3**2.0d0,F_0 )
	Integral = 2.0d0*pi**2.5d0/((a+b)*(c+d)*(a+b+c+d)**0.5d0)*exp(-a*b/(a+b)*r1**2.0d0 - c*d/(c+d)*r2**2.0d0)*F_0
  return
  end



!====================================================
!                       POLYGON
!====================================================

  Subroutine Polygon(n,a)
	Double Precision :: a,pi,r,x,x_new,y,y_new,z
	Integer :: i,n
	pi = 4.0d0*datan(1.0d0)

	Open(Unit = 4, File = 'Input.xyz')

!--------------Calculating the radius of circumscribed circle-------

	r = dsqrt( (a**2.0d0)/( 2.0d0 *(1-dcos(2.0d0*pi/n)) ))

	x = r
	y = 0.0d0
	z = 0.0d0

	write(4,*) 'H',x,y,z

	do i=1,n-1
	  x_new = x*dcos(2.0d0*pi/n)-y*dsin(2.0d0*pi/n)
	  y_new = x*dsin(2.0d0*pi/n)+y*dcos(2.0d0*pi/n)
	  z = 0.0d0
	  write(4,*) 'H',x_new,y_new,z
	  x = x_new
	  y = y_new
	end do
	close(4)

  return
  end

	!====================================================
	!                       GRID
	!====================================================

	  Subroutine Grid(m,n,r)
		Double Precision :: r,x,y,z
		Integer :: i,j,m,n

		Open(Unit = 4, File = 'Input.xyz')

		z = 0.0d0

		do i=0,m-1
			do j=0,n-1
				x = j*r
				y = i*r
				write(4,*) 'H',x,y,z
			end do
		end do


		close(4)

	  return
	  end

!====================================================
!                       SWAP
!====================================================

	  Subroutine Swap(a,b)
		Double Precision :: a,b,temp

			temp = a
			a = b
			b = temp
	  return
		end


!====================================================
!                  END OF THE PROGRAM
!====================================================
