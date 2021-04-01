!====================================================
!                      ALLOCATE
!====================================================

  Subroutine alloc_array
    use array


    write(*,*)
    write(*,*) "|--------------DMET Calculations-------------|"
    write(*,*)
    write(*,*) "Please specify the number of atoms:"
    read(*,*) natom
    write(*,*) "Please specify the number of electrons:"
    read(*,*) nelec
    write(*,*) "Please specify the number of atomic orbitals on each atom:"
    read(*,*) norbit

    pi = 4.0d0*datan(1.0d0)

    threshold_HF = 1.0d-11
    threshold_OEP = 1.0d-5
    E_thr = 1.0d-12

    max_it = 100
    max_bondlength = 2.1d0
    bond_incr = 0.1d0

    alpha1 = 1.0d0
    beta1 = 1.0d0

    atom_basis = 1
    n_fragment = 10
    mu = 1.0d6
    unlocalized = 4
    norbit_fci = 2
    nelec_A = 2
    bool = 1
    bath = 0

!====================================================
!                  READING BASIS SET
!====================================================

    norbit_t = norbit*natom

    Allocate(Sub_basis(norbit_t))
    Allocate(Z(natom))
    Allocate(alpha(norbit_t,30))
    Allocate(coeff(norbit_t,30))


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


    nprimtv_t = 0
    do i=1,norbit_t

    	read(3,*) atomname,Sub_basis(i)

    	do j=1,Sub_basis(i)
    	  read (3,*) alpha(i,j), coeff(i,j)
    	end do

    	nprimtv_t = nprimtv_t + Sub_basis(i)
    	end do

  	close(3)

    write(*,*) "total primitives=",nprimtv_t

!====================================================
!         SCALING of Coefficients & Exponents
!----------------------------------------------------

    do i=1,norbit_t
      do j=1,Sub_basis(i)
    !    alpha(i,j) = alpha(i,j)*(zeta(i)**2.0d0)
    	  coeff(i,j) = coeff(i,j)*((2.0d0*alpha(i,j)/pi)**0.75d0)
    	end do
    end do



  lwork = 3*norbit_t
  lwork1 = norbit_t
  lwork2 = norbit_t


!	call Combinatory(nelec,2,mat_dim_temp)
!	call Combinatory(2*norbit_fci-nelec,2,mat_dim)

  mat_dim = 2*norbit_fci-1 + (norbit_fci-1)*(norbit_fci-2)/2

  write(*,*) "FCI H Dimension:",mat_dim,"x",mat_dim

  LWORK_FCI = 3*mat_dim

  ninter = norbit_t*(norbit_t+1)/2
  ninter = norbit_t**2


  Allocate(S_M(norbit_t,norbit_t))
  Allocate(S_inv(norbit_t,norbit_t))

  Allocate(S_M_T(norbit_t,norbit_t))
  Allocate(S_A(norbit,norbit))
  Allocate(S_half(norbit_t,norbit_t))
  Allocate(S_half_inv(norbit_t,norbit_t))
  Allocate(S_half_temp(norbit_t,norbit_t))
  Allocate(S_half_diag(norbit_t,norbit_t))
  Allocate(S_prime(norbit,norbit_t))
  Allocate(S_prime_T(norbit_t,norbit))

  Allocate(U(norbit_t,norbit_t))
  Allocate(H_t(mat_dim,mat_dim))
  Allocate(H_core(norbit_t,norbit_t))
  Allocate(H_core_orth(norbit_t,norbit_t))
  Allocate(H_core_new(norbit_t,norbit_t))
  Allocate(KE_t(norbit_t,norbit_t))
  Allocate(KE_t_temp(norbit_t,norbit_t))
  Allocate(KE_t_prime(norbit_t,norbit_t))
  Allocate(KE_f(norbit_t,norbit_t))

  Allocate(P_new_orth(norbit_t,norbit_t))
  Allocate(P_orth(norbit_t,norbit_t))
  Allocate(P_temp(norbit_t,norbit_t))

  Allocate(KE_new(norbit_t,norbit_t))
  Allocate(PE_t(norbit_t,norbit_t))
  Allocate(PE_t_init(norbit_t,norbit_t))
  Allocate(PE_t_new(norbit_t,norbit_t))
  Allocate(PE_t_final(norbit_t,norbit_t))
  Allocate(PE_t_old(norbit_t,norbit_t))

  Allocate(P_new(norbit_t,norbit_t))
  Allocate(P_fci(norbit_t,norbit_t))
  Allocate(P_fci_A(norbit_t,norbit_t))
  Allocate(P_fci_A_mo(norbit_t,norbit_t))
  Allocate(P_fci_A_temp(norbit_t,norbit_t))
  Allocate(P_HF_A(norbit_t,norbit_t))

  Allocate(delta_den(norbit_t))
  Allocate(P_HF(norbit_t,norbit_t,nelec_A/2))


  Allocate(P_garnet(norbit_t,norbit_t))
  Allocate(P_garnet_temp(norbit_t,norbit_t))
  Allocate(P_garnet_nat(norbit_t,norbit_t))
  Allocate(eigen_C_nat(norbit_t,norbit_t))
  Allocate(eigen_C_unnat(norbit_t,norbit_t))
  Allocate(P_garnet_evalue_un(norbit_t))
  Allocate(P_garnet_evalue(norbit_t))
  Allocate(Project_garnet_temp(norbit_t,norbit_t))
  Allocate(Project_garnet(norbit_t,norbit_t))
  Allocate(Gamma_garnet(norbit_t,norbit_t))
  Allocate(eigen_c_nat_ao(norbit_t,norbit_t))
  Allocate(eigen_c_unnat_ao(norbit_t,norbit_t))
  Allocate(P_old(norbit_t,norbit_t))
  Allocate(Project_temp(norbit_t,norbit))

  Allocate(Project_temp_A(norbit_t,norbit))
  Allocate(Project_new(norbit_t,norbit_t))
  Allocate(Project_new_temp(norbit_t,norbit_t))

  Allocate(Gama(nelec/2,nelec/2))
  Allocate(Eigen_Gama(nelec/2,nelec/2))
  Allocate(Gama_diag(nelec/2,nelec/2))
  Allocate(Gama_temp(nelec/2,norbit_t))
  Allocate(Gama_B(norbit_t,norbit_t))
  Allocate(P1(norbit_t,norbit_t))
  Allocate(P2(norbit_t,norbit_t))

  Allocate(P_in(norbit_t,norbit_t))
  Allocate(P_check(norbit_t,norbit_t))
  Allocate(PS_check(norbit_t,norbit_t))

  Allocate(Slater_coeff(norbit_fci,norbit_fci))

  Allocate(direction(norbit_t,norbit_t))
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


  Allocate(del_P(norbit_t,norbit_t))
  Allocate(VP(norbit_t,norbit_t))
  Allocate(F_Prime(norbit_t,norbit_t))
  Allocate(F_A_Prime(norbit_t,norbit_t))
  Allocate(F(norbit_t,norbit_t))
  Allocate(F_A(norbit_t,norbit_t))
  Allocate(Eigen_F_A(norbit_t,norbit_t))
  Allocate(Eigen_F_A_inv(norbit_t,norbit_t))

  Allocate(Eigen_F_A_prime(norbit_t,norbit_t))
  Allocate(G(norbit_t,norbit_t))
  Allocate(X(norbit_t,norbit_t))
  Allocate(X_inv(norbit_t,norbit_t))
  Allocate(Eigen_C_Prime(norbit_t,norbit_t))
  Allocate(Eigen_C_H_t(mat_dim,mat_dim))
  Allocate(E_CI(mat_dim))
  Allocate(Eigen_C(norbit_t,norbit_t))
  Allocate(Eigen_FCI(norbit_t,norbit_fci))
  Allocate(eigen_c_occ(norbit_t,nelec/2))
  Allocate(eigen_c_occ_T(nelec/2,norbit_t))
  Allocate(eigen_c_occ_local(norbit_t,nelec/2))
  Allocate(M_2e(norbit_t,norbit_t,norbit_t,norbit_t))
  Allocate(M_2e_new(norbit_t,norbit_t,norbit_t,norbit_t))
  Allocate(Eigen_value(norbit_t))
  Allocate(Eigen_value_half(norbit_t))
  Allocate(E_value(nelec/2))
  Allocate(E_value_F_A(norbit_t))
  Allocate(Work(lwork))
  Allocate(Work1(lwork1))
  Allocate(Work2(lwork2))
  Allocate(Work_FCI(LWORK_FCI))
  Allocate(IPIV(norbit))
  Allocate(IPIV1(norbit_t))
  Allocate(IPIV2(norbit_t))

  Allocate(Zeta(natom))
  Allocate(R(natom,3))


  return
  end
