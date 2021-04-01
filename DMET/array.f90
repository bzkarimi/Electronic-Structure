module array

  Integer :: column, geo, i, itemp, j, jtemp, k, ktemp, l, ltemp,ninter
  Integer :: norbit_t, norbit_fci, norbit, n_fragment, nelec_A, nelec
  Integer :: m, n, o, p, nside, max_it, natom, nbas, ifragment
  Integer :: row, step, lwork, info,info1,info2, atom_basis, unlocalized,lwork1
  Integer :: LWORK_FCI, countenv, mat_dim, mat_dim_temp,nprimtv_t,lwork2
  Integer :: iatom,jatom,katom,latom,counti,countj,stp, bool, bath
  Integer,Allocatable,Dimension(:) :: IPIV,Sub_basis,IPIV1,IPIV2

  Double Precision :: a, alpha1, b, beta1, bond_incr,d, d1, d2, d3
  Double Precision :: E_corr,EN, EN_t,Integral, threshold_HF, threshold_OEP
  Double Precision :: Ws_new, Ws_old, xp, yp, zp, xq, yq, zq, KE
  Double Precision :: Tr_K_fci, Tr_P, Tr_K, stp_hf, delta_OEP, delta_HF
  Double Precision :: EN_old, EN_new, mu, alpha_inv, max_bondlength
  Double Precision :: N_elec,PE, pi,S, Tr_check, minimum, delta_E, E_thr

  Double Precision,Allocatable,Dimension(:) :: E_CI,dgrad,delta_den,dgra
  Double Precision,Allocatable,Dimension(:) :: Gradient_new,Gradient_old,hu,af
  Double Precision,Allocatable,Dimension(:) :: P_garnet_evalue_un,Eigen_value_half
  Double Precision,Allocatable,Dimension(:) :: P_garnet_evalue,Work_FCI,Z
  Double Precision,Allocatable,Dimension(:) :: E_value_F_A,Eigen_value,E_value
  Double Precision,Allocatable,Dimension(:) :: Zeta,Work,direction_1,Work1,work2

  Double Precision,Allocatable,Dimension(:,:) :: alpha,coeff,Eigen_C_H_t
  Double Precision,Allocatable,Dimension(:,:) :: eigen_c_occ_local,eigen_c_occ
  Double Precision,Allocatable,Dimension(:,:) :: eigen_c_occ_T,Eigen_C,Eigen_FCI
  Double Precision,Allocatable,Dimension(:,:) :: Eigen_C_Prime,Eigen_F_A
  Double Precision,Allocatable,Dimension(:,:) :: Eigen_F_A_prime, Eigen_Gama,KE_f
  Double Precision,Allocatable,Dimension(:,:) :: F_A,F,Project_new_temp,H_core_orth
  Double Precision,Allocatable,Dimension(:,:) :: del_P,G,Gama_B,Gama,S_inv
  Double Precision,Allocatable,Dimension(:,:) :: P_fci,P_fci_A,P_HF_A,KE_t_prime
  Double Precision,Allocatable,Dimension(:,:) :: Gama_temp,Gama_t, H_t,KE_t_temp
  Double Precision,Allocatable,Dimension(:,:) :: S_M,P_check,PS_check,R,X_inv
  Double Precision,Allocatable,Dimension(:,:) :: Project_temp_A,Project_temp
  Double Precision,Allocatable,Dimension(:,:) :: X,KE_new,U,PE_t,PE_t_new
  Double Precision,Allocatable,Dimension(:,:) :: H_core,S_A,F_Prime,H_core_new
  Double Precision,Allocatable,Dimension(:,:) :: VP,S_prime,P1,P2,S_M_T,KE_t
  Double Precision,Allocatable,Dimension(:,:) :: P_new,P_old,P_in,S_half_temp
  Double Precision,Allocatable,Dimension(:,:) :: direction,Project_new
  Double Precision,Allocatable,Dimension(:,:) :: P_fci_A_mo,P_fci_A_temp
  Double Precision,Allocatable,Dimension(:,:) :: Slater_coeff,P_garnet,F_A_Prime
  Double Precision,Allocatable,Dimension(:,:) :: Project_garnet_temp,Project_garnet
  Double Precision,Allocatable,Dimension(:,:) :: S_half,Gamma_garnet
  Double Precision,Allocatable,Dimension(:,:) :: P_garnet_temp,S_half_diag
  Double Precision,Allocatable,Dimension(:,:) :: Grad_new_column,direction_column
  Double Precision,Allocatable,Dimension(:,:) :: dgra_row, PE_t_old
  Double Precision,Allocatable,Dimension(:,:) :: direction_row,PE_t_final
  Double Precision,Allocatable,Dimension(:,:) :: PE_t_init,Gama_diag
  Double Precision,Allocatable,Dimension(:,:) :: dgra_column,eigen_C_nat
  Double Precision,Allocatable,Dimension(:,:) :: Hessian_inv,S_prime_T
  Double Precision,Allocatable,Dimension(:,:) :: eigen_C_unnat,eigen_c_unnat_ao
  Double Precision,Allocatable,Dimension(:,:) :: eigen_c_nat_ao,uhu,Eigen_F_A_inv
  Double Precision,Allocatable,Dimension(:,:) :: S_half_inv,P_garnet_nat
  Double Precision,Allocatable,Dimension(:,:) :: hu_row,alpha_bracket,hu_column
  Double Precision,Allocatable,Dimension(:,:) :: out_product_a,out_product_b
  Double Precision,Allocatable,Dimension(:,:) :: temp1,temp2,temp3
  Double Precision,Allocatable,Dimension(:,:) :: P_temp,P_orth,P_new_orth

  Double Precision,Allocatable,Dimension(:,:,:) :: P_HF

  Double Precision,Allocatable,Dimension(:,:,:,:) :: M_2e,M_2e_new

  Character(70) :: title,atomname


end module array
