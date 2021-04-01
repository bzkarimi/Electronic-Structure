  Subroutine all_integrals
    USE array

!=====================================================
! Integral calculations are done based on the book,
! "modern quantum chemistry" by Szabu and Ostland.
! to better undersand the notation please refer to chapters
! 2 & 3 of this book and appendices A and B.
!=====================================================

!====================================================
!                  OVERLAP INTEGRAL
!----------------------------------------------------

	open(Unit = 2,File ="Output.txt" )
	write(2,*) title

!-------------loop over atoms------------------------
	do i=1,norbit_t
		iatom = int((i-1)/norbit)+1

	  do j=1,norbit_t
      S_M(i,j) = 0.0d0

			jatom = int((j-1)/norbit)+1

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

	H_core = 0.0d0
!-------------loop over atoms------------------------
	do i=1,norbit_t
		iatom = int((i-1)/norbit)+1

	  do j=1,norbit_t
			KE_t(i,j) = 0.0d0
			jatom = int((j-1)/norbit)+1

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

	do i=1,norbit_t
		iatom = int((i-1)/norbit)+1

	  do j=1,norbit_t
			PE_t(i,j) = 0.0d0
			jatom = int((j-1)/norbit)+1


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
!----------------------------------------------------
			H_core(i,j) = H_core(i,j) + PE_t(i,j)

	  end do
	end do


!====================================================
!                 2e-Integrals
!----------------------------------------------------

  M_2e = 0.0d0

  do i=1,norbit_t
		iatom = int((i-1)/norbit)+1

    do j=1,norbit_t
			jatom = int((j-1)/norbit)+1

      do k=1,norbit_t
				katom = int((k-1)/norbit)+1

        do l=1,norbit_t
					latom = int((l-1)/norbit)+1


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



  !---------------------Overlap Matrix------------------
  		write(2,*) "Overlap Matrix S"
  		do i=1,norbit_t
  			do j=1,norbit_t
  			  write(2, "(f16.8)", advance = "No") S_M(i,j)
  			end do
  			write(2,*)
  		end do
  		write(2,*) "****************************************"
  		write(2,*) "  ************************************  "

  !---------------------Kinetic Matrix------------------
  		write(2,*) "Kinetic Matrix"

  		do i=1,norbit_t
  		  do j=1,norbit_t
  		    write(2, "(f16.8)", advance = "No") KE_t(i,j)
  		  end do
  		write(2,*)
  		end do
  		write(2,*) "****************************************"
  		write(2,*) "  ************************************  "

  !---------------------Potential Matrix-----------------

  		write(2,*) "Potential Matrix"
  		do i=1,norbit_t
  		  do j=1,norbit_t
  		    write(2, "(f16.8)", advance = "No") PE_t(i,j)
  			end do
  		write(2,*)
  	end do
  	write(2,*) "****************************************"
  		write(2,*) "  ************************************  "

!===========================================================================
!                   WRITING THE WHOLE 2e-MATRIX
!---------------------------------------------------------------------------

!	write(2,*) "All 2e Matrix"
!	do i=1,norbit_t
!	  do j=1,norbit_t
!	    do k=1,norbit_t
!	      do l=1,norbit_t
!	        write(2, "(A1,I3,I3,A1,I3,I3,A1,f16.8)") '(',i,j,'|',k,l,')',M_2e(i,j,k,l)
!	      end do
!	    end do
!	  end do
!	end do
!
!	write(2,*) "**********************************************"


!---------------------Core Hamiltonian-----------------

!	write(2,*) "H_core Matrix"
!	do i=1,norbit_t
!	  do j=1,norbit_t
!	    write(2, "(f16.8)", advance = "No") H_core(i,j)
!	  end do
!	write(2,*)
!	end do



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
