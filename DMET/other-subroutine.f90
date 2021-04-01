
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
!                  Kronecker Delta
!====================================================

  	Subroutine Kronecker(m,n,p)

  		Integer :: m,n,p

  		if (m.eq.n) then
  			p = 1
  		else
  			p = 0
  		end if

  	  return
  	  end
!====================================================
!                 Factorial
!====================================================

  		Subroutine Combinatory(m,n,p)

  			Integer :: m,n,p,i,ans

  	!-------------m!-------------------------
  			ans = 1
  			do i=1,m
  				ans = ans*i
  			end do
  			p = ans
  	!-------------n!-------------------------
  			ans = 1
  			do i=1,n
  				ans = ans*i
  			end do
  			p = p/ans
  	!-------------(m-n)!-------------------------
  			ans = 1
  			do i=1,m-n
  				ans = ans*i
  			end do
  			p = p/ans

  			return
  			end
