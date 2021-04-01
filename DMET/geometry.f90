  Subroutine input_geo
    USE array

!====================================================
!                  INPUT GEOMETRY
!====================================================
	40 continue

	write(*,*) 'Which Geometry are you interested in:'
	write(*,*) '1. Polygon'
	write(*,*) '2. Grid'
	write(*,*) '3. Line'
	read(*,*) geo

	if (geo.eq.1) then
		write(*,*) 'Please insert the number of Polygon side:'
		read(*,*) nside
		write(*,*) 'Please insert the Polygon side length(bond length):'
		read(*,*) a
	else if (geo.eq.2) then
		write(*,*) 'Please insert the number of rows:'
		read(*,*) row
		write(*,*) 'Please insert the number of columns:'
		read(*,*) column
		write(*,*) 'Please insert the bond length:'
		read(*,*) a
	else if (geo.eq.3) then
		write(*,*) 'Please insert the bond length(high-level calculation):'
		read(*,*) a
		write(*,*) 'Please insert the bond length(low-level calculation):'
		read(*,*) b
	else
		write(*,*) 'Input data is invalid!'
		goto 40
	end if


	if (geo.eq.1) then
  	call Polygon(nside,a)
  else if (geo.eq.2) then
		call Grid(row,column,a)
	else if (geo.eq.3) then
		call line(natom,a,b)
	end if

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
!                       LINE
!====================================================

		Subroutine Line(n,a,b)
		Double Precision :: a,b,x,y,z
		Integer :: i,n

		Open(Unit = 4, File = 'Input.xyz')

		x = 0.0d0
		y = 0.0d0
    z = 0.0d0

		write(4,*) 'H',x,y,z
		z = z + a
		write(4,*) 'H',x,y,z

		do i=1,n-2
			z = z + b
			write(4,*) 'H',x,y,z
		end do


		close(4)

		return
		end
