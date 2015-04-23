SUBROUTINE READ_RAMSES
	
	USE DIMENSION, ONLY : nparticles, nx, ny, nz
	USE VECTOR, only	: x, y, z, mass
	
	IMPLICIT NONE
	
	INTEGER :: unit_part, pos
	CHARACTER(LEN = 128) :: filename
	
	unit_part = 1
	
	filename = './output_00003/part_00003.out00001'
	
	open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	read(unit_part) 
	read(unit_part) 
	read(unit_part) nparticles
	
	write(*,'(a i10 a)') 'There are', nparticles, ' particles'
	ALLOCATE( x(nparticles) )	
	ALLOCATE( y(nparticles) )	
	ALLOCATE( z(nparticles) )
	ALLOCATE( mass(nparticles) )
	
	nx = nparticles ** (1.0/3.0)
	ny = nx
	nz = nx
	
	do pos = 1,5
    	read(unit_part) !skip
	end do
!     read(unit_part) ! Skip
!     read(unit_part) ! Skip
!     read(unit_part) ! Skip
!     read(unit_part) ! Skip
	
	read(unit_part) x
	read(unit_part) y
	read(unit_part) z
	
	do pos = 1,3
		read(unit_part) mass ! Skip velocity
	end do
	
	read(unit_part) mass ! 
	

	
	write(*,'(a i4 i4 i4)') 'Dimensions:', nx, ny, nz
	
	x = x * (nx - 1) + 1.0
	y = y * (ny - 1) + 1.0
	z = z * (nz - 1) + 1.0
	
	! Maybe you should also read the mass later	
	
	close(unit_part)
	
	do pos = 1, 10
		write(*,'(a F7.2 F7.2 F7.2 e15.6)') 'Position:', x(pos), y(pos), z(pos), mass(pos)
	end do
	
	
	
END SUBROUTINE READ_RAMSES