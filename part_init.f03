SUBROUTINE PART_INIT(nx, ny, nz, nparticles, x, y, z, mass)
	
	IMPLICIT NONE
	
	INTEGER :: nx, ny, nz, nparticles
	REAL*8, DIMENSION(nparticles) :: x, y, z, mass
	
	INTEGER :: unit_part, pos
	CHARACTER(LEN = 128) :: filename
	
	unit_part = 1
	
	filename = './output_00003/part_00003.out00001'
	
	open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	
	do pos = 1,5
    	read(unit_part) !skip
	end do
	
	read(unit_part) x
	read(unit_part) y
	read(unit_part) z
	
	do pos = 1,3
		read(unit_part) mass ! Skip velocity
	end do
	
	read(unit_part) mass !
	
! 	write(*,'(a i4 i4 i4)') 'Dimensions:', nx, ny, nz
		
	x = x * (nx - 1) + 1.0
	y = y * (ny - 1) + 1.0
	z = z * (nz - 1) + 1.0
	
	close(unit_part)
	
! 	do pos = 1, 10
! 		write(*,'(a F10.5 F10.5 F10.5 e15.6)') 'Position:', x(pos), y(pos), z(pos), mass(pos)
! 	end do
	
	
END SUBROUTINE PART_INIT