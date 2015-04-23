SUBROUTINE GET_DIM(nx, ny, nz, nparticles)
	
	IMPLICIT NONE
	
	INTEGER :: nx, ny, nz, nparticles, unit_part, pos
	CHARACTER(LEN = 128) :: filename
	
	unit_part = 1
	
	filename = './output_00003/part_00003.out00001'
	
	open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	read(unit_part) 
	read(unit_part) 
	read(unit_part) nparticles
	
! 	write(*,'(a i10 a)') 'There are', nparticles, ' particles'
	
	nx = nparticles ** (1.0/3.0)
	ny = nx
	nz = nx
	
! 	write(*,'(a i4 i4 i4)') 'Dimensions:', nx, ny, nz

! 	close(unit_part)
	
END SUBROUTINE GET_DIM