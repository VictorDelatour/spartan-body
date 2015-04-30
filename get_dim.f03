SUBROUTINE GET_DIM(nx, ny, nz, nparticles)
	
	IMPLICIT NONE
	
	INTEGER :: nx, ny, nz, nparticles, unit_info, min_level
	CHARACTER(LEN = 128) :: filename, buffer
	
	unit_info = 1
	
	filename = './output_00003/info_00003.txt'
	
	open(unit = unit_info, file = filename, status = 'old', form = 'formatted')

	read(unit_info, '(A13, I11)')
	read(unit_info, '(A13, I11)')
	read(unit_info, '(A13, I11)') buffer, min_level
	
	close(unit_info)

	nx = 2**min_level
	ny = nx
	nz = nx
	
END SUBROUTINE GET_DIM