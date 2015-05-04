SUBROUTINE GET_DIM(nx, ny, nz, nparticles)
	
	IMPLICIT NONE
	
	INTEGER :: nx, ny, nz, nparticles, unit_info, unit_part, min_level
	CHARACTER(LEN = 128) :: filename, buffer
	
	unit_info = 1
	unit_part = 10
	
	filename = './output_00003/info_00003.txt'
	
	open(unit = unit_info, file = filename, status = 'old', form = 'formatted')

	read(unit_info, '(A13, I11)')
	read(unit_info, '(A13, I11)')
	read(unit_info, '(A13, I11)') buffer, min_level
	
	close(unit_info)

	nx = 2**min_level
	ny = nx
	nz = nx
	
	filename = './output_00003/part_00003.out00001'
	
	open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	read(unit_part) 
	read(unit_part) 
	read(unit_part) nparticles
	close(unit_part)
	
END SUBROUTINE GET_DIM