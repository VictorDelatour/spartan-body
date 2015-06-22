SUBROUTINE GET_DIM(nx, ny, nz, nparticles, nproc)
	
	IMPLICIT NONE
	
	INTEGER :: nproc, pos, index
	INTEGER :: nx, ny, nz, nparticles, nparticles_file
	INTEGER :: unit_info, unit_part, min_level
	
	CHARACTER(LEN = 128) :: folder, filename
	CHARACTER(LEN = 128) :: buffer
	CHARACTER(LEN = 128) :: nx_string, pos_string, index_string
	
	unit_info = 1
	unit_part = 10	
	nparticles = 0
	
	write(nx_string, "(I0)") nx
	
	if (nx .eq. 512) then
		index = 6	
	else
		index = 3
	end if
	
	write(index_string, "(I0)") index
	
	
	
	folder = "./output_00003/test" // TRIM(nx_string)
	filename = TRIM(folder) // "/info_0000" // TRIM(index_string) // ".txt"
	
	open(unit = unit_info, file = filename, status = 'old', form = 'formatted')
	read(unit_info, '(A13, I11)') buffer, nproc
	close(unit_info)
	
	do pos = 1, nproc
		
		write(pos_string, "(I0)") pos
! 		filename = TRIM(folder) // "/part_00003.out0000" // TRIM(pos_string)

		if (pos < 10) then
			filename = TRIM(folder) // "/part_0000" // TRIM(index_string) // ".out0000" // TRIM(pos_string)
		else if (pos < 100) then
			filename = TRIM(folder) // "/part_0000" // TRIM(index_string) // ".out000" // TRIM(pos_string)
		else 
			filename = TRIM(folder) // "/part_0000" // TRIM(index_string) // ".out00" // TRIM(pos_string)
		end if
		
		open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	
		read(unit_part)
		read(unit_part)
		read(unit_part) nparticles_file
		
		nparticles = nparticles + nparticles_file
		
		close(unit_part)
		
	end do

	

	
END SUBROUTINE GET_DIM