SUBROUTINE PART_INIT(nx, ny, nz, nparticles, nproc, x, y, z, vx, vy, vz, mass)
	
	!
	! Reads all the files containing the particles to import their position, velocity and mass
	!
	
	IMPLICIT NONE
	
	INTEGER :: nx, ny, nz, nparticles, nproc, nparticles_file, total_npart
	INTEGER :: skip, i
	INTEGER :: first, last
	INTEGER, DIMENSION(:), ALLOCATABLE :: skip_indices
	REAL*8, DIMENSION(nparticles) :: x, y, z, vx, vy, vz, mass 
	REAL*8, DIMENSION(:), ALLOCATABLE :: buffer
	
	INTEGER :: unit_part, unit_info, pos, index
	CHARACTER(LEN = 128) :: folder, filename,  char_buffer
	CHARACTER(LEN = 128) :: nx_string, pos_string, index_string
		
	write(nx_string, "(I0)") nx
	
	if (nx .eq. 512) then
		index = 6	
	else
		index = 3
	end if
	
	write(index_string, "(I0)")  index
		
	
	folder = "./output_00003/test" // TRIM(nx_string)
	
	first = 1
	last = 0
	
	do pos = 1, nproc
		
		write(pos_string, "(I0)") pos
		
		if (pos < 10) then
			filename = TRIM(folder) // "/part_0000" // TRIM(index_string) // ".out0000" // TRIM(pos_string)
		else if (pos < 100) then
			filename = TRIM(folder) // "/part_0000" // TRIM(index_string) // ".out000" // TRIM(pos_string)
		else 
			filename = TRIM(folder) // "/part_0000" // TRIM(index_string) // ".out00" // TRIM(pos_string)
		end if
		
		write(*,*) filename
		
		open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	
		read(unit_part)
		read(unit_part)
		read(unit_part) nparticles_file
		
		do i = 1,5
			read(unit_part)
		end do
		
		last = first + nparticles_file - 1

		read(unit_part) x(first:last)
		read(unit_part) y(first:last)
		read(unit_part) z(first:last)
		
		read(unit_part) vx(first:last)
		read(unit_part) vy(first:last)
		read(unit_part) vz(first:last)
		
		read(unit_part) mass(first:last)
		
		close(unit_part)
		
		first = last + 1
		
	end do
	
	
	x = x * (128 - 1) + 1.0
	y = y * (128 - 1) + 1.0
	z = z * (128 - 1) + 1.0
	
	
END SUBROUTINE PART_INIT