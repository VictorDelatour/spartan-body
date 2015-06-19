SUBROUTINE PART_INIT(nx, ny, nz, nparticles, nproc, x, y, z, vx, vy, vz, mass)
	
	IMPLICIT NONE
	
	INTEGER :: nx, ny, nz, nparticles, nproc, nparticles_file, total_npart
	INTEGER :: skip, i
	INTEGER :: first, last
	INTEGER, DIMENSION(:), ALLOCATABLE :: skip_indices
	REAL*8, DIMENSION(nparticles) :: x, y, z, vx, vy, vz, mass 
	REAL*8, DIMENSION(:), ALLOCATABLE :: buffer
	
	INTEGER :: unit_part, unit_info, pos
	CHARACTER(LEN = 128) :: folder, filename,  char_buffer
	CHARACTER(LEN = 128) :: nx_string, pos_string
		
	write(nx_string, "(I0)") nx
	folder = "./output_00003/test" // TRIM(nx_string)
	
	first = 1
	last = 0
	
	do pos = 1, nproc
		
		write(pos_string, "(I0)") pos
		filename = TRIM(folder) // "/part_00003.out0000" // TRIM(pos_string)
		write(*,*) filename
		
		open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	
		read(unit_part)
		read(unit_part)
		read(unit_part) nparticles_file
		
		do i = 1,5
			read(unit_part)
		end do
		
		last = first + nparticles_file - 1
		
! 		write(*,*) first, last

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
	
! 	filename = './output_00003/part_00003.out00001'
!
! 	open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
!
! 	read(unit_part)
! 	read(unit_part)
! 	read(unit_part) nparticles_file
!
! 	write(*,'(a i10)') "Particles in file", nparticles_file
!
! 	do pos = 1,5
!     	read(unit_part) !skip
! 	end do
! 
! 	read(unit_part) x
! 	read(unit_part) y
! 	read(unit_part) z
!
! 	read(unit_part) vx
! 	read(unit_part) vy
! 	read(unit_part) vz
!
! 	read(unit_part) mass !
!
	x = x * (nx - 1) + 1.0
	y = y * (ny - 1) + 1.0
	z = z * (nz - 1) + 1.0
	
! 	close(unit_part)
!
	
	
END SUBROUTINE PART_INIT