SUBROUTINE PART_INIT(nx, ny, nz, nparticles, x, y, z, vx, vy, vz, mass)
	
	IMPLICIT NONE
	
	INTEGER :: nx, ny, nz, nparticles, nparticles_file
	INTEGER :: skip, i
	INTEGER, DIMENSION(:), ALLOCATABLE :: skip_indices
	REAL*8, DIMENSION(nparticles) :: x, y, z, vx, vy, vz, mass 
	REAL*8, DIMENSION(:), ALLOCATABLE :: buffer
	
	INTEGER :: unit_part, pos
	CHARACTER(LEN = 128) :: filename
	
	unit_part = 1
	
	filename = './output_00003/part_00003.out00001'
	
	open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	
	read(unit_part)
	read(unit_part)
	read(unit_part) nparticles_file
	
	do pos = 1,5
    	read(unit_part) !skip
	end do
	
	if (nparticles_file /= nparticles) then
		
		ALLOCATE( buffer(nparticles_file) )
		
		skip = floor( real(nparticles_file) / real(nparticles) )
		
		ALLOCATE( skip_indices(nparticles) )
		skip_indices = [(i, i=0,nparticles-1, 1)] * skip + 1

		
		
		read(unit_part) buffer
		x = buffer( skip_indices )
! 		x = buffer( 1:nparticles )
		
		read(unit_part) buffer
		y = buffer( skip_indices )
! 		y = buffer( 1:nparticles )
		
		read(unit_part) buffer
		z = buffer( skip_indices )
! 		z = buffer( 1:nparticles )
		
		
		read(unit_part) buffer
		vx = buffer( skip_indices )
! 		vx = buffer( 1:nparticles )
		
		read(unit_part) buffer
		vy = buffer( skip_indices )
! 		vy = buffer( 1:nparticles )
		
		read(unit_part) buffer
		vz = buffer( skip_indices )
! 		vz = buffer( 1:nparticles )
		
		
		read(unit_part) buffer
		mass = buffer( skip_indices )
! 		mass = buffer( 1:nparticles )
		
	else
		
		read(unit_part) x
		read(unit_part) y
		read(unit_part) z
	
		read(unit_part) vx
		read(unit_part) vy
		read(unit_part) vz

		read(unit_part) mass !
		
	end if
	

	
	
	
! 	write(*,'(a i4 i4 i4)') 'Dimensions:', nx, ny, nz
		
	x = x * (nx - 1) + 1.0
	y = y * (ny - 1) + 1.0
	z = z * (nz - 1) + 1.0
	
	close(unit_part)
	
! 	do pos = 1, 10
! 		write(*,'(a F10.5 F10.5 F10.5 e15.6)') 'Position:', x(pos), y(pos), z(pos), mass(pos)
! 	end do
	
	
END SUBROUTINE PART_INIT