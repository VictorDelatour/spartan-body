SUBROUTINE WRITE_DENSITY(nx, ny, nz, density)
	
	IMPLICIT NONE
	
	INTEGER :: ierror, bov_id, density_id
	INTEGER, INTENT(in) :: nx, ny, nz
	REAL*8, DIMENSION(nx, ny, nz), INTENT(in) :: density
	
	bov_id = 7
	density_id = 8
	
	open(unit = bov_id, file = "./data/density_agnew.bov", status = "replace", action = "readwrite")
	
	write(bov_id,*) "TIME: 1"
	write(bov_id,*) "DATA_FILE: density_agnew.dat"
	write(bov_id,*) "DATA_SIZE: ", nx, ny, nz
	write(bov_id,*) "DATA_FORMAT: DOUBLE"
	write(bov_id,*) "VARIABLE: density"
	write(bov_id,*) "DATA_ENDIAN: LITTLE"
	write(bov_id,*) "CENTERING: zonal"
	write(bov_id,*) "BRICK_ORIGIN: ", 1, 1, 1
	write(bov_id,*) "BRICK_SIZE: ", nx, ny, nz
	
	close(unit = bov_id)
	
	open(unit = density_id, file = "./data/density_agnew.dat", status = "replace", action = "readwrite", iostat = ierror, &
	    &  form = "unformatted", access = "direct", recl = nx*ny*nz*8)
		
	write(unit = density_id, rec = 1) density
	
	close(unit = density_id)
	
END SUBROUTINE