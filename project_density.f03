SUBROUTINE PROJECT_DENSITY(nx, ny, nz, nparticles, x, y, z, mass, density)
	
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: nx, ny, nz, nparticles
	REAL*8, DIMENSION(nparticles) :: x, y, z, mass
	REAL*8, DIMENSION(nx, ny, nz) :: density ! Maybe be careful with dimensions, fortran is column-major
	
	REAL*8 :: didx, didy, didz, didxp,  didyp,  didzp
	REAL*8 :: d1, d2, d3, d4, d5, d6, d7, d8
	
	INTEGER :: particle
	INTEGER :: idx, idy, idz, idxp, idyp, idzp
	INTEGER :: ierror
	
	mass(:) = 1.0
	
	do particle = 1, nparticles
		
		idx = int( x(particle) + 0.5 )
		idy = int( y(particle) + 0.5 )
		idz = int( z(particle) + 0.5 )
		
		if(idx .eq. 0) idx = nx
		if(idy .eq. 0) idy = ny
		if(idz .eq. 0) idz = nz
		
		idxp = idx + 1
		idyp = idy + 1
		idzp = idz + 1
		
		if(idx .eq. nx) idxp = 1
		if(idy .eq. ny) idyp = 1
		if(idz .eq. nz) idzp = 1
		
		! Computes the distance from the middle point of the cell?
		didx = x(particle) - float(idx) + 0.5
		didy = y(particle) - float(idy) + 0.5
		didz = z(particle) - float(idz) + 0.5
		
		didxp = 1.0 - didx
		didyp = 1.0 - didy
		didzp = 1.0 - didz
				
		! Computes the weights for trilinear interpolation
		d1 = didx 	* didy 	* didz
		d2 = didxp 	* didy 	* didz
		d3 = didx 	* didyp * didz
		d4 = didxp 	* didyp * didz
		d5 = didx 	* didy 	* didzp
		d6 = didxp 	* didy 	* didzp
		d7 = didx 	* didyp * didzp
		d8 = didxp 	* didyp * didzp
		
		
		density(idx,  idy,  idz)  = density(idx,  idy,  idz)  + d1 * mass(particle)
		density(idxp, idy,  idz)  = density(idxp, idy,  idz)  + d2 * mass(particle)
		density(idx,  idyp, idz)  = density(idx,  idyp, idz)  + d3 * mass(particle)
		density(idxp, idyp, idz)  = density(idxp, idyp, idz)  + d4 * mass(particle)
		density(idx,  idy,  idzp) = density(idx,  idy,  idzp) + d5 * mass(particle)
		density(idxp, idy,  idzp) = density(idxp, idy,  idzp) + d6 * mass(particle)
		density(idx,  idyp, idzp) = density(idx,  idyp, idzp) + d7 * mass(particle)
		density(idxp, idyp, idzp) = density(idxp, idyp, idzp) + d8 * mass(particle)
		
	end do
	
	CALL WRITE_DENSITY(nx, ny, nz, density)
	
END SUBROUTINE PROJECT_DENSITY