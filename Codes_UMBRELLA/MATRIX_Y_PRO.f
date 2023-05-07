	SUBROUTINE MATRIX_Y_SUB
	
	implicit none

	include 'NOOFPARTICLES.inc'
	include 'POSITIONS.inc'
	INCLUDE 'MATRIX_Y.inc'

	integer Jnb
	integer i1,i2
	double precision D,DX,DY
	double precision DXS,DYS,ChiRSq


	double precision DETY
	double precision MAT11,MAT12,MAT21,MAT22

C	MATRIX-Y
	DO I1=1,N
	    MATY(I1,1,1) =0.0d0
	    MATY(I1,1,2) =0.0d0
	    MATY(I1,2,1) =0.0d0
	    MATY(I1,2,2) =0.0d0
	ENDDO

c	Jnb=1 to 6 = 6 nearest neighbors
	DO Jnb=1,NBONDS
	  I1=Part1(Jnb)
	  I2=Part2(Jnb)
	  
	    D = (inX(I2)-inX(I1))
	    DX = D - (INBOXX)*ANINT(D/INBOXX)
	    D = (inY(I2)-inY(I1))
	    DY = D - INBOXY*ANINT(D/INBOXY)

	    ChiRSQ = DX**2 + DY**2

	    MATY(I1,1,1) = MATY(I1,1,1) + DX*DX
	    MATY(I1,1,2) = MATY(I1,1,2) + DX*DY
	    MATY(I1,2,1) = MATY(I1,2,1) + DY*DX
	    MATY(I1,2,2) = MATY(I1,2,2) + DY*DY

	    DX=-DX
	    DY=-DY
	    MATY(I2,1,1) = MATY(I2,1,1) + DX*DX
	    MATY(I2,1,2) = MATY(I2,1,2) + DX*DY
	    MATY(I2,2,1) = MATY(I2,2,1) + DY*DX
	    MATY(I2,2,2) = MATY(I2,2,2) + DY*DY
	ENDDO


C	INVERSE CALCULATION	
	DO I1=1,N
	  MAT11=MATY(I1,1,1)
	  MAT12=MATY(I1,1,2)
	  MAT21=MATY(I1,2,1)
	  MAT22=MATY(I1,2,2)

	  DETY = MATY(I1,1,1)*MATY(I1,2,2) -  MATY(I1,1,2)*MATY(I1,2,1)
	  DETY = 1.0D0/DETY

	    MATY(I1,1,1) = MAT22*DETY
	    MATY(I1,1,2) = -MAT12*DETY
	    MATY(I1,2,1) = -MAT21*DETY
	    MATY(I1,2,2) = MAT11*DETY
	ENDDO

	return 
	end
