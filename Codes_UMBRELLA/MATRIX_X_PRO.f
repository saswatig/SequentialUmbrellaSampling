	SUBROUTINE MATRIX_X_SUB(ENERGY_SW)
	IMPLICIT NONE 

	include 'NOOFPARTICLES.inc'	
	INCLUDE 'POSITIONS.inc'
	INCLUDE 'MATRIX_X.inc'

	INTEGER ENERGY_SW

	INTEGER I,J,I1,JNB,J1
	DOUBLE PRECISION RXIJ,RYIJ
	DOUBLE PRECISION DX,DY,DXS,DYS,rijsq,D

	IF(ENERGY_SW.EQ.0)THEN
	  DO I1=1,N
	  MATX(I1,1,1)=0.0D0
	  MATX(I1,1,2)=0.0D0
	  MATX(I1,2,1)=0.0D0
	  MATX(I1,2,2)=0.0D0
	  ENDDO
CCCCCCCCCCCCCCCCCCCCCCC
C	MATRIX X
cccccccccccccccc
c	LOOP OVER ALL MOLECULES
	DO 115 Jnb=1,nbonds
	  I=Part1(Jnb)
	  J=Part2(Jnb)
	  rxij=x(i)-x(j)
	  ryij=y(i)-y(j)

C	PERIODIC BOUNDARY CONDITION
	  rxij=rxij-boxx*anint(rxij/boxx)
          ryij=ryij-boxy*anint(ryij/boxy)	

	  rijsq=rxij*rxij+ryij*ryij

	    D = inX(I)-inX(J)
	    DX = D - INBOXX*ANINT(D/INBOXX)
	    D = inY(I)-inY(J)
	    DY = D - INBOXY*ANINT(D/INBOXY)


		MATX(I,1,1) = MATX(I,1,1) + rxij*DX
		MATX(I,1,2) = MATX(I,1,2) + rxij*DY
	    	MATX(I,2,1) = MATX(I,2,1) + ryij*DX
	    	MATX(I,2,2) = MATX(I,2,2) + ryij*DY

	DX=-DX
	DY=-DY
	DXS=-rxij
	DYS=-ryij

		MATX(J,1,1) = MATX(J,1,1) + DXS*DX
		MATX(J,1,2) = MATX(J,1,2) + DXS*DY
	    	MATX(J,2,1) = MATX(J,2,1) + DYS*DX
	    	MATX(J,2,2) = MATX(J,2,2) + DYS*DY
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
115	CONTINUE
	ELSE
	  DO I1=1,7
	  iF(I1.EQ.7)then
	  J1=ENERGY_SW
	  else
	  J1=NEBORS(ENERGY_SW,I1)
	  endif	  
	  MATX(J1,1,1)=0.0D0
	  MATX(J1,1,2)=0.0D0
	  MATX(J1,2,1)=0.0D0
	  MATX(J1,2,2)=0.0D0
	  ENDDO
CCCCCCCCCCCCCCCCCCCCCCC
C	MATRIX X
cccccccccccccccc
c	LOOP OVER ALL MOLECULES
	DO 116 Jnb=1,7
	  iF(Jnb.EQ.7)then
	  I=ENERGY_SW
	  else
	  I=nebors(ENERGY_SW,Jnb)
	  endif
	  do j1=1,6
	  J=nebors(I,j1)
	  rxij=x(i)-x(j)
	  ryij=y(i)-y(j)

C	PERIODIC BOUNDARY CONDITION
	  rxij=rxij-boxx*anint(rxij/boxx)
          ryij=ryij-boxy*anint(ryij/boxy)	

	  rijsq=rxij*rxij+ryij*ryij

	    D = inX(I)-inX(J)
	    DX = D - INBOXX*ANINT(D/INBOXX)
	    D = inY(I)-inY(J)
	    DY = D - INBOXY*ANINT(D/INBOXY)


		MATX(I,1,1) = MATX(I,1,1) + rxij*DX
		MATX(I,1,2) = MATX(I,1,2) + rxij*DY
	    	MATX(I,2,1) = MATX(I,2,1) + ryij*DX
	    	MATX(I,2,2) = MATX(I,2,2) + ryij*DY

	enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
116	CONTINUE
	ENDIF
	
	RETURN
	END
