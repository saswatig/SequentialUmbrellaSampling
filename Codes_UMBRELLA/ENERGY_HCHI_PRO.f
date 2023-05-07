	SUBROUTINE ENERGY_HCHI_SUB(ENERGY_SW,V_HCHI,P_HCHI)
	IMPLICIT NONE

	INCLUDE 'NOOFPARTICLES.inc'
	INCLUDE 'POSITIONS.inc'
	INCLUDE 'MATRIX_Y.inc'
	INCLUDE 'MATRIX_X.inc'
	INCLUDE 'ENERGY_HCHI.inc'

	INTEGER ENERGY_SW
	DOUBLE PRECISION V_HCHI,P_HCHI

	DOUBLE PRECISION HCHI,FSUM

	double precision RSSQ,ChiRSq
	integer i1,i2
	integer i,j,Jnb,j1,j2
	double precision D,DX,DY
	double precision DXS,DYS
        double precision fxhi,fyhi

	HCHI=0.D0
	fsum=0.0d0

	IF(ENERGY_SW.EQ.0)THEN

	CALL MATRIX_Y_SUB
	CALL MATRIX_X_SUB(ENERGY_SW)
C	CALCULATE STRAIN TENSOR
	DO I1=1,N
	    EP(I1,1,1) = MATX(I1,1,1)*MATY(I1,1,1)
     &	   	+ MATX(I1,1,2)*MATY(I1,2,1) - 1.0D0

	    EP(I1,1,2) = MATX(I1,1,1)*MATY(I1,1,2)
     &		+ MATX(I1,1,2)*MATY(I1,2,2) - 0.0D0

	    EP(I1,2,1) = MATX(I1,2,1)*MATY(I1,1,1)
     &		+ MATX(I1,2,2)*MATY(I1,2,1) - 0.0D0

	    EP(I1,2,2) = MATX(I1,2,1)*MATY(I1,1,2)
     &		+ MATX(I1,2,2)*MATY(I1,2,2) - 1.0D0
	
	  CHI(I1)=0.0D0
	  Press_hChi(i1)=0.0d0
	ENDDO

c	LOOP OVER ALL MOLECULES
	DO 112 Jnb=1,nbonds
	  I=Part1(Jnb)
	  J=Part2(Jnb)

	  DXS=x(i)-x(j)
	  DYS=y(i)-y(j)
	  DXS = DXS - BOXX*ANINT(DXS/BOXX)
	  DYS = DYS - BOXY*ANINT(DYS/BOXY)
	  RSSQ=DXS**2+DYS**2
	
	  DX=inx(i)-inx(j)
	  DY=iny(i)-iny(j)
	  DX = DX - INBOXX*ANINT(DX/INBOXX)
	  DY = DY - INBOXY*ANINT(DY/INBOXY)
 	  ChiRSQ = DX**2 + DY**2


		  CHI(i) = CHI(i) + (DXS - (1.0D0+EP(I,1,1))*DX
     &		- EP(I,1,2)*DY)**2 + (DYS - EP(I,2,1)*DX
     &		- (1.0D0+EP(I,2,2))*DY)**2

	HCHI=HCHI+(DXS - (1.0D0+EP(I,1,1))*DX
     &		- EP(I,1,2)*DY)**2 + (DYS - EP(I,2,1)*DX
     &		- (1.0D0+EP(I,2,2))*DY)**2

          fxhi=4.0d0*DXS-2.0d0*(1.0d0+Ep(j,1,1))*DX
     &          -2.0d0*Ep(j,1,2)*DY
          fyhi=4.0d0*DYS-2.0d0*Ep(j,2,1)*DX
     &          -2.0d0*(1.0d0+Ep(j,2,2))*DY

        fsum=fsum+fxhi*DXS+fyhi*DYS
	Press_hChi(i)=Press_hChi(i)+fxhi*DXS+fyhi*DYS

	DX=-DX
	DY=-DY
	DXS=-DXS
	DYS=-DYS


		  CHI(j) = CHI(j) + (DXS - (1.0D0+EP(j,1,1))*DX
     &		- EP(j,1,2)*DY)**2 + (DYS - EP(j,2,1)*DX
     &		- (1.0D0+EP(j,2,2))*DY)**2

	HCHI=HCHI+(DXS - (1.0D0+EP(j,1,1))*DX
     &		- EP(j,1,2)*DY)**2 + (DYS - EP(j,2,1)*DX
     &		- (1.0D0+EP(j,2,2))*DY)**2

          fxhi=4.0d0*DXS-2.0d0*(1.0d0+Ep(i,1,1))*DX
     &          -2.0d0*Ep(i,1,2)*DY
          fyhi=4.0d0*DYS-2.0d0*Ep(i,2,1)*DX
     &          -2.0d0*(1.0d0+Ep(i,2,2))*DY

        fsum=fsum+fxhi*DXS+fyhi*DYS
	Press_hChi(j)=Press_hChi(j)+fxhi*DXS+fyhi*DYS

112	CONTINUE
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ELSE
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	CALL MATRIX_X_SUB(ENERGY_SW)
C	CALCULATE STRAIN TENSOR
	DO J2=1,7
	  if(j2.eq.7)then
	  i1=ENERGY_SW	
	  else
	  i1=nebors(ENERGY_SW,j2)
	  endif
	    EP(I1,1,1) = MATX(I1,1,1)*MATY(I1,1,1)
     &	   	+ MATX(I1,1,2)*MATY(I1,2,1) - 1.0D0

	    EP(I1,1,2) = MATX(I1,1,1)*MATY(I1,1,2)
     &		+ MATX(I1,1,2)*MATY(I1,2,2) - 0.0D0

	    EP(I1,2,1) = MATX(I1,2,1)*MATY(I1,1,1)
     &		+ MATX(I1,2,2)*MATY(I1,2,1) - 0.0D0

	    EP(I1,2,2) = MATX(I1,2,1)*MATY(I1,1,2)
     &		+ MATX(I1,2,2)*MATY(I1,2,2) - 1.0D0
            
	    CHI(I1)=0.0D0
	    Press_hChi(i1)=0.0d0
	ENDDO

C	CALCULATE CHI 

c	LOOP OVER ALL MOLECULES
	DO 113 Jnb=1,7
	  if(Jnb.eq.7)then
	  I=ENERGY_SW
	  else
	  I=nebors(ENERGY_SW,Jnb)
	  endif
	  Do J1=1,6
	  J=nebors(I,J1)

	  DXS=x(i)-x(j)
	  DYS=y(i)-y(j)
	  DXS = DXS - BOXX*ANINT(DXS/BOXX)
	  DYS = DYS - BOXY*ANINT(DYS/BOXY)
	  RSSQ=DXS**2+DYS**2
	
	  DX=inx(i)-inx(j)
	  DY=iny(i)-iny(j)
	  DX = DX - INBOXX*ANINT(DX/INBOXX)
	  DY = DY - INBOXY*ANINT(DY/INBOXY)
 	  ChiRSQ = DX**2 + DY**2


		  CHI(i) = CHI(i) + (DXS - (1.0D0+EP(I,1,1))*DX
     &		- EP(I,1,2)*DY)**2 + (DYS - EP(I,2,1)*DX
     &		- (1.0D0+EP(I,2,2))*DY)**2


CCCCC NEED EP(J,?,?) TO CALCULATE fxhi AND fyhi, EP VALUES RETAINED IN THE COMMON STATEMENT IN 'PCHI.inc' CCCCCCCCCCCCCCCCCC
          fxhi=4.0d0*DXS-2.0d0*(1.0d0+Ep(j,1,1))*DX
     &          -2.0d0*Ep(j,1,2)*DY
          fyhi=4.0d0*DYS-2.0d0*Ep(j,2,1)*DX
     &          -2.0d0*(1.0d0+Ep(j,2,2))*DY

	Press_hChi(i)=Press_hChi(i)+fxhi*DXS+fyhi*DYS

	EndDo
113	CONTINUE

	Do i1=1,n
	hchi=hchi+chi(i1)
	fsum=fsum+Press_hChi(i1)
	enddo


cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	V_HCHI=HCHI
        P_hChi=(0.50d0/24.0d0)*fsum/(boxx*boxy)
	RETURN
	END
