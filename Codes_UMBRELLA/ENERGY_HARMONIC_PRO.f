	SUBROUTINE ENERGY_HARMONIC_SUB(ENERGY_SW,V_HARMONIC,P_HARMONIC)
	IMPLICIT NONE
	
	INCLUDE 'INITIAL_PARAMETERS.inc'	
	INCLUDE 'NOOFPARTICLES.inc'
	INCLUDE 'POSITIONS.inc'
	INCLUDE 'ENERGY_HARMONIC.inc'

	INTEGER ENERGY_SW
	DOUBLE PRECISION V_HARMONIC,P_HARMONIC

	INTEGER I,JNB,J,J1,J2
	DOUBLE PRECISION V,FSUM 
	DOUBLE PRECISION RXIJ,RYIJ,RIJSQ
	DOUBLE PRECISION D,DX,DY
	DOUBLE PRECISION VIJ,FIJ,FXIJ,FYIJ

	v=0.0d0
	fsum=0.0d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	IF(ENERGY_SW.eq.0)THEN
C-----------------------
       	  Do I=1,N	
	    Press_H(I)=0.0d0
	    POTEN_H(I)=0.0D0
	  EndDo
c	  LOOP OVER ALL MOLECULES
	  DO 111 Jnb=1,nbonds
	    I=Part1(Jnb)
	    J=Part2(Jnb)
	    rxij=x(i)-x(j)
	    ryij=y(i)-y(j)
C	    PERIODIC BOUNDARY CONDITION
	    rxij=rxij-boxx*anint(rxij/boxx)
            ryij=ryij-boxy*anint(ryij/boxy)	

	    rijsq=rxij*rxij+ryij*ryij

	    D = inX(I)-inX(J)
	    DX = D - BOXX*ANINT(D/BOXX)
	    D = inY(I)-inY(J)
	    DY = D - BOXY*ANINT(D/BOXY)

	    vij=.50d0*Econst*(abs(dsqrt(rijsq))-lp)**2
	    v=v+vij

	    POTEN_H(I)=POTEN_H(I)+VIJ
	    POTEN_H(J)=POTEN_H(J)+VIJ

	    fij=-Econst*(abs(dsqrt(rijsq))-lp)
	    fxij=fij*(rxij)
            fyij=fij*(ryij)	
	    fxij=fxij/abs(dsqrt(rijsq))
	    fyij=fyij/abs(dsqrt(rijsq))

	    fsum=fsum+fxij*rxij+fyij*ryij

	    Press_H(i)=Press_H(i)+fxij*rxij+fyij*ryij
            Press_H(j)=Press_H(j)+fxij*rxij+fyij*ryij
111	  CONTINUE
	  P_HARMONIC=0.50d0*fsum/(boxx*boxy)
	  V_HARMONIC=V
c-------------------------
	ELSE
C-------------------------
c	  LOOP OVER 7 PARTICLE
	  DO J2=1,7
	    if(j2.eq.7)then
	    i=ENERGY_SW	
	    else
	    i=nebors(ENERGY_SW,j2)
	    endif
	    Press_H(i)=0.0d0
	    POTEN_H(I)=0.0D0
	  ENDDO

	  DO 113 Jnb=1,7
	    if(Jnb.eq.7)then
	    I=ENERGY_SW
	    else
	    I=nebors(ENERGY_SW,Jnb)
	    endif
	    Do J1=1,6
	      j=nebors(i,j1)
	      rxij=x(i)-x(j)
	      ryij=y(i)-y(j)
C	      PERIODIC BOUNDARY CONDITION
	      rxij=rxij-boxx*anint(rxij/boxx)
              ryij=ryij-boxy*anint(ryij/boxy)	

	      rijsq=rxij*rxij+ryij*ryij

	      D = inX(I)-inX(J)
	      DX = D - BOXX*ANINT(D/BOXX)
	      D = inY(I)-inY(J)
	      DY = D - BOXY*ANINT(D/BOXY)

	      vij=.50d0*Econst*(abs(dsqrt(rijsq))-lp)**2

	      POTEN_H(I)=POTEN_H(I)+VIJ
	      
	      fij=-Econst*(abs(dsqrt(rijsq))-lp)
	      fxij=fij*(rxij)
	      fyij=fij*(ryij)	
	      fxij=fxij/abs(dsqrt(rijsq))
	      fyij=fyij/abs(dsqrt(rijsq))
	      Press_H(i)=Press_H(i)+fxij*rxij+fyij*ryij
	    EndDo
113	  CONTINUE

	  Do I=1,N
	    fsum=fsum+Press_H(i)
	    V=V+POTEN_H(I)
	  EndDo
	  P_HARMONIC=0.250d0*fsum/(boxx*boxy)
	  V_HARMONIC=0.50D0*V
C--------------------------------------------
	ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	RETURN
	END
