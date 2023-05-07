	SUBROUTINE ENERGY_WCA_SUB(ENERGY_SW,V_WCA,P_WCA)
	IMPLICIT NONE

	INCLUDE 'NOOFPARTICLES.inc'
	INCLUDE 'POSITIONS.inc'
	INCLUDE 'VERLETLIST.inc'
	INCLUDE 'ENERGY_WCA.inc'


	INTEGER ENERGY_SW
	DOUBLE PRECISION V_WCA,P_WCA

	DOUBLE PRECISION SIGSQ,RCUTSQ,RCUT6,ECUT
	DOUBLE PRECISION V,FSUM
	INTEGER I,JJ,J,J2
	DOUBLE PRECISION RXI,RYI,RXIJ,RYIJ,RIJSQ
	DOUBLE PRECISION SR2,SR6,VIJ,WIJ
	DOUBLE PRECISION FIJ,FXIJ,FYIJ

	DOUBLE PRECISION VER_R_Check

	sigsq=sigma**2
        rcutsq=rcut**2
	rcut6=(dble(1/rcut))**6
	ecut=1.0d0

	 v=0.0d0
	 fsum=0.0d0

	IF(ENERGY_SW.eq.0)THEN
	  CALL CELLLIST_MAPS_SUB
          CALL CELLLIST_LINKS_SUB
          CALL VERLETLIST_SUB
C---------------------------------
C	LOOP OVER ALL PARTICLES 
          DO I=1,N	
	    Press_WCA(I)=0.0d0	
	    POTEN_WCA(I)=0.0D0
              Do jj=1,NVLIST(I)	
                rxi=x(i)
	        ryi=y(i)
                j=Vlist(I,JJ)
                rxij=rxi-x(j)
	        ryij=ryi-y(j)
C	        PERIODIC BOUNDARY CONDITION
	        rxij=rxij-boxx*nint(rxij/boxx)
                ryij=ryij-boxy*nint(ryij/boxy)	

	        rijsq=rxij*rxij+ryij*ryij

	        if(rijsq.lt.rcutsq)then
	          sr2=sigsq/rijsq
		  sr6=sr2**3
		  vij=4.0d0*sr6*(sr6-1.0d0)+ecut
		  wij=sr6*(sr6-.50d0)
		  v=v+vij
		  POTEN_WCA(I)=POTEN_WCA(I)+VIJ

                  fij=wij/rijsq
                  fxij=fij*rxij
                  fyij=fij*ryij

                  fsum=fsum+fxij*rxij+fyij*ryij
	          Press_WCA(I)=Press_WCA(I)+fxij*rxij+fyij*ryij
	        endif
	      EndDo
	    ENDDO
C	  INCORPORATE ENERGY FACTORS
          P_WCA=0.50d0*24.0d0*fsum/(boxx*boxy)
	  V_WCA=0.50d0*v
c---------------------
	ELSE
c---------------------
C	  CHECK IF NEW VERLET LIST REQUIRED

          VER_R_Check=dsqrt((X(ENERGY_SW)-XVL(ENERGY_SW))**2
     &			   +(Y(ENERGY_SW)-YVL(ENERGY_SW))**2)
          VER_R_Check=abs(VER_R_Check)
          IF(VER_R_Check.gt.0.50d0*(Rcut_Verlet-Rcut))Then
          CALL CELLLIST_LINKS_SUB
          CALL VERLETLIST_SUB
          ENDIF

c	  LOOP OVER 7 PARTICLE
	  Do J2=1,NVLIST(ENERGY_SW)+1
	    if(j2.eq.NVLIST(ENERGY_SW)+1)then
	    i=ENERGY_SW
	    else
	    i=VLIST(ENERGY_SW,j2)
	    endif
	    Press_WCA(i)=0.0d0
	    POTEN_WCA(I)=0.0D0
	  EndDo

	  Do J2=1,NVLIST(ENERGY_SW)+1
	    if(j2.eq.NVLIST(ENERGY_SW)+1)then
	    i=ENERGY_SW
	    else
	    i=VLIST(ENERGY_SW,j2)
	    endif

            Do jj=1,NVLIST(I)
	      rxi=x(i)
	      ryi=y(i)

	      j=Vlist(I,JJ)

	      rxij=rxi-x(j)
	      ryij=ryi-y(j)

C	      PERIODIC BOUNDARY CONDITION
	      rxij=rxij-boxx*nint(rxij/boxx)
              ryij=ryij-boxy*nint(ryij/boxy)	

	      rijsq=rxij*rxij+ryij*ryij

	      if(rijsq.lt.rcutsq)then
		sr2=sigsq/rijsq
		sr6=sr2**3
		vij=4.0d0*sr6*(sr6-1.0d0)+ecut
		wij=sr6*(sr6-.50d0)
		 
		POTEN_WCA(I)=POTEN_WCA(I)+VIJ

                fij=wij/rijsq
                fxij=fij*rxij
                fyij=fij*ryij
	        Press_WCA(I)=Press_WCA(I)+fxij*rxij+fyij*ryij
	      endif
	    EndDo
	  EndDo

	    do I=1,N
	      fsum=fsum+Press_WCA(I)
              v=v+POTEN_WCA(I)
	    enddo

C	  INCORPORATE ENERGY FACTORS
          P_WCA=0.50d0*24.0d0*fsum/(boxx*boxy)
	  V_WCA=0.5d0*v
C------------------------------
	ENDIF

	RETURN
	END
