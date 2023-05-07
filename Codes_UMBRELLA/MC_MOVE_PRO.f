	SUBROUTINE MC_MOVE_SUB(ICYCLE,EN,PRESSURE,V_HCHI,V_WCA)
	IMPLICIT NONE

	INCLUDE 'INITIAL_PARAMETERS.inc'
	INCLUDE 'NOOFPARTICLES.inc'
	INCLUDE 'POSITIONS.inc'
	INCLUDE 'UMBRELLA_SAMPLING.inc'
	INCLUDE 'ENERGY_HARMONIC.inc'
	INCLUDE 'ENERGY_WCA.inc'
	INCLUDE 'MATRIX_X.inc'
	INCLUDE 'ENERGY_HCHI.inc'

	INTEGER ICYCLE

	INTEGER IDUM,O
	DOUBLE PRECISION RAN3

	DOUBLE PRECISION EN0,ENN,EN
	DOUBLE PRECISION PRESSURE0,PRESSUREN,PRESSURE
	DOUBLE PRECISION V_HCHI0,V_HCHIN,V_HCHI
	DOUBLE PRECISION V_WCA0,V_WCAN,V_WCA

	DOUBLE PRECISION X0,Y0

C	TO RESTORE THE ENERGY AND PRESSURE OF EACH PARTICLE FOR A REJECTED MC MOVE
	DOUBLE PRECISION ENPART_H0(7),PRESSPART_H0(7)
	DOUBLE PRECISION ENPART_WCA0(N),PRESSPART_WCA0(N)
	DOUBLE PRECISION ENPART_HCHI0(7),PRESSPART_HCHI0(7)
	DOUBLE PRECISION EPPART_HCHI0(7,2,2),MATXPART_HCHI0(7,2,2)

	INTEGER I,J,i1,i2

	DOUBLE PRECISION CHI_AV


	IDUM=98776
	O=INT(RAN3(IDUM)*N)+1

C====================================================================================C
C==================== UMBRELLA SAMPLING =============================================C
C====================================================================================C
	IF(ICYCLE.EQ.1)THEN
	CHECK_UM_CONF=1
	ENDIF

	IF(MOD((ICYCLE-1),RUNINBIN).eq.0)THEN
	  IF(CHECK_UM_CONF.EQ.0)THEN
C==========COLLECT DATA===========================C
	  BINNO=INT((ICYCLE-1)/RUNINBIN)+1
	  Do I=1,BINNO-1
	    LBND=LOWLIMUM+(I-1)*BINSZ
	    UBND=LOWLIMUM+(I+1)*BINSZ
	    write(24,'(1(1x,I15),4(1x,g16.9))')I,LBND,UBND,
     &		                    UHISTO(I,1),UHISTO(I,2)
	  EndDo
C==========COLLECT DATA===========================C
	  STOP 'UMBRELLA SAMPLING NEW WINDOW DID NOT GET NEW CONFIG'
	  ELSE
	  CHECK_UM_CONF=0
	  ENDIF

	  DO I=1,N
	  X(I)=XUM(I)
	  Y(I)=YUM(I)
	  ENDDO
	  BOXX=BOXXUM
	  BOXY=BOXYUM
	
	  LBND=LOWLIMUM+INT((ICYCLE-1)/RUNINBIN)*BINSZ
	  UBND=LOWLIMUM+(INT((ICYCLE-1)/RUNINBIN)+2)*BINSZ
          BINNO=INT((ICYCLE-1)/RUNINBIN)+1

	  CALL ENERGY_SUB(ICYCLE,0,EN0,PRESSURE0,V_HCHI0,V_WCA0)
	ELSE
	  CALL ENERGY_SUB(ICYCLE,O,EN0,PRESSURE0,V_HCHI0,V_WCA0)
	ENDIF

C====================================================================================C
C==================== UMBRELLA SAMPLING =============================================C
C====================================================================================C


C	TO RESTORE THE ENERGY AND PRESSURE OF EACH PARTICLE FOR A REJECTED MC MOVE
	DO J=1,7
	  if(j.eq.7)then
	  i=O
	  else
	  i=nebors(O,j)
	  endif
	  ENPART_H0(J)=POTEN_H(I)
	  ENPART_HCHI0(J)=CHI(I)

	  PRESSPART_H0(J)=PRESS_H(I)
	  PRESSPART_HCHI0(J)=PRESS_HCHI(I)

	  do i1=1,2
	  do i2=1,2
	  EPPART_HCHI0(J,i1,i2)=EP(I,i1,i2)
	  MATXPART_HCHI0(J,i1,i2)=MATX(I,i1,i2)
	  enddo
	  enddo
	ENDDO

	DO J=1,N
	  ENPART_WCA0(J)=POTEN_WCA(J)
	  PRESSPART_WCA0(J)=PRESS_WCA(J)
	ENDDO
C=====================================================================
		
	X0=X(O)
	Y0=Y(O)
	X(O)=X(O)+(RAN3(IDUM)-0.50D0)*MDELX*LP
	Y(O)=Y(O)+(RAN3(IDUM)-0.50D0)*MDELY*LP

        if(x(O).gt.+0.50d0*boxx)then
        x(O)=x(O)-boxx
        end if
        if(x(O).lt.-0.50d0*boxx)then
        x(O)=x(O)+boxx
        end if

        if(y(O).gt.+.50d0*boxy)then
        y(O)=y(O) - boxy
        endif
        if(y(O).lt.-0.50d0*boxy)then
        y(O)=y(O)+boxy
        endif
C========================================================================================C
C========= UMBRELLA SAMPLING ============================================================C
C========================================================================================C

	CALL ENERGY_SUB(ICYCLE,O,ENN,PRESSUREN,V_HCHIN,V_WCAN)

	CHI_AV=V_HCHIN/DFLOAT(N)

	IF(ICYCLE.EQ.1)THEN
	DO I1=1,UBINS
	  Do I2=1,2
	  UHISTO(I1,I2)=0.0D0
	  EndDo
	ENDDO
	ENDIF
	IF(CHI_AV.LT.LBND)THEN
	UHISTO(BINNO,1)=UHISTO(BINNO,1)+1.0D0
	ENDIF

	IF(CHI_AV.GT.UBND)THEN
	UHISTO(BINNO,2)=UHISTO(BINNO,2)+1.0D0
	ENDIF

	IF(CHI_AV.GE.LBND.AND.CHI_AV.LE.UBND)THEN
               If(CHI_AV.GT.(LBND+BINSZ).AND.CHI_AV.LE.UBND)Then
		  do I=1,N
		    XUM(I)=X(I)
	            YUM(I)=Y(I)
c		    write(86,*)i,xum(i),yum(i),lbnd,ubnd,chi_av,binno
		  enddo
c		    write(86,*)
c		    write(87,*)icycle,lbnd,ubnd,chi_av,binno
	 	    BOXXUM=BOXX
		    BOXYUM=BOXY
		    CHECK_UM_CONF=1
	       EndIf
CCCCCCC   ACCEPTENCE RULE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	    IF(RAN3(IDUM).LT.EXP(-BETA*(ENN-EN0)))THEN
	       EN=ENN
	       PRESSURE=PRESSUREN
	       V_HCHI=V_HCHIN
	       V_WCA=V_WCAN
               If(CHI_AV.GT.(LBND+BINSZ).AND.CHI_AV.LE.UBND)Then
	          UHISTO(BINNO,2)=UHISTO(BINNO,2)+1.0D0
	       Else
	          UHISTO(BINNO,1)=UHISTO(BINNO,1)+1.0D0
	       EndIf
	    ELSE
	       EN=EN0
	       PRESSURE=PRESSURE0
	       V_HCHI=V_HCHI0
	       V_WCA=V_WCA0
	       X(O)=X0
	       Y(O)=Y0
C	       TO RESTORE THE ENERGY AND PRESSURE OF EACH PARTICLE FOR A REJECTED MC MOVE
	       DO J=1,7
	         if(j.eq.7)then
	           i=O
	         else
	           i=nebors(O,j)
	         endif
	         POTEN_H(I)=ENPART_H0(J)
	         CHI(I)=ENPART_HCHI0(J)

	         PRESS_H(I)=PRESSPART_H0(J)
	         PRESS_HCHI(I)=PRESSPART_HCHI0(J)

	         do i1=1,2
	         do i2=1,2
	           EP(I,i1,i2)=EPPART_HCHI0(J,i1,i2)
	           MATX(I,i1,i2)=MATXPART_HCHI0(J,i1,i2)
	         enddo
	         enddo
	       ENDDO
	       
	       DO J=1,N	
		POTEN_WCA(J)=ENPART_WCA0(J)
	        PRESS_WCA(J)=PRESSPART_WCA0(J)
	       ENDDO
	     ENDIF
CCCCCCC   ACCEPTENCE RULE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	ELSE
	       EN=EN0
	       PRESSURE=PRESSURE0
	       V_HCHI=V_HCHI0
	       V_WCA=V_WCA0
	       X(O)=X0
	       Y(O)=Y0
C	       TO RESTORE THE ENERGY AND PRESSURE OF EACH PARTICLE FOR A REJECTED MC MOVE
	       DO J=1,7
	         if(j.eq.7)then
	           i=O
	         else
	           i=nebors(O,j)
	         endif
	         POTEN_H(I)=ENPART_H0(J)
	         CHI(I)=ENPART_HCHI0(J)

	         PRESS_H(I)=PRESSPART_H0(J)
	         PRESS_HCHI(I)=PRESSPART_HCHI0(J)

	         do i1=1,2
	         do i2=1,2
	           EP(I,i1,i2)=EPPART_HCHI0(J,i1,i2)
	           MATX(I,i1,i2)=MATXPART_HCHI0(J,i1,i2)
	         enddo
	         enddo
	       ENDDO

	       DO J=1,N	
		POTEN_WCA(J)=ENPART_WCA0(J)
	        PRESS_WCA(J)=PRESSPART_WCA0(J)
	       ENDDO
	ENDIF
C=========================================================================================C
C========= UMBRELLA SAMPLING =============================================================C
C=========================================================================================C

C==========COLLECT DATA===========================C
	IF (ICYCLE.EQ.NCYCLE)THEN
	  Do I=1,UBINS
	    LBND=LOWLIMUM+(I-1)*BINSZ
	    UBND=LOWLIMUM+(I+1)*BINSZ
	    write(24,'(1(1x,I15),4(1x,g16.9))')I,LBND,UBND,
     &		                    UHISTO(I,1),UHISTO(I,2)
	  EndDo
	ENDIF
C==========COLLECT DATA===========================C


	RETURN
	END

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C----------RANDOM NUMBER GENERTOR----------------

ccccccccccccccccccccccccccccccc
	FUNCTION ran3(idum)
	IMPLICIT NONE

              INTEGER idum
              INTEGER MBIG,MSEED,MZ

C       REAL MBIG,MSEED,MZ

        double precision ran3,FAC
        PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)

C       PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)

        INTEGER i,iff,ii,inext,inextp,k
        INTEGER mj,mk,ma(55)

C       REAL mj,mk,ma(55)

        SAVE iff,inext,inextp,ma

        DATA iff /0/

        if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
        endif
        inext=inext+1
        if(inext.eq.56)inext=1
        inextp=inextp+1
        if(inextp.eq.56)inextp=1
        mj=ma(inext)-ma(inextp)
        if(mj.lt.MZ)mj=mj+MBIG
        ma(inext)=mj
        ran3=mj*FAC

        return
         END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccc
	
	 


