	SUBROUTINE MC_VOL_SUB(ICYCLE,EN,PRESSURE,V_HCHI,V_WCA)
	IMPLICIT NONE

	INCLUDE 'INITIAL_PARAMETERS.inc'
	INCLUDE 'NOOFPARTICLES.inc'
	INCLUDE 'POSITIONS.inc'
	INCLUDE 'UMBRELLA_SAMPLING.inc'
	INCLUDE 'ENERGY_HARMONIC.inc'
	INCLUDE 'ENERGY_WCA.inc'
	INCLUDE 'MATRIX_Y.inc'
	INCLUDE 'MATRIX_X.inc'
	INCLUDE 'ENERGY_HCHI.inc'

	INTEGER ICYCLE

	DOUBLE PRECISION EN0,ENN,EN
	DOUBLE PRECISION PRESSURE0,PRESSUREN,PRESSURE
	DOUBLE PRECISION V_HCHI0,V_HCHIN,V_HCHI
	DOUBLE PRECISION V_WCA0,V_WCAN,V_WCA
	DOUBLE PRECISION X0(N),Y0(N)
	
	INTEGER IDUM
	DOUBLE PRECISION RAN5

	DOUBLE PRECISION BOXX0,BOXY0,VOL0,VOLN,VOL
	DOUBLE PRECISION TEMP1,TEMP2,ARG
	
	INTEGER I	

	DOUBLE PRECISION CHI_AV

	IDUM=875654
		
	CALL ENERGY_SUB(ICYCLE,0,EN0,PRESSURE0,V_HCHI0,V_WCA0)

CCCCCCC CHANGE VOLUME CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	BOXX0=BOXX
	BOXY0=BOXY
	VOL0=BOXX0*BOXY0
	TEMP1=RAN5(IDUM)-0.50D0
	TEMP2=RAN5(IDUM)-0.50D0
	BOXX=BOXX+TEMP1*VXMAX*BOXX
	BOXY=BOXY+TEMP2*VYMAX*BOXY
	VOLN=BOXX*BOXY

ccccccc TRIAL MOVE cccccccccccccccccccccccc
	DO I=1,N
	X0(I)=X(I)
	Y0(I)=Y(I)
C	INX0(I)=INX(I)
C	INY0(I)=INY(I)
	X(I)=X(I)*BOXX/BOXX0
	Y(I)=Y(I)*BOXY/BOXY0
C	INX(I)=INX(I)*BOXX/BOXX0
C	INY(I)=INY(I)*BOXY/BOXY0
	
        if(x(I).gt.+0.50d0*boxx)then
        x(I)=x(I)-boxx
        end if
        if(x(I).lt.-0.50d0*boxx)then
        x(I)=x(I)+boxx
        end if

        if(y(I).gt.+.50d0*boxy)then
        y(I)=y(I) - boxy
        endif
        if(y(I).lt.-0.50d0*boxy)then
        y(I)=y(I)+boxy
        endif

	ENDDO
ccccccc TRIAL MOVE ccccccccccccccccccccccccc

	CALL ENERGY_SUB(ICYCLE,0,ENN,PRESSUREN,V_HCHIN,V_WCAN)

CCCCCCC ACCEPTENCE RULE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	arg=-beta*((enn-en0)+constP*(voln-vol0)
     &	-(n+1)*log(voln/vol0)/beta)

C============================================================================C
C============UMBRELLA SAMPLING===============================================C
C============================================================================C
	
	CHI_AV=V_HCHIN/DFLOAT(N)

	IF(RAN5(IDUM).LT.EXP(arg).AND.CHI_AV.GE.LBND
     &.AND.CHI_AV.LE.UBND)THEN
	  VOL=VOLN
	  EN=ENN
	  V_HCHI=V_HCHIN
	  V_WCA=V_WCAN
	  PRESSURE=PRESSUREN
	ELSE
	  VOL=VOL0
	  EN=EN0
	  V_HCHI=V_HCHI0
	  V_WCA=V_WCA0
	  PRESSURE=PRESSURE0
	  Do I=1,N
	  X(I)=X0(I)
	  Y(I)=Y0(I)	
C	  INX(I)=INX0(I)
C	  INY(I)=INY0(I)
	  EndDO
	  BOXX=BOXX0
	  BOXY=BOXY0

	  CALL ENERGY_SUB(ICYCLE,0,EN,PRESSURE,V_HCHI,V_WCA)
	ENDIF
C============================================================================C
C============UMBRELLA SAMPLING===============================================C
C============================================================================C

C	WRITE(24,*)ICYCLE,EN,boxx*boxy,dfloat(N)/(BOXX*BOXY),
C     &	PRESSURE

	RETURN
	END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
ccccccccccccccccccccccccccccccc
	FUNCTION ran5(idum)
	IMPLICIT NONE

              INTEGER idum
              INTEGER MBIG,MSEED,MZ

C       REAL MBIG,MSEED,MZ

        double precision ran5,FAC
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
        ran5=mj*FAC

        return
         END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccc
	
	
	

