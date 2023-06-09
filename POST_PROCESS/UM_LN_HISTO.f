	IMPLICIT NONE

	INTEGER I,J,LAST_UM_BINNO,SBIN,UBINS
	PARAMETER (UBINS=1000)
	DOUBLE PRECISION LBND(UBINS),UBND(UBINS)
	DOUBLE PRECISION UHISTO(UBINS,2),LNHISTO(UBINS)
        CHARACTER*123 F1,F2

        OPEN(11,FILE='INPUT')
        READ(11,'(A)')F2
        F1=F2
        CLOSE(11)

	I=1
	OPEN(26,FILE=F1)
88	READ(26,'(1(1x,I15),4(1x,g16.9))',END=87)I,LBND(I),
     &	UBND(I),UHISTO(I,1),UHISTO(I,2)
	I=I+1
	GOTO 88
87	CLOSE(26)

	LAST_UM_BINNO=I-1	

	DO I=1,LAST_UM_BINNO
	LNHISTO(I)=0.0D0
	ENDDO
	SBIN=1
	DO I=1,LAST_UM_BINNO
	  if(UHISTO(I,1).EQ.0.0d0)then
            if(UHISTO(I+1,1).EQ.0.0d0)then
            SBIN=I
            else
            SBIN=I
            GOTO 988    
            endif

	  endif
	ENDDO
	
988	DO I=SBIN,LAST_UM_BINNO
	  Do J=SBIN,I
	  if(UHISTO(J,1).GT.0.0d0)then
	  LNHISTO(I)=LNHISTO(I)+LOG(UHISTO(J,2)/UHISTO(J,1))	
	  endif	
	  EndDo
	ENDDO

        OPEN(11,FILE='OUTPUT')
        READ(11,'(A)')F2
        F1=F2
        CLOSE(11)

        OPEN(28,FILE=F1)
	DO I=1,LAST_UM_BINNO
	write(28,'(1x,I15,2(1x,g16.9))')I,
     &	(LBND(I)+UBND(I))*0.50D0,LNHISTO(I)
	ENDDO
        CLOSE(28)
	STOP
	END

	
