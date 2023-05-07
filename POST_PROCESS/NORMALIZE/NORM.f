	INTEGER DIV,I1,HABI
	DOUBLE PRECISION CHI(1000),PHI(1000),NCONST
	DOUBLE PRECISION PCHI(1000),DIFF
        CHARACTER*123 F1,F2

        OPEN(11,FILE='INPUT')
        READ(11,'(A)')F2
        F1=F2
        CLOSE(11)

c        OPEN(12,FILE='h.dat')
c        READ(12,*)H
c        CLOSE(12)
	
	DIV=1
	OPEN(87,FILE=F1)
98	READ(87,'(1x,I15,2(1x,g16.9))',END=99)HABI,
     &	CHI(DIV),PHI(DIV)
	phi(div)=exp(phi(div))
	DIV=DIV+1
	GOTO 98
99	CLOSE(87)
	DIFF=CHI(2)-CHI(1)

	  NCONST=0.0D0	
	  DO I1=1,div-1
	  NCONST=NCONST+PHI(I1)*DIFF
	  ENDDO

        OPEN(11,FILE='OUTPUT')
        READ(11,'(A)')F2
        F1=F2
        CLOSE(11)

	open(25,file=F1)
	  DO I1=1,div-1
	  PCHI(I1)=PHI(I1)/(NCONST)
	  WRITE(25,'(3(1x,g16.9))')CHI(I1),PHi(I1),PCHI(I1)
	  ENDDO	
          write(25,*)
          write(25,*)
	close(25)
CCCCCCCCCCCCCCCCCCCC
	STOP
	END
