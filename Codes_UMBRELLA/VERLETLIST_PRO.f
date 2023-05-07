	SUBROUTINE VERLETLIST_SUB

	IMPLICIT NONE
	include 'NOOFPARTICLES.inc'
	include 'POSITIONS.inc'
	include 'CELLLIST.inc'
	include 'VERLETLIST.inc'

	double precision rxi,ryi,rcutsq
        double precision rijsq,rxij,ryij
        integer icell,jcello,jcell,nabor,i,j
	integer NVLISTI


        rcutsq=rcut_verlet**2

	DO I=1,N
	  NVLIST(I)=0
	  XVL(I)=X(I)
	  YVL(I)=Y(I)
	ENDDO

C	LOOP OVER ALL CELLS
	
	do 5000 icell=1,ncell
		
		i=head(icell)

c	LOOP OVER ALL MOLECULES IN THE CELL
1000	if (i.gt.0) then
	
          rxi=x(i)
	  ryi=y(i)

	  NVLISTI=NVLIST(i)
C	LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL
	j=list(i)

2000	if(j.gt.0)then

	  rxij=rxi-x(j)
	  ryij=ryi-y(j)

C	PERIODIC BOUNDARY CONDITION
	  rxij=rxij-boxx*nint(rxij/boxx)
          ryij=ryij-boxy*nint(ryij/boxy)	

	  rijsq=rxij*rxij+ryij*ryij

	  if(rijsq.lt.rcutsq)then
		NVLISTI=NVLISTI+1
		NVLIST(j)=NVLIST(j)+1
	        VLIST(i,NVLISTI)=j
	        VLIST(j,NVLIST(j))=i	        
	  endif

	  j=list(j)

	  goto 2000
	
	endif

C	LOOP OVER NEIGHBOURING CELLS
	jcello=4*(icell-1)

	do 4000 nabor=1,4

	  jcell=map(jcello+nabor)

C	LOOP OVER ALL MOLECULES IN NEIGHBOURING CELLS

	  j=head(jcell)

3000	  if(j.ne.0)then

		rxij=rxi-x(j)
		ryij=ryi-y(j)

		rxij=rxij-boxx*nint(rxij/boxx)
                ryij=ryij-boxy*nint(ryij/boxy)

		rijsq=rxij*rxij+ryij*ryij
		
		if(rijsq.lt.rcutsq)then
		  NVLISTI=NVLISTI+1
		  NVLIST(J)=NVLIST(J)+1
	          VLIST(i,NVLISTI)=j
	          VLIST(j,NVLIST(j))=i	  
	  	endif

		j=list(j)
		
		goto 3000
	  endif

4000	continue
	
	NVLIST(i)=NVLISTI
	i=list(i)

	goto 1000
	
	endif

5000	continue

	RETURN
	END	

	
