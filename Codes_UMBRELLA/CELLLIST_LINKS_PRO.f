	SUBROUTINE CELLLIST_LINKS_SUB

        implicit none
	
	include 'NOOFPARTICLES.inc'	
	include 'POSITIONS.inc'
	include 'CELLLIST.inc'
	include 'VERLETLIST.inc'
        
        integer a,b
        double precision celli,cellx,celly
        integer icell,i

        do icell=1,ncell
        head(icell)=0
c	write(*,*)ncell
        end do
c        do i=1,n
c        write(90,*)x(i),y(i)
c        end do

        celli=dble(m)
        cellx=boxx/celli
	celly=boxy/celli
c	write(*,*)cell,box,rcut
        if (cellx.lt.rcut_verlet.or.celly.lt.rcut_verlet)then
        stop 'cell size too small for cutoff'
        end if

        do i=1,n
        
        a=int((x(i)+.50d0*boxx)/cellx)
        b=int((y(i)+.50d0*boxy)/celly)
        icell=1+a+b*m
        
        list(i)=head(icell)
        head(icell)=i
c        write(101,*)icell,head(icell)
        end do
c	write(101,*)
c	do i=1,ncell
c	write(102,*)i,head(i)
c	enddo
c	write(102,*)
c        do i=1,n
c        write(106,*)x(i),y(i)
c        end do
        return
        end

