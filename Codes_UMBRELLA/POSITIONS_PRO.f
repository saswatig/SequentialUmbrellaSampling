	SUBROUTINE POSITIONS_SUB
	IMPLICIT NONE

	INCLUDE 'INITIAL_PARAMETERS.inc'
	INCLUDE 'NOOFPARTICLES.inc'
	INCLUDE 'POSITIONS.inc'	
	INCLUDE 'UMBRELLA_SAMPLING.inc'

	INTEGER NBCOUNT,IP,IB,J,I,NBCOUNTER(N),NPART,JP
	DOUBLE PRECISION DX,DY,DR,XC,YC

	LP=DSQRT((2.0D0)/(DSQRT(3.0D0)*RHO))
        boxx = dfloat(nu)*lp
        boxy = dfloat(nu)*0.50d0*dsqrt(3.0d0)*lp
	INBOXX=BOXX
	INBOXY=BOXY

	BOXXUM=BOXX
	BOXYUM=BOXY

	nbcount=0
        ip = 0
        ib = 0
        do j = 1,nu/2
           do i = 1,nu
              ip = ip+1
              x(ip)=dfloat(i-1)*lp
              y(ip)=dfloat(j-1)*dsqrt(3.0d0)*lp
              ip = ip+1
              x(ip)=x(ip-1)+0.50d0*lp
              y(ip)=y(ip-1)+dsqrt(3.0d0)*0.50d0*lp
           end do
        end do
        npart = ip

        write(20,*)npart,' atoms'
        write(20,*)3*npart,' bonds'
        write(20,*)
        write(20,*)' 1 atom types'
        write(20,*)' 1 bond types'
        write(20,*)
        write(20,*)'0.0 ',boxx,' xlo xhigh'
        write(20,*)'0.0 ',boxy,' ylo yhigh'
        write(20,*)'-0.001 0.001 zlo zhigh'
        write(20,*)
        write(20,*)'Atoms '
        write(20,*)
        
        do ip = 1,npart
           write(20,*)ip,x(ip),y(ip)
	nbcounter(ip)=0
        end do
        write(20,*)
        write(20,*)'Bonds '
        write(20,*)
        ib = 0 
        do ip = 1,npart-1
           do jp = ip+1,npart
              dx = x(jp) - x(ip)
              if(abs(dx).ge.boxx/2.0d0)then
                 if(dx.gt.0)then
                    dx=dx-boxx
                 else
                    dx=dx+boxx
                 endif
              endif
              dy = y(jp) - y(ip)
              if(abs(dy).ge.boxy/2.0d0)then
                 if(dy.gt.0)then
                    dy=dy-boxy
                 else
                    dy=dy+boxy
                 endif
              endif
              dr = dsqrt(dx*dx+dy*dy)
              if(dr.le.1.7*lp)then 
                 ib = ib+1
		nbcount=nbcount+1
		part1(nbcount)=ip
		part2(nbcount)=jp
		nbcounter(ip)=nbcounter(ip)+1
		nbcounter(jp)=nbcounter(jp)+1
		nebors(ip,nbcounter(ip))=jp
		nebors(jp,nbcounter(jp))=ip
                 write(20,*)ib,' 1 ',ip,jp
              endif
            end do
        end do
        write(20,*)
        write(20,*)'Bond Coeffs'
        write(20,*)
        write(20,*)'1 0.5', lp
        write(20,*)
        write(20,*)'Masses'
        write(20,*)
        write(20,*)'1 1.0'

	xc=0.0d0
	yc=0.0d0
	do ip = 1,n
	xc=xc+x(ip)
	yc=yc+y(ip)
        end do
	xc=xc/dfloat(n)
	yc=yc/dfloat(n)
	do ip=1,n
	x(ip)=x(ip)-xc
	y(ip)=y(ip)-yc
	inx(ip)=x(ip)
	iny(ip)=y(ip)
	XUM(IP)=X(IP)
	YUM(IP)=Y(IP)
	enddo

	RETURN
	END
