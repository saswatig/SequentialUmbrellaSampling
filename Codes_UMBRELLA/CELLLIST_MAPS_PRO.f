  	SUBROUTINE CELLLIST_MAPS_SUB

        implicit none
    
	include 'NOOFPARTICLES.inc'    
	include 'CELLLIST.inc'

        integer ix,iy,imap,icell

        icell(ix,iy)=1+mod(ix-1+m,m)+mod(iy-1+m,m)*m

        do iy=1,m
                do ix=1,m
                
                imap=(icell(ix,iy)-1)*4

                map(imap+1)=icell(ix+1,iy  )
                map(imap+2)=icell(ix+1,iy+1)
                map(imap+3)=icell(ix  ,iy+1)
                map(imap+4)=icell(ix-1,iy+1)
c        write(*,*)icell(ix,iy),imap,map(imap+1),map(imap+2),
c     & map(imap+3),map(imap+4)
                enddo
        enddo
        
        return
        end
