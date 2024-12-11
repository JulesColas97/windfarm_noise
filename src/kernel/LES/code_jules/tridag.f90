subroutine tridag
use sim_param, only : a,b,c,p=>tzz,RHS_col=>dummy4
use param,only : rp,lh,ld,ny,nz,coord,nproc,mpi_rp,up,down,localComm,status,ierr
use param,only : kmin2,kmax2 ! Defined in grid.f90
implicit none
!  RHS and p are interleaved as complex arrays
real(rp),dimension(lh,ny,nz+1)::gam
real(rp)::bet(lh,ny)
integer :: jx,jy,j,tag0,q,ii,ir
integer :: nchunks,cstart,cend
integer,parameter::chunksize=4

nchunks = ny/chunksize

if(coord==0) then
 do jy= 1,ny
  do jx= 1,lh-1
   ii=2*jx
   ir=ii-1
   bet(jx,jy)   =b(jx,jy,1)
   p(ir:ii,jy,0)=RHS_col(ir:ii,jy,0)/bet(jx,jy)
  end do
 end do
end if

do q=1,nchunks
 cstart=1+(q-1)*chunksize
 cend  =cstart+chunksize-1
 tag0=10*(q-1)
 if (coord/= 0) then
  call mpi_recv(bet(1,cstart)  ,lh*chunksize,mpi_rp,down,tag0+2,localComm,status,ierr)
  call mpi_recv(p  (1,cstart,0),ld*chunksize,mpi_rp,down,tag0+3,localComm,status,ierr)
 end if

 do j=2,kmax2
  do jy=cstart,cend
   do jx=1,lh-1
    gam(jx,jy,j)   =c(j-1)/bet(jx,jy)
    bet(jx,jy)     =b(jx,jy,j)-a(j)*gam(jx,jy,j)

    ii=2*jx
    ir=ii-1
!    p(ir:ii,jy,j-1)=(RHS_col(ir:ii,jy,j)-a(j)*p(ir:ii,jy,j-2))/bet(jx,jy)
    p(ir:ii,jy,j-1)=(RHS_col(ir:ii,jy,j-1)-a(j)*p(ir:ii,jy,j-2))/bet(jx,jy)
   enddo
  enddo 
 enddo

 if (coord/=nproc-1) then
  call mpi_send (bet(1,cstart)     ,lh*chunksize,mpi_rp,up,tag0+2,localComm,ierr)
  call mpi_send (p  (1,cstart,nz-1),ld*chunksize,mpi_rp,up,tag0+3,localComm,ierr)
 end if
end do

do q=1,nchunks!ny
 cstart=1+(q-1)*chunksize
 cend  =cstart+chunksize-1
 tag0=10*(q-1)
 if (coord/=nproc-1) then  
  call mpi_recv (p(1  ,cstart,nz)  ,ld*chunksize,mpi_rp,up,tag0+4,localComm,status,ierr)
  call mpi_recv (gam(1,cstart,nz+1),lh*chunksize,mpi_rp,up,tag0+5,localComm,status,ierr)
 end if

 do j=nz,kmin2,-1
  do jy=cstart,cend
   do jx=1,lh-1
    ii=2*jx
    ir=ii-1
    p(ir:ii,jy,j-1)=p(ir:ii,jy,j-1)-gam(jx,jy,j+1)*p(ir:ii,jy,j)
   enddo
  end do
 end do

! if (coord/= 0) then
 call mpi_send (p  (1,cstart,1),ld*chunksize,mpi_rp,down,tag0+4,localComm,ierr)
 call mpi_send (gam(1,cstart,2),lh*chunksize,mpi_rp,down,tag0+5,localComm,ierr)
! endif
end do

end subroutine tridag
