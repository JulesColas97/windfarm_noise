!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates and writes out time averaged velocity histograms
! Contains the following subroutines: 
!
! - pdf             : Driver routine, determines when to calculate and write out pdf
! - pdf_init        : Initilizes arrays used in computation of pdf
! - checkpoint_pdf  : Writes PDF to file
! - tavg_streamw_pdf: Collects velocity values per bin and saves it as histogram
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined(TAVG_PDF) 
subroutine pdf
use param, only : jt_total,total_time,dt,rp,coord
use param, only : checkpoint_nskip,tavg_nskip,tavg_nstart,tavg_start_time,tavg_end_time
use param, only : pdf_dt,pdf_total_time
use io, only : pdf_initialized
implicit none

! Store data both on numered as well as normal histogram files
if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint_pdf(2) !Numbered files
if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint_pdf(1) !Regular files

!  Determine if time summations are to be calculated
if (total_time >= tavg_start_time .and. total_time <= tavg_end_time) then
 pdf_dt = pdf_dt + dt

 ! Initialization:
 if(total_time >= tavg_start_time .and. .not.pdf_initialized) then
  tavg_nstart = jt_total
  if (coord == 0) then
   write(*,*) '-------------------------------'
   write(*,"(1a,i9,1a,i9)") 'Starting PDF calculations from ', tavg_nstart
   write(*,*) '-------------------------------'
  endif  ! coord==0
  call pdf_init
 endif ! end initialization

 if(mod(jt_total-tavg_nstart,tavg_nskip).eq.0.and.jt_total-tavg_nstart.gt.0) then
  call tavg_streamw_pdf
  pdf_total_time= pdf_total_time+pdf_dt ! Update total averaging time
  pdf_dt=0.0_rp ! To determine time interval with next call to "pdf"
 endif
endif 
end subroutine pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Checkpoint routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint_pdf (version)
use param, only : rp,path,nz,nbins
use io, only : pdf_initialized
use io, only : tavg_u_histogram,tavg_u_xy,tavg_u2_xy
use io, only : tavg_v_histogram,tavg_v_xy,tavg_v2_xy
use io, only : tavg_w_histogram,tavg_w_xy,tavg_w2_xy
implicit none
integer,intent(IN) :: version
character (64) :: folder

! Make backup of the histogram data files 
if( pdf_initialized ) then
 folder=trim(path)//'output'
 call write_PDF(version,tavg_u_histogram(1:nbins,1:nz-1),'tavg_u_histogram',folder,0,nz-1)
 call write_PDF(version,tavg_v_histogram(1:nbins,1:nz-1),'tavg_v_histogram',folder,0,nz-1)
 call write_PDF(version,tavg_w_histogram(1:nbins,1:nz-1),'tavg_w_histogram',folder,0,nz-1)
 call write_1Dfield(version,tavg_u_xy(1:nz-1),'tavg_u_xy',folder,0,nz-1)
 call write_1Dfield(version,tavg_u2_xy(1:nz-1),'tavg_u2_xy',folder,0,nz-1)
 call write_1Dfield(version,tavg_v_xy(1:nz-1),'tavg_v_xy',folder,0,nz-1)
 call write_1Dfield(version,tavg_v2_xy(1:nz-1),'tavg_v2_xy',folder,0,nz-1)
 call write_1Dfield(version,tavg_w_xy(1:nz-1),'tavg_w_xy',folder,0,nz-1)
 call write_1Dfield(version,tavg_w2_xy(1:nz-1),'tavg_w2_xy',folder,0,nz-1)
endif
end subroutine checkpoint_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the pdf variable to zero or read from file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pdf_init
use param, only : rp,nz,path,coord,mpi_rp,localComm,ierr
use param, only : pdf_total_time,pdf_dt
use param, only : nbins,bin_min_u,bin_max_u,bin_min_v,bin_max_v,bin_min_w,bin_max_w
use io, only : pdf_initialized
use io, only : binsize_u,startbin_u,binsize_v,startbin_v,binsize_w,startbin_w
use io, only : tavg_u_histogram,tavg_u_xy,tavg_u2_xy
use io, only : tavg_v_histogram,tavg_v_xy,tavg_v2_xy
use io, only : tavg_w_histogram,tavg_w_xy,tavg_w2_xy
implicit none
real(rp) :: test_sum_hist
character (64) :: folder
logical :: exst

allocate(tavg_u_histogram(nbins,1:nz-1))
allocate(tavg_u_xy(1:nz-1))
allocate(tavg_u2_xy(1:nz-1))
allocate(tavg_v_histogram(nbins,1:nz-1))
allocate(tavg_v_xy(1:nz-1))
allocate(tavg_v2_xy(1:nz-1))
allocate(tavg_w_histogram(nbins,1:nz-1))
allocate(tavg_w_xy(1:nz-1))
allocate(tavg_w2_xy(1:nz-1))

folder=trim(path)//'output'

!Initialize time averaging structures 
  pdf_dt           =0.0_rp
  pdf_total_time   =0.0_rp

  tavg_u_histogram(:,:) =0.0_rp
  tavg_u_xy(:)          =0.0_rp
  tavg_u2_xy(:)         =0.0_rp
  tavg_v_histogram(:,:) =0.0_rp
  tavg_v_xy(:)          =0.0_rp
  tavg_v2_xy(:)         =0.0_rp
  tavg_w_histogram(:,:) =0.0_rp
  tavg_w_xy(:)          =0.0_rp
  tavg_w2_xy(:)         =0.0_rp

  binsize_u = (bin_max_u - bin_min_u) / float(nbins)

  ! Create a "startbin" which will be added to ibinu 
  ! in subroutine tavg_pdf to make sure that ibinu never 
  ! gets <=0
  if(bin_min_u < 0.0_rp .and. bin_max_u > 0.0_rp) then
   startbin_u = -floor(bin_min_u/binsize_u) + 1
  else
   startbin_u = 1
  endif

  binsize_v = (bin_max_v - bin_min_v) / float(nbins)
  if(bin_min_v < 0.0_rp .and. bin_max_v > 0.0_rp) then
   startbin_v = -floor(bin_min_v/binsize_v) + 1
  else
   startbin_v = 1
  endif

  binsize_w = (bin_max_w - bin_min_w) / float(nbins)
  if(bin_min_w < 0.0_rp .and. bin_max_w > 0.0_rp) then
   startbin_w = -floor(bin_min_w/binsize_w) + 1
  else
   startbin_w = 1
  endif

inquire (file=trim(folder)//'/tavg_u_histogram.h5',exist=exst)
if (.not. exst) then ! Initialize when no previous data are present
  if (coord == 0) then
    write(*,*)'No previous PDF data - starting from scratch.'
  endif
else ! Read data from file
  if (coord == 0) then
   write(*,*) 'Read in PDF from previous run'
  endif
  !Send information from master to other tasks
  call mpi_bcast(pdf_total_time,1,mpi_rp,coord,localComm,ierr)

 call read_PDF(tavg_u_histogram(1:nbins,1:nz-1),'tavg_u_histogram',folder,0,nz-1)
 call read_PDF(tavg_v_histogram(1:nbins,1:nz-1),'tavg_v_histogram',folder,0,nz-1)
 call read_PDF(tavg_w_histogram(1:nbins,1:nz-1),'tavg_w_histogram',folder,0,nz-1)
 call read_1Dfield(tavg_u_xy(1:nz-1),'tavg_u_xy',folder,0,nz-1)
 call read_1Dfield(tavg_u2_xy(1:nz-1),'tavg_u2_xy',folder,0,nz-1)
 call read_1Dfield(tavg_v_xy(1:nz-1),'tavg_v_xy',folder,0,nz-1)
 call read_1Dfield(tavg_v2_xy(1:nz-1),'tavg_v2_xy',folder,0,nz-1)
 call read_1Dfield(tavg_w_xy(1:nz-1),'tavg_w_xy',folder,0,nz-1)
 call read_1Dfield(tavg_w2_xy(1:nz-1),'tavg_w2_xy',folder,0,nz-1)

  test_sum_hist = sum(tavg_u_histogram(1,1:nz-1)) + sum(tavg_u_histogram(nbins,1:nz-1)) &
              & + sum(tavg_v_histogram(1,1:nz-1)) + sum(tavg_v_histogram(nbins,1:nz-1)) &
              & + sum(tavg_w_histogram(1,1:nz-1)) + sum(tavg_w_histogram(nbins,1:nz-1))

  if (test_sum_hist > 0.01_rp) then
         stop "!! ADJUSTMENT of range of bin_min or bin_max needed. When done, & 
             &  DELETE old tavg_*_histogram.h5 files !!!"
  endif
endif
! Set global switch that tavg of pdf has been initialized:
pdf_initialized=.true.
end subroutine pdf_init

subroutine tavg_streamw_pdf
use sim_param,only : u,v,w
use param, only : nx,ny,nz,rp,inxny
use param, only : nbins,pdf_total_time,pdf_dt
use io, only: binsize_u,startbin_u
use io, only: binsize_v,startbin_v
use io, only: binsize_w,startbin_w
use io, only: tavg_u_histogram,tavg_u_xy,tavg_u2_xy
use io, only: tavg_v_histogram,tavg_v_xy,tavg_v2_xy
use io, only: tavg_w_histogram,tavg_w_xy,tavg_w2_xy
implicit none
real(rp),dimension(nbins,nz-1) :: temp_hist_u
real(rp),dimension(nbins,nz-1) :: temp_hist_v
real(rp),dimension(nbins,nz-1) :: temp_hist_w
real(rp) :: gg1,gg2
integer :: i,j,k,ibinu,ibinv,ibinw

! Used for time averaging:
gg1=pdf_total_time/(pdf_total_time+pdf_dt)
gg2=pdf_dt        /(pdf_total_time+pdf_dt)

! in tavg_u_histogram all u-values of one bin summed up,
! result is written to file in subroutine checkpoint using write_pdf.f90
! mean of u and u2 are calculated and saved for the normalization of the pdf 

do k=1,nz-1
 temp_hist_u(:,k) = 0._rp
 temp_hist_v(:,k) = 0._rp
 temp_hist_w(:,k) = 0._rp
 do j=1,ny
   do i = 1,nx
      ibinu = startbin_u + floor(u(i,j,k)/binsize_u)
      ibinv = startbin_v + floor(v(i,j,k)/binsize_v)
      ibinw = startbin_w + floor(w(i,j,k)/binsize_w)

      ! Insure that ibin does not extend range of histogram, in
      ! case it does, values are saved at end of histogram:
      ! (startbin insures that its never <=0, defined in tavg_init) 
      ibinu = min(ibinu,nbins)
      ibinv = min(ibinv,nbins)
      ibinw = min(ibinw,nbins)

      ibinu = max(ibinu,1)
      ibinv = max(ibinv,1)
      ibinw = max(ibinw,1)

      temp_hist_u(ibinu,k) = temp_hist_u(ibinu,k) + 1._rp
      temp_hist_v(ibinv,k) = temp_hist_v(ibinv,k) + 1._rp
      temp_hist_w(ibinw,k) = temp_hist_w(ibinw,k) + 1._rp
   enddo
 enddo
 tavg_u_histogram(:,k) = gg1 * tavg_u_histogram(:,k) + gg2 * temp_hist_u(:,k)
 tavg_v_histogram(:,k) = gg1 * tavg_v_histogram(:,k) + gg2 * temp_hist_v(:,k)
 tavg_w_histogram(:,k) = gg1 * tavg_w_histogram(:,k) + gg2 * temp_hist_w(:,k)

 tavg_u_xy(k)   = gg1 * tavg_u_xy(k)  + gg2 * inxny*sum(u(1:nx,1:ny,k))
 tavg_u2_xy(k)  = gg1 * tavg_u2_xy(k) + gg2 * inxny*sum(u(1:nx,1:ny,k)**2)
 tavg_v_xy(k)   = gg1 * tavg_v_xy(k)  + gg2 * inxny*sum(v(1:nx,1:ny,k))
 tavg_v2_xy(k)  = gg1 * tavg_v2_xy(k) + gg2 * inxny*sum(v(1:nx,1:ny,k)**2)
 tavg_w_xy(k)   = gg1 * tavg_w_xy(k)  + gg2 * inxny*sum(w(1:nx,1:ny,k))
 tavg_w2_xy(k)  = gg1 * tavg_w2_xy(k) + gg2 * inxny*sum(w(1:nx,1:ny,k)**2)
enddo

!tavg_u_histogram = tavg_u_histogram * inxny / binsize_u Can be done in
!postprocessing
!tavg_v_histogram = tavg_v_histogram * inxny / binsize_v
!tavg_w_histogram = tavg_w_histogram * inxny / binsize_w
end subroutine tavg_streamw_pdf
#endif
