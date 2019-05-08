program whatever

integer :: i, p, n
integer :: dx
integer :: imax
integer :: tmax

real :: t, dt, yr2sec, maxVolumeLoss
real :: tout
real :: dtout
real :: tprogress
real :: dtprogress

real, dimension (:,:,:), allocatable :: Deposit_fractions

real*16, dimension (:,:), allocatable :: h, hn
real*16, dimension (:), allocatable :: x, bedrock, totalHeight
real*16, dimension (3) :: k, q0
real*16 :: diffMagnitude, diffGradient, subsidence_rate

dx= 1.e3       ! grid spacing (m)
imax= 101      ! number of nodes
dt= 1          ! time step in (yr)
tmax= 300000     ! max number of years
dtout = 10000     ! nr of years between write output
tout = 0          ! threshold to write next output in years (increased every write action by dtout)
tprogress = 0     !
dtprogress = tmax/1000 !
k(1)= 3.2e-3      ! diffusivity (m2/s)
k(2)= 3.2e-3      ! diffusivity (m2/s)
k(3)= 3.2e-4      ! diffusivity (m2/s)
q0(1)= 0!3.e-8      ! sediment input (m2/s)
q0(2)= 1.e-10      ! sediment input (m2/s)
q0(3)= 2.e-10      ! sediment input (m2/s)
subsidence_rate = 0.!4e-4*dt/3 !m/yr
maxVolumeLoss = 0

allocate ( Deposit_fractions(0:imax,0:500,1:3) )  !imax (x_nodes, y_nodes, nr_of_grain_sizes)
allocate ( x(0:imax) )
allocate ( bedrock(0:imax), totalHeight(0:imax) )
allocate ( h(0:imax,3), hn(0:imax,3) )

yr2sec = 60*60*24*365.25

! dt = 0.9*(dx*dx)/(2.e0*maxval(k))
! print*, 'max dt = ', 0.9*(dx*dx)/(2.e0*maxval(k))


do i = 0,imax      ! initialise elevation
  do p = 1,3
      h(i,p)= 0.
!       h(i,3) = 400*sin(3.1415*0.5*i/imax)
      bedrock(i) = 0.
      totalHeight(i) = 0.
  end do
  x(i)=i*dx/1000
end do

t = 0.
dt = 1
do while (t < tmax)      ! main loop (a loop over time)
    
    if (t.ge.tout) then
        n= nint(tout/dtout)
        call writeOutput(n, imax, x, hn, bedrock, totalHeight)
        tout = tout + dtout
    endif
    
    dt = min( dt, 0.7*(dx*dx)/(2.e0*(maxval(k)+maxVolumeLoss+subsidence_rate))/yr2sec )
    
    do p=1,3
        hn(0,p)= hn(0,p) + q0(p)*dt*365.25*24*3600 !1000!-bedrock(0)!h(0,p) + q0(p)*dt*365.25*24*3600  ! boundary condition at left end of basin
!         hn(0,2)= bedrock(0)
        hn(0,1)= bedrock(0)
        hn(imax,p) = 0. ! boundary condition at right end of basin
    end do
    
    call diffusion(h, hn, bedrock, totalHeight, imax, dx, dt, k)
    
!         call abrasion (hn, dx, imax, dt, bedrock, maxVolumeLoss)
    call subsidence (hn, imax, dt, bedrock, subsidence_rate) !!! call subsedence must be before boundary conditions are applied !!!
    
    
    if (t > tprogress) then
        call progress (t,real(tmax),real(0))
        tprogress = tprogress + dtprogress
    end if
    
    do i = 0,imax
        do p=1,3
            h(i,p)= hn(i,p)          ! new values to become "old" in next step
        end do
        totalHeight(i) = hn(i,3)+hn(i,2)+hn(i,1)-2*bedrock(i) !!-3*bedrock for hn(1,2,3) and +1bedrock for totalHeight
    end do
!     print*, hn(5)-bedrock(5)
    t = t + dt
    
end do
! write(*,*) char(13) !This removes an artifact created by the progress bar


open (101,file= 'output.dat', status='replace', action='write')
do i= 0,imax
    write(101,*) x(i), h(i,1), h(i,2), h(i,3), bedrock(i), totalHeight(i)
end do
close (101)

!-----------------------------------------------

end ! end of main program

!-----------------------------------------------

subroutine diffusion(h, hn, bedrock, totalHeight, imax, dx, dt, k)

  integer :: imax, dx, i, p

  real*16 :: diffMagnitude, diffGradient
  real*16, dimension (0:imax) :: dh
  real*16, dimension (0:imax,3) :: h, hn
  real*16, dimension (0:imax) :: totalHeight
  real*16, dimension (0:imax) :: bedrock
  real*16, dimension (3) :: k
  real :: dt
  
  do i = 1,imax
      dh(i) = 0
  end do
  
  do i = 1,imax-1      ! find new elevation at all internal nodes
      do p=1,3            
!           hn(i,p)= h(i,p) + ((k(p)*dt*365.25*24*3600)/(dx*dx))*( h(i+1,p) - 2.*h(i,p) + h(i-1,p) )  !!Standard diffusion equation
          diffGradient = totalHeight(i+1) - 2*totalHeight(i) + totalHeight(i-1) 
          diffMagnitude = ((k(p)*dt*365.25*24*3600)/(dx*dx))
          
          if (totalHeight(i-1) > totalHeight(i+1)) then  !! If slope is negative
              hn(i,p)= h(i,p) + min(diffMagnitude*diffGradient, (h(i-1,p)-bedrock(i-1)) )
              hn(i-1,p)= hn(i-1,p) - min(diffMagnitude*diffGradient, (h(i-1,p)-bedrock(i-1)) )
          else   !! If slope is negative
              hn(i,p)= h(i,p) + min(diffMagnitude*diffGradient, (h(i+1,p)-bedrock(i+1)) )
              hn(i+1,p)= hn(i+1,p) - min(diffMagnitude*diffGradient, (h(i+1,p)-bedrock(i+1)) )
          endif
          dh(i) = dh(i) + hn(i,p) - h(i,p)
      end do
  end do
  
  return
end !!  End subroutine 'diffusion'

!-----------------------------------------------

subroutine abrasion(h, dx, imax, dt, bedrock, maxVolumeLoss)

  integer :: dx, imax, i, p
  real :: dt, maxVolumeLoss
  real*16, dimension (3) :: volumeLoss, massLoss
  real*16, dimension (0:imax,3) :: h
  real*16, dimension (0:imax) :: bedrock
  
!   do i=0,imax
!     h(i,3) = h(i,3) - 1.e-6*h(i,3)
!     h(i,2) = h(i,2) + 1.e-6*h(i,3) - 0.4e-5*h(i,2)
!     h(i,1) = h(i,1) + 0.4e-5*h(i,2)
!   end do
  maxVolumeLoss = 0
  do i = 1,imax       
    volumeLoss(3) = 1e-9*dx*(h(i,3)-bedrock(i))*dt     !gravel
!     volumeLoss(2) = 4.e-2*(h(i,2)-bedrock(i))*dt!5.e-6*h(i,2)*dt    !sand
    volumeLoss(1) = 0.!0.000001*h(i,1)*dt !silt
    
    if (maxval(volumeLoss) > maxVolumeLoss) maxVolumeLoss = maxval(volumeLoss)    !needed to determine max dt
    
!     if (i == 3) print*, h(i,2)
    h(i,3) = h(i,3)-volumeLoss(3)
    h(i,2) = h(i,2)-volumeLoss(2)+1.0*volumeLoss(3)
!     h(i,2) = h(i,2) + volumeLoss(3)-volumeLoss(2)!.05*(dt*4e-4/3)*(imax-i )
!     h(i,1) = h(i,1)+volumeLoss(2)!-volumeLoss(1)+0.1*volumeLoss(3)!+volumeLoss(2)
    h(i,1) = h(i,1)-volumeLoss(1)+0.0*volumeLoss(3)!+volumeLoss(2)
!     if (i == 3) print*, h(i,2)
  end do
  
  return
end !end abrasion subroutine

!-----------------------------------------------

subroutine subsidence (h, imax, dt, bedrock, subsidence_rate)

  integer :: imax, i, p
  real :: dt
  real*16 :: subsidence_rate
  real*16, dimension (0:imax,3) :: h
  real*16, dimension (0:imax) :: totalHeight
  real*16, dimension (0:imax) :: bedrock
  
  !note, the "/3" represents the three grain sizes being lowered an equal amount, totalling to subsidence. To be changed to relative proportions later
!   subsidence_rate = dt*4e-4/3 !m/yr
  
  do i = 0,imax
    bedrock(i) = bedrock(i)-subsidence_rate*(imax-i)
!     h(i,1) = bedrock(i)
!     h(i,2) = bedrock(i)
!     totalHeight(i) = totalHeight(i)-subsidence_rate*(imax-i )
    do p =1,3
      h(i,p) = max( h(i,p)-subsidence_rate*(imax-i) , bedrock(i))
!       h(i,p) = h(i,p)-subsidence_rate*(imax-i)
    end do 
  end do
!   do i=0,imax
!     h(i,1) = bedrock(i)
!     h(i,2) = bedrock(i)
!   end do 
  
  return
end !end subroutine "subsidence"

!-----------------------------------------------

subroutine writeOutput(n, imax, x, h, bedrock, totalHeight)
  implicit none
  
  character :: ext, ext2*2, ext3*3, ext4*4, ext5*5, ext6*6
  character*80 :: fna
  integer :: n, imax, i
  real*16, dimension (0:imax,3) :: h
  real*16, dimension (0:imax) :: totalHeight, bedrock, x
!   print*, n
  if (n.lt.10) then
    write (ext,'(i1)') n
    fna= 'output'//ext//'.xy'
  else if (n.ge.10.and.n.lt.100) then
    write (ext2,'(i2)') n
    fna= 'output'//ext2//'.xy'
  else if (n.ge.100.and.n.lt.1000) then
    write (ext3,'(i3)') n
    fna= 'output'//ext3//'.xy'
  else if (n.ge.1000.and.n.lt.10000) then
    write (ext4,'(i4)') n
    fna= 'output'//ext4//'.xy'
  else if (n.ge.10000.and.n.lt.100000) then
    write (ext5,'(i5)') n
    fna= 'output'//ext5//'.xy'
  else if (n.ge.100000.and.n.lt.1000000) then
    write (ext6,'(i6)') n
    fna= 'output'//ext6//'.xy'
  else
    write (*,*) 'error: file name cannot be constructed.'
    stop
  endif
  
  open (20,file= "outputData/"//fna)
  do i= 0,imax
    write (20,*) x(i), h(i,1), h(i,2), h(i,3), bedrock(i), totalHeight(i)
  end do
  close (20)
  
end
  
!-----------------------------------------------

subroutine progress(t,maximum,minimum)
  implicit none
  real ::maximum,minimum
  real :: t
  integer::k
  character(len=50)::bar="???% |                                        |     "
  write(unit=bar(1:3),fmt="(i3)") nint((100*(abs(t)-minimum)/(maximum-minimum)))
  do k=1, nint((40*(abs(t)-minimum)*(1/(maximum-minimum))))
    bar(6+k:6+k)="#"
  enddo
  ! print the progress bar.
  write(unit=6,fmt="(a1,a51)",advance="no") char(13),bar
  !print*, ((t-minimum)*(1/(maximum-minimum)))
  !(t-min)*(max/(max-min))
  return
end subroutine progress
