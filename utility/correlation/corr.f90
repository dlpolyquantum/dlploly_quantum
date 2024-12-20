
!*********************************************************************************
!
! DL_POLY Quantum additonal program for calculating velocity and dipole autocorrelation
! functions and resulting diffusion coefficients and absorption spectra, respectively
!
! Requires FFTW3 package for performing Fourier transform
!
! Authors - Nathan London and Dil Limbu
! 
!************************************************************************************

program correlation 
  use ISO_C_BINDING
  implicit none 

  !include 'fftw3.f03'
  include 'mpif.h'

  character(20) :: filenm
  character(5)  :: filenum
  character(8)  :: cortype
  logical       :: trimtraj
  integer       :: i,j,k,idnode,mxnode,ierr
  integer       :: traj0,traj1
  integer       :: wrtstep,wrtmol,wrttraj
  integer       :: nstep,ntraj,nfreq,nfreq2,nummols,keycorr,nfull,nstart,nend
  integer       :: ntotal,nskip
  real(8)      :: dt, domega,omegamax,omegaconv,pi,tau,domega2
  real(8), allocatable :: time(:),corr(:,:),corravg(:),error(:),dip(:,:,:),diff(:)
  real(8), allocatable :: dip_der(:,:,:)                         !omar
  real(8), allocatable ::spec(:,:),spectemp(:), window(:),windowfull(:)
  real(8), allocatable ::corravg0(:),spec0(:),corrsim(:),corrsim0(:),hold(:)
  real(8), allocatable ::arg(:),arg2(:),specInt(:),specInt2(:),timeSec(:)

  nskip=0

    
  open(1,file="CONTROL",status='old')

  read(1,*) ntotal, nskip
  read(1,*) nummols
  read(1,*) ntraj
  read(1,*) dt
  read(1,*) cortype
  read(1,*) tau
  read(1,*) trimtraj

  close(1)

  if(trim(cortype).eq.'velocity')then
    keycorr=1
  else
    keycorr=2
  endif
  
  nstep = ntotal-nskip

  dt=dt*1.d-12
  domega=1.d0/(dble(nstep+1)*dt)
  pi = 4.d0 * atan(1.d0)
  omegaconv = 29979245800.d0
  omegamax=5000 * omegaconv
  nfreq2 = 1000
  domega2 = omegamax/dble(nfreq2)

  nfreq=int(omegamax/domega)

  allocate(time(nstep+1))
  allocate(timeSec(nstep+1))
  allocate(corr(ntraj,nstep+1))
  allocate(error(nstep+1))
  allocate(spec(2,nfreq))
  allocate(windowfull(2*(nstep+1)+1))
  allocate(window(nstep+1))
  allocate(dip(nummols,nstep+1,3))
  allocate(dip_der(nummols,nstep+1,3))                     !omar
  allocate(spec0(nfreq))
  allocate(spectemp(nfreq))
  allocate(diff(nstep/2))
  allocate(hold(10))
  
  allocate(corravg(nstep+1))
  allocate(corravg0(nstep+1))
  allocate(corrsim(nstep+1))
  allocate(corrsim0(nstep+1))
  allocate(arg(nstep+1))
  allocate(specInt(nfreq2))
  allocate(arg2(nstep+1))
  allocate(specInt2(nfreq2))
  
  do j=1,nfreq
    spec(1,j)=(j-1)*domega/(29929245800.d0)
  enddo
  
  !tau = 40*1.d-12
  !tau = 1000*1d-12
  !windowfull(:)=1.d0
  tau = tau*1.d-12
  do j=1,nstep+1
    if((j*dt).le.tau) then
      window(j) = (cos(pi*(j*dt)/(2.0*tau)))**2
    else
      window(j) = 0.0
    endif
  !windowfull(j) = (cos(pi*(j*dt)/((nstep*dt))))**2
  
  enddo

  !window(:) = 0.0
  !window(:)=windowfull(0:(nstep+1)/2)
  corravg(:)=0.d0
  corrsim(:)=0.d0
  spec(2,:)=0.d0
  diff(:)=0.d0
 
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,idnode,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mxnode,ierr)
  
  traj0=(idnode*ntraj)/mxnode+1
  traj1=((idnode+1)*ntraj)/mxnode

  write(*,*) idnode,traj0
  write(*,*) idnode,traj1
  do i=traj0,traj1
   write(*,*)i
   write(filenum,'(i0.3)') i
    if(keycorr.eq.1)then
      filenm ='traj' // trim(filenum) // '/CORVEL'
    else if(keycorr.eq.2) then
      filenm ='traj' // trim(filenum) // '/CORDIP'
    endif
    write(*,*) idnode, filenm
    open((idnode+2)*mxnode,file=filenm,status='old')  
    read((idnode+2)*mxnode,*) cortype,hold(1)
    if(trimtraj) then
      do j=1,nskip
        read((idnode+2)*mxnode,'(1e14.6)') hold(1)
        do k=1,nummols
          read((idnode+2)*mxnode,'(3e14.6)') hold(2:4)
        enddo
      enddo
    endif
    do j=1,nstep+1
      read((idnode+2)*mxnode,'(1e14.6)') time(j)
      do k=1,nummols
        read((idnode+2)*mxnode,'(3e14.6)') dip(k,j,:)
      enddo
    enddo
    
    close((idnode+2)*mxnode)
  
    write(*,*) time(1) 
    if(trimtraj) then
        time(:) = time(:) - time(1)
    endif
    write(*,*) time(1) 
 
    if(keycorr.eq.1)then
      call velocity(nstep,nummols,i,idnode,time,dip,corravg,corrsim)
    else if(keycorr.eq.2)then
      dt = time(2) - time(1)                        ! time step assuming uniform time intervals
      do k=1,nummols                                ! dipole moment derivative calcultion using centered difference
          do j=2,nstep
              dip_der(k,j,:) = (dip(k,j+1,:) - dip(k,j-1,:)) / (2.0 * dt)
          enddo
          dip_der(k,1,:) = (dip(k,2,:) - dip(k,1,:)) / dt                ! forward difference for the first time step
          dip_der(k,nstep+1,:) = (dip(k,nstep+1,:) - dip(k,nstep,:)) / dt            ! backward difference for the last time step
      enddo
      call spectrum(idnode,mxnode,nstep,nfreq,nummols,i,time,window,dip_der,corravg,spec,corrsim)
    endif
  
  enddo

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  wrttraj=ntraj
  wrtmol=nummols
  wrtstep=nstep+1
  
  spectemp(:)=spec(2,:)

  if(keycorr.eq.1)then
    call MPI_REDUCE(corravg,corravg0,nstep+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(corrsim,corrsim0,nstep+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(idnode.eq.0)then
      corravg0(:)=corravg0(:)/dble(ntraj)
      corrsim0(:)=corrsim0(:)/dble(ntraj)/dble(nummols)
      do i=2,nstep/2
        call diffusion(i,time,corravg0,diff(i))
      enddo
    endif
  else if(keycorr.eq.2)then
    call MPI_REDUCE(corravg,corravg0,nstep+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(corrsim,corrsim0,nstep+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(spectemp,spec0,nfreq,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if(idnode.eq.0)then
      timeSec(:) = time(:)*1.d-12
      do i=1,nfreq2
        do j=1,nstep+1
          arg(j) = cos(domega2*(i-1)*timeSec(j)*2.0*pi)*corravg0(j)/dble(ntraj)
          arg2(j) =cos(domega2*(i-1)*timeSec(j)*2.0*pi)* &
            corrsim0(j)*window(j)/dble(ntraj)/dble(nummols)
        enddo
        call trapezoid_int(nstep+1,timeSec,arg,specInt(i))
        call trapezoid_int(nstep+1,timeSec,arg2,specInt2(i))
      enddo
      corravg0(:)=corravg0(:)/dble(ntraj)/(0.2081943**2)
      corrsim0(:)=corrsim0(:)/dble(ntraj)/dble(nummols)/(0.2081943**2)
      spec0(:)=spec0(:)/dble(ntraj)
    endif
  endif

  if(idnode.eq.0) open(3,file='CORAVG',status='replace')
  
  if(idnode.eq.0)then
    do i=1,nstep
      if(keycorr.eq.1)then
        write(3,'(500e14.6)') time(i),corravg0(i),corrsim0(i),window(i)
      else if(keycorr.eq.2)then
        write(3,'(500e14.6)') time(i),corravg0(i),corrsim0(i)*window(i),window(i)
      end if
    enddo
    close(3)
  endif

  if(idnode.eq.0.and.keycorr.eq.2)then
    open(4,file='SPEC',status='replace')

    do i=1,nfreq
      !write(4,'(2e14.6)') spec(1,i),spec0(i)*spec(1,i)**2
      write(4,'(2e14.6)') spec(1,i),spec0(i)
    enddo
    close(4)
  else if (idnode.eq.0.and.keycorr.eq.1)then
    open(4,file='DIFF',status='replace')

    do i=1,nstep/2
      write(4,'(2e14.6)') time(i),diff(i)
    enddo
    close(4)

  endif
  
  if(idnode.eq.0.and.keycorr.eq.2)then
    open(4,file='SPECint',status='replace')

    do i=1,nfreq2
      !write(4,'(3e14.6)') domega2*(i-1)/omegaconv,specInt(i)*(domega2*(i-1)/omegaconv)**2,&
      !  specInt2(i)*(domega2*(i-1)/omegaconv)**2
      write(4,'(3e14.6)') domega2*(i-1)/omegaconv,specInt(i),&
        specInt2(i)
    enddo
    close(4)
  endif

  deallocate(time,timeSec,corr,error,spec,windowfull,window,dip,spec0,spectemp,diff,hold)
  deallocate(corravg,corravg0,corrsim,corrsim0,arg2,specInt2)
  call mpi_finalize(ierr)

end program correlation

subroutine error_bar(nstep,ntraj,corr,corravg,error)
  implicit none
  
  integer, intent(in) :: nstep,ntraj
  integer :: i,j
  real(8), intent(in) :: corr(ntraj,nstep+1),corravg(nstep+1)
  real(8), intent(out) :: error(nstep+1)

  error(:)=0.d0

  do i=1,nstep+1
    do j=1,ntraj
      error(i)=error(i)+(corr(j,i)-corravg(i))**2
    enddo
    error(i)=sqrt(error(i))/dble(ntraj)
  enddo

end subroutine error_bar

subroutine diffusion(nstep,time,corravg,diff)
  implicit none
  
  integer,intent(in) :: nstep
  real(8),intent(in) :: time(nstep+1),corravg(nstep+1)
  real(8),intent(out) :: diff

  call trapezoid_int(nstep+1,time,corravg,diff)

  diff = diff/3.d0

end subroutine diffusion

subroutine trapezoid_int(n,x,func,integrand)
  implicit none

  integer :: i 
  integer :: n  !number of points
  real(8) :: x(n), func(n)  !the x coordinates and cooresponding function values
  real(8) :: integrand

  integrand = 0.0d0
  do i = 2, n
    integrand = integrand + (func(i-1)+func(i)) * abs(x(i)-x(i-1))/2.0d0
  enddo

end subroutine trapezoid_int

subroutine velocity(nstep,nummols,traj,idnode,time,dip,corravg,corrsim)

  integer :: j,k,l,nstep,nfreq,nummols,traj,idnode
  real(8) :: time(nstep),corravg(nstep),hold(nstep,3)
  real(8) :: dip(nummols,nstep,3),corrsim(nstep)
  real(8) :: pi, dt
  integer(8) :: planf1,planb,planf2
  complex(8) :: holdcfull(nstep,3),holdc(nstep,3),dipconv(nstep)
  character(20) :: filenm
  character(5)  :: filenum

  pi = 4.d0 * atan(1.d0)
 
      
  call dfftw_plan_dft_r2c_1d(planf1,nstep,hold(:,1),holdc(:,1),FFTW_MEASURE)
  call dfftw_plan_dft_c2r_1d(planb,nstep,holdc(:,1),hold(:,1),FFTW_MEASURE)

  hold(:,:)=0.d0
  holdc(:,:)=0.d0
  holdcfull(:,:)=0.d0

  dipconv(:)=0.d0

  do l=1,nummols
    hold(:,:)=dip(l,:,:)

    call dfftw_execute_dft_r2c(planf1,hold(:,1),holdc(:,1))
    call dfftw_execute_dft_r2c(planf1,hold(:,2),holdc(:,2))
    call dfftw_execute_dft_r2c(planf1,hold(:,3),holdc(:,3))

    do j=1,nstep
      do k=1,3
        dipconv(j) = dipconv(j) + holdc(j,k)*conjg(holdc(j,k))
        corrsim(j) = corrsim(j) + dip(l,j,k)*dip(l,1,k)
      enddo
    enddo
  
  dt=time(2)-time(1)
  
  dipconv(:)=dipconv(:)/(dble(nstep))
  holdc(:,1)=dipconv(:)

  call dfftw_execute(planb)
  hold(:,1) = hold(:,1)/dble(nstep)
  corravg(:) = corravg(:) + hold(:,1)/dble(nummols)
  corrsim(:) = corrsim(:)

  enddo
  write(filenum,'(i0.3)') traj
  filenm ='VELtraj' // trim(filenum)
  write(*,*) idnode, filenm
  open(idnode+5,file=filenm,status='replace')  
    do i=1,nstep
      write(idnode+5,'(2e14.6)') time(i),hold(i,1)
    enddo
  close(idnode+5)

  call dfftw_destroy_plan(planf1)  
  call dfftw_destroy_plan(planb)  

end subroutine velocity

subroutine spectrum(idnode,mxnode,nstep,nfreq,nummols,traj,time,window,dip,corravg,spec,corrsim)

  integer :: j,k,l,nstep,nfreq,nummols,traj,idnode,mxnode
  real(8) :: time(nstep+1),corravg(nstep+1),hold((nstep+1),3)
  real(8) :: dip(nummols,nstep+1,3),corrsim(nstep+1)
  real(8) :: pi, dt
  real(8) :: spec(2,nfreq),window(nstep+1)
  integer(8) :: planf1,planb,planf2
  complex(8) :: holdcfull((nstep+1),3),holdc((nstep+1),3),dipconv((nstep+1))
  character(20) :: filenm
  character(5)  :: filenum

  pi = 4.d0 * atan(1.d0)
 
  call dfftw_plan_dft_r2c_1d(planf1,(nstep+1),hold(:,1),holdc(:,1),FFTW_MEASURE)
  call dfftw_plan_dft_c2r_1d(planb,(nstep+1),holdc(:,1),hold(:,1),FFTW_MEASURE)
  call dfftw_plan_r2r_1d(planf2,(nstep+1),hold(:,1),hold(:,1),FFTW_R2HC,FFTW_MEASURE)
  !call dfftw_plan_dft_r2c_1d(planf1,nstep+1,hold(:,1),holdc(:,1),FFTW_MEASURE)
  !call dfftw_plan_dft_c2r_1d(planb,nstep+1,holdc(:,1),hold(:,1),FFTW_MEASURE)
  !call dfftw_plan_r2r_1d(planf2,(nstep+1),hold(:,1),hold(:,1),FFTW_R2HC,FFTW_MEASURE)

  hold(:,:)=0.d0
  holdc(:,:)=0.d0
  holdcfull(:,:)=0.d0

  dipconv(:)=0.d0

  do l=1,nummols
    hold(:,:)=dip(l,:,:)
    !write(*,*) hold(1:10,1)
    call dfftw_execute_dft_r2c(planf1,hold(:,1),holdc(:,1))
    call dfftw_execute_dft_r2c(planf1,hold(:,2),holdc(:,2))
    call dfftw_execute_dft_r2c(planf1,hold(:,3),holdc(:,3))

    do j=1,(nstep+1)
      do k=1,3
        dipconv(j) = dipconv(j) + holdc(j,k)*conjg(holdc(j,k))
        corrsim(j) = corrsim(j) + hold(j,k)*hold(1,k) 
      enddo
    enddo
  enddo
  
  dt=time(2)-time(1)
  
  dipconv(:)=dipconv(:)*dt/(2.d0*pi)/dble((nstep+1))/dble(nummols)
  holdc(:,1)=dipconv(:)

  call dfftw_execute(planb)
  hold(:,1) = hold(:,1)/dble((nstep+1))
  !  write(*,*) hold(1:10,1)
  !corravg(:) = corravg(:) + hold(:,1)/(dt/(2.d0*pi))
  !corrsim(:) = corrsim(:)

  do j=1,(nstep+1)
    hold(j,1)=hold(j,1)*window(j)
    corravg(j) = corravg(j) + hold(j,1)/(dt/(2.d0*pi))!*window(j)
    !corrsim(j) = corrsim(j)*window(j)
  enddo

  call dfftw_execute(planf2)
  !  write(*,*) hold(1:10,1)

  spec(2,:) = spec(2,:) + hold(1:nfreq,1)
    
  write(filenum,'(i0.3)') traj
  filenm ='traj' // trim(filenum) // '/SPECindiv'
  open((idnode+4)*mxnode,file=filenm,status='replace')  
    do i=1,nfreq
      write((idnode+4)*mxnode,'(2e14.6)') spec(1,i),hold(i,1)*spec(1,i)**2
    enddo
  close((idnode+4)*mxnode)

  call dfftw_destroy_plan(planf1)  
  call dfftw_destroy_plan(planb)  
  call dfftw_destroy_plan(planf2)  

end subroutine spectrum
