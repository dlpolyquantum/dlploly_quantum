      module correlation_module

c***********************************************************************
c     
c     dl_poly_quantum module for calculation real time correlation
c     functions
c
c     authors - Nathan London and Dil Limbu 2023
c
c***********************************************************************

      use setup_module, only: mspimd,nbeads,corr,mxbuff
      use config_module, only: xxx,yyy,zzz,vxx,vyy,vzz,buffer,cell,
     x                         rcell
      use site_module, only: chgsit,wgtsit,nummols,numsit
      use utility_module, only: sdot0,invert
      use error_module, only: error

      implicit none

      real(8), allocatable, save :: xoper(:),yoper(:),zoper(:)

      contains

      subroutine corr_init
     x  (idnode,mxnode,natms,ntpmls,itmols,keyens,keycorr,nummols,
     x  numsit)

c***********************************************************************
c     
c     dl_poly_quantum subroutine for initializing correlation function
c     calculation  
c
c     authors - Nathan London and Dil Limbu 2023
c
c***********************************************************************

      character*6 name
      integer, intent(in) :: idnode,mxnode,natms,ntpmls,itmols,keycorr
      integer, intent(in) :: keyens
      integer, intent(in) :: nummols(ntpmls),numsit(ntpmls)
      integer :: i

c     allocate operator value arrays

      allocate(xoper(nummols(itmols)))
      allocate(yoper(nummols(itmols)))
      allocate(zoper(nummols(itmols)))

      xoper=0.d0
      yoper=0.d0
      zoper=0.d0

c     open output file based on operator type

      if(idnode.eq.0)then
         
        if(keycorr.eq.1)then
          name='CORVEL'
          open(unit=corr,file=name,status='replace')
          write(corr,'(A8,I5)') "velocity",nummols(itmols)
        elseif(keycorr.eq.2)then
          name='CORDIP'
          open(unit=corr,file=name,status='replace')
          write(corr,'(A8,I5)') "dipole",nummols(itmols)
        endif
      endif

c     calculate initial value of operator

      if(keyens.eq.62)then
        call  calc_val_cent
     x    (idnode,mxnode,natms,ntpmls,itmols,keycorr,0,nummols,numsit)
      else 
        call calc_val_bead
     x    (idnode,mxnode,natms,ntpmls,itmols,keycorr,0,nummols,numsit)
      endif

c     write initial value to file

      if(idnode.eq.0)then
        
        write(corr,'(1e14.6)') 0.d0
c        do i=1,nummols(itmols)
c          write(corr,'(3e14.6)') xoper(i),yoper(i),zoper(i)
c        enddo
          write(corr,'(3e14.6)') sum(xoper(:)),sum(yoper(:)),
     x      sum(zoper(:))
         
      endif
      
      end subroutine corr_init
      
      subroutine correlation
     x  (idnode,mxnode,natms,ntpmls,itmols,keyens,keycorr,nstep,nummols,
     x  numsit,tstep)
c***********************************************************************
c     
c     dl_poly_quantum subroutine for doing correlation function
c     calculation during simulation
c
c     authors - Nathan London and Dil Limbu 2023
c
c***********************************************************************
      integer, intent(in) :: idnode,mxnode,natms,ntpmls,itmols,nstep
      integer, intent(in) :: keyens,keycorr
      integer, intent(in) :: nummols(ntpmls),numsit(ntpmls)
      real(8), intent(in) :: tstep
      integer :: i

c     calculate value of operator

      if(keyens.eq.62)then
        call  calc_val_cent
     x    (idnode,mxnode,natms,ntpmls,itmols,keycorr,nstep,nummols,
     x    numsit)
      else 
        call calc_val_bead
     x    (idnode,mxnode,natms,ntpmls,itmols,keycorr,nstep,nummols,
     x    numsit)
      endif
      
c     write value of operator to file
      if(idnode.eq.0)then
         
        write(corr,'(1e14.6)') nstep*tstep
c        do i=1,nummols(itmols)
c          write(corr,'(3e14.6)') xoper(i),yoper(i),zoper(i)
c        enddo
          write(corr,'(3e14.6)') sum(xoper(:)),sum(yoper(:)),
     x      sum(zoper(:))
         
      endif
      
      end subroutine correlation

      subroutine calc_val_cent
     x  (idnode,mxnode,natms,ntpmls,itmols,keycorr,nstep,nummols,numsit)
c***********************************************************************
c     
c     dl_poly_quantum subroutine for calculating value of operator
c     for PA-CMD simulations
c
c     authors - Nathan London and Dil Limbu 2023
c
c***********************************************************************

      integer, intent(in) :: idnode,mxnode,natms,ntpmls,itmols,nstep
      integer, intent(in) :: keycorr
      integer, intent(in) :: nummols(ntpmls),numsit(ntpmls)
      integer i,j,k,imol0,imol1,jsite,ksite
      real(8) datx(nbeads),daty(nbeads),datz(nbeads)
      real(8) mass(numsit(itmols)),charge(numsit(itmols))
      real(8) atmx(numsit(itmols)),atmy(numsit(itmols))
      real(8) atmz(numsit(itmols))
      real(8) molmass,avgx,avgy,avgz,xval,yval,zval
      real(8) comx,comy,comz

      imol0=(idnode*nummols(itmols))/mxnode+1
      imol1=((idnode+1)*nummols(itmols))/mxnode
    
      jsite=0
      ksite=0
      
c     get atom indices if specificied molecule is not first type
      if(itmols.gt.1)then
        do i=1,itmols-1
          jsite=jsite+nummols(i)*numsit(i)
          ksite=ksite+numsit(i)
        enddo
      endif
     
c     calculate molecular mass      
      molmass=0.d0
      do i=1,numsit(itmols)
        mass(i)=wgtsit(ksite+i)
        molmass=molmass+mass(i)
      enddo
      
c     get atom charges for dipole moment 
      if(keycorr.eq.2)then
        do i=1,numsit(itmols)
          charge(i)=chgsit(ksite+i)
        enddo
      endif

c     get velocities or postions depending on operator      
      do i=imol0,imol1
        xval=0.d0
        yval=0.d0
        zval=0.d0
        do j=1,numsit(itmols)
          do k=1,nbeads
            if(keycorr.eq.1)then 
              datx(k)=vxx((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              daty(k)=vyy((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              datz(k)=vzz((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
            elseif(keycorr.eq.2)then
              datx(k)=xxx((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              daty(k)=yyy((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              datz(k)=zzz((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
            endif
          enddo

c     bead average before operator
          call bead_average(datx,daty,datz,avgx,avgy,avgz)
          atmx(j)=avgx
          atmy(j)=avgy
          atmz(j)=avgz
        enddo

c     get molecule position with minimum image convention
        if(keycorr.eq.2) call mol_gather(numsit(itmols),atmx,atmy,atmz)

c     calculate center-of-mass
        do j=1,numsit(itmols)
          xval=xval+mass(j)*atmx(j)/molmass
          yval=yval+mass(j)*atmy(j)/molmass
          zval=zval+mass(j)*atmz(j)/molmass
        enddo

c     calculate dipole moment based on distance from center-of-mass        
        if(keycorr.eq.2)then
c          comx=xval
c          comy=yval
c          comz=zval
          comx=0.d0
          comy=0.d0
          comz=0.d0

          xval=0.d0
          yval=0.d0
          zval=0.d0

          do j=1,numsit(itmols)
            xval=xval+charge(j)*(atmx(j)-comx)
            yval=yval+charge(j)*(atmy(j)-comy)
            zval=zval+charge(j)*(atmz(j)-comz)
          enddo
        
        endif

c     store operater value for molecule        
          xoper(i)=xval
          yoper(i)=yval
          zoper(i)=zval
      enddo

c     sum over nodes
      if(mxnode.gt.1)then
         call merge
     x   (idnode,mxnode,nummols(itmols),mxbuff,xoper,yoper,zoper,buffer)
      endif
      
      end subroutine calc_val_cent

      subroutine calc_val_bead
     x  (idnode,mxnode,natms,ntpmls,itmols,keycorr,nstep,nummols,numsit)
c***********************************************************************
c     
c     dl_poly_quantum subroutine for calculating value of operator
c     for real-time dynamics simulations
c
c     authors - Nathan London and Dil Limbu 2023
c
c***********************************************************************

      integer, intent(in) :: idnode,mxnode,natms,ntpmls,itmols,nstep
      integer, intent(in) :: keycorr
      integer, intent(in) :: nummols(ntpmls),numsit(ntpmls)
      integer i,j,k,imol0,imol1,jsite,ksite
      real(8) datx(nbeads),daty(nbeads),datz(nbeads)
      real(8) mass(numsit(itmols)),charge(numsit(itmols))
      real(8) atmx(numsit(itmols)),atmy(numsit(itmols))
      real(8) atmz(numsit(itmols))
      real(8) molmass,avgx,avgy,avgz,xval,yval,zval
      real(8) comx,comy,comz

      imol0=(idnode*nummols(itmols))/mxnode+1
      imol1=((idnode+1)*nummols(itmols))/mxnode
      
      jsite=0
      ksite=0
c     get atom indices if specificied molecule is not first type
      if(itmols.gt.1)then
        do i=1,itmols-1
          jsite=jsite+nummols(i)*numsit(i)
          ksite=ksite+numsit(i)
        enddo
      endif
      
c     calculate molecular mass      
      molmass=0.d0
      do i=1,numsit(itmols)
        mass(i)=wgtsit(ksite+i)
        molmass=molmass+mass(i)
      enddo
      
c     get atom charges for dipole moment 
      if(keycorr.eq.2)then
        do i=1,numsit(itmols)
          charge(i)=chgsit(ksite+i)
        enddo
      endif
      
      do i=imol0,imol1
        xval=0.d0
        yval=0.d0
        zval=0.d0
        do k=1,nbeads
          datx(k)=0.d0
          daty(k)=0.d0
          datz(k)=0.d0
          do j=1,numsit(itmols)
            if(keycorr.eq.1)then 
              atmx(j)=vxx((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              atmy(j)=vyy((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              atmz(j)=vzz((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
            elseif(keycorr.eq.2)then
              atmx(j)=xxx((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              atmy(j)=yyy((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
              atmz(j)=zzz((k-1)*natms+(i-1)*numsit(itmols)+jsite+j)
            endif
          enddo

c     get molecule position with minimum image convention
          if(keycorr.eq.2)then
            call mol_gather(numsit(itmols),atmx,atmy,atmz)
          endif
          
c     get center-of-mass value          
c          if(keycorr.eq.1)then
          do j=1,numsit(itmols)
            datx(k)=datx(k)+mass(j)*atmx(j)/molmass
            daty(k)=daty(k)+mass(j)*atmy(j)/molmass
            datz(k)=datz(k)+mass(j)*atmz(j)/molmass
          enddo
c          endif

c     calculate dipole moment based on distance from center-of-mass        
          if(keycorr.eq.2)then
c            write(*,*) "atmx: ",atmx(:)
c            write(*,*) "charge: ",charge(:)
c          comx=xval
c          comy=yval
c          comz=zval
          comx=0.d0
          comy=0.d0
          comz=0.d0
c            write(*,*) "comx: ",comx
            datx(k)=0.d0
            daty(k)=0.d0
            datz(k)=0.d0

            do j=1,numsit(itmols)
              datx(k)=datx(k)+charge(j)*(atmx(j)-comx)
              daty(k)=daty(k)+charge(j)*(atmy(j)-comy)
              datz(k)=datz(k)+charge(j)*(atmz(j)-comz)
            enddo
c            write(*,*) "datx: ",datx
            comx=0.d0
            comy=0.d0
            comz=0.d0
c            write(*,*) "comx: ",comx
            datx(k)=0.d0
            daty(k)=0.d0
            datz(k)=0.d0

            do j=1,numsit(itmols)
              datx(k)=datx(k)+charge(j)*(atmx(j)-comx)
              daty(k)=daty(k)+charge(j)*(atmy(j)-comy)
              datz(k)=datz(k)+charge(j)*(atmz(j)-comz)
            enddo
c            write(*,*) "datx: ",datx
          
          endif

        enddo

c     calculate bead-averaged value of operator        
        call bead_average(datx,daty,datz,avgx,avgy,avgz)
        xval=avgx
        yval=avgy
        zval=avgz
          
c     store operater value for molecule        
        xoper(i)=xval
        yoper(i)=yval
        zoper(i)=zval
      enddo

c     sum over nodes      
      if(mxnode.gt.1)then
         call merge
     x   (idnode,mxnode,nummols(itmols),mxbuff,xoper,yoper,zoper,buffer)

      endif

      end subroutine calc_val_bead

      subroutine bead_average(datx,daty,datz,avgx,avgy,avgz)
c***********************************************************************
c     
c     dl_poly_quantum subroutine for calculating bead-average
c     of a quantity
c
c     authors - Nathan London and Dil Limbu 2023
c
c***********************************************************************

      implicit none
      integer i
      real(8), intent(in) :: datx(nbeads),daty(nbeads),datz(nbeads)
      real(8) :: avgx,avgy,avgz

      avgx=0.d0
      avgy=0.d0
      avgz=0.d0

      do i=1,nbeads
         avgx=avgx+datx(i)
         avgy=avgy+daty(i)
         avgz=avgz+datz(i)
      enddo
  
      avgx=avgx/dble(nbeads)
      avgy=avgy/dble(nbeads)
      avgz=avgz/dble(nbeads)

      end subroutine bead_average

      subroutine mol_gather(numsit,atmx,atmy,atmz)
c***********************************************************************
c     
c     dl_poly_quantum subroutine for returning unwrapped positions
c     of a molecule to ensure molecule is whole for calculating 
c     center-of-mass within minimum image convention   
c
c     authors - Nathan London and Dil Limbu 2023
c
c***********************************************************************
        
      integer, intent(in) :: numsit
      integer i
      real(8) :: atmx(numsit),atmy(numsit),atmz(numsit)
      real(8) :: dxx,dyy,dzz,sxx,syy,szz,det

      call invert(cell,rcell,det)

      sxx=atmx(1)*rcell(1)+atmy(1)*rcell(4)+atmz(1)*rcell(7)
      syy=atmx(1)*rcell(2)+atmy(1)*rcell(5)+atmz(1)*rcell(8)
      szz=atmx(1)*rcell(3)+atmy(1)*rcell(6)+atmz(1)*rcell(9)
      sxx=sxx-anint(sxx)
      syy=syy-anint(syy)
      szz=szz-anint(szz)
      atmx(1)=sxx*cell(1)+syy*cell(4)+szz*cell(7)
      atmy(1)=sxx*cell(2)+syy*cell(5)+szz*cell(8)
      atmz(1)=sxx*cell(3)+syy*cell(6)+szz*cell(9)

      do i=2,numsit
        
        dxx=atmx(i)-atmx(1)
        dyy=atmy(i)-atmy(1)
        dzz=atmz(i)-atmz(1)
        sxx=dxx*rcell(1)+dyy*rcell(4)+dzz*rcell(7)
        syy=dxx*rcell(2)+dyy*rcell(5)+dzz*rcell(8)
        szz=dxx*rcell(3)+dyy*rcell(6)+dzz*rcell(9)
        sxx=sxx-anint(sxx)
        syy=syy-anint(syy)
        szz=szz-anint(szz)
        dxx=sxx*cell(1)+syy*cell(4)+szz*cell(7)
        dyy=sxx*cell(2)+syy*cell(5)+szz*cell(8)
        dzz=sxx*cell(3)+syy*cell(6)+szz*cell(9)
        atmx(i)=atmx(1)+dxx
        atmy(i)=atmy(1)+dyy
        atmz(i)=atmz(1)+dzz

      enddo
       
      end subroutine mol_gather
      
      end module correlation_module

