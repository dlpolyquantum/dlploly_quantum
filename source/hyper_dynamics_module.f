      module hyper_dynamics_module
      
c***********************************************************************
c     
c     dl_poly module for defining hyperdynamics routines
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      use config_module
      use ensemble_tools_module
      use forces_module
      use nlist_builders_module
      use optimiser_module
      use property_module
      use setup_module
      use shake_module
      use temp_scalers_module
      use vv_motion_module
      
      implicit none
      
      integer, parameter :: mxtrn=10
      integer, parameter :: mxbsn=1000
      integer, parameter :: mxneb=8
      integer, parameter :: mxdiffs=300
      integer, parameter :: hyper_tag=35000      
      
      integer numbsn,numpro,numtrk,ndiff,maxtrk,numdark,home_bsn,numbpd
      integer nbsa(mxbsn),nbsb(mxbsn),ktrn(mxtrn)
      real(8) xtrn(mxtrn),ytrn(mxtrn)
      real(8) tstop,tkeres,timhyp,timres,tboost,boost,vbase
      
      integer, allocatable :: idabsn(:),keymin(:)
      real(8), allocatable :: xbas(:),ybas(:),zbas(:)
      real(8), allocatable :: xchk(:),ychk(:),zchk(:)
      real(8), allocatable :: xres(:),yres(:),zres(:)
      real(8), allocatable :: vxrs(:),vyrs(:),vzrs(:)
      real(8), allocatable :: fxrs(:),fyrs(:),fzrs(:)
      real(8), allocatable :: xhyp(:),yhyp(:),zhyp(:)
      real(8), allocatable :: vxhp(:),vyhp(:),vzhp(:)
      real(8), allocatable :: fxhp(:),fyhp(:),fzhp(:)
      real(8), allocatable :: xdiffs(:),ydiffs(:),zdiffs(:)
      real(8), allocatable :: celneb(:,:),path(:),optk(:,:)
      real(8), allocatable :: xneb(:),yneb(:),zneb(:),engneb(:)
      real(8), allocatable :: fxneb(:),fyneb(:),fzneb(:)
      real(8), allocatable :: hxneb(:),hyneb(:),hzneb(:)
      real(8), allocatable :: taux(:),tauy(:),tauz(:)
      real(8), allocatable :: track(:)

      integer bsn_1(maxneb),bsn_2(maxneb)
      real(8) strhyp(9),strres(9),engbsn(2)
      real(8) celbas(9),celhyp(9),celchk(9),celres(9)

      save numbsn,numtrk,numpro,ndiff,numdark,timres
      save xbas,ybas,zbas,xchk,ychk,zchk,timhyp,vbase
      save xres,yres,zres,vxrs,vyrs,vzrs,fxrs,fyrs,fzrs
      save xhyp,yhyp,zhyp,vxhp,vyhp,vzhp,fxhp,fyhp,fzhp
      save celbas,celhyp,celres,celchk,strhyp,strres
      save idabsn,nbsa,nbsb,xdiffs,ydiffs,zdiffs,tkeres
      save xneb,yneb,zneb,engneb,taux,tauy,tauz,keymin
      save fxneb,fyneb,fzneb,hxneb,hyneb,hzneb,path
      save optk,tstop,tboost,boost,numbpd
      
      contains
      
      subroutine alloc_hyper_arrays(idnode,mxnode)
      
c***********************************************************************
c     
c     dl_poly subroutine for defining hyperdynamics arrays and
c     initialising control variables
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=16
      
      logical safe
      integer i,fail,idnode,mxnode,nebmax
      dimension fail(nnn)

      safe=.true.
      nebmax=msatms*(mxneb+1)
      
c     initialise control variables
      
      numbpd=0
      numtrk=0
      numbsn=0
      numpro=0
      ndiff=0
      numdark=0
      home_bsn=0
      tkeres=0.d0
      tstop=1.d30
      timhyp=0.d0
      timres=0.d0
      boost=1.d0
      tboost=0.d0
      vbase=-huge(1.d0)
      do i=1,maxneb
        bsn_1(i)=0
        bsn_2(i)=0
      enddo
      
c     allocate working arrays

      fail(:)=0
      
      allocate (xbas(msatms),ybas(msatms),zbas(msatms),stat=fail(1))
      allocate (xchk(msatms),ychk(msatms),zchk(msatms),stat=fail(2))
      allocate (xres(msatms),yres(msatms),zres(msatms),stat=fail(3))
      allocate (vxrs(msatms),vyrs(msatms),vzrs(msatms),stat=fail(4))
      allocate (fxrs(msatms),fyrs(msatms),fzrs(msatms),stat=fail(5))
      allocate (xhyp(msatms),yhyp(msatms),zhyp(msatms),stat=fail(6))
      allocate (vxhp(msatms),vyhp(msatms),vzhp(msatms),stat=fail(7))
      allocate (fxhp(msatms),fyhp(msatms),fzhp(msatms),stat=fail(8))
      allocate (xdiffs(mxdiffs),ydiffs(mxdiffs),zdiffs(mxdiffs),
     x  stat=fail(9))
      allocate (idabsn(mxdiffs),keymin(0:mxneb),stat=fail(10))
      allocate (xneb(nebmax),yneb(nebmax),zneb(nebmax),stat=fail(11))
      allocate (fxneb(nebmax),fyneb(nebmax),fzneb(nebmax),stat=fail(12))
      allocate (taux(msatms),tauy(msatms),tauz(msatms),stat=fail(13))
      allocate (engneb(0:mxneb),celneb(9,0:mxneb),path(0:mxneb),
     x  stat=fail(14))
      allocate (hxneb(nebmax),hyneb(nebmax),hzneb(nebmax),stat=fail(15))
      allocate (optk(5,0:mxneb),stat=fail(16))
      
      if(any(fail.gt.0))safe=.false.      
      if(mxnode.gt.1)call gstate(safe)    
      if(.not.safe)call error(idnode,1115)
      
      end subroutine alloc_hyper_arrays
      
      subroutine hyper_start
     x  (ltad,lbpd,lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,idnode,
     x  imcon,keyfce,keyfld,keyshl,keytol,kmax1,kmax2,kmax3,multt,
     x  mxnode,natms,ngrp,nhko,nlatt,nneut,nospl,nscons,nstbgr,
     x  nstep,nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,
     x  ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,ntcons,ntrack,alpha,
     x  delr,dlrpot,drewd,elrc,virlrc,epsq,fmax,opttol,rctter,
     x  rcut,rcutfb,rcuttb,rprim,rvdw,temp,tstep,volm,sigma,
     x  hyp_units)
      
c***********************************************************************
c     
c     dl_poly subroutine for starting a hyperdynamics simulation
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2007
c     
c***********************************************************************

      implicit none
      
      logical ltad,lbpd,lfcap,lneut,lnsq,loglnk,lzeql,newlst,savflg
      integer nblock,idnode,imcon,keyfce,keyfld,keyshl,keytol
      integer kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt
      integer nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond
      integer ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,i,j
      integer ntpvdw,ntshl,ntteth,numblock,iatm0,iatm1,ntcons
      integer ktol,pass,fail,ntrack
      real(8) alpha,delr,dlrpot,drewd,elrc,epsq,fmax,opttol,rctter
      real(8) rcut,rcutfb,rcuttb,rprim,rvdw,temp,tstep,volm,engcfg
      real(8) virlrc,cvgerr,dum,otol,cgerr,sigma,hyp_units
      
c     allocate track array for BPD
      
      allocate (track(0:nblock/ntrack),stat=fail)
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     current block number
      
      numblock=nstep/nblock
      
c     open hyperdynamics journal file
      
      if(idnode.eq.0)
     x  open(nevnt,file="EVENTS",form='formatted',position='append')
      
c     set up hyperdynamics for simulation start
      
      if(nstep.eq.0)then
        
c     initialise bias potential boost factor
      
        numbpd=0
        boost=1.d0
        tboost=0.d0
      
c     set basin difference markers
        
        do i=1,mxbsn
          
          nbsa(i)=0
          nbsb(i)=0
          
        enddo
         
c     store the starting configuration
        
        savflg=.true.
        tkeres=sigma
        call store_config(savflg,idnode,mxnode,natms,strres,celres,
     x    xres,yres,zres,vxrs,vyrs,vzrs,fxrs,fyrs,fzrs)
        
c     minimise starting structure
        
        call define_minimum_state
     x    (lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,pass,
     x    idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons,
     x    kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,
     x    nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,
     x    ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,
     x    ntpvdw,ntshl,ntteth,alpha,delr,dlrpot,drewd,
     x    elrc,virlrc,epsq,fmax,opttol,rctter,rcut,rcutfb,
     x    rcuttb,rprim,rvdw,temp,tstep,volm,engcfg,cvgerr)
        
c     define zero energy for BPD dynamics mode
        
        vbase=engcfg
        
c     write events entry for minimisation
        
        if(idnode.eq.0)then
          
          write(nevnt,'("MIN",i10,3i6,1p,3e14.5)')
     x      nstep,pass,numblock,keytol,opttol/hyp_units,
     x      engcfg/hyp_units,cvgerr/hyp_units
          write(nrite,'(1x,"MIN",i10,3i6,1p,3e14.5)')
     x      nstep,pass,numblock,keytol,opttol/hyp_units,
     x      engcfg/hyp_units,cvgerr/hyp_units
          write(nrite,"(1x,120('-'))")
          
        endif
      
c     save minimised starting structure as basin file
        
        call write_reference_config
     x    ('CFGBSN','BASINS',nbsn,numbsn,natms,imcon,idnode,engcfg)
        
c     save details of starting home basin
        
        engbsn(1)=engcfg
        
        do i=1,9
          celbas(i)=cell(i)
        enddo
        
        call invert(cell,rcell,dum)
        
        j=0
        do i=iatm0,iatm1
          
          j=j+1
          xbas(j)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
          ybas(j)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
          zbas(j)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
          
        enddo
        
c     restore the starting configuration
        
        savflg=.false.
        call store_config(savflg,idnode,mxnode,natms,strres,celres,
     x    xres,yres,zres,vxrs,vyrs,vzrs,fxrs,fyrs,fzrs)
        
      else
        
c     restore previous data from hyperdynamics backup file
      
        call hyper_open(ltad,idnode,mxnode,natms,nsteql)
                
c     reset home basin for hyperdynamics (home basin is 0 for TAD)
        
        if(lbpd)home_bsn=numbsn-1
        
c     store the current configuration
        
        savflg=.true.
        call store_config(savflg,idnode,mxnode,natms,strhyp,celhyp,
     x    xhyp,yhyp,zhyp,vxhp,vyhp,vzhp,fxhp,fyhp,fzhp)
        
c     read minimised starting structure from home basin file
      
        call read_reference_config
     x    ('CFGBSN','BASINS',nbsn,home_bsn,natms,imcon,idnode,engcfg)
        
c     save details of current home basin
        
        engbsn(1)=engcfg
        
        do i=1,9
          celbas(i)=cell(i)
        enddo
        
        call invert(cell,rcell,dum)
        
        j=0
        do i=iatm0,iatm1
          
          j=j+1
          xbas(j)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
          ybas(j)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
          zbas(j)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
          
        enddo
        
c     restore the current configuration
        
        savflg=.false.
        call store_config(savflg,idnode,mxnode,natms,strhyp,celhyp,
     x    xhyp,yhyp,zhyp,vxhp,vyhp,vzhp,fxhp,fyhp,fzhp)
        
      endif
      
      return
      end subroutine hyper_start
      
      subroutine hyper_driver
     x  (seek,ltad,lbpd,recycle,lfcap,lneut,lnsq,loglnk,lzeql,
     x  newlst,prechk,tadall,nebgo,nblock,ntrack,idnode,imcon,
     x  keyfce,keyfld,keyshl,keytol,kmax1,kmax2,kmax3,multt,
     x  mxnode,natms,ngrp,ntcons,nhko,nlatt,nneut,nospl,nscons,
     x  nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,
     x  ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,blkout,
     x  alpha,delr,dlrpot,drewd,elrc,virlrc,epsq,fmax,
     x  opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,temp,
     x  tstep,volm,engcfg,catchrad,sprneb,deltad,tlow,engtke,
     x  tolnce,hyp_units,ebias,vmin)
      
c***********************************************************************
c     
c     dl_poly subroutine for implementing a hyperdynamics simulation
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2007
c     
c***********************************************************************
      
      implicit none

      character*8 seek
      logical lbpd,ltad,lfcap,lneut,lnsq,loglnk,lzeql,newlst,lneb
      logical lrefmin,same,savflg,recycle,scan,prechk,tadall,nebgo
      integer nblock,idnode,imcon,keyfce,keyfld,keyshl,keytol,ntrack
      integer kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt
      integer nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond
      integer ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw
      integer ntshl,ntteth,blkout,numblock,bsn1,bsn2,itrack
      integer nturn,ntcons,mdiff,newbsn,iatm0,iatm1,pass,i,j,itrk
      real(8) alpha,delr,dlrpot,drewd,elrc,epsq,fmax,opttol,rctter
      real(8) rcut,rcutfb,rcuttb,rprim,rvdw,temp,tstep,volm,catchrad
      real(8) cvgerr,estar,engcfg,cfgtmp,engcpe,engsrp,catch
      real(8) vircpe,engmet,virmet,virlrc,engtbp,virtbp,dum
      real(8) engfbp,virfbp,engter,virter,engbnd,virbnd,engang
      real(8) virang,engdih,virdih,enginv,virinv,engtet,virtet
      real(8) engshl,shlke,virshl,engfld,virfld,virsrp,sprneb
      real(8) deltad,deltal,tlow,timhop,timlow,engtke,tolnce,hyp_units
      real(8) ebias,vmin
      
      data bsn1,bsn2/0,0/
      
c     control variables
      
      lneb=.false.
      numblock=nstep/nblock
      maxtrk=nblock/ntrack
      if(numdark.eq.0)numdark=nsteql
      lrefmin=(mod(nstep,nblock).eq.0)
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     BPD/TAD simulation time
      
      timhyp=timhyp+tstep
      
c     track the tboost value
      
      if(mod(nstep,ntrack).eq.0)track(mod(numtrk,maxtrk))=tboost
      
c     provisional check for transition - compare current config with 
c     the reference state (not in dark period)
      
      same=.true.
      scan=.false.
      if(prechk.and.(mod(nstep,ntrack).eq.0).and.
     x  (lbpd.or.(ltad.and.(nstep.gt.numdark))))then
        
        catch=0.65d0*catchrad
        call check_for_transition
     x    (seek,same,scan,idnode,mxnode,natms,imcon,mdiff,nblock,
     x    catch)
        
        if(.not.same.and.idnode.eq.0)then
          
          write(nevnt,'("PRE",i10)')nstep
          write(nrite,'(1x,"PRE",i10)')nstep
          write(nrite,"(1x,120('-'))")
        
        endif
        
      endif
      
      if(.not.same.or.lrefmin)then
        
c     store the current configuration
        
        savflg=.true.
        call store_config(savflg,idnode,mxnode,natms,strhyp,celhyp,
     x    xhyp,yhyp,zhyp,vxhp,vyhp,vzhp,fxhp,fyhp,fzhp)
        
c     minimise current structure
        
        call define_minimum_state
     x    (lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,pass,
     x    idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons,
     x    kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,
     x    nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,
     x    ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,
     x    ntpvdw,ntshl,ntteth,alpha,delr,dlrpot,drewd,
     x    elrc,virlrc,epsq,fmax,opttol,rctter,rcut,rcutfb,
     x    rcuttb,rprim,rvdw,temp,tstep,volm,cfgtmp,cvgerr)
        
c     write events entry for minimisation
      
        if(idnode.eq.0)then
          
          write(nevnt,'("MIN",i10,3i6,1p,3e14.5)')
     x      nstep,pass,numblock,keytol,opttol/hyp_units,
     x      cfgtmp/hyp_units,cvgerr/hyp_units
          write(nrite,'(1x,"MIN",i10,3i6,1p,3e14.5)')
     x      nstep,pass,numblock,keytol,opttol/hyp_units,
     x      cfgtmp/hyp_units,cvgerr/hyp_units
          write(nrite,"(1x,120('-'))")
          
        endif
        
c     confirm any transition
        
        if(ltad)scan=.true.
        call check_for_transition
     x    (seek,same,scan,idnode,mxnode,natms,imcon,mdiff,nblock,
     x    catchrad)
        
c     transition detected - proceed with transition analysis
        
        if(.not.same)then
          
c     store new basin energy
          
          engbsn(2)=cfgtmp
          
c     save new minimised state (bias potential dynamics only)
          
          if(lbpd)then
            
            do i=1,9
              celres(i)=cell(i)
            enddo
            
            j=0
            do i=iatm0,iatm1
              
              j=j+1
              xres(j)=xxx(i)
              yres(j)=yyy(i)
              zres(j)=zzz(i)
              
            enddo

          endif
          
c     record transition (for TAD only if outside blackout period)
          
          if(lbpd.or.nstep.gt.numdark)then
            
c     check if transition results in unique new basin (TAD only)
            
            if(ltad)call check_basins(newbsn,mdiff,mxnode)
            
c     analysis of new basin
            
            if(lbpd.or.tadall.or.newbsn.eq.numbsn)then
              
c     set difference counters and pointers (TAD only)
              
              if(ltad)then
                
                if(numbsn.gt.mxbsn)call error(idnode,2330)
              
                ndiff=mdiff
                
                if(numbsn.gt.1)then
                  nbsa(numbsn)=nbsb(numbsn-1)+1
                else
                  nbsa(numbsn)=1
                endif
                
                nbsb(numbsn)=mdiff
                
              endif
              
c     save the basin file and store basin energy
              
              call write_reference_config
     x          ('CFGBSN','BASINS',nbsn,numbsn,natms,imcon,idnode,
     x          cfgtmp)
              
c     determine minimum (reaction) path and activation energy
              
              if(nebgo)call neb_driver
     x          (lfcap,lneut,lnsq,loglnk,lzeql,newlst,lneb,bsn1,
     x          bsn2,idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,
     x          keytol,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,ngrp,
     x          ntcons,ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,
     x          keyshl,ntfree,keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,
     x          ntshl,nscons,delr,dlrpot,engcpe,engsrp,epsq,rcut,
     x          rprim,rvdw,vircpe,virsrp,alpha,drewd,volm,
     x          engmet,virmet,elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,
     x          engfbp,virfbp,rctter,engter,virter,engbnd,virbnd,
     x          engang,virang,engdih,virdih,enginv,virinv,engtet,
     x          virtet,engshl,shlke,virshl,engfld,virfld,cfgtmp,fmax,
     x          temp,tstep,opttol,sprneb,hyp_units)
              
c     analyse the transition - determine e-star and destination state
              
              if(nebgo)call transition_properties
     x          (seek,ltad,lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,
     x          idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons,
     x          kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,
     x          nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,
     x          ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,
     x          ntpvdw,ntshl,ntteth,nturn,numbsn,alpha,delr,dlrpot,
     x          drewd,elrc,virlrc,epsq,fmax,opttol,rctter,rcut,rcutfb,
     x          rcuttb,rprim,rvdw,temp,tstep,volm,cfgtmp,cvgerr,estar,
     x          catchrad,hyp_units)
              
c     estimate time of transition from past trajectory
              
              call transition_time
     x          (seek,lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,
     x          idnode,imcon,keyfce,keyfld,keyshl,keytol,kmax1,
     x          kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,
     x          nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,
     x          ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,
     x          ntpter,ntrack,ntpvdw,ntshl,ntteth,ntcons,itrk,
     x          alpha,delr,dlrpot,drewd,elrc,virlrc,epsq,fmax,
     x          opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,temp,
     x          tstep,volm,cfgtmp,cvgerr,catchrad,timhop,hyp_units)
              
c     update TAD control variables
              
              if(ltad)then
                
c     update blackout period
                
                numdark=nblock*((nstep+blkout)/nblock+1)
                
c     calculate stopping time
                
                timlow=timhop*exp(-(estar/temp-estar/tlow)/boltz)
                tstop=min(tstop,deltad*(timlow/deltad)**(tlow/temp))
                
c     write transition data for TAD only
                
                if(idnode.eq.0)then
                  
                  write(nevnt,'("TRA",i10,3i6,1p,4e14.5)')
     x              nstep,home_bsn,numbsn-1,nturn,estar/hyp_units,
     x              timhop,timlow,tstop
                  write(nevnt,'("BLK",2i10)')nstep,numdark
                  write(nrite,'(1x,"TRA",i10,3i6,1p,4e14.5)')
     x              nstep,home_bsn,numbsn-1,nturn,estar/hyp_units,
     x              timhop,timlow,tstop
                  write(nrite,'(1x,"BLK",2i10)')nstep,numdark
                  write(nrite,"(1x,120('-'))")
                  
                endif
                
              elseif(nebgo.and.idnode.eq.0)then
                
c     write transition data for bias potential dynamics with NEB
                
                write(nevnt,'("TRA",i10,3i6,1p,3e14.5)')
     x            nstep,home_bsn,numbsn-1,nturn,estar/hyp_units,
     x            timhop,timhop*tboost
                write(nrite,'(1x,"TRA",i10,3i6,1p,3e14.5)')
     x            nstep,home_bsn,numbsn-1,nturn,estar/hyp_units,
     x            timhop,timhop*tboost
                write(nrite,"(1x,120('-'))")
                
                numbpd=0
                tboost=0.0
                
              elseif(idnode.eq.0)then
                
c     write transition data for bias potential dynamics without NEB
                
                write(nevnt,'("TRA",i10,2i6,1p,3e14.5)')
     x            nstep,home_bsn,numbsn-1,ebias/hyp_units,
     x            timhop,timhop*tboost
                write(nrite,'(1x,"TRA",i10,2i6,1p,3e14.5)')
     x            nstep,home_bsn,numbsn-1,ebias/hyp_units,
     x            timhop,timhop*tboost
                write(nrite,"(1x,120('-'))")
                
                numbpd=0
                tboost=0.0
                
              endif
              
            else
              
c     update blackout period when transition not unique (TAD only)
              
              numdark=nblock*((nstep+blkout)/nblock+1)
              
              if(idnode.eq.0)then
                
                write(nevnt,'("TRR",i10,3i6)')nstep,home_bsn,newbsn
                write(nevnt,'("BLK",2i10)')nstep,numdark
                write(nrite,'(1x,"TRR",i10,3i6)')nstep,home_bsn,newbsn
                write(nrite,'(1x,"BLK",2i10)')nstep,numdark
                write(nrite,"(1x,120('-'))")
                
              endif
              
            endif
            
c     actions when new basin has been visited before (TAD only)
            
          elseif(ltad)then
            
c     ignore dark transition and extend blackout period (TAD only)
            
            if(idnode.eq.0)then
              
              write(nevnt,'("TRI",i10)')nstep
              write(nrite,'(1x,"TRI",i10)')nstep
              write(nrite,"(1x,120('-'))")
              
            endif
            
            if(nstep.le.nsteql)then
              
              nsteql=nsteql+blkout
              numdark=nsteql+blkout
              if(idnode.eq.0)then
                
                write(nevnt,'("EQL",2i10)')nstep,nsteql
                write(nrite,'(1x,"EQL",2i10)')nstep,nsteql
                write(nrite,"(1x,120('-'))")
                
              endif
              
            else
              
              numdark=nblock*((nstep+blkout)/nblock+1)
              if(idnode.eq.0)then
                
                write(nevnt,'("BLK",3i10)')nstep,numdark
                write(nrite,'(1x,"BLK",3i10)')nstep,numdark
                write(nrite,"(1x,120('-'))")
                
              endif
              
            endif
            
          endif
          
          if(ltad)then
            
c     return to the block starting state after transition (TAD only)
            
            timhyp=timres
            savflg=.false.
            call store_config(savflg,idnode,mxnode,natms,strres,celres,
     x        xres,yres,zres,vxrs,vyrs,vzrs,fxrs,fyrs,fzrs)
            
c     scramble the velocities (and conserve system energy)
            
            call regauss(idnode,imcon,mxnode,natms,ngrp,nscons,ntcons,
     x        ntshl,keyshl,tkeres,temp,tolnce)
          
          elseif(lbpd)then
            
c     reset reference state to new basin (bias potential dynamics only)
            
            home_bsn=home_bsn+1
            engbsn(1)=engbsn(2)
      
            do i=1,9
              celbas(i)=celres(i)
            enddo
            
            call invert(celbas,rcell,dum)
            
            do i=1,iatm1-iatm0+1
              
              xbas(i)=rcell(1)*xres(i)+rcell(4)*yres(i)+rcell(7)*zres(i)
              ybas(i)=rcell(2)*xres(i)+rcell(5)*yres(i)+rcell(8)*zres(i)
              zbas(i)=rcell(3)*xres(i)+rcell(6)*yres(i)+rcell(9)*zres(i)
              
            enddo
            
c     restore current hyperdynamics configuration
            
            savflg=.false.
            call store_config(savflg,idnode,mxnode,natms,strhyp,celhyp,
     x        xhyp,yhyp,zhyp,vxhp,vyhp,vzhp,fxhp,fyhp,fzhp)
            
c     reset boost factor
            
            numbpd=0
            tboost=0.d0
            
          endif
          
c     no transition detected so restore current trajectory
          
        else
          
c     restore the current configuration
          
          savflg=.false.
          call store_config(savflg,idnode,mxnode,natms,strhyp,celhyp,
     x      xhyp,yhyp,zhyp,vxhp,vyhp,vzhp,fxhp,fyhp,fzhp)
          
c     save the block configuration as reset state (TAD only)
          
          if(ltad.and.lrefmin)then
            
            savflg=.true.
            tkeres=engtke
            timres=timhyp
            call store_config(savflg,idnode,mxnode,natms,strres,celres,
     x        xres,yres,zres,vxrs,vyrs,vzrs,fxrs,fyrs,fzrs)
            
          endif
          
        endif
        
      endif
      
c     close down if TAD stopping time reached
      
      if(ltad.and.tstop.lt.timhyp)recycle=.false.
            
c     write a tracking file
      
      if(mod(nstep,ntrack).eq.0)then
        
        itrack=mod(numtrk,maxtrk)
        call write_reference_config
     x    ('CFGTRK','TRACKS',ntrk,itrack,natms,imcon,idnode,engcfg)

        numtrk=numtrk+1
        
      endif
      
      return
      end subroutine hyper_driver
      
      subroutine define_minimum_state
     x  (lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,pass,
     x  idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons,
     x  kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,
     x  nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,
     x  ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,
     x  ntpvdw,ntshl,ntteth,alpha,delr,dlrpot,drewd,
     x  elrc,virlrc,epsq,fmax,opttol,rctter,rcut,rcutfb,
     x  rcuttb,rprim,rvdw,temp,tstep,volm,engcfg,cvgerr)

c***********************************************************************
c     
c     dl_poly subroutine for controlling subroutine calls for a
c     structural minimisation to define a minimum state for
c     hyperdynamics simulations
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c***********************************************************************
      
      implicit none

      logical lfcap,nogofr,lneut,lnsq,loglnk,stropt,lzeql
      logical newlst,ltad,lsolva,lfree,lexcite,lpimd,lfqcmd

      integer idnode,imcon,keyfce,keyfld,keyshl,keytol,nbeads
      integer keystr,kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp
      integer nhko,nlatt,nneut,nospl,nscons,nstbgr,nstep,nsteql
      integer ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet
      integer ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf
      integer nblock,pass,mxpass,mstep,ntcons,nsolva,isolva

      real(8) alpha,delr,dlrpot,drewd,elrc,engang,engbnd
      real(8) engcpe,engdih,engfbp,engfld,enginv,engshl,engsrp
      real(8) engtbp,engter,engtet,epsq,fmax,opttol,rctter
      real(8) rcut,rcutfb,rcuttb,rprim,rvdw,shlke,engcfg,temp
      real(8) tstep,virang,virbnd,vircpe,virdih,virfbp
      real(8) virfld,virinv,virlrc,virmet,virshl,virsrp
      real(8) virtbp,virter,virtet,volm,engmet,cfgold,cvgerr
      real(8) hnorm,grad0,grad1,ff1,sgn,engord,virord
      real(8) engrng,virrng,qmsbnd
      
      integer keyens

      data mxpass/1000/
      
c     control variables
      
      pass=0
      keystr=0
      numrdf=0
      ltad=.true.
      engcfg=1.d30
      stropt=.false.
      nogofr=.false.
      
c     dummy variables
      
      lpimd=.false.
      lfqcmd=.false.
      lsolva=.false.
      lfree=.false.
      lexcite=.false.
      nbeads=1
      nsolva=0
      isolva=1
      engord=0.d0
      virord=0.d0
      
c     relax the current structure
      
      do while(.not.stropt.and.pass.lt.mxpass)
        
        pass=pass+1
        cfgold=engcfg
        mstep=nstep+pass

c     construct verlet neighbour list
        
        call nlist_driver
     x    (newlst,lneut,lnsq,loglnk,ltad,natms,nbeads,idnode,mxnode,
     x    imcon,nneut,keyfce,rcut,delr,tstep)
        
c     calculate atomic forces
        
        call force_manager
     x    (newlst,lneut,lnsq,nogofr,lzeql,loglnk,lfcap,lsolva,lfree,
     x    lexcite,lpimd,idnode,mxnode,natms,imcon,mstep,nstbgr,nsteql,
     x    numrdf,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,
     x    ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,keyshl,
     x    keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,nsolva,
     x    isolva,nbeads,delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw,
     x    vircpe,virsrp,alpha,drewd,volm,engmet,virmet,elrc,virlrc,
     x    rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp,rctter,engter,
     x    virter,engbnd,virbnd,engang,virang,engdih,virdih,enginv,
     x    virinv,engtet,virtet,engshl,shlke,virshl,engfld,virfld,
     x    engcfg,fmax,temp,engord,virord,engrng,virrng,qmsbnd,keyens,
     x    lfqcmd)
        
c     frozen atoms option
        
        call freeze(natms)
        
c     structure optimisation
        
        call strucopt
     x    (stropt,keystr,keytol,idnode,mxnode,natms,ntcons,nscons,
     x    imcon,ngrp,ntfree,tstep,opttol,engcfg,hnorm,grad0,grad1,
     x    ff1,sgn)
        
        cvgerr=abs(engcfg-cfgold)
        
      enddo
            
      return
      end subroutine define_minimum_state
      
      subroutine write_reference_config
     x  (fnam,direct,nfil,nnn,natms,imcon,idnode,engcfg)
      
c***********************************************************************
c     
c     dlpoly utility to write a minimum structure file in CONFIG
c     format
c     
c     copyright daresbury laboratory
c     author      w.smith  june 2006
c     
c**********************************************************************
      
      implicit none
      
      character*6 fnam
      character*4 tail
      character*6 direct
      integer nfil,nnn,i,natms,imcon,idnode,levcfg
      real(8) engcfg
      
      levcfg=0
      
c     node zero handles i/o
      
      if(idnode.eq.0)then
        
c     write configuration data to new configuration file
        
        write(tail,'(i4.4)')nnn
        open(nfil,file=direct//'/'//fnam//tail,form='formatted')
        
        write(nfil,'(a10)')fnam//tail
        write(nfil,'(3i10,g20.12)') levcfg,imcon,natms,engcfg

        if(imcon.gt.0) write(nfil,'(3f20.12)') cell
        
        do i=1,natms
          
          write(nfil,'(a8,i10)') atmnam(i),i
          write(nfil,'(3g20.10)') xxx(i),yyy(i),zzz(i)
          
        enddo
        
        close (nfil)
        
      endif
      
      nnn=nnn+1
      
      return
      end subroutine write_reference_config
      
      subroutine read_reference_config
     x  (fnam,direct,nfil,nnn,natms,imcon,idnode,engcfg)
      
c***********************************************************************
c     
c     dlpoly utility to read a reference structure file in CONFIG
c     format
c     
c     copyright daresbury laboratory
c     author      w.smith  february 2007
c     
c**********************************************************************
      
      implicit none
      
      character*6 fnam
      character*4 tail
      character*6 direct
      integer nfil,nnn,i,natms,imcon,idnode,levcfg
      real(8) engcfg
      
c     node zero handles i/o
      
      if(idnode.eq.0)then
        
c     read configuration data from configuration file on proc 0
        
        write(tail,'(i4.4)')nnn
        open(nfil,file=direct//'/'//fnam//tail,form='formatted')
        
        read(nfil,*)
        read(nfil,'(3i10,g20.12)')levcfg,imcon,natms,engcfg
        buffer(1)=dble(levcfg)
        buffer(2)=dble(imcon)
        buffer(3)=dble(natms)
        buffer(4)=engcfg
        if(imcon.gt.0) read(nfil,'(3f20.12)') cell
        do i=1,9
          buffer(i+4)=cell(i)
        enddo
        call gdsum(buffer(1),13,buffer(14))
        
        do i=1,natms
          
          read(nfil,'(a8)') atmnam(i)
          read(nfil,'(3g20.10)') xxx(i),yyy(i),zzz(i)
          
        enddo
        
        close (nfil)
        
      else
        
c     gather data from configuration file on procs > 0
        
        do i=1,13
          buffer(i)=0.d0
        enddo
        call gdsum(buffer(1),13,buffer(14))

        levcfg=nint(buffer(1))
        imcon=nint(buffer(2))
        natms=nint(buffer(3))
        engcfg=buffer(4)
        do i=1,9
          cell(i)=buffer(i+4)
        enddo
        do i=1,natms
          xxx(i)=0.d0
          yyy(i)=0.d0
          zzz(i)=0.d0
        enddo
        
      endif
      
c     global gather of atomic coordinates
      
      call gdsum(xxx,natms,buffer)
      call gdsum(yyy,natms,buffer)
      call gdsum(zzz,natms,buffer)
      
      return
      end subroutine read_reference_config
      
      subroutine store_config(lsave,idnode,mxnode,natms,strold,celold,
     x  xold,yold,zold,vxold,vyold,vzold,fxold,fyold,fzold)

c***********************************************************************
c     
c     dlpoly hyperdynamics routine for storing the current 
c     configuration
c     
c     copyright daresbury laboratory
c     author      w.smith  sep  2006
c     
c**********************************************************************
      
      implicit none

      logical lsave
      integer idnode,mxnode,natms,iatm0,iatm1,i,j
      
      real(8) strold(9),celold(9)
      real(8) xold(msatms),yold(msatms),zold(msatms)
      real(8) vxold(msatms),vyold(msatms),vzold(msatms)
      real(8) fxold(msatms),fyold(msatms),fzold(msatms)
      
c     block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
      if(lsave)then
        
c     save cell and stress tensors
        
        do i=1,9
          
          celold(i)=cell(i)
          strold(i)=stress(i)
          
        enddo
        
c     save configuration
        
        j=0
        do i=iatm0,iatm1
          
          j=j+1
          xold(j)=xxx(i)
          yold(j)=yyy(i)
          zold(j)=zzz(i)
          vxold(j)=vxx(i)
          vyold(j)=vyy(i)
          vzold(j)=vzz(i)
          fxold(j)=fxx(i)
          fyold(j)=fyy(i)
          fzold(j)=fzz(i)
          
        enddo
        
      else
        
c     restore cell and stress tensors
        
        do i=1,9
          
          cell(i)=celold(i)
          stress(i)=strold(i)
          
        enddo
        
c     restore configuration
        
        j=0
        do i=iatm0,iatm1
          
          j=j+1
          xxx(i)=xold(j)
          yyy(i)=yold(j)
          zzz(i)=zold(j)
          vxx(i)=vxold(j)
          vyy(i)=vyold(j)
          vzz(i)=vzold(j)
          fxx(i)=fxold(j)
          fyy(i)=fyold(j)
          fzz(i)=fzold(j)
          
        enddo
        
c     replication of full configuration data
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
          call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
          call merge(idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
          
        endif
        
      endif
      
      return
      end subroutine store_config
      
      subroutine check_for_transition
     x  (seek,same,scan,idnode,mxnode,natms,imcon,mdiff,nblock,
     x  catchrad)
      
c***********************************************************************
c     
c     dlpoly hyperdynamics routine for checking when a transition
c     has occured in a configuration
c     
c     copyright daresbury laboratory
c     author      w.smith  sep  2006
c     
c**********************************************************************
      
      implicit none

      character*8 seek
      logical same,safe,scan,all
      integer idnode,mxnode,natms,imcon,nblock,mdiff
      integer iatm0,iatm1,i,j
      real(8) catchrad,catch2,rr2,dum,sxx,syy,szz,txx,tyy,tzz,pp2
      
      all=(seek.eq.'all     ')
      
c     flag for comparing structures
      
      same=.true.
      
c     block indices

      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     initialise search parameters
      
      catch2=catchrad**2
      
c     construct coordinate check arrays

      do i=1,9
        celchk(i)=cell(i)
      enddo

c     store structure in reduced coordinates (target atoms only)
      
      call invert(cell,rcell,dum)
      
      j=0
      do i=iatm0,iatm1
        
        j=j+1
        
        if(all.or.atmnam(i).eq.seek)then
          
          xchk(j)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
          ychk(j)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
          zchk(j)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
          
        endif
        
      enddo
      
c     compare current structure with reference basin
      
      j=0
      pp2=0.d0
      safe=.true.
      mdiff=ndiff
      do i=iatm0,iatm1
        
        j=j+1
        
        if(all.or.atmnam(i).eq.seek)then
          
c     calculate separations in reduced units
          
          sxx=xchk(j)-xbas(j)
          syy=ychk(j)-ybas(j)
          szz=zchk(j)-zbas(j)
          
c     calculate minimum image separations
          
          sxx=sxx-nint(sxx)
          syy=syy-nint(syy)
          szz=szz-nint(szz)
          
c     set trial structure at minimum displacements from reference
          
          xchk(j)=xbas(j)+sxx
          ychk(j)=ybas(j)+syy
          zchk(j)=zbas(j)+szz
          
c     calculate atomic separations in real coordinates
          
          txx=(celchk(1)*xchk(j)+celchk(4)*ychk(j)+celchk(7)*zchk(j))
     x      -(celbas(1)*xbas(j)+celbas(4)*ybas(j)+celbas(7)*zbas(j))
          tyy=(celchk(2)*xchk(j)+celchk(5)*ychk(j)+celchk(8)*zchk(j))
     x      -(celbas(2)*xbas(j)+celbas(5)*ybas(j)+celbas(8)*zbas(j))
          tzz=(celchk(3)*xchk(j)+celchk(6)*ychk(j)+celchk(9)*zchk(j))
     x      -(celbas(3)*xbas(j)+celbas(6)*ybas(j)+celbas(9)*zbas(j))
          
c     calculate total structural difference
          
          rr2=txx**2+tyy**2+tzz**2
          pp2=max(pp2,rr2)
          
c     identify and store the displaced atoms
          
          if(scan.and.rr2.ge.catch2)then
            
            mdiff=mdiff+1
            
            if(mdiff.le.mxdiffs)then
              
              idabsn(mdiff)=i
              xdiffs(mdiff)=txx
              ydiffs(mdiff)=tyy
              zdiffs(mdiff)=tzz
              
            else
              
              safe=.false.
              
            endif
            
          endif
          
        endif
        
      enddo
      
c     global check on diffs arrays
      
      if(scan)then
        
        if(mxnode.gt.1)call gstate(safe)
        if(.not.safe)then
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'number of current diffs',i10)")mdiff
          call error(idnode,2340)
          
        endif
        
      endif
      
c     global transition check
      
      same=(pp2.lt.catch2)
      if(mxnode.gt.1)call gstate(same)
      
      return
      end subroutine check_for_transition
      
      subroutine check_basins(newbsn,mdiff,mxnode)
      
c***********************************************************************
c     
c     dlpoly hyperdynamics routine for checking that a new basin is not
c     one already known
c     
c     copyright daresbury laboratory
c     author      w.smith  jan  2007
c     
c**********************************************************************
      
      implicit none
      
      logical same
      integer newbsn,ia,ib,ic,id,ibsn,i,j,k,mxnode,mdiff
      
      ibsn=1
      newbsn=0
      ib=mdiff
      ia=ndiff+1
      same=.false.
      do while(.not.same.and.ibsn.lt.numbsn)
        
        ic=nbsa(ibsn)
        id=nbsb(ibsn)
        
        if(ib-ia.eq.id-ic)then
          
          same=.true.
          
          do k=0,ib-ia
            
            i=ia+k
            j=ic+k
            
            if(.not.((idabsn(i).eq.idabsn(j)).and.
     x        (abs(xdiffs(i)-xdiffs(j)).lt.0.1d0).and.
     x        (abs(ydiffs(i)-ydiffs(j)).lt.0.1d0).and.
     x        (abs(zdiffs(i)-zdiffs(j)).lt.0.1d0)))same=.false.
            
          enddo
          
        endif
        
c     check if same on all processors
        
        if(mxnode.gt.1)call gstate(same)
        if(same)newbsn=ibsn
        
        ibsn=ibsn+1
        
      enddo
      
c     if not same - must be new basin!
      
      if(.not.same)newbsn=numbsn
      
      return
      end subroutine check_basins
      
      subroutine neb_driver
     x  (lfcap,lneut,lnsq,loglnk,lzeql,newlst,lneb,bsn1,
     x  bsn2,idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,
     x  keytol,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,ngrp,
     x  ntcons,ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,
     x  keyshl,ntfree,keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,
     x  ntshl,nscons,delr,dlrpot,engcpe,engsrp,epsq,rcut,
     x  rprim,rvdw,vircpe,virsrp,alpha,drewd,volm,
     x  engmet,virmet,elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,
     x  engfbp,virfbp,rctter,engter,virter,engbnd,virbnd,
     x  engang,virang,engdih,virdih,enginv,virinv,engtet,
     x  virtet,engshl,shlke,virshl,engfld,virfld,engcfg,fmax,
     x  temp,tstep,opttol,sprneb,hyp_units)
      
c***********************************************************************
c     
c     dl_poly subroutine for controlling a nudged elastic band 
c     calculation
c     
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2007
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lneb,newlst,lneut,lnsq,stropt
      logical lzeql,loglnk,lfcap
      integer idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,mstep
      integer keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw
      integer ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,keyshl
      integer keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,nscons
      integer ngrp,keytol,ntfree,iatm0,iatm1,iatm2,ntcons
      integer i,j,k,n,pass,mxpass,nkinks,bsn1,bsn2,itrk
      real(8) delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw
      real(8) vircpe,virsrp,alpha,drewd,volm,engmet,virmet
      real(8) elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp
      real(8) rctter,engter,virter,engbnd,virbnd,engang,virang
      real(8) engdih,virdih,enginv,virinv,engtet,virtet,engshl
      real(8) shlke,virshl,engfld,virfld,engcfg,fmax,temp,tstep
      real(8) sprneb,opttol,hyp_units,fac,xxn,yyn,zzn,tol,cvg,dum
      
      data mxpass/100/
      
c     control variables
      
      stropt=.false.
      do n=0,mxneb
        
        keymin(n)=0
        do i=1,5
          optk(i,n)=0.d0
        enddo
        
      enddo
      if(lneb)numpro=-(100*bsn1+bsn2)
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      iatm2=iatm1-iatm0+1
      
c     read in the required start and end point configurations
      
      if(lneb)then
        
c     read data for first reference structure
        
        call read_reference_config
     x    ('CFGBSN','BASINS',nbsn,bsn1,natms,imcon,idnode,engcfg)

        engbsn(1)=engcfg
        
        do i=1,9
          celbas(i)=cell(i)
        enddo
        call invert(cell,rcell,dum)
        
        j=0
        do i=iatm0,iatm1
          
          j=j+1
          xbas(j)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
          ybas(j)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
          zbas(j)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
          
        enddo
        
c     read data for second reference structure
        
        call read_reference_config
     x    ('CFGBSN','BASINS',nbsn,bsn2,natms,imcon,idnode,engcfg)

        engbsn(2)=engcfg
        
        do i=1,9
          celchk(i)=cell(i)
        enddo
        call invert(cell,rcell,dum)
        
        j=0
        do i=iatm0,iatm1
          
          j=j+1
          xchk(j)=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
          ychk(j)=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
          zchk(j)=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)
          
        enddo
        
      endif
      
c     construct initial configurations in chain
      
      k=0
      do n=0,mxneb
        
        j=0
        fac=dble(n)/dble(mxneb)
        
c     construct linear mix of cell vectors
        
        do i=1,9
          celneb(i,n)=(1.d0-fac)*celbas(i)+fac*celchk(i)
        enddo
        
c     construct configurations by linear interpolation
        
        do i=iatm0,iatm1
          
          j=j+1
          k=k+1

          xxn=xbas(j)+fac*((xchk(j)-xbas(j))-nint(xchk(j)-xbas(j)))
          yyn=ybas(j)+fac*((ychk(j)-ybas(j))-nint(ychk(j)-ybas(j)))
          zzn=zbas(j)+fac*((zchk(j)-zbas(j))-nint(zchk(j)-zbas(j)))
          xneb(k)=celneb(1,n)*xxn+celneb(4,n)*yyn+celneb(7,n)*zzn
          yneb(k)=celneb(2,n)*xxn+celneb(5,n)*yyn+celneb(8,n)*zzn
          zneb(k)=celneb(3,n)*xxn+celneb(6,n)*yyn+celneb(9,n)*zzn
  
        enddo
        
      enddo
      
c     start of NEB optimisation
      
      pass=0
      safe=.false.
      do while(.not.safe.and.pass.lt.mxpass)
        
        pass=pass+1
        safe=.true.
        mstep=nstep+pass
        
c     calculate system forces on all chain configurations
        
        call neb_system_forces
     x    (lfcap,lneut,lnsq,loglnk,lzeql,newlst,idnode,mxnode,
     x    natms,mstep,imcon,nstbgr,nsteql,keyfce,kmax1,kmax2,
     x    kmax3,nhko,nlatt,ntpvdw,ntpmet,nospl,multt,nneut,ntptbp,
     x    ntpfbp,ntpter,keyshl,keyfld,ntbond,ntangl,ntdihd,ntinv,
     x    ntteth,ntshl,delr,dlrpot,engcpe,engsrp,epsq,
     x    rcut,rprim,rvdw,vircpe,virsrp,alpha,drewd,volm,
     x    engmet,virmet,elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,
     x    engfbp,virfbp,rctter,engter,virter,engbnd,virbnd,
     x    engang,virang,engdih,virdih,enginv,virinv,engtet,
     x    virtet,engshl,shlke,virshl,engfld,virfld,engcfg,fmax,
     x    temp,tstep)
        
c     calculate spring forces on all chain configurations
        
        call neb_spring_forces(idnode,mxnode,natms,nkinks,sprneb)
        
c     energy minimisation of each chain configuration
        
        do n=0,mxneb
          
c     construct cell vectors for nth chain configuration
          
          do i=1,9
            cell(i)=celneb(i,n)
          enddo
          
c     construct coordinate and force arrays for nth chain configuration
          
          k=n*iatm2
          do i=iatm0,iatm1
            
            k=k+1
            xxx(i)=xneb(k)
            yyy(i)=yneb(k)
            zzz(i)=zneb(k)
            fxx(i)=fxneb(k)
            fyy(i)=fyneb(k)
            fzz(i)=fzneb(k)
            
          enddo
          
c     restore search direction vector if keymin > 0
          
          if(keymin(n).gt.0)then
            
            k=n*iatm2
            do i=iatm0,iatm1
              
              k=k+1
              hhx(i)=hxneb(k)
              hhy(i)=hyneb(k)
              hhz(i)=hzneb(k)
              
            enddo
            
          endif
          
c     form complete global arrays
          
          if(mxnode.gt.1)then
            
            call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
            call merge(idnode,mxnode,natms,mxbuff,fxx,fyy,fzz,buffer)
            if(keymin(n).gt.0)
     x        call merge(idnode,mxnode,natms,mxbuff,hhx,hhy,hhz,buffer)
            
          endif
          
c     structure optimisation
          
          call strucopt
     x      (stropt,keymin(n),keytol,idnode,mxnode,natms,ntcons,nscons,
     x      imcon,ngrp,ntfree,tstep,10.d0*opttol,engneb(n),optk(1,n),
     x      optk(2,n),optk(3,n),optk(4,n),optk(5,n))
          safe=safe.and.stropt
          stropt=.false.
          
c     update coordinate arrays for nth chain configuration
          
          k=n*iatm2
          do i=iatm0,iatm1
            
            k=k+1
            xneb(k)=xxx(i)
            yneb(k)=yyy(i)
            zneb(k)=zzz(i)
            hxneb(k)=hhx(i)
            hyneb(k)=hhy(i)
            hzneb(k)=hhz(i)
            
          enddo
          
        enddo
        
      enddo
      
c     convergence check
      
c$$$      if(.not.safe)then
c$$$        
c$$$        call error(idnode,2320)
c$$$        
c$$$      else
      
c     save neb profile
        
        call write_profile(idnode,mxnode,natms,hyp_units)
        
c     write neb summary
        
        if(idnode.eq.0)then
          
          if(lneb)then
            
            write(nrite,'(/,1x,"summary of NEB calculation",/)')
            write(nrite,'(1x,"path and energy for state",i4,
     x        " ---> state",i4," transition")')bsn1,bsn2
            write(nrite,'(1x,"convergence status :",l4)')safe
            write(nrite,'(1x,"obtained after ",i4," iterations",/)')pass
             
            do n=0,mxneb
              write(nrite,'(6x,1p,2e14.6)')path(n),engneb(n)/hyp_units
            enddo
           
          else
            
            write(nevnt,'("NEB",i10,3i6,1p,2e14.5)')
     x        nstep,pass,mxpass,mxneb+1,engbsn(1)/hyp_units,
     x        engbsn(2)/hyp_units
            write(nrite,'(1x,"NEB",i10,3i6,1p,2e14.5)')
     x        nstep,pass,mxpass,mxneb+1,engbsn(1)/hyp_units,
     x        engbsn(2)/hyp_units
            write(nrite,"(1x,120('-'))")
            
          endif
          
        endif
        
c$$$      endif
      
c     end of NEB optimisation
      
      return
      end subroutine neb_driver
      
      subroutine neb_system_forces
     x  (lfcap,lneut,lnsq,loglnk,lzeql,newlst,idnode,mxnode,
     x  natms,mstep,imcon,nstbgr,nsteql,keyfce,kmax1,kmax2,
     x  kmax3,nhko,nlatt,ntpvdw,ntpmet,nospl,multt,nneut,ntptbp,
     x  ntpfbp,ntpter,keyshl,keyfld,ntbond,ntangl,ntdihd,ntinv,
     x  ntteth,ntshl,delr,dlrpot,engcpe,engsrp,epsq,
     x  rcut,rprim,rvdw,vircpe,virsrp,alpha,drewd,volm,
     x  engmet,virmet,elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,
     x  engfbp,virfbp,rctter,engter,virter,engbnd,virbnd,
     x  engang,virang,engdih,virdih,enginv,virinv,engtet,
     x  virtet,engshl,shlke,virshl,engfld,virfld,engcfg,fmax,
     x  temp,tstep)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating system forces in a nudged 
c     elastic band calculation
c     
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2007
c     
c***********************************************************************
      
      implicit none
      
      logical newlst,lneut,lnsq,nogofr,lzeql,loglnk,lfcap,ltad
      logical lsolva,lfree,lexcite,lpimd,lfqcmd
      integer idnode,mxnode,natms,imcon,nstbgr,nsteql,mstep
      integer numrdf,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw
      integer ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,keyshl
      integer keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl
      integer iatm0,iatm1,iatm2,i,k,n,nsolva,isolva,nbeads
      real(8) delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw
      real(8) vircpe,virsrp,alpha,drewd,volm,engmet,virmet
      real(8) elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp
      real(8) rctter,engter,virter,engbnd,virbnd,engang,virang
      real(8) engdih,virdih,enginv,virinv,engtet,virtet,engshl
      real(8) shlke,virshl,engfld,virfld,engcfg,fmax,temp,tstep
      real(8) engord,virord,engrng,virrng,qmsbnd
      
      integer keyens

      numrdf=0
      ltad=.true.
      nogofr=.false.
      
c     dummy variables
      
      lpimd=.false.
      lfqcmd=.false.
      lsolva=.false.
      lfree=.false.
      lexcite=.false.
      nbeads=1
      nsolva=0
      isolva=1
      engord=0.d0
      virord=0.d0
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      iatm2=iatm1-iatm0+1
      
c     calculate system forces for all chain configurations
      
      do n=0,mxneb
        
c     construct cell vectors for one chain configuration
        
        do i=1,9
          cell(i)=celneb(i,n)
        enddo
        
c     construct coordinate array for one chain configuration
        
        k=n*iatm2
        do i=iatm0,iatm1
          
          k=k+1
          xxx(i)=xneb(k)
          yyy(i)=yneb(k)
          zzz(i)=zneb(k)
          
        enddo
        
c     form complete global arrays
        
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,natms,mxbuff,xxx,yyy,zzz,buffer)
        
c     construct verlet neighbour list
        
        call nlist_driver
     x    (newlst,lneut,lnsq,loglnk,ltad,natms,nbeads,idnode,mxnode,
     x    imcon,nneut,keyfce,rcut,delr,tstep)
        
c     calculate atomic forces for one chain configuration
        
        call force_manager
     x    (newlst,lneut,lnsq,nogofr,lzeql,loglnk,lfcap,lsolva,lfree,
     x    lexcite,lpimd,idnode,mxnode,natms,imcon,mstep,nstbgr,nsteql,
     x    numrdf,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,
     x    ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,keyshl,
     x    keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,nsolva,
     x    isolva,nbeads,delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw,
     x    vircpe,virsrp,alpha,drewd,volm,engmet,virmet,elrc,virlrc,
     x    rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp,rctter,engter,
     x    virter,engbnd,virbnd,engang,virang,engdih,virdih,enginv,
     x    virinv,engtet,virtet,engshl,shlke,virshl,engfld,virfld,
     x    engcfg,fmax,temp,engord,virord,engrng,virrng,qmsbnd,keyens,
     x    lfqcmd)
        
c     store configuration energy of chain configuration
        
        engneb(n)=engcfg
        
c     frozen atoms option
        
        call freeze(natms)
        
c     allocate forces to atoms of chain configuration
        
        k=n*iatm2
        do i=iatm0,iatm1
          
          k=k+1
          fxneb(k)=fxx(i)
          fyneb(k)=fyy(i)
          fzneb(k)=fzz(i)
          
        enddo
        
      enddo
      
      return
      end subroutine neb_system_forces
      
      subroutine neb_spring_forces(idnode,mxnode,natms,nkinks,sprneb)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating spring forces in a nudged 
c     elastic band calculation
c     
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2007
c     
c***********************************************************************
      
      implicit none
      
      integer i,j,k,n,kp,km,idnode,mxnode,natms,iatm0,iatm1,iatm2
      integer nkinks
      real(8) rp2,rm2,tau2,fpar,vv0,vp1,vm1,aaa,bbb,txx,tyy,tzz
      real(8) uxx,uyy,uzz,wxx,wyy,wzz,sxx,syy,szz,rxx,ryy,rzz
      real(8) sprneb,fac,kink,det
      real(8) rcella(9),rcellb(9),rcellc(9),cella(9),cellb(9),cellc(9)
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      iatm2=iatm1-iatm0+1
      
c     energies of first and last basins
      
      engneb(0)=engbsn(1)
      engneb(mxneb)=engbsn(2)
      
c     calculate spring tangents for all chain configurations
      
      nkinks=0
      do n=1,mxneb-1
        
        rp2=0.d0
        rm2=0.d0
        tau2=0.d0
        fpar=0.d0
        kink=0.d0
        vv0=engneb(n)
        vp1=engneb(n+1)
        vm1=engneb(n-1)
        
c     calculate tangent vector and scalar product with system force
        
        if(vp1.gt.vv0.and.vv0.gt.vm1)then
          
          aaa=1.d0
          bbb=0.d0
          
        else if(vp1.lt.vv0.and.vv0.lt.vm1)then
          
          aaa=0.d0
          bbb=1.d0
          
        else if(vp1.gt.vm1)then
          
          aaa=max(abs(vp1-vv0),abs(vv0-vm1))
          bbb=min(abs(vp1-vv0),abs(vv0-vm1))
          
        else
          
          aaa=min(abs(vp1-vv0),abs(vv0-vm1))
          bbb=max(abs(vp1-vv0),abs(vv0-vm1))
          
        endif
        
c     invert cell matrices
        
        do i=1,9
          cella(i)=celneb(i,n-1)
          cellb(i)=celneb(i,n)
          cellc(i)=celneb(i,n+1)
        enddo
        call invert(cella,rcella,det)
        call invert(cellb,rcellb,det)
        call invert(cellc,rcellc,det)

        j=0
        k=n*iatm2
        do i=iatm0,iatm1
          
          j=j+1
          k=k+1
          km=k-iatm2

c     calculate first spring vector (pbc corrected)

          sxx=rcellb(1)*xneb(k)+rcellb(4)*yneb(k)+rcellb(7)*zneb(k)
          syy=rcellb(2)*xneb(k)+rcellb(5)*yneb(k)+rcellb(8)*zneb(k)
          szz=rcellb(3)*xneb(k)+rcellb(6)*yneb(k)+rcellb(9)*zneb(k)
          rxx=rcella(1)*xneb(km)+rcella(4)*yneb(km)+rcella(7)*zneb(km)
          ryy=rcella(2)*xneb(km)+rcella(5)*yneb(km)+rcella(8)*zneb(km)
          rzz=rcella(3)*xneb(km)+rcella(6)*yneb(km)+rcella(9)*zneb(km)
          rxx=rxx-nint(rxx-sxx)
          ryy=ryy-nint(ryy-syy)
          rzz=rzz-nint(rzz-szz)
          txx=xneb(k)-
     x      (rxx*celneb(1,n-1)+ryy*celneb(4,n-1)+rzz*celneb(7,n-1))
          tyy=yneb(k)-
     x      (rxx*celneb(2,n-1)+ryy*celneb(5,n-1)+rzz*celneb(8,n-1))
          tzz=zneb(k)-
     x      (rxx*celneb(3,n-1)+ryy*celneb(6,n-1)+rzz*celneb(9,n-1))

c     calculate second spring vector (pbc corrected)

          kp=k+iatm2
          rxx=rcellc(1)*xneb(kp)+rcellc(4)*yneb(kp)+rcellc(7)*zneb(kp)
          ryy=rcellc(2)*xneb(kp)+rcellc(5)*yneb(kp)+rcellc(8)*zneb(kp)
          rzz=rcellc(3)*xneb(kp)+rcellc(6)*yneb(kp)+rcellc(9)*zneb(kp)
          rxx=rxx-nint(rxx-sxx)
          ryy=ryy-nint(ryy-syy)
          rzz=rzz-nint(rzz-szz)
          uxx=-xneb(k)+
     x      rxx*celneb(1,n+1)+ryy*celneb(4,n+1)+rzz*celneb(7,n+1)
          uyy=-yneb(k)+
     x      rxx*celneb(2,n+1)+ryy*celneb(5,n+1)+rzz*celneb(8,n+1)
          uzz=-zneb(k)+
     x      rxx*celneb(3,n+1)+ryy*celneb(6,n+1)+rzz*celneb(9,n+1)

          rp2=rp2+uxx*uxx+uyy*uyy+uzz*uzz
          rm2=rm2+txx*txx+tyy*tyy+tzz*tzz
          wxx=aaa*uxx+bbb*txx
          wyy=aaa*uyy+bbb*tyy
          wzz=aaa*uzz+bbb*tzz
          taux(j)=wxx
          tauy(j)=wyy
          tauz(j)=wzz
          tau2=tau2+wxx*wxx+wyy*wyy+wzz*wzz
          fpar=fpar+wxx*fxneb(k)+wyy*fyneb(k)+wzz*fzneb(k)
          kink=kink+txx*uxx+tyy*uyy+tzz*uzz
          
        enddo
        
        if(mxnode.gt.1)then
          
          buffer(1)=rp2
          buffer(2)=rm2
          buffer(3)=tau2
          buffer(4)=fpar
          buffer(5)=kink
          call gdsum(buffer(1),5,buffer(6))
          rp2=buffer(1)
          rm2=buffer(2)
          tau2=buffer(3)
          fpar=buffer(4)
          kink=buffer(5)
          
        endif
        
c     check for kinking of NEB
        
        kink=cos(kink/sqrt(rp2*rm2))
        if(kink.lt.0.5d0)nkinks=nkinks+1
        
c     calculate final forces
        
        j=0
        k=n*iatm2
        tau2=sqrt(tau2)
        fac=(sprneb*(sqrt(rp2)-sqrt(rm2))-fpar/tau2)/tau2
        do i=iatm0,iatm1
          
          j=j+1
          k=k+1
          fxneb(k)=fxneb(k)+fac*taux(j)
          fyneb(k)=fyneb(k)+fac*tauy(j)
          fzneb(k)=fzneb(k)+fac*tauz(j)
          
        enddo
        
      enddo
      
c     abort if kinks detected
      
      if(nkinks.gt.0)then
        
        if(idnode.eq.0)
     x    write(nrite,'(1x,"number of kinks detected ",i6)')nkinks
        call error(idnode,2350)
        
      endif
      
      return
      end subroutine neb_spring_forces
      
      subroutine transition_properties
     x  (seek,ltad,lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,
     x  idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons,kmax1,
     x  kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,nneut,nospl,
     x  nscons,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,ntfree,
     x  ntinv,ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,
     x  nturn,numbsn,alpha,delr,dlrpot,drewd,elrc,virlrc,epsq,
     x  fmax,opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,temp,
     x  tstep,volm,cfgtmp,cvgerr,estar,catchrad,hyp_units)
      
c***********************************************************************
c     
c     dl_poly subroutine for analysing the NEB path and determining
c     the destination state (if not end of chain).
c     
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2007
c     
c***********************************************************************
      
      implicit none
      
      character*8 seek
      logical lfcap,lneut,lnsq,loglnk,lzeql,newlst,ltad,scan,same
      integer nblock,idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons
      integer kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt
      integer nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond
      integer ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,pass
      integer ntpvdw,ntshl,ntteth,nturn,i,k,n,iatm0,iatm1,iatm2
      integer numblock,numbsn,mdiff
      real(8) alpha,delr,dlrpot,drewd,elrc,epsq,fmax,opttol,rctter
      real(8) rcut,rcutfb,rcuttb,rprim,rvdw,temp,tstep,volm,cfgtmp
      real(8) virlrc,cvgerr,estar,catchrad,hyp_units
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      iatm2=iatm1-iatm0+1
      
c     analyse neb profile
      
      call scan_profile(nturn,estar)
      
c     determine true new state from first maximum
      
      if(nturn.gt.1)then
        
        i=1
        do while(ktrn(i).ge.0)
          i=i+1
        enddo
        n=-ktrn(i)
        
c     construct cell vectors for nth chain configuration
        
        do i=1,9
          cell(i)=celneb(i,n)
        enddo
        
c     construct coordinate force arrays for nth chain configuration
        
        k=n*iatm2
        do i=iatm0,iatm1
          
          k=k+1
          xxx(i)=xneb(k)
          yyy(i)=yneb(k)
          zzz(i)=zneb(k)
          
        enddo
        
c     now minimise structure - this is correct new state
        
        call define_minimum_state
     x    (lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,pass,
     x    idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons,
     x    kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,
     x    nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,
     x    ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,
     x    ntpvdw,ntshl,ntteth,alpha,delr,dlrpot,drewd,
     x    elrc,virlrc,epsq,fmax,opttol,rctter,rcut,rcutfb,
     x    rcuttb,rprim,rvdw,temp,tstep,volm,cfgtmp,cvgerr)
        
c     write events entry  for minimisation
        
        numblock=nstep/nblock
        
        if(idnode.eq.0)then
          
          write(nevnt,'("MIN",i10,3i6,1p,3e14.5)')
     x      nstep,pass,numblock,keytol,opttol/hyp_units,
     x      cfgtmp/hyp_units,cvgerr/hyp_units
          write(nrite,'(1x,"MIN",i10,3i6,1p,3e14.5)')
     x      nstep,pass,numblock,keytol,opttol/hyp_units,
     x      cfgtmp/hyp_units,cvgerr/hyp_units
          write(nrite,"(1x,120('-'))")
          
        endif
        
        if(ltad)then
          
c     determine differences for new state (TAD only)
        
          scan=.true.
          call check_for_transition
     x      (seek,same,scan,idnode,mxnode,natms,imcon,mdiff,nblock,
     x      catchrad)
          
c     set difference counters and pointers
          
          if(numbsn.gt.mxbsn)call error(idnode,2330)
          
          ndiff=mdiff
          
          if(numbsn.gt.1)then
            nbsa(numbsn)=nbsb(numbsn-1)+1
          else
            nbsa(numbsn)=1
          endif
          
          nbsb(numbsn)=mdiff
          
c     save minimised starting structure as basin file
          
          call write_reference_config
     x      ('CFGBSN','BASINS',nbsn,numbsn,natms,imcon,idnode,cfgtmp)
          
        endif
        
      endif
      
      return
      end subroutine transition_properties
      
      subroutine write_profile(idnode,mxnode,natms,hyp_units)
      
c***********************************************************************
c     
c     dl_poly subroutine for writing profile file for NEB path
c     
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2007
c     
c***********************************************************************
      
      character*4 tail
      integer idnode,mxnode,natms,i,j,k,n,iatm0,iatm1,iatm2
      real(8) hyp_units
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      iatm2=iatm1-iatm0+1
      
      if(idnode.eq.0)then
        
c     open profile data file
      
        if(numpro.lt.0)then

          write(tail,'(i4.4)')iabs(numpro)
          open(npro,file='PROFILES'//'/'//'PRX'//tail//'.XY',
     x      form='formatted')
          
        else
          
          write(tail,'(i4.4)')numpro
          open(npro,file='PROFILES'//'/'//'PRO'//tail//'.XY',
     x      form='formatted')
          
        endif
        
      endif
      
c     calculate path
      
      path(0)=0.d0
      if(idnode.eq.0)write(npro,'(1p,2e14.6)')path(0),
     x  engneb(0)/hyp_units
      
      do n=1,mxneb
        
        path(n)=0.d0
        k=n*iatm2
        j=(n-1)*iatm2
        do i=iatm0,iatm1
          
          j=j+1
          k=k+1
          path(n)=(xneb(k)-xneb(j))**2+(yneb(k)-yneb(j))**2+
     x      (zneb(k)-zneb(j))**2+path(n)
          
        enddo
        
        if(mxnode.gt.1)call gdsum(path(n),1,buffer(1))
        
        path(n)=sqrt(path(n))+path(n-1)
        if(idnode.eq.0)write(npro,'(1p,2e14.6)')path(n),
     x    engneb(n)/hyp_units
        
      enddo
      
      numpro=numpro+1
      
      if(idnode.eq.0)close(npro)
      
      return
      end subroutine write_profile
      
      subroutine scan_profile(nturn,estar)
      
c*********************************************************************
c     
c     dl_poly  routine for analysing neb energy profile
c     
c     copyright - daresbury laboratory
c     author    - w.smith january 2007
c     
c*********************************************************************

      implicit none
      
      integer, parameter :: nscan=100

      integer i,np,n1,n2,npnts,fail,nturn
      real(8) di,dj,rpd,uu,vv,v0,ss,estar
      real(8), allocatable :: aa(:),dd(:),gg(:),zz(:)
      
c     allocate working arrays
      
      allocate (aa(0:mxneb),dd(0:mxneb),gg(0:mxneb),zz(0:mxneb),
     x  stat=fail)
      
      npnts=mxneb+1
      n1=npnts-1
      n2=npnts-2
      
c     calculate spline coefficients
      
      gg(0)=0.d0
      dd(0)=path(1)-path(0)
      
      do i=1,n1-1
        
        dd(i)=path(i+1)-path(i)
        gg(i)=2.d0*(path(i+1)-path(i-1))
        zz(i)=6.d0*((engneb(i+1)-engneb(i))/dd(i)-
     x    (engneb(i)-engneb(i-1))/dd(i-1))
        
      enddo
      
      gg(n1)=0.d0
      dd(n1)=0.d0
      aa(0)=0.d0
      aa(1)=dd(1)/gg(1)
      
      do i=2,n2-1
        
        gg(i)=gg(i)-dd(i-1)*aa(i-1)
        aa(i)=dd(i)/gg(i)
        
      enddo
      
      gg(n1-1)=gg(n1-1)-dd(n2-1)*aa(n2-1)
      gg(1)=zz(1)/gg(1)
      
      do i=2,n1-1
        gg(i)=(zz(i)-dd(i-1)*gg(i-1))/gg(i)
      enddo
      
      do i=1,n2-1
        gg(n1-i)=gg(n1-i)-aa(n1-i)*gg(npnts-i)
      enddo
      
c     now scan across the profile locating maxima and minima
            
      np=1
      nturn=0
      ss=1.d0
      v0=engneb(0)
      rpd=(path(npnts-1)-path(0))/dble(nscan)
      
      do i=2,nscan-1
        
        uu=rpd*dble(i)+path(0)
        
        do while(np.lt.npnts.and.uu.gt.path(np))
          np=np+1
        enddo

        di=uu-path(np-1)
        dj=path(np)-uu
        vv=(di*engneb(np)+dj*engneb(np-1)-di*dj*
     x    ((dd(np-1)+dj)*gg(np-1)+(dd(np-1)+di)*gg(np))/6.d0)/dd(np-1)
        
        if(ss.gt.0.d0.and.vv.le.v0)then
          
          nturn=nturn+1
          xtrn(nturn)=uu
          ytrn(nturn)=vv
          ktrn(nturn)=np
          
        else if(ss.lt.0.d0.and.vv.gt.v0)then
          
          nturn=nturn+1
          xtrn(nturn)=uu
          ytrn(nturn)=vv
          ktrn(nturn)=-np
          
        endif
        
        ss=sign(1.d0,vv-v0)
        v0=vv
        
      enddo
      
c     estimated activation energy
      
      i=1
      do while(ktrn(i).lt.0)
        i=i+1
      enddo
      estar=ytrn(i)-engbsn(1)
      
c     deallocate working arrays
      
      deallocate (aa,dd,gg,zz,stat=fail)
      
      return
      end subroutine scan_profile
      
      subroutine transition_time
     x  (seek,lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,idnode,
     x  imcon,keyfce,keyfld,keyshl,keytol,kmax1,kmax2,kmax3,
     x  multt,mxnode,natms,ngrp,nhko,nlatt,nneut,nospl,nscons,
     x  nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,
     x  ntpfbp,ntpmet,ntptbp,ntpter,ntrack,ntpvdw,ntshl,ntteth,
     x  ntcons,itrk,alpha,delr,dlrpot,drewd,elrc,virlrc,epsq,
     x  fmax,opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,temp,
     x  tstep,volm,cfgtmp,cvgerr,catchrad,timhop,hyp_units)

c*********************************************************************
c     
c     dl_poly  routine for estimating the time of a transition
c     from a backlog of previous configurations
c     
c     copyright - daresbury laboratory
c     author    - w.smith february 2007
c     
c*********************************************************************
      
      implicit none
      
      character*8 seek
      logical same,minflg,lfcap,lneut,lnsq,loglnk,scan
      logical lzeql,newlst
      integer nblock,idnode,imcon,keyfce,keyfld,keyshl,keytol
      integer kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt
      integer nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond
      integer ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter
      integer ntpvdw,ntshl,ntteth,ntcons,ntrack,itrk,mdiff,pass
      integer nback,ntry,numblock
      real(8) alpha,delr,dlrpot,drewd,elrc,virlrc,epsq,fmax,opttol
      real(8) rctter,rcut,rcutfb,rcuttb,rprim,rvdw,temp,tstep
      real(8) volm,cfgtmp,cvgerr,eng,tol,catchrad,timhop,hyp_units
      
c     determine starting tracking file
      
      nback=min(numtrk,maxtrk)
      
c     search track files for transition
        
      itrk=0
      same=.false.
      
      do while(.not.same.and.itrk.le.nback)
        
        itrk=itrk+1
        ntry=mod(numtrk-itrk,maxtrk)
        
        call read_reference_config
     x    ('CFGTRK','TRACKS',ntrk,ntry,natms,imcon,idnode,eng)
        
        call define_minimum_state
     x    (lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,pass,
     x    idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons,
     x    kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,
     x    nneut,nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,
     x    ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,
     x    ntpvdw,ntshl,ntteth,alpha,delr,dlrpot,drewd,
     x    elrc,virlrc,epsq,fmax,opttol,rctter,rcut,rcutfb,
     x    rcuttb,rprim,rvdw,temp,tstep,volm,cfgtmp,cvgerr)
        
c     write events entry for minimisation (normally deactivated)
        
c$$$        if(idnode.eq.0)then
c$$$          
c$$$          numblock=nstep/nblock
c$$$          write(nevnt,'("MIN",i10,3i6,1p,3e14.5)')
c$$$     x      nstep,pass,numblock,keytol,opttol/hyp_units,
c$$$     x      cfgtmp/hyp_units,cvgerr/hyp_units
c$$$          write(nrite,'(1x,"MIN",i10,3i6,1p,3e14.5)')
c$$$     x      nstep,pass,numblock,keytol,opttol/hyp_units,
c$$$     x      cfgtmp/hyp_units,cvgerr/hyp_units
c$$$          write(nrite,"(1x,120('-'))")
c$$$          
c$$$        endif
        
c     check if still in base state
        
        scan=.false.
        call check_for_transition
     x    (seek,same,scan,idnode,mxnode,natms,imcon,mdiff,nblock,
     x    catchrad)
        
      enddo
      
      timhop=timhyp-tstep*dble(ntrack)*(dble(itrk)-0.5d0)
      tboost=track(ntry)
      
      return
      end subroutine transition_time

      subroutine scramble_velocities(idnode,natms)
      
c***********************************************************************
c     
c     dlpoly hyperdynamics routine for randomising velocities after a 
c     transition has occured (use with identical species only)
c     
c     copyright daresbury laboratory
c     author      w.smith  jan  2007
c     
c**********************************************************************
            
      implicit none

      integer idnode,natms,i,j,k,m,n
      real(8) vvv

      do j=1,10
        
        do i=1,natms
          
          k=int(natms*duni())+1
          vvv=vxx(i)
          vxx(i)=vxx(k)
          vxx(k)=vvv
          m=int(natms*duni())+1
          vvv=vyy(i)
          vyy(i)=vyy(m)
          vyy(m)=vvv
          n=int(natms*duni())+1
          vvv=vzz(i)
          vzz(i)=vzz(n)
          vzz(n)=vvv
          
        enddo
        
      enddo
      
      return
      end subroutine scramble_velocities
      
      subroutine hyper_close(ltad,idnode,mxnode,natms,nsteql)
      
c***********************************************************************
c     
c     dlpoly routine for saving hyperdynamics restart data
c     
c     copyright daresbury laboratory
c     author      w.smith  dec  2007
c     
c**********************************************************************
      
      implicit none
      
      logical ltad
      integer idnode,mxnode,natms,nsteql
      integer iatm0,iatm1,i,j,k,n,last,ierr,netdif
      real(8) buff(2)
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     total number of atomic differences
      
      buff(1)=dble(ndiff)
      call gdsum(buff(1),1,buff(2))
      netdif=nint(buff(1))
      
      if(idnode.eq.0)then
        
c     open hyperdynamics restart file
      
        open(nhrs,file="HYPREV",form="unformatted")
      
c     write control variables
        
        write(nhrs)ltad,dble(numbsn),dble(numtrk),dble(numpro),
     x    dble(netdif),dble(numdark),dble(nsteql),dble(numbpd),
     x    timhyp,timres,tstop,tkeres,tboost,vbase,strres,celres
        write(nhrs)track
      endif
      
      if(ltad.and.numbsn.gt.1)then
        
c     load basin difference data
        
        n=0
        do i=1,numbsn-1
          
          do j=nbsa(i),nbsb(i)
            
            buffer(n+1)=dble(idabsn(j))
            buffer(n+2)=dble(i)
            buffer(n+3)=xdiffs(j)
            buffer(n+4)=ydiffs(j)
            buffer(n+5)=zdiffs(j)
            n=n+5
            
          enddo
          
        enddo
        last=n
        
c     write basin difference data
        
        do k=1,mxnode-1
          
          if(idnode.eq.0)then
            
            
            call csend(hyper_tag+k,buff,1,k,ierr)
            call crecv(2*hyper_tag+k,buff,1)
            call crecv(3*hyper_tag+k,buffer(last+1),nint(buff(1)))
            last=nint(buff(1))+last
            
          elseif(k.eq.idnode)then
            
            call crecv(hyper_tag+k,buff,1)
            buff(1)=dble(last)
            call csend(2*hyper_tag+k,buff,1,0,ierr)
            call csend(3*hyper_tag+k,buffer,last,0,ierr)
            
          endif
          
        enddo
        
        if(idnode.eq.0)write(nhrs)(buffer(i),i=1,last)
        call gsync()
        
      endif
      
c     load reference block configuration data
        
      j=0
      k=1
      do i=iatm0,iatm1
        
        buffer(j+1)=xres(k)
        buffer(j+2)=yres(k)
        buffer(j+3)=zres(k)
        j=j+3
        k=k+1
        
      enddo
      last=j
      
c     write reference block configuration data
      
      do k=1,mxnode-1
        
        if(idnode.eq.0)then
          
          call csend(hyper_tag+k,buff,1,k,ierr)
          call crecv(2*hyper_tag+k,buff,1)
          call crecv(3*hyper_tag+k,buffer(last+1),nint(buff(1)))
          last=nint(buff(1))+last
          
        elseif(k.eq.idnode)then
          
          call crecv(hyper_tag+k,buff,1)
          buff(1)=dble(last)
          call csend(2*hyper_tag+k,buff,1,0,ierr)
          call csend(3*hyper_tag+k,buffer,last,0,ierr)
          
        endif
        
      enddo
      
      if(idnode.eq.0)then
        
        write(nhrs)(buffer(i),i=1,last)
        close(nhrs)
        
      endif
      call gsync()
      
      return
      end subroutine hyper_close
      
      subroutine hyper_open(ltad,idnode,mxnode,natms,nsteql)
      
c***********************************************************************
c     
c     dlpoly routine for reading hyperdynamics restart data
c     
c     copyright daresbury laboratory
c     author      w.smith  dec  2007
c     
c**********************************************************************
      
      implicit none
      
      logical ltad,mtad
      integer idnode,mxnode,natms,nsteql
      integer iatm0,iatm1,i,j,k,n,last,netdif,ierr
      real(8) buff(1)
      
c     block indices
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
c     restore control variables
      
      mtad=.true.
      if(idnode.eq.0)then
        
c     open hyperdynamics restart file
      
        open(nhrs,file="HYPOLD",form="unformatted")
      
        read(nhrs)mtad,(buffer(i),i=1,30)
        read(nhrs)track
        
      else
        
        do i=1,30
          buffer(i)=0.d0
        enddo
        track(:)=0.d0
        
      endif
      
c     check restart file is tad compatible
      
      call gstate(mtad)
      if(ltad.and.(.not.mtad))call error(idnode,2341)

      call gdsum(buffer(1),31,buffer(32))
      
      numbsn=nint(buffer(1))
      numtrk=nint(buffer(2))
      numpro=nint(buffer(3))
      netdif=nint(buffer(4))
      numdark=nint(buffer(5))
      nsteql=nint(buffer(6))
      numbpd=nint(buffer(7))
      timhyp=buffer(8)
      timres=buffer(9)
      tstop=buffer(10)
      tkeres=buffer(11)
      tboost=buffer(12)
      vbase=buffer(13)
      do i=1,9
        
        strres(i)=buffer(i+13)
        celres(i)=buffer(i+22)

      enddo
      last=size(track)
      call gdsum(track,last,buffer)

      if(ltad.and.numbsn.gt.1)then
        
c     restore basin difference data
        
        last=5*netdif
        if(idnode.eq.0)read(nhrs)(buffer(i),i=1,last)
        
        do k=1,mxnode-1
          
          if(idnode.eq.0)then
            
            call csend(hyper_tag+k,buffer,last,k,ierr)
            
          elseif(k.eq.idnode)then
            
            call crecv(hyper_tag+k,buffer,last)
            
          endif
          
        enddo
        
c     reject nonlocal basin difference data
        
        j=0
        do i=1,last,5
          
          n=nint(buffer(i))
          if(n.ge.iatm0.and.n.le.iatm1)then
            
            buffer(j+1)=buffer(i)
            buffer(j+2)=buffer(i+1)
            buffer(j+3)=buffer(i+2)
            buffer(j+4)=buffer(i+3)
            buffer(j+5)=buffer(i+4)
            j=j+5
            
          endif
          
        enddo
        last=j
        
c     unload basin difference data
        
        n=0
        nbsa(1)=1
        do i=1,numbsn-1
          
          if(i.gt.1)nbsa(i)=n+1
          
          do j=1,last,5
            
            if(nint(buffer(j+1)).eq.i)then
              
              n=n+1
              idabsn(n)=nint(buffer(j))
              xdiffs(n)=buffer(j+2)
              ydiffs(n)=buffer(j+3)
              zdiffs(n)=buffer(j+4)
              
            endif
            
          enddo
          
          nbsb(i)=n
          
        enddo
        ndiff=n
        call gsync()
        
      endif
      
c     retrieve reference block configuration data

      last=3*natms
      if(idnode.eq.0)read(nhrs)(buffer(i),i=1,last)
      
      do k=1,mxnode-1
        
c     read reference block configuration data
        
        if(idnode.eq.0)then
          
          call csend(hyper_tag+k,buffer,last,k,ierr)
          
        elseif(k.eq.idnode)then
          
          call crecv(hyper_tag+k,buffer,last)
          
        endif
        
      enddo
      
c     unload reference block configuration data
        
      n=1
      j=3*(iatm0-1)
      do i=iatm0,iatm1
        
        xres(n)=buffer(j+1)
        yres(n)=buffer(j+2)
        zres(n)=buffer(j+3)
        j=j+3
        n=n+1
        
      enddo
      
      if(idnode.eq.0)close(nhrs)
      call gsync()
      
      return
      end subroutine hyper_open
      
      subroutine bpd_forces(natms,keybpd,vmin,ebias,temp,engcfg)
      
c***********************************************************************
c     
c     dl_poly subroutine for scaling forces in a bias potential dynamics 
c     simulation using hamelberg, mongan and mccammon factor
c     J. Chem. Phys. 120 (2004) 11919
c     
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2008
c     
c***********************************************************************
      
      integer i,natms,mynode,keybpd
      real(8) alpha,vmin,ebias,beta,temp,engcfg,eboost,hscale
      real(8) engtmp
      
      boost=1.d0
      numbpd=numbpd+1
      
c     reset potential energy wrt base level
      
      if(keybpd.eq.1)then
        engtmp=engcfg-vbase
      else
        engtmp=engcfg-engbsn(1)
      endif

      if(ebias.gt.engtmp)then
        
c     bias potental boost
        
        alpha=ebias*(ebias-vmin)/vmin
        beta=1.d0/(boltz*temp*dble(natms))
        eboost=(ebias-engtmp)**2/(alpha+ebias-engtmp)
        boost=exp(beta*eboost)
        
c     bias potential forces scaling factor
        
        hscale=(alpha/(alpha+ebias-engtmp))**2
        
c     scale forces
        
        do i=1,natms
          
          fxx(i)=fxx(i)*hscale
          fyy(i)=fyy(i)*hscale
          fzz(i)=fzz(i)*hscale
          
        enddo
        
      endif
      
c     accumulative average of boost factor
      
      tboost=boost/dble(numbpd)+dble(numbpd-1)*tboost/dble(numbpd)
      
      return
      end subroutine bpd_forces

      end module hyper_dynamics_module
