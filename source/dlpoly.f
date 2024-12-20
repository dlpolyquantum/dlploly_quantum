      program dlpoly_quantum
      
c***********************************************************************
c     
c     dl_poly classic is an stfc/ccp5 program package for the 
c     dynamical simulation of molecular systems.
c     
c     dl_poly is the copyright of the stfc daresbury laboratory, 
c     daresbury, warrington wa4 4ad. 
c     
c     neither the stfc, daresbury laboratory, ccp5 nor the authors
c     of this package claim that it is free from errors and do not
c     accept liability for any loss or damage that may arise from
c     its use. it is the users responsibility to verify that the 
c     package dl_poly is fit for the purpose the user intends for
c     it.
c     
c     users of this package are recommended to consult the dl_poly
c     user manual for the full description of its use and purpose.
c     
c     authors: w.smith and t.r.forester 1995
c     copyright daresbury laboratory 1995
c
c     M.R.Momeni & F.A.Shakib have added the following modules and
c     subroutines to the DL_POLY CLASSIC VERSION 1.10 along with
c     necessary changes in a variety of other modules. They do not 
c     claim that their modifications are free from errors and do 
c     not accept liability for any loss or damage that may arise 
c     from its use. it is the users responsibility to verify that 
c     the modifications are fit for the purpose the user intends for
c     it.
c
c     nhc_module.f
c     water_module.f
c     nptvv_nhc subroutine in vv_motion_module.f
c     nvtvv_nhc subroutine in vv_motion_module.f
c     nhc_part subroutine in ensemble_tools_module.f
c     nhc_baro subroutine in ensemble_tools_module.f
c     double-wall Lennard-Jones potential in external_field_module.f
c
c     The modifications in other modules are specified whenever 
c     occured.
c
c     authors: M.R.Momeni and F.A.Shakib
c     copyright M.R.Momeni and F.A.Shakib 2021
c     
c     Method Development and Materials Simulation Laboratory
c
c                         DL_POLY QUANTUM VERSION 1.0
c
c***********************************************************************
      
c     declare required modules
      
      use angles_module
      use bonds_module
      use config_module
      use core_shell_module
      use define_system_module
      use dihedral_module
      use driver_module
      use ewald_module
      use exclude_module
      use external_field_module
      use forces_module
      use four_body_module
      use hkewald_module
      use hyper_dynamics_module
      use integrator_module
      use inversion_module
      use metal_module
      use metafreeze_module
      use nlist_builders_module
      use pair_module
      use pimd_module
      use pmf_module
      use property_module
      use rigid_body_module
      use setup_module
      use shake_module
      use site_module
      use solvation_module
      use spme_module
      use temp_scalers_module
      use tersoff_module
      use tether_module
      use three_body_module
      use utility_module
      use vdw_module
      use vv_pimd_module
      use nhc_module
      use water_module
      use correlation_module
      use pimd_piglet_module
      use fqcmd_module

      implicit none
      
      character*1 hms,dec
      character*8 seek
      logical safe

      logical ltscal,lzeql,loptim,ltraj,lgofr,lpgr,lfcap,recycle
      logical newlst,lneut,loglnk,lnsq,lzden,lshmov,lcnb,ltad,lneb
      logical stropt,lzero,nolink,newgau,lminim,lminnow,lhit,lbpd
      logical prechk,tadall,lexcite,lsolva,lfree,lfrmas,lswitch
      logical lghost,llswitch,lnfic,nebgo,lpsoc,redirect,lpimd
      logical inhc,lmsite,lcorr,lfqcmd,lfcmd
      
      integer npage,lines,idnode,mxnode,memr,intsta,istraj,nsbzdn
      integer keyens,keyfce,keyres,keytrj,kmax1,kmax2,kmax3,multt
      integer nstack,nstbgr,nstbpo,nhko,nlatt,nstbts,nsteql,nstraj
      integer nstrun,nospl,keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw
      integer ntptbp,ntpmet,ntpfbp,nshels,imcon,levcfg,nneut,minstp
      integer ntangl,ntbond,ntcons,ntdihd,ntinv,ntpmf,nspmf,ntfree
      integer ntteth,ntshl,nstep,numacc,numrdf,nzden,nscons,i,k
      integer ntpter,keyshl,isw,keyver,keystr,keytol,numgau,khit
      integer nhit,keybpd,ntrack,nblock,blkout,numneb,nturn,mode
      integer natms2,ntghost,nsolva,isolva,nofic,iadd
      integer nrespa,iqt4,keycorr,molcorr,wrtcorr,wrtrdf

      real(8) alpha,delr,epsq,fmax,press,quattol,rcut,rprim,rvdw,taup
      real(8) taut,temp,timcls,timjob,tolnce,tstep,tzero,dlrpot,drewd
      real(8) engunit,rcuttb,rctter,rcutfb,degfre,degrot,chit,conint
      real(8) elrc,virlrc,engbnd,volm,degshl,chip,virbnd,engang,virang
      real(8) engdih,virdih,enginv,virinv,engtbp,virtbp,engter,virter
      real(8) engfbp,virfbp,engsrp,virsrp,engcpe,vircpe,vircon,vircom
      real(8) engfld,virfld,engshl,virshl,shlke,engtet,virtet,virpmf
      real(8) consv,engke,engrot,sigma,virtot,engcfg,prntim,simtim
      real(8) sigma_nhc,sigma_volm,alpha_volm,g_qt4f
      real(8) stpeng,stpeth,stpprs,stptmp,stpvir,stpvol,width,zlen
      real(8) timelp,engmet,virmet,pass0,pass1,pass2,rlxtol,opttol
      real(8) catchrad,sprneb,deltad,tlow,engtke,ehit,xhit,yhit,zhit
      real(8) ebias,vmin,hyp_units,estar,chit_shl,sigma_shl,engthe,chi
      real(8) engord,virord,engqpi,engqvr,engrng,virrng,qmsrgr,qmsbnd
      real(8) uuu(102),gaumom(0:5)
      real(8), allocatable :: tbuffer(:)
      
      data timelp/0.d0/,lminnow/.false./,ntrack/10/
      data npage,lines/8,0/,recycle/.true./
      data pass0/0.d0/,pass1/0.d0/,pass2/0.d0/
      data delr,epsq,press,quattol,rprim,rvdw/6*0.d0/
      data temp,timcls,timjob,tolnce,rlxtol/5*0.d0/
      data gaumom/6*0.d0/
      
c     set up the communications
      
      call initcomms()
      call gsync()
      
c     determine processor identities
      
      call machine(idnode,mxnode)
      
c     activate for limited-life executable
      
c$$$      call bomb(idnode,2020,6,26)
      
      allocate (tbuffer(10),stat=memr)
      
      call parset(redirect,idnode,mxnode,tbuffer)
      
c     open main printing file
      
      if(.not.redirect.and.idnode.eq.0)open(nrite,file='OUTPUT')
      if(idnode.eq.0) write (nrite,
     x  "(/,20x,'DL_POLY Quantum 2.0',
     x   /,/,20x,'Running on ',i4,' nodes',/,/)") mxnode

c     allocate arrays for each function
      
      call alloc_ang_arrays(idnode,mxnode)
      call alloc_bnd_arrays(idnode,mxnode)
      call alloc_config_arrays(idnode,mxnode)
      call alloc_csh_arrays(idnode,mxnode)
      call alloc_dih_arrays(idnode,mxnode)
      call alloc_ewald_arrays(idnode,mxnode)
      call alloc_exc_arrays(idnode,mxnode)
      call alloc_exi_arrays(idnode,mxnode)
      call alloc_fbp_arrays(idnode,mxnode)
      call alloc_fld_arrays(idnode,mxnode)
      call alloc_free_arrays(idnode,mxnode)
      call alloc_hke_arrays(idnode,mxnode)
      call alloc_hyper_arrays(idnode,mxnode)
      call alloc_inv_arrays(idnode,mxnode)
      call alloc_met_arrays(idnode,mxnode)
      call alloc_pair_arrays(idnode,mxnode)
      call alloc_pimd_arrays(idnode,mxnode)
      call alloc_pmf_arrays(idnode,mxnode)
      call alloc_prp_arrays(idnode,mxnode)
      call alloc_rgbdy_arrays(idnode,mxnode)
      call alloc_shake_arrays(idnode,mxnode)
      call alloc_site_arrays(idnode,mxnode)
      call alloc_sol_arrays(idnode,mxnode)
      call alloc_spme_arrays(idnode,mxnode)
      call alloc_tbp_arrays(idnode,mxnode)
      call alloc_ter_arrays(idnode,mxnode)
      call alloc_tet_arrays(idnode,mxnode)
      call alloc_vdw_arrays(idnode,mxnode)
      
c     start clock
      
      call timchk(0,tzero)
      
c     input the control parameters defining the simulation
      
      call simdef
     x  (seek,lfcap,lgofr,lnsq,loptim,lzero,lminim,lpgr,ltraj,ltscal,
     x  lzeql,lzden,nolink,newgau,lhit,lbpd,ltad,lneb,prechk,tadall,
     x  lsolva,lfree,lfrmas,lexcite,lswitch,lghost,lnfic,nebgo,lpsoc,
     x  lpimd,inhc,lmsite,lcorr,lfqcmd,lfcmd,idnode,minstp,intsta,
     x  istraj,keybpd,keyens,keyfce,keyres,keyver,keytrj,keycorr,
     x  molcorr,wrtcorr,kmax1,kmax2,
     x  kmax3,multt,nstack,nstbgr,nsbzdn,nstbpo,nhko,nlatt,nstbts,
     x  nsteql,nstraj,nstrun,nospl,keytol,numgau,khit,nhit,nblock,
     x  ntrack,blkout,numneb,mode,nsolva,isolva,nofic,nbeads,nchain,
     x  nrespa,g_qt4f,alpha,delr,epsq,fmax,press,quattol,rcut,rprim,
     x  rvdw,taup,taut,temp,timcls,timjob,tolnce,tstep,rlxtol,opttol,
     x  zlen,ehit,xhit,yhit,zhit,ebias,vmin,catchrad,sprneb,deltad,tlow,
     x  hyp_units,chi,nsp1,wrtrdf)

c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

c     allocate NHC thermostat arrays

      call alloc_nhc_arrays(idnode,mxnode,nchain)

c     allocate arrays related to qtip4p/f water model

      call alloc_water_arrays(idnode,mxnode)
 
c *******************************************************************

c     allocate arrays related to piglet thermostat

      call alloc_piglet_arrays(idnode,mxnode)

c     input the system force field
      
      call sysdef
     x  (lneut,lnsq,lsolva,lfree,lexcite,lswitch,lghost,lpimd,idnode,
     x  keyfce,keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw,ntptbp,ntpmet,
     x  ntpfbp,ntpter,nshels,keyshl,ntghost,keyver,dlrpot,engunit,
     x  rvdw,rcuttb,rctter,rcutfb)
      
      if(ntpmet.gt.0.and.multt.gt.1)call error(idnode,153)
      
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

c     give indices to each M-site of each water molecule
c     if water model qtip4p/f is requested

      if(lmsite)then
             
         call water_index(idnode,mxnode,nbeads,ntpmls,iqt4)
      endif   

c *******************************************************************      

c     construct initial configuration of system
      
      call sysgen
     x  (loglnk,lneut,nolink,lfree,lfrmas,lpimd,lmsite,idnode,imcon,
     x  keyens,keyfce,keyres,levcfg,multt,mxnode,ntpmls,nbeads,iqt4,
     x  g_qt4f,delr,rcut,temp,volm,uuu)
      
c     determine system dimensions
      
      call dcell(cell,celprp)
      width=min(celprp(7),celprp(8),celprp(9))
      
c     construct initial bookkeeping arrays
   
      call sysbook
     x  (loglnk,lneut,lshmov,lcnb,lsolva,lghost,lmsite,idnode,
     x  imcon,mxnode,natms,nneut,ngrp,nscons,ntangl,ntbond,ntcons,
     x  ntdihd,ntinv,ntpmls,ntpmf,nspmf,ntfree,ntteth,ntshl,
     x  ntghost,degfre,degrot)
      
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

c     amend number of free atoms for pimd system and remove
c     M-sites from number of DOFs in classical and pimd 
c     simulations

      if(lmsite)then
         ntfree=nbeads*(ntfree-nummols(ntpmls))
      else
         ntfree=nbeads*ntfree
      endif   

c *******************************************************************      
      
c     reset atom numbers for excitation simulation
      
      if(lghost)then
        natms2=natms-ntghost
      else
        natms2=natms
      endif
      
c     set initial system temperature
      
      call systemp
     x  (lpimd,inhc,idnode,imcon,keyres,mxnode,natms2,nbeads,ngrp,
     x  nscons,ntcons,ntfree,ntshl,levcfg,keyshl,keyens,degfre,degshl,
     x  nchain,degrot,engke,tolnce,temp,sigma,sigma_nhc,sigma_volm,
     x  alpha_volm,uuu,tstep)
      
c     read thermodynamic and structural data from restart file
      
      call sysinit
     x  (lgofr,lzden,lsolva,lfree,lghost,lpsoc,lpimd,idnode,imcon,
     x  keyfce,keyres,mxnode,natms,nbeads,ntshl,nstep,numacc,numrdf,
     x  ntpatm,ntpmet,ntpvdw,nzden,chip,chit,conint,elrc,engunit,
     x  virlrc,rvdw,volm,virtot,vircom,tboost,chit_shl,engthe,gaumom)
      
c     metadynamics by d. quigley
      
      if(lmetadyn) then
        
c     make copy of excluded atom list for use by metadynamics
        
        call exclude_copy_mtd(idnode)
        
c     initialise metadynamics, read order parameter definitions etc.
        
        call define_metadynamics(idnode,mxnode,natms,ntpatm,temp)   

      end if

c     synchronise LRC, SIC and system charge terms for switching
      
      llswitch=.false.
      if(lswitch)then
        
        if(nstep.ge.nswitch)then
          
          if(mod((nstep-nswitch)/niswitch,2).eq.0)then
            
            call switch_atm(lfrmas)
            call switch(elrc,virlrc)
            llswitch=.true.
            
          endif
          
        endif
        
      endif
      
c     zero long range component of stress
      
      stresl(:)=0.d0
      
c     zero constraint terms
      
      vircon=0.d0
      virpmf=0.d0
      if(lminim.or.loptim.or.ntcons.eq.0)strcns(:)=0.d0

c     metadynamics by d. quigley
      
      sigma_shl=boltz*degshl*0.5d0      
      
c     convert BPD parameters to internal units
      
      if(lbpd)then
        
        ebias=0.5d0*boltz*degfre*ebias
        vmin=0.5d0*boltz*degfre*vmin
        
      endif


c     time check

      call timchk(1,tzero)
      
c      numPot=3
c      numPoints=847
c      write(*,*) fxx(:)
c      write(*,*) "engbnd ", engbnd
c      call allocate_fqcmd_arrays(idnode,mxnode)
c      call read_pot_tables(idnode,mxnode)
c      write(*,*) "fqcmd test start"
c      nqcbnd=1
c      nqcang=1
      if(lfqcmd) then
        call fqcmd_setup(idnode,mxnode,natms,
     x  ntbond,ntangl,ntpatm,ntpmls,lpimd)
      endif
c      if(lfqcmd.and.nstep.gt.nsteql) then
c      if(lfqcmd) then
c        call  get_qcent_internal(idnode,mxnode,imcon,natms,
c     x  ntbond,ntangl,ntpmls,rcut)
c        call get_qcent_cart(idnode,mxnode,imcon,
c     x  ntpmls,natms,ntbond,rcut,lpimd)
c      endif
c     control variable for structure optimizer
      
      keystr=0
      stropt=.false.
      
      if(lminim)then
        
c     first step of minimisation programme

        if(idnode.eq.0)write(nrite,"(1x,120('-'))")
        
        call minimiser
     x    (lfcap,lneut,lnsq,loglnk,lzeql,newlst,idnode,imcon,keyfce,
     x    keyfld,keyshl,keytol,kmax1,kmax2,kmax3,multt,mxnode,natms,
     x    ngrp,nhko,nlatt,nneut,nospl,nscons,ntcons,nstbgr,nstep,
     x    nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet,
     x    ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf,alpha,delr,dlrpot,
     x    drewd,elrc,engang,engbnd,engcpe,engdih,engfbp,engfld,enginv,
     x    engshl,engsrp,engtbp,engter,engtet,epsq,fmax,opttol,rctter,
     x    rcut,rcutfb,rcuttb,rprim,rvdw,shlke,engcfg,temp,tstep,
     x    virang,virbnd,vircpe,virdih,virfbp,virfld,virinv,virlrc,
     x    virmet,virshl,virsrp,virtbp,virter,virtet,volm,engmet,
     x    virtot,sigma,tolnce,engunit,engord,virord,lfqcmd)
        
      elseif((lpimd.or.keyver.gt.0).and.nstep.eq.0)then
        
c     calculate initial conditions for velocity verlet
        
        call molecular_dynamics
     x    (lfcap,lgofr,lneut,lnsq,loglnk,loptim,lzeql,lzero,newlst,
     x    stropt,recycle,ltad,lsolva,lfree,lghost,lpimd,idnode,imcon,
     x    keyfce,keyfld,keyshl,keystr,keytol,kmax1,kmax2,kmax3,multt,
     x    mxnode,natms,ngrp,nhko,nlatt,nneut,nospl,nscons,nstbgr,nstep,
     x    nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet,
     x    ntptbp,ntpter,ntpvdw,ntshl,ntteth,ntcons,numrdf,nsolva,
     x    isolva,nbeads,alpha,delr,dlrpot,drewd,elrc,engang,engbnd,
     x    engcpe,engdih,engfbp,engfld,enginv,engshl,engsrp,engtbp,
     x    engter,engtet,epsq,fmax,opttol,rctter,rcut,rcutfb,rcuttb,
     x    rprim,rvdw,shlke,engcfg,temp,tstep,virang,virbnd,vircpe,
     x    virdih,virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
     x    virtbp,virter,virtet,volm,engmet,virtot,engord,virord,
     x    engrng,virrng,qmsbnd,keyens,lfqcmd)
       
c     bias potential dynamics option - reset forces
        
        if(lbpd)call bpd_forces(natms,keybpd,vmin,ebias,temp,engcfg)
        
      endif
      
c     stage initial forces for pimd
      
      if(lpimd)then

        if(keyens.le.42)then

          call stage_forces(lmsite,idnode,
     x               mxnode,natms,nbeads, ntpmls,g_qt4f)

        elseif(keyens.ge.43)then

          call force2norm(lmsite,idnode,
     x               mxnode,natms,nbeads,ntpmls,g_qt4f)
        endif

      elseif(lmsite)then

        call qt4_force_redist(idnode,mxnode,nbeads,ntpmls,g_qt4f) 

      endif
      
c      numPot=3
c      numPoints=847
c      write(*,*) fxx(:)
c      write(*,*) "engbnd ", engbnd
c      call allocate_fqcmd_arrays(idnode,mxnode)
c      call read_pot_tables(idnode,mxnode)
c      write(*,*) "fqcmd test start"
c      call fqcmd_setup(idnode,mxnode,ntpatm)
c      call test_interp(imcon,natms)
c      call fqcmd_correct_force(idnode,mxnode,imcon,natms,
c     x  engsrp,engcfg)
c      write(*,*) fxx(:)
c      write(*,*) "xxx ",xxx(:)
c      write(*,*) "engbnd ", engbnd

      if(ltad.or.(lbpd.and.keybpd.eq.2))then
        
c     construct the first reference state
        
        call hyper_start
     x    (ltad,lbpd,lfcap,lneut,lnsq,loglnk,lzeql,newlst,nblock,idnode,
     x    imcon,keyfce,keyfld,keyshl,keytol,kmax1,kmax2,kmax3,multt,
     x    mxnode,natms,ngrp,nhko,nlatt,nneut,nospl,nscons,nstbgr,
     x    nstep,nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,
     x    ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,ntcons,ntrack,alpha,
     x    delr,dlrpot,drewd,elrc,virlrc,epsq,fmax,opttol,rctter,
     x    rcut,rcutfb,rcuttb,rprim,rvdw,temp,tstep,volm,sigma,
     x    hyp_units)
        
      endif
      
c     perform selected NEB calculation
      
      if(lneb)then
        
        do i=1,numneb
          
          call neb_driver
     x      (lfcap,lneut,lnsq,loglnk,lzeql,newlst,lneb,bsn_1(i),
     x      bsn_2(i),idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,
     x      keytol,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,ngrp,
     x      ntcons,ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,
     x      keyshl,ntfree,keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,
     x      ntshl,nscons,delr,dlrpot,engcpe,engsrp,epsq,rcut,
     x      rprim,rvdw,vircpe,virsrp,alpha,drewd,volm,
     x      engmet,virmet,elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,
     x      engfbp,virfbp,rctter,engter,virter,engbnd,virbnd,
     x      engang,virang,engdih,virdih,enginv,virinv,engtet,
     x      virtet,engshl,shlke,virshl,engfld,virfld,engcfg,fmax,
     x      temp,tstep,opttol,sprneb,hyp_units)
          
          call scan_profile(nturn,estar)
          
          if(idnode.eq.0)then
            
            write(nrite,"(1x,120('-'))")
            write(nrite,'(1x,"TRA",3i6,1p,4e14.5)')
     x        bsn_1(i),bsn_2(i),nturn,estar/hyp_units
            write(nrite,"(1x,120('-'))")
            
          endif
          
        enddo
        
c     bypass the MD cycle for this option
        
        recycle=.false.
        
      endif
      
      if(lcorr)then
        call corr_init(idnode,mxnode,natms,ntpmls,molcorr,keyens,
     x    keycorr,nummols,numsit)
      endif
      
c***********************************************************************
c     start of molecular dynamics calculations
c***********************************************************************
      
      do while(recycle)
        
c     increase step counter
        
        nstep=nstep+1
        recycle=(nstep.lt.nstrun)
        
c     store velocities for free energy or solvation simulation
        
        if(keyver.eq.0)then
          
          if(lsolva)then
            
            vxo_sol(:)=vxx(:)
            vyo_sol(:)=vyy(:)
            vzo_sol(:)=vzz(:)
            
          elseif(lfree)then
            
            vxo_fre(:)=vxx(:)
            vyo_fre(:)=vyy(:)
            vzo_fre(:)=vzz(:)
            
          endif
          
        endif
        
c     molecular switching option for excitation
        
        if(lswitch)then
          
          if(nstep.ge.nswitch)then
            
            if(mod(nstep-nswitch,niswitch).eq.0)then
              
              call switch_atm(lfrmas)
              call switch(elrc,virlrc)
              llswitch=.not.llswitch
              
            endif
            
          endif
          
        endif
        
c     switch on the minimiser
        
        if(lminim)then
          
          lminnow=(mod(nstep,minstp).eq.0)
          
        endif
        
c     conserved quantity (other than K + U)
        
        consv=0.d0
        
c     energy accumulators
        
        if(.not.lminnow)then
          
          engke=0.d0
          engrot=0.d0
          
        endif

c     calculate volume of simulation cell
        
        if(imcon.ne.0.and.imcon.ne.6)then
          
          call dcell(cell,celprp)
          volm=celprp(10)
          if(imcon.eq.4)then
            
            volm=0.5d0*celprp(10)
            
          elseif(imcon.eq.5)then
            
            volm=0.5d0*celprp(10)
            
          elseif(imcon.eq.7)then
            
            volm=0.5d0*celprp(10)
            
          endif
          
        else
          
          volm=0.d0
          
        endif
        
c     reset sutton chen long range corrections (constant pressure only)
        
        if(ntpmet.gt.0.and.keyens.ge.4.and.keyens.le.7) then
          
          call lrcmetal
     x      (idnode,imcon,natms,ntpatm,engunit,rvdw,volm)
          
        endif
        
c     activate the impact option at designated time step
        
        if(lhit.and.nstep.eq.nhit)call impact
     x    (khit,natms,idnode,mxnode,ehit,xhit,yhit,zhit)
        
c     integrate equations of motion stage 1 of velocity verlet
        
        if(lpimd)then
          
          isw=1
c          call pimd_integrate
c     x      (lmsite,isw,idnode,mxnode,imcon,ntpmls,natms,keyens,nstep,
c     x      tstep,taut,g_qt4f,temp,engke,engthe,chi,uuu,gaumom)
          
          call pimd_integrate
     x      (lmsite,isw,idnode,mxnode,imcon,ntpmls,natms,keyens,nstep,
     x      tstep,g_qt4f,temp,engke,engthe,chi,uuu,gaumom,
     x      safe,nrespa,ntpatm,ntshl,keyshl,taut,taup,sigma,
     x      sigma_nhc,sigma_volm,alpha_volm,virtot,vircon,virlrc,virrng,
     x      press,volm,chit,consv,conint,elrc,chit_shl,sigma_shl)
      
        elseif(keyver.gt.0)then
          
          isw=1
          if(.not.loptim)then
            
            if(llswitch)call copy_force(idnode,mxnode)
            
            call vv_integrate
     x        (lcnb,lshmov,lnfic,lmsite,isw,idnode,mxnode,imcon,natms2,
     x        nstep,nchain,nrespa,ntpmls,ngrp,keyens,nscons,ntcons,
     x        ntpatm,ntfree,nspmf,ntpmf,mode,nofic,ntshl,keyshl,tstep,
     x        engke,engrot,tolnce,vircon,vircom,virtot,temp,press,volm,
     x        sigma,sigma_nhc,sigma_volm,alpha_volm,taut,taup,
     x        chit,chip,consv,conint,g_qt4f,elrc,virlrc,virpmf,
     x        chit_shl,sigma_shl,gaumom)
            
            if(lghost)call update_ghost(idnode,mxnode)
            
            if(keyens.ge.4.and.keyens.le.7)then
              
              if(lfree.or.lghost)
     x          call lrcorrect_fre(lfree,volm,elrc,virlrc)
              if(lsolva)call lrcorrect_sol(lghost,volm)
              
            endif
            
          endif

c     scale t=0 tether reference positions (constant pressure only)
          
          if(keyens.ge.4.and.keyens.le.7) then
            
            call xscale(idnode,mxnode,natms,keyens,imcon,tstep)
            
          endif
          
        endif
        
        if(lminnow)then
          
          call minimiser
     x      (lfcap,lneut,lnsq,loglnk,lzeql,newlst,idnode,imcon,keyfce,
     x      keyfld,keyshl,keytol,kmax1,kmax2,kmax3,multt,mxnode,natms,
     x      ngrp,nhko,nlatt,nneut,nospl,nscons,ntcons,nstbgr,nstep,
     x      nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet,
     x      ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf,alpha,delr,dlrpot,
     x      drewd,elrc,engang,engbnd,engcpe,engdih,engfbp,engfld,enginv,
     x      engshl,engsrp,engtbp,engter,engtet,epsq,fmax,opttol,rctter,
     x      rcut,rcutfb,rcuttb,rprim,rvdw,shlke,engcfg,temp,tstep,
     x      virang,virbnd,vircpe,virdih,virfbp,virfld,virinv,virlrc,
     x      virmet,virshl,virsrp,virtbp,virter,virtet,volm,engmet,
     x      virtot,sigma,tolnce,engunit,engord,virord,lfqcmd)
          
        elseif(loptim.or.keyshl.ne.2)then
          
          call molecular_dynamics
     x      (lfcap,lgofr,lneut,lnsq,loglnk,loptim,lzeql,lzero,
     x      newlst,stropt,recycle,ltad,lsolva,lfree,lghost,lpimd,
     x      idnode,imcon,keyfce,keyfld,keyshl,keystr,keytol,kmax1,
     x      kmax2,kmax3,multt,mxnode,natms,ngrp,nhko,nlatt,nneut,
     x      nospl,nscons,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,
     x      ntfree,ntinv,ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,
     x      ntteth,ntcons,numrdf,nsolva,isolva,nbeads,alpha,delr,
     x      dlrpot,drewd,elrc,engang,engbnd,engcpe,engdih,engfbp,
     x      engfld,enginv,engshl,engsrp,engtbp,engter,engtet,epsq,
     x      fmax,opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,shlke,
     x      engcfg,temp,tstep,virang,virbnd,vircpe,virdih,
     x      virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
     x      virtbp,virter,virtet,volm,engmet,virtot,engord,virord,
     x      engrng,virrng,qmsbnd,keyens,lfqcmd)
          
        else
          
          call shell_relaxation
     x      (lfcap,lgofr,lneut,lnsq,loglnk,lzeql,newlst,ltad,lsolva,
     x      lfree,lghost,idnode,imcon,keyfce,keyfld,keyshl,
     x      kmax1,kmax2,kmax3,multt,mxnode,natms,nhko,nlatt,nneut,
     x      nospl,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,ntinv,
     x      ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf,
     x      ntpmls,nsolva,isolva,alpha,delr,dlrpot,drewd,elrc,engang,
     x      engbnd,engcpe,engdih,engfbp,engfld,enginv,engshl,engsrp,
     x      engtbp,engter,engtet,epsq,fmax,rctter,rcut,rcutfb,rcuttb,
     x      rprim,rvdw,shlke,engcfg,temp,tstep,virang,virbnd,vircpe,
     x      virdih,virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
     x      virtbp,virter,virtet,volm,engmet,virtot,rlxtol,pass0,
     x      pass1,pass2,engord,virord,lfqcmd)
          
        endif
        
c     stage forces for pimd
        if(lpimd.and.(keyens.le.42))call stage_forces(lmsite,idnode,
     x               mxnode,natms,nbeads, ntpmls,g_qt4f)
        if(lpimd.and.(keyens.ge.43))call force2norm(lmsite,idnode,
     x               mxnode,natms,nbeads,ntpmls,g_qt4f)
c      call fqcmd_correct_force(idnode,mxnode,imcon,natms,
c     x  engsrp,engcfg)
c      write(*,*) "xxx ",xxx(:)
c      write(*,*) "engbnd ", engbnd
        
c     bias potential dynamics option - reset forces
        
        if(lbpd)call bpd_forces(natms,keybpd,vmin,ebias,temp,engcfg)
        
c     switching option for excitation simulation
        
        if(llswitch)call copy_force(idnode,mxnode)
        
c     integrate equations of motion
        
        if(lpimd)then
          
          isw=2
c          call pimd_integrate
c     x      (lmsite,isw,idnode,mxnode,imcon,ntpmls,natms,keyens,nstep,
c     x      tstep,taut,g_qt4f,temp,engke,engthe,chi,uuu,gaumom)
          
          call pimd_integrate
     x      (lmsite,isw,idnode,mxnode,imcon,ntpmls,natms,keyens,nstep,
     x      tstep,g_qt4f,temp,engke,engthe,chi,uuu,gaumom,
     x      safe,nrespa,ntpatm,ntshl,keyshl,taut,taup,sigma,
     x      sigma_nhc,sigma_volm,alpha_volm,virtot,vircon,virlrc,virrng,
     x      press,volm,chit,consv,conint,elrc,chit_shl,sigma_shl)
      
        elseif(keyver.eq.0)then
          
c     integrate equations of motion by leapfrog verlet
          
          if(.not.(loptim.or.lminnow))call lf_integrate
     x      (lcnb,lshmov,lnfic,idnode,mxnode,imcon,natms2,nstep,ngrp,
     x      keyens,nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nofic,
     x      tstep,engke,engrot,tolnce,quattol,vircon,vircom,virtot,
     x      temp,press,volm,sigma,taut,taup,chit,chip,consv,conint,
     x      elrc,virlrc,virpmf,gaumom)
          
        else if(keyver.gt.0)then
          
c     integrate equations of motion by velocity verlet (stage 2)
          
          isw=2
          if(.not.loptim)call vv_integrate
     x      (lcnb,lshmov,lnfic,lmsite,isw,idnode,mxnode,imcon,natms2,
     x      nstep,nchain,nrespa,ntpmls,ngrp,keyens,nscons,ntcons,
     x      ntpatm,ntfree,nspmf,ntpmf,mode,nofic,ntshl,keyshl,tstep,
     x      engke,engrot,tolnce,vircon,vircom,virtot,temp,press,volm,
     x      sigma,sigma_nhc,sigma_volm,alpha_volm,taut,taup,
     x      chit,chip,consv,conint,g_qt4f,elrc,virlrc,virpmf,
     x      chit_shl,sigma_shl,gaumom)
          
        endif
        
c     update the atomic positions for the ghost molecule
        
        if(lghost)call update_ghost(idnode,mxnode)
        
c     long range correction adjustment for free energy and solvation
        
        if(keyens.ge.4.and.keyens.le.7)then
          
          if(lfree.or.lghost)call lrcorrect_fre(lfree,volm,elrc,virlrc)
          if(lsolva)call lrcorrect_sol(lghost,volm)
          
        endif
        
c     application of transition analysis procedures
        
        if(ltad.or.(lbpd.and.keybpd.eq.2))then
          
          engtke=engke+engrot
          call hyper_driver
     x      (seek,ltad,lbpd,recycle,lfcap,lneut,lnsq,loglnk,lzeql,
     x      newlst,prechk,tadall,nebgo,nblock,ntrack,idnode,imcon,
     x      keyfce,keyfld,keyshl,keytol,kmax1,kmax2,kmax3,multt,
     x      mxnode,natms,ngrp,ntcons,nhko,nlatt,nneut,nospl,nscons,
     x      nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,
     x      ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,blkout,
     x      alpha,delr,dlrpot,drewd,elrc,virlrc,epsq,fmax,
     x      opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,temp,
     x      tstep,volm,engcfg,catchrad,sprneb,deltad,tlow,engtke,
     x      tolnce,hyp_units,ebias,vmin)
          
        endif
        
c     reset average boost factor in BPD during equilibration
        
        if(lbpd.and.keybpd.eq.1)then
          
          if(lzeql.and.nstep.le.nsteql)then
            
            numbpd=0
            tboost=0.d0
            
          endif
          
        endif
        
c     calculate shell kinetic energy
        
        if(keyshl.eq.1)then
          
          call corshl(idnode,mxnode,ntshl,shlke)
          
        endif
        
c     scale t=0 tether reference positions (constant pressure only)
        
        if(keyver.eq.0.and.keyens.ge.4.and.keyens.le.7) then
          
          call xscale(idnode,mxnode,natms,keyens,imcon,tstep)
          
        endif
        
c     apply temperature scaling
        
        if((ltscal.and.nstep.le.nsteql).and.
     x    mod(nstep-nsteql,nstbts).eq.0)then
          
          chit=0.d0
          chit_shl=0.d0
          chip=0.d0
          eta(:)=0.d0

          if(keyshl.eq.1) then
            
            do k=1,4
              call vscaleg(idnode,mxnode,imcon,natms2,ngrp,sigma)
              call shlqnch(idnode,mxnode,ntshl,temp)
            enddo
            
          else
            
            call vscaleg(idnode,mxnode,imcon,natms2*nbeads,ngrp,sigma)
            
          endif
          
        endif
        
c     reset atom velocities at intervals if required
        
        if(newgau)then
          if(mod(nstep,numgau).eq.0)then
            if(lpimd) then 
              call gauss(natms*nbeads,vxx,vyy,vzz)
              call momenta2norm(idnode,mxnode,natms)
              call reset_pimd_momenta(idnode,mxnode,natms,sigma)
            else
            call regauss
     x      (idnode,imcon,mxnode,natms2*nbeads,ngrp,nscons,ntcons,
     x      ntshl,keyshl,sigma,temp,tolnce)
            endif
          endif
        endif
        
c     calculate quantum energy
        
        if(lpimd)then

          if(keyens.ge.43)then 

            call quantum_energy_nm
     x        (idnode,mxnode,natms,temp,engke,engcfg,engrng,engqpi,
     x        engqvr,qmsrgr)

          else

            call quantum_energy
     x        (idnode,mxnode,natms,temp,engke,engcfg,engrng,engqpi,
     x        engqvr,qmsrgr)

            call ring_energy
     x        (idnode,mxnode,natms,temp,engrng,virrng,qmsbnd,stress)

          endif

          engcfg=engcfg+engrng
        endif

c     calculate correlation function        
        if(lcorr.and.(mod(nstep,wrtcorr).eq.0))then
          call correlation
     x      (idnode,mxnode,natms,ntpmls,molcorr,keyens,keycorr,nstep,
     x      nummols,numsit,tstep) 
        endif
      
      if((lfqcmd.and.nstep.gt.nsteql)
     x  .and.(mod(nstep,wrtrdf).eq.0)) then
        call  get_qcent_internal(idnode,mxnode,imcon,natms,
     x  ntbond,ntangl,ntpmls,rcut,lfcmd)
        call get_qcent_cart(idnode,mxnode,imcon,
     x  ntpmls,natms,ntbond,rcut,lpimd,lfcmd)
      endif

c     calculate physical quantities
        
        if(nstep.gt.0)call static
     x    (lbpd,lzeql,lpimd,idnode,intsta,imcon,keyens,natms,nstack,
     x    nstep,nsteql,ntpatm,numacc,mxnode,nblock,keybpd,numbpd,
     x    consv,degfre,degrot,engang,engbnd,engcpe,engdih,enginv,
     x    engke,engrot,engsrp,engunit,engcfg,stpeng,stpeth,stpprs,
     x    stptmp,stpvir,stpvol,tstep,virbnd,engfbp,vircom,vircon,
     x    vircpe,virsrp,engfld,virfld,engtbp,virtbp,virpmf,virshl,
     x    engshl,engtet,virtet,degshl,shlke,virang,width,engmet,
     x    virmet,engter,virter,boost,tboost,engqpi,engqvr,engrng,
     x    virrng,qmsrgr,qmsbnd,engthe)

c     reset moment accumulators in equilbration period
        
        if(nstep.eq.nsteql)gaumom(0)=0.d0
        
c     z density calculation
        
        if(lzden.and.((.not.lzeql).or.(nstep.gt.nsteql))) then
          
          call zden0(idnode,natms,mxnode,nzden,zlen)
          
        endif
        
c     terminate program if boundary conditions violated
        
        if(imcon.gt.0.and.rcut.gt.width)then
          
          levcfg=2
          call revive
     x      (lgofr,lzden,lpimd,idnode,imcon,mxnode,natms,nbeads,levcfg,
     x      nstep,nzden,numacc,numrdf,chip,chit,conint,tstep,engcfg,
     x      virtot,vircom,tboost,chit_shl,gaumom)
          
          if(lpimd)then
            
            call write_thermostats(idnode,mxnode,natms,temp)
            if(keyens.eq.41)call save_rnd_cfg(idnode,mxnode,uuu)
            
          endif
          
          call error(idnode,95)
          
        endif
        
c     line-printer output every nstbpo steps
        
        if(nstep.eq.1.or.(nstep.gt.1.and.mod(nstep,nstbpo).eq.0))then
          
          call timchk(0,timelp)
          if(idnode.eq.0)then
            
            call get_prntime(hms,timelp,prntim)
            call get_simtime(dec,nstep,tstep,simtim)
            if(mod(lines,npage).eq.0)then
              write(nrite,"(1x,120('-'),
     x        /,/,1x,'    step',5x,'eng_tot',4x,'temp_tot',5x,
     x        'eng_cfg',5x,'eng_vdw',5x,'eng_cou',5x,'eng_bnd',
     x        5x,'eng_ang',5x,'eng_dih',5x,'eng_tet',/,1x,
     x        'time    ',5x,' eng_pv',4x,'temp_rot',5x,'vir_cfg',
     x        5x,'vir_vdw',5x,'vir_cou',5x,'vir_bnd',5x,'vir_ang',
     x        5x,'vir_con',5x,'vir_tet',/,1x,'cpu time',6x,
     x        'volume',4x,'temp_shl',5x,'eng_shl',5x,'vir_shl',
     x        7x,'alpha',8x,'beta',7x,'gamma',5x,'vir_pmf',
     x        7x,'press')")
              if(lpimd)write(nrite,"(1x,'pimd    ',5x,'eng_qpi',5x,
     x        'eng_qvr',5x,'eng_rng',5x,'vir_rng',5x,'qms_rgr',5x,
     x        'qms_bnd',5x,'eng_the')")
              write(nrite,"(/,/,1x,120('-'))")
            endif
            write(nrite,"(1x,i8,1p,9e12.4,/,1x,0p,f7.3,a1,1p,9e12.4,
     x        /,1x,0p,f7.3,a1,1p,9e12.4)")
     x      nstep,(stpval(i),i=1,9),
     x      simtim,dec,(stpval(i),i=10,18),
     x      prntim,hms,(stpval(i),i=19,27)
          iadd=mxnstk-7
          if(lpimd)write(nrite,"(9x,1p,9e12.4)")
     x      (stpval(iadd+i),i=1,7)
          write(nrite,"(/,1x,' rolling',1p,9e12.4,/,1x,'averages',
     x        1p,9e12.4,/,9x,1p,9e12.4)") (ravval(i),i=1,27)
          if(lpimd)write(nrite,"(9x,1p,9e12.4)")
     x      (ravval(iadd+i),i=1,7)
          write(nrite,"(1x,120('-'))")
          
        endif
        
        lines=lines+1
        
      endif
      
c     report end of equilibration period
      
      if((.not.loptim).and.(.not.lzero).and.(nstep.ge.nsteql))then
        
        if((ltscal.and.idnode.eq.0).and.(nstep.eq.nsteql))
     x    write(nrite,"(/,/,1x,'switching off temperature ',
     x        'scaling at step ',i6,/,/,/,1x,120('-'))") nstep
        ltscal=.false.
        
      endif
      
c     write trajectory data
      
      if(ltraj.and.nstep.ge.nstraj) then
        if(idnode.eq.0.and.mod(nstep-nstraj,istraj).eq.0)then
          
          call traject
     x      (ltraj,idnode,imcon,istraj,keytrj,natms*nbeads,
     x      nstraj,nstep,tstep,lpimd)
          
        endif
        
      endif
      
c     write solvation energy file
      
      if(lsolva.and.nstep.ge.nsolva)then
        
        if(mod(nstep-nsolva,isolva).eq.0)then
          
          call solva_temp(idnode,mxnode,natms2,keyver)
          call solvation_write(lexcite,lswitch,idnode,natms,
     x      nstep,nsolva,isolva,tstep,engunit,elrc)
          
        endif
        
      endif
      
c     write free energy file
      
      if(lfree.and.nstep.ge.nfrn)then
        
        if(mod(nstep-nfrn,ifrn).eq.0)then
          
          call free_kinetic(lfrmas,idnode,mxnode,keyver)
          call free_energy_write(idnode,nstep,engunit)
          
        endif
        
      endif
      
c     save restart data in event of system crash
      
      if(mod(nstep,ndump).eq.0.and.nstep.ne.nstrun)then
        
        levcfg=2
        call revive
     x    (lgofr,lzden,lpimd,idnode,imcon,mxnode,natms,nbeads,levcfg,
     x    nstep,nzden,numacc,numrdf,chip,chit,conint,tstep,engcfg,
     x    virtot,vircom,tboost,chit_shl,gaumom)
        
        if(ltad.or.lbpd)
     x    call hyper_close(ltad,idnode,mxnode,natms,nsteql)
        
        if(lpimd)then
          
          call write_thermostats(idnode,mxnode,natms,temp)
          if(keyens.eq.41)call save_rnd_cfg(idnode,mxnode,uuu)
          
        endif
        
      endif
      
c     cycle time check
      
      call timchk(0,timelp)
      recycle=(recycle.and.timjob-timelp.gt.timcls)
      enddo
      
c***********************************************************************
c     end of molecular dynamics calculations
c***********************************************************************
      
c     last time check
      
      call timchk(0,timelp)
      call get_prntime(hms,timjob,prntim)
      if(idnode.eq.0)write(nrite,
     x  "(/,/,1x,'run terminating. elapsed cpu time = ',1p,e13.5,
     x   ', job time = ',0p,f7.3,a1,', close time = ',f7.2,'s',/)")
     x  timelp,prntim,hms,timcls
      
c     shell relaxation convergence statistics
      
      if(.not.loptim.and.keyshl.eq.2)then
        
        if(idnode.eq.0)write(nrite,
     x    "(/,/,1x,'shell relaxation statistics : average cycles = ',
     xf8.3,' maximum cycles = ',f8.3)")pass1,pass2
        
      endif
      
c     produce summary of simulation
      
      levcfg=2
      if(loptim)levcfg=0
      if(.not.lneb)call result
     x  (ltad,lbpd,lgofr,lpgr,lzden,lpimd,idnode,imcon,keyens,mxnode,
     x  natms,nbeads,levcfg,nzden,nstep,ntpatm,numacc,numrdf,keybpd,
     x  chip,chit,conint,rcut,tstep,engcfg,volm,virtot,vircom,zlen,
     x  tboost,chit_shl,gaumom)
      
c     write hyperdynamics restart file
      
      if(ltad.or.lbpd)
     x  call hyper_close(ltad,idnode,mxnode,natms,nsteql)
      
c     write pimd thermostats file
      
      if(lpimd)then
        
        call write_thermostats(idnode,mxnode,natms,temp)
        if(keyens.eq.41)call save_rnd_cfg(idnode,mxnode,uuu)
        
      endif
     
      if(lfqcmd) then 
      call write_fqcmd_rdfs(idnode,mxnode,nsteql,nstrun,rcut,volm,
     x  temp,lpimd,ntpatm,wrtrdf)
      endif
c     close output channels
      
      if(idnode.eq.0)then
        
        close (nrite)
        close (nstats)
        close (nhist)
        close (nevnt)
        close (ntherm)
        close (npuni)
        close (corr)  

      endif
      
c      call dealloc_fqcmd_arrays(idnode,mxnode)      
c     terminate job
      
      call exitcomms()
      
      end
