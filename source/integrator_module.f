      module integrator_module

c***********************************************************************
c     
c     dl_poly module for selecting verlet integration schemes
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c     nvtvv_nhc and nptvv_nhc subroutines are added for NVT and
c     NPT ensembles together with Nose-Hoover Chain thermostat/barostat
c
c     copyright - M.R.Momeni and F.A. Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c     Method Development and Materials Simulation Laboratory
c
c***********************************************************************
      
      use error_module
      use lf_motion_module
      use lf_rotation1_module
      use lf_rotation2_module
      use pmf_module
      use temp_scalers_module
      use vv_motion_module
      use vv_rotation1_module
      use vv_rotation2_module
      use vv_pimd_module
      
      logical, parameter :: gaussian_switch=.true.
      
      contains
      
      subroutine lf_integrate
     x  (lcnb,lshmov,lnfic,idnode,mxnode,imcon,natms,nstep,ngrp,
     x  keyens,nscons,ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nofic,
     x  tstep,engke,engrot,tolnce,quattol,vircon,vircom,virtot,
     x  temp,press,volm,sigma,taut,taup,chit,chip,consv,conint,
     x  elrc,virlrc,virpmf,gaumom)

c***********************************************************************
c     
c     dl_poly subroutine for selecting the integration algorithm
c     to solve the equations of motion. based on the leapfrog
c     verlet algorithm
c     
c     copyright - daresbury laboratory
c     author    - w. smith december 2005
c     
c***********************************************************************

      implicit none

      logical safe,safep,safeq,lcnb,lshmov,lnfic
      integer idnode,mxnode,imcon,natms,ngrp,keyens,nscons,nofic
      integer ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nstep
      real(8) tstep,engke,engrot,tolnce,quattol,vircon,vircom
      real(8) virtot,temp,press,volm,sigma,taut,taup,chit,chip
      real(8) consv,conint,elrc,virlrc,virpmf
      real(8) gaumom(0:5)
      
      safe=.true.
      safeq=.true.
      safep=.true.
      
      if(ngrp.eq.0) then
        
        if(keyens.eq.0) then

c     verlet leapfrog 

          call nve_1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,tolnce,tstep,vircon)
          
        else if(keyens.eq.1) then

c     Evans Gaussian Temperature constraints
          
          call nvt_e1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,tolnce,tstep,vircon)
          
        else if(keyens.eq.2) then

c     Berendsen thermostat
          
          call nvt_b1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      engke,taut,sigma,tolnce,tstep,vircon)
          
        else if(keyens.eq.3) then

c     Nose-Hoover thermostat
          
          call nvt_h1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,nscons,ntcons,
     x      chit,consv,conint,engke,taut,sigma,tolnce,tstep,vircon)
          
        elseif(keyens.eq.4) then

c     Berendsen thermostat and isotropic barostat

          call npt_b1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,elrc,engke,virlrc,press,taup,taut,sigma,tolnce,
     x      tstep,virtot,vircon,volm)

        else if(keyens.eq.5) then

c     Nose-Hoover thermostat and isotropic barostat 

          call npt_h1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,chip,chit,conint,consv,elrc,engke,virlrc,press,
     x      taup,taut,sigma,temp,tolnce,tstep,virtot,vircon,volm)

        else if(keyens.eq.6) then

c     Berendsen thermostat and barostat (cell shape varying)

          call nst_b1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,mode,elrc,engke,virlrc,press,taup,taut,sigma,
     x      tolnce,tstep,vircon,volm)

        else if(keyens.eq.7) then

c     Nose-Hoover thermostat and barostat (cell shape varying)
          
          call nst_h1
     x      (safe,lshmov,idnode,imcon,mxnode,natms,ntpatm,nscons,
     x      ntcons,mode,chit,conint,consv,elrc,engke,virlrc,press,
     x      taup,taut,sigma,temp,tolnce,tstep,vircon,volm)

        elseif(keyens.eq.8) then

c     Potential of mean force in NVE

            call pmflf
     x        (safe,safep,lshmov,idnode,imcon,mxnode,natms,nscons,
     x        ntcons,nspmf,ntpmf,engke,tolnce,tstep,vircon,virpmf)

        endif

      elseif(ngrp.gt.0) then

c     apply rigid body equations of motion
        
        if(keyens.eq.0) then
          
          if(.not.lcnb) then

            call nveq_1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x        vircon)

          else

            call nveq_2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,tolnce,tstep,vircom,
     x        vircon)

          endif

        elseif(keyens.eq.1) then

c     invalid option

          call error(idnode,430)
          
        elseif(keyens.eq.2) then
          
          if(.not.lcnb) then

            call nvtq_b1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,sigma,taut,tolnce,
     x        tstep,vircom,vircon)

          else

            call nvtq_b2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,quattol,sigma,taut,tolnce,
     x        tstep,vircom,vircon)
          
          endif

        elseif(keyens.eq.3) then
          
          if(.not.lcnb) then 

            call nvtq_h1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,chit,consv,conint,engke,engrot,quattol,
     x        sigma,taut,tolnce,tstep,vircom,vircon)

          else

            call nvtq_h2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,conint,consv,chit,engke,engrot,quattol,
     x        sigma,taut,tolnce,tstep,vircom,vircon)

          endif
            
        elseif(keyens.eq.4) then

          if(.not.lcnb) then

            call nptq_b1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,virtot,vircom,
     x        vircon,volm)
          
          else

            call nptq_b2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,
     x        virtot,volm)

          endif

        elseif(keyens.eq.5) then
          
          if(.not.lcnb) then 

            call nptq_h1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,virtot,vircom,vircon,volm)

          else

            call nptq_h2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,chip,chit,consv,conint,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,vircom,vircon,virtot,volm)

          endif
            
        elseif(keyens.eq.6) then

          if(.not.lcnb) then

            call nstq_b1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm)

          else

            call nstq_b2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,elrc,engke,engrot,virlrc,press,
     x        quattol,sigma,taup,taut,tolnce,tstep,vircom,vircon,volm)

          endif

        elseif(keyens.eq.7) then

          if(.not.lcnb) then

            call nstq_h1
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,chit,conint,consv,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,vircom,vircon,volm)

          else

            call nstq_h2
     x        (safe,safeq,lshmov,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,chit,conint,consv,elrc,engke,
     x        engrot,virlrc,press,quattol,sigma,taup,taut,temp,tolnce,
     x        tstep,vircom,vircon,volm)
            
          endif

        else

c     invalid option

          call error(idnode,430)

        endif

      endif
      
c     calculate gaussian moments of momenta
      
      call gaussian_moments(idnode,mxnode,natms,temp,gaumom)
      
c    check on convergence of pmf-shake
      
      if(ntpmf.gt.0) then
        
        if(mxnode.gt.1) call gstate(safep)
        if(.not.safep) call error(idnode,438)

      endif    

c    check on convergence of shake

      if(ntcons.gt.0) then

        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe) call error(idnode,105)

      endif    

c     check on convergence of quaternion algorithm

      if(ngrp.gt.0) then

        if(mxnode.gt.1) call gstate(safeq)
        if(.not.safeq) call error(idnode,321)

      endif

c     eliminate "flying ice cube" in long simulations (Berendsen)
      
      if(lnfic.and.(keyens.eq.2.or.keyens.eq.4.or.keyens.eq.6))then
        
        if(mod(nstep,nofic).eq.0)then
          
          call vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)
          
        endif
        
      endif
      
      return
      end subroutine lf_integrate

      subroutine vv_integrate
     x  (lcnb,lshmov,lnfic,lmsite,isw,idnode,mxnode,imcon,natms,nstep,
     x  nchain,nrespa,ntpmls,ngrp,keyens,nscons,ntcons,ntpatm,ntfree,
     x  nspmf,ntpmf,mode,nofic,ntshl,keyshl,tstep,engke,engrot,tolnce,
     x  vircon,vircom,virtot,temp,press,volm,sigma,sigma_nhc,sigma_volm,
     x  alpha_volm,taut,taup,chit,chip,consv,conint,g_qt4f,
     x  elrc,virlrc,virpmf,chit_shl,sigma_shl,gaumom)

c***********************************************************************
c     
c     dl_poly subroutine for selecting the integration algorithm
c     to solve the equations of motion. based on the velocity
c     verlet algorithm
c     
c     copyright - daresbury laboratory
c     author    - w. smith february 2005
c     
c***********************************************************************

      implicit none

      logical safe,safep,lcnb,lshmov,lnfic,lmsite
      integer isw,idnode,mxnode,imcon,natms,ngrp,keyens,nscons
      integer nchain,nrespa,ntpmls
      integer ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nstep,nofic
      integer ntshl,keyshl
      real(8) tstep,engke,engrot,tolnce,vircon,vircom
      real(8) virtot,temp,press,volm,sigma,taut,taup,chit,chip
      real(8) consv,conint,elrc,virlrc,virpmf,chit_shl,sigma_shl
      real(8) sigma_nhc,sigma_volm,alpha_volm,g_qt4f
      real(8) gaumom(0:5)

      if(ngrp.eq.0) then
        
        if(keyens.eq.0) then

c     velocity verlet

          call nvevv_1
     x      (safe,lshmov,lmsite,isw,idnode,mxnode,natms,imcon,
     x      nscons,ntcons,ntpmls,tstep,engke,tolnce,vircon,g_qt4f)
          
        else if(keyens.eq.1) then

c     Evans Gaussian Temperature constraints
          
          call nvtvv_e1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,tstep,engke,tolnce,vircon)
          
        else if(keyens.eq.2) then

c     Berendsen thermostat
          
          call nvtvv_b1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,tstep,taut,sigma,engke,tolnce,vircon,
     x      lmsite,g_qt4f,ntpmls)
          
        else if(keyens.eq.3) then

c     Nose-Hoover thermostat
          
          call nvtvv_h1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntshl,keyshl,tstep,taut,sigma,chit,consv,
     x      conint,engke,tolnce,vircon,chit_shl,sigma_shl,
     x      lmsite,g_qt4f,ntpmls)
          
        elseif(keyens.eq.4) then

c     Berendsen thermostat and isotropic barostat

          call nptvv_b1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntpatm,tstep,taut,taup,sigma,engke,press,elrc,
     x      virlrc,tolnce,virtot,vircon,volm,
     x      lmsite,ntpmls,g_qt4f)

        else if(keyens.eq.5) then

c     Nose-Hoover thermostat and isotropic barostat 

          call nptvv_h1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,temp,
     x      chip,chit,consv,conint,engke,elrc,tolnce,vircon,
     x      virtot,virlrc,volm,press,chit_shl,sigma_shl,
     x      lmsite,ntpmls,g_qt4f)

        else if(keyens.eq.6) then

c     Berendsen thermostat and barostat (cell shape varying)

          call nstvv_b1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntpatm,mode,tstep,taut,taup,sigma,engke,press,
     x      elrc,virlrc,tolnce,vircon,volm)

        else if(keyens.eq.7) then

c     Nose-Hoover thermostat and barostat (cell shape varying)
          
          call nstvv_h1
     x      (safe,lshmov,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,ntpatm,mode,ntshl,keyshl,tstep,taut,taup,sigma,
     x      temp,chit,consv,conint,engke,elrc,tolnce,vircon,
     x      virlrc,volm,press,chit_shl,sigma_shl)

        elseif(keyens.eq.8) then

c     Potential of mean force in NVE

          call pmfvv
     x      (safe,safep,lshmov,isw,idnode,mxnode,imcon,natms,nscons,
     x      ntcons,nspmf,ntpmf,engke,tolnce,tstep,vircon,virpmf)

        elseif(keyens.eq.9) then

c     Nose-Hoover Chain thermostat in NVT ensemble
        
c          call nvtvv_nhc
c     x      (safe,lmsite,isw,idnode,mxnode,natms,imcon,
c     x      nchain,nrespa,ntpmls,         
c     x      ntshl,keyshl,tstep,taut,sigma,sigma_nhc,g_qt4f,chit,
c     x      consv,conint,engke,vircon,chit_shl,sigma_shl)

           call nvtvv_nhc
     x      (safe,lshmov,lmsite,isw,idnode,mxnode,natms,imcon,nscons,
     x      ntcons,nchain,nrespa,ntpmls,         
     x      ntshl,keyshl,tstep,taut,sigma,sigma_nhc,g_qt4f,chit,
     x      consv,conint,engke,tolnce,vircon,chit_shl,sigma_shl)

        elseif(keyens.eq.10) then

c     Nose-Hoover Chain thermostat in NPT ensemble
        
          call nptvv_nhc
     x    (safe,lmsite,isw,idnode,mxnode,natms,imcon,
     x    nchain,nrespa,ntpmls,        
     x    ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,
     x    sigma_nhc,sigma_volm,alpha_volm,virtot,vircon,virlrc,
     x    g_qt4f,press,volm,chit,consv,
     x    conint,engke,elrc,chit_shl,sigma_shl)

        endif

      elseif(ngrp.gt.0) then

c     apply rigid body equations of motion
        
        if(keyens.eq.0) then
          
          if(.not.lcnb) then

            call nveqvv_1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,tolnce,tstep,vircom,vircon)

          else

            call nveqvv_2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,tolnce,tstep,vircom,vircon)
            
          endif

        elseif(keyens.eq.1) then

c     invalid option

          call error(idnode,430)
          
        elseif(keyens.eq.2) then
          
          if(.not.lcnb) then

            call nvtqvv_b1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,taut,sigma,tolnce,tstep,
     x        vircom,vircon)

          else

            call nvtqvv_b2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,engke,engrot,taut,sigma,tolnce,tstep,
     x        vircom,vircon)
            
          endif

        elseif(keyens.eq.3) then
          
          if(.not.lcnb) then 

            call nvtqvv_h1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntshl,keyshl,chit,consv,conint,engke,
     x        engrot,taut,sigma,tolnce,tstep,vircom,vircon,chit_shl,
     x        sigma_shl)

          else

            call nvtqvv_h2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntshl,keyshl,chit,consv,conint,engke,
     x        engrot,taut,sigma,tolnce,tstep,vircom,vircon,chit_shl,
     x        sigma_shl)

          endif
          
        elseif(keyens.eq.4) then

          if(.not.lcnb) then

            call nptqvv_b1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,engke,engrot,press,taut,taup,sigma,
     x        tolnce,tstep,vircom,vircon,elrc,virlrc,virtot,volm)
            
          else

            call nptqvv_b2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,engke,engrot,press,taut,taup,sigma,
     x        tolnce,tstep,vircom,vircon,elrc,virlrc,virtot,volm)

          endif

        elseif(keyens.eq.5) then
          
          if(.not.lcnb) then 

            call nptqvv_h1
     x        (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x        ntcons,ntpatm,ntfree,ntshl,keyshl,tstep,taut,taup,sigma,
     x        temp,chip,chit,consv,conint,engke,engrot,elrc,tolnce,
     x        vircon,virtot,virlrc,vircom,volm,press,chit_shl,
     x        sigma_shl)

          else
      
            call nptqvv_h2
     x        (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x        ntcons,ntpatm,ntfree,ntshl,keyshl,tstep,taut,taup,sigma,
     x        temp,chip,chit,consv,conint,engke,engrot,elrc,tolnce,
     x        vircom,vircon,virtot,virlrc,volm,press,chit_shl,
     x        sigma_shl)
      
          endif
          
        elseif(keyens.eq.6) then

          if(.not.lcnb) then

            call nstqvv_b1
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,engke,engrot,press,taut,taup,
     x        sigma,tolnce,tstep,vircom,vircon,elrc,virlrc,volm)

          else

            call nstqvv_b2
     x        (safe,lshmov,isw,imcon,idnode,mxnode,natms,ngrp,nscons,
     x        ntcons,ntfree,ntpatm,mode,engke,engrot,press,taut,taup,
     x        sigma,tolnce,tstep,vircom,vircon,elrc,virlrc,volm)

          endif

        elseif(keyens.eq.7) then

          if(.not.lcnb) then

            call nstqvv_h1
     x        (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x        ntcons,ntpatm,ntfree,mode,ntshl,keyshl,tstep,taut,taup,
     x        sigma,temp,chit,consv,conint,engke,engrot,elrc,tolnce,
     x        vircon,virlrc,vircom,volm,press,chit_shl,sigma_shl)

          else

            call nstqvv_h2
     x        (safe,lshmov,isw,idnode,mxnode,natms,imcon,ngrp,nscons,
     x        ntcons,ntpatm,ntfree,mode,ntshl,keyshl,tstep,taut,taup,
     x        sigma,temp,chit,consv,conint,engke,engrot,elrc,tolnce,
     x        vircom,vircon,virlrc,volm,press,chit_shl,sigma_shl)
            
          endif

        else

c     invalid option

          call error(idnode,430)

        endif

      endif
      
c     calculate gaussian moments of momenta
      
      if(isw.eq.2)then
        
        call gaussian_moments(idnode,mxnode,natms,temp,gaumom)

      endif
      
c     check on convergence of pmf-shake

      if(ntpmf.gt.0) then

        if(mxnode.gt.1) call gstate(safep)
        if(.not.safep) call error(idnode,438)

      endif    

c     check on convergence of shake

      if(ntcons.gt.0) then

        if(mxnode.gt.1) call gstate(safe)
        if(.not.safe) call error(idnode,105)

      endif    

c     eliminate "flying ice cube" in long simulations (Berendsen)
      
      if(lnfic.and.(keyens.eq.2.or.keyens.eq.4.or.keyens.eq.6))then
        
        if(mod(nstep,nofic).eq.0)then
          
          call vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)
          
        endif
        
      endif
      
      return
      end subroutine vv_integrate

c      subroutine pimd_integrate
c     x  (lmsite,isw,idnode,mxnode,imcon,ntpmls,natms,keyens,nstep,tstep,
c     x  taut,g_qt4f,temp,engke,engthe,chi,uuu,gaumom)
      
c***********************************************************************
c     
c     dl_poly subroutine for selecting the integration algorithm
c     to solve the pimd equations of motion. based on the velocity
c     verlet algorithm
c     
c     copyright - daresbury laboratory
c     author    - w. smith september 2016
c     
c***********************************************************************

c      implicit none

c      logical lmsite
c      integer isw,idnode,mxnode,imcon,ntpmls,natms,keyens,nstep
c      real(8) tstep,engke,engthe,temp,taut,chi,g_qt4f
c      real(8) uuu(102),gaumom(0:5)
      
       subroutine pimd_integrate
     x      (lmsite,isw,idnode,mxnode,imcon,ntpmls,natms,keyens,nstep,
     x      tstep,g_qt4f,temp,engke,engthe,chi,uuu,gaumom,
     x      safe,nrespa,ntpatm,ntshl,keyshl,taut,taup,sigma,
     x      sigma_nhc,sigma_volm,alpha_volm,virtot,vircon,virlrc,virrng,
     x      press,volm,chit,consv,conint,elrc,chit_shl,sigma_shl)
      
c***********************************************************************
c     
c     dl_poly subroutine for selecting the integration algorithm
c     to solve the equations of motion. based on the velocity
c     verlet algorithm
c     
c     copyright - daresbury laboratory
c     author    - w. smith february 2005
c     
c***********************************************************************

      implicit none

      logical safe,safep,lcnb,lshmov,lnfic,lmsite
      integer isw,idnode,mxnode,imcon,natms,ngrp,keyens,nscons
      integer nchain,nrespa,ntpmls
      integer ntcons,ntpatm,ntfree,nspmf,ntpmf,mode,nstep,nofic
      integer ntshl,keyshl
      real(8) tstep,engke,engrot,tolnce,vircon,vircom,virrng
      real(8) virtot,temp,press,volm,sigma,taut,taup,chit,chip
      real(8) consv,conint,elrc,virlrc,virpmf,chit_shl,sigma_shl
      real(8) sigma_nhc,sigma_volm,alpha_volm,g_qt4f
      real(8) gaumom(0:5)
      real(8) :: uuu(102),chi,engthe

      
      engthe=0.d0
      
      if(keyens.eq.40) then
        
c     nvt ensemble
        
        call pimd_nvt
     x  (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,temp,
     x  g_qt4f,engke,engthe)
        
      elseif(keyens.eq.41) then
        
c     nvt ensemble - gentle thermostat
        
        call pimd_nvt_gth1
     x  (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,temp,
     x  g_qt4f,engke,engthe,chi,uuu)
        
      else if(keyens.eq.42) then
        
c     nvt ensemble - nose-hoover chains
        
        call pimd_nvt_nhc
     x    (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,
     x    g_qt4f,temp,engke,engthe,nrespa)
        
      else if(keyens.eq.43) then
        
c     nvt ensemble - nose-hoover chains - normal mode
        
        call pimd_nvt_nhc_nm
     x    (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,
     x    g_qt4f,temp,engke,engthe,nrespa,sigma)
        
      else if(keyens.eq.44) then
        
c     nvt ensemble - PILE thermostat - normal mode
        
        call pimd_nvt_pile_nm
     x    (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,
     x    g_qt4f,temp,engke,engthe,uuu)

      else if(keyens.eq.45) then
        
c     nvt ensemble - piglet thermostat - normal mode
        
        call pimd_nvt_piglet
     x    (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,
     x    g_qt4f,temp,engke,engthe,uuu,virtot,vircon,virrng)

      else if(keyens.eq.51) then
        
c     npt ensemble - nhc barostat - normal mode
        
        call pimd_npt_nhc_nm
     x    (safe,lmsite,isw,idnode,mxnode,natms,imcon,
     x    nrespa,ntpmls,        
     x    ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,
     x    sigma_nhc,sigma_volm,alpha_volm,virtot,vircon,virlrc,
     x    g_qt4f,press,volm,chit,consv,
     x    conint,engke,elrc,chit_shl,sigma_shl,temp,engthe)
      
      else if(keyens.eq.52) then
        
c     npt ensemble - PILE thermostat - normal mode
        
        call pimd_npt_pile_nm
     x    (safe,lmsite,isw,idnode,mxnode,natms,imcon,
     x    ntpmls,ntpatm,ntshl,keyshl,tstep,taut,taup,sigma,
     x    sigma_volm,alpha_volm,virtot,vircon,virlrc,
     x    g_qt4f,press,volm,chit,consv,
     x    conint,engke,elrc,chit_shl,sigma_shl,temp,engthe,uuu)

      else if(keyens.eq.61) then
        
c     nve ensemble (RPMD) - normal mode

        call pimd_nve
     x    (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,
     x    g_qt4f,temp,engke)
      
      else if(keyens.eq.62) then

c     PACMD - 

        call pacmd
     x    (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,taut,
     x    g_qt4f,temp,engke,engthe,nrespa)
       
      else if(keyens.eq.63) then

c     TRPMD
        call trpmd
     x    (lmsite,isw,idnode,mxnode,natms,imcon,ntpmls,tstep,
     x    g_qt4f,temp,engke)

      else 
c     invalid option
        
        call error(idnode,430)
        
      endif

      if(isw.eq.2)then
        
        if(gaussian_switch)then
          
c     calculate gaussian moments of bead momentum
          
          call bead_moments
     x      (idnode,mxnode,natms,nbeads,temp,gaumom)
          
        else
          
c     calculate gaussian moments of thermostats
          if(keyens.ne.44)then
          call thermostat_moments
     x      (idnode,mxnode,natms,nbeads,nchain,temp,taut,gaumom)
          endif
        endif

      endif
      
      return
      end subroutine pimd_integrate

      subroutine gaussian_moments(idnode,mxnode,natms,temp,gaumom)

c***********************************************************************
c     
c     dl_poly module for calculating the gaussian moments of atomic
c     momenta
c     copyright - daresbury laboratory
c     author    - w. smith    dec 2016
c     
c***********************************************************************
      
      implicit none

      integer idnode,mxnode,natms,iatm0,iatm1,i
      real(8) ppp,ff1,ff2,temp,gaumom(0:5)
      
      iatm0=(idnode*natms)/mxnode+1
      iatm1=((idnode+1)*natms)/mxnode
      
      do i=iatm0,iatm1
        
        if(lstfrz(i).gt.0)cycle
        
        ff1=gaumom(0)+1.d0
        ff2=gaumom(0)/(gaumom(0)+1.d0)
        ppp=weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)/(boltz*temp)
        gaumom(0)=ff1
        gaumom(1)=ppp/(3.d0*ff1)+gaumom(1)*ff2
        gaumom(2)=ppp**2/(15.d0*ff1)+gaumom(2)*ff2
        gaumom(3)=ppp**3/(105.d0*ff1)+gaumom(3)*ff2
        gaumom(4)=ppp**2*ppp**2/(945.d0*ff1)+gaumom(4)*ff2
        gaumom(5)=ppp**3*ppp**2/(10395.d0*ff1)+gaumom(5)*ff2
        
      enddo
      
      end subroutine gaussian_moments

      subroutine bead_moments(idnode,mxnode,natms,nbeads,temp,gaumom)

c**********************************************************************
c     
c     dl_poly_classic subroutine for calculating the moments of the 
c     bead momentum distribution function in a path intergral molecular
c     dynamics simulation
c     
c     copyright - daresbury laboratory
c     author    - w.smith dec 2016
c     
c**********************************************************************
      
      implicit none
      
      integer, intent(in ) :: idnode,mxnode,natms,nbeads
      real(8), intent(in)  :: temp

      integer i,iatm0,iatm1
      real(8) ppp
      real(8) gaumom(0:5)
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)

      do i=1,iatm1-iatm0
        
        ppp=(pxx(i)**2+pyy(i)**2+pzz(i)**2)*rzmass(i)/(boltz*temp)
        gaumom(0)=gaumom(0)+1.d0
        gaumom(1)=ppp/(3.d0*gaumom(0))+
     x    gaumom(1)*(gaumom(0)-1.d0)/gaumom(0)
        gaumom(2)=ppp**2/(15.d0*gaumom(0))+
     x    gaumom(2)*(gaumom(0)-1.d0)/gaumom(0)
        gaumom(3)=ppp**3/(105.d0*gaumom(0))+
     x    gaumom(3)*(gaumom(0)-1.d0)/gaumom(0)
        gaumom(4)=ppp**2*ppp**2/(945.d0*gaumom(0))+
     x    gaumom(4)*(gaumom(0)-1.d0)/gaumom(0)
        gaumom(5)=ppp**3*ppp**2/(10395.d0*gaumom(0))+
     x    gaumom(5)*(gaumom(0)-1.d0)/gaumom(0)
        
      enddo

      end subroutine bead_moments

      subroutine thermostat_moments
     x  (idnode,mxnode,natms,nbeads,nchain,temp,taut,gaumom)

c**********************************************************************
c     
c     dl_poly_classic subroutine for calculating the moments of the 
c     thermostat momentum distribution function in a path intergral 
c     molecular dynamics simulation
c     
c     copyright - daresbury laboratory
c     author    - w.smith dec 2016
c     
c**********************************************************************
      
      implicit none
      
      integer, intent(in ) :: idnode,mxnode,natms,nbeads,nchain
      real(8), intent(in)  :: temp,taut

      integer i,j,k,iatm0,iatm1
      real(8) ppp,qqq,qqk
      real(8) gaumom(0:5)
      
      iatm0=nbeads*((idnode*natms)/mxnode)
      iatm1=nbeads*(((idnode+1)*natms)/mxnode)
      
      qqq=boltz*temp*taut**2
      qqk=hbar**2/(boltz*temp*dble(nbeads))

      
      do k=1,nbeads
        
        do i=k,iatm1-iatm0,nbeads
          
          do j=1,nchain
            
            ppp=(pcx(j,i)**2+pcy(j,i)**2+pcz(j,i)**2)/
     x        (qqq*boltz*temp)
            gaumom(0)=gaumom(0)+1.d0
            gaumom(1)=ppp/(3.d0*gaumom(0))+
     x        gaumom(1)*(gaumom(0)-1.d0)/gaumom(0)
            gaumom(2)=ppp**2/(15.d0*gaumom(0))+
     x        gaumom(2)*(gaumom(0)-1.d0)/gaumom(0)
            gaumom(3)=ppp**3/(105.d0*gaumom(0))+
     x        gaumom(3)*(gaumom(0)-1.d0)/gaumom(0)
            gaumom(4)=ppp**2*ppp**2/(945.d0*gaumom(0))+
     x        gaumom(4)*(gaumom(0)-1.d0)/gaumom(0)
            gaumom(5)=ppp**3*ppp**2/(10395.d0*gaumom(0))+
     x        gaumom(5)*(gaumom(0)-1.d0)/gaumom(0)
            
          enddo
          
        enddo
        
        qqq=qqk

      enddo

      end subroutine thermostat_moments
      
      end module integrator_module
