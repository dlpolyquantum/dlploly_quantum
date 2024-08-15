      module forces_module
      
c***********************************************************************
c     
c     dl_poly module for calculation of atomic forces
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     adapted   - d. quigley : metadynamics
c     
c***********************************************************************
      
      use config_module
      use coulomb_module
      use ensemble_tools_module
      use error_module
      use ewald_module
      use exclude_module
      use external_field_module
      use four_body_module
      use hkewald_module
      use metafreeze_module
      use metal_module
      use neu_coul_module
      use nlist_builders_module
      use pair_module
      use property_module
      use setup_module
      use solvation_module
      use spme_module
      use tersoff_module
      use three_body_module
      use utility_module
      use vdw_module
      use fqcmd_module

      contains
      
      subroutine force_manager
     x  (newlst,lneut,lnsq,lgofr,lzeql,loglnk,lfcap,lsolva,lfree,
     x  lghost,lpimd,idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,
     x  numrdf,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,ntpmet,
     x  nospl,multt,nneut,ntptbp,ntpfbp,ntpter,keyshl,keyfld,
     x  ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,nsolva,isolva,
     x  nbeads,delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw,
     x  vircpe,virsrp,alpha,drewd,volm,engmet,virmet,elrc,virlrc,
     x  rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp,rctter,engter,
     x  virter,engbnd,virbnd,engang,virang,engdih,virdih,enginv,
     x  virinv,engtet,virtet,engshl,shlke,virshl,engfld,virfld,
     x  engcfg,fmax,temp,engord,virord,engrng,virrng,qmsbnd,keyens,
     x  lfqcmd)
      
c*********************************************************************
c     
c     dl_poly subroutine to manage the calculation of the atomic forces
c     from all force field terms.
c     
c     copyright - daresbury laboratory
c     author    - w.smith
c     
c*********************************************************************
      
      use pimd_module, only : ring_forces,ring_energy

      implicit none
      
      logical newlst,lneut,lnsq,lgofr,lzeql,loglnk,lfcap,lsolva
      logical lfree,lghost,llsolva,lpimd,lfqcmd,safe
      
      integer idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,numrdf
      integer keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,ntpmet,ksatms
      integer i,j,k,m,n,nospl,multt,nneut,ntbond,ntangl,ntdihd,nsolva
      integer isolva,ntinv,ntteth,ntshl,ntptbp,ntpfbp,ntpter,keyshl
      integer keyfld,nbeads,nn,nsatm,ibase,nbase,numatm
      integer fail(1:8)
      integer keyens
      
      real(8) delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw
      real(8) vircpe,virsrp,alpha,drewd,volm,engmet,virmet,qfactor
      real(8) elrc,virlrc,rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp
      real(8) rctter,engter,virter,engbnd,virbnd,engang,virang,engdih
      real(8) virdih,enginv,virinv,engtet,virtet,engshl,virshl,engfld
      real(8) virfld,fmax,temp,shlke,engcfg,tmpeng,tmpvir,engord,virord
      real(8) engcpt,engsrt,vircpt,virsrt,engrng,virrng,qmsbnd,engtmp
      real(8) virtmp
      
      integer, allocatable :: lstbak(:,:),lenbak(:)
      real(8), allocatable :: xxxbak(:),yyybak(:),zzzbak(:)
      real(8), allocatable :: fxxbak(:),fyybak(:),fzzbak(:)
      
c     allocate storage arrays
      
      if(lpimd)then
        
        safe=.true.
        fail(:)=0
        allocate (xxxbak(1:natms),          stat=fail(1))
        allocate (yyybak(1:natms),          stat=fail(2))
        allocate (zzzbak(1:natms),          stat=fail(3))
        allocate (fxxbak(1:natms),          stat=fail(4))
        allocate (fyybak(1:natms),          stat=fail(5))
        allocate (fzzbak(1:natms),          stat=fail(6))
        allocate (lenbak(1:msatms),          stat=fail(7))
        allocate (lstbak(1:msatms,1:mxlist), stat=fail(8))
        if(any(fail.gt.0))safe=.false.
        if(mxnode.gt.1)call gstate(safe)
        if(.not.safe)call error(idnode,519)
        
      endif
      
      llsolva=.false.
      nsatm=(natms-idnode-1)/mxnode+1
      
c     pimd factor for bead interactions
      
      qfactor=1.d0/dble(nbeads)
      
c     initialize energy and virial accumulators
      
      engbnd=0.d0
      virbnd=0.d0
      engang=0.d0
      virang=0.d0
      engdih=0.d0
      virdih=0.d0
      enginv=0.d0
      virinv=0.d0
      engtbp=0.d0
      virtbp=0.d0
      engter=0.d0
      virter=0.d0
      engfbp=0.d0
      virfbp=0.d0
      engsrp=0.d0
      virsrp=0.d0
      engcpe=0.d0
      vircpe=0.d0
      engfld=0.d0
      virfld=0.d0
      engshl=0.d0
      virshl=0.d0
      shlke =0.d0
      engtet=0.d0
      virtet=0.d0
      engmet=0.d0
      virmet=0.d0
      virord=0.0d0
      engord=0.0d0
      engrng=0.d0
      virrng=0.d0
      
      if(lmetadyn)then
        
        eng_loc=0.0d0
        vir_loc=0.0d0
        fxx_loc(:)=0.0d0
        fyy_loc(:)=0.0d0
        fzz_loc(:)=0.0d0
        stress_loc(:)=0.0d0
        
      endif
      
c     initialise free energy accumulators
      
      if(lfree)then
        
        ang_fre=0.d0
        bnd_fre=0.d0
        dih_fre=0.d0
        inv_fre=0.d0
        tbp_fre=0.d0
        fbp_fre=0.d0
        cou_fre=0.d0
        vdw_fre=0.d0
        shl_fre=0.d0
        ang_vir=0.d0
        bnd_vir=0.d0
        dih_vir=0.d0
        inv_vir=0.d0
        tbp_vir=0.d0
        fbp_vir=0.d0
        cou_vir=0.d0
        vdw_vir=0.d0
        shl_vir=0.d0
        eng_cfg_fre=0.d0
        vir_cfg_fre=0.d0
        
      endif
      
c     initialise solvation and excitation arrays
      
      if(lsolva)then
        
        if(keyfce/2.gt.0)lcomp(6)=.true.
        if(mod(keyfce,2).eq.1)lcomp(7)=.true.
        if(mod(nstep-nsolva,isolva).eq.0)then
          
          llsolva=.true.
          cou_sol(:)=0.d0
          vdw_sol(:)=0.d0
          
          if(lghost)then
            
            cou_exc(:)=0.d0
            vdw_exc(:)=0.d0
            
          endif
          
        endif
        
      endif
      
c     zero stress tensor
      
      if(nstep.gt.0)then
        
        do i=1,9
          
          stress(i)=0.d0
          
        enddo
        
      endif
      
c     start of loop over pimd beads
      
      if(lpimd)then
        
c     pimd: store original configuration data
        
        do i=1,natms
          
          xxxbak(i)=xxx(i)
          yyybak(i)=yyy(i)
          zzzbak(i)=zzz(i)
          
        enddo
        
        do n=1,nsatm
          
          lenbak(n)=lentry(n)
          
        enddo
        
        do n=1,nsatm
          do j=1,lentry(n)
            
            lstbak(n,j)=list(n,j)
            
          enddo
        enddo
        
      endif
      
c     cycle through the beads of the path integral ring polymer
c     this applies only to interactions defined by atom types
      
      do k=1,nbeads
        
        engsrt=0.d0
        virsrt=0.d0
        engcpt=0.d0
        vircpt=0.d0
        
c     initialise the force arrays
        
        do i=1,natms
          
          fxx(i)=0.d0
          fyy(i)=0.d0
          fzz(i)=0.d0
          
        enddo
        
c     select beads with identity k
        
        if(lpimd.and.k.gt.1)then
          
          ibase=(k-1)*natms
          
          do i=1,natms
            
            xxx(i)=xxx(i+ibase)
            yyy(i)=yyy(i+ibase)
            zzz(i)=zzz(i+ibase)
            
          enddo
          
          nbase=(k-1)*nsatm
          
          do n=1,nsatm
            
            lentry(n)=lentry(n+nbase)
            
          enddo
            
          do n=1,nsatm
            do m=1,lentry(n)
              
              list(n,m)=list(n+nbase,m)-ibase
              
            enddo
          enddo
          
        endif
        
        if(keyfce.gt.0)then
          
c     calculate pair forces, including coulombic forces
          
          if(lnsq)then
            
c     multiple timestep - all-pairs
            
            call multiple_nsq
     x        (lnsq,lgofr,lzeql,newlst,lsolva,lfree,lghost,idnode,
     x        imcon,keyfce,multt,mxnode,natms,nstep,nstbgr,nsteql,
     x        numrdf,nsolva,isolva,delr,dlrpot,engcpt,engsrt,epsq,
     x        rcut,rprim,rvdw,vircpt,virsrt)
            
            engsrp=engsrp+engsrt*qfactor
            engcpe=engcpe+engcpt*qfactor
            virsrp=virsrp+virsrt*qfactor
            vircpe=vircpe+vircpt*qfactor
            
          elseif(.not.lneut)then         
            
c     single timestep
            
            if(multt.eq.1)then
              
              call forces
     x          (loglnk,lgofr,lzeql,lsolva,lfree,lghost,idnode,imcon,
     x          keyfce,kmax1,kmax2,kmax3,nhko,nlatt,mxnode,ntpvdw,
     x          ntpmet,natms,nstbgr,nstep,nsteql,numrdf,nospl,nsolva,
     x          isolva,alpha,dlrpot,drewd,engcpt,engsrt,epsq,rcut,rvdw,
     x          vircpt,virsrt,volm,engmet,virmet)
              
              engsrp=engsrp+engsrt*qfactor
              engcpe=engcpe+engcpt*qfactor
              virsrp=virsrp+virsrt*qfactor
              vircpe=vircpe+vircpt*qfactor
              
            else
              
              call multiple
     x          (loglnk,lgofr,lzeql,newlst,lsolva,lfree,lghost,idnode,
     x          imcon,keyfce,nlatt,kmax1,kmax2,kmax3,nhko,multt,
     x          mxnode,natms,nstep,nstbgr,nsteql,numrdf,nospl,nsolva,
     x          isolva,alpha,dlrpot,drewd,engcpt,engsrt,epsq,rcut,
     x          rprim,rvdw,vircpt,virsrt,volm)
              
              engsrp=engsrp+engsrt*qfactor
              engcpe=engcpe+engcpt*qfactor
              virsrp=virsrp+virsrt*qfactor
              vircpe=vircpe+vircpt*qfactor
              
            endif
            
          elseif(lneut)then
            
c     neutral groups
            
            if(multt.eq.1)then
              
              call forces_neu
     x          (lgofr,lzeql,lsolva,lfree,lghost,idnode,imcon,keyfce,
     x          mxnode,natms,nneut,nstbgr,nstep,nsteql,numrdf,nsolva,
     x          isolva,dlrpot,engcpt,engsrt,epsq,rcut,rvdw,alpha,
     x          vircpt,virsrt)
              
              engsrp=engsrp+engsrt
              engcpe=engcpe+engcpt
              virsrp=virsrp+virsrt
              vircpe=vircpe+vircpt
              
            else
              
              call multiple_neu
     x          (lgofr,lzeql,newlst,lsolva,lfree,lghost,idnode,imcon,
     x          keyfce,multt,mxnode,natms,nneut,nstbgr,nstep,nsteql,
     x          numrdf,nsolva,isolva,delr,dlrpot,engcpt,engsrt,epsq,
     x          rprim,rcut,rvdw,alpha,vircpt,virsrt)
              
              engsrp=engsrp+engsrt
              engcpe=engcpe+engcpt
              virsrp=virsrp+virsrt
              vircpe=vircpe+vircpt
              
            endif
            
          endif
          
        endif
        
c     calculate intramolecular forces
        call intra_forces
     x    (llsolva,lfree,lghost,idnode,mxnode,imcon,natms,nbeads,nstep,
     x    keyfce,keyshl,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,epsq,
     x    engbnd,virbnd,engang,virang,engcpe,vircpe,engsrp,virsrp,
     x    engdih,virdih,enginv,virinv,engtet,virtet,engshl,virshl,
     x    dlrpot,rcut,rvdw,alpha)
        
c     external field
        
        if(keyfld.gt.0)then
          
          call extnfld
     x      (idnode,imcon,keyfld,mxnode,natms,engtmp,virtmp)
          
          engfld=engtmp*qfactor
          virfld=virtmp*qfactor
          
        endif
        
c     calculate three body forces
        
        if(ntptbp.gt.0)then
          
          call thbfrc
     x      (llsolva,lfree,lghost,idnode,mxnode,natms,imcon,rcuttb,
     x      engtmp,virtmp)
          
          engtbp=engtbp+engtmp*qfactor
          virtbp=virtbp+virtmp*qfactor
          
        endif
        
c     calculate four body forces
        
        if(ntpfbp.gt.0)then
          
          call fbpfrc
     x      (llsolva,lfree,lghost,idnode,mxnode,natms,imcon,rcutfb,
     x      engfbp,virfbp)
          
          engfbp=engfbp+engtmp*qfactor
          virfbp=virfbp+virtmp*qfactor
          
        endif
        
c     calculate tersoff potential forces
        
        if(ntpter.gt.0)then
          
          call tersoff
     x      (idnode,mxnode,natms,imcon,rctter,engtmp,virtmp)
          
          engter=engter+engtmp*qfactor
          virter=virter+virtmp*qfactor
          
        endif
        
        if(lpimd)then
          
          ibase=(k-1)*natms
          
          if(k.eq.1)then
            
c     store forces for bead 1
            
            do i=1,natms
              
              fxxbak(i)=fxx(i)*qfactor
              fyybak(i)=fyy(i)*qfactor
              fzzbak(i)=fzz(i)*qfactor
              
            enddo
            
          else
            
c     store forces for bead k>1
            
            ibase=(k-1)*natms
            
            do i=1,natms
              
              fxx(i+ibase)=fxx(i)*qfactor
              fyy(i+ibase)=fyy(i)*qfactor
              fzz(i+ibase)=fzz(i)*qfactor
              
            enddo
            
          endif
          
        endif
        
c     end loop over pimd beads
        
      enddo
      
c     restore original configuration data
      
      if(lpimd)then
        
        do i=1,natms
          
          xxx(i)=xxxbak(i)
          yyy(i)=yyybak(i)
          zzz(i)=zzzbak(i)
          fxx(i)=fxxbak(i)
          fyy(i)=fyybak(i)
          fzz(i)=fzzbak(i)
          
        enddo
        
        do n=1,nsatm
          
          lentry(n)=lenbak(n)
          
        enddo
        
        do n=1,nsatm
          do j=1,lentry(n)
            
            list(n,j)=lstbak(n,j)
            
          enddo
        enddo
        
c     end of loop over pimd beads
        
      endif
      
c     add in long range corrections to energy and pressure
      
      engsrp=engsrp+elrc
      virsrp=virsrp+virlrc
      engmet=engmet+elrcm(0)
      virmet=virmet+vlrcm(0)
      if(lfree)then
        
        vdw_fre=vdw_fre+elrc_fre
        vdw_vir=vdw_vir+vlrc_fre
        
      endif
      
c     deallocate arrays
      
      if(lpimd)then

        fail(:)=0
        safe=.true.
        deallocate(xxxbak,yyybak,zzzbak,stat=fail(1))
        deallocate(fxxbak,fyybak,fzzbak,stat=fail(2))
        deallocate(lstbak,lenbak,stat=fail(3))
        if(any(fail.gt.0))safe=.false.
        if(mxnode.gt.1)call gstate(safe)
        if(.not.safe)call error(idnode,528)
        
      endif
      
c     pimd simulations: reduce stress tensor by bead number
      
      if(lpimd)stress(:)=stress(:)*qfactor
      
c     metadynamics option : use potential energy as order parameter
      
      if(lmetadyn)then
        
        tmpeng=engsrp+engcpe+engbnd+engang+engdih+engfld+
     x    engtbp+engfbp+engshl+enginv+engter+engmet
        
        tmpvir=vircpe+virsrp+virbnd+virtbp+virter+virfld+
     x    virang+virshl+virtet+virmet
        
        call metafreeze_driver
     x    (imcon,natms,temp,nstep,tmpeng,tmpvir,engord,virord)
        
      endif
      
c     calculate ring forces for pimd option
      
      if(lpimd)then
        
        if(keyens.le.42) then
                
          call ring_forces
     x      (idnode,mxnode,natms,temp,engrng,virrng,qmsbnd,stress)

        else

          call ring_energy
     x      (idnode,mxnode,natms,temp,engrng,virrng,qmsbnd,stress)

        endif
c        write(6,*)"engrng",engrng
      endif
     
      if(lfqcmd.and.(.not.lpimd))then
        call fqcmd_correct_force(idnode,mxnode,imcon,natms,
     x  engsrp)
      endif 
      
c     global summation of force arrays (basic replicated data strategy)
      
      numatm=nbeads*natms
      call global_sum_forces(numatm,mxnode,fxx,fyy,fzz)
c       write(*,*) "fxx after ", fxx(:)
      
c     global sum of stress arrays
      
      if(mxnode.gt.1)call gdsum(stress,9,buffer)
      
c     add long range correction to diagonal terms of stress tensor
      
      stress(1)=stress(1)-(virlrc+vlrcm(0))/3.d0
      stress(5)=stress(5)-(virlrc+vlrcm(0))/3.d0
      stress(9)=stress(9)-(virlrc+vlrcm(0))/3.d0
      
c     cap forces in equilibration mode
      
      if(nstep.le.nsteql.and.lfcap)
     x  call fcap(lfcap,natms,fmax,temp)
      
c     total configuration energy
      
      engcfg=engsrp+engcpe+engbnd+engang+engdih+engfld+engtbp+
     x  engfbp+engshl+enginv+engter+engmet+engtet
     
c      write(nrite,"('conf. enrergy',e12.4)")engcfg

c     total derivative of the configurational free energy
      
      if(lfree)then
        
        eng_cfg_fre=dlambda*(ang_fre+bnd_fre+dih_fre+inv_fre+
     x    tbp_fre+fbp_fre+cou_fre+vdw_fre+shl_fre)
        vir_cfg_fre=dlambda*(ang_vir+bnd_vir+dih_vir+inv_vir+
     x    tbp_vir+fbp_vir+cou_vir+vdw_vir+shl_vir)
        
      endif
      
c     sum solvation and excitation energies for pair forces
      
      if(mxnode.gt.1)then
        
        if(llsolva)then
          
          call gdsum(vdw_sol,mxtmls_sol2,buffer)
          call gdsum(cou_sol,mxtmls_sol2,buffer)
          
          if(lghost)then
            
            call gdsum(vdw_exc,mxtmls_exc2,buffer)
            call gdsum(cou_exc,mxtmls_exc2,buffer)
            
          endif
          
        endif
        
      endif
      
c     add long range corrections to solvation terms
      
      if(lsolva)then
        
        vdw_sol(:)=vdw_sol(:)+elrc_sol(:)
        if(lghost)vdw_exc(:)=vdw_exc(:)+elrc_exc(:)
        
      endif
      
      return
      end subroutine force_manager

      subroutine intra_forces
     x  (llsolva,lfree,lghost,idnode,mxnode,imcon,natms,nbeads,nstep,
     x  keyfce,keyshl,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,epsq,
     x  engbnd,virbnd,engang,virang,engcpe,vircpe,engsrp,virsrp,engdih,
     x  virdih,enginv,virinv,engtet,virtet,engshl,virshl,dlrpot,rcut,
     x  rvdw,alpha)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating intramolecular forces
c     including: bonds, angles, dihedrals, inversions, tethers
c     and core-shell forces
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory
c     author    - w. smith march 2016.
c     
c***********************************************************************

      implicit none
      
      logical llsolva,lfree,lghost
      integer idnode,mxnode,imcon,ntbond,ntangl,ntdihd,ntinv,ntteth
      integer ntshl,natms,nbeads,nstep,keyfce,keyshl
      real(8) epsq,qfactor,engtmp,virtmp,engbnd,virbnd,engang,virang
      real(8) engcpe,vircpe,engsrp,virsrp,engdih,virdih,engcpt,vircpt
      real(8) engsrt,virsrt,enginv,virinv,engtet,virtet,engshl,virshl
      real(8) dlrpot,rcut,rvdw,alpha

      qfactor=1.d0/dble(nbeads)
      
c     calculate bond forces
      
      if(ntbond.gt.0)then
        
        call bndfrc
     x    (llsolva,lfree,lghost,idnode,imcon,mxnode,ntbond,epsq,
     x    engtmp,virtmp)
        
        engbnd=engbnd+engtmp*qfactor
        virbnd=virbnd+virtmp*qfactor
c        write(6,*) "engbnd",engbnd
c        write(6,*) "virbnd",virbnd 
      endif
      
c     calculate valence angle forces
      
      if(ntangl.gt.0)then
        
        call angfrc
     x    (llsolva,lfree,lghost,idnode,imcon,mxnode,ntangl,engtmp,
     x    virtmp)
        
        engang=engang+engtmp*qfactor
        virang=virang+virtmp*qfactor
        
      endif
        
c     calculate dihedral forces
      
      if(ntdihd.gt.0)then
        
        call dihfrc
     x    (llsolva,lfree,lghost,idnode,imcon,mxnode,ntdihd,keyfce,
     x    dlrpot,epsq,engcpt,engtmp,engsrt,rcut,rvdw,alpha,vircpt,
     x    virtmp,virsrt)
        
        engcpe=engcpe+engcpt*qfactor
        vircpe=vircpe+vircpt*qfactor
        engsrp=engsrp+engsrt*qfactor
        virsrp=virsrp+virsrt*qfactor
        engdih=engdih+engtmp*qfactor
        virdih=virdih+virtmp*qfactor
        
      endif
      
c     calculate inversion forces
      
      if(ntinv.gt.0)then
        
        call invfrc
     x    (llsolva,lfree,lghost,idnode,imcon,mxnode,ntinv,engtmp,
     x    virtmp)
        
        enginv=enginv+engtmp*qfactor
        virinv=virinv+virtmp*qfactor
        
      endif
      
c     calculate tethered atom forces
      
      if(ntteth.gt.0)then
        
        call tethfrc
     x    (idnode,mxnode,imcon,natms,nstep,ntteth,engtmp,virtmp)
        
        engtet=engtet+engtmp*qfactor
        virtet=virtet+virtmp*qfactor
        
      endif
      
c     calculate shell model forces
      
      if(keyshl.gt.0)then
        
        call shlfrc
     x    (llsolva,lfree,lghost,idnode,imcon,mxnode,ntshl,engshl,virshl)
        
      endif
      
      return
      end subroutine intra_forces
      
      subroutine forces
     x  (loglnk,lgofr,lzeql,lsolva,lfree,lghost,idnode,imcon,keyfce,
     x  kmax1,kmax2,kmax3,nhko,nlatt,mxnode,ntpvdw,ntpmet,natms,
     x  nstbgr,nstep,nsteql,numrdf,nospl,nsolva,isolva,alpha,dlrpot,
     x  drewd,engcpe,engsrp,epsq,rcut,rvdw,vircpe,virsrp,volm,engmet,
     x  virmet)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating interatomic forces
c     using the verlet neighbour list
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     amended  - t. forester sept 1994
c     amended  - w. smith june 1995 for metal potentials
c     
c     key:
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     = 2,3  ------ Ewald sum                         : ewald1,2,3
c     = 4,5  ------ distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     = 8,9  ------ truncated and shifted coulombic   : coul4
c     = 10,11 ----- reaction field                    : coul3
c     = 12,13 ----- Smoothed Particle Mesh Ewald      : ewald[_spme,2,3]
c     = 14,15 ----- Hautman-Klein-Ewald               : hkewald1,2,3
c     
c***********************************************************************
      
      implicit none
      
      logical lgofr,lzeql,loglnk,lewald,lspme,lhke,newjob,lcshft,safe
      logical lsolva,lfree,lghost,llsolva
      
      integer idnode,imcon,keyfce,kmax1,kmax2,kmax3,nhko,nlatt
      integer mxnode,ntpvdw,natms,nstbgr,nstep,nsteql,numrdf
      integer ntpmet,nospl,nsolva,isolva,i,j,k,ii
      
      real(8) alpha,dlrpot,drewd,engcpe,engsrp,epsq,rcut,rvdw,eps
      real(8) vircpe,virsrp,volm,engacc,engac1,viracc,engmet,virmet
      
      save newjob
      
      data newjob/.true./
      
      safe=.true.
      llsolva=.false.
      if(lsolva)then
        llsolva=(mod(nstep-nsolva,isolva).eq.0)
      endif
      lhke=(keyfce/2.eq.7)
      lspme=(keyfce/2.eq.6)
      lewald=(keyfce/2.eq.1)
      lcshft=(keyfce/2.eq.4.or.keyfce/2.eq.5)
      
c     create ewald interpolation arrays
      
      if(newjob)then
        
        if(lhke)then
          
          call hkgen(idnode,nhko,nlatt,alpha,drewd,rcut)
          
        else if(lewald.or.lspme.or.lcshft)then
          
          call erfcgen(alpha,drewd,rcut)
          
        endif
        
        newjob=.false.
        
      endif
      
c     initialise force arrays
      
      do i=1,natms
        
        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0
        
      enddo
      
c     calculate local density in metals
      
      if(ntpmet.gt.0)then
        
        call metdens
     x    (idnode,imcon,mxnode,natms,engmet,virmet)
        
        stress(1)=stress(1)-virmet/3.d0
        stress(5)=stress(5)-virmet/3.d0
        stress(9)=stress(9)-virmet/3.d0
        
      endif
      
c     fourier contribution to coulombic forces in Ewald sum
      
      if(lewald)then
        
        call ewald1
     x    (lsolva,llsolva,lfree,lghost,idnode,mxnode,natms,imcon,
     x    kmax1,kmax2,kmax3,engac1,viracc,alpha,volm,epsq)
        
        engcpe=engcpe+engac1
        vircpe=vircpe+viracc
        
      endif
      
c     hautman-klein-ewald method
      
      if(lhke)then
        
c     fourier terms of hk-ewald
        
        call hkewald1
     x    (idnode,mxnode,natms,imcon,nhko,kmax1,kmax2,
     x    engacc,viracc,alpha,epsq)
        
        engcpe=engcpe+engacc
        vircpe=vircpe+viracc
        
c     real space terms of hk-ewald
        
        call hkewald2
     x    (idnode,mxnode,nhko,nlatt,imcon,natms,engacc,viracc,
     x    drewd,rcut,epsq)
        
        engcpe=engcpe+engacc
        vircpe=vircpe+viracc
        
      endif
      
c     smoothed particle mesh ewald
      
      if(lspme)then
        
        call ewald_spme
     x    (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,nospl,
     x    engac1,viracc,alpha,volm,epsq)
        
        engcpe=engcpe+engac1
        vircpe=vircpe+viracc
        
      endif
      
c     outer loop over atoms
      
      ii=0
      
      do i=idnode+1,natms,mxnode
        
        ii=ii+1
        
c     calculate interatomic distances
        
        do k=1,lentry(ii)
          
          j=list(ii,k)
          ilist(k)=j
          
          xdf(k)=xxx(i)-xxx(j)
          ydf(k)=yyy(i)-yyy(j)
          zdf(k)=zzz(i)-zzz(j)
          
        enddo
        
c     periodic boundary conditions
        
        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)
        
c     square of distances
        
        do k=1,lentry(ii)
          
          rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
          
        enddo
        
c     calculate metal forces and potential
        
        if(ntpmet.gt.0)then
          
          call metfrc(safe,i,lentry(ii),engacc,viracc)
          
          engmet=engmet+engacc
          virmet=virmet+viracc
          
        endif
        
c     calculate short range force and potential terms
        
        if(ntpvdw.gt.0.and.mod(keyfce,2).eq.1)then
          
          call srfrce
     x      (llsolva,lfree,lghost,i,lentry(ii),engacc,viracc,
     x      rvdw,dlrpot)
          
          engsrp=engsrp+engacc
          virsrp=virsrp+viracc
          
        endif
        
c     calculate coulombic force and potential terms
c     (real space contributions to ewald sum)
        
        if(lewald.or.lspme)then
          
          call ewald2(llsolva,lfree,lghost,i,lentry(ii),engacc,
     x      viracc,drewd,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.2)then
          
          call coul2
     x      (llsolva,lfree,lghost,i,lentry(ii),engacc,viracc,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.3)then
          
          call coul0
     x      (llsolva,lfree,lghost,i,lentry(ii),engacc,viracc,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.4)then
          
          call coul4
     x      (llsolva,lfree,lghost,i,lentry(ii),engacc,viracc,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.5)then
          
          call coul3
     x      (llsolva,lfree,lghost,i,lentry(ii),engacc,viracc,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        endif
        
c     accumulate radial distribution functions
        
        if(lgofr.and.((.not.lzeql).or.(nstep.gt.nsteql)))then
          
          if(mod(nstep,nstbgr).eq.0)then
            
            call rdf0(i,lentry(ii),rcut)
            
          endif
          
        endif
        
      enddo      
      
c     check metal interpolation is safe
      
      if(ntpmet.gt.0)then
        
        if(mxnode.gt.1)call gstate(safe)
        if(.not.safe)call error(idnode,142)
        
      endif
      
c     calculate corrections for intramolecular coulomb terms in
c     Ewald sum
      
      if(lewald.or.lspme.or.lhke)then
        
        eps=epsq
        if(loglnk)eps=eps*2.0d0
        
c     outer loop over atoms
        
        ii=0
        
        do i=idnode+1,natms,mxnode
          
          ii=ii+1
          
c     calculate interatomic distances
          
          do k=1,nexatm(ii)
            
            j=lexatm(ii,k)
            jlist(k)=j
            
            xdf(k)=xxx(i)-xxx(j)
            ydf(k)=yyy(i)-yyy(j)
            zdf(k)=zzz(i)-zzz(j)
            
          enddo
          
c     periodic boundary condition
          
          call images(imcon,0,1,nexatm(ii),cell,xdf,ydf,zdf)
          
c     calculate correction terms
          
          if(lhke)then
            
            call hkewald3(i,ii,engacc,viracc,eps)
            
          else
            
            call ewald3
     x        (llsolva,lfree,lghost,i,ii,engacc,viracc,alpha,eps)
            
          endif
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        enddo
        
      endif
      
c     counter for rdf statistics outside loop structure
      
      if(lgofr.and.((.not.lzeql).or.(nstep.gt.nsteql)))then
        
        if(mod(nstep,nstbgr).eq.0)then
          
          numrdf=numrdf+1
          
        endif
        
      endif
      
c     sum up contributions to short range and coulombic potential
      
      if(mxnode.gt.1)then
        
        buffer(1)=engsrp
        buffer(2)=virsrp
        buffer(3)=engcpe
        buffer(4)=vircpe
        buffer(5)=engmet
        buffer(6)=virmet
        buffer(7)=vdw_fre
        buffer(8)=cou_fre
        buffer(9)=vdw_vir
        buffer(10)=cou_vir
        call gdsum(buffer(1),10,buffer(11))
        engsrp=buffer(1)
        virsrp=buffer(2)
        engcpe=buffer(3)
        vircpe=buffer(4)
        engmet=buffer(5)
        virmet=buffer(6)
        vdw_fre=buffer(7)
        cou_fre=buffer(8)
        vdw_vir=buffer(9)
        cou_vir=buffer(10)
        
      endif

      return
      end subroutine forces
      
      subroutine forces_neu
     x  (lgofr,lzeql,lsolva,lfree,lghost,idnode,imcon,keyfce,
     x  mxnode,natms,nneut,nstbgr,nstep,nsteql,numrdf,nsolva,
     x  isolva,dlrpot,engcpe,engsrp,epsq,rcut,rvdw,alpha,
     x  vircpe,virsrp)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating interatomic forces
c     using the verlet neighbour list
c     neutral group implemenation - no Ewald sum option
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     
c     modified  - t. forester april 1993
c     key:
c     
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     = 2,3  ------ invalid
c     = 4,5  ------ distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     = 8,9  ------ invalid
c     
c***********************************************************************
      
      implicit none
      
      logical lgofr,lzeql,newlst,lchk,lsolva,lfree,lghost,llsolva
      
      integer idnode,imcon,keyfce,mxnode,natms,nneut,nstbgr
      integer nstep,nsteql,numrdf,i,fail,jneu,jj0,jj1,j
      integer ibig,ia,ineu,isn,ik,nsolva,isolva
      real(8) dlrpot,engcpe,engsrp,epsq,rcut,rvdw,vircpe
      real(8) virsrp,engacc,viracc,anorm,alpha
      
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      
      dimension fail(2)
      
      data fail/0,0/
      
c     allocate working arrays
      
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(1))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(2))
      do i=1,2
        if(fail(i).ne.0)call error(idnode,1820)
      enddo
      
      llsolva=.false.
      if(lsolva)then
        llsolva=(mod(nstep-nsolva,isolva).eq.0)
      endif
      
c     initialise force and stress arrays
      
      do i=1,natms
        
        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0
        
      enddo
      
      do i=1,9
        stress(i)=0.d0
      enddo
      
c     initialise energy and virial accumulators
      
      engcpe=0.d0
      engsrp=0.d0
      vircpe=0.d0
      virsrp=0.d0
      
c     intra group vectors com
      
      do jneu=1,nneut
        
        jj0=neulst(jneu)
        jj1=neulst(jneu+1)-1
        
c     loop over jneu sites
        
        do j=jj0,jj1
          
          txx(j)=xxx(j)-xxx(jj0)
          tyy(j)=yyy(j)-yyy(jj0)
          tzz(j)=zzz(j)-zzz(jj0)
          
        enddo
        
      enddo
      
      call images(imcon,0,1,natms,cell,txx,tyy,tzz)
      
      do jneu=1,nneut
        
        jj0=neulst(jneu)
        jj1=neulst(jneu+1)-1
        
c     loop over jneu sites
        
        do j=jj0,jj1
          
          xxx(j)=txx(j)+xxx(jj0)
          yyy(j)=tyy(j)+yyy(jj0)
          zzz(j)=tzz(j)+zzz(jj0)
          
        enddo
        
c     centre of molecule
        
        uxx(jneu)=0.d0
        uyy(jneu)=0.d0
        uzz(jneu)=0.d0
        anorm=1.d0/dble(jj1-jj0+1)
        
        do j=jj0,jj1
          
          uxx(jneu)=uxx(jneu)+xxx(j)*anorm
          uyy(jneu)=uyy(jneu)+yyy(j)*anorm
          uzz(jneu)=uzz(jneu)+zzz(j)*anorm
          
        enddo
        
c     vector from site to geometric centre
        
        do j=jj0,jj1
          
          txx(j)=xxx(j)-uxx(jneu)
          tyy(j)=yyy(j)-uyy(jneu)
          tzz(j)=zzz(j)-uzz(jneu)
          
        enddo
        
      enddo
      
c     outer loop over neutral groups
      
      lchk=.true.
      ibig=0
      ia=0
      
      do ineu=idnode+1,nneut,mxnode
        
        ia=ia+1
        
c     calculate interatomic distances
        
        newlst=.true.
        
        isn=1
        call neutlst
     x    (newlst,lchk,isn,imcon,idnode,ineu,ia,ik,
     x    txx,tyy,tzz,uxx,uyy,uzz)
        
c     trap possible array bound exception 
        
        ibig=max(ibig,ik)
        if(ik.gt.mxxdf)ik=0
        
c     calculate short range force and potential terms
        
        if(mod(keyfce,2).eq.1)then
          
          call srfrceneu
     x      (llsolva,lfree,lghost,ik,engacc,viracc,dlrpot,rvdw)
          
          engsrp=engsrp+engacc
          virsrp=virsrp+viracc
          
        endif
        
c     calculate coulombic force and potential terms
        
        if(keyfce/2.eq.2)then
          
          call coul2neu
     x      (llsolva,lfree,lghost,ik,engacc,viracc,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.3)then
          
          call coul0neu
     x      (llsolva,lfree,lghost,ik,engacc,viracc,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.4)then
          
          call error(idnode,250)
          
        elseif(keyfce/2.eq.5)then
          
          call coul3neu
     x      (llsolva,lfree,lghost,ik,engacc,viracc,epsq,rcut,alpha)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        endif
        
c     accumulate radial distribution functions
        
        if( ((.not.lzeql).or.(nstep.gt.nsteql)).and.(lgofr).and.
     x    mod(nstep,nstbgr).eq.0)then
          
          call rdf0neu(ik,rcut)
          
        endif
        
      enddo
      
c     check on validity of call to neutlst
      
      if(mxnode.gt.1)call gstate(lchk)
      if(.not.lchk)then 
        call gimax(ibig,1,i)
        if(idnode.eq.0)write(nrite,*)'mxxdf must be at least ',ibig
        if(idnode.eq.0)write(nrite,*)'mxxdf is currently ',mxxdf
        call  error(idnode,478)
      endif
      
      if(keyfce/2.eq.1.or.keyfce/2.eq.6)call error(idnode,250)
      
c     counter for rdf statistics outside loop structure
      
      if( ((.not.lzeql).or.(nstep.gt.nsteql)).and.(lgofr).and.
     x  mod(nstep,nstbgr).eq.0)numrdf=numrdf+1
      
c     sum up contributions to short range and coulombic potential
      
      if(mxnode.gt.1)then
        
        buffer(1)=engsrp
        buffer(2)=virsrp
        buffer(3)=engcpe
        buffer(4)=vircpe
        buffer(5)=vdw_fre
        buffer(6)=cou_fre
        buffer(7)=vdw_vir
        buffer(8)=cou_vir
        call gdsum(buffer(1),8,buffer(9))
        engsrp=buffer(1)
        virsrp=buffer(2)
        engcpe=buffer(3)
        vircpe=buffer(4)
        vdw_fre=buffer(5)
        cou_fre=buffer(6)
        vdw_vir=buffer(7)
        cou_vir=buffer(8)
        
      endif
      
c     deallocate work arrays
      
      deallocate (txx,tyy,tzz,stat=fail(1))
      deallocate (uxx,uyy,uzz,stat=fail(2))
      
      return
      end subroutine forces_neu
      
      subroutine multiple
     x  (loglnk,lgofr,lzeql,newlst,lsolva,lfree,lghost,idnode,
     x  imcon,keyfce,nlatt,kmax1,kmax2,kmax3,nhko,multt,
     x  mxnode,natms,nstep,nstbgr,nsteql,numrdf,nospl,nsolva,
     x  isolva,alpha,dlrpot,drewd,engcpe,engsrp,epsq,rcut,rprim,
     x  rvdw,vircpe,virsrp,volm)
c***************************************************************************
c     
c     dl_poly subroutine for multiple time step algorithm
c     reciprocal space calculated on long time steps.
c     
c     copyright daresbury laboratory
c     
c     author  t. forester,  may 1993
c     
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     = 2,3  ------ Ewald sum                         : ewald1,2,3,4
c     = 4,5  ------ distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     = 8,9  ------ truncated and shifted coulombic   : coul4
c     = 10,11 ----- reaction field                    : coul3
c     = 12,13 ----- Smoothed Particle Mesh Ewald      : ewald[_spme,2,3,4]
c     = 14,15 ----- Hautman-Klein-Ewald               : hkewald1,2,3,4
c     
c****************************************************************************
      
      implicit none
      
      integer, parameter :: nnn=5

      logical newplst,newlst,lgofr,lzeql,lgr,loglnk,lewald,lspme
      logical lhke,newjob,lcshft,lsolva,lfree,lghost,llsolva
      integer idnode,imcon,keyfce,nlatt,kmax1,kmax2,kmax3,nhko,multt
      integer mxnode,natms,nstep,nstbgr,nsteql,numrdf,nospl,fail
      integer numlsts,i,nstep0,nsolva,isolva,ii,k,j,ik
      real(8) alpha,dlrpot,drewd,engcpe,engsrp,epsq,rcut,rprim,rvdw
      real(8) vircpe,virsrp,volm,stresp,engcpl,engacc,viracc,engac1
      real(8) vircpl,eps,ann,engsr1,viracl,engsrl,virsrl,virac2,engcp1
      real(8) vircp1,engacl,engac2,virsr1
      
      real(8), allocatable :: fpx(:),fpy(:),fpz(:)
      real(8), allocatable :: vdw_sol_put(:),cou_sol_put(:)
      real(8), allocatable :: vdw_exc_put(:),cou_exc_put(:)
      
      dimension stresp(9),fail(nnn)
      
      save engcpl,engsrl,vircpl,virsrl,nstep0,numlsts,engcp1,vircp1
      save engsr1,virsr1,stresp,fpx,fpy,fpz,newjob
      save vdw_sol_put,cou_sol_put,vdw_exc_put,cou_exc_put
      
      data newjob/.true./
      data numlsts/-1/
      
      llsolva=.false.
      if(lsolva)then
        llsolva=(mod(nstep-nsolva,isolva).eq.0)
      endif
      lhke=(keyfce/2.eq.7)
      lspme=(keyfce/2.eq.6)
      lewald=(keyfce/2.eq.1)
      lcshft=(keyfce/2.eq.4.or.keyfce/2.eq.5)
      if(newlst)nstep0=nstep
      newplst=(newlst).or.(mod(nstep-nstep0,multt).eq.0)
      
c     allocate working arrays
      
      if(newjob)then
        
        do i=1,nnn
          fail(i)=0
        enddo
        allocate (fpx(mxatms),fpy(mxatms),fpz(mxatms),stat=fail(1))
        if(lsolva)then
          
          allocate (vdw_sol_put(mxtmls_sol2),stat=fail(2))
          allocate (cou_sol_put(mxtmls_sol2),stat=fail(3))
          if(lghost)then
            allocate (vdw_exc_put(mxtmls_exc2),stat=fail(4))
            allocate (cou_exc_put(mxtmls_exc2),stat=fail(5))
          endif
          
        endif
        do i=1,nnn
          if(fail(i).ne.0)call error(idnode,1840)
        enddo
        
      endif
      
c     create ewald interpolation arrays
      
      if(newjob)then
        
        if(lspme.or.lewald.or.lcshft)then
          
          call erfcgen(alpha,drewd,rcut)
          
        endif
        
      endif
      
      newjob=.false.
      
c     divide neighbour list into primary and secondary neighbours
      
      if(newplst)then        
        
        numlsts=numlsts+1
        call primlst(idnode,mxnode,natms,imcon,rprim)
        
      endif
      
c     flag for accumulating rdfs
      
      lgr=.false.
      if(nstbgr.gt.0)lgr=(mod(numlsts,nstbgr).eq.0)
      lgr=(lgr.and.(newplst.and.lgofr))
      lgr=(lgr.and.((.not.lzeql).or.(nstep-nsteql.gt.0)))
      
c     zero force and stress arrays
      
      do i=1,natms
        
        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0
        
      enddo
      
      do i=1,9
        stress(i)=0.d0
      enddo
      
c     ********************PROCESS SECONDARY NEIGHBOURS******************
      
      if(newplst.or.(mod(nstep-nstep0,multt).le.1))then
        
c     zero accumulators
        
        engcpl=0.d0
        vircpl=0.d0
        engsrl=0.d0
        virsrl=0.d0
        llsolva=lsolva
        if(lsolva)then
          
          vdw_sol(:)=0.d0
          cou_sol(:)=0.d0
          
          if(lghost)then
            
            vdw_exc(:)=0.d0
            cou_exc(:)=0.d0
            
          endif
          
        endif
        
c     calculate fourier contribution to secondary coulombic forces
        
        if(lewald.or.lspme.or.lhke)then
          
          if(lewald)then
            
            call ewald1
     x        (lsolva,llsolva,lfree,lghost,idnode,mxnode,natms,imcon,
     x        kmax1,kmax2,kmax3,engac1,viracc,alpha,volm,epsq)
            
c     hautman-klein-ewald method
            
          elseif(lhke)then
            
            call hkewald1
     x        (idnode,mxnode,natms,imcon,nhko,kmax1,kmax2,
     x        engac1,viracc,alpha,epsq)
            
c     real space terms of hk-ewald
            
            call hkewald2
     x        (idnode,mxnode,nhko,nlatt,imcon,natms,engac2,
     x        virac2,drewd,rcut,epsq)
            
            engac1=engac1+engac2
            viracc=viracc+virac2
            
          elseif(lspme)then
            
c     smoothed particle mesh ewald
            
            call ewald_spme
     x        (idnode,mxnode,natms,imcon,kmax1,kmax2,kmax3,nospl,
     x        engac1,viracc,alpha,volm,epsq)
            
          endif
          
          engcpl=engcpl+engac1
          vircpl=vircpl+viracc
          
c     calculate corrections for intramolecular coulomb terms in 
c     Ewald sum
c     note: if using link cells - have double counted excluded 
c     interactions use temporary adjustment to relative dielectric
c     constant
          
          eps=epsq
          if(loglnk)eps=epsq*2.0d0
          
c     calculate self interaction corrections for fourier contributions
          
          ii=0
          
          do i=idnode+1,natms,mxnode
            
            ii=ii+1
            
c     calculate interatomic distances
            
            do k=1,nexatm(ii)
              
              j=lexatm(ii,k)
              jlist(k)=j
              
              xdf(k)=xxx(i)-xxx(j)
              ydf(k)=yyy(i)-yyy(j)
              zdf(k)=zzz(i)-zzz(j)
              
            enddo
            
c     periodic boundary condition
            
            call images(imcon,0,1,nexatm(ii),cell,xdf,ydf,zdf)
            
c     calculate correction terms
            
            if(lhke)then
              
              call hkewald3(i,ii,engacc,viracc,eps)
              
            else
              
              call ewald3
     x          (lsolva,lfree,lghost,i,ii,engacc,viracc,alpha,eps)
              
            endif
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          enddo
          
        endif
        
c     calculate pair contributions to secondary neighbour forces
        
        ii=0
        
        do i=idnode+1,natms,mxnode
          
          ii=ii+1
          
c     calculate interatomic distances
          
          ik=0
          
          do k=1,lentry(ii)
            
            j=list(ii,k)
            
            if(j.gt.0)then
              
              ik=ik+1
              ilist(ik)=j
              xdf(ik)=xxx(i)-xxx(j)
              ydf(ik)=yyy(i)-yyy(j)
              zdf(ik)=zzz(i)-zzz(j)
              
            endif
            
          enddo
          
c     periodic boundary conditions
          
          call images(imcon,0,1,ik,cell,xdf,ydf,zdf)
          
c     square of distance
          
          do k=1,ik
            
            rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
            
          enddo
          
c     accumulate radial distribution functions
          
          if(lgr)call rdf0(i,ik,rcut)
          
c     calculate short range force and potential terms
          
          if(mod(keyfce,2).eq.1)then
            
            call srfrce
     x        (lsolva,lfree,lghost,i,ik,engacc,viracc,rvdw,dlrpot)
            
            engsrl=engsrl+engacc
            virsrl=virsrl+viracc
            
          endif
          
c     calculate coulombic force and potential terms
c     (real space contributions to ewald sum)
          
          if(lewald.or.lspme)then
            
            call ewald2
     x        (lsolva,lfree,lghost,i,ik,engacc,viracc,drewd,rcut,epsq)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif(keyfce/2.eq.2)then
            
            call coul2
     x        (lsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif(keyfce/2.eq.3)then
            
            call coul0
     x        (lsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif(keyfce/2.eq.4)then
            
            call coul4
     x        (lsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif(keyfce/2.eq.5)then
            
            call coul3
     x        (lsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          endif
          
        enddo
        
c     store secondary forces and stress tensor
        
        do i=1,natms
          
          flx(i)=fxx(i)
          fly(i)=fyy(i)
          flz(i)=fzz(i)
          fxx(i)=0.d0
          fyy(i)=0.d0
          fzz(i)=0.d0
          
        enddo
        
        do i=1,9
          stresl(i)=stress(i)
          stress(i)=0.d0
        enddo
        
c     store solvation and excitation arrays
        
        if(lsolva)then
          
          vdw_sol_lng(:)=vdw_sol(:)
          cou_sol_lng(:)=cou_sol(:)
          vdw_sol(:)=0.d0
          cou_sol(:)=0.d0
          
          if(lghost)then
            
            vdw_exc_lng(:)=vdw_exc(:)
            cou_exc_lng(:)=cou_exc(:)
            vdw_exc(:)=0.d0
            cou_exc(:)=0.d0
            
          endif
          
        endif
        
      endif
      
c     ****************END OF SECONDARY NEIGHBOUR PROCESSING*************
      
c     ********************PROCESS PRIMARY NEIGHBOURS********************
      
c     zero accumulators for total energies and virials
      
      engcpe=0.d0
      engsrp=0.d0
      vircpe=0.d0
      virsrp=0.d0
      
c     calculate pair force contributions
      
      ii=0
      
      do i=idnode+1,natms,mxnode
        
        ii=ii+1
        
c     calculate interatomic distances
        
        ik=0
        
        do k=1,lentry(ii)
          
          j=-list(ii,k)
          
          if(j.gt.0)then
            
            ik=ik+1
            ilist(ik)=j
            xdf(ik)=xxx(i)-xxx(j)
            ydf(ik)=yyy(i)-yyy(j)
            zdf(ik)=zzz(i)-zzz(j)
            
          endif
          
        enddo
        
c     periodic boundary conditions
        
        call images(imcon,0,1,ik,cell,xdf,ydf,zdf)
        
c     square of distance
        
        do k=1,ik
          
          rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
          
        enddo
        
c     accumulate radial distribution functions
        
        if(lgr)call rdf0(i,ik,rcut)
        
c     calculate short range force and potential terms
        
        if(mod(keyfce,2).eq.1)then
          
          call srfrce
     x      (llsolva,lfree,lghost,i,ik,engacc,viracc,rvdw,dlrpot)
          
          engsrp=engsrp+engacc
          virsrp=virsrp+viracc
          
        endif
        
c     calculate coulombic force and potential terms
c     (real space contributions to ewald sum)
        
        if(lewald.or.lspme.or.lhke)then
          
          if(newplst.or.
     x      (mod(nstep-nstep0,multt).le.1))then
            
            if(lhke)then
              
              call hkewald4(i,ik,engacc,viracc,engacl,viracl,rcut,epsq)
              
            else
              
              call ewald4
     x          (llsolva,lfree,lghost,i,ik,engacc,viracc,engacl,viracl,
     x          drewd,rcut,epsq)
              
            endif
            
            engcpe=engcpe+engacc
            vircpe=vircpe+viracc
            engcpl=engcpl+engacl
            vircpl=vircpl+viracl
            
          else
            
            call coul0
     x        (llsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
            
            engcpe=engcpe+engacc
            vircpe=vircpe+viracc
            
          endif
          
        elseif(keyfce/2.eq.2)then
          
          call coul2
     x      (llsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.3)then
          
          call coul0
     x      (llsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.4)then
          
          call coul4
     x      (llsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.5)then
          
          call coul3
     x      (llsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        endif
        
      enddo
      
c     **************END OF PRIMARY NEIGHBOUR PROCESSING*****************
      
c     counter for rdf statistics outside loop structure
      
      if(lgr)numrdf=numrdf+1
      
c     add in secondary neighbour contributions to force, energy etc
      
      if(newplst)then
        
        do i=1,natms
          
          fpx(i)=flx(i)
          fpy(i)=fly(i)
          fpz(i)=flz(i)
          
        enddo
        
        do i=1,9
          stresp(i)=stresl(i)
        enddo
        
        engsr1=engsrl
        virsr1=virsrl
        engcp1=engcpl
        vircp1=vircpl
        
c     store solvation and excitation arrays
        
        if(lsolva)then
          
          vdw_sol_put(:)=vdw_sol_lng(:)
          cou_sol_put(:)=cou_sol_lng(:)
          
          if(lghost)then
            
            vdw_exc_put(:)=vdw_exc_lng(:)
            cou_exc_put(:)=cou_exc_lng(:)
            
          endif
          
        endif
        
      endif
      
c     calculate force increments etc
      
      if(mod(nstep-nstep0,multt).eq.1)then
        
        do i=1,natms
          
          flx(i)=flx(i)-fpx(i)
          fly(i)=fly(i)-fpy(i)
          flz(i)=flz(i)-fpz(i)
          
        enddo
        
        do i=1,9
          stresl(i)=stresl(i)-stresp(i)
        enddo
        
        virsrl=virsrl-virsr1
        engsrl=engsrl-engsr1
        vircpl=vircpl-vircp1
        engcpl=engcpl-engcp1
        
c     solvation and excitation increments
        
        if(lsolva)then
          
          vdw_sol_lng(:)=vdw_sol_lng(:)-vdw_sol_put(:)
          cou_sol_lng(:)=cou_sol_lng(:)-cou_sol_put(:)
          
          if(lghost)then
            
            vdw_exc_lng(:)=vdw_exc_lng(:)-vdw_exc_put(:)
            cou_exc_lng(:)=cou_exc_lng(:)-cou_exc_put(:)
            
          endif
          
        endif
        
      endif
      
c     extrapolate long range terms
      
      ann=dble(mod(nstep-nstep0,multt))
      
      do i=1,natms
        
        fxx(i)=fpx(i)+flx(i)*ann+fxx(i)
        fyy(i)=fpy(i)+fly(i)*ann+fyy(i)
        fzz(i)=fpz(i)+flz(i)*ann+fzz(i)
        
      enddo
      
      do i=1,9
        stress(i)=stress(i)+stresp(i)+stresl(i)*ann
      enddo
      
      engsrp=engsr1+engsrl*ann+engsrp
      virsrp=virsr1+virsrl*ann+virsrp
      engcpe=engcp1+engcpl*ann+engcpe
      vircpe=vircp1+vircpl*ann+vircpe
      
c     solvation and excitation extrapolation
      
      if(llsolva)then
        
        vdw_sol(:)=vdw_sol_put(:)+vdw_sol_lng(:)*ann+vdw_sol(:)
        cou_sol(:)=cou_sol_put(:)+cou_sol_lng(:)*ann+cou_sol(:)
        
        if(lghost)then
          
          vdw_exc(:)=vdw_exc_put(:)+vdw_exc_lng(:)*ann+vdw_exc(:)
          cou_exc(:)=cou_exc_put(:)+cou_exc_lng(:)*ann+cou_exc(:)
          
        endif
        
      endif
      
c     sum up contributions to short range and coulombic potential
      
      if(mxnode.gt.1)then
        
        buffer(1)=engsrp
        buffer(2)=virsrp
        buffer(3)=engcpe
        buffer(4)=vircpe
        buffer(5)=vdw_fre
        buffer(6)=cou_fre
        buffer(7)=vdw_vir
        buffer(8)=cou_vir
        call gdsum(buffer(1),8,buffer(9))
        engsrp=buffer(1)
        virsrp=buffer(2)
        engcpe=buffer(3)
        vircpe=buffer(4)
        vdw_fre=buffer(5)
        cou_fre=buffer(6)
        vdw_vir=buffer(7)
        cou_vir=buffer(8)
        
      endif
      
      return
      end subroutine multiple
      
      subroutine multiple_neu
     x  (lgofr,lzeql,newlst,lsolva,lfree,lghost,idnode,imcon,
     x  keyfce,multt,mxnode,natms,nneut,nstbgr,nstep,nsteql,
     x  numrdf,nsolva,isolva,delr,dlrpot,engcpe,engsrp,epsq,
     x  rprim,rcut,rvdw,alpha,vircpe,virsrp)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating interatomic forces
c     using the verlet neighbour list
c     neutral group implemenation - no Ewald sum option
c     multiple timestep option
c     
c     parallel replicated data version
c     
c     fpx,fpy,fpz : forces from electrostatics fron rprim < r <= rcut
c     fxx,fyy,fzz : total force
c     
c     copyright daresbury laboratory april 1994
c     author  - t. forester april 1993
c     key:
c     
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     = 2,3  ------ invalid
c     = 4,5  ------ distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     = 8,9  ------ invalid
c     = 10,11 ----- reaction field                    : coul3
c     
c***********************************************************************
      
      implicit none
      
      integer, parameter :: nnn=7

      logical lgofr,lzeql,newlst,newplst,lgr,lchk,newjob,lsolva
      logical lfree,lghost,llsolva
      integer idnode,imcon,keyfce,multt,mxnode,natms,nneut,nstbgr
      integer nstep,nsteql,numrdf,fail,i,numlsts,jneu,jj0,j
      integer jj1,ineu,ia,isn,ibig,ik,nstep0,nsolva,isolva
      real(8) delr,dlrpot,engcpe,engsrp,epsq,rprim,rcut,rvdw,vircpe
      real(8) virsrp,engcpl,vircpl,engsrl,virsrl,anorm,ann,stresp
      real(8) engacc,viracc,engsr1,virsr1,engcp1,vircp1,alpha
      
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: uxx(:),uyy(:),uzz(:)
      real(8), allocatable :: fpx(:),fpy(:),fpz(:)
      real(8), allocatable :: vdw_sol_put(:),cou_sol_put(:)
      real(8), allocatable :: vdw_exc_put(:),cou_exc_put(:)
      
      dimension fail(nnn),stresp(9)
      
      save engcpl,engsrl,vircpl,virsrl,nstep0,numlsts,engcp1,vircp1
      save engsr1,virsr1,stresp,fpx,fpy,fpz
      save vdw_sol_put,cou_sol_put,vdw_exc_put,cou_exc_put
      
      data newjob/.true./
      data numlsts/-1/
      
c     allocate working arrays
      
      do i=1,nnn
        fail(i)=0
      enddo
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(1))
      allocate (uxx(mxatms),uyy(mxatms),uzz(mxatms),stat=fail(2))
      if(newjob)then
        
        allocate (fpx(mxatms),fpy(mxatms),fpz(mxatms),stat=fail(3))
        if(lsolva)then
          
          allocate (vdw_sol_put(mxtmls_sol2),stat=fail(4))
          allocate (cou_sol_put(mxtmls_sol2),stat=fail(5))
          
          if(lghost)then
            allocate (vdw_exc_put(mxtmls_exc2),stat=fail(6))
            allocate (cou_exc_put(mxtmls_exc2),stat=fail(7))
          endif
          
        endif
        newjob=.false.
        
      endif
      do i=1,nnn
        if(fail(i).ne.0)call error(idnode,1850)
      enddo
      
      llsolva=.false.
      if(lsolva)then
        llsolva=(mod(nstep-nsolva,isolva).eq.0)
      endif

c     error if ewald sum requested
      
      if(keyfce/2.eq.1.or.keyfce/2.eq.6)call error(idnode,250)
      
c     create list of primary and secondary neighbours
      
      if(newlst)nstep0=nstep
      newplst=(mod(nstep-nstep0,multt).eq.0)
      
      if(newplst)then        
        
        numlsts=numlsts+1
        call prneulst(newlst,imcon,idnode,mxnode,nneut,rprim)
        
      endif
      
c     zero accumulators for total energies and virials
      
      engcpe=0.d0
      engsrp=0.d0
      vircpe=0.d0
      virsrp=0.d0
      
c     zero force arrays
      
      do i=1,natms
        
        fxx(i)=0.d0
        fyy(i)=0.d0
        fzz(i)=0.d0
        
      enddo
      
c     zero stress arrays
      
      do i=1,9
        stress(i)=0.d0
      enddo
      
c     flag for accumulating rdfs
      
      lgr=.false.
      if(nstbgr.gt.0)lgr=(mod(numlsts,nstbgr).eq.0)
      lgr=(lgr.and.(newplst.and.lgofr))
      lgr=(lgr.and.((.not.lzeql).or.(nstep-nsteql.gt.0)))
      
c     intra group vectors com
      
      do jneu=1,nneut
        
        jj0=neulst(jneu)
        jj1=neulst(jneu+1)-1
        
c     loop over jneu sites
        
        do j=jj0,jj1
          
          txx(j)=xxx(j)-xxx(jj0)
          tyy(j)=yyy(j)-yyy(jj0)
          tzz(j)=zzz(j)-zzz(jj0)
          
        enddo
        
      enddo
      
      call images(imcon,0,1,natms,cell,txx,tyy,tzz)
      
      do jneu=1,nneut
        
        jj0=neulst(jneu)
        jj1=neulst(jneu+1)-1
        
c     loop over jneu sites
        
        do j=jj0,jj1
          
          xxx(j)=txx(j)+xxx(jj0)
          yyy(j)=tyy(j)+yyy(jj0)
          zzz(j)=tzz(j)+zzz(jj0)
          
        enddo
        
c     centre of molecule
        
        uxx(jneu)=0.d0
        uyy(jneu)=0.d0
        uzz(jneu)=0.d0
        anorm=1.d0/dble(jj1-jj0+1)
        
        do j=jj0,jj1
          
          uxx(jneu)=uxx(jneu)+xxx(j)*anorm
          uyy(jneu)=uyy(jneu)+yyy(j)*anorm
          uzz(jneu)=uzz(jneu)+zzz(j)*anorm
          
        enddo
        
c     vector from site to geometric centre
        
        do j=jj0,jj1
          
          txx(j)=xxx(j)-uxx(jneu)
          tyy(j)=yyy(j)-uyy(jneu)
          tzz(j)=zzz(j)-uzz(jneu)
          
        enddo
        
      enddo
      
c     ********************PROCESS SECONDARY NEIGHBOURS********************
        
      lchk=.true.
      ibig=0
      ia=0
      
      if(newplst.or.(mod(nstep-nstep0,multt).le.1))then
      
c     zero accumulators for secondary neighbour energies and virial
        
        engcpl=0.d0
        vircpl=0.d0
        engsrl=0.d0
        virsrl=0.d0
        
c     initialise solvation and excitation  arrays

        if(lsolva)then
          
          cou_sol(:)=0.d0
          vdw_sol(:)=0.d0
          
          if(lghost)then
            
            cou_exc(:)=0.d0
            vdw_exc(:)=0.d0
            
          endif
          
        endif
        
c     outer loop over neutral groups
        
        do ineu=idnode+1,nneut,mxnode
          
          ia=ia+1
          
c     calculate interatomic distances
          
          isn=-1
          call neutlst
     x      (.true.,lchk,isn,imcon,idnode,ineu,ia,ik,
     x      txx,tyy,tzz,uxx,uyy,uzz)
          
c     trap possible array bound exception 
          
          ibig=max(ibig,ik)
          if(ik.gt.mxxdf)ik=0
          
c     calculate short range force and potential terms
          
          if(mod(keyfce,2).eq.1.and.(rvdw.gt.rprim-delr))then
            
            call srfrceneu
     x        (lsolva,lfree,lghost,ik,engacc,viracc,dlrpot,rvdw)
            
            engsrl=engsrl+engacc
            virsrl=virsrl+viracc
            
          endif
          
c     calculate coulombic force and potential terms
          
          if(keyfce/2.eq.2)then
            
            call coul2neu
     x        (lsolva,lfree,lghost,ik,engacc,viracc,epsq)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif(keyfce/2.eq.3)then
            
            call coul0neu
     x        (lsolva,lfree,lghost,ik,engacc,viracc,epsq)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          elseif(keyfce/2.eq.4)then
            
            call error(idnode,250)
            
          elseif(keyfce/2.eq.5)then
            
            call coul3neu
     x        (lsolva,lfree,lghost,ik,engacc,viracc,epsq,rcut,alpha)
            
            engcpl=engcpl+engacc
            vircpl=vircpl+viracc
            
          else
            
            call error(idnode,250)
            
          endif
          
c     accumulate radial distribution functions out to rcut
          
          if(lgr)then
            call rdf0neu(ik,rcut)
          endif
          
        enddo
        
c     store secondary forces and stress tensor
        
        do i=1,natms
          
          flx(i)=fxx(i)
          fly(i)=fyy(i)
          flz(i)=fzz(i)
          fxx(i)=0.d0
          fyy(i)=0.d0
          fzz(i)=0.d0
          
        enddo
        
        do i=1,9
          
          stresl(i)=stress(i)
          stress(i)=0.d0
          
        enddo
        
c     store solvation and excitation arrays
        
        if(lsolva)then
          
          vdw_sol_lng(:)=vdw_sol(:)
          cou_sol_lng(:)=cou_sol(:)
          vdw_sol(:)=0.d0
          cou_sol(:)=0.d0
          
          if(lghost)then
            
            vdw_exc_lng(:)=vdw_exc(:)
            cou_exc_lng(:)=cou_exc(:)
            vdw_exc(:)=0.d0
            cou_exc(:)=0.d0
            
          endif
          
        endif
        
      endif
      
c     ****************END OF SECONDARY NEIGHBOUR PROCESSING*************
      
c     ********************PROCESS PRIMARY NEIGHBOURS********************
      
      ia=0
      do ineu=idnode+1,nneut,mxnode
        
        ia=ia+1
        
c     calculate interatomic distances
        
        isn=1        
        call neutlst
     x    (.true.,lchk,isn,imcon,idnode,ineu,ia,ik,
     x    txx,tyy,tzz,uxx,uyy,uzz)
        
c     trap possible array bound exception 
        
        ibig=max(ibig,ik)
        if(ik.gt.mxxdf)ik=0
        
c     calculate short range force and potential terms
        
        if(mod(keyfce,2).eq.1)then
          
          call srfrceneu
     x      (llsolva,lfree,lghost,ik,engacc,viracc,dlrpot,rvdw)
          
          engsrp=engsrp+engacc
          virsrp=virsrp+viracc
          
        endif
        
c     calculate coulombic force and potential terms
        
        if(keyfce/2.eq.2)then
          
          call coul2neu
     x      (llsolva,lfree,lghost,ik,engacc,viracc,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.3)then
          
          call coul0neu
     x      (llsolva,lfree,lghost,ik,engacc,viracc,epsq)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        elseif(keyfce/2.eq.4)then
          
          call error(idnode,250)
          
        elseif(keyfce/2.eq.5)then
          
          call coul3neu
     x      (llsolva,lfree,lghost,ik,engacc,viracc,epsq,rcut,alpha)
          
          engcpe=engcpe+engacc
          vircpe=vircpe+viracc
          
        endif
        
c     accumulate radial distribution functions out to rcut
        
        if(lgr)then
          call rdf0neu(ik,rcut)
        endif
        
      enddo
      
c     ******************END OF PRIMARY NEIGHBOUR PROCESSING*************
      
c     check on validity of call to neutlst
      
      if(mxnode.gt.1)call gstate(lchk)
      if(.not.lchk)then 
        call gimax(ibig,1,i)
        if(idnode.eq.0)write(nrite,*)'mxxdf must be at least ',ibig
        if(idnode.eq.0)write(nrite,*)'mxxdf is currently ',mxxdf
        call  error(idnode,479)
      endif
      
c     counter for rdf statistics outside loop structure
      
      if(lgr)numrdf=numrdf+1
      
c     add in secondary neighbour contributions to force, energy etc
      
      if(newplst)then
        
        do i=1,natms
          
          fpx(i)=flx(i)
          fpy(i)=fly(i)
          fpz(i)=flz(i)
          
        enddo
        
        do i=1,9
          stresp(i)=stresl(i)
        enddo
        
        engsr1=engsrl
        virsr1=virsrl
        engcp1=engcpl
        vircp1=vircpl
        
c     store solvation and excitation arrays
        
        if(lsolva)then
          
          vdw_sol_put(:)=vdw_sol_lng(:)
          cou_sol_put(:)=cou_sol_lng(:)
          
          if(lghost)then
            
            vdw_exc_put(:)=vdw_exc_lng(:)
            cou_exc_put(:)=cou_exc_lng(:)
            
          endif
          
        endif
        
      endif
      
c     calculate force increments etc
      
      if(mod(nstep-nstep0,multt).eq.1)then
        
        do i=1,natms
          
          flx(i)=flx(i)-fpx(i)
          fly(i)=fly(i)-fpy(i)
          flz(i)=flz(i)-fpz(i)
          
        enddo
        
        do i=1,9
          stresl(i)=stresl(i)-stresp(i)
        enddo
        
        virsrl=virsrl-virsr1
        engsrl=engsrl-engsr1
        vircpl=vircpl-vircp1
        engcpl=engcpl-engcp1
        
c     solvation and excitation increments
        
        if(lsolva)then
          
          vdw_sol_lng(:)=vdw_sol_lng(:)-vdw_sol_put(:)
          cou_sol_lng(:)=cou_sol_lng(:)-cou_sol_put(:)
          
          if(lghost)then
            
            vdw_exc_lng(:)=vdw_exc_lng(:)-vdw_exc_put(:)
            cou_exc_lng(:)=cou_exc_lng(:)-cou_exc_put(:)
            
          endif
          
        endif
        
      endif
      
c     extrapolate long range terms
      
      ann=dble(mod(nstep-nstep0,multt))
      
      do i=1,natms
        
        fxx(i)=fpx(i)+flx(i)*ann+fxx(i)
        fyy(i)=fpy(i)+fly(i)*ann+fyy(i)
        fzz(i)=fpz(i)+flz(i)*ann+fzz(i)
        
      enddo
      
      do i=1,9
        stress(i)=stress(i)+stresl(i)*ann+stresp(i)
      enddo
      
      engsrp=engsr1+engsrl*ann+engsrp
      virsrp=virsr1+virsrl*ann+virsrp
      engcpe=engcp1+engcpl*ann+engcpe
      vircpe=vircp1+vircpl*ann+vircpe
      
c     solvation and excitation extrapolation
      
      if(llsolva)then
        
        vdw_sol(:)=vdw_sol_put(:)+vdw_sol_lng(:)*ann+vdw_sol(:)
        cou_sol(:)=cou_sol_put(:)+cou_sol_lng(:)*ann+cou_sol(:)
        
        if(lghost)then
          
          vdw_exc(:)=vdw_exc_put(:)+vdw_exc_lng(:)*ann+vdw_exc(:)
          cou_exc(:)=cou_exc_put(:)+cou_exc_lng(:)*ann+cou_exc(:)
          
        endif
        
      endif
      
c     sum up contributions to short range and coulombic potential
      
      if(mxnode.gt.1)then
        
        buffer(1)=engsrp
        buffer(2)=virsrp
        buffer(3)=engcpe
        buffer(4)=vircpe
        buffer(5)=vdw_fre
        buffer(6)=cou_fre
        buffer(7)=vdw_vir
        buffer(8)=cou_vir
        call gdsum(buffer(1),8,buffer(9))
        engsrp=buffer(1)
        virsrp=buffer(2)
        engcpe=buffer(3)
        vircpe=buffer(4)
        vdw_fre=buffer(5)
        cou_fre=buffer(6)
        vdw_vir=buffer(7)
        cou_vir=buffer(8)
        
      endif
      
c     deallocate work arrays
      
      deallocate (txx,tyy,tzz,uxx,uyy,uzz,stat=fail(1))
      
      return
      end subroutine multiple_neu
      
      subroutine multiple_nsq
     x  (lnsq,lgofr,lzeql,newlst,lsolva,lfree,lghost,idnode,
     x  imcon,keyfce,multt,mxnode,natms,nstep,nstbgr,nsteql,
     x  numrdf,nsolva,isolva,delr,dlrpot,engcpe,engsrp,epsq,
     x  rcut,rprim,rvdw,vircpe,virsrp)
      
c***********************************************************************
c     
c     dl_poly subroutine for multiple time step algorithm 
c     to be used with all-pairs option
c     
c     flx,fly,flz : forces from electrostatics from r > rcut
c     fpx,fpy,fpz : forces from electrostatics from rprim < r <= rcut
c     fxx,fyy,fzz : total force
c     
c     copyright daresbury laboratory 1993
c     
c     author  t. forester,  may 1993
c     
c     keyfce = odd  ------ short range potentials calculated : srfrce
c     = 0,1  ------ no electrostatics
c     Ewald sum --- not used
c     = 4,5  ------ Distance dependent dielectric     : coul2
c     = 6,7  ------ coulombic                         : coul0
c     truncated and shifted coulombic -- not used
c     reaction field - not used
c     
c****************************************************************************
      
      implicit none
      
      integer, parameter :: nnn=5
      logical newplst,newlst,lgofr,lzeql,lgr,lnsq,newjob,lsolva
      logical lfree,lghost,llsolva
      integer idnode,imcon,keyfce,multt,mxnode,natms,nstep,nstbgr
      integer nsteql,numrdf,fail,nstep0,ii,ik,k,numlsts,nsolva
      integer isolva,i,j
      real(8) delr,dlrpot,engcpe,engsrp,engcp3,epsq,rcut
      real(8) rprim,rvdw,vircpe,virsrp,vircp3,rcut1,engcp2,vircp2
      real(8) engsr2,virsr2,stresp,engacc,viracc
      
      real(8), allocatable :: fpx(:),fpy(:),fpz(:)
      real(8), allocatable :: vdw_sol_put(:),cou_sol_put(:)
      real(8), allocatable :: vdw_exc_put(:),cou_exc_put(:)
      
      dimension stresp(9),fail(nnn)
      
      save engsr2,virsr2,engcp2,vircp2,nstep0,numlsts,stresp,fpx,fpy,fpz
      save vdw_sol_put,cou_sol_put,vdw_exc_put,cou_exc_put

      data numlsts/-1/
      data newjob/.true./
      
c     allocate work arrays
      
      if(newjob)then
        
        do i=1,nnn
          fail(i)=0
        enddo
        allocate (fpx(mxatms),fpy(mxatms),fpz(mxatms),stat=fail(1))
        if(lsolva)then
          
          allocate (vdw_sol_put(mxtmls_sol2),stat=fail(2))
          allocate (cou_sol_put(mxtmls_sol2),stat=fail(3))
          if(lghost)then
            allocate (vdw_exc_put(mxtmls_exc2),stat=fail(4))
            allocate (cou_exc_put(mxtmls_exc2),stat=fail(5))
          endif
          
        endif
        do i=1,nnn
          if(fail(i).ne.0)call error(idnode,1860)
        enddo
        newjob=.false.
        
      endif
      
      if(lnsq)then
        
        llsolva=.false.
        if(lsolva)then
          llsolva=(mod(nstep-nsolva,isolva).eq.0)
        endif
        
c     divide neighbour list into primary and secondary neighbours
        
        if(newplst)then        
          
          numlsts=numlsts+1
          call primlst(idnode,mxnode,natms,imcon,rprim)
          
        endif
        
c     flag for accumulating rdfs
        
        lgr=(lgofr.and.(.not.lzeql.or.(nstep-nsteql.gt.0)))
        lgr=(lgr.and.newplst.and.(mod(numlsts,nstbgr).eq.0))
        
c     set extended cutoff for electrostatics - secondary shell
        
        rcut1=rcut+delr
        
        if(newlst)nstep0=nstep
        newplst=(newlst.or.mod(nstep-nstep0,multt).eq.0)
        
c     ********************PROCESS TERTIARY NEIGHBOURS*********************
        
        if(newplst)then
          
          call coul_nsq
     x      (lsolva,lfree,lghost,idnode,mxnode,natms,imcon,epsq,rcut,
     x      engcp3,vircp3)
          
        endif
        
c     ****************END OF TERTIARY NEIGHBOUR PROCESSING**************
        
c     ********************PROCESS SECONDARY NEIGHBOURS********************
        
        if(newplst)then
          
c     zero accumulators for secondary neighbour energies and virial
          
          engcp2=0.d0
          vircp2=0.d0
          engsr2=0.d0
          virsr2=0.d0
          
c     zero secondary forces
          
          do i=1,natms
            
            fxx(i)=0.d0
            fyy(i)=0.d0
            fzz(i)=0.d0
            
          enddo
          
c     zero solvation and excitation arrays
        
          if(lsolva)then
            
            vdw_sol(:)=0.d0
            cou_sol(:)=0.d0
            
            if(lghost)then
              
              vdw_exc(:)=0.d0
              cou_exc(:)=0.d0
              
            endif
            
          endif

c     zero stress tensor
          
          do i=1,9
            stress(i)=0.d0
          enddo
          
          ii=0
          do i=idnode+1,natms,mxnode
            
            ii=ii+1
            
c     calculate interatomic vectors
            
            ik=0
            do k=1,lentry(ii)
              
              j=list(ii,k)
              
              if(j.gt.0)then
                
                ik=ik+1
                ilist(ik)=j
                xdf(ik)=xxx(i)-xxx(j)
                ydf(ik)=yyy(i)-yyy(j)
                zdf(ik)=zzz(i)-zzz(j)
                
              endif
              
            enddo
            
c     periodic boundary condition only for interactions > rprim
            
            call images(imcon,0,1,ik,cell,xdf,ydf,zdf)
            
c     square of interatomic distances
            
            do k=1,ik
              
              rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
              
            enddo
            
c     short range forces for secondary shell
            
            if((mod(keyfce,2).eq.1).and.(rvdw.gt.rprim-delr))then
              
              call srfrce
     x          (lsolva,lfree,lghost,i,ik,engacc,viracc,rvdw,dlrpot)
              
              engsr2=engsr2+engacc
              virsr2=virsr2+viracc
              
            endif
            
c     calculate coulombic force and potential terms
            
            if(keyfce/2.eq.1.or.keyfce/2.eq.6)then
              
              call error(idnode,424)
              
            elseif(keyfce/2.eq.2)then
              
c     distance dependent dielectric
              
              call coul2
     x          (lsolva,lfree,lghost,i,ik,engacc,viracc,rcut1,epsq)
              
              engcp2=engcp2+engacc
              vircp2=vircp2+viracc
              
            elseif(keyfce/2.eq.3)then
              
c     coulombic potential
              
              call coul0
     x          (lsolva,lfree,lghost,i,ik,engacc,viracc,rcut1,epsq)
              
              engcp2=engcp2+engacc
              vircp2=vircp2+viracc
              
            elseif(keyfce/2.eq.4)then
              
c     truncated shifted coulombic potential
              
              call error(idnode,424)
              
            endif
            
c     accumulate radial distribution functions : out to rcut
            
            if(lgr)call rdf0(i,ik,rcut)
            
          enddo
          
c     store secondary forces and stress tensor
          
          do i=1,natms
            
            fpx(i)=fxx(i)
            fpy(i)=fyy(i)
            fpz(i)=fzz(i)
            
          enddo
          
          do i=1,9
            stresp(i)=stress(i)
          enddo
          
c     store solvation and excitation arrays
          
          if(lsolva)then
            
            vdw_sol_put(:)=vdw_sol(:)
            cou_sol_put(:)=cou_sol(:)
            
            if(lghost)then
              
              vdw_exc_put(:)=vdw_exc(:)
              cou_exc_put(:)=cou_exc(:)
              
            endif
            
          endif
          
        endif
        
c     ****************END OF SECONDARY NEIGHBOUR PROCESSING*************
        
c     ********************PROCESS PRIMARY NEIGHBOURS********************
        
c     zero accumulators for total energies and virials
        
        engcpe=0.d0
        engsrp=0.d0
        vircpe=0.d0
        virsrp=0.d0
        
c     zero primary forces
        
        do i=1,natms
          
          fxx(i)=0.d0
          fyy(i)=0.d0
          fzz(i)=0.d0
          
        enddo
        
c     zero stress tensor
        
        do i=1,9
          stress(i)=0.d0
        enddo
        
c     zero solvation and excitation arrays
        
        if(llsolva)then
          
          vdw_sol(:)=0.d0
          cou_sol(:)=0.d0
          
          if(lghost)then
            
            vdw_exc(:)=0.d0
            cou_exc(:)=0.d0
            
          endif
          
        endif
        
c     calculate primary pair force contributions
        
        ii=0
        
        do i=idnode+1,natms,mxnode
          
          ii=ii+1
          
c     calculate interatomic distances
          
          ik=0
          
          do k=1,lentry(ii)
            
            j=-list(ii,k)
            
            if(j.gt.0)then
              
              ik=ik+1
              ilist(ik)=j
              xdf(ik)=xxx(i)-xxx(j)
              ydf(ik)=yyy(i)-yyy(j)
              zdf(ik)=zzz(i)-zzz(j)
              
            endif
            
          enddo
          
c     periodic boundary conditions
          
          call images(imcon,0,1,ik,cell,xdf,ydf,zdf)
          
c     square of interatomic distances
          
          do k=1,ik
            
            rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
            
          enddo
          
c     accumulate radial distribution functions : out to rcut
          
          if(lgr)call rdf0(i,ik,rcut)
          
c     calculate short range force and potential terms
          
          if(mod(keyfce,2).eq.1)then
            
            call srfrce
     x        (llsolva,lfree,lghost,i,ik,engacc,viracc,rvdw,dlrpot)
            
            engsrp=engsrp+engacc
            virsrp=virsrp+viracc
            
          endif
          
c     calculate coulombic force and potential terms
          
          if(keyfce/2.eq.1.or.keyfce/2.eq.6)then
            
            call error(idnode,424)
            
          elseif(keyfce/2.eq.2)then
            
c     distance dependent dielectric
            
            call coul2
     x        (llsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
            
            engcpe=engcpe+engacc
            vircpe=vircpe+viracc
            
          elseif(keyfce/2.eq.3)then
            
c     coulombic potential
            
            call coul0
     x        (llsolva,lfree,lghost,i,ik,engacc,viracc,rcut,epsq)
            
            engcpe=engcpe+engacc
            vircpe=vircpe+viracc
            
          elseif(keyfce/2.eq.4)then
            
            call error(idnode,424)
            
          endif
          
        enddo
        
c     **************END OF PRIMARY NEIGHBOUR PROCESSING*****************
        
c     counter for rdf statistics outside loop structure
        
        if(lgr)numrdf=numrdf+1
        
c     add in secondary and tertiary neighbour contributions to 
c     force, energy etc
        
        do i=1,natms
          
          fxx(i)=fxx(i)+fpx(i)+flx(i)
          fyy(i)=fyy(i)+fpy(i)+fly(i)
          fzz(i)=fzz(i)+fpz(i)+flz(i)
          
        enddo
        
        do i=1,9
          stress(i)=stress(i)+stresp(i)
        enddo
        
        engsrp=engsrp+engsr2
        virsrp=virsrp+virsr2
        
        engcpe=engcpe+engcp2+engcp3 
        vircpe=vircpe+vircp2+vircp3
        
c     calculate solvation and excitation arrays
        
        if(llsolva)then
          
          vdw_sol(:)=vdw_sol(:)+vdw_sol_put(:)+vdw_sol_lng(:)
          cou_sol(:)=cou_sol(:)+cou_sol_put(:)+cou_sol_lng(:)
          
          if(lghost)then
            
            vdw_exc(:)=vdw_exc(:)+vdw_exc_put(:)+vdw_exc_lng(:)
            cou_exc(:)=cou_exc(:)+cou_exc_put(:)+cou_exc_lng(:)
            
          endif
          
        endif
        
c     sum up contributions to short range and coulombic potential
        
        if(mxnode.gt.1)then 
          
          buffer(1)=engsrp
          buffer(2)=virsrp
          buffer(3)=engcpe
          buffer(4)=vircpe
          buffer(5)=vdw_fre
          buffer(6)=cou_fre
          buffer(7)=vdw_vir
          buffer(8)=cou_vir
          call gdsum(buffer(1),8,buffer(9))
          engsrp=buffer(1)
          virsrp=buffer(2)
          engcpe=buffer(3)
          vircpe=buffer(4)
          vdw_fre=buffer(5)
          cou_fre=buffer(6)
          vdw_vir=buffer(7)
          cou_vir=buffer(8)
          
        endif
        
      endif
      
      return
      end subroutine multiple_nsq
      
      subroutine neutlst
     x  (newlst,lchk,isn,imcon,idnode,ineu,ia,ll,
     x  txx,tyy,tzz,uxx,uyy,uzz)
      
c***********************************************************************
c     
c     dlpoly routine to create pair lists for neutral group
c     implementations.
c     loops over group ineu
c     
c     replicated data version
c     
c     copyright daresbury laboratory 1994
c     author t.forester march 1994
c     
c     isn = -1 => secondary neighbours
c     isn =  1 => primary neighbours - must contain excld interactions
c     
c***********************************************************************
      
      implicit none
      
      logical newlst,lchk,lexc
      integer isn,imcon,idnode,ineu,ia,ll,i,jj,jj0,jj1
      integer fail,ibig,keyexc,lenia,j,jneu,in0,in1
      real(8) txx,tyy,tzz,uxx,uyy,uzz
      
      dimension txx(mxatms),tyy(mxatms),tzz(mxatms)
      dimension uxx(mxatms),uyy(mxatms),uzz(mxatms)
      
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      
      data fail/0/
      
c     allocate work arrays
      
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail)
      if(fail.ne.0)call error(idnode,1830)
      
      if(newlst)then
        
        ibig=0
        
c     set cutoff radius
        
        ll=0
        
c     number of excludes found
        
        if(isn.lt.0)then
          keyexc=nexatm(ia)+2
        else
          keyexc=1
        endif
        
c     do centre - centre distances
        
        lenia=lentry(ia)
        
        do j=1,lenia
          
          jneu=abs(list(ia,j))
          xxt(j)=uxx(ineu)-uxx(jneu)
          yyt(j)=uyy(ineu)-uyy(jneu)
          zzt(j)=uzz(ineu)-uzz(jneu)
          
        enddo
        
        call images(imcon,0,1,lenia,cell,xxt,yyt,zzt)
        
c     working intragroup vectors of central group 
c     - for periodic boundaries
        
        in0=neulst(ineu)
        in1=neulst(ineu+1)-1
        
c     loop over neutral groups sites of a  
        
        
c     loop over groups in list
        
        do jj=1,lentry(ia)
          
          jneu=list(ia,jj)*isn
          
          if(jneu.gt.0)then
            
            do i=in0,in1
              
              jj0=neulst(jneu)
              jj1=neulst(jneu+1)-1
              
              if(ineu.eq.jneu)jj0=i+1
              
c     loop over jneu sites
              
              do j=jj0,jj1
                
c     reject atoms in excluded pair list
                
                lexc=.false.     
                
                if(keyexc.lt.nexatm(ia))then
                  
                  if(lexatm(ia,keyexc).eq.i)then
                    if(lexatm(ia,keyexc+1).eq.j)then
                      lexc=.true.
                      keyexc=keyexc+2
                    endif
                  endif   
                  
                endif
                
c     reject frozen atom pairs
                
                if(lstfrz(i).ne.0)then
                  if(lstfrz(j).ne.0)lexc=.true.
                endif
                
                if(.not.lexc)then
                  
                  ll=ll+1
                  if(ll.le.mxxdf)then
                    
                    xdf(ll)=txx(i)+xxt(jj)-txx(j)
                    ydf(ll)=tyy(i)+yyt(jj)-tyy(j)
                    zdf(ll)=tzz(i)+zzt(jj)-tzz(j)
                    rsqdf(ll)=xdf(ll)**2+ydf(ll)**2+zdf(ll)**2
                    ilist(ll)=i
                    jlist(ll)=j
                    
                  else
                    
                    lchk=.false.
                    ibig=max(ibig,ll)
                    
                  endif
                  
                endif
                
              enddo
              
            enddo
            
          endif
          
        enddo
        
      endif
      
c     deallocate work arrays
      
      deallocate (xxt,yyt,zzt,stat=fail)
      
      return
      end subroutine neutlst
      
      end module forces_module
      
