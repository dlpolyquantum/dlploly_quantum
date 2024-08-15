      module define_system_module
      
c***********************************************************************
c     
c     dl_poly module for utility subroutines and functions
c     copyright - daresbury laboratory
c     author    - w. smith     aug 2006
c     adapted   - p.-a. cazade oct 2007, solvation etc
c     adapted   - d. quigley   nov 2010, metadynamics
c     adapted   - w. smith     aug 2016, pimd
c     
c***********************************************************************
      
      use angles_module
      use bonds_module
      use config_module
      use core_shell_module
      use dihedral_module
      use ensemble_tools_module
      use error_module
      use ewald_module
      use exclude_module
      use external_field_module
      use four_body_module
      use hkewald_module
      use hyper_dynamics_module
      use inversion_module
      use metafreeze_module
      use metal_module
      use parse_module
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
      use nhc_module
      use water_module
      use fqcmd_module

      contains
      
      subroutine simdef
     x  (seek,lfcap,lgofr,lnsq,loptim,lzero,lminim,lpgr,ltraj,ltscal,
     x  lzeql,lzden,nolink,newgau,lhit,lbpd,ltad,lneb,prechk,tadall,
     x  lsolva,lfree,lfrmas,lexcite,lswitch,lghost,lnfic,nebgo,lpsoc,
     x  lpimd,inhc,lmsite,lcorr,lfqcmd,idnode,minstp,intsta,istraj,
     x  keybpd,keyens,keyfce,keyres,keyver,keytrj,keycorr,molcorr,
     x  wrtcorr,kmax1,kmax2,kmax3,multt,nstack,nstbgr,nsbzdn,nstbpo,
     x  nhko,nlatt,nstbts,nsteql,nstraj,nstrun,nospl,keytol,numgau,
     x  khit,nhit,nblock,ntrack,blkout,numneb,mode,nsolva,isolva,
     x  nofic,nbeads,nchain,nrespa,g_qt4f,alpha,delr,epsq,fmax,press,
     x  quattol,rcut,rprim,rvdw,taup,taut,temp,timcls,timjob,tolnce,
     x  tstep,rlxtol,opttol,zlen,ehit,xhit,yhit,zhit,ebias,vmin,
     x  catchrad,sprneb,deltad,tlow,hyp_units,chi,nsp1)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading the simulation control parameters
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith july 1992.
c     adapted   - p.-a. cazade oct 2007, solvation etc
c     
c     modified
c     author   - t.forester       may  1993
c     
c***********************************************************************
      
      implicit none
      
      character*8 cunit,seek
      character*1 hms
      character*1 directive(lenrec)
      logical lsolva,lfree,lfrmas,lexcite,lswitch,lghost,lnfic,lpsoc
      logical ltscal,lzeql,loptim,ltraj,lfcap,lgofr,lpgr,lpres,safe
      logical lstep,ltemp,lcut,ldelr,lprim,lrfce,lens,novdw,lrvdw,kill
      logical lnsq,lzden,lewald,lspme,lhke,loop,lzero,nolink,newgau
      logical lminim,lminopt,ltad,lneb,lhit,lbpd,prechk,tadall,nebgo
      logical lpimd,lver,inhc,lmsite,lcorr,lfqcmd
      integer idnode,intsta,istraj,keyens,keyfce,keyres,nstbpo,nsbzdn
      integer keytrj,keycorr,molcorr,wrtcorr,kmax1,kmax2,kmax3,multt
      integer nstack,nstbgr,khit
      integer nhit,nhko,nlatt,nstbts,nsteql,nstraj,nstrun,nospl,ntrack
      integer idum,imcon,keyver,keytol,nblock,blkout,numgau,nbeads
      integer minstp,numneb,i,keybpd,mode,nsolva,isolva,nofic,nchain
      integer nrespa
      real(8) alpha,delr,epsq,fmax,press,quattol,rcut,rprim,rvdw,taup
      real(8) taut,temp,timcls,timjob,tolnce,tstep,rlxtol,opttol
      real(8) eps,tol,fm,densvar,delrdf,delzdn,zlen,ehit,hyp_units
      real(8) catchrad,sprneb,deltad,tlow,xhit,yhit,zhit,ebias,vmin
      real(8) prntim,chi
      real(8) g_qt4f
      integer nsp1
      
CSGIC      real(8) dummy
CCRAY      real(8) dummy
CFFTWc     FFTW instruction codes
CFFTW
CFFTW      integer FFTW_FORWARD,FFTW_BACKWARD
CFFTW      parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
CFFTW
CFFTW      integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
CFFTW      parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)
CFFTW
CFFTW      integer FFTW_ESTIMATE,FFTW_MEASURE
CFFTW      parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
CFFTW
CFFTW      integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
CFFTW      parameter (FFTW_OUT_OF_PLACE=0)
CFFTW      parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
CFFTW
CFFTW      integer FFTW_THREADSAFE
CFFTW      parameter (FFTW_THREADSAFE=128)
CFFTW
      
c     intitialize system variables: temperature,pressure,ensemble key
c     force key, cutoff, primary cutoff, verlet shell width, relative
c     dielectric constant,timestep,temperature scaling flag, 
c     temp scaling interval
      
      mode=0
      nhko=0
      nlatt=0
      nsteql=0
      nstrun=0
      minstp=0
      keybpd=0
      keyres=0
      keyens=0
      keyver=0
      nstbts=0
      nstbgr=0
      nsbzdn=0
      nstbpo=100
      nstack=mxstak
      intsta=0
      nstraj=0
      istraj=1
      keytrj=0
      numgau=1
      kmax1=0
      kmax2=0
      kmax3=0
      nospl=min(8,mxspl)
      isolva=1
      nsolva=0
      niswitch=0
      nswitch=0
      nofic=1000
      nbeads=1
      nchain=1
      nrespa=1
      keyfce=0
      multt=1
      keytol=0
      alpha=0.d0
      taut=0.d0
      taup=0.d0
      fmax=1000.d0
      tstep=0.d0
      temp=0.d0
      press=0.d0
      rcut=0.d0
      rprim=0.d0
      rvdw=0.d0
      delr=0.d0
      epsq=1.d0
      rlxtol=1.d0
      opttol=1.d0
      tolnce=1.d-8
      quattol=1.d-8
      timjob=0.d0
      timcls=0.d0
      delrdf=0.d0
      delzdn=0.d0
      zlen=0.d0
      ehit=0.d0
      xhit=0.d0
      yhit=0.d0
      zhit=0.d0
      vmin=0.d0
      ebias=0.d0
      catchrad=0.d0
      pfree=0.d0
      chi=0.d0
      
      lhit=.false.
      lbpd=.false.
      ltad=.false.
      lneb=.false.
      loop=.true.
      lnfic=.false.
      lpsoc=.false.
      lzero=.false.
      ltscal=.false.
      lewald=.false.
      lspme=.false.
      lhke=.false.
      lgofr=.false.
      lpgr=.false.
      lzeql=.true.
      loptim=.false.
      lminim=.false.
      lminopt=.false.
      ltraj=.false.
      lfcap=.false.
      ltemp=.false.
      lstep=.false.
      lcut=.false.
      ldelr=.false.
      lprim=.false.
      lrfce=.false.
      lens=.false.
      novdw=.false.
      lrvdw=.false.
      lpres=.false.
      kill=.false.
      lnsq=.false.
      lzden=.false.
      nolink=.false.
      newgau=.false.
      prechk=.false.
      tadall=.false.
      lsolva=.false.
      lfree=.false.
      lfrmas=.false.
      lexcite=.false.
      lswitch=.false.
      lghost=.false.
      nebgo=.true.
      lpimd=.false.
      inhc=.false.
      lmsite=.false.
      lcorr=.false.
      lfqcmd=.false.
      lver=.false.
      seek='all     '
      
c     open the simulation input file
      
      if(idnode.eq.0)open(nread,file='CONTROL',status='old')
      
c     read job title
      
      call getrec(safe,idnode,nread)
      if(.not.safe)call abort_control_read(1,idnode,nread)
      
      call copystring(record,sysname,80)
      if(idnode.eq.0)then 
        
        write(nrite,"(3(1x,120('*'),/),1x,15('*'),5x,80a1,5x,15('*'),/,
     x    3(1x,120('*'),/),/,/,1x,'SIMULATION CONTROL PARAMETERS',/)")
     x    sysname
        
      endif
      
c     read and process directives from CONTROL file
      
      do while(loop)
        
        call getrec(safe,idnode,nread)
        if(.not.safe)call abort_control_read(1,idnode,nread)
        
c     convert to lowercase and strip out leading blanks
        
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        call copystring(record,directive,lenrec)
        
        if(record(1).eq.'#'.or.record(1).eq.' ')then
          
c     record is commented out
          cycle
          
        elseif(findstring('redirect',directive,idum))then
          
c     ignore this option in this context
          cycle
          
        elseif(findstring('steps',directive,idum))then
          
c     number of timesteps
          
          nstrun=intstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'selected number of timesteps',3x,i10)")nstrun
          
        elseif(findstring('integ',directive,idum))then
          
c     choice of integration algorithm
          
          if(findstring('leapfrog',directive,idum))then
            
            if(lver)then
              call error(idnode,-94)
              kill=.true.
            endif
            lver=.true.
            keyver=0
            if(idnode.eq.0)write(nrite,
     x        "(/,1x,'leapfrog verlet integration selected')")
            
          elseif(findstring('velocity',directive,idum))then
            
            if(lver)then
              call error(idnode,-94)
              kill=.true.
            endif
            lver=.true.
            keyver=1
            if(idnode.eq.0)write(nrite,
     x        "(/,1x,'velocity verlet integration selected')")
            
          endif
          
        elseif(findstring('no fic',directive,idum))then
          
c     cancel possible "flying ice cube" in Berendsen thermostats
          
          lnfic=.true.
          nofic=intstr(directive,lenrec,idum)
          
        elseif(findstring('shells',directive,idum).and.
     x      findstring('on',directive,idum).and.
     x      findstring('cores',directive,idum))then
          
c     put shells on cores at start - shell model only (else null)
          
          lpsoc=.true.
          
        elseif(findstring('densvar',directive,idum))then
          
c     specify allowed density variation
          
          densvar=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'allowed density variation   ',3x,1p,e12.4)")
     x      densvar
          
        elseif(findstring('no link',directive,idum))then
          
c     switch off link cell option
          
          nolink=.true.
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'link cells option deactivated')")
          
        elseif(findstring('equil',directive,idum))then
          
c     number of equilibration timesteps
          
          nsteql=intstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'equilibration period        ',3x,i10)")nsteql
          
        elseif(findstring('restart',directive,idum))then
          
c     restart control
          
          if(findstring('noscale',directive,idum))then
            
            keyres=3
            if(idnode.eq.0)write(nrite,
     x        "(/,1x,'noscale restart requested')")
            
          elseif(findstring('scale',directive,idum))then
            
            keyres=2
            if(idnode.eq.0)write(nrite,
     x        "(/,1x,'scaled restart requested')")
            
          else
            
            keyres=1
            if(idnode.eq.0)write(nrite,"(/,1x,'restart requested')")
            
          endif
          
        elseif(findstring('ensemble',directive,idum))then
          
c     ensemble selection
          
          call ensemble_selection(directive,lens,kill,inhc,idnode,
     x         keyens,mode,nchain,nrespa,taut,taup)
        
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

c     Assigning gamma and (1-gamma) values if water model qtip4p/f 
c     is requested in CONTROL file

        elseif(findstring('qtip4pf',directive,idum))then

          lmsite=.true.
          g_qt4f=0.73612d0
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'qtip4p/f water model is requested',
     x      /,1x,'gamma set to  ',3x,e12.4)")
     x      g_qt4f
 
c *******************************************************************      

        elseif(findstring('regauss',directive,idum))then
          
c     re-initialise velocities option (regaussing)
          
          newgau=.true.
          numgau=intstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'regaussing option activated',
     x      /,1x,'regaussing interval set to  ',3x,i10)")
     x      numgau
          
        elseif(findstring('scale',directive,idum))then
          
          nstbts=intstr(directive,lenrec,idum)
          if(nstbts.gt.0)then
            ltscal=.true.
            if(idnode.eq.0)write(nrite,
     x        "(/,1x,'temperature scaling on' 
     x        /,1x,'temperature scaling interval',3x,i10)")
     x        nstbts
            
          endif
          
        elseif(findstring('rdf',directive,idum))then
          
          if(findstring('print',directive,idum))then
            
            lpgr=.true.
            lpgr=(lgofr.and.lpgr)
            if(idnode.eq.0)write(nrite,
     x        "(/,1x,'g(r) printing option on      ')")
            
          else
            
            lgofr=.true.
            nstbgr=intstr(directive,lenrec,idum)
            delrdf=dblstr(directive,lenrec,idum)
            if(nstbgr.eq.0)nstbgr=10
            if(delrdf.lt.1.d-8)delrdf=0.05d0
            
            if(idnode.eq.0)then
              
              write(nrite,
     x          "(/,/,1x,'radial distribution functions on ',
     x          /,1x,'g(r) collection interval    ',3x,i10)")nstbgr
              write(nrite,
     x          "(1x,'g(r) bin width              ',3x,1p,e12.4)")
     x          delrdf
              
            endif
            
          endif
          
        elseif(findstring('zden',directive,idum))then
          
          lzden=.true.
          nsbzdn=intstr(directive,lenrec,idum)
          delzdn=dblstr(directive,lenrec,idum)
          zlen=dblstr(directive,lenrec,idum)
          if(nsbzdn.eq.0)nsbzdn=10
          if(delzdn.lt.1.d-8)then
            zlen=0.1d0*dble(mxzdn)
            delzdn=0.1d0
          elseif(zlen.lt.1.d-8)then
            zlen=delzdn*dble(mxzdn)
          endif
          if(idnode.eq.0)then
            
            write(nrite,
     x        "(/,/,1x,'Z density profile requested',
     x        /,1x,'zdensity collection interval',3x,i10)")nsbzdn
            write(nrite,
     x        "(1x,'zdensity bin width          ',3x,1p,e12.4)")
     x        delzdn
            write(nrite,
     x        "(1x,'zdensity range              ',3x,1p,e12.4)")
     x        zlen
            
          endif
          
        elseif(findstring('collect',directive,idum))then
          
          lzeql=.false.
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'equilibration included in overall averages')")
          
        elseif(findstring('neb',directive,idum))then
          
c     activate nudged elastic band option
          
          call neb_option(directive,lneb,lminopt,idnode,
     x      numneb,keytol,sprneb,opttol,hyp_units)
          
c     read path integral option
c*******************************************************************
c 
c   Additional PI methods added
c
c   copyright - Dil Limbu and Nathan London
c   authors - Dil Limbu and Nathan London 2023
c
c*******************************************************************            

        elseif(findstring('pimd',directive,idum))then
          
          if(lver)then
            call error(idnode,-94)
            kill=.true.
          endif
          if(lens)then
            call error(idnode,-414)
            kill =.true.
          endif
          lens=.true.
          lver=.true.
          lpimd=.true.
          keyver=2
          keyens=40
          if(findstring('npt',directive,idum))then
c barostat options            
            if(findstring('nhc',directive,idum))then
              keyens=51
              nbeads=intstr(directive,lenrec,idum)
              nrespa=intstr(directive,lenrec,idum)
              nchain=intstr(directive,lenrec,idum)
              taut=dblstr(directive,lenrec,idum)
              taup=dblstr(directive,lenrec,idum)
              if((nbeads.eq.0).or.(nrespa.eq.0).or.(nchain.eq.0).or.
     x          (taut.eq.0).or.(taup.eq.0).or.(nchain.eq.1))
     x          call error(idnode,4011)
              nrespa=max(nrespa,1)
              nchain=max(nchain,1)
              inhc=.true.         
            elseif(findstring('pile',directive,idum))then
              keyens=52
              nbeads=intstr(directive,lenrec,idum)
              taut=dblstr(directive,lenrec,idum)
              taup=dblstr(directive,lenrec,idum)
              if((nbeads.eq.0).or.(taut.eq.0).or.(taup.eq.0))
     x          call error(idnode,4012)
            endif
          elseif(findstring('nvt',directive,idum))then
            keyens=40
            nbeads=intstr(directive,lenrec,idum)
            taut=dblstr(directive,lenrec,idum)
            if((nbeads.eq.0).or.(taut.eq.0))
     x        call error(idnode,4013)
          elseif(findstring('gth',directive,idum))then
            keyens=41
            nbeads=intstr(directive,lenrec,idum)
            taut=dblstr(directive,lenrec,idum)
            chi=dblstr(directive,lenrec,idum)
            if((nbeads.eq.0).or.(taut.eq.0).or.(chi.eq.0))
     x        call error(idnode,4014)
          elseif(findstring('nhc',directive,idum))then
            keyens=42
            nbeads=intstr(directive,lenrec,idum)
            nrespa=intstr(directive,lenrec,idum)
            nchain=intstr(directive,lenrec,idum)
            taut=dblstr(directive,lenrec,idum)
            if((nbeads.eq.0).or.(nrespa.eq.0).or.(nchain.eq.0).or.
     x        (taut.eq.0).or.(nchain.eq.1))
     x        call error(idnode,4015)
            nrespa=max(nrespa,1)
            nchain=max(nchain,1)
          elseif(findstring('nm',directive,idum))then
            keyens=43
            nbeads=intstr(directive,lenrec,idum)
            nrespa=intstr(directive,lenrec,idum)
            nchain=intstr(directive,lenrec,idum)
            taut=dblstr(directive,lenrec,idum)
            if((nbeads.eq.0).or.(nrespa.eq.0).or.(nchain.eq.0).or.
     x        (taut.eq.0).or.(nchain.eq.1))
     x        call error(idnode,4016)
            nrespa=max(nrespa,1)
            nchain=max(nchain,1)
          elseif(findstring('pile',directive,idum))then
            keyens=44
            nbeads=intstr(directive,lenrec,idum)
            taut=dblstr(directive,lenrec,idum)
            if((nbeads.eq.0).or.(taut.eq.0))
     x        call error(idnode,4017)
          elseif(findstring('piglet',directive,idum))then
            keyens=45
            nbeads=intstr(directive,lenrec,idum)
            nsp1=intstr(directive,lenrec,idum)
            if((nbeads.eq.0).or.(nsp1.eq.0))
     x        call error(idnode,4018)
            nsp1=nsp1+1
          elseif(findstring('nve',directive,idum))then
            keyens=61
            nbeads=intstr(directive,lenrec,idum)
            if(nbeads.eq.0)
     x        call error(idnode,4019)
          elseif(findstring('pacmd',directive,idum))then
            keyens=62
            nbeads=intstr(directive,lenrec,idum)
            nrespa=intstr(directive,lenrec,idum)
            nchain=intstr(directive,lenrec,idum)
            taut=dblstr(directive,lenrec,idum)
            if((nbeads.eq.0).or.(nrespa.eq.0).or.(nchain.eq.0).or.
     x        (taut.eq.0).or.(nchain.eq.1))
     x        call error(idnode,4020)
            nrespa=max(nrespa,1)
            nchain=max(nchain,1)
          elseif(findstring('trpmd',directive,idum))then
            keyens=63
            nbeads=intstr(directive,lenrec,idum)
            if(nbeads.eq.1)
     x        call error(idnode,4021)
          else
c     default is nvt
            keyens=40
            nbeads=intstr(directive,lenrec,idum)
            taut=dblstr(directive,lenrec,idum)
            if((nbeads.eq.0).or.(taut.eq.0))
     x        call error(idnode,4022)
          endif
          if(nbeads.eq.0)nbeads=num_beads_default
          if(taut.le.1.d-6)taut=1.d0
          
          if(idnode.eq.0)then
            write(nrite,"(/,1x,'PIMD option selected')")
            write(nrite,"(1x,'Number of quantum beads/atom :',i5)")
     x        nbeads
            if(keyens.eq.40)then
              write(nrite,
     x          "(1x,'Canonical Ensemble with Nose-Hoover Thermostat')")
              write(nrite,"(1x,'Thermostat relaxation time (ps):',
     x          1p,e12.4)")taut
            elseif(keyens.eq.41)then
              write(nrite,
     x          "(1x,'Canonical Ensemble with Gentle Thermostat')")
              write(nrite,"(1x,'Thermostat relaxation time (ps):',
     x          1p,e12.4)")taut
              write(nrite,"(1x,'Stochastic force parameter:',
     x          1p,e12.4)")chi
            elseif(keyens.eq.42)then
              write(nrite,
     x          "(1x,'Canonical Ensemble with Nose-Hoover Chains')")
              write(nrite,"(1x,'Number of Nose-Hoover chains :',i5)")
     x          nchain
              write(nrite,"(1x,'Thermostat relaxation time (ps):',
     x          1p,e12.4)")taut
            elseif(keyens.eq.43)then
              write(nrite,
     x          "(1x,'Canonical Ensemble in normal modes with mNHC')")
              write(nrite,"(1x,'Number of  RESPA steps       :',i5)")
     x          nrespa
              write(nrite,"(1x,'Number of Nose-Hoover chains :',i5)")
     x          nchain
              write(nrite,"(1x,'Thermostat relaxation time (ps):',
     x          1p,e12.4)")taut
            elseif(keyens.eq.44)then
              write(nrite,
     x          "(1x,'Canonical Ensemble in normal modes with PILE')")
              write(nrite,"(1x,'Thermostat relaxation time (ps):',
     x          1p,e12.4)")taut
            elseif(keyens.eq.45)then
              write(nrite,
     x          "(1x,'Canonical Ensemble in normal modes with PIGLET')")
              write(nrite,"(1x,'Thermostat w/ extra no. of momenta of:',
     x          1p,i5)")nsp1-1
            elseif(keyens.eq.51)then
              write(nrite,
     x        "(/,1x,'NPT Simulation in normal mode with mNHC',
     x          /,1x,'thermostat relaxation time',1p,e12.4,
     x          /,1x,'barostat relaxation time  ',1p,e12.4,
     x          /,1x,'number of RESPA steps  ',1p,i6,
     x          /,1x,'number of chains       ',1p,i6)")
     x                taut,taup,nrespa,nchain
            elseif(keyens.eq.52)then
              write(nrite,
     x        "(/,1x,'NPT Simulation in normal mode with PILE',
     x          /,1x,'thermostat relaxation time',1p,e12.4,
     x          /,1x,'barostat relaxation time  ',1p,e12.4)")
     x                taut,taup
            elseif(keyens.eq.61)then
              write(nrite,
     x          "(1x,'RPMD in normal modes')")
            elseif(keyens.eq.62)then
              write(nrite,
     x          "(1x,'Partialy Adiabatic CMD')")
              write(nrite,"(1x,'Number of  RESPA steps       :',i5)")
     x          nrespa
              write(nrite,"(1x,'Number of Nose-Hoover chains :',i5)")
     x          nchain
              write(nrite,"(1x,'Thermostat relaxation time (ps):',
     x          1p,e12.4)")taut
            elseif(keyens.eq.63)then
              write(nrite,
     x          "(1x,'Thermostatted RPMD')")
            endif
          endif
         
c     correlation function option
        elseif(findstring('corr',directive,idum))then
          lcorr=.true.
          if(findstring('velocity',directive,idum))then
            keycorr=1
            molcorr=intstr(directive,lenrec,idum)
            wrtcorr=intstr(directive,lenrec,idum)
          elseif(findstring('dipole',directive,idum))then
            keycorr=2
            molcorr=intstr(directive,lenrec,idum)
            wrtcorr=intstr(directive,lenrec,idum)
c         default is velocity autocorrelation for molecule type 1
c         written every 10 steps
          else
            keycorr=1
            wrtcorr=10
          endif
          if(molcorr.eq.0) molcorr=1
          if(idnode.eq.0)then
            if(keycorr.eq.1)then
              write(nrite,"(/,1x,'Correlation function: velocity')")
              write(nrite,"(/,1x,'Correlation function molecule:',i5)")
     x          molcorr
              write(nrite,"(/,1x,'Correlation function steps:',i5)")
     x          wrtcorr
            elseif(keycorr.eq.2)then 
              write(nrite,"(/,1x,'Correlation function: dipole')")
              write(nrite,"(/,1x,'Correlation function molecule:',i5)")
     x          molcorr
              write(nrite,"(/,1x,'Correlation function steps:',i5)")
     x          wrtcorr
            endif
          endif
c     fast-qcmd option
        elseif(findstring('fqcmd',directive,idum))then
          lfqcmd=.true.
          numPot=intstr(directive,lenrec,idum)
          numPoints=intstr(directive,lenrec,idum)
          nrdfpts=intstr(directive,lenrec,idum)
          if(idnode.eq.0)then
            write(nrite,"(/,1x,'Performing f-QCMD')")
            write(nrite,"(/,1x,'Number of non-bonded tables:',i5)")
     x        numPot
            write(nrite,"(/,1x,'Number of points in tables:',i5)")
     x        numPoints
            write(nrite,"(/,1x,'Number of points in output rdfs:',i5)")
     x        nrdfpts
          endif
       elseif(findstring('impact',directive,idum))then
          
c     activate the impact option
          
          if(lhit)call error(idnode,516)
          lhit=.true.
          khit=intstr(directive,lenrec,idum)
          nhit=intstr(directive,lenrec,idum)
          ehit=dblstr(directive,lenrec,idum)
          xhit=dblstr(directive,lenrec,idum)
          yhit=dblstr(directive,lenrec,idum)
          zhit=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)then
            
            write(nrite,"(/,1x,'impact option selected')")
            write(nrite,"(
     x        /,1x,'identity of impact atom        ',i10,
     x        /,1x,'time step of impact            ',i10,
     x        /,1x,'impact recoil energy (keV)     ',1p,e12.4,
     x        /,1x,'impact direction (x component) ',1p,e12.4,
     x        /,1x,'impact direction (y component) ',1p,e12.4,
     x        /,1x,'impact direction (z component) ',1p,e12.4)")
     x        khit,nhit,ehit,xhit,yhit,zhit
            
          endif
          
c     convert impact energy from keV to internal units
          
          ehit=ehit*9648530.821d0
          
        elseif(findstring('bpd',directive,idum))then
          
c     activate the BPD option
          
          call bpd_option(directive,seek,lbpd,ltad,lminopt,prechk,
     x      nebgo,keybpd,idnode,nblock,ntrack,keytol,ebias,vmin,
     x      catchrad,sprneb,opttol,hyp_units)
          
        elseif(findstring('tad',directive,idum))then
          
c     activate temperature accelerated dynamics option
          
          call tad_option(directive,ltad,lbpd,lminopt,prechk,tadall,
     x      idnode,nblock,ntrack,blkout,keytol,catchrad,sprneb,tlow,
     x      deltad,opttol,hyp_units)
          
        elseif(findstring('minim',directive,idum))then
          
          if(lminopt)call error(idnode,225)
          if(findstring('forc',directive,idum))keytol=0
          if(findstring('ener',directive,idum))keytol=1
          if(findstring('posi',directive,idum))keytol=2
          hyp_units=energy_unit()
          minstp=intstr(directive,lenrec,idum)
          opttol=dblstr(directive,lenrec,idum)
          call getword(cunit,directive,8,lenrec)
          lminim=.true.
          loptim=.false.
          lzero=.false.
          ltscal=.false.
          lminopt=.true.
          
          if(idnode.eq.0)then
            
            write(nrite,
     x        "(/,1x,'minimisation programme requested')")
            write(nrite,
     x        "(1x,'structure minimisation interval ',
     x        3x,i10)")minstp
            call print_optim(keytol)
            write(nrite,
     x        "(1x,'structure minimisation tolerance',
     x        3x,1p,e12.4,1x,a8)")opttol,cunit
            
          endif
          if(keytol.lt.2)opttol=opttol*hyp_units
          
        elseif(findstring('optim',directive,idum))then
          
          if(lminopt)call error(idnode,225)
          if(findstring('forc',directive,idum))keytol=0
          if(findstring('ener',directive,idum))keytol=1
          if(findstring('posi',directive,idum))keytol=2
          hyp_units=energy_unit()
          opttol=dblstr(directive,lenrec,idum)
          call getword(cunit,directive,8,lenrec)
          loptim=.true.
          lminim=.false.
          lzero=.false.
          ltscal=.false.
          lminopt=.true.
          
          if(idnode.eq.0)then
            
            write(nrite,
     x        "(/,1x,'structure optimisation requested')")
            call print_optim(keytol)
            write(nrite,
     x        "(1x,'tolerance for structure optimisation ',
     x        3x,1p,e12.4,1x,a8)")opttol,cunit
            
          endif
          if(keytol.lt.2)opttol=opttol*hyp_units
          
        elseif(findstring('zero',directive,idum))then
          
          if(lminopt)call error(idnode,225)
          temp=1.d0
          lzero=.true.
          loptim=.false.
          lminim=.false.
          ltemp=.true.
          ltscal=.false.
          lminopt=.true.
          
          if(idnode.eq.0)then
            
            write(nrite,
     x        "(/,1x,'zero K optimisation requested')")
            write(nrite,
     x        "(' temperature reset to',1p,e12.4)")1.d0
            
          endif
          
        else if(findstring('solva',directive,idum))then
          
          call solvation_option
     x      (directive,lsolva,idnode,nsolva,isolva)
          
        else if(findstring('decomp',directive,idum))then
          
          call solvation_option
     x      (directive,lsolva,idnode,nsolva,isolva)
          
        elseif(findstring('metafreeze',directive,idum).or.
     x      findstring('metadyn',directive,idum))then
          
c     activate metadynamics option - d. quigley
          
          call metadyn_option
     x      (directive,lmetadyn,lstein,ltet,lglobpe,llocpe,idnode,
     x      ncolvar,nq4,nq6,ntet,hkey,meta_step_int,globpe_scale,
     x      locpe_scale,ref_W_aug,h_aug,wt_Dt)
          
        else if(findstring('free',directive,idum))then
          
          call free_energy_option(directive,lfree,lfrmas,idnode)
          
        else if(findstring('excite',directive,idum))then
          
          call excitation_option
     x      (directive,lsolva,lexcite,lghost,idnode,nsolva,isolva)
          
        else if(findstring('switch',directive,idum))then
          
          call switching_option
     x      (directive,lsolva,lswitch,lghost,idnode,nsolva,isolva)
          
        elseif(findstring('print',directive,idum))then
          
          nstbpo=intstr(directive,lenrec,idum)
          nstbpo=max(nstbpo,1)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'data printing interval      ',3x,i10)")nstbpo
          
        elseif(findstring('stack',directive,idum))then
          
          nstack=intstr(directive,lenrec,idum)
          
c     reset stack limit if too large
          
          nstack=min(nstack,mxstak)
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'data stacking interval      ',3x,i10)")nstack
          
        elseif(findstring('stats',directive,idum))then
          
          intsta=intstr(directive,lenrec,idum)
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'statistics file interval    ',3x,i10)")intsta
          
        elseif(findstring('traj',directive,idum))then
          
          ltraj=.true.
          nstraj=intstr(directive,lenrec,idum)
          istraj=max(intstr(directive,lenrec,idum),1)
          keytrj=intstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'trajectory file option on  ',
     x      /,1x,'trajectory file start       ',3x,i10,
     x      /,1x,'trajectory file interval    ',3x,i10
     x      /,1x,'trajectory file info key    ',3x,i10)")
     x      nstraj,istraj,keytrj
          
        elseif(findstring('ewald',directive,idum).or.
     x      findstring('spme',directive,idum).or.
     x      findstring('hke',directive,idum))then
          
c     read Ewald or HK-Ewald or SPM-Ewald sum parameters
          
          call ewald_selection(directive,lhke,lspme,lewald,lcut,
     x      lrfce,kill,idnode,keyfce,imcon,nhko,nlatt,kmax1,kmax2,
     x      kmax3,alpha,rcut)
          
        elseif(findstring('distan',directive,idum))then
          
          keyfce=4
          if(idnode.eq.0)write(nrite,
     x      "(/,/,1x,'Electrostatics : Distance dependent dielectric')")
          
          if(lrfce)then
            call  error(idnode,-416)
            kill=.true.
          endif
          
          lrfce=.true.
          
        elseif(findstring('coul',directive,idum))then
          
          keyfce=6
          if(idnode.eq.0)write(nrite,
     x      "(/,/,1x,'Electrostatics : Coulombic potential')")
          
          if(lrfce)then
            call  error(idnode,-416)
            kill=.true.
          endif
          
          lrfce=.true.
          
        elseif(findstring('shift',directive,idum))then
          
          keyfce=8
          alpha=0.d0
          
          if(idnode.eq.0)write(nrite,
     x      "(/,/,1x,'Electrostatics : Shifted Coulombic potential')")
          
          if(findstring('precision',directive,idum))then
            
            eps=dblstr(directive,lenrec,idum)
            if(.not.lcut)then
              call error(idnode,-435)
              kill=.true.
            else
              
              eps=min(abs(eps),0.5d0)
              tol=sqrt(abs(log(eps*rcut)))
              alpha=sqrt(abs(log(eps*rcut*tol)))/rcut
              if(idnode.eq.0)then
                
                write(nrite,
     x            "(1x,'Specified precision parameter : ',1p,
     x            e12.4)")eps
                write(nrite,
     x            "(1x,'Calculated damping parameter: ',1p,
     x            e12.4)")alpha
                
              endif
              
            endif
            
          elseif(findstring('damp',directive,idum))then
            
            alpha=dblstr(directive,lenrec,idum)
            if(idnode.eq.0)write(nrite,
     x        "(1x,'Specified damping parameter : ',1p,e12.4)")
     x        alpha
            
          endif
          
          if(lrfce)then
            call  error(idnode,-416)
            kill=.true.
          endif
          
          lrfce=.true.
          
        elseif(findstring('reaction',directive,idum))then
          
          keyfce=10
          alpha=0.d0
          
          if(idnode.eq.0)write(nrite,
     x      "(/,/,1x,'Electrostatics : reaction field')")
          
          if(findstring('precision',directive,idum))then
            
            eps=dblstr(directive,lenrec,idum)
            if(.not.lcut)then
              call error(idnode,-435)
              kill=.true.
            else
              
              eps=min(abs(eps),0.5d0)
              tol=sqrt(abs(log(eps*rcut)))
              alpha=sqrt(abs(log(eps*rcut*tol)))/rcut
              if(idnode.eq.0)then
                
                write(nrite,
     x            "(1x,'Specified precision parameter : ',1p,
     x            e12.4)")eps
                write(nrite,
     x            "(1x,'Calculated damping parameter: ',1p,
     x            e12.4)")alpha
                
              endif
              
            endif
            
          elseif(findstring('damp',directive,idum))then
            
            alpha=dblstr(directive,lenrec,idum)
            if(idnode.eq.0)write(nrite,
     x        "(1x,'Specified damping parameter : ',1p,e12.4)")
     x        alpha
            
          endif
          
          if(lrfce)then
            call  error(idnode,-416)
            kill=.true.
          endif
          
          lrfce=.true.
          
        elseif(findstring('cap',directive,idum))then
          
          lfcap=.true.
          fm=dblstr(directive,lenrec,idum)
          if(fm.gt.0.d0)fmax=fm
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'force capping :',16x,1p,e12.4,' kT/A')")fmax
          
        elseif(findstring('no vdw',directive,idum))then
          
          if(idnode.eq.0)write(nrite,
     x      "(/,/,1x,'short-range potential terms off')")
          novdw=.true.
          
        elseif(findstring('no elec',directive,idum))then
          
          keyfce=0
          if(idnode.eq.0)write(nrite,
     x      "(/,/,1x,'electrostatic potential terms off')")
          
          if(lrfce)then
            call  error(idnode,-416)
            kill=.true.
          endif
          
          lrfce=.true.
          
        elseif(findstring('mult',directive,idum))then
          
          multt=max(intstr(directive,lenrec,idum),1)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'multiple timestep interval  ',3x,i10)")multt
          
        elseif(findstring('timestep',directive,idum))then
          
          lstep=.true.
          tstep=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'simulation timestep         ',3x,1p,e12.4)")tstep
          
        elseif(findstring('temp',directive,idum))then
          
          ltemp=.true.
          temp=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'simulation temperature      ',3x,1p,e12.4)")temp
          
        elseif(findstring('pres',directive,idum))then
          
          press=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'simulation pressure (katm)  ',3x,1p,e12.4)")press
          
c     convert from katm to internal units of pressure
          
          press=press/prsunt
          lpres=.true.
          
        elseif(findstring('prim',directive,idum))then
          
c     primary cutoff
          
          lprim=.true.
          rprim=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'primary neighbour cut off   ',3x,1p,e12.4)")rprim
          
        elseif(findstring('rvdw',directive,idum))then
          
c     cutoff for short range potentials
          
          rvdw=dblstr(directive,lenrec,idum)
          lrvdw=.true.
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'real space cut off (vdw)    ',3x,1p,e12.4)")rvdw
          
        elseif(findstring('delr',directive,idum))then
          
c     Verlet shell width
          
          ldelr=.true.
          delr=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'border width of Verlet shell',3x,1p,e12.4)")delr
          
        elseif(findstring('cut',directive,idum))then
          
c     cutoff
          
          lcut=.true.
          rcut=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'real space cut off          ',3x,1p,e12.4)")rcut
          
        elseif(findstring('eps',directive,idum))then
          
c     relative dielectric constant
          
          epsq=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,/,1x,'relative dielectric constant',3x,1p,e12.4)")epsq
          
        elseif(findstring('rlxtol',directive,idum))then
          
c     force tolerance for shell relaxation
          
          rlxtol=max(rlxtol,dblstr(directive,lenrec,idum))
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'tolerance for shell relaxn. ',3x,1p,e12.4)")rlxtol
          
        elseif(findstring('shake',directive,idum))then
          
c     tolerance for shake
          
          tolnce=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'tolerance for SHAKE         ',3x,1p,e12.4)")tolnce
          
        elseif(findstring('quaternion',directive,idum))then
          
c     tolerance for quaternion integration
          
          quattol=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'tolerance for Quaternions   ',3x,1p,e12.4)")quattol
          
        elseif(findstring('job time',directive,idum))then
          
c     time for simulation (in seconds/minutes/hours/days or indefinite)
          
          if(findstring('indef',directive,idum))then
            timjob=1.0d6*365.25d0*24.d0*60.d0*60.d0
          else
            timjob=dblstr(directive,lenrec,idum)
            if(findstring('m',directive,idum))then
              timjob=6.0d1*timjob
            elseif(findstring('h',directive,idum))then
              timjob=3.6d3*timjob
            elseif(findstring('d',directive,idum))then
              timjob=8.64d4*timjob
            endif
          endif
          
          call get_prntime(hms,timjob,prntim)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'user allocated job time (',a1,') ',3x,f8.4)")
     x      hms,prntim
          
        elseif(findstring('close time',directive,idum))then
          
c     time for winding up a job (in seconds)
          
          timcls=dblstr(directive,lenrec,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'job closure time        (s) ',3x,f8.3)")timcls
          
        elseif(findstring('all pairs',directive,idum))then
          
c     full minimum image - N^2 interactions each timestep
          
          lnsq=.true.
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'All-pairs requested for electrostatics')")
          
        elseif(findstring('nospl',directive,idum))then
          
c     set ewald_spme interpolation order
          
          nospl=intstr(directive,lenrec,idum)
          
        elseif(findstring('finish',directive,idum))then
          
c     safe termination of reading CONTROL file
          
          loop=.false.
          
        else
          
c     unrecognised directive in control file
          
          kill=.true.
          if(idnode.eq.0)write(nrite,"(/,/,100a1)")record
          call error(idnode,-3)
          
        endif
        
      enddo
      
c     check on steps before temperature scaling
      
      if(nstbts.eq.0)nstbts=nstrun+1
      
c     conduct consistency checks on directives
      
      if(lminim)then
        
c     ensure final configuration follows minimisation
        
        if(minstp.eq.0)minstp=nstrun
        nstrun=minstp*(nstrun/minstp)
        
      endif
      
c     check force options
      
      if(.not.lrfce)then
        
c     check the long range force option has been properly specified
        
        if(.not.novdw)then
          
          kill=.true.
          call error(idnode,-383)

        endif
        
      else
        
c     check the short range force option has been properly specified
        
        if(novdw)then
          
          if(keyfce.eq.0)then
            
            lcut=.true.
            ldelr=.true.
            keyfce=keyfce+1
            
          endif
          
        else
          
          keyfce=keyfce+1
          
        endif
        
      endif
      
c     if tad selected use only leap frog
      
      if(ltad.and.keyver.eq.1)then
        
        if(idnode.eq.0)write(nrite,
     x    "(/,1x,'switching to leapfrog for TAD dynamics')")
        keyver=0
        
      endif
      
c     error checking 
      
      if(lmetadyn.and.keyens.ne.3.and.keyens.ne.5.and.keyens.ne.7)then
        
        kill=.true.
        call error(idnode,-2360)
        
      endif
      
      if(lsolva.or.lfree.or.lexcite.or.lswitch)then
        
        if(lspme)then
          
          kill=.true.
          call error(idnode,-601)
          
        endif
        
        if(lhke)then
          
          kill=.true.
          call error(idnode,-602)
          
        endif
        
      endif
      
      if(lghost.and.nstbgr.ne.isolva)then
        
        call warning(idnode,130,dble(isolva),0.d0,0.d0)
        nstbgr=isolva
        
      endif
      if(lfree.and.lgofr)then
        
        call warning(idnode,140,0.d0,0.d0,0.d0)
        lgofr=.false.
        lpgr=.false.
        
      endif
      if(loptim)then
        
        temp=0.d0
        
      elseif(.not.ltemp)then
        
        kill=.true.
        call error(idnode,-380)
        
      endif
      
      if(.not.lstep)then
        
        kill=.true.
        call error(idnode,-381)
        
      endif
      
      if(.not.lcut)then
        
        kill=.true.
        call error(idnode,-382)
        
      endif
      
c     check if van der Waals cutoff set
      
      if(.not.lrvdw.and.mod(keyfce,2).eq.1)then
        
        if(rcut.gt.0.d0)then
          
          rvdw=rcut
          
        else
          
          kill=.true.
          call error(idnode,-402)
          
        endif      
        
      endif
      
      if(.not.ldelr)then
        
        kill=.true.
        call error(idnode,-384)
        
      endif
      
      if(multt.gt.1)then
        
        if(.not.lprim)then
          
          kill=.true.
          call error(idnode,-385)
          
        elseif(rprim.gt.rcut)then
          
          kill=.true.
          call error(idnode,-386)
          
        endif
        
      endif
      
c     check settings in nvt ensemble
      
      if(keyens.ge.2.and.keyens.le.3)then
        
        if(taut.le.0.d0)then
          
          kill=.true.
          call error(idnode,-464)
          
        endif
        
      endif
      
c     check settings in npt ensemble
      
      if(keyens.ge.4.and.keyens.le.7)then
        
        if(.not.lpres)then
          
          kill=.true.
          call error(idnode,-387)
          
        endif
        
c     check barostat and thermostat rates non zero
        
        if(taut.le.0.d0)then
          
          kill=.true.
          call error(idnode,-464)
          
        endif
        if(taup.le.0.d0)then
          
          kill=.true.
          call error(idnode,-466)
          
        endif
        
      endif
      
c     check multiple timestep cutoffs are sensible
      
      if(multt.gt.1)then
        if(rcut-rprim.lt.delr)then
          
          kill=.true.
          call error(idnode,-398)
          
        endif
      endif
      
c     check rcut > rvdw (for verlet list constructor)
      
      if(rcut.lt.rvdw)then 
        
        kill=.true.
        call error(idnode,-400)
        
      endif
      
c     check spme is not being used with incorrect pbc
      
      if(lspme)then
        
        if(imcon.eq.0.or.imcon.eq.6)then
          
          kill=.true.
          call error(idnode,-513)
          
        endif
        
      endif
      
c     check on all-pairs calculation request
      
      if(lnsq)then
        
        if(multt.eq.1)then
          
          kill=.true.
          call error(idnode,-422)
          
        endif
        
        if(keyfce/2.lt.2.or.keyfce/2.gt.3)then
          
          kill=.true.
          call error(idnode,-424)
          
        endif
        
      endif
      
c     cancel rdf option if no vdw or coulombic forces
      
      if(lgofr.and.keyfce.eq.0)then
        
        lgofr=.false.
        call warning(idnode,120,0.d0,0.d0,0.d0)
        
      endif
      
      if(kill)call abort_control_read(2,idnode,nread)
      
c     close CONTROL file
      
      if(idnode.eq.0)close(nread)
      
      return
      end subroutine simdef
      
      subroutine sysdef
     x  (lneut,lnsq,lsolva,lfree,lexcite,lswitch,lghost,lpimd,idnode,
     x  keyfce,keyfld,natms,ngrp,ntpatm,ntpmls,ntpvdw,ntptbp,ntpmet,
     x  ntpfbp,ntpter,nshels,keyshl,ntghost,keyver,dlrpot,engunit,
     x  rvdw,rcuttb,rctter,rcutfb)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading in the molecular specifications
c     of the system to be simulated
c     version for rigid unit data and neutral groups
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith may 1992.
c     amended   - w.smith march 1994 
c     amended   - t.forester april 1994
c     amended   - w.smith  dec 1994 - getrec etc
c     amended   - a.smondyrev may 2000 - keydih=5 for 
c     ryckaert-bellemans potential in dihedrals
c     amended   - a.smondyrev may 2000 - keydih=6 for 
c     fluorinated ryckaert-bellemans potential in dihedrals
c     adapted   - p.-a. cazade oct 2007, solvation etc
c     
c***********************************************************************
      
      implicit none
      
      logical lunits,lmols,lneut,ltable,lnsq,lshl,safe,lpmf,lrig,lcon
      logical loop1,loop2,lsolva,lfree,lexcite,lswitch,lghost,lpimd
      logical lfreeze
      
      integer idnode,mxnode,keyfce,keyfld,natms,ngrp,ntpatm,ntpmls
      integer ntpvdw,ntptbp,ntpmet,ntpfbp,nshels,ksite
      integer nsite,nconst,nangle,ndihed,ninver,nbonds
      integer nvoids
      integer nteth,nspmf,itmols,i,idum,keyver
      integer ntpter,keyshl,iatm,natmsr,ntghost
      
      real(8) dlrpot,engunit,rvdw,rcuttb,rctter,rcutfb
      real(8) sumchg
      
      data loop1/.true./,loop2/.true./
      
c     initialise system counters: atomic site index, number of 
c     constraints, bond angles, dihedrals, inversions, chemical bonds,
c     unique atom types, total number of atoms,
c     total number of rigid groups, number of tethered atoms,
c     number of three body potentials
      
      nsite=0
      nconst=0
      nangle=0
      ndihed=0
      ninver=0
      nbonds=0
      nvoids=0
      ntpatm=0
      natms=0
      ngrp=0
      nteth=0
      ntptbp=0
      ntpter=0
      ntpmet=0
      ntpvdw=0
      ntpfbp=0
      nshels=0
      nspmf=0
      keyfld=0
      keyshl=0
      ntghost=0
      natmsr=0
      ntcons_ghost=0

      lunits=.false.
      lmols=.false.
      lneut=.false.
      ltable=.false.
      lmetab=.false.
      lshl=.false.
      lpmf=.false.
      lrig=.false.
      lcon=.false.
      lfreeze=.false.
      engunit=1.d0
      
      numbonds(:)=0
      numpmf(:)=0
      numcon(:)=0
      numdih(:)=0
      numinv(:)=0
      numgrp(:)=0
      numsit(:)=0
      numteth(:)=0
      numshl(:)=0
      npmf(:)=0
      indpmf(:)=0
      
c     open force field data file
      
      if(idnode.eq.0)open (nfield,file='FIELD',status='old')
      
      if(idnode.eq.0)
     x  write(nrite,"(/,/,'SYSTEM SPECIFICATION')")
      
      call getrec(safe,idnode,nfield)
      if(.not.safe)call abort_field_read(1,idnode,nfield)
      
c     read and process directives from field file
      
      do while(loop1)
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)call abort_field_read(1,idnode,nfield)
        
c     convert to lowercase and remove leading blanks
        
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        
        if(record(1).eq.'#'.or.record(1).eq.' ')then
          
c     record is commented out
          cycle
          
        elseif(findstring('units',record,idum))then
          
c     identify energy unit for input/output
          
          lunits=.true.
          call define_units(idnode,engunit)
          
c     neutral group control option
          
        elseif(findstring('neut',record,idum))then
          
          lneut=.true.
          if(idnode.eq.0)
     x      write(nrite,"(/,' neutral group implementation in use')")
          
c     can't have neutral groups with all-pairs
          
          if(lnsq)call error(idnode,426)
          
c     specify molecular species
          
        elseif(findstring('molecu',record,idum))then
          
c     number of molecular types
          
          if(lmols)call error(idnode,11)
          lmols=.true.
          ntpmls=intstr(record,lenrec,idum)
          
          if(idnode.eq.0)
     x      write(nrite,"(/,/,1x,'number of molecular types',6x,i10)")
     x      ntpmls
          
          if(ntpmls.gt.mxtmls)call error(idnode,10)
          
c     initialise total system charge
          
          sumchg=0.d0
          
c     read in molecular characteristics
          
          do itmols=1,ntpmls
            
            if(idnode.eq.0)
     x        write(nrite,"(/,1x,'molecular species type',9x,i10)")
     x        itmols
            
c     name of molecular species
            
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abort_field_read(1,idnode,nfield)
            
            call copystring(record,molnam(1,itmols),40)
            if(idnode.eq.0)
     x        write(nrite,"(/,/,1x,'name of species:',13x,40a1)")
     x        (molnam(i,itmols),i=1,40)
            
c     stop processing if energy unit has not been specified
            
            if(.not.lunits)call error(idnode,6)
            
c     read molecular data
            
            loop2=.true.
            
            do while(loop2)
              
              call getrec(safe,idnode,nfield)
              if(.not.safe)call abort_field_read(1,idnode,nfield)
              
              call lowcase(record,lenrec)
              call strip(record,lenrec)
              
              ksite=0
              
              if(findstring('nummol',record,idum))then
                
                nummols(itmols)=intstr(record,lenrec,idum)
                if(idnode.eq.0)
     x            write(nrite,"(/,1x,'number of molecules  ',
     x            10x,i10)")nummols(itmols)
                
              elseif(findstring('atoms',record,idum))then
                
c     read in atomic details
                
                call define_atoms
     x            (safe,lneut,lfreeze,idnode,itmols,nsite,ksite,ntpatm)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     read core - shell spring parameters
                
              elseif(findstring('shell',record,idum))then
                
                lshl=.true.
                call define_core_shell
     x            (safe,idnode,itmols,nshels,nsite,keyshl,
     x            engunit)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     read chemical bond force constant and bondlength
                
              elseif(findstring('bonds',record,idum))then
                
                call define_bonds
     x            (safe,idnode,itmols,nbonds,nsite,engunit)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

c     Assign void interactions to water model qtip4p/f to exclude
c     intramolecular interactions related to M-sites from Ewald sum
                
              elseif(findstring('voids',record,idum))then
                
                call define_voids
     x            (safe,idnode,mxnode,itmols,nvoids)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c *******************************************************************      

c     read bond atom indices and constraint bondlength
                
              elseif(findstring('constr',record,idum))then
                
                lcon=.true.
                call define_constraints
     x            (safe,lghost,idnode,itmols,nconst,nsite,natmsr)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     read pmf bond atom indices, weights and constraint bondlength
                
              elseif(findstring('pmf',record,idum))then
                
                if(lpmf)call error(idnode,484)
                lpmf=.true.
                call define_pmf(safe,idnode,itmols,nspmf)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     read intramolecular angular potential parameters
                
              elseif(findstring('angles',record,idum))then
                
                call define_angles
     x            (safe,idnode,itmols,nangle,nsite,engunit)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     read intramolecular dihedral potential parameters
                
              elseif(findstring('dihedr',record,idum))then
                
                call define_dihedrals
     x            (safe,idnode,itmols,ndihed,nsite,engunit)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     read intramolecular inversion potential parameters
                
              elseif(findstring('invers',record,idum))then
                
                call define_inversions
     x            (safe,idnode,itmols,ninver,nsite,engunit)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     read rigid body data
                
              elseif(findstring('rigid',record,idum))then
                
                lrig=.true.
                call define_rigid_body
     x            (safe,lghost,idnode,itmols,ngrp,natmsr)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     read tethered atom indices and tethering parameters
                
              elseif(findstring('teth',record,idum))then
                
                call define_tethers
     x            (safe,idnode,itmols,nteth,nsite,engunit)
                if(.not.safe)call abort_field_read(1,idnode,nfield)
                
c     finish of data for one molecular type
                
              elseif(findstring('finish',record,idum))then
                
c     running total of number of atoms in system
                
                natms=natms+nummols(itmols)*numsit(itmols)
                if(natms.gt.mxatms)call error(idnode,75)
                
c     check core-shell units are not both in same rigid body unit
                
                call check_shells(idnode,itmols,nshels,ngrp)
                
                loop2=.false.
                
              else
                
c     error exit for unidentified directive in molecular data
                
                if(idnode.eq.0)write(nrite,'(12x,a)')record
                call error(idnode,12)
                
              endif
              
            enddo
            
c     construction of atmolt table for solvation calculation
            
            if(lsolva)then
              
              do iatm=natmsr+1,natms
                atmolt(iatm)=itmols
              enddo
              natmsr=natms
              
            endif
            
          enddo
          
c     construction of atm_fre table for free energy or excitation
          
          if(lfree.or.lexcite.or.lswitch)then
            
            atm_fre(:)=0
            if((ind_fre(1).ne.0).and.(ind_fre(2).ne.0))then
              
              do iatm=ind_fre(1),ind_fre(2)
                atm_fre(iatm)=1
              enddo
              
            endif
            
            if((ind_fre(3).ne.0).and.(ind_fre(4).ne.0))then
              
              do iatm=ind_fre(3),ind_fre(4)
                
                atm_fre(iatm)=2
                if(lghost)ntghost=ntghost+1
                
              enddo
              
            endif
            
          endif
          
c     calculate system charge
          
          call check_syschg(idnode,ntpmls,sumchg)
          
c     read in the nonbonded potential energy parameters
          
        elseif(findstring('vdw',record,idum))then
          
          call define_van_der_waals
     x      (safe,ltable,lunits,lmols,idnode,ntpvdw,
     x      ntpatm,keyfce,dlrpot,rvdw,engunit)
          if(.not.safe)call abort_field_read(1,idnode,nfield)
          
c     read in the metal potential energy parameters
          
        elseif(findstring('met',record,idum))then
          
          call define_metals
     x      (safe,lunits,lmols,idnode,ntpmet,ntpatm,rvdw,engunit)
          if(.not.safe)call abort_field_read(1,idnode,nfield)
          
c     read the three body potential energy parameters
          
        elseif(findstring('tbp',record,idum))then
          
          call define_three_body
     x      (safe,lunits,lmols,idnode,ntptbp,ntpatm,rcuttb,engunit)
          if(.not.safe)call abort_field_read(1,idnode,nfield)
          
c     read the tersoff potential energy parameters
          
        elseif(findstring('tersoff',record,idum))then
          
          call define_tersoff
     x      (safe,lunits,lmols,idnode,ntpter,ntpatm,rctter,engunit)
          if(.not.safe)call abort_field_read(1,idnode,nfield)
          
c     read in the four body potential energy parameters
          
        elseif(findstring('fbp',record,idum))then
          
          call define_four_body
     x      (safe,lunits,lmols,idnode,ntpfbp,ntpatm,
     x      rcutfb,engunit)
          if(.not.safe)call abort_field_read(1,idnode,nfield)
          
c     read external field data
          
        elseif(findstring('extern',record,idum))then
          
          call define_external_field
     x      (safe,lunits,idnode,keyfld,engunit)
          if(.not.safe)call abort_field_read(1,idnode,nfield)
          
c     normal end of FIELD file
          
        elseif(findstring('close',record,idum))then
          
          loop1=.false.
          if(ntpvdw.eq.0.and.ntpmet.eq.0.and.
     x      mod(keyfce,2).eq.1)call error(idnode,145)
          
c     error exit for unidentified directive
          
        else
          
          if(idnode.eq.0)write(nrite,'(100a)')record
          call abort_field_read(2,idnode,nfield)
          
        endif
        
      enddo
      
c     close force field file
      
      if(idnode.eq.0)close (nfield)
      
      
      if(lshl.and.idnode.eq.0)then
        
        if(keyshl.eq.1)write(nrite,
     x    "(/,/,'adiabatic shell model in operation')")
        
        if(keyshl.eq.2)write(nrite,
     x    "(/,/,'relaxed shell model in operation')")
        
      endif
      
      if(lshl.and.keyshl.eq.0)call error(idnode,1951)
      
c     if metadynamics and shell selected use only velocity verlet
      
      if(lshl.and.lmetadyn.and.keyver.eq.0)then
        
        if(idnode.eq.0)write(nrite,
     x    "(/,1x,'switching to velocity verlet for metadynamics')")
        keyver=1
        
      endif
      
c     error exit if pimd incompatible options selected
      
      if(lpimd)then
        
         if(lshl)call error(idnode,-518)
         if(lcon)call error(idnode,-520)
         if(lrig)call error(idnode,-522)
         if(lneut)call error(idnode,-524)
         if(lmetadyn)call error(idnode,-525)
         if(lsolva)call error(idnode,-526)
         if(lfree)call error(idnode,-527)
         if(lfreeze)call error(idnode,-529)
         if(lshl.or.lcon.or.lrig.or.lneut.or.lmetadyn.or.lsolva.or.
     x     lfree.or.lfreeze)call error(idnode,0)
         
      endif
      
      return
      end subroutine sysdef
      
      subroutine sysgen
     x  (loglnk,lneut,nolink,lfree,lfrmas,lpimd,lmsite,idnode,imcon,
     x  keyens,keyfce,keyres,levcfg,multt,mxnode,ntpmls,nbeads,iqt4,
     x  g_qt4f,delr,rcut,temp,volm,uuu)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading the configuration data file
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith         july 1997
c     adapted   - p.-a. cazade oct 2007, solvation etc
c     
c***********************************************************************
      
      implicit none
      
      character*1 atname(8)
      
      logical loglnk,safe,lneut,nolink,lfree,lfrmas,lpimd,chkmat
      logical lmsite
      integer idnode,imcon,keyens,keyfce,keyres,levcfg,multt
      integer ntpmls,i,k,l,m,n,indatm,indnam,indneu,ilx,ily,ilz
      integer ncells,idum,mxnode,nbeads,numatm,mbeads,matms,iqt4
      real(8) g_qt4f
      real(8) delr,rcut,volm,xcoord,ycoord,zcoord,totmas,xveloc
      real(8) yveloc,zveloc,xforce,yforce,zforce,axx,rt3,xhi,yhi,zhi
      real(8) width,dum,dum1,dum2,test,disp,temp
      real(8) com(3),uuu(102)
      
c     open the system input file
      
      if(idnode.eq.0)open (nconf,file='CONFIG')
      
c     read the CONFIG file header
      
      call getrec(safe,idnode,nconf)
      if(.not.safe)call abort_config_read(1,idnode,nconf)
      
      call copystring(record,cfgname,80)
      if(idnode.eq.0)write(nrite,
     x  "(/,1x,'configuration file name: ',/,/,10x,80a1)")cfgname
      
      call getrec(safe,idnode,nconf)
      if(.not.safe)call abort_config_read(1,idnode,nconf)
      
      levcfg=intstr(record,lenrec,idum)
      imcon=intstr(record,lenrec,idum)
      
      if(idnode.eq.0)write(nrite,
     x  "(/,/,1x,'selected image convention',6x,i10)")imcon
      
c     check config file contents for consistent data
      
      if((imcon.eq.0.or.imcon.eq.6).and.
     x  (keyfce/2.eq.1.or.keyfce/2.eq.6))
     x  call error(idnode,180)
      
      if(imcon.eq.0.and.(.not.lneut).and.(keyfce.gt.1)
     x  .and.(multt.eq.1))call warning(idnode,30,0.d0,0.d0,0.d0)
      
      if(imcon.eq.0.and.(keyens.ge.4.and.keyens.le.7))
     x  call error(idnode,390)
      
      if(imcon.le.2.and.(keyens.eq.6.or.keyens.eq.7))imcon=3
      if(keyres.gt.0.and.levcfg.lt.1)call error(idnode,85)
      
c     specify molecular dynamics simulation cell
     
      if(imcon.eq.0)then
        
c     if no periodic boundaries - set zero values for cell 
c     vectors and cell volume
        
        cell(:)=0.d0
        volm=0.d0
        
      else
        
c     read cell vectors
        
        call getrec(safe,idnode,nconf)
        if(.not.safe)call abort_config_read(1,idnode,nconf)
        cell(1)=dblstr(record,lenrec,idum)
        cell(2)=dblstr(record,lenrec,idum)
        cell(3)=dblstr(record,lenrec,idum)
        call getrec(safe,idnode,nconf)
        if(.not.safe)call abort_config_read(1,idnode,nconf)
        cell(4)=dblstr(record,lenrec,idum)
        cell(5)=dblstr(record,lenrec,idum)
        cell(6)=dblstr(record,lenrec,idum)
        call getrec(safe,idnode,nconf)
        if(.not.safe)call abort_config_read(1,idnode,nconf)
        cell(7)=dblstr(record,lenrec,idum)
        cell(8)=dblstr(record,lenrec,idum)
        cell(9)=dblstr(record,lenrec,idum)
        
      endif
      
c     read the atomic coordinates
      
      indatm=0
      indnam=0
      indneu=0
      safe=.true.
      
c     site multiplicity factor for pimd

      numatm=nbeads*mxatms
c      write(6,*)"numatm",numatm 
c    restructure config read for pimd when keyres > 0

      if(keyres.eq.0.and.keyens.lt.61)then
        mbeads=1
        matms=mxatms
      else
        mbeads=nbeads
        matms=mxatms*nbeads
      endif

c     read atomic coordinates, velocities and forces
      
      do n=1,mbeads
        
        do k=1,ntpmls
          
          do l=1,nummols(k)
            
            do m=1,numsit(k)
              
              indatm=indatm+1
c              write(6,*)"indatm",indatm 
              if(indatm.gt.numatm)call error(idnode,45)
              
              xxx(indatm)=0.d0
              yyy(indatm)=0.d0
              zzz(indatm)=0.d0
              vxx(indatm)=0.d0
              vyy(indatm)=0.d0
              vzz(indatm)=0.d0
              fxx(indatm)=0.d0
              fyy(indatm)=0.d0
              fzz(indatm)=0.d0
              
              if(idnode.eq.0)then
                
                read(nconf,'(8a1)',end=100)atname
                read(nconf,'(3f20.0)',end=100)xcoord,ycoord,zcoord
                if(levcfg.gt.0)then
                  read(nconf,'(3f20.0)',end=100)xveloc,yveloc,zveloc
                endif
                if(levcfg.gt.1)then
                  read(nconf,'(3f20.0)',end=100)xforce,yforce,zforce
                endif
                
c     strip blanks off atom name
                
                call strip(atname,8)
                
                if(sitnam(indnam+m).eq.mkwd8(atname))then
                  
                  xxx(indatm)=xcoord
                  yyy(indatm)=ycoord
                  zzz(indatm)=zcoord

                  if(levcfg.gt.0)then
                    
                    vxx(indatm)=xveloc
                    vyy(indatm)=yveloc
                    vzz(indatm)=zveloc
                    
                  endif
                  
                  if(levcfg.gt.1)then
                    
                    fxx(indatm)=xforce
                    fyy(indatm)=yforce
                    fzz(indatm)=zforce
                    
                  endif
                  
                else
                  
                  write(nrite,"(/,/,'unidentified atom label :',8a1,
     x              ': atom number ',i5)")atname,indatm
                  safe=.false.
                  
                endif
                
              endif
              
              call gstate(safe)
              if(.not.safe)call error(idnode,25)
              
              ltype(indatm)=ltpsit(indnam+m)
              weight(indatm)=wgtsit(indnam+m)
              chge(indatm)=chgsit(indnam+m)
              atmnam(indatm)=sitnam(indnam+m)
              lstfrz(indatm)=lfzsit(indnam+m)
              if(lneut)lstneu(indatm)=nugrp(indnam+m)+indneu
              
c     reset atomic masses according to free energy definitions
              
              if(lfree)then
                
                weight_sav(indatm)=weight(indatm)
                
                if(lfrmas)then
                  
                  if(indatm.ge.ind_fre(1).and.indatm.le.ind_fre(2))then
                    weight(indatm)=lambda1*weight(indatm)
                  elseif(indatm.ge.ind_fre(3).and.indatm.le.ind_fre(4))
     x                then
                    weight(indatm)=lambda2*weight(indatm)
                  endif
                  
                endif
                
              endif
              
            enddo
            
            indneu=indneu+nugrp(indnam+numsit(k))
            
          enddo
          
          indnam=indnam+numsit(k)
          
        enddo

        indnam=0

      enddo

c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

c     Initialize the position of M-site based on O and H for 
c     qtip4p/f water model 

      if(lmsite)then

        call qtip4pf(idnode,mxnode,imcon,nbeads,ntpmls,g_qt4f)

      endif

c ******************************************************************  

      if(mxnode.gt.1)then
        
        call gdsum(xxx,matms,buffer)
        call gdsum(yyy,matms,buffer)
        call gdsum(zzz,matms,buffer)
        
        if(levcfg.gt.0)then
          
          call gdsum(vxx,matms,buffer)
          call gdsum(vyy,matms,buffer)
          call gdsum(vzz,matms,buffer)
          
        endif
        
        if(levcfg.gt.1)then
          
          call gdsum(fxx,matms,buffer)
          call gdsum(fyy,matms,buffer)
          call gdsum(fzz,matms,buffer)
          
        endif
        
      endif
c      write(6,*) "xxx", xxx(:) 
c     for pimd expand initial atomic system to quantum system
      
      if(lpimd.and.keyres.eq.0.and.keyens.lt.61)then
        
        m=indatm*nbeads+1
        
        do n=1,nbeads
          
          do i=indatm,1,-1
            
            m=m-1
            if(m.gt.numatm)call error(idnode,45)
            disp=1.d-3**(hbar/sqrt(boltz*temp*weight(i)))/dble(nbeads)
            
c     spread beads to fraction of de broglie wavelength
            
            xxx(m)=xxx(i)+disp*(2.d0*duni()-1.d0)
            yyy(m)=yyy(i)+disp*(2.d0*duni()-1.d0)
            zzz(m)=zzz(i)+disp*(2.d0*duni()-1.d0)
            
            fxx(m)=fxx(i)
            fyy(m)=fyy(i)
            fzz(m)=fzz(i)
            ltype(m)=ltype(i)
            weight(m)=weight(i)
            chge(m)=chge(i)
            atmnam(m)=atmnam(i)
            
          enddo
          
        enddo
        
      endif

c     initialise parallel random number sequence

      dum=puni(1,uuu)
      
c     check integrity of cell vectors : for cubic, TO and RD cases
c     ie. cell(1)=cell(5)=cell(9) (or cell(9)/sqrt(2) for RD)
      
      if((imcon.eq.1).or.(imcon.eq.4).or.(imcon.eq.5))then
        
        axx=(abs(cell(1))+abs(cell(5)))/2.d0
        test=1.d-8*axx
        if(abs(cell(1)-axx).gt.test)call error(idnode,410)
        if(abs(cell(5)-axx).gt.test)call error(idnode,410)
        if(imcon.eq.5)then
          if(abs(cell(9)-axx*sqrt(2.d0)).gt.test)
     x      call error(idnode,410)
        else
          if(abs(cell(9)-axx).gt.test)call error(idnode,410)
        endif
        
      endif
      
c     check integrity of hexagonal prism cell vectors
      
      if(imcon.eq.7)then
        
        rt3=sqrt(3.d0)
        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
     x    call error(idnode,410)
        
      endif
      
c     check 2D PBC for imcon=6
      
      if(imcon.eq.6)then
        
        if(abs(cell(3)).gt.1.d-10)call error(idnode,410)
        if(abs(cell(6)).gt.1.d-10)call error(idnode,410)
        if(abs(cell(7)).gt.1.d-10)call error(idnode,410)
        if(abs(cell(8)).gt.1.d-10)call error(idnode,410)
        
      endif
      
c     check for diagonal cell matrix if appropriate
      
      if((imcon.eq.1).or.(imcon.eq.2).or.(imcon.eq.4).or.
     x  (imcon.eq.5).or.(imcon.eq.7))then
        
        chkmat=.false.
        if(abs(cell(2)).gt.1.d-10)chkmat=.true.
        if(abs(cell(3)).gt.1.d-10)chkmat=.true.
        if(abs(cell(4)).gt.1.d-10)chkmat=.true.
        if(abs(cell(6)).gt.1.d-10)chkmat=.true.
        if(abs(cell(7)).gt.1.d-10)chkmat=.true.
        if(abs(cell(8)).gt.1.d-10)chkmat=.true.
        if(chkmat)call error(idnode,410)
        
      endif
      
c     put centre of mass at centre of coordinates if imcon=0
      
      if(imcon.eq.0)then
        
        totmas=getmass(numatm,idnode,mxnode)
        call getcom(numatm,idnode,mxnode,totmas,com)
        
        do i=1,numatm
          
          xxx(i)=xxx(i)-com(1)
          yyy(i)=yyy(i)-com(2)
          zzz(i)=zzz(i)-com(3)
          
        enddo
        
      endif
      
c     set widths if unset - needed for check on link cells below
      
      if(imcon.eq.0.or.imcon.eq.6)then
        
        xhi=abs(xxx(1))
        yhi=abs(yyy(1))
        zhi=abs(zzz(1))
        do i=2,numatm
          
          xhi=max(xhi,abs(xxx(i)))
          yhi=max(yhi,abs(yyy(i)))
          zhi=max(zhi,abs(zzz(i)))
          
        enddo
        if(imcon.eq.0)then
          
          cell(1)=max(2.d0*xhi+rcut+delr,3.d0*(rcut+delr))
          cell(5)=max(2.d0*yhi+rcut+delr,3.d0*(rcut+delr))
          cell(9)=max(2.d0*zhi+rcut+delr,3.d0*(rcut+delr))
          
        endif
        
        if(imcon.eq.6.and.cell(9).lt.1.d-6)then
          
          cell(9)=max(2.d0*zhi+rcut+delr,3.d0*(rcut+delr))
          
        endif
        
      endif
      
c     calculate dimensional properties of simulation cell
      
      call dcell(cell,celprp)
      
      if(imcon.eq.0)then
        
        volm=0.d0
        
      elseif(imcon.eq.4)then
        
        volm=0.5d0*celprp(10)
        
      elseif(imcon.eq.5)then
        
        volm=0.5d0*celprp(10)
        
      elseif(imcon.eq.7)then
        
        volm=0.5d0*celprp(10)
        
      else
        
        volm=celprp(10)
        
      endif
      
      if(idnode.eq.0)then
        
        write(nrite,"(/,/,1x,'simulation cell vectors'/,/)")
        write(nrite,"(21x,3f12.6)")cell
        
        write(nrite,
     x    "(/,/,1x,'system volume     ',2x,1p,g22.12)")volm
        
      endif
      
c     check value of cutoff and reset if necessary
      
      if(imcon.gt.0)then
        
        width=min(celprp(7),celprp(8),celprp(9))/2.d0
        if(imcon.eq.4)width=sqrt(3.d0)*cell(1)/4.d0
        if(imcon.eq.5)width=cell(1)/2.d0
        if(imcon.eq.6)width=min(celprp(7),celprp(8))/2.d0
        
c     halt program if potential cutoff exceeds cell width
        
        if(rcut.gt.width)call error(idnode,95)
        
      endif
      
c     decide on whether to use link cells for verlet list constructor
      
      if(nolink)then
        
        loglnk=.false.
        
      else
        
        loglnk=.true.
        
        ilx=int(celprp(7)/(rcut+delr))
        ily=int(celprp(8)/(rcut+delr))
        ilz=int(celprp(9)/(rcut+delr))
        if(ilx.lt.3.or.ily.lt.3.or.ilz.lt.3)loglnk=.false.
        ncells=ilx*ily*ilz
        if(lneut.and.ncells.le.36)loglnk=.false.
        if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)loglnk=.false.
        if(loglnk.and.ncells.gt.mxcell)then
          
          dum1=dble(ncells)
          dum2=dble(mxcell)
          call warning(idnode,90,dum1,dum2,dum2)
          loglnk=.false.
          
        endif
        
      endif
      
      if(loglnk.and.idnode.eq.0)
     x  write(nrite,"(/,/,' link cell algorithm in use')")
      
      if(idnode.eq.0)close (nconf)
      
c     ensure PBC compliance of starting structure
      
      if(keyres.eq.0.and.imcon.gt.0)then
        
        call images(imcon,idnode,mxnode,numatm,cell,xxx,yyy,zzz)
        if(mxnode.gt.1)
     x    call merge(idnode,mxnode,numatm,mxbuff,xxx,yyy,zzz,buffer)
        
      endif
      
      return
      
c     error exit for config file read
      
  100 call abort_config_read(2,idnode,nconf)
      
      end subroutine sysgen
      
      subroutine sysinit
     x  (lgofr,lzden,lsolva,lfree,lghost,lpsoc,lpimd,idnode,imcon,
     x  keyfce,keyres,mxnode,natms,nbeads,ntshl,nstep,numacc,numrdf,
     x  ntpatm,ntpmet,ntpvdw,nzden,chip,chit,conint,elrc,engunit,
     x  virlrc,rvdw,volm,virtot,vircom,tboost,chit_shl,engthe,gaumom)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading the REVIVE file data and 
c     defining the initial thermodynamic and structural accumulators.
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith         july 1997
c     adapted   - p.-a. cazade oct 2007, solvation etc
c     adapted   - d. quigley nov 2010, metadynamics
c     
c***********************************************************************
      
      implicit none
      
      logical lgofr,lzden,lfree,lsolva,lghost,lpsoc,lpimd
      integer idnode,imcon,keyfce,keyres,mxnode,natms,nbeads,nstep
      integer numacc,numrdf,ntpatm,nzden,i,j,k,ntpmet,ntshl,ntpvdw
      real(8) chip,chit,conint,elrc,engunit,virlrc,rvdw,volm
      real(8) dnumrd,dnstep,dnumac,dnzden,virtot,vircom,tboost
      real(8) chit_shl,engthe
      real(8) gaumom(0:5)
      
      engthe=0.d0
      
c     read or initialise accumulator arrays
      
      if(keyres.eq.1.and.idnode.eq.0)then
        
c     read accumulator data from dump file
        
        open(nrest,file='REVOLD',form='unformatted')
        
        read(nrest) dnstep,dnumac,dnumrd,chit,chip,conint,dnzden,
     x    tboost,chit_shl
        read(nrest) virtot,vircom,eta,strcns,strbod
        read(nrest) stpval
        read(nrest) sumval
        read(nrest) ssqval
        read(nrest) zumval
        read(nrest) ravval
        read(nrest) stkval
        read(nrest) xx0,yy0,zz0
        read(nrest) xxs,yys,zzs
        
        if(lgofr) read(nrest)rdf
        if(lzden) read(nrest)zdens
        read(nrest)gaumom
        
        nstep=nint(dnstep)
        numacc=nint(dnumac)
        numrdf=nint(dnumrd)
        nzden=nint(dnzden)
        close (nrest)
        
      else
         
c     initialise step counters
        
        nstep=0
        numacc=0
        numrdf=0
        nzden=0
        
c     initialise temperature and pressure coupling parameters
c     and integral for conserved quantity
        
        chit=0.d0
        chip=0.d0
        conint=0.d0
        virtot=0.d0
        vircom=0.d0
        chit_shl=0.d0
        do i=1,9
          
          eta(i)=0.d0
          strcns(i)=0.d0
          strbod(i)=0.d0
          
        enddo

c     initialise bias potential boost factor
        
        tboost=0.d0
        
c     initialise accumulator arrays
        
        do i=1,mxnstk
          
          stpval(i)=0.d0
          sumval(i)=0.d0
          ssqval(i)=0.d0
          zumval(i)=0.d0
          ravval(i)=0.d0
          
        enddo
        
        do i=1,mxatms
          
          xx0(i)=0.d0
          yy0(i)=0.d0
          zz0(i)=0.d0
          xxs(i)=0.d0
          yys(i)=0.d0
          zzs(i)=0.d0
          
        enddo
        
        do j=1,mxnstk
          do i=1,mxstak
            
            stkval(i,j)=0.d0
            
          enddo
        enddo
        
        if(lgofr)then
          
          do i=1,mxxtyp
            do j=1,mxrdf
              
              rdf(j,i)=0.d0
              
            enddo
          enddo
          
        endif
        
        if(lzden)then
          
          do i=1,mxatyp
            do j=1,mxzdn
              
              zdens(j,i)=0.d0
              
            enddo
          enddo
          
        endif

        do i=0,5
          
          gaumom(i)=0.d0
          
        enddo
        
      endif
      
c     put shells on cores at start
      
      if(lpsoc.and.keyres.ne.1.and.ntshl.gt.0)
     x  call put_shells_on_cores(idnode,mxnode,ntshl)
      
c     if restart then broadcast stored variables via a global sum
      
      if(keyres.eq.1.and.mxnode.gt.1)then
        
        if(mxbuff.lt.natms.or.mxbuff.lt.mxnstk*mxstak)
     x    call error(idnode,186)
        
        buffer(1)=chit
        buffer(2)=chip
        buffer(3)=conint
        buffer(4)=dble(nstep)
        buffer(5)=dble(numacc)
        buffer(6)=dble(numrdf)
        buffer(7)=dble(nzden)
        buffer(8)=tboost
        buffer(9)=virtot
        buffer(10)=vircom
        buffer(11)=chit_shl
        call gdsum(buffer(1),11,buffer(12))
        chit=buffer(1)
        chip=buffer(2)
        conint=buffer(3)
        nstep=nint(buffer(4))
        numacc=nint(buffer(5))
        numrdf=nint(buffer(6))
        nzden=nint(buffer(7))
        tboost=buffer(8)
        virtot=buffer(9)
        vircom=buffer(10)
        chit_shl=buffer(11)
        
        call gdsum(eta,9,buffer)
        call gdsum(strcns,9,buffer)
        call gdsum(strbod,9,buffer)
        call gdsum(stpval,mxnstk,buffer)
        call gdsum(sumval,mxnstk,buffer)
        call gdsum(ssqval,mxnstk,buffer)
        call gdsum(zumval,mxnstk,buffer)
        call gdsum(ravval,mxnstk,buffer)    
        call gdsum(stkval,mxnstk*mxstak,buffer)
        call gdsum(xx0,natms,buffer)
        call gdsum(yy0,natms,buffer)
        call gdsum(zz0,natms,buffer)
        call gdsum(xxs,natms,buffer)
        call gdsum(yys,natms,buffer)
        call gdsum(zzs,natms,buffer)
        
c     for rdf table - broadcast and normalise
        
        if(lgofr)then
          
          do k=1,mxxtyp
            call gdsum(rdf(1,k),mxrdf,buffer)
            do j=1,mxrdf
              
              rdf(j,k)=rdf(j,k)/dble(mxnode)
              
            enddo
          enddo
          
        endif
        
        if(lzden)then
          
          do k=1,mxatyp
            call gdsum(zdens(1,k),mxzdn,buffer)
            do j=1,mxzdn
              
              zdens(j,k)=zdens(j,k)/dble(mxnode)
              
            enddo
          enddo
          
        endif

c     broadcast and normalise gaussian moments
        
        call gdsum(gaumom(0),6,buffer(1))
        
        gaumom(0)=gaumom(0)*dble((((idnode+1)*natms)/mxnode)-
     x    (idnode*natms)/mxnode)/dble(natms)
        
      endif
      
c     number densities and long-range corrections
      
      elrc=0.d0       
      virlrc=0.d0
      
      if(imcon.eq.0.or.imcon.eq.6)volm=4.d0*pi/3.d0*rvdw**3
      
      call lrcorrect
     x  (lsolva,lfree,lghost,idnode,imcon,keyfce,natms,
     x  ntpatm,ntpvdw,elrc,engunit,virlrc,rvdw,volm)
      
      if(lmetab.or.ntpmet.eq.0)then
        
        elrcm(0)=0.d0
        vlrcm(0)=0.d0
        
      else 
        
        call lrcmetal
     x    (idnode,imcon,natms,ntpatm,engunit,rvdw,volm)
        
      endif
      
      if(imcon.eq.0.or.imcon.eq.6)volm=0.d0

      return
      end subroutine sysinit
      
      subroutine systemp
     x  (lpimd,inhc,idnode,imcon,keyres,mxnode,natms,nbeads,ngrp,nscons,
     x  ntcons,ntfree,ntshl,levcfg,keyshl,keyens,degfre,degshl,nchain,
     x  degrot,engke,tolnce,temp,sigma,sigma_nhc,sigma_volm,alpha_volm,
     x  uuu,tstep)
      
c***********************************************************************
c     
c     dl_poly subroutine for setting the initial system temperature
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith         july 1997
c     
c***********************************************************************
      
      implicit none

      logical lpimd,inhc
      integer idnode,imcon,keyres,mxnode,natms,nbeads,ngrp,nscons
      integer ntcons,ntfree,ntshl,levcfg,i,io,k,keyshl,keyens,nchain
      real(8) degfre,degshl,degrot,tolnce,temp,sigma,engke,rsq,tstep
      real(8) sigma_nhc,eta_nhc(nchain),peta(nchain)
      real(8) sigma_volm,alpha_volm,ksi(nchain),pksi(nchain)
      real(8) uuu(102)
      
c     number of degrees of freedom 
c     3 for com translation
c     3 for angular momentum about origin (non-periodic systems only)
     
      degfre=dble(3*(ntfree-ntshl)-3-ntcons)+degfre
      if(imcon.eq.0.or.imcon.eq.6)degfre=degfre-3.0d0
      if(imcon.eq.0.or.imcon.eq.6)degrot=max(0.d0,degrot-3.0d0)
      degshl=dble(3*ntshl)
      
c$$$     lose one degree of freedom if temperature constrained
c$$$     gaussian constraints
c$$$     if(keyens.eq.1)degfre=degfre-1.d0
      
      if(idnode.eq.0)
     x  write(nrite,"(/,/,' total degrees of freedom       ',f20.0,/,
     x  ' rotational degrees of freedom  ',f20.0,/,
     x  ' shell pseudo degrees of freedom',f20.0)")
     x  degfre,degrot,degshl
      if(degfre.lt.1.d0)call error(idnode,350)

c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

c     set temperature control parameter
      
      sigma=temp*boltz*degfre*0.5d0
      sigma_nhc=temp*boltz*0.5d0
      sigma_volm=temp*boltz*(degfre+1)
      alpha_volm=1.d0+(1.d0/degfre)

      if(inhc)then
        
c     initialise nhc variables
      
      call nhc_init(idnode,mxnode,nchain)
      
      endif

c *******************************************************************      
      
      if(lpimd)then
       
c     initialise pimd simulations
      
        call pimd_init
     x    (idnode,mxnode,natms,keyres,keyens,temp,sigma,engke,
     x    stress,uuu,tstep)

      else

c     initialise classical md simulations
        
        do i=1,natms*nbeads
          
          if(lstfrz(i).ne.0.or.weight(i).lt.1.d-6)then
            
            rmass(i)=0.d0
            weight(i)=0.d0
            
          else
            
            rmass(i)=1.d0/weight(i)
            
          endif
          
        enddo
        
c     generate starting velocities
        
        if(keyres.eq.0)then
          
          call gauss(natms*nbeads,vxx,vyy,vzz)
          
          do i=1,natms*nbeads
            
            rsq=sqrt(rmass(i))
            vxx(i)=vxx(i)*rsq
            vyy(i)=vyy(i)*rsq
            vzz(i)=vzz(i)*rsq
            
          enddo
          
          if(ntcons.gt.0)call quench
     x      (imcon,idnode,mxnode,natms,nscons,tolnce)
          
          if(ngrp.gt.0)call quatqnch(idnode,imcon,mxnode,natms,ngrp)
          
          if(keyshl.eq.1)then
            
            do k=1,4
              
              call vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)
              call shlqnch(idnode,mxnode,ntshl,temp)
              
            enddo
            
          else
            
            call vscaleg(idnode,mxnode,imcon,natms*nbeads,ngrp,sigma)
            
          endif
          
        elseif(keyres.eq.1.or.keyres.eq.3)then 
          
          if(ngrp.gt.0)call quatqnch(idnode,imcon,mxnode,natms,ngrp)
          
        elseif(keyres.eq.2)then
          
          if(ngrp.gt.0)then 
            
            call vscaleg
     x        (idnode,mxnode,imcon,natms,ngrp,sigma)
            
          elseif(keyshl.eq.1)then
            
            do k=1,4
              
              call vscaleg(idnode,mxnode,imcon,natms,ngrp,sigma)
              call shlqnch(idnode,mxnode,ntshl,temp)
              
            enddo
            
          else
            
            call vscaleg(idnode,mxnode,imcon,natms*nbeads,ngrp,sigma)
            
          endif
          
        endif

c     initial system kinetic energy and stress
        
        call kinstress(natms,idnode,mxnode,stress)
        engke=0.5d0*(stress(1)+stress(5)+stress(9))
        do i=1,9
          stress(i)=stress(i)/dble(mxnode)
        enddo
        
      endif
      
c     print out sample of initial configuration 
      
      if(idnode.eq.0)write(nrite,
     x  "(/,/,1x,'sample of starting configuration',/)")
      
      io=(natms*nbeads+19)/20
      if((levcfg.le.1).and.(idnode.eq.0))
     x  write(nrite,"(6x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)',
     x  7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',/,/)")
      if((levcfg.eq.2).and.(idnode.eq.0))
     x  write(nrite,"(6x,'i',7x,'x(i)',8x,'y(i)',8x,'z(i)',
     x  7x,'vx(i)',7x,'vy(i)',7x,'vz(i)',
     x  7x,'fx(i)',7x,'fy(i)',7x,'fz(i)',/,/)")
      
      do i=1,natms*nbeads,io
        
        if(levcfg.le.1)then
          
          if(idnode.eq.0)write(nrite,
     x      "(1x,i6,1p,3e12.4,3e12.4,3e12.4)")
     x      i,xxx(i),yyy(i),zzz(i),vxx(i),vyy(i),vzz(i)
          
        elseif(levcfg.eq.2)then
          
          if(idnode.eq.0)write(nrite,
     x      "(1x,i6,1p,3e12.4,3e12.4,3e12.4)")
     x      i,xxx(i),yyy(i),zzz(i),
     x      vxx(i),vyy(i),vzz(i),fxx(i),fyy(i),fzz(i)
          
        endif
        
      enddo

c     write out normal mode frequencies for testing - Nathan London
      if(keyens.ge.43)then
        if(idnode.eq.0)write(nrite,
     x    "(/,/,1x,'normal mode frequencies',/)")
        do i=1,nbeads
          if(idnode.eq.0)write(nrite,
     x      "(1e12.4)")nmfreq(i)
        enddo
      endif 
      return
      end subroutine systemp
      
      subroutine sysbook
     x  (loglnk,lneut,lshmov,lcnb,lsolva,lghost,lmsite,idnode,
     x  imcon,mxnode,natms,nneut,ngrp,nscons,ntangl,ntbond,ntcons,
     x  ntdihd,ntinv,ntpmls,ntpmf,nspmf,ntfree,ntteth,ntshl,
     x  ntghost,degfre,degrot)
      
c***********************************************************************
c     
c     dl_poly subroutine  defining global bookkeeping
c     arrays
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith         july 1997
c     adapted   - p.-a. cazade oct 2007, solvation etc
c     
c***********************************************************************
      
      implicit none
      
      logical loglnk,lneut,lshmov,lcnb,lsolva,lghost,lmsite
      integer idnode,imcon,mxnode,natms,nneut,ngrp,nscons,ntangl,i
      integer ntbond,ntcons,ntdihd,ntinv,ntpmls,ntpmf,nspmf,ntfree
      integer ntteth,ntshl,ii,isol,itmols,igsol,iggsol,natmsf,natmsl
      integer ntghost,natms2
      real(8) degfre,degrot
      
c     if excitation calculation, allow for ghost species
      
      if(lghost)then
        natms2=natms-ntghost
      else
        natms2=natms
      endif
      
c     neutral group bookkeeping
      
      if(lneut)call neutbook(lneut,idnode,natms,nneut)
      
c     rigid body bookkeeping 
      
      call quatbook
     x  (lsolva,idnode,imcon,mxnode,natms2,ngrp,ntpmls,
     x  ntfree,degfre,degrot)
      
c     if excitation calculation, allow for ghost species
      
      if(lghost)then
        
        numcon(mxtmls)=numcon(mxtmls)+ntcons_ghost
        numgrp(mxtmls)=numgrp(mxtmls)+ngrp_ghost
        
      endif
      
c     construct list of excluded pair interactions
      
      if(lneut)then
 
        call exclude(idnode,mxnode,natms,ntpmls,lmsite)
        call excludeneu(idnode,mxnode,nneut)
        
      elseif(.not.lneut)then
        
        call exclude(idnode,mxnode,natms,ntpmls,lmsite)
        
        if(loglnk)then
          
          call exclude_link(idnode,mxnode,ntpmls)
          
        else
          
          call exclude_atom(idnode,mxnode,natms,ntpmls)
          
        endif
        
      endif
      
c     if excitation calculation, allow for ghost species
      
      if(lghost)then
        
        numcon(mxtmls)=numcon(mxtmls)-ntcons_ghost
        numgrp(mxtmls)=numgrp(mxtmls)-ngrp_ghost
        
      endif
      
c     construct interaction lists for bonded forces
      
      call intlist
     x  (lshmov,lcnb,idnode,mxnode,natms2,nscons,ntangl,ntbond,
     x  ntcons,ntdihd,ntinv,ntpmls,ntteth,ntshl,ntpmf,nspmf,ngrp)
      
c     adaptations for solvation and excitation simulations
      
      if(lsolva.or.lghost)then
                
        natmsf=0
        natmsl=0
        natm_sol(:)=0
        const_sol(:)=numcon(:)*nummols(:)
        rigid_sol(:)=numgrp(:)*nummols(:)
        
        if(ngrp.eq.0)then
          
          do itmols=1,mxtmls
            
            natmsl=natmsl+numsit(itmols)*nummols(itmols)
            
            do isol=natmsf+1,natmsl
              
              if(lstfrz(isol).eq.0)then
                natm_sol(itmols)=natm_sol(itmols)+1
              endif
              
            enddo
            
            natmsf=natmsl
            
          enddo

        else
          
          ii=1
          
          do itmols=1,mxtmls
            
            natmsl=natmsl+numsit(itmols)*nummols(itmols)
            
            do isol=natmsf+1,natmsl
              
              if(lstgot_sol(ii).eq.isol)then
                ii=ii+1
              else
                
                if(lstfrz(isol).eq.0)then
                  natm_sol(itmols)=natm_sol(itmols)+1
                endif
                
              endif
              
            enddo
            
            natmsf=natmsl
            
          enddo
          
          degrot_sol(:)=degrot_sol(:)+dble(rigid_sol(:))*3.d0
          degfre_sol(:)=degrot_sol(:)+dble(rigid_sol(:))*3.d0
          
        endif
        
        if(lghost)natm_sol(mxtmls)=natm_sol(mxtmls)-ntghost
        degfre_sol(:)=dble(3*(natm_sol(:))-const_sol(:))+degfre_sol(:)
        
      endif
      
      return
      end subroutine sysbook
      
      subroutine define_units(idnode,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine for selecting energy units
c     
c     copyright - daresbury laboratory 
c     author    - w. smith august 2003
c     
c***********************************************************************
      
      implicit none
      
      integer idnode,idum,i
      real(8) engunit
      logical blank
      
      blank=.true.
      
      do i=6,lenrec
        if(record(i).ne.' ')blank=.false.
      enddo
      
      if(blank)then
        
        if(idnode.eq.0)
     x    write(nrite,"(/,' energy units=dl_poly internal ',
     x    'units ')")
        
      elseif(findstring('ev',record,idum))then
        
        engunit=9648.530821d0
        if(idnode.eq.0)
     x    write(nrite,"(/,' energy units=electron volts ')")
        
      elseif(findstring('kev',record,idum))then
        
        engunit=9648530.821d0
        if(idnode.eq.0)
     x    write(nrite,"(/,' energy units=kilo electron volts ')")
        
      elseif(findstring('kcal',record,idum))then
        
        engunit=418.4d0
        if(idnode.eq.0)
     x    write(nrite,"(/,' energy units=kcal/ mol ')")
        
      elseif(findstring('kj',record,idum))then
        
        engunit=1.d2
        if(idnode.eq.0)
     x    write(nrite,"(/,' energy units=kjoule/mol ')")
        
      elseif(findstring('k',record,idum))then
        
        engunit=boltz
        if(idnode.eq.0)
     x    write(nrite,"(/,' energy units=kelvin ')")
        
      elseif(findstring('internal',record,idum))then
        
        if(idnode.eq.0)
     x    write(nrite,"(/,' energy units=dl_poly internal',
     x    ' units ')")
        
      else
        
        if(idnode.eq.0)write(nrite,'(a)')record
        call error(idnode,5)
        
      endif
      
      return
      end subroutine define_units
      
      subroutine quatbook
     x  (lsolva,idnode,imcon,mxnode,natms,ngrp,ntpmls,ntfree,degfre,
     x  degrot)
      
c**************************************************************************
c     
c     dl_poly subroutine for setting up bookkeeping for rigid bodies
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1993
c     author      t.forester october 1993
c     amended     t.forester dec 1994 : block data
c     adapted   - p.-a. cazade oct 2007, solvation etc
c     
c*************************************************************************
      
      implicit none
      
      logical safe,pass1,pass2,linear,lsolva
      integer fail,idnode,imcon,mxnode,natms,ngrp,ntpmls,ntfree
      integer i,igrp,jgrp,kgrp,jr,jt,igrp1,igrp2,itmols,imols,lgrp,id
      integer ii,jj,isite,k,kk,ill,i1,i2,i3,j,ngp,ifre1,ifre2,ig,ij
      integer fngrp,lngrp
      real(8) degfre,degrot,dnorm,a1,rtall,rotall,rot,aa,rotinr,bb,rot1
      real(8) rsq,det,dettest,aq,bq,cq,dq,eq,fq,gq,hq,rnorm,tol,rotxyz
      real(8) rotlim,rrr
      
      integer, allocatable :: ind(:,:),lstgot(:)
      real(8), allocatable :: gaxs(:,:),rotmin(:),accum(:)
      real(8), allocatable :: txx(:),tyy(:),tzz(:)
      real(8), allocatable :: xxt(:),yyt(:),zzt(:)
      
      dimension rot(9),aa(9),rotinr(3,3),bb(9),rot1(3,3),fail(5)
      
      data fail/0,0,0,0,0/
      
c     allocate working arrays
      
      allocate (ind(mxgrp,3),lstgot(mxatms),stat=fail(1))
      allocate (txx(mxatms),tyy(mxatms),tzz(mxatms),stat=fail(2))
      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail(3))
      allocate (gaxs(mxungp,9),rotmin(mxungp),stat=fail(4))
      allocate (accum(mxungp),stat=fail(5))
      do i=1,5
        if(fail(i).ne.0)call error(idnode,1790)
      enddo
      
c     initialise bookkeeping indices
      
      igrp=0
      jgrp=0
      kgrp=0
      isite=0
      jr=0
      jt=0
      safe=.true.
      degfre=0.d0
      degrot=0.d0
      
c     rigid body identifier
      
      do i=1,natms
        lstbod(i)=0
      enddo
      
c     number of rigid groups in system
      
      ngrp=0
      do itmols=1,ntpmls
        ngrp=ngrp+nummols(itmols)*numgrp(itmols)
      enddo
      
c     block indices for groups
      
      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode
      
c     loop over molecule types
      
      do itmols=1,ntpmls
        
c     loop over molecules in system
        
        do imols=1,nummols(itmols)
          
c     construct rigid body site list: each processor has a different copy
          
          do lgrp=1,numgrp(itmols)
            
            igrp=igrp+1
            
            if(igrp.le.mxgrp)then
              
              lstgtp(igrp)=listyp(lgrp+kgrp)
              id=listyp(lgrp+kgrp)
              
              if((igrp.ge.igrp1).and.(igrp.le.igrp2))then
                
                jgrp=jgrp+1
                
                do jj=1,numgsit(id)
                  
                  jr=jr+1
                  jt=jt+1
                  
                  if(jr.le.mxatms.and.jt.le.mxatms)then
                    
                    lstrgd(jr)=lstgst(id,jj)+isite
                    lstgot(jt)=lstgst(id,jj)+isite
                    lstbod(lstgst(id,jj)+isite)=igrp
                    
                  else
                    
                    safe=.false.
                    
                  endif
                  
                enddo
                
              else
                
                do jj=1,numgsit(id)
                  
                  jt=jt+1
                  if(jt.le.mxatms)then
                    
                    lstgot(jt)=lstgst(id,jj)+isite
                    lstbod(lstgst(id,jj)+isite)=igrp
                    
                  else
                    
                    safe=.false.
                    
                  endif
                  
                enddo
                
              endif
              
            else
              
              safe=.false.
              
            endif
            
          enddo
          
          if(mxnode.gt.1)call gstate(safe)
          if(.not.safe)call error(idnode,304)
          isite=isite+numsit(itmols)
          
        enddo
        
        kgrp=kgrp+numgrp(itmols)
        
      enddo
      
      if(ngrp.eq.0)then
        
        j=0
        do i=1,natms
          
          if(lstfrz(i).eq.0)then
            
            j=j+1
            lstfre(j)=i
            
          endif
          
        enddo
        ntfree=j
        
      else
        
c     centre of mass of groups
c     assumes group dimensions are smaller than half box width
        
        do i=1,natms
          
          lstme(i)=0
          
        enddo
        
        do id=1,mxungp
          
          gmass(id)=0.d0
          
        enddo
        
        jr=0
        do ig=igrp1,igrp2
          
c     working com is first site in group
          
          i=lstrgd(jr+1)
          txx(ig)=xxx(i)
          tyy(ig)=yyy(i)
          tzz(ig)=zzz(i)
          
          id=lstgtp(ig)
          safe=.false.
          if(abs(gmass(id)).lt.1.d-10)safe=.true.
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            xxt(jr)=xxx(i)-txx(ig)
            yyt(jr)=yyy(i)-tyy(ig)
            zzt(jr)=zzz(i)-tzz(ig)
            if(safe)gmass(id)=gmass(id)+weight(i)
            
          enddo
          
        enddo
        
c     minimum image from working com
        
        call images(imcon,0,1,jr,cell,xxt,yyt,zzt)
        
        jr=0
        do ig=igrp1,igrp2
          
          gcmx(ig)=0.d0
          gcmy(ig)=0.d0
          gcmz(ig)=0.d0
          
          id=lstgtp(ig)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            gcmx(ig)=gcmx(ig)+weight(i)*xxt(jr)
            gcmy(ig)=gcmy(ig)+weight(i)*yyt(jr)
            gcmz(ig)=gcmz(ig)+weight(i)*zzt(jr)
            
          enddo
          
          gcmx(ig)=gcmx(ig)/gmass(id)+txx(ig)
          gcmy(ig)=gcmy(ig)/gmass(id)+tyy(ig)
          gcmz(ig)=gcmz(ig)/gmass(id)+tzz(ig)
          
        enddo
        
c     global communications
        
        if(mxnode.gt.1)then
          
          call merge(idnode,mxnode,ngrp,mxbuff,gcmx,gcmy,gcmz,buffer)
          
        endif
        
c     make sure all nodes have same copy of gmass
        
        if(mxnode.gt.1)then
          
          do id=1,mxungp
            
            accum(id)=0.d0
            if(gmass(id).gt.0.d0)accum(id)=1.d0
            
          enddo
          
          call gdsum(gmass(1),mxungp,buffer(1))
          call gdsum(accum(1),mxungp,buffer(1))
          
          do id=1,mxungp
            
            dnorm=max(1.d0,accum(id))
            gmass(id)=gmass(id)/dnorm
            
          enddo
          
        endif
        
c     find a group of each type on this node to 
c     find principal axis system of the group type
        
        do id=1,mxungp
          
          jr=0
          ij=0
          safe=.false.
          
          do while(.not.safe.and.ij.lt.ngrp)
            
            ij=ij+1
            jr=jr+numgsit(lstgtp(ij))
            if(lstgtp(ij).eq.id)safe=.true.
            
          enddo
          
          if(safe)then

c     rotational inertia accumulator
            
            do k=1,3
              
              do kk=1,3
                
                rotinr(k,kk)=0.d0
                
              enddo
              
            enddo
            
            jr=jr-numgsit(id)
            do j=1,numgsit(id)
              
              jr=jr+1
              i=lstgot(jr)
              
              xxt(jr)=xxx(i)-gcmx(ij)
              yyt(jr)=yyy(i)-gcmy(ij)
              zzt(jr)=zzz(i)-gcmz(ij)
              
              call images(imcon,0,1,1,cell,xxt(jr),yyt(jr),zzt(jr))
              
              rotinr(1,1)=rotinr(1,1)+weight(i)*(xxt(jr)**2)
              rotinr(1,2)=rotinr(1,2)+weight(i)*xxt(jr)*yyt(jr)
              rotinr(1,3)=rotinr(1,3)+weight(i)*xxt(jr)*zzt(jr)
              rotinr(2,2)=rotinr(2,2)+weight(i)*(yyt(jr)**2)
              rotinr(2,3)=rotinr(2,3)+weight(i)*yyt(jr)*zzt(jr)
              rotinr(3,3)=rotinr(3,3)+weight(i)*(zzt(jr)**2)
              
            enddo
            
            rotinr(2,1)=rotinr(1,2)
            rotinr(3,1)=rotinr(1,3)
            rotinr(3,2)=rotinr(2,3)
            
            call jacobi(rotinr,rot1,3)
            
            rot(1)=rot1(1,1)
            rot(4)=rot1(2,1)
            rot(7)=rot1(3,1)
            rot(2)=rot1(1,2)
            rot(5)=rot1(2,2)
            rot(8)=rot1(3,2)
            rot(3)=rot1(1,3)
            rot(6)=rot1(2,3)
            rot(9)=rot1(3,3)
            
c     rotational inertia accumulators
            
            rotinx(id,1)=0.d0
            rotiny(id,1)=0.d0
            rotinz(id,1)=0.d0
            
            jr=jr-numgsit(id)
            do j=1,numgsit(id)
              
              jr=jr+1
              i=lstgot(jr)
              
c     site positions in principal axis system
              
              gxx(id,j)=rot(1)*xxt(jr)+rot(4)*yyt(jr)+rot(7)*zzt(jr)
              gyy(id,j)=rot(2)*xxt(jr)+rot(5)*yyt(jr)+rot(8)*zzt(jr)
              gzz(id,j)=rot(3)*xxt(jr)+rot(6)*yyt(jr)+rot(9)*zzt(jr)
              
c     impose rounding 
              
              if(abs(gxx(id,j)).lt.1.d-8)gxx(id,j)=0.d0
              if(abs(gyy(id,j)).lt.1.d-8)gyy(id,j)=0.d0
              if(abs(gzz(id,j)).lt.1.d-8)gzz(id,j)=0.d0
              
c     rotational inertia tensor of group type
              
              rotinx(id,1)=rotinx(id,1)+
     x          weight(i)*(gyy(id,j)**2+gzz(id,j)**2)
              rotiny(id,1)=rotiny(id,1)+
     x          weight(i)*(gzz(id,j)**2+gxx(id,j)**2)
              rotinz(id,1)=rotinz(id,1)+
     x          weight(i)*(gxx(id,j)**2+gyy(id,j)**2)
              
            enddo
            
c     set axis system such that: Ixx >=Iyy >=Izz
            
            rotxyz=max(rotinx(id,1),rotiny(id,1),rotinz(id,1))
            
            if(rotxyz.ge.rotinx(id,1))then
              
              if(rotiny(id,1).ge.rotxyz)then
                
                do j=1,numgsit(id)
                  
                  a1=gxx(id,j)
                  gxx(id,j)=gyy(id,j)
                  gyy(id,j)=-a1
                  
                enddo
                
                rotiny(id,1)=rotinx(id,1)
                rotinx(id,1)=rotxyz
                
              elseif(rotinz(id,1).ge.rotxyz)then
                
                do j=1,numgsit(id)
                  
                  a1=gxx(id,j)
                  gxx(id,j)=gzz(id,j)
                  gzz(id,j)=-a1
                  
                enddo
                
                rotinz(id,1)=rotinx(id,1)
                rotinx(id,1)=rotxyz
                
              endif
              
            endif
            
            if(rotinz(id,1).gt.rotiny(id,1))then
              
              do j=1,numgsit(id)
                
                a1=gyy(id,j)
                gyy(id,j)=gzz(id,j)
                gzz(id,j)=-a1
                
              enddo
              
              a1=rotinz(id,1)
              rotinz(id,1)=rotiny(id,1)
              rotiny(id,1)=a1
              
            endif
            
c     set up principal axis system in terms of site positions
            
c     test for (near) linear unit
            
            ill=0
            rtall=(rotinx(id,1)+rotiny(id,1)+rotinz(id,1))
            
            if(rtall.gt.1.d-5)then
              rotall=rtall
            else
              rotall=1.d0
            endif
            
            rotmin(id)=min(rotinx(id,1),rotiny(id,1))
            rotmin(id)=min(rotmin(id),rotinz(id,1))/rotall
            
            if((rotinx(id,1)/rotall).lt.1.d-5)ill=ill+1
            if((rotiny(id,1)/rotall).lt.1.d-5)ill=ill+1
            if((rotinz(id,1)/rotall).lt.1.d-5)ill=ill+1
            
            if(ill.ge.2)then

c     point particle only
              
              ind(id,1)=1
              ind(id,2)=1
              ind(id,3)=1
              
              do jj=1,9
                gaxs(id,jj)=0.d0
              enddo
              
            elseif(ill.eq.1)then
              
c     linear molecule
              
              ind(id,1)=1
              ind(id,2)=2
              ind(id,3)=1
              
              aa(1)=gxx(id,1)-gxx(id,2)
              aa(4)=gyy(id,1)-gyy(id,2)
              aa(7)=gzz(id,1)-gzz(id,2)
              rsq=sqrt(aa(1)**2+aa(4)**2+aa(7)**2)
              
              if(abs(aa(7)/rsq).gt.0.5d0)then
                
                rsq=sqrt(aa(4)**2+aa(7)**2)
                aa(2)=0.d0
                aa(5)=aa(7)/rsq
                aa(8)=-aa(4)/rsq
                
              elseif(abs(aa(4)/rsq).gt.0.5d0)then
                
                rsq=sqrt(aa(4)**2+aa(1)**2)
                aa(2)=-aa(4)/rsq
                aa(5)=aa(1)/rsq
                aa(8)=0.d0
                
              elseif(abs(aa(1)/rsq).gt.0.5d0)then
                
                rsq=sqrt(aa(1)**2+aa(7)**2)
                aa(2)=-aa(7)/rsq
                aa(5)=0.d0
                aa(8)=aa(1)/rsq
                
              endif
              
              aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
              aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
              aa(9)=aa(1)*aa(5)-aa(4)*aa(2)
              
              call invert(aa,bb,det)
              
              if(abs(det).lt.1.d-5)call error(idnode,306)
              
              do j=1,9
                gaxs(id,j)=bb(j)
              enddo
              
            elseif(ill.eq.0)then
              
c     non-linear molecule
              
              i1=1
              i2=1
              i3=1
              pass1=.true.
              dettest=1.d-1
              
              do while(pass1.and.i2.lt.numgsit(id)-1)
                
                i2=i2+1
                i3=i2
                pass2=.true.
                
                do while(pass2.and.i3.lt.numgsit(id))
                  
                  i3=i3+1
                  
                  aa(1)=gxx(id,i1)-gxx(id,i2)
                  aa(4)=gyy(id,i1)-gyy(id,i2)
                  aa(7)=gzz(id,i1)-gzz(id,i2)
                  aa(2)=gxx(id,i1)-gxx(id,i3)
                  aa(5)=gyy(id,i1)-gyy(id,i3)
                  aa(8)=gzz(id,i1)-gzz(id,i3)
                  aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
                  aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
                  aa(9)=aa(1)*aa(5)-aa(4)*aa(2)
                  
c     invert matrix
                  
                  call invert(aa,bb,det)
                  
c     check on size of determinant - to see if the 3 sites are
c     too close to being linear for safety.
                  
                  pass2=abs(det).lt.dettest
                  
                enddo
                
                pass1=abs(det).lt.dettest
                
              enddo
              
              if(abs(det).lt.dettest)call error(idnode,306)
              
c     store indices used
              
              ind(id,1)=i1
              ind(id,2)=i2
              ind(id,3)=i3
              
c     store coefficients 
              
              do j=1,9
                
                gaxs(id,j)=bb(j)
                
              enddo
              
            endif
            
          endif
          
        enddo
        
c     check that rigid unit does not contain frozen atoms
        
        safe=.true.
        
        jr=0
        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            if(lstfrz(i).ne.0)safe=.false.
            
          enddo
          
        enddo
        
c     global check on error condition
        
        if(mxnode.gt.1)call gstate(safe)
        if(.not.safe)call error(idnode,360)
        
c     quaternions for all rigid groups in system
        
        jr=0
        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          i1=lstrgd(jr+ind(id,1))
          i2=lstrgd(jr+ind(id,2))
          i3=lstrgd(jr+ind(id,3))
          
          jr=jr+numgsit(id)
          
c     group basis vectors
          
          aa(1)=xxx(i1)-xxx(i2)
          aa(4)=yyy(i1)-yyy(i2)
          aa(7)=zzz(i1)-zzz(i2)
          
          call images(imcon,0,1,1,cell,aa(1),aa(4),aa(7))
          
          if(rotmin(id).gt.1.d-5)then
            
            aa(2)=xxx(i1)-xxx(i3)
            aa(5)=yyy(i1)-yyy(i3)
            aa(8)=zzz(i1)-zzz(i3)
            
          else
            
            rsq=sqrt(aa(1)**2+aa(4)**2+aa(7)**2)
            
            if(abs(aa(7)/rsq).gt.0.5d0)then
              
              rsq=sqrt(aa(4)**2+aa(7)**2)
              aa(2)=0.d0
              aa(5)=aa(7)/rsq
              aa(8)=-aa(4)/rsq
              
            elseif(abs(aa(4)/rsq).gt.0.5d0)then
              
              rsq=sqrt(aa(4)**2+aa(1)**2)
              aa(2)=-aa(4)/rsq
              aa(5)=aa(1)/rsq
              aa(8)=0.d0
              
            elseif(abs(aa(1)/rsq).gt.0.5d0)then
              
              rsq=sqrt(aa(1)**2+aa(7)**2)
              aa(2)=-aa(7)/rsq
              aa(5)=0.d0
              aa(8)=aa(1)/rsq
              
            endif
            
          endif
          
          call images(imcon,0,1,1,cell,aa(2),aa(5),aa(8))
          
          aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
          aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
          aa(9)=aa(1)*aa(5)-aa(4)*aa(2)
          
c     group rotational matrix
          
          rot(1)=gaxs(id,1)*aa(1)+gaxs(id,4)*aa(2)+gaxs(id,7)*aa(3)
          rot(2)=gaxs(id,2)*aa(1)+gaxs(id,5)*aa(2)+gaxs(id,8)*aa(3)
          rot(3)=gaxs(id,3)*aa(1)+gaxs(id,6)*aa(2)+gaxs(id,9)*aa(3)
          rot(4)=gaxs(id,1)*aa(4)+gaxs(id,4)*aa(5)+gaxs(id,7)*aa(6)
          rot(5)=gaxs(id,2)*aa(4)+gaxs(id,5)*aa(5)+gaxs(id,8)*aa(6)
          rot(6)=gaxs(id,3)*aa(4)+gaxs(id,6)*aa(5)+gaxs(id,9)*aa(6)
          rot(7)=gaxs(id,1)*aa(7)+gaxs(id,4)*aa(8)+gaxs(id,7)*aa(9)
          rot(8)=gaxs(id,2)*aa(7)+gaxs(id,5)*aa(8)+gaxs(id,8)*aa(9)
          rot(9)=gaxs(id,3)*aa(7)+gaxs(id,6)*aa(8)+gaxs(id,9)*aa(9)
          
c     determine quaternions from rotational matrix
          
          aq=rot(1)+rot(5)
          bq=rot(2)-rot(4)
          cq=rot(6)-rot(8)
          dq=rot(2)+rot(4)
          eq=rot(3)+rot(7)
          fq=rot(6)+rot(8)
          gq=rot(3)-rot(7)
          hq=rot(1)-rot(5)
          
          q0(ig)=0.5d0*sqrt(aq+sqrt(aq*aq+bq*bq))
          
          if(q0(ig).gt.1.d-4)then
            
            q1(ig)=-0.25d0*cq/q0(ig)
            q2(ig)=0.25d0*gq/q0(ig)
            q3(ig)=-0.25d0*bq/q0(ig)
            
          else
            
            q1(ig)=0.5d0*sqrt(hq+sqrt(hq*hq+dq*dq))
            
            if(q1(ig).gt.1.d-4)then
              
              q2(ig)=0.25d0*dq/q1(ig)
              q3(ig)=0.25d0*eq/q1(ig)
              
            else
              
              q2(ig)=0.5d0*sqrt(-hq+sqrt(hq*hq+dq*dq))
              
              if(q2(ig).gt.1.d-4)then
                
                q3(ig)=0.25d0*fq/q2(ig)
                
              else
                
                q3(ig)=1.d0
                
              endif
              
            endif
            
          endif
          
c     normalise quaternions
          
          rnorm=1.d0/sqrt(q0(ig)**2+q1(ig)**2+q2(ig)**2+q3(ig)**2)
          q0(ig)=rnorm*q0(ig)
          q1(ig)=rnorm*q1(ig)
          q2(ig)=rnorm*q2(ig)
          q3(ig)=rnorm*q3(ig)
          
        enddo
        
c     test for redundant degrees of freedom
c     and ensure rotational inertias are non-zero
        
        degrot=0.d0
        
        if(lsolva)then
          degrot_sol(:)=0.d0
        endif
        
        do ig=1,ngrp
          
          id=lstgtp(ig)
          rotall=1.d0/max(1.d-5,rotinx(id,1)+rotiny(id,1)+
     x      rotinz(id,1))
          
          if(rotall*rotinx(id,1).lt.1.d-5)then
            degrot=degrot-1.d0
          endif
          
          if(rotall*rotiny(id,1).lt.1.d-5)then
            degrot=degrot-1.d0
          endif
          
          if(rotall*rotinz(id,1).lt.1d-5)then
            degrot=degrot-1.d0
          endif
          
        enddo
        
c     rotational degrees of freedom and rigid body contribution
c     to total degrees of freedom
        
        degrot=degrot+dble(ngrp)*3.d0
        degfre=degrot+dble(ngrp)*3.d0
        
        if(lsolva)then
          
          fngrp=1
          lngrp=0
          
          do itmols=1,mxtmls
            
            lngrp=lngrp+nummols(itmols)*numgrp(itmols)
            
            do ig=fngrp,lngrp
              
              id=lstgtp(ig)
              rotall=1.d0/max(1.d-5,rotinx(id,1)+rotiny(id,1)+
     x          rotinz(id,1))
              
              if(rotall*rotinx(id,1).lt.1.d-5)then
                degrot_sol(itmols)=degrot_sol(itmols)-1.d0
              endif
              
              if(rotall*rotiny(id,1).lt.1.d-5)then
                degrot_sol(itmols)=degrot_sol(itmols)-1.d0
              endif
              
              if(rotall*rotinz(id,1).lt.1d-5)then
                degrot_sol(itmols)=degrot_sol(itmols)-1.d0
              endif
              
            enddo
            
            fngrp=lngrp+1
            
          enddo
          
        endif
        
c     summarise results
        
        if(idnode.eq.0)then
          
          if(gmass(1).gt.0.d0)then
            
            write(nrite,'(/,/,12x,a)')' summary of rigid body set up'
            
            do id=1,mxungp
              
              if(gmass(id).gt.0.d0)then
                
                write(nrite,'(/,a,i10)')' group of type ',id
                write(nrite,'(12x,a,f20.10)')' total mass    ',
     x            gmass(id)
                write(nrite,'(12x,a,3f20.10)')' rot. inertia  ',
     x            rotinx(id,1),rotiny(id,1),rotinz(id,1)
                write(nrite,'(/,12x,a,3(8x,a7))')' site','a coord',
     x            'b coord','c coord'
                do j=1,numgsit(id)
                  write(nrite,'(12x,i5,1p,3e15.5)')j,gxx(id,j),
     x              gyy(id,j),gzz(id,j)
                enddo
                
              endif
              
            enddo
            
          endif
          
        endif
        
c     find number of unique groups 
        
        ngp=0
        do ig=1,ngrp
          
          ngp=max(ngp,lstgtp(ig))
          
        enddo
        
c     calculate reciprocal of rotational inertias 
        
        do id=1,ngp
          
          rotlim=max(1.d-2,rotinx(id,1)+rotiny(id,1)+
     x      rotinz(id,1))*1.d-5
          
          if(rotinx(id,1).lt.rotlim)then
            rotinx(id,2)=0.d0
          else
            rotinx(id,2)=1.d0/rotinx(id,1)
          endif
          
          if(rotiny(id,1).lt.rotlim)then
            rotiny(id,2)=0.d0
          else
            rotiny(id,2)=1.d0/rotiny(id,1)
          endif
          
          if(rotinz(id,1).lt.rotlim)then
            rotinz(id,2)=0.d0
          else
            rotinz(id,2)=1.d0/rotinz(id,1)
          endif
          
        enddo
        
c     Check of quaternion set up with atomic positions
        
        jr=0
        do ig=igrp1,igrp2
          
c     group type
          
          id=lstgtp(ig)
          
c     new rotational matrix
          
          rot(1)=q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
          rot(2)=2.d0*(q1(ig)*q2(ig)-q0(ig)*q3(ig))
          rot(3)=2.d0*(q1(ig)*q3(ig)+q0(ig)*q2(ig))
          rot(4)=2.d0*(q1(ig)*q2(ig)+q0(ig)*q3(ig))
          rot(5)=q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
          rot(6)=2.d0*(q2(ig)*q3(ig)-q0(ig)*q1(ig))
          rot(7)=2.d0*(q1(ig)*q3(ig)-q0(ig)*q2(ig))
          rot(8)=2.d0*(q2(ig)*q3(ig)+q0(ig)*q1(ig))
          rot(9)=q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
          
          do j=1,numgsit(id)
            
            jr=jr+1
            i=lstrgd(jr)
            
            xxt(i)=rot(1)*gxx(id,j)+rot(2)*gyy(id,j)+
     x        rot(3)*gzz(id,j)+gcmx(ig)
            yyt(i)=rot(4)*gxx(id,j)+rot(5)*gyy(id,j)+
     x        rot(6)*gzz(id,j)+gcmy(ig)
            zzt(i)=rot(7)*gxx(id,j)+rot(8)*gyy(id,j)+
     x        rot(9)*gzz(id,j)+gcmz(ig)
            
            
            txx(jr)=xxx(i)-xxt(i)
            tyy(jr)=yyy(i)-yyt(i)
            tzz(jr)=zzz(i)-zzt(i)
            
          enddo
          
        enddo
        
        call images(imcon,0,1,jr,cell,txx,tyy,tzz)
        
c     set tolerance for testing quaternion setup.
        
        rsq=0.d0
        tol=1.d-2
        
        do i=1,jr
          
          rrr=txx(i)**2+tyy(i)**2+tzz(i)**2
          if(rrr.gt.tol)then 
            
            rsq=rrr
            
          endif
          
        enddo
        
c     exit if error in set up
        
        safe=.true.
        if(rsq.gt.tol)safe=.false.
        if(mxnode.gt.1)call gstate(safe)
        
        if(.not.safe)call  error(idnode,310)
        
c     sort lstgot into ascending order
        
        call shellsort(jt,lstgot)
        
c     check that no site is in more than 1 rigid group
        
        i=1
        safe=.true.
        do while(i.lt.jt)
          
          i=i+1
          linear=.true.
          do while(linear)
            
            linear=.false.
            
            if(lstgot(i).eq.lstgot(i-1))then
              
              linear=.true.
              safe=.false.
              jt=jt-1
              
              do j=i,jt
                lstgot(j)=lstgot(j+1)
              enddo
              
            endif
            
            if(i.ge.jt)linear=.false.
            
          enddo
          
        enddo
        
        if(.not.safe)call error(idnode,320)
        
c     list of 'free' sites
        
        ii=1
        jj=0
        do i=1,natms
          
          if(lstgot(ii).eq.i)then
            
            ii=ii+1
            
          else
            
            if(lstfrz(i).eq.0)then
              jj=jj+1
              lstfre(jj)=i
            endif
            
          endif
          
        enddo
        
c     number of free sites
        
        ntfree=jj
        
c     list of atoms integrated on this node
        
        jr=0
        do ig=igrp1,igrp2
          
          id=lstgtp(ig)
          jr=jr+numgsit(id)
          
        enddo
        
        do i=1,jr
          
          lstme(i)=lstrgd(i)
          
        enddo
        
c     block parameters for free atoms
        
        ifre1=(idnode*ntfree)/mxnode+1
        ifre2=((idnode+1)*ntfree)/mxnode
        
        do i=ifre1,ifre2
          
          jr=jr+1
          lstme(jr)=lstfre(i)
          
        enddo
        
c     exchange quaternion data with other nodes
        
        if(mxnode.gt.1)call merge4
     x    (idnode,mxnode,ngrp,mxbuff,q0,q1,q2,q3,buffer)
        
      endif
      
      if(lsolva)lstgot_sol(:)=lstgot(:)
      
c     deallocate work arrays
      
      deallocate (ind,lstgot,stat=fail(1))
      deallocate (txx,tyy,tzz,stat=fail(2))
      deallocate (xxt,yyt,zzt,stat=fail(3))
      deallocate (gaxs,rotmin,stat=fail(4))
      deallocate (accum,stat=fail(5))
      
      return
      end subroutine quatbook
      
      subroutine abort_field_read(kode,idnode,nfield)
      
c***********************************************************************
c     
c     dl_poly subroutine for aborting FIELD file read
c     
c     copyright - daresbury laboratory 
c     author    - w. smith    aug 2003
c     
c***********************************************************************
      
      implicit none
      
      integer kode,idnode,nfield
      
      if(idnode.eq.0)close (nfield)
      
      if(kode.eq.1)then
        
c     end of field file error exit
        
        call error(idnode,52)
        
      elseif(kode.eq.2)then
        
c     unrecognised directive in field file
        
        call error(idnode,4)
        
      endif
      
      return
      end subroutine abort_field_read
      
      subroutine abort_control_read(kode,idnode,nread)
      
c***********************************************************************
c     
c     dl_poly subroutine for aborting CONTROL file read
c     
c     copyright - daresbury laboratory 
c     author    - w. smith    aug 2003
c     
c***********************************************************************
      
      implicit none
      
      integer kode,idnode,nread
      
      if(idnode.eq.0)close (nread)
      
      if(kode.eq.1)then
        
c     end of control file error exit
        
        call error(idnode,53)
        
      elseif(kode.eq.2)then
        
c     general error exit from field file processing
        
        call error(idnode,0)
        
      endif
      
      return
      end subroutine abort_control_read
      
      subroutine abort_config_read(kode,idnode,nconf)
      
c***********************************************************************
c     
c     dl_poly subroutine for aborting CONTROL file read
c     
c     copyright - daresbury laboratory 
c     author    - w. smith    aug 2003
c     
c***********************************************************************
      
      implicit none
      
      integer kode,idnode,nconf
      
      if(idnode.eq.0)close (nconf)
      
      if(kode.eq.1)then
        
c     general error exit from field file processing
        
        call error(idnode,54)
        
      elseif(kode.eq.2)then
        
c     end of config file error exit
        
        call error(idnode,55)
        
      endif
      
      return
      end subroutine abort_config_read
      
      subroutine neutbook(lneut,idnode,natms,nneut)
      
c***********************************************************************
c     
c     dl_poly subroutine for neutral group bookkeeping
c     
c     copyright - daresbury laboratory
c     author    - w. smith    nov 2003
c     
c***********************************************************************
      
      implicit none
      
      logical lneut,safe
      integer idnode,natms,nneut,i
      
      safe=.true.
      
c     neutral group bookkeeping: sites must be listed consecutively
      
      if(lneut)then
        
        if(lstneu(1).ne.1)call error(idnode,230)
        
        neulst(1)=1
        nneut=1
        
        do i=2,natms
          
          safe=.false.
          if(lstneu(i).eq.lstneu(i-1))safe=.true.
          if(lstneu(i).eq.lstneu(i-1)+1)then
            
            safe=.true.
            nneut=nneut+1
            if(nneut.gt.mxneut)call error(idnode,220)
            neulst(nneut)=i
            
          endif
          
          if(.not.safe)call error(idnode,230)
          
        enddo
        
        neulst(nneut+1)=natms+1
        
      endif
      
      return
      
      end subroutine neutbook
      
      subroutine intlist
     x  (lshmov,lcnb,idnode,mxnode,natms,nscons,ntangl,ntbond,
     x  ntcons,ntdihd,ntinv,ntpmls,ntteth,ntshl,ntpmf,nspmf,ngrp)
      
c***********************************************************************
c     
c     dl_poly subroutine for constructing the interaction lists
c     for the entire simulated system
c     
c     parallel replicated dat version : block data
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith        july 1992
c     amended   - t.forester      oct 1993
c     amended   - t.forester      dec 1994 : block data
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lshmov,safe1,lcnb,lchk,lfail
      integer idnode,mxnode,natms,nscons,ntangl,ntbond,ntcons
      integer ntdihd,ntinv,ntpmls,ntteth,ntshl,ntpmf,nspmf
      integer ibonds,jbonds,kbonds,ipmf,jpmf,iangle,jangle,kangle
      integer idihed,jdihed,kdihed,iinver,jinver,kinver,iteths
      integer jteths,kteths,ishels,jshels,kshels,ntbon0,ntpmf0
      integer ntang0,ntdih0,ntinv0,nttet0,ntshl0,ntcon0,idum
      integer itmols,isite,iconst,jconst,kconst,ibnd1,ibnd2,ipmf1
      integer ipmf2,iang1,iang2,idih1,idih2,iinv1,iinv2,itet1
      integer itet2,ishl1,ishl2,imols,lbonds,lpmf,jj,nnn,langle
      integer ldihed,linver,lteths,lshels,i,ii,ntmp,klo,khi,ngrp
      integer klo0,ifail,iloop,nnode,nscons0,nscons1,icon,fail
      integer kcons,id,jdnode,lconst,itry,iatom,jatom,j,nfail
      real(8) tol
      
      integer, allocatable :: itest(:),index(:),kscons(:)
      integer, allocatable :: msite(:),mconst(:),listin(:)
      
      dimension fail(4)
      
      data fail/0,0,0,0/
      
c     allocate work arrays
      
      allocate (itest(mxtmls),index(mxtmls),stat=fail(1))
      allocate (msite(mxtmls),mconst(mxtmls),stat=fail(2))
      allocate (listin(mxatms),stat=fail(3))
      allocate (kscons(0:mxnode-1),stat=fail(4))
      do i=1,4
        if(fail(i).ne.0)call error(idnode,1800)
      enddo
      
c     initialise bookkeeping indices
      
      ibonds=0
      jbonds=0
      kbonds=0
      ipmf=0
      jpmf=0
      iangle=0
      jangle=0
      kangle=0
      idihed=0
      jdihed=0
      kdihed=0
      iinver=0
      jinver=0
      kinver=0
      iteths=0
      jteths=0
      kteths=0
      ishels=0
      jshels=0
      kshels=0
      safe=.true.
      safe1=.true.
      
c     find total number of bonds,pmf constraints,bond constraints,
c     angles,dihedrals,inversions, tethers,core-shells, in system 
c     - ignoring frozen atoms
      
      ntbon0=0
      ntpmf0=0
      ntcon0=0
      ntang0=0
      ntdih0=0
      ntinv0=0
      nttet0=0
      ntshl0=0
      nscons=0
      ntcons=0
      
      do itmols=1,ntpmls
        
        ntbon0=ntbon0+nummols(itmols)*numbonds(itmols)
        ntpmf0=ntpmf0+nummols(itmols)*numpmf(itmols)
        ntcon0=ntcon0+nummols(itmols)*numcon(itmols)
        ntang0=ntang0+nummols(itmols)*numang(itmols)
        ntdih0=ntdih0+nummols(itmols)*numdih(itmols)
        ntinv0=ntinv0+nummols(itmols)*numinv(itmols)
        nttet0=nttet0+nummols(itmols)*numteth(itmols)
        ntshl0=ntshl0+nummols(itmols)*numshl(itmols)
        
      enddo
      
      isite=0
      iconst=0
      jconst=0
      kconst=0
      
c     first and last index of bonds, angles etc for this node
      
      ibnd1=(idnode*ntbon0)/mxnode+1
      ibnd2=((idnode+1)*ntbon0)/mxnode
      
      ipmf1=(idnode*ntpmf0)/mxnode+1
      ipmf2=((idnode+1)*ntpmf0)/mxnode
      ntpmf=ntpmf0
      nspmf=ipmf2+1-ipmf1
      
      iang1=(idnode*ntang0)/mxnode+1
      iang2=((idnode+1)*ntang0)/mxnode
      
      idih1=(idnode*ntdih0)/mxnode+1
      idih2=((idnode+1)*ntdih0)/mxnode
      
      iinv1=(idnode*ntinv0)/mxnode+1
      iinv2=((idnode+1)*ntinv0)/mxnode
      
      itet1=(idnode*nttet0)/mxnode+1
      itet2=((idnode+1)*nttet0)/mxnode
      
      ishl1=(idnode*ntshl0)/mxnode+1
      ishl2=((idnode+1)*ntshl0)/mxnode
      
c     loop over molecule types
      
      do itmols=1,ntpmls
        
        
c     loop over molecules in system
        
        do imols=1,nummols(itmols)
          
c     construct bond constraint list later
c     construct chemical bond interaction list
          
          do lbonds=1,numbonds(itmols)
            
            ibonds=ibonds+1
            
            if(ibonds.ge.ibnd1.and.ibonds.le.ibnd2)then
              
              jbonds=jbonds+1
              if(jbonds.le.mxbond)then
                
                listbnd(jbonds,1)=lbonds+kbonds
                listbnd(jbonds,2)=lstbnd(lbonds+kbonds,1)
     x            +isite
                listbnd(jbonds,3)=lstbnd(lbonds+kbonds,2)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if(mxnode.gt.1)call gstate(safe)
          if(.not.safe)call error(idnode,31)
          
c     construct pmf site lists - no exclusions
          
          do lpmf=1,numpmf(itmols)
            
            ipmf=ipmf+1
            
            if(ipmf.ge.ipmf1.and.ipmf.le.ipmf2)then
              
              jpmf=jpmf+1
              if(jpmf.le.mspmf)then
                
                nnn=npmf(1)+npmf(2)
                if(nnn.le.mxspmf)then
                  
                  do jj=1,npmf(1)+npmf(2)
                    lstpmf(jj,jpmf)=indpmf(jj)+isite
                  enddo
                  
                else
                  
                  safe=.false.
                  
                endif
                
              else
                
                safe1=.false.
                
              endif
              
            endif
            
          enddo
          
          if(mxnode.gt.1)call gstate(safe1)
          if(.not.safe1)call error(idnode,458)
          
          if(mxnode.gt.1)call gstate(safe)
          if(.not.safe)call error(idnode,460)
          
c     construct valence angle interaction list
          
          do langle=1,numang(itmols)
            
            iangle=iangle+1
            
            if(iangle.ge.iang1.and.iangle.le.iang2)then
              
              jangle=jangle+1
              if(jangle.le.mxangl)then
                
                listang(jangle,1)=langle+kangle
                listang(jangle,2)=lstang(langle+kangle,1)
     x            +isite
                listang(jangle,3)=lstang(langle+kangle,2)
     x            +isite
                listang(jangle,4)=lstang(langle+kangle,3)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if(mxnode.gt.1)call gstate(safe)
          if(.not.safe)call error(idnode,51)
          
c     construct dihedral angle interaction list
          
          do ldihed=1,numdih(itmols)
            
            idihed=idihed+1
            
            if(idihed.ge.idih1.and.idihed.le.idih2)then
              
              jdihed=jdihed+1
              if(jdihed.le.mxdihd)then
                
                listdih(jdihed,1)=ldihed+kdihed
                listdih(jdihed,2)=lstdih(ldihed+kdihed,1)
     x            +isite
                listdih(jdihed,3)=lstdih(ldihed+kdihed,2)
     x            +isite
                listdih(jdihed,4)=lstdih(ldihed+kdihed,3)
     x            +isite
                listdih(jdihed,5)=lstdih(ldihed+kdihed,4)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if(mxnode.gt.1)call gstate(safe)
          if(.not.safe)call error(idnode,61)
          
c     construct inversion potential list
          
          do linver=1,numinv(itmols)
            
            iinver=iinver+1
            
            if(iinver.ge.iinv1.and.iinver.le.iinv2)then
              
              jinver=jinver+1
              if(jinver.le.mxinv)then
                
                listinv(jinver,1)=linver+kinver
                listinv(jinver,2)=lstinv(linver+kinver,1)
     x            +isite
                listinv(jinver,3)=lstinv(linver+kinver,2)
     x            +isite
                listinv(jinver,4)=lstinv(linver+kinver,3)
     x            +isite
                listinv(jinver,5)=lstinv(linver+kinver,4)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if(mxnode.gt.1)call gstate(safe)
          if(.not.safe)call error(idnode,77)
          
c     construct tethered atoms interaction list
          
          do lteths=1,numteth(itmols)
            
            iteths=iteths+1
            
            if(iteths.ge.itet1.and.iteths.le.itet2)then
              
              jteths=jteths+1
              if(jteths.le.msteth)then
                
                listtet(jteths,1)=lteths+kteths
                listtet(jteths,2)=lsttet(lteths+kteths)+isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if(mxnode.gt.1)call gstate(safe)
          if(.not.safe)call error(idnode,63)
          
c     construct core-shell list
          
          do lshels=1,numshl(itmols)
            
            ishels=ishels+1
            
            if(ishels.ge.ishl1.and.ishels.le.ishl2)then
              
              jshels=jshels+1
              if(jshels.le.mxshl)then
                
                listshl(jshels,1)=lshels+kshels
                listshl(jshels,2)=lstshl(lshels+kshels,1)
     x            +isite
                listshl(jshels,3)=lstshl(lshels+kshels,2)
     x            +isite
                
              else
                
                safe=.false.
                
              endif
              
            endif
            
          enddo
          
          if(mxnode.gt.1)call gstate(safe)
          if(.not.safe)call error(idnode,59)
          
          isite=isite+numsit(itmols)
          
        enddo
        
        kbonds=kbonds+numbonds(itmols)
        kangle=kangle+numang(itmols)
        kdihed=kdihed+numdih(itmols)
        kinver=kinver+numinv(itmols)
        kteths=kteths+numteth(itmols)
        kshels=kshels+numshl(itmols)
        
      enddo
      
c     store array counters for bookkeeping
      
      ntbond=ibonds
      ntangl=iangle
      ntdihd=idihed
      ntinv=iinver
      ntteth=iteths
      ntshl=ishels
      
c     pass bond constraint information to other nodes
      
      if(ntcon0.gt.0)then
        
        ntcons=ntcon0
        
c     find starting site no. and constraint no. for each molec. type
        
        msite(1)=0
        mconst(1)=0
        
        do itmols=2,ntpmls
          
          msite(itmols)=msite(itmols-1)+numsit(itmols-1)*
     x      nummols(itmols-1)
          mconst(itmols)=mconst(itmols-1)+numcon(itmols-1)
          
        enddo
        
c     sort molecules into ascending order of number of constraints
        
        do i=1,ntpmls
          
          itest(i)=numcon(i)
          index(i)=0
          
        enddo
        
        call shellsort(ntpmls,itest)
        
        do i=1,ntpmls
          
          lchk=.true.
          do j=1,ntpmls
            
            if(itest(i).eq.numcon(j))then
              
              if(lchk)then 
                index(i)=j
                lchk=.false.
                
              endif
              
              do ii=1,i-1
                if(index(ii).eq.j)lchk=.true.
              enddo
              
            endif
            
          enddo
          
        enddo
        
c     load balance to within 10%
        
        tol=1.0d0+(0.10d0)/2.d0
        kcons=(ntcons)/mxnode
        ntmp=0
        
c     find smallest constrained molecule to allocate to a node
        
        do i=1,ntpmls
          
          if(ntmp.le.mxnode)then
            
            if(numcon(index(i)).gt.0)then
              ntmp=ntmp+nummols(index(i))
              klo=max(0,kcons-numcon(index(i))/2)
              khi=klo+numcon(index(i))+1
            endif
            
          endif
          
        enddo
        
c     reset hi/lo limits if molecules contain too many constraints
        
        if(dble(khi)/dble(max(1,klo)).gt.tol)then
          klo=nint(dble(kcons)/tol)
          khi=nint(dble(kcons)*tol)+1
        endif
        
c     store lo value for later
        
        klo0=klo
        
c     begin assignment of constraints ----------------------------------
        
        ifail=-1
        lfail=.true.
        do while(lfail)
          
          ifail=ifail+1
          
          if(ifail.gt.ntpmls)then
            call error(idnode,432)
          endif
          
          iconst=0
          jconst=0
          kconst=0
          lconst=0
          
c     zero running totals of constraints on each processor
          
          do id=0,mxnode-1
            kscons(id)=0
          enddo
          
          iloop=0
          lfail=.false.
          iconst=0
          jconst=0
          nnode=0
          
c     assign difficult molecules in blocks
          
          if(ifail.gt.0)then
            
            nfail=0
            do i=1,ifail
              
              ii=ntpmls+1-i
              nfail=nfail+nummols(index(ii))*numcon(index(ii))
              
            enddo
            
c     decide on number of processors to split over
            
            nnode=int(dble(nfail)/dble(max(kcons,1))+1.d0/tol)
            nnode=max(2,nnode)
            nnode=min(nnode,mxnode)
            
c     assign to processors 0..nnode-1
            
            do id=0,nnode-1
              
              nscons0=(id*nfail)/nnode+1
              nscons1=((id+1)*nfail)/nnode
              
              kscons(id)=nscons1+1-nscons0
              
            enddo
            
c     this processors block
            
            nscons0=(idnode*nfail)/nnode+1
            nscons1=((idnode+1)*nfail)/nnode
            
c     assign in blocks
            
            do itmols=ntpmls,ntpmls-ifail+1,-1
              
              ii=index(itmols)
              icon=numcon(ii)
              kconst=mconst(ii)
              
              do imols=1,nummols(ii)
                
                isite=msite(ii)+(imols-1)*numsit(ii)
                
c     construct bond constraint list
                
                do lconst=1,numcon(ii)
                  
                  iconst=iconst+1
                  
                  if(iconst.ge.nscons0.and.iconst.le.nscons1)then
                    
                    jconst=jconst+1
                    
                    if(jconst.le.mxcons)then
                      
                      listcon(jconst,1)=lconst+kconst
                      iatom=lstcon(lconst+kconst,1)+isite
                      jatom=lstcon(lconst+kconst,2)+isite
                      
                      listcon(jconst,2)=iatom
                      listcon(jconst,3)=jatom
                      
                    else
                      
                      safe=.false.
                      
                    endif
                    
                  endif
                  
                enddo
                
              enddo
              
            enddo
            
          endif
          
c     assign non-problematic molecules
          
          jdnode=mod(nnode+1,mxnode)
          
          do itmols=ntpmls-ifail,1,-1
            
            ii=index(itmols)
            icon=numcon(ii)
            kconst=mconst(ii)
            
            do imols=1,nummols(ii)
              
              itry=0
              lchk=.true.
              do while(lchk)
                
                if(kscons(jdnode)+icon.le.klo)then
                  
                  if(jdnode.ne.idnode)then
                    kscons(jdnode)=kscons(jdnode)+icon
                    jdnode=mod(jdnode+1,mxnode)
                    lchk=.false.
                  else
                    
c     construct bond constraint list
                    
                    isite=msite(ii)+(imols-1)*numsit(ii)
                    do lconst=1,numcon(ii)
                      
                      jconst=jconst+1
                      
                      if(jconst.le.mxcons)then
                        
                        listcon(jconst,1)=lconst+kconst
                        iatom=lstcon(lconst+kconst,1)+isite
                        jatom=lstcon(lconst+kconst,2)+isite
                        listcon(jconst,2)=iatom
                        listcon(jconst,3)=jatom
                        
                      else
                        
                        safe=.false.
                        
                      endif
                      
                    enddo
                    
                    kscons(jdnode)=kscons(jdnode)+icon
                    jdnode=mod(jdnode+1,mxnode)
                    lchk=.false.
                    
                  endif
                  
                else
                  
                  jdnode=mod(jdnode+1,mxnode)
                  lchk=.true.
                  itry=itry+1
                  
                endif
                
                if(lchk.and.itry.gt.mxnode)then
                  
                  klo=kcons
                  kcons=khi
                  itry=0
                  iloop=iloop+1
                  
                endif
                
c     split molecule across nodes if have to
                
                if(iloop.gt.3)then
                  lfail=.true.
                  kcons=ntcons/mxnode
                  klo=klo0
                  lchk=.false.
                endif
                
              enddo
              
            enddo
            
          enddo
          
c     check no node has less than minimum number
          
          do id=0,mxnode-1
            if(kscons(id).lt.klo0)then 
              lfail=.true.
            endif
          enddo
          
        enddo
        
        if(mxnode.gt.1)call gstate(safe)
        if(.not.safe)then
          
          if(mxnode.gt.1)call gimax(jconst,1,idum)
          if(idnode.eq.0)write(nrite,'(a,i10,a,i10)')
     x      'Number of constraints found ',jconst,'Max allowed ',mxcons
          
          call error(idnode,41)
          
        endif
        
        nscons=kscons(idnode)
        
        call passcon
     x    (lshmov,idnode,mxnode,natms,nscons,lashap,lishap,listme,
     x    listin,listot,listcon,lstfrz)
        
      endif
      
      if(npmf(1).gt.0)then
        
        call passpmf
     x    (idnode,mxnode,natms,nspmf,listpm,listin,lstpmt,lstpmf,npmf)
        
      endif
      
c     pass rigid body data
      
      lcnb=.false.
      if(ntcons.gt.0.and.ngrp.gt.0)then
        
        call passquat
     x    (lcnb,idnode,mxnode,natms,ngrp,nscons,ntpmls,listin,
     x    listcon,lstrgd,lstout,lstcsit,lstgtp,nummols,numgrp,
     x    numgsit)
        
      endif
      
c     deallocate work arrays
      
      deallocate(itest,index,msite,stat=fail(1))
      deallocate(mconst,kscons,listin,stat=fail(2))
      
      return
      end subroutine intlist
      
      subroutine ensemble_selection
     x  (directive,lens,kill,inhc,idnode,keyens,mode,nchain,nrespa,taut,
     x  taup)
c***********************************************************************
c     
c     dl_poly subroutine for selecting the ensemble and reading 
c     the required parameters
c     copyright - daresbury laboratory
c     author    - w. smith    feb 2008
c     
c     NHC is added to NVT and NPT ensembles
c     copyright - M.R.Momeni and F.A.Shakib
c     authors   - M.R.Momeni and F.A.Shakib 2021
c
c***********************************************************************
      
      implicit none
      
      character*1 directive(lenrec)
      logical kill,lens,inhc
      integer keyens,idnode,idum,mode
      integer nrespa,nchain
      real(8) taut,taup
      
      if(findstring('nve',directive,idum))then

        keyens=0
        if(idnode.eq.0)write(nrite,
     x    "(/,1x,'microcanonical ensemble')")
        if(lens)then
          call error(idnode,-414)
          kill=.true.
        endif
        lens=.true.
        
      elseif(findstring('nvt',directive,idum))then
        
        if(findstring('evans',directive,idum))then
          
          keyens=1
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Evans Gaussian temperature constraints',
     x      ' in use')")
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          
        elseif(findstring('ber',directive,idum))then
          
          keyens=2
          taut=dblstr(directive,69,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Berendsen thermostat',
     x      /,1x,'thermostat relaxation time     ',1p,e12.4)")
     x      taut
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          
        elseif(findstring('hoover',directive,idum))then
          
          keyens=3
          taut=dblstr(directive,69,idum)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Nose-Hoover ',
     x      /,1x,'thermostat relaxation time     ',1p,e12.4)")
     x      taut
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

        elseif(findstring('nhc',directive,idum))then
         

          keyens=9
            taut=dblstr(directive,lenrec,idum)
            nrespa=intstr(directive,lenrec,idum)
            nchain=intstr(directive,lenrec,idum)
            if((taut.eq.0).or.(nrespa.eq.0).or.(nchain.eq.0))
     x        call error(idnode,4031)
            nrespa=max(nrespa,1)
            nchain=max(nchain,1)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Nose-Hoover Chain',
     x      /,1x,'thermostat relaxation time',1p,e12.4,
     x      /,1x,'number of RESPA steps  ',1p,i6,
     x      /,1x,'number of chains       ',1p,i6)")
     x      taut,nrespa,nchain
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          inhc=.true.         
        else

          kill=.true.
          if(idnode.eq.0)write(nrite,"(/,/,100a1)")record
          call error(idnode,-3)
          
        endif

c *******************************************************************      
        
      elseif(findstring('npt',directive,idum))then
        
        if(findstring('ber',directive,idum))then
          
          keyens=4
          taut=dblstr(directive,lenrec,idum)
          taup=dblstr(directive,lenrec,idum)
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Berendsen isotropic N-P-T',
     x      /,1x,'thermostat relaxation time     ',1p,e12.4,
     x      /,1x,'barostat relaxation time       ',1p,e12.4)")
     x      taut,taup
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          
        elseif(findstring('hoover',directive,idum))then
          
          keyens=5
          taut=dblstr(directive,lenrec,idum)
          taup=dblstr(directive,lenrec,idum)
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Nose-Hoover  (Melchionna) isotropic N-P-T ',
     x      /,1x,'thermostat relaxation time     ',1p,e12.4,
     x      /,1x,'barostat relaxation time       ',1p,e12.4)")
     x      taut,taup
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          
c *******************************************************************      
c     M.R.Momeni & F.A.Shakib
c     Method Development and Materials Simulation Laboratory

        elseif(findstring('nhc',directive,idum))then
         
          keyens=10
            taut=dblstr(directive,lenrec,idum)
            taup=dblstr(directive,lenrec,idum)
            nrespa=intstr(directive,lenrec,idum)
            nchain=intstr(directive,lenrec,idum)
            if((taut.eq.0).or.(taup.eq.0).or.(nrespa.eq.0).or.
     x         (nchain.eq.0))call error(idnode,4032)
            nrespa=max(nrespa,1)
            nchain=max(nchain,1)
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Nose-Hoover Chain',
     x      /,1x,'thermostat relaxation time',1p,e12.4,
     x      /,1x,'barostat relaxation time  ',1p,e12.4,
     x      /,1x,'number of RESPA steps ',1p,i6,
     x      /,1x,'number of chains      ',1p,i6)")
     x      taut,taup,nrespa,nchain
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          inhc=.true.         

c *******************************************************************      
        
      else
          
          kill=.true.
          if(idnode.eq.0)write(nrite,"(/,/,100a1)")record
          call error(idnode,-3)
          
        endif

      elseif(findstring('nst',directive,idum))then
        
        mode=0
        if(findstring('block',directive,idum))mode=1
        if(findstring('surf',directive,idum))mode=2
        if(findstring('slab',directive,idum))mode=3
        
        if(findstring('ber',directive,idum))then
          
          keyens=6
          taut=dblstr(directive,lenrec,idum)
          taup=dblstr(directive,lenrec,idum)
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Berendsen anisotropic N-P-T',
     x      /,1x,'thermostat relaxation time     ',1p,e12.4,
     x      /,1x,'barostat relaxation time       ',1p,e12.4)")
     x      taut,taup
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          
        elseif(findstring('hoover',directive,idum))then
          
          keyens=7
          taut=dblstr(directive,lenrec,idum)
          taup=dblstr(directive,lenrec,idum)
          
          if(idnode.eq.0)write(nrite,
     x      "(/,1x,'Nose-Hoover (Melchionna) anisotropic N-P-T ',
     x      /,1x,'thermostat relaxation time     ',1p,e12.4,
     x      /,1x,'barostat relaxation time       ',1p,e12.4)")
     x      taut,taup
          if(lens)then
            call error(idnode,-414)
            kill=.true.
          endif
          lens=.true.
          
        else
          
          kill=.true.
          if(idnode.eq.0)write(nrite,"(/,/,100a1)")record
          call error(idnode,-3)
          
        endif
        
        if(idnode.eq.0)then
          
          if(mode.eq.0)then
            write(nrite,"(/,1x,'NST mode 0 X<>Y<>Z')")
          elseif(mode.eq.1)then
            write(nrite,
     x        "(/,1x,'NST mode 1 X<>Y<>Z (rectangular block)')")
          elseif(mode.eq.2)then
            write(nrite,
     x        "(/,1x,'NST mode 2 X=Y<>Z (liquid surface)')")
          elseif(mode.eq.3)then
            write(nrite,
     x        "(/,1x,'NST mode 3 X<>Y<>Z (solid slab)')")
          endif
          
        endif
        
      elseif(findstring('pmf',directive,idum))then
        
        keyens=8
        if(idnode.eq.0)write(nrite,
     x    "(/,1x,'potential of mean force calculation (NVE)')")
        if(lens)then
          call error(idnode,-414)
          kill=.true.
        endif
        lens=.true.
        
      else
        
        call error(idnode,-436)
        kill=.true.
        
      endif
      
      return
      end subroutine ensemble_selection
      
      subroutine neb_option
     x  (directive,lneb,lminopt,idnode,numneb,keytol,sprneb,
     x  opttol,hyp_units)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading parameters for NEB option
c     copyright - daresbury laboratory
c     author    - w. smith    feb 2008
c     
c***********************************************************************
      
      implicit none
      
      character*8 cunit
      character*1 directive(lenrec)
      logical lneb,lminopt,endneb,safe
      integer numneb,idnode,keytol,i,idum
      real(8) sprneb,opttol,hyp_units
      
      if(lminopt)call error(idnode,225)
      lminopt=.true.
      lneb=.true.
      endneb=.false.
      numneb=intstr(directive,lenrec,idum)
      if(numneb.eq.0)numneb=1
      numneb=min(maxneb,numneb)
      
      hyp_units=1.d0
      do while(.not.endneb)
        
        call getrec(safe,idnode,nread)
        if(.not.safe)call abort_control_read(1,idnode,nread)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        call copystring(record,directive,lenrec)
        
        if(record(1).eq.'#'.or.record(1).eq.'&')then
c     information only - skip record
          cycle
        elseif(findstring('endneb',directive,idum))then
          endneb=.true.
        elseif(findstring('units',directive,idum))then
          hyp_units=energy_unit()
          call getword(cunit,directive,8,lenrec)
          call getword(cunit,directive,8,lenrec)
        elseif(findstring('basin_1',directive,idum))then
          call striptext(directive,lenrec,1)
          do i=1,numneb
            bsn_1(i)=intstr(directive,lenrec,idum)
          enddo
        elseif(findstring('basin_2',directive,idum))then
          call striptext(directive,lenrec,1)
          do i=1,numneb
            bsn_2(i)=intstr(directive,lenrec,idum)
          enddo
        elseif(findstring('neb_spring',directive,idum))then
          sprneb=dblstr(directive,lenrec,idum)
        elseif(findstring('forc',directive,idum))then
          keytol=0
          opttol=dblstr(directive,lenrec,idum)
        elseif(findstring('ener',directive,idum))then
          keytol=1
          opttol=dblstr(directive,lenrec,idum)
        elseif(findstring('posi',directive,idum))then
          keytol=2
          opttol=dblstr(directive,lenrec,idum)
        endif
        
      enddo
      
      if(idnode.eq.0)then
        
        write(nrite,"(/,1x,'NEB calculation controls')")
        write(nrite,"(/,1x,'identity of basin 1            ',
     x    10i10)")(bsn_1(i),i=1,numneb)
        write(nrite,"(1x,'identity of basin 2            ',
     x    10i10)")(bsn_2(i),i=1,numneb)
        write(nrite,
     x    "(1x,'NEB spring constant            ',e12.4,
     x    /,1x,'minimisation tolerance         ',e12.4,
     x    /,1x,'energy units                   ',2x,a8)")
     x    sprneb,opttol,cunit
        
        call print_optim(keytol)
        
      endif
      
c     units conversion
      
      sprneb=sprneb*hyp_units
      if(keytol.lt.2)opttol=opttol*hyp_units
      
      return
      end subroutine  neb_option
      
      subroutine bpd_option
     x  (directive,seek,lbpd,ltad,lminopt,prechk,nebgo,keybpd,idnode,
     x  nblock,ntrack,keytol,ebias,vmin,catchrad,sprneb,opttol,
     x  hyp_units)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading parameters for bias potential
c     dynamics option
c     copyright - daresbury laboratory
c     author    - w. smith    feb 2008
c     
c***********************************************************************
      
      implicit none
      
      character*8 cunit,seek
      character*1 directive(lenrec)
      logical lbpd,ltad,lminopt,prechk,endbpd,safe,nebgo
      integer keybpd,idnode,nblock,ntrack,keytol,idum
      real(8) ebias,vmin,catchrad,sprneb,opttol,hyp_units
      
      if(lminopt)call error(idnode,225)
      if(ltad)call error(idnode,2355)
      lminopt=.true.
      lbpd=.true.
      endbpd=.false.
      cunit=" dl_poly"
      if(idnode.eq.0)
     x  write(nrite,"(/,1x,'bias potential dynamics controls')")
      
      if(findstring('dyn',directive,idum))then
        
        keybpd=1
        hyp_units=energy_unit()
        ebias=dblstr(directive,lenrec,idum)
        vmin=dblstr(directive,lenrec,idum)
        call getword(cunit,directive,8,lenrec)
        if(idnode.eq.0)write(nrite,"(
     x    1x,'dynamics option selected       ',
     x    /,1x,'bias potential E_bias  (kelvin)',f10.4,
     x    /,1x,'bias potential V_min   (kelvin)',f10.4
     x    /,1x,'energy units                   ',2x,a8)")
     x    ebias,vmin,cunit
        
      elseif(findstring('path',directive,idum))then
        
        keybpd=2
        nebgo=.true.
        hyp_units=1.d0
        do while(.not.endbpd)
          
          call getrec(safe,idnode,nread)
          if(.not.safe)call abort_control_read(1,idnode,nread)
          call lowcase(record,lenrec)
          call strip(record,lenrec)
          call copystring(record,directive,lenrec)
          
          if(record(1).eq.'#'.or.record(1).eq.'&')then
c     information only - skip record
            cycle
          elseif(findstring('endbpd',directive,idum))then
            endbpd=.true.
          elseif(findstring('pre',directive,idum))then
            prechk=.true.
            if(findstring('false',directive,idum))prechk=.false.
          elseif(findstring('noneb',directive,idum))then
            nebgo=.false.
          elseif(findstring('target',directive,idum))then
            call getword(seek,directive,8,lenrec)
            call getword(seek,directive,8,lenrec)
          elseif(findstring('units',directive,idum))then
            hyp_units=energy_unit()
            call getword(cunit,directive,8,lenrec)
            call getword(cunit,directive,8,lenrec)
          elseif(findstring('ebias',directive,idum))then
            ebias=dblstr(directive,lenrec,idum)
          elseif(findstring('vmin',directive,idum))then
            vmin=dblstr(directive,lenrec,idum)
          elseif(findstring('num_block',directive,idum))then
            nblock=intstr(directive,lenrec,idum)
          elseif(findstring('num_track',directive,idum))then
            ntrack=intstr(directive,lenrec,idum)
          elseif(findstring('catch_radius',directive,idum))then
            catchrad=dblstr(directive,lenrec,idum)
          elseif(findstring('neb_spring',directive,idum))then
            sprneb=dblstr(directive,lenrec,idum)
          elseif(findstring('forc',directive,idum))then
            keytol=0
            opttol=dblstr(directive,lenrec,idum)
          elseif(findstring('ener',directive,idum))then
            keytol=1
            opttol=dblstr(directive,lenrec,idum)
          elseif(findstring('posi',directive,idum))then
            keytol=2
            opttol=dblstr(directive,lenrec,idum)
          endif
          
        enddo
        
        if(idnode.eq.0)then
          
          write(nrite,"(
     x      1x,'dynamics with path analysis selected',
     x      /,1x,'bias potential E_bias  (kelvin)',f10.4,
     x      /,1x,'bias potential V_min   (kelvin)',f10.4,
     x      /,1x,'steps per time block           ',i10,
     x      /,1x,'steps per tracking block       ',i10,
     x      /,1x,'configuration catch radius  (A)',f10.4,
     x      /,1x,'minimisation tolerance         ',e12.4,
     x      /,1x,'atom type to be tracked        ',2x,a8,
     x      /,1x,'energy units                   ',2x,a8)")
     x      ebias,vmin,nblock,ntrack,catchrad,opttol,seek,cunit
          if(nebgo)write(nrite,
     x      "(1x,'NEB spring constant            ',e12.4)")sprneb
          if(prechk)write(nrite,
     x      "(1x,'transition prechecking option selected')")
          call print_optim(keytol)
          
        endif
        
c     energy unit conversions
        
        sprneb=sprneb*hyp_units
        if(keytol.lt.2)opttol=opttol*hyp_units
        
      endif
      
      return
      end subroutine bpd_option
      
      subroutine tad_option
     x  (directive,ltad,lbpd,lminopt,prechk,tadall,idnode,nblock,
     x  ntrack,blkout,keytol,catchrad,sprneb,tlow,deltad,opttol,
     x  hyp_units)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading parameters for TAD option
c     copyright - daresbury laboratory
c     author    - w. smith    feb 2008
c     
c***********************************************************************
      
      implicit none
      
      character*8 cunit
      character*1 directive(lenrec)
      logical ltad,lbpd,lminopt,prechk,tadall,endtad,safe
      integer idnode,nblock,ntrack,blkout,keytol,idum
      real(8) catchrad,sprneb,deltad,tlow,opttol,hyp_units
      
      if(lminopt)call error(idnode,225)
      if(lbpd)call error(idnode,2355)
      lminopt=.true.
      ltad=.true.
      endtad=.false.
      hyp_units=1.d0
      
      do while(.not.endtad)
        
        call getrec(safe,idnode,nread)
        if(.not.safe)call abort_control_read(1,idnode,nread)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        call copystring(record,directive,lenrec)
        
        if(record(1).eq.'#'.or.record(1).eq.'&')then
c     information only - skip record
          cycle
        elseif(findstring('endtad',directive,idum))then
          endtad=.true.
        elseif(findstring('pre',directive,idum))then
          prechk=.true.
          if(findstring('false',directive,idum))prechk=.false.
        elseif(findstring('all',directive,idum))then
          tadall=.true.
          if(findstring('false',directive,idum))tadall=.false.
        elseif(findstring('units',directive,idum))then
          hyp_units=energy_unit()
          call getword(cunit,directive,8,lenrec)
          call getword(cunit,directive,8,lenrec)
        elseif(findstring('num_block',directive,idum))then
          nblock=intstr(directive,lenrec,idum)
        elseif(findstring('num_track',directive,idum))then
          ntrack=intstr(directive,lenrec,idum)
        elseif(findstring('blackout',directive,idum))then
          blkout=intstr(directive,lenrec,idum)
        elseif(findstring('catch_radius',directive,idum))then
          catchrad=dblstr(directive,lenrec,idum)
        elseif(findstring('neb_spring',directive,idum))then
          sprneb=dblstr(directive,lenrec,idum)
        elseif(findstring('deltad',directive,idum))then
          deltad=dblstr(directive,lenrec,idum)
        elseif(findstring('low_temp',directive,idum))then
          tlow=dblstr(directive,lenrec,idum)
        elseif(findstring('forc',directive,idum))then
          keytol=0
          opttol=dblstr(directive,lenrec,idum)
        elseif(findstring('ener',directive,idum))then
          keytol=1
          opttol=dblstr(directive,lenrec,idum)
        elseif(findstring('posi',directive,idum))then
          keytol=2
          opttol=dblstr(directive,lenrec,idum)
        endif
        
      enddo
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,1x,'TAD dynamics controls'
     x    /,1x,'steps per time block           ',i10,
     x    /,1x,'steps per tracking block       ',i10,
     x    /,1x,'steps in blackout periods      ',i10,
     x    /,1x,'configuration catch radius     ',1p,e12.4,
     x    /,1x,'NEB spring constant            ',e12.4,
     x    /,1x,'stopping parameter             ',e12.4,
     x    /,1x,'target low temperature         ',e12.4,
     x    /,1x,'minimisation tolerance         ',e12.4,
     x    /,1x,'energy units                   ',2x,a8)")
     x    nblock,ntrack,blkout,catchrad,sprneb,deltad,
     x    tlow,opttol,cunit
        if(prechk)write(nrite,
     x    "(1x,'transition prechecking option selected')")
        if(tadall)write(nrite,
     x    "(1x,'option for all basins analysis selected')")
        call print_optim(keytol)
        
      endif
      
c     energy unit conversions
      
      sprneb=sprneb*hyp_units
      if(keytol.lt.2)opttol=opttol*hyp_units
      
      return
      end subroutine tad_option

      subroutine metadyn_option
     x  (directive,lmetadyn,lstein,ltet,lglobpe,llocpe,idnode,
     x  ncolvar,nq4,nq6,ntet,hkey,meta_step_int,globpe_scale,
     x  locpe_scale,ref_W_aug,h_aug,wt_Dt)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading parameters for metadynamics option
c     copyright - daresbury laboratory
c     author    - w. smith    jan 2011
c     
c     note: default values are set in metafreeze_module
c     
c***********************************************************************
      
      implicit none
      
      character*1 directive(lenrec)
      logical lmetadyn,endmet,lstein,ltet,lglobpe,llocpe,safe
      integer idnode,idum,ncolvar,nq4,nq6,ntet,hkey,meta_step_int
      real(8) globpe_scale,locpe_scale,ref_W_aug,h_aug,wt_Dt
      
      lmetadyn=.true.
      endmet=.false.          
      
      do while(.not.endmet)
        
        call getrec(safe,idnode,nread)
        if(.not.safe)call abort_control_read(1,idnode,nread)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        call copystring(record,directive,lenrec)
        
        if(record(1).eq.'#'.or.record(1).eq.'&')then
c     information only - skip record
          cycle
        elseif(findstring('endmet',directive,idum))then
          endmet=.true.
        elseif(findstring('ncolvar',directive,idum))then
          ncolvar=intstr(directive,lenrec,idum)
        elseif(findstring('lstein',directive,idum))then
          lstein=.true.
          if(findstring('false',directive,idum))lstein=.false.
        elseif(findstring('ltet',directive,idum))then
          ltet=.true.
          if(findstring('false',directive,idum))ltet=.false.
        elseif(findstring('lglobpe',directive,idum))then
          lglobpe=.true.
          if(findstring('false',directive,idum))lglobpe=.false.
        elseif(findstring('llocpe',directive,idum))then
          llocpe=.true.
          if(findstring('false',directive,idum))llocpe=.false.
        elseif(findstring('globpe_scale',directive,idum))then
          globpe_scale=dblstr(directive,lenrec,idum)
        elseif(findstring('locpe_scale',directive,idum))then
          locpe_scale=dblstr(directive,lenrec,idum)
        elseif(findstring('nq4',directive,idum))then
          nq4=intstr(directive,lenrec,idum)
          nq4=intstr(directive,lenrec,idum) ! do twice - number in name!
        elseif(findstring('nq6',directive,idum))then
          nq6=intstr(directive,lenrec,idum)
          nq6=intstr(directive,lenrec,idum) ! do twice - number in name!
        elseif(findstring('ntet',directive,idum))then
          ntet=intstr(directive,lenrec,idum)
        elseif(findstring('meta_step_int',directive,idum))then
          meta_step_int=intstr(directive,lenrec,idum)
        elseif(findstring('ref_w_aug',directive,idum))then
          ref_W_aug=dblstr(directive,lenrec,idum)
        elseif(findstring('h_aug',directive,idum))then
          h_aug=dblstr(directive,lenrec,idum)
        elseif(findstring('hkey',directive,idum))then
          hkey=intstr(directive,lenrec,idum)
        elseif(findstring('wt_dt',directive,idum))then
          wt_dt=dblstr(directive,lenrec,idum)
        endif
        
      enddo
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,1x,'metadynamics controls'
     x    /,1x,'total number of collective variables',i10,
     x    /,1x,'steinhardt parameters option (Q4/Q6)',l10,
     x    /,1x,'tetrahedral parameters option (zeta)',l10,
     x    /,1x,'global potential parameter option   ',l10,
     x    /,1x,'local potential parameter option    ',l10,
     x    /,1x,'global potential param. scale factor',e12.4,
     x    /,1x,'local potential param. scale factor ',e12.4)")
     x    ncolvar,lstein,ltet,lglobpe,llocpe,globpe_scale,locpe_scale
        
        write(nrite,
     x    "(  1x,'number of Q4 atom pair types        ',i10,
     x      /,1x,'number of Q6 atom pair types        ',i10,
     x      /,1x,'number of zeta atom triplet types   ',i10)")
     x    nq4,nq6,ntet
        
        write(nrite,
     x    "(  1x,'gaussian deposition interval        ',i10,
     x      /,1x,'reference gaussian height           ',e12.4,
     x      /,1x,'gaussian width parameter            ',e12.4,
     x      /,1x,'height control key                  ',i10,
     x      /,1x,'well-tempered control parameter     ',e12.4)")
     x    meta_step_int,ref_W_aug,h_aug,hkey,wt_Dt
        
      endif
            
      return
      end subroutine metadyn_option
      
      subroutine ewald_selection
     x  (directive,lhke,lspme,lewald,lcut,lrfce,kill,idnode,keyfce,
     x  imcon,nhko,nlatt,kmax1,kmax2,kmax3,alpha,rcut)
      
c***********************************************************************
c     
c     dl_poly subroutine for selecting the ewald method and reading 
c     the required parameters
c     copyright - daresbury laboratory
c     author    - w. smith    feb 2008
c     
c***********************************************************************
      
      implicit none
      
      character*1 directive(lenrec)
      logical lhke,lspme,lewald,lcut,lrfce,kill,safe
      integer idnode,keyfce,imcon,nhko,nlatt,kmax1,kmax2,kmax3,idum
      integer kmaxpow2
      real(8) alpha,rcut,eps,tol,fac,tol1
      
      lhke=findstring('hke',directive,idum)
      lspme=findstring('spme',directive,idum)
      lewald=findstring('ewald',directive,idum)
      if(lewald)keyfce=2
      if(lspme)keyfce=12
      if(lhke)keyfce=14
      if(idnode.eq.0)open(nconf,file='CONFIG')
      call getrec(safe,idnode,nconf)
      call getrec(safe,idnode,nconf)
      imcon=intstr(record,lenrec,idum)
      imcon=intstr(record,lenrec,idum)
      if(.not.lhke.and.(imcon.eq.0.or.imcon.eq.6))then
        
        call error(idnode,-180)
        kill=.true.
        
      endif
      
      if(findstring('precision',directive,idum))then
        
        eps=dblstr(directive,lenrec,idum)
        if(idnode.eq.0)write(nrite,
     x    "(/,1x,'Ewald sum  precision    ',7x,1p,e12.4)")eps
        
        if(lhke)then
          
          nhko=min(intstr(directive,lenrec,idum),3)
          nlatt=min(intstr(directive,lenrec,idum),2)
          if(nlatt.eq.0)nlatt=1
          if(nhko.eq.0)nhko=1
          
        endif
        
        if(.not.lcut)then
          call error(idnode,-433)
          kill=.true.
        else
          
c     retreive cell vectors
          
          call getrec(safe,idnode,nconf)
          cell(1)=dblstr(record,lenrec,idum)
          cell(2)=dblstr(record,lenrec,idum)
          cell(3)=dblstr(record,lenrec,idum)
          call getrec(safe,idnode,nconf)
          cell(4)=dblstr(record,lenrec,idum)
          cell(5)=dblstr(record,lenrec,idum)
          cell(6)=dblstr(record,lenrec,idum)
          call getrec(safe,idnode,nconf)
          cell(7)=dblstr(record,lenrec,idum)
          cell(8)=dblstr(record,lenrec,idum)
          cell(9)=dblstr(record,lenrec,idum)
          
c     compute alpha and the kmax
          
          if(lewald.or.lspme)then
            
            call dcell(cell,celprp)
            eps=min(abs(eps),0.5d0)
            tol=sqrt(abs(log(eps*rcut)))
            alpha=sqrt(abs(log(eps*rcut*tol)))/rcut
            tol1=sqrt(-log(eps*rcut*(2.d0*tol*alpha)**2))
            fac=1.d0
            if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)
     x        fac=2.d0**(1.d0/3.d0)
            kmax1=nint(0.25d0+fac*celprp(1)*alpha*tol1/pi)
            kmax2=nint(0.25d0+fac*celprp(2)*alpha*tol1/pi)
            kmax3=nint(0.25d0+fac*celprp(3)*alpha*tol1/pi)
            
          elseif(lhke)then
            
            if(nhko.eq.0)then
              if(eps.le.1.d-6)then
                alpha=3.46d0/rcut
              elseif(eps.le.1.d-5)then
                alpha=3.14d0/rcut
              else
                alpha=2.76d0/rcut
              endif
            elseif(nhko.eq.1)then
              if(eps.le.1.d-6)then
                alpha=4.37d0/rcut
              elseif(eps.le.1.d-5)then
                alpha=4.08d0/rcut
              else
                alpha=3.75d0/rcut
              endif                
            elseif(nhko.eq.2)then
              if(eps.le.1.d-6)then
                alpha=5.01d0/rcut
              elseif(eps.le.1.d-5)then
                alpha=4.74d0/rcut
              else
                alpha=4.44d0/rcut
              endif
            elseif(nhko.eq.3)then
              if(eps.le.1.d-6)then
                alpha=5.55d0/rcut
              elseif(eps.le.1.d-5)then
                alpha=5.28d0/rcut
              else
                alpha=5.00d0/rcut
              endif
            endif
            alpha=alpha/dble(2*nlatt+1)
            if(abs(cell(9)).lt.1.d-8)cell(9)=1.d0
            call dcell(cell,celprp)
            tol=2.d0*alpha*sqrt(abs(log(eps*alpha)))
            tol1=2.d0*alpha*sqrt(abs(log(eps*alpha*tol)))
            kmax1=nint(0.25d0+0.5d0*celprp(1)*tol1/pi)
            kmax2=nint(0.25d0+0.5d0*celprp(2)*tol1/pi)
            kmax3=1
            
          endif
          
        endif
        
      else
        
        alpha=dblstr(directive,lenrec,idum)
        kmax1=intstr(directive,lenrec,idum)
        kmax2=intstr(directive,lenrec,idum)
        
        if(lhke)then
          
          kmax3=1
          nhko=min(intstr(directive,lenrec,idum),3)
          nlatt=min(intstr(directive,lenrec,idum),2)
          
        else
          
          kmax3=intstr(directive,lenrec,idum)
          
        endif
        
      endif
      
c     if spme double kmax and set to next power of 2, with current upper
c     limit of 512.
      
      if(lspme)then
        
        kmaxpow2=1
        do while (kmax1.gt.kmaxpow2.and.kmaxpow2.lt.256)
          kmaxpow2=kmaxpow2 * 2
        enddo
        kmax1=2 * kmaxpow2
        
        kmaxpow2=1
        do while (kmax2.gt.kmaxpow2.and.kmaxpow2.lt.256)
          kmaxpow2=kmaxpow2 * 2
        enddo
        kmax2=2 * kmaxpow2
        
        kmaxpow2=1
        do while (kmax3.gt.kmaxpow2.and.kmaxpow2.lt.256)
          kmaxpow2=kmaxpow2 * 2
        enddo
        kmax3=2 * kmaxpow2
        
      endif
      
      if(idnode.eq.0)then
        
        close(nconf)
        
        if(lspme)then
          
          write(nrite,
     x      "(/,1x,'Electrostatics : SPME  ')")
          
          write(nrite,
     x      "(/,1x,'Ewald convergence parameter    ',1p,e12.4,
     x      /,1x,'Ewald kmax1 kmax2 kmax3     ',3i5)")
     x      alpha,kmax1/2,kmax2/2,kmax3/2
          
        elseif(lhke)then
          
          write(nrite,
     x      "(/,1x,'Electrostatics : Hautman-Klein-Ewald sum  ')")
          
          write(nrite,
     x      "(/,1x,'Ewald convergence parameter    ',1p,e12.4,
     x      /,1x,'Ewald kmax1 kmax2              ',2i5)")
     x      alpha,kmax1,kmax2
          
          write(nrite,
     x      "(1x,'HKE expansion order     ',7x,i10,
     x      /,1x,'HKE lattice control     ',7x,i10)")nhko,nlatt
          
        else
          
          write(nrite,
     x      "(/,1x,'Electrostatics : Ewald sum  ')")
          
          write(nrite,
     x      "(/,1x,'Ewald convergence parameter    ',1p,e12.4,
     x      /,1x,'Ewald kmax1 kmax2 kmax3     ',3i5)")
     x      alpha,kmax1,kmax2,kmax3
          
        endif
        
      endif
      
      if(lspme)then
        
c     Initialize fft tables
        
CFFTW             call fftw3d_f77_create_plan
CFFTW     x           (fplan,kmaxd,kmaxe,kmaxf,
CFFTW     x            FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
CFFTW
CFFTW             call fftw3d_f77_create_plan
CFFTW     x           (bplan,kmaxd,kmaxe,kmaxf,
CFFTW     x            FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE)
        
CSGIC             call zzfft3d( 0,kmaxd,kmaxe,kmaxf,1.d0,dummy,1,1,
CSGIC     x                     dummy,1,1,ffttable,dummy,dummy )
        
CCRAY             call ccfft3d( 0,kmaxd,kmaxe,kmaxf,1.d0,dummy,1,1,
CCRAY     x                     dummy,1,1,ffttable,dummy,dummy )
        
      endif
      
      if(lspme)then
        
        if(kmax1.gt.kmaxd.or.kmax2.gt.kmaxe.or.kmax3.gt.kmaxf)then
          
          kill=.true.
          call error(idnode,-185)
          
        endif
        
      elseif(lhke)then
        
        if(kmax2.gt.kmaxb)then
          
          kill=.true.
          call error(idnode,-185)
          
        endif
        
      else
        
        if(kmax2.gt.kmaxb.or.kmax3.gt.kmaxc)then
          
          kill=.true.
          call error(idnode,-185)
          
        endif
        
      endif
      
      if(lrfce)then
        call  error(idnode,-416)
        kill=.true.
      endif
      lrfce=.true.
      
      return
      end subroutine ewald_selection
      
      subroutine print_optim(keytol)
      
c***********************************************************************
c     
c     dl_poly subroutine for printing the optimisation option
c     the required parameters
c     copyright - daresbury laboratory
c     author    - w. smith    feb 2008
c     
c***********************************************************************
      
      implicit none
      
      integer keytol
      
      if(keytol.eq.0)then
        write(nrite,
     x    "(1x,'convergence to minimum force selected')")
      elseif(keytol.eq.1)then
        write(nrite,
     x    "(1x,'convergence to minimum energy selected')")
      else
        write(nrite,
     x    "(1x,'convergence to minimum position selected')")
      endif
      
      return
      end subroutine print_optim
      
      function energy_unit()
      
c***********************************************************************
c     
c     dl_poly subroutine for assigning energy conversion factors
c     copyright - daresbury laboratory
c     author    - w. smith    feb 2008
c     
c***********************************************************************
      
      implicit none
      
      integer idum
      real(8) energy_unit
      
      energy_unit=1.d0
      if(findstring('ev',record,idum))then
        energy_unit=9648.530821d0
      elseif(findstring('kev',record,idum))then
        energy_unit=9648530.821d0
      elseif(findstring('kcal',record,idum))then
        energy_unit=418.4d0
      elseif(findstring('kj',record,idum))then
        energy_unit=1.d2
      elseif(findstring('k',record,idum))then
        energy_unit=boltz
      endif
      
      return
      end function energy_unit
      
      subroutine solvation_option
     x  (directive,lsolva,idnode,nsolva,isolva)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading parameters for solvation option
c     copyright - daresbury laboratory
c     authors   - w. smith and p.-a. cazade jul 2008
c     
c***********************************************************************
      
      implicit none
      
      character*1 directive(lenrec)
      logical lsolva,endsol,safe
      integer idnode,nsolva,isolva,idum
      
      lsolva=.true.
      endsol=.false.
      
      nsolva=intstr(directive,lenrec,idum)
      isolva=intstr(directive,lenrec,idum)
      
      if(nsolva.eq.0.and.isolva.eq.0)then
        
        do while(.not.endsol)
          
          call getrec(safe,idnode,nread)
          if(.not.safe)call abort_control_read(1,idnode,nread)
          call lowcase(record,lenrec)
          call strip(record,lenrec)
          call copystring(record,directive,lenrec)
          
          if(findstring('endsol',directive,idum))then
            endsol=.true.
          elseif(findstring('enddec',directive,idum))then
            endsol=.true.
          elseif(findstring('start',directive,idum))then
            nsolva=intstr(directive,lenrec,idum)
          elseif(findstring('inter',directive,idum))then
            isolva=max(intstr(directive,lenrec,idum),1)
          endif
          
        enddo
        
      endif
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,1x,'solvation calculation selected',
     x    /,1x,'start of solvation calculation ',i10,
     x    /,1x,'solvation calculation interval ',i10)")
     x    nsolva,isolva
        
      endif
      
      return
      end subroutine solvation_option
      
      subroutine free_energy_option(directive,lfree,lfrmas,idnode)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading parameters for free energy option
c     copyright - daresbury laboratory
c     authors   - w. smith and p.-a. cazade jul 2008
c     
c***********************************************************************
      
      implicit none
      
      character*1 directive(lenrec)
      logical lfree,lfrmas,endfre,safe
      integer idnode,idum
      
      mfree=1
      kfree=1
      lfree=.true.
      lfrmas=.false.
      endfre=.false.
      
      do while(.not.endfre)
        
        call getrec(safe,idnode,nread)
        if(.not.safe)call abort_control_read(1,idnode,nread)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        call copystring(record,directive,lenrec)
        
        if(findstring('endfre',directive,idum))then
          endfre=.true.
        elseif(findstring('start',directive,idum))then
          nfrn=intstr(directive,lenrec,idum)
        elseif(findstring('interval',directive,idum))then
          ifrn=intstr(directive,lenrec,idum)
        elseif(findstring('lambda',directive,idum))then
          pfree=dblstr(directive,lenrec,idum)
        elseif(findstring('mix',directive,idum))then
          mfree=intstr(directive,lenrec,idum)
        elseif(findstring('expo',directive,idum))then
          kfree=intstr(directive,lenrec,idum)
        elseif(findstring('reset_mass',directive,idum))then
          lfrmas=.true.
          if(findstring('false',directive,idum))lfrmas=.false.
        elseif(findstring('system_a',directive,idum))then
          ind_fre(1)=intstr(directive,lenrec,idum)
          ind_fre(2)=intstr(directive,lenrec,idum)
        elseif(findstring('system_b',directive,idum))then
          ind_fre(3)=intstr(directive,lenrec,idum)
          ind_fre(4)=intstr(directive,lenrec,idum)
        endif
        
      enddo
      
      if(mfree.eq.1)kfree=1
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,1x,'free energy option selected',
     x    /,1x,'start of free energy calculation ',i10,
     x    /,1x,'sampling interval                ',i10,
     x    /,1x,'free energy parameter (lambda)   ',f10.3,
     x    /,1x,'mixing rule selected             ',i10,
     x    /,1x,'mixing rule exponent             ',i10,
     x    /,1x,'system A first atom              ',i10,
     x    /,1x,'system A last atom               ',i10,
     x    /,1x,'system B first atom              ',i10,
     x    /,1x,'system B last atom               ',i10,
     x    /,1x,'mass scaling option              ',l10)")
     x    nfrn,ifrn,pfree,mfree,kfree,ind_fre,lfrmas
        
      endif
      
c     define free energy scaling parameters

      call freegen()
      
      return
      end subroutine free_energy_option
      
      subroutine excitation_option
     x  (directive,lsolva,lexcite,lghost,idnode,nsolva,isolva)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading parameters for excitation option
c     copyright - daresbury laboratory
c     authors   - w. smith and p.-a. cazade jul 2008
c     
c***********************************************************************
      
      implicit none
      
      character*1 directive(lenrec)
      logical lsolva,lexcite,lghost,endexc,safe
      integer idnode,nsolva,isolva,idum
      
      lsolva=.true.
      lghost=.true.
      lexcite=.true.
      endexc=.false.
      
      do while(.not.endexc)
        
        call getrec(safe,idnode,nread)
        if(.not.safe)call abort_control_read(1,idnode,nread)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        call copystring(record,directive,lenrec)
        
        if(findstring('endexc',directive,idum))then
          endexc=.true.
        elseif(findstring('start',directive,idum))then
          nsolva=intstr(directive,lenrec,idum)
        elseif(findstring('inter',directive,idum))then
          isolva=intstr(directive,lenrec,idum)
        elseif(findstring('system_a',directive,idum))then
          ind_fre(1)=intstr(directive,lenrec,idum)
          ind_fre(2)=intstr(directive,lenrec,idum)
        elseif(findstring('system_b',directive,idum))then
          ind_fre(3)=intstr(directive,lenrec,idum)
          ind_fre(4)=intstr(directive,lenrec,idum)
        endif
        
      enddo
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,1x,'excitation option selected',
     x    /,1x,'energy decomposition start       ',i10,
     x    /,1x,'energy decomposition interval    ',i10,
     x    /,1x,'system A first atom              ',i10,
     x    /,1x,'system A last atom               ',i10,
     x    /,1x,'system B first atom              ',i10,
     x    /,1x,'system B last atom               ',i10)")
     x   nsolva,isolva,ind_fre
        
      endif
      
      return
      end subroutine excitation_option
      
      subroutine switching_option
     x  (directive,lsolva,lswitch,lghost,idnode,nsolva,isolva)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading parameters for switching option
c     copyright - daresbury laboratory
c     authors   - w. smith and p.-a. cazade jul 2008
c     
c***********************************************************************
      
      implicit none
      
      character*1 directive(lenrec)
      logical lsolva,lswitch,lghost,endswi,safe
      integer idnode,nsolva,isolva,idum
      
      lsolva=.true.
      lghost=.true.
      lswitch=.true.
      endswi=.false.
      niswitch=0
      
      do while(.not.endswi)
        
        call getrec(safe,idnode,nread)
        if(.not.safe)call abort_control_read(1,idnode,nread)
        call lowcase(record,lenrec)
        call strip(record,lenrec)
        call copystring(record,directive,lenrec)
        
        if(findstring('endswi',directive,idum))then
          endswi=.true.
        elseif(findstring('start',directive,idum))then
          nsolva=intstr(directive,lenrec,idum)
        elseif(findstring('inter',directive,idum))then
          isolva=intstr(directive,lenrec,idum)
        elseif(findstring('period',directive,idum))then
          niswitch=max(intstr(directive,lenrec,idum),2)
        elseif(findstring('system_a',directive,idum))then
          ind_fre(1)=intstr(directive,lenrec,idum)
          ind_fre(2)=intstr(directive,lenrec,idum)
        elseif(findstring('system_b',directive,idum))then
          ind_fre(3)=intstr(directive,lenrec,idum)
          ind_fre(4)=intstr(directive,lenrec,idum)
        endif
        
      enddo
      
      if(niswitch.eq.0)niswitch=nsolva
      nswitch=nsolva
      
      if(idnode.eq.0)then
        
        write(nrite,
     x    "(/,1x,'switching option selected',
     x    /,1x,'energy decomposition start       ',i10,
     x    /,1x,'energy decomposition interval    ',i10,
     x    /,1x,'switching period                 ',i10,
     x    /,1x,'system A first atom              ',i10,
     x    /,1x,'system A last atom               ',i10,
     x    /,1x,'system B first atom              ',i10,
     x    /,1x,'system B last atom               ',i10)")
     x   nsolva,isolva,niswitch,ind_fre
        
      endif
      
      return
      end subroutine switching_option
      
      end module define_system_module
      
