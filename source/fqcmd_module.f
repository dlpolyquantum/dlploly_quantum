      module fqcmd_module

c**********************************************************
c        
c     dl_poly_quantum module for fast qausi-centroid MD
c
c     refs: 
c       Fletcher, T. et al. J. Chem. Phys. 2021,
c         155,231101
c       Lawrence, J. E. et al. J. Phys. Chem. B 2023,
c         127,9172-9180
c
c     Author - Nathan London 2024
c        
c
c**********************************************************

      use setup_module, only: mspimd,nbeads,corr,mxbuff,mxvdw,
     x  mxatms,mxbond,mxangl,mxtmls,msbad,mxsite,boltz
      use config_module, only: xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,
     x  buffer,cell,rcell,ltype,list,lentry
      use site_module, only: nummols,numsit,unqatm,wgtsit
      use error_module, only: error
      use parse_module, only: getrec, getword, lenrec, record
      use utility_module, only: loc2,images,global_sum_forces
      use pair_module, only: ilist,xdf,ydf,zdf
      use bonds_module, only: numbonds,listbnd
      use angles_module, only: numang,listang
      use correlation_module, only: mol_gather

      implicit none

      real(8), allocatable, save :: potTables(:,:,:)
      real(8), allocatable, save :: qCent(:,:),qcrdf(:,:)
      real(8), allocatable, save :: qcx(:),qcy(:),qcz(:)
      integer, allocatable, save :: potTypes(:)

      integer, save:: numPot, numPoints, totMol,nqcbnd,nqcang
      integer, save:: nrdfpts
      contains

      subroutine allocate_fqcmd_arrays(idnode,mxnode,natms,
     x  ntbond,ntangl)
c*************************************************************
c
c     dl_poly quantum subroutine to allocate arrays needed
c        for f-QCMD calculations
c        
c     Authors: Nathan London 2024
c        
c*************************************************************
      implicit none

      logical safe
      integer, intent(in) :: idnode,mxnode,natms,ntbond,ntangl
      integer :: i
      integer, dimension(1:7) :: fail

      safe=.true.      
      
      fail(:)=0

c     get total number of molecules in system      
      totMol = 0
c     loop over molecule types
      do i=1,mxtmls
        totMol = totMol + nummols(i)
      enddo
c     allocate arrays      
      allocate (potTables(numPot,3,numPoints+4),stat=fail(1))
      allocate (potTypes(mxvdw),stat=fail(2))
      allocate (qCent(totMol,ntbond+ntangl),stat=fail(3))
      allocate (qcx(mxatms),stat=fail(4))
      allocate (qcy(mxatms),stat=fail(5))
      allocate (qcz(mxatms),stat=fail(6))
      allocate (qcrdf(2*numPot+nqcbnd+nqcang,nrdfpts),stat=fail(7))

c     initialize RDF array to be zero at the start     
      qcrdf(:,:)=0.d0
      if(any(fail.gt.0)) safe=.false.
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,540)
      
      end subroutine allocate_fqcmd_arrays

      subroutine dealloc_fqcmd_arrays(idnode,mxnode)

      implicit none

      logical safe

      integer, intent(in) :: idnode,mxnode
      integer, dimension(1:2) :: fail

      fail(:) = 0
      safe = .true.

      deallocate(potTables,stat=fail(1))
      deallocate(potTypes,stat=fail(2))

      if(any(fail.gt.0)) safe=.false.
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,535)

      end subroutine dealloc_fqcmd_arrays
      
      subroutine read_pot_tables(idnode,mxnode,ntpatm)
c********************************************************
c
c     dl_poly quantum subroutine to read in f-QCMD 
c       correction potential tables
c
c     Authors: Nathan London 2024
c
c********************************************************
      implicit none

      integer, intent(in) :: idnode,mxnode,ntpatm
      
      logical safe
      integer :: i,j,hold,katom1,katom2,jtpatm,keypot
      real(8) :: rconv, econv
      character(20) :: filename
      character(5) :: filenm
      character(8) :: atom1,atom2
      
      potTypes(:)=0
      potTables(:,:,:)=0.d0

c     left over conversion from au to angstrom
      rconv = 1.d0
c     convert energy from kcal/mol to internal energy units      
      econv = 418.4d0 

c     loop over pair interactions
      do i=1,numPot
c       get pair table filename
        write(filenm,'(i1)') i
        filename = 'Pair_' // trim(filenm) // '_pot.table'
        if(idnode.eq.0) open(48,file=filename,status='old')
c       get interaction atom types
        call getrec(safe,idnode,48)
        call getword(atom1,record,8,lenrec)
        call getword(atom2,record,8,lenrec)
       
c       define pair interaction table index based on atom
c         types       
        do jtpatm=1,ntpatm
          if(atom1.eq.unqatm(jtpatm)) katom1=jtpatm
          if(atom2.eq.unqatm(jtpatm)) katom2=jtpatm
        enddo
        keypot = loc2(katom1,katom2)
        potTypes(keypot)=i
        
c       read in table        
        if(idnode.eq.0) then
          do j=1,numPoints
            read(48,*) hold, potTables(i,:,j)
          enddo
          close(48)
        endif
      enddo
c     convert to internal units      
      if(idnode.eq.0) then
        potTables(:,1,:) = potTables(:,1,:)*rconv
        potTables(:,2,:) = potTables(:,2,:)*econv
        potTables(:,3,:) = potTables(:,3,:)*econv/rconv
      endif

c     send potential tables to all nodes      
      if(mxnode.gt.1) then
        do i=1,numPot
          do j=1,3
            call gdsum(potTables(i,j,:),numPoints,buffer)
          enddo
        enddo
      endif
      end subroutine read_pot_tables

      subroutine fqcmd_setup(idnode,mxnode,natms,ntbond,
     x  ntangl,ntpatm,lpimd)
c********************************************************
c
c     dl_poly quantum subroutine to setup f-QCMD simulation 
c       at the start of a run
c
c     Authors: Nathan London 2024
c
c********************************************************

      implicit none

      logical, intent(in) :: lpimd
      integer, intent(in) :: idnode,mxnode,natms,ntbond
      integer, intent(in) :: ntangl,ntpatm

      call allocate_fqcmd_arrays(idnode,mxnode,natms,ntbond,
     x  ntangl)

      call read_pot_tables(idnode,mxnode,ntpatm)

      end subroutine fqcmd_setup

      subroutine fqcmd_correct_force(idnode,mxnode,imcon,natms,
     x  engsrp)
c********************************************************
c
c     dl_poly quantum subroutine to add the correction
c       potential and force for nonbonded interations
c       in an f-QCMD simulation through interpolation 
c       of correction potential/force tables
c
c     Authors: Nathan London 2024
c
c********************************************************
      
      integer, intent(in) :: idnode,mxnode,imcon,natms
      real(8), intent(inout) :: engsrp

      integer :: interpJ, iatm,jatm,j,k,l,m,ii
      real(8) :: engfqc,engacc,gamma,ai,aj,ab,rab,rrab,rcut,fi(3)
      real(8) :: fx,fy,fz
      real(8) :: rdr,ppp,vk0,vk1,vk2,gk0,gk1,gk2,t1,t2

      rcut = potTables(1,1,numPoints)
      rdr = 1.d0/(potTables(1,1,2)-potTables(1,1,1))
      engfqc = 0.d0
      ii=0
      
c     loop over atams      
      do iatm=idnode+1,natms,mxnode
        ii=ii+1
        fi(1) = fxx(iatm)
        fi(2) = fyy(iatm)
        fi(3) = fzz(iatm)
c       loop over atoms in neighbor list and calc distance
        do k=1,lentry(ii)
          j=list(ii,k)
          ilist(k) = j
          xdf(k) = xxx(iatm) - xxx(j) 
          ydf(k) = yyy(iatm) - yyy(j) 
          zdf(k) = zzz(iatm) - zzz(j) 
        enddo
c       PBC correction to distances
        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)

c       determine table to use based on atom types        
        ai = dble(ltype(iatm))
        do m=1,lentry(ii)
          jatm=ilist(m)
          aj = dble(ltype(jatm))

          if(ai.gt.aj) then
            ab = ai*(ai-1.d0)*0.5d0+aj+0.5d0
          else
            ab = aj*(aj-1.d0)*0.5d0+ai+0.5d0
          endif

          k=potTypes(int(ab))
     
c         calcuate length of distance vector         
          if(k.ne.0) then
            rrab = 0.d0
            rab = sqrt(xdf(m)**2+ydf(m)**2+zdf(m)**2)
            if(rab.gt.1d-6) rrab=1.d0/rab

c           if distance within cutoff calculate potential
c             and force through interpolation of tables            
            if(rab.lt.rcut) then
              l = int(rab*rdr)
              ppp = rab*rdr-dble(l)

              if(l.eq.0) then
                engacc = potTables(k,2,1)
                gamma = potTables(k,3,1)
              else
                vk0 = potTables(k,2,l)
                vk1 = potTables(k,2,l+1)
                vk2 = potTables(k,2,l+2)
                t1 = vk0+(vk1-vk0)*ppp
                t2 = vk1+(vk2-vk1)*(ppp-1.d0)
                engacc = t1+(t2-t1)*ppp*0.5d0
                
                gk0 = potTables(k,3,l)
                gk1 = potTables(k,3,l+1)
                gk2 = potTables(k,3,l+2)
                t1 = gk0+(gk1-gk0)*ppp
                t2 = gk1+(gk2-gk1)*(ppp-1.d0)
                gamma = (t1+(t2-t1)*ppp*0.5d0)
              endif
            else
              gamma = 0.d0
              engacc = 0.d0
            endif
          else 
            gamma = 0.d0
            engacc = 0.d0
          endif

c         get vector forces          
          gamma = gamma*rrab
          fx = gamma*xdf(m)
          fy = gamma*ydf(m)
          fz = gamma*zdf(m)
          
c         add correction forces           
          fi(1) = fi(1) + fx
          fi(2) = fi(2) + fy
          fi(3) = fi(3) + fz

          fxx(jatm) = fxx(jatm) - fx
          fyy(jatm) = fyy(jatm) - fy
          fzz(jatm) = fzz(jatm) - fz

          engfqc = engfqc + engacc
        enddo

        fxx(iatm) = fi(1)
        fyy(iatm) = fi(2)
        fzz(iatm) = fi(3)
       enddo
      
      end subroutine fqcmd_correct_force
     
      subroutine get_qcent_internal(idnode,mxnode,imcon,natms,
     x  ntbond,ntangl,ntpmls,rcut)
c********************************************************
c
c     dl_poly quantum subroutine to calculate the 
c       curvilinear coordinates of the quasi-centroids
c       of each molecule
c     (currently limited to bonds and angles)
c
c     Authors: Nathan London 2024
c
c********************************************************

      implicit none

      integer, intent(in) :: idnode,mxnode,imcon,natms
      integer, intent(in) :: ntbond,ntangl,ntpmls
      real(8), intent(in) :: rcut

      logical :: safe
      integer :: ii,itmol,imol,imol0,imolbnd,ibnd1,ibnd2,ibase
      integer :: imolang,iang1,iang2,ibnd,iang
      integer :: i,j,k, ia,ib,ic,fail(2)
  
      real(8) :: rab,rbc,xab,xbc,yab,ybc,zab,zbc
      real(8) :: cost,theta,pi,rdr,rtheta
      real(8), allocatable :: xdab(:),ydab(:),zdab(:)
      real(8), allocatable :: xdbc(:),ydbc(:),zdbc(:)
      real(8) :: bondVal(ntbond), angVal(ntangl)

c     allocate temporary arrays for atom distances      
      safe = .true.
      fail(:) = 0

      allocate(xdab(msbad),ydab(msbad),zdab(msbad),stat=fail(1))
      allocate(xdbc(msbad),ydbc(msbad),zdbc(msbad),stat=fail(2))
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,540)

c     get indices for parallel calculations      
      ibnd1 = (idnode*ntbond)/mxnode+1
      ibnd2 = ((idnode+1)*ntbond)/mxnode
      iang1 = (idnode*ntangl)/mxnode+1
      iang2 = ((idnode+1)*ntangl)/mxnode

      imol0 = (idnode*totmol)/mxnode+1

c     define histogram bin spacing      
      rdr = (1.3225d0)/dble(nrdfpts)
      pi = 4.0*atan(1.0)
      rtheta = pi/dble(nrdfpts)

      qCent(:,:) = 0.d0
      bondVal(:) = 0.d0
      angVal(:) = 0.d0

c     loop over beads to get bond distances      
      do k=1,nbeads
        ibase = (k-1)*natms
       
        ii=0

        do i=ibnd1,ibnd2
          ii = ii + 1

          ia = listbnd(ii,2) + ibase
          ib = listbnd(ii,3) + ibase
        
          xdab(ii)=xxx(ia)-xxx(ib)
          ydab(ii)=yyy(ia)-yyy(ib)
          zdab(ii)=zzz(ia)-zzz(ib)
        enddo

        call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
        
        ii=0

        do i=ibnd1,ibnd2
          ii = ii + 1

          rab = sqrt((xdab(ii)**2+ydab(ii)**2)+zdab(ii)**2)
c       store bond distances in temporary array
          bondVal(i) = bondVal(i) + rab/dble(nbeads)
        enddo
 
c       calculate angle values        
        ii=0

        do i=iang1,iang2
          ii = ii + 1

          ia = listang(ii,2) + ibase
          ib = listang(ii,3) + ibase
          ic = listang(ii,4) + ibase 

          xdab(ii)=xxx(ia)-xxx(ib)
          ydab(ii)=yyy(ia)-yyy(ib)
          zdab(ii)=zzz(ia)-zzz(ib)
          xdbc(ii)=xxx(ic)-xxx(ib)
          ydbc(ii)=yyy(ic)-yyy(ib)
          zdbc(ii)=zzz(ic)-zzz(ib)
        enddo

        call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
        call images(imcon,0,1,ii,cell,xdbc,ydbc,zdbc)
        
        ii=0

        do i=iang1,iang2
          ii = ii + 1

          rab = sqrt((xdab(ii)**2+ydab(ii)**2)+zdab(ii)**2)
          xab = xdab(ii)/rab
          yab = ydab(ii)/rab
          zab = zdab(ii)/rab
          
          rbc = sqrt((xdbc(ii)**2+ydbc(ii)**2)+zdbc(ii)**2)
          xbc = xdbc(ii)/rbc
          ybc = ydbc(ii)/rbc
          zbc = zdbc(ii)/rbc

          cost = xab*xbc + yab*ybc + zab*zbc
          theta = acos(cost)
c       store angle values in temporary array
          angVal(i) = angVal(i) + theta/dble(nbeads)
        enddo
      enddo

c     merge temporary arrays over processors      
      if(mxnode.gt.1) then
        call gdsum(bondVal,ntbond,buffer)
        call gdsum(angVal,ntangl,buffer)
      endif

c     store values in molecule partitioned QCMD arrays        
      imol = 0
      ibnd=0
      iang=0
      do i=1,ntpmls
        do j=1,nummols(i)
          imol = imol+1
          imolbnd=0
          imolang=0
          do k=1,numbonds(i)
            imolbnd = imolbnd + 1
            ibnd = ibnd +1
            qCent(imol,imolbnd) = bondVal(ibnd)
          enddo
          do k=1,numang(i)
            imolang = imolang + 1
            iang = iang +1
            qCent(imol,ntbond+imolang) = angVal(iang)
          enddo
        enddo
      enddo

c     histogram internal coordinates
c     (hardcoded for water currently)      
      do i=1,nummols(1)
        call qchistogram(1,qCent(i,1),0.d0,rdr)
        call qchistogram(1,qCent(i,2),0.d0,rdr)
        call qchistogram(2,qCent(i,ntbond+1),0.d0,rtheta)
      enddo
      safe = .true.
      fail = 0

      deallocate(xdab,ydab,zdab,stat=fail(1))
      deallocate(xdbc,ydbc,zdbc,stat=fail(2))
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,540)
      end subroutine get_qcent_internal

      subroutine get_qcent_cart(idnode,mxnode,imcon,
     x  ntpmls,natms,ntbond,rcut,lpimd)
c**********************************************************
c        
c     dl_poly_quantum subroutine for calculating the 
c       cartesian coordinates of the quasi-centroids through
c       throught eckart-like frame rotation with respect
c       to the cartesian centroids        
c
c     refs: 
c       Lawrence, J. E.; et al. J. Phys. Chem. B 2023,
c         127,9172-9180.
c       Trenins, G.; Haggard, C.; Althorpe, S. C. J. Chem.
c         Phys. 2022, 151, 174108.           
c       Krasnoshchekov, S. V.; Isayeva, E. V.; Stpanov, N.F.   
c         J. Chem. Phys. 2014, 140, 154104.
c           
c     Author - Nathan London 2024
c        
c
c**********************************************************

      implicit none

      logical, intent(in) :: lpimd
      integer, intent(in) :: idnode,mxnode,imcon,ntpmls
      integer, intent(in) :: natms,ntbond
      real(8), intent(in) :: rcut

      logical :: safe
      integer :: fail(2),ierr
      integer :: itmols,isite,imol,imol0,imol1
      integer :: ii,i,j,k,m,jsite,ksite,copystart,copyend
      integer :: ai,aj,iatm,jatm,ab
      real(8) :: molmass,centcomx,centcomy,centcomz
      real(8) :: qcomx,qcomy,qcomz,rad,theta,rdr,rab  
      real(8) :: mass(mxsite),centx(mxsite),centy(mxsite)
      real(8) :: centz(mxsite),cmat(10),rotmat(3,3),perdif(3)
      real(8) atmx(mxsite),atmy(mxsite)
      real(8) atmz(mxsite)
      real(8), allocatable :: xtmp(:),ytmp(:),ztmp(:)

c      write(*,*) "get cart"
      safe = .true.
      fail(:) = 0 
      allocate(xtmp(mxsite),ytmp(mxsite),stat=fail(1))
      allocate(ztmp(mxsite),stat=fail(2))
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,540)
      
      jsite = 0
      ksite = 0

c     define histrogram bin width
      rdr = rcut/dble(nrdfpts)

c     perform eckart rotation for pimd with more than 1 bead      
      if(lpimd.and.(nbeads.gt.1)) then

c       loop over molecule types        
        do itmols=1,ntpmls

c         get parallization indiced
          imol0 = (idnode*nummols(itmols))/mxnode+1
          imol1 = ((idnode+1)*nummols(itmols))/mxnode
          
          if(itmols.gt.1) then
            jsite = jsite + nummols(itmols-1)*
     x      numsit(itmols-1)
            ksite = ksite + numsit(itmols-1)
          endif

c         get molecule mass for COM calc          
          molmass=0.d0
          do isite=1,numsit(itmols)
            mass(isite) = wgtsit(ksite+isite)
            molmass = molmass + mass(isite)
          enddo
          
          qcx(:) = 0.d0
          qcy(:) = 0.d0
          qcz(:) = 0.d0
          do imol=imol0,imol1
            centx(:) = 0.d0
            centy(:) = 0.d0
            centz(:) = 0.d0
            centcomx = 0.d0
            centcomy = 0.d0
            centcomz = 0.d0
            qcomx = 0.d0
            qcomy = 0.d0
            qcomz = 0.d0

c           loop over atoms in molecule            
            do j=1,numsit(itmols)
c             loop over beads to get cartesian centroid            
              do k=1,nbeads
                centx(j) = centx(j) +
     x          xxx((k-1)*natms+(imol-1)
     x          *numsit(itmols)+jsite+j)
                centy(j) = centy(j) +
     x          yyy((k-1)*natms+(imol-1)
     x          *numsit(itmols)+jsite+j)
                centz(j) = centz(j) +
     x          zzz((k-1)*natms+(imol-1)
     x          *numsit(itmols)+jsite+j)
              enddo
              centx(j) = centx(j)/dble(nbeads)
              centy(j) = centy(j)/dble(nbeads)
              centz(j) = centz(j)/dble(nbeads)
c             make sure molecule is whole for PBC              
              call mol_gather(numsit(itmols),
     x            centx,centy,centz)

c             get centroid COM
              centcomx = centcomx + mass(j)*centx(j)
              centcomy = centcomy + mass(j)*centy(j)
              centcomz = centcomz + mass(j)*centz(j)

c             define temporay quasi-centroid cartesian
c             coordinates based on internal coordinates
c             (works for 3 atom molecules only)              
              if(j.eq.1)then
                xtmp(j) = 0.d0
                ytmp(j) = 0.d0
                ztmp(j) = 0.d0
              else if (j.eq.2) then
                xtmp(j) = 0.d0
                ytmp(j) = 0.d0
                ztmp(j) = qCent(imol,1)
              else if (j.eq.3) then
                rad = qCent(imol,2)
                theta = qCent(imol,ntbond+1)
                xtmp(j) = rad*sin(theta)
                ytmp(j) = 0.d0
                ztmp(j) = rad*cos(theta)
              else
                xtmp(j) = 0.d0
                ytmp(j) = 0.d0
                ztmp(j) = 0.d0
              endif

c             calculate quasi-centroid COM              
              qcomx = qcomx + mass(j)*xtmp(j)
              qcomy = qcomy + mass(j)*ytmp(j)
              qcomz = qcomz + mass(j)*ztmp(j)

            enddo
            
            centcomx = centcomx/molmass
            centcomy = centcomy/molmass
            centcomz = centcomz/molmass
            qcomx = qcomx/molmass
            qcomy = qcomy/molmass
            qcomz = qcomz/molmass

c           shift both coordinates so COM at origin            
            centx(:) = centx(:) - centcomx
            centy(:) = centy(:) - centcomy
            centz(:) = centz(:) - centcomz
           
            xtmp(:) = xtmp(:) - qcomx   
            ytmp(:) = ytmp(:) - qcomy   
            ztmp(:) = ztmp(:) - qcomz    
        
c           rotate the quasi-centroid coordinates with
c             eckart-frame rotation      
            call get_cmat(numsit(itmols),centx,centy,centz,
     x      xtmp,ytmp,ztmp,mass,cmat)
            call get_rotmat(cmat,rotmat) 
            call rotate(numsit(itmols),rotmat,xtmp,ytmp,ztmp)

c           shift both coordinates so COM is at original
c           centroid COM            
            xtmp(:) = xtmp(:) + centcomx   
            ytmp(:) = ytmp(:) + centcomy   
            ztmp(:) = ztmp(:) + centcomz    
            centx(:) = centx(:) + centcomx
            centy(:) = centy(:) + centcomy
            centz(:) = centz(:) + centcomz

c           copy quasi-centroid coords to main arrays            
            copystart = (imol-1)*numsit(itmols)+jsite+1
            copyend = (imol-1)*numsit(itmols)+jsite+numsit(itmols)
            qcx(copystart:copyend) = xtmp(1:numsit(itmols))
            qcy(copystart:copyend) = ytmp(1:numsit(itmols))
            qcz(copystart:copyend) = ztmp(1:numsit(itmols))
          enddo
        enddo

c       merge coords over nodes        
        if(mxnode.gt.1) then
          call gsync()
           call gdsum(qcx,natms,buffer)
           call gdsum(qcy,natms,buffer)
           call gdsum(qcz,natms,buffer)
        endif
 
c     for f-QCMD (not PIMD) simulations just use cartesian coords      
      else
        qcx(:) = xxx(:)
        qcy(:) = yyy(:)
        qcz(:) = zzz(:)
      endif

c     calculate non-bonded distances for RDFs      
      ii=0
      do iatm=idnode+1,natms,mxnode
        ii=ii+1
        do k=1,lentry(ii)
          j=list(ii,k)
          ilist(k) = j
          xdf(k) = qcx(iatm) - qcx(j) 
          ydf(k) = qcy(iatm) - qcy(j) 
          zdf(k) = qcz(iatm) - qcz(j) 
        enddo

        call images(imcon,0,1,lentry(ii),cell,xdf,ydf,zdf)
        ai = dble(ltype(iatm))
        do m=1,lentry(ii)
          jatm=ilist(m)
          aj = dble(ltype(jatm))

          if(ai.gt.aj) then
            ab = ai*(ai-1.d0)*0.5d0+aj+0.5d0
          else
            ab = aj*(aj-1.d0)*0.5d0+ai+0.5d0
          endif

          k=potTypes(int(ab))

c         call histogram subroutines
c           for PIMD use standard binning
c           for f-QCMD use force matching          
          if(k.ne.0) then
            rab = sqrt(xdf(m)**2+ydf(m)**2+zdf(m)**2)
            if(lpimd.and.(nbeads.gt.1)) then
              call qchistogram(k+nqcbnd+nqcang,rab,0.d0,rdr)
            else
              call qc_force_hist(nqcbnd+nqcang+numPot+k,iatm,
     x          jatm,m,rab,0.d0,rdr)
            endif
          endif
        enddo
      enddo

      safe = .true.
      fail(:) = 0
      deallocate(xtmp,ytmp,stat=fail(1))
      deallocate(ztmp,stat=fail(2))
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,540)
      end subroutine get_qcent_cart
      
      subroutine get_cmat(numsit,centx,centy,centz,
     x  xtmp,ytmp,ztmp,mass,cmat)
c**********************************************************
c        
c     dl_poly_quantum subroutine for calculating the 
c       matrix of mass weighed difference between 
c       quasi-centroid and centroid positions
c
c     refs: 
c       Krasnoshchekov, S. V.; Isayeva, E. V.; Stpanov, N.F.   
c         J. Chem. Phys. 2014, 140, 154104.
c           
c     Author - Nathan London 2024
c        
c
c**********************************************************

        implicit none

        integer, intent(in) :: numsit
        real(8), intent(in) :: mass(mxsite)
        real(8), intent(in) :: centx(mxsite)
        real(8), intent(in) :: centy(mxsite)
        real(8), intent(in) :: centz(mxsite)
        real(8), intent(in) :: xtmp(mxsite)
        real(8), intent(in) :: ytmp(mxsite)
        real(8), intent(in) :: ztmp(mxsite)
        real(8), intent(inout) :: cmat(10)
        
        integer :: i
        real(8) :: xplus,xminus,yplus,yminus
        real(8) :: zplus,zminus
        
        cmat(:) = 0.d0
c       calulate sums and differnces  
        do i=1,numsit
          xplus = centx(i) + xtmp(i)
          xminus = centx(i) - xtmp(i)
          yplus = centy(i) + ytmp(i)
          yminus = centy(i) - ytmp(i)
          zplus = centz(i) + ztmp(i)
          zminus = centz(i) - ztmp(i)

c       construct c matrix          
          cmat(1) = cmat(1) + mass(i)*
     x      (xminus**2+yminus**2+zminus**2)
          cmat(3) = cmat(3) + mass(i)*
     x      (xminus**2+yplus**2+zplus**2)
          cmat(6) = cmat(6) + mass(i)*
     x      (xplus**2+yminus**2+zplus**2)
          cmat(10) = cmat(10) + mass(i)*
     x      (xplus**2+yplus**2+zminus**2)

          cmat(2) = cmat(2) + mass(i)*
     x      (yplus*zminus - yminus*zplus)
          cmat(4) = cmat(4) + mass(i)*
     x      (xminus*zplus - xplus*zminus)
          cmat(7) = cmat(7) + mass(i)*
     x      (xplus*yminus - xminus*yplus)

          cmat(5) = cmat(5) + mass(i)*
     x      (xminus*yminus - xplus*yplus)
          cmat(8) = cmat(8) + mass(i)*
     x      (xminus*zminus - xplus*zplus)
          cmat(9) = cmat(9) + mass(i)*
     x      (zminus*yminus - zplus*yplus)
        enddo

      end subroutine get_cmat
     
      subroutine get_rotmat(cmat,rotmat)
c**********************************************************
c        
c     dl_poly_quantum subroutine for calculating the 
c       rotation matrix for the eckart-like frame 
c
c     refs: 
c       Krasnoshchekov, S. V.; Isayeva, E. V.; Stpanov, N.F.   
c         J. Chem. Phys. 2014, 140, 154104.
c           
c     Author - Nathan London 2024
c        
c
c**********************************************************

        implicit none

        real(8), intent(in) :: cmat(10)
        real(8), intent(inout) :: rotmat(3,3)

        integer :: info,i,j
        real(8) :: quat(4),eigval(4),eigvec(4,4),work(12)

c       get eigenvalues/vectors of c matrix        
        call DSPEV('V','U',4,cmat,eigval,eigvec,4,work,info)

c       store eigenvector for lowest eigenvalue as quaternion        
        if(info.eq.0) then
          quat(:) = eigvec(:,1)
        else
          call error(1,540)
        endif

c       construct rotation matrix from quaternion        
        rotmat(1,1) = quat(1)**2+quat(2)**2-quat(3)**2
     x    -quat(4)**2
        rotmat(2,2) = quat(1)**2-quat(2)**2+quat(3)**2
     x    -quat(4)**2
        rotmat(3,3) = quat(1)**2-quat(2)**2-quat(3)**2
     x    +quat(4)**2

        rotmat(2,1) = 2.d0*(quat(2)*quat(3)+quat(1)*quat(4))
        rotmat(3,1) = 2.d0*(quat(2)*quat(4)-quat(1)*quat(3))
        rotmat(1,2) = 2.d0*(quat(2)*quat(3)-quat(1)*quat(4))
        rotmat(3,2) = 2.d0*(quat(3)*quat(4)+quat(1)*quat(2))
        rotmat(1,3) = 2.d0*(quat(2)*quat(4)+quat(1)*quat(3))
        rotmat(2,3) = 2.d0*(quat(3)*quat(4)-quat(1)*quat(2))

      end subroutine get_rotmat

      subroutine rotate(numsit,rotmat,xtmp,ytmp,ztmp)
c**********************************************************
c        
c     dl_poly_quantum subroutine for rotating molecules 
c       into eckart-like frame
c
c     Author - Nathan London 2024
c        
c
c**********************************************************

        implicit none
        integer, intent(in) :: numsit
        real(8), intent(in) :: rotmat(3,3)
        real(8), intent(inout) :: xtmp(mxsite)
        real(8), intent(inout) :: ytmp(mxsite)
        real(8), intent(inout) :: ztmp(mxsite)

        integer :: i
        real(8) :: x,y,z

        do i=1,numsit
          x = rotmat(1,1)*xtmp(i)+rotmat(2,1)*ytmp(i)
     x      + rotmat(3,1)*ztmp(i)
          y = rotmat(1,2)*xtmp(i)+rotmat(2,2)*ytmp(i)
     x      + rotmat(3,2)*ztmp(i)
          z = rotmat(1,3)*xtmp(i)+rotmat(2,3)*ytmp(i)
     x      + rotmat(3,3)*ztmp(i)
          xtmp(i) = x
          ytmp(i) = y
          ztmp(i) = z
        enddo
      end subroutine rotate

      subroutine qchistogram(rdf,val,xMin,dBins)
c**********************************************************
c        
c     dl_poly_quantum subroutine for histograming using
c       standard binning
c
c     Author - Nathan London 2024
c        
c
c**********************************************************

      implicit none

      integer, intent(in) :: rdf
      real(8), intent(in) :: val,xMin,dBins

      integer :: bin

      bin = int((val-xMin)/dBins+1)
      if(bin.eq.1) then
        write(*,*) "rVal", val
      endif
      if(bin.gt.0.and.bin.le.nrdfpts) then
        qcrdf(rdf,bin) = qcrdf(rdf,bin) + 1.d0/dBins
      endif

      end subroutine qchistogram

      subroutine qc_force_hist(rdf,iatm,jatm,m,val,xMin,dBins)
c**********************************************************
c        
c     dl_poly_quantum subroutine for calculating the 
c       radial distribution function using force matching 
c
c     refs: 
c       Rotenberg, B. J. Chem. Phys. 2020, 155, 150902   
c           
c     Author - Nathan London 2024
c        
c
c**********************************************************

      implicit none

      integer, intent(in) :: rdf,iatm,jatm,m
      real(8), intent(in) :: val,xMin,dBins

      integer :: i,bin
      real(8) :: fdx,fdy,fdz,factor,rconv,econv

c     calculate force difference between atom pair      
      fdx = (fxx(iatm)-fxx(jatm))
      fdy = (fyy(iatm)-fyy(jatm))
      fdz = (fzz(iatm)-fzz(jatm))

      factor = 0.5d0*(fdx*xdf(m)+fdy*ydf(m)+fdz*zdf(m))
     x  /(val**3)

c     add value to bins according to heaviside function      
      bin = int((val-xMin)/dBins+1)+1
      if(bin.gt.0.and.bin.le.nrdfpts) then
        qcrdf(rdf,bin:nrdfpts) = qcrdf(rdf,bin:nrdfpts)
     x    + factor
      endif

      end subroutine qc_force_hist
  
      subroutine write_fqcmd_rdfs(idnode,mxnode,nsteql,
     x  nstrun,rcut,volm,temp,lpimd)
c**********************************************************
c        
c     dl_poly_quantum subroutine for writing out distribution 
c       fucntions to match with format of IBI code
c
c     refs: 
c       Rotenberg, B. J. Chem. Phys. 2020, 155, 150902   
c       Lawrence, J. E. et al. J. Phys. Chem. B 2023,
c         127,9172-9180
c
c     Author - Nathan London 2024
c
c**********************************************************

      implicit none

      logical, intent(in) :: lpimd
      integer, intent(in) :: idnode,mxnode,nsteql,nstrun
      real(8), intent(in) :: rcut,volm,temp   

      integer :: i,j
      real(8) :: ncount,pi,rdrbnd,rdrvdw,rtheta,rconv,r
      real(8) :: dens,factor(2*numPot),econv

      rdrbnd = (1.3225d0)/dble(nrdfpts)
      rdrvdw = rcut/dble(nrdfpts)
      pi = 4.0*atan(1.0)
      rtheta = pi/dble(nrdfpts)
      rconv = 0.529177210903
      econv = 262550.02  
      dens = dble(nummols(1))/volm

c     merge RDFs over nodes      
      if(mxnode.gt.1) then
        do i=1,2*numPot+nqcbnd+nqcang
          call gdsum(qcrdf(i,:),nrdfpts,buffer)
        enddo
      endif

c     calculate normalization factors depending on RDF
c       type      
      factor(:) = 4.d0*pi*dble(nummols(1))*dens
     x  *dble(nstrun-nsteql)
      do i=1,2*numPot
        if (i.eq.1) then
          factor(i) = factor(i)/2.d0
        else if (i.eq.2.or.i.eq.3) then
          factor(i) = factor(i)*2.d0
        else if (i.eq.4) then
          factor(i) = 2.0d0/(boltz*temp*4.d0*pi)*volm/
     x      (dble(nummols(1))*dble(nummols(1)))
     x      /dble(nstrun-nsteql)
        else if (i.eq.5) then
          factor(i) = 1.d0/(boltz*temp*4.d0*pi)*volm/
     x      (dble(nummols(1))*2.d0*dble(nummols(1)))
     x      /dble(nstrun-nsteql)
        else if (i.eq.6) then
          factor(i) = 2.0d0/(boltz*temp*4.d0*pi)*volm/
     x      (2.d0*dble(nummols(1))*2.d0*dble(nummols(1)))
     x      /dble(nstrun-nsteql)
        endif
      enddo

c     normalize intramolecular distributions to have unit
c       area      
      call trapezoid_int(nrdfpts,rdrbnd,qcrdf(1,:),ncount)
      qcrdf(1,:) = qcrdf(1,:)/ncount
      call trapezoid_int(nrdfpts,rtheta,qcrdf(2,:),ncount)
      qcrdf(2,:) = qcrdf(2,:)/ncount
     
c     write out to files      
      if(idnode.eq.0) then
        open(48,file='Average_Intra.d',status='replace')
        do i=1,nrdfpts
          write(48,'(2e14.6)') (dble(i-0.5)*
     x      rdrbnd),qcrdf(1,i)
        enddo
        close(48)
        
        open(48,file='Average_Angle.d',status='replace')
        do i=1,nrdfpts
          write(48,'(2e14.6)') (dble(i-0.5)*rtheta),qcrdf(2,i)
        enddo
        close(48)
        
        open(48,file='Average_RDF.d',status='replace')
        if(lpimd)then
          do i=1,nrdfpts
            r = (dble(i-0.5)*rdrvdw)
            write(48,'(1e14.6)',advance='no') r
            do j=1,numPot
              write(48,'(1e14.6)',advance='no')
     x          qcrdf(nqcbnd+nqcang+j,i)/factor(j)/
     x          (r**2)
            enddo
            write(48,*) " "
          enddo
        else
          do i=1,nrdfpts
            r = (dble(i-.5)*rdrvdw)
            write(48,'(1e14.6)',advance='no') r
            do j=numPot+1,2*numPot
              write(48,'(1e14.6)',advance='no')
     x          qcrdf(nqcbnd+nqcang+j,i)*factor(j)
            enddo
            write(48,*) " "
          enddo
        endif
        close(48)

      endif
      
      end subroutine write_fqcmd_rdfs
    
      subroutine trapezoid_int(n,dx,func,integrand)
        implicit none
c**********************************************************
c        
c     dl_poly_quantum subroutine for calculating numerical
c       integral using trapezoid method 
c       
c     Author - Nathan London 2024
c
c**********************************************************

        integer :: i 
        integer :: n  
        real(8) :: dx, func(n)  
        real(8) :: integrand

        integrand = 0.0d0
        do i = 2, n
          integrand = integrand + (func(i-1)+func(i)) 
     x    * dx/2.0d0
        enddo

      end subroutine trapezoid_int
      
      end module fqcmd_module
