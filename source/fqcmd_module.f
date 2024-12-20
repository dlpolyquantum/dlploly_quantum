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
     x  mxatms,mxbond,mxangl,mxtmls,msbad,mxsite,boltz,nrite,nfield,
     x  pi
      use config_module, only: xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz,
     x  buffer,cell,rcell,ltype,list,lentry
      use site_module, only: nummols,numsit,unqatm,wgtsit,dens
      use error_module, only: error
      use parse_module, only: getrec, getword, lenrec, record
      use utility_module, only: loc2,images,global_sum_forces
      use pair_module, only: ilist,xdf,ydf,zdf
      use bonds_module, only: numbonds,listbnd,prmbnd
      use angles_module, only: numang,listang,prmang
      use correlation_module, only: mol_gather
      use parse_module
c      use define_system_module, only: abort_field_read

      implicit none

      real(8), allocatable, save :: potTables(:,:,:)
      real(8), allocatable, save :: qCent(:,:),qcrdf(:,:)
      real(8), allocatable, save :: qcx(:),qcy(:),qcz(:)
      real(8), allocatable, save :: ctdx(:),ctdy(:),ctdz(:)
      integer, allocatable, save :: potTypes(:),bndtyp(:,:)
      integer, allocatable, save :: angtyp(:,:)

      integer, save:: numPot, numPoints, totMol,nqcbnd,nqcang
      integer, save:: nrdfpts
      real(8), save:: rmax
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
      integer, dimension(1:20) :: fail

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
      allocate (ctdx(natms),stat=fail(8))
      allocate (ctdy(natms),stat=fail(9))
      allocate (ctdz(natms),stat=fail(10))
      allocate (bndtyp(mxtmls,ntbond),stat=fail(11))
      allocate (angtyp(mxtmls,ntangl),stat=fail(12))

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
     
      subroutine read_fqcmd_files(idnode,mxnode,ntbond,
     x  ntangl,ntpmls,lpimd)
c********************************************************
c
c     dl_poly quantum subroutine to read in f-QCMD 
c       details files
c
c     Authors: Nathan London 2024
c
c********************************************************

      implicit none 

      integer, intent(in) :: idnode,mxnode,ntbond,ntangl
      integer, intent(in) :: ntpmls
      logical, intent(in) :: lpimd

      character*1 :: molnam(40)
      character*8 keyword
      character*1 message(80)
      logical :: safe,loop1,loop2
      integer :: i,idum,itmols

      integer :: ibond,iang,bnd,bndt,ang,angt
      integer :: nbonds,nbt,nangls,nat
      real(8) engunit,parpot(6)
     
      write(*,*) "mxbond", mxbond 
      engunit =  418.4d0 
      data loop1/.true./,loop2/.true./
c     open force field data file
      
      if(idnode.eq.0)open (nfield,file='FQCMD',status='old')
      
      if(idnode.eq.0)
     x  write(nrite,"(/,/,'f-QCMD SPECIFICATION')")
      
      bndt=0
      angt=0
      bndtyp(:,:)=0
      call getrec(safe,idnode,nfield)
      
c     loop over molecules          
          do itmols=1,ntpmls
            
            bnd=1
            ang=1
            if(idnode.eq.0)
     x        write(nrite,"(/,1x,'molecular species type',9x,i10)")
     x        itmols
            
c     name of molecular species
            
            call getrec(safe,idnode,nfield)
            
            call copystring(record,molnam(1),40)
            if(idnode.eq.0)
     x        write(nrite,"(/,/,1x,'name of species:',13x,40a1)")
     x        (molnam(i),i=1,40)
            
            
c     read molecular data
            
            loop2=.true.
            
            do while(loop2)
              
              call getrec(safe,idnode,nfield)
c              if(.not.safe)call abort_field_read(1,idnode,nfield)
              
              call lowcase(record,lenrec)
              call strip(record,lenrec)
              
              
c     read chemical bond force constant and bondlength
                
              if(findstring('bonds',record,idum))then
                write(nrite,*) "bond"  
                nbt=intstr(record,lenrec,idum)
       
                write(nrite,*) "number of bond types: ",nbt 
                do ibond=1,nbt
                  bndt = bndt+1
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)return
              
                  call copystring(record,message,80)
                  call lowcase(record,4)
                  call getword(keyword,record,4,lenrec)

                  nbonds=intstr(record,lenrec,idum)
                  bndtyp(itmols,bnd:bnd+nbonds-1)=bndt
                  bnd=bnd+nbonds
                  write(nrite,*) "number of occurances: ",nbonds
                enddo
                
                
c     read intramolecular angular potential parameters
                
              elseif(findstring('angles',record,idum))then
                write(nrite,*) "angle"  
                nat=intstr(record,lenrec,idum)
       
                write(nrite,*) "number of angle types: ",nat 
                do ibond=1,nat

                  angt = angt+1
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)return
              
                  call copystring(record,message,80)
                  call lowcase(record,4)
                  call getword(keyword,record,4,lenrec)

                  nangls=intstr(record,lenrec,idum)
                  angtyp(itmols,ang:ang+nangls-1)=angt
                  ang=ang+nangls
                  write(nrite,*) "number of occurances: ",nangls 
                enddo
                
c     read intramolecular dihedral potential parameters
                
              elseif(findstring('dihedr',record,idum))then
                
c     finish of data for one molecular type
                
              elseif(findstring('finish',record,idum))then
                
c     running total of number of atoms in system
                
                loop2=.false.
                
              else
                
c     error exit for unidentified directive in molecular data
                
                if(idnode.eq.0)write(nrite,'(12x,a)')record
                call error(idnode,12)
                
              endif
              
            enddo
            
          
c     normal end of FIELD file
          
c        elseif(findstring('close',record,idum))then
          
          loop1=.false.
c          if(ntpvdw.eq.0.and.ntpmet.eq.0.and.
c     x      mod(keyfce,2).eq.1)call error(idnode,145)
          
c     error exit for unidentified directive
          
c        else
          
c          if(idnode.eq.0)write(nrite,'(100a)')record
c          call abort_field_read(2,idnode,nfield)
          
c        endif
        
      enddo

      end subroutine read_fqcmd_files

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
        if(i.lt.10)then
          write(filenm,'(i1)') i
        else
          write(filenm,'(i2)') i
        endif
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
c        write(*,*) potTables(i,:,numPoints)
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
     x  ntangl,ntpatm,ntpmls,lpimd)
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
      integer, intent(in) :: ntangl,ntpatm,ntpmls

      call allocate_fqcmd_arrays(idnode,mxnode,natms,ntbond,
     x  ntangl)

      call read_fqcmd_files(idnode,mxnode,ntbond,
     x  ntangl,ntpmls,lpimd)
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
     x  ntbond,ntangl,ntpmls,rcut,lfcmd)
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

      logical, intent(in) :: lfcmd
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
      real(8) :: bondVal2(ntbond), angVal2(ntangl)

c     allocate temporary arrays for atom distances      
      safe = .true.
      fail(:) = 0

      allocate(xdab(msbad),ydab(msbad),zdab(msbad),stat=fail(1))
      allocate(xdbc(msbad),ydbc(msbad),zdbc(msbad),stat=fail(2))
      if(any(fail.gt.0))safe=.false.
      if(mxnode.gt.1) call gstate(safe)
      if(.not.safe) call error(idnode,540)

c     calculate centroid values
      call get_centroid(natms)

c     get indices for parallel calculations      
      ibnd1 = (idnode*ntbond)/mxnode+1
      ibnd2 = ((idnode+1)*ntbond)/mxnode
      iang1 = (idnode*ntangl)/mxnode+1
      iang2 = ((idnode+1)*ntangl)/mxnode

      imol0 = (idnode*totmol)/mxnode+1

c     define histogram bin spacing      
      rdr = rmax/dble(nrdfpts)
      pi = 4.0*atan(1.0)
      rtheta = pi/dble(nrdfpts)

      qCent(:,:) = 0.d0
      bondVal(:) = 0.d0
      angVal(:) = 0.d0

c     Calculate bonds and angles from curvilinear coords      
      if(nbeads.gt.1) then
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

     
      endif 

c     Calculate bond and angle from centroid values      
      bondVal2(:) = 0.d0
      angVal2(:) = 0.d0
c     loop over beads to get bond distances      
       
        ii=0

        do i=ibnd1,ibnd2
          ii = ii + 1

          ia = listbnd(ii,2) 
          ib = listbnd(ii,3) 
           
          xdab(ii)=ctdx(ia)-ctdx(ib)
          ydab(ii)=ctdy(ia)-ctdy(ib)
          zdab(ii)=ctdz(ia)-ctdz(ib)
        enddo

        call images(imcon,0,1,ii,cell,xdab,ydab,zdab)
        
        ii=0

        do i=ibnd1,ibnd2
          ii = ii + 1

          rab = sqrt((xdab(ii)**2+ydab(ii)**2)+zdab(ii)**2)
c       store bond distances in temporary array
          bondVal2(i) = rab
        enddo
 
c       calculate angle values        
        ii=0

        do i=iang1,iang2
          ii = ii + 1

          ia = listang(ii,2)
          ib = listang(ii,3)
          ic = listang(ii,4) 

          xdab(ii)=ctdx(ia)-ctdx(ib)
          ydab(ii)=ctdy(ia)-ctdy(ib)
          zdab(ii)=ctdz(ia)-ctdz(ib)
          xdbc(ii)=ctdx(ic)-ctdx(ib)
          ydbc(ii)=ctdy(ic)-ctdy(ib)
          zdbc(ii)=ctdz(ic)-ctdz(ib)
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
          angVal2(i) = theta
        enddo

c     merge temporary arrays over processors      
      if(mxnode.gt.1) then
        call gdsum(bondVal,ntbond,buffer)
        call gdsum(angVal,ntangl,buffer)
        call gdsum(bondVal2,ntbond,buffer)
        call gdsum(angVal2,ntangl,buffer)
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
c           Use curvilinear coord definintion          
            if(i.eq.ntpmls.and.nbeads.gt.1.and.(.not.lfcmd))then
              qCent(imol,imolbnd) = bondVal(ibnd)
c           Use cartesian coord definintion          
            else
              qCent(imol,imolbnd) = bondVal2(ibnd)
            endif
            call qchistogram(bndtyp(i,imolbnd),qCent(imol,imolbnd),
     x        0.d0,rdr)
          enddo
          if(i.eq.ntpmls)then
        endif
          do k=1,numang(i)
            imolang = imolang + 1
            iang = iang +1
            if(i.eq.ntpmls.and.nbeads.gt.1.and.(.not.lfcmd))then
              qCent(imol,ntbond+imolang) = angVal(iang)
            else
              qCent(imol,ntbond+imolang) = angVal2(iang)
            endif
c     histogram internal coordinates
            call qchistogram(nqcbnd+angtyp(i,imolang),
     x        qCent(imol,ntbond+imolang),0.d0,rtheta)
          enddo
        enddo
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
     x  ntpmls,natms,ntbond,rcut,lpimd,lfcmd)
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

      logical, intent(in) :: lpimd,lfcmd
      integer, intent(in) :: idnode,mxnode,imcon,ntpmls
      integer, intent(in) :: natms,ntbond
      real(8), intent(in) :: rcut

      logical :: safe
      integer :: fail(2),ierr
      integer :: itmols,isite,imol,imol0,imol1,molacc 
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
      molacc = 0

c     define histrogram bin width
      rdr = rcut/dble(nrdfpts)

      qcx(:) = 0.d0
      qcy(:) = 0.d0
      qcz(:) = 0.d0
c       loop over molecule types        
      do itmols=1,ntpmls
c     perform eckart rotation for pimd with more than 1 bead
c       and only water in mixed system       
      if(itmols.gt.1) then
        jsite = jsite + nummols(itmols-1)*
     x      numsit(itmols-1)
        ksite = ksite + numsit(itmols-1)
        molacc = molacc + nummols(itmols-1)
      endif
        if(lpimd.and.(nbeads.gt.1).and.(.not.lfcmd).and.
     x      itmols.eq.ntpmls) then

c         get parallization indiced
          imol0 = (idnode*nummols(itmols))/mxnode+1
          imol1 = ((idnode+1)*nummols(itmols))/mxnode
          
c         get molecule mass for COM calc          
          molmass=0.d0
          do isite=1,numsit(itmols)
            mass(isite) = wgtsit(ksite+isite)
            molmass = molmass + mass(isite)
          enddo
c         copy over centroid
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
            copystart = (imol-1)*numsit(itmols)+jsite+1
            copyend = (imol-1)*numsit(itmols)+jsite+numsit(itmols)
            centx(1:numsit(itmols))=ctdx(copystart:copyend)
            centy(1:numsit(itmols))=ctdy(copystart:copyend)
            centz(1:numsit(itmols))=ctdz(copystart:copyend)
            call mol_gather(numsit(itmols),
     x            centx,centy,centz)
            do j=1,numsit(itmols)

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
                ztmp(j) = qCent(imol+molacc,1)
              else if (j.eq.3) then
                rad = qCent(imol+molacc,2)
                theta = qCent(imol+molacc,ntbond+1)
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
            qcx(copystart:copyend) = xtmp(1:numsit(itmols))
            qcy(copystart:copyend) = ytmp(1:numsit(itmols))
            qcz(copystart:copyend) = ztmp(1:numsit(itmols))
          enddo

c     for f-CMD use cartesian centooid          
        else if ((lfcmd.and.nbeads.gt.1).or.
     x    (lpimd.and.nbeads.gt.1.and.itmols.lt.ntpmls)) then
          imol0 = (idnode*nummols(itmols))/mxnode+1
          imol1 = ((idnode+1)*nummols(itmols))/mxnode
          
          do imol=imol0,imol1
            copystart = (imol-1)*numsit(itmols)+jsite+1
            copyend = (imol-1)*numsit(itmols)+jsite+numsit(itmols)
            qcx(copystart:copyend) = ctdx(copystart:copyend)
            qcy(copystart:copyend) = ctdy(copystart:copyend)
            qcz(copystart:copyend) = ctdz(copystart:copyend)
c            write(*,*) "qcx", qcx(copystart:copyend)
          enddo
c     for f-QCMD (not PIMD) simulations just use cartesian coords      
        else
          imol0 = (idnode*nummols(itmols))/mxnode+1
          imol1 = ((idnode+1)*nummols(itmols))/mxnode
          
          do imol=imol0,imol1
            copystart = (imol-1)*numsit(itmols)+jsite+1
            copyend = (imol-1)*numsit(itmols)+jsite+numsit(itmols)
            qcx(copystart:copyend) = xxx(copystart:copyend)
            qcy(copystart:copyend) = yyy(copystart:copyend)
            qcz(copystart:copyend) = zzz(copystart:copyend)
          enddo
        endif
      enddo
c       merge coords over nodes        
      if(mxnode.gt.1) then
        call gsync()
         call gdsum(qcx,natms,buffer)
         call gdsum(qcy,natms,buffer)
         call gdsum(qcz,natms,buffer)
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
            if(lpimd) then
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

      if(val.eq.0)then
        write(*,*) 'val is 0'
      endif
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
     x  nstrun,rcut,volm,temp,lpimd,ntpatm,wrtrdf)
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
      integer, intent(in) :: ntpatm,wrtrdf
      real(8), intent(in) :: rcut,volm,temp   

      integer :: i,j,k,ia,ib,ab
      real(8) :: ncount,pi,rdrbnd,rdrvdw,rtheta,rconv,r
      real(8) :: factor(2*numPot),econv

      rdrbnd = rmax/dble(nrdfpts)
      rdrvdw = rcut/dble(nrdfpts)
      pi = 4.0*atan(1.0)
      rtheta = pi/dble(nrdfpts)
      rconv = 0.529177210903
      econv = 262550.02  

c     merge RDFs over nodes      
      if(mxnode.gt.1) then
        do i=1,2*numPot+nqcbnd+nqcang
          call gdsum(qcrdf(i,:),nrdfpts,buffer)
        enddo
      endif

c     calculate normalization factors depending on RDF
      do ia=1,ntpatm
        do ib=ia,ntpatm
          ab=(ib*(ib-1))/2+ia
          k=potTypes(ab)
          if(k.ne.0) then
            factor(k) = volm*dens(ia)*dens(ib)*
     x        dble((nstrun-nsteql)/wrtrdf)*4.d0*pi
            if((ia.eq.ib).and.(volm*dens(ia).gt.1.d0))then
              factor(k) = factor(k)*0.5d0
            endif
            factor(numPot+k) = factor(k)*boltz*temp
          endif
        enddo
      enddo  

c     normalize intramolecular distributions to have unit
c       area      
      do i=1,nqcbnd
        call trapezoid_int(nrdfpts,rdrbnd,qcrdf(i,:),ncount)
        qcrdf(i,:) = qcrdf(i,:)/ncount
      enddo
      do i=1,nqcang
        call trapezoid_int(nrdfpts,rtheta,qcrdf(nqcbnd+i,:),ncount)
        qcrdf(nqcbnd+i,:) = qcrdf(nqcbnd+i,:)/ncount
      enddo
     
c     write out to files      
      if(idnode.eq.0) then
        open(48,file='Average_Intra.d',status='replace')
        do i=1,nrdfpts
          write(48,'(1e14.6)',advance='no') (dble(i-0.5)*
     x      rdrbnd)
          do j=1,nqcbnd
            write(48,'(1e14.6)',advance='no')qcrdf(j,i)
          enddo
          write(48,*) " "
        enddo
        close(48)
        
        open(48,file='Average_Angle.d',status='replace')
        do i=1,nrdfpts
          write(48,'(1e14.6)',advance='no') (dble(i-0.5)*
     x      rtheta)
          do j=1,nqcang
            write(48,'(1e14.6)',advance='no')
     x        qcrdf(nqcbnd+j,i)
          enddo
          write(48,*) " "
        enddo
        close(48)
        
        open(48,file='Average_RDF.d',status='replace')
c        if(lpimd.and.nbeads.gt.1)then
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
     x          qcrdf(nqcbnd+nqcang+j,i)/factor(j)
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
     
      subroutine get_centroid(natms)

      integer, intent(in) :: natms
      integer :: i

      do i=1,natms
        ctdx(i) = sum(xxx(i:natms*nbeads:natms))/dble(nbeads)
        ctdy(i) = sum(yyy(i:natms*nbeads:natms))/dble(nbeads)
        ctdz(i) = sum(zzz(i:natms*nbeads:natms))/dble(nbeads)
      enddo

      end subroutine get_centroid 
      end module fqcmd_module
