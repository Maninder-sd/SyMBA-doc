c converts binary file to ascii file without graze collisions

      include '../swift.inc'

      real*8 mass(NTPMAX)
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

      integer npart,nbodm,nbod,nbod0,ierr,ifol,istep
      integer iflgchk,iu,i,id,nleft,ium,k
      integer io_read_hdr,io_read_line,io_read_mass,io_read_mass_r
      integer io_read_hdr_r,io_read_line_r

      real*8 tloss,tadd,told
      integer iwhy,itype,inew,itest,iwhyold
      integer nlost,ndebris,ntotal,ninitial

      real*8 t0,tstop,dt,dtout,dtdump
      real*8 t,tnext

      real*8 tm
      real*8 rmin,rmax,rmaxu,qmin,rpl(NTPMAX),rhill(NTPMAX)
      logical*2 lclose
      real*8 a,e,inc,capom,omega,capm,j2rp2,j4rp4
      real*8 peri,apo,tcoll,obar
      real*8 dsurv,twall,tout

      integer isign, i2sign, i3sign, iLR, iSLR, itarg, iproj, ideb

      real*8 dr, coll(NTPMAX)

      integer iloss,iaccrete

      integer nexist, nbodn,n

      integer plist(NTPMAX),plist_inv(NTPMAX)
      integer slist(NTPMAX),slist_inv(NTPMAX)

c jen added
      real*8 mtotn(NTPMAX),mcoren(NTPMAX),mmantn(NTPMAX)
      real*8 LRMassn(NTPMAX),CMFi,CMFt,CMFp,CMFLR,CMFSLR
      real*8 x,y,z,vx,vy,vz,gm,gmsum,ap,ep,pinc,pcapom,pomega,pcapm
      integer ialpha
      integer ilossn(NTPMAX),iaccreten(NTPMAX),ntotaln(NTPMAX)

      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
      REAL(KIND=DP) :: mtot(NTPMAX),mcore(NTPMAX),mmant(NTPMAX)
      REAL(KIND=DP) :: munit,Fm,Mratio,mu,mearth,mdeb,mdcore,mdmant
      REAL(KIND=DP) :: tmass, pmass, LRMass, SLRMass, massi(NTPMAX)
      character*80 fdebris, fdiscard
      logical*2 dloss


c for subroutines
      real*8 tmtot,pmtot,tcore,pcore,tmant,pmant
      real*8 LRc,LRm,SLRc,SLRm


      character*80 outfile,inparfile,inplfile,fopenstat,inplfilen
      character*80 cname,c2,c3,c4,c5,c6,c7
      character*80 filenamef,filenamec,inpmfile

      dr = PI/180.0
      Fm = 7.51/7.58
c in kg
      mu=1.660538921e-27
      mearth=5.972e24
      Mratio=Fm*(55.845*mu)/((55.845*mu)+(24.305*mu)+(28.0855*mu)
     &   +(3.*15.999*mu))
      munit=5.0428958e28
      n=0

c     Read in param file for collision treatments
      write(*,*) 'Enter name of composition parameter file:'
      read(*,999) inpmfile
c inpmfile='ctypes.in'
      call io_init_cparam(inpmfile,cname,c2,c3,c4,c5,c6,c7,CMFi,tout,
     &  dloss,fdebris,fdiscard)


c     Prompt and read name of planet data file
c     write(*,*) ' '
c     write(*,*) 'Enter name of planet data file : '
c     read(*,*) inplfile
      inplfile='dump_pl.dat'
      call io_init_pl_symba(inplfile,lclose,iflgchk,nbod,npart,mass,
     &     xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

c     Get data for the run and the test particles
c      write(*,*) 'Enter name of parameter data file : '
c      read(*,999) inparfile
      inparfile='dump_param.dat'

      if (dloss) then
         call io_init_param_dloss(inparfile,t0,tstop,dt,dtout,dtdump,
     &    twall,dsurv,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,
     &         fopenstat)
      else
         call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,twall,
     &     iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,
     &     fopenstat)
      endif


      open(40,file=fdiscard,status='old',iostat=ierr)


      if (ierr.ne.0) then
         write(*,*) 'Could not open discard_mass.out'
         stop
      endif

      open(50,file=fdebris,status='old',iostat=ierr)
      if (ierr.ne.0) then
         write(*,*) 'Could not open add_debris.out 1'
         stop
      endif


      nlost=0
      ndebris=0
      do
         read(40,*,iostat=ierr) tloss,iwhy
         if (ierr.gt.0) then
            write(*,*) 'Could not read discard_mass.out'
            stop
         else if (ierr.lt.0) then
            exit
         else if (tloss.gt.t0) then
            write(*,*) 'reached last dump time'
            goto 3
         else
            if (iwhy.eq.2) then
c read in the target, projectile and total centre of mass values
               read(40,*) isign, itarg, tmass
               read(40,*)
               read(40,*)
               read(40,*) i2sign, iproj, pmass
               read(40,*)
               read(40,*)
               read(40,*)
               read(40,*)
               read(40,*)

               read(50,*,iostat=ierr) tadd,itype,inew
               if (ierr.gt.0) then
                  write(*,*) 'Could not read add_debris.out 2'
                  stop
               endif
               if (abs(tloss-tadd)/tloss.gt.1E-10) then
                  write(*,*) 'Time mismatch between loss and add'
                  write(*,*) tloss,tadd
                  stop
               endif

c     itype is the collision type
c     1 = Super-catastrophic disruption
c     2 = Catastrophic disruption
c     3 = Erosion
c     4 = Partial accretion
c     5 = Hit'n'spray
c     6 = Hit'n'run
c     7 = Bounce
c     8 = Perfect merger
c     9 = Graze'n'merge

               if (itype.eq.8.or.itype.eq.9) then
                  nlost=nlost+1
                  read(50,*)iLR,LRMass

               elseif (itype.eq.7) then
                  read(50,*)iLR,LRMass
                  read(50,*)iSLR,SLRMass

               elseif (itype.eq.6) then
                  read(50,*)iLR, LRMass
                  read(50,*)iSLR, SLRMass
                  ndebris=ndebris+inew

                  if (inew.gt.0) then
                     read(50,*)ideb,mdeb
                     do i=2,1+inew-1
                        read(50,*)
                     enddo
                  endif

               elseif (itype.eq.1) then
                  nlost=nlost+2
                  ndebris=ndebris+inew
                  if (inew.gt.0) then
                     read(50,*)ideb,mdeb
                     do i=2,inew
                        read(50,*)
                     enddo
                  endif

               else
                  nlost=nlost+1
                  ndebris=ndebris+inew
                  read(50,*) iLR,LRMass
                  if (inew.gt.0) then
                     read(50,*)ideb,mdeb
                     do i=2,inew
                        read(50,*)
                     enddo
                  endif
               endif


            else
               nlost=nlost+1
               do i=1,3
                  read(40,*)
               enddo

            endif
         endif
      enddo

 3    continue

      close(40)
      close(50)

 1000 format(a100,i6)
      ninitial=nbod+nlost-ndebris
      write(*,1000)"Number of initial bodies (embryos, planetesimals and
     & debris): ",ninitial
      write(*,1000)"Number of debris created: ",ndebris
      ntotal=nbod+nlost
      write(*,1000)"Number of total bodies: ",ntotal
      write(*,1000)"Number of bodies removed: ",nlost
      write(*,1000)"Number of surviving bodies and debris: ",nbod

      nexist = ninitial
c maintain this variable for after I read in planet file
      nbodn=nbod

c     Get data for the run and the test particles
c     write(*,*) 'Enter name of parameter data file : '
c     read(*,999) inparfile
c      inparfile='dump_param.dat'
c      call io_init_param(inparfile,t0,tstop,dt,dtout,dtdump,
c     &     iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,fopenstat)

c      open(unit=60,file='follow.all')

      open(unit=70,file='debris.origin')

c need to be able to specify/modify these names so you can do multiple runs with different parameters and not have them write over each other

c      write(filenamec,fmt='(a,a,a)')'pl.',cname,'compositionsf'
c      write(*,*)filenamec
      open(unit=80,file='pl.'//trim(cname)//'compositions-nograzefa')

c      write(filenamef,fmt='(a,a,a)')'follow.',cname,'collisionsf'
      open(unit=90,file='follow.'//trim(cname)//'collisions-nograzefa')

c     plist takes in initial index and returns current index or -1 if lost, debris start with an initial index
c     slist takes in initial index and returns current index or index of accreting body, debris start with an initial index
      do i=1,ntotal
            plist(i) = i
            slist(i) = i
      enddo

      do k=2,ntotal
         i=plist(k)
         plist_inv(i)=k
         i=slist(k)
         slist_inv(i)=k
      enddo

      if(btest(iflgchk,0)) then
         write(*,*) ' Reading an FXDR binary file '
      else if(btest(iflgchk,1)) then
         write(*,*) ' Reading an real*4 binary file '
      else
         write(*,*) ' ERROR: no binary file format specified '
         write(*,*) '        in param file '
         stop
      endif

      iu = 20
      ium = 30
      if(btest(iflgchk,0))  then ! bit 0 is set
         call io_open_fxdr(outfile, 'r', .true., iu, ierr)
         call io_open_fxdr('mass.'//outfile, 'r', .true.,
     &        ium, ierr)
      else
         call io_open(iu,outfile,'old','UNFORMATTED',ierr)
         call io_open(ium,'mass.'//outfile,'old','UNFORMATTED',ierr)
      endif

      open(40,file=fdiscard,status='old',iostat=ierr)
      if (ierr.ne.0) then
         write(*,*) 'Could not open discard_mass.out'
         stop
      endif

      open(50,file=fdebris,status='old',iostat=ierr)
      if (ierr.ne.0) then
         write(*,*) 'Could not open add_debris.out 3'
         stop
      endif

c     if tout < dtout then defaults to dtout
c     write(*,*)'output timestep'
c     read(*,*)tout
c      tout=1.0E5


      tnext=0.0

      read(40,*,iostat=ierr) tloss,iwhy
      if (ierr.gt.0) then
         write(*,*) 'Could not read discard_mass.out'
         stop
      else if (ierr.lt.0) then
         stop
      endif

c need to read in initial body masses from planet file
c Prompt and read name of planet data file
      write(*,*)''
      write(*,*)'Enter name of planet data file:'
      read(*,999) inplfilen
 999  format(a)
      call io_init_pl_symba(inplfilen,lclose,iflgchk,nbod,npart,mass,
     &     xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

c reset nbod
      massi = mass
      write(*,*)'at t =0 nexist is, nbod',nexist,nbod
      if(nexist.ne.nbod) then
         stop
      endif
      do i=nbod+1,ntotal
         massi(i)=-99.
      enddo
      mcore = massi*CMFi
      mmant = massi-mcore
      mtot = massi
      nbod=nbodn

 1    continue


         if(btest(iflgchk,0))  then ! bit 0 is set
            ierr = io_read_hdr(iu,t,nbod,nleft)
            if(ierr.ne.0) goto 2
            ierr = io_read_mass(tm,nbodm,mass,ium)
            if(ierr.ne.0) goto 2
         else
            ierr = io_read_hdr_r(iu,t,nbod,nleft)
            if(ierr.ne.0) goto 2
            ierr = io_read_mass_r(tm,nbodm,mass,ium)
            if(ierr.ne.0) goto 2
         endif

         if(tm.ne.t) then
            write(*,*) ' Mismatch between time:',t,tm
            goto 2
         endif

         if(nbodm.ne.nbod) then
            write(*,*) ' Mismatch between number of bodies at time tm:'
     &      ,tm,nbod,nbodm
            goto 2
         endif
c use this loop to calculate composition

      do while (tloss.le.t)

         if (iwhy.eq.2) then
            read(40,*)i,iaccrete,tmass
            read(40,*)x,y,z
            read(40,*)vx,vy,vz
c            do i=1,2
c                read(40,*)x,y,z
c            enddo
            read(40,*)i,iloss,pmass
            do i=1,5
               read(40,*)
            enddo

c calculate semi-major axis of collision
         gm=TWOPI*TWOPI
         gmsum=gm+tmass
         call orbel_xv2el(x,y,z,vx,vy,vz,gmsum,
     &     ialpha,ap,ep,pinc,pcapom,pomega,pcapm)
c         rinc=rinc/pi*180.
c         capom=capom/pi*180.
c         omega=omega/pi*180.
c         capm=capm/pi*180.
c         write(13,*)gm/gm0*328000, a,e,rinc,capom,omega,capm

            read(50,*,iostat=ierr) tadd,itype,inew
            if (ierr.gt.0) then
               write(*,*) 'Could not read add_debris.out 4'
               stop
            endif
            if (abs(tloss-tadd)/tloss.gt.1E-10) then
               write(*,*) 'Time mismatch between loss and add:'
               write(*,*) tloss,tadd
               stop
            endif
            if(inew.gt.0) then
               do i=1,inew
                  nexist=nexist+1
                  write(70,4000)tloss,plist_inv(nexist),
     &            plist_inv(iaccrete)
               enddo
            endif
c     itype is the collision type
c     1 = Super-catastrophic disruption
c     2 = Catastrophic disruption
c     3 = Erosion
c     4 = Partial accretion
c     5 = Hit'n'spray
c     6 = Hit'n'run
c     7 = Bounce
c     8 = Perfect merger
c     9 = Graze'n'merge
c            write(*,*)ntotal
            if (itype.eq.8.or.itype.eq.9) then
               read(50,*)iLR,LRMass !read in LRMass

c              check that you're not calling a mass that's been destroyed
               if(mtot(iaccrete).eq.-99.)then
                  write(*,*)'mtot = -99 for',iaccrete,tadd
               endif
               if(mtot(iloss).eq.-99.)then
                  write(*,*)'mtot = -99 for',iloss,tadd
               endif
c              check that core + mantle = total mass
c               if((mcore(iaccrete)+mmant(iaccrete)).gt.
c     &         tmass+1e-10)then
c                  write(*,*)'problem with target mass',mtot(iaccrete)
c                  write(*,*)tmass,mcore(iaccrete),mmant(iaccrete)
c               endif

c              calculate compositions
               call merge(mtot(iaccrete),mcore(iaccrete),
     & mmant(iaccrete),mtot(iloss),mcore(iloss),mmant(iloss),
     & LRMass,LRc,LRm)


               CMFt = mcore(iaccrete)/mtot(iaccrete)
               CMFp = mcore(iloss)/mtot(iloss)
               CMFLR = LRc/(LRc+LRm)

               write(90,5000)tloss,ap,iaccrete,plist_inv(iaccrete),
     &   tmass,CMFt,iloss,plist_inv(iloss),pmass,CMFp,itype,
     &   iLR,LRMass,CMFLR,0,0e0,0.,0,0,0e0

c              set new values for mass
               mcore(iLR)=LRc
               mmant(iLR)=LRm
               mtot(iLR)=mcore(iLR)+mmant(iLR)

c              check calculations are correct
               if(mtot(iLR).gt.(LRMass+1e-10))then
                  write(*,*)'mtot greater than LRMass at',tadd
                  write(*,*) mcore(iLR),mmant(iLR),mtot(iLR)
                  write(*,*) pmass,mcore(iloss),mmant(iloss),
     &         tmass,mcore(iaccrete),mmant(iaccrete),LRMass
               endif
               if(mmant(iLR).lt.0)then
                  write(*,*)'mant=negative for',iLR,itype
               endif

               told=tloss
               read(40,*,iostat=ierr) tloss,iwhy
               if (ierr.gt.0) then
                  write(*,*) 'Could not read discard_mass.out'
                  stop
               else if (ierr.lt.0) then
                  exit
               endif
c              check that iLR = iaccrete
               if(iLR.ne.iaccrete) then
                  write(*,*)'iacc =\= iLR',iLR,iaccrete
               endif

               if (tloss.eq.told) then
c              keep track of values and go to next loop
                  n=n+1
                  ilossn(n)=iloss
                  ntotaln(n)=ntotal
                  mcoren(n)=LRc
                  mmantn(n)=LRm
                  mtotn(n)=LRc+LRm
                  iaccreten(n)=iaccrete
                  LRMassn(n)=LRMass
                  goto 4

               elseif(tloss.ne.told.and.n.ne.0) then

c                 call plist the correct number of times
                  do i=1,n
                     nexist=nexist-1
                     call plist_loss(ilossn(n),ntotaln(n),plist,
     &               mcore,mmant,mtot,massi,nexist)
                     call slist_collision(ilossn(n),iaccreten(n),
     &               ntotaln(n),slist)
                  enddo
                  write(*,*)'tloss eq told finished for n, t=',n,told
                  n=0

               endif

c              remove the lost body 
               nexist=nexist-1
               call plist_loss(iloss,ntotal,plist,mcore,mmant,
     &         mtot,massi,nexist)
               call slist_collision(iloss,iaccrete,ntotal,slist)

            elseif (itype.eq.7) then
c              bounce collision
               read(50,*)iLR,LRMass
               read(50,*)iSLR,SLRMass

c              check that you're not calling a mass that's been destroyed
               if(mtot(iaccrete).eq.-99.)then
                  write(*,*)'mtot = -99 for',iaccrete,tadd,mtot
               endif
               if(mtot(iloss).eq.-99.)then
                  write(*,*)'mtot = -99 for',iloss,tadd,mtot
               endif
   
c              calculate compositions
               call bounce(c7,mtot(iaccrete),mcore(iaccrete),
     &  mmant(iaccrete),mtot(iloss),mcore(iloss),mmant(iloss),
     & LRMass,SLRMass,LRc,LRm,SLRc,SLRm)

               CMFt = mcore(iaccrete)/mtot(iaccrete)
               CMFp = mcore(iloss)/mtot(iloss)
               CMFLR = LRc/(LRc+LRm)
               CMFSLR = SLRc/(SLRc+SLRm)

c              write to follow.collisions
               write(90,5000)tloss,ap,iaccrete,plist_inv(iaccrete),
     &         tmass,CMFt,iloss,plist_inv(iloss),pmass,CMFp,itype,
     &         iLR,LRMass,CMFLR,iSLR,SLRMass,CMFSLR,0,0,0e0

c              assign new values
               mcore(iLR)=LRc
               mmant(iLR)=LRm
               mtot(iLR)=mcore(iLR)+mmant(iLR)
               mcore(iSLR)=SLRc
               mmant(iSLR)=SLRm
               mtot(iSLR)=mcore(iSLR)+mmant(iSLR)

c              check that values are positive
               if(mmant(iLR).lt.0)then
                  write(*,*)'mant=negative for',iLR,itype
               endif
               if(mmant(iSLR).lt.0)then
                  write(*,*)'mant=negative for',iSLR,itype
               endif
               if(mtot(iLR).le.1e-11)then
                  write(*,*)'mtot=0 for',itype
               endif
               if(mtot(iSLR).le.1e-11)then
                  write(*,*)'mtot=0 for',itype
               endif
               if(mtot(iLR).gt.(LRMass+1e-10)) then
                  write(*,*)'mtot greater than LRMass for colltype7'
               endif
               if(mtot(iSLR).gt.(SLRMass+1e-10))then
                  write(*,*)'mtot greater than SLRMass for colltype7'
               endif

               told=tloss
               read(40,*,iostat=ierr) tloss,iwhy
               if (ierr.gt.0) then
                  write(*,*) 'Could not read discard_mass.out'
                  stop
               else if (ierr.lt.0) then
                  exit
               endif

            elseif (itype.eq.6) then
c              hit and run
               read(50,*)iLR,LRMass
               read(50,*)iSLR,SLRMass

c              read in ideb and mdeb
               if (inew.gt.0) then
                  read(50,*)ideb,mdeb
               elseif (inew.le.0) then
                  ideb = 0
                  mdeb = 0.
               endif

c              check that you're not calling mass that's been destroyed
               if(mtot(iaccrete).eq.-99.)then
                  write(*,*)'mtot = -99 for',iaccrete,tadd,mtot
               endif
               if(mtot(iloss).eq.-99.)then
                  write(*,*)'mtot = -99 for',iloss,tadd,mtot
               endif

c              calculate compositions
               call hitnrun(c6,mtot(iaccrete),mcore(iaccrete),
     & mmant(iaccrete),mtot(iloss),mcore(iloss),mmant(iloss),
     & LRMass,SLRMass,LRc,LRm,SLRc,SLRm,mdcore,mdmant)


               CMFt = mcore(iaccrete)/mtot(iaccrete)
               CMFp = mcore(iloss)/mtot(iloss)
               CMFLR = LRc/(LRc+LRm)
               CMFSLR = SLRc/(SLRc+SLRm)

c              write to follow.collisions
               write(90,5000)tloss,ap,iaccrete,plist_inv(iaccrete),
     &         tmass,CMFt,iloss,plist_inv(iloss),pmass,CMFp,itype,iLR,
     &         LRMass,CMFLR,iSLR,SLRMass,CMFSLR,inew-1,ideb,mdeb

c              set new values
               mcore(iLR)=LRc
               mcore(iSLR)=SLRc
               mmant(iLR)=LRm
               mmant(iSLR)=SLRm
               mtot(iLR)=mcore(iLR)+mmant(iLR)
               mtot(iSLR)=mcore(iSLR)+mmant(iSLR)

c              check that values are positive
               if(mmant(iLR).lt.0)then
                  write(*,*)'mant=negative for',iLR,itype
               endif
               if(mmant(iSLR).lt.0)then
                  write(*,*)'mant=negative for',iSLR,itype
               endif
               if(mtot(iLR).le.1e-12) then
                  write(*,*)'mtot=0 for',itype
               endif
               if(mtot(iSLR).le.1e-12) then
                  write(*,*)'mtot=0 for',itype
               endif
               if(mtot(iLR).gt.(LRMass+1e-10))then
                  write(*,*)'mtot greater than LRMass for colltype6'
               endif
               if(mtot(iSLR).gt.(SLRMass+1e-10))then
                  write(*,*)'mtot greater than SLRMass for colltype6'
               endif

c              set debris values
               mcore(ideb)=mdcore*mdeb
               mmant(ideb)=mdmant*mdeb
               mtot(ideb)=mdeb
               if (inew.gt.0) then
                  do i=2,inew
                     read(50,*)ideb,mdeb
                     mcore(ideb)=mdcore*mdeb
                     mmant(ideb)=mdmant*mdeb
                     mtot(ideb)=mdeb

c                    check debris values are valid
                     if(mtot(ideb).le.0)then
                        write(*,*)'debris mtot=0 for',itype
                     endif

                  enddo
               endif

               told=tloss
               read(40,*,iostat=ierr) tloss,iwhy
               if (ierr.gt.0) then
                  write(*,*) 'Could not read discard_mass.out'
                  stop
               else if (ierr.lt.0) then
                  exit
               endif

            elseif (itype.ge.2.and.itype.le.3) then
c              disruption, erosion and partial accretion
               read(50,*)iLR,LRmass

               if (inew.gt.0) then
                  read(50,*)ideb,mdeb
               elseif (inew.le.0) then
                  ideb = 0
                  mdeb = 0.
               endif

c              check you're not reading in particles that have been delected
               if(mtot(iaccrete).eq.-99.)then
                  write(*,*)'mtot = -99 for',iaccrete,tadd,mtot
               endif
               if(mtot(iloss).eq.-99)then
                  write(*,*)'mtot = -99 for',iloss,tadd,mtot
               endif

c              calculate composition
               call erosion(c3,mtot(iaccrete),mcore(iaccrete),
     & mmant(iaccrete),mtot(iloss),mcore(iloss),mmant(iloss),
     & LRMass,LRc,LRm,mdcore,mdmant)

               CMFt = mcore(iaccrete)/mtot(iaccrete)
               CMFp = mcore(iloss)/mtot(iloss)
               CMFLR = LRc/(LRc+LRm)

c              write to follow.collisions
               write(90,5000)tloss,ap,iaccrete,plist_inv(iaccrete),
     &         tmass,CMFt,iloss,plist_inv(iloss),pmass,CMFp,itype,iLR,
     &         LRMass,CMFLR,0,0e0,0.,inew,ideb,mdeb

c              set new values
               mcore(iLR)=LRc
               mmant(iLR)=LRm
               mtot(iLR)=mcore(iLR)+mmant(iLR)
c               if(mdcore.lt.1e-10)then
c                     write(*,*)'-core',mcore(iLR),LRMass,mdcore,mdmant
c                     write(*,*)mtot(iloss),mtot(iaccrete),
c     &            mcore(iaccrete)

c              check that values make sense
               if(mdmant.lt.0)then
                     write(*,*)'-mantle',mcore(iLR),LRMass,mdmant,mdcore
                     write(*,*)mtot(iloss),mtot(iaccrete),
     &            mcore(iaccrete)
               endif
               if(mmant(iLR).lt.0)then
                  write(*,*)'mant=negative for',iLR,itype
               endif
               if(mtot(iLR).le.1e-11)then
                  write(*,*)'mtot=0 for',itype
               endif
               if(mtot(iLR).gt.(LRMass+1e-10))then
                  write(*,*)'mtot greater than LRMass for colltype2-4'
               endif

c              set new debris values
               mcore(ideb)=mdcore*mdeb
               mmant(ideb)=mdmant*mdeb
               mtot(ideb)=mdeb
               if (inew.gt.0) then
                  do i=2,inew
                     read(50,*)ideb,mdeb
                     mcore(ideb)=mdeb*mdcore
                     mmant(ideb)=mdeb*mdmant
                     mtot(ideb)=mdeb

c                    check values are valid
                     if(mtot(ideb).le.0)then
                        write(*,*)'debris mtot=0 for',itype
                     endif
c                  if(mmant(ideb).lt.0)then
c                     write(*,*)'mmant -ve for',ideb,itype
c                  endif
                  enddo
               endif

               nexist=nexist-1
               call plist_loss(iloss,ntotal,plist,mcore,mmant,
     &         mtot,massi,nexist)
               call slist_collision(iloss,iaccrete,ntotal,slist)

               told=tloss
               read(40,*,iostat=ierr) tloss,iwhy
               if (ierr.gt.0) then
                  write(*,*) 'Could not read discard_mass.out'
                  stop
               else if (ierr.lt.0) then
                  exit
               endif

            elseif (itype.eq.4) then
c              partial accretion
               read(50,*)iLR,LRmass
               if(inew.gt.0) then
                  read(50,*)ideb,mdeb
               elseif(inew.le.0) then
                  ideb = 0
                  mdeb = 0.
               endif

c              check you're not reading in particles that have been delected
               if(mtot(iaccrete).eq.-99.)then
                  write(*,*)'mtot = -99 for',iaccrete,tadd,mtot
               endif
               if(mtot(iloss).eq.-99)then
                  write(*,*)'mtot = -99 for',iloss,tadd,mtot
               endif

c              calculate composition
               call partacc(c4,mtot(iaccrete),mcore(iaccrete),
     & mmant(iaccrete),mtot(iloss),mcore(iloss),mmant(iloss),
     & LRMass,LRc,LRm,mdcore,mdmant)

               CMFt = mcore(iaccrete)/mtot(iaccrete)
               CMFp = mcore(iloss)/mtot(iloss)
               CMFLR = LRc/(LRc+LRm)

c              write to follow.collisions
               write(90,5000)tloss,ap,iaccrete,plist_inv(iaccrete),
     &         tmass,CMFt,iloss,plist_inv(iloss),pmass,CMFp,itype,iLR,
     &         LRMass,CMFLR,0,0e0,0.,inew-1,ideb,mdeb

c              set new values
               mcore(iLR)=LRc
               mmant(iLR)=LRm
               mtot(iLR)=mcore(iLR)+mmant(iLR)
c               if(mdcore.lt.1e-10)then
c                     write(*,*)'-core',mcore(iLR),LRMass,mdcore,mdmant
c                     write(*,*)mtot(iloss),mtot(iaccrete),
c     &            mcore(iaccrete)

c              check that values make sense
               if(mdmant.lt.0)then
                     write(*,*)'-mantle',mcore(iLR),LRMass,mdmant,mdcore
                     write(*,*)mtot(iloss),mtot(iaccrete),
     &            mcore(iaccrete)
               endif
               if(mmant(iLR).lt.0)then
                  write(*,*)'mant=negative for',iLR,itype
               endif
               if(mtot(iLR).le.1e-12)then
                  write(*,*)'mtot=0 for',itype
               endif
               if(mtot(iLR).gt.(LRMass+1e-10))then
                  write(*,*)'mtot greater than LRMass for colltype2-4'
               endif

c              set new debris values
               mcore(ideb)=mdcore*mdeb
               mmant(ideb)=mdmant*mdeb
               mtot(ideb)=mdeb
               if (inew.gt.0) then
                  do i=2,inew
                     read(50,*)ideb,mdeb
                     mcore(ideb)=mdeb*mdcore
                     mmant(ideb)=mdeb*mdmant
                     mtot(ideb)=mdeb

c                    check values are valid
                     if(mtot(ideb).le.0)then
                        write(*,*)'debris mtot=0 for',itype
                     endif
c                  if(mmant(ideb).lt.0)then
c                     write(*,*)'mmant -ve for',ideb,itype
c                  endif
                  enddo
               endif

               nexist=nexist-1
               call plist_loss(iloss,ntotal,plist,mcore,mmant,
     &         mtot,massi,nexist)
               call slist_collision(iloss,iaccrete,ntotal,slist)

               told=tloss
               read(40,*,iostat=ierr) tloss,iwhy
               if (ierr.gt.0) then
                  write(*,*) 'Could not read discard_mass.out'
                  stop
               else if (ierr.lt.0) then
                  exit
               endif

            elseif(itype.eq.5) then
c              hit'n'spray
               read(50,*)iLR,LRMass
               if (inew.gt.0) then
                  read(50,*)ideb,mdeb
               elseif (inew.le.0) then
                  ideb = 0
                  mdeb = 0.
               endif

c              check that values are valid
               if(mtot(iaccrete).eq.-99.)then
                  write(*,*)'mtot = -99 for',iaccrete,tadd,mtot
               endif
               if(mtot(iloss).eq.-99.)then
                  write(*,*)'mtot = -99 for',iloss,tadd,mtot
               endif

c              calculate composition
               call hitnspray(c5,mtot(iaccrete),mcore(iaccrete),
     & mmant(iaccrete),mtot(iloss),mcore(iloss),mmant(iloss),
     & LRMass,LRc,LRm,mdcore,mdmant)

               CMFt = mcore(iaccrete)/mtot(iaccrete)
               CMFp = mcore(iloss)/mtot(iloss)
               CMFLR = LRc/(LRc+LRm)

c              write to follow.collisions
               write(90,5000)tloss,ap,iaccrete,plist_inv(iaccrete),
     &         tmass,CMFt,iloss,plist_inv(iloss),pmass,CMFp,itype,iLR,
     &         LRMass,CMFLR,0,0e0,0.,inew,ideb,mdeb

c              set new values
               mcore(iLR)=LRc
               mmant(iLR)=LRm
               mtot(iLR)=mcore(iLR)+mmant(iLR)

c              check values are valid
               if(mmant(iLR).lt.0)then
                  write(*,*)'mant=negative for',iLR,itype
               endif
               if(mtot(iLR).le.1e-11)then
                  write(*,*)'mtot=0 for',itype
               endif
               if(mtot(iLR).gt.(LRMass+1e-10))then
                  write(*,*)'mtot greater than LRMass for colltype5'
               endif

c              set new debris values
               if (inew.gt.0) then
                  mcore(ideb)=mdcore*mdeb
                  mmant(ideb)=mdmant*mdeb
                  mtot(ideb)=mdeb
                  do i=3,1+inew
                     read(50,*)ideb,mdeb
                     mcore(ideb)=mdcore*mdeb
                     mmant(ideb)=mdmant*mdeb
                     mtot(ideb)=mdeb

c                    check debris values are valid
                     if(mtot(ideb).le.0)then
                        write(*,*)'debris mtot=0 for',ideb,itype
                     endif
                     if(mmant(ideb).lt.0)then
                        write(*,*)'mmant -ve for',ideb,itype
                     endif
                  enddo
               endif

               nexist=nexist-1
               call plist_loss(iloss,ntotal,plist,mcore,mmant,
     &         mtot,massi,nexist)
               call slist_collision(iloss,iaccrete,ntotal,slist)

               told=tloss
               read(40,*,iostat=ierr) tloss,iwhy
               if (ierr.gt.0) then
                  write(*,*) 'Could not read discard_mass.out'
                  stop
               else if (ierr.lt.0) then
                  exit
               endif

            elseif (itype.eq.1) then
c              super-catastrophic
               mdcore=(mcore(iaccrete)+mcore(iloss))
     &    /(mtot(iaccrete)+mtot(iloss))
               mdmant=1.-mdcore

               if(inew.gt.0) then
                  read(50,*) ideb,mdeb
c                 set new values
                  mcore(ideb) = mdcore*mdeb
                  mmant(ideb)=mdmant*mdeb
                  mtot(ideb)=mdeb
               elseif(inew.le.0) then
                  ideb = 0
                  mdeb = 0.
               endif

               CMFt = mcore(iaccrete)/mtot(iaccrete)
               CMFp = mcore(iloss)/mtot(iloss)

c              write to follow.collisions
               write(90,5000)tloss,ap,iaccrete,plist_inv(iaccrete),
     &  tmass,CMFt,iloss,plist_inv(iloss),pmass,CMFp,
     &  itype,0,0e0,0.,0,0e0,0.,inew,ideb,mdeb

c              set new values
               if (inew.gt.0) then
                  do i=2,inew
                     read(50,*)ideb,mdeb
                     mcore(ideb)=mdcore*mdeb
                     mmant(ideb)=mdmant*mdeb
                     mtot(ideb)=mdeb

c                    check that values are valid
                     if(mtot(ideb).le.0)then
                        write(*,*)'debris mtot=0 for',itype
                     endif
                     if(mmant(ideb).lt.0)then
                        write(*,*)'mmant -ve for',ideb,itype
                     endif
                  enddo
               endif

               nexist=nexist-2
c               if(nexist.le.280) then
c                  write(*,*) 'nexist in coll 1 =',nexist,tloss
c               endif
               if(iloss.gt.iaccrete) then
                  call plist_loss(iloss,ntotal,plist,mcore,mmant,
     &            mtot,massi,(nexist+1))
c                  write(*,*) 'nexist in coll 1=',nexist+1,tloss
                  call plist_loss(iaccrete,ntotal,plist,mcore,mmant,
     &            mtot,massi,nexist)
c                  write(*,*) 'nexist in coll 1=',nexist,tloss
                  call slist_loss(iloss,ntotal,slist)
                  call slist_loss(iaccrete,ntotal,slist)
               else
                  call plist_loss(iaccrete,ntotal,plist,mcore,mmant,
     &            mtot,massi,(nexist+1))
c                  write(*,*) 'nexist in coll 1=',nexist+1,tloss
                  call plist_loss(iloss,ntotal,plist,mcore,mmant,
     &            mtot,massi,nexist)
c                  write(*,*) 'nexist in coll 1=',nexist,tloss
                  call slist_loss(iaccrete,ntotal,slist)
                  call slist_loss(iloss,ntotal,slist)
               endif

               told=tloss
               read(40,*,iostat=ierr) tloss,iwhy
               if (ierr.gt.0) then
                  write(*,*) 'Could not read discard_mass.out'
                  stop
               else if (ierr.lt.0) then
                  exit
               endif

            else
c              leaving this just in case I messed something up
               read(50,*)ideb,mdeb
               write(*,*)'no collision type'
               do i=2,1+inew
                  read(50,*)
               enddo
               nexist=nexist-1
               call plist_loss(iloss,ntotal,plist,mcore,mmant,
     &         mtot,massi,nexist)
               call slist_collision(iloss,iaccrete,ntotal,slist)

               told=tloss
               read(40,*,iostat=ierr) tloss,iwhy
               if (ierr.gt.0) then
                  write(*,*) 'Could not read discard_mass.out'
                  stop
               else if (ierr.lt.0) then
                  exit
               endif

            endif
         else
            read(40,*) i,iloss
            if (ierr.gt.0) then
               write(*,*) 'Could not read discard_mass.out',ierr,tloss
               stop
            else if (ierr.lt.0) then
               exit
            endif
            nexist=nexist-1
            call plist_loss(iloss,ntotal,plist,mcore,mmant,
     &      mtot,massi,nexist)
            call slist_loss(iloss,ntotal,slist)
            do i=1,2
               read(40,*)
            enddo

            told=tloss
            iwhyold=iwhy
            read(40,*,iostat=ierr) tloss,iwhy
            if (ierr.gt.0) then
               write(*,*) 'Could not read discard_mass.out'
               stop
            else if (ierr.lt.0) then
               exit
            endif

         endif

c        check that you have the right number of particles
         if(mtot(nexist).eq.-99.) then
            write(*,*)'nexist = ',nexist
            write(*,*)'mtot(nexist-1)',mtot(nexist-1)
            write(*,*)'time',tloss
            goto 2
         endif

c        leave this here to double check that nothing is happening with
c        the other collision types
         if (told.eq.tloss.and.iwhy.eq.2) then
            if (itype.eq.1) then
               write(*,*) 'fuck',itype,tloss

            elseif(itype.eq.8.or.itype.eq.9) then
               write(*,*) 'type 8 or 9',tloss,iwhyold

            elseif(itype.ge.2.and.itype.le.5) then
               write(*,*) 'fuck',itype,tloss

            elseif(itype.eq.6.or.itype.eq.7) then
               write(*,*) 'fuck',itype,tloss
            endif
         endif

         do k=2,ntotal
            i=plist(k)
            plist_inv(i)=k
            i=slist(k)
            slist_inv(i)=k
         enddo

 4       continue

      enddo

         do i=2,nbod
            if(btest(iflgchk,0))  then ! bit 0 is set
               ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm)
            else
               ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm)
            endif
            if(ierr.ne.0) goto 2

            if(t.ge.tnext)then
c               write(60,3000)t,plist_inv(i),a,e,inc/dr,
c     &      mass(i)*8442.69
               write(80,6000)t,plist_inv(i),a,e,inc/dr,
     &  mass(i)*munit/mearth,i,mcore(i)*munit/mearth,
     &  mmant(i)*munit/mearth,mtot(i)*munit/mearth
            endif
         enddo
 2000    format(es10.4,1x,9(i6,1x))
 3000    format(es10.4,1x,i6,1x,f8.3,1x,f8.4,1x,f8.3,1x,es10.3)
 4000    format(es10.4,1x,i6,1x,i6)
 5000    format(es13.7,1x,f8.3,1x,2(i6,1x,i6,1x,es10.3,1x,f8.4,1x),
     &   i6,1x,2(i6,1x,es10.3,1x,f8.4,1x),i6,1x,i6,1x,es10.3)
 7000    format(es10.4,1x,i6,1x,es10.3,1x,i6,1x,es10.3,1x,
     &   i6,1x,i6,1x,es10.3,1x,i6,1x,es10.3,1x,i6,1x,i6,1x,es10.3)
 6000    format(es13.7,1x,i6,1x,f8.3,1x,f8.4,1x,f8.3,1x,es10.3,1x,i6,1x,
     &        es10.3,1x,es10.3,1x,es10.3,1x)

         if(tout.lt.1e3.and.t.ge.1e5) tout=1e3
         if (t.ge.tnext) tnext=tnext+tout

         goto 1

 2    continue

      close(20)
      close(30)
      close(40)
      close(50)
c      close(60)
      close(70)
      close(80)
      close(90)

      stop
      end

c-------------------------------------------
      subroutine plist_loss(iloss,ntotal,plist,mcore,mmant,
     & mtot,massi,nexist)

      include 'swift.inc'

      integer i,iloss,ntotal,plist(ntotal),nexist

      real*8 mcore(NTPMAX),mmant(NTPMAX),mtot(NTPMAX)
      real*8 massi(NTPMAX)

      do i=1,ntotal
         if(plist(i).eq.iloss) then
               plist(i) = -1
         endif
      enddo

      do i=1,ntotal
         if(plist(i).gt.iloss) then
            plist(i) = plist(i) - 1
         endif
      enddo

      do i=1,nexist
         if(i.ge.iloss) then
            mcore(i) = mcore(i+1)
            mmant(i) = mmant(i+1)
            mtot(i) = mtot(i+1)
            massi(i) = massi(i+1)
         endif
         if(mtot(i).eq.-99.) then
            write(*,*) 'mtot=-99',i,nexist,iloss
         endif
      enddo
      
      mcore(nexist+1) = -99.
      mmant(nexist+1) = -99.
      mtot(nexist+1) = -99.
      massi(nexist+1) = -99.

      if (mtot(nexist+2).ne.-99.) then
         write(*,*) 'mtot(nexist+2)=',mtot(nexist+2),nexist
      endif

      return
      end

c-------------------------------------------
      subroutine slist_loss(iloss,ntotal,slist)

      include 'swift.inc'

      integer i,iloss,ntotal,slist(ntotal)

      do i=1,ntotal
         if(slist(i).eq.iloss) then
               slist(i) = -1
         endif
      enddo

      do i=1,ntotal
         if(slist(i).gt.iloss) then
            slist(i) = slist(i) - 1
         endif
      enddo

      return
      end

c-------------------------------------------
      subroutine slist_collision(iloss,iaccrete,ntotal,slist)

      include 'swift.inc'

      integer i,iloss,iaccrete,ntotal,slist(ntotal)

      do i=1,ntotal
         if(slist(i).eq.iloss) then
               slist(i) = iaccrete
         endif
      enddo

      do i=1,ntotal
         if(slist(i).gt.iloss) then
            slist(i) = slist(i) - 1
         endif
      enddo

      return
      end

c-------------------------------------------
      subroutine tlist_loss(iloss,ntotal,tlist)

      include 'swift.inc'

      integer i,iloss,ntotal,tlist(ntotal)

      do i=1,ntotal
         if(tlist(i).eq.iloss) then
               tlist(i) = -1
         endif
      enddo

      do i=1,ntotal
         if(tlist(i).gt.iloss) then
            tlist(i) = tlist(i) - 1
         endif
      enddo

      return
      end

c-------------------------------------------
      subroutine merge(tmtot,tcore,tmant,pmtot,pcore,pmant,
     & LRMass,LRc,LRm)

      include 'swift.inc'

      real*8 tmtot,pmtot,tcore,pcore,tmant,pmant
      real*8 LRc,LRm
      real*8 mdcore,mdmant,LRMass

c merges cores, debris mass (if any) is mantle
      LRc=tcore+pcore
      LRm=LRMass-LRc
      return
      end

c-------------------------------------------
      subroutine bounce(c7,tmtot,tcore,tmant,pmtot,pcore,pmant,
     & LRMass,SLRMass,LRc,LRm,SLRc,SLRm)

      include 'swift.inc'

      real*8 tmtot,pmtot,tcore,pcore,tmant,pmant
      real*8 LRc,LRm,SLRc,SLRm
      real*8 mdcore,mdmant,LRMass,SLRMass,ba
c parameters
      character*80 c7
c internal
      real*8 ra,rl,l
      real*8 rca,rcl,rfa,rfl,cmfa,cmfl

      LRc=tcore
      SLRc=pcore
      LRm=LRMass-LRc
      SLRm=SLRMass-SLRc

      return
      end

c-------------------------------------------
      subroutine erosion(c3,tmtot,tcore,tmant,pmtot,pcore,
     & pmant,LRMass,LRc,LRm,mdcore,mdmant)
c for now I'll have this for collision types 2, 3 and 4, but should
c make separate functions for each

      include 'swift.inc'

      real*8 tmtot,pmtot,tcore,pcore,tmant,pmant
      real*8 LRc,LRm
      real*8 mdcore,mdmant,LRMass
c parameters
      character*80 c3

      if(LRMass.ge.(tcore+pcore)) then
         if(c3.eq.'mincore'.and.LRMass.gt.tmtot) then
            LRc=tcore+pcore
            LRm=LRMass-LRc
            mdcore=0.
            mdmant=1.
         elseif(c3.eq.'mincore'.and.LRMass.le.tmtot) then
            LRc=tcore
            LRm=LRMass-LRc
            mdcore=pcore/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         elseif(c3.eq.'maxcore') then
            LRc=tcore+pcore
            LRm=LRMass-LRc
            mdcore=0.
            mdmant=1.
         endif
      elseif(LRMass.le.pcore) then
         LRc=LRMass
         LRm=0.
         mdcore=(tcore+pcore-LRMass)/(pmtot+tmtot-LRMass)
         mdmant=1.-mdcore
      elseif(LRMass.ge.tcore) then
         if(c3.eq.'mincore') then
            LRc=tcore
            LRm=LRMass-LRc
            mdcore=pcore/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         elseif(c3.eq.'maxcore') then
            LRc=LRMass
            LRm=0.
            mdcore=(tcore+pcore-LRMass)/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         endif
      elseif(LRMass.gt.pcore) then
         if(c3.eq.'mincore') then
            LRc=pcore
            LRm=LRMass-pcore
            mdcore=tcore/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         elseif(c3.eq.'maxcore') then
            LRc=LRMass
            LRm=0.
            mdcore=(tcore+pcore-LRMass)/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         endif
      endif

      return
      end

c-------------------------------------------
      subroutine partacc(c4,tmtot,tcore,tmant,pmtot,pcore,
     & pmant,LRMass,LRc,LRm,mdcore,mdmant)

      include 'swift.inc'

      real*8 tmtot,pmtot,tcore,pcore,tmant,pmant
      real*8 LRc,LRm
      real*8 mdcore,mdmant,LRMass
      character*80 c4
c
      if(LRMass.ge.(tcore+pcore)) then
         if(c4.eq.'mincore'.and.LRMass.gt.tmtot) then
            LRc=tcore+pcore
            LRm=LRMass-LRc
            mdcore=0.
            mdmant=1.
         elseif(c4.eq.'mincore'.and.LRMass.le.tmtot) then
            LRc=tcore
            LRm=LRMass-LRc
            mdcore=pcore
            mdmant=1.-mdcore
         elseif(c4.eq.'maxcore') then
            LRc=tcore+pcore
            LRm=LRMass-LRc
            mdcore=0.
            mdmant=1.
            write(*,*) 'cores merge partial accretion'
         endif
      elseif(LRMass.le.pcore) then
         LRc=LRMass
         LRm=0.
         mdcore=(tcore+pcore-LRMass)/(pmtot+tmtot-LRMass)
         mdmant=1.-mdcore
      elseif(LRMass.ge.tcore) then
         if(c4.eq.'mincore') then
            LRc=tcore
            LRm=LRMass-pcore
            mdcore=tcore/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         elseif(c4.eq.'maxcore') then
            LRc=LRMass
            LRm=0.
            mdcore=(pcore+tcore-LRMass)/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
            write(*,*) 'cores partially merge (type 4)'
         endif
      elseif(LRMass.gt.pcore) then
         if(c4.eq.'mincore') then
            LRc=pcore
            LRm=LRMass-pcore
            mdcore=tcore/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         elseif(c4.eq.'maxcore') then
            LRc=LRMass
            LRm=0.
            mdcore=(tcore+pcore-LRMass)/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         endif
      endif

      return
      end
c-------------------------------------------
      subroutine hitnspray(c5,tmtot,tcore,tmant,pmtot,pcore,
     & pmant,LRMass,LRc,LRm,mdcore,mdmant)

      include 'swift.inc'

      real*8 tmtot,pmtot,tcore,pcore,tmant,pmant
      real*8 LRc,LRm,SLRc,SLRm,smllcore,lrgcore
      real*8 mdcore,mdmant,LRMass,SLRMass,mdeb
c parameters
      character*80 c5

      mdeb=pmtot+tmtot-LRMass
      if(LRMass.gt.tmtot) then
         if(LRMass.ge.(pcore+tcore)) then
c cores only merge if core is too big to be part of debris
            LRc=tcore+pcore
            LRm=LRMass-LRc
            mdcore=0.
            mdmant=1.
            write(*,*) 'cores merge for hit and spray'
            if(LRm.lt.0) then
               write(*,*) 'mdcore is -ve for 1.2'
            endif
c         elseif(c5.eq.'maxcore') then
c            LRc=tcore+pcore
c            LRm=LRMass-LRc
c            mdcore=0.
c            mdmant=1.
         elseif(LRMass.lt.(pcore+tcore)) then
            LRc=LRMass
            LRm=0.
            mdcore=(pcore+tcore-LRMass)/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         endif
      elseif(LRMass.ge.tcore) then
c     ensure that core remains same unless debris is too small for pcore
         if(pcore.lt.mdeb) then
            LRc=tcore
            LRm=LRMass-LRc
            mdcore=pcore/(pmtot+tmtot-LRMass)
            mdmant=1.-mdcore
         elseif(pcore.gt.mdeb) then
            LRc=tcore+(pcore-mdeb)
            LRm=LRMass-LRc
            mdcore=1.
            mdmant=0.
         endif
      elseif(LRMass.lt.tcore) then
c         write(*,*) 'hitnspray LRMass too small'
         LRc=LRMass
         LRm=0.
         mdcore=(pcore+tcore-LRMass)/(pmtot+tmtot-LRMass)
         mdmant=1.-mdcore
      endif

      return
      end


c-------------------------------------------
      subroutine hitnrun(c6,tmtot,tcore,tmant,pmtot,pcore,pmant,
     & LRMass,SLRMass,LRc,LRm,SLRc,SLRm,mdcore,mdmant)

      include 'swift.inc'

      real*8 tmtot,pmtot,tcore,pcore,tmant,pmant
      real*8 LRc,LRm,SLRc,SLRm
      real*8 mdcore,mdmant,LRMass,SLRMass
c parameters
      character*80 c6
c internal
      real*8 corerem
c should I include this?
      real*8 c6cf !this is fraction of projectile core accreted

c cores are the same
c         elseif(c6.eq.'medcore') then
c only minimal core change
c            LRc=tcore+pcore*c6cf
c         endif
      if(SLRMass.ge.pcore.and.
     &         LRMass.ge.tcore) then
         LRc=tcore
         SLRc=pcore
         LRm=LRMass-LRc
         SLRm=SLRMass-SLRc
         mdcore=0.
         mdmant=1.
      elseif(SLRMass.lt.pcore) then
         SLRc=SLRMass
         SLRm=0.
         corerem=pcore-SLRMass
         if(c6.eq.'maxcore'.and.
     &   LRMass.ge.(tcore+corerem)) then
            LRc=tcore+corerem
            LRm=LRMass-LRc
            mdcore=0.
            mdmant=1.
         elseif(c6.eq.'maxcore'.and.
     &   LRMass.lt.(tcore+corerem)) then
            LRc=tcore
            LRm=LRMass-LRc
            mdcore=corerem/(pmtot+tmtot-LRMass-SLRMass)
            mdmant=1.-mdcore
         elseif(c6.eq.'mincore'.and.LRMass.ge.tcore) then
            LRc=tcore
            LRm=LRMass-LRc
            mdcore=corerem/(pmtot+tmtot-LRMass-SLRMass)
            mdmant=1.-mdcore
         elseif(LRMass.lt.tcore) then
            LRc=LRMass
            LRm=0.
            mdcore=(LRMass-LRc+corerem)/(pmtot+tmtot-
     &     LRMass-SLRMass)
            mdmant=1.-mdcore
         endif
      elseif(LRMass.lt.tcore) then
         LRc=LRMass
         LRm=0.
         SLRc=pcore
         SLRm=SLRMass-pcore
         mdcore=(tcore-LRMass)/(pmtot+tmtot-LRMass-SLRMass)
         mdmant=1.-mdcore
      endif

      return
      end

c-------------------------------------------
      subroutine GetRadius(m,r)
      real*8 m,r
      if (m.le.1.18554e-6) then ! 0.01 Earth masses
         r=1.01854e-3*m**0.333333
      else if(m.le.2.37108e-4) then ! 2.0 Earth masses
         r=6.34540e-4*m**0.298657
      else if(m.le.6.82761e-4) then ! 5.8 Earth masses
         r=4.31273e-4*m**0.252394
      else if(m.le.1.08189e-2) then ! 91 Earth masses
         r=1.55290e-2*m**0.744031
      else if(m.le.2.2943) then ! 20000 Earth masses
         r=4.49484e-4*m**(-0.038559)
      else ! stars
         r=2.18062e-4*m**0.832465
      endif
      return
      end

c-------------------------------------------
      subroutine getcorerad(m,cmf,crfrac)
      real*8 m,cmf,crfrac

c     function parameters
      real*8 a,ba,c,d,f

      a=0.300854247  !calculated in superearth-results
      ba=1.12704335
      c=1.23554869e-3
      d=0.43016866
      f=0.293045443

      if(cmf.eq.0) then   !artificially set boundaries
         crfrac = 0.
      else if(cmf.eq.1) then
         crfrac = 1.
      else
         crfrac=(a*cmf**ba/(c*m+d))+f
      endif

      return
      end


c************************************************************************
c     IO_INIT_CPARAM.F
c************************************************************************

c     Input for composition calculations
c     Parameters: options are 'mincore' or 'maxcore'
c     Format:
c     c2, c3, c4
c     c5, c6, c7

      subroutine io_init_cparam(infile,cname,c2,c3,c4,c5,c6,c7,CMFi,
     &   tout,dloss,fdebris,fdiscard)

      include '../swift.inc'
      include '../io/io.inc'

c.. Input
      character*(*) infile
      character*80 cname,c2,c3,c4,c5,c6,c7,fdebris,fdiscard
      real*8 CMFi, tout
      logical*2 dloss

      integer ierr

      write(*,*) 'Composition parameter file is ',infile
      call io_open(7,infile,'old','formatted',ierr)

      read(7,*) cname
      read(7,*) c2,c3,c4
      write(*,*)'c2, c3, c4:',c2,c3,c4
      read(7,*) c5,c6,c7
      write(*,*) 'c5, c6, c7:',c5,c6,c7
      read(7,*) CMFi,tout
      write(*,*) 'CMFi, tout:',CMFi,tout
      read(7,*) dloss
      write(*,*) 'dloss: ',dloss
      read(7,*) fdebris, fdiscard
      write(*,*) 'fdebris, fdiscard: ',fdebris,fdiscard

      close(unit=7)

      return
      end       !io_init_cparam.f
c************************************************************************
c     IO_INIT_PL_SYMBA.F
c************************************************************************
c     IO_INIT_PL_SYMBA reads in the data for the Sun and planets for
c     symba routines
c
c     Input:
c     infile        ==> File name to read from (character*80)
c     lclose        ==> .true. --> discard particle if it gets
c     too close to a planet. Read in that
c     distance in io_init_pl_symba
c     (logical*2 scalar)
c     iflgchk        ==>  bit 5 set ==>  include J2 and J4 terms
c
c     Output:
c     nbod          ==>  number of massive bodies (int scalar)
c     mass          ==>  mass of bodies (real array)
c     xh,yh,zh      ==>  initial position in Helio coord
c     (real arrays)
c     vxh,vyh,vzh   ==>  initial position in Helio coord
c     (real arrays)
c     rpl           ==>  physical size of planet
c     (real array)
c     rhill         ==>  size of planet's hills sphere
c     (real array)
c     j2rp2,j4rp4   ==>  J2*radii_pl^2 and  J4*radii_pl^4
c     (real scalars)
c
c     Remarks: Based on io_init_pl
c     Authors:  Hal Levison
c     Date:    11/21/96
c     Last revision: 1/10/97

      subroutine io_init_pl_symba(infile,lclose,iflgchk,nbod,npart,mass,
     &     xh,yh,zh,vxh,vyh,vzh,rpl,rhill,j2rp2,j4rp4)

      include '../swift.inc'
      include '../io/io.inc'

c...  Input
      character*(*) infile
      integer iflgchk
      logical*2 lclose

c...  Output
      real*8 mass(NTPMAX),rpl(NTPMAX),j2rp2,j4rp4
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX),rhill(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)
      integer nbod,npart

c...  Internal
      integer j,ierr,ibad
      real*8 r2hill(NTPMAX),rhrat

c-----
c...  Executable code
      iflgchk = 0
      write(*,*) 'Planet data file is ',infile
      call io_open(7,infile,'old','formatted',ierr)

c     Read number of planets
      read(7,*) nbod,npart

      if(nbod.gt.NTPMAX) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   The number of massive bodies,',nbod,','
         write(*,*) '   is too large, it must be less than',NTPMAX
         call util_exit(1)
      endif

c     write(*,23) nbod,nbod-npart
c23   format(/,'Number of bodies (incl. the Sun) is ',i4,' with ',i4
c    &,' debris')

c     For each planet read mass,
c     and helioc. position and vel .
      if(btest(iflgchk,5))  then ! bit 5 is set
         read(7,*) mass(1),j2rp2,j4rp4
      else
         read(7,*) mass(1)
         j2rp2 = 0.0d0
         j4rp4 = 0.0d0
      endif
      read(7,*) xh(1),yh(1),zh(1)
      read(7,*) vxh(1),vyh(1),vzh(1)
      rpl(1) = 0.0d0
      rhill(1) = 0.0d0

      if(  (xh(1).ne.0.0d0) .or.
     &     (yh(1).ne.0.0d0) .or.
     &     (zh(1).ne.0.0d0) .or.
     &     (vxh(1).ne.0.0d0) .or.
     &     (vyh(1).ne.0.0d0) .or.
     &     (vzh(1).ne.0.0d0) ) then
         write(*,*) ' SWIFT ERROR: in io_init_pl_symba: '
         write(*,*) '   Input MUST be in heliocentric coordinates '
         write(*,*) '   Position and Vel. of Massive body 1 .ne. 0'
         call util_exit(1)
      endif

      do j=2,nbod
         if(lclose) then
            read(7,*) mass(j),rhill(j),rpl(j)
         else
            read(7,*) mass(j),rhill(j)
         endif
         read(7,*) xh(j),yh(j),zh(j)
         read(7,*) vxh(j),vyh(j),vzh(j)
      enddo

      close(unit = 7)

c...  check to see if the hills spheres are ok
      call util_hills(nbod,mass,xh,yh,zh,vxh,vyh,vzh,r2hill)
      ibad = 0
      do j=2,nbod
         rhill(j) = sqrt(r2hill(j))
         rhrat = rhill(j)/sqrt(r2hill(j))
         if( (rhrat.gt.2.0) .or. (rhrat.lt.0.5) ) then
            ibad = ibad + 1
         endif
      enddo

      if(ibad.ne.0) then
         write(*,*) 'Warning in io_init_pl_symba:'
         write(*,*) '   Hill''s spheres are not consistent on ',
     &        ibad,' objects'
      endif

      return
      end                       ! io_init_pl_symba.f



c************************************************************************
c                          IO_INIT_PARAM.F
c************************************************************************
c INIT_PARAM reads in the parameters for the integration.
c
c      Input:
c            infile   ==> File name to read from (character*80)
c
c      Output:
c            t0       ==> Initial time (real scalar)
c            tstop    ==> final time (real scalar)
c            dt       ==> time step  (real scalar)
c            dtout    ==> time between binary outputs (real scalar)
c            dtdump   ==> time between dumps  (real scalar)
c            dsurv    ==>  fraction of debris that are gone
c            iflgchk  ==>  =0 don't run diagnostic routines
c                          bit 0 set ==>  write int*2 binary data file
c                          bit 1 set ==>  write real*4 binary file
c                          bit 2 set ==>  calc energy of system wrt time
c                          bit 3 set ==>  calc jacobi of the test particles
c                          bit 4 set ==>  check if particles are removed
c                          bit 5 set ==>  include J2 and J4 terms
c      rmin,rmax      ==>  maximum and min distance from Sun
c                                if <0  then don't check
c                                    (real scalar)
c      rmaxu          ==>  maximum distance from Sun in not bound
c                                 if <0  then don't check
c                                      (real scalar)
c       qmin          ==> Smallest perihelion distance
c                                 if <0  then don't check
c                                      (real scalar)
c       lclose        ==> .true. --> discard particle if it gets
c                                    too close to a planet. Read in that
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c       outfile       ==>  Name of binary output file (character*80)
c       fopenstat     ==>  The status flag for the open statements of the
c                          output files.  Must be one of the following:
c                                 new      (die if the file exists)
c                                 append   (add to what is there)
c                                 unknown  (just write over what is there)
c                                 (character*80)
c
c
c Remarks:
c Authors:  Martin Duncan
c Date:    3/2/93
c Last revision:  5/10/94  HFL

      subroutine io_init_param(infile,t0,tstop,dt,dtout,dtdump,twall,
     &         iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,
     &         fopenstat)

      include '../swift.inc'
      include '../io/io.inc'

c...    Input
      character*(*) infile

c...  Outputs:
      integer iflgchk
      real*8 t0,tstop,dt,dsurv
      real*8 dtout,dtdump,twall
      real*8 rmin,rmax,rmaxu,qmin
      logical*2 lclose
      character*80 outfile,fopenstat

c...  Internals
      logical*1 lflg(0:IO_NBITS-1)
      integer i,ierr

c-----
c...  Executable code

c  Jen 2018 added in twall to the param file
      write(*,*) 'Parameter data file is ',infile
      call io_open(7,infile,'old','formatted',ierr)
      read(7,*) t0,tstop,dt
      write(*,*) 't0,tstop,dt : ',t0,tstop,dt
      read(7,*) dtout,dtdump,twall
      write(*,*) 'dtout,dtdump,twall : ',dtout,dtdump,twall
c      read(7,*) dsurv
c      write(*,*) 'dsurv :',dsurv
      read(7,*) (lflg(i),i=IO_NBITS-1,0,-1)

      iflgchk=0
      do i=0,IO_NBITS-1
         if(lflg(i)) then
            iflgchk = ibset(iflgchk,i)
         endif
      enddo

      write(*,*) (lflg(i),i=IO_NBITS-1,0,-1),' = ',iflgchk

      if(btest(iflgchk,0) .and. btest(iflgchk,1))  then
         write(*,*) ' SWIFT ERROR: in io_init_param:'
         write(*,*) '    Invalid logical flags '
         write(*,*) '    You cannot request that both a real and ',
     &                '       an integer binary file be written '
         call util_exit(1)
      endif

      if(btest(iflgchk,4))  then ! bit 4 is set
         read(7,*) rmin,rmax,rmaxu,qmin,lclose
         write(*,*) 'rmin,rmax,rmaxu,qmin,lclose :',rmin,
     &          rmax,rmaxu,qmin,lclose
      else
         rmin = -1.0
         rmax = -1.0
         rmaxu = -1.0
         qmin = -1.0
         lclose = .false.
      endif

      if(btest(iflgchk,0) .or. btest(iflgchk,1))  then
         read(7,999) outfile
 999     format(a)
         write(*,*) 'outfile : ', outfile
         write(*,*) ' '
      endif

      read(7,999) fopenstat
      if(  (fopenstat(1:3).ne.'new') .and.
     &       (fopenstat(1:3).ne.'NEW') .and.
     &       (fopenstat(1:7).ne.'unknown') .and.
     &       (fopenstat(1:7).ne.'UNKNOWN') .and.
     &       (fopenstat(1:6).ne.'append') .and.
     &       (fopenstat(1:6).ne.'APPEND') ) then
         write(*,*) ' SWIFT ERROR: in io_init_param:'
         write(*,*) '    Invalid status flag:',fopenstat(1:7),':'
         call util_exit(1)
      endif

      close(unit = 7)

      return
      end     ! io_init_param

c************************************************************************
c                          IO_INIT_PARAM_DLOSS.F
c************************************************************************
c INIT_PARAM reads in the parameters for the integration.
c
c      Input:
c            infile   ==> File name to read from (character*80)
c
c      Output:
c            t0       ==> Initial time (real scalar)
c            tstop    ==> final time (real scalar)
c            dt       ==> time step  (real scalar)
c            dtout    ==> time between binary outputs (real scalar)
c            dtdump   ==> time between dumps  (real scalar)
c            dsurv    ==>  fraction of debris that survives the collision
c            iflgchk  ==>  =0 don't run diagnostic routines
c                          bit 0 set ==>  write int*2 binary data file
c                          bit 1 set ==>  write real*4 binary file
c                          bit 2 set ==>  calc energy of system wrt time
c                          bit 3 set ==>  calc jacobi of the test particles
c                          bit 4 set ==>  check if particles are removed
c                          bit 5 set ==>  include J2 and J4 terms
c      rmin,rmax      ==>  maximum and min distance from Sun
c                                if <0  then don't check
c                                    (real scalar)
c      rmaxu          ==>  maximum distance from Sun in not bound
c                                 if <0  then don't check
c                                      (real scalar)
c       qmin          ==> Smallest perihelion distance
c                                 if <0  then don't check
c                                      (real scalar)
c       lclose        ==> .true. --> discard particle if it gets
c                                    too close to a planet. Read in that
c                                    distance in io_init_pl
c                                      (logical*2 scalar)
c       outfile       ==>  Name of binary output file (character*80)
c       fopenstat     ==>  The status flag for the open statements of the
c                          output files.  Must be one of the following:
c                                 new      (die if the file exists)
c                                 append   (add to what is there)
c                                 unknown  (just write over what is there)
c                                 (character*80)
c
c
c Remarks:
c Authors:  Martin Duncan
c Date:    3/2/93
c Last revision:  5/10/94  HFL

      subroutine io_init_param_dloss(infile,t0,tstop,dt,dtout,dtdump,
     &    twall,dsurv,iflgchk,rmin,rmax,rmaxu,qmin,lclose,outfile,
     &         fopenstat)

      include '../swift.inc'
      include '../io/io.inc'

c...    Input
      character*(*) infile

c...  Outputs:
      integer iflgchk
      real*8 t0,tstop,dt,dsurv
      real*8 dtout,dtdump,twall
      real*8 rmin,rmax,rmaxu,qmin
      logical*2 lclose
      character*80 outfile,fopenstat

c...  Internals
      logical*1 lflg(0:IO_NBITS-1)
      integer i,ierr

c-----
c...  Executable code

c  Jen 2018 added in twall to the param file
      write(*,*) 'Parameter data file is ',infile
      call io_open(7,infile,'old','formatted',ierr)
      read(7,*) t0,tstop,dt
      write(*,*) 't0,tstop,dt : ',t0,tstop,dt
      read(7,*) dtout,dtdump,twall
      write(*,*) 'dtout,dtdump,twall : ',dtout,dtdump,twall
      read(7,*) dsurv
      write(*,*) 'dsurv :',dsurv
      read(7,*) (lflg(i),i=IO_NBITS-1,0,-1)

      iflgchk=0
      do i=0,IO_NBITS-1
         if(lflg(i)) then
            iflgchk = ibset(iflgchk,i)
         endif
      enddo

      write(*,*) (lflg(i),i=IO_NBITS-1,0,-1),' = ',iflgchk

      if(btest(iflgchk,0) .and. btest(iflgchk,1))  then
         write(*,*) ' SWIFT ERROR: in io_init_param:'
         write(*,*) '    Invalid logical flags '
         write(*,*) '    You cannot request that both a real and ',
     &                '       an integer binary file be written '
         call util_exit(1)
      endif

      if(btest(iflgchk,4))  then ! bit 4 is set
         read(7,*) rmin,rmax,rmaxu,qmin,lclose
         write(*,*) 'rmin,rmax,rmaxu,qmin,lclose :',rmin,
     &          rmax,rmaxu,qmin,lclose
      else
        rmin = -1.0
        rmax = -1.0
        rmaxu = -1.0
        qmin = -1.0
        lclose = .false.
      endif

      if(btest(iflgchk,0) .or. btest(iflgchk,1))  then
         read(7,999) outfile
 999     format(a)
         write(*,*) 'outfile : ', outfile
         write(*,*) ' '
      endif

      read(7,999) fopenstat
      if(  (fopenstat(1:3).ne.'new') .and.
     &       (fopenstat(1:3).ne.'NEW') .and.
     &       (fopenstat(1:7).ne.'unknown') .and.
     &       (fopenstat(1:7).ne.'UNKNOWN') .and.
     &       (fopenstat(1:6).ne.'append') .and.
     &       (fopenstat(1:6).ne.'APPEND') ) then
         write(*,*) ' SWIFT ERROR: in io_init_param:'
         write(*,*) '    Invalid status flag:',fopenstat(1:7),':'
         call util_exit(1)
      endif

      close(unit = 7)

      return
      end     ! io_init_param


