c   Gererate distributed initial position and velocity for embryos
c   and planetesimals with power law
c   modified by Jen 2019
      include 'swift.inc'

      integer NPLR
c      parameter(NPLR=9)
      real*8 SMASSYR
      parameter(SMASSYR=TWOPI*TWOPI)

      real*8 mass(NTPMAX),rpl(NTPMAX),rhill(NTPMAX)
      real*8 xh(NTPMAX),yh(NTPMAX),zh(NTPMAX)
      real*8 vxh(NTPMAX),vyh(NTPMAX),vzh(NTPMAX)

      real*8 xb(NTPMAX),yb(NTPMAX),zb(NTPMAX)
      real*8 vxb(NTPMAX),vyb(NTPMAX),vzb(NTPMAX)

      real*8 xcm,ycm,zcm
      real*8 vxcm,vycm,vzcm

      real*8 masso(NTPMAX),rplo(NTPMAX),j2rp2,j4rp4
      real*8 xho(NTPMAX),yho(NTPMAX),zho(NTPMAX),rhillo(NTPMAX)
      real*8 vxho(NTPMAX),vyho(NTPMAX),vzho(NTPMAX)

      real*8 xbo(NTPMAX),ybo(NTPMAX),zbo(NTPMAX)
      real*8 vxbo(NTPMAX),vybo(NTPMAX),vzbo(NTPMAX)

      real*8 gm,a,e,inc,capom,omega,capm,dr,au
      real*8 arand,erand,incrand,caporand,omegarand,capmrand
c      real*8 rmass(NPLR),epoch,year,msys,rphy(NPLR),munit
      real*8 epoch,year,msys,munit,num,mt,mearth
      real*8 fachill,r2hill(NTPMAX),dum
      logical*2 lclose, fjup
      integer iflgchk,id
      character*80 inplfile,infile
      real*8 totmass,density,spac,pmass,vd,a0,a1,a2,minemb
      real*8 tot,tot2,emass,fac,Miso,rh,alp,g,b,am,x,emb,embm
      integer iseed,nembryo,j,k,ninit,nfin,np
      real*8 meane, meani, sigmae, sigmai

      integer nbod,ialpha,i,ip1,ip2,iuflg,icflg,irflg

      data dr/1.7453292519943294d-2/ 
      data epoch/2449101.0d0/
      data year/365.2422d0/
      data au/1.495978707d8/    ! km
      data munit/5.0428958d31/ !in g
      data mearth/5.972d27/ !in g 

 1    write(*,*) ' Units Menu: '
      write(*,*) '       0 ==> Solar masses, and AU '
      write(*,*) '       1 ==> AU, and Years '
      read(*,*) iuflg
      if( (iuflg.ne.0) .and. (iuflg.ne.1) ) goto 1

 2    write(*,*) ' Coordinate Menu: '
      write(*,*) '       0 ==> ecliptic '
      write(*,*) '       1 ==> invariable plane '
      read(*,*) icflg
      if( (icflg.ne.0) .and. (icflg.ne.1) ) goto 2

 3    write(*,*) ' Planetary Radius Menu: '
      write(*,*) '       0 ==> Do not include a radius '
      write(*,*) '       1 ==> Include the physical radius of planet '
      read(*,*) irflg
      if( (irflg.ne.0) .and. (irflg.ne.1) ) goto 3

c      nbod = NPLR + 1

      if(iuflg.ne.1) then
         mass(1) = 1.0d0
      else
         mass(1) = SMASSYR
      endif

      write(*,*) ' Mass of the Sun is ',mass(1)

      xh(1) = 0.0
      yh(1) = 0.0
      zh(1) = 0.0
      vxh(1) = 0.0
      vyh(1) = 0.0
      vzh(1) = 0.0

c     Read in param file for initial values
      infile = 'init_params_ex.in'
      call init_iparam(infile,totmass,alp,np,a0,a2,meane,meani)
      nbod=1
      density=2.01d8 !density of material is 3g/cm

c equal to 25Me- 859 = 68.4*4*pi
c for 50Me embryo7s - 1570 = 125*4*pi
c for 75Me (37.5 embryos) = 102.5*4*pi = 1288
c for alp = 0, 25 Me, const= 678
c      totmass=(am*4.*PI*(1.d0**(2.-alp)-0.05d0**(2.-alp)))*(au*1d5)**2
c      totmass=3.75d1


      write(*,*) totmass
      
c convert to munits if iuflg = 1
      if(iuflg.eq.1) then 
         totmass=totmass*mearth/munit
         write(*,*) 'totmass in munit',totmass
      endif

      if(irflg.eq.0) then
         lclose = .false.
      else if(irflg.eq.1) then
         lclose = .true.
      endif

      ip1 = 2

      call random_seed()

      pmass=totmass/real(np)
      write(*,*) pmass
c      sd=7*(a/(1d4*au))**(-3./2.) !in g/cm^2
      call random_seed()

      sigmae = sqrt(2/PI)*meane
      sigmai = sqrt(2/PI)*meani
      
      do i=ip1,ip1+np
         call random_number(arand)
c draw from power law distribution (surface density)
         g=1.-alp
c         a=arand*(a2-a0)+a0
         a=(a0**g+(a2**g-a0**g)*arand)**(1./g)
         mass(i)=pmass
         rpl(i)=(3.d0*mass(i)/(4.d0*PI*density))**(1./3.)
         rhill(i)=a*(mass(i)/(mass(1)*3d0))**(1./3.)
         call random_number(erand)
         call rayleigh(sigmae,erand,e)
         do while (e.gt.0.8) !make sure e < 0.8
            call random_number(erand)
            call rayleigh(sigmae,erand,e)
         enddo
         call random_number(incrand)
         call rayleigh(sigmai,incrand,inc)
         call random_number(capmrand)
         capm = TWOPI*capmrand
         call random_number(omegarand)
         omega = TWOPI*omegarand
         call random_number(caporand)
         capom = TWOPI*caporand
         gm = mass(1)+mass(i)
         ialpha = -1
         call orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,
     &      xh(i),yh(i),zh(i),vxh(i),vyh(i),vzh(i))

      enddo
      
      nbod=nbod+np
      write(*,*) 'nbod is',nbod

c not sure if I need this at all
      call coord_h2b(nbod,mass,xh,yh,zh,vxh,vyh,vzh,
     &     xb,yb,zb,vxb,vyb,vzb,msys)
      
      if(icflg.eq.1) then
        call invar(nbod,mass,xb,yb,zb,vxb,vyb,vzb)
      endif


      write(*,*) ' '
      write(*,*) 'Enter name of planet data file : '
      read(*,999) inplfile
 999  format(a)
      iflgchk = 0
      j2rp2 = 0.0d0
      j4rp4 = 0.0d0
      call io_dump_pl_symba(inplfile,nbod,mass,xh,yh,zh,
     &     vxh,vyh,vzh,lclose,iflgchk,rpl,rhill,j2rp2,j4rp4)

      stop
      end

c--------------------------------------------------------------------
      subroutine rayleigh(sigma,x,new)

      real*8 sigma, x, new

      new = sqrt(-2*sigma**2*log(1-x))

      return
      end

c--------------------------------------------------------------------
      subroutine invar(nbod,mass,x,y,z,vx,vy,vz)

      include 'swift.inc'

      integer i,nbod

      real*8 mass(nbod)
      real*8 x(nbod),y(nbod),z(nbod)
      real*8 vx(nbod),vy(nbod),vz(nbod)
      real*8 xt,yt,zt
      real*8 vxt,vyt,vzt

      real*8 c1,c2,c3,c,bot
      real*8 a(3,3)

      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0

      do i=1,nbod
         c1 = c1 + mass(i)*( y(i)*vz(i) - z(i)*vy(i) )
         c2 = c2 + mass(i)*( z(i)*vx(i) - x(i)*vz(i) )
         c3 = c3 + mass(i)*( x(i)*vy(i) - y(i)*vx(i) )
      enddo

      write(*,*) ' c1,c2,c3 ',c1,c2,c3

      c = sqrt( c1*c1 + c2*c2 + c3*c3 )
      c1 = c1/c
      c2 = c2/c
      c3 = c3/c
      bot = 1.0d0/(1.0d0 + c3)

      a(1,1) = 1.0d0 - c1*c1*bot
      a(1,2) = -1.0d0*c1*c2*bot
      a(1,3) = -1.0d0*c1

      a(2,1) = -1.0d0*c1*c2*bot
      a(2,2) = 1.0d0 - c2*c2*bot
      a(2,3) = -1.0d0*c2

      a(3,1) = c1
      a(3,2) = c2
      a(3,3) = c3

      do i=1,nbod
          xt = a(1,1)*x(i) + a(1,2)*y(i) + a(1,3)*z(i) 
          yt = a(2,1)*x(i) + a(2,2)*y(i) + a(2,3)*z(i) 
          zt = a(3,1)*x(i) + a(3,2)*y(i) + a(3,3)*z(i) 
          x(i) = xt
          y(i) = yt
          z(i) = zt
          vxt = a(1,1)*vx(i) + a(1,2)*vy(i) + a(1,3)*vz(i) 
          vzt = a(3,1)*vx(i) + a(3,2)*vy(i) + a(3,3)*vz(i) 
          vx(i) = vxt
          vy(i) = vyt
          vz(i) = vzt
       enddo

c..    check

      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0

      do i=1,nbod
         c1 = c1 + mass(i)*( y(i)*vz(i) - z(i)*vy(i) )
         c2 = c2 + mass(i)*( z(i)*vx(i) - x(i)*vz(i) )
         c3 = c3 + mass(i)*( x(i)*vy(i) - y(i)*vx(i) )
      enddo

      write(*,*) ' c1,c2,c3 ',c1,c2,c3

      return
      end

****************************************************
c inputs are:
c
c totmass: mass of disk in Earth masses
c     alp: slope of disk
c      am: coefficient for Miso formula
c     emb: maximum embryo mass
c  minemb: minimum embryo mass
c nembryo: number of embryos
c      a0: inner disk limit
c      a2: outer disk limit
c       b: hill radius spacing
c      np: number of planetesimals
c
****************************************************
      subroutine init_iparam(infile,totmass,alp,np,a0,a2,meane,meani)

      include '../swift.inc'
      include '../io/io.inc'


      character*(*) infile
      real*8 totmass,alp,am,emb,minemb,a0,a2,b
      real*8 meane, meani
      integer np,nembryo
      integer ierr

      write(*,*) 
      call io_open(7,infile,'old','formatted',ierr)

      read(7,*) totmass,alp,np
      write(*,*) 'totmass, alp, np:',totmass,alp,np
      read(7,*) a0, a2
      write(*,*) 'amin, amax',a0,a2
      read(7,*) meane,meani
      write(*,*) 'meane, meani',meane,meani

      close(unit=7)

      return
      end    !init_iparam



