c**************************************************************************
c  Property code written by R K Choudhury and B K Sahoo, IIA,             *
c  Bangalore. Code is modified by Narendra Nath Dutta, IIT, KGP in the    *
c  calcualtions of E1 transition(len), E1 transition(vel), E2 transition  *
c  (len), M1 transition, hyperfine A and B constants.                     *
c                                                                         *
c  PNC part is incomplete                                                 *
c                                                                         *
c               Last modified on 05-02-2013                               *
c**************************************************************************
c  cpp
c
c**************************************************************************
c     Sample input.dbar for neutral K                                     *
c                                                                         *
c      7  14 12 12 11 11 10 10                                            *
c         3 2 2 0 0 0 0                                                   *
c         1 0 1                                                           *
c         3 1 39 19                                                       *
c         1 2 1 1 1 0 1 0 0 0 0 1                                         *
c                                                                         *
c**************************************************************************
c
c dimension variables
#define      MN     750
#define   MNBAS     150
#define   MNOCC     100
#define   MNEXC     100
#define   MNSYM      11
#define     MNS       6
#define     MXV      11
#define    MDIM     18000000
c#define   N2ORB    MNBAS*(MNBAS+1)/2
c#define   N2BAS    MNBAS*MNBAS

c constants
#define    STDIN     5
#define    STDOUT    6
#define    WFNIN     7
#define    LINEAR    15
#define    NTFILE    16
#define    DIPOLE    17
#define    DIPOLEV   18
#define    HYPER     19
#define    QUADPULE  20
#define    QUADPULEV 21
#define    QUADB     22
#define    M1PROP    23
#define    M2PROP    24
#define    GFACTOR   25
#define    GFACTORC  26
#define    OVERLAP   27
#define    OCTC      28
#define    NSI       29
#define    NSD       30



c*******************************************************************
c     This programme calculates the value of the diagrams
c*******************************************************************
      
      implicit real*8(a-h,o-z)

      dimension nbas(MNSYM),iocc(MNSYM)
      dimension ti(MDIM),d(MNBAS,MNBAS),dc(MNBAS,MNBAS)
      dimension tf(MDIM)

      call daopen
      call symmcal
      call setgrd
      call readinp(nsym,nbas,iocc,iopt,nprint,idip,ivel,iquad,
     :ihyp,iquadB,ioctC,im1p,im2p,igfact,iovlp,insi,insd,inms)
      call setup(nsym,nbasis,nbas,iocc)
      call grasprd(nbasis)
      call setqic

       
c--------------------------------------------------------------
c     T1 and T2 cluster operators are read
c--------------------------------------------------------------                                                                                 
      rewind NTFILE
      read(NTFILE)n,(ti(i),i=1,n)                                        
      call symm(iopt,nprint,ti,tf)
      print*,'I am over1'

c---------------------------------------------------------------
c      T+ T
c     e  e  diagrams are calculated
c---------------------------------------------------------------

      if(iovlp.eq.1)then
      open(OVERLAP,file='overlap.dat',form='unformatted',
     :         status='unknown')
      call cinit(d)
      call overlap(d,ti,tf)
      write(OVERLAP)((d(i,j),i=1,nbasis),j=1,nbasis)
      endif
      print*,'I am over2'


c---------------------------------------------------------------
c              O = Electric dipole transition operator (Len Gg)
c       T+   T
c      e  O e   digrams are calculated
c---------------------------------------------------------------

      if(idip.eq.1)then
      if(ivel.eq.0)goto 1
      if(ivel.eq.2)goto 1
          goto 2
  1   continue
      call cinit(d)
      call Calce1mat(nbasis)
      open(DIPOLE,file='dipole.dat',form='unformatted',status='unknown')
      call dipolehh(d,ti,tf)
      call dipolehp(d,ti,tf)
      call dipolepp(d,ti,tf)
      write(DIPOLE)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(DIPOLE)
  2   continue


c----------------------------------------------------------------
c              O = Electric dipole transition operator (Vel Gg)
c       T+   T
c      e  O e   digrams are calculated
c----------------------------------------------------------------


      if(ivel.eq.1)goto 3
      if(ivel.eq.2)goto 3
          goto 4
  3   continue
      call cinit(d)
      call Calce1vel(nbasis)
      open(DIPOLEV,file='dipolev.dat',form='unformatted',
     :status='unknown')
      call dipvelhh(d,ti,tf)
      call dipvelhp(d,ti,tf)
      call dipvelpp(d,ti,tf)
      write(DIPOLEV)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(DIPOLEV)
  4   continue

      endif

      print*,'I am over3'


c-------------------------------------------------------------------
c              O = Electric quadrupole transition operator (Len Gg)
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

      if(iquad.eq.1)then
      if(ivel.eq.0)goto 5
      if(ivel.eq.2)goto 5
          goto 6
  5   continue
      call cinit(d)
      call Calquad(nbasis)
      open(QUADPULE,file='quadpule.dat',form='unformatted',
     :status='unknown')
      call quadpulehh(d,ti,tf)
      call quadpulehp(d,ti,tf)
      call quadpulepp(d,ti,tf)
      write(QUADPULE)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(QUADPULE)
  6   continue

c-------------------------------------------------------------------
c              O = Electric quadrupole transition operator (Vel Gg)
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

      if(ivel.eq.1)goto 7
      if(ivel.eq.2)goto 7
          goto 8
  7   continue
      call cinit(d)
      call Calquadvel(nbasis)
      open(QUADPULEV,file='quadpulev.dat',form='unformatted',
     :status='unknown')
      call quadruvelhh(d,ti,tf)
      call quadruvelhp(d,ti,tf)
c      call quadruvelph(d,ti,tf)
      call quadruvelpp(d,ti,tf)
      write(QUADPULEV)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(QUADPULEV)
  8   continue
      endif

      print*,'I am over4'

c-------------------------------------------------------------------
c                 O = Magnetic dipole hyperfine operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

      if(ihyp.eq.1)then
      call cinit(d)
      call Calchyp(nbasis)
      open(HYPER,file='hyper.dat',form='unformatted',status='unknown')
      call hyperhh(d,ti,tf)
      call hyperhp(d,ti,tf)
      call hyperpp(d,ti,tf)
      write(HYPER)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(HYPER)
      endif

      print*,'I am over5'

c-------------------------------------------------------------------
c                 O = Electric quadrupole hyperfine operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

      if(iquadB.eq.1)then
      call cinit(d)
      call CalcB(nbasis)
      open(QUADB,file='quadb.dat',form='unformatted',status='unknown')
      call quadbhh(d,ti,tf)
      call quadbhp(d,ti,tf)
      call quadbpp(d,ti,tf)
      write(QUADB)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(QUADB)
      endif

      print*,'I am over6'

c-------------------------------------------------------------------
c                 O = Magnetic octupole hyperfine operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

      if(ioctC.eq.1)then
      call cinit(d)
      call CalcC(nbasis)
      open(OCTC,file='octc.dat',form='unformatted',status='unknown')
      call octchh(d,ti,tf)
      call octchp(d,ti,tf)
      call octcpp(d,ti,tf)
      write(OCTC)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(OCTC)
      endif

      print*,'I am over7'

c-------------------------------------------------------------------
c                 O = Magnetic dipole transition operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

      if(im1p.eq.1)then
      call cinit(d)
      call Calcm1prop(nbasis)
      open(M1PROP,file='m1prop.dat',form='unformatted',status='unknown')
      call m1prophh(d,ti,tf)
      call m1prophp(d,ti,tf)
      call m1proppp(d,ti,tf)
      write(M1PROP)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(M1PROP)
      endif

      print*,'I am over8'

c-------------------------------------------------------------------
c                 O = Magnetic octupole transition operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

      if(im2p.eq.1)then
      call cinit(d)
      call Calcm2prop(nbasis)
      open(M2PROP,file='m2prop.dat',form='unformatted',status='unknown')
      call m2prophh(d,ti,tf)
      call m2prophp(d,ti,tf)
      call m2proppp(d,ti,tf)
      write(M2PROP)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(M2PROP)
      endif

      print*,'I am over9'

c-------------------------------------------------------------------
c                 O = G factor operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------


      if(igfact.eq.1)then
      call cinit(d)
      call Calcgfact(nbasis)
      open(GFACTOR,file='gfact.dat',form='unformatted',status='unknown')
      open(GFACTORC,file='gfactc.dat',form='unformatted',
     :status='unknown')
      call gfacthh(d,dc,ti,tf)
      call gfacthp(d,dc,ti,tf)
      call gfactpp(d,dc,ti,tf)
      write(GFACTOR)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(GFACTOR)
      write(GFACTORC)((dc(i,j),i=1,nbasis),j=1,nbasis)
      close(GFACTORC)
      endif

      print*,'I am over10'

c-------------------------------------------------------------------
c                 O = Nuclear spin independent (NSI) operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

       if(insi.eq.1)then
       call cinit(d)
       call Calcnsira(nbasis)
       open(NSI,file='nsi.dat',form='unformatted',status='unknown')
       call nsihh(d,ti,tf)
       call nsihp(d,ti,tf)
       call nsipp(d,ti,tf)
       write(NSI)((d(i,j),i=1,nbasis),j=1,nbasis)
       close(NSI)
       endif

       print*,'I am over11'


c-------------------------------------------------------------------
c                 O = Nuclear spin dependent (NSD) operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

       if(insd.eq.1)then
       call cinit(d)
       call Calcnsdra(nbasis)
       open(NSD,file='nsd.dat',form='unformatted',status='unknown')
       call nsdhh(d,ti,tf)
       call nsdhp(d,ti,tf)
       call nsdpp(d,ti,tf)
       write(NSD)((d(i,j),i=1,nbasis),j=1,nbasis)
       close(NSD)
       endif

       print*,'I am over12'
c-------------------------------------------------------------------
c                 O = Magnetic dipole hyperfine operator
c       T+   T
c      e  O e   digrams are calculated
c-------------------------------------------------------------------

      if(inms.eq.1)then
      call cinit(d)
      open(NMS,file='nms.dat',form='unformatted',status='unknown')
      call nmshh(d,ti,tf)
      call nmshp(d,ti,tf)
      call nmspp(d,ti,tf)
      write(NMS)((d(i,j),i=1,nbasis),j=1,nbasis)
      close(NMS)
      endif

       call daclose

       stop
       end

c************************************************************************

      subroutine findk2(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)
c     this subroutine finds the common k value
      implicit real*8 (a-h,o-z)
      real*8 jj
      dimension ja(MXV),jb(MXV),
     :          jp(MXV),jq(MXV),jeven(MXV)

      do i=1,MXV
      jeven(i)=0
      enddo

      iab=ia*ib
      ipq=ip*iq

      jminab=idint(dmax1(orba,orbb)-dmin1(orba,orbb))
      jminpq=idint(dmax1(orbp,orbq)-dmin1(orbp,orbq))

      jmaxab=idint(orba+orbb)
      jmaxpq=idint(orbp+orbq)

      ll=0
      kk=0
      do i=jminab,jmaxab
      jtot=jmaxab+i
      jj=float(jtot)-2.0*(float(jtot/2))
      if(jj.ne.0.and.iab.gt.0)then
      ll=ll+1
      ja(ll)=i
      endif
      if(jj.eq.0.and.iab.lt.0)then
      kk=kk+1
      jb(kk)=i
      endif
      enddo

      mm=0
      nn=0
      do i=jminpq,jmaxpq
      jtot=jmaxpq+i
      jj=float(jtot)-2.0*(float(jtot/2))
      if(jj.ne.0.and.ipq.gt.0)then
      mm=mm+1
      jp(mm)=i
      endif
      if(jj.eq.0.and.ipq.lt.0)then
      nn=nn+1
      jq(nn)=i
      endif
      enddo


      kkk=0
      do i=1,ll
      do j=1,mm
      if(ja(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
      endif
      enddo
      enddo

      do i=1,ll
      do j=1,nn
      if(ja(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
      endif
      enddo
      enddo

      do i=1,kk
      do j=1,mm
      if(jb(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
      endif
      enddo
      enddo

      do i=1,kk
      do j=1,nn
      if(jb(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
      endif
      enddo
      enddo
      return
      end
 
      subroutine init(a,n)
c     initialization
      real*8 a(n)
      do i=1,n
      a(i)=0
      enddo
      return
      end

      subroutine iloc(i,j,k,l,ij,kl)
c     pack two index to one index
      implicit integer (a-z)
      ij=j*1000+i
      kl=l*1000+k
      return
      end

        subroutine sixj(rj1,rj2,rj3,rl1,rl2,rl3,s6j)
        implicit real*8 (a-h,o-z)
        common/factor/fct(0:50),mfd
!comment function s6j(j1,j2,j3,l1,l2,l3) calculates the 6j symbol
!        of its arguments.  it returns an error message if the
!        sum of angular momenta involved in a triangular
!        relationship is not an integer;
!for all j1,j2,j3,l1,l2,l3 such that
!        fixp(2*j1) and fixp(2*j2) and fixp(2*j3) and
!        fixp(2*l1) and fixp(2*l2) and fixp(2*l3)
! let s6j(j1,j2,j3,l1,l2,l3) =
!    begin scalar ws6j,delprod,zmin,zmax,tria,trib,tric,trid,s;

        if(mfd.eq.33)goto 704
        mfd=33
        fct(0)=1.
        fct(1)=1.
        do 20 i=2,33
        fct(i)=fct(i-1)*i
20      continue

704    tria=tri6j(rj1,rj2,rj3)
       trib=tri6j(rl1,rj2,rl3)
       tric=tri6j(rj1,rl2,rl3)
       trid=tri6j(rl1,rl2,rj3)
c       type *,tria,trib,tric,trid
       if(tria.eq.-1..or.trib.eq.-1..or.tric.eq.-1..or.trid.eq.-1.)then
       print *,'inconsistent arguments to 6j symbol'
       print *,tria,trib,tric,trid
       print *,rj1,rj2,rj3
       print *,rl1,rl2,rl3
        return
        endif
       if(tria.eq.0..or.trib.eq.0..or.tric.eq.0..or.trid.eq.0.)then
c         type *,'triangle condition'
          s6j=0
        return
        endif
!    comment apply formula in edmonds [(6.3.7) p.99]

       delprod=tridel(rj1,rj2,rj3)*tridel(rj1,rl2,rl3)*
     *         tridel(rl1,rj2,rl3)*tridel(rl1,rl2,rj3)
c       type *,'delprod=',delprod
c       type *,'args to min',(rj1+rj2+rl1+rl2),(rj2+rj3+rl2+rl3),
c     *                       (rj3+rj1+rl3+rl1)
       izmin=dmax1((rj1+rj2+rj3),(rj1+rl2+rl3),(rl1+rj2+rl3)
     *           ,(rl1+rl2+rj3))
       izmax=dmin1((rj1+rj2+rl1+rl2),(rj2+rj3+rl2+rl3)
     *           ,(rj3+rj1+rl3+rl1))
        ws6j=0
c       type *,'izmin,izmax=',izmin,izmax
        do 10 iz=izmin,izmax
        z=iz
         ws6j=ws6j+((-1)**iz)*fct(int(z+1.))/(fct(z-rj1-rj2-rj3)*
     *           fct(int(z-rj1-rl2-rl3))*
     *           fct(int(z-rl1-rj2-rl3))*
     *           fct(int(z-rl1-rl2-rj3))*
     *           fct(int(rj1+rj2+rl1+rl2-z))*
     *           fct(int(rj2+rj3+rl2+rl3-z))*
     *           fct(int(rj3+rj1+rl3+rl1-z)))
10      continue
       s6j=delprod*ws6j
       return
       end

        function tridel(rm1,rm2,rm3)
        implicit real*8 (a-h,o-z)
        common/factor/fct(0:50),mfd

c       type *,'+++',rm1,rm2,rm3
!comment function tridel evaluates the delta symbol defined on
!        p.99 of edmonds;
!for all m1,m2,m3 such that fixp(2*m1) and fixp(2*m2) and fixp(2*m3)
c       type *,m1+m2-m3,m1-m2+m3,-m1+m2+m3,m1+m2+m3+1

        i1=rm1+rm2-rm3
        i2=rm1-rm2+rm3
        i3=-rm1+rm2+rm3
        i4=rm1+rm2+rm3+1

        tridel=sqrt(fct(i1)*fct(i2)*fct(i3)/fct(i4))

c       type *,'***',rm1,rm2,rm3

        return
        end
        function tri6j(rj1,rj2,rj3)
        implicit real*8 (a-h,o-z)
        sumj=rj1+rj2+rj3
        isumj=int(sumj)
        if(abs(sumj-isumj).gt.1e-3)then !if j1+j2+j3 is not an integer
          tri6j=-1                      !then tri6j=-1
        goto 999
        endif
        if((2.*abs(rj1-rj2).le.2.*rj3).and.
     *        (2.*(rj1+rj2).ge.2.*rj3))then
          tri6j=1                               !triangle condition is satisfied
        goto 999
        endif

        tri6j=0                         !j1+j2+j3 is an integer
c       type *,rj1,rj2,rj3
999     return
        end
      
* ------------------------------------------------------------------
*    programe to calculate the 3j symbols
*-----------------------------------------------------------------
      function dr(cj1,cj2,j3,cm1,cm2,m3)
      implicit real*8(a-h,o-z)
      real*8 m3,j3
      ddr=0.d0
      if(abs(cm1+cm2+m3).gt.0.1d0)go to 80
      if(abs(cm1).gt.cj1+0.1d0)go to 80
      if(abs(cm2).gt.cj2+0.1d0)go to 80
      if(abs(cm1+cm2).gt.j3+0.1d0)go to 80
      if((cj1+cj2).lt.j3-0.1d0)go to 80
      if(abs(cj1-cj2).gt.j3+0.1d0)go to 80
      xnum=fact(j3+cj1-cj2,j3+cj2+cj1+1)*fact(j3-cj1+cj2,cj1-cm1)*
     1fact(cj1+cj2-j3,cj1+cm1)*fact(j3+m3,cj2-cm2)*
     2fact(j3-m3,cj2+cm2)
      xnum=sqrt(xnum)
      a=-cj1+cm1
      b=cj2-cj1-m3
      c=cj2+j3+cm1
      d=-cj1+cj2+j3
      e=j3-m3
      numin=0
      if(b.gt.0.d0)numin=b+0.1d0
      numax=c+0.1d0
      if(c.gt.d+0.1d0)numax=d+0.1d0
      f=numax
      if(f.gt.e+0.1d0)numax=e+0.1d0
      cj=cj1+cm2-m3
      kmi=numin+1
      kma=numax+1
      do 1 k=kmi,kma
      ck=k-1
      term=fact(1.d0,d-ck)*fact(c-ck,e-ck)*fact(ck-a,ck-b)*fact(1.d0,ck)
      ci=dmod(cj+ck,2.d0)
      if(dabs(ci).gt.0.5d0)term=-term
      ddr=ddr+term
   1  continue
      ddr=ddr*xnum
  80  dr=ddr
      return
      end
*************************************************************************
      double precision function fact(cj1,cj2)
      implicit real*8(a-h,o-z)
      n=abs(cj1-cj2)
      r=1.d0
      if(n.eq.0)go to 4
      k=1
      if(cj1.lt.cj2)k=-1
      cjx=dmax1(cj1,cj2)
      cjy=dmin1(cj1,cj2)
      hr=cjy
      do 1 i=1,n
      hr=hr+1.d0
      r=r*hr
  1   continue
      if(k.lt.0)r=1.d0/r
  4   fact=r
      return
      end
      
************************************************************************
*                                                                      *
      subroutine quad (result)
*                                                                      *
*   the argument result is an approximation  to the integral of f(r)   *
*   from  zero  to  infinity,  where the values of rp(i)*f(r(i)) are   *
*   tabulated in the array  ta(i). the integral in the interval zero   *
*   to r(2) is computed by use of an analytical fit                    *
*                                                                      *
*                                sigma                                 *
*                      f(r) = a r                                      *
*                                                                      *
*   a five-point  closed  newton-cotes  formula (cf. f b hildebrand,   *
*   introduction to numerical analysis,  edition, mcgraw-hill,   *
*   new york, 1974, p 93)  is  used  to  compute the integral in the   *
*   interval  r(2:mtp).  the  contribution  from  the  tail  of  the   *
*   function beyond the last  tabular  point  (mtp) is assumed to be   *
*   negligible. the method uses  mtp+3  tabulation points. array  ta   *
*   should therefore be dimensioned to at least  n+4 .                 *
*                                                                      *
*   no subroutines called.                                             *
*                                                                      *
*   written by farid a parpia, at oxford   last updated: 10 nov 1989   *
*                                                                      *
************************************************************************
*
      implicit doubleprecision (a-h, o-z)

      common/cons/zero,half,tenth,one,two,three,ten
     :      /def4/accy,nscf,nsic,nsolv
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
     :      /ncc/c1,c2,c3,c4

*   find first values that will permit computation of exponent
      mtpm1 = mtp-1
      do 3 i = 2,mtpm1
         tai = ta(i)
         if (abs (tai) .gt. zero) then
            ip1   = i+1
            taip1 = ta(ip1)
            quott = taip1/tai
            if (quott .gt. zero) then
*   exponent from fit
               frip1 = taip1/rp(ip1)
               fri   = tai  /rp(i  )
               ratio = frip1/fri
               rip1  = r (ip1)
               ri    = r (i  )
               sigma = log (ratio)/log (rip1/ri)
*   analytical integration and error estimate for interval r(1:i)
               fri    = ri*fri
               result = fri/(sigma+one)
*   set the tail to zero
               do 1 loc = 1,3
                  ta(mtp+loc) = zero
    1          continue
*   newton-cotes quadature for the remainder
               result = result+c1*tai
               do 2 loc = ip1,mtp,4
                  result = result+c2*(ta(loc  )+ta(loc+2))
     :                           +c3* ta(loc+1)
     :                           +c4* ta(loc+3)
    2          continue
               if (mod (mtp-i,4) .eq. 0) result = result-c1*ta(mtp)
               goto 4
            endif
         endif
    3 continue
      result = zero
    4 return
      end

************************************************************************
*                                                                      *
      subroutine yzk(k,i,j)
*                                                                      *
*   ----------------   section 09   subprogram 32   ----------------   *
*                                                                      *
*   this subroutine evaluates hartree y- and z-functions:              *
*                                                                      *
*               (k)            (k)           (k)                       *
*              y   (i,j;r) =  z   (i,j;r) + w   (i,j;r)                *
*                                                                      *
*   where                                                              *
*                                                                      *
*    (k)                                                               *
*   z   (i,j;r) =  i ( (s/r)   (p (s)*p (s) + q (s)*q (s)) ; 0 - r )   *
*                                i     j       i     j                 *
*                                                                      *
*   where                                                              *
*                                                                      *
*    (k)                    k+1                                        *
*   w   (i,j;r) =  i ( (r/s)   (p (s)*p (s) + q (s)*q (s)) ; r -       *
*                                i     j       i     j    infinity )   *
*                                                                      *
*   where  i ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
*   range  in  s .  the y-function is tabulated in  common/tatb/  in   *
*   array  tb , the z-function in array ta .                           *
*                                                                      *
*   subroutines called: draw                                           *
*                                                                      *
*   written by farid a parpia, at oxford   last updated: 13 nov 1989   *
*                                                                      *
************************************************************************
*
      implicit doubleprecision (a-h, o-z)

      dimension rhop(MN),rttk(MN),wk(MN),temp(MN),
     :          yk(MN),zk(MN)
*
      common/cnc5/cnc5c(2:5,2:4)
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /def4/accy,nscf,nsic,nsolv
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /ncc/c1,c2,c3,c4
     :      /tatb/ta(MN),tb(MN),mtp
     :      /wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :                                      mf(MNBAS)
*
      equivalence (ta(1),zk(1)),(tb(1),yk(1))
*                         k
*   for  k > 0  compute  r   and store in  rttk
      if (k .gt. 0) then
         do 1 ii = 2,n+4
            rttk(ii) = r(ii)**k
    1    continue
      endif
*   determine maximum tabulation point as location beyond which
*   rhop  (see comment statements below) would be zero; determine
*   other important locations
      mtpp1 = mtp+1
      mtpp3 = mtp+3
      mtpp4 = mtp+4
*   compute rp(s)*(p (s)*p (s)+q (s)*q (s)) and store in rhop
*                   i     j     i     j
      do 2 ii = 2,mtp
      rhop(ii) = rp(ii)*(pf(ii,i)*pf(ii,j)+qf(ii,i)*qf(ii,j))
    2 continue
*   fill array temp with r**k * rhop
      temp(1) = zero
      if (k .eq. 0) then
         do 3 ii = 2,mtp
            temp(ii) = rhop(ii)
    3    continue
      else
         do 4 ii = 2,mtp
            temp(ii) = rttk(ii)*rhop(ii)
    4    continue
      endif
*   set an additional four points to zero
      do 5 ii = mtpp1,mtpp4
         temp(ii) = zero
    5 continue
*                                     k
*   compute the first few values of  r  * zk  using semi-open
*   newton-cotes formulae
      zk(1) = zero
      do 7 ii = 2,4
         sum = zero
         do 6 kk = 2,5
            sum = sum+cnc5c(kk,ii)*temp(kk)
    6    continue
         zk(ii) = sum
    7 continue
*                         k
*   compute remainder of r  * zk: march out to mtp+3; use closed
*   newton-cotes formula
      do 8 ii = 5,mtpp3
         zk(ii) = zk(ii-4)+c1*(temp(ii-4)+temp(ii  ))
     :                    +c2*(temp(ii-3)+temp(ii-1))
     :                    +c3* temp(ii-2)
    8 continue
*                                       k   (k)
*   determine the asymptotic value of  r * z
*                   (0)
*   correction to  z   : in the manner of  c froese fischer,
*   the hartree-fock method for atoms, john wiley & sons,
*   new york, 1977, p 235.
      if (k .eq. 0) then
         if (i .eq. j) then
            zklim = one
         else
            zklim = zero
         endif
         do 10 kk = mtpp3,mtp,-1
            dif = zk(kk)-zklim
            if (abs (dif) .gt. accy) then
               do 9 ii = kk,2,-4
                  zk(ii) = zk(ii)-dif
    9          continue
            endif
   10    continue
      else
         zklim = zk(mtpp3)
      endif
*   tabulate  zk  for entire internal grid
      if (k .eq. 0) then
         do 11 ii = mtpp4,n
            zk(ii) = zklim
   11    continue
      else
         do 12 ii = 2,mtpp3
            zk(ii) = zk(ii)/rttk(ii)
   12    continue
         do 13 ii = mtpp4,n
            zk(ii) = zklim/rttk(ii)
   13    continue
      endif
*   start array wk / r**(k+1)
      np4 = n+4
      do 14 ii = np4,mtpp1,-1
         wk(ii) = zero
   14 continue
*             k+1
*   compute  r       and store in rttk
      if (k .gt. 0) then
         do 15 ii = 2,n
            rttk(ii) = rttk(ii)*r(ii)
   15    continue
      endif
*   fill array temp with rhop / r**(k+1) ; set temp(1) = zero
*   to avoid 0/0 case
      temp(1) = zero
      if (k .eq. 0) then
         do 16 ii = 2,mtp
            temp(ii) = rhop(ii)/r(ii)
   16    continue
      else
         do 17 ii = 2,mtp
            temp(ii) = rhop(ii)/rttk(ii)
   17    continue
      endif
*   compute remainder of wk / r**(k+1): march in to the origin
      do 18 ii = mtp,2,-1
         wk(ii) = wk(ii+4)+c1*(temp(ii  )+temp(ii+4))
     :                    +c2*(temp(ii+1)+temp(ii+3))
     :                    +c3*(temp(ii+2))
   18 continue
      wk(1) = zero
*   compute wk
      if (k .eq. 0) then
         do 19 ii = 2,mtp
            wk(ii) = wk(ii)*r(ii)
   19    continue
      else
         do 20 ii = 2,mtp
            wk(ii) = wk(ii)*rttk(ii)
   20    continue
      endif
*   assemble solution
      yk(1) = zero
      do 21 ii = 2,n
         yk(ii) = zk(ii)+wk(ii)
   21 continue
      return
      end
      subroutine ninej(rj11,rj12,rj13,rj21,rj22,rj23,rj31,rj32,
     : rj33,r9j)
      implicit real*8(a-h,o-z)
      common/factor/fct(0:50),mfd
      common/cons/zero,half,tenth,one,two,three,ten

      if(mfd.eq.33)goto 101
      mfd=33                          !calculate factorials if this
      fct(0)=1
      fct(1)=1                        !hasn't already been done
      do 20 i=2,mfd
      fct(i)=fct(i-1)*i
20    continue

101   tria=tri1(rj11,rj12,rj13)
      trib=tri1(rj21,rj22,rj23)
      tric=tri1(rj31,rj32,rj33)
      trid=tri1(rj11,rj21,rj31)
      trie=tri1(rj12,rj22,rj32)
      trif=tri1(rj13,rj23,rj33)

      if(tria.eq.(-one).or.trib.eq.(-one).or.tric.eq.(-one).or.
     *  trid.eq.(-one).or.trie.eq.(-one).or.trif.eq.(-one))then
      print *,'inconsistent arguments to 9j symbol'
      endif

      if(tria.eq.zero.or.trib.eq.zero.or.tric.eq.zero.or.
     *  trid.eq.zero.or.trie.eq.zero.or.trif.eq.zero)then
      r9j=0
      return
      endif

!    comment apply formula in edmonds [(6.4.3) p.101]
      rkmin=rmax3(abs(2*(rj11-rj33)),abs(2*(rj32-rj21)),
     *          abs(2*(rj12-rj23)))/2
      rkmax=rmin3(2*(rj11+rj33),2*(rj32+rj21),2*(rj12+rj23))/2
      if(rkmax.lt.rkmin)print *,'something wrong with k'
      sum=0
      rk=rkmin
99    if(rk.gt.rkmax)goto 22

      call sixj(rj11,rj21,rj31,rj32,rj33,rk,sja6)
      call sixj(rj12,rj22,rj32,rj21,rk,rj23,sjb6)
      call sixj(rj13,rj23,rj33,rk,rj11,rj12,sjc6)
      sum=sum+((-1)**(int(2*rk)))*(2*rk+1)*
     x sja6*sjb6*sjc6
      rk=rk+1
      goto 99
22    r9j=sum
      return
      end

      function tri1(rj1,rj2,rj3)
      implicit real*8 (a-h,o-z)
      sumj=rj1+rj2+rj3
      isumj=int(sumj)
      if(sumj.ne.isumj)then           !if j1+j2+j3 is not an integer
      tri1=-1                       !then tri1=-1
      goto 999
      endif
      if((2*abs(rj1-rj2).le.2*rj3).and.
     *        (2*(rj1+rj2).ge.2*rj3))then
      tri1=1                        !triangle condition is satisfied
      goto 999                      !so tri=1
      endif

      tri1=0                          !j1+j2+j3 is an integer

999   return
      end
      function rmax3(x,y,z)
      implicit real*8 (a-h,o-z)
      rmax3=x
      if(x.gt.y.and.x.gt.z)rmax3=x
      if(y.gt.x.and.y.gt.z)rmax3=y
      if(z.gt.x.and.z.gt.y)rmax3=z
      return
      end

      function rmin3(x,y,z)
      implicit real*8 (a-h,o-z)
      rmin3=x
      if(x.lt.y.and.x.lt.z)rmin3=x
      if(y.lt.x.and.y.lt.z)rmin3=y
      if(z.lt.x.and.z.lt.y)rmin3=z
      return
      end

      subroutine symmcal
      implicit real*8 (a-h,o-z)
      common/symmetry/mtbl(MNS,MNS)
c     1 stands for even and 2 stands for odd parity
c     even parity
      mtbl(1,1)=1
      mtbl(2,2)=1
      mtbl(3,3)=1
      mtbl(4,4)=1
      mtbl(5,5)=1
      mtbl(6,6)=1

      mtbl(1,3)=1
      mtbl(3,1)=1
      mtbl(3,5)=1
      mtbl(5,3)=1
      mtbl(1,5)=1
      mtbl(5,1)=1
      mtbl(2,6)=1
      mtbl(6,2)=1
      mtbl(4,6)=1
      mtbl(6,4)=1
      mtbl(2,4)=1
      mtbl(4,2)=1

c     odd parity

      mtbl(1,2)=2
      mtbl(2,1)=2

      mtbl(1,4)=2
      mtbl(4,1)=2

      mtbl(1,6)=2
      mtbl(6,1)=2

      mtbl(2,3)=2
      mtbl(3,2)=2

      mtbl(2,5)=2
      mtbl(5,2)=2

      mtbl(3,4)=2
      mtbl(4,3)=2

      mtbl(3,6)=2
      mtbl(6,3)=2

      mtbl(5,4)=2
      mtbl(4,5)=2

      mtbl(5,6)=2
      mtbl(6,5)=2

      return
      end

      subroutine setgrd
      implicit real*8 (a-h,o-z)
      common/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
      n = MN-10
      mtp=n
      h = 6.25d-02
      rnt = 2.0d-06
      r(1) = 0.
      rp(1) = rnt
      eph = exp (h)
      ett = 1.
*   set up the arrays r, rp, rpor
      do i = 2,MN
      ett = eph*ett
      ettm1 = ett-1.
      r(i) = rnt*ettm1
      rp(i) = rnt*ett
      rpor(i) = ett/ettm1
      enddo
      return
      end
      
      subroutine daopen
      implicit real*8 (a-h,o-z)

      open(STDIN,file='input.dbar',form='formatted',status='old')
      open(STDOUT,file='dbar.out',form='formatted',status='unknown')
      open(WFNIN,file='wfn.dat',form='unformatted',status='old')
      open(NTFILE,file='to.dat',form='unformatted',status='old')
      return
      end

c**********************************************************************

c***********************************************************************
c     This subroutine constructs the inpute file of the property       *
c     code                                                             *
c***********************************************************************


      subroutine readinp(nsym,nbas,iocc,iopt,nprint,idip,ivel,iquad,
     :ihyp,iquadB,ioctC,im1p,im2p,igfact,iovlp,insi,insd,inms)

      implicit real*8(a-h,o-z)
      common/iallparity/iall
      common/inparm/nparm
      common/pncp/amass,qw
      dimension nbas(MNSYM),iocc(MNSYM)  
    
c----------------------------------------------------------------------
c     read number of symmetry
c----------------------------------------------------------------------

      read(STDIN,*)nsym,(nbas(i),i=1,nsym)

c----------------------------------------------------------------------
c     occupancy of each orbital
c----------------------------------------------------------------------

      read(STDIN,*)(iocc(i),i=1,nsym)

c----------------------------------------------------------------------
c     truncation scheme and option flag
c----------------------------------------------------------------------

      read(STDIN,*)iopt,nprint,iall

c----------------------------------------------------------------------
c     nuclear property
c----------------------------------------------------------------------

      read(STDIN,*)nparm,qw,amass

c----------------------------------------------------------------------
c     property related optional flags
c---------------------------------------------------------------------

      read(STDIN,*)idip,ivel,iquad,ihyp,iquadB,ioctC,im1p,im2p,
     :insi,insd,igfact,iovlp
      read(STDIN,*)inms

      return
      end

c*********************************************************************

      subroutine grasprd(nbasis)
      implicit real*8 (a-h,o-z)
     
      common/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :       mf(MNBAS)
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
      common/orbital_energy/eorb(MNBAS)
c     local dimension
     
      rewind (WFNIN)
      read(WFNIN)h,n
      mtp=n
      read(WFNIN)(r(i),i = 1,n),(rp(i),i = 1,n),(rpor(i),i=1,n)
      do i=1,nbasis
      read(WFNIN)eorb(i)
      mf(i)=n
      mfj=mf(i)
      read(WFNIN)(pf(j,i),j=1,mfj),(qf(j,i),j=1,mfj)
      enddo
      rewind WFNIN

      write(STDOUT,9999)
      write(STDOUT,9998)(eorb(i),i=1,nbasis)

      mtp=n

9999  format(//,10x,17h ORBITAL ENERGIES,/)
9998  format(4d20.11)
      return
      end
     
      subroutine setup(nsym,nbasis,nbas,iocc)
      implicit real*8 (a-h,o-z)
      character*1 orbtyp,orbsym
      character*4 qj1,qj
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)
      common/maxj/jmax
      common/kappa/nka(MNBAS)

c     local dimension
      dimension kap(MNSYM),ia(MNSYM),iorbsm(MNSYM)
      dimension orbtyp(MNS),nkj(MNBAS)
      dimension orb(MNSYM),np(MNBAS)
      dimension qj(MNBAS),qj1(MNS)
      dimension nprimc(MNOCC),nkprimc(MNOCC)
      dimension nprime(MNEXC),nkprime(MNEXC)
      dimension nmc(MNBAS),nme(MNBAS),nm(MNBAS)
      dimension nbas(MNSYM),iocc(MNSYM)
      dimension jorbsm(MNSYM)
      data kap/-1, 1,-2, 2,-3, 3,-4, 4,-5,5,-6/
      data orb/0.5,0.5,1.5,1.5,2.5,2.5,3.5,3.5,4.5,4.5,5.5/
      data ia /1,-1,1,-1,1,-1,1,-1,1,-1,1/
      data iorbsm/1,2,2,3,3,4,4,5,5,6,6/
      data jorbsm/1,2,2,1,1,2,2,1,1,2,2/
      data qj1/'1/2','3/2','5/2','7/2','9/2','11/2'/
      data orbtyp/'s','p','d','f','g','h'/

      nocc=0
      nexcit=0
      nbasis=0
      jmax=2*orb(nsym)-1
      ii=0
      do isym=1,nsym
      ibasis=nbas(isym)
      if(ibasis.ne.0) then
      nbasis = nbasis + ibasis
      do jbas=1,ibasis
      ii = ii + 1
      nka(ii)=kap(isym)
      orbj(ii)=orb(isym)
      np(ii)=jbas
      iparity(ii)=jorbsm(isym)
      if(nka(ii).eq.kap(isym))iiq(ii)=ia(isym)
      if(nka(ii).eq.kap(isym).and.jbas.le.iocc(isym))then
      nocc=nocc+1
      kc(nocc)=ii
      nkprimc(nocc)=kap(isym)
      nprimc(nocc)=np(ii)
      orbc(nocc)=orbj(ii)
      nmc(nocc)=idint(2*(orbc(nocc))) +1
      if(nkprimc(nocc).eq.kap(isym))isymc(nocc)=iorbsm(isym)
      if(nkprimc(nocc).eq.kap(isym))iqc(nocc)=ia(isym)
      else
      nexcit=nexcit+1
      ke(nexcit)=ii
      nprime(nexcit)=np(ii)
      nkprime(nexcit)=nka(ii)
      orbe(nexcit)=orb(isym)
      nme(nexcit)=idint(2*(orbj(ii))) +1
      if(nkprime(nexcit).eq.kap(isym))isyme(nexcit)=iorbsm(isym)
      if(nkprime(nexcit).eq.kap(isym))iqe(nexcit)=ia(isym)
      nup=nup+1
      endif
      enddo
      endif
      enddo

c     write occupied and unoccupied orbitals
      if(iproc.eq.0)then
      write(STDOUT,101)
      write(STDOUT,100)
      endif
      do i=1,nocc
      if(orbc(i).eq.0.5)qj(i)=qj1(1)
      if(orbc(i).eq.1.5)qj(i)=qj1(2)
      if(orbc(i).eq.2.5)qj(i)=qj1(3)
      if(orbc(i).eq.3.5)qj(i)=qj1(4)
      if(orbc(i).eq.4.5)qj(i)=qj1(5)
      if(orbc(i).eq.5.5)qj(i)=qj1(6)
      if(isymc(i).eq.1)orbsym=orbtyp(1)
      if(isymc(i).eq.2)orbsym=orbtyp(2)
      if(isymc(i).eq.3)orbsym=orbtyp(3)
      if(isymc(i).eq.4)orbsym=orbtyp(4)
      if(isymc(i).eq.5)orbsym=orbtyp(5)
      if(isymc(i).eq.6)orbsym=orbtyp(6)
      if(isymc(i).eq.3)isymc(i)=1
      if(isymc(i).eq.4)isymc(i)=2
      if(isymc(i).eq.5)isymc(i)=1
      if(isymc(i).eq.6)isymc(i)=2
      write(STDOUT,200)orbsym,nprimc(i),nmc(i),kc(i),qj(i)
      enddo

      write(STDOUT,102)
      write(STDOUT,100)
      do i=1,nexcit
      if(orbe(i).eq.0.5)qj(i)=qj1(1)
      if(orbe(i).eq.1.5)qj(i)=qj1(2)
      if(orbe(i).eq.2.5)qj(i)=qj1(3)
      if(orbe(i).eq.3.5)qj(i)=qj1(4)
      if(orbe(i).eq.4.5)qj(i)=qj1(5)
      if(orbe(i).eq.5.5)qj(i)=qj1(6)
      if(isyme(i).eq.1)orbsym=orbtyp(1)
      if(isyme(i).eq.2)orbsym=orbtyp(2)
      if(isyme(i).eq.3)orbsym=orbtyp(3)
      if(isyme(i).eq.4)orbsym=orbtyp(4)
      if(isyme(i).eq.5)orbsym=orbtyp(5)
      if(isyme(i).eq.6)orbsym=orbtyp(6)
      if(isyme(i).eq.3)isyme(i)=1
      if(isyme(i).eq.4)isyme(i)=2
      if(isyme(i).eq.5)isyme(i)=1
      if(isyme(i).eq.6)isyme(i)=2
      write(STDOUT,200)orbsym,nprime(i),nme(i),ke(i),qj(i)
      enddo
 100  format(/,1x,'TYPE',2x,' N ',2x,'2J+1',2x,'ORB. NO.',
     $3x,' J ',/)
 101  format(/,10x,' OCCUPIED ORBITALS',/)
 102  format(/,10x,' EXCITED ORBITALS',/)
 200  format(2x,a,2x,3(i4,2x),6x,a3)
      return
      end



      subroutine setqic
      implicit real*8 (a-h,o-z)
      logical first
*
      dimension b13(13,13),cg(6),c5num(1:5,2:5),c6num(1:6,2:6)
*
      common/cnc5/cnc5c(2:5,2:4)
     :      /cons/zero,half,tenth,one,two,three,ten
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /cnc6/cnc6c(1:6,2:6)
     :      /lic13/a(13,13)
     :      /ncc/c1,c2,c3,c4
     :      /sbc/c(6)
*
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*   THIRTEEN-POINT  LAGRANGE  INTERPOLATION  COEFFICIENTS FOR FIRST    *
*   DERIVATIVE                                                         *
*
      DATA (B13( 1,I),I = 1,13) /            -1486442880.0D 00,
     :            5748019200.0D 00, -15807052800.0D 00,
     :           35126784000.0D 00, -59276448000.0D 00,
     :           75873853440.0D 00, -73766246400.0D 00,
     :           54195609600.0D 00, -29638224000.0D 00,
     :           11708928000.0D 00,  -3161410560.0D 00,
     :             522547200.0D 00,    -39916800.0D 00/
      DATA (B13( 2,I),I = 1,13) /              -39916800.0D 00,
     :            -967524480.0D 00,   2634508800.0D 00,
     :           -4390848000.0D 00,   6586272000.0D 00,
     :           -7903526400.0D 00,   7376624640.0D 00,
     :           -5269017600.0D 00,   2822688000.0D 00,
     :           -1097712000.0D 00,    292723200.0D 00,
     :             -47900160.0D 00,      3628800.0D 00/
      DATA (B13( 3,I),I = 1,13) /                3628800.0D 00,
     :             -87091200.0D 00,   -684478080.0D 00,
     :            1596672000.0D 00,  -1796256000.0D 00,
     :            1916006400.0D 00,  -1676505600.0D 00,
     :            1149603840.0D 00,   -598752000.0D 00,
     :             228096000.0D 00,    -59875200.0D 00,
     :               9676800.0D 00,      -725760.0D 00/
      DATA (B13( 4,I),I = 1,13) /                -725760.0D 00,
     :              13063680.0D 00,   -143700480.0D 00,
     :            -476910720.0D 00,   1077753600.0D 00,
     :            -862202880.0D 00,    670602240.0D 00,
     :            -431101440.0D 00,    215550720.0D 00,
     :             -79833600.0D 00,     20528640.0D 00,
     :              -3265920.0D 00,       241920.0D 00/
      DATA (B13( 5,I),I = 1,13) /                 241920.0D 00,
     :              -3870720.0D 00,     31933440.0D 00,
     :            -212889600.0D 00,   -303937920.0D 00,
     :             766402560.0D 00,   -447068160.0D 00,
     :             255467520.0D 00,   -119750400.0D 00,
     :              42577920.0D 00,    -10644480.0D 00,
     :               1658880.0D 00,      -120960.0D 00/
      DATA (B13( 6,I),I = 1,13) /                -120960.0D 00,
     :               1814400.0D 00,    -13305600.0D 00,
     :              66528000.0D 00,   -299376000.0D 00,
     :            -148262400.0D 00,    558835200.0D 00,
     :            -239500800.0D 00,     99792000.0D 00,
     :             -33264000.0D 00,      7983360.0D 00,
     :              -1209600.0D 00,        86400.0D 00/
      DATA (B13( 7,I),I = 1,13) /                  86400.0D 00,
     :              -1244160.0D 00,      8553600.0D 00,
     :             -38016000.0D 00,    128304000.0D 00,
     :            -410572800.0D 00,            0.0D 00,
     :             410572800.0D 00,   -128304000.0D 00,
     :              38016000.0D 00,     -8553600.0D 00,
     :               1244160.0D 00,       -86400.0D 00/
      DATA (B13( 8,I),I = 1,13) /                 -86400.0D 00,
     :               1209600.0D 00,     -7983360.0D 00,
     :              33264000.0D 00,    -99792000.0D 00,
     :             239500800.0D 00,   -558835200.0D 00,
     :             148262400.0D 00,    299376000.0D 00,
     :             -66528000.0D 00,     13305600.0D 00,
     :              -1814400.0D 00,       120960.0D 00/
      DATA (B13( 9,I),I = 1,13) /                 120960.0D 00,
     :              -1658880.0D 00,     10644480.0D 00,
     :             -42577920.0D 00,    119750400.0D 00,
     :            -255467520.0D 00,    447068160.0D 00,
     :            -766402560.0D 00,    303937920.0D 00,
     :             212889600.0D 00,    -31933440.0D 00,
     :               3870720.0D 00,      -241920.0D 00/
      DATA (B13(10,I),I = 1,13) /                -241920.0D 00,
     :               3265920.0D 00,    -20528640.0D 00,
     :              79833600.0D 00,   -215550720.0D 00,
     :             431101440.0D 00,   -670602240.0D 00,
     :             862202880.0D 00,  -1077753600.0D 00,
     :             476910720.0D 00,    143700480.0D 00,
     :             -13063680.0D 00,       725760.0D 00/
      DATA (B13(11,I),I = 1,13) /                 725760.0D 00,
     :              -9676800.0D 00,     59875200.0D 00,
     :            -228096000.0D 00,    598752000.0D 00,
     :           -1149603840.0D 00,   1676505600.0D 00,
     :           -1916006400.0D 00,   1796256000.0D 00,
     :           -1596672000.0D 00,    684478080.0D 00,
     :              87091200.0D 00,     -3628800.0D 00/
      DATA (B13(12,I),I = 1,13) /               -3628800.0D 00,
     :              47900160.0D 00,   -292723200.0D 00,
     :            1097712000.0D 00,  -2822688000.0D 00,
     :            5269017600.0D 00,  -7376624640.0D 00,
     :            7903526400.0D 00,  -6586272000.0D 00,
     :            4390848000.0D 00,  -2634508800.0D 00,
     :             967524480.0D 00,     39916800.0D 00/
      DATA (B13(13,I),I = 1,13) /               39916800.0D 00,
     :            -522547200.0D 00,   3161410560.0D 00,
     :          -11708928000.0D 00,  29638224000.0D 00,
     :          -54195609600.0D 00,  73766246400.0D 00,
     :          -75873853440.0D 00,  59276448000.0D 00,
     :          -35126784000.0D 00,  15807052800.0D 00,
     :           -5748019200.0D 00,   1486442880.0D 00/
*
      DATA B13DEN/479001600.0D 00/
*
*----------------------------------------------------------------------*
*   COEFFICIENTS FOR SIENKIEWICZ-BAYLIS ALGORITHM                      *
*
      DATA CG / 1771.0D 00,
     :          9235.0D 00,
     :          5890.0D 00,
     :          4610.0D 00,
     :            35.0D 00,
     :            59.0D 00/
*
      DATA DENOM /5760.0D 00/
*
*----------------------------------------------------------------------*
*                                                                      *
*   FIVE-POINT NEWTON-COTES COEFFICIENTS FOR CLOSED INTEGRATION. EX-   *
*   PRESSED AS RATIONAL NUMBERS                                        *
*
      DATA (C5NUM(I,2),I = 1,5)/ 251.0D 00, 646.0D 00,
     :       -264.0D 00, 106.0D 00, -19.0D 00/
      DATA (C5NUM(I,3),I = 1,5)/ 232.0D 00, 992.0D 00,
     :        192.0D 00,  32.0D 00,  -8.0D 00/
      DATA (C5NUM(I,4),I = 1,5)/ 243.0D 00, 918.0D 00,
     :        648.0D 00, 378.0D 00, -27.0D 00/
      DATA (C5NUM(I,5),I = 1,5)/ 224.0D 00,1024.0D 00,
     :        384.0D 00,1024.0D 00, 224.0D 00/
*
      DATA C5DEN/ 720.0D 00/
*
*----------------------------------------------------------------------*
*                                                                      *
*   SIX-POINT NEWTON-COTES COEFFICIENTS FOR CLOSED INTEGRATION.  EX-   *
*   PRESSED AS RATIONAL NUMBERS                                        *
*
      DATA (C6NUM(I,2),I = 1,6)/ 475.0D 00,1427.0D 00,
     :                          -798.0D 00, 482.0D 00,
     :                          -173.0D 00,  27.0D 00/
      DATA (C6NUM(I,3),I = 1,6)/ 448.0D 00,2064.0D 00,
     :                           224.0D 00, 224.0D 00,
     :                           -96.0D 00,  16.0D 00/
      DATA (C6NUM(I,4),I = 1,6)/ 459.0D 00,1971.0D 00,
     :                          1026.0D 00,1026.0D 00,
     :                          -189.0D 00,  27.0D 00/
      DATA (C6NUM(I,5),I = 1,6)/ 448.0D 00,2048.0D 00,
     :                           768.0D 00,2048.0D 00,
     :                           448.0D 00,   0.0D 00/
      DATA (C6NUM(I,6),I = 1,6)/ 475.0D 00,1875.0D 00,
     :                          1250.0D 00,1250.0D 00,
     :                          1875.0D 00, 475.0D 00/
*
      DATA C6DEN/1440.0D 00/
*
*----------------------------------------------------------------------*
*
c      DATA FIRST /.TRUE./
*

*   LAGRANGE INTERPOLATION COEFFICIENTS
*
*   DO THIS INITIALIZATION ONCE PER RUN ONLY
*
c      IF (FIRST) THEN
*
*   THIRTEEN-POINT COEFFICIENTS FOR DPBDT
*
         FACTOR = ONE/B13DEN
         DO 2 J = 1,13
            DO 1 I = 1,13
               a(I,J) = B13(I,J)*FACTOR
    1       CONTINUE
    2    CONTINUE
*
c         FIRST = .FALSE.
*
c      ENDIF
*
*   SIENKIEWICZ-BAYLIS COEFFICIENTS FOR SBSTEP
*
      C(1) = CG(1)/DENOM
      FACTOR = H/DENOM
      DO 3 I = 2,6
         C(I) = CG(I)*FACTOR
    3 CONTINUE
*
*   NEWTON-COTES COEFFICIENTS FOR YZK AND QUAD
*
      FACTOR = H/C5DEN
      DO 5 J = 2,4
         DO 4 I = 2,5
            CNC5C(I,J) = FACTOR*C5NUM(I,J)
    4    CONTINUE
    5 CONTINUE
*
      C1 = FACTOR*C5NUM(1,5)
      C2 = FACTOR*C5NUM(2,5)
      C3 = FACTOR*C5NUM(3,5)
      C4 = C1+C1
*
*   NEWTON-COTES COEFFICIENTS FOR START
*
      FACTOR = H/C6DEN
      DO 7 J = 2,6
         DO 6 I = 1,6
            CNC6C(I,J) = FACTOR*C6NUM(I,J)
    6    CONTINUE
    7 CONTINUE
 
C     NUMERICAL CONSTANTS

      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      THREE=3.0D0
      TEN=10.0D0
      TENTH=ONE/TEN
      HALF=ONE/TWO

      RETURN
      END

      subroutine daclose
      close (STDIN)
      close (STDOUT)
      close (WFNIN)
      return
      end


c***************************************************************************
c                                                                          *
c      SUBROUTINE   SYMM                                                   *
c                                                                          *
c     THIS SUBROUTINE SETS UP THE T1 AND T2 EQUATION INDEX                 *
c                                                                          *
c***************************************************************************

      subroutine symm(iopt,nprint,ti,tf)

      implicit real*8 (a-h,o-z)
      common/iallparity/iall
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/maxj/jmax
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/skip/nnskip(0:MXV,2),n2skip,nsing
      dimension ncheck(MDIM)

c     local dimension
      dimension kgot(MXV+1),ieven(0:MXV),iodd(0:MXV)
      dimension indx(MDIM),indx1(MDIM)
      dimension ti(MDIM),tf(MDIM)

      nsing=0

      do i=1, MDIM
         ncheck(i)=0
      enddo

c---------------------------------------------------------------------------
c     set the equation index for T1
c---------------------------------------------------------------------------

      if(iopt.ne.0)then
      write(STDOUT,9999)
      write(STDOUT,9997)

      do 1 ia=1,nocc
      iasym=isymc(ia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)

      iapsym=mtbl(ipsym,iasym)

      orba=orbc(ia)
      orbp=orbe(ip)

      iaa=iqc(ia)
      iap=iqe(ip)

      iia=kc(ia)
      iip=ke(ip)

      if(iapsym.ne.1)go to 1
      if(orba.ne.orbp)go to 1

      nsing=nsing+1

      idpa(iip,iia)=nsing

      indx(nsing)=nsing
      indx1(nsing)=nsing
      
      tf(nsing)=ti(nsing)

   1  continue

      endif

c----------------------------------------------------------------------------
c    set the skip information
c----------------------------------------------------------------------------

      do 2 i=0,MXV
      ieven(i)=0
      iodd(i)=0
      nnskip(i,1)=0
      nnskip(i,2)=0
  2   continue

      do 3 ia=1,nocc
      iasym=isymc(ia)

      do 3 ip=1,nexcit
      ipsym=isyme(ip)

      iapsym=mtbl(ipsym,iasym)

      orba=orbc(ia)
      orbp=orbe(ip)

      iaa=iqc(ia)
      iap=iqe(ip)

      iia=kc(ia)
      iip=ke(ip)

      call findk(orba,orbp,orba,orbp,iaa,iap,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 3
      do 4 nloop=1,kmax
      kk=kgot(nloop)
      if(iapsym.eq.1)then
           ieven(kk)=ieven(kk)+1
           idra(iip,iia,kk)=ieven(kk)
      else
           iodd(kk)=iodd(kk)+1
           idra(iip,iia,kk)=iodd(kk)
      endif

  4   continue
  3   continue

      n2skip=0

      if(iall.eq.0)then

      do 5 i=2,MXV,2
      nnskip(i,2)=nnskip(i-2,2)+ieven(i-2)*(ieven(i-2)+1)/2
  5   continue

      do  i=0,JMAX,2
      n2skip=n2skip+ieven(i)*(ieven(i)+1)/2
      enddo

      do 6 i=3,MXV-2,2
      nnskip(i,1)=nnskip(i-2,1)+iodd(i-2)*(iodd(i-2)+1)/2
  6   continue

       else

      do 9 i=1,MXV
      nnskip(i,2)=nnskip(i-1,2)+ieven(i-1)*(ieven(i-1)+1)/2
  9   continue

      do  i=0,MXV
      n2skip=n2skip+ieven(i)*(ieven(i)+1)/2
      enddo

      do 10 i=1,MXV
      nnskip(i,1)=nnskip(i-1,1)+iodd(i-1)*(iodd(i-1)+1)/2
 10   continue

      endif

      nn=nsing

      write(STDOUT,9998)
      write(STDOUT,9996)

c--------------------------------------------------------------------------
c    set the equation index for T2
c--------------------------------------------------------------------------

      do 15 ia=1,nocc
      iasym=isymc(ia)

      do 20 ip=1,nexcit
      ipsym=isyme(ip)

      iapsym=mtbl(ipsym,iasym)

      do 30 ib=1,nocc
      ibsym=isymc(ib)

      do 40 iq=1,nexcit
      iqsym=isyme(iq)

      ibqsym=mtbl(ibsym,iqsym)

      if(ibqsym.ne.iapsym)go to 40

      orba=orbc(ia)
      orbb=orbc(ib)
      orbp=orbe(ip)
      orbq=orbe(iq)

      iaa=iqc(ia)
      iab=iqc(ib)
      iap=iqe(ip)
      iaq=iqe(iq)

      iia=kc(ia)
      iib=kc(ib)
      iip=ke(ip)
      iiq=ke(iq)

      call findk(orba,orbp,orbb,orbq,iaa,iap,iab,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 45 nloop=1,kmax
      kk=kgot(nloop)
      nap=idra(iip,iia,kk)
      nbq=idra(iiq,iib,kk)
      call iloc(iip,iia,iiq,iib,ipa,iqb)
      if(ipa.lt.iqb)go to 45
      nn=nn+1
      neqn=ntloc(nap,nbq,iapsym,kk)

      indx(nn)=nn
      indx1(nn)=neqn

      if (ncheck(neqn).eq.1) then
         write(STDOUT,*)'double count in symm =', neqn
      stop
      endif

      ncheck(neqn)=1

      tf(neqn)=ti(neqn)*(-1)**(orbp+orba+orbq+orbb)

  45  continue
  40  continue
  30  continue
  20  continue
  15  continue

      ntmax=nn

9999  format(//,'SINGLE EXCITATIONS',/)
9998  format(//,'DOUBLE EXCITATIONS',/)
9997  format(3x,'IN','    IA ','  IP ','J(IA)',' J(IP)',/)
9996  format(3x,'IN','    IA ','  IB ','  IP ',' IQ' ,3x,'JA',
     $         4x,'JB',5x,'JP',5x,'JQ',5x,'ASM',2x,'BSM',2x,
     $         'PSM',2x,'QSM',2x,'KMAX',/)
9995  format(3i5,2(f5.1,2x))
9994  format(5i5,4(f5.1,2x),8(i4,1x))

      return
      end

c***********************************************************************

      function ntloc(iap,ibq,iapsym,kk)
      implicit real*8 (a-h,o-z)
c     common block
      common/skip/nnskip(0:MXV,2),n2skip,nsing
      mm=max0(iap,ibq)
      nn=min0(iap,ibq)
      if(iapsym.eq.1)then
      ntloc=mm*(mm-1)/2+nn+nnskip(kk,2)+nsing
      else
      ntloc=mm*(mm-1)/2+nn+nnskip(kk,1)+n2skip+nsing
      endif
      return
      end

c********************************************************************
  
      subroutine cinit(d)
      real*8 d(MNBAS,MNBAS)
      do 1 i=1,MNBAS
      do 1 j=1,MNBAS
 1    d(i,j)=0.0d0
      return
      end

c*****************************************************************
c     This subroutine calculated the hyperfine O bar diagram of  *
c     of the type  HH                                            *
c*************************************************************** *

      subroutine hyperhh(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)

      if(ibsym.ne.iasym)go to 1

c------------------------------------------------------------------
c     diagram 1  = <b|d|a>
c------------------------------------------------------------------

      call hyper(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      d(iib,iia)=d(iib,iia)+dipole

      sum=zero

      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

c----------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c----------------------------------------------------------------

      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.eq.ibsym)then

      ini=idpa(iip,iia)

      call hyper(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)

      endif
      endif
      endif

c----------------------------------------------------------------
c     diagram 3 =  <b|t1|p><p|d|a>
c----------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.eq.iasym)then

      inf=idpa(iip,iib)

      call hyper(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)

      endif
      endif
      endif

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.eq.ipbsym) then

c----------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c----------------------------------------------------------------

      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.iasym)then

      call hyper(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      if(orbb.eq.orbp)then
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 5 =  <bg|t2|rp><p|d|a><r|t1|g>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.eq.igsym)then

      call hyper(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c--------------------------------------------------------------------

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npb=idra(iip,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(npb,nrg,irgsym,1)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      ini=idpa(iip,iia)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then

c------------------------------------------------------------------
c     diagrams 8 = <g|d|r><b|t1|p><pr|t2|ag>
c------------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then

      call hyper(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      inf=idpa(iip,iib)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.ibsym)then

      call hyper(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)

      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      if(orbp.eq.orba)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
      if(orba.eq.orbp)then
c     sum=sum-dipole*tf(inf)*tf(inf)*fact
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(ipsym.ne.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c-------------------------------------------------------------------
c     diagram 12 = <b|t1|r><t|d|p><p|t1|a>
c-------------------------------------------------------------------

      call hyper(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

      inf=idpa(iir,iib)
      ini=idpa(iip,iia)

      sum=sum+dipole*tf(inf)*ti(ini)

 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

c-----------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c-----------------------------------------------------------------

      if(igsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call hyper(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ini=idpa(iip,iia)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c-------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c-------------------------------------------------------------------

      if(igsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call hyper(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iip,iib)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 8

      call hyper(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+1+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 10   continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,one,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11   continue
  9   continue
  8   continue

      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(iqsym.ne.irsym)go to 13

      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13

      call hyper(orbq,orbr,iiq,iir,orbq,orbr,kapq,kapr,dipole)

c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,one,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue

c-------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(one,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+one+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo

 14   continue
 13   continue

      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.ne.irsym)go to 16  
       
      call hyper(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)


c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      call sixj(one,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+one)
      fact=s6j*isign/(2*al1+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,one,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+one)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact

  18  continue
  17  continue
  16  continue
  12  continue

      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(idsym.ne.ibsym)go to 21

      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22

      call hyper(orbb,orbd,iib,iid,orbb,orbd,kapb,kapd,dipole)

c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

c----------------------------------------------------------------- 
c     diagram 22= <gd|t2|pr><b|d|d><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 24   continue
 23   continue
 22   continue
 21   continue    
  
      if(idsym.ne.iasym)go to 25

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26

      call hyper(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(npd,nrg,irgsym,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------------

      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

 28   continue
 27   continue
 26   continue
 25   continue      
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+1)

 1    continue   
  
      return
      end

c*****************************************************************
c     This subroutine calculated the hyperfine O bar diagram of  *
c     of the type  HP                                            *
c*************************************************************** *

      subroutine hyperhp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)

      if(ipsym.ne.iasym)go to 1

c-------------------------------------------------------------------
c     diagram 1  = <p|d|a>
c-------------------------------------------------------------------

      call hyper(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      d(iip,iia)=d(iip,iia)+dipole

      sum=zero
      
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 3

      call hyper(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.irgsym)go to 401

c---------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c---------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,ipasym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(2*one+1)
      sum=sum+dipole*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then   
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dipole*ti(ini)*fact

 4    continue

      endif

 401  continue

c--------------------------------------------------------------------------
c     diagram 4 = <r|t1|a><r|d|g><p|t1|g>
c--------------------------------------------------------------------------

      if(igsym.ne.irsym)go to 3
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3

      call hyper(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ini=idpa(iir,iia)
      inj=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*ti(inj)

 3    continue
  
c-------------------------------------------------------------------------
c     diagram 5 =  <r|t1|a><g|t1|r><p|d|g>
c-------------------------------------------------------------------------

      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5

      call hyper(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dipole)

      ini=idpa(iir,iia)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 5    continue

c-------------------------------------------------------------------
c     diagram 6 = <r|d|a><g|t|r><p|t|g>
c-------------------------------------------------------------------

      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6

      call hyper(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 6    continue

      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7

      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym)go to 8

      call hyper(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.iqasym)go to 8

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 7 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      if(orba.eq.orbq)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))

      endif
      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 8 =  <g|t1|r><p|d|q><qr|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 9    continue
 8    continue


      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      if(iasym.ne.ibsym)go to 10

      call hyper(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)

      if(irgsym.ne.ipbsym)go to 10

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 9 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbp.eq.orbb)then
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 10 = <g|t1|r><b|d|a><pr|t2|gb>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue
  7   continue


      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqasym) goto 13
      if(ipsym.ne.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13

      call hyper(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dipole)

c-----------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orba,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue
 14   continue
 13   continue


      if(orbp.ne.orbq)go to 12

      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12

      call hyper(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dipole)

c--------------------------------------------------------------------
c     diagram 13 =  <gb|t2|rq><pr|t2|bg><q|d|a>
c--------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 17   continue
 16   continue
 12   continue


      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2

      inf=idpa(iir,iig)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)

      if(irsym.ne.iqsym)go to 20

      call hyper(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipasym.ne.iqgsym)go to 20

c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------
  
      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(one,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 21   continue
 20   continue

      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)

      if(ibsym.ne.igsym)go to 22

      call hyper(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)

      if(ipasym.ne.irbsym)go to 22

c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nrb=idra(iir,iib,1)
      ini=ntloc(npa,nrb,ipasym,1)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*1+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(one,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 23   continue
 22   continue
  2   continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------

      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24

      ini=idpa(iiq,iia)

      call hyper(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      sum=sum+ti(ini)*dipole

 24   continue

c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------

      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25

      ini=idpa(iip,iib)

      call hyper(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      sum=sum-ti(ini)*dipole

 25   continue

      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 26

      call hyper(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*1+1)*(2*1+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

 27   continue

c------------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg> 
c------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)

      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(one,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*1+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 29   continue

      endif
      endif
      enddo
      endif

 28   continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(one,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact

 31   continue

      endif
      endif
      enddo
      endif

 30   continue

c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32

      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(one,orbr,orbg,al1,orbq,orbb,s6j1)
      call sixj(one,orba,orbp,al2,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 34   continue
 33   continue
 32   continue
 26   continue

 
      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+1)

 1    continue

      return
      end

c*****************************************************************
c     This subroutine calculates the hyperfine O bar diagram of  *
c     of the type  PP                                            *
c*************************************************************** *
     
      subroutine hyperpp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
    
      if(ipsym.ne.iqsym) go to 1

      call hyper(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

c-----------------------------------------------------------------
c     diagram 1 = <p|d|q>
c-----------------------------------------------------------------

      d(iip,iiq)=d(iip,iiq)+dipole

      sum=0.0d0
        
      sum6=0
      sum8=0
      
      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      if(iasym.eq.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then

      call hyper(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      ini=idpa(iip,iia)

c----------------------------------------------------------------
c     diagram 2 = <p|t|a><a|d|q>
c----------------------------------------------------------------

      sum=sum-dipole*ti(ini)

      endif
      endif
      endif

      if(iasym.eq.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then

      call hyper(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

      sum=sum-dipole*tf(inf)

      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.eq.ipsym)then

      call hyper(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orbq.eq.orba)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c-----------------------------------------------------------------

      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)then
      sum=sum+dipole*ti(ini)*tf(inf)*fact
 
      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.eq.irsym)then

      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then

      call hyper(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      nqa=idra(iiq,iia,1)
      inf=ntloc(nrg,nqa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c----------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.eq.irsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call hyper(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------  
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      npa=idra(iip,iia,1)
      ini=ntloc(nrg,npa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.eq.iqsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call hyper(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orba.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)then
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 
      endif
      enddo
      endif
      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

c-----------------------------------------------------------------
c     diagram 12 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

      if(irsym.eq.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then

      call hyper(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------
c     diagram 13 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------

      if(irsym.eq.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then

      call hyper(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  5   continue

c------------------------------------------------------------------
c     diagram 14 = <q|t1|a><a|d|g><g|t1|p
c------------------------------------------------------------------

      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      if(igsym.eq.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then

      call hyper(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)

      sum=sum+dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue     
 4    continue 
    
      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 9

      call hyper(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue

c-----------------------------------------------------------------
c     diagram 16 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)

      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,one,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 14   continue
 13   continue
 12   continue
  9   continue

      if(iasym.ne.idsym)go to 8

      call hyper(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)

      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8

c-----------------------------------------------------------------
c     digram 17 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(one,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+one+l1)/(2*al1+1)   
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+one+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,one,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*ti(ini)*tf(inf)*fact

 16   continue
 15   continue
  8   continue

c------------------------------------------------------------------
c     diagram 19 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------

      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nka(iis)

      if(issym.ne.irsym)go to 18

      call hyper(orbs,orbr,iis,iir,orbs,orbr,kaps,kapr,dipole)

      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18

      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,one,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2      
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 20   continue

c------------------------------------------------------------------
c     diagram 20 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,one,orbg,al2,orbp,s9j)
c      fact=(-1)**(orbq+orbg+l2+one+2*orbp)*s9j           
      fact=(-1)**(orba+orbg+orbr+orbs+l1+l2)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 21   continue
 19   continue
 18   continue

c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)

      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.ne.issym)go to 22

      call hyper(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,dipole)

      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)   
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 24   continue
 23   continue
 22   continue

c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)

      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.ne.issym)go to 25

      call hyper(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,dipole)

      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax 
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 27   continue
 26   continue
 25   continue
 17   continue
  7   continue
 
      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+1)

  1   continue

      return
      end

c*********************************************************
 

c*********************************************************
c  This subroutine calculates the reduced                *
c  matrix element  of the hyperfine operator             *
c*********************************************************

      subroutine hyper(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp)

      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2

      common/cons/zero,half,tenth,one,two,three,ten
      common/inthyp/rint(MNBAS,MNBAS)

      factor=(-1)**(j1+0.5)
      w3j1=dr(j1,one,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      hyp=-result*e1*(kap1+kap2)

      return
      end

c********************************************************************

c********************************************************************
c                                                                   *
c    This subroutine calculates the radial part of the reduced      *
c    matrix element of the hyperfine oparator                       *
c                                                                   *
c********************************************************************

      subroutine Calchyp(nbasis)

      implicit real*8 (a-h,o-z)

      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/inthyp/rint(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      mtp=n


      do i1 = 1, nbasis
      do i2 = 1, i1

      rint(i1,i2) = zero

      if(iparity(i1).eq.iparity(i2))then

      ta(1)=zero

      do 1 l=2,mtp
      ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)/(r(l)*r(l))
    1 continue
      call quad(result)

      endif

      rint(i1,i2) = result
      rint(i2,i1) = rint(i1,i2)

      enddo
      enddo

      return
      end

c**********************************************************************
c                                                                     *
c     This subroutine calculates the hh, pp and hp type of diagrams   *
c                                     T+ T                            *
c     coming from the contraction of e  e to calculate the normaliza  *
c     -tion constant in the property code                             *
c                                                                     *
c**********************************************************************


      subroutine overlap(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kapp/nka(MNBAS)

      dimension kgot(MXV),igot(MXV),jgot(MXV)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=1,MXV
      kgot(i)=0
      igot(i)=0
      enddo

c---------------------------------------------------------------
c                    PP type diagram
c---------------------------------------------------------------

      do 10 ip=1,nexcit
      orbp=orbe(ip)
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)

      do 10 iq=1,nexcit
      orbq=orbe(iq)
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)

      if(ipsym.ne.iqsym)go to 10
      if(orbp.ne.orbq)go to 10

      sum=0

c-----------------------------------------------------
c     diagram 1   <q|t1|a><a|t1|p>
c-----------------------------------------------------

      do 15 ia=1,nocc
      iia=kc(ia)
      orba=orbc(ia)
      iasym=isymc(ia)
      iaa=iqc(ia)

      if(orba.ne.orbp)go to 15
      if(ipsym.ne.iasym)go to 15

      ini=idpa(iip,iia)
      inf=idpa(iiq,iia)

      sum=sum-ti(ini)*tf(inf)

 15   continue

      do 20 ia=1,nocc
      iia=kc(ia)
      orba=orbc(ia)
      iasym=isymc(ia)
      iaa=iqc(ia)

      do 20 ib=1,nocc
      iib=kc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      iab=iqc(ib)

      do 20 ir=1,nexcit
      iir=ke(ir)
      orbr=orbe(ir)
      irsym=isyme(ir)
      iar=iqe(ir)

      iabsym=mtbl(iasym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      iprsym=mtbl(ipsym,irsym)

      if(iabsym.ne.iprsym)go to 20

c----------------------------------------------------------------
c     diagram 2 = <qr|t2|ab><pr|t2|ab>
c----------------------------------------------------------------

      call findk(orbr,orbb,orba,orbq,iar,iab,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 20
      do 25 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrb=idra(iir,iib,l1)
      nqa=idra(iiq,iia,l1)
      npa=idra(iip,iia,l1)
      inf=ntloc(nrb,nqa,irbsym,l1)
      ini=ntloc(nrb,npa,irbsym,l1)
      isign=(-1)**(orba+orbb+orbp+orbr)
      fact=isign/((2*l1+1)*(2*orbp+1))
      sum=sum-tf(inf)*ti(ini)*fact

c-----------------------------------------------------------------
c     diagram 3 = <qr|t2|ab><pr|t2|ba>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,jgot,jmax)
      if(jmax.eq.0)go to 25
      do 30 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npb=idra(iip,iib,l2)
      ini=ntloc(nra,npb,irasym,l2)
      isign=(-1)**(orba+orbb+orbp+orbr+2*orbr)
      call sixj(al1,orbr,orbb,al2,orbp,orba,s6j)
      fact=isign*s6j/(2*orbp+1)
      sum=sum+tf(inf)*ti(ini)*fact

 30   continue
 25   continue
 20   continue

      d(iip,iiq)=sum

 10   continue

c---------------------------------------------------------------
c            HH type diagram
c---------------------------------------------------------------
 
      do 40 ia=1,nocc
      iia=kc(ia)
      orba=orbc(ia)
      iasym=isymc(ia)
      iaa=iqc(ia)

      do 40 ib=1,nocc
      iib=kc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      iab=iqc(ib)

      if(orba.ne.orbb)go to 40
      if(iasym.ne.ibsym)go to 40

      sum=0.0d0

c-------------------------------------------------------------
c     diagram 1= <b|t1|p><p|t1|a>
c------------------------------------------------------------- 
      
      do 50 ip=1,nexcit
      iip=ke(ip)
      orbp=orbe(ip)
      ipsym=isyme(ip)
      iap=iqe(ip)

      if(ipsym.ne.iasym)go to 50
      if(orbp.ne.orba)go to 50
      if(ipsym.ne.ibsym)go to 50
      if(orbp.ne.orbb)go to 50

      inf=idpa(iip,iib)
      ini=idpa(iip,iia)

      sum=sum+tf(inf)*ti(ini)

 50   continue

      do 60 ip=1,nexcit
      orbp=orbe(ip)
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)

      do 60 iq=1,nexcit
      orbq=orbe(iq)
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)

      do 60 ig=1,nocc
      iig=kc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      iag=iqc(ig)

      iqgsym=mtbl(iqsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(iqgsym.ne.ipbsym)go to 60

c-----------------------------------------------------------------
c     diagram 2 = <bg|t2|pq><ag|t2|pq>
c-----------------------------------------------------------------

      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 60
      do 70 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npb=idra(iip,iib,l1)
      npa=idra(iip,iia,l1)
      inf=ntloc(nqg,npb,iqgsym,l1)
      ini=ntloc(nqg,npa,iqgsym,l1)
      isign=(-1)**(orbg+orbq+orbp+orbb)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum+tf(inf)*ti(ini)*fact

c-----------------------------------------------------------------
c     diagram 3 = <bg|t2|pq><ag|t2|qp>
c-----------------------------------------------------------------

      call findk(orbg,orbp,orbq,orba,iag,iap,iaq,iaa,jgot,jmax)
      if(jmax.eq.0)go to 70
      do 80 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqa=idra(iiq,iia,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nqa,npg,iqasym,l2)
      call sixj(al1,orbp,orbb,al2,orbq,orbg,s6j)
      fact=(-1)**(orbg+orbp+orba+orbq+2*orbq)*s6j/(2*orba+1)
      sum=sum-tf(inf)*ti(ini)*fact

 80   continue
 70   continue
 60   continue

      d(iib,iia)=sum

 40   continue
      
c-----------------------------------------------------------
c       diagram type PH or HP
c-----------------------------------------------------------

      
      do 100 ia=1,nocc
      iia=kc(ia)
      orba=orbc(ia)
      iasym=isymc(ia)
      iaa=iqc(ia)

      do 100 ip=1,nexcit
      iip=ke(ip)
      orbp=orbe(ip)
      ipsym=isyme(ip)
      iap=iqe(ip)

      if(orbp.ne.orba)go to 100
      if(ipsym.ne.iasym)go to 100

      in=idpa(iip,iia)

c-------------------------------------------------------------
c     diagram 1 = <p|t1|a>
c-------------------------------------------------------------

      d(iip,iia)=d(iip,iia)+ti(in)

      do 110 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iib=kc(ib)
      iab=iqc(ib)

      do 120 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iiq=ke(iq)
      iaq=iqe(iq)

      if(orbq.ne.orbb)go to 120
      if(iqsym.ne.ibsym)go to 120

c--------------------------------------------------------------
c     diagram 2 = <pq|t2|ab><q|t1|b>
c--------------------------------------------------------------

      npa=idra(iip,iia,0)
      nqb=idra(iiq,iib,0)
      iqbsym=mtbl(iqsym,ibsym)
      ini=ntloc(npa,nqb,iqbsym,0)
      inf=idpa(iiq,iib)
      fact=dsqrt((2*orbb+1)/(2*orba+1))
      d(iip,iia)=d(iip,iia)+tf(inf)*ti(ini)*fact

c--------------------------------------------------------------
c     diagram 3 = <pq|t2|ba><q|t1|b>
c--------------------------------------------------------------

      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.eq.0)go to 120
      do 130 iloop=1,imax
      l1=igot(iloop)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npb,nqa,iqasym,l1)
      fact=(-1)**(orba+orbb+2*orbb+l1)/(2*orba+1)
      d(iip,iia)=d(iip,iia)-tf(inf)*ti(ini)*fact

 130  continue
 120  continue
 110  continue

      d(iia,iip)=d(iip,iia)

 100  continue

      return
      end
c------------------------------------------------------------
c*****************************************************************
c     This subroutine calculated the NMS O bar diagrams of  *
c     of the type  HH                                            *
c*************************************************************** *

      subroutine nmshh(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)

      if(ibsym.ne.iasym)go to 1

c------------------------------------------------------------------
c     diagram 1  = <b|d|a>
c------------------------------------------------------------------

      call nms(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      d(iib,iia)=d(iib,iia)+dipole

      sum=zero

      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      
c----------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c----------------------------------------------------------------
      
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      
      ini=idpa(iip,iia)

      call nms(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)

      endif
      endif
      endif
      
c----------------------------------------------------------------
c     diagram 3 =  <b|t1|p><p|d|a>
c----------------------------------------------------------------
      
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.eq.iasym)then

      inf=idpa(iip,iib)

      call nms(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)

      endif
      endif
      endif

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.eq.ipbsym) then

c----------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c----------------------------------------------------------------
      
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.iasym)then
      
      call nms(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)
      
      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      if(orbb.eq.orbp)then
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      
      endif
      endif
      enddo
      endif
c--------------------------------------------------------------------
c     diagrams 5 =  <bg|t2|rp><p|d|a><r|t1|g>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.eq.igsym)then

      call nms(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c--------------------------------------------------------------------
      
      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npb=idra(iip,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(npb,nrg,irgsym,1)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      
      endif
      enddo
      endif
      
c--------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c--------------------------------------------------------------------
      
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      ini=idpa(iip,iia)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then

c------------------------------------------------------------------
c     diagrams 8 = <g|d|r><b|t1|p><pr|t2|ag>
c------------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then
      
      call nms(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      
      inf=idpa(iip,iib)
      
      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      
      endif
      enddo
      endif
      
c-------------------------------------------------------------------
c     diagrams 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      inf=idpa(iip,iib)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.ibsym)then

      call nms(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)

      inf=idpa(iir,iig)
      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      if(orbp.eq.orba)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      
      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------
      
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
      if(orba.eq.orbp)then
c     sum=sum-dipole*tf(inf)*tf(inf)*fact
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(ipsym.ne.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c-------------------------------------------------------------------
c     diagram 12 = <b|t1|r><t|d|p><p|t1|a>
c-------------------------------------------------------------------
      
      call nms(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)
 
      inf=idpa(iir,iib)
      ini=idpa(iip,iia)

      sum=sum+dipole*tf(inf)*ti(ini)
      
 5    continue
      
      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

c-----------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c-----------------------------------------------------------------

      if(igsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call nms(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ini=idpa(iip,iia)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c-------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c-------------------------------------------------------------------
      
      if(igsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call nms(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)
      
      ini=idpa(iip,iig)
      inf=idpa(iip,iib)
      
      sum=sum-dipole*ti(ini)*tf(inf)
      
      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue
      
      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 8

      call nms(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      
      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8
      
      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+1+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 10   continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,one,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11   continue
  9   continue
  8   continue

      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
 
      if(iqsym.ne.irsym)go to 13

      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      
      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13
      
      call nms(orbq,orbr,iiq,iir,orbq,orbr,kapq,kapr,dipole)
      
c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------
      
      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)
      
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,one,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue

c-------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(one,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+one+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      
      enddo
      
 14   continue
 13   continue
      
      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.ne.irsym)go to 16  
       
      call nms(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

      
c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------
      
      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      call sixj(one,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+one)
      fact=s6j*isign/(2*al1+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,one,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+one)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact

  18  continue
  17  continue
  16  continue
  12  continue

      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(idsym.ne.ibsym)go to 21

      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22
      
      call nms(orbb,orbd,iib,iid,orbb,orbd,kapb,kapd,dipole)
      
c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)
      
      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 22= <gd|t2|pr><b|d|d><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 24   continue
 23   continue
 22   continue
 21   continue
      
      if(idsym.ne.iasym)go to 25
      
      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26
      
      call nms(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)
      
c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------
      
      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(npd,nrg,irgsym,l1)
 
      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------------

      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact
      
 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue
      
      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+1)

 1    continue   

      return
      end
      
c*****************************************************************
c     This subroutine calculated the NMS O bar diagram of  *
c     of the type  HP                                            *
c*************************************************************** *
      
      subroutine nmshp(d,ti,tf)


      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      
      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)
      
      if(ipsym.ne.iasym)go to 1
      
c-------------------------------------------------------------------
c     diagram 1  = <p|d|a>
c-------------------------------------------------------------------
      
      call nms(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      
      d(iip,iia)=d(iip,iia)+dipole
      
      sum=zero
      
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      
      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 3

      call nms(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.irgsym)go to 401

c---------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c---------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,ipasym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(2*one+1)
      sum=sum+dipole*ti(ini)*fact
      
      endif
      enddo
      endif
      
c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c---------------------------------------------------------------------
      
      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then   
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dipole*ti(ini)*fact
      
 4    continue
      
      endif
      
 401  continue
      
c--------------------------------------------------------------------------
c     diagram 4 = <r|t1|a><r|d|g><p|t1|g>
c--------------------------------------------------------------------------

      if(igsym.ne.irsym)go to 3
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3

      call nms(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ini=idpa(iir,iia)
      inj=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*ti(inj)

 3    continue

c-------------------------------------------------------------------------
c     diagram 5 =  <r|t1|a><g|t1|r><p|d|g>
c-------------------------------------------------------------------------

      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5

      call nms(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dipole)

      ini=idpa(iir,iia)
      inf=idpa(iir,iig)
      
      sum=sum-dipole*ti(ini)*tf(inf)
      
 5    continue
      
c-------------------------------------------------------------------
c     diagram 6 = <r|d|a><g|t|r><p|t|g>
c-------------------------------------------------------------------
      
      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6

      call nms(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iip,iig)
      
      sum=sum-dipole*ti(ini)*tf(inf)
      
 6    continue
      
      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7

      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym)go to 8

      call nms(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.iqasym)go to 8

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 7 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      if(orba.eq.orbq)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))

      endif
      endif
      enddo
      endif
      
c---------------------------------------------------------------------
c     diagram 8 =  <g|t1|r><p|d|q><qr|t2|ga>
c---------------------------------------------------------------------
      
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 9    continue
 8    continue


      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      if(iasym.ne.ibsym)go to 10

      call nms(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)

      if(irgsym.ne.ipbsym)go to 10

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 9 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbp.eq.orbb)then
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))
      
      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 10 = <g|t1|r><b|d|a><pr|t2|gb>
c-----------------------------------------------------------------
      
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10 
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      
 11   continue
 10   continue
  7   continue


      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqasym) goto 13
      if(ipsym.ne.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13

      call nms(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dipole)
 
c-----------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)
      
      call findk(orbr,orbg,orbq,orba,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue
 14   continue
 13   continue


      if(orbp.ne.orbq)go to 12

      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12
      
      call nms(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dipole)
      
c--------------------------------------------------------------------
c     diagram 13 =  <gb|t2|rq><pr|t2|bg><q|d|a>
c--------------------------------------------------------------------
      
      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)
 
      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 17   continue
 16   continue
 12   continue


      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2
      
      inf=idpa(iir,iig)
      
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)

      if(irsym.ne.iqsym)go to 20
      
      call nms(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)
      
      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      
      if(ipasym.ne.iqgsym)go to 20
      
c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------
  
      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(one,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 21   continue
 20   continue

      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)
      
      if(ibsym.ne.igsym)go to 22
      
      call nms(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      
      if(ipasym.ne.irbsym)go to 22
      
c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------
      
      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nrb=idra(iir,iib,1)
      ini=ntloc(npa,nrb,ipasym,1)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*1+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(one,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 23   continue
 22   continue
  2   continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------
      
      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24
      
      ini=idpa(iiq,iia)
      
      call nms(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
      
      sum=sum+ti(ini)*dipole
      
 24   continue
      
c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------
 
      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25

      ini=idpa(iip,iib)

      call nms(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      sum=sum-ti(ini)*dipole

 25   continue

      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
 
      if(igsym.ne.irsym)go to 26
      
      call nms(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)
      
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27
      
c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*1+1)*(2*1+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

 27   continue

c------------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg>
c------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      
      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(one,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*1+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 29   continue

      endif
      endif
      enddo
      endif
      
 28   continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(one,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact
      
 31   continue
      
      endif
      endif
      enddo
      endif

 30   continue
      
c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------
      
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      
      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32
      
      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(one,orbr,orbg,al1,orbq,orbb,s6j1)
      call sixj(one,orba,orbp,al2,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 34   continue
 33   continue
 32   continue
 26   continue


      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+1)

 1    continue
      
      return
      end
      
c*****************************************************************
c     This subroutine calculates the NMS O bar diagram of  *
c     of the type  PP                                            *
c*************************************************************** *
      
      subroutine nmspp(d,ti,tf)
      
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)
 
      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)
      
      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym) go to 1

      call nms(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

c-----------------------------------------------------------------
c     diagram 1 = <p|d|q>
c-----------------------------------------------------------------

      d(iip,iiq)=d(iip,iiq)+dipole
      
      sum=0.0d0
        
      sum6=0
      sum8=0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      if(iasym.eq.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then
      
      call nms(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)
      
      ini=idpa(iip,iia)
      
c----------------------------------------------------------------
c     diagram 2 = <p|t|a><a|d|q>
c----------------------------------------------------------------

      sum=sum-dipole*ti(ini)

      endif
      endif
      endif

      if(iasym.eq.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then

      call nms(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

      sum=sum-dipole*tf(inf)

      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.eq.ipsym)then

      call nms(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orbq.eq.orba)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif
      
c-----------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c-----------------------------------------------------------------
      
      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)then
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      
      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.eq.irsym)then

      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then

      call nms(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      nqa=idra(iiq,iia,1)
      inf=ntloc(nrg,nqa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      
      endif
      enddo
      endif
      
c----------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c----------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      
      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.eq.irsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call nms(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      npa=idra(iip,iia,1)
      ini=ntloc(nrg,npa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      
      endif
      enddo
      endif
      
c----------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      
      enddo
      endif
      endif
      endif
      endif
      endif

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.eq.iqsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call nms(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orba.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif
      
c-----------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c-----------------------------------------------------------------
      
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)then
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      
      endif
      enddo
      endif
      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

c-----------------------------------------------------------------
c     diagram 12 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

      if(irsym.eq.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then

      call nms(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)
      
      endif
      endif
      endif
      endif
      endif
      
c---------------------------------------------------------------
c     diagram 13 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------
      
      if(irsym.eq.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then

      call nms(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,dipole)
      
      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)
      
      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  5   continue

c------------------------------------------------------------------
c     diagram 14 = <q|t1|a><a|d|g><g|t1|p
c------------------------------------------------------------------

      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      if(igsym.eq.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then

      call nms(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)
      
      sum=sum+dipole*ti(ini)*tf(inf)
      
      endif
      endif
      endif
      endif
      endif

 6    continue     
 4    continue 

      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      
      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 9

      call nms(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      
 11   continue
 10   continue

c-----------------------------------------------------------------
c     diagram 16 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)

      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,one,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      
 14   continue
 13   continue
 12   continue
  9   continue
      
      if(iasym.ne.idsym)go to 8 

      call nms(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)
      
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)
      
      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8
      
c-----------------------------------------------------------------
c     digram 17 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------
      
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(one,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+one+l1)/(2*al1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+one+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,one,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      
 16   continue
 15   continue
  8   continue
      
c------------------------------------------------------------------
c     diagram 19 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------
      
      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is) 
      kaps=nka(iis)

      if(issym.ne.irsym)go to 18
      
      call nms(orbs,orbr,iis,iir,orbs,orbr,kaps,kapr,dipole)
      
      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18

      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,one,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 20   continue

c------------------------------------------------------------------
c     diagram 20 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------
      
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax 
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,one,orbg,al2,orbp,s9j)
c      fact=(-1)**(orbq+orbg+l2+one+2*orbp)*s9j           
      fact=(-1)**(orba+orbg+orbr+orbs+l1+l2)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      
 21   continue
 19   continue
 18   continue
      
c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------
      
      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)

      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.ne.issym)go to 22

      call nms(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,dipole)

      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)   
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      
 24   continue
 23   continue
 22   continue
      
c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)

      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.ne.issym)go to 25

      call nms(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,dipole)

      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      
 27   continue
 26   continue
 25   continue
 17   continue
  7   continue

      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+1)

  1   continue

      return
      end

c*********************************************************
c*********************************************************
c  This subroutine calculates reduced matrix element of  *
c  the NMS (isotopic shift) operator                                *
c*********************************************************

      subroutine nms(j1,j2,i1,i2,m1,m2,kap1,kap2,val)

      implicit real*8(a-h,o-z)
      real*8 j1,j2,m1,m2

      dimension nmsr(MNBAS,MNBAS)
      common/cons/zero,half,tenth,one,two,three,ten

       call Calcnms(i1,kap1,nmsr)
       factor=(two*j2+one)
       w3j1=dr(j1,one,j2,half,zero,-half)
       e1=factor*(w3j1**2)
       result=nmsr(i1,i2)
       val=result*e1
       return
       end


  
c********************************************************************
      
c*********************************************************
c  This subroutine calculates the radial integral of the *
c  reduced matrix element of the NMS operator      *
c*********************************************************

      subroutine Calcnms(i1,kap1,nmsr)

      implicit double precision(a-h,o-z)

       dimension gem1(MN),gem2(MN),gem3(MN),gem4(MN)
     :  ,dpf(MN,MNBAS),dqf(MN,MNBAS),nmsr(MNBAS,MNBAS)
cc     :,dppf(MN,MNBAS),dqqf(MN,MNBAS)
cc     :  ,gemb1(MN),gemb2(MN)

      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tptq/tpf(MN),tqf(MN)
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
cc      common/radial/nmsr(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)


       mtp=n
        i2=i1
         ka=kap1

       if (ka.gt.0) then
         la = ka
       else
         la = -ka-1
       endif

        ja = iabs(ka)*2-1

        la_til = ja-la

         call diff1(i1)
       do i=1,mtp
         dpf(i,i1)=tpf(i)
         dqf(i,i1)=tqf(i)
       enddo


       do i=1,mtp

         gem1(i)=zero
         gem2(i)=zero
         gem3(i)=zero
         gem4(i)=zero

       enddo

         kk1=dble(la*(la+1))
         kk2=dble(la_til*(la_til+1))


cc       do i1=1,nbasis
cc       do i2=1,i1

         
       nmsr(i1,i2) = zero

cc       if(iparity(i1).eq.iparity(i2))then

cc       ta(1)=zero

       do i=2,mtp
       
         gem1(i)=((kk1)/(r(i)*r(i)))*pf(i,i1)*pf(i,i2)*rp(i)
         
         gem2(i)=(dqf(i,i1)*dqf(i,i2))/rp(i)

         gem3(i)=((kk2)/(r(i)*r(i)))*qf(i,i1)*qf(i,i2)*rp(i)
       
         gem4(i)=(dpf(i,i1)*dpf(i,i2))/rp(i)
         
cc------------------------------------------------------------------
cc       ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)/(r(l)*r(l))
cc-------------------------------------------------------------------
        enddo
       
         call quad_new (gem1,rll1)
         
         call quad_new (gem2,rdss2)
         rdss2=rdss2/(h*h)

         call quad_new (gem3,rss3)

         call quad_new (gem4,rdll4)
         rdll4=rdll4/(h*h)

cc       call quad(result)
         result=(rll1+rdss2+rss3+rdll4)*half

       nmsr(i1,i2)=result
       nmsr(i2,i1)=result



cc       enddo
cc       enddo

       return
       end

c************************************************************

c**********************************************************************
!---------------------------------------------------------
      subroutine quad_new(aa,result)
!--------------------------------------------------------
      implicit real*8(a-h, o-z)
      dimension aa(MN)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     : /ncc/c1,c2,c3,c4
!------------------------------------------------------------------------------
!   find first values that will permit computation of exponent
!------------------------------------------------------------------------------
      result = zero
      mtp=n
      mtpm1 = mtp - 1
      do 3 i = 2,mtpm1
        tai = aa(i)
        if (dabs (tai) .gt. zero) then
          ip1   = i+1
          taip1 = aa(ip1)
          quott = taip1/tai
          if (quott .gt. zero) then
!------------------------------------------------------------------------------
!                          exponent from fit
!------------------------------------------------------------------------------
               frip1 = taip1/rp(ip1)
               fri   = tai  /rp(i  )
               ratio = frip1/fri
               rip1  = r (ip1)
               ri    = r (i  )
               sigma = dlog (ratio)/dlog (rip1/ri)
!-----------------------------------------------------------------------------
!       analytical integration and error estimate for interval r(1:i)
!-----------------------------------------------------------------------------
               fri    = ri*fri
               result = fri/(sigma+one)
!------------------------------------------------------------------------------
!                       set the tail to zero
!------------------------------------------------------------------------------
               do 1 loc = 1,3
                  aa(mtp+loc) = zero
    1          continue
!------------------------------------------------------------------------------
!              newton-cotes quadature for the remainder
!------------------------------------------------------------------------------
               result = result+c1*tai
               do 2 loc = ip1,mtp,4
                  result = result+c2*(aa(loc  )+aa(loc+2))
     :                           +c3* aa(loc+1)
     :                           +c4* aa(loc+3)
    2          continue
               if (mod (mtp-i,4) .eq. 0) result = result-c1*aa(mtp)
             goto 4
            endif
         endif
    3 continue
!------------------------------------------------------------------------------
!          no value which will permit computation of exponent
!------------------------------------------------------------------------------
      result = zero
    4 return
      end
c**********************************************************
c              differentiation    
c*********************************************************
!-----------------------------------------------------------
      subroutine diff1(J)
!-----------------------------------------------------------
      implicit real*8 (a-h, o-z)

      COMMON/cons/zero,half,tenth,one,two,three,ten
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /lic13/a(13,13)
     :      /tptq/tpf(MN),tqf(MN)
     :      /wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),mf(MNBAS)

      equivalence (a1,a(7,1)),(a2,a(7,2)),(a3,a(7,3)),
     :            (a4,a(7,4)),(a5,a(7,5)),(a6,a(7,6))

!   COMPUTE DERIVATIVE IN THREE SEPARATE REGIONS

!   FIRST, POINTS 1 TO 6

      do i = 1,6
         hdpbdt = 0.0d 00
         hdqbdt = 0.0d 00
         do k = 1,13
            aik = a(i,k)
            hdpbdt = hdpbdt+aik*pf(k,J)
            hdqbdt = hdqbdt+aik*qf(k,J)
         enddo
         tpf(i) = hdpbdt
         tqf(i) = hdqbdt
      enddo
      do i = 7,n-6
         tpf(i) =  a1*(pf(i-6,J)-pf(i+6,J))+a2*(pf(i-5,J)-pf(i+5,J))
     :           +a3*(pf(i-4,J)-pf(i+4,J))+a4*(pf(i-3,J)-pf(i+3,J))
     :           +a5*(pf(i-2,J)-pf(i+2,J))+a6*(pf(i-1,J)-pf(i+1,J))
         tqf(i) =  a1*(qf(i-6,J)-qf(i+6,J))+a2*(qf(i-5,J)-qf(i+5,J))
     :           +a3*(qf(i-4,J)-qf(i+4,J))+a4*(qf(i-3,J)-qf(i+3,J))
     :           +a5*(qf(i-2,J)-qf(i+2,J))+a6*(qf(i-1,J)-qf(i+1,J))
       enddo

!   LAST, POINTS N-5 TO N

      do i = n-5,n
         irow = i-n+13
         hdpbdt = 0.0d 00
         hdqbdt = 0.0d 00
         do k = 1,13
            aik = a(irow,k)
            loc = n-13+k
            hdpbdt = hdpbdt+aik*pf(loc,J)
            hdqbdt = hdqbdt+aik*qf(loc,J)
         enddo
         tpf(i) = hdpbdt
         tqf(i) = hdqbdt
       enddo
       return
       end
!---------------------------------------------------------


c============================================================   
      subroutine findk(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)
c=============================================================
c     this subroutine finds the common k value
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'

      common/iallparity/iall

      dimension jeven(MXV+1)

       if(iall.eq.0)then

      call findk2(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,kkk)
       return

       else

      if(iall.eq.3)then
      call findk_br(orba,orbb,orbp,orbq,ia,ib,ip,iq,jeven,kkk)
      return

       else

      eps=1.0d-10

      kkk=0
      do k=0, MXV
         if(idelta(orba,orbb,k).eq.1.and.
     &      idelta(orbp,orbq,k).eq.1) then
            kkk=kkk+1
            jeven(kkk)=k
         endif
      enddo

         endif
         endif

      return
      end


      subroutine findk_br(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)
C     THIS SUBROUTINE FINDS THE ALLOWED MULTIPOLE MOMENT
C     VARIBALES "ORBA, ORBB ETC. ARE J-VALUE OF THE A, B
C     ORBITAL, IA, IB ETC. ARE 'IA" VALUE.
C
C     ON EXIT KKK--> NO. OF ALLOWED MOMENT AND JEVEN
C     STORES THE MOMENT VALUE
      implicit real*8 (a-h,o-z)
c      include 'mpif.h'
      real*8 jj
      dimension ja(MXV),jb(MXV),jp(MXV),jq(MXV)
      dimension jeven(MXV+1)
      do i=1,MXV
      jeven(i)=0
      enddo
      iab=ia*ib
      ipq=ip*iq
      jminab=idint(dmax1(orba,orbb)-dmin1(orba,orbb))
      jminpq=idint(dmax1(orbp,orbq)-dmin1(orbp,orbq))
      jmaxab=idint(orba+orbb)
      jmaxpq=idint(orbp+orbq)
      ll=0
      kk=0
      do i=jminab,jmaxab
      jtot=jmaxab+i
      jj=dfloat(jtot)-2.0*(dfloat(jtot/2))
      if(jj.ne.0.and.iab.lt.0)then
      ll=ll+1
      ja(ll)=i
      endif
      if(jj.eq.0.and.iab.gt.0)then
      kk=kk+1
      jb(kk)=i
      endif
      enddo
      mm=0
      nn=0
      do i=jminpq,jmaxpq
      jtot=jmaxpq+i
      jj=dfloat(jtot)-2.0*(dfloat(jtot/2))
      if(jj.ne.0.and.ipq.lt.0)then
      mm=mm+1
      jp(mm)=i
      endif
      if(jj.eq.0.and.ipq.gt.0)then
      nn=nn+1
      jq(nn)=i
      endif
      enddo
      kkk=0
      do i=1,ll
      do j=1,mm
      if(ja(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
      endif
      enddo
      enddo
      do i=1,ll
      do j=1,nn
      if(ja(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=ja(i)
      endif
      enddo
      enddo
      do i=1,kk
      do j=1,mm
      if(jb(i).eq.jp(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
      endif
      enddo
      enddo
      do i=1,kk
      do j=1,nn
      if(jb(i).eq.jq(j))then
      kkk=kkk+1
      jeven(kkk)=jb(i)
      endif
      enddo
      enddo

      return
      end

      integer function idelta(ja,jp,k)
      implicit none
      real*8 ja,jp,eps
      integer k,kmax,kmin
      eps=1.0d-10
      kmax = int(ja + jp + eps)
      kmin = int(abs(ja - jp)+eps)
      if (k.le.kmax.and.k.ge.kmin) then
         idelta = 1   ! triangle condition satisfied
      else
         idelta = 0   ! not satisfied
      endif
      return
      end

c**********************************************************************


c****************************************************************
c         This subroutine calculates the reduced matrix element *
c         of E1 operator in velocity gauge                      *
c****************************************************************


      subroutine e1vel(i1,i2,jorb1,jorb2,kap1,kap2,e1m)

      implicit double precision (a-h,o-z)
      double precision j1,j2,m1,m2,jorb1,jorb2
      common/cons/zero,half,tenth,one,two,three,ten
      common/inte1vel/rint1(MNBAS,MNBAS),rint2(MNBAS,MNBAS)
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)

      common/info/nocc,nexcit
      common/orbital_energy/eorb(MNBAS)

      nbasis=nocc+nexcit

      alpha=137.0359895
      e1m=0.0

      fact=eorb(i2)-eorb(i1)
      tj1=idint(2.0*jorb1)
      tj2=idint(2.0*jorb2)
      factor=(-1)**(jorb1+0.5)
      j1=tj1/two
      j2=tj2/two
      w3j1=dr(j1,one,j2,half,zero,-half)
      e1rad=rint1(i1,i2)
      e2rad=rint2(i1,i2)
      e1=dsqrt(tj1+one)*dsqrt(tj2+one)*factor*w3j1
     :*(e1rad - (kap1-kap2)*e2rad)
      e1m=e1*alpha/fact

      return
      end

c****************************************************************************

c***********************************************************************
c                                                                      *
c          This subroutine calculate the radial part of the E1 matrix  *
c          element in velocity gauge                                   *
c                                                                      *
c***********************************************************************



      subroutine Calce1vel(nbasis)

      implicit real*8 (a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/inte1vel/rint1(MNBAS,MNBAS),rint2(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      mtp = n

c------------------------anty symmetric part------------------------

      do i1 = 1, nbasis
      do i2 = 1, nbasis

      rint1(i1,i2) = zero

      ta(1)=zero
      do 1 l=2,mtp
      ta(l)=(pf(l,i1)*qf(l,i2)-pf(l,i2)*qf(l,i1))*rp(l)
    1 continue
      call quad(result)

      rint1(i1,i2) = result

      enddo
      enddo

c----------------------symmetric part-------------------------------

      do i1 = 1, nbasis
      do i2 = 1, i1

      rint2(i1,i2) = zero

      ta(1)=zero
      do 2 l=2,mtp
      ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)
    2 continue
      call quad(result1)

      rint2(i1,i2) = result1
      rint2(i2,i1) = rint2(i1,i2)

      enddo
      enddo

      return
      end

c**********************************************************************

c*****************************************************************
c     This subroutine calculates the E2 transition O bar diagram *
c     of the type  HH                                            *
c*************************************************************** *


      subroutine quadpulehh(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)

      if(ibsym.ne.iasym)go to 1

c------------------------------------------------------------------
c     diagram 1  = <b|d|a>
c------------------------------------------------------------------

      call quadpul(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      d(iib,iia)=d(iib,iia)+dipole

      sum=zero

      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

c----------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c----------------------------------------------------------------

      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.eq.ibsym)then

      ini=idpa(iip,iia)

      call quadpul(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)

      endif
      endif
      endif

c----------------------------------------------------------------
c     diagram 3 =  <b|t1|p><p|d|a>
c----------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.eq.iasym)then

      inf=idpa(iip,iib)

      call quadpul(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)

      endif
      endif
      endif

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.eq.ipbsym) then

c----------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c----------------------------------------------------------------

      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.iasym)then

      call quadpul(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      if(orbb.eq.orbp)then
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 5 =  <bg|t2|rp><p|d|a><r|t1|g>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.eq.igsym)then

      call quadpul(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c--------------------------------------------------------------------

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npb=idra(iip,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(npb,nrg,irgsym,2)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(two*two+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+2)*s6j
      ini=idpa(iip,iia)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then

c------------------------------------------------------------------
c     diagrams 8 = <g|d|r><b|t1|p><pr|t2|ag>
c------------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then

      call quadpul(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,irgsym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(two*two+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+2)*s6j
      inf=idpa(iip,iib)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.ibsym)then

      call quadpul(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)

      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      if(orbp.eq.orba)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
      if(orba.eq.orbp)then
c     sum=sum-dipole*tf(inf)*tf(inf)*fact
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(ipsym.ne.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c-------------------------------------------------------------------
c     diagram 12 = <b|t1|r><t|d|p><p|t1|a>
c-------------------------------------------------------------------

      call quadpul(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

      inf=idpa(iir,iib)
      ini=idpa(iip,iia)

      sum=sum+dipole*tf(inf)*ti(ini)

 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

c-----------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c-----------------------------------------------------------------

      if(igsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call quadpul(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ini=idpa(iip,iia)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c-------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c-------------------------------------------------------------------

      if(igsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call quadpul(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iip,iib)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 8

      call quadpul(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,two,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,two,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+2+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 10   continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,two,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11   continue
  9   continue
  8   continue


      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(iqsym.ne.irsym)go to 13

      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13

      call quadpul(orbq,orbr,iiq,iir,orbq,orbr,kapq,kapr,dipole)

c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,two,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,two,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue

c-------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(two,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+two+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo

 14   continue
 13   continue

      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.ne.irsym)go to 16

      call quadpul(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      call sixj(two,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+two)
      fact=s6j*isign/(2*al1+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,two,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+two)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact

  18  continue
  17  continue
  16  continue
  12  continue

      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(idsym.ne.ibsym)go to 21

      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22

      call quadpul(orbb,orbd,iib,iid,orbb,orbd,kapb,kapd,dipole)

c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 22= <gd|t2|pr><b|d|d><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 24   continue
 23   continue
 22   continue
 21   continue

      if(idsym.ne.iasym)go to 25

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26

      call quadpul(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(npd,nrg,irgsym,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------------

      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

       d(iib,iia)=d(iib,iia)+sum
       d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+one)

 1    continue

      return
      end

c*****************************************************************
c     This subroutine calculates the E2 transition O bar diagram *
c     of the type  HP                                            *
c*************************************************************** *

      subroutine quadpulehp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)

      if(ipsym.ne.iasym)go to 1

c-------------------------------------------------------------------
c     diagram 1  = <p|d|a>
c-------------------------------------------------------------------

      call quadpul(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      d(iip,iia)=d(iip,iia)+dipole

      sum=zero

      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 3

      call quadpul(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.irgsym)go to 401

c---------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c---------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,ipasym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(2*two+1)
      sum=sum+dipole*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum-dipole*ti(ini)*fact

 4    continue

      endif

 401  continue

c--------------------------------------------------------------------------
c     diagram 4 = <r|t1|a><r|d|g><p|t1|g>
c--------------------------------------------------------------------------

      if(igsym.ne.irsym)go to 3
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3

      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
    
      call quadpul(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      sum=sum-dipole*ti(ini)*ti(inj)

 3    continue

c-------------------------------------------------------------------------
c     diagram 5 =  <r|t1|a><g|t1|r><p|d|g>
c-------------------------------------------------------------------------

      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5

      call quadpul(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dipole)

      ini=idpa(iir,iia)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 5    continue

c-------------------------------------------------------------------
c     diagram 6 = <r|d|a><g|t|r><p|t|g>
c-------------------------------------------------------------------

      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6

      call quadpul(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dipole)

      inf=idpa(iir,iig)
      ini=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 6    continue

      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7

      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym)go to 8

      call quadpul(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.iqasym)go to 8

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 7 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      if(orba.eq.orbq)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))

      endif
      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 8 =  <g|t1|r><p|d|q><qr|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 9    continue
 8    continue


      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      if(iasym.ne.ibsym)go to 10

      call quadpul(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)

      if(irgsym.ne.ipbsym)go to 10

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 9 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbp.eq.orbb)then
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbp+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 10 = <g|t1|r><b|d|a><pr|t2|gb>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue
  7   continue


      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqasym) goto 13
      if(ipsym.ne.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13

      call quadpul(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dipole)

c-----------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orba,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue
 14   continue
 13   continue

      if(orbp.ne.orbq)go to 12

      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12

      call quadpul(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dipole)

c--------------------------------------------------------------------
c     diagram 13 =  <gb|t2|rq><pr|t2|bg><q|d|a>
c--------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 17   continue
 16   continue
 12   continue

      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2

      inf=idpa(iir,iig)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)

      if(irsym.ne.iqsym)go to 20

      call quadpul(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipasym.ne.iqgsym)go to 20

c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nqg=idra(iiq,iig,2)
      ini=ntloc(npa,nqg,ipasym,2)
      fact=(-1)**(orbg+orbq+2*orbg+2)/(2*2+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(two,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+2)*s6j
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 21   continue
 20   continue

      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)

      if(ibsym.ne.igsym)go to 22

      call quadpul(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)

      if(ipasym.ne.irbsym)go to 22

c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nrb=idra(iir,iib,2)
      ini=ntloc(npa,nrb,ipasym,2)
      fact=(-1)**(orbb+orbr+2*orbb+2)/(2*2+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(two,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+2)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 23   continue
 22   continue
  2   continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------

      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24

      ini=idpa(iiq,iia)

      call quadpul(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      sum=sum+ti(ini)*dipole

 24   continue

c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------

      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25

      ini=idpa(iip,iib)

      call quadpul(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      sum=sum-ti(ini)*dipole

 25   continue

      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 26

      call quadpul(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------


      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nqb,ipasym,2)
      inf=ntloc(nrg,nqb,iqbsym,2)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*2+1)*(2*2+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

 27   continue

c------------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg>
c------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(nrg,nqb,iqbsym,2)

      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(two,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*2+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 29   continue

      endif
      endif
      enddo
      endif

 28   continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      ini=ntloc(npa,nqb,iqbsym,2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(two,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2*2+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*2+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact

 31   continue

      endif
      endif
      enddo
      endif

 30   continue

c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32

      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(two,orbr,orbg,al1,orbq,orbb,s6j1)
      call sixj(two,orba,orbp,al2,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 34   continue
 33   continue
 32   continue
 26   continue

      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+one)
       
 1    continue
      return
      end


      subroutine quadpuleph(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      nbasis=nocc+nexcit

c     P-H type <p|d|a>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nak(iip)
      if(ipsym.ne.iasym)go to 1
c     diagram 1  = <p|d|a>
      call quadpul(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      d(iip,iia)=d(iip,iia)+dipole
      sum=zero

c     Diagrams that are linear in T
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      if(igsym.ne.irsym)go to 3
      call quadpul(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.irgsym)go to 401
c     diagram 2 = sum(r,g) <pr|t2|ag><g|d|r>

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,ipasym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(2*two+1)
      sum=sum+dipole*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 3 = -sum(r,g)<g|d|r><pr|t2|ga>
      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.eq.0)go to 401
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum-dipole*ti(ini)*fact
 4    continue
 401  continue

c     Non-linear contribution
c     diagram 8 = -sum(r,g) <r|t1|a><g|d|r><p|t1|g>
      if(ipsym.ne.igsym)go to 3
      if(irsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3
      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
      call quadpul(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      sum=sum-dipole*ti(ini)*ti(inj)
 3    continue
c     diagram 9 = -sum(r,g) <r|t1|a><g|t1|r><p|d|g>
      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5
      call quadpul(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dipole)
      ini=idpa(iir,iia)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 5    continue
c     diagram 10 = -sum(r,g) <r|d|a><g|t|r><p|t|g>
      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6
      ini=idpa(iir,iig)
      inf=idpa(iip,iig)
      call quadpul(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dipole)
      sum=sum-dipole*ti(ini)*tf(inf)
 6    continue

      if(igsym.ne.irsym)go to 7
      if(orbg.ne.orbr)go to 7
c     statement valid for diagrams 4,5,6,7
      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      if(ipsym.ne.iqsym)go to 8
c     statement valid for diagrams 4,5
      call quadpul(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.iqasym)go to 8
      inf=idpa(iir,iig)

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     diagram 4 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ag>
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbp+one))
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(2*orbq+1))
      endif
      enddo
      endif

c     diagram 5 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ga>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 9    continue
 8    continue

      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)
      if(iasym.ne.ibsym)go to 10
      call quadpul(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)
      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)
      if(irgsym.ne.ipbsym)go to 10
      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
c     diagram 6 = -sum(b,r,g) <g|t1|r><b|d|a><pr|t2|bg>
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbq+one))
      endif
      enddo
      endif

c     diagram 7 = sum(b,r,g) <g|t1|r><b|d|a><pr|t2|gb>
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      if(orbp.ne.orbb)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
  7   continue

c     diagram 11,12, 13 and 14
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)
      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)
      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      if(ipsym.ne.ibsym)go to 13
c     statement valid for diagram 11, 12
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqasym) goto 13
      call quadpul(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dipole)
c     diagram 11 = -sum (q,g,r,b)<gb|t2|rq><qr|t2|ag><p|d|b>
      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
c     diagram 12 = sum (q,g,r,b)<gb|t2|rq><qr|t2|ga><p|d|b>
      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
 14   continue
 13   continue

      if(orbp.ne.orbq)go to 12
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12
c     diagram 13 = -sum(q,r,g,b) <gb|t2|rq><pr|t2|bg><q|d|a>
      call quadpul(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dipole)
      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      npb=idra(iip,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
c     diagram 14 = sum(q,r,g,b) <gb|t2|rq><pr|t2|gb><q|d|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 17   continue
 16   continue
 12   continue

c     diagram 15 to 18
      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2
      inf=idpa(iir,iig)
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nak(iiq)
      if(irsym.ne.iqsym)go to 20
      call quadpul(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(ipasym.ne.iqgsym)go to 20

      call findk(orba,orbp,orbq,orbg,iaa,iap,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nqg=idra(iiq,iig,2)
      ini=ntloc(npa,nqg,ipasym,2)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*2+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(two,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+2)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact
 21   continue
 20   continue

c     diagram 17 and 18
      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)
      if(ibsym.ne.igsym)go to 22
      call quadpul(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      if(ipasym.ne.irbsym)go to 22

      call findk(orba,orbp,orbr,orbb,iaa,iap,iar,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrb=idra(iir,iib,2)
      ini=ntloc(npa,nrb,ipasym,2)
      fact=(-1)**(orbb+orbr+2*orbb+2)/(2*2+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(two,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+2)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 23   continue
 22   continue
  2   continue

c     diagram 19
      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)
      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24
      ini=idpa(iiq,iia)
      call quadpul(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
      sum=sum+ti(ini)*dipole
 24   continue
c     diagram 20
      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)
      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25
      ini=idpa(iip,iib)
      call quadpul(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)
      sum=sum-ti(ini)*dipole
 25   continue

c     diagram 21-24
      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)
      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)
      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)
      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)
      if(igsym.ne.irsym)go to 26
c     statement is valid for diagram 21 to 25
      call quadpul(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 21
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

      call findk(orba,orbp,orbq,orbb,iaa,iap,iaq,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nqb,ipasym,2)
      inf=ntloc(nrg,nqb,iqbsym,2)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2.0*2+1)*(2.0*2+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif

 27   continue
c     diagram 22
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(nrg,nqb,iqbsym,2)
      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(2.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2.0*2+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact
 29   continue
      endif
      endif
      enddo
      endif
 28   continue

c     diagram 23
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 30
      do jloop=1,jmax
      l2=jgot(jloop)
      if (l2.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      ini=ntloc(npa,nqb,iqbsym,2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(2.0d0,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2.0*2+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2.0*2+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact
 31   continue
      endif
      endif
      enddo
 30   continue
c     diagram 24
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32
      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(2.0d0,orbr,orbg,al2,orbq,orbb,s6j1)
      call sixj(2.0d0,orba,orbp,al1,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 34   continue
 33   continue
 32   continue
 26   continue

      d(iip,iia)=d(iip,iia)+sum
 1    continue
      return
      end

c*****************************************************************
c     This subroutine calculates the E2 transition O bar diagram *
c     of the type  PP                                            *
c*************************************************************** *

      subroutine quadpulepp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym) go to 1

      call quadpul(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

c-----------------------------------------------------------------
c     diagram 1 = <p|d|q>
c-----------------------------------------------------------------

      d(iip,iiq)=d(iip,iiq)+dipole

      sum=0.0d0

      sum6=0
      sum8=0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      if(iasym.eq.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then

      call quadpul(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      ini=idpa(iip,iia)

c----------------------------------------------------------------
c     diagram 2 = <p|t|a><a|d|q>
c----------------------------------------------------------------

      sum=sum-dipole*ti(ini)

      endif
      endif
      endif

      if(iasym.eq.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then

      call quadpul(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

      sum=sum-dipole*tf(inf)

      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.eq.ipsym)then

      call quadpul(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orbq.eq.orba)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c-----------------------------------------------------------------

      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.eq.irsym)then

      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then

      call quadpul(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      nrg=idra(iir,iig,2)
      nqa=idra(iiq,iia,2)
      inf=ntloc(nrg,nqa,irgsym,2)
      fact=(-1)**(orbg+orbr+2*orbr+2)/(2*two+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c----------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(two,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.eq.irsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call quadpul(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      nrg=idra(iir,iig,2)
      npa=idra(iip,iia,2)
      ini=ntloc(nrg,npa,irgsym,2)
      fact=(-1)**(orbg+orbr+2*orbr+2)/(2*two+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(two,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.eq.iqsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call quadpul(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orba.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)sum=sum+dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

c-----------------------------------------------------------------
c     diagram 12 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

      if(irsym.eq.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then

      call quadpul(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------
c     diagram 13 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------

      if(irsym.eq.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then

      call quadpul(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  5   continue


c------------------------------------------------------------------
c     diagram 14 = <q|t1|a><a|d|g><g|t1|p
c------------------------------------------------------------------

      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      if(igsym.eq.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then

      call quadpul(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)

      sum=sum+dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue


      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 9

      call quadpul(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,two,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,two,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue

c-----------------------------------------------------------------
c     diagram 16 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)

      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,two,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+two)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 14   continue
 13   continue
 12   continue
  9   continue

      if(iasym.ne.idsym)go to 8

      call quadpul(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)

      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8

c-----------------------------------------------------------------
c     digram 17 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(two,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+two+l1)/(2*al1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+two+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,two,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*ti(ini)*tf(inf)*fact

 16   continue
 15   continue
  8   continue

c------------------------------------------------------------------
c     diagram 19 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------

      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nka(iis)

      if(issym.ne.irsym)go to 18

      call quadpul(orbs,orbr,iis,iir,orbs,orbr,kaps,kapr,dipole)

      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18

      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,two,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,two,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+two+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 20   continue

c------------------------------------------------------------------
c     diagram 20 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,two,orbg,al2,orbp,s9j)
c      fact=(-1)**(orbq+orbg+l2+two+2*orbp)*s9j
      fact=(-1)**(orba+orbg+orbr+orbs+l1+l2)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 21   continue
 19   continue
 18   continue

c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)

      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.ne.issym)go to 22

      call quadpul(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,dipole)

      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 24   continue
 23   continue
 22   continue

c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)

      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.ne.issym)go to 25

      call quadpul(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,dipole)

      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 27   continue
 26   continue
 25   continue
 17   continue
  7   continue

       d(iip,iiq)=d(iip,iiq)+sum
       d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+one)

  1   continue

      return
      end

c*************************************************************************

c*************************************************************************
c                                                                        *
c    This subroutine calculates the radial part of the reduced matrix    *
c    element for E2 transition operator                                  *
c                                                                        *
c*************************************************************************

      subroutine Calquad(nbasis)


      implicit real*8 (a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/intquad/rint(MNBAS,MNBAS),kap1,kap2
c      common/intquad/rint(MNBAS,MNBAS),orb1,orb2
      common/orbital_energy/eorb(MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      alpha=137.0359895

      mtp=n

      do i1 = 1, nbasis
      do i2 = 1, nbasis

      rint(i1,i2) = zero

      if(iparity(i1).eq.iparity(i2))then

      ta(1)=zero

      do 1 l=2,mtp
      ta(l)=r(l)**2*rp(l)*((pf(l,i1)*pf(l,i2)+qf(l,i2)*qf(l,i1))
     : + r(l)*(eorb(i1)-eorb(i2))/(alpha*7)*(((kap1-kap2)/3*
     :(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1)))+
     :(pf(l,i1)*qf(l,i2)-qf(l,i1)*pf(l,i2))))
    1 continue
      call quad(result)

      endif

      rint(i1,i2) = result


      enddo
      enddo

      return
      end

c***************************************************************************

      subroutine Calquadvel(nbasis)
c--------------------------------------------------------------------------
c     Fill the array rint()
c--------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/intquadvel/rint(MNBAS,MNBAS),kap1,kap2
      common/orbital_energy/eorb(MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      alpha=137.0359895

      mtp=n

      do i1 = 1, nbasis
c      do i2 = 1, nbasis
      do i2 = 1, i1

      rint(i1,i2) = zero

      if(iparity(i1).eq.iparity(i2))then
      ta(1)=zero
      do 1 l=2,mtp
      ta(l)=r(l)*rp(l)*(2*(pf(l,i1)*qf(l,i2)-qf(l,i1)*pf(l,i2))/
     :((eorb(i2)-eorb(i1))/alpha)-
     :(kap1-kap2)*((r(l)**2/21*(eorb(i1)-eorb(i2)))+(r(l)*alpha/
     :(eorb(i2)-eorb(i1))*
     :(pf(l,i1)*qf(l,i2)+qf(l,i1)*pf(l,i2)))))
    1 continue

      call quad(result)

      endif

      rint(i1,i2) = result
      rint(i2,i1) = rint(i1,i2)

      enddo
      enddo

      return
      end

c***************************************************************************

c***************************************************************************
c                                                                          *
c       This subrourine calculates the radial integral of the reduced      *
c       matrix element of the E2 hyperfine operator                        *
c                                                                          *
c*************************************************************************** 

      subroutine CalcB(nbasis)


      implicit real*8 (a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/intB/rint(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      mtp=n

      do i1 = 1, nbasis
      do i2 = 1, i1

      rint(i1,i2) = zero

      if(iparity(i1).eq.iparity(i2))then

      ta(1)=zero
      do 1 l=2,mtp
      ta(l)=(pf(l,i1)*pf(l,i2)+qf(l,i2)*qf(l,i1))*rp(l)/r(l)**3
    1 continue

      call quad(result)

      endif

      rint(i1,i2) = result
      rint(i2,i1) = rint(i1,i2)

      enddo
      enddo

      return
      end

c*******************************************************************

c********************************************************************
c                                                                   *
c    This subroutine calculates the radial part of the reduced      *
c    matrix elements for the magnetic dipole transition operator    *
c                                                                   *
c********************************************************************


      subroutine Calcm1prop(nbasis)

      implicit real*8 (a-h,o-z)

      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/intm1/rint(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      mtp=n
   

      do i1 = 1, nbasis
      do i2 = 1, i1

      rint(i1,i2) = zero

      if(iparity(i1).eq.iparity(i2))then

      ta(1)=zero

      do 1 l=2,mtp
      ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)*r(l)
    1 continue
      call quad(result)
     
      endif

      rint(i1,i2) = result
      rint(i2,i1) = rint(i1,i2)

      enddo
      enddo

      return
      end

c**************************************************************************

      subroutine Calcm2prop(nbasis)
c--------------------------------------------------------------------------
c     Fill the array rint()
c--------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/intm2/rint(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      mtp=n

      do i1 = 1, nbasis
      do i2 = 1, i1

      rint(i1,i2) = zero
      if(iparity(i1).ne.iparity(i2))then
      ta(1)=zero
      do 1 l=2,mtp
      ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)*r(l)**2
    1 continue

      call quad(result)

      endif

      rint(i1,i2) = result
      rint(i2,i1) = rint(i1,i2)

      enddo
      enddo

      return
      end

      subroutine Calcgfact(nbasis)
c--------------------------------------------------------------------------
c     Fill the array rint()
c--------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/intg/rint(MNBAS,MNBAS)
      common/intgc/rintc(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      mtp=n

      do i1 = 1, nbasis
      do i2 = 1, i1

      rint(i1,i2) = zero
      if(iparity(i1).eq.iparity(i2))then
      ta(1)=zero
      do 1 l=2,mtp
      ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)*r(l)
    1 continue

      call quad(result)

      endif

      rint(i1,i2) = result
      rint(i2,i1) = rint(i1,i2)

      rintc(i1,i2) = zero
      if(iparity(i1).eq.iparity(i2))then
      ta(1)=zero
      do 2 l=2,mtp
      ta(l)=(pf(l,i1)*pf(l,i2)+qf(l,i2)*qf(l,i1))*rp(l)
    2 continue

      call quad(resultc)

      endif

      rintc(i1,i2) = resultc
      rintc(i2,i1) = rintc(i1,i2)

      enddo
      enddo

      return
      end

      subroutine CalcC(nbasis)
c--------------------------------------------------------------------------
c     Fill the array rint()
c--------------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/intoctc/rint(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      mtp=n

      do i1 = 1, nbasis
      do i2 = 1, i1

      rint(i1,i2) = zero
      if(iparity(i1).eq.iparity(i2))then
      ta(1)=zero
      do 1 l=2,mtp
      ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)/(r(l)**4)
    1 continue

      call quad(result)

      endif

      rint(i1,i2) = result
      rint(i2,i1) = rint(i1,i2)

      enddo
      enddo

      return
      end

c***********************************************************************
c     This subroutine calculated the electric dipole O bar diagram of  *
c     of the type HH                                                   *
c***********************************************************************

      subroutine dipolehh(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      jgot(i)=0
      igot(i)=0
      enddo

      nbasis=nocc+nexcit

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)

      if(ibsym.eq.iasym)go to 1

c-----------------------------------------------------------------
c     diagram 1  = <b|d|a>
c-----------------------------------------------------------------

      call e1mat(iib,iia,orbb,orba,kapb,kapa,dipole)

      d(iib,iia)=d(iib,iia)+dipole

      sum=zero

      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

c------------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c------------------------------------------------------------------

      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.ne.ibsym)then

      ini=idpa(iip,iia)

      call e1mat(iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)

      endif
      endif
      endif

c------------------------------------------------------------------
c     diagram 3 = <b|t1|p><p|d|a>
c------------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.ne.iasym)then

      inf=idpa(iip,iib)

      call e1mat(iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)

      endif
      endif
      endif

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.eq.ipbsym) then

c---------------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c---------------------------------------------------------------------

      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.ne.iasym)then

      call e1mat(iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

      call findk(orbb,orbp,orbg,orbr,iab,iap,iag,iar,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      if(orbb.eq.orbp)then
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 5 = <bg|t2|rp><p|d|a><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym.eq.irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.ne.igsym)then

      call e1mat(iir,iig,orbr,orbg,kapr,kapg,dipole)

c-------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c-------------------------------------------------------------------

      ini=idpa(iip,iia)

      call findk(orbb,orbp,orbg,orbr,iab,iap,iag,iar,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npb=idra(iip,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(npb,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c---------------------------------------------------------------------

      ini=idpa(iip,iia)

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c----------------------------------------------------------------
c     diagrams 8 = <<g|d|r><b|t1|p><pr|t2|ag>
c----------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)

      if(ipasym .eq. irgsym) then
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.ne.igsym)then

      call e1mat(iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      call sixj(one,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      inf=idpa(iip,iib)
      ini=ntloc(npg,nra,ipgsym,l1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

      if(ipasym.eq.irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.ne.ibsym)then

      call e1mat(iib,iip,orbb,orbp,kapb,kapp,dipole)

      inf=idpa(iir,iig)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      if(orbp.eq.orba)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      if(orba.eq.orbp) then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif
      endif

 3    continue
 2    continue

      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      if(ipsym.eq.irsym)go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c------------------------------------------------------------------
c     diagram 12 = <b|t1|r><t|d|p><p|t1|a>
c------------------------------------------------------------------

      call e1mat(iir,iip,orbr,orbp,kapr,kapp,dipole)

      inf=idpa(iir,iib)
      ini=idpa(iip,iia)

      sum=sum+dipole*tf(inf)*ti(ini)

 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

c------------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c------------------------------------------------------------------

      if(igsym.ne.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call e1mat(iib,iig,orbb,orbg,kapb,kapg,dipole)

      ini=idpa(iip,iia)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c---------------------------------------------------------------------

      if(igsym.ne.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call e1mat(iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iip,iib)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nak(iid)

      if(igsym.eq.idsym)go to 8

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

      call e1mat(iid,iig,orbd,orbg,kapd,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+1+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 10   continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,one,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11   continue
  9   continue
  8   continue

      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      if(iqsym.eq.irsym)go to 13

      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13

      call e1mat(iiq,iir,orbq,orbr,kapq,kapr,dipole)

c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,one,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue

c----------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c----------------------------------------------------------------------

      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(one,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+one+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo

 14   continue
 13   continue

      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.eq.irsym)go to 16

      call e1mat(iir,iip,orbr,orbp,kapr,kapp,dipole)

c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      ini=ntloc(nqg,npa,ipasym,l1)     
      call sixj(one,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+one)
      fact=s6j*isign/(2*al1+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,one,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+one)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact

  18  continue
  17  continue
  16  continue
  12  continue

      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nak(iid)

      if(idsym.eq.ibsym)go to 21

      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22

      call e1mat(iib,iid,orbb,orbd,kapb,kapd,dipole)


c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 22=   <gd|t2|pr><b|d|d><pr|t2|ag>
c---------------------------------------------------------------------

      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 24   continue
 23   continue
 22   continue
 21   continue

      if(idsym.eq.iasym)go to 25

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26

      call e1mat(iid,iia,orbd,orba,kapd,kapa,dipole)

c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(npd,nrg,irgsym,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l2)
      ini=ntloc(npd,nrg,irgsym,l1)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+1)

 1    continue

      return
      end

c*****************************************************************
c     This subroutine calculated the hyperfine O bar diagram of  *
c     of the type  HP                                            *
c*************************************************************** *


      subroutine dipolehp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nak(iip)

      if(ipsym.eq.iasym)go to 1

c--------------------------------------------------------------
c     diagram 1  = <p|d|a>
c--------------------------------------------------------------

      call e1mat(iip,iia,orbp,orba,kapp,kapa,dipole)

      d(iip,iia)=d(iip,iia)+dipole

      sum=zero


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      if(igsym.eq.irsym)go to 3

      call e1mat(iig,iir,orbg,orbr,kapg,kapr,dipole)

      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.irgsym)go to 401

c--------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c--------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,ipasym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(2*one+1)
      sum=sum+dipole*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c----------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.eq.0)go to 401
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dipole*ti(ini)*fact

 4    continue
 401  continue

c--------------------------------------------------------------------------
c     diagram 4 = <r|t1|a><r|d|g><p|t1|g>
c--------------------------------------------------------------------------

      if(igsym.eq.irsym)go to 3
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3

      call e1mat(iig,iir,orbg,orbr,kapg,kapr,dipole)

      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
  
      sum=sum-dipole*ti(ini)*ti(inj)

 3    continue

c------------------------------------------------------------------------
c     diagram 5 =  <r|t1|a><g|t1|r><p|d|g>
c-------------------------------------------------------------------------

      if(ipsym.eq.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5

      call e1mat(iip,iig,orbp,orbg,kapp,kapg,dipole)

      ini=idpa(iir,iia)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 5    continue

c-------------------------------------------------------------------
c     diagram 6 = <r|d|a><g|t|r><p|t|g>
c-------------------------------------------------------------------

      if(irsym.eq.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6

      ini=idpa(iir,iig)
      inf=idpa(iip,iig)

      call e1mat(iir,iia,orbr,orba,kapr,kapa,dipole)

      sum=sum-dipole*ti(ini)*tf(inf)

 6    continue

      if(igsym.ne.irsym)go to 7
      if(orbg.ne.orbr)go to 7

      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      if(ipsym.eq.iqsym)go to 8

      call e1mat(iip,iiq,orbp,orbq,kapp,kapq,dipole)

      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.iqasym)go to 8

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 7 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      if(orba.eq.orbq)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(2*orbq+1))

      endif
      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 8 =  <g|t1|r><p|d|q><qr|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 9    continue
 8    continue

      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)

      if(iasym.eq.ibsym)go to 10

      call e1mat(iib,iia,orbb,orba,kapb,kapa,dipole)

      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)

      if(irgsym.ne.ipbsym)go to 10

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 9 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbp.eq.orbb)then
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 10 = <g|t1|r><b|d|a><pr|t2|gb>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)then
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 
      endif
 11   continue
 10   continue
  7   continue


      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)

      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipsym.eq.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      if(irgsym.ne.iqasym) goto 13

      call e1mat(iip,iib,orbp,orbb,kapp,kapb,dipole)

c---------------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c---------------------------------------------------------------------

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orba,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue
 14   continue
 13   continue

      if(orbp.ne.orbq)go to 12

      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.eq.iasym)go to 12

c-------------------------------------------------------------------
c     diagram 13 = <gb|t2|rq><pr|t2|bg><q|d|a>
c-------------------------------------------------------------------

      call e1mat(iiq,iia,orbq,orba,kapq,kapa,dipole)

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif


c-------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c-------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 17   continue
 16   continue
 12   continue

      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2

      inf=idpa(iir,iig)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nak(iiq)

      if(irsym.eq.iqsym)go to 20

      call e1mat(iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipasym.ne.iqgsym)go to 20

c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orba,orbp,orbq,orbg,iaa,iap,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(one,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact

 21   continue
 20   continue


      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)

      if(ibsym.eq.igsym)go to 22

      call e1mat(iib,iig,orbb,orbg,kapb,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)

      if(ipasym.ne.irbsym)go to 22

c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbb,iaa,iap,iar,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrb=idra(iir,iib,1)
      ini=ntloc(npa,nrb,ipasym,1)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*1+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(one,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 23   continue
 22   continue
  2   continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------

      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)

      if(iqsym.eq.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24

      ini=idpa(iiq,iia)

      call e1mat(iip,iiq,orbp,orbq,kapp,kapq,dipole)

      sum=sum+ti(ini)*dipole

 24   continue

c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------

      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)

      if(ibsym.eq.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25

      ini=idpa(iip,iib)

      call e1mat(iib,iia,orbb,orba,kapb,kapa,dipole)

      sum=sum-ti(ini)*dipole

 25   continue

      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)

      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)

      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      if(igsym.eq.irsym)go to 26

      call e1mat(iir,iig,orbr,orbg,kapr,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------

      call findk(orba,orbp,orbq,orbb,iaa,iap,iaq,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*1+1)*(2*1+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

 27   continue

c----------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg>
c------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)

      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(one,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*1+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 29   continue

      endif
      endif
      enddo
      endif

 28   continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 30
      do jloop=1,jmax
      l2=jgot(jloop)
      if (l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(one,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2.0*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact

 31   continue

      endif
      endif
      enddo
      
 30   continue

c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32

      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(one,orbr,orbg,al1,orbq,orbb,s6j1)
      call sixj(one,orba,orbp,al2,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 34   continue
 33   continue
 32   continue
 26   continue

      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+1)

 1    continue

      return
      end

c*****************************************************************
c     This subroutine calculated the hyperfine O bar diagram of  *
c     of the type  PP                                            *
c*************************************************************** *

      subroutine dipolepp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nak(iip)

      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)

      if(ipsym.eq.iqsym) go to 1

      call e1mat(iip,iiq,orbp,orbq,kapp,kapq,dipole)

c--------------------------------------------------------------
c     diagram 1 = <p|d|q>
c--------------------------------------------------------------

      d(iip,iiq)=d(iip,iiq)+dipole

      sum=0.0d0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      if(iasym.ne.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then

      call e1mat(iia,iiq,orba,orbq,kapa,kapq,dipole)

      ini=idpa(iip,iia)

c----------------------------------------------------------------
c     diagram 2 = <p|t|a><a|d|q>
c----------------------------------------------------------------

      sum=sum-dipole*ti(ini)

      endif
      endif
      endif

      if(iasym.ne.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then

      call e1mat(iip,iia,orbp,orba,kapp,kapa,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

      sum=sum-dipole*tf(inf)

      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.ne.ipsym)then

      call e1mat(iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
c     sum=sum-dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      if(orbq.eq.orba)then
      sum=sum-dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orba+1))

      endif
      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c---------------------------------------------------------------------

      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba) then
      sum=sum+dipole*ti(ini)*tf(inf)*fact
 
      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.ne.irsym)then

      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then

      call e1mat(iir,iig,orbr,orbg,kapr,kapg,dipole)

      ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      nqa=idra(iiq,iia,1)
      inf=ntloc(nrg,nqa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.ne.irsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call e1mat(iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      npa=idra(iip,iia,1)
      ini=ntloc(nrg,npa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.ne.iqsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call e1mat(iia,iiq,orba,orbq,kapa,kapq,dipole)

      inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orba.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(two*orba+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp) then
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

c------------------------------------------------------------------
c     diagram 12 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

      if(irsym.ne.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then

      call e1mat(iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------
c     diagram 13 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------

      if(irsym.ne.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then

      call e1mat(iip,iir,orbp,orbr,kapp,kapr,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  5   continue

c------------------------------------------------------------------
c     diagram 14 = <q|t1|a><a|d|g><g|t1|p
c------------------------------------------------------------------

      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nak(iia)

      if(igsym.ne.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then

      call e1mat(iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)

      sum=sum+dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nak(iia)

      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nak(iid)

      if(igsym.eq.idsym)go to 9

      call e1mat(iid,iig,orbd,orbg,kapd,kapg,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue

c------------------------------------------------------------------
c     diagram 16 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)

      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,one,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 14   continue
 13   continue
 12   continue
  9   continue

      if(iasym.eq.idsym)go to 8

      call e1mat(iid,iia,orbd,orba,kapd,kapa,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)

      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8

c-----------------------------------------------------------------
c     digram 17 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(one,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+one+l1)/(2*al1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+one+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,one,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 16   continue
 15   continue
  8   continue

c------------------------------------------------------------------
c     diagram 19 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------

      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nak(iis)

      if(issym.eq.irsym)go to 18

      call e1mat(iis,iir,orbs,orbr,kaps,kapr,dipole)

      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18

      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,one,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 20   continue

c------------------------------------------------------------------
c     diagram 20 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,one,orbg,al2,orbp,s9j)
c      fact=(-1)**(orbq+orbg+l2+one+2*orbp)*s9j
      fact=(-1)**(orba+orbg+orbr+orbs+l1+l2)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 21   continue
 19   continue
 18   continue

c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)

      if (irgsym .ne. iassym) goto 22
      if (irgsym .ne. iqasym) goto 22
      if(ipsym.eq.issym)go to 22

      call e1mat(iip,iis,orbp,orbs,kapp,kaps,dipole)

      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 24   continue
 23   continue
 22   continue

c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)

      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.eq.issym)go to 25

      call e1mat(iis,iiq,orbs,orbq,kaps,kapq,dipole)

      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 27   continue
 26   continue
 25   continue
 17   continue
  7   continue

      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+1)

  1   continue

      return
      end


c*********************************************************************


c*********************************************************************
c      The reduced matrix element of the electric dipole operator is *
c      calculated                                                    *
c*********************************************************************

      subroutine e1mat(i1,i2,jorb1,jorb2,kap1,kap2,e1m)

      implicit double precision (a-h,o-z)
      double precision j1,j2,jorb1,jorb2

      common/cons/zero,half,tenth,one,two,three,ten
      common/orbital_energy/eorb(MNBAS)
      common/inte1mat/rint1(MNBAS,MNBAS),rint2(MNBAS,MNBAS)
     :,rint3(MNBAS,MNBAS)

      tj1=idint(2.0*jorb1)
      tj2=idint(2.0*jorb2)
      factor=(-1)**(jorb1+0.5)
      j1=tj1/two
      j2=tj2/two
      w3j1=dr(j1,one,j2,half,zero,-half)
      e1rad=rint1(i1,i2)

      e1=dsqrt(tj1+one)*dsqrt(tj2+one)*factor*w3j1

      e1m=e1*e1rad

      return
      end

c*******************************************************************

c*******************************************************************
c                                                                  *
c  This subroutine calculates the radial part of the reduced       *
c  metrix elements of the electric dipole operator                 *
c                                                                  *
c*******************************************************************

      subroutine Calce1mat(nbasis)


      implicit real*8 (a-h,o-z)

      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/inte1mat/rint1(MNBAS,MNBAS),rint2(MNBAS,MNBAS)
     :,rint3(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

      do i1 = 1, nbasis
      do i2 = 1, i1

      rint1(i1,i2) = 0.0d0

      if(iparity(i1).ne.iparity(i2))then
      
      rint1(i1,i2) = frint(i1,i2,1)
      rint1(i2,i1) = rint1(i1,i2)

      endif
      enddo
      enddo

      return
      end


c***********************************************************************


c***********************************************************************
c                                                                      *
c                                                                      *
c                                                                      *
c          The value of this  frint is an approximation to:            *
c                                                                      *
c              k                                                       *
c         i ( r  *  ( p (r)*p (r) + q (r)*q (r) ; 0 to infinity)       *
c                      i     j       i     j                           *
c                                                                      *
c   where   i ( g(r) ; range )  denotes  the  integral  of g(r) over   *
c   range.                                                             *
c                                                                      *
c   subroutines called: quad.                                          *
c                                                                      *
c   written by farid a parpia, at oxford   last updated: 18 may 1989   *
c                                                                      *
c***********************************************************************


      double precision function frint (i,j,k)

      implicit real*8 (a-h,o-z)

      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)

      mtp = n

      ta(1) = zero

      do 1 l = 2,mtp
      ta(l)=(r(l)**k)*(pf(l,i)*pf(l,j)+qf(l,i)*qf(l,j))*rp(l)
    1 continue
      call quad (result)

      frint = result

      return
      end

c*******************************************************************

      
c***********************************************************************
c     This subroutine calculates the NSI PNC O bar diagram of          *
c     of the type HH                                                   *
c***********************************************************************

      
       subroutine nsihh(d,ti,tf)

       implicit real*8 (a-h,o-z)
       common/jvalue/orbc(MNOCC),orbe(MNEXC)
       common/index/ke(MNEXC),kc(MNOCC)
       common/info/nocc,nexcit
       common/parity/isymc(MNOCC),isyme(MNEXC)
       common/symmetry/mtbl(MNS,MNS)
       common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
       common/orbital_energy/eorb(MNBAS)
       common/idra1/idra(MNBAS,MNBAS,0:MXV)
       common/idpa1/idpa(MNBAS,MNBAS)
       common/kappa/nak(MNBAS)
       common/inparm/nparm
       common/pncp/amass,qw


       dimension kgot(MXV),igot(MXV),jgot(MXV)
       dimension ti(MDIM),d(MNBAS,MNBAS),tf(MDIM)

       data zero,half,one,two/0.0,0.5,1.0,2.0/

       do i=1,MXV
       kgot(i)=0
       jgot(i)=0
       igot(i)=0
       enddo

       nbasis=nocc+nexcit

       do 1 ia=1,nocc
       iasym=isymc(ia)
       iia=kc(ia)
       orba=orbc(ia)
       iaa=iqc(ia)
       kapa=nak(iia)

       do 1 ib=1,ia
       ibsym=isymc(ib)
       iib=kc(ib)
       iab=iqc(ib)
       orbb=orbc(ib)
       kapb=nak(iib)

       if(ibsym.eq.iasym)go to 1

c-----------------------------------------------------------------
c     diagram 1  = <b|d|a>
c-----------------------------------------------------------------

       call nsi(iib,iia,orbb,orba,kapb,kapa,dipole)
       d(iib,iia)=d(iib,iia)+dipole

       sum=zero

       do 2 ip=1,nexcit
       ipsym=isyme(ip)
       iip=ke(ip)
       orbp=orbe(ip)
       iap=iqe(ip)
       kapp=nak(iip)

c------------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c------------------------------------------------------------------

       if(ipsym.eq.iasym)then
       if(orbp.eq.orba)then
       if(ipsym.ne.ibsym)then


       ini=idpa(iip,iia)

       call nsi(iib,iip,orbb,orbp,kapb,kapp,dipole)
       sum=sum+dipole*ti(ini)

       endif
       endif
       endif

c------------------------------------------------------------------
c     diagram 3 = <b|t1|p><p|d|a>
c------------------------------------------------------------------

       if(orbp.eq.orbb)then
       if(ipsym.eq.ibsym)then
       if(ipsym.ne.iasym)then

       inf=idpa(iip,iib)
       call nsi(iip,iia,orbp,orba,kapp,kapa,dipole)
       sum=sum+dipole*tf(inf)

       endif
       endif
       endif

       do 3 ir=1,nexcit
       irsym=isyme(ir)
       iir=ke(ir)
       orbr=orbe(ir)
       iar=iqe(ir)
       kapr=nak(iir)

       do 3 ig=1,nocc
       igsym=isymc(ig)
       iig=kc(ig)
       orbg=orbc(ig)
       iag=iqc(ig)
       kapg=nak(iig)

       irgsym=mtbl(irsym,igsym)
       ipbsym=mtbl(ipsym,ibsym)

c---------------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c---------------------------------------------------------------------

       if(irgsym.eq.ipbsym) then
       if(orbg.eq.orbr)then
       if(irsym.eq.igsym)then
       if(ipsym.ne.iasym)then

       ini=idpa(iir,iig)

       call findk(orbb,orbp,orbg,orbr,iab,iap,iag,iar,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       npb=idra(iip,iib,0)
       nrg=idra(iir,iig,0)
       inf=ntloc(npb,nrg,irgsym,0)
       call nsi(iip,iia,orbp,orba,kapp,kapa,dipole)
c      sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
       sum=sum+dipole*ti(ini)*tf(inf)*
     : dsqrt((two*orbg+1)/(two*orbp+1))

       endif
       enddo
       endif


c-------------------------------------------------------------------
c     diagrams 5 = <bg|t2|rp><p|d|a><r|t1|g>
c-------------------------------------------------------------------

       call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
       if(imax.ne.0)then
       if(orbb.eq.orbp)then
       do iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npg=idra(iip,iig,l1)
       nrb=idra(iir,iib,l1)
       ipgsym=mtbl(ipsym,igsym)
       inf=ntloc(npg,nrb,ipgsym,l1)
       fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
       sum=sum-dipole*ti(ini)*tf(inf)*fact

       enddo
       endif
       endif

       endif
       endif
       endif
       endif

       if(ipbsym.eq.irgsym) then
       if(orbp.eq.orba)then
       if(ipsym.eq.iasym)then
       if(irsym.ne.igsym)then

       call nsi(iir,iig,orbr,orbg,kapr,kapg,dipole)

c-------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c-------------------------------------------------------------------

       ini=idpa(iip,iia)

       call findk(orbb,orbp,orbg,orbr,iab,iap,iag,iar,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       npb=idra(iip,iib,0)
       nrg=idra(iir,iig,0)
       inf=ntloc(npb,nrg,irgsym,0)
       fact=(-1)**(orbg+orbr+0+2*orbr)/(two*zero+one)
       sum=sum+dipole*ti(ini)*tf(inf)*fact

       endif
       enddo
       endif

c---------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c---------------------------------------------------------------------

       ini=idpa(iip,iia)

       call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npg=idra(iip,iig,l1)
       nrb=idra(iir,iib,l1)
       ipgsym=mtbl(ipsym,igsym)
       inf=ntloc(npg,nrb,ipgsym,l1)
       call sixj(zero,orbr,orbg,al1,orba,orbb,s6j)
       fact=(-1)**(orbr+orbg+0)*s6j
       sum=sum-dipole*ti(ini)*tf(inf)*fact

       enddo
       endif

       endif
       endif
       endif
       endif

c----------------------------------------------------------------
c     diagram 8 =  <g|d|r><b|t1|p><pr|t2|ag>
c----------------------------------------------------------------

       ipasym=mtbl(ipsym,iasym)

       if(ipasym .eq. irgsym) then
       if(orbp.eq.orbb)then
       if(ipsym.eq.ibsym)then
       if(irsym.ne.igsym)then

       call nsi(iig,iir,orbg,orbr,kapg,kapr,dipole)

       inf=idpa(iip,iib)

       call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       npa=idra(iip,iia,0)
       nrg=idra(iir,iig,0)
       ini=ntloc(npa,nrg,irgsym,0)
       fact=(-1)**(orbg+orbr+0+2*orbr)/(two*zero+one)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif


c-------------------------------------------------------------------
c     diagram 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

       call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npg=idra(iip,iig,l1)
       nra=idra(iir,iia,l1)
       ipgsym=mtbl(ipsym,igsym)
       call sixj(zero,orbr,orbg,al1,orbb,orba,s6j)
       fact=(-1)**(orbr+orbg+0)*s6j
       inf=idpa(iip,iib)
       ini=ntloc(npg,nra,ipgsym,l1)
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       enddo
       endif

       endif
       endif
       endif
       endif

c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

       if(ipasym.eq.irgsym) then
       if(orbg.eq.orbr)then
       if(irsym.eq.igsym)then
       if(ipsym.ne.ibsym)then

       call nsi(iib,iip,orbb,orbp,kapb,kapp,dipole)

       inf=idpa(iir,iig)

       call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       npa=idra(iip,iia,0)
       nrg=idra(iir,iig,0)
       ini=ntloc(npa,nrg,irgsym,0)
c      sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
       sum=sum+dipole*tf(inf)*ti(ini)*
     : dsqrt((two*orbg+1)/(two*orbp+1))

       endif
       enddo
       endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------

       call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
       if(imax.ne.0)then
       if(orba.eq.orbp) then
       do iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npg=idra(iip,iig,l1)
       nra=idra(iir,iia,l1)
       ipgsym=mtbl(ipsym,igsym)
       ini=ntloc(npg,nra,ipgsym,l1)
       inf=idpa(iir,iig)
       fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       enddo
       endif
       endif

       endif
       endif
       endif
       endif

 3     continue
 2     continue

       do 4 ip=1,nexcit
       ipsym=isyme(ip)
       iip=ke(ip)
       orbp=orbe(ip)
       iap=iqe(ip)
       kapp=nak(iip)

       do 5 ir=1,nexcit
       irsym=isyme(ir)
       iir=ke(ir)
       orbr=orbe(ir)
       iar=iqe(ir)
       kapr=nak(iir)

       if(ipsym.eq.irsym)go to 5
       if(ipsym.ne.iasym)go to 5
       if(orbp.ne.orba)go to 5
       if(irsym.ne.ibsym)go to 5
       if(orbr.ne.orbb)go to 5

c------------------------------------------------------------------
c     diagram 12 = <b|t1|r><r|d|p><p|t1|a>
c------------------------------------------------------------------

       call nsi(iir,iip,orbr,orbp,kapr,kapp,dipole)

       inf=idpa(iir,iib)
       ini=idpa(iip,iia)
       sum=sum+dipole*tf(inf)*ti(ini)

 5     continue

       do 6 ig=1,nocc
       igsym=isymc(ig)
       iig=kc(ig)
       orbg=orbc(ig)
       iag=iqc(ig)
       kapg=nak(iig)

c------------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c------------------------------------------------------------------

       if(igsym.ne.ibsym)then
       if(ipsym.eq.iasym)then
       if(orbp.eq.orba)then
       if(ipsym.eq.igsym)then
       if(orbp.eq.orbg)then

       call nsi(iib,iig,orbb,orbg,kapb,kapg,dipole)

       ini=idpa(iip,iia)
       inf=idpa(iip,iig)

       sum=sum-dipole*ti(ini)*tf(inf)

       endif
       endif
       endif
       endif
       endif

c---------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c---------------------------------------------------------------------

       if(igsym.ne.iasym)then
       if(ipsym.eq.ibsym)then
       if(orbp.eq.orbb)then
       if(ipsym.eq.igsym)then
       if(orbp.eq.orbg)then

       call nsi(iig,iia,orbg,orba,kapg,kapa,dipole)

       ini=idpa(iip,iig)
       inf=idpa(iip,iib)

       sum=sum-dipole*ti(ini)*tf(inf)

       endif
       endif
       endif
       endif
       endif

 6     continue
 4     continue

       do 7 ip=1,nexcit
       ipsym=isyme(ip)
       iip=ke(ip)
       orbp=orbe(ip)
       iap=iqe(ip)
       kapp=nak(iip)

       do 7 ir=1,nexcit
       irsym=isyme(ir)
       iir=ke(ir)
       orbr=orbe(ir)
       iar=iqe(ir)
       kapr=nak(iir)

       do 7 ig=1,nocc
       igsym=isymc(ig)
       iig=kc(ig)
       orbg=orbc(ig)
       iag=iqc(ig)
       kapg=nak(iig)

       do 8 id=1,nocc
       idsym=isymc(id)
       iid=kc(id)
       orbd=orbc(id)
       iad=iqc(id)
       kapd=nak(iid)

       if(igsym.eq.idsym)go to 8

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

       call nsi(iid,iig,orbd,orbg,kapd,kapg,dipole)

       ipasym=mtbl(ipsym,iasym)
       ipbsym=mtbl(ipsym,ibsym)
       irgsym=mtbl(irsym,igsym)
       irdsym=mtbl(irsym,idsym)

       if(ipasym.ne.irdsym)go to 8
       if(ipbsym.ne.irgsym)go to 8

       call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
       if(imax.eq.0)go to 8
       do 9 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npb=idra(iip,iib,l1)
       nrg=idra(iir,iig,l1)
       inf=ntloc(npb,nrg,ipbsym,l1)

       call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
       if(jmax.eq.0)go to 9
       do 10 jloop=1,jmax
       l2=jgot(jloop)
       al2=float(l2)
       npa=idra(iip,iia,l2)
       nrd=idra(iir,iid,l2)
       ini=ntloc(npa,nrd,ipasym,l2)
       call sixj(orbg,orbd,zero,al2,al1,orbr,s6j1)
       call sixj(orbb,orba,zero,al2,al1,orbp,s6j2)
       fact=(-1)**(orba+orbb+0+2*orbg)*s6j1*s6j2
       sum=sum-dipole*tf(inf)*ti(ini)*fact

 10    continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

       call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
       if(jmax.eq.0)go to 9
       do 11 jloop=1,jmax
       l2=jgot(jloop)
       al2=float(l2)
       npd=idra(iip,iid,l2)
       nra=idra(iir,iia,l2)
       ipdsym=mtbl(ipsym,idsym)
       ini=ntloc(npd,nra,ipdsym,l2)
       call ninej(al1,orbp,orbb,orbg,orbd,zero,orbr,al2,orba,s9j)
       fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
       sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11    continue
  9    continue
  8    continue

       do 12 iq=1,nexcit
       iqsym=isyme(iq)
       iiq=ke(iq)
       orbq=orbe(iq)
       iaq=iqe(iq)
       kapq=nak(iiq)

       if(iqsym.eq.irsym)go to 13

       iqgsym=mtbl(iqsym,igsym)
       irgsym=mtbl(irsym,igsym)
       ipasym=mtbl(ipsym,iasym)
       ipbsym=mtbl(ipsym,ibsym)

       if(irgsym.ne.ipasym)go to 13
       if(iqgsym.ne.ipbsym)go to 13

       call nsi(iiq,iir,orbq,orbr,kapq,kapr,dipole)

c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------

       call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
       if(imax.eq.0)go to 13
       do 14 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npb=idra(iip,iib,l1)
       nqg=idra(iiq,iig,l1)
       inf=ntloc(npb,nqg,ipbsym,l1)

       call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
       if(jmax.eq.0)go to 14
       do 15 jloop=1,jmax
       l2=jgot(jloop)
       al2=float(l2)
       nrg=idra(iir,iig,l2)
       npa=idra(iip,iia,l2)
       ini=ntloc(npa,nrg,ipasym,l2)
       call sixj(orbq,orbr,zero,al2,al1,orbg,s6j1)
       call sixj(orbb,orba,zero,al2,al1,orbp,s6j2)
       fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
       sum=sum+dipole*tf(inf)*ti(ini)*fact

 15    continue

c----------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c----------------------------------------------------------------------

       call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
       if(jmax.eq.0)go to 14
       do jloop=1,jmax
       l2=jgot(jloop)
       al2=float(l2)
       nra=idra(iir,iia,l2)
       npg=idra(iip,iig,l2)
       ipgsym=mtbl(ipsym,igsym)
       ini=ntloc(npg,nra,ipgsym,l2)
       call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
       call sixj(zero,orbq,orbr,al2,orba,orbb,s6j2)
       fact=(-1)**(orbr+orbq+zero+orbg+orbp+l2)*s6j1*s6j2
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       enddo

 14    continue
 13    continue

       iqgsym=mtbl(iqsym,igsym)
       irbsym=mtbl(irsym,ibsym)
       ipasym=mtbl(ipsym,iasym)

       if(ipasym.ne.iqgsym)go to 16
       if(irbsym.ne.iqgsym)go to 16
       if(ipsym.eq.irsym)go to 16

       call nsi(iir,iip,orbr,orbp,kapr,kapq,dipole)

c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------

       call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
       if(jmax.eq.0)go to 16
       do 17 jloop=1,jmax
       l1=jgot(jloop)
       al1=float(l1)
       call sixj(zero,orbp,orbr,al1,orbb,orba,s6j)
       isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+zero)
       fact=s6j*isign/(2*al1+1)
       nqg=idra(iiq,iig,l1)
       nrb=idra(iir,iib,l1)

       call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       k1=igot(iloop)
       if(k1.eq.l1)then
       nqg=idra(iiq,iig,l1)
       npa=idra(iip,iia,l1)
       inf=ntloc(nqg,nrb,irbsym,l1)
       ini=ntloc(nqg,npa,ipasym,l1)
       sum=sum+tf(inf)*ti(ini)*dipole*fact

       endif
       enddo
       endif


c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

       call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
       if(imax.eq.0)go to 17
       do 18 iloop=1,imax
       l2=igot(iloop)
       al2=float(l2)
       npg=idra(iip,iig,l2)
       nqa=idra(iiq,iia,l2)
       iqasym=mtbl(iqsym,iasym)
       ini=ntloc(npg,nqa,iqasym,l2)
       call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
       call sixj(orba,orbb,zero,orbr,orbp,al1,s6j2)
       isign=(-1)**(orbg+orbq+al1+orbp+orbr+zero)
       fact=isign*s6j1*s6j2
       sum=sum-tf(inf)*ti(ini)*dipole*fact

 18    continue
 17    continue
 16    continue
 12    continue

       do 20 id=1,nocc
       idsym=isymc(id)
       iid=kc(id)
       orbd=orbc(id)
       iad=iqc(id)
       kapd=nak(iid)

       if(idsym.eq.ibsym)go to 21

       irgsym=mtbl(irsym,igsym)
       ipasym=mtbl(ipsym,iasym)
       irdsym=mtbl(irsym,idsym)
       ipdsym=mtbl(ipsym,idsym)

       if(irgsym.ne.ipasym)go to 22
       if(irgsym.ne.ipdsym)go to 22
       if(orba.ne.orbd)go to 22

       call nsi(iib,iid,orbb,orbd,kapb,kapd,dipole)

c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

       call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
       if(imax.eq.0)go to 22
       do 23 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npa=idra(iip,iia,l1)
       nrg=idra(iir,iig,l1)
       npd=idra(iip,iid,l1)
       ini=ntloc(npa,nrg,ipasym,l1)

       call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
       if(jmax.ne.0)then
       do jloop=1,jmax
       k1=jgot(jloop)
       if(k1.eq.l1)then
       nrg=idra(iir,iig,l1)
       npd=idra(iip,iid,l1)
       inf=ntloc(npd,nrg,irgsym,l1)
       isign=(-1)**(orbg+orbr+orbd+orbp)
       fact=isign/((2*l1+1)*(2*orba+1))
       sum=sum-ti(ini)*tf(inf)*dipole*fact

       endif
       enddo
       endif

c---------------------------------------------------------------------
c     diagram 22=   <gd|t2|pr><b|d|d><pr|t2|ag>
c---------------------------------------------------------------------

       call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
       if(jmax.eq.0)go to 23
       do 24 jloop=1,jmax
       l2=jgot(jloop)
       al2=float(l2)
       npa=idra(iip,iia,l1)
       nrg=idra(iir,iig,l1)
       nrd=idra(iir,iid,l2)
       npg=idra(iip,iig,l2)
       ipgsym=mtbl(ipsym,igsym)
       ini=ntloc(npa,nrg,ipasym,l1)
       inf=ntloc(npg,nrd,ipgsym,l2)
       call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
       isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
       fact=isign*s6j/(2*orba+1)
       sum=sum+ti(ini)*tf(inf)*dipole*fact

 24    continue
 23    continue
 22    continue
 21    continue

       if(idsym.eq.iasym)go to 25

       irgsym=mtbl(irsym,igsym)
       ipbsym=mtbl(ipsym,ibsym)
       ipdsym=mtbl(ipsym,idsym)

       if(irgsym.ne.ipdsym)go to 26
       if(irgsym.ne.ipbsym)go to 26
       if(orbb.ne.orbd)go to 26

       call nsi(iid,iia,orbd,orba,kapd,kapa,dipole)

c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------

       call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
       if(imax.eq.0)go to 26
       do 27 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npb=idra(iip,iib,l1)
       nrg=idra(iir,iig,l1)
       npd=idra(iip,iid,l1)

       call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
       if(jmax.ne.0)then
       do jloop=1,jmax
       k1=jgot(jloop)
       if(k1.eq.l1)then
       npb=idra(iip,iib,l1)
       nrg=idra(iir,iig,l1)
       inf=ntloc(npb,nrg,ipbsym,l1)
       ini=ntloc(npd,nrg,irgsym,l1)
       isign=(-1)**(orbg+orbr+orbb+orbp)
       fact=isign/((2*l1+1)*(2*orbb+1))
       sum=sum-tf(inf)*ti(ini)*dipole*fact

       endif
       enddo
       endif

c------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------

       call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
       if(jmax.eq.0)go to 27
       do 28 jloop=1,jmax
       l2=jgot(jloop)
       al2=float(l2)
       npg=idra(iip,iig,l2)
       nrb=idra(iir,iib,l2)
       nrg=idra(iir,iig,l1)
       npd=idra(iip,iid,l1)
       ipgsym=mtbl(ipsym,igsym)
       inf=ntloc(npg,nrb,ipgsym,l2)
       ini=ntloc(npd,nrg,irgsym,l1)
       call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
       isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
       fact=isign*s6j/(2*orbb+1)
       sum=sum+tf(inf)*ti(ini)*dipole*fact

 28    continue
 27    continue
 26    continue
 25    continue
 20    continue
  7    continue

       d(iib,iia)=d(iib,iia)+sum
       d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+1)

  1    continue

       return
       end


c*****************************************************************
c     This subroutine calculates the NSI PNC O bar diagram of    *
c     of the type  HP                                            *
c*****************************************************************



       subroutine nsihp(d,ti,tf)

       implicit real*8 (a-h,o-z)
       common/jvalue/orbc(MNOCC),orbe(MNEXC)
       common/index/ke(MNEXC),kc(MNOCC)
       common/info/nocc,nexcit
       common/parity/isymc(MNOCC),isyme(MNEXC)
       common/symmetry/mtbl(MNS,MNS)
       common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
       common/orbital_energy/eorb(MNBAS)
       common/idra1/idra(MNBAS,MNBAS,0:MXV)
       common/idpa1/idpa(MNBAS,MNBAS)
       common/kappa/nak(MNBAS)
       common/inparm/nparm
       common/pncp/amass,qw


       dimension kgot(MXV),igot(MXV),jgot(MXV)
       dimension d(MNBAS,MNBAS)
       dimension ti(MDIM),tf(MDIM)

       data zero,half,one,two/0.0,0.5,1.0,2.0/

       do i=1,MXV
       kgot(i)=0
       igot(i)=0
       jgot(i)=0
       enddo

       do 1 ia=1,nocc
       iasym=isymc(ia)
       iia=kc(ia)
       orba=orbc(ia)
       iaa=iqc(ia)
       kapa=nak(iia)

       do 1 ip=1,nexcit
       ipsym=isyme(ip)
       iip=ke(ip)
       iap=iqe(ip)
       orbp=orbe(ip)
       kapp=nak(iip)

       if(ipsym.eq.iasym)go to 1

c--------------------------------------------------------------
c     diagram 1  = <p|d|a>
c--------------------------------------------------------------

       call nsi(iip,iia,orbp,orba,kapp,kapa,dipole)

       d(iip,iia)=d(iip,iia)+dipole

       sum=zero


       do 2 ir=1,nexcit
       irsym=isyme(ir)
       iir=ke(ir)
       orbr=orbe(ir)
       iar=iqe(ir)
       kapr=nak(iir)

       do 2 ig=1,nocc
       igsym=isymc(ig)
       iig=kc(ig)
       orbg=orbc(ig)
       iag=iqc(ig)
       kapg=nak(iig)

       if(igsym.eq.irsym)go to 3

       call nsi(iig,iir,orbg,orbr,kapg,kapr,dipole)

       ipasym=mtbl(ipsym,iasym)
       irgsym=mtbl(irsym,igsym)

       if(ipasym.ne.irgsym)go to 401

c--------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c--------------------------------------------------------------------

       call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       npa=idra(iip,iia,0)
       nrg=idra(iir,iig,0)
       ini=ntloc(npa,nrg,ipasym,0)
       fact=(-1)**(orbg+orbr+0+2*orbr)/(2*zero+1)
       sum=sum+dipole*ti(ini)*fact

       endif
       enddo
       endif

c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c----------------------------------------------------------------------

       call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
       if(imax.eq.0)go to 401
       do 4 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npg=idra(iip,iig,l1)
       nra=idra(iir,iia,l1)
       ipgsym=mtbl(ipsym,igsym)
       ini=ntloc(npg,nra,ipgsym,l1)
       call sixj(zero,orbr,orbg,al1,orbp,orba,s6j)
       fact=(-1)**(orbg+orbr+zero)*s6j
       sum=sum-dipole*ti(ini)*fact

  4    continue
 401   continue

c---------------------------------------------------------------------
c     diagram 8 = <r|t1|a><g|d|r><p|t1|g>
c---------------------------------------------------------------------

       if(ipsym.ne.igsym)go to 3
       if(orbp.ne.orbg)go to 3
       if(irsym.ne.iasym)go to 3
       if(orbr.ne.orba)go to 3
       ini=idpa(iir,iia)
       inj=idpa(iip,iig)
       call nsi(iig,iir,orbg,orbr,kapg,kapr,dipole)
       sum=sum-dipole*ti(ini)*ti(inj)

 3     continue

c----------------------------------------------------------------------
c     diagram 9 = <r|t1|a><g|t1|r><p|d|g>
c----------------------------------------------------------------------

       if(ipsym.eq.igsym)go to 5
       if(irsym.ne.iasym)go to 5
       if(orbr.ne.orba)go to 5
       if(igsym.ne.irsym)go to 5
       if(orbr.ne.orbg)go to 5

       call nsi(iip,iig,orbp,orbg,kapp,kapg,dipole)

       ini=idpa(iir,iia)
       inf=idpa(iir,iig)

       sum=sum-dipole*ti(ini)*tf(inf)

 5     continue

c---------------------------------------------------------------------
c     diagram 10 = <r|d|a><g|t|r><p|t|g>
c---------------------------------------------------------------------

       if(irsym.eq.iasym)go to 6
       if(irsym.ne.igsym)go to 6
       if(orbr.ne.orbg)go to 6
       if(igsym.ne.ipsym)go to 6
       if(orbg.ne.orbp)go to 6

       ini=idpa(iir,iig)
       inf=idpa(iip,iig)

       call nsi(iir,iia,orbr,orba,kapr,kapa,dipole)

       sum=sum-dipole*ti(ini)*tf(inf)

 6     continue

       if(igsym.ne.irsym)go to 7
       if(orbg.ne.orbr)go to 7

       do 8 iq=1,nexcit
       iqsym=isyme(iq)
       iiq=ke(iq)
       orbq=orbe(iq)
       iaq=iqe(iq)
       kapq=nak(iiq)

       if(ipsym.eq.iqsym)go to 8

       call nsi(iip,iiq,orbp,orbq,kapp,kapq,dipole)

       iqasym=mtbl(iqsym,iasym)
       irasym=mtbl(irsym,iasym)
       irgsym=mtbl(irsym,igsym)
       if(irgsym.ne.iqasym)go to 8

       inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 4 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

       call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       nrg=idra(iir,iig,0)
       nqa=idra(iiq,iia,0)
       ini=ntloc(nrg,nqa,iqasym,0)
c      sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
       sum=sum+dipole*tf(inf)*ti(ini)*
     : dsqrt((2*orbg+1)/(2*orbq+1))

       endif
       enddo
       endif

c---------------------------------------------------------------------
c     diagram 5 = <g|t1|r><p|d|q><qr|t2|ga>
c--------------------------------------------------------------------

       call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
       if (imax.eq.0)go to 8
       do 9 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nqg=idra(iiq,iig,l1)
       nra=idra(iir,iia,l1)
       ini=ntloc(nqg,nra,irasym,l1)
       if(orba.ne.orbq)go to 9
       fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
       sum=sum-dipole*tf(inf)*ti(ini)*fact

 9     continue
 8     continue

       do 10 ib=1,nocc
       ibsym=isymc(ib)
       iib=kc(ib)
       orbb=orbc(ib)
       iab=iqc(ib)
       kapb=nak(iib)

       if(iasym.eq.ibsym)go to 10

       call nsi(iib,iia,orbb,orba,kapb,kapa,dipole)
       ipbsym=mtbl(ipsym,ibsym)
       ipgsym=mtbl(ipsym,igsym)
       if(irgsym.ne.ipbsym)go to 10
       inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 6 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

       call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       nrg=idra(iir,iig,0)
       npb=idra(iip,iib,0)
       ini=ntloc(nrg,npb,irgsym,0)
       if(orbp.eq.orbb)then
       sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))

       endif
       endif
       enddo
       endif

c-------------------------------------------------------------------
c     diagram 7 = <g|t1|r><b|d|a><pr|t2|gb>
c-------------------------------------------------------------------

       call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
       if (imax.eq.0)go to 10
       if(orbp.ne.orbb)go to 10
       do 11 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npg=idra(iip,iig,l1)
       nrb=idra(iir,iib,l1)
       ini=ntloc(npg,nrb,ipgsym,l1)
       fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

 11    continue
 10    continue
  7    continue



       do 12 iq=1,nexcit
       iqsym=isyme(iq)
       iiq=ke(iq)
       orbq=orbe(iq)
       iaq=iqe(iq)
       kapq=nak(iiq)

       do 12 ib=1,nocc
       ibsym=isymc(ib)
       iib=kc(ib)
       orbb=orbc(ib)
       iab=iqc(ib)
       kapb=nak(iib)

       irgsym=mtbl(irsym,igsym)
       iqbsym=mtbl(iqsym,ibsym)
       iqasym=mtbl(iqsym,iasym)

       if(ipsym.eq.ibsym)go to 13
       if(orbb.ne.orba)go to 13
       if(irgsym.ne.iqbsym)go to 13
       if(irgsym.ne.iqasym) goto 13

       call nsi(iip,iib,orbp,orbb,kapp,kapb,dipole)

c---------------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c---------------------------------------------------------------------

       call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
       if(imax.eq.0)go to 13
       do 14 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nrg=idra(iir,iig,l1)
       nqb=idra(iiq,iib,l1)
       inf=ntloc(nrg,nqb,iqbsym,l1)

       call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iaa,jgot,jmax)
       if(jmax.ne.0)then
       do jloop=1,jmax
       k1=jgot(jloop)
       if(k1.eq.l1)then
       nrg=idra(iir,iig,l1)
       nqa=idra(iiq,iia,l1)
       ini=ntloc(nrg,nqa,irgsym,l1)
       fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif

c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

       call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
       if(kmax.eq.0)go to 14
       do 15 jloop=1,kmax
       l2=kgot(jloop)
       al2=float(l2)
       nra=idra(iir,iia,l2)
       nqg=idra(iiq,iig,l2)
       iqgsym=mtbl(iqsym,igsym)
       ini=ntloc(nra,nqg,iqgsym,l2)
       call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
       fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

 15    continue
 14    continue
 13    continue

       if(orbp.ne.orbq)go to 12

       iqbsym=mtbl(iqsym,ibsym)
       irgsym=mtbl(irsym,igsym)
       ipgsym=mtbl(ipsym,igsym)
       ipbsym=mtbl(ipsym,ibsym)

       if(irgsym.ne.ipbsym) goto 12
       if(irgsym.ne.iqbsym)go to 12
       if(iqsym.eq.iasym)go to 12

c-------------------------------------------------------------------
c     diagram 13 = <gb|t2|rq><pr|t2|bg><q|d|a>
c-------------------------------------------------------------------

       call nsi(iiq,iia,orbq,orba,kapq,kapa,dipole)

       call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
       if(imax.eq.0)go to 12
       do 16 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nrg=idra(iir,iig,l1)
       nqb=idra(iiq,iib,l1)
       npb=idra(iip,iib,l1)
       inf=ntloc(nrg,nqb,irgsym,l1)

       call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
       if(jmax.ne.0)then
       do jloop=1,jmax
       k1=jgot(jloop)
       if(k1.eq.l1)then
       nrg=idra(iir,iig,l1)
       npb=idra(iip,iib,l1)
       ini=ntloc(nrg,npb,irgsym,l1)
       fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif

c-------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c-------------------------------------------------------------------

       call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
       if(kmax.eq.0)go to 16
       do 17 jloop=1,kmax
       l2=kgot(jloop)
       al2=float(l2)
       nrb=idra(iir,iib,l2)
       npg=idra(iip,iig,l2)
       ini=ntloc(nrb,npg,ipgsym,l2)
       call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
       fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

 17    continue
 16    continue
 12    continue

       if(orbg.ne.orbr)go to 2
       if(igsym.ne.irsym)go to 2

       inf=idpa(iir,iig)

       do 20 iq=1,nexcit
       iqsym=isyme(iq)
       iiq=ke(iq)
       iaq=iqe(iq)
       orbq=orbe(iq)
       kapq=nak(iiq)

       if(irsym.eq.iqsym)go to 20

       call nsi(iir,iiq,orbr,orbq,kapr,kapq,dipole)

       ipasym=mtbl(ipsym,iasym)
       iqgsym=mtbl(iqsym,igsym)
       iqasym=mtbl(iqsym,iasym)

       if(ipasym.ne.iqgsym)go to 20

c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------

       call findk(orba,orbp,orbq,orbg,iaa,iap,iaq,iag,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       npa=idra(iip,iia,0)
       nqg=idra(iiq,iig,0)
       ini=ntloc(npa,nqg,ipasym,0)
       fact=(-1)**(orbg+orbq+2*orbg+0)/(2*zero+1)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

       call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
       if(imax.eq.0)go to 20
       do 21 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nqa=idra(iiq,iia,l1)
       npg=idra(iip,iig,l1)
       ini=ntloc(nqa,npg,iqasym,l1)
       call sixj(zero,orbq,orbr,al1,orbp,orba,s6j)
       fact=(-1)**(orbq+orbr+0)*s6j
       sum=sum-dipole*ti(ini)*tf(inf)*fact

 21    continue
 20    continue


       do 22 ib=1,nocc
       ibsym=isymc(ib)
       iib=kc(ib)
       iab=iqc(ib)
       orbb=orbc(ib)
       kapb=nak(iib)

       if(ibsym.eq.igsym)go to 22

       call nsi(iib,iig,orbb,orbg,kapb,kapg,dipole)

       ipasym=mtbl(ipsym,iasym)
       irbsym=mtbl(irsym,ibsym)
       irasym=mtbl(irsym,iasym)

       if(ipasym.ne.irbsym)go to 22

c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------

       call findk(orba,orbp,orbr,orbb,iaa,iap,iar,iab,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       npa=idra(iip,iia,0)
       nrb=idra(iir,iib,0)
       ini=ntloc(npa,nrb,ipasym,0)
       fact=(-1)**(orbb+orbr+2*orbb+0)/(2*zero+1)
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

       call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
       if(imax.eq.0)go to 22
       do 23 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nra=idra(iir,iia,l1)
       npb=idra(iip,iib,l1)
       ini=ntloc(nra,npb,irasym,l1)
       call sixj(zero,orbr,orbb,al1,orbp,orba,s6j)
       fact=(-1)**(orbb+orbr+0)*s6j
       sum=sum+dipole*tf(inf)*ti(ini)*fact

 23    continue
 22    continue
  2    continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------

       do 24 iq=1,nexcit
       iqsym=isyme(iq)
       orbq=orbe(iq)
       iaq=iqe(iq)
       iiq=ke(iq)
       kapq=nak(iiq)

       if(iqsym.eq.ipsym)go to 24
       if(iqsym.ne.iasym)go to 24
       if(orbq.ne.orba)go to 24

       ini=idpa(iiq,iia)

       call nsi(iip,iiq,orbp,orbq,kapp,kapq,dipole)

       sum=sum+ti(ini)*dipole

 24    continue

c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------

       do 25 ib=1,nocc
       ibsym=isymc(ib)
       orbb=orbc(ib)
       iab=iqc(ib)
       iib=kc(ib)
       kapb=nak(iib)

       if(ibsym.eq.iasym)go to 25
       if(ibsym.ne.ipsym)go to 25
       if(orbp.ne.orbb)go to 25

       ini=idpa(iip,iib)

       call nsi(iib,iia,orbb,orba,kapb,kapa,dipole)

       sum=sum-ti(ini)*dipole

 25    continue

       do 26 iq=1,nexcit
       iqsym=isyme(iq)
       orbq=orbe(iq)
       iaq=iqe(iq)
       iiq=ke(iq)
       kapq=nak(iiq)

       do 26 ir=1,nexcit
       irsym=isyme(ir)
       orbr=orbe(ir)
       iar=iqe(ir)
       iir=ke(ir)
       kapr=nak(iir)

       do 26 ib=1,nocc
       ibsym=isymc(ib)
       orbb=orbc(ib)
       iab=iqc(ib)
       iib=kc(ib)
       kapb=nak(iib)

       do 26 ig=1,nocc
       igsym=isymc(ig)
       orbg=orbc(ig)
       iag=iqc(ig)
       iig=kc(ig)
       kapg=nak(iig)

       if(igsym.eq.irsym)go to 26

       call nsi(iir,iig,orbr,orbg,kapr,kapg,dipole)

       ipasym=mtbl(ipsym,iasym)
       iqbsym=mtbl(iqsym,ibsym)
       irgsym=mtbl(irsym,igsym)

       if(ipasym.ne.iqbsym)go to 27
       if(irgsym.ne.iqbsym)go to 27

c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------

       call findk(orba,orbp,orbq,orbb,iaa,iap,iaq,iab,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       npa=idra(iip,iia,0)
       nqb=idra(iiq,iib,0)
       nrg=idra(iir,iig,0)
       ini=ntloc(npa,nqb,ipasym,0)
       inf=ntloc(nrg,nqb,iqbsym,0)
       fact=(-1)**(orbg+orbb+orbr+orbq)/((2.0*0+1)*(2.0*0+1))
       sum=sum+ti(ini)*tf(inf)*dipole*fact

       endif
       enddo
       endif

 27    continue

c----------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg>
c------------------------------------------------------------------

       ipbsym=mtbl(ipsym,ibsym)
       iqasym=mtbl(iqsym,iasym)

       if(irgsym.ne.iqbsym)go to 28
       if(ipbsym.ne.iqasym)go to 28

       call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
       if(jmax.ne.0)then
       do jloop=1,jmax
       l2=jgot(jloop)
       if(l2.eq.0)then
       nqb=idra(iiq,iib,0)
       nrg=idra(iir,iig,0)
       inf=ntloc(nrg,nqb,iqbsym,0)

       call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
       if(imax.ne.0)then
       do 29 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npb=idra(iip,iib,l1)
       nqa=idra(iiq,iia,l1)
       ini=ntloc(npb,nqa,iqasym,l1)
       call sixj(zero,orbq,orbb,al1,orbp,orba,s6j)
       fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2.0*0+1)
       sum=sum-tf(inf)*ti(ini)*dipole*fact

 29    continue

       endif
       endif
       enddo
       endif

 28    continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

       ipasym=mtbl(ipsym,iasym)
       iqbsym=mtbl(iqsym,ibsym)
       irbsym=mtbl(irsym,ibsym)
       iqgsym=mtbl(iqsym,igsym)

       if(iqgsym.ne.irbsym)go to 30
       if(ipasym.ne.iqbsym)go to 30

       call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
       if(jmax.eq.0)go to 30
       do jloop=1,jmax
       l2=jgot(jloop)
       if (l2.eq.0)then
       npa=idra(iip,iia,0)
       nqb=idra(iiq,iib,0)
       ini=ntloc(npa,nqb,iqbsym,0)

       call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
       if(imax.ne.0)then
       do 31 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nqg=idra(iiq,iig,l1)
       nrb=idra(iir,iib,l1)
       inf=ntloc(nqg,nrb,iqgsym,l1)
       call sixj(zero,orbb,orbq,al1,orbg,orbr,s6j)
c      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2.0*1+1)
       fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2.0*0+1)
       sum=sum-ti(ini)*tf(inf)*dipole*fact

 31    continue

       endif
       endif
       enddo

 30    continue

c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

       ipbsym=mtbl(ipsym,ibsym)
       iqasym=mtbl(iqsym,iasym)
       irbsym=mtbl(irsym,ibsym)
       iqgsym=mtbl(iqsym,igsym)

       if(iqgsym.ne.irbsym)go to 32
       if(ipbsym.ne.iqasym)go to 32
 
       call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
       if(kmax.eq.0)go to 32
       do 33 kloop=1,kmax
       l2=kgot(kloop)
       al2=float(l2)
       npb=idra(iip,iib,l2)
       nqa=idra(iiq,iia,l2)
       ini=ntloc(npb,nqa,ipbsym,l2)

       call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
       if(imax.eq.0)go to 33
       do 34 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nqg=idra(iiq,iig,l1)
       nrb=idra(iir,iib,l1)
       inf=ntloc(nqg,nrb,iqgsym,l1)
       call sixj(zero,orbr,orbg,al2,orbq,orbb,s6j1)
       call sixj(zero,orba,orbp,al1,orbb,orbq,s6j2)
       fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
       sum=sum+ti(ini)*tf(inf)*dipole*fact

 34    continue
 33    continue
 32    continue
 26    continue

       d(iip,iia)=d(iip,iia)+sum
       d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+1)

  1    continue

       return
       end

c*****************************************************************
c     This subroutine calculated the hyperfine O bar diagram of  *
c     of the type  PP                                            *
c*************************************************************** *

       subroutine nsipp(d,ti,tf)

       implicit real*8 (a-h,o-z)
       common/jvalue/orbc(MNOCC),orbe(MNEXC)
       common/index/ke(MNEXC),kc(MNOCC)
       common/info/nocc,nexcit
       common/parity/isymc(MNOCC),isyme(MNEXC)
       common/symmetry/mtbl(MNS,MNS)
       common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
       common/orbital_energy/eorb(MNBAS)
       common/idra1/idra(MNBAS,MNBAS,0:MXV)
       common/idpa1/idpa(MNBAS,MNBAS)
       common/kappa/nak(MNBAS)
       common/inparm/nparm
       common/pncp/amass,qw


       dimension kgot(MXV),igot(MXV),jgot(MXV)
       dimension d(MNBAS,MNBAS)
       dimension ti(MDIM),tf(MDIM)

       data zero,half,one,two/0.0,0.5,1.0,2.0/

       do i=1,MXV
       kgot(i)=0
       igot(i)=0
       jgot(i)=0
       enddo


       do 1 ip=1,nexcit
       ipsym=isyme(ip)
       orbp=orbe(ip)
       iap=iqe(ip)
       iip=ke(ip)
       kapp=nak(iip)

       do 1 iq=1,ip
       iqsym=isyme(iq)
       orbq=orbe(iq)
       iaq=iqe(iq)
       iiq=ke(iq)
       kapq=nak(iiq)

       if(ipsym.eq.iqsym) go to 1

       call nsi(iip,iiq,orbp,orbq,kapp,kapq,dipole)

c--------------------------------------------------------------
c     diagram 1 = <p|d|q>
c--------------------------------------------------------------

       d(iip,iiq)=d(iip,iiq)+dipole

       sum=0.0d0

       do 2 ia=1,nocc
       iasym=isymc(ia)
       orba=orbc(ia)
       iia=kc(ia)
       iaa=iqc(ia)
       kapa=nak(iia)

       if(iasym.ne.iqsym)then
       if(iasym.eq.ipsym)then
       if(orbp.eq.orba)then

       call nsi(iia,iiq,orba,orbq,kapa,kapq,dipole)

       ini=idpa(iip,iia)

c----------------------------------------------------------------
c     diagram 2 = <q|d|a><a|t1|p>
c----------------------------------------------------------------

       sum=sum-dipole*ti(ini)

       endif
       endif
       endif

       if(iasym.ne.ipsym)then
       if(iasym.eq.iqsym)then
       if(orbq.eq.orba)then

       call nsi(iip,iia,orbp,orba,kapp,kapa,dipole)

       inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

       sum=sum-dipole*tf(inf)

       endif
       endif
       endif

       do 3 ig=1,nocc
       igsym=isymc(ig)
       orbg=orbc(ig)
       iag=iqc(ig)
       iig=kc(ig)
       kapg=nak(iig)

       do 3 ir=1,nexcit
       irsym=isyme(ir)
       orbr=orbe(ir)
       iar=iqe(ir)
       iir=ke(ir)
       kapr=nak(iir)

       irgsym=mtbl(irsym,igsym)
       iqasym=mtbl(iqsym,iasym)

       if(irgsym.eq.iqasym)then
       if(igsym.eq.irsym)then
       if(orbg.eq.orbr)then
       if(iasym.ne.ipsym)then

       ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

       call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       nrg=idra(iir,iig,0)
       nqa=idra(iiq,iia,0)
       inf=ntloc(nrg,nqa,irgsym,0)
       call nsi(iip,iia,orbp,orba,kapp,kapa,dipole)
c      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
       sum=sum-dipole*tf(inf)*ti(ini)*
     : dsqrt((2*orbg+1)/(two*orba+1))

       endif
       enddo
       endif

c---------------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c---------------------------------------------------------------------

       call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
       if(imax.ne.0)then
       if(orbq.eq.orba) then
       do iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nqg=idra(iiq,iig,l1)
       nra=idra(iir,iia,l1)
       irasym=mtbl(irsym,iasym)
       inf=ntloc(nqg,nra,irasym,l1)
       fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
       sum=sum+dipole*ti(ini)*tf(inf)*fact

       enddo
       endif
       endif
       endif
       endif
       endif
       endif

       if(orba.eq.orbp)then
       if(iasym.eq.ipsym)then
       if(igsym.ne.irsym)then

       iqasym=mtbl(iqsym,iasym)

       if(irgsym.eq.iqasym)then

       call nsi(iir,iig,orbr,orbg,kapr,kapg,dipole)

       ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

       call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       nrg=idra(iir,iig,0)
       nqa=idra(iiq,iia,0)
       inf=ntloc(nrg,nqa,irgsym,0)
       fact=(-1)**(orbg+orbr+2*orbr+0)/(2*zero+1)
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif

c-------------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c-------------------------------------------------------------------

       call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nqg=idra(iiq,iig,l1)
       nra=idra(iir,iia,l1)
       irasym=mtbl(irsym,iasym)
       inf=ntloc(nqg,nra,irasym,l1)
       call sixj(zero,orbg,orbr,al1,orbp,orbq,s6j)
       fact=(-1)**(orbg+orbr+zero)*s6j
       sum=sum+dipole*ti(ini)*tf(inf)*fact

       enddo
       endif
       endif
       endif
       endif
       endif


       if(orba.eq.orbq)then
       if(iasym.eq.iqsym)then
       if(igsym.ne.irsym)then

       ipasym=mtbl(ipsym,iasym)

       if(irgsym.eq.ipasym)then

       call nsi(iig,iir,orbg,orbr,kapg,kapr,dipole)

       inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

       call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       nrg=idra(iir,iig,0)
       npa=idra(iip,iia,0)
       ini=ntloc(nrg,npa,irgsym,0)
       fact=(-1)**(orbg+orbr+2*orbr+0)/(2*zero+1)
       sum=sum-dipole*ti(ini)*tf(inf)*fact

       endif
       enddo
       endif

c------------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c------------------------------------------------------------------

       call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npg=idra(iip,iig,l1)
       nra=idra(iir,iia,l1)
       irasym=mtbl(irsym,iasym)
       ini=ntloc(npg,nra,irasym,l1)
       call sixj(zero,orbg,orbr,al1,orbq,orbp,s6j)
       fact=(-1)**(orbg+orbr+zero)*s6j
       sum=sum+dipole*tf(inf)*ti(ini)*fact

       enddo
       endif
       endif
       endif
       endif
       endif

       if(orbr.eq.orbg)then
       if(irsym.eq.igsym)then
       if(iasym.ne.iqsym)then

       ipasym=mtbl(ipsym,iasym)

       if(irgsym.eq.ipasym)then

       call nsi(iia,iiq,orba,orbq,kapa,kapq,dipole)

       inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

       call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
       if(imax.ne.0)then
       do iloop=1,imax
       l1=igot(iloop)
       if(l1.eq.0)then
       nrg=idra(iir,iig,0)
       npa=idra(iip,iia,0)
       ini=ntloc(nrg,npa,irgsym,0)
c      sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
       sum=sum-dipole*ti(ini)*tf(inf)*
     : dsqrt((2*orbg+1)/(two*orba+1))

       endif
       enddo
       endif

c--------------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c--------------------------------------------------------------------

       call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
       if(imax.ne.0)then
       if(orba.eq.orbp) then
       do iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       npg=idra(iip,iig,l1)
       nra=idra(iir,iia,l1)
       irasym=mtbl(irsym,iasym)
       ini=ntloc(npg,nra,irasym,l1)
       fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

       enddo
       endif
       endif
       endif
       endif
       endif
       endif

  3    continue
  2    continue

       do 4 ig=1,nocc
       igsym=isymc(ig)
       orbg=orbc(ig)
       iag=iqc(ig)
       iig=kc(ig)
       kapg=nak(iig)

       do 5 ir=1,nexcit
       irsym=isyme(ir)
       orbr=orbe(ir)
       iar=iqe(ir)
       iir=ke(ir)
       kapr=nak(iir)

c------------------------------------------------------------------
c     diagram 13 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

       if(irsym.ne.iqsym)then
       if(igsym.eq.ipsym)then
       if(orbg.eq.orbp)then
       if(igsym.eq.irsym)then
       if(orbg.eq.orbr)then

       call nsi(iir,iiq,orbr,orbq,kapr,kapq,dipole)

       ini=idpa(iip,iig)
       inf=idpa(iir,iig)

       sum=sum-dipole*ti(ini)*tf(inf)

       endif
       endif
       endif
       endif
       endif

c---------------------------------------------------------------
c     diagram 14 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------

       if(irsym.ne.ipsym)then
       if(igsym.eq.irsym)then
       if(orbg.eq.orbr)then
       if(igsym.eq.iqsym)then
       if(orbg.eq.orbq)then

       call nsi(iip,iir,orbp,orbr,kapp,kapr,dipole)

       ini=idpa(iir,iig)
       inf=idpa(iiq,iig)

       sum=sum-dipole*ti(ini)*tf(inf)

       endif
       endif
       endif
       endif
       endif

  5    continue

c------------------------------------------------------------------
c     diagram 12 = <q|t1|a><a|d|g><g|t1|p>
c------------------------------------------------------------------

       do 6 ia=1,nocc
       iasym=isymc(ia)
       orba=orbc(ia)
       iaa=iqc(ia)
       iia=kc(ia)
       kapa=nak(iia)

       if(igsym.ne.iasym)then
       if(igsym.eq.ipsym)then
       if(orbg.eq.orbp)then
       if(iasym.eq.iqsym)then
       if(orba.eq.orbq)then

       call nsi(iig,iia,orbg,orba,kapg,kapa,dipole)

       ini=idpa(iip,iig)
       inf=idpa(iiq,iia)

       sum=sum+dipole*ti(ini)*tf(inf)

       endif
       endif
       endif
       endif
       endif

  6    continue
  4    continue

       do 7 ig=1,nocc
       igsym=isymc(ig)
       orbg=orbc(ig)
       iag=iqc(ig)
       iig=kc(ig)
       kapg=nak(iig)

       do 7 ir=1,nexcit
       irsym=isyme(ir)
       orbr=orbe(ir)
       iar=iqe(ir)
       iir=ke(ir)
       kapr=nak(iir)

       do 7 ia=1,nocc
       iasym=isymc(ia)
       orba=orbc(ia)
       iaa=iqc(ia)
       iia=kc(ia)
       kapa=nak(iia)

       do 8 id=1,nocc
       idsym=isymc(id)
       orbd=orbc(id)
       iad=iqc(id)
       iid=kc(id)
       kapd=nak(iid)

       if(igsym.eq.idsym)go to 9

       call nsi(iid,iig,orbd,orbg,kapd,kapg,dipole)

       irgsym=mtbl(irsym,igsym)
       iqasym=mtbl(iqsym,iasym)
       irdsym=mtbl(irsym,idsym)
       ipasym=mtbl(ipsym,iasym)

       if(ipasym.ne.irdsym)go to 9
       if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

       call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
       if(imax.eq.0)go to 9
       do 10 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nrg=idra(iir,iig,l1)
       nqa=idra(iiq,iia,l1)
       inf=ntloc(nrg,nqa,irgsym,l1)

       call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
       if(kmax.eq.0)go to 10
       do 11 kloop=1,kmax
       l2=kgot(kloop)
       al2=float(l2)
       nrd=idra(iir,iid,l2)
       npa=idra(iip,iia,l2)
       ini=ntloc(nrd,npa,irdsym,l2)
       isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
       call sixj(orbg,orbd,zero,al2,al1,orbr,s6j1)
       call sixj(orbq,orbp,zero,al2,al1,orba,s6j2)
       fact=isign*s6j1*s6j2
       sum=sum+dipole*tf(inf)*ti(ini)*fact

 11    continue
 10    continue

c------------------------------------------------------------------
c     diagram 20 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------

       irgsym=mtbl(irsym,igsym)
       iqasym=mtbl(iqsym,iasym)
       ipdsym=mtbl(ipsym,idsym)
       irasym=mtbl(irsym,iasym)

       if(irgsym.ne.iqasym)go to 12
       if(ipdsym.ne.irasym)go to 12

       call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
       if(imax.eq.0)go to 12
       do 13 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nrg=idra(iir,iig,l1)
       nqa=idra(iiq,iia,l1)
       ini=ntloc(nrg,nqa,irgsym,l1)

       call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
       if(kmax.eq.0)go to 13
       do 14 kloop=1,kmax
       l2=kgot(kloop)
       al2=float(l2)
       nra=idra(iir,iia,l2)
       npd=idra(iip,iid,l2)
       inf=ntloc(nra,npd,irasym,l2)
       call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
       call sixj(orbp,al2,orbd,orbg,zero,orbq,s6j2)
       fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
       sum=sum-dipole*tf(inf)*ti(ini)*fact

 14    continue
 13    continue
 12    continue
  9    continue



       if(iasym.eq.idsym)go to 8

       call nsi(iid,iia,orbd,orba,kapd,kapa,dipole)

       irgsym=mtbl(irsym,igsym)
       iqasym=mtbl(iqsym,iasym)
       irdsym=mtbl(irsym,idsym)
       ipgsym=mtbl(ipsym,igsym)

       if(ipgsym.ne.irdsym)go to 8
       if(iqasym.ne.irgsym)go to 8

c-----------------------------------------------------------------
c     digram 19 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------

       call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
       if(imax.eq.0)go to 8
       do 15 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nrg=idra(iir,iig,l1)
       nqa=idra(iiq,iia,l1)
       inf=ntloc(nrg,nqa,irgsym,l1)

       call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
       if(jmax.ne.0)then
       do jloop=1,jmax
       k1=jgot(jloop)
       if(k1.eq.l1)then
       nrg=idra(iir,iig,l1)
       npd=idra(iip,iid,l1)
       ini=ntloc(nrg,npd,irgsym,l1)
       call sixj(zero,orba,orbd,al1,orbp,orbq,s6j)
       fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+zero+l1)
     &  /(2*al1+1)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif

 42    continue

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

       call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
       if(kmax.eq.0)go to 15
       do 16 kloop=1,kmax
       l2=kgot(kloop)
       al2=float(l2)
       nrd=idra(iir,iid,l2)
       npg=idra(iip,iig,l2)
       ini=ntloc(nrd,npg,irdsym,l2)
       isign=(-1)**(orba+orbd+zero+orbg+orbr+l1)
       call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
       call sixj(orbq,al1,orba,orbd,zero,orbp,s6j2)
       fact=isign*s6j1*s6j2
       sum=sum-dipole*tf(inf)*ti(ini)*fact

 16    continue
 15    continue
  8    continue


       do 17 is=1,nexcit
       issym=isyme(is)
       orbs=orbe(is)
       ias=iqe(is)
       iis=ke(is)
       kaps=nak(iis)

       if(issym.eq.irsym)go to 18

       call nsi(iis,iir,orbs,orbr,kaps,kapr,dipole)

       isgsym=mtbl(issym,igsym)
       ipasym=mtbl(ipsym,iasym)
       iqasym=mtbl(iqsym,iasym)
       irgsym=mtbl(irsym,igsym)

       if(irgsym.ne.ipasym)go to 18
       if(isgsym.ne.iqasym)go to 18

c------------------------------------------------------------------
c     diagram 16 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------

       call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
       if(imax.eq.0)go to 18
       do 19 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nsg=idra(iis,iig,l1)
       nqa=idra(iiq,iia,l1)
       inf=ntloc(nsg,nqa,isgsym,l1)

       call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
       if(kmax.eq.0)go to 19
       do 20 kloop=1,kmax
       l2=kgot(kloop)
       al2=float(l2)
       nrg=idra(iir,iig,l2)
       npa=idra(iip,iia,l2)
       ini=ntloc(nrg,npa,irgsym,l2)
       call sixj(orbs,zero,orbr,al2,orbg,al1,s6j1)
       call sixj(orbq,orbp,zero,al2,al1,orba,s6j2)
       fact=(-1)**(orbp+orbq+zero+2*orbq)*s6j1*s6j2
       sum=sum-dipole*tf(inf)*ti(ini)*fact

 20    continue

c------------------------------------------------------------------
c     diagram 17 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------

       call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
       if(kmax.eq.0)go to 19
       do 21 kloop=1,kmax
       l2=kgot(kloop)
       al2=float(l2)
       nra=idra(iir,iia,l2)
       npg=idra(iip,iig,l2)
       irasym=mtbl(irsym,iasym)
       ini=ntloc(nra,npg,irasym,l2)
       call ninej(al1,orba,orbq,orbs,orbr,zero,orbg,al2,orbp,s9j)
       fact=(-1)**(orbq+orbg+l2+zero+2*orbp)*s9j
       sum=sum+dipole*tf(inf)*ti(ini)*fact

 21    continue
 19    continue
 18    continue

c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------

       irgsym=mtbl(irsym,igsym)
       iassym=mtbl(iasym,issym)
       iqasym=mtbl(iqsym,iasym)

       if (irgsym .ne. iassym) goto 22
       if (irgsym .ne. iqasym) goto 22
       if(ipsym.eq.issym)go to 22

       call nsi(iip,iis,orbp,orbs,kapp,kaps,dipole)

       if(orbs.ne.orbq)go to 22

       call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
       if(imax.eq.0)go to 22
       do 23 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nrg=idra(iir,iig,l1)
       nqa=idra(iiq,iia,l1)
       irgsym=mtbl(irsym,igsym)
       inf=ntloc(nrg,nqa,irgsym,l1)

       call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,jgot,jmax)
       if(jmax.ne.0)then
       do jloop=1,jmax
       k1=jgot(jloop)
       if(k1.eq.l1)then
       nrg=idra(iir,iig,l1)
       nsa=idra(iis,iia,l1)
       ini=ntloc(nrg,nsa,irgsym,l1)
       fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

       call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
       if(kmax.eq.0)go to 23
       do 24 kloop=1,kmax
       l2=kgot(kloop)
       al2=float(l2)
       nsg=idra(iis,iig,l2)
       nra=idra(iir,iia,l2)
       irasym=mtbl(irsym,iasym)
       ini=ntloc(nsg,nra,irasym,l2)
       call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
       fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

 24    continue
 23    continue
 22    continue

c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

       irgsym=mtbl(irsym,igsym)
       iapsym=mtbl(iasym,ipsym)
       iassym=mtbl(iasym,issym)

       if(irgsym .ne. iassym) goto 25
       if(irgsym .ne. iapsym) goto 25
       if(iqsym.eq.issym)go to 25

       call nsi(iis,iiq,orbs,orbq,kaps,kapq,dipole)

       if(orbs.ne.orbp)go to 25

       call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
       if(imax.eq.0)go to 25
       do 26 iloop=1,imax
       l1=igot(iloop)
       al1=float(l1)
       nrg=idra(iir,iig,l1)
       nsa=idra(iis,iia,l1)
       irgsym=mtbl(irsym,igsym)
       inf=ntloc(nrg,nsa,irgsym,l1)

       call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
       if(jmax.ne.0)then
       do jloop=1,jmax
       k1=jgot(jloop)
       if(k1.eq.l1)then
       npa=idra(iip,iia,l1)
       nrg=idra(iir,iig,l1)
       ini=ntloc(nrg,npa,irgsym,l1)
       fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
       sum=sum-dipole*tf(inf)*ti(ini)*fact

       endif
       enddo
       endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

       call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
       if(kmax.eq.0)go to 26
       do 27 kloop=1,kmax
       l2=kgot(kloop)
       al2=float(l2)
       npg=idra(iip,iig,l2)
       nra=idra(iir,iia,l2)
       irasym=mtbl(irsym,iasym)
       ini=ntloc(npg,nra,irasym,l2)
       call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
       fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
       sum=sum+dipole*tf(inf)*ti(ini)*fact

  27   continue
  26   continue
  25   continue
  17   continue
   7   continue

       d(iip,iiq)=d(iip,iiq)+sum
       d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+1)

   1   continue

       return
       end


c********************************************************************



c*******************************************************************
c                                                                  *
c  This subroutine calculates the radial part of the reduced       *
c  metrix elements of the NSI PNC operator                         *
c                                                                  *
c*******************************************************************


c********************************************************************


        subroutine Calcnsira(nbasis)


        implicit real*8 (a-h,o-z)

        common/cons/zero,half,tenth,one,two,three,ten
     :  /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /tatb/ta(MN),tb(MN),mtp
     :  /wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :  mf(MNBAS)
     :  /disn/aa(MN)

        common/intnsimat/rintnsi(MNBAS,MNBAS)
        common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)
        common/kappa/nak(MNBAS)
        common/inparm/nparm
        common/pncp/amass,qw


        do i1 = 1, nbasis
        do i2 = 1, i1

        rintnsi(i1,i2) = 0.0d0

        if(iparity(i1).ne.iparity(i2))then
        if(nak(i1).ne.-nak(i2))then

        ta(1)=zero

        do 1 l=2,mtp
        call nuc
        ta(l)=((qf(l,i1)*pf(l,i2)-pf(l,i1)*qf(l,i2))*aa(l)*rp(l))
    1   continue
        call quad(result)

        endif
        endif

        rintnsi(i1,i2) = result
        rintnsi(i2,i1) = -rintnsi(i1,i2)

        enddo
        enddo

        return
        end

c*******************************************************************

c*********************************************************************
c      The reduced matrix element of the NSI PNC matrix element is   *
c      calculated                                                    *
c*********************************************************************

        subroutine nsi(i1,i2,jorb1,jorb2,kap1,kap2,ansim)

        implicit double precision (a-h,o-z)
        double precision j1,j2,jorb1,jorb2

        common/cons/zero,half,tenth,one,two,three,ten
        common/orbital_energy/eorb(MNBAS)
        common/intnsimat/rintnsi(MNBAS,MNBAS)
        common/inparm/nparm
        common/pncp/amass,qw


        ansim= zero      

        f=2*1.41421356237
        g=2.2225e-14
        factor=dble(qw*g/f)

        ans1=rintnsi(i1,i2)

        ansim=factor*ans1

        return
        end

c***********************************************************************


c***********************************************************************
c       This subroutine calculates the nuclear distribution function   *
c       as a point, uniform sphere or fermi distribution function      *
c***********************************************************************

              
       subroutine nuc

       implicit doubleprecision (a-h, o-z) 

       common /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :        /cons/zero,half,tenth,one,two,three,ten
     :        /disn/aa(MN)
       common/inparm/nparm
       common/pncp/amass,qw

c---------------------------------------------------------------------
c    Set the parameters for nuclear distribution function
c---------------------------------------------------------------------

       pi = 3.141592653589790d 000
       tfm=2.3d0
       ainfcm = 0.529177249d-08
       fmtoau = 1.0d-13/ainfcm
       t=tfm*fmtoau
       const= 4.0d0*log(3.0d0)
       parma=t/const
       rrmsfm=0.836d0*amass**(1.d0/3.d0)+0.57d0
       rrms = rrmsfm*fmtoau
       factor = rrms**2-(7.0d 00/5.0d 00)
     :                    *(pi**2)*(parma**2)
       parmc = sqrt(5.0d 00/3.0d 00)*dsqrt(factor)

c---------------------------------------------------------------------
c                  For point nucleus
c---------------------------------------------------------------------

       if (nparm.eq.o) then

       do 1 i=1,MN
           aa(i)=1
    1  continue

       endif

c---------------------------------------------------------------------
c                  For uniform spherical distribution
c---------------------------------------------------------------------

       if (nparm.eq.1) then

       rn2=(5.0d 00/3.0d 00)*(rrms**2)
       rn=dsqrt(rn2)

       do 2 i=1,MN
       fac=4*pi*rn**3
       aa(i)= dble(3/fac)
    2  continue

       endif

c---------------------------------------------------------------------
c                 For Fermi type distribution
c                 Taking the analytical expansion of Mohanty 
c                 and Perpia PRA, 46, 3735, (1992)
c--------------------------------------------------------------------

       if (nparm.eq.2) then

       c=parmc
       a=parma

       do 3 i=1, MN

       ri=r(i)
       rmc=ri-c
       rmcba=rmc/a
       rbc=ri/c

       if (rbc.le.one) then

       call es1(rmcba, s2rcba)
       aa(i)= 1+s2rcba
      
       else

       call es2(-rmcba, s3rcba)
       aa(i)=s3rcba

       endif

    3  continue

       endif

c---------------------------------------------------------------------
c                 For Fermi type distribution
c                 Without any analytical expansion
c--------------------------------------------------------------------

       if (nparm.eq.3) then

       c=parmc
       a=parma

       do 4 i=1,MN
          ri=r(i)
          rmc=ri-c
          rmcba=rmc/a
          enf=1+dexp(rmcba)
          aa(i)=dble(1/enf)
    4  continue

       endif

       return
       end

c**********************************************************************

c**********************************************************************
c                        evaluate the sum of the series               *
c                                                                     *
c          k(f)=+sum_(1)^(infinity)(-1)^n  exp(n*f)                   *
c                                                                     *
c       for k = 2, 3 to machine precision. this is a utility          *
c**********************************************************************

       subroutine es1 (f,s2f)

       implicit doubleprecision (a-h, o-z)
       common /cons/zero,half,tenth,one,two,three,ten

       n=0
       s2f= zero
       fase=one
    1  n=n+1
        en= dble(n)
        fase=-fase
        enf=dexp(en*f)
        term2=fase*enf
        s2last=s2f
        s2f=s2f+term2
        if (dabs (s2f) .ne. dabs (s2last)) goto 1

        return
        end

c**********************************************************************

c**********************************************************************
c                        evaluate the sum of the series               *
c                                                                     *
c          k(f)=-sum_(1)^(infinity)(-1)^n  exp(n*f)                   *
c                                                                     *
c       for k = 2, 3 to machine precision. this is a utility          *
c**********************************************************************


       subroutine es2 (f,s2f)

       implicit doubleprecision (a-h, o-z)
       common /cons/zero,half,tenth,one,two,three,ten

       n=0
       s2f= zero
       fase=-one
    1  n=n+1
        en= dble(n)
        fase=-fase
        enf=dexp(en*f)
        term2=fase*enf
        s2last=s2f
        s2f=s2f+term2
        if (dabs (s2f) .ne. dabs (s2last)) goto 1

        return
        end

c**********************************************************************


c***********************************************************************
c     This subroutine calculates the NSD PNC operator diagram of       *
c     of the type HH                                                   *
c***********************************************************************

      subroutine nsdhh(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)
      common/inparm/nparm
      common/pncp/amass,qw


      dimension kgot(MXV),igot(MXV),jgot(MXV)
      dimension ti(MDIM),d(MNBAS,MNBAS),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=1,MXV
      kgot(i)=0
      jgot(i)=0
      igot(i)=0
      enddo

      nbasis=nocc+nexcit

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)

      if(ibsym.eq.iasym)go to 1
 
c-----------------------------------------------------------------
c     diagram 1  = <b|d|a>
c-----------------------------------------------------------------

      call nsdmat(iib,iia,orbb,orba,kapb,kapa,dipole)
      d(iib,iia)=d(iib,iia)+dipole

      sum=zero

      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

c------------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c------------------------------------------------------------------

      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.ne.ibsym)then

      ini=idpa(iip,iia)

      call nsdmat(iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)

      endif
      endif
      endif

c------------------------------------------------------------------
c     diagram 3 = <b|t1|p><p|d|a>
c------------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.ne.iasym)then

      inf=idpa(iip,iib)
      call nsdmat(iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)

      endif
      endif
      endif

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.eq.ipbsym) then

c---------------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c---------------------------------------------------------------------

      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.ne.iasym)then

      ini=idpa(iir,iig)

      call findk(orbb,orbp,orbg,orbr,iab,iap,iag,iar,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
      call nsdmat(iip,iia,orbp,orba,kapp,kapa,dipole)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 5 = <bg|t2|rp><p|d|a><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      if(orbb.eq.orbp)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif
      endif

      if(ipbsym.eq.irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.ne.igsym)then

      call nsdmat(iir,iig,orbr,orbg,kapr,kapg,dipole)

c-------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c-------------------------------------------------------------------

      ini=idpa(iip,iia)

      call findk(orbb,orbp,orbg,orbr,iab,iap,iag,iar,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npb=idra(iip,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(npb,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c---------------------------------------------------------------------

      ini=idpa(iip,iia)

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c----------------------------------------------------------------
c     diagrams 8 = <g|d|r><b|t1|p><pr|t2|ag>
c----------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.ne.igsym)then

      call nsdmat(iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      call sixj(one,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      inf=idpa(iip,iib)
      ini=ntloc(npg,nra,ipgsym,l1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

      if(ipasym.eq.irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.ne.ibsym)then

      call nsdmat(iib,iip,orbb,orbp,kapb,kapp,dipole)

      inf=idpa(iir,iig)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      if(orba.eq.orbp) then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      if(ipsym.eq.irsym)go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c------------------------------------------------------------------
c     diagram 12 = <b|t1|r><r|d|p><p|t1|a>
c------------------------------------------------------------------

      call nsdmat(iir,iip,orbr,orbp,kapr,kapp,dipole)
      inf=idpa(iir,iib)
      ini=idpa(iip,iia)
      sum=sum+dipole*tf(inf)*ti(ini)

 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

c------------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c------------------------------------------------------------------

      if(igsym.ne.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call nsdmat(iib,iig,orbb,orbg,kapb,kapg,dipole)

      ini=idpa(iip,iia)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c---------------------------------------------------------------------

      if(igsym.ne.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call nsdmat(iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iip,iib)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nak(iid)

      if(igsym.eq.idsym)go to 8

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

      call nsdmat(iid,iig,orbd,orbg,kapd,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+1+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 10   continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,one,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11   continue
  9   continue
  8   continue

      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      if(iqsym.eq.irsym)go to 13

      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13

      call nsdmat(iiq,iir,orbq,orbr,kapq,kapr,dipole)

c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,one,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue

c----------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c----------------------------------------------------------------------

      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(one,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+one+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo

 14   continue
 13   continue

      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.eq.irsym)go to 16

      call nsdmat(iir,iip,orbr,orbp,kapr,kapq,dipole)

c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      call sixj(one,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+one)
      fact=s6j*isign/(2*al1+1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,one,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+one)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 18   continue
 17   continue
 16   continue
 12   continue

      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nak(iid)

      if(idsym.eq.ibsym)go to 21

      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22

      call nsdmat(iib,iid,orbb,orbd,kapb,kapd,dipole)

c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 22=   <gd|t2|pr><b|d|d><pr|t2|ag>
c---------------------------------------------------------------------

      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 24   continue
 23   continue
 22   continue
 21   continue


      if(idsym.eq.iasym)go to 25

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26

      call nsdmat(iid,iia,orbd,orba,kapd,kapa,dipole)

c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      ini=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l2)
      ini=ntloc(npd,nrg,irgsym,l1)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+1)

 1    continue

      return
      end


c*****************************************************************
c     This subroutine calculates NSD PNC operator diagram of     *
c     of the type  HP                                            *
c*************************************************************** *


      subroutine nsdhp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)
      common/inparm/nparm
      common/pncp/amass,qw


      dimension kgot(MXV),igot(MXV),jgot(MXV)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=1,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo


      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nak(iip)

      if(ipsym.eq.iasym)go to 1

c--------------------------------------------------------------
c     diagram 1  = <p|d|a>
c--------------------------------------------------------------

      call nsdmat(iip,iia,orbp,orba,kapp,kapa,dipole)

      d(iip,iia)=d(iip,iia)+dipole

      sum=zero


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      if(igsym.eq.irsym)go to 3

      call nsdmat(iig,iir,orbg,orbr,kapg,kapr,dipole)

      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.irgsym)go to 401

c--------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c--------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,ipasym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(2*one+1)
      sum=sum+dipole*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c----------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.eq.0)go to 401
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dipole*ti(ini)*fact

 4    continue
 401  continue


c---------------------------------------------------------------------
c     diagram 8 = <r|t1|a><g|d|r><p|t1|g>
c---------------------------------------------------------------------

      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3
      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
      call nsdmat(iig,iir,orbg,orbr,kapg,kapr,dipole)
      sum=sum-dipole*ti(ini)*ti(inj)

 3    continue

c----------------------------------------------------------------------
c     diagram 9 = <r|t1|a><g|t1|r><p|d|g>
c----------------------------------------------------------------------

      if(ipsym.eq.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5

      call nsdmat(iip,iig,orbp,orbg,kapp,kapg,dipole)

      ini=idpa(iir,iia)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 5    continue

c---------------------------------------------------------------------
c     diagram 10 = <r|d|a><g|t|r><p|t|g>
c---------------------------------------------------------------------

      if(irsym.eq.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6

      ini=idpa(iir,iig)
      inf=idpa(iip,iig)

      call nsdmat(iir,iia,orbr,orba,kapr,kapa,dipole)

      sum=sum-dipole*ti(ini)*tf(inf)

 6    continue

      if(igsym.ne.irsym)go to 7
      if(orbg.ne.orbr)go to 7

      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      if(ipsym.eq.iqsym)go to 8

      call nsdmat(iip,iiq,orbp,orbq,kapp,kapq,dipole)

      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.iqasym)go to 8
      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 4 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(2*orbq+1))

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 5 = <g|t1|r><p|d|q><qr|t2|ga>
c--------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 9    continue
 8    continue

      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)

      if(iasym.eq.ibsym)go to 10

      call nsdmat(iib,iia,orbb,orba,kapb,kapa,dipole)

      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)

      if(irgsym.ne.ipbsym)go to 10

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 6 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbp.eq.orbb)then
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))

      endif
      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagram 7 = <g|t1|r><b|d|a><pr|t2|gb>
c-------------------------------------------------------------------

      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      if(orbp.ne.orbb)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue
  7   continue


      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)

      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipsym.eq.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      if(irgsym.ne.iqasym) goto 13

      call nsdmat(iip,iib,orbp,orbb,kapp,kapb,dipole)

c---------------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c---------------------------------------------------------------------

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue
 14   continue
 13   continue

      if(orbp.ne.orbq)go to 12

      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.eq.iasym)go to 12

c-------------------------------------------------------------------
c     diagram 13 = <gb|t2|rq><pr|t2|bg><q|d|a>
c-------------------------------------------------------------------

      call nsdmat(iiq,iia,orbq,orba,kapq,kapa,dipole)

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      npb=idra(iip,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif


c-------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c-------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 17   continue
 16   continue
 12   continue

      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2

      inf=idpa(iir,iig)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nak(iiq)

      if(irsym.eq.iqsym)go to 20

      call nsdmat(iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipasym.ne.iqgsym)go to 20

c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orba,orbp,orbq,orbg,iaa,iap,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(one,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact

 21   continue
 20   continue

      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)

      if(ibsym.eq.igsym)go to 22

      call nsdmat(iib,iig,orbb,orbg,kapb,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)

      if(ipasym.ne.irbsym)go to 22

c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbb,iaa,iap,iar,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrb=idra(iir,iib,1)
      ini=ntloc(npa,nrb,ipasym,1)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*1+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(one,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 23   continue
 22   continue
  2   continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------

      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)

      if(iqsym.eq.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24

      ini=idpa(iiq,iia)

      call nsdmat(iip,iiq,orbp,orbq,kapp,kapq,dipole)

      sum=sum+ti(ini)*dipole

 24   continue

c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------

      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)

      if(ibsym.eq.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25

      ini=idpa(iip,iib)

      call nsdmat(iib,iia,orbb,orba,kapb,kapa,dipole)

      sum=sum-ti(ini)*dipole

 25   continue

      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)

      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)

      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      if(igsym.eq.irsym)go to 26

      call nsdmat(iir,iig,orbr,orbg,kapr,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------

      call findk(orba,orbp,orbq,orbb,iaa,iap,iaq,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2.0*1+1)*(2.0*1+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

 27   continue

c----------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg>
c------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)

      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(1.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2.0*1+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 29   continue

      endif
      endif
      enddo
      endif

 28   continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 30
      do jloop=1,jmax
      l2=jgot(jloop)
      if (l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(1.0d0,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2.0*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2.0*1+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact

 31   continue

      endif
      endif
      enddo

 30   continue

c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32

      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(1.0d0,orbr,orbg,al2,orbq,orbb,s6j1)
      call sixj(1.0d0,orba,orbp,al1,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 34   continue
 33   continue
 32   continue
 26   continue

      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+1)

  1   continue

      return
      end



c*****************************************************************
c     This subroutine calculates NSD PNC operator diagram of     *
c     of the type  PP                                            *
c*************************************************************** *


      subroutine nsdpp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)
      common/inparm/nparm
      common/pncp/amass,qw


      dimension kgot(MXV),igot(MXV),jgot(MXV)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=1,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nak(iip)

      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)

      if(ipsym.eq.iqsym) go to 1

      call nsdmat(iip,iiq,orbp,orbq,kapp,kapq,dipole)

c--------------------------------------------------------------
c     diagram 1 = <p|d|q>
c--------------------------------------------------------------

      d(iip,iiq)=d(iip,iiq)+dipole

      sum=0.0d0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      if(iasym.ne.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then

      call nsdmat(iia,iiq,orba,orbq,kapa,kapq,dipole)

      ini=idpa(iip,iia)

c----------------------------------------------------------------
c     diagram 2 = <p|t|a><a|d|q>
c----------------------------------------------------------------

      sum=sum-dipole*ti(ini)

      endif
      endif
      endif

      if(iasym.ne.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then

      call nsdmat(iip,iia,orbp,orba,kapp,kapa,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

      sum=sum-dipole*tf(inf)

      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.ne.ipsym)then

      ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
      call nsdmat(iip,iia,orbp,orba,kapp,kapa,dipole)
c     sum=sum-dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      sum=sum-dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orba+1))

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c---------------------------------------------------------------------

      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      if(orbq.eq.orba) then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif
      endif


      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.ne.irsym)then

      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then

      call nsdmat(iir,iig,orbr,orbg,kapr,kapg,dipole)

      ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      nqa=idra(iiq,iia,1)
      inf=ntloc(nrg,nqa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.ne.irsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call nsdmat(iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      npa=idra(iip,iia,1)
      ini=ntloc(nrg,npa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.ne.iqsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call nsdmat(iia,iiq,orba,orbq,kapa,kapq,dipole)

      inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(two*orba+1))

      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      if(orba.eq.orbp) then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

c------------------------------------------------------------------
c     diagram 13 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

      if(irsym.ne.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then

      call nsdmat(iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------
c     diagram 14 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------

      if(irsym.ne.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then

      call nsdmat(iip,iir,orbp,orbr,kapp,kapr,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  5   continue

c------------------------------------------------------------------
c     diagram 12 = <q|t1|a><a|d|g><g|t1|p
c------------------------------------------------------------------

      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nak(iia)

      if(igsym.ne.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then

      call nsdmat(iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)

      sum=sum+dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  6   continue
  4   continue

      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nak(iia)

      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nak(iid)

      if(igsym.eq.idsym)go to 9

      call nsdmat(iid,iig,orbd,orbg,kapd,kapg,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue

c------------------------------------------------------------------
c     diagram 20 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)

      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      inf=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,one,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 14   continue
 13   continue
 12   continue
  9   continue


      if(iasym.eq.idsym)go to 8

      call nsdmat(iid,iia,orbd,orba,kapd,kapa,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)

      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8

c-----------------------------------------------------------------
c     digram 19 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(one,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+one+l1)
     & /(2*al1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

 42   continue

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+one+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,one,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 16   continue
 15   continue
  8   continue

c------------------------------------------------------------------
c     diagram 16 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------

      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nak(iis)

      if(issym.eq.irsym)go to 18

      call nsdmat(iis,iir,orbs,orbr,kaps,kapr,dipole)

      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18

      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,one,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 20   continue

c------------------------------------------------------------------
c     diagram 17 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,one,orbg,al2,orbp,s9j)
      fact=(-1)**(orbq+orbg+l2+one+2*orbp)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 21   continue
 19   continue
 18   continue

c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)

      if (irgsym .ne. iassym) goto 22
      if (irgsym .ne. iqasym) goto 22
      if(ipsym.eq.issym)go to 22

      call nsdmat(iip,iis,orbp,orbs,kapp,kaps,dipole)

      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 24   continue
 23   continue
 22   continue

c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)

      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.eq.issym)go to 25

      call nsdmat(iis,iiq,orbs,orbq,kaps,kapq,dipole)

      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 27   continue
 26   continue
 25   continue
 17   continue
  7   continue

      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+1)

  1   continue

      return
      end

c*********************************************************************


c*******************************************************************
c                                                                  *
c  This subroutine calculates the radial part of the reduced       *
c  metrix elements of the NSD PNC operator                         *
c                                                                  *
c*******************************************************************


c********************************************************************


        subroutine Calcnsdra(nbasis)


        implicit real*8 (a-h,o-z)

        common/cons/zero,half,tenth,one,two,three,ten
     :  /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /tatb/ta(MN),tb(MN),mtp
     :  /wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :  mf(MNBAS)
     :  /disn/aa(MN)
        common/intnsdmat/rintnsd(MNBAS,MNBAS)
        common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)
        common/kappa/nak(MNBAS)
        common/inparm/nparm
        common/pncp/amass,qw

        do i1 = 1, nbasis
        do i2 = 1, i1
 
        orbi1=orbj(i1)
        orbi2=orbj(i2)
        
        if(nak(i1).lt.0)then

        dli1=orbi1-half
        dlbi1=orbi1+half

        else

        dli1=orbi1+half
        dlbi1=orbi1-half

        endif

        if(nak(i2).lt.0)then

        dli2=orbi2-half
        dlbi2=orbi2+half

        else

        dli2=orbi2+half
        dlbi2=orbi2-half

        endif

        isign1=(-1)**(orbi1+dlbi1-half)
        isign2=(-1)**(orbi1+dli1-half)

        call sixj(orbi2,orbi1,one,half,half,dlbi1,s6j1)
        call sixj(orbi2,orbi1,one,half,half,dli1,s6j2)

        fact1=isign1*sqrt(6*(2*orbi1+1)*(2*orbi2+1))*s6j1
        fact2=isign2*sqrt(6*(2*orbi1+1)*(2*orbi2+1))*s6j2

        rintnsd(i1,i2) = 0.0d0

        if(iparity(i1).ne.iparity(i2))then
        if(dlbi1.eq.dli2)then

        ta(1)=zero

        do 1 l=2,mtp
        call nuc
        ta(l)=((qf(l,i1)*pf(l,i2)*fact1-
     :  pf(l,i1)*qf(l,i2)*fact2)*aa(l)*rp(l))
    1   continue
        call quad(result)

        endif
        endif

        rintnsd(i1,i2) = result
        rintnsd(i2,i1) = -rintnsd(i1,i2)

        enddo
        enddo

        return
        end

c**********************************************************************



c*********************************************************************
c      The reduced matrix element of the NSD PNC matrix element is   *
c      calculated                                                    *
c*********************************************************************

        subroutine nsdmat(i1,i2,jorb1,jorb2,kap1,kap2,ansdm)

        implicit double precision (a-h,o-z)
        double precision j1,j2,jorb1,jorb2

        common/cons/zero,half,tenth,one,two,three,ten
        common/orbital_energy/eorb(MNBAS)
        common/intnsdmat/rintnsd(MNBAS,MNBAS)
        common/inparm/nparm
        common/pncp/amass,qw

        ansdm= zero

        f=1.41421356237
        g=2.2225e-14
        factor=dble(g/f)

        ansd1=rintnsd(i1,i2)

        ansdm=factor*ansd1

        return
        end


c***********************************************************************


c***********************************************************************
c     This subroutine calculated the electric dipole O bar diagram of  *
c     of the type HH                                                   *
c***********************************************************************


      subroutine dipvelhh(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      jgot(i)=0
      igot(i)=0
      enddo

      nbasis=nocc+nexcit

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)

      if(ibsym.eq.iasym)go to 1

c-----------------------------------------------------------------
c     diagram 1  = <b|d|a>
c-----------------------------------------------------------------

      call e1vel(iib,iia,orbb,orba,kapb,kapa,dipole)
      d(iib,iia)=d(iib,iia)+dipole

      sum=zero


      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

c------------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c------------------------------------------------------------------

      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.ne.ibsym)then

      ini=idpa(iip,iia)

      call e1vel(iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)

      endif
      endif
      endif

c------------------------------------------------------------------
c     diagram 3 = <b|t1|p><p|d|a>
c------------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.ne.iasym)then

      inf=idpa(iip,iib)

      call e1vel(iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)

      endif
      endif
      endif

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.eq.ipbsym) then

c---------------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c---------------------------------------------------------------------

      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.ne.iasym)then

      call e1vel(iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

      call findk(orbb,orbp,orbg,orbr,iab,iap,iag,iar,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      if(orbb.eq.orbp)then
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 5 = <bg|t2|rp><p|d|a><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif     
      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym.eq.irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.ne.igsym)then

      call e1vel(iir,iig,orbr,orbg,kapr,kapg,dipole)

c-------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c-------------------------------------------------------------------

      ini=idpa(iip,iia)

      call findk(orbb,orbp,orbg,orbr,iab,iap,iag,iar,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npb=idra(iip,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(npb,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c---------------------------------------------------------------------

      ini=idpa(iip,iia)

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c----------------------------------------------------------------
c     diagrams 8 = <<g|d|r><b|t1|p><pr|t2|ag>
c----------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)

      if(ipasym .eq. irgsym) then
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.ne.igsym)then

      call e1vel(iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      call sixj(one,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      inf=idpa(iip,iib)
      ini=ntloc(npg,nra,ipgsym,l1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

      if(ipasym.eq.irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.ne.ibsym)then

      call e1vel(iib,iip,orbb,orbp,kapb,kapp,dipole)

      inf=idpa(iir,iig)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      if(orbp.eq.orba)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
      if(orba.eq.orbp) then
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      if(ipsym.eq.irsym)go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c------------------------------------------------------------------
c     diagram 12 = <b|t1|r><t|d|p><p|t1|a>
c------------------------------------------------------------------

      call e1vel(iir,iip,orbr,orbp,kapr,kapp,dipole)

      inf=idpa(iir,iib)
      ini=idpa(iip,iia)

      sum=sum+dipole*tf(inf)*ti(ini)

 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

c------------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c------------------------------------------------------------------

      if(igsym.ne.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call e1vel(iib,iig,orbb,orbg,kapb,kapg,dipole)

      ini=idpa(iip,iia)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c---------------------------------------------------------------------


      if(igsym.ne.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call e1vel(iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iip,iib)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nak(iip)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nak(iid)

      if(igsym.eq.idsym)go to 8

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

      call e1vel(iid,iig,orbd,orbg,kapd,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+1+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 10   continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,one,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11   continue
  9   continue
  8   continue

      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      if(iqsym.eq.irsym)go to 13

      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13

      call e1vel(iiq,iir,orbq,orbr,kapq,kapr,dipole)

c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,one,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue

c----------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c----------------------------------------------------------------------

      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(one,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+one+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo

 14   continue
 13   continue

      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.eq.irsym)go to 16

      call e1vel(iir,iip,orbr,orbp,kapr,kapp,dipole)

c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then      
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      call sixj(one,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+one)
      fact=s6j*isign/(2*al1+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,one,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+one)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact

  18  continue
  17  continue
  16  continue
  12  continue


      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nak(iid)

      if(idsym.eq.ibsym)go to 21

      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22

      call e1vel(iib,iid,orbb,orbd,kapb,kapd,dipole)

c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 22=   <gd|t2|pr><b|d|d><pr|t2|ag>
c---------------------------------------------------------------------

      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 24   continue
 23   continue
 22   continue
 21   continue

      if(idsym.eq.iasym)go to 25

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26

      call e1vel(iid,iia,orbd,orba,kapd,kapa,dipole)

c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(npd,nrg,irgsym,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l2)
      ini=ntloc(npd,nrg,irgsym,l1)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+one)

 1    continue

      return
      end

c*****************************************************************
c     This subroutine calculated the hyperfine O bar diagram of  *
c     of the type  HP                                            *
c*************************************************************** *

      subroutine dipvelhp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      nbasis=nocc+nexcit

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nak(iip)

      if(ipsym.eq.iasym)go to 1

c--------------------------------------------------------------
c     diagram 1  = <p|d|a>
c--------------------------------------------------------------

      call e1vel(iip,iia,orbp,orba,kapp,kapa,dipole)

      d(iip,iia)=d(iip,iia)+dipole

      sum=zero

      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      if(igsym.eq.irsym)go to 3

      call e1vel(iig,iir,orbg,orbr,kapg,kapr,dipole)

      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.irgsym)go to 401

c--------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c--------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,ipasym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(2*one+1)
      sum=sum+dipole*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c----------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.eq.0)go to 401
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dipole*ti(ini)*fact

 4    continue
 401  continue

c--------------------------------------------------------------------------
c     diagram 4 = <r|t1|a><r|d|g><p|t1|g>
c--------------------------------------------------------------------------

      if(igsym.eq.irsym)go to 3
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3

      call e1vel(iig,iir,orbg,orbr,kapg,kapr,dipole)

      ini=idpa(iir,iia)
      inj=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*ti(inj)

 3    continue

c------------------------------------------------------------------------
c     diagram 5 =  <r|t1|a><g|t1|r><p|d|g>
c-------------------------------------------------------------------------

      if(ipsym.eq.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5

      call e1vel(iip,iig,orbp,orbg,kapp,kapg,dipole)

      ini=idpa(iir,iia)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 5    continue

c-------------------------------------------------------------------
c     diagram 6 = <r|d|a><g|t|r><p|t|g>
c-------------------------------------------------------------------

      if(irsym.eq.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6

      ini=idpa(iir,iig)
      inf=idpa(iip,iig)

      call e1vel(iir,iia,orbr,orba,kapr,kapa,dipole)

      sum=sum-dipole*ti(ini)*tf(inf)

 6    continue

      if(igsym.ne.irsym)go to 7
      if(orbg.ne.orbr)go to 7

      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      if(ipsym.eq.iqsym)go to 8

      call e1vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)

      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.iqasym)go to 8

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 7 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     sum=sum+dipole*ti(inf)*tf(ini)*dsqrt(2*orbg+1)
      if(orba.eq.orbq)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(2*orbq+1))
      
      endif
      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 8 =  <g|t1|r><p|d|q><qr|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 9    continue
 8    continue

      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)

      if(iasym.eq.ibsym)go to 10

      call e1vel(iib,iia,orbb,orba,kapb,kapa,dipole)

      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)

      if(irgsym.ne.ipbsym)go to 10

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 9 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbb.eq.orbp)then
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbp+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 10 = <g|t1|r><b|d|a><pr|t2|gb>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)then
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
 11   continue
 10   continue
  7   continue

      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)

      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)

      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipsym.eq.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      if(irgsym.ne.iqasym) goto 13

      call e1vel(iip,iib,orbp,orbb,kapp,kapb,dipole)

c---------------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c---------------------------------------------------------------------

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orba,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif


c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue
 14   continue
 13   continue

      if(orbp.ne.orbq)go to 12

      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.eq.iasym)go to 12

c-------------------------------------------------------------------
c     diagram 13 = <gb|t2|rq><pr|t2|bg><q|d|a>
c-------------------------------------------------------------------

      call e1vel(iiq,iia,orbq,orba,kapq,kapa,dipole)

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif


c-------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c-------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 17   continue
 16   continue
 12   continue


      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2

      inf=idpa(iir,iig)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nak(iiq)

      if(irsym.eq.iqsym)go to 20

      call e1vel(iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipasym.ne.iqgsym)go to 20

c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orba,orbp,orbq,orbg,iaa,iap,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(one,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact

 21   continue
 20   continue

      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)

      if(ibsym.eq.igsym)go to 22

      call e1vel(iib,iig,orbb,orbg,kapb,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)

      if(ipasym.ne.irbsym)go to 22

c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbb,iaa,iap,iar,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrb=idra(iir,iib,1)
      ini=ntloc(npa,nrb,ipasym,1)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*1+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(one,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 23   continue
 22   continue
  2   continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------

      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)

      if(iqsym.eq.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24

      ini=idpa(iiq,iia)

      call e1vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)

      sum=sum+ti(ini)*dipole

 24   continue

c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------

      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)

      if(ibsym.eq.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25

      ini=idpa(iip,iib)

      call e1vel(iib,iia,orbb,orba,kapb,kapa,dipole)

      sum=sum-ti(ini)*dipole

 25   continue

      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)

      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)

      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      if(igsym.eq.irsym)go to 26

      call e1vel(iir,iig,orbr,orbg,kapr,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------

      call findk(orba,orbp,orbq,orbb,iaa,iap,iaq,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*1+1)*(2*1+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

 27   continue

c----------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg>
c------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)

      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(1.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2.0*1+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 29   continue

      endif
      endif
      enddo
      endif

 28   continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 30
      do jloop=1,jmax
      l2=jgot(jloop)
      if (l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(one,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2.0*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact

 31   continue

      endif
      endif
      enddo

 30   continue

c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32

      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(one,orbr,orbg,al1,orbq,orbb,s6j1)
      call sixj(one,orba,orbp,al2,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 34   continue
 33   continue
 32   continue
 26   continue

      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+one)

 1    continue

      return
      end

      subroutine dipvelph(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      nbasis=nocc+nexcit

c     P-H type <p|d|a>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)
      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nak(iip)
      if(ipsym.eq.iasym)go to 1
c     diagram 1  = <p|d|a>
      call e1vel(iip,iia,orbp,orba,kapp,kapa,dipole)
      d(iip,iia)=d(iip,iia)+dipole
      sum=zero


c     Diagrams that are linear in T
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)
      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)
      if(igsym.eq.irsym)go to 3
      call e1vel(iig,iir,orbg,orbr,kapg,kapr,dipole)
      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.irgsym)go to 401
c     diagram 2 = sum(r,g) <pr|t2|ag><g|d|r>

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,ipasym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(2*one+1)
      sum=sum+dipole*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 3 = -sum(r,g)<g|d|r><pr|t2|ga>
      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.eq.0)go to 401
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dipole*ti(ini)*fact
 4    continue
 401  continue
c     Non-linear contribution
c     diagram 8 = -sum(r,g) <r|t1|a><g|d|r><p|t1|g>
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(irsym.eq.igsym)go to 3
      if(orbr.ne.orba)go to 3
      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
      call e1vel(iig,iir,orbg,orbr,kapg,kapr,dipole)
      sum=sum-dipole*ti(ini)*ti(inj)
 3    continue
c     diagram 9 = -sum(r,g) <r|t1|a><g|t1|r><p|d|g>
      if(ipsym.eq.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5
      call e1vel(iip,iig,orbp,orbg,kapp,kapg,dipole)
      ini=idpa(iir,iia)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 5    continue
c     diagram 10 = -sum(r,g) <r|d|a><g|t|r><p|t|g>
      if(irsym.eq.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6
      ini=idpa(iir,iig)
      inf=idpa(iip,iig)
      call e1vel(iir,iia,orbr,orba,kapr,kapa,dipole)
      sum=sum-dipole*ti(ini)*tf(inf)
 6    continue

      if(igsym.ne.irsym)go to 7
      if(orbg.ne.orbr)go to 7
c     statement valid for diagrams 4,5,6,7
      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      if(ipsym.eq.iqsym)go to 8
c     statement valid for diagrams 4,5
      call e1vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)
      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.iqasym)go to 8
      inf=idpa(iir,iig)

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     diagram 4 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ag>
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbp+one))
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(2*orbq+1))
      endif
      enddo
      endif

c     diagram 5 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ga>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 9    continue
 8    continue

      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)
      if(iasym.eq.ibsym)go to 10
      call e1vel(iib,iia,orbb,orba,kapb,kapa,dipole)
      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)
      if(irgsym.ne.ipbsym)go to 10
      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
c     diagram 6 = -sum(b,r,g) <g|t1|r><b|d|a><pr|t2|bg>
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbq+one))
      endif
      enddo
      endif

c     diagram 7 = sum(b,r,g) <g|t1|r><b|d|a><pr|t2|gb>
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      if(orbp.ne.orbb)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
  7   continue


c     diagram 11,12, 13 and 14
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)
      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)
      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      if(ipsym.eq.ibsym)go to 13
c     statement valid for diagram 11, 12
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqasym) goto 13
      call e1vel(iip,iib,orbp,orbb,kapp,kapb,dipole)
c     diagram 11 = -sum (q,g,r,b)<gb|t2|rq><qr|t2|ag><p|d|b>
      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
c     diagram 12 = sum (q,g,r,b)<gb|t2|rq><qr|t2|ga><p|d|b>
      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
 14   continue
 13   continue

      if(orbp.ne.orbq)go to 12
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.eq.iasym)go to 12
c     diagram 13 = -sum(q,r,g,b) <gb|t2|rq><pr|t2|bg><q|d|a>
      call e1vel(iiq,iia,orbq,orba,kapq,kapa,dipole)
      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      npb=idra(iip,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
c     diagram 14 = sum(q,r,g,b) <gb|t2|rq><pr|t2|gb><q|d|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 17   continue
 16   continue
 12   continue

c     diagram 15 to 18
      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2
      inf=idpa(iir,iig)
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nak(iiq)
      if(irsym.eq.iqsym)go to 20
      call e1vel(iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(ipasym.ne.iqgsym)go to 20

      call findk(orba,orbp,orbq,orbg,iaa,iap,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(one,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact
 21   continue
 20   continue

c     diagram 17 and 18
      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)
      if(ibsym.eq.igsym)go to 22
      call e1vel(iib,iig,orbb,orbg,kapb,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      if(ipasym.ne.irbsym)go to 22

      call findk(orba,orbp,orbr,orbb,iaa,iap,iar,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrb=idra(iir,iib,1)
      ini=ntloc(npa,nrb,ipasym,1)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*1+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(one,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 23   continue
 22   continue
  2   continue

c     diagram 19
      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)
      if(iqsym.eq.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24
      ini=idpa(iiq,iia)
      call e1vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)
      sum=sum+ti(ini)*dipole
 24   continue
c     diagram 20
      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)
      if(ibsym.eq.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25
      ini=idpa(iip,iib)
      call e1vel(iib,iia,orbb,orba,kapb,kapa,dipole)
      sum=sum-ti(ini)*dipole
 25   continue

c     diagram 21-24
      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)
      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)
      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)
      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)
      if(igsym.eq.irsym)go to 26
c     statement is valid for diagram 21 to 25
      call e1vel(iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 21
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

      call findk(orba,orbp,orbq,orbb,iaa,iap,iaq,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2.0*1+1)*(2.0*1+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif

 27   continue
c     diagram 22
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(1.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2.0*1+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact
 29   continue
      endif
      endif
      enddo
      endif
 28   continue

c     diagram 23
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 30
      do jloop=1,jmax
      l2=jgot(jloop)
      if (l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(1.0d0,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2.0*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2.0*1+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact
 31   continue
      endif
      endif
      enddo
 30   continue
c     diagram 24
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32
      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(1.0d0,orbr,orbg,al2,orbq,orbb,s6j1)
      call sixj(1.0d0,orba,orbp,al1,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 34   continue
 33   continue
 32   continue
 26   continue

      d(iip,iia)=d(iip,iia)+sum
 1    continue
      return
      end


c*****************************************************************
c     This subroutine calculated the hyperfine O bar diagram of  *
c     of the type  PP                                            *
c*************************************************************** *

      subroutine dipvelpp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      nbasis=nocc+nexcit

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nak(iip)

      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)

      if(ipsym.eq.iqsym) go to 1

      call e1vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)

c--------------------------------------------------------------
c     diagram 1 = <p|d|q>
c--------------------------------------------------------------

      d(iip,iiq)=d(iip,iiq)+dipole

      sum=0.0d0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      if(iasym.ne.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then

      call e1vel(iia,iiq,orba,orbq,kapa,kapq,dipole)

      ini=idpa(iip,iia)

c----------------------------------------------------------------
c     diagram 2 = <p|t|a><a|d|q>
c----------------------------------------------------------------

      sum=sum-dipole*ti(ini)

      endif
      endif
      endif

      if(iasym.ne.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then

      call e1vel(iip,iia,orbp,orba,kapp,kapa,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

      sum=sum-dipole*tf(inf)

      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.ne.ipsym)then

      ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
      call e1vel(iip,iia,orbp,orba,kapp,kapa,dipole)
c     diagram 4 = -sum (a,g,r) <ag|t2|qr><r|t|g><p|d|a>
c     sum=sum-dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      if(orbq.eq.orba)then
      sum=sum-dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orba+1))

      endif
      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c---------------------------------------------------------------------

      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba) then
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.ne.irsym)then

      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then

      call e1vel(iir,iig,orbr,orbg,kapr,kapg,dipole)

      ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      nqa=idra(iiq,iia,1)
      inf=ntloc(nrg,nqa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.ne.irsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call e1vel(iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      npa=idra(iip,iia,1)
      ini=ntloc(nrg,npa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.ne.iqsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call e1vel(iia,iiq,orba,orbq,kapa,kapq,dipole)

      inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orba.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(two*orba+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp) then
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

c------------------------------------------------------------------
c     diagram 12 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

      if(irsym.ne.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then

      call e1vel(iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------
c     diagram 13 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------

      if(irsym.ne.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then

      call e1vel(iip,iir,orbp,orbr,kapp,kapr,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  5   continue

c------------------------------------------------------------------
c     diagram 14 = <q|t1|a><a|d|g><g|t1|p
c------------------------------------------------------------------

      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nak(iia)

      if(igsym.ne.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then

      call e1vel(iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)

      sum=sum+dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)

      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nak(iia)

      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nak(iid)

      if(igsym.eq.idsym)go to 9

      call e1vel(iid,iig,orbd,orbg,kapd,kapg,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue

c------------------------------------------------------------------
c     diagram 16 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)

      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,one,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 14   continue
 13   continue
 12   continue
  9   continue

      if(iasym.eq.idsym)go to 8

      call e1vel(iid,iia,orbd,orba,kapd,kapa,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)

      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8

c-----------------------------------------------------------------
c     digram 17 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(one,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+one+l1)/(2*al1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+one+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,one,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 16   continue
 15   continue
  8   continue

c------------------------------------------------------------------
c     diagram 19 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------

      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nak(iis)

      if(issym.eq.irsym)go to 18

      call e1vel(iis,iir,orbs,orbr,kaps,kapr,dipole)

      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18

      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,one,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 20   continue

c------------------------------------------------------------------
c     diagram 20 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,one,orbg,al2,orbp,s9j)
c      fact=(-1)**(orbq+orbg+l2+one+2*orbp)*s9j
      fact=(-1)**(orba+orbg+orbr+orbs+l1+l2)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 21   continue
 19   continue
 18   continue

c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)

      if (irgsym .ne. iassym) goto 22
      if (irgsym .ne. iqasym) goto 22
      if(ipsym.eq.issym)go to 22

      call e1vel(iip,iis,orbp,orbs,kapp,kaps,dipole)

      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 24   continue
 23   continue
 22   continue

c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)

      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.eq.issym)go to 25

      call e1vel(iis,iiq,orbs,orbq,kaps,kapq,dipole)

      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 27   continue
 26   continue
 25   continue
 17   continue
  7   continue

      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+one)

  1   continue

      return
      end


      subroutine quadruvelhh(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/


      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

c     H-H type <b|d|a>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,nocc
c      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)
      if(ibsym.ne.iasym)go to 1
c     diagram 1  = <b|d|a>
      call e2vel(iib,iia,orbb,orba,kapb,kapa,dipole)
      d(iib,iia)=d(iib,iia)+dipole
      sum=zero

c     Diagrams that are linear in T
      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
c     diagram 2 = sum(p) <b|d|p><p|t1|a>
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      ini=idpa(iip,iia)
      call e2vel(iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)
      endif
      endif
      endif
c     diagram 3 = sum(p) <b|t1|p><p|d|a>
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      inf=idpa(iip,iib)
      call e2vel(iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)
      endif
      endif
      endif


c     diagrams nonlinear in T

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.eq.ipbsym) then

c     diagrams 4 = sum(p,r,g) <bg|t2|pr><p|d|a><r|t1|g>
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.iasym)then

      ini=idpa(iir,iig)
      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
      call e2vel(iip,iia,orbp,orba,kapp,kapa,dipole)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      endif
      enddo
      endif


c     diagrams 5 = -sum(p,r,g) <bg|t2|rp><p|d|a><r|t1|g>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)

      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif

      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.eq.igsym)then
      call e2vel(iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 6 = sum(p,g,r)<bg|t2|pr><r|d|g><p|t1|a>

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npb=idra(iip,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(npb,nrg,irgsym,2)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(two*two+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      endif
      enddo
      endif


c     diagram 7 = -sum(p,g,r)<bg|t2|rp><r|d|g><p|t1|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      ini=idpa(iip,iia)
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      enddo
      endif
      endif
      endif
      endif
      endif


      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then
c     diagrams 8 = sum(p,r,g) <<g|d|r><b|t1|p><pr|t2|ag>
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then
      call e2vel(iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)
      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,irgsym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(two*two+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif


c     diagrams 9 = -sum(p,r,g) <<g|d|r><b|t1|p><pr|t2|ga>
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+2)*s6j
      inf=idpa(iip,iib)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      enddo
      endif
      endif
      endif
      endif
      endif


c     diagrams 10 = (p,r,g) <g|t1|r><b|d|p><pr|t2|ag>
      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.ibsym)then
      call e2vel(iib,iip,orbb,orbp,kapb,kapp,dipole)
      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      endif
      enddo
      endif


c     diagrams 11 =-(p,r,g) <g|t1|r><b|d|p><pr|t2|ga>
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)

      if(orba.eq.orbp)then
c     sum=sum-dipole*tf(inf)*tf(inf)*fact
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif

      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      if(ipsym.ne.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c     diagram 11 = sum <p,r) <b|t1|r><t|d|p><p|t1|a>
      call e2vel(iir,iip,orbr,orbp,kapr,kapp,dipole)
      inf=idpa(iir,iib)
      ini=idpa(iip,iia)
      sum=sum+dipole*tf(inf)*ti(ini)
 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
c     diagram 12 = - sum (p,g) <g|t1|p><b|d|g><p|t1|a>
      if(igsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      call e2vel(iib,iig,orbb,orbg,kapb,kapg,dipole)
      ini=idpa(iip,iia)
      inf=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
c     diagram 13 = - sum (p,g) <g|d|a><b|t1|p><p|t1|g>
      if(igsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      call e2vel(iig,iia,orbg,orba,kapg,kapa,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iip,iib)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
 6    continue
 4    continue


c     diagrams 15 and 15 (a)
      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)
      if(igsym.ne.idsym)go to 8
c     diagram 15 = - sum (p,r,g,d) <gb|t2|rp><d|d|g><pr|t2|ad>
      call e2vel(iid,iig,orbd,orbg,kapd,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8
      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,two,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,two,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+2+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 10   continue
c     diagram 15(a) = sum (p,r,g,d) <gb|t2|rp><d|d|g><pr|t2|da>
      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,two,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact
 11   continue
  9   continue
  8   continue

c     diagrams 16 and 17
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      if(iqsym.ne.irsym)go to 13
      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13
      call e2vel(iiq,iir,orbq,orbr,kapq,kapr,dipole)
c     diagram 16 = sum(p,q,r,g) <bg|t2|pq><q|d|r><pr|t2|ag>
      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,two,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,two,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
c     diagram 17 = -sum (p,q,r,g) <bg|t2|pq><q|d|r><pr|t2|ga>
      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(two,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+two+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      enddo
 14   continue
 13   continue
c     diagrams 18 and 19
      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.ne.irsym)go to 16
      call e2vel(iir,iip,orbr,orbp,kapr,kapp,dipole)

c     diagram 18= sum (p,r,q,g) <bg|t2|rq><r|d|p><pq|t2|ag>
      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      call sixj(two,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+two)
      fact=s6j*isign/(2*al1+1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact
      endif
      enddo
      endif

c     diagram 19= -sum (p,r,q,g) <bg|t2|rq><r|d|p><pq|t2|ga>
      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,two,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+two)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact
  18  continue
  17  continue
  16  continue
  12  continue
c     diagrams 20,21,22,23
      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)
      if(idsym.ne.ibsym)go to 21
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22
      call e2vel(iib,iid,orbb,orbd,kapb,kapd,dipole)

c     diagram 20= -sum (p,r,d,g) <gd|t2|rp><b|d|d><pr|t2|ag>
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif

c     diagram 21=  sum (p,r,d,g) <gd|t2|pr><b|d|d><pr|t2|ag>
      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 24   continue
 23   continue
 22   continue
 21   continue
      if(idsym.ne.iasym)go to 25
      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26
      call e2vel(iid,iia,orbd,orba,kapd,kapa,dipole)

c     diagram 22= -sum (p,r,d,g) <bg|t2|pr><d|d|a><pr|t2|dg>
      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact
      endif
      enddo
      endif

c     diagram 23=  sum (p,r,d,g) <gb|t2|pr><d|d|a><pr|t2|dg>
      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact
 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
c     d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+one)
 1    continue
      return
      end

      subroutine quadruvelhp(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo


c     H-P/P-H type <p|d|a>/<a|d|p>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)
      if(ipsym.ne.iasym)go to 1
c     diagram 1  = <p|d|a>
      call e2vel(iip,iia,orbp,orba,kapp,kapa,dipole)
      d(iip,iia)=d(iip,iia)+dipole
      sum=zero


c     Diagrams that are linear in T
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
      if(igsym.ne.irsym)go to 3
      call e2vel(iig,iir,orbg,orbr,kapg,kapr,dipole)
      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.irgsym)go to 401
c     diagram 2 = sum(r,g) <pr|t2|ag><g|d|r>

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,ipasym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(2*two+1)
      sum=sum+dipole*ti(ini)*fact
      endif
      enddo
      endif


c     diagram 3 = -sum(r,g)<g|d|r><pr|t2|ga>
      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum-dipole*ti(ini)*fact
 4    continue
      endif
 401  continue

c     Non-linear contribution
c     diagram 8 = -sum(r,g) <r|t1|a><r|d|g><p|t1|g>
      call e2vel(iig,iir,orbg,orbr,kapg,kapr,dipole)
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3
      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*ti(inj)
 3    continue

c     diagram 9 = -sum(r,g) <r|t1|a><g|t1|r><p|d|g>
      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5
      call e2vel(iip,iig,orbp,orbg,kapp,kapg,dipole)
      ini=idpa(iir,iia)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 5    continue

c     diagram 10 = -sum(r,g) <r|d|a><g|t|r><p|t|g>
      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6
      call e2vel(iir,iia,orbr,orba,kapr,kapa,dipole)
      inf=idpa(iir,iig)
      ini=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 6    continue


      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7
      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      if(ipsym.ne.iqsym)go to 8
      call e2vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)
      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      if(irgsym.ne.iqasym)go to 8
      inf=idpa(iir,iig)

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     diagram 4 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ag>
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))
      endif
      enddo
      endif

c     diagram 5 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ga>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 9    continue
 8    continue


      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      if(iasym.ne.ibsym)go to 10
      call e2vel(iib,iia,orbb,orba,kapb,kapa,dipole)
      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)
      if(irgsym.ne.ipbsym)go to 10
      inf=idpa(iir,iig)

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbp.eq.orbb)then
c     diagram 6 = -sum(b,r,g) <g|t1|r><b|d|a><pr|t2|bg>
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbp+1))
      endif
      endif
      enddo
      endif

c     diagram 7 = sum(b,r,g) <g|t1|r><b|d|a><pr|t2|gb>
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
  7   continue

c     diagram 11,12, 13 and 14
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqasym) goto 13
      if(ipsym.ne.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      call e2vel(iip,iib,orbp,orbb,kapp,kapb,dipole)

c     diagram 11 = -sum (q,g,r,b)<gb|t2|rq><qr|t2|ag><p|d|b>
      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 12 = sum (q,g,r,b)<gb|t2|rq><qr|t2|ga><p|d|b>
      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
 14   continue
 13   continue


      if(orbp.ne.orbq)go to 12
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12
c     diagram 13 = -sum(q,r,g,b) <gb|t2|rq><pr|t2|bg><q|d|a>
      call e2vel(iiq,iia,orbq,orba,kapq,kapa,dipole)

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 14 = sum(q,r,g,b) <gb|t2|rq><pr|t2|gb><q|d|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 17   continue
 16   continue
 12   continue


c     diagram 15 to 18
      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2
      inf=idpa(iir,iig)
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)
      if(irsym.ne.iqsym)go to 20
      call e2vel(iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(ipasym.ne.iqgsym)go to 20

      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+2)/(2*2+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(two,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+2)*s6j
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 21   continue
 20   continue
      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)
      if(ibsym.ne.igsym)go to 22
      call e2vel(iib,iig,orbb,orbg,kapb,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      if(ipasym.ne.irbsym)go to 22

      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(iloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nrb=idra(iir,iib,2)
      ini=ntloc(npa,nrb,ipasym,2)
      fact=(-1)**(orbb+orbr+2*orbb+2)/(2*2+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(two,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+2)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 23   continue
 22   continue
  2   continue


c     diagram 19
      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24
      ini=idpa(iiq,iia)
      call e2vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)
      sum=sum+ti(ini)*dipole
 24   continue
c     diagram 20
      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)
      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25
      ini=idpa(iip,iib)
      call e2vel(iib,iia,orbb,orba,kapb,kapa,dipole)
      sum=sum-ti(ini)*dipole
 25   continue
c     diagram 21-24
      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      if(igsym.ne.irsym)go to 26
      call e2vel(iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 21
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*1+1)*(2*1+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif
 27   continue

c     diagram 22
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(2.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*2+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact
 29   continue
      endif
      endif
      enddo
      endif
 28   continue
c     diagram 23
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      ini=ntloc(npa,nqb,iqbsym,2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(2.0d0,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*2+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact
 31   continue
      endif
      endif
      enddo
      endif
 30   continue
c     diagram 24
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32
      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(2.0d0,orbr,orbg,al2,orbq,orbb,s6j1)
      call sixj(2.0d0,orba,orbp,al1,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 34   continue
 33   continue
 32   continue
 26   continue

      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orbp+orba+one)
 1    continue
      return
      end


      subroutine quadruvelph(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nak(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      nbasis=nocc+nexcit

c     P-H type <p|d|a>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nak(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nak(iip)
      if(ipsym.ne.iasym)go to 1
c     diagram 1  = <p|d|a>
      call e2vel(iip,iia,orbp,orba,kapp,kapa,dipole)
      d(iip,iia)=d(iip,iia)+dipole
      sum=zero

c     Diagrams that are linear in T
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nak(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nak(iig)

      if(igsym.ne.irsym)go to 3
      call e2vel(iig,iir,orbg,orbr,kapg,kapr,dipole)
      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.irgsym)go to 401
c     diagram 2 = sum(r,g) <pr|t2|ag><g|d|r>

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,ipasym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(2*two+1)
      sum=sum+dipole*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 3 = -sum(r,g)<g|d|r><pr|t2|ga>
      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.eq.0)go to 401
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum-dipole*ti(ini)*fact
 4    continue
 401  continue

c     Non-linear contribution
c     diagram 8 = -sum(r,g) <r|t1|a><g|d|r><p|t1|g>
      if(ipsym.ne.igsym)go to 3
      if(irsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3
      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
      call e2vel(iig,iir,orbg,orbr,kapg,kapr,dipole)
      sum=sum-dipole*ti(ini)*ti(inj)
 3    continue
c     diagram 9 = -sum(r,g) <r|t1|a><g|t1|r><p|d|g>
      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5
      call e2vel(iip,iig,orbp,orbg,kapp,kapg,dipole)
      ini=idpa(iir,iia)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 5    continue
c     diagram 10 = -sum(r,g) <r|d|a><g|t|r><p|t|g>
      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6
      ini=idpa(iir,iig)
      inf=idpa(iip,iig)
      call e2vel(iir,iia,orbr,orba,kapr,kapa,dipole)
      sum=sum-dipole*ti(ini)*tf(inf)
 6    continue

      if(igsym.ne.irsym)go to 7
      if(orbg.ne.orbr)go to 7
c     statement valid for diagrams 4,5,6,7
      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      if(ipsym.ne.iqsym)go to 8
c     statement valid for diagrams 4,5
      call e2vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)
      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.iqasym)go to 8
      inf=idpa(iir,iig)

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     diagram 4 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ag>
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbp+one))
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(2*orbq+1))
      endif
      enddo
      endif

c     diagram 5 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ga>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 9    continue
 8    continue

      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)
      if(iasym.ne.ibsym)go to 10
      call e2vel(iib,iia,orbb,orba,kapb,kapa,dipole)
      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)
      if(irgsym.ne.ipbsym)go to 10
      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
c     diagram 6 = -sum(b,r,g) <g|t1|r><b|d|a><pr|t2|bg>
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbq+one))
      endif
      enddo
      endif

c     diagram 7 = sum(b,r,g) <g|t1|r><b|d|a><pr|t2|gb>
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      if(orbp.ne.orbb)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
  7   continue

c     diagram 11,12, 13 and 14
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nak(iiq)
      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nak(iib)
      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      if(ipsym.ne.ibsym)go to 13
c     statement valid for diagram 11, 12
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqasym) goto 13
      call e2vel(iip,iib,orbp,orbb,kapp,kapb,dipole)
c     diagram 11 = -sum (q,g,r,b)<gb|t2|rq><qr|t2|ag><p|d|b>
      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
c     diagram 12 = sum (q,g,r,b)<gb|t2|rq><qr|t2|ga><p|d|b>
      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
 14   continue
 13   continue

      if(orbp.ne.orbq)go to 12
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12
c     diagram 13 = -sum(q,r,g,b) <gb|t2|rq><pr|t2|bg><q|d|a>
      call e2vel(iiq,iia,orbq,orba,kapq,kapa,dipole)
      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      npb=idra(iip,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
c     diagram 14 = sum(q,r,g,b) <gb|t2|rq><pr|t2|gb><q|d|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 17   continue
 16   continue
 12   continue

c     diagram 15 to 18
      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2
      inf=idpa(iir,iig)
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nak(iiq)
      if(irsym.ne.iqsym)go to 20
      call e2vel(iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(ipasym.ne.iqgsym)go to 20

      call findk(orba,orbp,orbq,orbg,iaa,iap,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nqg=idra(iiq,iig,2)
      ini=ntloc(npa,nqg,ipasym,2)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*2+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(two,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+2)*s6j
      sum=sum-dipole*ti(ini)*tf(inf)*fact
 21   continue
 20   continue

c     diagram 17 and 18
      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nak(iib)
      if(ibsym.ne.igsym)go to 22
      call e2vel(iib,iig,orbb,orbg,kapb,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      if(ipasym.ne.irbsym)go to 22

      call findk(orba,orbp,orbr,orbb,iaa,iap,iar,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrb=idra(iir,iib,2)
      ini=ntloc(npa,nrb,ipasym,2)
      fact=(-1)**(orbb+orbr+2*orbb+2)/(2*2+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(two,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+2)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 23   continue
 22   continue
  2   continue

c     diagram 19
      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)
      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24
      ini=idpa(iiq,iia)
      call e2vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)
      sum=sum+ti(ini)*dipole
 24   continue
c     diagram 20
      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)
      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25
      ini=idpa(iip,iib)
      call e2vel(iib,iia,orbb,orba,kapb,kapa,dipole)
      sum=sum-ti(ini)*dipole
 25   continue

c     diagram 21-24
      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nak(iiq)
      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nak(iir)
      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nak(iib)
      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nak(iig)
      if(igsym.ne.irsym)go to 26
c     statement is valid for diagram 21 to 25
      call e2vel(iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 21
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

      call findk(orba,orbp,orbq,orbb,iaa,iap,iaq,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nqb,ipasym,2)
      inf=ntloc(nrg,nqb,iqbsym,2)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2.0*2+1)*(2.0*2+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif

 27   continue
c     diagram 22
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(nrg,nqb,iqbsym,2)
      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(2.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2.0*2+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact
 29   continue
      endif
      endif
      enddo
      endif
 28   continue

c     diagram 23
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 30
      do jloop=1,jmax
      l2=jgot(jloop)
      if (l2.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      ini=ntloc(npa,nqb,iqbsym,2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(2.0d0,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2.0*2+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2.0*2+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact
 31   continue
      endif
      endif
      enddo
 30   continue
c     diagram 24
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32
      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(2.0d0,orbr,orbg,al2,orbq,orbb,s6j1)
      call sixj(2.0d0,orba,orbp,al1,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 34   continue
 33   continue
 32   continue
 26   continue

      d(iip,iia)=d(iip,iia)+sum
 1    continue
      return
      end

      subroutine quadruvelpp(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

c      do 1 iq=1,nexcit
      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym) go to 1
      call e2vel(iip,iiq,orbp,orbq,kapp,kapq,dipole)
c     diagram 1 =<p|d|q>
      d(iip,iiq)=d(iip,iiq)+dipole
      sum=0.0d0

      sum6=0
      sum8=0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      if(iasym.eq.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then
      call e2vel(iia,iiq,orba,orbq,kapa,kapq,dipole)
      ini=idpa(iip,iia)
c     diagram 2 = -sum (a) <p|t|a><a|d|q>
      sum=sum-dipole*ti(ini)
      endif
      endif
      endif

      if(iasym.eq.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then
      call e2vel(iip,iia,orbp,orba,kapp,kapa,dipole)
      inf=idpa(iiq,iia)
c     diagram 3 = -sum (a) <p|d|a><a|t|q>
      sum=sum-dipole*tf(inf)
      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.eq.ipsym)then
      ini=idpa(iir,iig)

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
      call e2vel(iip,iia,orbp,orba,kapp,kapa,dipole)
c     diagram 4 = -sum (a,g,r) <ag|t2|qr><r|t|g><p|d|a>
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      endif
      enddo
      endif

c     diagram 5 =  sum (a,g,r) <ag|t2|rq><r|t|g><p|d|a>
      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 6 and 7
      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.eq.irsym)then
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.eq.iqasym)then
      call e2vel(iir,iig,orbr,orbg,kapr,kapg,dipole)
      ini=idpa(iip,iia)

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      nrg=idra(iir,iig,2)
      nqa=idra(iiq,iia,2)
      inf=ntloc(nrg,nqa,irgsym,2)
      fact=(-1)**(orbg+orbr+2*orbr+2)/(2*two+1)
c     diagram 6 = -sum (a,g,r) <ag|t2|qr><r|d|g><p|t|a>
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 7 =  sum (a,g,r) <ag|t2|rq><r|d|g><p|t|a>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(two,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 8 and 9
      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.eq.irsym)then
      ipasym=mtbl(ipsym,iasym)
      if(irgsym.eq.ipasym)then
      call e2vel(iig,iir,orbg,orbr,kapg,kapr,dipole)
      inf=idpa(iiq,iia)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      nrg=idra(iir,iig,2)
      npa=idra(iip,iia,2)
      ini=ntloc(nrg,npa,irgsym,2)
      fact=(-1)**(orbg+orbr+2*orbr+2)/(2*two+1)
c     diagram 8 = -sum (a,g,r) <g|d|r><a|t|q><pr|t2|ag>
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif
      enddo
      endif

c     diagram 9 =  sum (a,g,r) <g|d|r><a|t|q><pr|t2|ga>
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(two,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 10 and 11

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.eq.iqsym)then
      ipasym=mtbl(ipsym,iasym)
      if(irgsym.eq.ipasym)then
      call e2vel(iia,iiq,orba,orbq,kapa,kapq,dipole)
      inf=idpa(iir,iig)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     diagram 10 = -sum (a,g,r) <g|t|r><a|d|q><pr|t2|ag>
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      endif
      enddo
      endif

c     diagram 11 = -sum (a,g,r) <g|t|r><a|d|q><pr|t2|ga>
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)sum=sum+dipole*tf(inf)*ti(ini)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
c     diagram 13
      if(irsym.eq.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      call e2vel(iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif

c     diagram 14
      if(irsym.eq.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then
      call e2vel(iip,iir,orbp,orbr,kapp,kapr,dipole)
      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif

  5   continue
c     diagram 12
      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      if(igsym.eq.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then
      call e2vel(iig,iia,orbg,orba,kapg,kapa,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)
      sum=sum+dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
 6    continue
 4    continue
c     diagram 15 and 18
      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)
      if(igsym.ne.idsym)go to 9
      call e2vel(iid,iig,orbd,orbg,kapd,kapg,dipole)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,two,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,two,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
c     diagram 20
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)
      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,two,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+two)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 14   continue
 13   continue
 12   continue
  9   continue
      if(iasym.ne.idsym)go to 8
      call e2vel(iid,iia,orbd,orba,kapd,kapa,dipole)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)
      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(two,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+two+l1)
     & /(2*al1+1)
c     diagram 19
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 18
      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+two+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,two,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*ti(ini)*tf(inf)*fact
 16   continue
 15   continue
  8   continue
c     diagram 16
      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nka(iis)
      if(issym.ne.irsym)go to 18
      call e2vel(iis,iir,orbs,orbr,kaps,kapr,dipole)
      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18
      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)
      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,two,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,two,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 20   continue
c     diagram 17
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,two,orbg,al2,orbp,s9j)
      fact=(-1)**(orbq+orbg+l2+two+2*orbp)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 21   continue
 19   continue
 18   continue
c     diagram 21 and 22
      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)
      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.ne.issym)go to 22
      call e2vel(iip,iis,orbp,orbs,kapp,kaps,dipole)
      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 24   continue
 23   continue
 22   continue
c     diagram 23 and 24
      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)
      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.ne.issym)go to 25
      call e2vel(iis,iiq,orbs,orbq,kaps,kapq,dipole)
      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 27   continue
 26   continue
 25   continue
 17   continue
  7   continue

      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+one)
  1   continue
      return
      end

c*****************************************************************
c     This subroutine calculates the E2 hyperfine O bar diagram  *
c     of the type  HH                                            *
c*************************************************************** *


      subroutine quadbhh(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)

      if(ibsym.ne.iasym)go to 1

c------------------------------------------------------------------
c     diagram 1  = <b|d|a>
c------------------------------------------------------------------

      call quadB(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      d(iib,iia)=d(iib,iia)+dipole

      sum=zero

      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

c----------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c----------------------------------------------------------------

      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.eq.ibsym)then

      ini=idpa(iip,iia)

      call quadB(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)

      endif
      endif
      endif

c----------------------------------------------------------------
c     diagram 3 =  <b|t1|p><p|d|a>
c----------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.eq.iasym)then

      inf=idpa(iip,iib)

      call quadB(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)

      endif
      endif
      endif


      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.eq.ipbsym) then

c----------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c----------------------------------------------------------------

      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.iasym)then

      call quadB(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)    

      ini=idpa(iir,iig)

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      if(orbb.eq.orbp)then
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 5 =  <bg|t2|rp><p|d|a><r|t1|g>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.eq.igsym)then

      call quadB(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c--------------------------------------------------------------------

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npb=idra(iip,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(npb,nrg,irgsym,2)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(two*two+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+2)*s6j
      ini=idpa(iip,iia)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then

c------------------------------------------------------------------
c     diagrams 8 = <g|d|r><b|t1|p><pr|t2|ag>
c------------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then

      call quadB(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,irgsym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(two*two+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+2)*s6j
      inf=idpa(iip,iib)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.ibsym)then

      call quadB(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)

      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      if(orbp.eq.orba)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
      if(orba.eq.orbp)then
c     sum=sum-dipole*tf(inf)*tf(inf)*fact
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(ipsym.ne.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c-------------------------------------------------------------------
c     diagram 12 = <b|t1|r><t|d|p><p|t1|a>
c-------------------------------------------------------------------

      call quadB(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

      inf=idpa(iir,iib)
      ini=idpa(iip,iia)

      sum=sum+dipole*tf(inf)*ti(ini)

 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

c-----------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c-----------------------------------------------------------------

      if(igsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call quadB(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ini=idpa(iip,iia)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c-------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c-------------------------------------------------------------------

      if(igsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call quadB(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iip,iib)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue


      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 8

      call quadB(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,two,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,two,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+2+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 10   continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,two,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11   continue
  9   continue
  8   continue

      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(iqsym.ne.irsym)go to 13

      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13

      call quadB(orbq,orbr,iiq,iir,orbq,orbr,kapq,kapr,dipole)

c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,two,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,two,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue

c-------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(two,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+two+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo

 14   continue
 13   continue

      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.ne.irsym)go to 16

      call quadB(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      call sixj(two,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+two)
      fact=s6j*isign/(2*al1+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,two,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+two)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact

  18  continue
  17  continue
  16  continue
  12  continue

      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(idsym.ne.ibsym)go to 21

      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22

      call quadB(orbb,orbd,iib,iid,orbb,orbd,kapb,kapd,dipole)

c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 22= <gd|t2|pr><b|d|d><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 24   continue
 23   continue
 22   continue
 21   continue

      if(idsym.ne.iasym)go to 25

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26

      call quadB(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(npd,nrg,irgsym,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------------

      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+one)

 1    continue

      return
      end


c*****************************************************************
c     This subroutine calculates the E2 hyperfine O bar diagram  *
c     of the type  HP                                            *
c*****************************************************************


      subroutine quadbhp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo


      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)

      if(ipsym.ne.iasym)go to 1

c-------------------------------------------------------------------
c     diagram 1  = <p|d|a>
c-------------------------------------------------------------------

      call quadB(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      d(iip,iia)=d(iip,iia)+dipole

      sum=zero

      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 3

      call quadB(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.irgsym)go to 401

c---------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c---------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,ipasym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(2*two+1)
      sum=sum+dipole*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum-dipole*ti(ini)*fact

 4    continue

      endif

 401  continue

c--------------------------------------------------------------------------
c     diagram 4 = <r|t1|a><r|d|g><p|t1|g>
c--------------------------------------------------------------------------

      if(igsym.ne.irsym)go to 3
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3

      call quadB(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ini=idpa(iir,iia)
      inj=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*ti(inj)

 3    continue

c-------------------------------------------------------------------------
c     diagram 5 =  <r|t1|a><g|t1|r><p|d|g>
c-------------------------------------------------------------------------

      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5

      call quadB(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dipole)

      ini=idpa(iir,iia)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 5    continue

c-------------------------------------------------------------------
c     diagram 6 = <r|d|a><g|t|r><p|t|g>
c-------------------------------------------------------------------

      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6

      call quadB(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dipole)

      inf=idpa(iir,iig)
      ini=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 6    continue


      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7

      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym)go to 8

      call quadB(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.iqasym)go to 8

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 7 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      if(orba.eq.orbq)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))

      endif
      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 8 =  <g|t1|r><p|d|q><qr|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 9    continue
 8    continue


      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      if(iasym.ne.ibsym)go to 10

      call quadB(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)

      if(irgsym.ne.ipbsym)go to 10

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 9 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbp.eq.orbb)then
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 10 = <g|t1|r><b|d|a><pr|t2|gb>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue
  7   continue


      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqasym) goto 13
      if(ipsym.ne.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13

      call quadB(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dipole)

c-----------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orba,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue
 14   continue
 13   continue


      if(orbp.ne.orbq)go to 12

      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12

      call quadB(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dipole)

c--------------------------------------------------------------------
c     diagram 13 =  <gb|t2|rq><pr|t2|bg><q|d|a>
c--------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 17   continue
 16   continue
 12   continue

      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2

      inf=idpa(iir,iig)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)

      if(irsym.ne.iqsym)go to 20

      call quadB(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipasym.ne.iqgsym)go to 20

c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nqg=idra(iiq,iig,2)
      ini=ntloc(npa,nqg,ipasym,2)
      fact=(-1)**(orbg+orbq+2*orbg+2)/(2*2+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(two,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+2)*s6j
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 21   continue
 20   continue

      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)

      if(ibsym.ne.igsym)go to 22

      call quadB(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)

      if(ipasym.ne.irbsym)go to 22

c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nrb=idra(iir,iib,2)
      ini=ntloc(npa,nrb,ipasym,2)
      fact=(-1)**(orbb+orbr+2*orbb+2)/(2*2+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(two,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+2)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 23   continue
 22   continue
  2   continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------

      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24

      ini=idpa(iiq,iia)

      call quadB(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      sum=sum+ti(ini)*dipole

 24   continue

c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------

      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25

      ini=idpa(iip,iib)

      call quadB(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      sum=sum-ti(ini)*dipole

 25   continue

      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 26

      call quadB(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nqb,ipasym,2)
      inf=ntloc(nrg,nqb,iqbsym,2)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*2+1)*(2*2+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

 27   continue

c------------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg>
c------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28


      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(nrg,nqb,iqbsym,2)

      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(two,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*2+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 29   continue

      endif
      endif
      enddo
      endif

 28   continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      ini=ntloc(npa,nqb,iqbsym,2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(two,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2*2+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*2+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact

 31   continue

      endif
      endif
      enddo
      endif

 30   continue

c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32

      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(two,orbr,orbg,al1,orbq,orbb,s6j1)
      call sixj(two,orba,orbp,al2,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 34   continue
 33   continue
 32   continue
 26   continue


      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+one)

 1    continue

      return
      end

c*****************************************************************
c     This subroutine calculates the E2 hyperfine O bar diagram  *
c     of the type  PP                                            *
c*************************************************************** *


      subroutine quadbpp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym) go to 1

      call quadB(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

c-----------------------------------------------------------------
c     diagram 1 = <p|d|q>
c-----------------------------------------------------------------

      d(iip,iiq)=d(iip,iiq)+dipole

      sum=0.0d0

      sum6=0
      sum8=0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      if(iasym.eq.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then

      call quadB(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      ini=idpa(iip,iia)

c----------------------------------------------------------------
c     diagram 2 = <p|t|a><a|d|q>
c----------------------------------------------------------------

      sum=sum-dipole*ti(ini)

      endif
      endif
      endif

      if(iasym.eq.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then

      call quadB(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

      sum=sum-dipole*tf(inf)

      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.eq.ipsym)then

      call quadB(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orbq.eq.orba)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c-----------------------------------------------------------------

      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.eq.irsym)then

      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then

      call quadB(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      nrg=idra(iir,iig,2)
      nqa=idra(iiq,iia,2)
      inf=ntloc(nrg,nqa,irgsym,2)
      fact=(-1)**(orbg+orbr+2*orbr+2)/(2*two+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c----------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(two,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.eq.irsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call quadB(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      nrg=idra(iir,iig,2)
      npa=idra(iip,iia,2)
      ini=ntloc(nrg,npa,irgsym,2)
      fact=(-1)**(orbg+orbr+2*orbr+2)/(2*two+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(two,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.eq.iqsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call quadB(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orba.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)sum=sum+dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

c-----------------------------------------------------------------
c     diagram 12 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

      if(irsym.eq.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then

      call quadB(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------
c     diagram 13 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------

      if(irsym.eq.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then

      call quadB(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  5   continue

c------------------------------------------------------------------
c     diagram 14 = <q|t1|a><a|d|g><g|t1|p>
c------------------------------------------------------------------

      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      if(igsym.eq.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then

      call quadB(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)

      sum=sum+dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 9

      call quadB(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,two,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,two,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue

c-----------------------------------------------------------------
c     diagram 16 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)

      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,two,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+two)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 14   continue
 13   continue
 12   continue
  9   continue

      if(iasym.ne.idsym)go to 8

      call quadB(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)

      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8

c-----------------------------------------------------------------
c     digram 17 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(two,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+two+l1)/(2*al1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+two+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,two,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*ti(ini)*tf(inf)*fact

 16   continue
 15   continue
  8   continue

c------------------------------------------------------------------
c     diagram 19 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------

      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nka(iis)

      if(issym.ne.irsym)go to 18

      call quadB(orbs,orbr,iis,iir,orbs,orbr,kaps,kapr,dipole)

      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18

      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,two,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,two,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+two+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 20   continue

c------------------------------------------------------------------
c     diagram 20 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,two,orbg,al2,orbp,s9j)
c      fact=(-1)**(orbq+orbg+l2+two+2*orbp)*s9j
      fact=(-1)**(orba+orbg+orbr+orbs+l1+l2)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 21   continue
 19   continue
 18   continue

c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)

      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.ne.issym)go to 22

      call quadB(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,dipole)

      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 24   continue
 23   continue
 22   continue

c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)

      if(irgsym.ne.iassym) goto 25
      if(irgsym.ne.iapsym) goto 25
      if(iqsym.ne.issym)go to 25

      call quadB(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,dipole)

      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 27   continue
 26   continue
 25   continue
 17   continue
  7   continue

      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbq+orbp+one)

  1   continue

      return
      end

      subroutine octchh(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV),igot(MXV),jgot(MXV)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two,three/0.0,0.5,1.0,2.0,3.0/

      do i=1,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

c     H-H type <b|d|a>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)
      if(ibsym.ne.iasym)go to 1
c     diagram 1  = <b|d|a>
      call hyperc(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)
      d(iib,iia)=d(iib,iia)+dipole
      sum=zero


c     Diagrams that are linear in T
      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
c     diagram 2 = sum(p) <b|d|p><p|t1|a>
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      ini=idpa(iip,iia)
      call hyperc(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)
      endif
      endif
      endif
c     diagram 3 = sum(p) <b|t1|p><p|d|a>
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      inf=idpa(iip,iib)
      call hyperc(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)
      endif
      endif
      endif


c     diagrams nonlinear in T

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.eq.ipbsym) then

c     diagrams 4 = sum(p,r,g) <bg|t2|pr><p|d|a><r|t1|g>
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.iasym)then

      ini=idpa(iir,iig)
      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
      call hyperc(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      endif
      enddo
      endif

c     diagrams 5 = -sum(p,r,g) <bg|t2|rp><p|d|a><r|t1|g>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)

      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif

      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.eq.igsym)then
      call hyperc(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 6 = sum(p,g,r)<bg|t2|pr><r|d|g><p|t1|a>

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.3)then
      npb=idra(iip,iib,3)
      nrg=idra(iir,iig,3)
      inf=ntloc(npb,nrg,irgsym,3)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+one+2*orbr)/(two*three+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      endif
      enddo
      endif


c     diagram 7 = -sum(p,g,r)<bg|t2|rp><r|d|g><p|t1|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(three,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      ini=idpa(iip,iia)
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      enddo
      endif
      endif
      endif
      endif
      endif


      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then
c     diagrams 8 = sum(p,r,g) <<g|d|r><b|t1|p><pr|t2|ag>
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then
      call hyperc(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)
      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.3)then
      npa=idra(iip,iia,3)
      nrg=idra(iir,iig,3)
      ini=ntloc(npa,nrg,irgsym,3)
      fact=(-1)**(orbg+orbr+one+2*orbr)/(two*three+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif


c     diagrams 9 = -sum(p,r,g) <<g|d|r><b|t1|p><pr|t2|ga>
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(three,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+one)*s6j
      inf=idpa(iip,iib)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      enddo
      endif
      endif
      endif
      endif
      endif


c     diagrams 10 = (p,r,g) <g|t1|r><b|d|p><pr|t2|ag>
      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.ibsym)then
      call hyperc(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      endif
      enddo
      endif


c     diagrams 11 =-(p,r,g) <g|t1|r><b|d|p><pr|t2|ga>
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)

      if(orba.eq.orbp)then
c     sum=sum-dipole*tf(inf)*tf(inf)*fact
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif

      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      if(ipsym.ne.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c     diagram 11 = sum <p,r) <b|t1|r><t|d|p><p|t1|a>
      call hyperc(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)
      inf=idpa(iir,iib)
      ini=idpa(iip,iia)
      sum=sum+dipole*tf(inf)*ti(ini)
 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
c     diagram 12 = - sum (p,g) <g|t1|p><b|d|g><p|t1|a>
      if(igsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      call hyperc(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)
      ini=idpa(iip,iia)
      inf=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
c     diagram 13 = - sum (p,g) <g|d|a><b|t1|p><p|t1|g>
      if(igsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      call hyperc(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iip,iib)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
 6    continue
 4    continue


c     diagrams 15 and 15 (a)
      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)
      if(igsym.ne.idsym)go to 8
c     diagram 15 = - sum (p,r,g,d) <gb|t2|rp><d|d|g><pr|t2|ad>
      call hyperc(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8
      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,three,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,three,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+1+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 10   continue
c     diagram 15(a) = sum (p,r,g,d) <gb|t2|rp><d|d|g><pr|t2|da>
      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,three,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact
 11   continue
  9   continue
  8   continue

c     diagrams 16 and 17
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      if(iqsym.ne.irsym)go to 13
      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13
      call hyperc(orbq,orbr,iiq,iir,orbq,orbr,kapq,kapr,dipole)
c     diagram 16 = sum(p,q,r,g) <bg|t2|pq><q|d|r><pr|t2|ag>
      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,three,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,three,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
c     diagram 17 = -sum (p,q,r,g) <bg|t2|pq><q|d|r><pr|t2|ga>
      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(one,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+one+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      enddo
 14   continue
 13   continue
c     diagrams 18 and 19
      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.ne.irsym)go to 16
      call hyperc(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

c     diagram 18= sum (p,r,q,g) <bg|t2|rq><r|d|p><pq|t2|ag>
      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      call sixj(one,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+one)
      fact=s6j*isign/(2*al1+1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact
      endif
      enddo
      endif

c     diagram 19= -sum (p,r,q,g) <bg|t2|rq><r|d|p><pq|t2|ga>
      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,three,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+one)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact
  18  continue
  17  continue
  16  continue
  12  continue
c     diagrams 20,21,22,23
      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)
      if(idsym.ne.ibsym)go to 21
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22
      call hyperc(orbb,orbd,iib,iid,orbb,orbd,kapb,kapd,dipole)

c     diagram 20= -sum (p,r,d,g) <gd|t2|rp><b|d|d><pr|t2|ag>
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif

c     diagram 21=  sum (p,r,d,g) <gd|t2|pr><b|d|d><pr|t2|ag>
      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 24   continue
 23   continue
 22   continue
 21   continue
      if(idsym.ne.iasym)go to 25
      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipdsym)go to 26

      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26
      call hyperc(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

c     diagram 22= -sum (p,r,d,g) <bg|t2|pr><d|d|a><pr|t2|dg>
      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact
      endif
      enddo
      endif

c     diagram 23=  sum (p,r,d,g) <gb|t2|pr><d|d|a><pr|t2|dg>
      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact
 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(Orba+Orbb+1)
 1    continue
      return
      end


      subroutine octchp(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV),igot(MXV),jgot(MXV)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two,three/0.0,0.5,1.0,2.0,3.0/

      do i=1,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo


c     H-P/P-H type <p|d|a>/<a|d|p>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)
      if(ipsym.ne.iasym)go to 1
c     diagram 1  = <p|d|a>
      call hyperc(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      d(iip,iia)=d(iip,iia)+dipole
      sum=zero

c     Diagrams that are linear in T
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
      if(igsym.ne.irsym)go to 3
      call hyperc(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.irgsym)go to 401
c     diagram 2 = sum(r,g) <pr|t2|ag><g|d|r>

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.3)then
      npa=idra(iip,iia,3)
      nrg=idra(iir,iig,3)
      ini=ntloc(npa,nrg,ipasym,3)
      fact=(-1)**(orbg+orbr+one+2*orbr)/(2*three+1)
      sum=sum+dipole*ti(ini)*fact
      endif
      enddo
      endif


c     diagram 3 = -sum(r,g)<g|d|r><pr|t2|ga>
      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(three,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dipole*ti(ini)*fact
 4    continue
      endif
 401  continue

c     Non-linear contribution
c     diagram 8 = -sum(r,g) <r|t1|a><r|d|g><p|t1|g>
      call hyperc(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3
      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*ti(inj)
 3    continue

c     diagram 9 = -sum(r,g) <r|t1|a><g|t1|r><p|d|g>
      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5
      call hyperc(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dipole)
      ini=idpa(iir,iia)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 5    continue

c     diagram 10 = -sum(r,g) <r|d|a><g|t|r><p|t|g>
      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6
      call hyperc(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dipole)
      inf=idpa(iir,iig)
      ini=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 6    continue

      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7
      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      if(ipsym.ne.iqsym)go to 8
      call hyperc(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      if(irgsym.ne.iqasym)go to 8
      inf=idpa(iir,iig)

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     diagram 4 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ag>
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))
      endif
      enddo
      endif

c     diagram 5 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ga>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 9    continue
 8    continue

      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      if(iasym.ne.ibsym)go to 10
      call hyperc(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)
      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)
      if(irgsym.ne.ipbsym)go to 10
      inf=idpa(iir,iig)

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbb.eq.orbp)then
c     diagram 6 = -sum(b,r,g) <g|t1|r><b|d|a><pr|t2|bg>
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))
      endif
      endif
      enddo
      endif

c     diagram 7 = sum(b,r,g) <g|t1|r><b|d|a><pr|t2|gb>
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
  7   continue


c     diagram 11,12, 13 and 14
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqasym) goto 13
      if(ipsym.ne.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      call hyperc(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dipole)

c     diagram 11 = -sum (q,g,r,b)<gb|t2|rq><qr|t2|ag><p|d|b>
      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 12 = sum (q,g,r,b)<gb|t2|rq><qr|t2|ga><p|d|b>
      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
 14   continue
 13   continue


      if(orbp.ne.orbq)go to 12
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12
c     diagram 13 = -sum(q,r,g,b) <gb|t2|rq><pr|t2|bg><q|d|a>
      call hyperc(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dipole)

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 14 = sum(q,r,g,b) <gb|t2|rq><pr|t2|gb><q|d|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 17   continue
 16   continue
 12   continue


c     diagram 15 to 18
      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2
      inf=idpa(iir,iig)
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)
      if(irsym.ne.iqsym)go to 20
      call hyperc(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(ipasym.ne.iqgsym)go to 20

      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.3)then
      npa=idra(iip,iia,3)
      nqg=idra(iiq,iig,3)
      ini=ntloc(npa,nqg,ipasym,3)
      fact=(-1)**(orbg+orbq+2*orbg+3)/(2*3+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(three,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 21   continue
 20   continue
      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)
      if(ibsym.ne.igsym)go to 22
      call hyperc(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      if(ipasym.ne.irbsym)go to 22

      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(iloop)
      if(l2.eq.3)then
      npa=idra(iip,iia,3)
      nrb=idra(iir,iib,3)
      ini=ntloc(npa,nrb,ipasym,3)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*3+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(three,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 23   continue
 22   continue
  2   continue


c     diagram 19
      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24
      ini=idpa(iiq,iia)
      call hyperc(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
      sum=sum+ti(ini)*dipole
 24   continue
c     diagram 20
      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)
      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25
      ini=idpa(iip,iib)
      call hyperc(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)
      sum=sum-ti(ini)*dipole
 25   continue
c     diagram 21-24
      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      if(igsym.ne.irsym)go to 26
      call hyperc(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 21
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.3)then
      npa=idra(iip,iia,3)
      nqb=idra(iiq,iib,3)
      nrg=idra(iir,iig,3)
      ini=ntloc(npa,nqb,ipasym,3)
      inf=ntloc(nrg,nqb,iqbsym,3)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*3+1)*(2*3+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif
 27   continue

c     diagram 22
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.3)then
      nqb=idra(iiq,iib,3)
      nrg=idra(iir,iig,3)
      inf=ntloc(nrg,nqb,iqbsym,3)
      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(3.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*1+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact
 29   continue
      endif
      endif
      enddo
      endif
 28   continue
c     diagram 23
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.3)then
      npa=idra(iip,iia,3)
      nqb=idra(iiq,iib,3)
      ini=ntloc(npa,nqb,iqbsym,3)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(3.0d0,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*3+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact
 31   continue
      endif
      endif
      enddo
      endif
 30   continue
c     diagram 24
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32
      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(3.0d0,orbr,orbg,al2,orbq,orbb,s6j1)
      call sixj(3.0d0,orba,orbp,al1,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 34   continue
 33   continue
 32   continue
 26   continue


      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(Orba+Orbp+1)
 1    continue
      return
      end


      subroutine octcpp(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV),igot(MXV),jgot(MXV)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two,three/0.0,0.5,1.0,2.0,3.0/

      do i=1,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)
      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym) go to 1
      call hyperc(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
c     diagram 1 =<p|d|q>
      d(iip,iiq)=d(iip,iiq)+dipole
      sum=0.0d0

      sum6=0
      sum8=0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      if(iasym.eq.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then
      call hyperc(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)
      ini=idpa(iip,iia)
c     diagram 2 = -sum (a) <p|t|a><a|d|q>
      sum=sum-dipole*ti(ini)
      endif
      endif
      endif

      if(iasym.eq.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then
      call hyperc(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      inf=idpa(iiq,iia)
c     diagram 3 = -sum (a) <p|d|a><a|t|q>
      sum=sum-dipole*tf(inf)
      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.eq.ipsym)then
      ini=idpa(iir,iig)

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
      call hyperc(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
c     diagram 4 = -sum (a,g,r) <ag|t2|qr><r|t|g><p|d|a>
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      endif
      enddo
      endif

c     diagram 5 =  sum (a,g,r) <ag|t2|rq><r|t|g><p|d|a>
      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 6 and 7
      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.eq.irsym)then
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.eq.iqasym)then
      call hyperc(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)
      ini=idpa(iip,iia)

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.3)then
      nrg=idra(iir,iig,3)
      nqa=idra(iiq,iia,3)
      inf=ntloc(nrg,nqa,irgsym,3)
      fact=(-1)**(orbg+orbr+2*orbr+3)/(2*three+1)
c     diagram 6 = -sum (a,g,r) <ag|t2|qr><r|d|g><p|t|a>
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 7 =  sum (a,g,r) <ag|t2|rq><r|d|g><p|t|a>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(three,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 8 and 9
      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.eq.irsym)then
      ipasym=mtbl(ipsym,iasym)
      if(irgsym.eq.ipasym)then
      call hyperc(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      inf=idpa(iiq,iia)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.3)then
      nrg=idra(iir,iig,3)
      npa=idra(iip,iia,3)
      ini=ntloc(nrg,npa,irgsym,3)
      fact=(-1)**(orbg+orbr+2*orbr+one)/(2*three+1)
c     diagram 8 = -sum (a,g,r) <g|d|r><a|t|q><pr|t2|ag>
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif
      enddo
      endif

c     diagram 9 =  sum (a,g,r) <g|d|r><a|t|q><pr|t2|ga>
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(three,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 10 and 11

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.eq.iqsym)then
      ipasym=mtbl(ipsym,iasym)
      if(irgsym.eq.ipasym)then
      call hyperc(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)
      inf=idpa(iir,iig)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     diagram 10 = -sum (a,g,r) <g|t|r><a|d|q><pr|t2|ag>
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      endif
      enddo
      endif

c     diagram 11 = -sum (a,g,r) <g|t|r><a|d|q><pr|t2|ga>
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)sum=sum+dipole*tf(inf)*ti(ini)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
c     diagram 13
      if(irsym.eq.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      call hyperc(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif

c     diagram 14
      if(irsym.eq.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then
      call hyperc(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,dipole)
      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif

  5   continue
c     diagram 12
      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      if(igsym.eq.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then
      call hyperc(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)
      sum=sum+dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
 6    continue
 4    continue
c     diagram 15 and 18
      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)
      if(igsym.ne.idsym)go to 9
      call hyperc(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,three,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,three,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
c     diagram 20
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)
      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,three,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 14   continue
 13   continue
 12   continue
  9   continue
      if(iasym.ne.idsym)go to 8
      call hyperc(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)
      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(three,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+one+l1)
     & /(2*al1+1)
c     diagram 19
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 18
      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+one+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,three,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*ti(ini)*tf(inf)*fact
 16   continue
 15   continue
  8   continue
c     diagram 16
      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nka(iis)
      if(issym.ne.irsym)go to 18
      call hyperc(orbs,orbr,iis,iir,orbs,orbr,kaps,kapr,dipole)
      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18
      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)
      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,three,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,three,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 20   continue
c     diagram 17
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,three,orbg,al2,orbp,s9j)
      fact=(-1)**(orbq+orbg+l2+one+2*orbp)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 21   continue
 19   continue
 18   continue
c     diagram 21 and 22
      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)
      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.ne.issym)go to 22
      call hyperc(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,dipole)
      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 24   continue
 23   continue
 22   continue
c     diagram 23 and 24
      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)
      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.ne.issym)go to 25
      call hyperc(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,dipole)
      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 27   continue
 26   continue
 25   continue
 17   continue
  7   continue


      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+1)
  1   continue
      return
      end

c*********************************************************
      subroutine hyperc(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp)
c*********************************************************
c  this subroutine calculates the edm t-pt reduced
c  matrix elements between orbitals with j=1/2
c*********************************************************
      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
      common/intoctc/rint(MNBAS,MNBAS)

      factor=(-1)**(j1+0.5)
      w3j1=dr(j1,three,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      hyp=-result*e1*(kap1+kap2)/three

      return
      end

c*****************************************************************

c*****************************************************************
c       This subroutine calculates the reduced                   *
c       matrix elements of the m1 transition operator            *
c*****************************************************************

      subroutine m1prop(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp)

      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
      common/intm1/rint(MNBAS,MNBAS)

      data alpha/137.0359895/

      factor=(-1)**(j1+0.5)
      w3j1=dr(j1,one,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      hyp=-result*e1*(kap1+kap2)*alpha

      return
      end

c*********************************************************
      subroutine m2prop(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp)
c*********************************************************
c  this subroutine calculates the edm t-pt reduced
c  matrix elements between orbitals with j=1/2
c*********************************************************
      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
      common/intm2/rint(MNBAS,MNBAS)

      data alpha/137.0359895/

      factor=(-1)**(j1+0.5)*two/three
      w3j1=dr(j1,two,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      hyp=-result*e1*(kap1+kap2)*alpha

      return
      end

c*********************************************************
      subroutine gfact(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp,hypc)
c*********************************************************
c  this subroutine calculates the edm t-pt reduced
c  matrix elements between orbitals with j=1/2
c*********************************************************
      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
      common/intg/rint(MNBAS,MNBAS)
      common/intgc/rintc(MNBAS,MNBAS)
      data VL/137.03599976D+00/
      data CORR/0.00116D+00/

      factor=(-1)**(j1+0.5)
      w3j1=dr(j1,one,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      resultc=rintc(i1,i2)
      hyp=-result*e1*(kap1+kap2)*VL
      hypc=result*e1*(kap1+kap2-one)*CORR

      return
      end


c*********************************************************

c*********************************************************
c                                                        *
c  This subroutine calculates the reduced matrix element *
c  of the E2 hyperfine operator                          *
c                                                        *
c*********************************************************


      subroutine quadB(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp)


      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
      common/intB/rint(MNBAS,MNBAS)

      factor=(-1)**(j1+0.5)
      w3j1=dr(j1,two,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      hyp=-result*e1

      return
      end

c*****************************************************************

c*****************************************************************
c     This subroutine calculates the M1 transition O bar diagram *
c     of the type  HH                                            *
c*************************************************************** *


      subroutine m1prophh(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)

      if(ibsym.ne.iasym)go to 1

c------------------------------------------------------------------
c     diagram 1  = <b|d|a>
c------------------------------------------------------------------

      call m1prop(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      d(iib,iia)=d(iib,iia)+dipole

      sum=zero

      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

c----------------------------------------------------------------
c     diagram 2 = <b|d|p><p|t1|a>
c----------------------------------------------------------------

      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.eq.ibsym)then

      ini=idpa(iip,iia)

      call m1prop(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)

      endif
      endif
      endif

c----------------------------------------------------------------
c     diagram 3 =  <b|t1|p><p|d|a>
c----------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.eq.iasym)then

      inf=idpa(iip,iib)

      call m1prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)

      endif
      endif
      endif

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.eq.ipbsym) then

c----------------------------------------------------------------
c     diagrams 4 = <bg|t2|pr><p|d|a><r|t1|g>
c----------------------------------------------------------------

      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.iasym)then

      call m1prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      if(orbb.eq.orbp)then
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 5 =  <bg|t2|rp><p|d|a><r|t1|g>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)
      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.eq.igsym)then

      call m1prop(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 6 = <bg|t2|pr><r|d|g><p|t1|a>
c--------------------------------------------------------------------

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npb=idra(iip,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(npb,nrg,irgsym,1)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif


c--------------------------------------------------------------------
c     diagram 7 = <bg|t2|rp><r|d|g><p|t1|a>
c--------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      ini=idpa(iip,iia)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then


c------------------------------------------------------------------
c     diagrams 8 = <g|d|r><b|t1|p><pr|t2|ag>
c------------------------------------------------------------------

      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then

      call m1prop(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     diagrams 9 = <g|d|r><b|t1|p><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      inf=idpa(iip,iib)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


c------------------------------------------------------------------
c     diagrams 10 = <g|t1|r><b|d|p><pr|t2|ag>
c------------------------------------------------------------------

      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.ibsym)then

      call m1prop(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)

      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      if(orbp.eq.orba)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))

      endif
      endif
      enddo
      endif

c--------------------------------------------------------------------
c     diagrams 11 = <g|t1|r><b|d|p><pr|t2|ga>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)
      if(orba.eq.orbp)then
c     sum=sum-dipole*tf(inf)*tf(inf)*fact
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(ipsym.ne.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c-------------------------------------------------------------------
c     diagram 12 = <b|t1|r><t|d|p><p|t1|a>
c-------------------------------------------------------------------

      call m1prop(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

      inf=idpa(iir,iib)
      ini=idpa(iip,iia)

      sum=sum+dipole*tf(inf)*ti(ini)

 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

c-----------------------------------------------------------------
c     diagram 13 = <g|t1|p><b|d|g><p|t1|a>
c-----------------------------------------------------------------

      if(igsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call m1prop(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ini=idpa(iip,iia)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c-------------------------------------------------------------------
c     diagram 14 = <g|d|a><b|t1|p><p|t1|g>
c-------------------------------------------------------------------

      if(igsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then

      call m1prop(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iip,iib)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue


      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 8

      call m1prop(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

c--------------------------------------------------------------------
c     diagram 15 = <gb|t2|rp><d|d|g><pr|t2|ad>
c--------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+1+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 10   continue

c---------------------------------------------------------------------
c     diagram 16 = <gb|t2|rp><d|d|g><pr|t2|da>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,one,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact

 11   continue
  9   continue
  8   continue


      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(iqsym.ne.irsym)go to 13

      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13

      call m1prop(orbq,orbr,iiq,iir,orbq,orbr,kapq,kapr,dipole)

c-------------------------------------------------------------------
c     diagram 17 = <bg|t2|pq><q|d|r><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,one,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue

c-------------------------------------------------------------------
c     diagram 18 = <bg|t2|pq><q|d|r><pr|t2|ga>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(one,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+one+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      enddo

 14   continue
 13   continue

      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.ne.irsym)go to 16

      call m1prop(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

c------------------------------------------------------------------
c     diagram 19= <bg|t2|rq><r|d|p><pq|t2|ag>
c------------------------------------------------------------------

      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      call sixj(one,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+one)
      fact=s6j*isign/(2*al1+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 20= <bg|t2|rq><r|d|p><pq|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,one,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+one)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact

  18  continue
  17  continue
  16  continue
  12  continue

      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)

      if(idsym.ne.ibsym)go to 21

      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22

      call m1prop(orbb,orbd,iib,iid,orbb,orbd,kapb,kapd,dipole)

c-------------------------------------------------------------------
c     diagram 21= <gd|t2|rp><b|d|d><pr|t2|ag>
c-------------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 22= <gd|t2|pr><b|d|d><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 24   continue
 23   continue
 22   continue
 21   continue

      if(idsym.ne.iasym)go to 25

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)

      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26

      call m1prop(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

c-----------------------------------------------------------------
c     diagram 23= <bg|t2|pr><d|d|a><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(npd,nrg,irgsym,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------------
c     diagram 24= <gb|t2|pr><d|d|a><pr|t2|dg>
c------------------------------------------------------------------------

      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact

 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(orba+orbb+one)

 1    continue

      return
      end


c*****************************************************************
c     This subroutine calculated the M1 transition O bar diagram *
c     of the type  HP                                            *
c*************************************************************** *

      subroutine m1prophp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)

      if(ipsym.ne.iasym)go to 1

c-------------------------------------------------------------------
c     diagram 1  = <p|d|a>
c-------------------------------------------------------------------

      call m1prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      d(iip,iia)=d(iip,iia)+dipole

      sum=zero

      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 3

      call m1prop(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.irgsym)go to 401

c---------------------------------------------------------------------
c     diagram 2 = <pr|t2|ag><g|d|r>
c---------------------------------------------------------------------

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,ipasym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(2*one+1)
      sum=sum+dipole*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------------
c     diagram 3 = <g|d|r><pr|t2|ga>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dipole*ti(ini)*fact

 4    continue

      endif

 401  continue

c--------------------------------------------------------------------------
c     diagram 4 = <r|t1|a><r|d|g><p|t1|g>
c--------------------------------------------------------------------------

      if(igsym.ne.irsym)go to 3
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3

      call m1prop(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      ini=idpa(iir,iia)
      inj=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*ti(inj)

 3    continue

c-------------------------------------------------------------------------
c     diagram 5 =  <r|t1|a><g|t1|r><p|d|g>
c-------------------------------------------------------------------------

      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5

      call m1prop(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dipole)

      ini=idpa(iir,iia)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 5    continue

c-------------------------------------------------------------------
c     diagram 6 = <r|d|a><g|t|r><p|t|g>
c-------------------------------------------------------------------

      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6

      call m1prop(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iip,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

 6    continue


      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7

      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym)go to 8

      call m1prop(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.iqasym)go to 8

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 7 =  <g|t1|r><p|d|q><qr|t2|ag>
c----------------------------------------------------------------

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      if(orba.eq.orbq)then
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))

      endif
      endif
      enddo
      endif

c---------------------------------------------------------------------
c     diagram 8 =  <g|t1|r><p|d|q><qr|t2|ga>
c---------------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 9    continue
 8    continue


      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      if(iasym.ne.ibsym)go to 10

      call m1prop(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)

      if(irgsym.ne.ipbsym)go to 10

      inf=idpa(iir,iig)

c----------------------------------------------------------------
c     diagram 9 = <g|t1|r><b|d|a><pr|t2|bg>
c----------------------------------------------------------------

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbb.eq.orbp)then
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbp+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 10 = <g|t1|r><b|d|a><pr|t2|gb>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue
  7   continue

      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)

      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)

      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqasym) goto 13
      if(ipsym.ne.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13

      call m1prop(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dipole)

c-----------------------------------------------------------------
c     diagram 11 = <gb|t2|rq><qr|t2|ag><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orba,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 12 = <gb|t2|rq><qr|t2|ga><p|d|b>
c-----------------------------------------------------------------

      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 15   continue
 14   continue
 13   continue


      if(orbp.ne.orbq)go to 12

      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)

      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12

      call m1prop(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dipole)

c--------------------------------------------------------------------
c     diagram 13 =  <gb|t2|rq><pr|t2|bg><q|d|a>
c--------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 14 = <gb|t2|rq><pr|t2|gb><q|d|a>
c------------------------------------------------------------------

      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 17   continue
 16   continue
 12   continue


      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2

      inf=idpa(iir,iig)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)

      if(irsym.ne.iqsym)go to 20

      call m1prop(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(ipasym.ne.iqgsym)go to 20

c-------------------------------------------------------------------
c     diagram 15= <pq|t2|ag><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 16 = <pq|t2|ga><q|d|r><r|t1|g>
c------------------------------------------------------------------

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(one,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 21   continue
 20   continue

      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)

      if(ibsym.ne.igsym)go to 22

      call m1prop(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)

      if(ipasym.ne.irbsym)go to 22

c-------------------------------------------------------------------
c     diagram 17= <pr|t2|ab><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nrb=idra(iir,iib,1)
      ini=ntloc(npa,nrb,ipasym,1)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*1+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c    diagram 18 = <pr|t2|ba><b|d|g><r|t1|g>
c-------------------------------------------------------------------

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(one,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 23   continue
 22   continue
  2   continue

c--------------------------------------------------------------
c     diagram 19 = <a|t1|q><q|d|p>
c--------------------------------------------------------------

      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24

      ini=idpa(iiq,iia)

      call m1prop(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

      sum=sum+ti(ini)*dipole

 24   continue

c------------------------------------------------------------------
c     diagram 20 = <p|t1|b><b|d|a>
c------------------------------------------------------------------

      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25

      ini=idpa(iip,iib)

      call m1prop(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)

      sum=sum-ti(ini)*dipole

 25   continue

      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)

      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(igsym.ne.irsym)go to 26

      call m1prop(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)

      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

c-------------------------------------------------------------------
c     diagram 21 = <pq|t2|ab><r|d|g><qr|t2|bg>
c-------------------------------------------------------------------

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*1+1)*(2*1+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact

      endif
      enddo
      endif

 27   continue

c------------------------------------------------------------------
c     diagram 22 =  <pq|t2|ba><r|d|g><qr|t2|bg>
c------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)

      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(1.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*1+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact

 29   continue

      endif
      endif
      enddo
      endif

 28   continue

c-------------------------------------------------------------------
c     diagram 23 = <pq|t2|ab><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(1.0d0,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact

 31   continue

      endif
      endif
      enddo
      endif

 30   continue

c-------------------------------------------------------------------
c     diagram 24 = <pq|t2|ba><r|d|g><qr|t2|gb>
c-------------------------------------------------------------------

      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)

      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32

      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)

      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(1.0d0,orbr,orbg,al1,orbq,orbb,s6j1)
      call sixj(1.0d0,orba,orbp,al2,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact

 34   continue
 33   continue
 32   continue
 26   continue


      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(orba+orbp+one)

 1    continue

      return
      end

c*****************************************************************
c     This subroutine calculates the M1 transition O bar diagram *
c     of the type  PP                                            *
c*************************************************************** *

      subroutine m1proppp(d,ti,tf)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym) go to 1

      call m1prop(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)

c-----------------------------------------------------------------
c     diagram 1 = <p|d|q>
c-----------------------------------------------------------------

      d(iip,iiq)=d(iip,iiq)+dipole

      sum=0.0d0

      sum6=0
      sum8=0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      if(iasym.eq.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then

      call m1prop(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      ini=idpa(iip,iia)

c----------------------------------------------------------------
c     diagram 2 = <p|t|a><a|d|q>
c----------------------------------------------------------------

      sum=sum-dipole*ti(ini)

      endif
      endif
      endif

      if(iasym.eq.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then

      call m1prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 3 = <p|d|a><a|t|q>
c----------------------------------------------------------------

      sum=sum-dipole*tf(inf)

      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.eq.ipsym)then

      call m1prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)

      ini=idpa(iir,iig)

c------------------------------------------------------------------
c     diagram 4 = <ag|t2|qr><r|t|g><p|d|a>
c------------------------------------------------------------------

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orbq.eq.orba)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 5 = <ag|t2|rq><r|t|g><p|d|a>
c-----------------------------------------------------------------

      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)then
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.eq.irsym)then

      iqasym=mtbl(iqsym,iasym)

      if(irgsym.eq.iqasym)then

      call m1prop(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)

      ini=idpa(iip,iia)

c-----------------------------------------------------------------
c     diagram 6 = <ag|t2|qr><r|d|g><p|t|a>
c-----------------------------------------------------------------

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      nqa=idra(iiq,iia,1)
      inf=ntloc(nrg,nqa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------
c     diagram 7 = <ag|t2|rq><r|d|g><p|t|a>
c----------------------------------------------------------------

      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.eq.irsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call m1prop(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iiq,iia)

c----------------------------------------------------------------
c     diagram 8 = <g|d|r><a|t|q><pr|t2|ag>
c----------------------------------------------------------------

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      npa=idra(iip,iia,1)
      ini=ntloc(nrg,npa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
      sum=sum-dipole*ti(ini)*tf(inf)*fact

      endif
      enddo
      endif

c----------------------------------------------------------------
c     diagram 9 = <g|d|r><a|t|q><pr|t2|ga>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.eq.iqsym)then

      ipasym=mtbl(ipsym,iasym)

      if(irgsym.eq.ipasym)then

      call m1prop(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)

      inf=idpa(iir,iig)

c-----------------------------------------------------------------
c     diagram 10 = <g|t|r><a|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      if(orba.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))

      endif
      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 11 = <g|t|r><a|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)then
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif
      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

c-----------------------------------------------------------------
c     diagram 12 = <r|t1|g><r|d|q><g|t1|p>
c-----------------------------------------------------------------

      if(irsym.eq.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then

      call m1prop(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iir,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

c---------------------------------------------------------------
c     diagram 13 = <r|t1|g><p|d|r><g|t1|q>
c---------------------------------------------------------------

      if(irsym.eq.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then

      call m1prop(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,dipole)

      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)

      sum=sum-dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

  5   continue

      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      if(igsym.eq.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then

      call m1prop(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)

      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)

      sum=sum+dipole*ti(ini)*tf(inf)

      endif
      endif
      endif
      endif
      endif

 6    continue
 4    continue

      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)

      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)

      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      if(igsym.ne.idsym)go to 9

      call m1prop(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)

      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9

c------------------------------------------------------------------
c     diagram 15 = <qr|t2|ag><g|d|d><pr|t2|ad>
c------------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 11   continue
 10   continue

c-----------------------------------------------------------------
c     diagram 16 = <qr|t2|ag><g|d|d><pr|t2|da>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)

      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,one,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 14   continue
 13   continue
 12   continue
  9   continue

      if(iasym.ne.idsym)go to 8

      call m1prop(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)

      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8

c-----------------------------------------------------------------
c     digram 17 = <qr|t2|ag><a|d|d><pr|t2|dg>
c-----------------------------------------------------------------

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(one,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+one+l1)/(2*al1+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 18 = <qr|t2|ag><a|d|d><pr|t2|gd>
c-----------------------------------------------------------------

      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+one+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,one,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*ti(ini)*tf(inf)*fact

 16   continue
 15   continue
  8   continue

c------------------------------------------------------------------
c     diagram 19 = <qs|t2|ag><s|d|r><pr|t2|ag>
c------------------------------------------------------------------

      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nka(iis)

      if(issym.ne.irsym)go to 18

      call m1prop(orbs,orbr,iis,iir,orbs,orbr,kaps,kapr,dipole)

      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18

      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,one,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact

 20   continue

c------------------------------------------------------------------
c     diagram 20 = <qs|t2|ag><s|d|r><pr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,one,orbg,al2,orbp,s9j)
c      fact=(-1)**(orbq+orbg+l2+one+2*orbp)*s9j
      fact=(-1)**(orba+orbg+orbr+orbs+l1+l2)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 21   continue
 19   continue
 18   continue

c------------------------------------------------------------------
c     diagram 21 = <qr|t2|ag><p|d|s><sr|t2|ag>
c------------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)

      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.ne.issym)go to 22

      call m1prop(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,dipole)

      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     diagram 22 = <qr|t2|ag><p|d|s><sr|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 24   continue
 23   continue
 22   continue

c-----------------------------------------------------------------
c     diagram 23 = <sr|t2|ag><s|d|q><pr|t2|ag>
c-----------------------------------------------------------------

      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)

      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.ne.issym)go to 25

      call m1prop(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,dipole)

      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     diagram 24 = <sr|t2|ag><s|d|q><pr|t2|ga>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact

 27   continue
 26   continue
 25   continue
 17   continue
  7   continue

      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+1)

  1   continue

      return
      end

      subroutine m2prophh(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

c     H-H type <b|d|a>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)
      if(ibsym.eq.iasym)go to 1
c     diagram 1  = <b|d|a>
      call m2prop(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)
      d(iib,iia)=d(iib,iia)+dipole
      sum=zero


c     Diagrams that are linear in T
      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
c     diagram 2 = sum(p) <b|d|p><p|t1|a>
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.ne.ibsym)then
      ini=idpa(iip,iia)
      call m2prop(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      sum=sum+dipole*ti(ini)
      endif
      endif
      endif
c     diagram 3 = sum(p) <b|t1|p><p|d|a>
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.ne.iasym)then
      inf=idpa(iip,iib)
      call m2prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      sum=sum+dipole*tf(inf)
      endif
      endif
      endif


c     diagrams nonlinear in T

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.eq.ipbsym) then

c     diagrams 4 = sum(p,r,g) <bg|t2|pr><p|d|a><r|t1|g>
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.ne.iasym)then

      ini=idpa(iir,iig)
      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
      call m2prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
c     sum=sum+dipole*ti(ini)*tf(inf)*dsqrt(two*orbg+1)
      sum=sum+dipole*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      endif
      enddo
      endif


c     diagrams 5 = -sum(p,r,g) <bg|t2|rp><p|d|a><r|t1|g>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)

      if(orbb.eq.orbp)then
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif

      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.ne.igsym)then
      call m2prop(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 6 = sum(p,g,r)<bg|t2|pr><r|d|g><p|t1|a>

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npb=idra(iip,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(npb,nrg,irgsym,2)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(two*two+one)
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      endif
      enddo
      endif


c     diagram 7 = -sum(p,g,r)<bg|t2|rp><r|d|g><p|t1|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+2)*s6j
      ini=idpa(iip,iia)
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      enddo
      endif
      endif
      endif
      endif
      endif


      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then
c     diagrams 8 = sum(p,r,g) <<g|d|r><b|t1|p><pr|t2|ag>
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.ne.igsym)then
      call m2prop(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)

      inf=idpa(iip,iib)
      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,irgsym,2)

      fact=(-1)**(orbg+orbr+2+2*orbr)/(two*two+one)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif


c     diagrams 9 = -sum(p,r,g) <<g|d|r><b|t1|p><pr|t2|ga>
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+2)*s6j
      inf=idpa(iip,iib)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      enddo
      endif
      endif
      endif
      endif
      endif


c     diagrams 10 = (p,r,g) <g|t1|r><b|d|p><pr|t2|ag>
      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.ne.ibsym)then
      call m2prop(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dipole)
      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(two*orbg+1)
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      endif
      enddo
      endif


c     diagrams 11 =-(p,r,g) <g|t1|r><b|d|p><pr|t2|ga>
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)

      if(orba.eq.orbp)then
c     sum=sum-dipole*tf(inf)*tf(inf)*fact
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif

      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      if(ipsym.eq.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c     diagram 11 = sum <p,r) <b|t1|r><t|d|p><p|t1|a>
      call m2prop(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)
      inf=idpa(iir,iib)
      ini=idpa(iip,iia)
      sum=sum+dipole*tf(inf)*ti(ini)
 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
c     diagram 12 = - sum (p,g) <g|t1|p><b|d|g><p|t1|a>
      if(igsym.ne.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      call m2prop(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)
      ini=idpa(iip,iia)
      inf=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
c     diagram 13 = - sum (p,g) <g|d|a><b|t1|p><p|t1|g>
      if(igsym.ne.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      call m2prop(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iip,iib)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
 6    continue
 4    continue


c     diagrams 15 and 15 (a)
      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)
      if(igsym.eq.idsym)go to 8
c     diagram 15 = - sum (p,r,g,d) <gb|t2|rp><d|d|g><pr|t2|ad>
      call m2prop(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8
      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,two,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,two,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+2+2*orbg)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 10   continue
c     diagram 15(a) = sum (p,r,g,d) <gb|t2|rp><d|d|g><pr|t2|da>
      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,two,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dipole*tf(inf)*ti(ini)*s9j*fact
 11   continue
  9   continue
  8   continue

c     diagrams 16 and 17
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      if(iqsym.eq.irsym)go to 13
      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13
      call m2prop(orbq,orbr,iiq,iir,orbq,orbr,kapq,kapr,dipole)
c     diagram 16 = sum(p,q,r,g) <bg|t2|pq><q|d|r><pr|t2|ag>
      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,two,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,two,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
c     diagram 17 = -sum (p,q,r,g) <bg|t2|pq><q|d|r><pr|t2|ga>
      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(two,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+two+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      enddo
 14   continue
 13   continue
c     diagrams 18 and 19
      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.eq.irsym)go to 16
      call m2prop(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dipole)

c     diagram 18= sum (p,r,q,g) <bg|t2|rq><r|d|p><pq|t2|ag>
      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      call sixj(two,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+two)
      fact=s6j*isign/(2*al1+1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact
      endif
      enddo
      endif

c     diagram 19= -sum (p,r,q,g) <bg|t2|rq><r|d|p><pq|t2|ga>
      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,two,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+two)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dipole*fact
  18  continue
  17  continue
  16  continue
  12  continue
c     diagrams 20,21,22,23
      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)
      if(idsym.eq.ibsym)go to 21
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22
      call m2prop(orbb,orbd,iib,iid,orbb,orbd,kapb,kapd,dipole)

c     diagram 20= -sum (p,r,d,g) <gd|t2|rp><b|d|d><pr|t2|ag>
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif

c     diagram 21=  sum (p,r,d,g) <gd|t2|pr><b|d|d><pr|t2|ag>
      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 24   continue
 23   continue
 22   continue
 21   continue
      if(idsym.eq.iasym)go to 25
      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26
      call m2prop(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)

c     diagram 22= -sum (p,r,d,g) <bg|t2|pr><d|d|a><pr|t2|dg>
      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dipole*fact
      endif
      enddo
      endif

c     diagram 23=  sum (p,r,d,g) <gb|t2|pr><d|d|a><pr|t2|dg>
      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dipole*fact
 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(Orba+Orbb+one)
 1    continue
      return
      end

      subroutine m2prophp(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)


      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo


c     H-P/P-H type <p|d|a>/<a|d|p>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)
      if(ipsym.eq.iasym)go to 1
c     diagram 1  = <p|d|a>
      call m2prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      d(iip,iia)=d(iip,iia)+dipole
      sum=zero


c     Diagrams that are linear in T
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
      if(igsym.eq.irsym)go to 3
      call m2prop(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.irgsym)go to 401
c     diagram 2 = sum(r,g) <pr|t2|ag><g|d|r>

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nrg,ipasym,2)
      fact=(-1)**(orbg+orbr+2+2*orbr)/(2*two+1)
      sum=sum+dipole*ti(ini)*fact
      endif
      enddo
      endif


c     diagram 3 = -sum(r,g)<g|d|r><pr|t2|ga>
      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(two,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum-dipole*ti(ini)*fact
 4    continue
      endif
 401  continue

c     Non-linear contribution
c     diagram 8 = -sum(r,g) <r|t1|a><r|d|g><p|t1|g>
      call m2prop(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3
      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*ti(inj)
 3    continue


c     diagram 9 = -sum(r,g) <r|t1|a><g|t1|r><p|d|g>
      if(ipsym.eq.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5
      call m2prop(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dipole)
      ini=idpa(iir,iia)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 5    continue

c     diagram 10 = -sum(r,g) <r|d|a><g|t|r><p|t|g>
      if(irsym.eq.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6
      call m2prop(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dipole)
      inf=idpa(iir,iig)
      ini=idpa(iip,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
 6    continue


      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7
      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      if(ipsym.eq.iqsym)go to 8
      call m2prop(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      if(irgsym.ne.iqasym)go to 8
      inf=idpa(iir,iig)

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     diagram 4 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ag>
c     sum=sum+dipole*tf(inf)*ti(ini)*dsqrt(2*orbg+1)
      sum=sum+dipole*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))
      endif
      enddo
      endif

c     diagram 5 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ga>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 9    continue
 8    continue


      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      if(iasym.eq.ibsym)go to 10
      call m2prop(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)
      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)
      if(irgsym.ne.ipbsym)go to 10
      inf=idpa(iir,iig)

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbp.eq.orbb)then
c     diagram 6 = -sum(b,r,g) <g|t1|r><b|d|a><pr|t2|bg>
      sum=sum-dipole*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))
      endif
      endif
      enddo
      endif

c     diagram 7 = sum(b,r,g) <g|t1|r><b|d|a><pr|t2|gb>
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
  7   continue

c     diagram 11,12, 13 and 14
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqasym) goto 13
      if(ipsym.eq.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      call m2prop(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dipole)

c     diagram 11 = -sum (q,g,r,b)<gb|t2|rq><qr|t2|ag><p|d|b>
      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 12 = sum (q,g,r,b)<gb|t2|rq><qr|t2|ga><p|d|b>
      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 15   continue
 14   continue
 13   continue


      if(orbp.ne.orbq)go to 12
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.eq.iasym)go to 12
c     diagram 13 = -sum(q,r,g,b) <gb|t2|rq><pr|t2|bg><q|d|a>
      call m2prop(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dipole)

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 14 = sum(q,r,g,b) <gb|t2|rq><pr|t2|gb><q|d|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 17   continue
 16   continue
 12   continue


c     diagram 15 to 18
      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2
      inf=idpa(iir,iig)
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)
      if(irsym.eq.iqsym)go to 20
      call m2prop(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(ipasym.ne.iqgsym)go to 20

      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nqg=idra(iiq,iig,2)
      ini=ntloc(npa,nqg,ipasym,2)
      fact=(-1)**(orbg+orbq+2*orbg+2)/(2*2+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(two,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+2)*s6j
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 21   continue
 20   continue
      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)
      if(ibsym.eq.igsym)go to 22
      call m2prop(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dipole)
      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      if(ipasym.ne.irbsym)go to 22

      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(iloop)
      if(l2.eq.2)then
      npa=idra(iip,iia,2)
      nrb=idra(iir,iib,2)
      ini=ntloc(npa,nrb,ipasym,2)
      fact=(-1)**(orbb+orbr+2*orbb+2)/(2*2+1)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(two,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+2)*s6j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 23   continue
 22   continue
  2   continue


c     diagram 19
      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
      if(iqsym.eq.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24
      ini=idpa(iiq,iia)
      call m2prop(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
      sum=sum+ti(ini)*dipole
 24   continue
c     diagram 20
      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)
      if(ibsym.eq.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25
      ini=idpa(iip,iib)
      call m2prop(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dipole)
      sum=sum-ti(ini)*dipole
 25   continue
c     diagram 21-24
      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      if(igsym.eq.irsym)go to 26
      call m2prop(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)
c     diagram 21
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      npa=idra(iip,iia,2)
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      ini=ntloc(npa,nqb,ipasym,2)
      inf=ntloc(nrg,nqb,iqbsym,2)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*2+1)*(2*2+1))
      sum=sum+ti(ini)*tf(inf)*dipole*fact
      endif
      enddo
      endif
 27   continue

c     diagram 22
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.2)then
      nqb=idra(iiq,iib,2)
      nrg=idra(iir,iig,2)
      inf=ntloc(nrg,nqb,iqbsym,2)
      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(2.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*2+1)
      sum=sum-tf(inf)*ti(ini)*dipole*fact
 29   continue
      endif
      endif
      enddo
      endif
 28   continue
c     diagram 23
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(2.0d0,orbb,orbq,al1,orbg,orbr,s6j)
c     fact=(-1)**(orbb+orbr+orbq+orbg)*s6j/(2*2+1)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*2+1)
      sum=sum-ti(ini)*tf(inf)*dipole*fact
 31   continue
      endif
      endif
      enddo
      endif
 30   continue
c     diagram 24
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32
      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(2.0d0,orbr,orbg,al2,orbq,orbb,s6j1)
      call sixj(2.0d0,orba,orbp,al1,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dipole*fact
 34   continue
 33   continue
 32   continue
 26   continue


      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(Orba+Orbp+one)
 1    continue
      return
      end

      subroutine m2proppp(d,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)
      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(ipsym.eq.iqsym) go to 1
      call m2prop(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
c     diagram 1 =<p|d|q>
      d(iip,iiq)=d(iip,iiq)+dipole
      sum=0.0d0

      sum6=0
      sum8=0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      if(iasym.ne.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then
      call m2prop(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)
      ini=idpa(iip,iia)
c     diagram 2 = -sum (a) <p|t|a><a|d|q>
      sum=sum-dipole*ti(ini)
      endif
      endif
      endif

      if(iasym.ne.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then
      call m2prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
      inf=idpa(iiq,iia)
c     diagram 3 = -sum (a) <p|d|a><a|t|q>
      sum=sum-dipole*tf(inf)
      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.ne.ipsym)then
      ini=idpa(iir,iig)

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
      call m2prop(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dipole)
c     diagram 4 = -sum (a,g,r) <ag|t2|qr><r|t|g><p|d|a>
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      endif
      enddo
      endif

c     diagram 5 =  sum (a,g,r) <ag|t2|rq><r|t|g><p|d|a>
      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 6 and 7
      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.ne.irsym)then
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.eq.iqasym)then
      call m2prop(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dipole)
      ini=idpa(iip,iia)

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      nrg=idra(iir,iig,2)
      nqa=idra(iiq,iia,2)
      inf=ntloc(nrg,nqa,irgsym,2)
      fact=(-1)**(orbg+orbr+2*orbr+2)/(2*two+1)
c     diagram 6 = -sum (a,g,r) <ag|t2|qr><r|d|g><p|t|a>
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 7 =  sum (a,g,r) <ag|t2|rq><r|d|g><p|t|a>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(two,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 8 and 9
      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.ne.irsym)then
      ipasym=mtbl(ipsym,iasym)
      if(irgsym.eq.ipasym)then
      call m2prop(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dipole)
      inf=idpa(iiq,iia)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.2)then
      nrg=idra(iir,iig,2)
      npa=idra(iip,iia,2)
      ini=ntloc(nrg,npa,irgsym,2)
      fact=(-1)**(orbg+orbr+2*orbr+2)/(2*two+1)
c     diagram 8 = -sum (a,g,r) <g|d|r><a|t|q><pr|t2|ag>
      sum=sum-dipole*ti(ini)*tf(inf)*fact
      endif
      enddo
      endif

c     diagram 9 =  sum (a,g,r) <g|d|r><a|t|q><pr|t2|ga>
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(two,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+two)*s6j
      sum=sum+dipole*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 10 and 11

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.ne.iqsym)then
      ipasym=mtbl(ipsym,iasym)
      if(irgsym.eq.ipasym)then
      call m2prop(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dipole)
      inf=idpa(iir,iig)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     diagram 10 = -sum (a,g,r) <g|t|r><a|d|q><pr|t2|ag>
c     sum=sum-dipole*ti(ini)*tf(inf)*dsqrt(2*orbg+1)
      sum=sum-dipole*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      endif
      enddo
      endif

c     diagram 11 = -sum (a,g,r) <g|t|r><a|d|q><pr|t2|ga>
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)sum=sum+dipole*tf(inf)*ti(ini)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
c     diagram 13
      if(irsym.ne.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      call m2prop(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iir,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif

c     diagram 14
      if(irsym.ne.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then
      call m2prop(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,dipole)
      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)
      sum=sum-dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif

  5   continue
c     diagram 12
      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      if(igsym.ne.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then
      call m2prop(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dipole)
      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)
      sum=sum+dipole*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
 6    continue
 4    continue
c     diagram 15 and 18
      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)
      if(igsym.eq.idsym)go to 9
      call m2prop(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dipole)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,two,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,two,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
c     diagram 20
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)
      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,two,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+two)
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 14   continue
 13   continue
 12   continue
  9   continue
      if(iasym.eq.idsym)go to 8
      call m2prop(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dipole)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)
      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(two,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+two+l1)
     & /(2*al1+1)
c     diagram 19
      sum=sum+dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 18
      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+two+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,two,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dipole*ti(ini)*tf(inf)*fact
 16   continue
 15   continue
  8   continue
c     diagram 16
      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nka(iis)
      if(issym.eq.irsym)go to 18
      call m2prop(orbs,orbr,iis,iir,orbs,orbr,kaps,kapr,dipole)
      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18
      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)
      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,two,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,two,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+two+2*orbq)*s6j1*s6j2
      sum=sum-dipole*tf(inf)*ti(ini)*fact
 20   continue
c     diagram 17
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,two,orbg,al2,orbp,s9j)
      fact=(-1)**(orbq+orbg+l2+two+2*orbp)*s9j
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 21   continue
 19   continue
 18   continue
c     diagram 21 and 22
      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)
      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.eq.issym)go to 22
      call m2prop(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,dipole)
      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 24   continue
 23   continue
 22   continue
c     diagram 23 and 24
      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)
      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.eq.issym)go to 25
      call m2prop(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,dipole)
      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dipole*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dipole*tf(inf)*ti(ini)*fact
 27   continue
 26   continue
 25   continue
 17   continue
  7   continue


      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+one)
  1   continue
      return
      end

      subroutine gfacthh(d,dc,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS),dc(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

c     H-H type <b|d|a>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)

      do 1 ib=1,ia
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      orbb=orbc(ib)
      if(ibsym.ne.iasym)go to 1
c     diagram 1  = <b|d|a>
      call gfact(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dip,dipc)
      d(iib,iia)=d(iib,iia)+dip
      dc(iib,iia)=dc(iib,iia)+dipc
      sum=zero
      sumc=zero


c     Diagrams that are linear in T
      do 2 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
c     diagram 2 = sum(p) <b|d|p><p|t1|a>
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      ini=idpa(iip,iia)
      call gfact(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dip,dipc)
      sum=sum+dip*ti(ini)
      sumc=sumc+dipc*ti(ini)
      endif
      endif
      endif
c     diagram 3 = sum(p) <b|t1|p><p|d|a>
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      inf=idpa(iip,iib)
      call gfact(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dip,dipc)
      sum=sum+dip*tf(inf)
      sumc=sumc+dipc*tf(inf)
      endif
      endif
      endif


c     diagrams nonlinear in T

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 3 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)

      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.eq.ipbsym) then

c     diagrams 4 = sum(p,r,g) <bg|t2|pr><p|d|a><r|t1|g>
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.iasym)then

      ini=idpa(iir,iig)
      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if (l1.eq.0)then
      npb=idra(iip,iib,0)
      nrg=idra(iir,iig,0)
      inf=ntloc(npb,nrg,irgsym,0)
      call gfact(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dip,dipc)
      sum=sum+dip*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      sumc=sumc+dipc*ti(ini)*tf(inf)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      endif
      enddo
      endif


c     diagrams 5 = -sum(p,r,g) <bg|t2|rp><p|d|a><r|t1|g>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      fact=(-1)**(orbp+orbg+l1+2*orbg)/(two*orbp+1)

      if(orbb.eq.orbp)then
      sum=sum-dip*ti(ini)*tf(inf)*fact
      sumc=sumc-dipc*ti(ini)*tf(inf)*fact
      endif

      enddo
      endif
      endif
      endif
      endif
      endif

      if(ipbsym .eq. irgsym) then
      if(orbp.eq.orba)then
      if(ipsym.eq.iasym)then
      if(irsym.eq.igsym)then
      call gfact(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dip,dipc)
c     diagram 6 = sum(p,g,r)<bg|t2|pr><r|d|g><p|t1|a>

      call findk(orbg,orbr,orbp,orbb,iag,iar,iap,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npb=idra(iip,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(npb,nrg,irgsym,1)
      ini=idpa(iip,iia)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dip*ti(ini)*tf(inf)*fact
      sumc=sumc+dipc*ti(ini)*tf(inf)*fact
      endif
      enddo
      endif


c     diagram 7 = -sum(p,g,r)<bg|t2|rp><r|d|g><p|t1|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ipgsym=mtbl(ipsym,igsym)
      inf=ntloc(npg,nrb,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orba,orbb,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      ini=idpa(iip,iia)
      sum=sum-dip*ti(ini)*tf(inf)*fact
      sumc=sumc-dipc*ti(ini)*tf(inf)*fact
      enddo
      endif
      endif
      endif
      endif
      endif


      ipasym=mtbl(ipsym,iasym)
      if(ipasym .eq. irgsym) then
c     diagrams 8 = sum(p,r,g) <<g|d|r><b|t1|p><pr|t2|ag>
      if(orbp.eq.orbb)then
      if(ipsym.eq.ibsym)then
      if(irsym.eq.igsym)then
      call gfact(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dip,dipc)

      inf=idpa(iip,iib)
      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,irgsym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(two*one+one)
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif


c     diagrams 9 = -sum(p,r,g) <<g|d|r><b|t1|p><pr|t2|ga>
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbb,orba,s6j)
      fact=(-1)**(orbr+orbg+1)*s6j
      inf=idpa(iip,iib)
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
      enddo
      endif
      endif
      endif
      endif
      endif


c     diagrams 10 = (p,r,g) <g|t1|r><b|d|p><pr|t2|ag>
      if(ipasym .eq. irgsym) then
      if(orbg.eq.orbr)then
      if(irsym.eq.igsym)then
      if(ipsym.eq.ibsym)then
      call gfact(orbb,orbp,iib,iip,orbb,orbp,kapb,kapp,dip,dipc)
      inf=idpa(iir,iig)

      call findk(orbg,orbr,orbp,orba,iag,iar,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      npa=idra(iip,iia,0)
      nrg=idra(iir,iig,0)
      ini=ntloc(npa,nrg,irgsym,0)
      sum=sum+dip*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      sumc=sumc+dipc*tf(inf)*ti(ini)*
     :dsqrt((two*orbg+1)/(two*orbp+1))
      endif
      enddo
      endif


c     diagrams 11 =-(p,r,g) <g|t1|r><b|d|p><pr|t2|ga>
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      inf=idpa(iir,iig)
      fact=(-1)**(orba+orbr+l1+2*orbr)/(two*orbp+1)

      if(orba.eq.orbp)then
      sum=sum-dip*ti(ini)*tf(inf)*fact
      sumc=sumc-dipc*ti(ini)*tf(inf)*fact
      endif

      enddo
      endif

      endif
      endif
      endif
      endif

 3    continue
 2    continue


      do 4 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      do 5 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      if(ipsym.ne.irsym) go to 5
      if(ipsym.ne.iasym)go to 5
      if(orbp.ne.orba)go to 5
      if(irsym.ne.ibsym)go to 5
      if(orbr.ne.orbb)go to 5

c     diagram 11 = sum <p,r) <b|t1|r><t|d|p><p|t1|a>
      call gfact(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dip,dipc)
      inf=idpa(iir,iib)
      ini=idpa(iip,iia)
      sum=sum+dip*tf(inf)*ti(ini)
      sumc=sumc+dipc*tf(inf)*ti(ini)
 5    continue

      do 6 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
c     diagram 12 = - sum (p,g) <g|t1|p><b|d|g><p|t1|a>
      if(igsym.eq.ibsym)then
      if(ipsym.eq.iasym)then
      if(orbp.eq.orba)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      call gfact(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dip,dipc)
      ini=idpa(iip,iia)
      inf=idpa(iip,iig)
      sum=sum-dip*ti(ini)*tf(inf)
      sumc=sumc-dipc*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
c     diagram 13 = - sum (p,g) <g|d|a><b|t1|p><p|t1|g>
      if(igsym.eq.iasym)then
      if(ipsym.eq.ibsym)then
      if(orbp.eq.orbb)then
      if(ipsym.eq.igsym)then
      if(orbp.eq.orbg)then
      call gfact(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dip,dipc)
      ini=idpa(iip,iig)
      inf=idpa(iip,iib)
      sum=sum-dip*ti(ini)*tf(inf)
      sumc=sumc-dipc*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
 6    continue
 4    continue


c     diagrams 15 and 15 (a)
      do 7 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      kapp=nka(iip)
      do 7 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 7 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
      do 8 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)
      if(igsym.ne.idsym)go to 8
c     diagram 15 = - sum (p,r,g,d) <gb|t2|rp><d|d|g><pr|t2|ad>
      call gfact(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dip,dipc)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      if(ipasym.ne.irdsym)go to 8
      if(ipbsym.ne.irgsym)go to 8
      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,igot,imax)
      if(imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 10 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l2)
      nrd=idra(iir,iid,l2)
      ini=ntloc(npa,nrd,ipasym,l2)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(orba+orbb+1+2*orbg)*s6j1*s6j2
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
 10   continue
c     diagram 15(a) = sum (p,r,g,d) <gb|t2|rp><d|d|g><pr|t2|da>
      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,jgot,jmax)
      if(jmax.eq.0)go to 9
      do 11 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npd=idra(iip,iid,l2)
      nra=idra(iir,iia,l2)
      ipdsym=mtbl(ipsym,idsym)
      ini=ntloc(npd,nra,ipdsym,l2)
      call ninej(al1,orbp,orbb,orbg,orbd,one,orbr,al2,orba,s9j)
      fact=(-1)**(l1+l2+orbr+orbg+orbp+orbd)
      sum=sum+dip*tf(inf)*ti(ini)*s9j*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*s9j*fact
 11   continue
  9   continue
  8   continue

c     diagrams 16 and 17
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      if(iqsym.ne.irsym)go to 13
      iqgsym=mtbl(iqsym,igsym)
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipasym)go to 13
      if(iqgsym.ne.ipbsym)go to 13
      call gfact(orbq,orbr,iiq,iir,orbq,orbr,kapq,kapr,dip,dipc)
c     diagram 16 = sum(p,q,r,g) <bg|t2|pq><q|d|r><pr|t2|ag>
      call findk(orbq,orbg,orbp,orbb,iaq,iag,iap,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqg=idra(iiq,iig,l1)
      inf=ntloc(npb,nqg,ipbsym,l1)
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,jgot,jmax)
      if(jmax.eq.0)go to 14
      do 15 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(npa,nrg,ipasym,l2)
      call sixj(orbq,orbr,one,al2,al1,orbg,s6j1)
      call sixj(orbb,orba,one,al2,al1,orbp,s6j2)
      fact=(-1)**(l1+l2+orba+orbb+2*orbq)*s6j1*s6j2
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 15   continue
c     diagram 17 = -sum (p,q,r,g) <bg|t2|pq><q|d|r><pr|t2|ga>
      call findk(orbr,orba,orbp,orbg,iar,iaa,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 14
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l2)
      call sixj(al1,orbg,orbq,al2,orbb,orbp,s6j1)
      call sixj(one,orbq,orbr,al2,orba,orbb,s6j2)
      fact=(-1)**(orbr+orbq+one+orbg+orbp+l2)*s6j1*s6j2
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
      enddo
 14   continue
 13   continue
c     diagrams 18 and 19
      iqgsym=mtbl(iqsym,igsym)
      irbsym=mtbl(irsym,ibsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.iqgsym)go to 16
      if(irbsym.ne.iqgsym)go to 16
      if(ipsym.ne.irsym)go to 16
      call gfact(orbr,orbp,iir,iip,orbr,orbp,kapr,kapp,dip,dipc)

c     diagram 18= sum (p,r,q,g) <bg|t2|rq><r|d|p><pq|t2|ag>
      call findk(orbr,orbb,orbq,orbg,iar,iab,iaq,iag,jgot,jmax)
      if(jmax.eq.0)go to 16
      do 17 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      call sixj(one,orbp,orbr,al1,orbb,orba,s6j)
      isign=(-1)**(orbg+orbq+2*orbg+l1+orbr+orbp+one)
      fact=s6j*isign/(2*al1+1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)

      call findk(orbp,orba,orbq,orbg,iap,iaa,iaq,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      k1=igot(iloop)
      if(k1.eq.l1)then
      nqg=idra(iiq,iig,l1)
      npa=idra(iip,iia,l1)
      inf=ntloc(nqg,nrb,irbsym,l1)
      ini=ntloc(nqg,npa,ipasym,l1)
      sum=sum+tf(inf)*ti(ini)*dip*fact
      sumc=sumc+tf(inf)*ti(ini)*dipc*fact
      endif
      enddo
      endif

c     diagram 19= -sum (p,r,q,g) <bg|t2|rq><r|d|p><pq|t2|ga>
      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nqa=idra(iiq,iia,l2)
      iqasym=mtbl(iqsym,iasym)
      ini=ntloc(npg,nqa,iqasym,l2)
      call sixj(orbq,orba,al2,orbp,orbg,al1,s6j1)
      call sixj(orba,orbb,one,orbr,orbp,al1,s6j2)
      isign=(-1)**(orbg+orbq+al1+orbp+orbr+one)
      fact=isign*s6j1*s6j2
      sum=sum-tf(inf)*ti(ini)*dip*fact
      sumc=sumc-tf(inf)*ti(ini)*dipc*fact
  18  continue
  17  continue
  16  continue
  12  continue
c     diagrams 20,21,22,23
      do 20 id=1,nocc
      idsym=isymc(id)
      iid=kc(id)
      orbd=orbc(id)
      iad=iqc(id)
      kapd=nka(iid)
      if(idsym.ne.ibsym)go to 21
      irgsym=mtbl(irsym,igsym)
      ipasym=mtbl(ipsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipasym)go to 22
      if(irgsym.ne.ipdsym)go to 22
      if(orba.ne.orbd)go to 22
      call gfact(orbb,orbd,iib,iid,orbb,orbd,kapb,kapd,dip,dipc)

c     diagram 20= -sum (p,r,d,g) <gd|t2|rp><b|d|d><pr|t2|ag>
      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npa,nrg,ipasym,l1)

      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      inf=ntloc(npd,nrg,irgsym,l1)
      isign=(-1)**(orbg+orbr+orbd+orbp)
      fact=isign/((2*l1+1)*(2*orba+1))
      sum=sum-ti(ini)*tf(inf)*dip*fact
      sumc=sumc-ti(ini)*tf(inf)*dipc*fact
      endif
      enddo
      endif

c     diagram 21=  sum (p,r,d,g) <gd|t2|pr><b|d|d><pr|t2|ag>
      call findk(orbr,orbd,orbp,orbg,iar,iad,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 23
      do 24 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npa,nrg,ipasym,l1)
      inf=ntloc(npg,nrd,ipgsym,l2)
      call sixj(al2,orbr,orba,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orba+orbp+2*orbr)
      fact=isign*s6j/(2*orba+1)
      sum=sum+ti(ini)*tf(inf)*dip*fact
      sumc=sumc+ti(ini)*tf(inf)*dipc*fact
 24   continue
 23   continue
 22   continue
 21   continue
      if(idsym.ne.iasym)go to 25
      irgsym=mtbl(irsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      ipdsym=mtbl(ipsym,idsym)
      if(irgsym.ne.ipdsym)go to 26
      if(irgsym.ne.ipbsym)go to 26
      if(orbb.ne.orbd)go to 26
      call gfact(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dip,dipc)

c     diagram 22= -sum (p,r,d,g) <bg|t2|pr><d|d|a><pr|t2|dg>
      call findk(orbr,orbg,orbp,orbd,iar,iag,iap,iad,igot,imax)
      if(imax.eq.0)go to 26
      do 27 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)

      call findk(orbr,orbg,orbp,orbb,iar,iag,iap,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npb=idra(iip,iib,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npb,nrg,ipbsym,l1)
      isign=(-1)**(orbg+orbr+orbb+orbp)
      fact=isign/((2*l1+1)*(2*orbb+1))
      sum=sum-tf(inf)*ti(ini)*dip*fact
      sumc=sumc-tf(inf)*ti(ini)*dipc*fact
      endif
      enddo
      endif

c     diagram 23=  sum (p,r,d,g) <gb|t2|pr><d|d|a><pr|t2|dg>
      call findk(orbr,orbb,orbp,orbg,iar,iab,iap,iag,jgot,jmax)
      if(jmax.eq.0)go to 27
      do 28 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nrb=idra(iir,iib,l2)
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npd,nrg,irgsym,l1)
      inf=ntloc(npg,nrb,ipgsym,l2)
      call sixj(al2,orbr,orbb,al1,orbp,orbg,s6j)
      isign=(-1)**(orbg+orbr+orbb+orbp+2*orbb)
      fact=isign*s6j/(2*orbb+1)
      sum=sum+tf(inf)*ti(ini)*dip*fact
      sumc=sumc+tf(inf)*ti(ini)*dipc*fact
 28   continue
 27   continue
 26   continue
 25   continue
 20   continue
  7   continue

      d(iib,iia)=d(iib,iia)+sum
      d(iia,iib)=d(iib,iia)*(-1)**(Orba+Orbb+one)
 1    continue
      return
      end

      subroutine gfacthp(d,dc,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension ti(MDIM),d(MNBAS,MNBAS),dc(MNBAS,MNBAS)
      dimension tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo


c     H-P/P-H type <p|d|a>/<a|d|p>
      do 1 ia=1,nocc
      iasym=isymc(ia)
      iia=kc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      iap=iqe(ip)
      orbp=orbe(ip)
      kapp=nka(iip)
      if(ipsym.ne.iasym)go to 1
c     diagram 1  = <p|d|a>
      call gfact(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dip,dipc)
      d(iip,iia)=d(iip,iia)+dip
      dc(iip,iia)=dc(iip,iia)+dipc
      sum=zero
      sumc=zero


c     Diagrams that are linear in T
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      iir=ke(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 2 ig=1,nocc
      igsym=isymc(ig)
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      kapg=nka(iig)
      if(igsym.ne.irsym)go to 3
      call gfact(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dip,dipc)
      ipasym=mtbl(ipsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.irgsym)go to 401
c     diagram 2 = sum(r,g) <pr|t2|ag><g|d|r>

      call findk(orba,orbp,orbr,orbg,iaa,iap,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nrg,ipasym,1)
      fact=(-1)**(orbg+orbr+1+2*orbr)/(2*one+1)
      sum=sum+dip*ti(ini)*fact
      sumc=sumc+dipc*ti(ini)*fact
      endif
      enddo
      endif


c     diagram 3 = -sum(r,g)<g|d|r><pr|t2|ga>
      call findk(orba,orbr,orbp,orbg,iaa,iar,iap,iag,igot,imax)
      if(imax.ne.0) then
      do 4 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      ipgsym=mtbl(ipsym,igsym)
      ini=ntloc(npg,nra,ipgsym,l1)
      call sixj(one,orbr,orbg,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum-dip*ti(ini)*fact
      sumc=sumc-dipc*ti(ini)*fact
 4    continue
      endif
 401  continue

c     Non-linear contribution
c     diagram 8 = -sum(r,g) <r|t1|a><r|d|g><p|t1|g>
      call gfact(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dip,dipc)
      if(ipsym.ne.igsym)go to 3
      if(orbp.ne.orbg)go to 3
      if(irsym.ne.iasym)go to 3
      if(orbr.ne.orba)go to 3
      ini=idpa(iir,iia)
      inj=idpa(iip,iig)
      sum=sum-dip*ti(ini)*ti(inj)
      sumc=sumc-dipc*ti(ini)*ti(inj)
 3    continue


c     diagram 9 = -sum(r,g) <r|t1|a><g|t1|r><p|d|g>
      if(ipsym.ne.igsym)go to 5
      if(irsym.ne.iasym)go to 5
      if(orbr.ne.orba)go to 5
      if(igsym.ne.irsym)go to 5
      if(orbr.ne.orbg)go to 5
      call gfact(orbp,orbg,iip,iig,orbp,orbg,kapp,kapg,dip,dipc)
      ini=idpa(iir,iia)
      inf=idpa(iir,iig)
      sum=sum-dip*ti(ini)*tf(inf)
      sumc=sumc-dipc*ti(ini)*tf(inf)
 5    continue

c     diagram 10 = -sum(r,g) <r|d|a><g|t|r><p|t|g>
      if(irsym.ne.iasym)go to 6
      if(irsym.ne.igsym)go to 6
      if(orbr.ne.orbg)go to 6
      if(igsym.ne.ipsym)go to 6
      if(orbg.ne.orbp)go to 6
      call gfact(orbr,orba,iir,iia,orbr,orba,kapr,kapa,dip,dipc)
      inf=idpa(iir,iig)
      ini=idpa(iip,iig)
      sum=sum-dip*ti(ini)*tf(inf)
      sumc=sumc-dipc*ti(ini)*tf(inf)
 6    continue


      if(igsym.ne.irsym) go to 7
      if(orbg.ne.orbr) go to 7
      do 8 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      if(ipsym.ne.iqsym)go to 8
      call gfact(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dip,dipc)
      iqasym=mtbl(iqsym,iasym)
      irasym=mtbl(irsym,iasym)
      if(irgsym.ne.iqasym)go to 8
      inf=idpa(iir,iig)

      call findk(orba,orbq,orbr,orbg,iaa,iaq,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      ini=ntloc(nrg,nqa,iqasym,0)
c     diagram 4 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ag>
      sum=sum+dip*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))
      sumc=sumc+dipc*tf(inf)*ti(ini)*
     :dsqrt((2*orbg+1)/(two*orbq+1))
      endif
      enddo
      endif

c     diagram 5 = sum(q,r,g) <g|t1|r><p|d|q><qr|t2|ga>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if (imax.eq.0)go to 8
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      ini=ntloc(nqg,nra,irasym,l1)
      if(orba.ne.orbq)go to 9
      fact=(-1)**(orba+orbr+al1+2*orba)/(2*orba+1)
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
 9    continue
 8    continue


      do 10 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      if(iasym.ne.ibsym)go to 10
      call gfact(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dip,dipc)
      ipbsym=mtbl(ipsym,ibsym)
      ipgsym=mtbl(ipsym,igsym)
      if(irgsym.ne.ipbsym)go to 10
      inf=idpa(iir,iig)

      call findk(orbp,orbb,orbr,orbg,iap,iab,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npb=idra(iip,iib,0)
      ini=ntloc(nrg,npb,irgsym,0)
      if(orbb.eq.orbp)then
c     diagram 6 = -sum(b,r,g) <g|t1|r><b|d|a><pr|t2|bg>
      sum=sum-dip*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))
      sumc=sumc-dipc*tf(inf)*ti(ini)*dsqrt((2*orbg+1)/(2*orbb+1))
      endif
      endif
      enddo
      endif

c     diagram 7 = sum(b,r,g) <g|t1|r><b|d|a><pr|t2|gb>
      call findk(orbp,orbg,orbr,orbb,iap,iag,iar,iab,igot,imax)
      if (imax.eq.0)go to 10
      do 11 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nrb=idra(iir,iib,l1)
      ini=ntloc(npg,nrb,ipgsym,l1)
      if(orbp.ne.orbb)go to 11
      fact=(-1)**(orbp+orbg+al1+2*orbg)/(2*orbb+1)
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
  7   continue

c     diagram 11,12, 13 and 14
      do 12 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      kapq=nka(iiq)
      do 12 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      kapb=nka(iib)
      irgsym=mtbl(irsym,igsym)
      iqbsym=mtbl(iqsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqasym) goto 13
      if(ipsym.ne.ibsym)go to 13
      if(orbb.ne.orba)go to 13
      if(irgsym.ne.iqbsym)go to 13
      call gfact(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,dip,dipc)

c     diagram 11 = -sum (q,g,r,b)<gb|t2|rq><qr|t2|ag><p|d|b>
      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iab,igot,imax)
      if(imax.eq.0)go to 13
      do 14 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,iqbsym,l1)

      call findk(orbr,orbg,orbq,orbb,iar,iag,iaq,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(nrg,nqa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbq+orba)/((2*l1+1)*(2*orba+1))
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 12 = sum (q,g,r,b)<gb|t2|rq><qr|t2|ga><p|d|b>
      call findk(orbr,orba,orbq,orbg,iar,iaa,iaq,iag,kgot,kmax)
      if(kmax.eq.0)go to 14
      do 15 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqg=idra(iiq,iig,l2)
      iqgsym=mtbl(iqsym,igsym)
c     inf=ntloc(nra,nqg,iqgsym,l2)
      ini=ntloc(nra,nqg,iqgsym,l2)
      call sixj(al1,orbq,orba,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbq+orbr+orba+2*orba)*s6j/(2*orba+1)
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 15   continue
 14   continue
 13   continue


      if(orbp.ne.orbq)go to 12
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(irsym,igsym)
      ipgsym=mtbl(ipsym,igsym)
      ipbsym=mtbl(ipsym,ibsym)
      if(irgsym.ne.ipbsym) goto 12
      if(irgsym.ne.iqbsym)go to 12
      if(iqsym.ne.iasym)go to 12
c     diagram 13 = -sum(q,r,g,b) <gb|t2|rq><pr|t2|bg><q|d|a>
      call gfact(orbq,orba,iiq,iia,orbq,orba,kapq,kapa,dip,dipc)

      call findk(orbg,orbr,orbb,orbq,iag,iar,iab,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 16 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqb=idra(iiq,iib,l1)
      inf=ntloc(nrg,nqb,irgsym,l1)

      call findk(orbg,orbr,orbb,orbp,iag,iar,iab,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nrg,npb,irgsym,l1)
      fact=(-1)**(orbg+orbr+orbb+orbp)/((2*l1+1)*(2*orbp+1))
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 14 = sum(q,r,g,b) <gb|t2|rq><pr|t2|gb><q|d|a>
      call findk(orbb,orbr,orbg,orbp,iab,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 16
      do 17 jloop=1,kmax
      l2=kgot(jloop)
      al2=float(l2)
      nrb=idra(iir,iib,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrb,npg,ipgsym,l2)
      call sixj(al1,orbp,orbb,al2,orbr,orbg,s6j)
      fact=(-1)**(orbg+orbr+orbp+orbb+2*orbp)*s6j/(2*orbp+1)
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 17   continue
 16   continue
 12   continue


c     diagram 15 to 18
      if(orbg.ne.orbr)go to 2
      if(igsym.ne.irsym)go to 2
      inf=idpa(iir,iig)
      do 20 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      iaq=iqe(iq)
      orbq=orbe(iq)
      kapq=nka(iiq)
      if(irsym.ne.iqsym)go to 20
      call gfact(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dip,dipc)
      ipasym=mtbl(ipsym,iasym)
      iqgsym=mtbl(iqsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(ipasym.ne.iqgsym)go to 20

      call findk(orbq,orbg,orbp,orba,iaq,iag,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqg=idra(iiq,iig,1)
      ini=ntloc(npa,nqg,ipasym,1)
      fact=(-1)**(orbg+orbq+2*orbg+1)/(2*1+1)
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbq,orba,orbp,orbg,iaq,iaa,iap,iag,igot,imax)
      if(imax.eq.0)go to 20
      do 21 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqa=idra(iiq,iia,l1)
      npg=idra(iip,iig,l1)
      ini=ntloc(nqa,npg,iqasym,l1)
      call sixj(one,orbq,orbr,al1,orbp,orba,s6j)
      fact=(-1)**(orbq+orbr+1)*s6j
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
 21   continue
 20   continue
      do 22 ib=1,nocc
      ibsym=isymc(ib)
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      kapb=nka(iib)
      if(ibsym.ne.igsym)go to 22
      call gfact(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,dip,dipc)
      ipasym=mtbl(ipsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      irasym=mtbl(irsym,iasym)
      if(ipasym.ne.irbsym)go to 22

      call findk(orbr,orbb,orba,orbp,iar,iab,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(iloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nrb=idra(iir,iib,1)
      ini=ntloc(npa,nrb,ipasym,1)
      fact=(-1)**(orbb+orbr+2*orbb+1)/(2*1+1)
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orbr,orba,orbb,orbp,iar,iaa,iab,iap,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      npb=idra(iip,iib,l1)
      ini=ntloc(nra,npb,irasym,l1)
      call sixj(one,orbr,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbb+orbr+1)*s6j
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 23   continue
 22   continue
  2   continue


c     diagram 19
      do 24 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
      if(iqsym.ne.ipsym)go to 24
      if(iqsym.ne.iasym)go to 24
      if(orbq.ne.orba)go to 24
      ini=idpa(iiq,iia)
      call gfact(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dip,dipc)
      sum=sum+ti(ini)*dip
      sumc=sumc+ti(ini)*dipc
 24   continue
c     diagram 20
      do 25 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      kapb=nka(iib)
      if(ibsym.ne.iasym)go to 25
      if(ibsym.ne.ipsym)go to 25
      if(orbp.ne.orbb)go to 25
      ini=idpa(iip,iib)
      call gfact(orbb,orba,iib,iia,orbb,orba,kapb,kapa,dip,dipc)
      sum=sum-ti(ini)*dip
      sumc=sumc-ti(ini)*dipc
 25   continue
c     diagram 21-24
      do 26 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      do 26 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 26 ib=1,nocc
      ibsym=isymc(ib)
      orbb=orbc(ib)
      iab=iqc(ib)
      iib=kc(ib)
      do 26 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      if(igsym.ne.irsym)go to 26
      call gfact(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dip,dipc)
c     diagram 21
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irgsym=mtbl(irsym,igsym)
      if(ipasym.ne.iqbsym)go to 27
      if(irgsym.ne.iqbsym)go to 27

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      ini=ntloc(npa,nqb,ipasym,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      fact=(-1)**(orbg+orbb+orbr+orbq)/((2*1+1)*(2*1+1))
      sum=sum+ti(ini)*tf(inf)*dip*fact
      sumc=sumc+ti(ini)*tf(inf)*dipc*fact
      endif
      enddo
      endif
 27   continue

c     diagram 22
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.ne.iqbsym)go to 28
      if(ipbsym.ne.iqasym)go to 28

      call findk(orbq,orbb,orbr,orbg,iaq,iab,iar,iag,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      nqb=idra(iiq,iib,1)
      nrg=idra(iir,iig,1)
      inf=ntloc(nrg,nqb,iqbsym,1)
      call findk(orbp,orbb,orbq,orba,iap,iab,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do 29 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npb=idra(iip,iib,l1)
      nqa=idra(iiq,iia,l1)
      ini=ntloc(npb,nqa,iqasym,l1)
      call sixj(1.0d0,orbq,orbb,al1,orbp,orba,s6j)
      fact=(-1)**(orbg+orbr+2*orbg+orbq+orbb)*s6j/(2*1+1)
      sum=sum-tf(inf)*ti(ini)*dip*fact
      sumc=sumc-tf(inf)*ti(ini)*dipc*fact
 29   continue
      endif
      endif
      enddo
      endif
 28   continue
c     diagram 23
      ipasym=mtbl(ipsym,iasym)
      iqbsym=mtbl(iqsym,ibsym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 30
      if(ipasym.ne.iqbsym)go to 30

      call findk(orbq,orbb,orbp,orba,iaq,iab,iap,iaa,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      if(l2.eq.1)then
      npa=idra(iip,iia,1)
      nqb=idra(iiq,iib,1)
      ini=ntloc(npa,nqb,iqbsym,1)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.ne.0)then
      do 31 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(1.0d0,orbb,orbq,al1,orbg,orbr,s6j)
      fact=(-1)**(2*orbb+orbb+orbr+orbq+orbg)*s6j/(2*1+1)
      sum=sum-ti(ini)*tf(inf)*dip*fact
      sumc=sumc-ti(ini)*tf(inf)*dipc*fact
 31   continue
      endif
      endif
      enddo
      endif
 30   continue
c     diagram 24
      ipbsym=mtbl(ipsym,ibsym)
      iqasym=mtbl(iqsym,iasym)
      irbsym=mtbl(irsym,ibsym)
      iqgsym=mtbl(iqsym,igsym)
      if(iqgsym.ne.irbsym)go to 32
      if(ipbsym.ne.iqasym)go to 32
      call findk(orbq,orba,orbp,orbb,iaq,iaa,iap,iab,kgot,kmax)
      if(kmax.eq.0)go to 32
      do 33 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npb=idra(iip,iib,l2)
      nqa=idra(iiq,iia,l2)
      ini=ntloc(npb,nqa,ipbsym,l2)
      call findk(orbq,orbg,orbr,orbb,iaq,iag,iar,iab,igot,imax)
      if(imax.eq.0)go to 33
      do 34 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nrb=idra(iir,iib,l1)
      inf=ntloc(nqg,nrb,iqgsym,l1)
      call sixj(1.0d0,orbr,orbg,al2,orbq,orbb,s6j1)
      call sixj(1.0d0,orba,orbp,al1,orbb,orbq,s6j2)
      fact=(-1)**(orbb+orbr+orbq+orbg)*s6j1*s6j2
      sum=sum+ti(ini)*tf(inf)*dip*fact
      sumc=sumc+ti(ini)*tf(inf)*dipc*fact
 34   continue
 33   continue
 32   continue
 26   continue


      d(iip,iia)=d(iip,iia)+sum
      d(iia,iip)=d(iip,iia)*(-1)**(Orba+Orbp+one)
 1    continue
      return
      end

      subroutine gfactpp(d,dc,ti,tf)
      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS),dc(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)
      do 1 iq=1,ip
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(ipsym.ne.iqsym) go to 1
      call gfact(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dip,dipc)
c     diagram 1 =<p|d|q>
      d(iip,iiq)=d(iip,iiq)+dip
      dc(iip,iiq)=dc(iip,iiq)+dipc
      sum=0.0d0
      sumc=0.0d0

      sum6=0
      sum8=0

      do 2 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      iaa=iqc(ia)
      kapa=nka(iia)
      if(iasym.eq.iqsym)then
      if(iasym.eq.ipsym)then
      if(orbp.eq.orba)then
      call gfact(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dip,dipc)
      ini=idpa(iip,iia)
c     diagram 2 = -sum (a) <p|t|a><a|d|q>
      sum=sum-dip*ti(ini)
      sumc=sumc-dipc*ti(ini)
      endif
      endif
      endif

      if(iasym.eq.ipsym)then
      if(iasym.eq.iqsym)then
      if(orbq.eq.orba)then
      call gfact(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dip,dipc)
      inf=idpa(iiq,iia)
c     diagram 3 = -sum (a) <p|d|a><a|t|q>
      sum=sum-dip*tf(inf)
      sumc=sumc-dipc*tf(inf)
      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.eq.iqasym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(iasym.eq.ipsym)then
      ini=idpa(iir,iig)

      call findk(orbg,orbr,orbq,orba,iag,iar,iaq,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      nqa=idra(iiq,iia,0)
      inf=ntloc(nrg,nqa,irgsym,0)
      call gfact(orbp,orba,iip,iia,orbp,orba,kapp,kapa,dip,dipc)
c     diagram 4 = -sum (a,g,r) <ag|t2|qr><r|t|g><p|d|a>
      sum=sum-dip*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      sumc=sumc-dipc*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      endif
      enddo
      endif

c     diagram 5 =  sum (a,g,r) <ag|t2|rq><r|t|g><p|d|a>
      call findk(orbg,orbq,orbr,orba,iag,iaq,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      fact=(-1)**(orbq+orbg+al1+2*orbp)/(2*orbq+1)
      if(orbq.eq.orba)sum=sum+dip*ti(ini)*tf(inf)*fact
      if(orbq.eq.orba)sumc=sumc+dipc*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 6 and 7
      if(orba.eq.orbp)then
      if(iasym.eq.ipsym)then
      if(igsym.eq.irsym)then
      iqasym=mtbl(iqsym,iasym)
      if(irgsym.eq.iqasym)then
      call gfact(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,dip,dipc)
      ini=idpa(iip,iia)

      call findk(orbq,orba,orbr,orbg,iaq,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      nqa=idra(iiq,iia,1)
      inf=ntloc(nrg,nqa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
c     diagram 6 = -sum (a,g,r) <ag|t2|qr><r|d|g><p|t|a>
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 7 =  sum (a,g,r) <ag|t2|rq><r|d|g><p|t|a>
      call findk(orbq,orbg,orbr,orba,iaq,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      inf=ntloc(nqg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dip*ti(ini)*tf(inf)*fact
      sumc=sumc+dipc*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 8 and 9
      if(orba.eq.orbq)then
      if(iasym.eq.iqsym)then
      if(igsym.eq.irsym)then
      ipasym=mtbl(ipsym,iasym)
      if(irgsym.eq.ipasym)then
      call gfact(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,dip,dipc)
      inf=idpa(iiq,iia)

      call findk(orbp,orba,orbr,orbg,iap,iaa,iar,iag,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.1)then
      nrg=idra(iir,iig,1)
      npa=idra(iip,iia,1)
      ini=ntloc(nrg,npa,irgsym,1)
      fact=(-1)**(orbg+orbr+2*orbr+1)/(2*one+1)
c     diagram 8 = -sum (a,g,r) <g|d|r><a|t|q><pr|t2|ag>
      sum=sum-dip*ti(ini)*tf(inf)*fact
      sumc=sumc-dipc*ti(ini)*tf(inf)*fact
      endif
      enddo
      endif

c     diagram 9 =  sum (a,g,r) <g|d|r><a|t|q><pr|t2|ga>
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      call sixj(one,orbg,orbr,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+one)*s6j
      sum=sum+dip*ti(ini)*tf(inf)*fact
      sumc=sumc+dipc*ti(ini)*tf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 10 and 11

      if(orbr.eq.orbg)then
      if(irsym.eq.igsym)then
      if(iasym.eq.iqsym)then
      ipasym=mtbl(ipsym,iasym)
      if(irgsym.eq.ipasym)then
      call gfact(orba,orbq,iia,iiq,orba,orbq,kapa,kapq,dip,dipc)
      inf=idpa(iir,iig)

      call findk(orbr,orbg,orbp,orba,iar,iag,iap,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nrg=idra(iir,iig,0)
      npa=idra(iip,iia,0)
      ini=ntloc(nrg,npa,irgsym,0)
c     diagram 10 = -sum (a,g,r) <g|t|r><a|d|q><pr|t2|ag>
      sum=sum-dip*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      sumc=sumc-dipc*ti(ini)*tf(inf)*
     :dsqrt((2*orbg+1)/(2*orba+1))
      endif
      enddo
      endif

c     diagram 11 = -sum (a,g,r) <g|t|r><a|d|q><pr|t2|ga>
      call findk(orbp,orbg,orbr,orba,iap,iag,iar,iaa,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nra=idra(iir,iia,l1)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l1)
      fact=(-1)**(orbp+orbg+2*orba+al1)/(2*orba+1)
      if(orba.eq.orbp)sum=sum+dip*tf(inf)*ti(ini)*fact
      if(orba.eq.orbp)sumc=sum+dipc*tf(inf)*ti(ini)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

  3   continue
  2   continue

      do 4 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 5 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
c     diagram 13
      if(irsym.eq.iqsym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      call gfact(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,dip,dipc)
      ini=idpa(iip,iig)
      inf=idpa(iir,iig)
      sum=sum-dip*ti(ini)*tf(inf)
      sumc=sumc-dipc*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif

c     diagram 14
      if(irsym.eq.ipsym)then
      if(igsym.eq.irsym)then
      if(orbg.eq.orbr)then
      if(igsym.eq.iqsym)then
      if(orbg.eq.orbq)then
      call gfact(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,dip,dipc)
      ini=idpa(iir,iig)
      inf=idpa(iiq,iig)
      sum=sum-dip*ti(ini)*tf(inf)
      sumc=sumc-dipc*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif

  5   continue
c     diagram 12
      do 6 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      if(igsym.eq.iasym)then
      if(igsym.eq.ipsym)then
      if(orbg.eq.orbp)then
      if(iasym.eq.iqsym)then
      if(orba.eq.orbq)then
      call gfact(orbg,orba,iig,iia,orbg,orba,kapg,kapa,dip,dipc)
      ini=idpa(iip,iig)
      inf=idpa(iiq,iia)
      sum=sum+dip*ti(ini)*tf(inf)
      sumc=sumc+dipc*ti(ini)*tf(inf)
      endif
      endif
      endif
      endif
      endif
 6    continue
 4    continue
c     diagram 15 and 18
      do 7 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      do 7 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      do 7 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iaa=iqc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      do 8 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)
      if(igsym.ne.idsym)go to 9
      call gfact(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,dip,dipc)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipasym=mtbl(ipsym,iasym)
      if(ipasym.ne.irdsym)go to 9
      if(iqasym.ne.irgsym)go to 9
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 9
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      call findk(orbd,orbr,orba,orbp,iad,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 10
      do 11 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrd,npa,irdsym,l2)
      isign=(-1)**(orbp+orbq+l1+l2+2*orbg)
      call sixj(orbg,orbd,one,al2,al1,orbr,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 11   continue
 10   continue
c     diagram 20
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      ipdsym=mtbl(ipsym,idsym)
      irasym=mtbl(irsym,iasym)
      if(irgsym.ne.iqasym)go to 12
      if(ipdsym.ne.irasym)go to 12
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 12
      do 13 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)
      call findk(orba,orbr,orbd,orbp,iaa,iar,iad,iap,kgot,kmax)
      if(kmax.eq.0)go to 13
      do 14 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npd=idra(iip,iid,l2)
      ini=ntloc(nra,npd,irasym,l2)
      call sixj(orba,orbr,al2,orbg,orbq,al1,s6j1)
      call sixj(orbp,al2,orbd,orbg,one,orbq,s6j2)
      fact=s6j1*s6j2*(-1)**(orba+orbr+l2+orbg+orbd+one)
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
 14   continue
 13   continue
 12   continue
  9   continue
      if(iasym.ne.idsym)go to 8
      call gfact(orbd,orba,iid,iia,orbd,orba,kapd,kapa,dip,dipc)
      irgsym=mtbl(irsym,igsym)
      iqasym=mtbl(iqsym,iasym)
      irdsym=mtbl(irsym,idsym)
      ipgsym=mtbl(ipsym,igsym)
      if(ipgsym.ne.irdsym)go to 8
      if(iqasym.ne.irgsym)go to 8
      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 8
      do 15 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orbp,orbd,iag,iar,iap,iad,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      npd=idra(iip,iid,l1)
      ini=ntloc(nrg,npd,irgsym,l1)
      call sixj(one,orba,orbd,al1,orbp,orbq,s6j)
      fact=s6j*(-1)**(orbg+orbr+2*orbr+orba+orbd+one+l1)
     & /(2*al1+1)
c     diagram 19
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

c     diagram 18
      call findk(orbd,orbr,orbg,orbp,iad,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 15
      do 16 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrd=idra(iir,iid,l2)
      npg=idra(iip,iig,l2)
      ini=ntloc(nrd,npg,irdsym,l2)
      isign=(-1)**(orba+orbd+one+orbg+orbr+l1)
      call sixj(orbg,al1,orbr,orbd,al2,orbp,s6j1)
      call sixj(orbq,al1,orba,orbd,one,orbp,s6j2)
      fact=isign*s6j1*s6j2
      sum=sum-dip*ti(ini)*tf(inf)*fact
      sumc=sumc-dipc*ti(ini)*tf(inf)*fact
 16   continue
 15   continue
  8   continue
c     diagram 16
      do 17 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)
      kaps=nka(iis)
      if(issym.ne.irsym)go to 18
      call gfact(orbs,orbr,iis,iir,orbs,orbr,kaps,kapr,dip,dipc)
      isgsym=mtbl(issym,igsym)
      ipasym=mtbl(ipsym,iasym)
      iqasym=mtbl(iqsym,iasym)
      irgsym=mtbl(irsym,igsym)
      if(irgsym.ne.ipasym)go to 18
      if(isgsym.ne.iqasym)go to 18
      call findk(orbg,orbs,orba,orbq,iag,ias,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 18
      do 19 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqa=idra(iiq,iia,l1)
      inf=ntloc(nsg,nqa,isgsym,l1)
      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 20 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nrg=idra(iir,iig,l2)
      npa=idra(iip,iia,l2)
      ini=ntloc(nrg,npa,irgsym,l2)
      call sixj(orbs,one,orbr,al2,orbg,al1,s6j1)
      call sixj(orbq,orbp,one,al2,al1,orba,s6j2)
      fact=(-1)**(orbp+orbq+one+2*orbq)*s6j1*s6j2
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
 20   continue
c     diagram 17
      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 19
      do 21 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      npg=idra(iip,iig,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nra,npg,irasym,l2)
      call ninej(al1,orba,orbq,orbs,orbr,one,orbg,al2,orbp,s9j)
      fact=(-1)**(orbq+orbg+l2+one+2*orbp)*s9j
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 21   continue
 19   continue
 18   continue
c     diagram 21 and 22
      irgsym=mtbl(irsym,igsym)
      iassym=mtbl(iasym,issym)
      iqasym=mtbl(iqsym,iasym)
      if (irgsym.ne.iassym) goto 22
      if (irgsym.ne.iqasym) goto 22
      if(ipsym.ne.issym)go to 22
      call gfact(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,dip,dipc)
      if(orbs.ne.orbq)go to 22

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,igot,imax)
      if(imax.eq.0)go to 22
      do 23 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqa=idra(iiq,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nqa,irgsym,l1)

      call findk(orbg,orbr,orba,orbq,iag,iar,iaa,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,nsa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbq)/((2*orbq+1)*(2*l1+1))
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orba,orbr,orbg,orbs,iaa,iar,iag,ias,kgot,kmax)
      if(kmax.eq.0)go to 23
      do 24 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      nsg=idra(iis,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(nsg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbq,orba,s6j)
      fact=(-1)**(orba+orbq+orbr+orbg+2*orbq)*s6j/(2*orbq+1)
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 24   continue
 23   continue
 22   continue
c     diagram 23 and 24
      irgsym=mtbl(irsym,igsym)
      iapsym=mtbl(iasym,ipsym)
      iassym=mtbl(iasym,issym)
      if(irgsym .ne. iassym) goto 25
      if(irgsym .ne. iapsym) goto 25
      if(iqsym.ne.issym)go to 25
      call gfact(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,dip,dipc)
      if(orbs.ne.orbp)go to 25

      call findk(orbg,orbr,orba,orbs,iag,iar,iaa,ias,igot,imax)
      if(imax.eq.0)go to 25
      do 26 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsa=idra(iis,iia,l1)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nrg,nsa,irgsym,l1)

      call findk(orbg,orbr,orba,orbp,iag,iar,iaa,iap,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      npa=idra(iip,iia,l1)
      nrg=idra(iir,iig,l1)
      ini=ntloc(nrg,npa,irgsym,l1)
      fact=(-1)**(orbg+orbr+orba+orbs)/((2*orbp+1)*(2*l1+1))
      sum=sum-dip*tf(inf)*ti(ini)*fact
      sumc=sumc-dipc*tf(inf)*ti(ini)*fact
      endif
      enddo
      endif

      call findk(orba,orbr,orbg,orbp,iaa,iar,iag,iap,kgot,kmax)
      if(kmax.eq.0)go to 26
      do 27 kloop=1,kmax
      l2=kgot(kloop)
      al2=float(l2)
      npg=idra(iip,iig,l2)
      nra=idra(iir,iia,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(npg,nra,irasym,l2)
      call sixj(al1,orbr,orbg,al2,orbp,orba,s6j)
      fact=(-1)**(orba+orbs+orbr+orbg+2*orbg)*s6j/(2*orbp+1)
      sum=sum+dip*tf(inf)*ti(ini)*fact
      sumc=sumc+dipc*tf(inf)*ti(ini)*fact
 27   continue
 26   continue
 25   continue
 17   continue
  7   continue


      d(iip,iiq)=d(iip,iiq)+sum
      d(iiq,iip)=d(iip,iiq)*(-1)**(orbp+orbq+one)
  1   continue
      return
      end

c************************************************************************


c************************************************************************
c                                                                       *
c     This subroutine calculates the reduced  matrix element of the     *
c     E2 transition operator                                            *
c                                                                       *
c************************************************************************


      subroutine quadpul(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp)


      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
      common/intquad/rint(MNBAS,MNBAS),kapp1,kapp2
c      common/intquad/rint(MNBAS,MNBAS),orb1,orb2

      orb1=j1
      orb2=j2

      kapp1=kap1
      kapp2=kap2

      factor=(-1)**(j1+0.5)
      w3j1=dr(j1,two,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      hyp=result*e1

      return
      end

c*********************************************************************

      subroutine e2vel(i1,i2,j1,j2,kap1,kap2,hyp)
      implicit double precision(a-h,o-z)
      real*8 j1,j2
      common/cons/zero,half,tenth,one,two,three,ten
      common/intquadvel/rint(MNBAS,MNBAS),kapp1,kapp2

      kapp1=kap1
      kapp2=kap2
      factor=(-1)**(j1+0.5)
      w3j1=dr(j1,two,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      hyp=result*e1

      return
      end

