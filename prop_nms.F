c*************************************************************************
c  Property code written by R K Choudhury and B K Sahoo, IIA,             *
c  Bangalore. Code is modified by Narendra Nath Dutta, IIT, KGP in the    *
c  calcualtions of E1 transition(len), E2 transition                      *
c  (len), M1 transition, hyperfine A and B constants.                     *
c                                                                         *
c  First 46 diagrams are checked and considered.                          *
c                                                                         *
c  PNC part is incomplete                                                 *
c                                                                         *
c               Last modified on 05-02-2013                               *
c**************************************************************************

c**************************************************************************
c     Sample input.prop for neutral K                                     *
c                                                                         *
c     7  14 12 12 11 11 10 10                                             *
c        3 2 2 0 0 0 0                                                    *
c        1 0 1                                                            *
c        3 1 39 19                                                        *
c        4 17 0 0 0 0.5 0.5                                               *
c        1 2 1 1 1 1 1 1 1 1 1                                            *
c     so4s                                                                *
c     so4p1                                                               *
c                                                                         *
c**************************************************************************


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC                                                                  CC
CC               IALL = 0   for only couloumb interaction           CC
CC                    = 1   for all parity channel                  CC
CC                    = 3   for only breit interaction              CC
CC                                                                  CC
CC                                                                  CC
CC               IOPT = 0   for CCD                                 CC
CC                    = 1   for CCSD                                CC
CC                                                                  CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC                                                                  CC
CC    INPUT FILE NAME:   input.prop                                 CC
CC                                                                  CC
CC     READ:                                                        CC
CC    1. a) No. of symm.  b)No.of basis for each symm.              CC
CC    2. a) No. of occupied orbitals for each symm.                 CC
CC    3. a) iopt  b) nprint  c) iall                                CC
CC    4. a) iprop b) initial orbital state index c) final orbital   CC
CC              state index d) g_I  e) j-value f) kappa value       CC
CC    5. a) Open-shell for initial state                            CC
CC    6. a) Open-shell for final   state                            CC
CC                                                                  CC
CC                                                                  CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c  cpp
c
c dimension variables
#define     MN       750
#define     MNBAS    150
#define     MNOCC    100
#define     MNEXC    100
c #define     MNEXC    80
#define     MNSYM    11
#define     MNS      6
#define     MXV      11
#define     MDIM     2800000
c#define    N2ORB    MNBAS*(MNBAS+1)/2 
c#define    N2BAS    MNBAS*MNBAS
c#define    N2ORB    2701
c#define    N2BAS    5329


c constants
#define    STDIN      5
#define    STDOUT     6
#define    WFNIN      7
#define    OCTC       16
#define    DIPOLE     17
#define    DIPOLEV    18
#define    HYPER      19
#define    QUADPULE   20
#define    QUADPULEV  21
#define    QUADB      22
#define    M1PROP     23
#define    M2PROP     24
#define    GFACTOR    25
#define    GFACTORC   26
#define    NIFILE     27
#define    NFFILE     28
#define    OVERLAP    29
#define    NTFILE     30

c*******************************************************************
c     this programme calculates the value of the diagrams
c*******************************************************************
      
      implicit real*8(a-h,o-z)
 
      real*8 jinit,jfin
      dimension nbas(MNSYM),iocc(MNSYM)
      dimension si(MDIM),d(MNBAS,MNBAS),sf(MDIM),dc(MNBAS,MNBAS)
      dimension ti(MDIM),tf(MDIM),sj(MDIM)

      call daopen
      call symmcal
      call setgrd
      call readinp(nsym,nbas,iocc,iopt,nprint,idip,ivel,iquad,ihyp,
     :iquadB,ioctc,im1p,im2p,igfact,ival,jval,gi,qm,ocm,jinit,jfin,
     :insi,insd,inms)
      print*,'inms =',inms
      call setup(nsym,nbasis,nbas,iocc)
      call flopn(ival,jval)
      call grasprd(nbasis)
      call setqic

c--------------------------------------------------------------
c     T1 and T2 cluster operators are read
c--------------------------------------------------------------       

      rewind(NTFILE)
      read(NTFILE)ntmax,(ti(i),i=1,ntmax)

      call symm(iopt,nprint,ti,tf)

c-------------------------------------------------------------------
c     Reading initial so files containting 
c     initial state information
c-------------------------------------------------------------------

      call sloci(iopt,ival,ns1max,nprint)

      rewind(NIFILE)
      read(NIFILE)ns1max,(si(i),i=1,ns1max)

c--------------------------------------------------------------
c     Set the conditions such that if initial and final
c     states are same than no need to calculate E1, M1, M2
c     matrix elements
c--------------------------------------------------------------

      if(idip.eq.1)goto 1111
c      if(iquad.eq.1)goto 1111
      if(im1p.eq.1)goto 1111
      if(im2p.eq.1)goto 1111
      if(ival.ne.jval)goto 1111

      do i=1,ns1max
      sj(i)=si(i)
      enddo

      goto 1112

 1111 continue

c-------------------------------------------------------------------
c     Reading final so files containting
c     final state information
c-------------------------------------------------------------------
     
      rewind(NFFILE)
      read(NFFILE)ns2max,(sj(i),i=1,ns2max)

 1112 continue

      call slocf(iopt,jval,nsmax,nprint,sj,sf)

c---------------------------------------------------------------
c             T+ T     
c     {1+Sf+}e  e{1+Si}  diagrams are calculated
c---------------------------------------------------------------
      
      rewind(OVERLAP)
      read(OVERLAP)((d(i,j),i=1,nbasis),j=1,nbasis)
      call cnormi(d,si,ival)
      call cnormf(d,sf,jval)

c---------------------------------------------------------------
c      O = Electric dipole transition operator (Len Gg)                      
c
c             T+   T   
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------

      if(idip.eq.1)then
      write(*,*)'------------------------------------------------'
      write(*,9989)

      if(ivel.eq.0) goto 1
      if(ivel.eq.2) goto 1
          goto 2

   1  continue

      call cinit(d)
      call calce1mat(nbasis)
      open(DIPOLE,file='dipole.dat',form='unformatted',status='unknown')
      rewind(DIPOLE)
      read(DIPOLE)((d(i,j),i=1,nbasis),j=1,nbasis)
      call dipoleprop(d,si,sf,ti,tf,ival,jval)
      call caldip

   2   continue

c---------------------------------------------------------------
c      O = Electric dipole transition operator (Vel Gg)
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(ivel.eq.1) goto 3
      if(ivel.eq.2) goto 3
          goto 4

   3  continue

      call cinit(d)
      call calce1vel(nbasis)
      open(DIPOLEV,file='dipolev.dat',form='unformatted',
     :status='unknown')
      rewind(DIPOLEV)
      read(DIPOLEV)((d(i,j),i=1,nbasis),j=1,nbasis)
      call dipvelprop(d,si,sf,ti,tf,ival,jval)
      call caldipvel

   4  continue

      endif

c---------------------------------------------------------------
c      O = Electric quadrupole transition operator (Len Gg)
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------

      if(iquad.eq.1)then
      write(*,*)'**********************************************'
      write(*,9990)

      if(ivel.eq.0) goto 5
      if(ivel.eq.2) goto 5
          goto 6

   5  continue

      call cinit(d)
      call calcquad(nbasis)
      open(QUADPULE,file='quadpule.dat',form='unformatted',
     :status='unknown')
      rewind(QUADPULE)
      read(QUADPULE)((d(i,j),i=1,nbasis),j=1,nbasis)
      call quadpuleprop(d,si,sf,ti,tf,ival,jval)
      call calquad

   6  continue

c---------------------------------------------------------------
c      O = Electric quadrupole transition operator (Vel Gg)
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(ivel.eq.1) goto 7
      if(ivel.eq.2) goto 7
          goto 8

   7  continue

      call cinit(d)
      call calcquadvel(nbasis)
      open(QUADPULEV,file='quadpulev.dat',form='unformatted',
     :status='unknown')
      rewind(QUADPULEV)
      read(QUADPULEV)((d(i,j),i=1,nbasis),j=1,nbasis)
      call quadvelprop(d,si,sf,ti,tf,ival,jval)
      call calquadvel

  8   continue

      endif

c---------------------------------------------------------------
c      O = Magnetic dipole hyperfine operator 
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(ihyp.eq.1)then
      write(*,*)'============================================='
      write(*,9991)
      call cinit(d)
      open(HYPER,file='hyper.dat',form='unformatted',status='unknown')
      rewind(HYPER)
      read(HYPER)((d(i,j),i=1,nbasis),j=1,nbasis)
      call calch(nbasis)
      call hyperprop(d,si,sf,ti,tf,ival,jval,gi,jinit,jfin)
      call calhyp
      endif

c---------------------------------------------------------------
c      O = Electric quadrupole hyperfine operator
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(iquadB.eq.1)then
      write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,9992)
      call cinit(d)
      open(QUADB,file='quadb.dat',form='unformatted',status='unknown')
      rewind(QUADB)
      read(QUADB)((d(i,j),i=1,nbasis),j=1,nbasis)
      call calcB(nbasis)
      call quadBprop(d,si,sf,ti,tf,ival,jval,qm,jinit,jfin)
      call calquadB
      endif


c---------------------------------------------------------------
c      O = Magnetic octupole hyperfine operator
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(ioctc.eq.1)then
      write(*,*)'______________________________________________'
      write(*,9993)
      call cinit(d)
      open(OCTC,file='octc.dat',form='unformatted',status='unknown')
      rewind(OCTC)
      read(OCTC)((d(i,j),i=1,nbasis),j=1,nbasis)
      call calcC(nbasis)
      call octCprop(d,si,sf,ti,tf,ival,jval,ocm,jinit,jfin)
      call calC
      endif
c----------------------------------------------------------------
c      o= normal mass shift operator
c---------------------------------------------------------------

      if(inms.eq.1)then
      write(*,*)'============================================='
      write(STDOUT,*)'sourav111'
      write(*,9988)
      call cinit(d)
      open(HYPER,file='hyper.dat',form='unformatted',status='unknown')
      rewind(HYPER)
      read(HYPER)((d(i,j),i=1,nbasis),j=1,nbasis)
      write(STDOUT,*)'sourav222'
cc      call calch(nbasis)
      call nmsp(d,si,sf,ti,tf,ival,jval,gi,jinit,jfin)
      call calhyp
      endif

c---------------------------------------------------------------
c      O = Magnetic dipole transition operator
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------

      if(im1p.eq.1)then
      write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      write(*,9994)
      call cinit(d)
      call calcm1prop(nbasis)
      open(M1PROP,file='m1prop.dat',form='unformatted',status='unknown')
      rewind(M1PROP)
      read(M1PROP)((d(i,j),i=1,nbasis),j=1,nbasis)
      call m1transprop(d,si,sf,ti,tf,ival,jval,jinit,jfin)
      call calm1
      endif

c---------------------------------------------------------------
c      O = Magnetic octupole transition operator
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(im2p.eq.1)then
      write(*,*)'##############################################'
      write(*,9995)
      call cinit(d)
      call calcm2prop(nbasis)
      open(M2PROP,file='m2prop.dat',form='unformatted',status='unknown')
      rewind(M2PROP)
      read(M2PROP)((d(i,j),i=1,nbasis),j=1,nbasis)
      call m2transprop(d,si,sf,ti,tf,ival,jval,jinit,jfin)
      call calm2
      endif

c---------------------------------------------------------------
c      O = G-Factor
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(igfact.eq.1)then
      write(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
      write(*,9996)
      call cinit(d)
      call calcgfact(nbasis)
      open(GFACTOR,file='gfact.dat',form='unformatted',status='unknown')
      open(GFACTORC,file='gfactc.dat',form='unformatted',
     :status='unknown')
      rewind(GFACTOR)
      read(GFACTOR)((d(i,j),i=1,nbasis),j=1,nbasis)
      rewind(GFACTORC)
      read(GFACTORC)((dc(i,j),i=1,nbasis),j=1,nbasis)
      call gfactprop(d,dc,si,sf,ti,tf,ival,jval,jinit,jfin)
      call calgfact
      endif


c---------------------------------------------------------------
c          O = Nuclear spin independent (NSI) PNC operator
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(insi.eq.1)then
      write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,9997)
      call cinit(d)
      open(NSI,file='nsi.dat',form='unformatted',status='unknown')
      rewind(NSI)
      read(NSI)((d(i,j),i=1,nbasis),j=1,nbasis)
      call calcnsira(nbasis)
      call nsiprop(d,si,sf,ti,tf,ival,jval)
      call calnsi
      endif

c---------------------------------------------------------------
c         O = Nuclear spin dependent (NSD) PNC operator
c
c             T+   T
c     {1+Sf+}e  O e{1+Si}diagrams are calculated
c---------------------------------------------------------------


      if(insd.eq.1)then
      write(*,*)'""""""""""""""""""""""""""""""""""""""""""""""'
      write(*,9998)
      call cinit(d)
      open(NSD,file='nsd.dat',form='unformatted',status='unknown')
      rewind(NSD)
      read(NSD)((d(i,j),i=1,nbasis),j=1,nbasis)
      call calcnsdra(nbasis)
      call nsdprop(d,si,sf,ti,tf,ival,jval)
      call calnsd
      endif



      write(*,*)'&^&^&^&^&^&^&^&@It is over@&^&^&^&^&^&^&^&^&^&^&^&'

      call daclose
 9988 format(//,1x,' NMS matrix elements',//)
 9989 format(//,1x,' E1 transition ',//)
 9990 format(//,1x,' E2 transition  ',//)
 9991 format(//,1x,' Hyperfine constants(A) ',//)
 9992 format(//,1x,' Hyperconstants(B) ',//)
 9993 format(//,1x,' Hyperconstants(C) ',//)
 9994 format(//,1x,' M1 transition  ',//)
 9995 format(//,1x,' M2 transition  ',//)
 9996 format(//,1x,' g-fact calculation  ',//)
 9997 format(//,1x,' NSI PNC matrix element ',//)
 9998 format(//,1x,' NSD PNC matrix element ',//)

      stop
      end

c***********************************************************************

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
      if(mfd.eq.33)goto 704
      mfd=33
      fct(0)=1.
      fct(1)=1.
      do 20 i=2,33
      fct(i)=fct(i-1)*i
20    continue

704   tria=tri6j(rj1,rj2,rj3)
      trib=tri6j(rl1,rj2,rl3)
      tric=tri6j(rj1,rl2,rl3)
      trid=tri6j(rl1,rl2,rj3)
c     type *,tria,trib,tric,trid
      if(tria.eq.-1..or.trib.eq.-1..or.tric.eq.-1..or.trid.eq.-1.)then
      print *,'inconsistent arguments to 6j symbol'
      print *,tria,trib,tric,trid
      print *,rj1,rj2,rj3
      print *,rl1,rl2,rl3
      return
      endif
      if(tria.eq.0..or.trib.eq.0..or.tric.eq.0..or.trid.eq.0.)then
c     type *,'triangle condition'
      s6j=0
      return
      endif
!    comment apply formula in edmonds [(6.3.7) p.99]
      delprod=tridel(rj1,rj2,rj3)*tridel(rj1,rl2,rl3)*
     *        tridel(rl1,rj2,rl3)*tridel(rl1,rl2,rj3)
c     type *,'delprod=',delprod
c     type *,'args to min',(rj1+rj2+rl1+rl2),(rj2+rj3+rl2+rl3),
c     *                       (rj3+rj1+rl3+rl1)
      izmin=dmax1((rj1+rj2+rj3),(rj1+rl2+rl3),(rl1+rj2+rl3)
     *          ,(rl1+rl2+rj3))
      izmax=dmin1((rj1+rj2+rl1+rl2),(rj2+rj3+rl2+rl3)
     *          ,(rj3+rj1+rl3+rl1))
      ws6j=0
c     type *,'izmin,izmax=',izmin,izmax
      do 10 iz=izmin,izmax
      z=iz
      ws6j=ws6j+((-1)**iz)*fct(int(z+1.))/(fct(z-rj1-rj2-rj3)*
     *         fct(int(z-rj1-rl2-rl3))*
     *         fct(int(z-rl1-rj2-rl3))*
     *         fct(int(z-rl1-rl2-rj3))*
     *         fct(int(rj1+rj2+rl1+rl2-z))*
     *         fct(int(rj2+rj3+rl2+rl3-z))*
     *         fct(int(rj3+rj1+rl3+rl1-z)))
10    continue
      s6j=delprod*ws6j
      return
      end

      function tridel(rm1,rm2,rm3)
      implicit real*8 (a-h,o-z)
      common/factor/fct(0:50),mfd

      i1=rm1+rm2-rm3
      i2=rm1-rm2+rm3
      i3=-rm1+rm2+rm3
      i4=rm1+rm2+rm3+1

      tridel=sqrt(fct(i1)*fct(i2)*fct(i3)/fct(i4))
c     type *,'***',rm1,rm2,rm3

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
     *   (2.*(rj1+rj2).ge.2.*rj3))then
      tri6j=1                               !triangle condition is satisfied
      goto 999
      endif

      tri6j=0                         !j1+j2+j3 is an integer
c     type *,rj1,rj2,rj3
999   return
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

      subroutine ninej(rj11,rj12,rj13,rj21,rj22,rj23,rj31,rj32,
     : rj33,r9j)
      implicit real*8(a-h,o-z)
      common/factor/fct(0:50),mfd

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

      if(tria.eq.-1.or.trib.eq.-1.or.tric.eq.-1.or.
     *  trid.eq.-1.or.trie.eq.-1.or.trif.eq.-1)then
      print *,'inconsistent arguments to 9j symbol'
      endif

      if(tria.eq.0.or.trib.eq.0.or.tric.eq.0.or.
     *  trid.eq.0.or.trie.eq.0.or.trif.eq.0)then
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

      mtbl(2,4)=1
      mtbl(4,2)=1

      mtbl(1,5)=1
      mtbl(5,1)=1

      mtbl(3,5)=1
      mtbl(5,3)=1

      mtbl(2,6)=1
      mtbl(6,2)=1
      mtbl(4,6)=1
      mtbl(6,4)=1


c     odd parity

      mtbl(1,2)=2
      mtbl(2,1)=2

      mtbl(1,4)=2
      mtbl(4,1)=2

      mtbl(2,3)=2
      mtbl(3,2)=2
      mtbl(1,6)=2
      mtbl(6,1)=2
      mtbl(3,6)=2
      mtbl(6,3)=2

      mtbl(5,6)=2
      mtbl(6,5)=2

      mtbl(3,4)=2
      mtbl(4,3)=2

      mtbl(2,5)=2
      mtbl(5,2)=2

      mtbl(5,4)=2
      mtbl(4,5)=2

      return
      end

      subroutine setgrd
      implicit real*8 (a-h,o-z)
      common/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
      n = MN-10
      mtp=n
      h = 5.0d-02
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

      open(STDIN,file='input.prop',form='formatted',status='old')
      open(STDOUT,file='prop_opt.out',form='formatted',status='unknown')
      open(WFNIN,file='wfn.dat',form='unformatted',status='old')
      open(NTFILE,file='to.dat',form='unformatted',status='unknown')
      open(OVERLAP,file='overlap.dat',form='unformatted',
     :      status='unknown')
      return
      end

c***********************************************************************

c***********************************************************************
c     This subroutine constructs the inpute file of the property       *
c     code                                                             *
c***********************************************************************

      subroutine readinp(nsym,nbas,iocc,iopt,nprint,idip,ivel,iquad,
     :ihyp,iquadB,ioctc,im1p,im2p,igfact,ival,jval,gi,qm,ocm,jinit,
     :jfin,insi,insd,inms)

      implicit real*8(a-h,o-z)
      real*8 jinit,jfin
      common/iallparity/iall
      common/inparm/nparm
      common/pncp/amass,qw
      dimension nbas(MNSYM),iocc(MNSYM)

c------------------------------------------------------------------------
c     read number of symmetry
c------------------------------------------------------------------------

      read(STDIN,*)nsym,(nbas(i),i=1,nsym)

c-------------------------------------------------------------------------
c     occupancy of each orbital
c-------------------------------------------------------------------------

      read(STDIN,*)(iocc(i),i=1,nsym)

c--------------------------------------------------------------------------
c     truncation scheme and option flag
c--------------------------------------------------------------------------

      read(STDIN,*)iopt,nprint,iall

c---------------------------------------------------------------------------
c     nuclear property
c---------------------------------------------------------------------------

      read(STDIN,*)nparm,qw,amass

c---------------------------------------------------------------------------
c     information of active orbital
c---------------------------------------------------------------------------

      read(STDIN,*)ival,jval,gi,qm,ocm,jinit,jfin
      write(STDOUT,*)'Initial and Final orbitals = ', ival,jval

c---------------------------------------------------------------------------
c     property related optional flags
c---------------------------------------------------------------------------

      read(STDIN,*)idip,ivel,iquad,ihyp,iquadB,ioctc,im1p,im2p
     :insi,insd,igfact
      read(STDIN,*)inms

c----------------------------------------------------------------------------
c     information of active orbital
c----------------------------------------------------------------------------

      return
      end

c******************************************************************************

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
      dimension orbtyp(MNS)
      dimension orb(MNSYM),np(MNBAS)
      dimension qj(MNBAS),qj1(MNS)
      dimension nprimc(MNOCC),nkprimc(MNOCC)
      dimension nprime(MNEXC),nkprime(MNEXC)
      dimension nmc(MNBAS),nme(MNBAS)
      dimension nbas(MNSYM),iocc(MNSYM)

      data kap/-1, 1,-2, 2,-3, 3,-4,4,-5,5,-6/
      data orb/0.5,0.5,1.5,1.5,2.5,2.5,3.5,3.5,4.5,4.5,5.5/
      data ia /1,-1,1,-1,1,-1,1,-1,1,-1,1/
      data iorbsm/1,2,2,3,3,4,4,5,5,6,6/
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
      iparity(ii)=iorbsm(isym)
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
      DATA FIRST /.TRUE./
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
      return
      end

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
      common/skip1/npbps(0:MXV,2)

c     local dimension
      dimension kgot(MXV+1),ieven(0:MXV),iodd(0:MXV)
      dimension indx(MDIM),indx1(MDIM)
      dimension ti(MDIM),tf(MDIM),ncheck(MDIM)
      
      nsing=0

      do i=1, MDIM
         ncheck(i)=0
      enddo

c------------------ set the equation index for T1-----------------------------

      if(iopt.eq.1)then

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

c------------------------Set the skip information-------------------------------

      do 2 i=0,MXV
      ieven(i)=0
      iodd(i)=0
      npbps(i,1)=0
      npbps(i,2)=0
      nnskip(i,1)=0
   2  nnskip(i,2)=0

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

      call findk(orba,orbp,orba,orbp,iaa,iap,iaa,iap,kgot,
     $kmax)
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

c-----------------------iall=0 for only even rank----------------------------

      do 5 i=2,MXV,2
      nnskip(i,2)=nnskip(i-2,2)+ieven(i-2)*(ieven(i-2)+1)/2
  5   continue

      do  i=0,JMAX,2
      n2skip=n2skip+ieven(i)*(ieven(i)+1)/2
      enddo

      do 6 i=3,MXV-2,2
      nnskip(i,1)=nnskip(i-2,1)+iodd(i-2)*(iodd(i-2)+1)/2
  6   continue

      do 7 i=0,MXV,2
      npbps(i,2)=ieven(i)
  7   continue

      do 8 i=1,MXV,2
      npbps(i,1)=iodd(i)
  8   continue

      else

      do 13 i=1,MXV
      nnskip(i,2)=nnskip(i-1,2)+ieven(i-1)*(ieven(i-1)+1)/2
 13   continue

      do  i=0,MXV
      n2skip=n2skip+ieven(i)*(ieven(i)+1)/2
      enddo

      do 14 i=1,MXV
      nnskip(i,1)=nnskip(i-1,1)+iodd(i-1)*(iodd(i-1)+1)/2
 14   continue

      do 15 i=0,MXV
      npbps(i,2)=ieven(i)
 15   continue

      do 16 i=0,MXV
      npbps(i,1)=iodd(i)
 16   continue

      endif

      nn=nsing

c------------------------ set the equation index for T2-------------------------

      do 17 ia=1,nocc
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

      if (ncheck(neqn).eq.1) then
         print*,'double count symm =', neqn
      stop
      endif

      ncheck(neqn)=1

      indx(nn)=nn
      indx1(nn)=neqn

      tf(neqn)=ti(neqn)*(-1)**(orbp+orba+orbq+orbb)

  45  continue
  40  continue
  30  continue
  20  continue
  17  continue

      ntmax=nn

      do i=1,nn
      do j=1,nn-1
      if (indx1(j).gt.indx1(j+1))then
      ind=indx1(j)
      indx1(j)=indx1(j+1)
      indx1(j+1)=ind
      endif
      enddo
      enddo


      do i=1,nn
      ii=indx1(i)-indx(i)
      if(ii.ne.0)then
      write(*,*)indx(i),indx1(i),indx1(i)-indx(i)
      endif
      enddo

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

c*********************************************************************************


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

c**********************************************************************
c    This subroitine calculates the electric dipole property diagram  *
c**********************************************************************

      subroutine dipoleprop(d,si,sf,ti,tf,ival,jval)

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
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)
 
      common/valued/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,sum45,sum46
 
      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS),si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo
   
      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0
      sum45=0
      sum46=0

      irank=1
      rank=1.0d0

c--------p is the final and q is the initial state-----------

      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iiq.ne.ival)go to 1
      if(ipsym.eq.iqsym) go to 1

c-------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c-------------------------------------------------------
               
      sum=sum+d(iip,iiq)   
      call e1mat(iip,iiq,orbp,orbq,kapp,kapq,dirac)    
      sum1=sum1+d(iip,iiq)


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

c     if(iir.ne.ival)then
      if(irsym.ne.ipsym)then
      if(orbq.eq.orbr)then
      if(irsym.eq.iqsym)then
      if(iir.ne.iiq)then

      in=jdpa(iiq,iir)

c------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c------------------------------------------------------

      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)
  
      endif
      endif
      endif
      endif

      if(irsym.ne.iqsym)then
      if(orbr.eq.orbp)then

      if(irsym.eq.ipsym)then
      if(iip.ne.iir)then

      in=kdpa(iip,iir)

c-------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-------------------------------------------------------

      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)
   
      endif
      endif
      endif
      endif
c     endif
      
      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.eq.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact

  6   continue

      endif

c----------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact
  
  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      in=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(in)*fact
      sum6=sum6+d(iir,iig)*sf(in)*fact
 
  9   continue

      endif

c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      in=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(in)*fact     
      sum7=sum7-d(iir,iig)*sf(in)*fact
  
 10   continue
  8   continue
  4   continue
  3   continue

      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c--------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c--------------------------------------------------------------

      if(issym.eq.irsym)go to 11
      if(orbp.ne.orbr)go to 11
      if(ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs)go to 11
      if(iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)
 
 11   continue
      
      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.eq.igsym)go to 12

      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      
      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(iqssym.ne.irgsym)go to 14

c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------
      
      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact
  
 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)              
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact
 
 18   continue
 17   continue
 14   continue

      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.and.iqsym.ne.issym)go to 20

      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c-----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c-----------------------------------------------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact 

 22   continue
 21   continue

c-----------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c-----------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)          
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue
  2   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.eq.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact
 
      endif
      enddo
      endif

  42  continue

c---------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c---------------------------------------------------------------

      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact
 
 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)      
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact 

 48   continue
 47   continue
 46   continue

c-----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c-----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact
  
 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)

      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52
      if(igsym.eq.idsym)go to 52

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)

c------------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c------------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      irdsym=mtbl(irsym,idsym)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)      
      fact=(-1)**(orbp+orbq+irank+2*orbg)*s6j1*s6j2
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact
   
  54  continue

c-----------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact

  55  continue
  53  continue
  52  continue
  30  continue

c-------Two-body effective operator (computed directly)------------

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

      hyp1=0

      if(ipsym.ne.issym)then
      call e1mat(iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif

      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then

      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)

      if(iqssym.eq.irasym)then

c---------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c---------------------------------------------------------------

      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact

      endif
      enddo
      endif

c---------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c---------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif

      hyp2=0

      if(issym.ne.iqsym)then
      call e1mat(iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      if(ipssym.eq.irasym)then

c-----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif      
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      if(irasym.eq.iqssym)then
      
      hyp1=0

      if(itsym.ne.issym)then
      call e1mat(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif

      endif
      endif
      endif
      endif
      endif


      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      irpsym=mtbl(irsym,ipsym)
      itasym=mtbl(itsym,iasym)
      if(irpsym.eq.itasym)then
      
      hyp1=0

      if(itsym.ne.issym)then
      call e1mat(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif
    
      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)

c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif

      endif
      endif
      endif
      endif
      endif

      enddo
      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)
     
      hyp1=0

      if(ipsym.ne.ibsym)then
      call e1mat(iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.iqssym)then
      if(irasym.eq.isbsym)then
      if(orbq.eq.orbb)then

c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      hyp2=0

      if(iqsym.ne.ibsym)then
      call e1mat(iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.ipssym)then
      if(irasym.eq.isbsym)then
      if(orbp.eq.orbb)then

c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      
      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      hyp=0

      if(iasym.eq.igsym)go to 610

      call e1mat(iig,iia,orbg,orba,kapg,kapa,hyp)
   
      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
c     if(orbb.ne.orbq)go to 620

      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

c--------------------------------------------------------------------
c     t2dt1s2     diagram 27 = <rs|t2|ab><a|d|g><b|t1|p><sg|s2|qr>
c--------------------------------------------------------------------
     
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+irank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact

 640  continue
      endif

c------------------------------------------------------------------
c     t2dt1s2     diagram 28 = <rs|t2|ba><a|d|g><b|t1|p><sg|s2|qr>
c------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact

  650 continue

      endif

  630 continue
  620 continue

      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
c     if(orbb.ne.orbp)go to 660

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2dt1s2     diagram 29 = <pa|s2|sr><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+irank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact

  680 continue
      endif

c-------------------------------------------------------------------
c     s2t1dt2     diagram 30 = <pa|s2|rs><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------
      
      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact

  690 continue
      endif

  670 continue
  660 continue
  610 continue

      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.ne.isgsym)go to 605
      if(irasym.ne.iqssym)go to 605
      
      if(ibsym.eq.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
      if(orbp.ne.orbb)go to 605

      in=idpa(iip,iib)

      hyp=0

      call e1mat(iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700

c---------------------------------------------------------------------
c     t2dt1s2     diagram 31 = <rs|t2|ag><g|d|b><p|t1|b><rs|s2|aq>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

c----------------------------------------------------------------------
c     t2dt1s2     diagram 32 = <rs|t2|ga><g|d|b><p|t1|b><rs|s2|aq>
c----------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact

  720 continue
  710 continue
  605 continue
  700 continue

      hyp=0

      if(ibsym.ne.igsym)then
      call e1mat(iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)

      if(ipssym.ne.irasym)go to 600
      if(isgsym.ne.irasym)go to 600
      
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 33 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ag>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

c------------------------------------------------------------------
c     s2t1dt2     diagram 34 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbp+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact

  750 continue
  740 continue
  730 continue
  600 continue
     
      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      
      hyp=0

      if(itsym.ne.issym)then
      call e1mat(iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      isbsym=mtbl(issym,ibsym)
      irasym=mtbl(irsym,iasym)
      iqtsym=mtbl(iqsym,itsym)

      if(isbsym.ne.irasym)go to 815
      if(iqtsym.ne.irasym)go to 815
      
      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 35 = <rs|t2|ab><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)

      call findk(orbs,orbb,orbr,orba,ias,iab,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 812
      do 811 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      if(l2.eq.l1)then
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+irank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

      endif

  811 continue
  812 continue

c-------------------------------------------------------------------
c     t2dt1s2     diagram 45 = <rs|t2|ba><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbb,ias,iaa,iar,iab,kgot,kmax)
      if(kmax.eq.0)go to 814
      do 813 kloop=1,kmax
      l1=kgot(kloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrb=idra(iir,iib,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrb,isasym,l1)
      call sixj(orbp,rank,orbq,orbt,al2,orbs,s6j1)
      call sixj(orbp,al2,orbs,orba,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbt+orbs+l2+orbr+irank+orba)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum45=sum45+hyp*tf(inf)*si(ini)*ti(in)*fact

  813 continue
  814 continue  
  830 continue
  820 continue
  815 continue

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      itbsym=mtbl(itsym,ibsym)

      if(ipssym.ne.irasym)go to 900
      if(itbsym.ne.irasym)go to 900

      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 36 = <pr|s2|sa><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------


      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)

      call findk(orbp,orbs,orbr,orba,iap,ias,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 912
      do 911 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      if(l2.eq.l1)then
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+irank+orba+orbr+al2)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact

      endif

  911 continue
  912 continue

c-------------------------------------------------------------------
c     s2t1dt2     diagram 46 = <ps|s2|ra><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------

      call findk(orbp,orbr,orbs,orba,iap,iar,ias,iaa,kgot,kmax)
      if(kmax.eq.0)go to 914
      do 913 kloop=1,kmax
      l1=kgot(kloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      nsa=idra(iis,iia,l1)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l1)
      call sixj(orbp,rank,orbq,orbt,al2,orbs,s6j1)
      call sixj(orbp,al2,orbs,orba,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbt+orbs+orbr+orba+l2+irank)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum46=sum46+hyp*ti(ini)*sf(inf)*tf(in)*fact

  913 continue
  914 continue
  930 continue
  900 continue
  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)

      hyp=0

      if(itsym.ne.ibsym)then
      call e1mat(iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      isqsym=mtbl(issym,iqsym)

      if(irasym.eq.isbsym)then
      if(irasym.eq.isqsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

c-------------------------------------------------------------------
c     s1t2ds2     diagram 41 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|aq>
c-------------------------------------------------------------------
      
      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nra=idra(iir,iia,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t2ds2     diagram 42 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|qa>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif

      enddo
      endif

      endif
      endif

      endif
      endif
      endif
      endif


      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then

      in=jdpa(iiq,iit)

      hyp=0

      if(itsym.ne.ibsym)then
      call e1mat(iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      ispsym=mtbl(issym,ipsym)

      if(irasym.eq.isbsym)then
      if(irasym.eq.ispsym)then

c-------------------------------------------------------------------
c     s2dt1s1     diagram 43 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ba>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

c--------------------------------------------------------------------
c     s2dt1s1     diagram 44 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ab>
c--------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue

      print*,'Dirac Fock Value for E1 =', dirac

   1  continue
 100  continue

      return
      end


c*****************************************************************

c******************************************************************
c     This subroutine calculates E1 reduced matrix element        *
c******************************************************************

      subroutine e1mat(i1,i2,jorb1,jorb2,kap1,kap2,e1m)

      implicit double precision (a-h,o-z)
      double precision j1,j2,jorb1,jorb2
      common/cons/zero,half,tenth,one,two,three,ten
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

c*********************************************************************


c*********************************************************************
c      This subroutine calculate the radial integral part for E1     *
c      matrix element                                                *
c*********************************************************************

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
c   The value of this  frint is an approximation to:                   *
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

c*********************************************************************


c**********************************************************************
c    This subroitine calculates the NSI PNC property diagram          *
c**********************************************************************

      subroutine nsiprop(d,si,sf,ti,tf,ival,jval)

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
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)
      common/inparm/nparm
      common/pncp/amass,qw


      common/valuensi/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS),si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0

      irank=0
      rank=0.0d0

c--------p is the final and q is the initial state-----------

      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iiq.ne.ival)go to 1
      if(ipsym.eq.iqsym) go to 1

c-------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c-------------------------------------------------------

      sum=sum+d(iip,iiq)
      call nsimat(iip,iiq,orbp,orbq,kapp,kapq,dirac)
      sum1=sum1+d(iip,iiq)


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(iir.ne.ival)then
      if(irsym.ne.ipsym)then
      if(orbq.eq.orbr)then
      if(irsym.eq.iqsym)then
      if(iir.ne.iiq)then

      in=jdpa(iiq,iir)

c------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c------------------------------------------------------

      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)

      endif
      endif
      endif
      endif

      if(irsym.ne.iqsym)then
      if(orbr.eq.orbp)then
      if(irsym.eq.ipsym)then
      if(iip.ne.iir)then

      in=kdpa(iip,iir)

c-------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-------------------------------------------------------

      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)

      endif
      endif
      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.eq.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact

  6   continue

      endif

c----------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact

  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      in=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(in)*fact
      sum6=sum6+d(iir,iig)*sf(in)*fact

  9   continue

      endif

c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      in=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(in)*fact
      sum7=sum7-d(iir,iig)*sf(in)*fact

 10   continue
  8   continue
  4   continue
  3   continue

      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c--------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c--------------------------------------------------------------

      if(issym.eq.irsym)go to 11
      if(orbp.ne.orbr)go to 11
      if(ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs)go to 11
      if(iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      in2=kdpa(iip,iir)
      in1=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(in1)*sf(in2)
      sum8=sum8+d(iir,iis)*si(in1)*sf(in2)

 11   continue

      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.eq.igsym)go to 12

      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(iqssym.ne.irgsym)go to 14

c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact

 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact

 18   continue
 17   continue
 14   continue

      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.or.iqsym.ne.issym)go to 20

      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c-----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c-----------------------------------------------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact

 22   continue
 21   continue

c-----------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c-----------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue
  2   continue


      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.eq.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+rank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

  42  continue

c---------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c---------------------------------------------------------------

      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(al2,orbp,orbt,al1,orbg,orbr,s6j1)
      call sixj(al2,orbt,orbp,rank,orbq,orbs,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+rank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact

 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(al1,al2,rank,orbq,orbp,orbr,s6j1)
      call sixj(al1,al2,rank,orbs,orbt,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact

 48   continue
 47   continue
 46   continue

c-----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c-----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+rank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact

 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)

      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52
      if(igsym.eq.idsym)go to 52

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)

c------------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c------------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      irdsym=mtbl(irsym,idsym)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(al2,al1,rank,orbp,orbq,orbs,s6j1)
      call sixj(al1,al2,rank,orbd,orbg,orbr,s6j2)
      fact=(-1)**(orbp+orbq+rank+2*orbg)*s6j1*s6j2
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact

  54  continue

c-----------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2+2*(orbd+orbr))*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact

  55  continue
  53  continue
  52  continue
  30  continue

c----------------------------------------------------------------
c     Two-body effective operator (computed directly) 
c----------------------------------------------------------------

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

      hyp1=0

      if(ipsym.ne.issym)then
      call nsimat(iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif

      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then
      iqssym=mtbl(iqsym,issym)

c---------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c---------------------------------------------------------------

      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact

      endif
      enddo
      endif

c---------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c---------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      endif

      hyp2=0

      if(issym.ne.iqsym)then
      call nsimat(iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)
      ipssym=mtbl(ipsym,issym)

c-----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

      hyp1=0

      if(itsym.ne.issym)then
      call nsimat(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0

      if(itsym.ne.issym)then
      call nsimat(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)

c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      enddo
      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

      hyp1=0

      if(ipsym.ne.ibsym)then
      call nsimat(iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.eq.iqssym)then
      if(orbq.eq.orbb)then

c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      hyp2=0

      if(iqsym.ne.ibsym)then
      call nsimat(iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)
      if(irasym.eq.ipssym)then
      if(orbp.eq.orbb)then

c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.eq.igsym)go to 610

      hyp=0

      call nsimat(iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620

      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

c--------------------------------------------------------------------
c     t2dt1s2     diagram 27 = <rs|t2|ab><a|d|g><b|t1|p><sg|s2|qr>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact

 640  continue

      endif

c------------------------------------------------------------------
c     t2dt1s2     diagram 28 = <rs|t2|ba><a|d|g><b|t1|p><sg|s2|qr>
c------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact

  650 continue

      endif

  630 continue
  620 continue

      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2dt1s2     diagram 29 = <pa|s2|sr><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact

  680 continue

      endif

c-------------------------------------------------------------------
c     s2t1dt2     diagram 30 = <pa|s2|rs><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact

  690 continue

      endif

  670 continue
  660 continue
  610 continue


      if(ibsym.eq.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
      if(orbp.ne.orbb)go to 605

      in=idpa(iip,iib)

      hyp=0

      call nsimat(iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700

c---------------------------------------------------------------------
c     t2dt1s2     diagram 31 = <rs|t2|ag><g|d|b><p|t1|b><rs|s2|aq>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

c----------------------------------------------------------------------
c     t2dt1s2     diagram 32 = <rs|t2|ga><g|d|b><p|t1|b><rs|s2|aq>
c----------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact

  720 continue
  710 continue
  605 continue
  700 continue

      hyp=0

      if(ibsym.ne.igsym)then
      call nsimat(iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 33 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ag>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

c------------------------------------------------------------------
c     s2t1dt2     diagram 34 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbp,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact

  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0

      if(itsym.ne.issym)then
      call nsimat(iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 35 = <rs|t2|ab><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)

      call findk(orbs,orbb,orbr,orba,ias,iab,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 812
      do 811 jloop=1,jmax
      l3=jgot(jloop)
      al3=float(l3)
      if(l2.eq.l3)then      
      nsb=idra(iis,iib,l2)      
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

      endif

  811 continue
  812 continue 
  830 continue
  820 continue
  815 continue


      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 36 = <pr|s2|sa><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)

      call findk(orbp,orbs,orbr,orba,iap,ias,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 912
      do 911 jloop=1,jmax
      l3=jgot(jloop)
      al3=float(l3)
      if(l2.eq.l3)then
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact

      endif

  911 continue
  912 continue
  930 continue
  900 continue
  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)

      hyp=0

      if(itsym.ne.ibsym)then
      call nsimat(iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

c-------------------------------------------------------------------
c     s1t2ds2     diagram 41 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nra=idra(iir,iia,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t2ds2     diagram 42 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|qa>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif
      enddo
      endif
      endif

      endif
      endif
      endif
      endif


      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then

      in=jdpa(iiq,iit)

      hyp=0

      if(itsym.ne.ibsym)then
      call nsimat(iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.isbsym)then

c-------------------------------------------------------------------
c     s2dt1s1     diagram 43 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ba>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

c--------------------------------------------------------------------
c     s2dt1s1     diagram 44 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ab>
c--------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue

      print*,'Dirac Fock Value for NSI =', dirac

   1  continue
 100  continue

      return
      end

c*****************************************************************

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

c*******************************************************************

c*********************************************************************
c      The reduced matrix element of the NSI PNC matrix element is   *
c      calculated                                                    *
c*********************************************************************

        subroutine nsimat(i1,i2,jorb1,jorb2,kap1,kap2,ansim)

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

c*********************************************************************


c**********************************************************************
c    This subroitine calculates the NSD PNC property diagram          *
c**********************************************************************

      subroutine nsdprop(d,si,sf,ti,tf,ival,jval)

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
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)
      common/inparm/nparm
      common/pncp/amass,qw

      common/valuensd/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS),si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0


      irank=1
      rank=1.0d0

c--------p is the final and q is the initial state-----------

      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iiq.ne.ival)go to 1
      if(ipsym.eq.iqsym) go to 1

c-------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c-------------------------------------------------------

      sum=sum+d(iip,iiq)
      call nsdmat(iip,iiq,orbp,orbq,kapp,kapq,dirac)
      sum1=sum1+d(iip,iiq)


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(iir.ne.ival)then
      if(irsym.ne.ipsym)then
      if(orbq.eq.orbr)then
      if(irsym.eq.iqsym)then
      if(iir.ne.iiq)then

      in=jdpa(iiq,iir)

c------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c------------------------------------------------------

      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)

      endif
      endif
      endif
      endif

      if(irsym.ne.iqsym)then
      if(orbr.eq.orbp)then
      if(irsym.eq.ipsym)then
      if(iip.ne.iir)then

      in=kdpa(iip,iir)

c-------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-------------------------------------------------------

      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)

      endif
      endif
      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.eq.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact

  6   continue

      endif

c----------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact

  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      in=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(in)*fact
      sum6=sum6+d(iir,iig)*sf(in)*fact

  9   continue

      endif

c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      in=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(in)*fact
      sum7=sum7-d(iir,iig)*sf(in)*fact

 10   continue
  8   continue
  4   continue
  3   continue

      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c--------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c--------------------------------------------------------------

      if(issym.eq.irsym)go to 11
      if(orbp.ne.orbr)go to 11
      if(ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs)go to 11
      if(iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      in2=kdpa(iip,iir)
      in1=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(in1)*sf(in2)
      sum8=sum8+d(iir,iis)*si(in1)*sf(in2)

 11   continue

      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.eq.igsym)go to 12

      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(iqssym.ne.irgsym)go to 14

c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact

 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact

 18   continue
 17   continue
 14   continue

      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.or.iqsym.ne.issym)go to 20

      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c-----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c-----------------------------------------------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact

 22   continue
 21   continue

c-----------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c-----------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue
  2   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.eq.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+rank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

  42  continue

c---------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c---------------------------------------------------------------

      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(al2,orbp,orbt,al1,orbg,orbr,s6j1)
      call sixj(al2,orbt,orbp,rank,orbq,orbs,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+rank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact

 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(al1,al2,rank,orbq,orbp,orbr,s6j1)
      call sixj(al1,al2,rank,orbs,orbt,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact

 48   continue
 47   continue
 46   continue

c-----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c-----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+rank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact

 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)

      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52
      if(igsym.eq.idsym)go to 52

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)


c------------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c------------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      irdsym=mtbl(irsym,idsym)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(al2,al1,rank,orbp,orbq,orbs,s6j1)
      call sixj(al1,al2,rank,orbd,orbg,orbr,s6j2)
      fact=(-1)**(orbp+orbq+rank+2*orbg)*s6j1*s6j2
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact

  54  continue

c-----------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2+2*(orbd+orbr))*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact

  55  continue
  53  continue
  52  continue
  30  continue

c----------------------------------------------------------------
c     Two-body effective operator (computed directly)
c----------------------------------------------------------------

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

      hyp1=0

      if(ipsym.ne.issym)then
      call nsdmat(iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif

      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then
      iqssym=mtbl(iqsym,issym)

c---------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c---------------------------------------------------------------

      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact

      endif
      enddo
      endif

c---------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c---------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      endif


      hyp2=0

      if(issym.ne.iqsym)then
      call nsdmat(iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)
      ipssym=mtbl(ipsym,issym)

c-----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

      hyp1=0

      if(itsym.ne.issym)then
      call nsdmat(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0

      if(itsym.ne.issym)then
      call nsdmat(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)


c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      enddo
      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

      hyp1=0

      if(ipsym.ne.ibsym)then
      call nsdmat(iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.eq.iqssym)then
      if(orbq.eq.orbb)then

c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      hyp2=0

      if(iqsym.ne.ibsym)then
      call nsdmat(iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)
      if(irasym.eq.ipssym)then
      if(orbp.eq.orbb)then


c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.eq.igsym)go to 610

      hyp=0

      call nsdmat(iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620

      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

c--------------------------------------------------------------------
c     t2dt1s2     diagram 27 = <rs|t2|ab><a|d|g><b|t1|p><sg|s2|qr>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact

 640  continue

      endif

c------------------------------------------------------------------
c     t2dt1s2     diagram 28 = <rs|t2|ba><a|d|g><b|t1|p><sg|s2|qr>
c------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact

  650 continue

      endif

  630 continue
  620 continue


      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2dt1s2     diagram 29 = <pa|s2|sr><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact

  680 continue

      endif

c-------------------------------------------------------------------
c     s2t1dt2     diagram 30 = <pa|s2|rs><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact

  690 continue

      endif

  670 continue
  660 continue
  610 continue


      if(ibsym.eq.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
      if(orbp.ne.orbb)go to 605

      in=idpa(iip,iib)

      hyp=0

      call nsdmat(iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700

c---------------------------------------------------------------------
c     t2dt1s2     diagram 31 = <rs|t2|ag><g|d|b><p|t1|b><rs|s2|aq>
c---------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

c----------------------------------------------------------------------
c     t2dt1s2     diagram 32 = <rs|t2|ga><g|d|b><p|t1|b><rs|s2|aq>
c----------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact

  720 continue
  710 continue
  605 continue
  700 continue


      hyp=0

      if(ibsym.ne.igsym)then
      call nsdmat(iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 33 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ag>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

c------------------------------------------------------------------
c     s2t1dt2     diagram 34 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ga>
c------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbp,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact

  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0

      if(itsym.ne.issym)then
      call nsdmat(iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 35 = <rs|t2|ab><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)

      call findk(orbs,orbb,orbr,orba,ias,iab,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 812
      do 811 jloop=1,jmax
      l3=jgot(jloop)
      if(l2.eq.l3)then
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

      endif

  811 continue
  812 continue
  830 continue
  820 continue
  815 continue


      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 36 = <pr|s2|sa><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------


      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)

      call findk(orbp,orbs,orbr,orba,iap,ias,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 912
      do 911 jloop=1,jmax
      l3=jgot(jloop)
      al3=float(l3)
      if(l2.eq.l3)then
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact

      endif

  911 continue
  912 continue
  930 continue
  900 continue
  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)

      hyp=0

      if(itsym.ne.ibsym)then
      call nsdmat(iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

c-------------------------------------------------------------------
c     s1t2ds2     diagram 41 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nra=idra(iir,iia,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t2ds2     diagram 42 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|qa>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif
      enddo
      endif
      endif

      endif
      endif
      endif
      endif


      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then

      in=jdpa(iiq,iit)

      hyp=0

      if(itsym.ne.ibsym)then
      call nsdmat(iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.isbsym)then

c-------------------------------------------------------------------
c     s2dt1s1     diagram 43 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ba>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

c--------------------------------------------------------------------
c     s2dt1s1     diagram 44 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ab>
c--------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif


  800 continue
  300 continue
  200 continue

      print*,'Dirac Fock Value for NSD =', dirac

   1  continue
 100  continue

      return
      end

c*****************************************************************

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


c****************************************************************

c****************************************************************
c    This subroitine calculates the hyperfine property diagram  *
c****************************************************************

      subroutine hyperprop(d,si,sf,ti,tf,
     :             ival,jval,gi,jinit,jfin)


      implicit real*8 (a-h,o-z)
      real*8 jinit,jfin

      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      common/value/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,sum45,sum46
     
      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/
      data VL/137.03598950D+00/
      data EMF/1836.152701D+00/
      data AUEV/27.2116D+00/
      data EVMHZ/2.417988D+08/

      BETAE=ONE/(TWO*VL)
      BETAN=BETAE/EMF
      CONVF=AUEV*EVMHZ
      CON=GI*BETAN*CONVF

      print*,'CON',con

      d1=dr(jinit,1.0d0,jfin,-jinit,0.0d0,jfin)
      con=con*d1

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0
      sum45=0
      sum46=0

      irank=1
      rank=1.0d0

c-----p is the final and q is the initial state----------------

      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      if(jfin.ne.orbp)go to 100
      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(jinit.ne.orbq)go to 1
      if(iiq.ne.ival)go to 1

      if(ipsym.ne.iqsym) go to 1

c---------------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c---------------------------------------------------------------

      sum=sum+d(iip,iiq)

      call hyper(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dirac)
      sum1=sum1+d(iip,iiq)


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(iir.ne.iiq)then
      if(irsym.eq.ipsym)then
      if(irsym.eq.iqsym.and.orbq.eq.orbr)then

      in=jdpa(iiq,iir)

c------------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c------------------------------------------------------------

      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)

      endif
      endif
      endif

      if(irsym.eq.iqsym)then
      if(irsym.eq.ipsym.and.orbr.eq.orbp)then
      if(iir.ne.iip)then

      in=kdpa(iip,iir)

c-----------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-----------------------------------------------------------

      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)

      endif
      endif
      endif

 2    continue
      
      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact

  6   continue

      endif

c------------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c------------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact

  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      jn=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(jn)*fact
      sum6=sum6+d(iir,iig)*sf(jn)*fact

  9   continue

      endif

c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      jn=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(jn)*fact
      sum7=sum7-d(iir,iig)*sf(jn)*fact

 10   continue
  8   continue
  4   continue
  3   continue

      do 11 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c-------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c-------------------------------------------------------------

      if(irsym.ne.issym)go to 11
      if(orbp.ne.orbr.and.ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs.and.iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)

 11   continue


      do 12 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 12

      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(iqssym.ne.irgsym)go to 14

c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact

 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)     
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact

 18   continue
 17   continue
 14   continue

      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.and.iqsym.ne.issym)go to 20

      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c----------------------------------------------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact

 22   continue
 21   continue

c------------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)        
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.ne.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

  42  continue

c-----------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c-----------------------------------------------------------------

      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact

 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact

 48   continue
 47   continue
 46   continue

c----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact

 51   continue
 50   continue
 49   continue
 41   continue
 40   continue


      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)

      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52
      if(igsym.ne.idsym)go to 52

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)

c-----------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c-----------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      irdsym=mtbl(irsym,idsym)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbp+orbq+irank+2*orbg)
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact

  54  continue

c------------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c------------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact

  55  continue
  53  continue
  52  continue
  30  continue

c----------------------------------------------------------------
c     Two-body effective operator (computed directly)
c----------------------------------------------------------------

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)


      hyp1=0

      if(ipsym.eq.issym)then
      call hyper(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif


      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      if(iqssym.eq.irasym)then

c-----------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c-----------------------------------------------------------------


      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif

      hyp2=0

      if(issym.eq.iqsym)then
      call hyper(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      if(ipssym.eq.irasym)then

c----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif      
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      if(irasym.eq.iqssym)then

      hyp1=0

      if(itsym.eq.issym)then
      call hyper(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif

      endif
      endif
      endif
      endif
      endif


      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      irpsym=mtbl(irsym,ipsym)
      itasym=mtbl(itsym,iasym)
      if(irpsym.eq.itasym)then

      hyp1=0

      if(itsym.eq.issym)then
      call hyper(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif
    
      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)

c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------


      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif

      endif
      endif
      endif
      endif
      endif

      enddo
      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

      hyp1=0

      if(ipsym.eq.ibsym)then
      call hyper(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.iqssym)then
      if(irasym.eq.isbsym)then
      if(orbq.eq.orbb)then

c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c-----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      hyp2=0

      if(iqsym.eq.ibsym)then
      call hyper(orbb,orbq,iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.ipssym)then
      if(irasym.eq.isbsym)then
      if(orbp.eq.orbb)then

c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      hyp=0

      if(iasym.ne.igsym)go to 610

      call hyper(orbg,orba,iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
c     if(orbb.ne.orbq)go to 620

      in=idpa(iip,iib)


      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 27 = <rs|t2|ab><a|d|g><b|t1|p><sg|s2|qr>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+irank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact

 640  continue
      endif

c-------------------------------------------------------------------
c     t2dt1s2     diagram 28 = <rs|t2|ba><a|d|g><b|t1|p><sg|s2|qr>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact

  650 continue
      endif

  630 continue
  620 continue

      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
c     if(orbb.ne.orbp)go to 660

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2dt1s2     diagram 29 = <pa|s2|sr><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+irank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact

  680 continue
      endif

c-------------------------------------------------------------------
c     s2t1dt2     diagram 30 = <pa|s2|rs><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact

  690 continue
      endif

  670 continue
  660 continue
  610 continue

      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.ne.isgsym)go to 605
      if(irasym.ne.iqssym)go to 605

      if(ibsym.ne.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
      if(orbp.ne.orbb)go to 605

      in=idpa(iip,iib)

      hyp=0

      call hyper(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700

c-------------------------------------------------------------------
c     t2dt1s2     diagram 31 = <rs|t2|ag><g|d|b><p|t1|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

c-------------------------------------------------------------------
c     t2dt1s2     diagram 32 = <rs|t2|ga><g|d|b><p|t1|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact

  720 continue
  710 continue
  605 continue
  700 continue

      hyp=0

      if(ibsym.eq.igsym)then
      call hyper(orbg,orbb,iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)

      if(ipssym.ne.irasym)go to 600
      if(isgsym.ne.irasym)go to 600

      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 33 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ag>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

c-------------------------------------------------------------------
c     s2t1dt2     diagram 34 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbp,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbp+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact

  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0

      if(itsym.eq.issym)then
      call hyper(orbs,orbt,iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      isbsym=mtbl(issym,ibsym)
      irasym=mtbl(irsym,iasym)
      iqtsym=mtbl(iqsym,itsym)

      if(isbsym.ne.irasym)go to 815
      if(iqtsym.ne.irasym)go to 815


      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 35 = <rs|t2|ab><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)

      call findk(orbs,orbb,orbr,orba,ias,iab,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 812
      do 811 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      if(l2.eq.l1)then
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+irank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

      endif

  811 continue
  812 continue

c-------------------------------------------------------------------
c     t2dt1s2     diagram 45 = <rs|t2|ba><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbb,ias,iaa,iar,iab,kgot,kmax)
      if(kmax.eq.0)go to 814
      do 813 kloop=1,kmax
      l1=kgot(kloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrb=idra(iir,iib,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrb,isasym,l1)
      call sixj(orbp,rank,orbq,orbt,al2,orbs,s6j1)
      call sixj(orbp,al2,orbs,orba,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbt+orbs+l2+orbr+irank+orba)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum45=sum45+hyp*tf(inf)*si(ini)*ti(in)*fact

c     print*, 'OK',hyp,tf(inf),si(ini),ti(ini),fact, sum45

  813 continue
  814 continue
  830 continue
  820 continue
  815 continue

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      itbsym=mtbl(itsym,ibsym)

      if(ipssym.ne.irasym)go to 900
      if(itbsym.ne.irasym)go to 900

      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 36 = <pr|s2|sa><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)

      call findk(orbp,orbs,orbr,orba,iap,ias,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 912
      do 911 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      if(l2.eq.l1)then
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+irank+orba+orbr+al2)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact

      endif

  911 continue
  912 continue

c-------------------------------------------------------------------
c     s2t1dt2     diagram 46 = <ps|s2|ra><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------

      call findk(orbp,orbr,orbs,orba,iap,iar,ias,iaa,kgot,kmax)
      if(kmax.eq.0)go to 914
      do 913 kloop=1,kmax
      l1=kgot(kloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      nsa=idra(iis,iia,l1)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l1)
      call sixj(orbp,rank,orbq,orbt,al2,orbs,s6j1)
      call sixj(orbp,al2,orbs,orba,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbt+orbs+orbr+orba+l2+irank)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum46=sum46+hyp*ti(ini)*sf(inf)*tf(in)*fact

  913 continue
  914 continue
  930 continue
  900 continue
  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)

      hyp=0

      if(itsym.eq.ibsym)then
      call hyper(orbt,orbb,iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      isqsym=mtbl(issym,iqsym)

      if(irasym.eq.isbsym)then
      if(irasym.eq.isqsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

c-------------------------------------------------------------------
c     s1t2ds2     diagram 41 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nqs=jdra(iiq,iis,l1)
      nra=idra(iir,iia,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t2ds2     diagram 42 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|qa>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbq,orbr,iaa,ias,iaq,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif

      enddo
      endif

      endif
      endif

      endif
      endif
      endif
      endif


      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then

      in=jdpa(iiq,iit)

      hyp=0

      if(itsym.eq.ibsym)then
      call hyper(orbb,orbt,iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      ispsym=mtbl(issym,ipsym)

      if(irasym.eq.isbsym)then
      if(irasym.eq.ispsym)then

c--------------------------------------------------------------------
c     s2dt1s1     diagram 43 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ba>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

c-------------------------------------------------------------------
c     s2dt1s1     diagram 44 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ab>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue
   1  continue
 100  continue

      sum1=sum1*con/jinit
      sum2=sum2*con/jinit
      sum3=sum3*con/jinit
      sum4=sum4*con/jinit
      sum5=sum5*con/jinit
      sum6=sum6*con/jinit
      sum7=sum7*con/jinit
      sum8=sum8*con/jinit
      sum9=sum9*con/jinit
      sum10=sum10*con/jinit
      sum11=sum11*con/jinit
      sum12=sum12*con/jinit
      sum13=sum13*con/jinit
      sum14=sum14*con/jinit
      sum15=sum15*con/jinit
      sum16=sum16*con/jinit
      sum17=sum17*con/jinit
      sum18=sum18*con/jinit
      sum19=sum19*con/jinit
      sum20=sum20*con/jinit
      sum21=sum21*con/jinit
      sum22=sum22*con/jinit
      sum23=sum23*con/jinit
      sum24=sum24*con/jinit
      sum25=sum25*con/jinit
      sum26=sum26*con/jinit
      sum27=sum27*con/jinit
      sum28=sum28*con/jinit
      sum29=sum29*con/jinit
      sum30=sum30*con/jinit
      sum31=sum31*con/jinit
      sum32=sum32*con/jinit
      sum33=sum33*con/jinit
      sum34=sum34*con/jinit
      sum35=sum35*con/jinit
      sum36=sum36*con/jinit
      sum37=sum37*con/jinit
      sum38=sum38*con/jinit
      sum39=sum39*con/jinit
      sum40=sum40*con/jinit
      sum41=sum41*con/jinit
      sum42=sum42*con/jinit
      sum43=sum43*con/jinit
      sum44=sum44*con/jinit
      sum45=sum45*con/jinit
      sum46=sum46*con/jinit

      print*,'Dirac - Fock value hyp(A) ',dirac*con/jinit

      return
      end

c*********************************************************


c*********************************************************
c  This subroutine calculates reduced matrix element of  *
c  the hyperfine operator                                *
c*********************************************************

      subroutine hyper(j1,j2,i1,i2,m1,m2,kap1,kap2,val)

      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2

      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
       common/radial/hyp(MNBAS,MNBAS)

       factor=(-1)**(j1+0.5)
       w3j1=dr(j1,one,j2,half,zero,-half)
       e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
       result=hyp(i1,i2)
       val=-result*e1*(kap1+kap2)

       return
       end


c**********************************************************


c*********************************************************
c  This subroutine calculates the radial integral of the *
c  reduced matrix element of the hyperfine operator      *
c*********************************************************

      subroutine Calch(Nbasis)

      implicit double precision(a-h,o-z)

      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/radial/hyp(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

       mtp=n

       do i1=1,nbasis
       do i2=1,i1

       hyp(i1,i2) = zero

c       if(iparity(i1).eq.iparity(i2))then

       ta(1)=zero

       do 1 l=2,mtp
       ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)/(r(l)*r(l))
    1  continue
       call quad(result)

       hyp(i1,i2)=result
       hyp(i2,i1)=result

c       endif

       enddo
       enddo

       return
       end

c************************************************************

      subroutine cinit(d)
      real*8 d(MNBAS,MNBAS)
      do 1 i=1,MNBAS
      do 1 j=1,MNBAS
 1    d(i,j)=0.0d0
      return
      end

c************************************************************
c****************************************************************
c This subroitine calculates the nms isotopic shift property diagrams
c****************************************************************

      subroutine nmsp(d,si,sf,ti,tf,
     :             ival,jval,gi,jinit,jfin)


      implicit real*8 (a-h,o-z)
      real*8 jinit,jfin

      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      common/value/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,sum45,sum46
      
      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)
      
      data zero,half,one,two/0.0,0.5,1.0,2.0/
      data VL/137.03598950D+00/
      data EMF/1836.152701D+00/
      data AUEV/27.2116D+00/
      data EVMHZ/2.417988D+08/
      
      BETAE=ONE/(TWO*VL)
      BETAN=BETAE/EMF
      CONVF=AUEV*EVMHZ
      CON=GI*BETAN*CONVF
      write(STDOUT,*),'sourav333'
      
      print*,'CON',con
      
      d1=dr(jinit,1.0d0,jfin,-jinit,0.0d0,jfin)
      con=con*d1
      
      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0
      sum45=0
      sum46=0
      
      irank=1
      rank=1.0d0
      
c-----p is the final and q is the initial state----------------
      
      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)
      
      if(jfin.ne.orbp)go to 100
      if(iip.ne.jval)go to 100
      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(jinit.ne.orbq)go to 1
      if(iiq.ne.ival)go to 1

      if(ipsym.ne.iqsym) go to 1

c---------------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c---------------------------------------------------------------

      sum=sum+d(iip,iiq)

      call nms(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dirac)
      sum1=sum1+d(iip,iiq)

      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      if(iir.ne.iiq)then
      if(irsym.eq.ipsym)then
      if(irsym.eq.iqsym.and.orbq.eq.orbr)then
      
      in=jdpa(iiq,iir)
      
c------------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c------------------------------------------------------------
      
      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)

      endif
      endif  
      endif

      if(irsym.eq.iqsym)then
      if(irsym.eq.ipsym.and.orbr.eq.orbp)then
      if(iir.ne.iip)then
      
      in=kdpa(iip,iir)

c-----------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-----------------------------------------------------------
      
      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)
      endif
      endif
      endif

 2    continue

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact
      
  6   continue
      
      endif
      
c------------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c------------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact

  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      jn=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(jn)*fact
      sum6=sum6+d(iir,iig)*sf(jn)*fact

  9   continue

      endif
c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------
      
      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      jn=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(jn)*fact
      sum7=sum7-d(iir,iig)*sf(jn)*fact
      
 10   continue
  8   continue
  4   continue
  3   continue
      
      do 11 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c-------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c-------------------------------------------------------------

      if(irsym.ne.issym)go to 11
      if(orbp.ne.orbr.and.ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs.and.iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)

 11   continue


      do 12 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      
      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 12
      
      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      
      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)
      
      if(iqssym.ne.irgsym)go to 14
      
c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact

 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact
      
 18   continue
 17   continue
 14   continue
      
      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.and.iqsym.ne.issym)go to 20
      
      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c----------------------------------------------------------------
      
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact

 22   continue
 21   continue

c------------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.ne.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

  42  continue  

c-----------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c-----------------------------------------------------------------
      
      itgsym=mtbl(itsym,igsym)
      
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      
      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact

 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact

 48   continue  
 47   continue
 46   continue
      
c----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c----------------------------------------------------------------
      
      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)
      
      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact

 51   continue
 50   continue
 49   continue
 41   continue
 40   continue
      do 52 id=1,nocc 
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)
      
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)
      
      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52
      if(igsym.ne.idsym)go to 52
      
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)
 
c-----------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c-----------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      irdsym=mtbl(irsym,idsym)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbp+orbq+irank+2*orbg)
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact

  54  continue

c------------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c------------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact
      
  55  continue
  53  continue
  52  continue
  30  continue
      
c----------------------------------------------------------------
c     Two-body effective operator (computed directly)
c----------------------------------------------------------------
      
      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir) 
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)
      
      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)


      hyp1=0

      if(ipsym.eq.issym)then
      call nms(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif


      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      if(iqssym.eq.irasym)then

c-----------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c-----------------------------------------------------------------


      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact
      
      endif
      enddo
      endif
      
c-----------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c-----------------------------------------------------------------
      
      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1) 
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact
      enddo
      endif
      endif
      endif

      hyp2=0

      if(issym.eq.iqsym)then
      call nms(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      if(ipssym.eq.irasym)then

c----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif
      
c-----------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact
      
      enddo
      endif      
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      if(irasym.eq.iqssym)then

      hyp1=0

      if(itsym.eq.issym)then
      call nms(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------
      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-----------------------------------------------------------------
      
      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif

      endif
      endif
      endif
      endif
      endif


      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      irpsym=mtbl(irsym,ipsym)
      itasym=mtbl(itsym,iasym)
      if(irpsym.eq.itasym)then

      hyp1=0

      if(itsym.eq.issym)then
      call nms(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)
      
c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------
      

      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact
      
      endif
      enddo
      endif
      
c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------
      
      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif

      endif
      endif
      endif
      endif
      endif

      enddo
      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)
      
      hyp1=0
      
      if(ipsym.eq.ibsym)then
      call nms(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif
      
      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.iqssym)then
      if(irasym.eq.isbsym)then
      if(orbq.eq.orbb)then
      
c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------
      
      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c-----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      hyp2=0
      
      if(iqsym.eq.ibsym)then
      call nms(orbb,orbq,iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif     

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      
      if(irasym.eq.ipssym)then
      if(irasym.eq.isbsym)then
      if(orbp.eq.orbb)then
      
c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------
      
      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      hyp=0
      
      if(iasym.ne.igsym)go to 610
      
      call nms(orbg,orba,iig,iia,orbg,orba,kapg,kapa,hyp)
      
      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)
      
      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
c     if(orbb.ne.orbq)go to 620
      
      in=idpa(iip,iib)

      
      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 27 = <rs|t2|ab><a|d|g><b|t1|p><sg|s2|qr>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+irank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact

 640  continue
      endif

c-------------------------------------------------------------------
c     t2dt1s2     diagram 28 = <rs|t2|ba><a|d|g><b|t1|p><sg|s2|qr>
c-------------------------------------------------------------------
      
      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact
      
  650 continue
      endif
      
  630 continue
  620 continue
      
      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)
      
      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
c     if(orbb.ne.orbp)go to 660

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2dt1s2     diagram 29 = <pa|s2|sr><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+irank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact

  680 continue
      endif

c-------------------------------------------------------------------
c     s2t1dt2     diagram 30 = <pa|s2|rs><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------
      
      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact
      
  690 continue
      endif
      
  670 continue
  660 continue
  610 continue
      
      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.ne.isgsym)go to 605
      if(irasym.ne.iqssym)go to 605

      if(ibsym.ne.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
      if(orbp.ne.orbb)go to 605

      in=idpa(iip,iib)

      hyp=0

      call nms(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700

c-------------------------------------------------------------------
c     t2dt1s2     diagram 31 = <rs|t2|ag><g|d|b><p|t1|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

c-------------------------------------------------------------------
c     t2dt1s2     diagram 32 = <rs|t2|ga><g|d|b><p|t1|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact
      
  720 continue
  710 continue
  605 continue
  700 continue
      
      hyp=0

      if(ibsym.eq.igsym)then
      call nms(orbg,orbb,iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)

      if(ipssym.ne.irasym)go to 600
      if(isgsym.ne.irasym)go to 600

      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 33 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ag>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact
      
c-------------------------------------------------------------------
c     s2t1dt2     diagram 34 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ga>
c-------------------------------------------------------------------
      
      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbp,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbp+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact
      
  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0

      if(itsym.eq.issym)then
      call nms(orbs,orbt,iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      isbsym=mtbl(issym,ibsym)
      irasym=mtbl(irsym,iasym)
      iqtsym=mtbl(iqsym,itsym)

      if(isbsym.ne.irasym)go to 815
      if(iqtsym.ne.irasym)go to 815


      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 35 = <rs|t2|ab><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------
  
      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)

      call findk(orbs,orbb,orbr,orba,ias,iab,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 812
      do 811 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      if(l2.eq.l1)then
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+irank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact
      
      endif
      
  811 continue
  812 continue

c-------------------------------------------------------------------
c     t2dt1s2     diagram 45 = <rs|t2|ba><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbb,ias,iaa,iar,iab,kgot,kmax)
      if(kmax.eq.0)go to 814
      do 813 kloop=1,kmax
      l1=kgot(kloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrb=idra(iir,iib,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrb,isasym,l1)
      call sixj(orbp,rank,orbq,orbt,al2,orbs,s6j1)
      call sixj(orbp,al2,orbs,orba,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbt+orbs+l2+orbr+irank+orba)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum45=sum45+hyp*tf(inf)*si(ini)*ti(in)*fact

c     print*, 'OK',hyp,tf(inf),si(ini),ti(ini),fact, sum45

  813 continue
  814 continue
  830 continue
  820 continue
  815 continue

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      itbsym=mtbl(itsym,ibsym)

      if(ipssym.ne.irasym)go to 900
      if(itbsym.ne.irasym)go to 900

      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900
      
      in=idpa(iiq,iib)
      
c-------------------------------------------------------------------
c     s2t1dt2     diagram 36 = <pr|s2|sa><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------
      
      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)
  
      call findk(orbp,orbs,orbr,orba,iap,ias,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 912
      do 911 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      if(l2.eq.l1)then
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+irank+orba+orbr+al2)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact

      endif

  911 continue
  912 continue

c-------------------------------------------------------------------
c     s2t1dt2     diagram 46 = <ps|s2|ra><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------

      call findk(orbp,orbr,orbs,orba,iap,iar,ias,iaa,kgot,kmax)
      if(kmax.eq.0)go to 914
      do 913 kloop=1,kmax
      l1=kgot(kloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      nsa=idra(iis,iia,l1)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l1)
      call sixj(orbp,rank,orbq,orbt,al2,orbs,s6j1)
      call sixj(orbp,al2,orbs,orba,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbt+orbs+orbr+orba+l2+irank)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum46=sum46+hyp*ti(ini)*sf(inf)*tf(in)*fact

  913 continue
  914 continue
  930 continue
  900 continue
  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)

      hyp=0

      if(itsym.eq.ibsym)then
      call nms(orbt,orbb,iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif
      
      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      isqsym=mtbl(issym,iqsym)
      
      if(irasym.eq.isbsym)then
      if(irasym.eq.isqsym)then
      
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

c-------------------------------------------------------------------
c     s1t2ds2     diagram 41 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nqs=jdra(iiq,iis,l1)
      nra=idra(iir,iia,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t2ds2     diagram 42 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|qa>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbq,orbr,iaa,ias,iaq,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym) 
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact
      
      enddo
      endif
      
      enddo
      endif
      
      endif
      endif
      
      endif
      endif
      endif
      endif


      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then

      in=jdpa(iiq,iit)

      hyp=0

      if(itsym.eq.ibsym)then
      call nms(orbb,orbt,iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      ispsym=mtbl(issym,ipsym)

      if(irasym.eq.isbsym)then
      if(irasym.eq.ispsym)then

c--------------------------------------------------------------------
c     s2dt1s1     diagram 43 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ba>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact
      
c-------------------------------------------------------------------
c     s2dt1s1     diagram 44 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ab>
c-------------------------------------------------------------------
      
      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact
      
      enddo
      endif
      enddo
      endif
      
      endif
      endif

      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue
   1  continue
 100  continue

      sum1=sum1*con/jinit
      sum2=sum2*con/jinit
      sum3=sum3*con/jinit
      sum4=sum4*con/jinit
      sum5=sum5*con/jinit
      sum6=sum6*con/jinit
      sum7=sum7*con/jinit
      sum8=sum8*con/jinit
      sum9=sum9*con/jinit
      sum10=sum10*con/jinit
      sum11=sum11*con/jinit
      sum12=sum12*con/jinit
      sum13=sum13*con/jinit
      sum14=sum14*con/jinit
      sum15=sum15*con/jinit
      sum16=sum16*con/jinit
      sum17=sum17*con/jinit
      sum18=sum18*con/jinit
      sum19=sum19*con/jinit
      sum20=sum20*con/jinit
      sum21=sum21*con/jinit
      sum22=sum22*con/jinit
      sum23=sum23*con/jinit
      sum24=sum24*con/jinit
      sum25=sum25*con/jinit
      sum26=sum26*con/jinit
      sum27=sum27*con/jinit
      sum28=sum28*con/jinit
      sum29=sum29*con/jinit
      sum30=sum30*con/jinit
      sum31=sum31*con/jinit
      sum32=sum32*con/jinit
      sum33=sum33*con/jinit
      sum34=sum34*con/jinit
      sum35=sum35*con/jinit
      sum36=sum36*con/jinit
      sum37=sum37*con/jinit
      sum38=sum38*con/jinit
      sum39=sum39*con/jinit
      sum40=sum40*con/jinit
      sum41=sum41*con/jinit
      sum42=sum42*con/jinit
      sum43=sum43*con/jinit
      sum44=sum44*con/jinit
      sum45=sum45*con/jinit
      sum46=sum46*con/jinit
      
      print*,'Dirac - Fock value hyp(A) ',dirac*con/jinit

      return
      end

c*********************************************************


c*********************************************************
c  This subroutine calculates reduced matrix element of  *
c  the isotopic shift operator                                *
c*********************************************************

      subroutine nms(j1,j2,i1,i2,m1,m2,kap1,kap2,val)

      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2

      dimension nmsr(MNBAS,MNBAS)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
cc     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
cc       common/radial/nmsr(MNBAS,MNBAS)
       
c       write(STDOUT,*)'sourav444'
       call Calcnms(i1,kap1,nmsr)
       factor=(two*j2+one)
       w3j1=dr(j1,one,j2,half,zero,-half)
       e1=factor*(w3j1**2)
       result=nmsr(i1,i2)
       val=result*e1
       return
       end
      

c**********************************************************


c*********************************************************
c  This subroutine calculates the radial integral of the *
c  reduced matrix element of the hyperfine operator      *
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

c       write(STDOUT,*)'sourav555'       
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


      subroutine octCprop(d,si,sf,ti,tf,
     :             ival,jval,ocm,jinit,jfin)

      implicit real*8 (a-h,o-z)
      real*8 jinit,jfin
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)

      common/kappa/nka(MNBAS)

      common/valueC/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two,three/0.0,0.5,1.0,2.0,3.0/
      data VL/137.03598950D+00/
      data EMF/1836.152701D+00/
      data AUEV/27.2116D+00/
      data EVMHZ/2.417988D+08/

      BETAE=ONE/(TWO*VL)
      BETAN=BETAE/EMF
      CONVF=AUEV*EVMHZ
      CONC=OCM*CONVF*BETAN*3.5711D-08

      print*,'CONC',conc

      d11=dr(jinit,3.0d0,jfin,-jinit,0.0d0,jfin)
      con=conc*d11

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0

      irank=3
      rank=3.0d0

c     p is the final and q is the initial state
      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)
      if(jfin.ne.orbp)go to 100
      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
      if(jinit.ne.orbq)go to 1
      if(iiq.ne.ival)go to 1
      if(ipsym.ne.iqsym) go to 1
c     diagram 1 =<f|d|i>
      sum=sum+d(iip,iiq)
      call hyperc(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dirac)
c      sum1=sum1+d(iip,iiq)
      sum1=sum1+dirac
c     diagram < | ds|>

      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      if(iir.ne.ival)then
      if(iir.ne.iiq)then
      if(irsym.eq.ipsym)then
      if(irsym.eq.iqsym.and.orbq.eq.orbr)then
      in=jdpa(iiq,iir)
c     diagram 2(a)=  sum (r) <p|d|r><r|s|q>
      call hyperc(orbp,orbr,iip,iir,orbp,orbr,kapp,kapr,oct)
      sum=sum+d(iip,iir)*si(in)
c      sum2=sum2+d(iip,iir)*si(in)
      sum2=sum2+oct*si(in)
      endif
      endif
      endif

      if(irsym.eq.iqsym)then
      if(irsym.eq.ipsym.and.orbr.eq.orbp)then
      if(iir.ne.iip)then
      in=kdpa(iip,iir)
c     diagram 3(a) =  sum (r) <p|s|r><r|d|q>
      call hyperc(orbr,orbq,iir,iiq,orbr,orbq,kapr,kapq,oct)
      sum=sum+d(iir,iiq)*sf(in)
c      sum3=sum3+d(iir,iiq)*sf(in)
      sum3=sum3+oct*sf(in)
      endif
      endif
      endif
      endif
 2    continue

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 4
      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)
      if(ipqsym.ne.irgsym)go to 4

c     diagram 2 (b) = sum(g,r) <g|d|r><pr|s|qg>
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+al1)/(2*al1+1)
      call hyperc(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,oct)
      sum=sum+d(iig,iir)*si(in)*fact
c      sum4=sum4+d(iig,iir)*si(in)*fact
      sum4=sum4+oct*si(in)*fact
  6   continue
      endif

c     diagram 2 (c) =-sum(g,r) <g|d|r><pr|s|gq>
      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(three,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      call hyperc(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,oct)
      sum=sum-d(iig,iir)*si(in)*fact
c      sum5=sum5-d(iig,iir)*si(in)*fact
      sum5=sum5-oct*si(in)*fact
  7   continue
  5   continue

c     diagram 3 (b)
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      jn=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+al1)/(2*al1+1)
      call hyperc(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,oct)
      sum=sum+d(iir,iig)*sf(jn)*fact
c      sum6=sum6+d(iir,iig)*sf(jn)*fact
      sum6=sum6+oct*sf(jn)*fact
  9   continue
      endif

c     diagram 3 (c)
      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      jn=nslocf(npr,nqg,iqgsym,l1)
      call sixj(three,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      call hyperc(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,oct)
      sum=sum-d(iir,iig)*sf(jn)*fact
c      sum7=sum7-d(iir,iig)*sf(jn)*fact
      sum7=sum7-oct*sf(jn)*fact
 10   continue
  8   continue
  4   continue
  3   continue

c     diagram 4e
      do 11 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      if(orbr.ne.orbs)go to 11
      if(orbp.ne.orbr.or.ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs.or.iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11
      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)
      call hyperc(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,oct)
      sum=sum+d(iir,iis)*si(ini)*sf(inf)
c      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+oct*si(ini)*sf(inf)
 11   continue

c     diagram 4a
      do 12 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)
      if(irsym.ne.igsym)go to 12
      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14
      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)
      if(iqssym.ne.irgsym)go to 14
      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      call hyperc(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,oct)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
c      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+oct*si(ini)*sf(inf)*fact
 16   continue
 15   continue
c     diagram 4b
      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      call hyperc(orbg,orbr,iig,iir,orbg,orbr,kapg,kapr,oct)
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
c      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-oct*si(ini)*sf(inf)*fact
 18   continue
 17   continue
 14   continue
      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.or.iqsym.ne.issym)go to 20
      ini=jdpa(iiq,iis)
      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)
      if(ipssym.ne.irgsym)go to 20
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      call hyperc(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,oct)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
c      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+oct*si(ini)*sf(inf)*fact
 22   continue
 21   continue
c     diagram 4d
      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      call hyperc(orbr,orbg,iir,iig,orbr,orbg,kapr,kapg,oct)
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
c      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-oct*si(ini)*sf(inf)*fact
 24   continue
 23   continue
 20   continue
 13   continue
 12   continue

          goto 999
      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.ne.itsym)go to 40
      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)
      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41
      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
c     diagram 4f
      call hyperc(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,oct)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
c      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+oct*si(ini)*sf(inf)*fact
      endif
      enddo
      endif

  42  continue
c     diagram 4(g)
      itgsym=mtbl(itsym,igsym)
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      call hyperc(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,oct)
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
c      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-oct*si(ini)*sf(inf)*fact
 45   continue
 44   continue
 43   continue
c     diagram 4(h)
      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      call hyperc(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,oct)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
c      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+oct*si(ini)*sf(inf)*fact
 48   continue
 47   continue
 46   continue
c     diagram 4(i)
      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)
      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      call hyperc(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,oct)
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
c      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-oct*si(ini)*sf(inf)*fact
 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

c     diagram 4(j)
      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)
      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52

      if(igsym.ne.idsym)go to 52
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)
      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbp+orbq+irank+2*orbg)
      call hyperc(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,oct)
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
c      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-oct*si(ini)*sf(inf)*fact
  54  continue
      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      call hyperc(orbd,orbg,iid,iig,orbd,orbg,kapd,kapg,oct)
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
c      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+oct*si(ini)*sf(inf)*fact
  55  continue
  53  continue
  52  continue
  30  continue

C     Two-body effective operator (computed directly)

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)
      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)
      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

c     diagrams  a1-a4 (sum19-sum22)

      hyp1=0
      if(ipsym.eq.issym)then
      call hyperc(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif


      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact
      enddo
      endif

      endif

      hyp2=0
      if(issym.eq.iqsym)then
      call hyperc(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)
      ipssym=mtbl(ipsym,issym)

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact
      enddo
      endif

      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

c     diagram 43
      hyp1=0
      if(itsym.eq.issym)then
      call hyperc(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 45
      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0
      if(itsym.eq.issym)then
      call hyperc(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)
      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact
      endif
      enddo
      endif


      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact
      enddo
      endif


      endif
      endif
      endif
      endif

      enddo

      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

c     diagrams  b1-b4 (sum19-sum22)
      hyp1=0
      if(ipsym.eq.ibsym)then
      call hyperc(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.eq.iqssym)then
      if(orbq.eq.orbb)then

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)

      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)

      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      hyp2=0
      if(iqsym.eq.ibsym)then
      call hyperc(orbb,orbq,iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)
      if(irasym.eq.ipssym)then
      if(orbp.eq.orbb)then

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.ne.igsym)go to 610
      hyp=0
c     c1-c2
      call hyperc(orbg,orba,iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620
      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact
 640  continue
      endif

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact
  650 continue
      endif

  630 continue
  620 continue

c     d1-d2
      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660
      in=idpa(iiq,iib)

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact
  680 continue
      endif

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact
  690 continue
      endif

  670 continue
  660 continue
  610 continue


      if(ibsym.ne.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
c     if(orbp.ne.orbg)go to 605
      if(orbp.ne.orbb)go to 605
      in=idpa(iip,iib)
c     c7-c8
      hyp=0
      call hyperc(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700
      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)


      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)

      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)


      if(irasym.eq.isgsym)then

      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact
      endif


      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)

      isasym=mtbl(issym,iasym)
      irgsym=mtbl(irsym,igsym)

      inf=ntloc(nsa,nrg,isasym,l1)

      if(irgsym.eq.isasym)then

      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact
      endif

  720 continue
  710 continue
  700 continue
  605 continue
c     d7-d8

      hyp=0
      if(ibsym.eq.igsym)then
      call hyperc(orbg,orbb,iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600
      in=idpa(iiq,iib)
      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact
  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0
      if(itsym.eq.issym)then
      call hyperc(orbs,orbt,iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)
c     c3
      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

  830 continue
  820 continue
  815 continue

c     d3
      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900
      in=idpa(iiq,iib)

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact
  930 continue
  900 continue

  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)
      hyp=0
      if(itsym.eq.ibsym)then
      call hyperc(orbt,orbb,iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nqs=jdra(iiq,iis,l1)
      nra=idra(iir,iia,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))

      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact
      endif
      enddo
      endif

      call findk(orba,orbs,orbq,orbr,iaa,ias,iaq,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)

      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif
      enddo
      endif

      endif

      endif
      endif
      endif
      endif


      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then
      in=jdpa(iiq,iit)

      hyp=0
      if(itsym.eq.ibsym)then
      call hyperc(orbb,orbt,iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact
      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue

  999 continue
   1  continue
 100  continue

      sum1=sum1*con
      sum2=sum2*con
      sum3=sum3*con
      sum4=sum4*con
      sum5=sum5*con
      sum6=sum6*con
      sum7=sum7*con
      sum8=sum8*con
      sum9=sum9*con
      sum10=sum10*con
      sum11=sum11*con
      sum12=sum12*con
      sum13=sum13*con
      sum14=sum14*con
      sum15=sum15*con
      sum16=sum16*con
      sum17=sum17*con
      sum18=sum18*con
      sum19=sum19*con
      sum20=sum20*con
      sum21=sum21*con
      sum22=sum22*con
      sum23=sum23*con
      sum24=sum24*con
      sum25=sum25*con
      sum26=sum26*con
      sum27=sum27*con
      sum28=sum28*con
      sum29=sum29*con
      sum30=sum30*con
      sum31=sum31*con
      sum32=sum32*con
      sum33=sum33*con
      sum34=sum34*con
      sum35=sum35*con
      sum36=sum36*con
      sum37=sum37*con
      sum38=sum38*con
      sum39=sum39*con
      sum40=sum40*con
      sum41=sum41*con
      sum42=sum42*con
      sum43=sum43*con
      sum44=sum44*con
      print*,'Dirac - Fock value hyp(C) ',dirac*con
      return
      end

c*********************************************************
  

      subroutine hyperc(j1,j2,i1,i2,m1,m2,kap1,kap2,val)
c*********************************************************
c  this subrout calculates the edm t-pt reduced
c  matrix elements between orbitals with j=1/2
c*********************************************************
      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
       common/radialc/hyp(MNBAS,MNBAS)
       mtp=n
c       ta(1)=zero
c       do 1 l=2,mtp
c       ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)/(r(l)**4)
c    1  continue
c       call quad(result)

       factor=(-1)**(j1+0.5)
       w3j1=dr(j1,three,j2,half,zero,-half)
       e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
       result=hyp(i1,i2)
       val=-result*e1*(kap1+kap2)/three
       return
       end

c***********************************************************************
c              SUBROUTINE   SLOCI                                      *
c                                                                      *
c     This subroutine sets the equation  index for the open shell      *
c      cluster amplitudes                                              *
c                                                                      *
c***********************************************************************

      subroutine sloci(iopt,iorb,nsmax,nprint)

      implicit real*8 (a-h,o-z)
      common/iallparity/iall
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
      common/maxj/jmax
      common/skip1/npbps(0:MXV,2)
      common/skip2/mmskip(0:MXV,2),m2skip
      common/skip3/ncis

c     local dimension
      dimension kgot(MXV+1),ieven(0:MXV),iodd(0:MXV)
      dimension npqbps(0:MXV,2),ncheck(MDIM)
      dimension indx(MDIM),indx1(MDIM)
      data zero/0.0d0/

      nbasis=nocc+nexcit

      do i=1,MDIM
         ncheck(i)=0
      enddo


      do i=0,MXV
       kgot(i)=0
       ieven(i)=0
       iodd(i)=0
       npqbps(i,1)=0
       npqbps(i,2)=0
       mmskip(i,1)=0
       mmskip(i,2)=0
      enddo

c--------------------------Set S1 equation index--------------------------------

      ncis=0

      if(iopt.eq.1)then

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)

      do 2 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)

      if(iip.ne.iorb)go to 2
      if(iip.eq.iiq)go to 2
      if(orbp.ne.orbq)go to 2
      if(ipsym.ne.iqsym)go to 2

      ncis=ncis+1
      jdpa(iip,iiq)=ncis
      indx(ncis)=ncis
      indx1(ncis)=ncis

  2   continue
  1   continue

      endif

c----------------------Set the skip information--------------------------------

      do 3 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      do 4 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)

      if(iip.ne.iorb)go to 4

      ipqsym=mtbl(ipsym,iqsym)

      call findk(orbp,orbq,orbp,orbq,iap,iaq,iap,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 4
      do 5 kloop=1,kmax
      kk=kgot(kloop)
      if(ipqsym.eq.1)then
      ieven(kk)=ieven(kk)+1
      jdra(iip,iiq,kk)=ieven(kk)
      else
      iodd(kk)=iodd(kk)+1
      jdra(iip,iiq,kk)=iodd(kk)
      endif

  5   continue
  4   continue
  3   continue

      m2skip=0

      if(iall.eq.0)then

      do 6 i=2,MXV,2
      npqbps(i,2)=npqbps(i-2,2)+ieven(i-2)*npbps(i-2,2)
      mmskip(i,2)=npqbps(i,2)
  6   continue

      do  i=0,JMAX,2
      m2skip=m2skip+ieven(i)*npbps(i,2)
      enddo

      do 7 i=3,MXV,2
      npqbps(i,1)=npqbps(i-2,1)+iodd(i-2)*npbps(i-2,1)
      mmskip(i,1)=npqbps(i,1)
  7   continue

      else

      do 10 i=1,MXV
      npqbps(i,2)=npqbps(i-1,2)+ieven(i-1)*npbps(i-1,2)
      mmskip(i,2)=npqbps(i,2)
 10   continue

      do  i=0,MXV
      m2skip=m2skip+ieven(i)*npbps(i,2)
      enddo

      do 11 i=1,MXV
      npqbps(i,1)=npqbps(i-1,1)+iodd(i-1)*npbps(i-1,1)
      mmskip(i,1)=npqbps(i,1)
 11   continue

      endif
   
      m2skip=m2skip+ncis

c-----------------------Set S2 equation index----------------------------------

      nn=ncis

      do 15 ip=1,nexcit
      ipsym=isyme(ip)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)

      ipqsym=mtbl(ipsym,iqsym)

      do 30 ir=1,nexcit
      irsym=isyme(ir)

      ipqrsm=mtbl(ipqsym,irsym)

      do 40 ig=1,nocc
      igsym=isymc(ig)

      if(ipqrsm.ne.igsym)go to 40

      orbg=orbc(ig)
      orbp=orbe(ip)
      orbq=orbe(iq)
      orbr=orbe(ir)

      iap=iqe(ip)
      iag=iqc(ig)
      iaq=iqe(iq)
      iar=iqe(ir)

      iip=ke(ip)
      iiq=ke(iq)
      iir=ke(ir)
      iig=kc(ig)

      if(iip.ne.iorb)go to 40

      call findk(orbp,orbq,orbr,orbg,iap,iaq,iar,iag,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 50 nloop=1,kmax
      kk=kgot(nloop)
      npq=jdra(iip,iiq,kk)
      nrg=idra(iir,iig,kk)
      nn=nn+1
      neqn=nsloci(npq,nrg,ipqsym,kk)

      if (ncheck(neqn).eq.1) then
        print*,'double count in sloci =', neqn
      stop
      endif

      ncheck(neqn)=1

      indx(nn)=nn
      indx1(nn)=neqn

  50  continue
  40  continue
  30  continue
  20  continue
  15  continue

      nsmax=nn

      do i=1,nn
      do j=1,nn-1
      if (indx1(j).gt.indx1(j+1))then
      ind=indx1(j)
      indx1(j)=indx1(j+1)
      indx1(j+1)=ind
      endif
      enddo
      enddo

      do i=1,nn
      ii=indx1(i)-indx(i)
      if(ii.ne.0)then
      write(*,*)indx(i),indx1(i),indx1(i)-indx(i)
      endif
      enddo
 
9999  format(//,'DOUBLE EXCITATIONS',/)
9998  format(3x,'IN','   IIP ',' IIQ ',' IIS ',' IIA',
     $1x,'JP',3X,'JQ',3x,'JR',3x,'JA',3x,'PSM',3x,'QSM',3x,
     $'RSM',3x,'ASM',3x,'KMAX',/)
9997  format(5i5,4f5.1,6i6)
9996  format(/,2x,'COUPLED CLUSTER A MATRIX (IP)',/)
9995  format(4f20.11)
9991  format(//,'SINGLE EXCITATIONS',/)
9990  format(3x,'IN','   IIP ',' IIQ',' J(IIP)',' J(IIQ)'/)
9989  format(3i5,2(f5.1,2x))

      return
      end

c******************************************************************************

      function nsloci(ipq,irg,ipar,kk)
      implicit real*8 (a-h,o-z)
c     common block
      common/skip1/npbps(0:MXV,2)
      common/skip2/mmskip(0:MXV,2),m2skip
      common/skip3/nsing

      if(ipar.eq.1)then
      nsloci=(ipq-1)*npbps(kk,2)+irg+mmskip(kk,2)+nsing
      else
      nsloci=(ipq-1)*npbps(kk,1)+irg+mmskip(kk,1)+m2skip
      endif
      return
      end


c***********************************************************************
c              SUBROUTINE   SLOCF                                      *
c                                                                      *
c     This subroutine sets the equation  index for the open shell      *
c      cluster amplitudes                                              *
c                                                                      *
c***********************************************************************

      subroutine slocf(iopt,iorb,nsmax,nprint,sj,sf)

      implicit real*8 (a-h,o-z)
      common/iallparity/iall
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
      common/idra3/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/jdpa(MNBAS,MNBAS)
      common/maxj/jmax
      common/skip1/npbps(0:MXV,2)
      common/skip6/mmskip(0:MXV,2),m2skip
      common/skip7/ncis

c     local dimension
      dimension kgot(MXV+1),ieven(0:MXV),iodd(0:MXV)
      dimension npqbps(0:MXV,2),ncheck(MDIM)
      dimension indx(MDIM),indx1(MDIM)
      dimension sf(MDIM),sj(MDIM)
      data zero/0.0d0/

      nbasis=nocc+nexcit
      do i=1,MDIM
         ncheck(i)=0
      enddo


      do i=0,MXV
       kgot(i)=0
       ieven(i)=0
       iodd(i)=0
       npqbps(i,1)=0
       npqbps(i,2)=0
       mmskip(i,1)=0
       mmskip(i,2)=0
      enddo

c--------------------------Set S1 equation index--------------------------------

      ncis=0

      if(iopt.eq.1)then

      do 1 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)

      do 2 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)

      if(iip.ne.iorb)go to 2
      if(iip.eq.iiq)go to 2
      if(orbp.ne.orbq)go to 2
      if(ipsym.ne.iqsym)go to 2

      ncis=ncis+1
      jdpa(iip,iiq)=ncis
      sf(ncis)=sj(ncis)
      indx(ncis)=ncis
      indx1(ncis)=ncis

  2   continue
  1   continue

      endif

c----------------------Set the skip information--------------------------------

      do 3 ip=1,nexcit
      ipsym=isyme(ip)
      iip=ke(ip)
      orbp=orbe(ip)
      iap=iqe(ip)

      do 4 iq=1,nexcit
      iqsym=isyme(iq)
      iiq=ke(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)

      if(iip.ne.iorb)go to 4

      ipqsym=mtbl(ipsym,iqsym)

      call findk(orbp,orbq,orbp,orbq,iap,iaq,iap,iaq,kgot,
     $kmax)
      if(kmax.eq.0)go to 4
      do 5 kloop=1,kmax
      kk=kgot(kloop)
      if(ipqsym.eq.1)then
      ieven(kk)=ieven(kk)+1
      jdra(iip,iiq,kk)=ieven(kk)
      else
      iodd(kk)=iodd(kk)+1
      jdra(iip,iiq,kk)=iodd(kk)
      endif

  5   continue
  4   continue
  3   continue

      m2skip=0

      if(iall.eq.0)then

      do 6 i=2,MXV,2
      npqbps(i,2)=npqbps(i-2,2)+ieven(i-2)*npbps(i-2,2)
      mmskip(i,2)=npqbps(i,2)
  6   continue

      do  i=0,JMAX,2
      m2skip=m2skip+ieven(i)*npbps(i,2)
      enddo

      do 7 i=3,MXV,2
      npqbps(i,1)=npqbps(i-2,1)+iodd(i-2)*npbps(i-2,1)
      mmskip(i,1)=npqbps(i,1)
  7   continue

      else

      do 10 i=1,MXV
      npqbps(i,2)=npqbps(i-1,2)+ieven(i-1)*npbps(i-1,2)
      mmskip(i,2)=npqbps(i,2)
 10   continue

      do  i=0,MXV
      m2skip=m2skip+ieven(i)*npbps(i,2)
      enddo

      do 11 i=1,MXV
      npqbps(i,1)=npqbps(i-1,1)+iodd(i-1)*npbps(i-1,1)
      mmskip(i,1)=npqbps(i,1)
 11   continue

      endif
 
      m2skip=m2skip+ncis

c-----------------------Set S2 equation index----------------------------------

      nn=ncis

      do 15 ip=1,nexcit
      ipsym=isyme(ip)

      do 20 iq=1,nexcit
      iqsym=isyme(iq)

      ipqsym=mtbl(ipsym,iqsym)

      do 30 ir=1,nexcit
      irsym=isyme(ir)

      ipqrsm=mtbl(ipqsym,irsym)

      do 40 ig=1,nocc
      igsym=isymc(ig)

      if(ipqrsm.ne.igsym)go to 40

      orbg=orbc(ig)
      orbp=orbe(ip)
      orbq=orbe(iq)
      orbr=orbe(ir)

      iap=iqe(ip)
      iag=iqc(ig)
      iaq=iqe(iq)
      iar=iqe(ir)

      iip=ke(ip)
      iiq=ke(iq)
      iir=ke(ir)
      iig=kc(ig)

      if(iip.ne.iorb)go to 40

      call findk(orbp,orbq,orbr,orbg,iap,iaq,iar,iag,kgot,
     $kmax)
      if(kmax.eq.0)go to 40
      do 50 nloop=1,kmax
      kk=kgot(nloop)
      npq=jdra(iip,iiq,kk)
      nrg=idra(iir,iig,kk)
      nn=nn+1
      neqn=nslocf(npq,nrg,ipqsym,kk)

      if(ncheck(neqn).eq.1)then
        print*,'double count in slocf =', neqn
      stop
      endif

      ncheck(neqn)=1

      indx(nn)=nn
      indx1(nn)=neqn

      sf(neqn)=sj(neqn)*(-1)**(orbp+orbq+orbr+orbg)

  50  continue
  40  continue
  30  continue
  20  continue
  15  continue

      nsmax=nn

9999  format(//,'DOUBLE EXCITATIONS',/)
9998  format(3x,'IN','   IIP ',' IIQ ',' IIS ',' IIA',
     $1x,'JP',3X,'JQ',3x,'JR',3x,'JA',3x,'PSM',3x,'QSM',3x,
     $'RSM',3x,'ASM',3x,'KMAX',/)
9997  format(5i5,4f5.1,6i6)
9996  format(/,2x,'COUPLED CLUSTER A MATRIX (IP)',/)
9995  format(4f20.11)
9991  format(//,'SINGLE EXCITATIONS',/)
9990  format(3x,'IN','   IIP ',' IIQ',' J(IIP)',' J(IIQ)'/)
9989  format(3i5,2(f5.1,2x))

      return
      end

c***************************************************************************8

      function nslocf(ipq,irg,ipar,kk)
      implicit real*8 (a-h,o-z)
c     common block
      common/skip1/npbps(0:MXV,2)
      common/skip6/mmskip(0:MXV,2),m2skip
      common/skip7/nsing

      if(ipar.eq.1)then
      nslocf=(ipq-1)*npbps(kk,2)+irg+mmskip(kk,2)+nsing
      else
      nslocf=(ipq-1)*npbps(kk,1)+irg+mmskip(kk,1)+m2skip
      endif
      return
      end

c***************************************************************
c     This subroutine  normalize the initial state             *
c***************************************************************

      subroutine cnormi(d,t,ival)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/norm/cnorm1,cnorm2
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
c     local dimension
 
      dimension t(MDIM),d(MNBAS,MNBAS)
      dimension igot(MXV+1),jgot(MXV+1),kgot(MXV+1)
      
      do i=0,MXV
      igot(i)=0
      kgot(i)=0
      jgot(i)=0
      enddo

c-----------------------------------------------------------------
c       diagram 1 = <p|d|p>
c-----------------------------------------------------------------

      cnorm=d(ival,ival)

      do 10 ip=1,nexcit
      iip=ke(ip)
      iap=iqe(ip)
      ipsym=isyme(ip)
      orbp=orbe(ip)

      if(iip.ne.ival)go to 10
     
      sum=0

      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0

c----------------------------------------------------------------
c     diagram 2,3 = <p|d|r><r|s1|p> + <p|s1|r><r|d|p>
c----------------------------------------------------------------

      do 20 ir=1,nexcit
      iir=ke(ir)
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)

      if(iir.eq.ival)go to 21
      if(irsym.ne.ipsym)go to 21
      if(orbp.ne.orbr)go to 21

      in=jdpa(iip,iir)

      sum=sum+2*d(iip,iir)*t(in)
      sum1=sum1+2*d(iip,iir)*t(in)
  
  21  continue

      do 30 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      iag=iqc(ig)

      if(igsym.ne.irsym)go to 30
      if(orbg.ne.orbr)go to 30

c------------------------------------------------------------
c     diagram 4,5 = <pg|s2|pr><r|d|g> + <r|d|g><pg|s2|pr>
c------------------------------------------------------------

      call findk(orbp,orbp,orbg,orbr,iap,iap,iag,iar,igot,imax)
      if(imax.eq.0)go to 30
      do 40 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.0)go to 40
      npp=jdra(iip,iip,l1)
      irgsym=mtbl(irsym,igsym)
      nrg=idra(iir,iig,l1)
      in=nsloci(npp,nrg,irgsym,l1)
      fact=dsqrt((2*orbg+1)/(2*orbp+1))
      sum=sum+2*d(iir,iig)*t(in)*fact
      sum2=sum2+2*d(iir,iig)*t(in)*fact
  
 40   continue

c--------------------------------------------------------------
c      diagram 6,7 = <pg|s2|rp><r|d|g> + <r|d|g><pg|s2|rp>
c--------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbp,iap,iar,iag,iap,igot,imax)
      if(imax.eq.0)go to 30
      do 50 iloop=1,imax
      l1=igot(iloop)
      npr=jdra(iip,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      npg=idra(iip,iig,l1)
      in=nsloci(npr,npg,ipgsym,l1)
      fact=(-1)**(orbp+orbr+2*orbr+l1)/(2*orbp+1)
      sum=sum-2*d(iir,iig)*t(in)*fact
      sum3=sum3-2*d(iir,iig)*t(in)*fact
 
 50   continue
 30   continue

      do 60 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      if(orbs.ne.orbr)go to 61
      if(issym.ne.irsym)go to 61
      if(iis.eq.iip)go to 61
      if(iir.eq.iip)go to 61
      if(orbp.ne.orbs)go to 61
      if(ipsym.ne.issym)go to 61
      if(orbp.ne.orbr)go to 61
      if(ipsym.ne.irsym)go to 61

      in1=jdpa(iip,iir)
      in2=jdpa(iip,iis)

c---------------------------------------------------------------
c         diagram 8 = <p|s1|r><r|d|s><s|s1|p>
c---------------------------------------------------------------

      sum4=sum4+d(iir,iis)*t(in1)*t(in2)
      sum=sum+d(iir,iis)*t(in1)*t(in2)  
 
  61  continue

      do 70 ig=1,nocc
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      igsym=isymc(ig)


      if(igsym.ne.irsym)go to 75
      if(orbg.ne.orbr)go to 75
      if(iis.eq.ival)go to 75
      if(orbs.ne.orbp)go to 75
      if(issym.ne.ipsym)go to 75

      irgsym=mtbl(irsym,igsym)
      ipssym=mtbl(ipsym,issym)

      if(irgsym.ne.ipssym)go to 75


c------------------------------------------------------------------
c          diagram 9,10 =  <pg|s2|sr><r|d|g><s|s1|p>
c                        + <p|s1|s><g|d|r><pg|s2|sr>
c-------------------------------------------------------------------


      in1=jdpa(iip,iis)
      nrg=idra(iir,iig,0)
      nps=jdra(iip,iis,0)
      irgsym=mtbl(irsym,igsym)
      in2=nsloci(nps,nrg,irgsym,0)
      fact=dsqrt((2*orbg+1)/(2*orbp+1))
      sum=sum+2*d(iig,iir)*t(in1)*t(in2)*fact
      sum5=sum5+2*d(iig,iir)*t(in1)*t(in2)*fact

c-----------------------------------------------------------------
c        diagram 11,12  =   <pg|s2|rs><r|d|g><s|s1|p>
c                         + <s|s1|p><g|d|r><pg|s2|rs>
c-----------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,igot,imax)
      if(imax.eq.0)go to 75
      do 80 iloop=1,imax
      l1=igot(iloop)
      npr=jdra(iip,iir,l1)
      nsg=idra(iis,iig,l1)
      isgsym=mtbl(issym,igsym)
      in2=nsloci(npr,nsg,isgsym,l1)
      fact=(-1)**(orbp+orbr+2*orbr+l1)/(2*orbp+1)
      sum=sum-2*d(iig,iir)*t(in1)*t(in2)*fact
      sum6=sum6-2*d(iig,iir)*t(in1)*t(in2)*fact 

 80   continue
 75   continue
 70   continue
 60   continue
 20   continue

      do 90 ir=1,nexcit
      iir=ke(ir)
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)

      do 90 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      iag=iqc(ig)

      do 90 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      do 100 it=1,nexcit
      itsym=isyme(it)
      orbt=orbe(it)
      iat=iqe(it)
      iit=ke(it)

      iptsym=mtbl(ipsym,itsym)

      if(issym.ne.itsym)go to 100
      if(orbs.ne.orbt)go to 100
      if(irgsym.ne.ipssym)go to 100
      if(irgsym.ne.iptsym)go to 100

c----------------------------------------------------------------
c     diagram 13 = <pg|s2|sr><s|d|t><pg|s2|tr>
c----------------------------------------------------------------

      call findk(orbp,orbs,orbr,orbg,iap,ias,iar,iag,igot,imax)
      if(imax.eq.0)go to 110
      do 120 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=jdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      irgsym=mtbl(irsym,igsym)
      npt=jdra(iip,iit,l1)
      in1=nsloci(nps,nrg,irgsym,l1)
      in2=nsloci(npt,nrg,irgsym,l1)
      fact=1/((2*l1+1)*(2*orbp+1))
      sum=sum+d(iis,iit)*t(in1)*t(in2)*fact
      sum7=sum7+d(iis,iit)*t(in1)*t(in2)*fact
 
 120  continue
 110  continue

c----------------------------------------------------------------
c     diagram 14 = <pg|s2|rs><s|d|t><pg|s2|tr>
c----------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,igot,imax)
      if(imax.eq.0)go to 130
      do 140 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=jdra(iip,iir,l1)
      nsg=idra(iis,iig,l1)
      isgsym=mtbl(issym,igsym)
      in1=nsloci(npr,nsg,isgsym,l1)

      call findk(orbp,orbt,orbr,orbg,iap,iat,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 140
      do 150 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npt=jdra(iip,iit,l2)
      nrg=idra(iir,iig,l2)
      irgsym=mtbl(irsym,igsym)
      in2=nsloci(npt,nrg,irgsym,l2)
      call sixj(al1,orbt,orbg,al2,orbr,orbp,s6j)
      fact=(-1)**(2*orbr)*s6j/(2*orbp+1)
      sum=sum-d(iis,iit)*t(in1)*t(in2)*fact
      sum8=sum8-d(iis,iit)*t(in1)*t(in2)*fact
 
  150 continue
  140 continue
  130 continue

c-----------------------------------------------------------------
c     diagram 15 = <pg|s2|rs><s|d|t><pg|s2|rt>
c-----------------------------------------------------------------
 
      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,igot,imax)
      if(imax.eq.0)go to 160
      do 170 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=jdra(iip,iir,l1)
      nsg=idra(iis,iig,l1)
      iprsym=mtbl(ipsym,irsym)
      ntg=idra(iit,iig,l1)    
      in1=nsloci(npr,nsg,iprsym,l1)
      in2=nsloci(npr,ntg,iprsym,l1)
      fact=1/((2*l1+1)*(2*orbp+1))
      sum=sum+d(iis,iit)*t(in1)*t(in2)*fact
      sum9=sum9+d(iis,iit)*t(in1)*t(in2)*fact
 
 170  continue
 160  continue

c-----------------------------------------------------------------
c     diagram 16 = <pg|s2|sr><s|d|t><pg|s2|rt>
c-----------------------------------------------------------------

      call findk(orbp,orbs,orbr,orbg,iap,ias,iar,iag,igot,imax)
      if(imax.eq.0)go to 180
      do 190 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=jdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      irgsym=mtbl(irsym,igsym)
      in1=nsloci(nps,nrg,irgsym,l1)

      call findk(orbp,orbr,orbt,orbg,iap,iar,iat,iag,jgot,jmax)
      if(jmax.eq.0)go to 200
      do 210 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npr=jdra(iip,iir,l2)
      ntg=idra(iit,iig,l2)
      itgsym=mtbl(itsym,igsym)
      in2=nsloci(npr,ntg,itgsym,l2)
      call sixj(al1,orbr,orbg,al2,orbt,orbp,s6j)
      fact=(-1)**(2*orbr)*s6j/(2*orbp+1)
      sum=sum-d(iis,iit)*t(in1)*t(in2)*fact
      sum10=sum10-d(iis,iit)*t(in1)*t(in2)*fact
 
  210 continue
  200 continue
  190 continue
  180 continue
  100 continue

      do 220 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(orbg.ne.orbd)go to 220
      if(igsym.ne.idsym)go to 220
      if(ipssym.ne.irgsym)go to 220
      if(ipssym.ne.irdsym)go to 220

c---------------------------------------------------------------
c     diagram 17 = <pg|s2|sr><g|d|d><pd|s2|sr>
c---------------------------------------------------------------

      call findk(orbp,orbs,orbr,orbg,iap,ias,iar,iag,igot,imax)
      if(imax.eq.0)go to 220
      do 230 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=jdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      in1=nsloci(nps,nrg,irgsym,l1)
      nrd=idra(iir,iid,l1)
      in2=nsloci(nps,nrd,irdsym,l1)
      fact=1/((2*l1+1)*(2*orbp+1))
      sum=sum-d(iid,iig)*t(in1)*t(in2)*fact
      sum11=sum11-d(iid,iig)*t(in1)*t(in2)*fact

c---------------------------------------------------------------
c     diagram 18 = <pg|s2|sr><g|d|d><pd|s2|rs>
c---------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbd,iap,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 230
      do 240 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npr=jdra(iip,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      in2=nsloci(npr,nsd,isdsym,l2)
      call sixj(al1,orbr,orbg,al2,orbs,orbp,s6j)
      fact=s6j*(-1)**(2*orbr)/(2*orbp+1)
      sum=sum+d(iid,iig)*t(in1)*t(in2)*fact
      sum12=sum12+d(iid,iig)*t(in1)*t(in2)*fact
 
 240  continue
 230  continue
 220  continue
  90  continue

      do 250 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      if(orbp.ne.orbs)go to 250
      if(iip.eq.iis)go to 250
      if(ipsym.ne.issym)go to 250

      in1=jdpa(iip,iis)

c--------------------------------------------------------------
c       diagram 19 = <p|s1|s><s|s1|p>
c--------------------------------------------------------------

      sum=sum+t(in1)*t(in1)
      sum13=sum13+t(in1)*t(in1)

 250  continue

      do 260 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      do 260 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)

      do 260 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

      irssym=mtbl(irsym,issym)
      ipgsym=mtbl(ipsym,igsym)
      irgsym=mtbl(irsym,igsym)

      if(ipgsym.ne.irssym)go to 260

c---------------------------------------------------------------
c        diagram 20 = <pg|s2|sr><pg|s2|sr>
c---------------------------------------------------------------

      call findk(orbp,orbs,orbr,orbg,iap,ias,iar,iag,igot,imax)
      if(imax.eq.0)go to 260
      do 270 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=jdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      in1=nsloci(nps,nrg,irgsym,l1)
      fact=1/((2*l1+1)*(2*orbp+1))
      sum14=sum14+fact*t(in1)*t(in1)
      sum=sum+fact*t(in1)*t(in1)

c---------------------------------------------------------------
c        diagram 21 = <pg|s2|sr><pg|s2|rs>
c---------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 270
      do 280 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npr=jdra(iip,iir,l2)
      nsg=idra(iis,iig,l2)
      isgsym=mtbl(issym,igsym)
      in2=nsloci(npr,nsg,isgsym,l2)
      call sixj(al1,orbr,orbg,al2,orbs,orbp,s6j)
      fact=s6j/(2*orbp+1)
      sum=sum+fact*t(in1)*t(in2)
      sum15=sum15+fact*t(in1)*t(in2)

 280  continue
 270  continue
 260  continue
  10  continue

      write(*,9000)dsqrt(1.0d0+sum+cnorm)
      cnorm1=dsqrt(1.0d0+sum+cnorm)

 9000 format(1x,'Normalization of the Initial state w.f=',f20.11)

      return
      end

c***************************************************************
c     This subroutine  normalize the final state               *
c***************************************************************


      subroutine cnormf(d,t,ival)

      implicit real*8 (a-h,o-z)
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/norm/cnorm1,cnorm2
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
c     local dimension
 
      dimension t(MDIM),d(MNBAS,MNBAS)
      dimension igot(MXV+1),jgot(MXV+1),kgot(MXV+1)
      
      do i=0,MXV
      igot(i)=0
      kgot(i)=0
      jgot(i)=0
      enddo

      cnorm=d(ival,ival)

c-----------------------------------------------------------------
c       diagram 1 = <p|d|p>
c-----------------------------------------------------------------

      do 10 ip=1,nexcit
      iip=ke(ip)
      ipsym=isyme(ip)
      iap=iqe(ip)
      ipsym=isyme(ip)
      orbp=orbe(ip)

      if(iip.ne.ival)go to 10
     
      sum=0

      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0

c----------------------------------------------------------------
c     diagram 2,3 = <p|d|r><r|s1|p> + <p|s1|r><r|d|p>
c----------------------------------------------------------------

      do 20 ir=1,nexcit
      iir=ke(ir)
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)

      if(iir.eq.ival)go to 21
      if(irsym.ne.ipsym)go to 21
      if(orbp.ne.orbr)go to 21

      in=kdpa(iip,iir)

      sum=sum+2*d(iip,iir)*t(in)
      sum1=sum1+2*d(iip,iir)*t(in) 
 
  21  continue

      do 30 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      iag=iqc(ig)

      if(igsym.ne.irsym)go to 30
      if(orbg.ne.orbr)go to 30

c------------------------------------------------------------
c     diagram 4,5 = <pg|s2|pr><r|d|g> + <r|d|g><pg|s2|pr>
c------------------------------------------------------------

      call findk(orbp,orbp,orbg,orbr,iap,iap,iag,iar,igot,imax)
      if(imax.eq.0)go to 30
      do 40 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.0)go to 40
      npp=kdra(iip,iip,l1)
      irgsym=mtbl(irsym,igsym)
      nrg=idra(iir,iig,l1)
      in=nslocf(npp,nrg,irgsym,l1)
      fact=dsqrt((2*orbg+1)/(2*orbp+1))
      sum=sum+2*d(iir,iig)*t(in)*fact
      sum2=sum2+2*d(iir,iig)*t(in)*fact
  
 40   continue

c--------------------------------------------------------------
c      diagram 6,7 = <pg|s2|rp><r|d|g> + <r|d|g><pg|s2|rp>
c--------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbp,iap,iar,iag,iap,igot,imax)
      if(imax.eq.0)go to 30
      do 50 iloop=1,imax
      l1=igot(iloop)
      npr=kdra(iip,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      npg=idra(iip,iig,l1)
      in=nslocf(npr,npg,ipgsym,l1)
      fact=(-1)**(orbp+orbr+2*orbr+l1)/(2*orbp+1)
      sum=sum-2*d(iir,iig)*t(in)*fact
      sum3=sum3-2*d(iir,iig)*t(in)*fact
 
 50   continue
 30   continue

      do 60 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      if(orbs.ne.orbr)go to 61
      if(issym.ne.irsym)go to 61
      if(iis.eq.iip)go to 61
      if(iir.eq.iip)go to 61
      if(orbp.ne.orbs)go to 61
      if(ipsym.ne.issym)go to 61
      if(orbp.ne.orbr)go to 61
      if(ipsym.ne.irsym)go to 61

      in1=kdpa(iip,iir)
      in2=kdpa(iip,iis)

c---------------------------------------------------------------
c         diagram 8 = <p|s1|r><r|d|s><s|s1|p>
c---------------------------------------------------------------

      sum4=sum4+d(iir,iis)*t(in1)*t(in2)
      sum=sum+d(iir,iis)*t(in1)*t(in2)
   
  61  continue

      do 70 ig=1,nocc
      iig=kc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      igsym=isymc(ig)

      if(igsym.ne.irsym)go to 75
      if(orbg.ne.orbr)go to 75
      if(iis.eq.ival)go to 75
      if(orbs.ne.orbp)go to 75
      if(issym.ne.ipsym)go to 75

      irgsym=mtbl(irsym,igsym)
      ipssym=mtbl(ipsym,issym)

      if(irgsym.ne.ipssym)go to 75
      
c------------------------------------------------------------------
c          diagram 9,10 =  <pg|s2|sr><r|d|g><s|s1|p>
c                        + <p|s1|s><g|d|r><pg|s2|sr>
c-------------------------------------------------------------------


      in1=kdpa(iip,iis)
      nrg=idra(iir,iig,0)
      nps=kdra(iip,iis,0)
      irgsym=mtbl(irsym,igsym)
      in2=nslocf(nps,nrg,irgsym,0)
      fact=dsqrt((2*orbg+1)/(2*orbp+1))
      sum=sum+2*d(iig,iir)*t(in1)*t(in2)*fact
      sum5=sum5+2*d(iig,iir)*t(in1)*t(in2)*fact
 
c-----------------------------------------------------------------
c        diagram 11,12  =   <pg|s2|rs><r|d|g><s|s1|p>
c                         + <s|s1|p><g|d|r><pg|s2|rs>
c-----------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,igot,imax)
      if(imax.eq.0)go to 75
      do 80 iloop=1,imax
      l1=igot(iloop)
      npr=kdra(iip,iir,l1)
      nsg=idra(iis,iig,l1)
      isgsym=mtbl(issym,igsym)
      in2=nslocf(npr,nsg,isgsym,l1)
      fact=(-1)**(orbp+orbr+2*orbr+l1)/(2*orbp+1)
      sum=sum-2*d(iig,iir)*t(in1)*t(in2)*fact
      sum6=sum6-2*d(iig,iir)*t(in1)*t(in2)*fact 

 80   continue
 75   continue
 70   continue
 60   continue
 20   continue

      do 90 ir=1,nexcit
      iir=ke(ir)
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)

      do 90 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      iag=iqc(ig)

      do 90 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      do 100 it=1,nexcit
      itsym=isyme(it)
      orbt=orbe(it)
      iat=iqe(it)
      iit=ke(it)

      iptsym=mtbl(ipsym,itsym)

      if(issym.ne.itsym)go to 100
      if(orbs.ne.orbt)go to 100
      if(irgsym.ne.ipssym)go to 100
      if(irgsym.ne.iptsym)go to 100

c----------------------------------------------------------------
c     diagram 13 = <pg|s2|sr><s|d|t><pg|s2|tr>
c----------------------------------------------------------------

      call findk(orbp,orbs,orbr,orbg,iap,ias,iar,iag,igot,imax)
      if(imax.eq.0)go to 110
      do 120 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      irgsym=mtbl(irsym,igsym)
      npt=kdra(iip,iit,l1)
      in1=nslocf(nps,nrg,irgsym,l1)
      in2=nslocf(npt,nrg,irgsym,l1)
      fact=1/((2*l1+1)*(2*orbp+1))
      sum=sum+d(iis,iit)*t(in1)*t(in2)*fact
      sum7=sum7+d(iis,iit)*t(in1)*t(in2)*fact
 
 120  continue
 110  continue

c----------------------------------------------------------------
c     diagram 14 = <pg|s2|rs><s|d|t><pg|s2|tr>
c----------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,igot,imax)
      if(imax.eq.0)go to 130
      do 140 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      nsg=idra(iis,iig,l1)
      isgsym=mtbl(issym,igsym)
      in1=nslocf(npr,nsg,isgsym,l1)

      call findk(orbp,orbt,orbr,orbg,iap,iat,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 140
      do 150 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npt=kdra(iip,iit,l2)
      nrg=idra(iir,iig,l2)
      irgsym=mtbl(irsym,igsym)
      in2=nslocf(npt,nrg,irgsym,l2)
      call sixj(al1,orbt,orbg,al2,orbr,orbp,s6j)
      fact=(-1)**(2*orbr)*s6j/(2*orbp+1)
      sum=sum-d(iis,iit)*t(in1)*t(in2)*fact
      sum8=sum8-d(iis,iit)*t(in1)*t(in2)*fact
 
  150 continue
  140 continue
  130 continue

c-----------------------------------------------------------------
c     diagram 15 = <pg|s2|rs><s|d|t><pg|s2|rt>
c-----------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,igot,imax)
      if(imax.eq.0)go to 160
      do 170 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      nsg=idra(iis,iig,l1)
      iprsym=mtbl(ipsym,irsym)
      ntg=idra(iit,iig,l1)    
      in1=nslocf(npr,nsg,iprsym,l1)
      in2=nslocf(npr,ntg,iprsym,l1)
      fact=1/((2*l1+1)*(2*orbp+1))
      sum=sum+d(iis,iit)*t(in1)*t(in2)*fact
      sum9=sum9+d(iis,iit)*t(in1)*t(in2)*fact
 
 170  continue
 160  continue

c-----------------------------------------------------------------
c     diagram 16 = <pg|s2|sr><s|d|t><pg|s2|rt>
c-----------------------------------------------------------------

      call findk(orbp,orbs,orbr,orbg,iap,ias,iar,iag,igot,imax)
      if(imax.eq.0)go to 180
      do 190 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      irgsym=mtbl(irsym,igsym)
      in1=nslocf(nps,nrg,irgsym,l1)

      call findk(orbp,orbr,orbt,orbg,iap,iar,iat,iag,jgot,jmax)
      if(jmax.eq.0)go to 200
      do 210 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npr=kdra(iip,iir,l2)
      ntg=idra(iit,iig,l2)
      itgsym=mtbl(itsym,igsym)
      in2=nslocf(npr,ntg,itgsym,l2)
      call sixj(al1,orbr,orbg,al2,orbt,orbp,s6j)
      fact=(-1)**(2*orbr)*s6j/(2*orbp+1)
      sum=sum-d(iis,iit)*t(in1)*t(in2)*fact
      sum10=sum10-d(iis,iit)*t(in1)*t(in2)*fact
 
  210 continue
  200 continue
  190 continue
  180 continue
  100 continue

      do 220 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)

      if(orbg.ne.orbd)go to 220
      if(igsym.ne.idsym)go to 220
      if(ipssym.ne.irgsym)go to 220
      if(ipssym.ne.irdsym)go to 220

c---------------------------------------------------------------
c     diagram 17 = <pg|s2|sr><g|d|d><pd|s2|sr>
c---------------------------------------------------------------
      
      call findk(orbp,orbs,orbr,orbg,iap,ias,iar,iag,igot,imax)
      if(imax.eq.0)go to 220
      do 230 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      in1=nslocf(nps,nrg,irgsym,l1)
      nrd=idra(iir,iid,l1)
      in2=nslocf(nps,nrd,irdsym,l1)
      fact=1/((2*l1+1)*(2*orbp+1))
      sum=sum-d(iid,iig)*t(in1)*t(in2)*fact
      sum11=sum11-d(iid,iig)*t(in1)*t(in2)*fact

c---------------------------------------------------------------
c     diagram 18 = <pg|s2|sr><g|d|d><pd|s2|rs>
c---------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbd,iap,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 230
      do 240 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npr=kdra(iip,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      in2=nslocf(npr,nsd,isdsym,l2)
      call sixj(al1,orbr,orbg,al2,orbs,orbp,s6j)
      fact=s6j*(-1)**(2*orbr)/(2*orbp+1)
      sum=sum+d(iid,iig)*t(in1)*t(in2)*fact
      sum12=sum12+d(iid,iig)*t(in1)*t(in2)*fact 

 240  continue
 230  continue
 220  continue
  90  continue

      do 250 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      if(orbp.ne.orbs)go to 250
      if(iip.eq.iis)go to 250
      if(ipsym.ne.issym)go to 250

      in1=kdpa(iip,iis)

c--------------------------------------------------------------
c       diagram 19 = <p|s1|s><s|s1|p>
c--------------------------------------------------------------

      sum=sum+t(in1)*t(in1)
      sum13=sum13+t(in1)*t(in1)

 250  continue


      do 260 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      ias=iqe(is)
      iis=ke(is)

      do 260 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iar=iqe(ir)
      iir=ke(ir)

      do 260 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)

      irssym=mtbl(irsym,issym)
      ipgsym=mtbl(ipsym,igsym)
      irgsym=mtbl(irsym,igsym)

      if(ipgsym.ne.irssym)go to 260

c---------------------------------------------------------------
c        diagram 20 = <pg|s2|sr><pg|s2|sr>
c---------------------------------------------------------------

      call findk(orbp,orbs,orbr,orbg,iap,ias,iar,iag,igot,imax)
      if(imax.eq.0)go to 260
      do 270 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      in1=nslocf(nps,nrg,irgsym,l1)
      fact=1/((2*l1+1)*(2*orbp+1))
      sum14=sum14+fact*t(in1)*t(in1)
      sum=sum+fact*t(in1)*t(in1)

c---------------------------------------------------------------
c        diagram 21 = <pg|s2|sr><pg|s2|rs>
c---------------------------------------------------------------

      call findk(orbp,orbr,orbs,orbg,iap,iar,ias,iag,jgot,jmax)
      if(jmax.eq.0)go to 270
      do 280 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      npr=kdra(iip,iir,l2)
      nsg=idra(iis,iig,l2)
      isgsym=mtbl(issym,igsym)
      in2=nslocf(npr,nsg,isgsym,l2)
      call sixj(al1,orbr,orbg,al2,orbs,orbp,s6j)
      fact=s6j/(2*orbp+1)
      sum=sum+fact*t(in1)*t(in2)
      sum15=sum15+fact*t(in1)*t(in2)

 280  continue
 270  continue
 260  continue
  10  continue

      write(*,9000)dsqrt(1.0d0+sum+cnorm)
      cnorm2=dsqrt(1.0d0+sum+cnorm)

 9000 format(1x,'Normalization of the Final state w.f=',f20.11)

      return
      end

c*********************************************************************

      subroutine flopn(ival,jval)
      implicit real*8 (a-h,o-z)
      character*8 siampl,sfampl
      read(STDIN,'(A)')siampl
      read(STDIN,'(A)')sfampl
      write(STDOUT,'(A)')siampl
      write(STDOUT,'(A)')sfampl
      open(unit=NIFILE,file=siampl,form='unformatted',
     : access='sequential',status='old')
c     if(ival.ne.jval)then
      open(unit=NFFILE,file=sfampl,form='unformatted',
     : access='sequential',status='old')
c     endif
      return
      end

c******************************************************************

c******************************************************************
c     This subroutine calculate the finel result with and without *
c     normalization correction to the magnetic dipole hyperfine   *
c     matrix element                                              *
c******************************************************************
     
      subroutine calhyp

      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/value/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,sum45,sum46

      cnorm=cnorm1*cnorm2
       write(STDOUT,*)'sourav666'

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44
      write(STDOUT,27)sum45
      write(STDOUT,28)sum46

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26
     :    +sum27+sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40+sum41+sum42
     :    +sum43+sum44+sum45+sum46

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)
   27 format(/,1x,'contribution from diagram 45',f20.11)
   28 format(/,1x,'contribution from diagram 46',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

c********************************************************************

      subroutine calC
      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valueC/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26+sum27
     :    +sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40+sum41+sum42
     :    +sum43+sum44

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

c******************************************************************

c******************************************************************
c     This subroutine calculate the finel result with and without *
c     normalization correction to the electric dipole transition  *
c     matrix element                                              *
c******************************************************************

      subroutine caldip

      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valued/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,sum45,sum46

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44
      write(STDOUT,27)sum45
      write(STDOUT,28)sum46

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26
     :    +sum27+sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40+sum41+sum42
     :    +sum43+sum44+sum45+sum46

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)
   27 format(/,1x,'contribution from diagram 45',f20.11)
   28 format(/,1x,'contribution from diagram 46',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

c**********************************************************************
c                                                                     *
c      This subroutine calculates the E2 transition property diagrams *
c      with the total results                                         *
c                                                                     *
c**********************************************************************


      subroutine calquad

      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valuequad/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,d3j

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1*d3j
      write(STDOUT,9981)sum2*d3j
      write(STDOUT,9982)sum3*d3j
      write(STDOUT,9972)sum4*d3j
      write(STDOUT,9983)sum5*d3j
      write(STDOUT,9973)sum6*d3j
      write(STDOUT,9984)sum7*d3j
      write(STDOUT,9985)sum8*d3j
      write(STDOUT,9975)sum9*d3j
      write(STDOUT,9986)sum10*d3j
      write(STDOUT,9976)sum11*d3j
      write(STDOUT,9987)sum12*d3j
      write(STDOUT,9988)sum13*d3j
      write(STDOUT,9989)sum14*d3j
      write(STDOUT,9990)sum15*d3j
      write(STDOUT,9991)sum16*d3j
      write(STDOUT,9992)sum17*d3j
      write(STDOUT,9993)sum18*d3j
      write(STDOUT,1)sum19*d3j
      write(STDOUT,2)sum20*d3j
      write(STDOUT,3)sum21*d3j
      write(STDOUT,4)sum22*d3j
      write(STDOUT,5)sum23*d3j
      write(STDOUT,6)sum24*d3j
      write(STDOUT,7)sum25*d3j
      write(STDOUT,8)sum26*d3j
      write(STDOUT,9)sum27*d3j
      write(STDOUT,10)sum28*d3j
      write(STDOUT,11)sum29*d3j
      write(STDOUT,12)sum30*d3j
      write(STDOUT,13)sum31*d3j
      write(STDOUT,14)sum32*d3j
      write(STDOUT,15)sum33*d3j
      write(STDOUT,16)sum34*d3j
      write(STDOUT,17)sum35*d3j
      write(STDOUT,18)sum36*d3j
      write(STDOUT,19)sum37*d3j
      write(STDOUT,20)sum38*d3j
      write(STDOUT,21)sum39*d3j
      write(STDOUT,22)sum40*d3j
      write(STDOUT,23)sum41*d3j
      write(STDOUT,24)sum42*d3j
      write(STDOUT,25)sum43*d3j
      write(STDOUT,26)sum44*d3j

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26
     :    +sum27+sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40
     :    +sum43+sum44

      corr=sum*fact*d3j
      write(*,9994)corr
      write(*,9995)sum*d3j
      write(*,9996)sum*d3j/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

c*********************************************************************

      subroutine calquadvel
      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valuequadvel/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26
c    :    +sum27+sum28+sum29+sum30+sum31+sum32
c    :    +sum33+sum34+sum35+sum36+sum37
c    :    +sum38+sum39+sum40
c    :    +sum43+sum44

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end


      subroutine grasprd(nbasis)
      implicit real*8 (a-h,o-z)

      common/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :       mf(MNBAS)
     :      /grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /tatb/ta(MN),tb(MN),mtp
      common/orbital_energy/eorb(MNBAS)
c     local dimension
      write(STDOUT,*)'sourav000000000'
    
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
      write(STDOUT,*)'99999999999999999'

      mtp=n

9999  format(//,10x,17h ORBITAL ENERGIES,/)
9998  format(4d20.11)
      return
      end

c*********************************************************
    
      


      subroutine CalcC(Nbasis)
c*********************************************************
c  this subrouti calculates the edm t-pt reduced
c  matrix elements between orbitals with j=1/2
c*********************************************************
      implicit double precision(a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)

       common/radialc/hyp(MNBAS,MNBAS)

       mtp=n

       do i1=1,nbasis
       do i2=1,i1
       ta(1)=zero
       do 1 l=2,mtp
       ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*
     : rp(l)/(r(l)**2*r(l)**2)
    1  continue
       call quad(result)
       hyp(i1,i2)=result
       hyp(i2,i1)=result
       enddo
       enddo

       return
       end


c****************************************************************

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


      subroutine findk(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     $jeven,kkk)
c     this subroutine finds the common k value
      implicit real*8 (a-h,o-z)
c     include 'mpif.h'
      common/iallparity/iall
      dimension jeven(MXV+1)

            if(iall.eq.0)then

      call findk2(orba,orbb,orbp,orbq,ia,ib,ip,iq,
     & jeven,kkk)
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

c**********************************************************************
c    This subroitine calculates the electric dipole property diagram  *
c**********************************************************************
 
      subroutine dipvelprop(d,si,sf,ti,tf,ival,jval)

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
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)
c
      common/valuedvel/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS),si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0

      irank=1
      rank=1.0d0

c--------p is the final and q is the initial state-----------

      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iiq.ne.ival)go to 1
      if(ipsym.eq.iqsym) go to 1

c-------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c-------------------------------------------------------

      call e1vel(iip,iiq,orbp,orbq,kapp,kapq,hyp1)
      print*,'Dirac-Fock value for E1_vel ',hyp1

      sum=sum+d(iip,iiq)
      sum1=sum1+d(iip,iiq)


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

c     if(iir.ne.ival)then
      if(irsym.ne.ipsym)then
      if(orbq.eq.orbr)then
      if(irsym.eq.iqsym)then
      if(iir.ne.iiq)then

      in=jdpa(iiq,iir)

c------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c------------------------------------------------------

      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)

      endif
      endif
      endif
      endif

      if(irsym.ne.iqsym)then
      if(orbr.eq.orbp)then
      if(irsym.eq.ipsym)then
      if(iip.ne.iir)then

      in=kdpa(iip,iir)

c-------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-------------------------------------------------------

      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)

      endif
      endif
      endif
      endif
c     endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.eq.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact

  6   continue

      endif

c----------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c----------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact

  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      in=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(in)*fact
      sum6=sum6+d(iir,iig)*sf(in)*fact

  9   continue

      endif

c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      in=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(in)*fact
      sum7=sum7-d(iir,iig)*sf(in)*fact

 10   continue
  8   continue
  4   continue
  3   continue

      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c--------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c--------------------------------------------------------------

      if(issym.eq.irsym)go to 11
      if(orbp.ne.orbr)go to 11
      if(ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs)go to 11
      if(iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)

 11   continue

      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.eq.igsym)go to 12

      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(iqssym.ne.irgsym)go to 14

c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact

 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact

 18   continue
 17   continue
 14   continue

      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.and.iqsym.ne.issym)go to 20

      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c-----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c-----------------------------------------------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact

 22   continue
 21   continue

c-----------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c-----------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue
  2   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.eq.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

  42  continue

c---------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c---------------------------------------------------------------

      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact

 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact

 48   continue
 47   continue
 46   continue

c-----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c-----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2) 
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact

 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)

      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52
      if(igsym.eq.idsym)go to 52

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)

c------------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c------------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      irdsym=mtbl(irsym,idsym)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=(-1)**(orbp+orbq+irank+2*orbg)*s6j1*s6j2
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact

  54  continue

c-----------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact

  55  continue
  53  continue
  52  continue
  30  continue

c------Two-body effective operator (computed directly)--------------

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

      hyp1=0

      if(ipsym.ne.issym)then
      call e1vel(iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif


      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then

      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)

      if(iqssym.eq.irasym)then

c---------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c---------------------------------------------------------------
      
      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact

      endif
      enddo
      endif

c---------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c---------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif

      hyp2=0

      if(issym.ne.iqsym)then
      call e1vel(iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      if(ipssym.eq.irasym)then

c-----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

      hyp1=0

      if(itsym.ne.issym)then
      call e1vel(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0

      if(itsym.ne.issym)then
      call e1vel(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)

c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      enddo
      endif
      endif


      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

      hyp1=0

      if(ipsym.ne.ibsym)then
      call e1vel(iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.iqssym)then
      if(irasym.eq.iqssym)then
      if(orbq.eq.orbb)then

c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      hyp2=0

      if(iqsym.ne.ibsym)then
      call e1vel(iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.ipssym)then
      if(irasym.eq.isbsym)then
      if(orbp.eq.orbb)then

c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.eq.igsym)go to 610
      hyp=0
c     c1-c2
      call e1vel(iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620
      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)


      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact
 640  continue
      endif

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact
  650 continue
      endif

  630 continue
  620 continue

c     d1-d2
      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660
      in=idpa(iiq,iib)

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact
  680 continue
      endif

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact
  690 continue
      endif

  670 continue
  660 continue
  610 continue


      if(ibsym.eq.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
c     if(orbp.ne.orbg)go to 605
      if(orbp.ne.orbb)go to 605
      in=idpa(iip,iib)
c     c7-c8
      hyp=0
      call e1vel(iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700
      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact

  720 continue
  710 continue
  605 continue
  700 continue
c     d7-d8

      hyp=0
      if(ibsym.ne.igsym)then
      call e1vel(iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600
      in=idpa(iiq,iib)
      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact
  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0
      if(itsym.ne.issym)then
      call e1vel(iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)
c     c3
      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

  830 continue
  820 continue
  815 continue

c     d3
      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900
      in=idpa(iiq,iib)

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact
  930 continue
  900 continue

  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)
      hyp=0
      if(itsym.ne.ibsym)then
      call e1vel(iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nra=idra(iir,iia,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact
      endif
      enddo
      endif

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact
      enddo
      endif
      enddo
      endif

      endif

      endif
      endif
      endif
      endif

      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then
      in=jdpa(iiq,iit)

      hyp=0
      if(itsym.ne.ibsym)then
      call e1vel(iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact
      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue

   1  continue
 100  continue
      return
      end

c********************************************************************

c********************************************************************
c                                                                   *
c    This subroitine calculates the E2 transition property diagram  *
c                                                                   *
c********************************************************************

      subroutine quadpuleprop(d,si,sf,ti,tf,ival,jval)

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
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)


      common/valuequad/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,d3j

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0

      irank=2
      rank=2.0d0

c-----p is the final and q is the initial state----------------

      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(iiq.ne.ival)go to 1

      if(ipsym.ne.iqsym) go to 1

      if(ival.eq.jval)then
      d3j=dr(orbp,two,orbq,-orbp,zero,orbq)
      else
      d3j=1.0d0
      endif

c---------------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c---------------------------------------------------------------

      sum=sum+d(iip,iiq)

      call quadpul(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dirac)
      print*,'Dirac-Fock value E2 ',dirac*d3j
      sum1=sum1+d(iip,iiq)


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(iir.ne.iiq)then
      if(irsym.eq.ipsym)then
      if(irsym.eq.iqsym.and.orbq.eq.orbr)then

      in=jdpa(iiq,iir)

c------------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c------------------------------------------------------------

      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)

      endif
      endif
      endif

      if(irsym.eq.iqsym)then
      if(irsym.eq.ipsym.and.orbr.eq.orbp)then
      if(iir.ne.iip)then

      in=kdpa(iip,iir)

c-----------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-----------------------------------------------------------

      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)

      endif
      endif
      endif


  2    continue

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact

  6   continue

      endif

c------------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c------------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact

  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      jn=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(jn)*fact
      sum6=sum6+d(iir,iig)*sf(jn)*fact

  9   continue

      endif

c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      jn=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(jn)*fact
      sum7=sum7-d(iir,iig)*sf(jn)*fact

 10   continue
  8   continue
  4   continue
  3   continue

      do 11 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c-------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c-------------------------------------------------------------

      if(irsym.ne.issym)go to 11
      if(orbp.ne.orbr.and.ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs.and.iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)

 11   continue

      do 12 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 12

      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(iqssym.ne.irgsym)go to 14

c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact

 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact

 18   continue
 17   continue
 14   continue

      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.and.iqsym.ne.issym)go to 20

      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c----------------------------------------------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact

 22   continue
 21   continue

c------------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.ne.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

  42  continue

c-----------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c-----------------------------------------------------------------

      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact

 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact

 48   continue
 47   continue
 46   continue

c----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact

 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
 
      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)

      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52

      if(igsym.ne.idsym)go to 52

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)

c-----------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c-----------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbp+orbq+irank+2*orbg)
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact

  54  continue

c------------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c------------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact

  55  continue
  53  continue
  52  continue
  30  continue
   
      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

      hyp1=0

      if(ipsym.eq.issym)then
      call quadpul(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif

      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then
        
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      if(iqssym.eq.irasym)then

c-----------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif

      hyp2=0

      if(issym.eq.iqsym)then
      call quadpul(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      if(ipssym.eq.irasym)then

c----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then


      hyp1=0

      if(itsym.eq.issym)then
      call quadpul(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif


      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0

      if(itsym.eq.issym)then
      call quadpul(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)

c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------


      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif
      endiF
      endif

      enddo
      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)


      hyp1=0

      if(ipsym.eq.ibsym)then
      call quadpul(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.iqssym)then
      if(irasym.eq.isbsym)then
      if(orbq.eq.orbb)then

c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c-----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      hyp2=0

      if(iqsym.eq.ibsym)then
      call quadpul(orbb,orbq,iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)

      if(irasym.eq.ipssym)then
      if(irasym.eq.isbsym)then
      if(orbp.eq.orbb)then

c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.ne.igsym)go to 610
      hyp=0
c     c1-c2
      call quadpul(orbg,orba,iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620
      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact
 640  continue
      endif

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact
  650 continue
      endif

  630 continue
  620 continue

c     d1-d2
      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660
      in=idpa(iiq,iib)

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact
  680 continue
      endif

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact
  690 continue
      endif

  670 continue
  660 continue
  610 continue

      if(ibsym.ne.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
c     if(orbp.ne.orbg)go to 605
      if(orbp.ne.orbb)go to 605
      in=idpa(iip,iib)
c     c7-c8
      hyp=0
      call quadpul(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700
      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)

      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)

      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)


      if(irasym.eq.isgsym)then

      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact
      endif


      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)

      isasym=mtbl(issym,iasym)
      irgsym=mtbl(irsym,igsym)

      inf=ntloc(nsa,nrg,isasym,l1)

      if(irgsym.eq.isasym)then

      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact
      endif

  720 continue
  710 continue
  700 continue
  605 continue
c     d7-d8

      hyp=0
      if(ibsym.eq.igsym)then
      call quadpul(orbg,orbb,iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600
      in=idpa(iiq,iib)
      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact
  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0
      if(itsym.eq.issym)then
      call quadpul(orbs,orbt,iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)
c     c3
      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

  830 continue
  820 continue
  815 continue

c     d3
      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900
      in=idpa(iiq,iib)

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact
  930 continue
  900 continue

  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)
      hyp=0
      if(itsym.eq.ibsym)then
      call quadpul(orbt,orbb,iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nqs=jdra(iiq,iis,l1)
      nra=idra(iir,iia,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))

      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact
      endif
      enddo
      endif

      call findk(orba,orbs,orbq,orbr,iaa,ias,iaq,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)

      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif
      enddo
      endif

      endif

      endif
      endif
      endif
      endif

      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then
      in=jdpa(iiq,iit)

      hyp=0
      if(itsym.eq.ibsym)then
      call quadpul(orbb,orbt,iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue
   1  continue
 100  continue

c      print*,'Dirac-Fock value E2 ',dirac*d3j

      return
      end

     
      subroutine quadvelprop(d,si,sf,ti,tf,ival,jval)
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
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)


      common/valuequadvel/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0

      irank=2
      rank=2.0d0

c     p is the final and q is the initial state
      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)
      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
      if(iiq.ne.ival)go to 1
      if(ipsym.ne.iqsym) go to 1
c     diagram 1 =<f|d|i>
      sum=sum+d(iip,iiq)
      call e2vel(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dipole)
      print*,'Dirac-Fock Value E2_vel ',dipole
      sum1=sum1+d(iip,iiq)
c     diagram < | ds|>

      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(iir.ne.iiq)then
      if(irsym.eq.ipsym)then
      if(irsym.eq.iqsym.and.orbq.eq.orbr)then
      in=jdpa(iiq,iir)
c     diagram 2(a)=  sum (r) <p|d|r><r|s|q>
      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)
      endif
      endif
      endif

      if(irsym.eq.iqsym)then
      if(irsym.eq.ipsym.and.orbr.eq.orbp)then
      if(iir.ne.iip)then
      in=kdpa(iip,iir)
c     diagram 3(a) =  sum (r) <p|s|r><r|d|q>
      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)
      endif
      endif
      endif

 2    continue

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      if(irsym.ne.igsym)go to 4
      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)
      if(ipqsym.ne.irgsym)go to 4
c     diagram 2 (b) = sum(g,r) <g|d|r><pr|s|qg>
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
c     diagram 2 (c) =-sum(g,r) <g|d|r><pr|s|gq>
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact
  6   continue
      endif
c     diagram 2 (c)
      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact
  7   continue
  5   continue
c     diagram 3 (b)
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      jn=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(jn)*fact
      sum6=sum6+d(iir,iig)*sf(jn)*fact
  9   continue
      endif
c     diagram 3 (c)
      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      jn=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(jn)*fact
      sum7=sum7-d(iir,iig)*sf(jn)*fact
 10   continue
  8   continue
  4   continue
  3   continue
c     diagram 4e
      do 11 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      if(orbr.ne.orbs)go to 11
      if(orbp.ne.orbr.or.ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs.or.iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11
      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)
      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)
 11   continue

c     diagram 4a
      do 12 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)
      if(irsym.ne.igsym)go to 12
      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14
      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)
      if(iqssym.ne.irgsym)go to 14
      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact
 16   continue
 15   continue
c     diagram 4b
      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact
 18   continue
 17   continue
 14   continue
      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.or.iqsym.ne.issym)go to 20
      ini=jdpa(iiq,iis)
      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)
      if(ipssym.ne.irgsym)go to 20
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact
 22   continue
 21   continue
c     diagram 4d
      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact
 24   continue
 23   continue
 20   continue
 13   continue
 12   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.ne.itsym)go to 40
      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)
      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41
      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
c     diagram 4f
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact
      endif
      enddo
      endif

  42  continue
c     diagram 4(g)
      itgsym=mtbl(itsym,igsym)
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact
 45   continue
 44   continue
 43   continue
c     diagram 4(h)
      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact
 48   continue
 47   continue
 46   continue
c     diagram 4(i)
      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)
      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact
 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

c     diagram 4(j)
      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)
      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52

      if(igsym.ne.idsym)go to 52
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)
      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbp+orbq+irank+2*orbg)
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact
  54  continue
      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact
  55  continue
  53  continue
  52  continue
  30  continue

C     Two-body effective operator (computed directly)

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)
      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)
      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

c     diagrams  a1-a4 (sum19-sum22)

      hyp1=0
      if(ipsym.eq.issym)then
      call e2vel(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif

      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact
      enddo
      endif

      endif

      hyp2=0
      if(issym.eq.iqsym)then
      call e2vel(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)
      ipssym=mtbl(ipsym,issym)

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact
      enddo
      endif

      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

c     diagram 43
      hyp1=0
      if(itsym.eq.issym)then
      call e2vel(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 45
      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0
      if(itsym.eq.issym)then
      call e2vel(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)
      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact
      endif
      enddo
      endif

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

      enddo

      endif
      endif


      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

c     diagrams  b1-b4 (sum19-sum22)
      hyp1=0
      if(ipsym.eq.ibsym)then
      call e2vel(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.eq.iqssym)then
      if(orbq.eq.orbb)then

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)

      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)

      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      hyp2=0
      if(iqsym.eq.ibsym)then
      call e2vel(orbb,orbq,iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)
      if(irasym.eq.ipssym)then
      if(orbp.eq.orbb)then

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.ne.igsym)go to 610
      hyp=0
c     c1-c2
      call e2vel(orbg,orba,iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620
      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)


      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact
 640  continue
      endif

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact
  650 continue
      endif

  630 continue
  620 continue

c     d1-d2
      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660
      in=idpa(iiq,iib)

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact
  680 continue
      endif

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact
  690 continue
      endif

  670 continue
  660 continue
  610 continue


      if(ibsym.ne.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
c     if(orbp.ne.orbg)go to 605
      if(orbp.ne.orbb)go to 605
      in=idpa(iip,iib)
c     c7-c8
      hyp=0
      call e2vel(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700
      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)

      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)

      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)


      if(irasym.eq.isgsym)then

      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact
      endif


      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)

      isasym=mtbl(issym,iasym)
      irgsym=mtbl(irsym,igsym)

      inf=ntloc(nsa,nrg,isasym,l1)

      if(irgsym.eq.isasym)then

      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact
      endif

  720 continue
  710 continue
  700 continue
  605 continue
c     d7-d8

      hyp=0
      if(ibsym.eq.igsym)then
      call e2vel(orbg,orbb,iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600
      in=idpa(iiq,iib)
      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact
  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0
      if(itsym.eq.issym)then
      call e2vel(orbs,orbt,iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)
c     c3
      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

  830 continue
  820 continue
  815 continue

c     d3
      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900
      in=idpa(iiq,iib)

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact
  930 continue
  900 continue

  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)
      hyp=0
      if(itsym.eq.ibsym)then
      call e2vel(orbt,orbb,iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nqs=jdra(iiq,iis,l1)
      nra=idra(iir,iia,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))

      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact
      endif
      enddo
      endif

      call findk(orba,orbs,orbq,orbr,iaa,ias,iaq,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)

      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif
      enddo
      endif

      endif

      endif
      endif
      endif
      endif

      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then
      in=jdpa(iiq,iit)

      hyp=0
      if(itsym.eq.ibsym)then
      call e2vel(orbb,orbt,iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact
      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue

 999  continue
   1  continue
 100  continue

      return
      end

c****************************************************************

c*******************************************************************
c    This subroitine calculates the E2 hyperfine property diagram  *
c*******************************************************************


      subroutine quadBprop(d,si,sf,ti,tf,ival,jval,qm,jinit,jfin)

      implicit real*8 (a-h,o-z)
      real*8 jorbv,jinit,jfin

      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      common/valueB/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,sum45,sum46

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      data VL/137.03598950D+00/
      data EMF/1836.152701D+00/
      data AUEV/27.2116D+00/
      data EVMHZ/2.417988D+08/

      CONVF=AUEV*EVMHZ
      CON=2*QM*CONVF*3.5711D-08

      d1=dr(jinit,2.0d0,jfin,-jfin,0.0d0,jfin)
      CON=d1*CON

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0
      sum45=0
      sum46=0

      irank=2
      rank=2.0d0

c-----p is the final and q is the initial state----------------

      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      if(jfin.ne.orbp)go to 100
      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(jinit.ne.orbq)go to 1
      if(iiq.ne.ival)go to 1

      if(ipsym.ne.iqsym) go to 1

c---------------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c---------------------------------------------------------------

      sum=sum+d(iip,iiq)

      call quadB(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dirac)
      sum1=sum1+d(iip,iiq)


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(iir.ne.iiq)then
      if(irsym.eq.ipsym)then
      if(irsym.eq.iqsym.and.orbq.eq.orbr)then

      in=jdpa(iiq,iir)

c------------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c------------------------------------------------------------

      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)

      endif
      endif
      endif

      if(irsym.eq.iqsym)then
      if(irsym.eq.ipsym.and.orbr.eq.orbp)then
      if(iir.ne.iip)then

      in=kdpa(iip,iir)

c-----------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-----------------------------------------------------------

      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)

      endif
      endif
      endif

 2    continue

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+al1)/(2*al1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact

  6   continue

      endif

c------------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c------------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact

  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      jn=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+al1)/(2*al1+1)
      sum=sum+d(iir,iig)*sf(jn)*fact
      sum6=sum6+d(iir,iig)*sf(jn)*fact

  9   continue

      endif

c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      jn=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(jn)*fact
      sum7=sum7-d(iir,iig)*sf(jn)*fact

 10   continue
  8   continue
  4   continue
  3   continue

      do 11 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c-------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c-------------------------------------------------------------

      if(irsym.ne.issym)go to 11
      if(orbp.ne.orbr.and.ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs.and.iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)

 11   continue

      do 12 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 12

      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(iqssym.ne.irgsym)go to 14

c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact

 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact

 18   continue
 17   continue
 14   continue

      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.and.iqsym.ne.issym)go to 20

      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c----------------------------------------------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact

 22   continue
 21   continue

c------------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.ne.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

  42  continue

c-----------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c-----------------------------------------------------------------

      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact

 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact

 48   continue
 47   continue
 46   continue

c----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact

 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)
      

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)

      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52

      if(igsym.ne.idsym)go to 52

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)

c-----------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c-----------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbp+orbq+irank+2*orbg)
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact

  54  continue

c------------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c------------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact

  55  continue
  53  continue
  52  continue
  30  continue


      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)


      hyp1=0

      if(ipsym.eq.issym)then
      call quadB(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif

      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      if(iqssym.eq.irasym)then

c-----------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c-----------------------------------------------------------------


      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif

      hyp2=0

      if(issym.eq.iqsym)then
      call quadB(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)
      
      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      if(ipssym.eq.irasym)then

c----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c-----------------------------------------------------------------


      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      if(irasym.eq.iqssym)then
      
      hyp1=0

      if(itsym.eq.issym)then
      call quadB(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif

      endif
      endif
      endif
      endif
      endif


      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      irpsym=mtbl(irsym,ipsym)
      itasym=mtbl(itsym,iasym)
      if(irpsym.eq.itasym)then

      hyp1=0

      if(itsym.eq.issym)then
      call quadB(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)

c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif

      endif
      endif
      endif
      endif
      endif

      enddo
      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

      hyp1=0

      if(ipsym.eq.ibsym)then
      call quadB(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.iqssym)then
      if(irasym.eq.isbsym)then
      if(orbq.eq.orbb)then

c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c-----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      hyp2=0

      if(iqsym.eq.ibsym)then
      call quadB(orbb,orbq,iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)

      if(irasym.eq.ipssym)then
      if(irasym.eq.isbsym)then
      if(orbp.eq.orbb)then

c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      hyp=0

      if(iasym.ne.igsym)go to 610
     
      call quadB(orbg,orba,iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
c     if(orbb.ne.orbq)go to 620

      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 27 = <rs|t2|ab><a|d|g><b|t1|p><sg|s2|qr>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+irank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact

 640  continue
      endif

c-------------------------------------------------------------------
c     t2dt1s2     diagram 28 = <rs|t2|ba><a|d|g><b|t1|p><sg|s2|qr>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact

  650 continue
      endif

  630 continue
  620 continue

      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
c     if(orbb.ne.orbp)go to 660

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2dt1s2     diagram 29 = <pa|s2|sr><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+irank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact

  680 continue
      endif

c-------------------------------------------------------------------
c     s2t1dt2     diagram 30 = <pa|s2|rs><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------


      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact

  690 continue
      endif

  670 continue
  660 continue
  610 continue

      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.ne.isgsym)go to 605
      if(irasym.ne.iqssym)go to 605

      if(ibsym.ne.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
      if(orbp.ne.orbb)go to 605

      in=idpa(iip,iib)

      hyp=0

      call quadB(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700

c-------------------------------------------------------------------
c     t2dt1s2     diagram 31 = <rs|t2|ag><g|d|b><p|t1|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

c-------------------------------------------------------------------
c     t2dt1s2     diagram 32 = <rs|t2|ga><g|d|b><p|t1|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact

  720 continue
  710 continue
  700 continue
  605 continue

      hyp=0

      if(ibsym.eq.igsym)then
      call quadB(orbg,orbb,iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)

      if(ipssym.ne.irasym)go to 600
      if(isgsym.ne.irasym)go to 600

      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 33 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ag>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

c-------------------------------------------------------------------
c     s2t1dt2     diagram 34 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbp,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbp+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact

  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0

      if(itsym.eq.issym)then
      call quadB(orbs,orbt,iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      isbsym=mtbl(issym,ibsym)
      irasym=mtbl(irsym,iasym)
      iqtsym=mtbl(iqsym,itsym)

      if(isbsym.ne.irasym)go to 815
      if(iqtsym.ne.irasym)go to 815

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 35 = <rs|t2|ab><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)

      call findk(orbs,orbb,orbr,orba,ias,iab,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 812
      do 811 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      if(l2.eq.l1)then
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+irank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

      endif

  811 continue
  812 continue

c-------------------------------------------------------------------
c     t2dt1s2     diagram 45 = <rs|t2|ba><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbs,orba,orbr,orbb,ias,iaa,iar,iab,kgot,kmax)
      if(kmax.eq.0)go to 814
      do 813 kloop=1,kmax
      l1=kgot(kloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrb=idra(iir,iib,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrb,isasym,l1)
      call sixj(orbp,rank,orbq,orbt,al2,orbs,s6j1)
      call sixj(orbp,al2,orbs,orba,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbt+orbs+l2+orbr+irank+orba)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum45=sum45+hyp*tf(inf)*si(ini)*ti(in)*fact

  813 continue
  814 continue 
  830 continue
  820 continue
  815 continue

      ipssym=mtbl(ipsym,issym)
      irasym=mtbl(irsym,iasym)
      itbsym=mtbl(itsym,ibsym)

      if(ipssym.ne.irasym)go to 900
      if(itbsym.ne.irasym)go to 900

      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 36 = <pr|s2|sa><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)

      call findk(orbp,orbs,orbr,orba,iap,ias,iar,iaa,jgot,jmax)
      if(jmax.eq.0)go to 912
      do 911 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      if(l2.eq.l1)then
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+irank+orba+orbr+al2)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact

      endif

  911 continue
  912 continue

c-------------------------------------------------------------------
c     s2t1dt2     diagram 46 = <ps|s2|ra><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------

      call findk(orbp,orbr,orbs,orba,iap,iar,ias,iaa,kgot,kmax)
      if(kmax.eq.0)go to 914
      do 913 kloop=1,kmax
      l1=kgot(kloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      nsa=idra(iis,iia,l1)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l1)
      call sixj(orbp,rank,orbq,orbt,al2,orbs,s6j1)
      call sixj(orbp,al2,orbs,orba,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbt+orbs+orbr+orba+l2+irank)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum46=sum46+hyp*ti(ini)*sf(inf)*tf(in)*fact

  913 continue
  914 continue      
  930 continue
  900 continue
  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)

      hyp=0

      if(itsym.eq.ibsym)then
      call quadB(orbt,orbb,iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      isqsym=mtbl(issym,iqsym)

      if(irasym.eq.isbsym)then
      if(irasym.eq.isqsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

c-------------------------------------------------------------------
c     s1t2ds2     diagram 41 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nqs=jdra(iiq,iis,l1)
      nra=idra(iir,iia,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t2ds2     diagram 42 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|qa>
c-------------------------------------------------------------------


      call findk(orba,orbs,orbq,orbr,iaa,ias,iaq,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif

      enddo
      endif

      endif
      endif

      endif
      endif
      endif
      endif


      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then
      in=jdpa(iiq,iit)

      hyp=0
      if(itsym.eq.ibsym)then
      call quadB(orbb,orbt,iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      ispsym=mtbl(issym,ipsym)

      if(irasym.eq.isbsym)then
      if(irasym.eq.ispsym)then

c--------------------------------------------------------------------
c     s2dt1s1     diagram 43 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ba>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

c-------------------------------------------------------------------
c     s2dt1s1     diagram 44 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ab>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue

   1  continue
 100  continue

      sum1=sum1*con
      sum2=sum2*con
      sum3=sum3*con
      sum4=sum4*con
      sum5=sum5*con
      sum6=sum6*con
      sum7=sum7*con
      sum8=sum8*con
      sum9=sum9*con
      sum10=sum10*con
      sum11=sum11*con
      sum12=sum12*con
      sum13=sum13*con
      sum14=sum14*con
      sum15=sum15*con
      sum16=sum16*con
      sum17=sum17*con
      sum18=sum18*con
      sum19=sum19*con
      sum20=sum20*con
      sum21=sum21*con
      sum22=sum22*con
      sum23=sum23*con
      sum24=sum24*con
      sum25=sum25*con
      sum26=sum26*con
      sum27=sum27*con
      sum28=sum28*con
      sum29=sum29*con
      sum30=sum30*con
      sum31=sum31*con
      sum32=sum32*con
      sum33=sum33*con
      sum34=sum34*con
      sum35=sum35*con
      sum36=sum36*con
      sum37=sum37*con
      sum38=sum38*con
      sum39=sum39*con
      sum40=sum40*con
      sum41=sum41*con
      sum42=sum42*con
      sum43=sum43*con
      sum44=sum44*con
      sum45=sum45*con
      sum46=sum46*con

      print*,'Dirac - Fock value hyp(B) ',dirac*con

      return
      end

c********************************************************************
c    This subroitine calculates the M1 transition property diagram  *
c********************************************************************


      subroutine m1transprop(d,si,sf,ti,tf,ival,jval,jinit,jfin)

      implicit real*8 (a-h,o-z)
      real*8 jinit,jfin

      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)

      common/kappa/nka(MNBAS)

      common/valuem1/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS)
      dimension si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0

      irank=1
      rank=1.0d0

c-----p is the final and q is the initial state----------------

      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)

      if(jfin.ne.orbp)go to 100
      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)

      if(jinit.ne.orbq)go to 1
      if(iiq.ne.ival)go to 1

      if(ipsym.ne.iqsym) go to 1

c---------------------------------------------------------------
c     d      diagram 1 = <p|d|q>
c---------------------------------------------------------------

      sum=sum+d(iip,iiq)

      call m1prop(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,dirac)
      sum1=sum1+d(iip,iiq)


      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(iir.ne.iiq)then
      if(irsym.eq.ipsym)then
      if(irsym.eq.iqsym.and.orbq.eq.orbr)then

      in=jdpa(iiq,iir)

c------------------------------------------------------------
c     ds1     diagram 2= <p|d|r><r|s1|q>
c-----------------------------------------------------------

      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)

      endif
      endif
      endif

      if(irsym.eq.iqsym)then
      if(irsym.eq.ipsym.and.orbr.eq.orbp)then
      if(iir.ne.iip)then

      in=kdpa(iip,iir)

c-----------------------------------------------------------
c     s1d     diagram 3 = <p|s1|r><r|d|q>
c-----------------------------------------------------------

      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)

      endif
      endif
      endif

 2    continue

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 4

      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)

      if(ipqsym.ne.irgsym)go to 4

c-----------------------------------------------------------------
c     ds2     diagram 4 = <g|d|r><pr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact

  6   continue

      endif

c------------------------------------------------------------------
c     ds2     diagram 5 = <g|d|r><pr|s2|gq>
c------------------------------------------------------------------

      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact

  7   continue
  5   continue

c-----------------------------------------------------------------
c     s2d     diagram 6 = <r|d|g><pg|s2|qr>
c-----------------------------------------------------------------

      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      jn=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(jn)*fact
      sum6=sum6+d(iir,iig)*sf(jn)*fact

  9   continue

      endif

c------------------------------------------------------------------
c     s2d     diagram 7 = <r|d|g><pg|s2|rq>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      jn=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(jn)*fact
      sum7=sum7-d(iir,iig)*sf(jn)*fact

 10   continue
  8   continue
  4   continue
  3   continue

      do 11 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

c-------------------------------------------------------------
c     s1ds1     diagram 8 = <p|s1|r><r|d|s><s|s1|q>
c-------------------------------------------------------------

      if(irsym.ne.issym)go to 11
      if(orbp.ne.orbr.and.ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs.and.iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11

      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)

      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)

 11   continue

      do 12 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      if(irsym.ne.igsym)go to 12

      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14

      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(iqssym.ne.irgsym)go to 14

c-----------------------------------------------------------------
c     s1ds2     diagram 9 = <p|s1|s><g|d|r><sr|s2|qg>
c-----------------------------------------------------------------

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact

 16   continue
 15   continue

c-----------------------------------------------------------------
c     s1ds2     diagram 10 = <p|s1|s><g|d|r><sr|s2|gq>
c-----------------------------------------------------------------

      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact

 18   continue
 17   continue
 14   continue

      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.and.iqsym.ne.issym)go to 20

      ini=jdpa(iiq,iis)

      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)

      if(ipssym.ne.irgsym)go to 20

c----------------------------------------------------------------
c     s2ds1     diagram 11 = <pg|s2|sr><r|d|g><s|s1|q>
c----------------------------------------------------------------

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact

 22   continue
 21   continue

c------------------------------------------------------------------
c     s2ds1     diagram 12 = <pg|s2|rs><r|d|g><s|t1|q>
c------------------------------------------------------------------

      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.ne.itsym)go to 40

      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)

      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

c----------------------------------------------------------------
c     s2ds2     diagram 13 = <pg|s2|tr><t|d|s><sr|s2|qg>
c----------------------------------------------------------------

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

  42  continue

c-----------------------------------------------------------------
c     s2ds2     diagram 14 = <pg|s2|rt><t|d|s><sr|s2|qg>
c-----------------------------------------------------------------

      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact

 45   continue
 44   continue
 43   continue

c----------------------------------------------------------------
c     s2ds2     diagram 15 = <pg|s2|rt><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)

      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact

 48   continue
 47   continue
 46   continue

c----------------------------------------------------------------
c     s2ds2     diagram 16 = <pg|s2|tr><t|d|s><rg|s2|qs>
c----------------------------------------------------------------

      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact

 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)
      kapd=nka(iid)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)

      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52
      if(igsym.ne.idsym)go to 52

      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)

c-----------------------------------------------------------------
c     s2ds2     diagram 17 = <pr|s2|sg><g|d|d><sr|s2|qd>
c-----------------------------------------------------------------

      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbp+orbq+irank+2*orbg)
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact

  54  continue

c------------------------------------------------------------------
c     s2ds2     diagram 18 = <pr|s2|sg><g|d|d><sr|s2|dq>
c------------------------------------------------------------------

      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact

  55  continue
  53  continue
  52  continue
  30  continue

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)


      hyp1=0

      if(ipsym.eq.issym)then
      call m1prop(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif


      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      if(iqssym.eq.irasym)then

c-----------------------------------------------------------------
c     t1ds2     diagram 19 = <r|t1|a><p|d|s><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     t1ds2     diagram 20 = <r|t1|a><p|d|s><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif

      hyp2=0

      if(issym.eq.iqsym)then
      call m1prop(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)

      irasym=mtbl(irsym,iasym)
      ipssym=mtbl(ipsym,issym)
      if(ipssym.eq.irasym)then

c----------------------------------------------------------------
c     s2dt1     diagram 21 = <pa|s2|sr><s|d|q><r|t1|a>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s2dt1     diagram 22 = <pa|s2|rs><s|d|q><r|t1|a>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

      hyp1=0

      if(itsym.eq.issym)then
      call m1prop(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

c-----------------------------------------------------------------
c     s1t1ds2     diagram 37 = <p|s1|t><t|d|s><r|t1|a><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact

      endif
      enddo
      endif

c-----------------------------------------------------------------
c     s1t1ds2     diagram 38 = <p|s1|t><t|d|s><r|t1|a><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0

      if(itsym.eq.issym)then
      call m1prop(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)

c-----------------------------------------------------------------
c     s2dt1s1     diagram 39 = <pa|s2|tr><t|d|s><r|t1|a><s|s1|q>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact

      endif
      enddo
      endif

c------------------------------------------------------------------
c     s2dt1s1      diagram 40 = <pa|s2|rt><t|d|s><r|t1|a><s|s1|q>
c------------------------------------------------------------------


      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact

      enddo
      endif
      endif
      endif
      endif
      endif

      enddo
      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

      hyp1=0

      if(ipsym.eq.ibsym)then
      call m1prop(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)
      isbsym=mtbl(issym,ibsym)

      if(irasym.eq.iqssym)then
      if(irasym.eq.isbsym)then
      if(orbq.eq.orbb)then

c-----------------------------------------------------------------
c     t2ds2      diagram 23 = <br|t2|sa><p|d|b><sr|s2|qa>
c-----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact

c-----------------------------------------------------------------
c     t2ds2     diagram 24 = <br|s2|sa><p|d|b><sr|s2|aq>
c-----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      hyp2=0

      if(iqsym.eq.ibsym)then
      call m1prop(orbb,orbq,iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)

      if(irasym.eq.ipssym)then
      if(irasym.eq.isbsym)then
      if(orbp.eq.orbb)then

c----------------------------------------------------------------
c     s2dt2     diagram 25 = <pr|s2|sa><b|d|q><rs|t2|ab>
c----------------------------------------------------------------

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact

c----------------------------------------------------------------
c     s2dt2     diagram 26 = <pr|s2|sa><b|d|q><rs|t2|ba>
c----------------------------------------------------------------

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.ne.igsym)go to 610

      hyp=0

      call m1prop(orbg,orba,iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)


      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620
      in=idpa(iip,iib)


      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 27 = <rs|t2|ab><a|d|g><b|t1|p><sg|s2|qr>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact

 640  continue

      endif

c-------------------------------------------------------------------
c     t2dt1s2     diagram 28 = <rs|t2|ba><a|d|g><b|t1|p><sg|s2|qr>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact

  650 continue

      endif

  630 continue
  620 continue

      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2dt1s2     diagram 29 = <pa|s2|sr><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact

  680 continue

      endif

c-------------------------------------------------------------------
c     s2t1dt2     diagram 30 = <pa|s2|rs><b|t1|q><a|d|g><rs|t2|gb>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact

  690 continue

      endif

  670 continue
  660 continue
  610 continue


      if(ibsym.ne.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
      if(orbp.ne.orbb)go to 605

      in=idpa(iip,iib)

      hyp=0

      call m1prop(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700

c-------------------------------------------------------------------
c     t2dt1s2     diagram 31 = <rs|t2|ag><g|d|b><p|t1|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      if(irasym.eq.isgsym)then
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

      endif

c-------------------------------------------------------------------
c     t2dt1s2     diagram 32 = <rs|t2|ga><g|d|b><p|t1|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      irgsym=mtbl(irsym,igsym)
      inf=ntloc(nsa,nrg,isasym,l1)
      if(irgsym.eq.isasym)then
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact

      endif

  720 continue
  710 continue
  700 continue
  605 continue

      hyp=0

      if(ibsym.eq.igsym)then
      call m1prop(orbg,orbb,iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600

      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 33 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ag>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

c-------------------------------------------------------------------
c     s2t1dt2     diagram 34 = <pr|s2|sa><b|t1|q><b|d|g><rs|t2|ga>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact

  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0

      if(itsym.eq.issym)then
      call m1prop(orbs,orbt,iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)

c-------------------------------------------------------------------
c     t2dt1s2     diagram 35 = <rs|t2|ab><s|d|t><p|t1|b><rt|s2|aq>
c-------------------------------------------------------------------

      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

  830 continue
  820 continue
  815 continue

      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900
      in=idpa(iiq,iib)

c-------------------------------------------------------------------
c     s2t1dt2     diagram 36 = <pr|s2|sa><q|t1|b><s|d|t><rt|t2|ab>
c-------------------------------------------------------------------

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact

  930 continue
  900 continue
  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)

      hyp=0

      if(itsym.eq.ibsym)then
      call m1prop(orbt,orbb,iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

c-------------------------------------------------------------------
c     s1t2ds2     diagram 41 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|aq>
c-------------------------------------------------------------------

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nqs=jdra(iiq,iis,l1)
      nra=idra(iir,iia,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact

      endif
      enddo
      endif

c-------------------------------------------------------------------
c     s1t2ds2     diagram 42 = <p|s1|t><rs|t2|ab><t|d|b><rs|s2|qa>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbq,orbr,iaa,ias,iaq,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif
      enddo
      endif
      endif

      endif
      endif
      endif
      endif


      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then

      in=jdpa(iiq,iit)

      hyp=0

      if(itsym.eq.ibsym)then
      call m1prop(orbb,orbt,iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then

c--------------------------------------------------------------------
c     s2dt1s1     diagram 43 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ba>
c--------------------------------------------------------------------

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

c-------------------------------------------------------------------
c     s2dt1s1     diagram 44 = <pr|s2|sa><t|d|b><t|s1|q><sr|t2|ab>
c-------------------------------------------------------------------

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue
   1  continue
 100  continue

      print*,'Dirac - Fock value M1 trns. amplitude ',dirac

      return
      end

      subroutine m2transprop(d,si,sf,ti,tf,ival,jval,jinit,jfin)
      implicit real*8 (a-h,o-z)
      real*8 jinit,jfin,mj
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

c
      common/valuem2/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

c     local dimension
c
      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS),si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      if(jinit.le.jfin)then
      mj=jinit
      else
      mj=jfin
      endif
      d11=dr(jfin,two,jinit,-mj,0.0d0,mj)

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0

      irank=2
      rank=2.0d0

c     p is the final and q is the initial state
      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)
      if(iip.ne.jval)go to 100
      if(orbp.ne.jfin) goto 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
      if(orbq.ne.jinit) goto 1
      if(iiq.ne.ival)go to 1
      if(ipsym.eq.iqsym) go to 1

c     diagram 1 =<f|d|i>
      call m2prop(iip,iiq,orbp,orbq,kapp,kapq,hyp1)
      print*,'Dirac-Fock Value M2', hyp1
      sum=sum+d(iip,iiq)
      sum1=sum1+d(iip,iiq)
c     diagram < | ds|>
      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      if(iir.ne.ival)then
      if(irsym.ne.ipsym)then
      if(orbq.eq.orbr)then
      if(irsym.eq.iqsym)then
      if(iir.ne.iiq)then
      in=jdpa(iiq,iir)
c     diagram 2(a)=  sum (r) <p|d|r><r|s|q>
      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)
      endif
      endif
      endif
      endif

      if(irsym.ne.iqsym)then
      if(orbr.eq.orbp)then
      if(irsym.eq.ipsym)then
      if(iip.ne.iir)then
      in=kdpa(iip,iir)
c     diagram 3(a) =  sum (r) <p|s|r><r|d|q>
      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)
      endif
      endif
      endif
      endif
      endif

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      if(irsym.eq.igsym)go to 4
      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)
      if(ipqsym.ne.irgsym)go to 4
c     diagram 2 (b) = sum(g,r) <g|d|r><pr|s|qg>
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact
  6   continue
      endif
c     diagram 2 (c) =-sum(g,r) <g|d|r><pr|s|gq>
      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact
  7   continue
  5   continue
c     diagram 3 (b)
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      in=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(in)*fact
      sum6=sum6+d(iir,iig)*sf(in)*fact
  9   continue
      endif
c     diagram 3 (c)
      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      in=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(in)*fact
      sum7=sum7-d(iir,iig)*sf(in)*fact
 10   continue
  8   continue
  4   continue
  3   continue
c     diagram 4e
      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      if(issym.eq.irsym)go to 11
      if(orbp.ne.orbr)go to 11
      if(ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs)go to 11
      if(iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11
      in2=kdpa(iip,iir)
      in1=jdpa(iiq,iis)
      sum=sum+d(iir,iis)*si(in1)*sf(in2)
      sum8=sum8+d(iir,iis)*si(in1)*sf(in2)
 11   continue

c     diagram 4a
      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)
      if(irsym.eq.igsym)go to 12
      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14
      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)
      if(iqssym.ne.irgsym)go to 14

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact
 16   continue
 15   continue
c     diagram 4b
      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact
 18   continue
 17   continue
 14   continue
c     diagram 4c
      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.or.iqsym.ne.issym)go to 20
      ini=jdpa(iiq,iis)
      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)
      if(ipssym.ne.irgsym)go to 20
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact
 22   continue
 21   continue
c     diagram 4d
      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact

 24   continue
 23   continue
 20   continue
 13   continue
 12   continue
  2   continue
      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)
      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)
      if(issym.eq.itsym)go to 40
      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)
      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41

      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+rank)
      fact=s6j*isign/(2*l1+1)
c     diagram 4f
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact
      endif
      enddo
      endif
  42  continue

c     diagram 4(g)
      itgsym=mtbl(itsym,igsym)
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(al2,orbp,orbt,al1,orbg,orbr,s6j1)
      call sixj(al2,orbt,orbp,rank,orbq,orbs,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+rank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact
 45   continue
 44   continue
 43   continue
c     diagram 4(h)
      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(al1,al2,rank,orbq,orbp,orbr,s6j1)
      call sixj(al1,al2,rank,orbs,orbt,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact
 48   continue
 47   continue
 46   continue
c     diagram 4(i)
      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)
      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+rank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact
 51   continue
 50   continue
 49   continue
 41   continue
 40   continue
c     diagram 4j and 4k
      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)
      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52

      if(igsym.eq.idsym)go to 52
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)
      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      irdsym=mtbl(irsym,idsym)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(al2,al1,rank,orbp,orbq,orbs,s6j1)
      call sixj(al1,al2,rank,orbd,orbg,orbr,s6j2)
      fact=(-1)**(orbp+orbq+rank+2*orbg)*s6j1*s6j2
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact
  54  continue
      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2+2*(orbd+orbr))*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact
  55  continue
  53  continue
  52  continue
  30  continue

C     Two-body effective operator (computed directly)

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)

      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)

      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

c     diagrams  a1-a4 (sum19-sum22)

      hyp1=0
      if(ipsym.ne.issym)then
      call m2prop(iip,iis,orbp,orbs,kapp,kaps,hyp1)
      endif

      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact
      enddo
      endif

      endif

      hyp2=0
      if(issym.ne.iqsym)then
      call m2prop(iis,iiq,orbs,orbq,kaps,kapq,hyp2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)
      ipssym=mtbl(ipsym,issym)

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact
      enddo
      endif

      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

c     diagram 43
      hyp1=0
      if(itsym.ne.issym)then
      call m2prop(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 45
      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0
      if(itsym.ne.issym)then
      call m2prop(iit,iis,orbt,orbs,kapt,kaps,hyp1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)
      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact
      endif
      enddo
      endif

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

      enddo

      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

c     diagrams  b1-b4 (sum19-sum22)
      hyp1=0
      if(ipsym.ne.ibsym)then
      call m2prop(iip,iib,orbp,orbb,kapp,kapb,hyp1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.eq.iqssym)then
      if(orbq.eq.orbb)then

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact
      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact
      enddo
      endif
      enddo
      endif

      endif
      endif

      hyp2=0
      if(iqsym.ne.ibsym)then
      call m2prop(iib,iiq,orbb,orbq,kapb,kapq,hyp2)
      endif

      ipssym=mtbl(ipsym,issym)
      if(irasym.eq.ipssym)then
      if(orbp.eq.orbb)then

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact
      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact
      enddo
      endif
      enddo
      endif

      endif
      endif

      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.eq.igsym)go to 610
      hyp=0
c     c1-c2
      call m2prop(iig,iia,orbg,orba,kapg,kapa,hyp)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620
      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact
 640  continue
      endif

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact
  650 continue
      endif

  630 continue
  620 continue

c     d1-d2
      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660
      in=idpa(iiq,iib)

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact
  680 continue
      endif

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact
  690 continue
      endif

  670 continue
  660 continue
  610 continue

      if(ibsym.eq.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
c     if(orbp.ne.orbg)go to 605
      if(orbp.ne.orbb)go to 605
      in=idpa(iip,iib)
c     c7-c8
      hyp=0
      call m2prop(iib,iig,orbb,orbg,kapb,kapg,hyp)
      if(orbq.ne.orbg)go to 700
      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact

  720 continue
  710 continue
  605 continue
  700 continue
c     d7-d8

      hyp=0
      if(ibsym.ne.igsym)then
      call m2prop(iig,iib,orbg,orbb,kapg,kapb,hyp)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600
      in=idpa(iiq,iib)
      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact
  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0
      if(itsym.ne.issym)then
      call m2prop(iis,iit,orbs,orbt,kaps,kapt,hyp)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)
c     c3
      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact

  830 continue
  820 continue
  815 continue

c     d3
      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900
      in=idpa(iiq,iib)

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact
  930 continue
  900 continue

  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)
      hyp=0
      if(itsym.ne.ibsym)then
      call m2prop(iit,iib,orbt,orbb,kapt,kapb,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nra=idra(iir,iia,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact
      endif
      enddo
      endif

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)
      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact
      enddo
      endif
      enddo
      endif

      endif

      endif
      endif
      endif
      endif

      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then
      in=jdpa(iiq,iit)

      hyp=0
      if(itsym.ne.ibsym)then
      call m2prop(iib,iit,orbb,orbt,kapb,kapt,hyp)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact
      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue

   1  continue
 100  continue
      return
      end

      subroutine gfactprop(d,dc,si,sf,ti,tf,ival,jval,jinit,jfin)
      implicit real*8 (a-h,o-z)
      real*8 jorbv,jfin,jinit
      common/jvalue/orbc(MNOCC),orbe(MNEXC)
      common/index/ke(MNEXC),kc(MNOCC)
      common/info/nocc,nexcit
      common/parity/isymc(MNOCC),isyme(MNEXC)
      common/symmetry/mtbl(MNS,MNS)
      common/kpavlaue/iqc(MNOCC),iqe(MNEXC)
      common/orbital_energy/eorb(MNBAS)
      common/idra1/idra(MNBAS,MNBAS,0:MXV)
      common/idpa1/idpa(MNBAS,MNBAS)
c     for initial state
      common/idra2/jdra(MNBAS,MNBAS,0:MXV)
      common/idpa2/jdpa(MNBAS,MNBAS)
c     for final state
      common/idra3/kdra(MNBAS,MNBAS,0:MXV)
      common/idpa3/kdpa(MNBAS,MNBAS)
      common/kappa/nka(MNBAS)

      common/valueg/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44
      common/valuegc/sumc1,sumc2,sumc3,sumc4,sumc5,sumc6,sumc7,
     :sumc8,sumc9,sumc10,sumc11,sumc12,sumc13,sumc14,sumc15,sumc16,
     :sumc17,sumc18,sumc19,sumc20,sumc21,sumc22,sumc23,sumc24,sumc25,
     :sumc26,sumc27,sumc28,sumc29,sumc30,sumc31,sumc32,sumc33,sumc34,
     :sumc35,sumc36,sumc37,sumc38,sumc39,sumc40,sumc41,sumc42,sumc43,
     :sumc44

      dimension kgot(MXV+1),igot(MXV+1),jgot(MXV+1)
      dimension d(MNBAS,MNBAS),dc(MNBAS,MNBAS)
      dimension si(MDIM),sf(MDIM)
      dimension ti(MDIM),tf(MDIM)

      data zero,half,one,two/0.0,0.5,1.0,2.0/
      data VL/137.03598950D+00/

      jorbv=jinit
      d1=dr(jorbv,1.0d0,jorbv,-jorbv,0.0d0,jorbv)

      con=1.0d0
      con=con*d1

      do i=0,MXV
      kgot(i)=0
      igot(i)=0
      jgot(i)=0
      enddo

      sum=0
      sum1=0
      sum2=0
      sum3=0
      sum4=0
      sum5=0
      sum6=0
      sum7=0
      sum8=0
      sum9=0
      sum10=0
      sum11=0
      sum12=0
      sum13=0
      sum14=0
      sum15=0
      sum16=0
      sum17=0
      sum18=0
      sum19=0
      sum20=0
      sum21=0
      sum22=0
      sum23=0
      sum24=0
      sum25=0
      sum26=0
      sum27=0
      sum28=0
      sum29=0
      sum30=0
      sum31=0
      sum32=0
      sum33=0
      sum34=0
      sum35=0
      sum36=0
      sum37=0
      sum38=0
      sum39=0
      sum40=0
      sum41=0
      sum42=0
      sum43=0
      sum44=0

      sumc=0
      sumc1=0
      sumc2=0
      sumc3=0
      sumc4=0
      sumc5=0
      sumc6=0
      sumc7=0
      sumc8=0
      sumc9=0
      sumc10=0
      sumc11=0
      sumc12=0
      sumc13=0
      sumc14=0
      sumc15=0
      sumc16=0
      sumc17=0
      sumc18=0
      sumc19=0
      sumc20=0
      sumc21=0
      sumc22=0
      sumc23=0
      sumc24=0
      sumc25=0
      sumc26=0
      sumc27=0
      sumc28=0
      sumc29=0
      sumc30=0
      sumc31=0
      sumc32=0
      sumc33=0
      sumc34=0
      sumc35=0
      sumc36=0
      sumc37=0
      sumc38=0
      sumc39=0
      sumc40=0
      sumc41=0
      sumc42=0
      sumc43=0
      sumc44=0

      irank=1
      rank=1.0d0

c     p is the final and q is the initial state
      do 100 ip=1,nexcit
      ipsym=isyme(ip)
      orbp=orbe(ip)
      iap=iqe(ip)
      iip=ke(ip)
      kapp=nka(iip)
      if(jfin.ne.orbp)go to 100
      if(iip.ne.jval)go to 100

      do 1 iq=1,nexcit
      iqsym=isyme(iq)
      orbq=orbe(iq)
      iaq=iqe(iq)
      iiq=ke(iq)
      kapq=nka(iiq)
      if(jorbv.ne.orbq)go to 1
      if(iiq.ne.ival)go to 1
      if(ipsym.ne.iqsym) go to 1
c     diagram 1 =<f|d|i>
      sum=sum+d(iip,iiq)
      call gfact(orbp,orbq,iip,iiq,orbp,orbq,kapp,kapq,
     :dirac,diracc)
      sum1=sum1+d(iip,iiq)
      sumc1=sumc1+dc(iip,iiq)
c     diagram < | ds|>

      do 2 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      if(iir.ne.iiq)then
      if(irsym.eq.ipsym)then
      if(irsym.eq.iqsym.and.orbq.eq.orbr)then
      in=jdpa(iiq,iir)
c     diagram 2(a)=  sum (r) <p|d|r><r|s|q>
      sum=sum+d(iip,iir)*si(in)
      sum2=sum2+d(iip,iir)*si(in)
      sumc2=sumc2+dc(iip,iir)*si(in)
      endif
      endif
      endif

      if(irsym.eq.iqsym)then
      if(irsym.eq.ipsym.and.orbr.eq.orbp)then
      if(iir.ne.iip)then
      in=kdpa(iip,iir)
c     diagram 3(a) =  sum (r) <p|s|r><r|d|q>
      sum=sum+d(iir,iiq)*sf(in)
      sum3=sum3+d(iir,iiq)*sf(in)
      sumc3=sumc3+dc(iir,iiq)*sf(in)
      endif
      endif
      endif

 2    continue

      do 3 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)

      do 3 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iag=iqc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      if(irsym.ne.igsym)go to 4
      ipqsym=mtbl(ipsym,iqsym)
      irgsym=mtbl(irsym,igsym)
      if(ipqsym.ne.irgsym)go to 4
c     diagram 2 (b) = sum(g,r) <g|d|r><pr|s|qg>
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
c     diagram 2 (c) =-sum(g,r) <g|d|r><pr|s|gq>
      if(imax.ne.0)then
      do 6 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 6
      nrg=idra(iir,iig,l1)
      nqp=jdra(iiq,iip,l1)
      in=nsloci(nqp,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iig,iir)*si(in)*fact
      sum4=sum4+d(iig,iir)*si(in)*fact
      sumc4=sumc4+dc(iig,iir)*si(in)*fact
  6   continue
      endif
c     diagram 2 (c)
      call findk(orbp,orbg,orbq,orbr,iap,iag,iaq,iar,igot,imax)
      if(imax.eq.0)go to 5
      do 7 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npg=idra(iip,iig,l1)
      nqr=jdra(iiq,iir,l1)
      ipgsym=mtbl(ipsym,igsym)
      in=nsloci(nqr,npg,ipgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(in)*fact
      sum5=sum5-d(iig,iir)*si(in)*fact
      sumc5=sumc5-dc(iig,iir)*si(in)*fact
  7   continue
  5   continue
c     diagram 3 (b)
      call findk(orbp,orbq,orbg,orbr,iap,iaq,iag,iar,igot,imax)
      if(imax.ne.0)then
      do 9 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      if(l1.ne.irank)go to 9
      nrg=idra(iir,iig,l1)
      npq=kdra(iip,iiq,l1)
      jn=nslocf(npq,nrg,ipqsym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+l1)/(2*l1+1)
      sum=sum+d(iir,iig)*sf(jn)*fact
      sum6=sum6+d(iir,iig)*sf(jn)*fact
      sumc6=sumc6+dc(iir,iig)*sf(jn)*fact
  9   continue
      endif
c     diagram 3 (c)
      call findk(orbp,orbr,orbq,orbg,iap,iar,iaq,iag,igot,imax)
      if(imax.eq.0)go to 8
      do 10 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nqg=idra(iiq,iig,l1)
      npr=kdra(iip,iir,l1)
      iqgsym=mtbl(iqsym,igsym)
      jn=nslocf(npr,nqg,iqgsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*sf(jn)*fact
      sum7=sum7-d(iir,iig)*sf(jn)*fact
      sumc7=sumc7-dc(iir,iig)*sf(jn)*fact
 10   continue
  8   continue
  4   continue
  3   continue
c     diagram 4e
      do 11 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 11 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      if(orbr.ne.orbs)go to 11
      if(orbp.ne.orbr.or.ipsym.ne.irsym)go to 11
      if(orbq.ne.orbs.or.iqsym.ne.issym)go to 11
      if(iir.eq.iip)go to 11
      if(iis.eq.iiq)go to 11
      inf=kdpa(iip,iir)
      ini=jdpa(iiq,iis)
      sum=sum+d(iir,iis)*si(ini)*sf(inf)
      sum8=sum8+d(iir,iis)*si(ini)*sf(inf)
      sumc8=sumc8+dc(iir,iis)*si(ini)*sf(inf)
 11   continue

c     diagram 4a
      do 12 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 12 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)
      if(irsym.ne.igsym)go to 12
      do 13 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)
      if(iis.eq.iip)go to 14
      if(orbs.ne.orbp)go to 14
      if(ipsym.ne.issym)go to 14
      inf=kdpa(iip,iis)
      iqssym=mtbl(iqsym,issym)
      irgsym=mtbl(irsym,igsym)
      if(iqssym.ne.irgsym)go to 14
      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 15
      do 16 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 16
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iig,iir)*si(ini)*sf(inf)*fact
      sum9=sum9+d(iig,iir)*si(ini)*sf(inf)*fact
      sumc9=sumc9+dc(iig,iir)*si(ini)*sf(inf)*fact
 16   continue
 15   continue
c     diagram 4b
      call findk(orbq,orbr,orbg,orbs,iaq,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 17
      do 18 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsg,iqrsym,l1)
      call sixj(rank,orbr,orbg,al1,orbp,orbq,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iig,iir)*si(ini)*sf(inf)*fact
      sum10=sum10-d(iig,iir)*si(ini)*sf(inf)*fact
      sumc10=sumc10-dc(iig,iir)*si(ini)*sf(inf)*fact
 18   continue
 17   continue
 14   continue
      if(iis.eq.iiq)go to 20
      if(orbs.ne.orbq.or.iqsym.ne.issym)go to 20
      ini=jdpa(iiq,iis)
      ipssym=mtbl(ipsym,issym)
      irgsym=mtbl(irsym,igsym)
      if(ipssym.ne.irgsym)go to 20
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 21
      do 22 iloop=1,imax
      l1=igot(iloop)
      if(l1.ne.irank)go to 22
      nrg=idra(iir,iig,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nrg,ipssym,l1)
      fact=(-1)**(orbg+orbr+2*orbg+irank)/(2*l1+1)
      sum=sum+d(iir,iig)*si(ini)*sf(inf)*fact
      sum11=sum11+d(iir,iig)*si(ini)*sf(inf)*fact
      sumc11=sumc11+dc(iir,iig)*si(ini)*sf(inf)*fact
 22   continue
 21   continue
c     diagram 4d
      call findk(orbp,orbr,orbg,orbs,iap,iar,iag,ias,igot,imax)
      if(imax.eq.0)go to 23
      do 24 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsg=idra(iis,iig,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsg,iprsym,l1)
      call sixj(rank,orbr,orbg,al1,orbq,orbp,s6j)
      fact=(-1)**(orbg+orbr+irank)*s6j
      sum=sum-d(iir,iig)*si(ini)*sf(inf)*fact
      sum12=sum12-d(iir,iig)*si(ini)*sf(inf)*fact
      sumc12=sumc12-dc(iir,iig)*si(ini)*sf(inf)*fact
 24   continue
 23   continue
 20   continue
 13   continue
 12   continue

      do 30 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      iar=iqe(ir)
      kapr=nka(iir)
      do 30 ig=1,nocc
      iig=kc(ig)
      iag=iqc(ig)
      orbg=orbc(ig)
      igsym=isymc(ig)
      kapg=nka(iig)

      do 30 is=1,nexcit
      iis=ke(is)
      ias=iqe(is)
      orbs=orbe(is)
      issym=isyme(is)
      kaps=nka(iis)

      do 40 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      itsym=isyme(it)
      kapt=nka(iit)

      if(issym.ne.itsym)go to 40
      irgsym=mtbl(irsym,igsym)
      iptsym=mtbl(ipsym,itsym)
      iqssym=mtbl(iqsym,issym)
      if(iptsym.ne.irgsym)go to 41
      if(iqssym.ne.irgsym)go to 41
      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 41
      do 42 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)

      call findk(orbq,orbs,orbg,orbr,iaq,ias,iag,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      k1=jgot(jloop)
      if(k1.eq.l1)then
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,irgsym,l1)
      call sixj(rank,orbt,orbs,al1,orbq,orbp,s6j)
      isign=(-1)**(orbg+orbr+2*orbg+al1+orbt+orbs+irank)
      fact=s6j*isign/(2*l1+1)
c     diagram 4f
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum13=sum13+d(iit,iis)*si(ini)*sf(inf)*fact
      sumc13=sumc13+dc(iit,iis)*si(ini)*sf(inf)*fact
      endif
      enddo
      endif

  42  continue
c     diagram 4(g)
      itgsym=mtbl(itsym,igsym)
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 43
      do 44 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      call findk(orbs,orbq,orbr,orbg,ias,iaq,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 44
      do 45 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrg=idra(iir,iig,l2)
      ini=nsloci(nqs,nrg,irgsym,l2)
      call sixj(orbg,orbr,al2,orbp,orbt,al1,s6j1)
      call sixj(orbs,orbq,al2,orbp,orbt,rank,s6j2)
      isign=(-1)**(orbg+orbr+al2+orbt+orbs+irank)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum14=sum14-d(iit,iis)*si(ini)*sf(inf)*fact
      sumc14=sumc14-dc(iit,iis)*si(ini)*sf(inf)*fact
 45   continue
 44   continue
 43   continue
c     diagram 4(h)
      isgsym=mtbl(issym,igsym)
      itgsym=mtbl(itsym,igsym)
      call findk(orbp,orbr,orbg,orbt,iap,iar,iag,iat,igot,imax)
      if(imax.eq.0)go to 46
      do 47 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npr=kdra(iip,iir,l1)
      ntg=idra(iit,iig,l1)
      inf=nslocf(npr,ntg,itgsym,l1)
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 47
      do 48 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbr,s6j1)
      call sixj(orbt,orbs,rank,al2,al1,orbg,s6j2)
      isign=(-1)**(orbp+orbq+al1+al2+2*orbt)
      fact=s6j1*s6j2*isign
      sum=sum+d(iit,iis)*si(ini)*sf(inf)*fact
      sum15=sum15+d(iit,iis)*si(ini)*sf(inf)*fact
      sumc15=sumc15+dc(iit,iis)*si(ini)*sf(inf)*fact
 48   continue
 47   continue
 46   continue
c     diagram 4(i)
      isgsym=mtbl(issym,igsym)
      irgsym=mtbl(irsym,igsym)
      call findk(orbp,orbt,orbg,orbr,iap,iat,iag,iar,igot,imax)
      if(imax.eq.0)go to 49
      do 50 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      npt=kdra(iip,iit,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(npt,nrg,irgsym,l1)
      call findk(orbs,orbg,orbr,orbq,ias,iag,iar,iaq,jgot,jmax)
      if(jmax.eq.0)go to 50
      do 51 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsg=idra(iis,iig,l2)
      ini=nsloci(nqr,nsg,isgsym,l2)
      call sixj(orbt,orbp,al1,orbq,orbs,rank,s6j1)
      call sixj(orbg,orbr,al1,orbq,orbs,al2,s6j2)
      isign=(-1)**(orbt+orbs+irank+orbg+orbr+al1)
      fact=s6j1*s6j2*isign
      sum=sum-d(iit,iis)*si(ini)*sf(inf)*fact
      sum16=sum16-d(iit,iis)*si(ini)*sf(inf)*fact
      sumc16=sumc16-dc(iit,iis)*si(ini)*sf(inf)*fact
 51   continue
 50   continue
 49   continue
 41   continue
 40   continue

c     diagram 4(j)
      do 52 id=1,nocc
      idsym=isymc(id)
      orbd=orbc(id)
      iad=iqc(id)
      iid=kc(id)

      irgsym=mtbl(irsym,igsym)
      irdsym=mtbl(irsym,idsym)
      ipssym=mtbl(ipsym,issym)
      iqssym=mtbl(iqsym,issym)
      if(ipssym.ne.irgsym)go to 52
      if(iqssym.ne.irdsym)go to 52

      if(igsym.ne.idsym)go to 52
      call findk(orbp,orbs,orbg,orbr,iap,ias,iag,iar,igot,imax)
      if(imax.eq.0)go to 52
      do 53 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nps=kdra(iip,iis,l1)
      nrg=idra(iir,iig,l1)
      inf=nslocf(nps,nrg,irgsym,l1)
      call findk(orbs,orbq,orbd,orbr,ias,iaq,iad,iar,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 54 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqs=jdra(iiq,iis,l2)
      nrd=idra(iir,iid,l2)
      ini=nsloci(nqs,nrd,irdsym,l2)
      call sixj(orbp,orbq,rank,al2,al1,orbs,s6j1)
      call sixj(orbg,orbd,rank,al2,al1,orbr,s6j2)
      fact=s6j1*s6j2*(-1)**(orbp+orbq+irank+2*orbg)
      sum=sum-d(iid,iig)*si(ini)*sf(inf)*fact
      sum17=sum17-d(iid,iig)*si(ini)*sf(inf)*fact
      sumc17=sumc17-dc(iid,iig)*si(ini)*sf(inf)*fact
  54  continue
      call findk(orbq,orbr,orbs,orbd,iaq,iar,ias,iad,jgot,jmax)
      if(jmax.eq.0)go to 53
      do 55 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nqr=jdra(iiq,iir,l2)
      nsd=idra(iis,iid,l2)
      isdsym=mtbl(issym,idsym)
      ini=nsloci(nqr,nsd,isdsym,l2)
      call ninej(al1,orbr,orbg,orbp,orbq,rank,orbs,al2,orbd,s9j)
      fact=(-1)**(orbg+orbr+orbs+orbd+l1+l2)*s9j
      sum=sum+d(iid,iig)*si(ini)*sf(inf)*fact
      sum18=sum18+d(iid,iig)*si(ini)*sf(inf)*fact
      sumc18=sumc18+dc(iid,iig)*si(ini)*sf(inf)*fact
  55  continue
  53  continue
  52  continue
  30  continue

C     Two-body effective operator (computed directly)

      do 200 ia=1,nocc
      iasym=isymc(ia)
      orba=orbc(ia)
      iia=kc(ia)
      kapa=nka(iia)
      iaa=iqc(ia)
      do 200 ir=1,nexcit
      irsym=isyme(ir)
      orbr=orbe(ir)
      iir=ke(ir)
      kapr=nka(iir)
      iar=iqe(ir)
      do 200 is=1,nexcit
      issym=isyme(is)
      orbs=orbe(is)
      iis=ke(is)
      kaps=nka(iis)
      ias=iqe(is)

c     diagrams  a1-a4 (sum19-sum22)

      hyp1=0
      if(ipsym.eq.issym)then
      call gfact(orbp,orbs,iip,iis,orbp,orbs,kapp,kaps,hyp1,hypc1)
      endif

      if(orba.eq.orbr)then
      if(irsym.eq.iasym)then
      inf=idpa(iir,iia)

      if(orbq.eq.orbs)then

      iqssym=mtbl(iqsym,issym)
      call findk(orba,orbr,orbs,orbq,iaa,iar,ias,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum19=sum19+hyp1*tf(inf)*si(ini)*fact
      sumc19=sumc19+hypc1*tf(inf)*si(ini)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(inf)*si(ini)*fact
      sum20=sum20+hyp1*tf(inf)*si(ini)*fact
      sumc20=sumc20+hypc1*tf(inf)*si(ini)*fact
      enddo
      endif

      endif

      hyp2=0
      if(issym.eq.iqsym)then
      call gfact(orbs,orbq,iis,iiq,orbs,orbq,kaps,kapq,hyp2,hypc2)
      endif

      if(orbs.eq.orbp)then

      ini=idpa(iir,iia)
      ipssym=mtbl(ipsym,issym)

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      kk=igot(iloop)
      if(kk.eq.0)then
      nra=idra(iir,iia,0)
      nps=kdra(iip,iis,0)
      inf=nslocf(nps,nra,ipssym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum21=sum21+hyp2*ti(ini)*sf(inf)*fact
      sumc21=sumc21+hypc2*ti(ini)*sf(inf)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbp,ias,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nsa,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp2*ti(ini)*sf(inf)*fact
      sum22=sum22+hyp2*ti(ini)*sf(inf)*fact
      sumc22=sumc22+hypc2*ti(ini)*sf(inf)*fact
      enddo
      endif

      endif

      do it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      kapt=nka(iit)
      itsym=isyme(it)
      orbt=orbe(it)

      if(itsym.eq.ipsym)then
      if(orbt.eq.orbp)then
      if(orbs.eq.orbq)then
      if(iit.ne.iip)then

c     diagram 43
      hyp1=0
      if(itsym.eq.issym)then
      call gfact(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1,hypc1)
      endif

      in=idpa(iir,iia)
      inf=kdpa(iip,iit)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      nqs=jdra(iiq,iis,0)
      iqssym=mtbl(iqsym,issym)
      ini=nsloci(nqs,nra,iqssym,0)
      fact=dsqrt((2*orba+1)/(2*orbq+1))
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum37=sum37+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sumc37=sumc37+hypc1*tf(in)*si(ini)*sf(inf)*fact
      endif
      enddo
      endif

      call findk(orbs,orba,orbr,orbq,ias,iaa,iar,iaq,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nqr=jdra(iiq,iir,l1)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l1)
      fact=(-1)**(orba+orbq+l1)/(2*orbq+1)
      sum=sum+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sum38=sum38+hyp1*tf(in)*si(ini)*sf(inf)*fact
      sumc38=sumc38+hypc1*tf(in)*si(ini)*sf(inf)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

c     diagram 45
      if(iqsym.eq.issym)then
      if(orbq.eq.orbs)then
      if(orbp.eq.orbt)then
      if(iis.ne.iiq)then

      hyp1=0
      if(itsym.eq.issym)then
      call gfact(orbt,orbs,iit,iis,orbt,orbs,kapt,kaps,hyp1,hypc1)
      endif

      in=idpa(iir,iia)
      ini=jdpa(iiq,iis)
      call findk(orba,orbr,orbp,orbt,iaa,iar,iap,iat,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      if(l1.eq.0)then
      nra=idra(iir,iia,0)
      npt=kdra(iip,iit,0)
      iptsym=mtbl(ipsym,itsym)
      inf=nslocf(npt,nra,iptsym,0)
      fact=dsqrt((2*orba+1)/(2*orbp+1))
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum39=sum39+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sumc39=sumc39+hypc1*ti(in)*sf(inf)*si(ini)*fact
      endif
      enddo
      endif

      call findk(orbt,orba,orbr,orbp,iat,iaa,iar,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nta=idra(iit,iia,l1)
      npr=kdra(iip,iir,l1)
      iprsym=mtbl(ipsym,irsym)
      inf=nslocf(npr,nta,iprsym,l1)
      fact=(-1)**(orba+orbp+l1)/(2*orbp+1)
      sum=sum+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sum40=sum40+hyp1*ti(in)*sf(inf)*si(ini)*fact
      sumc40=sumc40+hypc1*ti(in)*sf(inf)*si(ini)*fact
      enddo
      endif

      endif
      endif
      endif
      endif

      enddo

      endif
      endif

      do 300 ib=1,nocc
      iib=kc(ib)
      iab=iqc(ib)
      orbb=orbc(ib)
      ibsym=isymc(ib)
      kapb=nka(iib)

c     diagrams  b1-b4 (sum19-sum22)
      hyp1=0
      if(ipsym.eq.ibsym)then
      call gfact(orbp,orbb,iip,iib,orbp,orbb,kapp,kapb,hyp1,hypc1)
      endif

      irasym=mtbl(irsym,iasym)
      iqssym=mtbl(iqsym,issym)

      if(irasym.eq.iqssym)then
      if(orbq.eq.orbb)then

      call findk(orba,orbr,orbs,orbb,iaa,iar,ias,iab,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nra,irasym,l1)

      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*l1+1))
      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum23=sum23-hyp1*tf(inf)*si(ini)*fact
      sumc23=sumc23-hypc1*tf(inf)*si(ini)*fact

      call findk(orba,orbs,orbr,orbq,iaa,ias,iar,iaq,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      iqrsym=mtbl(iqsym,irsym)
      ini=nsloci(nqr,nsa,iqrsym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)

      sum=sum-hyp1*tf(inf)*si(ini)*fact
      sum24=sum24-hyp1*tf(inf)*si(ini)*fact
      sumc24=sumc24-hypc1*tf(inf)*si(ini)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif

      hyp2=0
      if(iqsym.eq.ibsym)then
      call gfact(orbb,orbq,iib,iiq,orbb,orbq,kapb,kapq,hyp2,hypc2)
      endif

      ipssym=mtbl(ipsym,issym)
      if(irasym.eq.ipssym)then
      if(orbp.eq.orbb)then

      call findk(orba,orbr,orbs,orbp,iaa,iar,ias,iap,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*l1+1))
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum25=sum25-hyp2*ti(ini)*sf(inf)*fact
      sumc25=sumc25-hypc2*ti(ini)*sf(inf)*fact

      call findk(orba,orbs,orbr,orbb,iaa,ias,iar,iab,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      irbsym=mtbl(irsym,ibsym)
      ini=ntloc(nrb,nsa,irbsym,l2)
      call sixj(al1,orbs,orbp,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp2*ti(ini)*sf(inf)*fact
      sum26=sum26-hyp2*ti(ini)*sf(inf)*fact
      sumc26=sumc26-hypc2*ti(ini)*sf(inf)*fact

      enddo
      endif
      enddo
      endif

      endif
      endif


      do 600 ig=1,nocc
      igsym=isymc(ig)
      orbg=orbc(ig)
      iig=kc(ig)
      kapg=nka(iig)
      iag=iqc(ig)

      if(iasym.ne.igsym)go to 610
      hyp=0
c     c1-c2
      call gfact(orbg,orba,iig,iia,orbg,orba,kapg,kapa,hyp,hypc)

      isbsym=mtbl(issym,ibsym)
      iqssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irasym.ne.isbsym)go to 620
      if(irgsym.ne.iqssym)go to 620
      if(ipsym.ne.ibsym)go to 620
      if(orbp.ne.orbb)go to 620
      if(orbb.ne.orbq)go to 620
      in=idpa(iip,iib)

      call findk(orbg,orbr,orbq,orbs,iag,iar,iaq,ias,jgot,jmax)
      if(jmax.eq.0) go to 620
      do 630 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nqs=jdra(iiq,iis,l1)
      ini=nsloci(nqs,nrg,iqssym,l1)

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do 640 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nra,nsb,isbsym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*si(ini)*tf(inf)*ti(in)*fact
      sum27=sum27-hyp*si(ini)*tf(inf)*ti(in)*fact
      sumc27=sumc27-hypc*si(ini)*tf(inf)*ti(in)*fact
 640  continue
      endif

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,igot,imax)
      if(imax.ne.0)then
      do 650 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      inf=ntloc(nrb,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*tf(inf)*si(ini)*ti(in)*fact
      sum28=sum28-hyp*tf(inf)*si(ini)*ti(in)*fact
      sumc28=sumc28-hypc*tf(inf)*si(ini)*ti(in)*fact
  650 continue
      endif

  630 continue
  620 continue

c     d1-d2
      isbsym=mtbl(issym,ibsym)
      ipssym=mtbl(iqsym,issym)
      irasym=mtbl(irsym,iasym)
      irgsym=mtbl(irsym,igsym)

      if(irgsym.ne.isbsym)go to 660
      if(irasym.ne.ipssym)go to 660
      if(orbq.ne.orbb)go to 660
      if(iqsym.ne.ibsym)go to 660
      if(orbb.ne.orbp)go to 660
      in=idpa(iiq,iib)

      call findk(orbg,orbr,orbb,orbs,iag,iar,iab,ias,igot,imax)
      if(imax.eq.0)go to 660
      do 670 iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nrg=idra(iir,iig,l1)
      nsb=idra(iis,iib,l1)
      isbsym=mtbl(issym,ibsym)
      ini=ntloc(nrg,nsb,isbsym,l1)

      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,jgot,jmax)
      if(jmax.ne.0)then
      do 680 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,ipssym,l2)
      call sixj(rank,al1,al2,orbr,orba,orbg,s6j1)
      call sixj(rank,al1,al2,orbs,orbp,orbq,s6j2)
      fact=(-1)**(orbp+orbq+rank)*s6j1*s6j2
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum29=sum29-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sumc29=sumc29-hypc*ti(ini)*sf(inf)*tf(in)*fact
  680 continue
      endif

      call findk(orba,orbs,orbp,orbr,iaa,ias,iap,iar,jgot,jmax)
      if(jmax.ne.0)then
      do 690 jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      npr=kdra(iip,iir,l2)
      isasym=mtbl(issym,iasym)
      inf=nslocf(npr,nsa,isasym,l2)
      call ninej(al2,orbr,orbp,orba,orbg,rank,orbs,al1,orbq,s9j)
      fact=(-1)**(orba+orbs+orbr+orbg+l1+l2)*s9j
      sum=sum-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum30=sum30-hyp*ti(ini)*sf(inf)*tf(in)*fact
      sumc30=sumc30-hypc*ti(ini)*sf(inf)*tf(in)*fact
  690 continue
      endif

  670 continue
  660 continue
  610 continue


      if(ibsym.ne.igsym)go to 605
      if(ipsym.ne.ibsym)go to 605
c     if(orbp.ne.orbg)go to 605
      if(orbp.ne.orbb)go to 605
      in=idpa(iip,iib)
c     c7-c8
      hyp=0
      call gfact(orbb,orbg,iib,iig,orbb,orbg,kapb,kapg,hyp,hypc)
      if(orbq.ne.orbg)go to 700
      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,igot,imax)
      if (imax.eq.0)go to 700
      do 710 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqs=jdra(iiq,iis,l2)

      irasym=mtbl(irsym,iasym)
      isgsym=mtbl(issym,igsym)

      ini=nsloci(nqs,nra,irasym,l2)
      nsg=idra(iis,iig,l2)


      if(irasym.eq.isgsym)then

      inf=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbq+1))
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum31=sum31+hyp*tf(inf)*si(ini)*ti(in)*fact
      sumc31=sumc31+hypc*tf(inf)*si(ini)*ti(in)*fact
      endif

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 710
      do 720 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)

      isasym=mtbl(issym,iasym)
      irgsym=mtbl(irsym,igsym)

      inf=ntloc(nsa,nrg,isasym,l1)

      if(irgsym.eq.isasym)then

      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum32=sum32+hyp*tf(inf)*si(ini)*ti(in)*fact
      sumc32=sumc32+hypc*tf(inf)*si(ini)*ti(in)*fact
      endif

  720 continue
  710 continue
  700 continue
  605 continue
c     d7-d8

      hyp=0
      if(ibsym.eq.igsym)then
      call gfact(orbg,orbb,iig,iib,orbg,orbb,kapg,kapb,hyp,hypc)
      endif

c     if(orbq.ne.orbg)go to 600
      if(orbp.ne.orbg)go to 600
      if(orbb.ne.orbq)go to 600
      if(iqsym.ne.ibsym)go to 600
      in=idpa(iiq,iib)
      call findk(orba,orbr,orbp,orbs,iaa,iar,iap,ias,igot,imax)
      if (imax.eq.0)go to 600
      do 740 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nps=kdra(iip,iis,l2)
      irasym=mtbl(irsym,iasym)
      inf=nslocf(nps,nra,irasym,l2)
      nsg=idra(iis,iig,l2)
      ini=ntloc(nra,nsg,irasym,l2)
      fact=(-1)*(orba+orbr+orbg+orbs)/((2*l2+1)*(2*orbp+1))
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum33=sum33+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sumc33=sumc33+hypc*ti(ini)*sf(inf)*tf(in)*fact

      call findk(orba,orbs,orbr,orbg,iaa,ias,iar,iag,jgot,jmax)
      if(jmax.eq.0)go to 740
      do 750 jloop=1,jmax
      l1=jgot(jloop)
      al1=float(l1)
      nsa=idra(iis,iia,l1)
      nrg=idra(iir,iig,l1)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nsa,nrg,isasym,l1)
      call sixj(al1,orbr,orbq,al2,orbs,orba,s6j)
      fact=(-1)**(orba+orbr+orbg+orbs)*s6j/(2*orbq+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum34=sum34+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sumc34=sumc34+hypc*ti(ini)*sf(inf)*tf(in)*fact
  750 continue
  740 continue
  730 continue
  600 continue

      do 800 it=1,nexcit
      iit=ke(it)
      iat=iqe(it)
      orbt=orbe(it)
      kapt=nka(iit)
      itsym=isyme(it)

      hyp=0
      if(itsym.eq.issym)then
      call gfact(orbs,orbt,iis,iit,orbs,orbt,kaps,kapt,hyp,hypc)
      endif

      if(orbp.ne.orbb)go to 815
      if(ipsym.ne.ibsym)go to 815
      in=idpa(iip,iib)
c     c3
      call findk(orbq,orbt,orbr,orba,iaq,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 820
      do 830 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      nqt=jdra(iiq,iit,l2)
      irasym=mtbl(irsym,iasym)
      ini=nsloci(nqt,nra,irasym,l2)
      nsb=idra(iis,iib,l2)
      inf=ntloc(nsb,nra,irasym,l2)
      call sixj(rank,orbs,orbt,al2,orbq,orbp,s6j)
      fact=(-1)**(orba+orbr+l2+orbs+orbt+rank)*s6j/(2*al2+1)
      sum=sum+hyp*tf(inf)*si(ini)*ti(in)*fact
      sum35=sum35+hyp*tf(inf)*si(ini)*ti(in)*fact
      sumc35=sumc35+hypc*tf(inf)*si(ini)*ti(in)*fact

  830 continue
  820 continue
  815 continue

c     d3
      if(iqsym.ne.ibsym)go to 900
      if(orbb.ne.orbq)go to 900
      in=idpa(iiq,iib)

      call findk(orbb,orbt,orbr,orba,iab,iat,iar,iaa,igot,imax)
      if(imax.eq.0)go to 900
      do 930 iloop=1,imax
      l2=igot(iloop)
      al2=float(l2)
      nra=idra(iir,iia,l2)
      ntb=idra(iit,iib,l2)
      irasym=mtbl(irsym,iasym)
      ini=ntloc(ntb,nra,irasym,l2)
      nps=kdra(iip,iis,l2)
      inf=nslocf(nps,nra,irasym,l2)
      call sixj(rank,orbt,orbs,al2,orbp,orbq,s6j)
      fact=(-1)**(orbs+orbt+rank+orba+orbr)*s6j/(2*al2+1)
      sum=sum+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sum36=sum36+hyp*ti(ini)*sf(inf)*tf(in)*fact
      sumc36=sumc36+hypc*ti(ini)*sf(inf)*tf(in)*fact
  930 continue
  900 continue

  810 continue

      if(orbq.eq.orbb)then
      if(orbp.eq.orbt)then
      if(ipsym.eq.itsym)then
      if(iip.ne.iit)then

      in=kdpa(iip,iit)
      hyp=0
      if(itsym.eq.ibsym)then
      call gfact(orbt,orbb,iit,iib,orbt,orbb,kapt,kapb,hyp,hypc)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then

      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      inf=ntloc(nra,nsb,irasym,l1)

      call findk(orba,orbr,orbq,orbs,iaa,iar,iaq,ias,kgot,kmax)
      if(kmax.ne.0)then
      do kloop=1,kmax
      k1=kgot(kloop)
      if(k1.eq.l1)then
      nqs=jdra(iiq,iis,l1)
      nra=idra(iir,iia,l1)
      ini=nsloci(nqs,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbq+1)*(2*al1+1))

      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum41=sum41-hyp*si(ini)*tf(inf)*sf(in)*fact
      sumc41=sumc41-hypc*si(ini)*tf(inf)*sf(in)*fact
      endif
      enddo
      endif

      call findk(orba,orbs,orbq,orbr,iaa,ias,iaq,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nqr=jdra(iiq,iir,l2)
      isasym=mtbl(issym,iasym)
      ini=nsloci(nqr,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbq+1)

      sum=sum-hyp*si(ini)*tf(inf)*sf(in)*fact
      sum42=sum42-hyp*si(ini)*tf(inf)*sf(in)*fact
      sumc42=sumc42-hypc*si(ini)*tf(inf)*sf(in)*fact

      enddo
      endif
      enddo
      endif

      endif

      endif
      endif
      endif
      endif

      if(orbp.eq.orbb)then
      if(orbq.eq.orbt)then
      if(iqsym.eq.itsym)then
      if(iiq.ne.iit)then
      in=jdpa(iiq,iit)

      hyp=0
      if(itsym.eq.ibsym)then
      call gfact(orbb,orbt,iib,iit,orbb,orbt,kapb,kapt,hyp,hypc)
      endif

      irasym=mtbl(irsym,iasym)
      isbsym=mtbl(issym,ibsym)
      if(irasym.eq.isbsym)then
      call findk(orba,orbr,orbb,orbs,iaa,iar,iab,ias,igot,imax)
      if(imax.ne.0)then
      do iloop=1,imax
      l1=igot(iloop)
      al1=float(l1)
      nra=idra(iir,iia,l1)
      nsb=idra(iis,iib,l1)
      ini=ntloc(nra,nsb,irasym,l1)
      nps=kdra(iip,iis,l1)
      inf=nslocf(nps,nra,irasym,l1)
      fact=(-1)**(orba+orbr+orbs+orbb)/((2*orbp+1)*(2*al1+1))
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum43=sum43-hyp*sf(inf)*ti(ini)*si(in)*fact
      sumc43=sumc43-hypc*sf(inf)*ti(ini)*si(in)*fact

      call findk(orba,orbs,orbb,orbr,iaa,ias,iab,iar,jgot,jmax)
      if(jmax.ne.0)then
      do jloop=1,jmax
      l2=jgot(jloop)
      al2=float(l2)
      nsa=idra(iis,iia,l2)
      nrb=idra(iir,iib,l2)
      isasym=mtbl(issym,iasym)
      ini=ntloc(nrb,nsa,isasym,l2)
      call sixj(al1,orbs,orbq,al2,orbr,orba,s6j)
      fact=s6j*(-1)**(orba+orbr+orbs+orbb)/(2*orbp+1)
      sum=sum-hyp*sf(inf)*ti(ini)*si(in)*fact
      sum44=sum44-hyp*sf(inf)*ti(ini)*si(in)*fact
      sumc44=sumc44-hypc*sf(inf)*ti(ini)*si(in)*fact
      enddo
      endif
      enddo
      endif

      endif
      endif
      endif
      endif
      endif

  800 continue
  300 continue
  200 continue

   1  continue
 100  continue

      sum1=sum1*con/jorbv
      sum2=sum2*con/jorbv
      sum3=sum3*con/jorbv
      sum4=sum4*con/jorbv
      sum5=sum5*con/jorbv
      sum6=sum6*con/jorbv
      sum7=sum7*con/jorbv
      sum8=sum8*con/jorbv
      sum9=sum9*con/jorbv
      sum10=sum10*con/jorbv
      sum11=sum11*con/jorbv
      sum12=sum12*con/jorbv
      sum13=sum13*con/jorbv
      sum14=sum14*con/jorbv
      sum15=sum15*con/jorbv
      sum16=sum16*con/jorbv
      sum17=sum17*con/jorbv
      sum18=sum18*con/jorbv
      sum19=sum19*con/jorbv
      sum20=sum20*con/jorbv
      sum21=sum21*con/jorbv
      sum22=sum22*con/jorbv
      sum23=sum23*con/jorbv
      sum24=sum24*con/jorbv
      sum25=sum25*con/jorbv
      sum26=sum26*con/jorbv
      sum27=sum27*con/jorbv
      sum28=sum28*con/jorbv
      sum29=sum29*con/jorbv
      sum30=sum30*con/jorbv
      sum31=sum31*con/jorbv
      sum32=sum32*con/jorbv
      sum33=sum33*con/jorbv
      sum34=sum34*con/jorbv
      sum35=sum35*con/jorbv
      sum36=sum36*con/jorbv
      sum37=sum37*con/jorbv
      sum38=sum38*con/jorbv
      sum39=sum39*con/jorbv
      sum40=sum40*con/jorbv
      sum41=sum41*con/jorbv
      sum42=sum42*con/jorbv
      sum43=sum43*con/jorbv
      sum44=sum44*con/jorbv

      sumc1=sumc1*d1/jorbv
      sumc2=sumc2*d1/jorbv
      sumc3=sumc3*d1/jorbv
      sumc4=sumc4*d1/jorbv
      sumc5=sumc5*d1/jorbv
      sumc6=sumc6*d1/jorbv
      sumc7=sumc7*d1/jorbv
      sumc8=sumc8*d1/jorbv
      sumc9=sumc9*d1/jorbv
      sumc10=sumc10*d1/jorbv
      sumc11=sumc11*d1/jorbv
      sumc12=sumc12*d1/jorbv
      sumc13=sumc13*d1/jorbv
      sumc14=sumc14*d1/jorbv
      sumc15=sumc15*d1/jorbv
      sumc16=sumc16*d1/jorbv
      sumc17=sumc17*d1/jorbv
      sumc18=sumc18*d1/jorbv
      sumc19=sumc19*d1/jorbv
      sumc20=sumc20*d1/jorbv
      sumc21=sumc21*d1/jorbv
      sumc22=sumc22*d1/jorbv
      sumc23=sumc23*d1/jorbv
      sumc24=sumc24*d1/jorbv
      sumc25=sumc25*d1/jorbv
      sumc26=sumc26*d1/jorbv
      sumc27=sumc27*d1/jorbv
      sumc28=sumc28*d1/jorbv
      sumc29=sumc29*d1/jorbv
      sumc30=sumc30*d1/jorbv
      sumc31=sumc31*d1/jorbv
      sumc32=sumc32*d1/jorbv
      sumc33=sumc33*d1/jorbv
      sumc34=sumc34*d1/jorbv
      sumc35=sumc35*d1/jorbv
      sumc36=sumc36*d1/jorbv
      sumc37=sumc37*d1/jorbv
      sumc38=sumc38*d1/jorbv
      sumc39=sumc39*d1/jorbv
      sumc40=sumc40*d1/jorbv
      sumc41=sumc41*d1/jorbv
      sumc42=sumc42*d1/jorbv
      sumc43=sumc43*d1/jorbv
      sumc44=sumc44*d1/jorbv
      print*,'Dirac-Fock value ',(dirac-diracc)*d1/jorbv

      return
      end

c****************************************************************************


c**********************************************************************
c       This subroutine calculates the reduced matrix element of E1   *
c       operator in velocity gauge                                    *
c**********************************************************************

      subroutine e1vel(i1,i2,jorb1,jorb2,kap1,kap2,e1m)

      implicit double precision (a-h,o-z)
      double precision j1,j2,m1,m2,jorb1,jorb2
      common/cons/zero,half,tenth,one,two,three,ten
      common/inte1vel/rint1(MNBAS,MNBAS),rint2(MNBAS,MNBAS)
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/orbital_energy/eorb(MNBAS)

      alpha=137.0359895
      e1m=0.0

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
      e1m=e1*alpha/(eorb(i2)-eorb(i1))

      return
      end

c**************************************************************************

c***************************************************************************
c       This subroutine calculates the radial part of the reduced matrix   *
c       element of E1 transition in velocity gauge                         *
c***************************************************************************

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

c--------------anty symmetric part----------------------------------------

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

c--------------symmetric part--------------------------------------------

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

c***************************************************************************

      subroutine e2vel(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp)
      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
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

c************************************************************************

c*************************************************************************
c                                                                        *
c     This subroutine calculates the radial part of the reduced matrix   *
c     element for E2 transition operator                                 *
c                                                                        *
c*************************************************************************


      subroutine Calcquad(nbasis)


      implicit double precision(a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/orbital_energy/eorb(MNBAS)
      common/kappa/nka(MNBAS)       
      common/intquad/rint(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

       alpha=137.0359895

       mtp=n
       
       do i1=1,nbasis
       do i2=1,nbasis

       kap1=nka(i1)
       kap2=nka(i2)

       rint(i1,i2)=zero

c       if(iparity(i1).eq.iparity(i2))then

       ta(1)=zero

       do 1 l=2,mtp
      ta(l)=r(l)**2*rp(l)*((pf(l,i1)*pf(l,i2)+qf(l,i2)*qf(l,i1))
     : + r(l)*(eorb(i1)-eorb(i2))/(alpha*7)*(((kap1-kap2)/3*
     :(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1)))+
     :(pf(l,i1)*qf(l,i2)-qf(l,i1)*pf(l,i2))))
    1  continue
       call quad(result)

       rint(i1,i2)=result

c       endif

       enddo
       enddo

      return
      end

c*************************************************************************


      subroutine Calcquadvel(nbasis)

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
      rint(i2,i1) = result

      enddo
      enddo

      return
      end

c*************************************************************


c*************************************************************
c     This subroutine calculates the reduced matrix element  *
c     for E2 transition operator                             *
c                                                            *
c*************************************************************



      subroutine quadpul(j1,j2,i1,i2,m1,m2,kap1,kap2,hyp)


      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
      common/intquad/rint(MNBAS,MNBAS)

      kapp1=kap1
      kapp2=kap2

      factor=(-1)**(j1+0.5)
      w3j1=dr(j1,two,j2,half,zero,-half)
      e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
      result=rint(i1,i2)
      hyp=result*e1

      if(i1.eq.25.or.i2.eq.4)then
c      print*, 'H', j1,j2,i1,i2,m1,m2,kap1,kap2,result,e1,hyp
      endif
               
      return
      end

c*************************************************************************

c**************************************************************************
c                                                                         *
c    This subroutine calculates the reduced matrix element of the E2      *
c    hyperfine operator                                                   *      
c                                                                         *
c**************************************************************************


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

c**************************************************************************

c**************************************************************************
c                                                                         *
c      This subroutine calculates the radial integral of the reduced      *
c      matrix element of the E2 hyperfine operator                        *
c                                                                         *
c************************************************************************** 

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

c      if(iparity(i1).eq.iparity(i2))then

      ta(1)=zero
      do 1 l=2,mtp
      ta(l)=(pf(l,i1)*pf(l,i2)+qf(l,i2)*qf(l,i1))*rp(l)/r(l)**3
    1 continue

      call quad(result)

c      endif

      rint(i1,i2) = result
      rint(i2,i1) = rint(i1,i2)

      enddo
      enddo

      return
      end

c***************************************************************************

c***************************************************************************
c                                                                          *
c       This subroutine calculates the radial part of the reduced          *
c       matrix element of the M1 transition operator                       *
c                                                                          *
c***************************************************************************

      subroutine Calcm1prop(nbasis)

      implicit double precision(a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/radialm1/hyp(MNBAS,MNBAS)
      common/syminfo/orbj(MNBAS),iiq(MNBAS),iparity(MNBAS)

       mtp=n

       do i1=1,nbasis
       do i2=1,i1
   
        hyp(i1,i2)=0

c       if(iparity(i1).eq.iparity(i2))then

       ta(1)=zero

       do 1 l=2,mtp
       ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)*r(l)
    1  continue
       call quad(result)

       hyp(i1,i2)=result
       hyp(i2,i1)=result

c       endif

       enddo
       enddo

       return
       end

c****************************************************************************



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
      implicit double precision(a-h,o-z)
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
      common/intg/rint(MNBAS,MNBAS)
      common/intgc/rintc(MNBAS,MNBAS)

       mtp=n

       do i1=1,nbasis
       do i2=1,i1

       ta(1)=zero
       do 1 l=2,mtp
       ta(l)=(pf(l,i1)*qf(l,i2)+pf(l,i2)*qf(l,i1))*rp(l)*r(l)
    1  continue
       call quad(result)
       rint(i1,i2)=result
       rint(i2,i1)=result

       ta(1)=zero
       do 2 l=2,mtp
       ta(l)=(pf(l,i1)*pf(l,i2)+qf(l,i2)*qf(l,i1))*rp(l)
    2  continue
       call quad(resultc)
       rintc(i1,i2)=resultc
       rintc(i2,i1)=resultc

       enddo
       enddo

      return
      end

c*********************************************************************

c*********************************************************************
c                                                                    *
c    This subroutine calculates the reduced matrix element for the   *
c     M1 transition operator                                         * 
c                                                                    *
c********************************************************************


      subroutine m1prop(j1,j2,i1,i2,m1,m2,kap1,kap2,val)

      implicit double precision(a-h,o-z)
      real*8 j1,j2,m1,m2
      common/cons/zero,half,tenth,one,two,three,ten
     :/grid/r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :/tatb/ta(MN),tb(MN),mtp
     :/wave/pz(MNBAS),pf(MN,MNBAS),qf(MN,MNBAS),
     :mf(MNBAS)
       common/radialm1/hyp(MNBAS,MNBAS)

       data alpha/137.0359895/

       factor=(-1)**(j1+0.5)
       w3j1=dr(j1,one,j2,half,zero,-half)
       e1=dsqrt(2*j1+one)*dsqrt(2*j2+one)*factor*w3j1
       result=hyp(i1,i2)
       val=-result*e1*(kap1+kap2)*alpha

       return
       end
       

c*********************************************************
      subroutine m2prop(i1,i2,j1,j2,kap1,kap2,hyp)
c*********************************************************
c  this subrouti calculates the edm t-pt reduced
c  matrix elements between orbitals with j=1/2
c*********************************************************
      implicit double precision(a-h,o-z)
      real*8 j1,j2
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
c  this subrouti calculates the edm t-pt reduced
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

c**********************************************************************
c                                                                     *
c      This subroutine calculates the E2 hyperfine diagrams with the  *
c      total results                                                  *
c                                                                     *
c**********************************************************************


      subroutine calquadB

      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valueB/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44,sum45,sum46

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44
      write(STDOUT,27)sum45
      write(STDOUT,28)sum46

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26
     :    +sum27+sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40+sum41+sum42
     :    +sum43+sum44+sum45+sum46

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)
   27 format(/,1x,'contribution from diagram 45',f20.11)
   28 format(/,1x,'contribution from diagram 46',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

c***********************************************************************

      subroutine calgfact
      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valueg/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44
      common/valuegc/sumc1,sumc2,sumc3,sumc4,sumc5,sumc6,sumc7,
     :sumc8,sumc9,sumc10,sumc11,sumc12,sumc13,sumc14,sumc15,sumc16,
     :sumc17,sumc18,sumc19,sumc20,sumc21,sumc22,sumc23,sumc24,sumc25,
     :sumc26,sumc27,sumc28,sumc29,sumc30,sumc31,sumc32,sumc33,sumc34,
     :sumc35,sumc36,sumc37,sumc38,sumc39,sumc40,sumc41,sumc42,sumc43,
     :sumc44

      data VL/137.03598950D+00/

c      cnorm=cnorm1*cnorm2
      cnorm=cnorm2*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26+sum27
     :    +sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40+sum41+sum42
     :    +sum43+sum44

      sumc=sumc1+sumc2+sumc3+sumc4+sumc5+sumc6+sumc7
     :    +sumc8+sumc9+sumc10+sumc11+sumc12+sumc13
     :    +sumc14+sumc15+sumc16+sumc17+sumc18
     :    +sumc19+sumc20+sumc21+sumc22+sumc23
     :    +sumc24+sumc25+sumc26+sumc27
     :    +sumc28+sumc29+sumc30+sumc31+sumc32
     :    +sumc33+sumc34+sumc35+sumc36+sumc37
     :    +sumc38+sumc39+sumc40+sumc41+sumc42
     :    +sumc43+sumc44

      corr=sum*fact
      corrQED=sumc*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm
      write(*,9997)corrQED
      write(*,9998)sumc
      write(*,9999)sumc/cnorm
      write(*,2000)cnorm
      write(*,2001)cnorm1
      write(*,2002)(sum+sumc)/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)
 9997 format(/,1x,'Norm QED Correction',f20.11)
 9998 Format(/,1x,'Without norm QED Correction',f20.11)
 9999 format(/,1x,'Approx QED Correction',f20.11)
 2000 format(/,1x,'Normalisation Factor',f20.11)
 2001 format(/,1x,'Normalisation Factor of Init',f20.11)
 2002 format(/,1x,'The Net g_j factor:         ',f20.11)

      return
      end

c**********************************************************************
c                                                                     *
c      This subroutine calculates M1 transition property diagrams     *
c      with the total results                                         *
c                                                                     *
c**********************************************************************

      subroutine calm1

      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valuem1/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26
     :    +sum27+sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40+sum41+sum42
     :    +sum43+sum44

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

c**********************************************************************

      subroutine calm2
      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valuem2/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26+sum27
     :    +sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40+sum41+sum42
     :    +sum43+sum44

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)
      return
      end

c******************************************************************
c     This subroutine calculate the finel result with and without *
c     normalization correction to the electric dipole transition  *
c     matrix element in velocity gauge                            *
c******************************************************************

      subroutine caldipvel

      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valuedvel/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26
     :    +sum27+sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40
     :    +sum43+sum44

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

c*********************************************************************

c******************************************************************
c     This subroutine calculate the finel result with and without *
c     normalization correction to the Nuclear spin independent    *
c     matrix element                                              *
c******************************************************************


      subroutine calnsi

      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valuensi/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26+sum27
     :    +sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40
     :    +sum43+sum44

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

**************************************************************************


c******************************************************************
c     This subroutine calculate the finel result with and without *
c     normalization correction to the nuclear spin dependent      *
c     matrix element                                              *
c******************************************************************


      subroutine calnsd

      implicit real*8 (a-h,o-z)

      common/norm/cnorm1,cnorm2
      common/valuensd/sum1,sum2,sum3,sum4,sum5,sum6,sum7,
     :sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15,sum16,
     :sum17,sum18,sum19,sum20,sum21,sum22,sum23,sum24,sum25,
     :sum26,sum27,sum28,sum29,sum30,sum31,sum32,sum33,sum34,
     :sum35,sum36,sum37,sum38,sum39,sum40,sum41,sum42,sum43,
     :sum44

      cnorm=cnorm1*cnorm2

      write(STDOUT,9980)sum1
      write(STDOUT,9981)sum2
      write(STDOUT,9982)sum3
      write(STDOUT,9972)sum4
      write(STDOUT,9983)sum5
      write(STDOUT,9973)sum6
      write(STDOUT,9984)sum7
      write(STDOUT,9985)sum8
      write(STDOUT,9975)sum9
      write(STDOUT,9986)sum10
      write(STDOUT,9976)sum11
      write(STDOUT,9987)sum12
      write(STDOUT,9988)sum13
      write(STDOUT,9989)sum14
      write(STDOUT,9990)sum15
      write(STDOUT,9991)sum16
      write(STDOUT,9992)sum17
      write(STDOUT,9993)sum18
      write(STDOUT,1)sum19
      write(STDOUT,2)sum20
      write(STDOUT,3)sum21
      write(STDOUT,4)sum22
      write(STDOUT,5)sum23
      write(STDOUT,6)sum24
      write(STDOUT,7)sum25
      write(STDOUT,8)sum26
      write(STDOUT,9)sum27
      write(STDOUT,10)sum28
      write(STDOUT,11)sum29
      write(STDOUT,12)sum30
      write(STDOUT,13)sum31
      write(STDOUT,14)sum32
      write(STDOUT,15)sum33
      write(STDOUT,16)sum34
      write(STDOUT,17)sum35
      write(STDOUT,18)sum36
      write(STDOUT,19)sum37
      write(STDOUT,20)sum38
      write(STDOUT,21)sum39
      write(STDOUT,22)sum40
      write(STDOUT,23)sum41
      write(STDOUT,24)sum42
      write(STDOUT,25)sum43
      write(STDOUT,26)sum44

      fact=(1.0d0-cnorm)/cnorm

      sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7
     :    +sum8+sum9+sum10+sum11+sum12+sum13
     :    +sum14+sum15+sum16+sum17+sum18
     :    +sum19+sum20+sum21+sum22+sum23
     :    +sum24+sum25+sum26+sum27
     :    +sum28+sum29+sum30+sum31+sum32
     :    +sum33+sum34+sum35+sum36+sum37
     :    +sum38+sum39+sum40
     :    +sum43+sum44

      corr=sum*fact
      write(*,9994)corr
      write(*,9995)sum
      write(*,9996)sum/cnorm

 9980 format(/,1x,'contribution from diagram  1',f20.11)
 9981 format(/,1x,'contribution from diagram  2',f20.11)
 9982 format(/,1x,'contribution from diagram  3',f20.11)
 9972 format(/,1x,'contribution from diagram  4',f20.11)
 9983 format(/,1x,'contribution from diagram  5',f20.11)
 9973 format(/,1x,'contribution from diagram  6',f20.11)
 9984 format(/,1x,'contribution from diagram  7',f20.11)
 9985 format(/,1x,'contribution from diagram  8',f20.11)
 9975 format(/,1x,'contribution from diagram  9',f20.11)
 9986 format(/,1x,'contribution from diagram 10',f20.11)
 9976 format(/,1x,'contribution from diagram 11',f20.11)
 9987 format(/,1x,'contribution from diagram 12',f20.11)
 9988 format(/,1x,'contribution from diagram 13',f20.11)
 9989 format(/,1x,'contribution from diagram 14',f20.11)
 9990 format(/,1x,'contribution from diagram 15',f20.11)
 9991 format(/,1x,'contribution from diagram 16',f20.11)
 9992 format(/,1x,'contribution from diagram 17',f20.11)
 9993 format(/,1x,'contribution from diagram 18',f20.11)
    1 format(/,1x,'contribution from diagram 19',f20.11)
    2 format(/,1x,'contribution from diagram 20',f20.11)
    3 format(/,1x,'contribution from diagram 21',f20.11)
    4 format(/,1x,'contribution from diagram 22',f20.11)
    5 format(/,1x,'contribution from diagram 23',f20.11)
    6 format(/,1x,'contribution from diagram 24',f20.11)
    7 format(/,1x,'contribution from diagram 25',f20.11)
    8 format(/,1x,'contribution from diagram 26',f20.11)
    9 format(/,1x,'contribution from diagram 27',f20.11)
   10 format(/,1x,'contribution from diagram 28',f20.11)
   11 format(/,1x,'contribution from diagram 29',f20.11)
   12 format(/,1x,'contribution from diagram 30',f20.11)
   13 format(/,1x,'contribution from diagram 31',f20.11)
   14 format(/,1x,'contribution from diagram 32',f20.11)
   15 format(/,1x,'contribution from diagram 33',f20.11)
   16 format(/,1x,'contribution from diagram 34',f20.11)
   17 format(/,1x,'contribution from diagram 35',f20.11)
   18 format(/,1x,'contribution from diagram 36',f20.11)
   19 format(/,1x,'contribution from diagram 37',f20.11)
   20 format(/,1x,'contribution from diagram 38',f20.11)
   21 format(/,1x,'contribution from diagram 39',f20.11)
   22 format(/,1x,'contribution from diagram 40',f20.11)
   23 format(/,1x,'contribution from diagram 41',f20.11)
   24 format(/,1x,'contribution from diagram 42',f20.11)
   25 format(/,1x,'contribution from diagram 43',f20.11)
   26 format(/,1x,'contribution from diagram 44',f20.11)

 9994 format(/,1x,'normalization     correction',f20.11)
 9995 format(/,1x,'sum (excluding correction)  ',f20.11)
 9996 format(/,1x,'sum (excluding correction)/normalization',f20.11)

      return
      end

c*********************************************************************

c*********************************************************************

      FUNCTION X12J(J1,J2,J3,J4,L1,L2,L3,L4,K1,K2,K3,K4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     CALCULATION OF 12J SYMBOL OF 2ND KIND AS DEFINED BY
C       YUTSIS ET AL., EQ. (19.3), P. 63
C     CALLS J6J (MOLSCAT SHULTEN-GORDON ALGORITHM) FOR 6J'S
C     VERSION 2. CHECKS FOR REPEATED 6J INDICES (SG, 22 NOV 94)
C
      PARAMETER (MXDIM=100)
      DIMENSION X6J1(MXDIM),X6J2(MXDIM),X6J3(MXDIM),X6J4(MXDIM)
      DATA TENTH/0.1D0/
C
C     CHECK FOR K4=0.  THIS WILL ARISE FROM LOU MONCHICK'S GHM EQS.
C     ONE MIGHT WANT TO CHECK FOR OTHER ZEROS AND REARRANGE
      IF (K4.EQ.0) THEN
C       USE YUTSIS EQ. (19.10)
        IF (J4.NE.L4.OR.J3.NE.L2) THEN
           X12J=0.D0
           RETURN
        ENDIF
        X12J=S6XJ(DBLE(L1),DBLE(L2),DBLE(K2),
     : DBLE(K3),DBLE(K1),DBLE(J1))*S6XJ(DBLE(K2),
     : DBLE(L3),DBLE(L4),DBLE(J2),DBLE(K3),DBLE(K1))/
     1       SQRT(DBLE((L2+L2+1)*(L4+L4+1)))
        LP=J1+L1-J2-L3
        IF (LP-2*(LP/2).NE.0) X12J=-X12J
        RETURN
      ENDIF
C
C     GENERAL CASE:  GET MIN, MAX VALUES OF SUMMATION INDEX
      IXMIN=MAX(ABS(K1-K2),ABS(K3-K4),ABS(J1-J3),ABS(J2-J4))
      IXMAX=MIN(K1+K2,K3+K4,J1+J3,J2+J4)
      IF (IXMIN.GT.IXMAX) THEN
C       12J MUST VANISH IN THIS CASE
        X12J=0.D0
        RETURN
      ENDIF
C
C     CONVERT EVERYONE TO DOUBLE PRECISION FOR INPUT TO J6J
      XJ1=J1
      XJ2=J2
      XJ3=J3
      XJ4=J4
      XL1=L1
      XL2=L2
      XL3=L3
      XL4=L4
      XK1=K1
      XK2=K2
      XK3=K3
      XK4=K4
C
C     GET 4 SETS OF 6J'S  -- CHECK TO SEE IF SOME ARE SAME
      IVAL1=MXDIM
      CALL J6J(XK1,XK2,XL1,XJ3,XJ1,IVAL1,XMIN1,X6J1)
      IF (K3.EQ.K1.AND.K4.EQ.K2.AND.L2.EQ.L1) THEN
        IVAL2=IVAL1
        XMIN2=XMIN1
        DO 2001 II=1,IVAL2
 2001   X6J2(II)=X6J1(II)
      ELSE
        IVAL2=MXDIM
        CALL J6J(XK3,XK4,XL2,XJ3,XJ1,IVAL2,XMIN2,X6J2)
      ENDIF
      IF (L3.EQ.L1.AND.J4.EQ.J3.AND.J2.EQ.J1) THEN
        IVAL3=IVAL1
        XMIN3=XMIN1
        DO 2002 II=1,IVAL3
 2002   X6J3(II)=X6J1(II)
      ELSEIF (J4.EQ.J3.AND.J2.EQ.J1.AND.L3.EQ.L2.AND.
     1        K1.EQ.K3.AND.K2.EQ.K4) THEN
        IVAL3=IVAL2
        XMIN3=XMIN2
        DO 2003 II=1,IVAL3
 2003   X6J3(II)=X6J2(II)
      ELSE
        IVAL3=MXDIM
        CALL J6J(XK1,XK2,XL3,XJ4,XJ2,IVAL3,XMIN3,X6J3)
      ENDIF
      IF (L4.EQ.L2.AND.J4.EQ.J3.AND.J2.EQ.J1) THEN
        IVAL4=IVAL2
        XMIN4=XMIN2
        DO 2004 II=1,IVAL4
 2004   X6J4(II)=X6J2(II)
      ELSEIF (K3.EQ.K1.AND.K4.EQ.K2.AND.L4.EQ.L3) THEN
        IVAL4=IVAL3
        XMIN4=XMIN3
        DO 2005 II=1,IVAL4
 2005   X6J4(II)=X6J3(II)
      ELSEIF (K3.EQ.K1.AND.K4.EQ.K2.AND.L4.EQ.L1.AND.
     1        J4.EQ.J3.AND.J2.EQ.J1) THEN
        IVAL4=IVAL1
        XMIN4=XMIN1
        DO 2006 II=1,IVAL4
 2006   X6J4(II)=X6J1(II)
      ELSE
        IVAL4=MXDIM
        CALL J6J(XK3,XK4,XL4,XJ4,XJ2,IVAL4,XMIN4,X6J4)
      ENDIF
C
C     SUM OVER IX=IXMIN,IXMAX
      SUM=0.D0
      DO 1000 IX=IXMIN,IXMAX
      IND1=1+IX-INT(XMIN1+TENTH)
      IF (IND1.LT.1.OR.IND1.GT.IVAL1) GO TO 1000
      IND2=1+IX-INT(XMIN2+TENTH)
      IF (IND2.LT.1.OR.IND2.GT.IVAL2) GO TO 1000
      IND3=1+IX-INT(XMIN3+TENTH)
      IF (IND3.LT.1.OR.IND3.GT.IVAL3) GO TO 1000
      IND4=1+IX-INT(XMIN4+TENTH)
      IF (IND4.LT.1.OR.IND4.GT.IVAL4) GO TO 1000
      SUM=SUM+X6J1(IND1)*X6J2(IND2)*X6J3(IND3)*X6J4(IND4)*DBLE(IX+IX+1)
 1000 CONTINUE
C     GET PARITY(L1-L2-L3+L4)
      LP=L1-L2-L3+L4
      IF (LP-2*(LP/2).NE.0) SUM=-SUM
      X12J=SUM
      RETURN
      END

      FUNCTION S6XJ(AJ1,AJ2,AJ3,AJ4,AJ5,AJ6)
C
C     CALCULATES 6-J SYMBOL:   _(J1 J2 J3 )_
C                               (J4 J5 J6 )
C     INTERFACE TO J6J ROUTINE.
C     MODIFIED BY S. GREEN 20 AUG 93; PASS DIMENSION OF XJ6J FOR CHECKING
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(MXDIM=200)
      DIMENSION XJ6J(MXDIM)
      IVAL=MXDIM
      CALL J6J(         AJ2,AJ3,
     1         AJ4,AJ5,AJ6,
     3         IVAL,XJ1MIN,XJ6J)
      IND=1+INT(AJ1)-INT(XJ1MIN+0.1D0)
      S6XJ=0.D0
      IF(IND.GE.1 .AND. IND.LE.IVAL) S6XJ=XJ6J(IND)
      RETURN
      END

      SUBROUTINE J6J(J2,J3,L1,L2,L3,IVAL,J1MIN,D6J)
      IMPLICIT DOUBLE PRECISION (A-H,J-Z)
      DIMENSION D6J(2)
      DATA ZERO/0.D0/,TENTH/0.1D0/,HALF/0.5D0/,ONE/1.D0/,TWO/2.D0/,
     $ CONST/1.0D-12/
      E(J1S)=SQRT((J1S-MJ23S)*(J23P1S-J1S)*(J1S-ML23S)*(L23P1S-J1S))
      F(J1,JJP1)=(TWO*J1+ONE)*(JJP1*(FACT-JJP1-TWO*LLP1)+FACT2)
C
C  THIS ROUTINE CALCULATES THE 6-J COEFFICIENTS FOR ALL PERMISSIBLE
C  VALUES OF J1 FOR FIXED VALUES OF J2, J3, L1, L2, AND L3 USING THE
C  RECURSIVE ALGORITHM OF K. SCHULTEN AND R. G. GORDON, J. MATH. PHYS.
C  VOL. 16, P. 1961, (1975).
C  PROGRAMMED BY D. E. FITZ, 10/22/79
C  MODIFIED BY S. GREEN 20 AUG 93 TO TEST DIMENSION ON D6J
C
      MXDIM=IVAL
      JJP2=J2*(J2+ONE)
      JJP3=J3*(J3+ONE)
      LLP1=L1*(L1+ONE)
      LLP2=L2*(L2+ONE)
      LLP3=L3*(L3+ONE)
      MJ23S=(J2-J3)**2
      ML23S=(L2-L3)**2
      J23P1S=(J2+J3+ONE)**2
      L23P1S=(L2+L3+ONE)**2
      FACT2=(LLP2-LLP3)*(JJP2-JJP3)
      FACT=JJP2+JJP3+LLP2+LLP3
      J1MIN=MAX(ABS(J2-J3),ABS(L2-L3))
      J1MAX=MIN(J2+J3,L2+L3)
      IVAL=INT(J1MAX-J1MIN+ONE+TENTH)
      IF (IVAL.GT.MXDIM) THEN
        WRITE(6,*) 'J6J: ARRAY D6J TOO SMALL. NEEDS ',IVAL,' BUT ONLY ',
     1   MXDIM,' SUPPLIED'
        STOP
      ENDIF
C
C  TEST FOR OTHER TRIANGULAR INEQUALITES.
C
      IL1=INT(TWO*L1+TENTH)
      IL2=INT(TWO*L2+TENTH)
      IL3=INT(TWO*L3+TENTH)
      IJ2=INT(TWO*J2+TENTH)
      IJ3=INT(TWO*J3+TENTH)
      IF((IJ2.LE.IL1+IL3.AND.IJ2.GE.IABS(IL1-IL3)).AND.
     $ (IJ3.LE.IL1+IL2.AND.IJ3.GE.IABS(IL1-IL2))) GO TO 11
      DO 12 I=1,IVAL
   12 D6J(I)=ZERO
      RETURN
C
   11 INMID=(IVAL+3)/2
      SGNV=J2+J3+L2+L3
      SGN=ONE
      ISIGN=INT(SGNV+TENTH)
      IF(MOD(ISIGN,2).NE.0) SGN=-ONE
      D6J(1)=HALF
C
C  UPWARD RECURSION.
C
      IF(IVAL.EQ.1) GO TO 40
      JJP1=J1MIN*(J1MIN+ONE)
      F1=F(J1MIN,JJP1)
      J1=J1MIN+ONE
      J1S=J1*J1
      E2=E(J1S)
      IF(J1MIN.LT.TENTH) GO TO 15
      D6J(2)=-F1*D6J(1)/(E2*J1MIN)
      GO TO 16
   15 D6J(2)=-HALF*(LLP2+JJP2-LLP1)*D6J(1)/SQRT(JJP2*LLP2)
   16 SCALE=D6J(2)
      IF(IVAL.EQ.2) GO TO 40
      DO 21 IJ2=3,INMID
      JJP1=J1*(J1+ONE)
      F1=F(J1,JJP1)
      J1=J1+ONE
      E1=E2
      J1S=J1*J1
      E2=E(J1S)
   21 D6J(IJ2)=-(F1*D6J(IJ2-1)+J1*E1*D6J(IJ2-2))/(E2*(J1-ONE))
      SCALE=D6J(INMID)
      IEXC=5
      IF(ABS(SCALE).GT.CONST) GO TO 18
      INMID=INMID-1
      SCALE=D6J(INMID)
      IEXC=3
      GO TO 30
   18 IF(IVAL.EQ.3) GO TO 40
C
C  DOWNWARD RECURSION.
C
   30 D6J(IVAL)=HALF
      J1=J1MAX
      J1S=J1*J1
      JJP1=J1*(J1+ONE)
      F1=F(J1,JJP1)
      E1=E(J1S)
      D6J(IVAL-1)=-F1*D6J(IVAL)/(E1*(J1+ONE))
      IEND=IVAL-INMID
      IF(IVAL.LE.IEXC) GO TO 31
      DO 32 IJ2=2,IEND
      J1=J1-ONE
      E2=E1
      J1S=J1*J1
      JJP1=J1*(J1+ONE)
      E1=E(J1S)
      F1=F(J1,JJP1)
   32 D6J(IVAL-IJ2)=-(J1*E2*D6J(IVAL-IJ2+2)+F1*D6J(IVAL-IJ2+1))/
     $ (E1*(J1+ONE))
C
C  MATCH UPWARD AND DOWNWARD RECURSIVE RESULTS BY SCALING.
C
   31 SCALE=SCALE/D6J(INMID)
      DO 33 IJ2=INMID,IVAL
   33 D6J(IJ2)=SCALE*D6J(IJ2)
C
C  NORMALIZE RESULTS AND SET PHASE.
C
   40 SUM=ZERO
      DO 41 IJ2=1,IVAL
      J1=J1MIN+DBLE(IJ2-1)
   41 SUM=SUM+(TWO*J1+ONE)*D6J(IJ2)**2
      RNORM=ONE/SQRT(SUM*(TWO*L1+ONE))
      IF((SGN*D6J(IVAL)).LT.ZERO) RNORM=-RNORM
      DO 42 IJ2=1,IVAL
   42 D6J(IJ2)=D6J(IJ2)*RNORM
      RETURN
      END

      FUNCTION THREEJ (J1,J2,J3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     COMPUTATION OF SPECIAL WIGNER 3J COEFFICIENT WITH
C     VANISHING PROJECTIONS.  SEE EDMONDS, P. 50.
C
C     THIS VERSION EVALUATES BINOM AND PARITY IN-LINE
C       SHOULD IMPROVE EFFICIENCY, ESPECIALLY ON CRAY;
C       ALSO GIVES IMPROVEMENT ON AMDAHL  (SG: 20 DEC 92)
C
C     STATEMENT FUNCTION FOR DELTA ASSOCIATED W/ RACAH AND SIXJ SYMBOLS
C     DELTA(I,J,K)= SQRT(1.D0/ ( BINOM(I+J+K+1,I+J-K) *
C    1                 BINOM(K+K+1,I-J+K) * DBLE(K+J-I+1) )  )
C
      I1=J1+J2+J3
      IF (I1-2*(I1/2).NE.0)    GO TO 8
    1 I2=J1-J2+J3
      IF (I2) 8,2,2
    2 I3=J1+J2-J3
      IF (I3) 8,3,3
    3 I4=-J1+J2+J3
      IF (I4) 8,4,4
    4 I5=I1/2
      I6=I2/2
      SIGN=1.D0
      IF (I5-2*(I5/2).NE.0) SIGN=-SIGN
C   7 THREEJ=SIGN*DELTA(J1,J2,J3)*BINOM(I5,J1)*BINOM(J1,I6)
C     B1,B2 ARE BINOM ASSOCIATED W/ DELTA
      N=J1+J2+J3+1
      M=J1+J2-J3
      NM = N-M
      MNM = MIN(NM,M)
      IF(MNM.LE.0) THEN
        B1=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 101 I = 1,MNM
        F = F+1.D0
        C = (FN-F)*B
  101   B = C/F
        B1 = B
      ENDIF
      N=J3+J3+1
      M=J1-J2+J3
      NM = N-M
      MNM = MIN(NM,M)
      IF(MNM.LE.0) THEN
        B2=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 102 I = 1,MNM
        F = F+1.D0
        C = (FN-F)*B
  102   B = C/F
        B2 = B
      ENDIF
      DELTA=SQRT(1.D0/(B1*B2*(J3+J2-J1+1)))
C     B3=BINOM(I5,J1),  B4=BINOM(J1,I6)
      N=I5
      M=J1
      NM = N-M
      MNM = MIN(NM,M)
      IF(MNM.LE.0) THEN
        B3=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 103 I = 1,MNM
        F = F+1.D0
        C = (FN-F)*B
  103   B = C/F
        B3 = B
      ENDIF
      N=J1
      M=I6
      NM = N-M
      MNM = MIN(NM,M)
      IF(MNM.LE.0) THEN
        B4=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 104 I = 1,MNM
        F = F+1.D0
        C = (FN-F)*B
  104   B = C/F
        B4 = B
      ENDIF
      THREEJ=SIGN*DELTA*B3*B4
      RETURN
    8 THREEJ=0.D0
      RETURN
      END

