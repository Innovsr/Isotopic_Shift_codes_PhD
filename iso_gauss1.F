#define   MNS   11
#define MAXVAL  10
#define MNOCC   30
#define MNBAS  300
#define   MNB  700
#define    MN  750
#define MFACT  500

      program main
      implicit real*8(a-h,o-z)
      call readinp
      call setgrd
      call setqic
      call gaussread
      call cal
      stop
      end

!-----------------------------------------------------
       subroutine readinp
!-----------------------------------------------------
       
      implicit real*8 (a-h,o-z)
      common /input/nn,ll,nb(MNS),nsy
     :  /symm/jsy,nbas(MNBAS)
     :  /atom/iz 

!1s,2s,3s,2p-,3p-,2p,3p,3d-,3d,4s write the grasp input file in this way
      open(unit=29,file='input_iso',form='formatted',status='old')
      read(29,*)jsy,(nbas(i),i=1,jsy)
      read(29,*)nsy,(nb(i),i=1,nsy)
! principle quantum no of valence
      read(29,*)nn
! orbital ang q num. of valence
      read(29,*)ll

!atomic num. of the system
      read(29,*)iz
        return
        end

!--------------------------------------------------
       subroutine gaussread
!-----------------------------------------------------
      implicit real*8(a-h,o-z) 
       dimension eorb(MNBAS),pz(MNBAS)
c     ,pf(MN,MNBAS),qf(MN,MNBAS)
      common /wave/ppf(MN,MNBAS),qqf(MN,MNBAS),mf(MNBAS)
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /input/nn,ll,nb(MNS),nsy
     :  /symm/jsy,nbas(MNBAS)
c     :  /nba/nbas(MNS)

      open(unit=7,file='wfn.dat',form='unformatted',status='old')
      kk=0
      nbasis=0
      rewind(7)
      read(7)h,n
      read(7)(r(i),i = 1,n),(rp(i),i = 1,n),(rpor(i),i=1,n)
      do j=1,jsy
      nbasis=nbasis+nbas(j)
      enddo
c      print*,nbasis
      do i=1,nbasis
      read(7)orben
      read(7)(ppf(ii,i),ii=1,n),(qqf(ii,i),ii=1,n)
      enddo

      rewind(7)

 
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine cal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            implicit real*8(a-h,o-z)
!       dimension eorb(MNBAS),pz(MNBAS)
      common /wave/ppf(MN,MNBAS),qqf(MN,MNBAS),mf(MNBAS)
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /input/nn,ll,nb(MNS),nsy
     :  /symm/jsy,nbas(MNBAS)



       nbass=0
      do j=1,jsy
      print*,nbas(j)
      enddo
      

c      if(nn.eq.2) then
c      do i=1+ll,2
c      nbass=nbass+nbas(i)
c      ia=nbass+1
c      ic=1
c      call gen_iso(ic,ia,breit,sms)
c      print*,breit
c       print*,sms
c      enddo
c      endif


      if(nn.ge.2) then
      do i=1+ll,1+(nn-1)*2
      

      breit3=0.0
      breit2=0.0
       breit4=0.0

      sms3=0.0
      sms2=0.0
      sms4=0.0


      if(i.eq.1.and.nsy.ge.2) then
        print*,'S',nn
        ia=nn
        ka=-1
        kc=1
        icc=nbas(1)+1
        iccc=nbas(1)+nb(2)
      do ic=icc,iccc
        call gen_iso(ic,ia,ka,kc,breit,sms)
        breit2=breit+breit2
        sms2=sms+sms2
      enddo
        ia=nn
      do j=1,2
        nbass=nbass+nbas(j)
      enddo
        icc=nbas(1)+nbas(2)+1
        iccc=nbas(1)+nbas(2)+nb(3)
        kc=-2
      do ic=icc,iccc
        call gen_iso(ic,ia,ka,kc,breit,sms)
        breit3=breit+breit3
        sms3=sms+sms3
      enddo
        breit1=breit2+breit3
        breit1=breit1*3606.49
        sms1=sms2+sms3
        sms1=sms1*3606.49
        print*,'breit =',breit1
        print*,'sms =',sms1
      endif

cccccccccccc end of S ccccccccccccccc     
 
      if(i.eq.2) then
        print*,'P-',nn
        ia=nbas(1)+(nn-1)
        iccc=nb(1)
      do ic=1,iccc
        ka=1
        kc=-1
        call gen_iso(ic,ia,ka,kc,breit,sms)
        breit2=breit+breit2
        sms2=sms+sms2
      enddo
      if(nsy.ge.4)then
      do j=1,3
        nbass=nbass+nbas(j)
      enddo
        icc=nbas(1)+nbas(2)+nbas(3)+1
        iccc=nbas(1)+nbas(2)+nbas(3)+nb(4)
      do ic=icc,iccc
        ka=1
        kc=2
        call gen_iso(ic,ia,ka,kc,breit,sms)
        breit3=breit+breit3
        sms3=sms+sms3
      enddo
       
      endif
        breit1=breit2+breit3
        breit1=breit1*3609.49
        sms1=sms2+sms3
        sms1=sms1*3609.49
        print*,'breit =',breit1
        print*,'sms =',sms1
      endif
cccccccccccccccccccccccccccccccc
      if(i.eq.3) then
      print*,'P+',nn
      do j=1,2
        nbass=nbass+nbas(j)
      enddo
        ia=nbas(1)+nbas(2)+(nn-1)
        iccc=nb(1)
      do ic=1,iccc
        ka=-2
        kc=-1
        call gen_iso(ic,ia,ka,kc,breit,sms)
        breit2=breit+breit2
        sms2=sms2+sms
      enddo
      if(nsy.ge.4)then
      do j=1,3
        nbass=nbass+nbas(j)
      enddo
        
        icc=nbas(1)+nbas(2)+nbas(3)+1
        iccc=nbas(1)+nbas(2)+nbas(3)+nb(4)
      do ic=icc,iccc
        kc=2
        ka=-2
        call gen_iso(ic,ia,ka,kc,breit,sms)
        breit3=breit+breit3
        sms3=sms+sms3
      enddo

      do j=1,4
        nbass=nbass+nbas(j)
      enddo
        icc=nbas(1)+nbas(2)+nbas(3)+nbas(4)+1
        iccc=nbas(1)+nbas(2)+nbas(3)+nbas(4)+nb(5)
      do ic=icc,iccc
        kc=-3
        ka=-2
        call gen_iso(ic,ia,ka,kc,breit,sms)
        breit4=breit+breit4
        sms4=sms+sms4
      enddo

      endif
        breit1=breit2+breit3+breit4
        breit1=breit1*3609.49
        sms1=sms2+sms3+sms4
        sms1=sms1*3609.49
        print*,'breit =',breit1
        print*,'sms =',sms1
      endif
ccccccccccccccccccccccccccccccc      
      if(i.eq.4) then
       print*,'D-',nn
      do j=1,3
         nbass=nbass+nbas(j)
      enddo
         ia=nbas(1)+nbas(2)+nbas(3)+(nn-2)
         icc=nbas(1)+1
         iccc=nbas(1)+nb(2)
      do ic=icc,iccc
         ka=2
         kc=1
         call gen_iso(ic,ia,ka,kc,breit,sms)
         breit2=breit+breit2
         sms2=sms2+sms
      enddo
      do j=1,2
         nbass=nbass+nbas(j)
      enddo
         icc=nbas(1)+nbas(2)+1
         iccc=nbas(1)+nbas(2)+nb(3)
      do ic=icc,iccc
         kc=-2
         ka=2
         call gen_iso(ic,ia,ka,kc,breit,sms)
         breit3=breit+breit3
         sms3=sms+sms3
      enddo
      if (nsy.ge.6) then
      do j=1,5
         nbass=nbass+nbas(j)
      enddo
         icc=nbas(1)+nbas(2)+nbas(3)+nbas(4)+nbas(5)+1
         iccc=nbas(1)+nbas(2)+nbas(3)+nbas(4)+nbas(5)+nb(6)
      do ic=icc,iccc
         ka=2
         kc=3
         call gen_iso(ic,ia,ka,kc,breit,sms)
         breit4=breit+breit4
         sms4=sms+sms4
      enddo
      endif
         breit1=breit3+breit2+breit4
         breit1=breit1*3609.49
         sms1=sms3+sms2+sms4
         sms1=sms1*3609.49
         print*,'breit =',breit1
         print*,'sms =',sms1
      endif
cccccccccccccccccccccccccccccccccccccccccccc      
       if(i.eq.5) then
       print*,'D+',nn
       do j=1,4
          nbass=nbass+nbas(j)
       enddo

          ia=nbas(1)+nbas(2)+nbas(3)+nbas(4)+1
       do j=1,2
          nbass=nbass+nbas(j)
       enddo
          icc=nbas(1)+nbas(2)+1
          iccc=nbas(1)+nbas(2)+nb(3)
       
       do ic=icc,iccc
          ka=-3
          kc=-2
          call gen_iso(ic,ia,ka,kc,breit,sms)
          breit2=breit+breit2
          sms2=sms2+sms
       enddo
       

       if(nsy.ge.6) then
       do j=1,5
          nbass=nbass+nbas(j)
       enddo
          icc=nbas(1)+nbas(2)+nbas(3)+nbas(4)+nbas(5)+1
          iccc=nbas(1)+nbas(2)+nbas(3)+nbas(4)+nbas(5)+nb(6)
       do ic=icc,iccc
          ka=-3
          kc=3
          call gen_iso(ic,ia,ka,kc,breit,sms)
          breit3=breit+breit3
          sms3=sms+sms3
       enddo
       do j=1,6
           nbass=nbass+nbas(j)
       enddo
          icc=nbas(1)+nbas(2)+nbas(3)+nbas(4)+nbas(5)+nbas(6)+1
          iccc=nbas(1)+nbas(2)+nbas(3)+nbas(4)+nbas(5)+nbas(6)+nb(7)
       do ic=icc,iccc
          ka=-3
          kc=-4
          call gen_iso(ic,ia,ka,kc,breit,sms)
          breit4=breit+breit4
          sms4=sms+sms4
       enddo
       endif


          breit1=breit3+breit2+breit4
          breit1=breit1*3609.49
          sms1=sms3+sms2+sms4
          sms1=sms1*3609.49
          print*,'breit =',breit1
          print*,'sms =',sms1
       endif


       enddo
       endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setgrd
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit real*8 (a-h,o-z)

      common /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
      n = MN-10
      mtp=n
      h = 5.0d-02
      rnt = 2.0d-06
      r(1) = 0.0d0
      rp(1) = rnt
      eph = exp (h)
      ett = 1.d0
!   set up the arrays r, rp, rpor
      do i = 2,MN
      ett = eph*ett
      ettm1 = ett-1.0d0
      r(i) = rnt*ettm1
      rp(i) = rnt*ett
      rpor(i) = ett/ettm1
      enddo
      return
      end
!------------------------------------------------------
      subroutine setqic
!--------------------------------------------------------------------------
! this  subroutine sets up the coefficients for subroutines
! quad, rinti, yzk
!--------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension b13(13,13)
      common/cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /lic13/a(13,13)
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
!------------------------------------------------------------------------
!   five-point newton-cotes coefficients for closed integration.
!   expressed as  rational numbers
!-------------------------------------------------------------------------
      data (c5num(i,2),i = 1,5)/ 251.0d 00, 646.0d 00,
     :       -264.0d 00, 106.0d 00, -19.0d 00/
      data (c5num(i,3),i = 1,5)/ 232.0d 00, 992.0d 00,
     :        192.0d 00,  32.0d 00,  -8.0d 00/
      data (c5num(i,4),i = 1,5)/ 243.0d 00, 918.0d 00,
     :        648.0d 00, 378.0d 00, -27.0d 00/
      data (c5num(i,5),i = 1,5)/ 224.0d 00,1024.0d 00,
     :        384.0d 00,1024.0d 00, 224.0d 00/
      data c5den/ 720.0d 00/
!--------------------------------------------------------------------------
      data (b13( 1,i),i = 1,13) /            -1486442880.0d 00,
     :            5748019200.0d 00, -15807052800.0d 00,
     :           35126784000.0d 00, -59276448000.0d 00,
     :           75873853440.0d 00, -73766246400.0d 00,
     :           54195609600.0d 00, -29638224000.0d 00,
     :           11708928000.0d 00,  -3161410560.0d 00,
     :             522547200.0d 00,    -39916800.0d 00/
      data (b13( 2,i),i = 1,13) /              -39916800.0d 00,
     :            -967524480.0d 00,   2634508800.0d 00,
     :           -4390848000.0d 00,   6586272000.0d 00,
     :           -7903526400.0d 00,   7376624640.0d 00,
     :           -5269017600.0d 00,   2822688000.0d 00,
     :           -1097712000.0d 00,    292723200.0d 00,
     :             -47900160.0d 00,      3628800.0d 00/
      data (b13( 3,i),i = 1,13) /                3628800.0d 00,
     :             -87091200.0d 00,   -684478080.0d 00,
     :            1596672000.0d 00,  -1796256000.0d 00,
     :            1916006400.0d 00,  -1676505600.0d 00,
     :            1149603840.0d 00,   -598752000.0d 00,
     :             228096000.0d 00,    -59875200.0d 00,
     :               9676800.0d 00,      -725760.0d 00/
      data (b13( 4,i),i = 1,13) /                -725760.0d 00,
     :              13063680.0d 00,   -143700480.0d 00,
     :            -476910720.0d 00,   1077753600.0d 00,
     :            -862202880.0d 00,    670602240.0d 00,
     :            -431101440.0d 00,    215550720.0d 00,
     :             -79833600.0d 00,     20528640.0d 00,
     :              -3265920.0d 00,       241920.0d 00/
      data (b13( 5,i),i = 1,13) /                 241920.0d 00,
     :              -3870720.0d 00,     31933440.0d 00,
     :            -212889600.0d 00,   -303937920.0d 00,
     :             766402560.0d 00,   -447068160.0d 00,
     :             255467520.0d 00,   -119750400.0d 00,
     :              42577920.0d 00,    -10644480.0d 00,
     :               1658880.0d 00,      -120960.0d 00/
      data (b13( 6,i),i = 1,13) /                -120960.0d 00,
     :               1814400.0d 00,    -13305600.0d 00,
     :              66528000.0d 00,   -299376000.0d 00,
     :            -148262400.0d 00,    558835200.0d 00,
     :            -239500800.0d 00,     99792000.0d 00,
     :             -33264000.0d 00,      7983360.0d 00,
     :              -1209600.0d 00,        86400.0d 00/
      data (b13( 7,i),i = 1,13) /                  86400.0d 00,
     :              -1244160.0d 00,      8553600.0d 00,
     :             -38016000.0d 00,    128304000.0d 00,
     :            -410572800.0d 00,            0.0d 00,
     :             410572800.0d 00,   -128304000.0d 00,
     :              38016000.0d 00,     -8553600.0d 00,
     :               1244160.0d 00,       -86400.0d 00/
      data (b13( 8,i),i = 1,13) /                 -86400.0d 00,
     :               1209600.0d 00,     -7983360.0d 00,
     :              33264000.0d 00,    -99792000.0d 00,
     :             239500800.0d 00,   -558835200.0d 00,
     :             148262400.0d 00,    299376000.0d 00,
     :             -66528000.0d 00,     13305600.0d 00,
     :              -1814400.0d 00,       120960.0d 00/
      data (b13( 9,i),i = 1,13) /                 120960.0d 00,
     :              -1658880.0d 00,     10644480.0d 00,
     :             -42577920.0d 00,    119750400.0d 00,
     :            -255467520.0d 00,    447068160.0d 00,
     :            -766402560.0d 00,    303937920.0d 00,
     :             212889600.0d 00,    -31933440.0d 00,
     :               3870720.0d 00,      -241920.0d 00/
      data (b13(10,i),i = 1,13) /                -241920.0d 00,
     :               3265920.0d 00,    -20528640.0d 00,
     :              79833600.0d 00,   -215550720.0d 00,
     :             431101440.0d 00,   -670602240.0d 00,
     :             862202880.0d 00,  -1077753600.0d 00,
     :             476910720.0d 00,    143700480.0d 00,
     :             -13063680.0d 00,       725760.0d 00/
      data (b13(11,i),i = 1,13) /                 725760.0d 00,
     :              -9676800.0d 00,     59875200.0d 00,
     :            -228096000.0d 00,    598752000.0d 00,
     :           -1149603840.0d 00,   1676505600.0d 00,
     :           -1916006400.0d 00,   1796256000.0d 00,
     :           -1596672000.0d 00,    684478080.0d 00,
     :              87091200.0d 00,     -3628800.0d 00/
      data (b13(12,i),i = 1,13) /               -3628800.0d 00,
     :              47900160.0d 00,   -292723200.0d 00,
     :            1097712000.0d 00,  -2822688000.0d 00,
     :            5269017600.0d 00,  -7376624640.0d 00,
     :            7903526400.0d 00,  -6586272000.0d 00,
     :            4390848000.0d 00,  -2634508800.0d 00,
     :             967524480.0d 00,     39916800.0d 00/
      data (b13(13,i),i = 1,13) /               39916800.0d 00,
     :            -522547200.0d 00,   3161410560.0d 00,
     :          -11708928000.0d 00,  29638224000.0d 00,
     :          -54195609600.0d 00,  73766246400.0d 00,
     :          -75873853440.0d 00,  59276448000.0d 00,
     :          -35126784000.0d 00,  15807052800.0d 00,
     :           -5748019200.0d 00,   1486442880.0d 00/

      data b13den/479001600.0d 00/

!------------------------------------------------------------------------
!   THIRTEEN-POINT COEFFICIENTS FOR DPBDT
!-----------------------------------------------------------------------
         factor = 1/b13den
         do j = 1,13
            do i = 1,13
               a(i,j) = b13(i,j)*factor
!            print*,'i,j,a(i,j)',i,j,a(i,j)
           enddo
          enddo
!------------------------------------------------------------------------
!                  newton-cotes coefficients for yzk and quad
!-----------------------------------------------------------------------
      factor = h/c5den
      do j = 2,4
         do i = 2,5
            cnc5c(i,j) = factor*c5num(i,j)
         enddo
       enddo
      c1 = factor*c5num(1,5)
      c2 = factor*c5num(2,5)
      c3 = factor*c5num(3,5)
      c4 = c1+c1

      do j = 2,5
         do i = 2,5
            c5num(i,j) = c5num(i,j)/c5den
         enddo
      enddo
      pi = 3.141592653589790d 000
      c = 1.3703598950000d 002
!     c = 1.00000d006
      zero=0.0d 00
      half=0.5d 00
      tenth=0.1d 00
      one=1.0d 00
      two=2.0d 00
      three=3.0d 00
      ten=10.0d 00
      return
      end
!----------------------------------------------------------------------
!      double precision function rint (i,j,k)

!   subroutines called: quad.                                          *
!************************************************************************

!      implicit real*8 (a-h,o-z)
!      common/cons/zero,half,tenth,one,two,three,ten
!     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
!     :  /wave/ppf(MN,MNBAS),qqf(MN,MNBAS),mf(MNBAS)
!      dimension taa(MN)

!      mtp = n

!   tabulate integrand as required for subroutine quad
!   value at first tabulation point is arbitrary

!      taa(1) = zero
!      do 1 l = 2,mtp
!      taa(l)=(r(l)**k)*(pf(l,i)*pf(l,j)+qf(l,i)*qf(l,j))*rp(l)
!    1 continue

!   perform integration

!      call quad (taa,result)
!      rint = result

!      return
!      end

!----------------------------------------------------------------------
       subroutine gen_iso(ic,ia,kka,kkc,br_tot,sms_tot)
!***********************************************************************
       implicit real*8(a-h, o-z)
       dimension gemat1(MN),gemat2(MN),gemat3(MN),gemat4(MN)
     :  ,dppf(MN,MNBAS),dqqf(MN,MNBAS),dpf(MN,MNBAS),dqf(MN,MNBAS),
     :gem1(MN),gem2(MN),gem3(MN),gem4(MN)
       common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /wave/ppf(MN,MNBAS),qqf(MN,MNBAS),mf(MNBAS)
     :  /tatb/tpf(MN),tqf(MN)
     :  /kakc/ka,kc,orba,orbc

        ka=kka
        kc=kkc
       
c       print*,'ic,ia,kc,ka',ic,ia,kc,ka
!       read(*,*)ic
!       read(*,*)ia   
 
       if (kc.eq.-1.or.kc.eq.1) orbc=0.5
       if (kc.eq.-2.or.kc.eq.2) orbc=1.5
       if (kc.eq.-3.or.kc.eq.3) orbc=2.5
       if (kc.eq.-4.or.kc.eq.4) orbc=3.5
       if (kc.eq.-5.or.kc.eq.5) orbc=4.5

       if (ka.eq.-1.or.ka.eq.1) orba=0.5
       if (ka.eq.-2.or.ka.eq.2) orba=1.5
       if (ka.eq.-3.or.ka.eq.3) orba=2.5
       if (ka.eq.-4.or.ka.eq.4) orba=3.5
       if (ka.eq.-5.or.ka.eq.5) orba=4.5
c        print*,'ia,ic,ka,kc,orba,orbc',ia,ic,ka,kc,orba,orbc

       ib=ic
       id=ia
!       print*,'ib,id',ib,id
       mtp=n
       call diff1 (ic)

       do i=2,mtp
       dpf(i,ic)=tpf(i)
       dqf(i,ic)=tqf(i)
       enddo
       call diff1 (id)
       do i=2,mtp
       dppf(i,id)=tpf(i)
       dqqf(i,id)=tqf(i)
       enddo
       do i=1,mtp
       gemat1(i)=zero
       gem1(i)=zero
       gemat2(i)=zero
       gem2(i)=zero
       gemat3(i)=zero
       gem3(i)=zero
       gemat4(i)=zero
       gem4(i)=zero
       enddo
       kk=1
!       print*,'give the value of orbc (j val of ic index) then orba
!     : (j value of ia index)'
!       read*,orbc
!       read*,orba
!       print*,'give the value kappa(ka)and kappa(kc)'      
!       read*,ka
!       read*,kc

       fact=-(2.d0*orbc+1.d0)*
     : clrx (ka,kk,kc)**2
!  for <3s,2p3 I I 2p3,3s> for Na

        kkc=(((ka)*(ka+1))-((kc)*(kc+1)))*half

        kkkc=(((-ka)*(-ka+1))+(kc*(-kc+1)))*half

        kkd=(((kc)*(kc+1))-((ka)*(ka+1)))*half

        kkkd=(((-kc)*(-kc+1))+((ka)*(-ka+1)))*half

c        print*,'kkc,kkkc,kkd,kkkd',kkc,kkkc,kkd,kkkd        

! for <3s,2p1 I I 2p1,3s> for Na
        
       do i=2,mtp
!***********************************
! gemat* are the integran with 1st derivative of PF and QF
!**********************************
       gemat1(i)=ppf(i,ia)*dpf(i,ic)
c       print*,pf(i,ia)

       gem1(i)=(kkc/r(i))*ppf(i,ia)*ppf(i,ic)*rp(i)

       gemat2(i)=qqf(i,ia)*dqf(i,ic)

       gem2(i)=(kkkc/r(i))*qqf(i,ia)*qqf(i,ic)*rp(i)

       gemat3(i)=ppf(i,ic)*dppf(i,ia)

       gem3(i)=(kkd/r(i))*ppf(i,ic)*ppf(i,ia)*rp(i)

       gemat4(i)=qqf(i,ic)*dqqf(i,ia)

       gem4(i)=(kkkd/r(i))*qqf(i,ic)*qqf(i,ia)*rp(i)

       enddo

       call quad (gemat1,rll1)
        rlll1=rll1/h
 
       call quad (gem1,rrl1)

       call quad (gemat2,rss2)
        rsss2=rss2/h

       call quad (gem2,rrs2)

       call quad (gemat3,rll3)
        rlll3=rll3/h

       call quad (gem3,rrl3)

       call quad (gemat4,rss4)
        rsss4=rss4/h

       call quad (gem4,rrs4)

      rllll1=rlll1-rrl1
      rssss2=rsss2-rrs2
      rllll3=rlll3-rrl3
      rssss4=rsss4-rrs4

       rac=(-rrl1-rrs2+rlll1+rsss2)**2
       rca=(-rrl3-rrs4+rlll3+rsss4)**2

       rllc=rllll1*rllll3
       rssc=rssss2*rssss4
       rlsc=rllll1*rssss4
       rslc=rssss2*rllll3

       vint=rllll3+rssss4
       ang1=((-1)**(orbc+0.5))*(SQRT(2.d0*orbc+1.d0))*
     : (SQRT(2.d0*orba+1.d0))*clrx (ka,kk,kc)     

       ang2=fact
       call breit(vint,ang1,ang2,ia,ic,smsres)
       br_tot=smsres
!       print*,'the total breit contribution is ',br_tot
 
       rllr=rllll1**2
       rssr=rssss2**2
       rlsr=rllll1*rssss2

!       print*,'rllc,rssc,rlsc,rslc,fact,ia,ic',rllc,rssc,rlsc,rslc,
!     : fact,ia,ic
!        print*,'rac,rca',rac,rca
       sms_tot=(rllc+rssc+rlsc+rslc)*fact
c       print*,'total cc value is =',sms_tot,fact
       totalr=(rllr+rssr+2*rlsr)*fact
!        print*,'total of rmbpt',totalr

!      r0a=rint(ia,ia,0)
!      r1a=rint(ia,ia,1)
!      r2a=rint(ia,ia,2)
!      r3a=rint(ia,ia,-1)
!      print*,'<AaIBa>,<AaIrIBa>,<AaIr**2IBa>,<AaI1/rIBa>',r0a,r1a
!     :,r2a,r3a

!      r0c=rint(ic,ic,0)
!      r1c=rint(ic,ic,1)
!      r2c=rint(ic,ic,2)
!      r3c=rint(ic,ic,-1)
!      print*,'<AcIBc>,<AcIrIBc>,<AcIr**2IBc>,<AcI1/rIBc>',r0c,r1c
!     :,r2c,r3c

       r4=vinti(ia,ic,ka,kc)
       r5=vinti(ic,ia,kc,ka) 
       r6=-fact*r4*r5
!       print*,'specific grasp, fact',r6,fact
       

       return
       end

!---------------------------------------------------------
        subroutine breit(vintbd,aang1,aang2,iaa,icc,smsbr)
!--------------------------------------------------------

       implicit real*8(a-h, o-z)
       dimension br1(MN),br2(MN)
       common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
     :  /wave/ppf(MN,MNBAS),qqf(MN,MNBAS),mf(MNBAS)
     :  /tatb/tpf(MN),tqf(MN)
     :  /kakc/ka,kc,orba,orbc
     :  /atom/iz



          mtp=n
          call sigma_1 (2,sig1)
          call sigma_2 (2,sig2)
!         print*,'sig1,sig2',sig1,sig2
!         print*,'write the Z value of the atom'
!          read*,z
         z=iz
        do i=2,mtp
          br1(i)=0
        enddo
        do i=2,mtp
      br1(i)=rpor(i)*(-qqf(i,iaa)*ppf(i,icc)*sig1+qqf(i,icc)
     :  *ppf(i,iaa)
     : *sig2)
        enddo
        call quad(br1,resbr1)
        smsbr1=-aang1*vintbd*resbr1*(-z/c)
!       print*,'aang1*vintbd*resbr1,br1',aang1,vintbd,resbr1,br1
        do i=2,mtp
        br2(i)=0
         enddo
         do i=2,mtp
       br2(i)=rpor(i)*(qqf(i,iaa)*ppf(i,icc)-qqf(i,icc)*ppf(i,iaa))
         enddo
        call quad(br2,resbr2)
        smsbr2=-aang2*vintbd*resbr2*(-z/c)
         
        smsbr= smsbr1+smsbr2
       
!       print*,'breit contribution first',smsbr1
!       print*,'breit contribution second',smsbr2

        return
        end
!----------------------------------------------------------
      FUNCTION VINTI (k,l,mk,nk)
!----------------------------------------------------------
             implicit real*8(a-h, o-z)
!      LOGICAL LDBPR
!      CHARACTER*2 NH


      dimension ttpf(MN),ttqf(MN)
      common /cons/zero,half,tenth,one,two,three,ten
     :  /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :  /wave/ppf(MN,MNBAS),qqf(MN,MNBAS),mf(MNBAS)
     :  /tatb/tpf(MN),tqf(MN)


!      mtp = MAX (mf(ia),mf(ic))

!   Piece involving derivatives

      mtp=n
      call diff1(l)
      do i=1,mtp
      ttpf(i) = 0
      enddo
      DO 1 i = 2,mtp
         ttpf(i) = ppf(i,k)*tpf(i)+qqf(i,k)*tqf(i)
    1 CONTINUE
      call quad(ttpf,piece1)
      piece1 = piece1/h
!      print*,'piece1',piece1

!   Pieces not involving derivatives

      kpj = mk
      kpk = nk
      ifact1 = kpj*(kpj+1)-kpk*(kpk+1)
      fact1 =half*(ifact1)
      ifact2 = -kpj*(-kpj+1)+kpk*(-kpk+1)
      fact2 =half*(ifact2)
     
!       print*,'fact1,fact2',mk,nk,fact1,fact2
       
      do i=1,mtp
      ttqf(i) = 0
      enddo
      DO 2 i = 2,mtp
         ttqf(i) = rpor(i)*( fact1*ppf(i,k)*ppf(i,l)
     :                    +fact2*qqf(i,k)*qqf(i,l))
    2 CONTINUE
      call quad (ttqf,piece2)
!      print*,'piece2',piece2

      vinti = piece1-piece2
!      print*,'vinti',vinti

      RETURN
      END
!************************************************************************
      subroutine sigma_1 (ipar,apart)
!                                                                      *
!   This  routine computes                                             *
!   [-kappa_a || sigma^(1) ||kappa_b]  if IPAR = 1                     *
!   <-kappa_a || sigma^(1) ||kappa_b>  if IPAR = 2                     *
!                                                                      *
!   Written by G. Gaigalas and E. Gaidamauskas                         *
!                                         Last revision: 16 Doc 2009   *
!                                                                      *
!***********************************************************************

      implicit real*8 (a-h,o-z)

      common/kakc/ka,kc,orba,orbc

!   Set KAP1 and KAP2

      kap1 = ka
      kap2 = kc
!      print*,'ka,kc',ka,kc

!   Determine the l quantum numbers

      if (kap1 .gt. 0) then
         l1 =  kap1
      else
         l1 = -kap1-1
      endif

      if (kap2 .gt. 0) then
         l2 =  kap2
      else
         l2 = -kap2-1
      endif

!   Determine the j quantum numbers and l_1 tilda

      j1 = iabs(ka)*2-1
      j2 = iabs(kc)*2-1
      l1_tilda = j1-l1
      if (l1_tilda.eq.l2) then

!   Determine the Racah W coefficients.

         call racah (1,j1,1,j2,2*l1_tilda,2,rac)
!         print*,'rac',rac
         if (mod(1+j1+1+j2,4).eq.0) rac=-rac
         apart =  rac*sqrt(6.0*(j2+1))
      else
         apart = 0.0
      end if
        
!   Determine the phase factori
        
      if(mod(2*l1_tilda+1+j1+2,4).eq.0) apart=-apart
      if(ipar.eq.2) then
         apart = dsqrt(dble(j1+1))*apart
      end if
      return
      end

!************************************************************************
!                                                                      *
      subroutine sigma_2 (ipar,apart)
!                                                                      *
!   This  routine computes                                             *
!   [-kappa_a || sigma^(2) ||kappa_b] if IPAR = 1                      *
!   <-kappa_a || sigma^(2) ||kappa_b> if IPAR = 2                      *
!                                                                      *
!                                                                      *
!   Written by G. Gaigalas and E. Gaidamauskas                         *
!                                         Last revision: 16 Doc 2009   *
!                                                                      *
!************************************************************************

      implicit real*8 (a-h,o-z)
!
      common/kakc/ka,kc,orba,orbc
!
!   Set KAP1 and KAP2

      kap1 = ka
      kap2 = kc

!   Determine the l quantum numbers

!      if (kap1 .gt. 0) then
      if (kap1 .gt. 0) then
         l1 =  kap1
      else
         l1 = -kap1-1
      endif
     
      if (kap2 .gt. 0) then
         l2 =  kap2
      else
         l2 = -kap2-1
      endif
   
!   Determine the j quantum numbers and l_2 tilda
     
      j1 = iabs(ka)*2-1
      j2 = iabs(kc)*2-1
      l2_tilda = j2 - l2
      if (l2_tilda .eq. l1) then
   
!   Determine the Racah W coefficients.
        
         call racah (1,j1,1,j2,2*l1,2,rac)
!         print*,'rac2',rac
         if (mod(1+j1+1+j2,4).eq.0) rac=-rac
         apart =  rac*sqrt(6.0*(j2+1))
      else
         apart = 0.0
      end if

!   Determine the phase factori

      if(mod(2*l1+1+j1+2,4).eq.0) apart=-apart
      if(ipar .eq. 2) then
         apart = dsqrt(dble(j1+1))*apart
      end if
      return
      end


!-----------------------------------------------------------
      subroutine diff1(J)
!-----------------------------------------------------------
      implicit real*8 (a-h, o-z)

      COMMON/cons/zero,half,tenth,one,two,three,ten
     :      /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     :      /lic13/a(13,13)
     :      /tatb/tpf(MN),tqf(MN)
     :      /wave/ppf(MN,MNBAS),qqf(MN,MNBAS),mf(MNBAS)

      equivalence (a1,a(7,1)),(a2,a(7,2)),(a3,a(7,3)),
     :            (a4,a(7,4)),(a5,a(7,5)),(a6,a(7,6))

!   COMPUTE DERIVATIVE IN THREE SEPARATE REGIONS

!   FIRST, POINTS 1 TO 6

      do i = 1,6
         hdpbdt = 0.0d 00
         hdqbdt = 0.0d 00
         do k = 1,13
            aik = a(i,k)
            hdpbdt = hdpbdt+aik*ppf(k,J)
            hdqbdt = hdqbdt+aik*qqf(k,J)
         enddo
         tpf(i) = hdpbdt
         tqf(i) = hdqbdt
      enddo
      do i = 7,n-6
      tpf(i)= a1*(ppf(i-6,J)-ppf(i+6,J))+a2*(ppf(i-5,J)-ppf(i+5,J))
     :       +a3*(ppf(i-4,J)-ppf(i+4,J))+a4*(ppf(i-3,J)-ppf(i+3,J))
     :       +a5*(ppf(i-2,J)-ppf(i+2,J))+a6*(ppf(i-1,J)-ppf(i+1,J))
      tqf(i)= a1*(qqf(i-6,J)-qqf(i+6,J))+a2*(qqf(i-5,J)-qqf(i+5,J))
     :      +a3*(qqf(i-4,J)-qqf(i+4,J))+a4*(qqf(i-3,J)-qqf(i+3,J))
     :      +a5*(qqf(i-2,J)-qqf(i+2,J))+a6*(qqf(i-1,J)-qqf(i+1,J))
       enddo

!   LAST, POINTS N-5 TO N

      do i = n-5,n
         irow = i-n+13
         hdpbdt = 0.0d 00
         hdqbdt = 0.0d 00
         do k = 1,13
            aik = a(irow,k)
            loc = n-13+k
            hdpbdt = hdpbdt+aik*ppf(loc,J)
            hdqbdt = hdqbdt+aik*qqf(loc,J)
         enddo
         tpf(i) = hdpbdt
         tqf(i) = hdqbdt
       enddo
       return
      end

!*********************************************************************
       subroutine racah (i,j,k,l,m,n,rac)
!                                                                      *
!   SUBROUTINE  to calculate Racah coefficients. The arguments I, J,   *
!   K, L, M, N should be twice their actual value. Works for integer   *
!   and  half-integer  values of  angular momenta. The routine makes   *
!   use of the GAM  array, thus  SUBROUTINE FACTT must be called be-   *
!   fore this routine is used.                                         *
!                                                                      *
!   Written by N S Scott                    Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************

      implicit real*8(a-h,o-z)

!      PARAMETER (MFACT = 500)

      common/facts/gam(MFACT)
              


      call factt(gam)
      j1 = i+j+m
      j2 = k+l+m
      j3 = i+k+n
      j4 = j+l+n
      if (((2*max(i,j,m)-j1).gt.0).or.(mod(j1,2).ne.0)) goto 2
      if (((2*max(k,l,m)-j2).gt.0).or.(mod(j2,2).ne.0)) goto 2
      if (((2*max(i,k,n)-j3).gt.0).or.(mod(j3,2).ne.0)) goto 2
      if (((2*max(j,l,n)-j4).gt.0).or.(mod(j4,2).ne.0)) goto 2
      goto 1
   2  rac = 0.0
!      print*,'11111111111111111'
      return

   1  continue
!      print*,'222222222222222222'
      j1 = j1/2
      j2 = j2/2
      j3 = j3/2
      j4 = j4/2
      j5 = (i+j+k+l)/2
      j6 = (i+l+m+n)/2
      j7 = (j+k+m+n)/2
      numin = max (j1,j2,j3,j4)+1
      numax = min (j5,j6,j7)+1
      rac = 1.0
      icount = 0

      if (numin.EQ.numax) goto 4
      numin = numin+1

      do kk = numin,numax
         ki = numax-icount
         rac = 1.0
     :       -(rac*dble(ki*(j5-ki+2)*(j6-ki+2)*(j7-ki+2))/
     :      dble((ki-1-j1)*(ki-1-j2)*(ki-1-j3)*(ki-1-j4)))
         icount = icount+1
        enddo

      numin = numin-1
   4  rac = rac*((-1.0)**(j5+numin+1))
     : *exp((gam(numin+1)-gam(numin-j1)
     : -gam(numin-j2)-gam(numin-j3)-gam(numin-j4)-gam(j5+2-numin)
     : -gam(j6+2-numin)-gam(j7+2-numin))+((gam(j1+1-i)+gam(j1+1-j)
     : +gam(j1+1-m)-gam(j1+2)+gam(j2+1-k)+gam(j2+1-l)+gam(J2+1-m)
     : -gam(j2+2)+gam(j3+1-i)+gam(j3+1-k)+gam(j3+1-n)-gam(j3+2)
     : +gam(j4+1-j)+gam(j4+1-l)+gam(j4+1-n)-gam(j4+2))
     : *0.5))

      return
      end

!---------------------------------------------------------
      subroutine quad(aa,result)
!--------------------------------------------------------
      implicit real*8(a-h, o-z)
      dimension aa(MN)
      common/cons/zero,half,tenth,one,two,three,ten
     : /grid/steplin,r(MN),rp(MN),rpor(MN),rnt,h,hp,n
     : /ncc/c1,c2,c3,c4,cnc5c(5,5),c5num(5,5),pi,c
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
!------------------------------------------------------------------------
      real*8 function clrx (kappaa,k,kappab)
!-----------------------------------------------------------------------
      implicit real*8(a-h, o-z)
      dimension gam(500)
      common/cons/zero,half,tenth,one,two,three,ten

       call factt(gam)
!----------------------------------------------------------------------
!               determine the absolute values of the kappas
!---------------------------------------------------------------------- 
        ka =iabs (kappaa)      
        kb =iabs (kappab)
!----------------------------------------------------------------------
!                   perform the triangularity check
!----------------------------------------------------------------------
         if ((iabs(ka-kb) .le. k) .and. (ka+kb-1 .ge. k)) then
!----------------------------------------------------------------------
!   triangularity satisfied; compute the 3j coefficient
!   begin with the logarithm of the square of the leading term
!----------------------------------------------------------------------
         exptrm = -dlog (dble (ka*kb))
         kapkb = ka+kb
         kabkp = kapkb+k
         kamkb = ka-kb
         kbmka = kb-ka
         exptrm = half*(exptrm+gam(kapkb-k  )+gam(kamkb+k+1)
     :                        +gam(kbmka+k+1)-gam(kabkp  +1) )
!----------------------------------------------------------------------
!   the remainder depends on the parity of ka+kb+k
!----------------------------------------------------------------------
         if (mod (kabkp,2) .eq. 0) then
!----------------------------------------------------------------------
!   computation for even parity case
!   include the phase factor: a minus sign if necessary
!----------------------------------------------------------------------
            if (mod (3*kabkp/2,2) .eq. 0) then
               clrx =  one
            else
               clrx = -one
            endif
!----------------------------------------------------------------------
!   include the contribution from the factorials
!----------------------------------------------------------------------
            exptrm = exptrm+gam((kabkp  +2)/2)-gam((kapkb-k  )/2)
     :                     -gam((kamkb+k+2)/2)-gam((kbmka+k+2)/2)

         else
            if (mod ((3*kabkp-1)/2,2) .eq. 0) then
               clrx =  one
            else
               clrx = -one
            endif
!----------------------------------------------------------------------
!   include the contribution from the factorials
!----------------------------------------------------------------------
            exptrm = exptrm+gam((kabkp  +1)/2)-gam((kapkb-k+1)/2)
     :                     -gam((kamkb+k+1)/2)-gam((kbmka+k+1)/2)

         endif
!----------------------------------------------------------------------
!   final assembly
!----------------------------------------------------------------------
         clrx = clrx*dexp (exptrm)

          else
!----------------------------------------------------------------------
!   triangularity violated; set the coefficient to zero
!----------------------------------------------------------------------
         clrx = zero

        endif

        return
        end
!-------------------------------------------------------------
      subroutine factt(gam)
!---------------------------------------------------------------------
! calculates the logs  of factorials required by the racah coefficient
! routine  dracah. written by n.s. scott.
!---------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension gam(500)
      common /cons/zero,half,tenth,one,two,three,ten
      data thirty/3.0d 01/

      gam(1) = one
      gam(2) = one
      x = two
      do i = 3,30
         gam(i) = gam(i-1)*x
         x = x+one
      enddo
      do  i = 1,30
         gam(i) = dlog(gam(i))
      enddo
      x = thirty
      do i = 31,500
         gam(i) = gam(i-1)+dlog(x)
         x = x+one
      enddo
      return
       end
