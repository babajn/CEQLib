      program spectra

      implicit none
      integer nmaxspec,nmaxlambda,nmaxtimes,nmaxtimesimpr
      integer nmaxlines,nmaxZ,nmaxlambdaext,i912
      integer nmaxtimesSSP,nmaxparam,nmaxCM
      parameter(nmaxparam=99)
      parameter(nmaxtimesSSP=600)
      parameter(nmaxspec=4500,nmaxCM=100,nmaxlambda=1300)
      parameter(nmaxtimes=20002,nmaxtimesimpr=100)
      parameter(nmaxlines=61,nmaxlambdaext=300)
      parameter(nmaxZ=7)
      character*72 fileSSPZ(nmaxZ),filespectra,a,fileexch
      character*72 fileSSPs,filescenarios
      character*100 header(nmaxZ,4)
      character*72 fileSFR,answerinfall
      character*72 answerneb,answerext,answerwind,answerZ
      character*72 char_fSNIa
      integer p,infall,intgr,i0,i01,i02,i03,i04,j0
      integer i,j,k,ntimes,nlambdacont,nlines,nZstellib
      integer nspec,nlambda,nlambdaext,typeSFR,nCM
      integer ntimesimpr,jimpr,timeSSP(nmaxtimesSSP)
      integer timeimpr(nmaxtimesimpr)
      integer invtime(nmaxtimes),q,error
      real time(nmaxtimes),twind
      real ejecta(nmaxZ,nmaxtimes),ZSSP(nmaxZ),Zexch
      real ejectaZ(nmaxZ,nmaxtimes)
      real fluxspec(-nmaxCM:nmaxspec),inclin,fluxbolZ,Zbol
      real SFR(nmaxtimes),dsigmagas,SFRlum(nmaxtimes)
      double precision dejecta,dejectaZ
      real sigmaBHNS(nmaxtimes),agestars(nmaxtimes)
      double precision dsigmaBHNS,dsigmaWD
      real massBHNS(nmaxZ,nmaxtimes),massBHNSint
      real sigmaWD(nmaxtimes),agebol,SFRmax,fluxbolage
      real massWD(nmaxZ,nmaxtimes),massWDint,fsub
      real massBHNS_IS,massBHNS_CB,massWD_IS,massWD_CB
      real massalive(nmaxZ,nmaxtimes)
      real ejecta_IS,ejecta_CB,ejectaZ_IS,ejectaZ_CB
      real nSNIaint,nSNII_IS,nSNII_CB
      real fluxstel(-nmaxCM:nmaxspec,nmaxlambda)
      real fluxgal(nmaxlambda)
      real x,sigmastars(nmaxtimes)
      real lambda(nmaxlambda)
      double precision NLymtot,NLym(nmaxZ,nmaxtimes),dble
      double precision NHItot,NHItotprev
      real ejectaint
      real g1,g2,g3,g4,gam1(100),gam2(100),gam3(100),gam4(100)
      real gline(nmaxlines),lambdaline(nmaxlines),lambdacont(100)
      real c,HeI,HeII,alphaTe,fline(nmaxlines),fneb(nmaxlambda)
      real sigmagas(nmaxtimes),sigmasub(nmaxtimes)
      real corextcont(nmaxlambda),corextinclincont(nmaxlambda)
      real corextline(nmaxlines),corextinclinline(nmaxlines)
      real Zgas(nmaxtimes),ZSFR(nmaxtimes)
      real flinetot(nmaxlines),f
      real fluxbol,nSNIItot,fluxext
      real nSNII(nmaxZ,nmaxtimes)
      real Zsol,Lsol,sigmaZ,dsigmaZ,ejectaZint,nSNIIint
      real tau,fSNIa,Zinfall
      real tinfall,nSNIa(nmaxZ,nmaxtimes),Zstars(nmaxtimes)
      real nSNIatot
      real Mgal(nmaxtimes)
      integer nused,iused(nmaxspec),nZ
      real Zext(5),coeffext(5,2),frac(5)
      real albedo,asym
      real lambdaext(nmaxlambdaext),tauext(nmaxlambdaext,2)
      real albedoext(nmaxlambdaext,2),asymext(nmaxlambdaext,2)
      real fluxSSP(nmaxZ,nmaxtimesSSP,nmaxlambda)
      real beta(nmaxtimes),alpha(nmaxtimes)
      real fluxbolSSP(nmaxZ,nmaxtimesSSP)
      integer iZinf(nmaxtimes),iZsup(nmaxtimes)
      integer nparam,iparam
      real param1,param2,param(nmaxparam)
      logical excessSFR
      
      integer i1,i2,i3,i4

      integer nslab1,nslab2,nslab3,nslab4
      parameter(nslab1=17,nslab2=9,nslab3=6,nslab4=10)
      real tauslab(nslab1),albslab(nslab2)
      real asymslab(nslab3),inclinslab(nslab4)
      real emergslab(nslab1,nslab2,nslab3)
      real emerginclinslab(nslab1,nslab2,nslab3,nslab4)

      integer nK1,nK2,nK3
      parameter(nK1=12,nK2=9,nK3=11)
      real tauKing(nK1),albKing(nK2),asymKing(nK3)
      real emergKing(nK1,nK2,nK3)

      integer iZinfk,iZsupk
      double precision agestarsi,sigmastarsi,Zstarsi
      double precision alphak,SFRk,SFRalphak,SFR1malphak
      integer invtimeq,invtimeqp1
      double precision betaq,betaqalphak,betaq1malphak
      double precision onembetaqalphak,onembetaq1malphak
      double precision fluxgalj,ffnebj
*     Lsol: erg/s; c: A/s.
      parameter(Zsol=0.02,Lsol=3.826e33,c=2.99792458e18)

      real interplinlog
      integer nSpitzer
      real taudustSpitzer(20),ySpitzer(20),y,taudust

      real NH,col_dens,mH,coeffextZ(2),tauV
      parameter(mH=1.67e-24)

      integer ntimesfile
      real timefile(nmaxtimes),SFRfile(nmaxtimes),Zfile(nmaxtimes)

      integer codeext,codeZ

***** Opening of the file containing the scenarios.

      error=1
      write(*,*) 'Name of the file containing the parameters', 
     $     ' of the star formation scenarios?'
      read(*,'(a)') filescenarios
      do while (error.ne.0) 
         open(50,file=filescenarios,status='old',iostat=error)
         if (error.ne.0) then
            write(*,*) 'This file does not exist!'
            write(*,*) 'Name of the file containing the',
     $           ' parameters of the star formation scenarios?'
            read(*,'(a)') filescenarios
         end if
      end do

***** Reading of the names of the SSPs files.

      read(50,'(a)') fileSSPs
      i=index(fileSSPs,':')
      i=i+1
      do while (fileSSPs(i:i).eq.' ') 
         i=i+1
      end do
      fileSSPs=fileSSPs(i:)
      open(80,file=fileSSPs,status='old')
      nZ=0
      do i=1,nmaxZ
         read(80,'(a)',end=999) fileSSPZ(i)
         open(20,file=fileSSPZ(i),status='old')
         read(20,'(a)') a
         read(20,'(29x,f5.4)') ZSSP(i)
         close(20)
         nZ=nZ+1
      end do
 999  close(80)     
* The files are ordered according to increasing metallicity.
      do i=1,nZ-1 
         do j=i+1,nZ
            if (ZSSP(j).lt.ZSSP(i)) then
               Zexch=ZSSP(i)
               ZSSP(i)=ZSSP(j)
               ZSSP(j)=Zexch
               fileexch=fileSSPZ(i)
               fileSSPZ(i)=fileSSPZ(j)
               fileSSPZ(j)=fileexch
            end if
         end do
      end do

      read(50,'(a)') char_fSNIa
      i=index(char_fSNIa,':')
      read(char_fSNIa(i+1:),*) fSNIa

***** Reading of the times for which the spectra are to be computed.

      open(30,status='old',file='ages.dat')
      ntimesimpr=0
      error=0
      do while (error.eq.0)
         read(30,*,iostat=error) timeimpr(ntimesimpr+1)
         if (error.eq.0) ntimesimpr=ntimesimpr+1
      end do
      close(30)

***** Reading of the stellar spectra.

      open(10,status='old',file='stellibLCBcor.dat')
      read(10,*) nZstellib,nlambda
      nspec=0
      do i=1,nZstellib
         read(10,*) j         
         do k=1,j
            read(10,'(a)') a
         end do
         nspec=nspec+j
      end do
      read(10,*) (lambda(k),k=1,nlambda)
      do i=1,nspec
         read(10,*) (fluxstel(i,k),k=1,nlambda)
      end do
      close(10)

      open(10,status='old',file='stellibCM.dat')
      read(10,*) nCM
      do i=1,nCM
         read(10,'(a)') a         
      enddo
      read(10,*) (lambda(k),k=1,nlambda)
      do i=1,nCM
         read(10,*) (fluxstel(-i,k),k=1,nlambda)
      end do
      close(10)

***** Position of the Lyman break: lambda(i912)<912 Angstroems<=lambda(i912+1)

      i912=1
      do while (lambda(i912+1).lt.912.)
         i912=i912+1
      end do

***** Calculations for the nebular continuum.

      open(60,status='old',file='HII.dat')
      read(60,*) HeI,HeII,alphaTe
      read(60,*) nlambdacont   
      do i=1,nlambdacont
         read(60,*) lambdacont(i),gam1(i),gam2(i),
     $        gam3(i),gam4(i)
      end do      
      j=1
      do i=1,nlambda
         do while (((lambdacont(j)-lambda(i))
     $        *(lambdacont(j+1)-lambda(i)).gt.0.)
     $        .and.(j.lt.nlambdacont))
            j=j+1
         end do
         if ((lambdacont(j)-lambda(i))
     $        *(lambdacont(j+1)-lambda(i)).le.0.) then
            g1=interplinlog(1./lambdacont(j),1./lambdacont(j+1),
     $           gam1(j),gam1(j+1),1./lambda(i))
            g2=interplinlog(1./lambdacont(j),1./lambdacont(j+1),
     $           gam2(j),gam2(j+1),1./lambda(i))
            g3=interplinlog(1./lambdacont(j),1./lambdacont(j+1),
     $           gam3(j),gam3(j+1),1./lambda(i))
            g4=interplinlog(1./lambdacont(j),1./lambdacont(j+1),
     $           gam4(j),gam4(j+1),1./lambda(i))
         endif
*     -1.9: gaunt factor ~ nu**(-0.1) at long wavelength
         if (lambda(i).ge.lambdacont(nlambdacont)) then
            g1=gam1(nlambdacont)
     $           *(lambda(i)/lambdacont(nlambdacont))**0.1
            g2=gam2(nlambdacont)
     $           *(lambda(i)/lambdacont(nlambdacont))**0.1
            g3=gam3(nlambdacont)
     $           *(lambda(i)/lambdacont(nlambdacont))**0.1
            g4=gam4(nlambdacont)
     $           *(lambda(i)/lambdacont(nlambdacont))**0.1
         endif
         fneb(i)=1.d-40*(g1+g2+HeI*g3+HeII*g4)*c/alphaTe
     $        /lambda(i)**2
      enddo

***** Calculations for the nebular lines.

      read(60,*) nlines
      do i=1,nlines
*     gline: line strength relative to Hbeta (first line).
         read(60,*) lambdaline(i),gline(i)
         fline(i)=1.244d-25*gline(i)/alphaTe
      end do
      close(60)

***** Reading of the extinction parameters for the slab model

      open(70,file='slab.dat',status='old')
      do i1=1,nslab1
         do i2=1,nslab2
            read(70,*) tauslab(i1),albslab(i2),
     $           (asymslab(i3),i3=1,nslab3)
            read(70,*) (emergslab(i1,i2,i3),i3=1,nslab3)
            do i3=1,nslab3
               emergslab(i1,i2,i3)=log(emergslab(i1,i2,i3))
            end do
            do i4=1,nslab4
               read(70,*) inclinslab(i4),(emerginclinslab(i1,i2,i3,i4),
     $              i3=1,nslab3)
               do i3=1,nslab3
                  emerginclinslab(i1,i2,i3,i4)=
     $                 log(emerginclinslab(i1,i2,i3,i4))
               end do
            end do
         end do
      end do
      close(70)

***** Reading of the extinction parameters for the King model

      open(10,file='King.dat',status='old')
      do i1=1,nK1
         do i2=1,nK2
            do i3=1,nK3
               read(10,*) tauKing(i1),albKing(i2),asymKing(i3),
     $              emergKing(i1,i2,i3)
               emergKing(i1,i2,i3)=log(emergKing(i1,i2,i3))
            end do
         end do
      end do
      close(10)

***** Reading of the dust properties (absorption, albedo, asymmetry)
***** for graphites and silicates. Different C/Si ratios are used depending
***** on the metallicity. 

      open(90,file='dust.dat',status='old')
      read(90,'(a)') a
      read(90,'(a)') a
      do i=2,4
         read(90,*) Zext(i),(coeffext(i,j),j=1,2)
      end do
      Zext(1)=0.
      coeffext(1,1)=coeffext(2,1)
      coeffext(1,2)=coeffext(2,2)
      Zext(5)=1.
      coeffext(5,1)=coeffext(4,1)
      coeffext(5,2)=coeffext(4,2)
      do i=1,5
         frac(i)=coeffext(i,1)/(coeffext(i,1)+coeffext(i,2))
      enddo
      read(90,'(a)') a
      read(90,*) nlambdaext
      do i=1,nlambdaext
         read(90,*) lambdaext(i),(tauext(i,j),albedoext(i,j),
     $        asymext(i,j),j=1,2)
      end do
      close(90)

      open(10,file='Spitzer.dat',status='old')
      read(10,*) nSpitzer
      do i=1,nSpitzer
         read(10,*) taudustSpitzer(i),ySpitzer(i)
      enddo
      close(10)

***** Reading of the flux emitted for an instantaneous burst in
***** the stellar spectra.

      do p=1,nZ
         open(20,status='old',file=fileSSPZ(p))
         do i=1,4
            read(20,'(a)') header(p,i)
         end do
         read(20,*) ntimes
         read(20,*) timeSSP(1),nused,x,fluxbolSSP(p,1),NHItotprev
         read(20,'(a)') a
         read(20,'(a)') a
         read(20,*) (iused(k),fluxspec(iused(k)),k=1,nused)
         do j=1,nlambda
            dble=0.d0
            do k=1,nused
               intgr=iused(k)
               dble=dble
     $              +fluxspec(intgr)*fluxstel(intgr,j)
            end do
            fluxSSP(p,1,j)=dble
         end do
         massalive(p,1)=1.
         do i=1,ntimes-1
            read(20,*) timeSSP(i+1),nused,x,fluxbolSSP(p,i+1),NHItot
            read(20,*) nSNII_IS,nSNII_CB,nSNIaint,
     $           massBHNS_IS,massBHNS_CB,massWD_IS,massWD_CB
            read(20,*) 
     $           ejecta_IS,ejecta_CB,ejectaZ_IS,ejectaZ_CB
            nSNIIint=(1.-fSNIa)*nSNII_IS+fSNIa*nSNII_CB
            nSNIaint=fSNIa*nSNIaint
            massBHNSint=(1.-fSNIa)*massBHNS_IS+fSNIa*massBHNS_CB
            massWDint=(1.-fSNIa)*massWD_IS+fSNIa*massWD_CB
            ejectaint=(1.-fSNIa)*ejecta_IS+fSNIa*ejecta_CB
            ejectaZint=(1.-fSNIa)*ejectaZ_IS+fSNIa*ejectaZ_CB
            read(20,*) (iused(k),fluxspec(iused(k)),k=1,nused)            
            do j=1,nlambda
               dble=0.d0
               do k=1,nused
                  intgr=iused(k)
                  dble=dble
     $                 +fluxspec(intgr)*fluxstel(intgr,j)
               end do
               fluxSSP(p,i+1,j)=dble
            end do
            do k=timeSSP(i)+1,timeSSP(i+1)    
               invtime(k)=i
               beta(k)=(timeSSP(i+1)-k+1.)/(timeSSP(i+1)-timeSSP(i))
               NLym(p,k)=NHItotprev+(NHItot-NHItotprev)
     $              *(k-timeSSP(i)-1.)/(timeSSP(i+1)-timeSSP(i))
               nSNII(p,k)=nSNIIint
               nSNIa(p,k)=nSNIaint
               ejecta(p,k)=ejectaint
               ejectaZ(p,k)=ejectaZint
               massBHNS(p,k)=massBHNSint
               massWD(p,k)=massWDint
*     Mass in stars still alive.
               massalive(p,k+1)=massalive(p,k) 
     $              -ejecta(p,k)-massBHNS(p,k)-massWD(p,k)
            end do
            NHItotprev=NHItot
         end do
         close(20)
      end do
      invtime(timeSSP(ntimes)+1)=ntimes-1
      beta(timeSSP(ntimes)+1)=0.

***** Reading of the star formation scenarios.
      
      read(50,'(a)',iostat=error) a
      do while (error.eq.0)
         read(50,'(a)') filespectra
         i=index(filespectra,':')
         i=i+1
         do while (filespectra(i:i).eq.' ') 
            i=i+1
         end do
         filespectra=filespectra(i:)
         open(40,status='new',file=filespectra,iostat=error)        
         if (error.ne.0) then
            write(*,*) filespectra(1:index(filespectra,' ')),
     $           'already existed;'
            do while (error.ne.0)
               filespectra=filespectra(1:index(filespectra,' ')-1)//'+'
               open(40,status='new',file=filespectra,iostat=error)
               if (error.eq.0) then
                  write(*,*) 'a new file named ',
     $                 filespectra(1:index(filespectra,' ')),
     $                 'has been created.'
                  write(*,*) ' '
               endif
            end do
         endif
         
         write(40,'(a)') header(1,3)
         write(40,'(a)') header(1,4)
         do p=1,nZ
            do i=1,2
               write(40,'(a)') header(p,i)
            end do
         end do
         write(40,'(a)') char_fSNIa

         call read_write_real(50,40,Zgas(1))

         read(50,'(a)') answerinfall
         write(40,'(a)') answerinfall
         if (answerinfall.eq.'Infall') then
            infall=1
            call read_write_real(50,40,tinfall)
            call read_write_real(50,40,Zinfall)
         else
            tinfall=1.
            Zinfall=0.
            infall=0
         end if

         call read_write_integer(50,40,typeSFR)
         if (typeSFR.ge.1.and.typeSFR.le.3) then            
            call read_write_real(50,40,param1)
            call read_write_real(50,40,param2)
         end if
         if (typeSFR.ge.10) then
            call read_write_integer(50,40,nparam)
            do iparam=1,nparam
               call read_write_real(50,40,param(iparam))
            end do
         end if
         if (typeSFR.le.-1) then
            call read_write_char(50,40,fileSFR)
         end if

         if (typeSFR.ge.-1) then
            read(50,'(a)') answerZ
            write(40,'(a)') answerZ
            if (index(answerZ,'No').ne.0) then
               codeZ=0
            else
               codeZ=1
            endif
            if (codeZ.eq.0) then
               call read_write_real(50,40,ZSFR(1))            
            else
               ZSFR(1)=Zgas(1)
            end if
         endif 

         call read_write_real(50,40,fsub)
         read(50,'(a)') answerwind
         write(40,'(a)') answerwind
         if (answerwind.eq.'Galactic winds') then
            call read_write_real(50,40,twind)
         else
            twind=nmaxtimes+1.
         end if         

         read(50,'(a)') answerneb
         write(40,'(a)') answerneb
         read(50,'(a)') answerext
         write(40,'(a)') answerext
         if (index(answerext,'No').ne.0) codeext=0
         if (index(answerext,'spheroidal').ne.0) codeext=1
         if (index(answerext,'averaged').ne.0) codeext=2
         if (index(answerext,'specific').ne.0) codeext=3
         if (codeext.eq.3) call read_write_real(50,40,inclin)

****** Initialization.

         Mgal(1)=(1.-infall)
         sigmagas(1)=Mgal(1)       
         sigmaZ=sigmagas(1)*Zgas(1)
         sigmaBHNS(1)=0.
         sigmaWD(1)=0.
         excessSFR=.false.
         sigmasub(1)=0.

***** Reading of the SFR (and possibly Z also) in fileSFR.         

         if (typeSFR.eq.-1) then
            open(10,file=fileSFR,status='old')
            ntimesfile=0
            error=0
            do while (error.eq.0)
               read(10,*,iostat=error) timefile(ntimesfile+1),
     $              SFRfile(ntimesfile+1)
               if (error.eq.0) ntimesfile=ntimesfile+1
            end do
            close(10)
         endif
         if (typeSFR.eq.-2) then
            open(10,file=fileSFR,status='old')
            ntimesfile=0
            error=0
            do while (error.eq.0)
               read(10,*,iostat=error) timefile(ntimesfile+1),
     $              SFRfile(ntimesfile+1),Zfile(ntimesfile+1)
               if (error.eq.0) ntimesfile=ntimesfile+1
            end do
            close(10)
         endif

***** Convolution of the SSPs with the star formation rate.

         i0=0
         do i=1,nmaxtimes-1            
            time(i)=i-1.
*
            if (typeSFR.le.-1) then
               call bracket(ntimesfile,timefile,time(i),i0)
               call Steffen(ntimesfile,timefile,SFRfile,.false.,
     $              time(i),SFR(i),i0)
               if (typeSFR.eq.-2) then
                  call Steffen(ntimesfile,timefile,Zfile,.false.,
     $                 time(i),ZSFR(i),i0)
               endif
            endif
*
            call weightZ(ZSFR(i),ZSSP,nZ,alpha(i),
     $           iZinf(i),iZsup(i))

*     Instantaneous burst

            if (typeSFR.eq.0) then 
               SFR(i)=0.        
               SFR(1)=1.
            end if

*     Constant star formation rate

            if (typeSFR.eq.1) then
               if (time(i).le.param2) then
                  SFR(i)=param1
               else
                  SFR(i)=0.
               end if
            end if

*     Exponentially decreasing or increasing star formation rate

            if (typeSFR.eq.2) SFR(i)=param2*exp(-time(i)/param1)/param1

*     SFR proportional to a power of the mass of gas

            if (typeSFR.eq.3) SFR(i)=sigmagas(i)**param1/param2

c            if (typeSFR.eq.n>=10) then
c               SFR(i)=your SFR law (note that i = time in Myr + 1)
c            end if

            if (time(i).ge.twind) SFR(i)=0.

*     Maximal SFR

            SFRmax=sigmagas(i)
            if (SFR(i).gt.SFRmax) then
               SFR(i)=SFRmax
               if (excessSFR.eqv..false.) then
                  write(40,'(2a,i5,a)') 'WARNING: the SFR has',
     $                 ' exceeded the maximal possible SFR at ',
     $                 i-1,' Myr!'
                  write(*,*) 'WARNING: ',filespectra,
     $                 ': the SFR has exceeded the',
     $                 ' maximal possible SFR at ',i-1,' Myr!'
                  excessSFR=.true.
               end if
            end if
            sigmagas(i)=sigmagas(i)-SFR(i)
            sigmaZ=sigmaZ-SFR(i)*Zgas(i)
            SFRlum(i)=(1.d0-fsub)*SFR(i)
            dejecta=0.d0
            dejectaZ=0.d0
            dsigmaBHNS=0.d0
            dsigmaWD=0.d0
            agestarsi=0.d0
            sigmastarsi=0.d0
            Zstarsi=0.d0
            do k=1,i
               iZinfk=iZinf(k)
               iZsupk=iZsup(k)
               q=i+1-k
               alphak=alpha(k)
               SFRk=SFRlum(k)
               SFRalphak=SFRk*alphak
               SFR1malphak=SFRk*(1.-alphak)
               dejecta=dejecta+SFRalphak*ejecta(iZinfk,q)
     $              +SFR1malphak*ejecta(iZsupk,q)
               dejectaZ=dejectaZ+SFRalphak*ejectaZ(iZinfk,q)
     $              +SFR1malphak*ejectaZ(iZsupk,q)
     $              +Zgas(k)*(SFRalphak*ejecta(iZinfk,q)
     $              +SFR1malphak*ejecta(iZsupk,q))
               dsigmaBHNS=dsigmaBHNS
     $              +SFRalphak*massBHNS(iZinfk,q)
     $              +SFR1malphak*massBHNS(iZsupk,q)
               dsigmaWD=dsigmaWD
     $              +SFRalphak*massWD(iZinfk,q)
     $              +SFR1malphak*massWD(iZsupk,q)
               sigmastarsi=sigmastarsi+SFRalphak*massalive(iZinfk,q)
     $              +SFR1malphak*massalive(iZsupk,q)
               Zstarsi=Zstarsi+ZSFR(k)
     $              *(SFRalphak*massalive(iZinfk,q)
     $              +SFR1malphak*massalive(iZsupk,q))
               agestarsi=agestarsi+(i-k)
     $              *(SFRalphak*massalive(iZinfk,q)
     $              +SFR1malphak*massalive(iZsupk,q))
            end do
            agestars(i)=agestarsi
            sigmastars(i)=sigmastarsi
            Zstars(i)=Zstarsi
            dsigmagas=dejecta+infall*exp(-time(i)/tinfall)/tinfall
            dsigmaZ=dejectaZ
     $           +infall*exp(-time(i)/tinfall)/tinfall*Zinfall
            if (sigmastars(i).gt.0.) then
               agestars(i)=agestars(i)/sigmastars(i)
               Zstars(i)=Zstars(i)/sigmastars(i)
            else
               agestars(i)=0.
               Zstars(i)=0.
            end if
            sigmasub(i+1)=sigmasub(i)+fsub*SFR(i)
            if (time(i).ge.twind) then
               Mgal(i+1)=Mgal(i)-sigmagas(i)-dejecta
               sigmaBHNS(i+1)=sigmaBHNS(i)+dsigmaBHNS
               sigmaWD(i+1)=sigmaWD(i)+dsigmaWD
               sigmagas(i+1)=0.
               sigmaZ=0.
               Zgas(i+1)=0.               
            else
               Mgal(i+1)=Mgal(i)+infall
     $              *exp(-time(i)/tinfall)/tinfall
               sigmaBHNS(i+1)=sigmaBHNS(i)+dsigmaBHNS
               sigmaWD(i+1)=sigmaWD(i)+dsigmaWD
               sigmagas(i+1)=max(sigmagas(i)+dsigmagas,0.)
               sigmaZ=max(sigmaZ+dsigmaZ,0.)
               if (sigmagas(i+1).le.0.) then 
                  Zgas(i+1)=0.
               else
                  Zgas(i+1)=sigmaZ/sigmagas(i+1)
               end if
            end if
            if (codeZ.eq.1) then
               ZSFR(i+1)=Zgas(i+1)
            else
               ZSFR(i+1)=ZSFR(1)
            end if
         end do

***** Synthetic spectra.

         write(40,'(60a1)') ('*',i=1,60)
         
         write(40,*) ntimesimpr,nlambda,nlines
         write(40,'(5(f9.0,2x))') (lambda(i),i=1,nlambda)
         write(40,'(5(f9.2,2x))') (lambdaline(i),i=1,nlines)
         i0=0
         do jimpr=1,ntimesimpr
            i=timeimpr(jimpr)+1
            if (answerneb.eq.'Nebular emission') then
               if (codeext.ne.0) then
                  taudust=0.5*Zgas(i)/Zsol
                  call bracket(nSpitzer,taudustSpitzer,taudust,i0)
                  call Steffen(nSpitzer,taudustSpitzer,ySpitzer,.false.,
     $                 taudust,y,i0)
                  f=y**3
               else
                  f=1.
               endif
            else
               f=0.
            endif
            do j=1,nlambda
               fluxgalj=0.d0
               ffnebj=f*fneb(j)
               do k=1,i
                  q=i+1-k
                  iZinfk=iZinf(k)
                  iZsupk=iZsup(k)
                  invtimeq=invtime(q)
                  invtimeqp1=invtimeq+1
                  betaq=beta(q)
                  alphak=alpha(k)
                  betaqalphak=betaq*alphak
                  betaq1malphak=betaq-betaqalphak
                  onembetaqalphak=alphak-betaqalphak
                  onembetaq1malphak=1.-alphak-betaq+betaqalphak
                  fluxgalj=fluxgalj+SFRlum(k)
     $                 *(betaqalphak*(fluxSSP(iZinfk,invtimeq,j)
     $                 +NLym(iZinfk,invtimeq)*ffnebj)
     $                 +betaq1malphak*(fluxSSP(iZsupk,invtimeq,j)
     $                 +NLym(iZsupk,invtimeq)*ffnebj)
     $                 +onembetaqalphak*(fluxSSP(iZinfk,invtimeqp1,j)
     $                 +NLym(iZinfk,invtimeqp1)*ffnebj)
     $                 +onembetaq1malphak
     $                 *(fluxSSP(iZsupk,invtimeqp1,j)
     $                 +NLym(iZsupk,invtimeqp1)*ffnebj))
               end do
               fluxgal(j)=fluxgalj
            end do
            do j=1,nlines
               flinetot(j)=0.
               do k=max(i-50,1),i
                  q=i+1-k
                  flinetot(j)=flinetot(j)+SFRlum(k)*(beta(q)
     $                 *(alpha(k)*f*NLym(iZinf(k),invtime(q))*fline(j)
     $                 +(1.-alpha(k))*f
     $                 *NLym(iZsup(k),invtime(q))*fline(j))+(1.-beta(q))
     $                 *(alpha(k)*f*NLym(iZinf(k),invtime(q)+1)*fline(j)
     $                 +(1.-alpha(k))
     $                 *f*NLym(iZsup(k),invtime(q)+1)*fline(j)))
               end do
            end do
            fluxbol=0.         
            fluxbolZ=0.
            fluxbolage=0.
            do k=1,i
               q=i+1-k
               fluxbol=fluxbol+SFRlum(k)*(beta(q)
     $              *(alpha(k)*fluxbolSSP(iZinf(k),invtime(q))
     $              +(1.-alpha(k))*fluxbolSSP(iZsup(k),invtime(q)))
     $              +(1.-beta(q))
     $              *(alpha(k)*fluxbolSSP(iZinf(k),invtime(q)+1)
     $              +(1.-alpha(k))*fluxbolSSP(iZsup(k),invtime(q)+1)))
               fluxbolZ=fluxbolZ+SFRlum(k)*ZSFR(k)*(beta(q)
     $              *(alpha(k)*fluxbolSSP(iZinf(k),invtime(q))
     $              +(1.-alpha(k))*fluxbolSSP(iZsup(k),invtime(q)))
     $              +(1.-beta(q))
     $              *(alpha(k)*fluxbolSSP(iZinf(k),invtime(q)+1)
     $              +(1.-alpha(k))*fluxbolSSP(iZsup(k),invtime(q)+1)))
               fluxbolage=fluxbolage+SFRlum(k)*(i-k)*(beta(q)
     $              *(alpha(k)*fluxbolSSP(iZinf(k),invtime(q))
     $              +(1.-alpha(k))*fluxbolSSP(iZsup(k),invtime(q)))
     $              +(1.-beta(q))
     $              *(alpha(k)*fluxbolSSP(iZinf(k),invtime(q)+1)
     $              +(1.-alpha(k))*fluxbolSSP(iZsup(k),invtime(q)+1)))
            end do
*     Calculation of the number of Lyman continuum photons. 
            NLymtot=0.
            do k=max(1,i-50),i  
               NLymtot=NLymtot
     $              +SFRlum(k)*(alpha(k)*NLym(iZinf(k),i+1-k)
     $              +(1.-alpha(k))*NLym(iZsup(k),i+1-k))
            end do
*     Calculation of the number of SNII.
            nSNIItot=0.
            do k=max(1,i-100),i 
               nSNIItot=nSNIItot
     $              +SFRlum(k)*(alpha(k)*nSNII(iZinf(k),i+1-k)
     $              +(1.-alpha(k))*nSNII(iZsup(k),i+1-k))
            end do
*     Calculation of the number of SNIa.
            nSNIatot=0.
            do k=1,i            
               nSNIatot=nSNIatot
     $              +SFRlum(k)*(alpha(k)*nSNIa(iZinf(k),i+1-k)
     $              +(1.-alpha(k))*nSNIa(iZsup(k),i+1-k))
            end do
            
********The Lyman continuum photons used to produce the nebular emission
********are removed from the spectrum.

            do j=1,i912
               fluxgal(j)=(1.-f)*fluxgal(j)
            end do
            
*******Extinction.
            
            fluxext=0.         
            tauV=0.
            if (codeext.ne.0) then
               call dust_composition(Zext,frac,
     $              Zgas(i),coeffextZ)

***********If there is some extinction, the Lyman continuum photons
***********not absorbed by the gas are absorbed by the dust.         
   
               if (Zgas(i).gt.0.) then
                  do j=1,i912-1
                     fluxext=fluxext+(fluxgal(j)+fluxgal(j+1))
     $                    *(lambda(j+1)-lambda(j))/2.
                     fluxgal(j)=0.
                  end do
                  fluxext=fluxext+(912.-lambda(i912))*fluxgal(i912)
                  fluxgal(i912)=0.
               
***********Idem for Lyman alpha.

                  fluxext=fluxext+flinetot(30)
                  flinetot(30)=0.     
               end if

***********Ellipticals.            
               
               if (codeext.eq.1) then
                  j0=0
                  NH=5.3e22*sigmagas(i)
                  col_dens=0.35*Zgas(i)*NH*mH*1.4
                  call ext(nlambdaext,lambdaext,tauext,
     $                 albedoext,asymext,5500.,coeffextZ,
     $                 tauV,albedo,asym,j0)
                  tauV=tauV*col_dens
                  
**************Continuum.
                  
                  i01=0
                  i02=0
                  i03=0
                  j0=0
                  do j=1,nlambda
                     call ext(nlambdaext,lambdaext,tauext,
     $                    albedoext,asymext,lambda(j),coeffextZ,
     $                    tau,albedo,asym,j0)
                     tau=tau*col_dens
                     call Steffen3D(tauKing,albKing,asymKing,emergKing,
     $                    nK1,nK2,nK3,tau,albedo,asym,corextcont(j),
     $                    i01,i02,i03)
                     corextcont(j)=exp(corextcont(j))
                  end do
                  do j=1,nlambda-1
                     fluxext=fluxext+((1.-corextcont(j))*fluxgal(j)
     $                    +(1.-corextcont(j+1))*fluxgal(j+1))/2.
     $                    *(lambda(j+1)-lambda(j)) 
                  end do
                  do j=1,nlambda
                     fluxgal(j)=fluxgal(j)*corextcont(j)
                  end do
                  
**************Lines.
                  
                  do j=1,nlines
                     i01=0
                     i02=0
                     i03=0
                     j0=0
                     call ext(nlambdaext,lambdaext,tauext,
     $                    albedoext,asymext,lambdaline(j),
     $                    coeffextZ,tau,albedo,asym,j0)
                     tau=tau*col_dens
                     call Steffen3D(tauKing,albKing,asymKing,emergKing,
     $                    nK1,nK2,nK3,tau,albedo,asym,corextline(j),
     $                    i01,i02,i03)
                     corextline(j)=exp(corextline(j))
                     fluxext=fluxext+(1.-corextline(j))*flinetot(j)
                     flinetot(j)=flinetot(j)*corextline(j)
                  end do               
               end if
               
***********Disk galaxies.
               
               if (codeext.ge.2) then
                  NH=6.8e21*sigmagas(i)
                  col_dens=0.35*Zgas(i)*NH*mH*1.4
                  j0=0
                  call ext(nlambdaext,lambdaext,tauext,
     $                 albedoext,asymext,5500.,coeffextZ,
     $                 tauV,albedo,asym,j0)
                  tauV=tauV*col_dens

**************Continuuum.

                  i01=0
                  i02=0
                  i03=0
                  i04=0
                  j0=0
                  do j=1,nlambda
                     call ext(nlambdaext,lambdaext,tauext,
     $                    albedoext,asymext,lambda(j),coeffextZ,
     $                    tau,albedo,asym,j0)
                     tau=tau*col_dens
                     call Steffen3D(tauslab,albslab,asymslab,emergslab,
     $                    nslab1,nslab2,nslab3,tau,albedo,asym,
     $                    corextcont(j),i01,i02,i03)
                     corextcont(j)=exp(corextcont(j))
                     if (codeext.eq.3) then
                        call Steffen4D(tauslab,albslab,asymslab,
     $                       inclinslab,emerginclinslab,nslab1,nslab2,
     $                       nslab3,nslab4,tau,albedo,asym,inclin,
     $                       corextinclincont(j),i01,i02,i03,i04)
                        corextinclincont(j)=exp(corextinclincont(j))
                     else
                        corextinclincont(j)=corextcont(j)
                     endif
                  end do
                  do j=1,nlambda-1
                     fluxext=fluxext+((1.-corextcont(j))*fluxgal(j)
     $                    +(1.-corextcont(j+1))*fluxgal(j+1))/2.
     $                    *(lambda(j+1)-lambda(j))
                  end do
                  do j=1,nlambda
                     fluxgal(j)=fluxgal(j)*corextinclincont(j)
                  end do

**************Lines.
                  
                  do j=1,nlines
                     i01=0
                     i02=0
                     i03=0
                     j0=0
                     call ext(nlambdaext,lambdaext,tauext,
     $                    albedoext,asymext,lambdaline(j),coeffextZ,
     $                    tau,albedo,asym,j0)
                     tau=tau*col_dens
                     call Steffen3D(tauslab,albslab,asymslab,emergslab,
     $                    nslab1,nslab2,nslab3,tau,albedo,asym,
     $                    corextline(j),i01,i02,i03)
                     corextline(j)=exp(corextline(j))
                     if (codeext.eq.3) then
                        call Steffen4D(tauslab,albslab,asymslab,
     $                       inclinslab,emerginclinslab,nslab1,nslab2,
     $                       nslab3,nslab4,tau,albedo,asym,inclin,
     $                       corextinclinline(j),i01,i02,i03,i04)
                        corextinclinline(j)=exp(corextinclinline(j))
                     else
                        corextinclinline(j)=corextline(j)
                     endif
                     fluxext=fluxext+(1.-corextline(j))*flinetot(j)
                     flinetot(j)=flinetot(j)*corextinclinline(j)
                  end do
               end if
            end if
            
            if (fluxbol.gt.0.) then
               fluxext=fluxext/fluxbol
               Zbol=fluxbolZ/fluxbol
               agebol=fluxbolage/fluxbol
            else
               fluxext=0.
               Zbol=0.
               agebol=0.
            end if
            write(40,'(i6,9(1x,e9.3))') i-1,  
     $           Mgal(i),sigmastars(i),
     $           sigmaWD(i),
     $           sigmaBHNS(i),
     $           sigmasub(i),
     $           sigmagas(i),Zgas(i),
     $           Zstars(i),Zbol
            write(40,'(6x,9(1x,e9.3))') fluxbol*Lsol,tauV,fluxext,
     $           SFR(i),NLymtot*Lsol,nSNIItot,nSNIatot,
     $           agestars(i),agebol
            write(40,'(5(e9.4,1x))') (Lsol*fluxgal(j),j=1,nlambda)
            write(40,'(5(e9.4,1x))') (Lsol*flinetot(j),j=1,nlines)
         end do
         close(40)

         read(50,'(a)',iostat=error) a
      end do
      close(50)

      end

***************

      subroutine weightZ(Z,ZSSP,nZ,alpha,iZinf,iZsup)

      implicit none
      integer nmaxZ,nmaxtimes
      parameter(nmaxZ=7)
      parameter(nmaxtimes=20002)
      integer j,nZ
      real Z,ZSSP(nmaxZ)
      real alpha
      integer iZinf,iZsup

      if (nZ.gt.1) then
         if (Z.le.ZSSP(1)) then 
            iZinf=1
            iZsup=1
            alpha=1.
         else
            if (Z.ge.ZSSP(nZ)) then
               iZinf=nZ
               iZsup=nZ
               alpha=1.
            else
               j=1
               do while (((Z-ZSSP(j))*(Z-ZSSP(j+1)).gt.0.)
     $              .and.(j.lt.nZ))
                  j=j+1
               end do
               iZinf=j
               iZsup=j+1
               alpha=log(ZSSP(j+1)/Z)/
     $              log(ZSSP(j+1)/ZSSP(j))
            end if
         end if
      else
         iZinf=1
         iZsup=1
         alpha=1.
      end if

      end

*****

      subroutine Steffen(n,x,y,condition,t,z,i) 
*     Steffen M. (1990), A&A 443, 450                            
      implicit none
      integer n,i
      real x(n),y(n)
      real a,b,c,d,t,z
      real h_i,h_im1,h_ip1
      real s_i,s_im1,s_ip1
      real p_i,p_ip1,y1_i,y1_ip1
      real tmx_i
      logical condition

      if (i.eq.1) then
         h_i=x(i+1)-x(i)
         s_i=(y(i+1)-y(i))/h_i
         h_ip1=x(i+2)-x(i+1)
         s_ip1=(y(i+2)-y(i+1))/h_ip1
         p_ip1=(s_i*h_ip1+s_ip1*h_i)/(h_i+h_ip1)
         y1_ip1=(sign(1.,s_i)+sign(1.,s_ip1))
     $        *min(abs(s_i),abs(s_ip1),0.5*abs(p_ip1))
         if (condition) then    
*     1st derivative is null in x=0
            y1_i=(6.*s_i*x(1)*(x(1)+h_i)
     $           -y1_ip1*x(1)*(3.*x(1)+2.*h_i))
     $           /(3.*x(1)**2+4.*x(1)*h_i+h_i**2)
         else
            y1_i=1.5*s_i-0.5*y1_ip1
         end if
      else
         if (i.eq.n-1) then
            h_i=x(i+1)-x(i)
            s_i=(y(i+1)-y(i))/h_i
            h_im1=x(i)-x(i-1)
            s_im1=(y(i)-y(i-1))/h_im1
            p_i=(s_im1*h_i+s_i*h_im1)/(h_im1+h_i)
            y1_i=(sign(1.,s_im1)+sign(1.,s_i))
     $        *min(abs(s_im1),abs(s_i),0.5*abs(p_i))
            y1_ip1=1.5*s_i-0.5*y1_i
         else
            h_i=x(i+1)-x(i)
            s_i=(y(i+1)-y(i))/h_i
            h_im1=x(i)-x(i-1)
            s_im1=(y(i)-y(i-1))/h_im1
            h_ip1=x(i+2)-x(i+1)
            s_ip1=(y(i+2)-y(i+1))/h_ip1
            p_i=(s_im1*h_i+s_i*h_im1)/(h_im1+h_i)
            p_ip1=(s_i*h_ip1+s_ip1*h_i)/(h_i+h_ip1)
            y1_i=(sign(1.,s_im1)+sign(1.,s_i))
     $        *min(abs(s_im1),abs(s_i),0.5*abs(p_i))
            y1_ip1=(sign(1.,s_i)+sign(1.,s_ip1))
     $        *min(abs(s_i),abs(s_ip1),0.5*abs(p_ip1))
         end if
      end if
      if (t.gt.x(n)) then
         z=y(n)+y1_ip1*(t-x(n))
      else
         tmx_i=t-x(i)
         a=(y1_i+y1_ip1-2.*s_i)/h_i**2
         b=(3.*s_i-2.*y1_i-y1_ip1)/h_i
         c=y1_i
         d=y(i)
         z=((a*tmx_i+b)*tmx_i+c)*tmx_i+d
      endif

      end

********************

      subroutine Steffen3D(x1,x2,x3,y,n1,n2,n3,
     $     t1,t2,t3,z,i01,i02,i03)
      
      implicit none
      integer n1,n2,n3,i1,i2,i3,n,i01,i02,i03
      parameter(n=100)
      real t1,t2,t3,z,x1(n1),x2(n2),x3(n3)
      real y(n1,n2,n3),y1(n),y2(n),y3(n)

      call bracket(n1,x1,t1,i01)
      call bracket(n2,x2,t2,i02)
      call bracket(n3,x3,t3,i03)
      do i1=max(1,i01-1),min(n1,i01+2)
         do i2=max(1,i02-1),min(n2,i02+2)
            do i3=max(1,i03-1),min(n3,i03+2)
              y3(i3)=y(i1,i2,i3)
           end do
           call Steffen(n3,x3,y3,.false.,t3,y2(i2),i03)
        end do
        call Steffen(n2,x2,y2,.false.,t2,y1(i1),i02)
      end do
      call Steffen(n1,x1,y1,.false.,t1,z,i01)

      end

********************

      subroutine Steffen4D(x1,x2,x3,x4,y,n1,n2,n3,n4,
     $     t1,t2,t3,t4,z,i01,i02,i03,i04)
      
      implicit none
      integer n1,n2,n3,n4,i1,i2,i3,i4,n,i01,i02,i03,i04
      parameter(n=100)
      real t1,t2,t3,t4,z,x1(n1),x2(n2),x3(n3),x4(n4)
      real y(n1,n2,n3,n4),y1(n),y2(n),y3(n),y4(n)

      call bracket(n1,x1,t1,i01)
      call bracket(n2,x2,t2,i02)
      call bracket(n3,x3,t3,i03)
      call bracket(n4,x4,t4,i04)
      do i1=max(1,i01-1),min(n1,i01+2)
         do i2=max(1,i02-1),min(n2,i02+2)
            do i3=max(1,i03-1),min(n3,i03+2)
               do i4=max(1,i04-1),min(n4,i04+2)
                 y4(i4)=y(i1,i2,i3,i4)
              end do              
              call Steffen(n4,x4,y4,.true.,t4,y3(i3),i04)
           end do
           call Steffen(n3,x3,y3,.false.,t3,y2(i2),i03)
        end do
        call Steffen(n2,x2,y2,.false.,t2,y1(i1),i02)
      end do
      call Steffen(n1,x1,y1,.false.,t1,z,i01)

      end

********

      real function interplinlog(x1,x2,y1,y2,x)

      implicit none
      
      real x1,x2,y1,y2,x,epsilon

      epsilon=1.e-37
      if (abs(x2-x1).lt.epsilon) then
         interplinlog=y1
      else
         if (y1.le.epsilon.or.y2.le.epsilon) then
            interplinlog=0.
         else
            interplinlog=y1*(y2/y1)**((x-x1)/(x2-x1))
         endif
      endif

      end

******

      subroutine dust_composition(Zext,frac,Z,coeffextZ)

      implicit none
      real Zext(5),frac(5),Z,coeffextZ(2)
      integer i

      i=0
      call bracket(5,Zext,Z,i)
      call Steffen(5,Zext,frac,.false.,Z,coeffextZ(1),i)
      coeffextZ(2)=1.-coeffextZ(1)

      end

***********

      subroutine ext(nlambdaext,lambdaext,tauext,
     $     albedoext,asymext,lambda,coeffextZ,tau,albedo,asym,j)

      implicit none
      integer nmaxlambdaext
      parameter(nmaxlambdaext=300)     
      integer nlambdaext,j
      real coeffextZ(2)
      real lambdaext(nmaxlambdaext),tauext(nmaxlambdaext,2)
      real albedoext(nmaxlambdaext,2),asymext(nmaxlambdaext,2)
      real lambda,tau,albedo,asym
      real tau1,tau2,albedo1,albedo2,asym1,asym2,alpha

      call bracket(nlambdaext,lambdaext,lambda,j)
      alpha=(1./lambda-1./lambdaext(j))/
     $     (1./lambdaext(j+1)-1./lambdaext(j))

      tau1=tauext(j,1)+alpha*(tauext(j+1,1)-tauext(j,1))
      tau2=tauext(j,2)+alpha*(tauext(j+1,2)-tauext(j,2))
      tau=coeffextZ(1)*tau1+coeffextZ(2)*tau2
      if (tau.gt.0.) then
         albedo1=albedoext(j,1)+alpha*(albedoext(j+1,1)-albedoext(j,1))
         albedo2=albedoext(j,2)+alpha*(albedoext(j+1,2)-albedoext(j,2))
         albedo=(coeffextZ(1)*tau1*albedo1
     $        +coeffextZ(2)*tau2*albedo2)/tau
         if (albedo.gt.1.) albedo=1.
         if (albedo.gt.0.) then
            asym1=asymext(j,1)+alpha*(asymext(j+1,1)-asymext(j,1))
            asym2=asymext(j,2)+alpha*(asymext(j+1,2)-asymext(j,2))
            asym=(coeffextZ(1)*tau1*albedo1*asym1
     $           +coeffextZ(2)*tau2*albedo2*asym2)/tau/albedo
            if (asym.gt.1.) asym=1.
            if (asym.lt.-1.) asym=-1.
         else
            albedo=0.
            asym=0.
         endif
      else
         tau=0.
         albedo=0.
         asym=0.
      endif

      end

******

      subroutine read_write_real(ur,uw,x)

      implicit none
      character*72 a
      integer ur,uw
      real x

      read(ur,'(a)') a
      write(uw,'(a)') a
      a=a(index(a,':')+1:)
      read(a,*) x

      end

******

      subroutine read_write_char(ur,uz,a)
      
      implicit none
      character*72 a
      integer ur,uz,i

      read(ur,'(a)') a
      write(uz,'(a)') a
      i=index(a,':')
      i=i+1
      do while (a(i:i).eq.' ') 
         i=i+1
      end do
      a=a(i:)

      end

******

      subroutine read_write_integer(ur,uw,i)

      implicit none
      character*72 a
      integer ur,uw,i

      read(ur,'(a)') a
      write(uw,'(a)') a
      a=a(index(a,':')+1:)
      read(a,*) i

      end

***** 

      subroutine bracket(n,x,t,i)

      implicit none
      integer n,iinf,isup,imed,i,di,niter
      real x(n),t,xinf,xsup
      
      if (i.eq.0) then
         niter=0
         iinf=1
         isup=n
         if (t.gt.x(iinf)) then
            if (t.ge.x(isup)) then
               iinf=isup-1
            else
               do while (iinf+1.lt.isup)
                  niter=niter+1
                  imed=(iinf+isup)/2
                  if (t.le.x(imed)) then
                     isup=imed
                  else
                     iinf=imed
                  end if
               end do
            end if
         end if         
      else
         niter=0
         iinf=1
         isup=n
         if (t.gt.x(iinf)) then
            if (t.ge.x(isup)) then
               iinf=isup-1
            else
               di=1
               if (t.ge.x(i)) then
                  iinf=i
                  isup=min(n,i+di)
                  xsup=x(isup)
                  do while(t.gt.xsup) 
                     niter=niter+1
                     di=2*di
                     isup=min(n,i+di)
                     xsup=x(isup)
                  enddo
               else
                  isup=i
                  iinf=max(1,i-di)
                  xinf=x(iinf)
                  do while(t.lt.xinf) 
                     niter=niter+1
                     di=2*di
                     iinf=max(1,i-di)
                     xinf=x(iinf)
                  enddo
               endif
               do while (iinf+1.lt.isup)
                  niter=niter+1
                  imed=(iinf+isup)/2
                  if (t.le.x(imed)) then
                     isup=imed
                  else
                     iinf=imed
                  end if
               end do
            end if
         end if
      end if
      i=iinf

      end
