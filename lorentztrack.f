c************************************************************************
c Track a particle trajectory obeying the Lorentz force, from an initial 
c position and energy.
      subroutine lorentztrack(nt,xt,vt,tt,dt,i,icase,nmax)
      implicit none
      integer nt,i,icase,nmax
      real xt(nt,3),vt(nt,3),tt(nt),dt
c xt,vt is the phase-space position (x1,x2,x3), (v1,v2,v3) at time tt.
c ith element contains initial conditions on entry. dt is time step.
c nt is dim of arrays. 
c On entry i is initiat step number.
c          nmax is max step to integrate to.
      integer nd,j,ierr,k
      real Etest(3),Btest(3)
      real y1(6),y2(6),time
      real dth,cond,thecond
      real small
      data small/1.e-4/
      data nd/6/
      data thecond/.005/

c  Set up starting point. From input data.
      write(*,'(a,i2,a,f8.4)') 'icase=',icase,'  dt=',dt
      write(*,1002)'x0=(',(xt(i,k),k=1,3),')   v0=(',(vt(i,k),k=1,3),')'
 1002 format(a,3f8.3,a,3f8.3,a)
      if(nmax.gt.nt)nmax=nt
      do k=1,3
         y2(k)=xt(i,k)
         y2(k+3)=vt(i,k)
      enddo
      time=tt(i)
      dth=dt
      do i=1,nmax
c     Iterate till we leave the mesh or run out of track length.
c         write(*,*)'y2=',y2(1),y2(2),y2(3),y2(4)
         do j=1,6
            y1(j)=y2(j)
         enddo
         do k=1,3
            xt(i,k)=y1(k)
            vt(i,k)=y1(k+3)
         enddo
c         write(*,'(a,6f10.4)')'y1:',y1
         tt(i)=time
 1003    continue
c decide how big dt must be if dynamically sized step.
         call RKADVC(dth,tt(i),ND,Y1,icase,Y2,cond,IERR)
         if(ierr.ne.0)goto 9
         if(thecond.ne.0)then
c Step length adaptation:
c Test whether this step was reasonable, if not, adjust step.
            if(cond.gt.thecond .and. dth.gt.dt/1024.)then
               dth=dth/2.
            elseif(dth.lt.dt .and. cond.lt.thecond/5.)then
               dth=dth*2.
            else
               goto 1004
            endif
            write(*,'(a,f10.7,a,i7,a,f10.6)')'Adjusting dt',dth,' step'
     $           ,i,' condition',cond
            goto 1003
         endif
 1004    continue
         time=tt(i)+dt
      enddo
      i=nmax-1
 9    continue
      end
c****************************************************************************
      subroutine RKFUN(yp,y,t,icase,ierr)
      real yp(6),y(6),t
      integer icase,ierr
C Lorentz RK advance function. Evaluates derivative of y into yp.
c x=y1-3, v=y4-6. 
c dv/dt = q/m(  E + vxB).   dx/dt = v
      real E(3), B(3)
      real qm
      parameter (qm=1.)
C Evaluate E and B, q/m.
      call getfields(y,t,icase,E,B)
C Enter derivatives.
      yp(1)=y(4)
      yp(2)=y(5)
      yp(3)=y(6)
      yp(4)=qm*(E(1) + y(5)*B(3) - y(6)*B(2))
      yp(5)=qm*(E(2) + y(6)*B(1) - y(4)*B(3))
      yp(6)=qm*(E(3) + y(4)*B(2) - y(5)*B(1))
c      write(*,*)yp
      ierr=0
      end
c****************************************************************************
      subroutine getfields(y,t,icase,E,B)
      real y(6),t,E(3),B(3)
      integer icase
      common /lorentzcom/par1,par2
      if(icase.eq.0)then
         E(1)=0.
         E(2)=0.
         E(3)=0.
         B(1)=0.
         B(2)=0.
         B(3)=1.
      elseif(icase.eq.1)then
c Uniform B and E
         E(1)=0.
         E(2)=0.05
         E(3)=0.
         B(1)=0.
         B(2)=0.
         B(3)=1.
      elseif(icase.eq.2)then
c Grad B
         E(1)=0.
         E(2)=0.
         E(3)=0.
         B(1)=0.
         B(2)=0.
         B(3)=1.+.1*y(1)
      elseif(icase.eq.3)then
c Curvature: B = \hat{e_\phi}
         E(1)=0.
         E(2)=0.
         E(3)=0.
         yr=sqrt(y(1)**2+ y(2)**2)
         if(yr .ne.0.) then
            B(1)=-y(2)/yr
            B(2)=y(1)/yr
         else
            B(1)=0.
            B(2)=0.
         endif
         B(3)=0.
      elseif(icase.eq.4)then
c Vacuum toroidal field
         E(1)=0.
         E(2)=0.
         E(3)=0.
         yr=(y(1)**2+ y(2)**2)*.2
         if(yr .ne.0.) then
            B(1)=-y(2)/yr
            B(2)=y(1)/yr
         else
            B(1)=0.
            B(2)=0.
         endif
         B(3)=0.
      elseif(icase.eq.5)then
c Dipole field 3(p.r) r /r^5 - p/r^3).  p in z direction. 
         E(1)=0.
         E(2)=0.
         E(3)=0.
         r=sqrt(y(1)**2+ y(2)**2 + y(3)**2)
         r5=r**5
         p=100.
         if(r .ne.0.) then
            B(1)=3.*y(3)*p*y(1)/r5
            B(2)=3.*y(3)*p*y(2)/r5
            B(3)=(3.*y(3)*p*y(3) - p*r*r) /r5
         else
            B(1)=0.
            B(2)=0.
            B(3)=0.
         endif
      elseif(icase.eq.6)then
c Tokamak with Ro=5, uniform current (B_theta\propto r)
c Bt=1. B_theta(1)=.080 q=rBt/RB_theta=2.5
         Ro=5.
         Bth0=0.080
         E(1)=0.
         E(2)=0.
         E(3)=0.
         yr=(y(1)**2+ y(2)**2)*.2
         if(yr .ne.0.) then
            B(1)=-y(2)/yr
            B(2)=y(1)/yr
            B(1)=B(1)-Bth0*y(3)*y(1)/yr
            B(2)=B(2)-Bth0*y(3)*y(2)/yr
            B(3)=Bth0*(yr-Ro)
         else
            B(1)=0.
            B(2)=0.
            B(3)=0.
         endif
      elseif(icase.eq.7)then
c Combination of point particle (Coulomb) and uniform E-field
         yr=y(1)**2+y(2)**2+y(3)**2+1.e-20
         do i=1,3
            B(i)=0.
            E(i)=-y(i)/yr**1.5
            if(i.eq.3)E(i)=E(i)+par1
         enddo
      endif
      end
c****************************************************************************
c Return the description.
c And draw any arrows needed.
      subroutine descfields(icase,desc)
      integer icase
      character*(*) desc
      real arrow(3,2)
      if(icase.eq.0)then
         desc='Uniform Bz'
c Draw arrow denoting Bz
         call nxyz2wxyz(-.3,-.3,-.2,arrow(1,1),arrow(2,1),arrow(3,1))
         call nxyz2wxyz(-.3,-.3,.2,arrow(1,2),arrow(2,2),arrow(3,2))
c         call vec3w(arrow(1,1),arrow(2,1),arrow(3,1),0)
c         call vec3w(arrow(1,2),arrow(2,2),arrow(3,2),1)
         bbl=(arrow(3,2)-arrow(3,1))*.15
c         write(*,*)'bbl=',bbl
         call accisgradinit(-40000,-40000,22000,64000,64000,64000)
         call arrow3path(arrow,2,8,0,.1*bbl,bbl,.3*bbl,1)
      elseif(icase.eq.1)then
c Draw arrow denoting Bz
         call nxyz2wxyz(-.3,-.3,-.2,arrow(1,1),arrow(2,1),arrow(3,1))
         call nxyz2wxyz(-.3,-.3,.2,arrow(1,2),arrow(2,2),arrow(3,2))
         bbl=(arrow(3,2)-arrow(3,1))*.15
         call accisgradinit(-40000,-40000,22000,64000,64000,64000)
         call arrow3path(arrow,2,6,0,.1*bbl,bbl,.3*bbl,1)
         call drcstr('  B!dz!d')
c Draw arrow denoting Ey
         call nxyz2wxyz(-.3,-.2,-.2,arrow(1,1),arrow(2,1),arrow(3,1))
         call nxyz2wxyz(-.3,.2,-.2,arrow(1,2),arrow(2,2),arrow(3,2))
         bbl=(arrow(2,2)-arrow(2,1))*.15
         call accisgradinit(22000,-40000,-40000,64000,64000,64000)
         call arrow3path(arrow,2,8,0,.1*bbl,bbl,.3*bbl,1)
         call drcstr('  E!dy!d')
         desc='Uniform Bz, Ey' 
      elseif(icase.eq.2)then
c Draw arrows denoting Bz
         zh=.1
         call accisgradinit(-40000,-40000,22000,64000,64000,64000)
         call nxyz2wxyz(-.2,-.3,-zh,arrow(1,1),arrow(2,1),arrow(3,1))
         x1=arrow(1,1)
         z1=-zh*(1.+0.1*x1)
         z2=zh*(1.+0.1*x1)
         call nxyz2wxyz(-.2,-.3,z1,arrow(1,1),arrow(2,1),arrow(3,1))
         call nxyz2wxyz(-.2,-.3,z2,arrow(1,2),arrow(2,2),arrow(3,2))
         bbl=(arrow(3,2)-arrow(3,1))*.15
c         bbl=(z2-z1)*.15
c         write(*,*)'x1,z1,z2',x1,z1,z2,bbl
         call arrow3path(arrow,2,8,0,.1*bbl,bbl,.3*bbl,1)
         call nxyz2wxyz(.2,-.3,-zh,arrow(1,1),arrow(2,1),arrow(3,1))
         x2=arrow(1,1)
         z1=-zh*(1.+0.1*x2)
         z2=zh*(1.+0.1*x2)
         call nxyz2wxyz(.2,-.3,z1,arrow(1,1),arrow(2,1),arrow(3,1))
         call nxyz2wxyz(.2,-.3,z2,arrow(1,2),arrow(2,2),arrow(3,2))
         bbl=(arrow(3,2)-arrow(3,1))*.15
         call arrow3path(arrow,2,8,0,.1*bbl,bbl,.3*bbl,1)
         desc='Gradient: Bz = 1 + 0.1 x' 
      elseif(icase.eq.3)then
         desc='Curvature: Constant |B|, azimuthal.'   
      elseif(icase.eq.4)then
         desc='Vacuum azimuthal field (B!A&!@1/r)'   
      elseif(icase.eq.5)then
         desc='Dipole B field'
      elseif(icase.eq.6)then
         desc='Tokamak (q=2.5)'
      elseif(icase.eq.7)then
         desc='Point charge in E-field'
      endif
      end
c****************************************************************************
        SUBROUTINE RKADVC(DX,X0,ND,Y1,icase,Y2,cond,IERR)
C 
C A fortran routine for ADVANCING an initial-value ODE integration
C by a fourth-order Runge-Kutta step. 
c Based upon NR but (1) a case switch, and (2) reporting conditioning.
C
C The equation to be solved is to be given in the form
C               dY/dX = F(Y,X)
C where Y is a vector of dimension equal to the order of the system.
C       INPUTS:
C       DX      The point spacing in the independent"var.
C       X0      The initial point.
C       ND      The order of the system. Must be <= 20
C       Y1      The initial vector, dimension ND (<=20)
c       icase   The case to be used in the function.
c  Outputs
C       Y2      The advanced vector, dimension ND.
c       COND    The condition of the advance defined as 
c                The magnitude difference in y-change 
c                divided by magnitude y-change, all squared.
c       ierr
C
C       RKFUN(YP,Y,X,icase,IERR)  A user-supplied subroutine which evaluates
C               the right hand side of the equation into vector YP.
C                The IERR is a parameter that may be set by FUN on a
C               condition that requires the integration to terminate.
C
C
C
        REAL YW,YK1,YK2,YK3,YK4,Y1,Y2
        INTEGER  ND,icase
        DIMENSION YW(20),YK1(20),YK2(20),YK3(20),YK4(20),Y1(*),Y2(*)
        real yc(20),ycmag2,ydcmag2
C
C                       TYPE INITIAL VALUES.
C       WRITE(*,50)X0,(Y1(J),J=1,ND)
        DX2=DX/2.
C
C                       ADVANCE THE INTEGRATION.
        XHALF=X0+DX2
        XPLUS=X0+DX
        CALL RKFUN(YK1,Y1,X0,icase,IERR)
        ycmag2=0.
        do J=1,ND
           YK1(J)=YK1(J)*DX
           yc(j)=yk1(j)
           ycmag2=ycmag2+yc(j)**2
           YW(J)=Y1(J)+YK1(J)/2.
        enddo
        CALL RKFUN(YK2,YW,XHALF,icase,IERR)
        DO J=1,ND
           YK2(J)=YK2(J)*DX
           YW(J)=Y1(J)+YK2(J)/2.
        enddo
        CALL RKFUN(YK3,YW,XHALF,icase,IERR)
        DO J=1,ND
           YK3(J)=YK3(J)*DX
           YW(J)=Y1(J)+YK3(J)
        enddo
        CALL RKFUN(YK4,YW,XPLUS,icase,IERR)
        ydcmag2=0.
        DO J=1,ND
           YK4(J)=YK4(J)*DX
           ydcmag2=ydcmag2+(yk4(j)-yc(j))**2
           Y2(J)=Y1(J)+(YK1(J)+2.*YK2(J)+2.*YK3(J)+YK4(J))/6.
        enddo
        cond=ydcmag2/ycmag2
C
C       WRITE(*,50)X0,(Y(J),J=1,ND)
50      FORMAT(8G10.3)
C
        RETURN
        END
c***********************************************************************
      subroutine polproj(nt,xt,rp,zp,nmax)
      integer nt,nmax
      real xt(nt,3)
      real rp(nt),zp(nt)
      do i=1,nmax
         rp(i)=sqrt(xt(i,1)**2+xt(i,2)**2)
         zp(i)=xt(i,3)
      enddo
      end
c***********************************************************************
c Main
c***********************************************************************
      program mainlorentz
      integer nt,i,icase,us
      parameter (nt=20000)
      real xt(nt,3),vt(nt,3),tt(nt),dt
      real rp(nt),zp(nt)
      integer ipc
      parameter (ipc=4)
      integer ipgf
      parameter (ipgf=50)
      character*50 filename
      character*32 firstline
      character*100 desc
      logical linset
      common /lorentzcom/par1,par2

      ips=0
      par1=1
      par2=1
c Input.
      linset=.false.
      call getarg(1,filename)
      if(filename(1:1) .eq. ' ') then
         filename='case.dat'
      else
         write(*,*)' Case file:',filename
      endif
c Read case, position, velocity, and dt.
      open(unit=4,file=filename,status='old',err=99)
      read(4,'(a)')firstline
      read(4,*,err=97,end=96)
     $     icase,xt(1,1),xt(1,2),xt(1,3),
     $     vt(1,1),vt(1,2),vt(1,3),dt,nmax,tanimin
      tanim=tanimin
      read(4,*,end=98,err=98)par1,par2
      close(4)
      goto 98
 96   write(*,*) 'File end encountered.'
 97   write(*,*) '**** Error in input file ',filename
      write(*,*)'After text first line, format must be I,3F,3F,F,I,F'
      stop ''
 99   continue
c Defaults.
      write(*,*) 'No file ',filename,'Using default initial conditions.'
      icase=1
      xt(1,1)=1.
      xt(1,2)=0.
      xt(1,3)=0.0001
      vt(1,1)=1.
      vt(1,2)=0.
      vt(1,3)=0.1
      dt=.1
      nmax=500
      tanim=0.3
 98   continue
      tt(1)=0.
      i=1
      call lorentztrack(nt,xt,vt,tt,dt,i,icase,nmax)
c      write(*,'(3f10.4)')((xt(i,j),j=1,3),i=1,100)
c Projected insert code
      call polproj(nt,xt,rp,zp,nmax)
      call minmax(rp,nmax,rmin,rmax)
      call minmax(zp,nmax,zmin,zmax)
      if(icase.gt.2 .and. icase.lt.7)then
         linset=.true.
         rscale=.2
         roff=.82
         zoff=.5
         scale=rscale/max((rmax-rmin),(zmax-zmin))
         write(*,*)'scale=',scale
         do i=1,nmax
            rp(i)=roff+scale*(rp(i)-rmin)
            zp(i)=zoff+scale*(zp(i)-zmin)
c     write(*,*)'rp,zp=',rp(i),zp(i)
         enddo
      endif
 101  format(4f10.4)
c plotting.
      write(*,*)'Plot control: a= toggle animation, r,e= rotate'
     $     ,' p= toggle print'
c      call pfset(3)
c Make the cube really a cube (default is .25,.25,.2,.5,.4).
      call setcube(.2,.2,.2,.5,.4)
 52   call auto3init(xt(1,1),xt(1,2),xt(1,3),nmax)
 51   continue
      call axident3()
c      call axproj(igetcorner())
      call descfields(icase,desc)
      call jdrwstr(.03,.03,desc,1.)
      if(linset)then
         call color(ipc)
         call jdrwstr(roff,zoff+rscale+.01,
     $        'Poloidal Projection',0.5)
         call jdrwstr(roff,zoff-.02,
     $        '!A_!BR!@=!A)!@(!Bx!@!u2!u+!By!@!u2!u)',1.)
         call jdrwstr(roff,zoff,'!Bz!A}!@',-.75)
         call color(15)
      endif
      call accisflush()
      if(tanim.ne.0)then
c Seems you can't sleep a short time.
c Shorter than about 1.e4 does not work with usleep. 
         us=(1000000.*tanim)/nmax
         istep=max(1,10000/us)
         us=us*istep
c         write(*,'(a,i6,a,i6,a,i4,a,f6.2)')'nmax=',nmax,'  us=',us,
c     $        ' istep=',istep,' tanim=',tanim
         if(linset)then
c            call color(ipc)
            call vecn(roff-scale*rmin,zoff,0)
            call vecn(roff-scale*rmin,zoff+scale*(zmax-zmin),1)
            call vec3w(0.,0.,zmin,0)
            call vec3w(0.,0.,zmax,1)
c            call color(15)
         endif
         call vec3w(xt(1,1),xt(1,2),xt(1,3),0)
         iinc=nmax/ipgf
c There's a problem with the vecglx driver. It can't flush properly.
         do i=2,nmax
            if(linset) then
               call color(ipc)
               call vecn(rp(i-1),zp(i-1),0)
               call vecn(rp(i),zp(i),1)
               call color(15)
            endif
            call vec3w(xt(i-1,1),xt(i-1,2),xt(i-1,3),0)
            call vec3w(xt(i,1),xt(i,2),xt(i,3),1)
            if(mod(i,istep).eq.0)then
c               write(*,*)'Calling second accisflush',i,us
               call accisflush()
c asleep puts spurious draws into the ps file (deliberately) 
c that makes it slow and ponderous for command line use.
               if(ips.eq.3.and.mod(i,iinc).eq.0) call asleep(us)
               if(mod(i,istep).eq.0)call usleep(500+us)
            endif
c        call poly3line(xt(1,1),xt(1,2),xt(1,3),i)
         enddo
      else
         call poly3line(xt(1,1),xt(1,2),xt(1,3),nmax)
      endif
      call accisflush()
c Limit to 30frames/s and ensure plotted.
c      write(*,*)'Calling eye3d'
      call usleep(33000)
      call eye3d(ival)
      call rotatezoom(ival)
      if(ival.eq.ichar('a'))then
         if(tanim.ne.0)then
            tanim=0
         else
            tanim=tanimin
         endif
      endif
      if(ival.eq.ichar('p'))then
         if(ips.eq.0) then
            call pfset(3)
            ips=3
            goto 52
         else
            call pfset(0)
            ips=0
         endif
      endif
      if(ival.ne.0 .and. ival.ne.ichar('q'))then
c Using accisclear concatenates all the drawing into one PS file.
         call accisclear()
c accisinit works with the cubeproj calls, producing lots of files
c         call accisinit()
         call cubeproj(icorner)
         call axproj(igetcorner())
         goto 51
      endif
c End plotting tidily without waiting.
      call prtend(' ')
      end
