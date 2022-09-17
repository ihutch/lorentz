c************************************************************************
c Track a particle trajectory obeying the Lorentz force, from an initial 
c position and energy.
      subroutine lorentztrack(nt,xt,vt,tt,dt,i)
      implicit none
      integer nt,i
      real xt(nt,3),vt(nt,3),tt(nt),dt
c xt,vt is the phase-space position (x1,x2,x3), (v1,v2,v3) at time tt.
c first element contains initial conditions on entry. dt is time step.
c nt is maximum step number (dim of arrays). i is working step number.

      integer nd,j,icase,ierr,k

      real y1(6),y2(6),scl,spm
      real small,ymstep2
      data small/1.e-4/
      data nd/4/

c  Set up starting point. From input data.
      do k=1,3
         y2(k)=xt(i,k)
         y2(k+3)=vt(i,k)
      enddo
      do i=1,nt-1
c     Iterate till we leave the mesh or run out of track length.
c         write(*,*)'y2=',y2(1),y2(2),y2(3),y2(4)
         do j=1,6
            y1(j)=y2(j)
         enddo
         do k=1,3
            xt(i,k)=y1(k)
            vt(i,k)=y1(k+3)
         enddo
c decide how big dt must be if dynamically sized step.
         call RKADVC(dt,tt(i),ND,Y1,Y2,icase,IERR)
c If we have exited the mesh.
c         if(y2(1).ge.r(maxi) .or. y2(1).le.r(1) .or.
c     $        y2(2).ge.z(maxj) .or. y2(2).le.z(1)) goto 9
c If a field-line has stopped advancing.
c         if(lfline .and.
c     $        ((y2(1)-y1(1))**2+(y2(2)-y1(2))**2).lt.ymstep2)goto 9
c If an iteration has terminated, probably we have left the mesh.
         if(ierr.ne.0)goto 9
         tt(i+1)=tt(i)+dt
      enddo
      i=nt-1
c      goto 9
c 8    write(*,*)'I, ierr returned..', i, ierr
 9    continue
c      write(*,*)'Leaving track:',i,' points.',' End:',y2(1),y2(2)
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
      yp(4)=qm*(E(3) + y(4)*B(2) - y(5)*B(1))
      end
c****************************************************************************
      subroutine getfields(y,t,icase,E,B)
      real y(6),t,E(3),B(3)
      integer icase
      if(icase.eq.1)then
c Uniform magnetic field in 3-direction.
         E(1)=0.
         E(2)=0.
         E(3)=0.
         B(1)=0.
         B(2)=0.
         B(3)=1.
      elseif(icase.eq.2)then
      endif
      end
c****************************************************************************
	SUBROUTINE RKADVC(DX,X0,ND,Y1,Y2,icase,IERR)
C 
C A fortran routine for ADVANCING an initial-value ODE integration
C by a fourth-order Runge-Kutta step. With a case switch.
C
C The equation to be solved is to be given in the form
C		dY/dX = F(Y,X)
C where Y is a vector of dimension equal to the order of the system.
C 	INPUTS:
C 	DX	The point spacing in the independent"var.
C 	X0	The initial point.
C 	ND	The order of the system. Must be <= 20
C	Y2	The advanced vector, dimension ND.
C	Y1	The initial vector, dimension ND (<=20)
c       icase   The case to be used in the function.
C
C	RKFUN(YP,Y,X,icase,IERR)  A user-supplied subroutine which evaluates
C		the right hand side of the equation into vector YP.
C		 The IERR is a parameter that may be set by FUN on a
C		condition that requires the integration to terminate.
C
C
C
	REAL YW,YK1,YK2,YK3,YK4,Y1,Y2
	INTEGER  ND,icase
	DIMENSION YW(20),YK1(20),YK2(20),YK3(20),YK4(20),Y1(20),Y2(20)
C
C			TYPE INITIAL VALUES.
C	WRITE(*,50)X0,(Y1(J),J=1,ND)
	DX2=DX/2.
C
C			ADVANCE THE INTEGRATION.
	XHALF=X0+DX2
	XPLUS=X0+DX
	CALL RKFUN(YK1,Y1,X0,icase,IERR)
	DO 210 J=1,ND
	YK1(J)=YK1(J)*DX
210	YW(J)=Y1(J)+YK1(J)/2.
	CALL RKFUN(YK2,YW,XHALF,icase,IERR)
	DO 220 J=1,ND
	YK2(J)=YK2(J)*DX
220	YW(J)=Y1(J)+YK2(J)/2.
	CALL RKFUN(YK3,YW,XHALF,icase,IERR)
	DO 230 J=1,ND
	YK3(J)=YK3(J)*DX
230	YW(J)=Y1(J)+YK3(J)
	CALL RKFUN(YK4,YW,XPLUS,icase,IERR)
	DO 240 J=1,ND
	YK4(J)=YK4(J)*DX
240	Y2(J)=Y1(J)+(YK1(J)+2.*YK2(J)+2.*YK3(J)+YK4(J))/6.
C
C	WRITE(*,50)X0,(Y(J),J=1,ND)
50	FORMAT(8G10.3)
C
	RETURN
	END
c***********************************************************************
c Test main only
      program testlorentz
      integer nt,i
      parameter (nt=100)
      real xt(nt,3),vt(nt,3),tt(nt),dt
      xt(1,1)=0.
      xt(1,2)=0.
      xt(1,3)=0.
      vt(1,1)=1.
      vt(1,1)=0.
      vt(1,1)=0.
      dt=.1
      tt(1)=0.
      call lorentztrack(nt,xt,vt,tt,dt,i)
      write(*,101)((xt(i,j),j=1,3),i=1,nt)
 101  format(3f10.4)
      end
c************************************************************************
