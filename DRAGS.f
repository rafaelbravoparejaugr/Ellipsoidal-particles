
c     rem ******************** 14000
c     rem comput drag froces
c     rem ******************** called by main

      subroutine drags (ji, rhof, rug, nu, us)


      implicit     none

      integer      i,nsteps

      real*8      a,al,alpha,ap,at,aux
      real*8      b,beta1t,beta2
      real*8      c,cd,cd0,cd90
      real*8      d1t,d2,delta1,delta2,delta3,deq
      real*8      eps
      real*8      fd,fi
      real*8      dv,dy 
      real*8      ji
      real*8      p,pi 
      real*8      nu
      real*8      Rey,rho2,rhof,rug 
      real*8      s,seq,sep,sph,sphc,sphl 
      real*8      u,us
      real*8      vol
      real*8      y1t,yc1,yc2,yc3
      real*8      x1t,xc1,xc2,xc3
      real*8      z,z0

      include     'param.h'
      include     'zd.h'
      include     'zf.h'
      include     'zrd.h'
      include     'zxyt23.h'

c     constant pi

      pi = acos(-1.d0)

c     assign a,b,c
      a = d(3,3)
      b = d(4,3)
      c = b

      xc1 = d(1,1) 
      yc1 = d(2,1)

      xc2 = d(1,2)
      yc2 = d(2,2)

      xc3 = d(1,3)
      yc3 = d(2,3)

      alpha = d(5,3)

c      xt23 =
c      yt23 = 

c     volume, surface, longitudinal and tranversal projected areas
      vol = 4.0d0/3.0d0*pi*a*b*c
      p   = 1.6075d0
      s   = 4.0d0*pi*(((a**p)*(b**p) + (a**p)*(c**p) 
     &+ (b**p)*(c**p))/3)**(1.0d0/p)
      al  = pi*a*b
      at  = pi*b*c

c     diametro y secciones equivalentes
      deq  = (6.0d0*vol/pi)**(1.0d0/3.0d0)
      seq  = pi*deq**2
      sep  = pi/4.0d0*deq**2
      sph  = seq/s
      sphc = sep/at
      sphl = sep/(s/2.0d0 - al)

c     distancia vertical proyectada por la elipse

      fi   = atan(b/a*tan(pi/2.0d0-alpha))
      dy   = abs(a*cos(fi)*sin(alpha)+b*sin(fi)*cos(alpha))
      dv   = abs(2*dy)

      write(*,*) fi, dv

c     velocidad media

      beta1t = atan(b*sin(fi),a*cos(fi))

      x1t =  a*cos(fi)*cos(alpha) - b*sin(fi)*sin(alpha)
      y1t =  a*cos(fi)*sin(alpha) + b*sin(fi)*cos(alpha)
      d1t =  sqrt(x1t**2 + y1t**2)



      beta2 = atan(yt23-yc3,xt23-xc3)

      rho2 =  atan(yc2-yt23,xc2-xt23)
 
      d2   = sqrt((xt23-xc3)**2+(yt23-yc3)**2)

      write(*,*) xt23, yt23, x1t, y1t, d1t
      write(*,*) beta2, rho2, d2
      write(*,*) beta1t

      delta1 = 2*ji*a
      delta2 = a*(1+sin(rho2+alpha))
      delta3 = abs(abs(d1t*sin(beta1t + alpha)) 
     &- abs(d2*sin(beta2+alpha)))
c      eps    = delta1 - delta2 - delta3

      write(*,*) delta1, delta2, delta3, eps

      eps = (yc3-abs(dy))-a*(1-2*ji)

      write(*,*) 'eps again', eps, yc3

c      pause

c     flujo laminar

      nsteps=1000 !numero de pasos para hacer la integracion numerica

      write(*,*) 'us*deq*rhof/nu', us*deq*rhof/nu

      if (us*deq*rhof/nu<4) then

c     area proyectada ap
      ap = pi*b*dv/2.0d0

c     velocidad media del fluido

      aux = 0.0d0
      do i = 1,nsteps
      z = eps + dv/nsteps*(0.5d0+(i-1))
      aux = aux + dv/nsteps*(2.0d0*sqrt(b**2*(1.0d0
     &-(z-(dy+eps))**2/dy**2))*z)
      end do !i
      u = rhof/(nu*ap)*us**2*aux


      write(*,*) 'laminar'
      write(*,*) 'ap', ap, 'dy', dy

      else ! turbulento

      z0 = rug*a*2;

c     area proyectada ap

      if (z0.gt.eps) then

      aux = 0.0d0
      do i = 1,nsteps
      z = z0 + (eps+dv-z0)/nsteps*(0.5d0+(i-1))
      aux = aux + (eps+dv-z0)/nsteps*(2*sqrt(b**2*(1
     &-(z-(dy+eps))**2/dy**2)))
      end do !i
      ap = aux  

c     velocidad media del fluido

      aux = 0.0d0
      do i = 1,nsteps
      z = z0 + (eps+dv-z0)/nsteps*(0.5d0+(i-1))
      aux = aux + (eps+dv-z0)/nsteps*(2.5*log(z/z0)
     &*2*sqrt(b**2*(1-(z-(dy+eps))**2/dy**2)))
      end do !i
      u = us/ap*aux
    
      else

      aux = 0.0d0
      do i = 1,nsteps
      z = eps + dv/nsteps*(0.5d0+(i-1))
      aux = aux + dv/nsteps*(2*sqrt(b**2*(1-(z-(dy+eps))**2/dy**2)))
      end do !i
      ap = aux  

c    velocidad media del fluido

      aux = 0.0d0
      do i = 1,nsteps
      z = eps + dv/nsteps*(0.5d0+(i-1))
      aux = aux + dv/nsteps*(2.5*Log(z/z0)*2*sqrt(b**2*(1
     &-(z-(dy+eps))**2/dy**2)))
      end do !ic
      u = us/ap*aux
     
      endif ! eps<z0 or eps >z0
      endif !laminar or turbulent

c     numero de reynolds

      Rey = rhof*u*deq/nu

      write(*,*) 'u', u

c     drag para inclinacion 0 grados

      cd0 = 8.0d0/Rey*1.0d0/sqrt(sphl) + 16.0d0/Rey*1.0d0/sqrt(sph) 
     &+  3.0d0/sqrt(Rey)*1.0d0/sph**(3.0d0/4.0d0) + 0.4210d0**(0.4d0*(
     &-log(sph))**0.2d0)*1.0d0/sphc

c     drag para inclinacion 90 grados

      cd90 = 8.0d0/Rey*1.0d0/sqrt(sphc) + 16.0d0/Rey*1.0d0/sqrt(sph) 
     &+ 3.0d0/sqrt(Rey)*1.0d0/sph**(3.0d0/4.0d0) + 0.4210d0**(0.4d0*(
     &-log(sph))**0.2d0)*1.0d0/sphl

c     drag para cualquier inclinacion

      cd = cd0 + (cd90 - cd0)*Sin(alpha)**3

c     fuerza de drag fd
      
      fd = cd/2*sep*rhof*u**2

      write(*,*) 'fd', fd
      

      f(1,3) = f(1,3) + fd
     
      end
