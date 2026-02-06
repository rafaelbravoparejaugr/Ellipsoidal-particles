
c     rem ******************** 900x
c     rem submatrix of no-rotation points
c     rem ******************** called by main

      subroutine xrotation (g0)

c     inputs -- n1,n1a,n1b,g(),k1(),k4(,1-2),g0,a(),sp(),f()
c     create -- t(),e(),a(),s(),a(),f()
c     output -- a(),f()

      implicit     none

      include     'param.h'
      include     'za.h'
      include     'zdz.h'
      include     'ze.h'
      include     'zf.h'
      include     'zg.h'
      include     'zh0.h'
      include     'zk.h'
      include     'zk1.h'
      include     'zk4.h'
      include     'zn.h'
      include     'zs.h'
      include     'zsp.h'
      include     'zt.h'

      real*8       g0, peng0

      integer      i,i0,i0c,i3
      real*8       x0,y0,dl2

c     mp: maximum number of input points/boundary vertices
c     md: maximum number of disks

c     ip(6,1)  : number of no rotation points
c     g(,1)/g(,2)/g(,3): coordinates/belonged disk no. of given points
c     t()      : displacement function matrix of point (x,y) in disk i0
c     e()      : temporary use for saving stiffness terms
c     k1()     : permutation matrix -  k1(old dk no)=new dk no
c     k4(i,1-2): starting number of k(),number of k() for i-th disk
c     a()      : stiffness matrix at start of current step
c     sp(,3)   : theta of no-rotation points
c                fixed angular spring stiffness: 1000*g0*r
c     s()      : temporary use for saving free terms
c     f()      : free terms at start of current step

c     8/14/96  add i0c as cluster number of i0 to be
c              consistent with rest of code

      call timer('xrotatio',33,1)

      do i = ip(5,2)+1, ip(6,2)

	    peng0 = g0*sp(5,i)
        i0    = g(3,i)
        if (i0.le.nd) then
          i0c = k1(1,i0)
          x0  = dz(1,i0c)
          y0  = dz(2,i0c)
        else
          i0c = i0 - nd + nd0
          x0  = h0(2,i0c)/h0(1,i0c)
          y0  = h0(3,i0c)/h0(1,i0c)
        endif
        dl2 = (x0 - g(1,i))**2 + (y0 - g(2,i))**2

c       Check if bolt fails at the start of current step

        if (peng0*dl2*sp(3,i).le.sp(6,i)) then

          i3 = k4(1,i0c)

c         Add to rotation terms

          if(cex.gt.0.0d0) a(3,3,i3) = a(3,3,i3) + peng0*dl2

c         Find free term, f1 (moment)

          if (sp(3,i).ne.0.0d0) then
            f(3,i0c) = f(3,i0c) + peng0*dl2*sp(3,i)
          endif
        else
          write (4,*) 'no-rot. const.:',i,' fail! m=',sp(3,i)*peng0
        endif
      end do ! i

      call timer('xrotatio',33,2)

      end
