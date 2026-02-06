
c     rem ******************** 14000
c     rem mass force
c     rem ******************** called by main

      subroutine mass

c     inputs -- o3,o4,n0,h0(,2),f()
c     create -- f()
c     output -- f()

      implicit     none

      integer      i

      include     'param.h'
      include     'zf.h'
      include     'zh0.h'
      include     'zn.h'
      include     'zp.h'

c     md: maximum number of disks

c     o3: unit mass force in x dirction
c     o4: unit mass force in y dirction
c     h0(,2): disk mass
c     n0: number of disks
c     f(): free terms before step iterations

      call timer('mass',15,1)

      do i = 1, nd0
        f(1,i) = f(1,i) + o3*h0(2,i)
        f(2,i) = f(2,i) + o4*h0(2,i)
      end do ! i
      do i = nd0+1, n0
        f(1,i) = f(1,i) + o0(i)*o3*h0(1,i)
        f(2,i) = f(2,i) + o0(i)*o4*h0(1,i)
      end do ! i

      call timer('mass',15,2)

      end
