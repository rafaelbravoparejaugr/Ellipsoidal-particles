
c     rem ******************** 25000
c     rem disk displacement function matrix
c     rem ******************** called by fixed,ptload,pun1,judgea,ratio,result

      subroutine displac (i0,x,y)

c     inputs -- x,y,i0,d()
c     create -- t()
c     output -- t()

      implicit     none

      integer      i0, i0c,ix
      real*8       x,y

      include     'param.h'
      include     'zdz.h'
      include     'zh0.h'
      include     'zk1.h'
      include     'zn.h'
      include     'zt.h'

c     md: maximum number of disks

c     id: 1 for x,y as cartesian coordinates; 2 for x,y as polar coordinates
c     i0: old disk number
c     d(,1-3): coordinates of center, radius of the specified disk
c     t(): displacement function matrix of point (x,y) in disk i0

c     displace asks for disk # but sends back cluster # (for disks)
c     8/14/96  correct above, do not send back cluster # (don't change i0)

c     call timer('displa',6,1)

      t(1,1) = 1.d0
      t(2,1) = 0.d0
      t(1,2) = 0.d0
      t(2,2) = 1.d0

      if (i0.le.nd) then
        i0c    = k1(1,i0)
        t(1,3) = dz(2,i0c) - y
        t(2,3) = x - dz(1,i0c)

c     This assumes that what is sent in for polys is i0 = nd+poly #
c     therefore ix = no # which is correct for h()

      else
        ix     = i0 - nd + nd0
        t(1,3) = h0(3,ix)/h0(1,ix) - y
        t(2,3) = x - h0(2,ix)/h0(1,ix)
      endif

c     call timer('displa',6,2)
      end
