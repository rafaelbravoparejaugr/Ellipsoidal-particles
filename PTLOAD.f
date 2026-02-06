
c     rem ******************** 12000
c     rem point loading
c     rem ******************** called by main

      subroutine ptload

c     inputs -- n6,n1,n1a,n1b,n1c,n7,g(),k1(),f()
c     create -- t(),s(),f()
c     output -- f()

      implicit     none

      include     'param.h'
      include     'zf.h'
      include     'zg.h'
      include     'zk1.h'
      include     'zn.h'
      include     'zs.h'
      include     'zsp.h'
      include     'zt.h'

      integer      i,i0,i0c, j
      real*8       x,y

c     mp: maximum number of input points/boundary vertices
c     md: maximum number of disks

c     ip(7,1): number of load points
c     g(,1)/g(,2)/g(,3): coordinates/belonged disk n0. of given points
c     sp(,3)/sp(,4): u v or point load in x,y dirction for given points
c     t(): displacement function matrix of point (x,y) in disk i0
c     s(): temporary use for saving free terms
c     k1(): permutation matrix -  k1(old dk no)=new dk no
c     f(): free terms before step iterations
c     8/16/96 add i0c as cluster number because displac no longer
c             sends back cluster # in i0 parameter

      call timer('ptload',21,1)

      do i = ip(6,2)+1, ip(7,2)

        i0 = g(3,i)

c       Assign n0# to i0c for use in f()

        if (i0.le.nd) then
          i0c = k1(1,i0)
        else
          i0c = i0-nd+nd0
        endif

        x = g(1,i)
        y = g(2,i)

        call displac (i0,x,y)

        do j = 1, 3
          s(j) = t(1,j)*sp(3,i)+t(2,j)*sp(4,i)
        end do ! j

c       Add to free terms

        do j = 1, 3
          f(j,i0c) = f(j,i0c)+s(j)
        end do ! j
      end do ! i

      call timer('ptload',21,2)

      end
