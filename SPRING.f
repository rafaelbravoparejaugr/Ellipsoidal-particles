
c     rem ******************** 90ab
c     rem submatrix of 1-d constraint/spring points
c     rem ******************** called by main

      subroutine spring (g0)

c     inputs -- n1,n1a,n1b,g(),k1(),k4(,1-2),g0,a(),sp(),f()
c     create -- t(),e(),a(),a(),f()
c     output -- a(),f()

      implicit     none

      include     'param.h'
      include     'za.h'
      include     'ze.h'
      include     'zf.h'
      include     'zg.h'
      include     'zh.h'
      include     'zk1.h'
      include     'zk4.h'
      include     'zn.h'
      include     'zsp.h'
      include     'zt.h'

      logical      iflag
      real*8       g0

      integer      i,i0,i0c,i3, j, l
      real*8       x,y, peng0

c     mp: maximum number of input points/boundary vertices
c     md: maximum number of disks

c     ip(2,1): number of 1-d constraint points
c     ip(3,1): number of spring points
c     g0: penalty (normal contact stiffness)
c     g(,1)/g(,2)/g(,3): coordinates/belonged disk n0. of given points
c     t(): displacement function matrix of point (x,y) in disk i0
c     e(): temporary use for saving stiffness terms
c     k1(): permutation matrix -  k1(old dk no)=new dk no
c     k4(i,1-2): starting number of k(),number of k() for i-th disk
c     a(): stiffness matrix at start of current step
c     sp(): lx,ly,del,k(ratio) of 1-d contraint/spring points/etc
c     f(): free terms at start of current step
c
c     8/16/96 add i0c as cluster number because displac no longer
c             sends back cluster # in i0 parameter

      call timer('spring',28,1)

      iflag = cex.gt.0.0d0
      do i = ip(1,2) + 1, ip(3,2)

        i0 = g(3,i)

c       Assign no# to i0c for use in f(), a(),k4()

        if (i0.le.nd) then
          i0c = k1(1,i0)
        else
          i0c = i0 - nd + nd0
        endif
        x = g(1,i)
        y = g(2,i)

        call displac (i0,x,y)

c       Find its stiffness matrix contribution

        do j = 1, 3
          t(3,j) = sp(1,i)*t(1,j) + sp(2,i)*t(2,j)
        end do ! j

        do j=1, 3
          do l=1, 3
            e(j,l) = t(3,j)*t(3,l)
          end do ! l
        end do ! j

        peng0 = g0*sp(5,i)
        if(iflag) then
          i3    = k4(1,i0c)
          do j = 1, 3
            do l = 1, 3
              a(l,j,i3) = a(l,j,i3) + peng0*e(j,l)
            end do ! l
          end do ! j
        endif

c       Find its free term contribution

        do j = 1, 3
          f(j,i0c) = f(j,i0c) + peng0*sp(3,i)*t(3,j)
        end do ! j

c       Damping terms

c        do l = 1,3
c          peng0 = 0.001*sp(5,i)*g0*h(l,i0c)
c          do j = 1, 3
c            f(j,i0c) = f(j,i0c) + peng0*e(j,l)
c          end do ! l
c        end do ! j

      end do ! i

      call timer('spring',28,2)

      end
