
c     rem ******************** 9000
c     rem submatrix of fixed points
c     rem ******************** called by main

      subroutine fixed (g0)

c     in -- n1,g(),k1(),k4(,1-2),g0,a(),f(),c()
c     create -- t(),e(),a(),a(),f()
c     out -- a(),f()

      implicit     none

      logical      iflag
      real*8       g0, peng0

      integer      i,i0,i0c,i3, j, l
      real*8       dxx,dyy, x,y

      include     'param.h'
      include     'za.h'
      include     'zc.h'
      include     'ze.h'
      include     'zf.h'
      include     'zg.h'
      include     'zk1.h'
      include     'zk4.h'
      include     'zn.h'
      include     'zsp.h'
      include     'zt.h'

c     mp: maximum number of input points/boundary vertices
c     md: maximum number of disks
c     ma: maximum number of elemts in [a]

c     ip(1,1): number of fixed points
c     g0: penalty (normal contact stiffness)
c     g(,1)/g(,2)/g(,3): coordinates/belonged disk of given points
c     c(): displacement errors for fixed/spring/measured points
c     t(): displacement function matrix of point (x,y) in disk i0
c     e(): temporary use for saving stiffness terms
c     k1(): permutation matrix -  k1(old dk no)=new dk no
c     k4(i,1-2): starting number of k(),number of k() for i-th disk
c     a(): stiffness matrix at start of current step
c     f(): free terms at start of current step
c     8/14/96 add i0c as cluster# because displac no longer sends back
c             cluster # as i0 term

      call timer('fixed',9,1)

	  iflag = cex.gt.0.0d0

      do i = 1, ip(1,2)
        i0 = g(3,i)

c       Assign no# to i0c for use in f() and k4()

        if (i0.le.nd) then
          i0c = k1(1,i0)
        else
          i0c = i0 - nd + nd0
        endif
        x = g(1,i)
        y = g(2,i)

        call displac (i0,x,y)

c       Find stiffness matrix contribution

        peng0 = sp(5,i)*g0
        do j = 1, 3
          do l = 1, 3
            e(j,l) = t(1,j)*t(1,l) + t(2,j)*t(2,l)
          end do ! l
        end do ! j
        if(iflag) then
          i3 = k4(1,i0c)
          do j = 1, 3
            do l = 1, 3
              a(l,j,i3) = a(l,j,i3) + peng0*e(j,l)
            end do ! l
          end do ! j
        endif

c       Find free term contribution

        dxx = c(1,i) + sp(3,i)
        dyy = c(2,i) + sp(4,i)
        if (dxx.ne.0.0d0 .or. dyy.ne.0.0d0) then
          do j = 1, 3
            f(j,i0c) = f(j,i0c) + peng0*(dxx*t(1,j) + dyy*t(2,j))
          end do ! j
        endif
      end do ! i

      call timer('fixed',9,2)
      end
