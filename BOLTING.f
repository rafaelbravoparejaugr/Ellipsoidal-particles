c     rem ******************** 90bt
c     rem submatrix of bolting points
c     rem ******************** called by main

      subroutine bolting (g0)

c     inputs -- n1,n1a,n1b,g(),k1(),k4(,1-2),g0,a(),sp(),f()
c     create -- t(),e(),a(),a(),f()
c     output -- a(),f()

      implicit     none

      integer      i,i1,i2,i3,i1c,i2c, j, l
      real*8       g0, x,y, peng0

      integer      indf

      include     'param.h'
      include     'za.h'
      include     'ze.h'
      include     'zf.h'
      include     'zg.h'
      include     'zk.h'
      include     'zk1.h'
      include     'zk4.h'
      include     'zn.h'
      include     'zsp.h'
      include     'zt.h'

c     mp: maximum number of input points/boundary vertices
c     md: maximum number of disks

c     ip(5,1): number of bolting points
c     g(,1)/g(,2)/g(,3): coordinates/belonged disk no. of given points
c     t(): displacement function matrix of point (x,y) in disk i0
c     e(): temporary use for saving stiffness terms
c     k1(): permutation matrix -  k1(old dk no) = new dk no
c     k4(i,1-2): starting number of k(),number of k() for i-th disk
c     a(): stiffness matrix at start of current step
c     sp(): lx,ly,del,k of bolting points
c     f(): free terms at start of current step

      call timer('boltin',2,1)

      do i = ip(4,2) + 1, ip(5,2), 2

        peng0 = g0*sp(5,i)

c       1st bolting point

        i1 = g(3,i)
        if (i1.le.nd) then
          i1c = k1(1,i1)
        else
          i1c = i1 - nd + nd0
        endif
        x = g(1,i)
        y = g(2,i)

        call displac (i1,x,y)

        do j = 1, 3
          e(4,j) = sp(1,i)*t(1,j) + sp(2,i)*t(2,j)
        end do ! j

c       2nd bolting point

        i2 = g(3,i+1)
        if (i2.le.nd) then
          i2c = k1(1,i2)
        else
          i2c = i2 - nd + nd0
        endif
        x = g(1,i+1)
        y = g(2,i+1)

        call displac (i2,x,y)

        do j = 1, 3
          e(5,j) = sp(1,i)*t(1,j) + sp(2,i)*t(2,j)
        end do ! j

c       Check if bolt fails at the start of current step

        if (sp(3,i)*peng0.le.sp(6,i)) then

          if(cex.gt.0.0d0) then

c           Find stiffness matrix, k11

            do j = 1, 3
              do l = 1, 3
                e(j,l) = e(4,j)*e(4,l)
              end do ! l
            end do ! j
            i3 = k4(1,i1c)
            do j = 1, 3
              do l = 1, 3
                a(l,j,i3) = a(l,j,i3) + peng0*e(j,l)
              end do ! l
            end do ! j

c           Find stiffness matrix, k22

            do j = 1, 3
              do l = 1, 3
                e(j,l) = e(5,j)*e(5,l)
              end do ! l
            end do ! j
            i3 = k4(1,i2c)
            do j = 1, 3
              do l = 1, 3
                a(l,j,i3) = a(l,j,i3) + peng0*e(j,l)
              end do ! l
            end do ! j

c           Find stiffness matrix, k12 (upper triangle)

            if (i1c.lt.i2c) then
              i3 = indf(i1c,i2c)
              do j = 1, 3
                do l = 1, 3
                  e(j,l) = e(4,j)*e(5,l)
                end do ! l
              end do ! j
              do j = 1, 3
                do l = 1, 3
                  a(l,j,i3) = a(l,j,i3) - peng0*e(j,l)
                end do ! l
              end do ! j
            else
              i3 = indf(i2c,i1c)
              do j = 1, 3
                do l = 1, 3
                e(j,l) = e(5,j)*e(4,l)
                end do ! l
              end do ! j
              do j = 1, 3
                do l = 1, 3
                  a(l,j,i3) = a(l,j,i3) - peng0*e(j,l)
                end do ! l
              end do ! j
            endif
          endif

c         Find free term, f1, f2

          if (sp(3,i).ne.0.d0) then
            peng0 = peng0*sp(3,i)
            do j = 1, 3
              f(j,i1c) = f(j,i1c) + peng0*e(4,j)
              f(j,i2c) = f(j,i2c) - peng0*e(5,j)
            end do ! j
          endif
        endif

      end do ! i

      call timer('boltin',2,2)

      end
