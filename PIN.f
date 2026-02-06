
c     rem ******************** 900c - new
c     rem submatrix of pin constraint points as jointing
c     rem ******************** called by main

      subroutine pin (g0)

c     in -- n1,n1a,n1b,g(),k1(),k4(,1-2),g0,a(),sp(),f()
c     create -- t(),e(),a(),f()
c     out -- a(),f()

      implicit     none

      include     'param.h'

      logical      iflag
      real*8       g0, peng0

      integer      i,i1,i1c,i2,i2c,i3, j, l
      real*8       x,y, dxx,dyy,dll, dn,ds, fx,fy

      integer      indf

      include     'za.h'
      include     'zc.h'
      include     'zd.h'
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

c     ip(4,1): number of pin constraint points
c     g0: penalty (normal contact stiffness)
c     g(,1)/g(,2)/g(,3): coordinates/belonged disk no. of given points
c     t(): displacement function matrix of point (x,y) in disk i0
c     e(): temporary use for saving stiffness terms
c     k1(): permutation matrix -  k1(old dk no)=new dk no
c     k4(i,1-2): starting number of k(),number of k() for i-th disk
c     a(): stiffness matrix at start of current step
c     sp(): lx,ly,del,k(ratio) of 1-d contraint/spring points/etc
c     f(): free terms at start of current step
c     8/14/96  add i1c and i2c as cluster #'s of disksi1&i2 because
c              displac no longer sends back cluster#

c     This routine assumes pin connections between disks only not polys

      call timer('pin',19,1)
      iflag = cex.gt.0.0d0
      do i = ip(3,2)+1, ip(4,2), 2

        peng0 = sp(5,i)*g0

        i1    = g(3,i)
        i1c   = k1(1,i1)
        i2    = g(3,i+1)
        i2c   = k1(1,i2)

c       Check if pins of disks fail at the start of current step

        if (i1.le.nd.and.i2.le.nd) then

          fx = (-c(1,i+1) + c(1,i) + sp(3,i))*peng0
          fy = (-c(2,i+1) + c(2,i) + sp(4,i))*peng0

          dxx = d(1,i2) - d(1,i1)
          dyy = d(2,i2) - d(2,i1)
          dll = sqrt(dxx*dxx + dyy*dyy)

          dn  = ( dxx*fx + dyy*fy)/dll
          ds  = (-dyy*fx + dxx*fy)/dll

          if (dn.gt.sp(1,i)) then
            goto 100

          elseif (abs(ds).gt.-dn*min(d(4,i1),d(4,i2))+sp(2,i)) then
            goto 100
          endif
        endif

c       1st point of pin constraint points

        x = g(1,i)
        y = g(2,i)

        call displac (i1,x,y)

        do j = 1, 3
          e(4,j) = t(1,j)
          e(5,j) = t(2,j)
        end do ! j

c       2nd point of pin constraint points

        x = g(1,i+1)
        y = g(2,i+1)
        call displac (i2,x,y)

        if(iflag) then

c         Find stiffness matrix, k11
          do j = 1, 3
            do l = 1, 3
              e(j,l) = e(4,j)*e(4,l) + e(5,j)*e(5,l)
            end do ! l
          end do ! j

          i3 = k4(1,i1c)
          do j = 1, 3
            do l = 1, 3
              a(l,j,i3) = a(l,j,i3) + peng0*e(j,l)
            end do ! l
          end do ! j

c         Find stiffness matrix, k22

          do j = 1, 3
            do l = 1, 3
              e(j,l) = t(1,j)*t(1,l) + t(2,j)*t(2,l)
            end do ! l
          end do ! j

          i3 = k4(1,i2c)
          do j = 1, 3
            do l = 1, 3
              a(l,j,i3) = a(l,j,i3) + peng0*e(j,l)
            end do ! l
          end do ! j

c         Find stiffness matrix, k12 (upper triangle)

          if (i1c.lt.i2c) then
            i3  =  indf(i1c,i2c)
            do j = 1, 3
              do l = 1, 3
                e(j,l) = e(4,j)*t(1,l) + e(5,j)*t(2,l)
              end do ! l
            end do ! j
            do j = 1, 3
              do l = 1, 3
                a(l,j,i3) = a(l,j,i3) - peng0*e(j,l)
              end do ! l
            end do ! j
          else

            i3  =  indf(i2c,i1c)
            do j = 1, 3
              do l = 1, 3
                e(j,l) = t(1,j)*e(4,l) + t(2,j)*e(5,l)
              end do ! l
            end do ! j

            do j = 1, 3
              do l = 1, 3
                a(l,j,i3) = a(l,j,i3) - peng0*e(j,l)
              end do ! l
            end do ! j
          endif

        endif ! iflag

c       Find free term, f1, f2

        dxx = c(1,i+1) - c(1,i) - sp(3,i)
        dyy = c(2,i+1) - c(2,i) - sp(4,i)
        if (dxx.ne.0.d0 .or. dyy.ne.0.d0) then

          do j = 1, 3
            f(j,i1c) = f(j,i1c) - peng0*(dxx*e(4,j) + dyy*e(5,j))
            f(j,i2c) = f(j,i2c) + peng0*(dxx*t(1,j) + dyy*t(2,j))
          end do ! j

        endif

100     continue

      end do ! i

      call timer('pin',19,2)
      end
