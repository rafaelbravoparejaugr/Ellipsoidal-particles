
c     rem ******************** 8000
c     rem inertia terms
c     rem ******************** called by main

      subroutine inertia (n3,n2b,i00,g3,NBeta,NGamma)

c     inputs -- n0,k4(,1-2),n2b,g3,o0(),i00,k1(),h0(),h()
c     create -- n3,a(),f(),m0(,1),g4,g5,q(),s01-s02,q(),a(),f()
c     output -- n3,m0(,1),a(),f()

      implicit     none

      integer      n3,n2b,i00
      real*8       g3

      integer      i,i2, j, l
      real*8       g4,g5,g6, x0,y0, s01,s02,NGamma,NBeta

      include     'param.h'
      include     'za.h'
      include     'zf.h'
      include     'zh.h'
      include     'zh0.h'
      include     'zk4.h'
      include     'zm0.h'
      include     'zn.h'
      include     'zp.h'
      include     'zq.h'

c     md: maximum number of disks
c     mi: maximum number of invasions

c     n0: number of disks
c     n3: total number of non-zero elements in a()=total no. of terms in k()
c     n2b: n2 + number of d-b contacts at start of current step
c     i00: static/dynamic control code - 0: for static (g5=0);
c                                        1: for implicit dynamic
c                                        2: for explicit dynamic
c     g3: current time step used
c     g4: initial inertia coefficients for k[]
c     g5: initial inertia coefficients for f[]
c     k4(i,1-2): starting number of k(),number of k() for i-th disk
c     a(): stiffness matrix at start of current step
c     f(): free terms at start of current step
c     m0(,1): contact indices before step iterations
c     q(): temporary use for save geometry terms of inertia
c     k1(): permutation matrix -  k1(old dk no)=new dk no
c     h0(,2): disk mass
c     d(,3): radius of the specified disk
c     h() : velocity vector of each disk at start of current step
c     ha(): acceleration vector of each disk at start of current step
c
c     Set initial values of a(),f(),m0(,1) to zero

      call timer('inertia',11,1)

      n3 = k4(1,n0) + k4(2,n0)-1
      do i = 1, n3
        do j = 1, 3
          do l = 1,3
            a(l,j,i) = 0.d0
          end do ! l
        end do ! j
      end do ! i

      do i = 1, n0
        do j = 1, 3
          f(j,i) = 0.d0
        end do ! j
      end do ! i

      do i = 1, n2b
        m0(1,i) = 0
      end do ! i

c     Set initial inertia coefficients

      if(i00.eq.1) then
	g6  = (1-2*NBeta)/(2*NBeta)
        g5  = 1/(g3*NBeta)
        g4  = 1/(NBeta*g3**2)
        cex = 1.0d0
      elseif(i00.eq.2) then
	g4  = 1.0d0
        g5  = 1.0d0
	cex = 0.0d0
      endif

c     Compute inertia terms

      do i = 1, n0
        i2 = k4(1,i)
        do l = 1, 3
          do j = 1, 3
            q(j,l) = 0.d0
          end do ! j
        end do ! l

        if (i.le.nd0) then
          q(1,1) = h0(2,i)
          q(2,2) = h0(2,i)
          q(3,3) = h0(3,i)

        else

          x0     = h0(2,i)/h0(1,i)
          y0     = h0(3,i)/h0(1,i)
          s01    = h0(4,i)-x0*h0(2,i)
          s02    = h0(5,i)-y0*h0(3,i)
          q(1,1) = h0(1,i)*o0(i)
          q(2,2) = h0(1,i)*o0(i)
          q(3,3) = (s01+s02)*o0(i)

        endif

        do j = 1, 3
          do l = 1, 3
            a(l,j,i2) = a(l,j,i2) + g4*q(j,l)
          end do ! l
        end do ! j

        if (i00.eq.0) then
          do j = 1, 3
            h (j,i) = 0.d0
            ha(j,i) = 0.d0
          end do ! j
        elseif(i00.eq.1) then
          do j = 1, 3
            do l = 1, 3
              f(j,i) = f(j,i) + g5*q(j,l)*h(l,i) + g6*q(j,l)*ha(l,i)
            end do ! l
          end do ! j
        elseif(i00.eq.2) then
          do j = 1, 3
            do l = 1, 3
              f(j,i) = f(j,i) - g5*q(j,l)*ha(l,i)
            end do ! l
          end do ! j
        endif

      end do ! i

      call timer('inertia',11,2)
      end
