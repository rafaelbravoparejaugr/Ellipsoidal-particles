
c     rem ******************** 5000
c     rem contact transfer
c     rem ******************** called by main

      subroutine transf (n2,n2b,i9,m3,m3b)

c     inputs -- n0,n2,n2b,m(),i9,m3,m3b,a(),m0(,3),k1(),f()
c     create -- c0(),a(,1-9),o(,1-2),m0(,2),m0(,3),o(,3)
c     output -- c0(),o(,1-3),f0(),m0(,2-3)

      implicit      none

      include      'param.h'
      include      'za.h'
      include      'zb.h'
      include      'zc0.h'
      include      'zdm.h'
      include      'zf0.h'
      include      'zk1.h'
      include      'zm0.h'
      include      'zn.h'
      include      'znp0.h'
      include      'zo.h'

      integer       n2,n2b,i9,m3,m3b

      integer       i,i1, j
      real*8        a1

      real*8        sgn

c     md: maximum number of disks
c     mi: maximum number of invasions

c     n2: number of possible contacts at start of current step
c     n2b: n2 + number of d-p + d-b contacts at start of current step
c     i9: current time step number
c     m3: number of d-d contacts at end of previous step
c     m3b: m3 + number of d-p + d-b contacts at end of previous step
c     d(): coordinates of disk center, their radia/friction angle/accumulated
c          rotation/radia*
c     c0(,1-2): starting and ending invading position no. for ball i
c     a(1-5): save dm(,1-5) of previous step
c     a(6-9): m0(,2), o(i,1-2), f0() of previous step
c     dm(,1-5): contact ball number/boundary no.,dips,contact friction angle
c               at end and start of previous step
c     o(,1-3): distance of penetration, sliding distance w.r.t reference line,
c              and (sign*o(,2)*-o(,1)) at start of current step
c     f0(): lock position transfered from previous step
c     m0(,3): save as m0(,2)
c     m0(,2): contact index - o=open;1=closed/sliding;2=locked

c     Compute c0(,1) and c0(,2) - distributing n2 to each disk

c     write(*,*) 'N2,N2B',n2,n2b
c     write(*,*) 'M3,M3B',m3,m3b

      call timer('transf',30,1)
      do i = 1, nd0
        c0(1,i) = 0.d0
        c0(2,i) = 0.d0
      end do ! i

      do i = 1, n2
        i1 = min( k1(1,int(dm(1,i))), k1(1,int(dm(2,i))))
        c0(1,i1) = c0(1,i1) + 1.d0
      end do ! i

      c0(2,1) = c0(1,1)
      do i = 2, nd0
        c0(2,i) = c0(1,i) + c0(2,i-1)
      end do ! i

      do i = 1, nd0
        c0(1,i) = c0(2,i) - c0(1,i) + 1.d0
      end do ! i

      if (i9.ne.1) then

c       Transfer contacts

c       disk-disk contacts

        do j = 1, n2
          do i = 1, m3
  
c           Contact transfer

            if (int(a(1,1,i)).eq.int(dm(1,j)).and.
     &          int(a(2,1,i)).eq.int(dm(2,j))) then
              m0(2,j) = int(a(3,2,i))
              o(1,j)  = a(1,3,i)
              o(2,j)  = a(2,3,i)
              f0(j)   = a(3,3,i)
              go to 100
            endif
          end do ! i
100       continue
        end do ! j

c       disk-boundary contacts

        do j = n2+1, n2b
          do i = m3+1, m3b

c           Contact transfer

            if (int(a(1,1,i)).eq.int(dm(1,j)).and.
     &          int(a(2,1,i)).eq.int(dm(2,j))) then
              m0(2,j) = int(a(3,2,i))
              o(1,j)  = a(1,3,i)
              o(2,j)  = a(2,3,i)
              f0(j)   = a(3,3,i)
            elseif (abs(m0(2,j)).eq.2.and.int(dm(2,j)).eq.0) then
              dm(3,j) = a(3,1,i) + b(3,i)
            endif
          end do ! i
        end do ! j

c       Save state of contacts at start of current step

        do i = 1, n2b
          m0(3,i) = m0(2,i)
          a1      = o(1,i)
          if (a1.ge.0.d0) a1 = 0.d0
          o(3,i)  = -a1*sgn(o(2,i))
        end do ! i
      endif
      call timer('transf',30,2)

      end
