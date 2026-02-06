
c     rem ******************** 6000
c     rem index matrix of non-zero storage,disk order;no graph t. req. for sor
c     rem ******************** called by main

      subroutine sparse (n2,n2p)

c     input  -- n0,c0(),dm(),b0()
c     create -- k4(,1-3),k(),j9,k1(),k2(),m1(),k4(,3),k3(),k()
c     output -- k4(,1-2),k(),k1(),k2()

      implicit      none

      include      'param.h'
      include      'zc0.h'
      include      'zdm.h'
      include      'zg.h'
      include      'zk.h'
      include      'zk1.h'
      include      'zk4.h'
      include      'zn.h'
      include      'znp0.h'

      integer       n2,n2p

      integer       i,i0,i1,i2,i3, j,j0,j2, l, nn

c     md       : maximum number of disks
c     mi       : maximum number of invasions
c     mk       : maximum number of elements in a[]
c     me       : maximum number of contacts per disk
c     mp       : maximum number of points

c     n0       : number of disks+polygons
c     j9       : number of produced connections
c     k4(i,1-3): starting number of k(), number of k(), and limit number of k()
c                for i-th disk
c     k()      : disk contact matrix
c     c0(,1-2) : starting and ending invading position no. for disk i
c     dm(,1-5) : contact ball number/boundary no.,dips,contact friction angle
c                at end and start of previous step
c     k1()     : permutation matrix -  k1(old dk no)=new dk no
c     k2()     : permutation matrix -  k2(new dk no)=old dk no
c     m1()     : temporary use for (1) k() for exchaged disk i (2) number
c                of non-zero terms of k() for each disk
c     k3()     : temporary use to save nonzero terms of k()

c     Set initial values of k4(,1-3) and k() for all nd disks and np polys

      call timer('sparse',27,1)
      do i = 1, n0
        k4(1,i) = 10*(i-1)+1
        k4(2,i) = 1
        k4(3,i) = 10
        i1      = k4(1,i)
        k(i1)   = i
      end do ! i

c     Register disk contacts from lower ordering disks;upper triangle

      do i = 1, nd0
        do j = 1, n2p
          if (j.le.n2.and.c0(1,i).gt.c0(2,i)) goto 210
          if (j.le.n2.and.(j.lt.c0(1,i).or.j.gt.c0(2,i))) goto 210
          i2 = dm(1,j)
          i2 = k1(1,i2)
          j2 = dm(2,j)
          if (j.le.n2) j2 = k1(1,j2)
          if (j.gt.n2) j2 = np0(abs(dm(2,j)))
          if (j.gt.n2.and.i.ne.i2.and.i.ne.j2) goto 210
          if (i2.eq.i) then
            j0 = j2
          else
            j0 = i2
          endif
          i0 = k4(1,i)+k4(2,i)-1
          if (k(i0).eq.j0) goto 210
          k4(2,i) = k4(2,i)+1
          if (k4(2,i).le.k4(3,i)) goto 218

c         Shift k() and k4(,1) when k4(,2) > k4(,3)

          k4(3,i) = k4(3,i)+10
          if (i.eq.n0) goto 218
          do l = k4(1,n0)+k4(2,n0)-1, k4(1,i+1), -1
            k(l+10) = k(l)
          end do ! l
          do l = i+1, n0
            k4(1,l) = k4(1,l)+10
          end do ! l

218       i3 = k4(1,i)+k4(2,i)-1
          k(i3) = j0
210     continue
        end do ! j
      end do ! i

c     Add pin/bolting connections

      do i = ip(3,2)+1, ip(5,2)-1, 2
        i2 = g(3,i)
        if (i2.le.nd) i2 = k1(1,i2)
        j2 = g(3,i+1)
        if (j2.le.nd) j2 = k1(1,j2)
        i0 = i2
        j0 = j2
        if (i2.gt.j2) then
          i0 = j2
          j0 = i2
        endif
        do j  =  k4(1,i0)+1, k4(1,i0)+k4(2,i0)-1
          if (k(j).eq.j0) goto 250
        end do ! j
        k4(2,i0) = k4(2,i0)+1
        if (k4(2,i0).le.k4(3,i0)) goto 258
        k4(3,i0) = k4(3,i0)+10
        if (i0.eq.n0) goto 258

        do j = k4(1,n0)+k4(2,n0)-1, k4(1,i0+1), -1
          k(j+10) = k(j)
        end do ! j

        do j = i0+1, n0
          k4(1,j) = k4(1,j)+10
        end do ! j

258     i3 = k4(1,i0)+k4(2,i0)-1
        k(i3) = j0

250     continue
      end do ! i

      call exceed (k4(1,n0)+k4(3,n0)-1,mk,'mk','in sparse......')

c     Delete zero-elements in k()

      nn = 0
      do i = 1, n0
        do j = k4(1,i), k4(1,i)+k4(2,i)-1
          nn = nn+1
          k(nn) = k(j)
        end do ! j
      end do ! i
      k4(1,1) = 1
      do i = 2, n0
        k4(1,i) = k4(1,i-1)+k4(2,i-1)
      end do ! i
      call timer('sparse',27,2)

      end
