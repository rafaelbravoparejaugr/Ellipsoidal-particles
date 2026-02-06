
c     rem ******************** 29000
c     rem ratio displacements and iteration drawing
c     rem ******************** called by main

      subroutine ratio (g2,q01,q03)

c     inputs -- n0,k0b(),d(),f(),g2,w(),m9,n3,b(),r()
c     create -- u(),q01,a(),f(),u()
c     output -- q01,a(),f(),u()

      implicit      none

      include      'param.h'
      include      'zb.h'
      include      'zbd.h'
      include      'zdz.h'
      include      'zk0b.h'
      include      'zn.h'
      include      'zt.h'
      include      'zu.h'
      include      'zw.h'
      include      'zxbd.h'

      real*8        g2,q01,q03

      integer       i,il,i1,i2, j, l, nn
      real*8        a1,a3, x,y, dx,dy, cc1,ss

c     md: maximum number of disks

c     n0: number of disks
c     m9: current iteration number (-1 for iteration funished)
c     n3: total number of non-zero terms in a()
c     g2: given maximum vertex displacement ratio
c     a1: maximum vertex displacement of current iteration for current step
c     i0: 1 for close->open; -1 for open->close
c     k0(): initial argument number of disk i in d()
c     d(): coordinates of disk center, their radia/friction angle/accumulated
c          rotation/radia*
c     k1(): permutation matrix -  k1(old bk no)=new bk no
c     f(): (1) displacement solutions; (2) recovered to be r()
c     u(): (1) disk center coord. due to current iteration; (2)
c          new coord. of disk center due to current iteration
c     w(): window limits
c     q01: ratio of max. vertex displac. to given value per step
c     b(): stiffness matrix before triangular decomposition
c     a(): recovered to be b()
c     r(): free term before triangular decomposition
c     ***note: in this routine b(,1-3) are the sol. displac. terms
c              a() and r() are not used.  above 3 lines probably
c              from another subroutine.

c     Compute maximum ratio of translation displacements at this iteration
c       /rotation

      call timer('ratio',24,1)

      a1 = 0.d0
      a3 = 0.d0

c     For disk groups

      do i = 1, nd0
        do j = 1, 3
          u(j,i) = b(j,i)
        end do ! j

        a1 = max(a1,sqrt(u(1,i)*u(1,i) + u(2,i)*u(2,i)))
        a3 = max(a3,abs(b(3,i)))

      end do ! i

c     For polygons

      il = nd0 + 1
      nn = nd - nd0
      do i = il, n0
        i1 = k0b(1,i-nd0)
        i2 = bd(1,i1)
        do j = i1+1, i1+i2+1

          x = bd(1,j)
          y = bd(2,j)

          call displac (i+nn,x,y)

          xbd(1,j) = 0.d0
          xbd(2,j) = 0.d0

          il = 2

c         Linearized updating

          if (il.eq.1) then
            do l = 1, 3
              xbd(1,j) = xbd(1,j) + t(1,l)*b(l,i)
              xbd(2,j) = xbd(2,j) + t(2,l)*b(l,i)
            end do ! l

c         True updating

          else
            xbd(1,j) = xbd(1,j) + b(1,i)
            xbd(2,j) = xbd(2,j) + b(2,i)
            dx = t(2,3)
            dy = -t(1,3)
            cc1 = cos(b(3,i)) - 1.d0
            ss  = sin(b(3,i))
            xbd(1,j) = xbd(1,j) + dx*cc1 - dy*ss
            xbd(2,j) = xbd(2,j) + dx*ss  + dy*cc1
          endif
          a1 = max(a1,sqrt(xbd(1,j)*xbd(1,j) + xbd(2,j)*xbd(2,j)))
        end do ! j
      end do ! i
      q01 = a1/(g2*(w(4) - w(2)))
      q03 = a3

c     Compute new coord. of disk group centers
c     *** u(,1-2) on right hand side below are displacements
c     *** u(,1-2) on left hand side are new coord. of clusters

      do i = 1, nd0
        u(1,i) = dz(1,i) + u(1,i)
        u(2,i) = dz(2,i) + u(2,i)
	u(3,i) = dz(5,i) + u(3,i)
      end do ! i

      call timer('ratio',24,2)
      end
