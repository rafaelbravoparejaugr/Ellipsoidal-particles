
c     rem ******************** 21000
c     rem compute displacement and errors
c     rem ******************** called by main

c     subroutine result (g3,i9,q02,irstep,i00)
      subroutine result (g2,g3,i9,q01,q02,q03,irstep,i00,NBeta,NGamma)

c     inputs -- n0,k0(),d(),u(),n1,n7,n6,g(),k1(),f(),sp(),c(),e0()
c               v0(),c11(),c22(),c12(),h()
c     create -- d(),t(),g(),sp(),c(),c11(),c22(),c12(),h()
c     output -- d(),g(),sp(),c(),c11(),c22(),c12(),h()

      implicit     none

      include     'param.h'
      include     'zb.h'
      include     'zbd.h'
      include     'zc.h'
      include     'zd.h'
      include     'zdz.h'
      include     'zg.h'
      include     'zh.h'
      include     'zk0b.h'
      include     'zk1.h'
      include     'zk2.h'
      include     'zn.h'
      include     'zrotat.h'
      include     'zsp.h'
      include     'zt.h'
      include     'zu.h'
      include     'zxbd.h'

      integer      i9,irstep,i00
c     real*8       g3,q02
      real*8       g2,g3,q01,q02,q03
      real*8       g3bar,NGamma,NBeta,haold

      integer      i,i0,i0c,i1,i1c,i2,il, j,ju
      real*8       x,y, x1,y1, dx,dy, dxx,dyy,dll, cc1,ss

c     mp: maximum number of input points/boundary vertices
c     md: maximum number of disks

c     n0      : Number of disks
c     n1      : Number of fixed points
c     n1a     : Number of 1-d constraint points
c     n1b     : Number of spring points
c     n1c     : Number of one-direction moving constraint points
c     nbt     : Number of bolting points
c     n1x     : Number of no rotation points
c     n7      : Number of measured points
c     n6      : Number of load points
c     g3      : Actual time step used
c     g(,1-3) : Coordinates/belonged disk no. of given points
c     sp()    : lx,ly,del,k(ratio) of 1-d constraint/spring points/etc
c     c()     : Displacement errors for fixed/measured points
c     k1()    : Permutation matrix -  k1(old dk no)=new dk no
c     t()     : Displacement function matrix of point (x,y) in disk i0
c     u()     : New positions of vertices due to current converging iteration
c     d()     : Coordinates of disk center, radia/friction angle/accumulated
c                 rotation/radia*
c     f()     : Displacement solutions
c     h()     : Velocity vector of disks
c     ha()    : Acceleration vector of disks
c     accrot(): Accumulated rotation for each particle

c     8/13/96 correction: b() array expecting cluster # not disk #
c             added i0c for use in b() array for points

c     Compute new positions/errors of points

      call timer('result',26,1)

c     Update disk group and polygon velocity vectors for use in next step

      if(i00.eq.1) then

c	 write(*,*) NBeta, NGamma
c	 pause

        do i = 1, n0
          do j = 1, 3

		  haold = ha(j,i)

	      ha(j,i) = 1/(NBeta*g3**2)*(b(j,i)-g3*h(j,i))
     &               -(1-2*NBeta)/(2*NBeta)*haold 

		  h(j,i)  = NGamma/(NBeta*g3)*b(j,i) + (1-NGamma/NBeta)*h(j,i)
     &              + (1-NGamma/(2*NBeta))*g3*haold


          end do ! j
        end do ! i

      elseif(i00.eq.2) then

        g3bar = 0.5d0*(g3 + g3old)
        do i = 1, n0
          do j = 1, 3

            h(j,i)  = h(j,i) + g3bar*ha(j,i)
            b(j,i)  = g3*h(j,i)


          end do ! j
        end do ! i
		call ratio(g2,q01,q03)
      endif

      do i = 1, nt

        i0 = int(g(3,i))

c       Assign no# to i0c for use in b()

        if (i0.le.nd) then
          i0c = k1(1,i0)
        else
          i0c = i0 - nd + nd0
        endif
        x = g(1,i)
        y = g(2,i)
        call displac (i0,x,y)
        x1 = 0.d0
        y1 = 0.d0
        il = 2
        if (il.eq.1) then

c         Linearized updating

          do j = 1, 3
            x1 = x1 + t(1,j)*b(j,i0c)
            y1 = y1 + t(2,j)*b(j,i0c)
          end do ! j
        else

c         True updating

          x1  =  x1 + b(1,i0c)
          y1  =  y1 + b(2,i0c)
          cc1 =  cos(b(3,i0c)) - 1.d0
          ss  =  sin(b(3,i0c))
          x1  =  x1 + t(2,3)*cc1 + t(1,3)*ss
          y1  =  y1 + t(2,3)*ss  - t(1,3)*cc1
        endif

c       Compute new del for 1-d const./spring/bolting/no-rot. points

        if (i.le.ip(1,2)) goto 213
        if (i.gt.ip(3,2)) goto 211

        sp(3,i) = sp(3,i) - sp(1,i)*x1 - sp(2,i)*y1

        goto 213

211     if (i.le.ip(4,2).or.i.gt.ip(5,2)) goto 212

        if (mod(i-ip(4,2),2).eq.1) then
          sp(3,i) = sp(3,i) - sp(1,i)*x1 - sp(2,i)*y1
        else
          sp(3,i-1) = sp(3,i-1) + sp(1,i-1)*x1 + sp(2,i-1)*y1
        endif

        goto 213

212     if (i.gt.ip(6,2)) goto 213
        sp(3,i) = sp(3,i) - b(3,i0c)

213     g(1,i) = x + x1
        g(2,i) = y + y1
        c(1,i) = c(1,i) - x1
        c(2,i) = c(2,i) - y1
      end do ! i

c     Update (lx,ly) of bolting lines for next step

      do i = ip(4,2)+1, ip(5,2), 2
        dxx = g(1,i+1) - g(1,i)
        dyy = g(2,i+1) - g(2,i)
        dll = 1.d0/sqrt(dxx*dxx + dyy*dyy)
        sp(1,i) = dxx*dll
        sp(2,i) = dyy*dll
      end do ! i

c     Update individual disk geometry after one step calculation

      do i = 1, nd
        i1c = k1(1,i)
        if (k2(2,i1c).eq.1) then
          d(1,i) = u(1,i1c)
          d(2,i) = u(2,i1c)
          d(5,i) = u(3,i1c)
        else
          dx     = d(1,i) - dz(1,i1c)
          dy     = d(2,i) - dz(2,i1c)
          cc1    = cos(b(3,i1c)) - 1.d0
          ss     = sin(b(3,i1c))
          d(1,i) = d(1,i) + b(1,i1c) + dx*cc1 - dy*ss
          d(2,i) = d(2,i) + b(2,i1c) + dx*ss  + dy*cc1
        endif
      end do ! i

c     Update disk group geometry after one step calculation

      do i = 1, nd0
        dz(1,i) = u(1,i)
        dz(2,i) = u(2,i)
        dz(5,i) = u(3,i)
      end do ! i

c     Track accumulated rotations of particles

      do i = 1,nd0
        if(i9.eq.1) then
          accrot(i) = b(3,i)
        else
          accrot(i) = accrot(i) + d(5,i)
        endif
      end do ! i

c     Output state

      if (irstep.ne.0) then
        if (mod(i9,irstep).eq.0) then
          do i = 1,nd0
            write (17,2000) q02,u(1,i),u(2,i),accrot(i)
          end do ! i
        endif
      endif

c     Update coordiates of polygon/boundary vertices

      do i = 1, np+nb
        i1 = k0b(1,i)
        i2 = bd(1,i1)
        ju = i1 + i2
        if (i.le.np) ju = i1 + i2 + 1
        do j = i1+1, ju
          bd(1,j) = bd(1,j) + xbd(1,j)
          bd(2,j) = bd(2,j) + xbd(2,j)
        end do ! j
      end do ! i

      call timer('result',26,2)
c     Formats

2000  format(1x,4e20.8)

      end
