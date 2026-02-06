
c     rem ******************** 13000
c     rem volume force
c     rem ******************** called by main

      subroutine volumn

c     inputs -- o1,o2,n0,k1(),h0(,1),f()
c     create -- f()
c     output -- f()

      implicit     none

      include     'param.h'
      include     'zf.h'
      include     'zh0.h'
      include     'zn.h'
      include     'zp.h'

      integer      i

c     md: maximum number of disks

c     o1: unit volume force in x dirction
c     o2: unit volume force in y dirction
c     n0: number of disks
c     k1(): permutation matrix -  k1(old dk no)=new dk no
c     f(): free terms before step iterations
c     h0(,1): ball area

      call timer('volumn',32,1)

      do i=1, n0
        f(1,i)=f(1,i)+o1*h0(1,i)
        f(2,i)=f(2,i)+o2*h0(1,i)
      end do ! i

      call timer('volumn',32,2)

      end
