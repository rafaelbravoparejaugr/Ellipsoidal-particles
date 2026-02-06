
c     rem statics control

      subroutine static (ist,e1,ev2,i9,rmax,q02)

      implicit         none

      include         'param.h'
      include         'zh.h'
      include         'zh0.h'
      include         'zn.h'
      include         'zp.h'
      include         'zstatic.h'

      integer          ist,i9
      real*8           e1,ev2,rmax,q02

      integer          i
      real*8           ev2p,grad,tkn,tma,tv2

c     md: maximum number of disks

c     e1  : statics criterion, when tkn/tma < e1, set ist=-1 -> stop running
c     tkn : sum of kinetic energy, first excluding angular contribution
c     tma : sum of disk mass
c     q02p: previous time
c     grad: gradient of avg. kinetic energy per mass

      call timer('static',29,1)
      ev2p = ev2
      tkn  = 0.d0
      tma  = 0.d0
      grad = 0.d0
      do i = 1, n0
        tv2 = h(1,i)*h(1,i) + h(2,i)*h(2,i)
        if (i.le.nd0) then
          tkn = tkn + .5d0*h0(2,i)*tv2
          tma = tma + h0(2,i)
        else
          tkn = tkn + .5d0*h0(1,i)*o0(i)*tv2
          tma = tma + h0(1,i)*o0(i)
        endif
      end do ! i

      ev2 = tkn/tma/rmax/rmax
      if (i9.ne.1) then
         grad = (ev2p-ev2)/(q02-q02p)
      endif
      q02p = q02

      if (ist.eq.1) then
        if (mod(i9,100).eq.0) then
          write (8,2000) q02,ev2,grad
          write(*,*) q02,ev2,grad
        endif

        if (i9.gt.1.and.grad.gt.0.d0 .and. grad.lt.e1) ist = -1
        if (ist.eq.-1) write (8,*) 'tkn,tma,ev2,grad:',tkn,tma,ev2,grad

      endif
      call timer('static',29,2)

2000  format (3e12.5)

      end
