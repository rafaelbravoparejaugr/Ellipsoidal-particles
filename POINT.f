
c     rem ******************** 600x
c     rem determination of point loadings at the start of present step
c     rem   using linear interpolation along the time scale
c     rem ******************** called by main

      subroutine point (q02)

      implicit     none

      real*8       q02

      integer      i,i1,i2, j
      real*8       a1,dnew

      include     'param.h'
      include     'zlp.h'
      include     'zn.h'
      include     'zpl.h'
      include     'zsp.h'

      call timer('point',20,1)

      do i = 1, ip(7,2)

        if (i.gt.ip(3,2).and.i.le.ip(5,2)) then
          i2 = mod(i-ip(3,2),2)
          if (i2.eq.0) goto 100
        endif

        do j = lp(1,i), lp(2,i)
          i1 = j
          if (pl(1,j).le.q02.and.pl(1,j+1).gt.q02) goto 201
        end do ! j

201     a1 = (q02 - pl(1,i1))/(pl(1,i1+1) - pl(1,i1))
        if (q02.ge.pl(1,lp(2,i))) a1 = 0.d0

        if (i.le.ip(1,2) .or.
     &     (i.gt.ip(3,2) .and. i.le.ip(4,2)) .or. i.gt.ip(6,2)) then

          sp(3,i) = pl(2,i1) + a1*(pl(2,i1+1) - pl(2,i1))
          sp(4,i) = pl(3,i1) + a1*(pl(3,i1+1) - pl(3,i1))

        else

          dnew = pl(2,i1) + a1*(pl(2,i1+1) - pl(2,i1))

          if (q02.eq.0.0d0) then
            sp(3,i) = dnew
          else
            sp(3,i) = sp(3,i) + dnew - sp(4,i)
          endif

          sp(4,i) = dnew

        endif

100   continue

      end do ! i

      call timer('point',20,2)
      end
