
c     rem ******************** integ
c     compute orientations, lengths of boundary segments
c       s,sx,sy of polygons?
c     rem ******************** called by main

      subroutine integ (i9,pi)

      implicit      none

      include      'param.h'
      include      'zbd.h'
      include      'zh0.h'
      include      'zk0b.h'
      include      'zn.h'

      integer       i9
      real*8        pi

      integer       i,i1,i2, j,ju
      real*8        e1,e2,e4, x2,y2,x3,y3, f1

c     Update s,sx,sy of rigid polygons

      call timer('integ',12,1)
      do i = nd0+1, n0
        i1 = k0b(1,i-nd0)
        i2 = bd(1,i1)
        h0(1,i) = 0.d0
        h0(2,i) = 0.d0
        h0(3,i) = 0.d0
        h0(4,i) = 0.d0
        h0(5,i) = 0.d0
        do j = i1+1, i1+i2
          x2      = bd(1,j)
          y2      = bd(2,j)
          x3      = bd(1,j+1)
          y3      = bd(2,j+1)
          f1      = (x2*y3 - x3*y2)
          h0(1,i) = h0(1,i) - f1*0.5d0
          h0(2,i) = h0(2,i) - f1*(x2+x3)/6.d0
          h0(3,i) = h0(3,i) - f1*(y2+y3)/6.d0
          h0(4,i) = h0(4,i) - f1*(x2*x2+x3*x3+x2*x3)/12.d0
          h0(5,i) = h0(5,i) - f1*(y2*y2+y3*y3+y2*y3)/12.d0
        end do ! j
      end do ! i

c     Update length and orientation of rigid boundaries and polygons

      do i = 1, np+nb
        if (i9.le.0.or.k0b(2,i).ne.1) then
          if (i.le.np) then
            ju = k0b(1,i) + bd(1,k0b(1,i))
          else
            ju = k0b(1,i) + bd(1,k0b(1,i)) - 1
          endif
          do j = k0b(1,i)+1, ju
            e1 = bd(1,j+1) - bd(1,j)
            e2 = bd(2,j+1) - bd(2,j)
            if (i9.eq.0) bd(3,j) = sqrt(e1*e1+e2*e2)

            if (e1.eq.0.d0) then
              if(e2.gt.0.d0) then
                e4 =  pi*0.5d0
              elseif(e2.lt.0.d0) then
                e4 = -pi*0.5d0
              endif
            else
              e4 = atan(e2/abs(e1))
              if (e1.lt.0.d0) e4 = pi - e4
              if (e4.lt.0.d0) e4 = 2.d0*pi + e4
            endif

            bd(4,j) = e4
          end do ! j
        endif
      end do ! i

      call timer('integ',12,2)
      end
