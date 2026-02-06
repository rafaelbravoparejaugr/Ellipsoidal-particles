      integer function indf(i1c,i2c)

      implicit     none

      include     'param.h'
      include     'zk.h'
      include     'zk4.h'

      integer      j,i1c,i2c

      do j = k4(1,i1c)+1, k4(1,i1c)+k4(2,i1c)-1
        indf = j
        if (k(j).eq.i2c) return
      end do ! j

      end
