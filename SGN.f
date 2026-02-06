
      real*8 function sgn(x)

      implicit none
      real*8   x

      if (x.gt.0.d0) then
        sgn =  1.d0
      else if (x.eq.0.d0) then
        sgn =  0.d0
      else
        sgn = -1.d0
      endif

      end
