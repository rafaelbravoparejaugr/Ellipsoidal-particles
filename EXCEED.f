
c     rem ********************
c     rem warning when declared dimension is exceeded
c     rem ********************

      subroutine exceed (i,j,s1,s2)

c     input  -- i,j,s1,s2
c     create -- none
c     output -- none

      implicit  none

      character s1*3,s2*15
      integer   i,j

      if (i.gt.j) then
        write (*,10) s1,s2
        write ( *,*) s1,' <',i
        write (2,10) s1,s2
        stop
      endif

c     Format

10    format (' Warning! ',a2,' is exceeded ',a15)

      end
