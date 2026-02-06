c$Id:$
      subroutine addext(fnam,fext,ifnam,ifext)
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose:       Add extender to file name for disk i/o
 
c      Inputs:
c        fname(*)   - File name without extender
c        fext(*)    - Extender
c        ifnam      - Length of 'fname'
c        ifext      - Length of 'fext'
 
c      Outputs:
c        fname(*)   - File name with added  '.' and 'fext'
c-----[--+---------+---------+---------+---------+---------+---------+-]
 
      implicit none
 
      integer iposl, iposx, ifnam,ifext
      character fnam*(*),fext*(*)
 
      integer ipos
      call timer('addext',1,1)
 
      iposl = ipos(fnam,ifnam) + 1
      iposx = ipos(fext,ifext)
 
      fnam(iposl:ifnam)         = '. '
      fnam(iposl+1:iposl+iposx) = fext(1:iposx)
 
      call timer('addext',1,2)
      end
