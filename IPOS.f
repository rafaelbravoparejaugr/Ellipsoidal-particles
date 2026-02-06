c$Id:$
      integer function ipos(file,nn)
c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Locate last character in character array
 
c      Inputs:
c         file(*) - Array to search
c         nn      - Length of array
 
c      Outputs:
c         ipos   - Position of last character
c-----[--+---------+---------+---------+---------+---------+---------+-]
 
      implicit none
 
      integer n,nn
      character*1 file(nn)
 
      save
 
      do n = nn,1,-1
        if(file(n).ne.' ') go to 100
      end do
      n    = 0
 
100   ipos = n
 
      end
