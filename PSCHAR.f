c
c   obtain character strings for restart files names

      subroutine pschar

      implicit      none

      include      'param.h'
      include      'zchl.h'
      include      'zncl.h'

      character*1   lab1
      character*2   lab2
      character*3   lab3
      character*4   lab4	  
      integer       i,k1

      do i=1,mxrsf
	  
	   write(*,*) i
        if (i.lt.10) then
          write(lab1,'(i1)')i
          k1=1
          chl(i) = lab1
          ncl(i) = k1
		  write(*,*) lab1
        else
        if (i.lt.100) then		
          write(lab2,'(i2)')i
          k1=2
		  chl(i) = lab2
          ncl(i) = k1
		  write(*,*) lab2
		else
        if (i.lt.1000) then		
          write(lab3,'(i3)')i
          k1=3
		  chl(i) = lab3
          ncl(i) = k1
		  write(*,*) lab3	
		else
        if (i.lt.10000) then		
          write(lab4,'(i4)')i
          k1=4
		  chl(i) = lab4
          ncl(i) = k1
		  write(*,*) lab4
		endif  		  	  
		endif
		endif 
        endif

      end do ! i

      end
