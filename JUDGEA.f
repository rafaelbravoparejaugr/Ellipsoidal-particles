
c     rem ******************** 
c     rem invasion judge after iteration
c     rem ******************** called by main

      subroutine judgea (n2,n2p,n2b,g0,s0,dd,m9)

c     inputs -- n2,n2b,m0(,2),dm(),d(),k1(),f(),g0,s0,w(),f0(),o(,2),
c               phi(),dd
c     create -- m0(,1),t(),q(),p(,1-2),d1,b1,p(,3),o(,1),m0(,2),f0(),
c                 o(,2),i7,i8,m9
c     output -- m0(,1-2),o(,1-2),f0(),m9

      implicit     none

      include     'param.h'

      integer      n2,n2p,n2b,m9
      real*8       g0,s0,dd

      integer      i,i0,i1,i2,i7,i8, j,j1,j2, j1a,j2a, l,l1
      real*8       a1,a2,a3,b1,b2,d1,q1,s1,t1,g00
      real*8       x0,y0, x1,y1,x2,y2,x3,y3, x01,y01,x02,y02,x03,y03

      real*8       p(3,3)

      real*8       sgn

      real*8	   xm,ym,rm1,rm2,thetam,phim,xs,ys,rs1,rs2
      real*8       thetas,phis,xcs,ycs

      include      'zb.h'
      include      'zbd.h'
      include      'zd.h'
      include      'zdm.h'
      include      'zf0.h'
      include      'zk1.h'
      include      'zm0.h'
      include      'zn.h'	  
      include      'znp0.h'
      include      'zo.h'
      include      'zq.h'
      include      'zrd.h'
      include      'zt.h'
      include      'zw.h'
      include      'zxbd.h'


c     md: maximum number of disks
c     mi: maximum number of invasions

c     n2: number of d-d contacts at start of current step
c     n2b: n2 + number of d-p + d-b contacts at start of current step
c     i7: changes in open->close during current iteration
c     i*8: changes in close->open during current iteration
c     n9: total iteration number from beginning
c     m9: current iteration number of current step (-1 for iteration funished)
c     g0: penalty (normal contact stiffness)
c     s0: closed-opened criterion
c     dd: conversion factor between radian and degree
c     m0(,1-2): contact indices at start and after current iteration
c     dm(,1-5): contact ball number/boundary no.,dips,contact friction angle
c               at end and start of previous step
c     rd(,3,2): coordinates of 3 contact-involved points
c     d(): coordinates of disk center, their radia/friction angle/accumulated
c          rotation/radia*
c     t(): displacement function matrix of point (x,y) in disk i0
c     q(): temporary use saved as coordinates of contact-involved vertices
c     p(): temporary use - p(,1-2) saved as displacements of 2 contact-involved
c          points
c     k1(): permutation matrix -  k1(old dk no)=new dk no
c     f(): displacement solutions
c     o(,1-2): distance of penetration, sliding distance w.r.t reference line
c              of each invasion at start of current iteration
c     w(): window limits
c     f0(i): lock position of invasion i

c  ****8/13/96 correction: b() expecting cluster number not disk #
c  ****8/14/96 made j2a and j1a cluster # of disks i1,i2 because
c             displac no longer sending back cluster #

c     Transfer open-close record before iteration to m0(,1)

      call timer('judgea',13,1)
      do i = 1, n2b
        m0(1,i) = m0(2,i)
      end do ! i

c     Invasion judge after iteration

      do i = 1, n2b

! Se calcula la penetración después de la iteración

        if (i.le.n2c) then
		
        i1 = dm(1,i)
        i2 = dm(2,i)


         rd(1,1,1) = rdd(1,1,i)+b(1,k1(1,i1))	 
         rd(1,2,1) = rdd(1,2,i)+b(2,k1(1,i1)) 
         rd(1,3,1) = rdd(1,3,i) 
         rd(1,4,1) = rdd(1,4,i) 		  
         rd(1,5,1) = rdd(1,5,i)+b(3,k1(1,i1))  	

         rd(2,1,1) = rdd(2,1,i)+b(1,k1(1,i2))	
         rd(2,2,1) = rdd(2,2,i)+b(2,k1(1,i2))	 
         rd(2,3,1) = rdd(2,3,i) 
         rd(2,4,1) = rdd(2,4,i) 		  
         rd(2,5,1) = rdd(2,5,i)+b(3,k1(1,i2))	
		 
	 cn = 1 		

		
		else

        i1 = dm(1,i)
        i2 = dm(2,i)
		
         rd(1,1,1) = rdd(1,1,i)+b(1,k1(1,i1))	
         rd(1,2,1) = rdd(1,2,i)+b(2,k1(1,i1))	 
         rd(1,3,1) = rdd(1,3,i) 
         rd(1,4,1) = rdd(1,4,i) 		  
         rd(1,5,1) = rdd(1,5,i)+b(3,k1(1,i1))	  	

         rd(2,1,1) = rdd(2,1,i)
         rd(2,2,1) = rdd(2,2,i) 
         rd(2,3,1) = rdd(2,3,i) 
         rd(2,4,1) = rdd(2,4,i) 		  
         rd(2,5,1) = rdd(2,5,i)
		 
	 cn = 0
		
	 endif ! if (i.le.n2c) then

C***************************************************************************************

        if (i.le.n2c) then
				 
          xm     = rdd(1,1,i)
	  ym     = rdd(1,2,i)
          rm1    = rdd(1,3,i)
	  rm2    = rdd(1,4,i)
	  thetam = rdd(1,5,i)
	  phim   = rdd(1,6,i) 

          xs     = rdd(2,1,i)
	  ys     = rdd(2,2,i)
          rs1    = rdd(2,3,i)
	  rs2    = rdd(2,4,i)
	  thetas = rdd(2,5,i)
	  phis   = rdd(2,6,i)

	x1 = (rm1*cos(phim)*cos(thetam) - rm2*sin(phim)*sin(thetam)) + xm
	y1 = (rm1*cos(phim)*sin(thetam) + rm2*sin(phim)*cos(thetam)) + ym

	xcs = (rs1*cos(phis)*cos(thetas) - rs2*sin(phis)*sin(thetas)) + xs
	ycs = (rs1*cos(phis)*sin(thetas) + rs2*sin(phis)*cos(thetas)) + ys
	
	
	x2 = xcs - rs1*sin(phis)*cos(thetas) - rs2*cos(phis)*sin(thetas) 
	y2 = ycs - rs1*sin(phis)*sin(thetas) + rs2*cos(phis)*cos(thetas) 
	x3 = xcs + rs1*sin(phis)*cos(thetas) + rs2*cos(phis)*sin(thetas) 
	y3 = ycs + rs1*sin(phis)*sin(thetas) - rs2*cos(phis)*cos(thetas)	

		 
        cn = 1
		
		else
		 
	xm     = rdd(1,1,i)
	ym     = rdd(1,2,i)
	rm1    = rdd(1,3,i)
	rm2    = rdd(1,4,i)
	thetam = rdd(1,5,i)
	phim   = rdd(1,6,i) 

	x2     = rdd(2,1,i)
	y2     = rdd(2,2,i)
	x3     = rdd(2,3,i)
	y3     = rdd(2,4,i)
	

	x1 = (rm1*cos(phim)*cos(thetam) - rm2*sin(phim)*sin(thetam)) + xm
	y1 = (rm1*cos(phim)*sin(thetam) + rm2*sin(phim)*cos(thetam)) + ym		 
		 
        cn = 0
		
	endif ! if (i.le.n2c) then

C***************************************************************************************
		  
		  
C***************************************************************************************
        do j = 1, 3
          do j1 = 1, 2
            p(j,j1) = 0.
          end do ! j1
        end do ! j
		
c       write(*,*) 'judgeA',x1,y1,0.5d0*(x2+x3),0.5d0*(y2+y3)		
c       write(*,*) 'judgeAcentros',rd(1,1,1),rd(1,2,1),rd(2,1,1),rd(2,2,1)

        if (i.le.n2p) then

          j1 = i1
          j2 = i2

          if (i.gt.n2) then

            if (i2.gt.0) then

c             'disk-poly side contact'

              j1a = k1(1,j1)
              j2  = np0(i2)
              j2a = j2 - nd + nd0
            else

c             'disk-poly vertex contact'

              j1  = np0(-i2)
              j1a = j1 - nd + nd0
              j2  = i1
              j2a = k1(1,j2)
            endif
          else

            j1a = k1(1,j1)
            j2a = k1(1,j2)
          endif

          call displac (j1,x1,y1)

          do l = 1, 2
            do l1 = 1, 3
              p(1,l) = p(1,l) + t(l,l1)*b(l1,j1a)
            end do ! l1
          end do ! l

          call displac (j2,x2,y2)

          do l = 1, 2
            do l1 = 1, 3
              p(2,l) = p(2,l) + t(l,l1)*b(l1,j2a)
            end do ! l1
          end do ! l

          call displac (j2,x3,y3)

          do l = 1, 2
            do l1 = 1, 3
              p(3,l) = p(3,l) + t(l,l1)*b(l1,j2a)
            end do ! l1
          end do ! l

        else

          i0 = i1
          if (i2.gt.0) then

            call displac (i0,x1,y1)

            do l = 1, 2
              do l1 = 1, 3
                p(1,l) = p(1,l) + t(l,l1)*b(l1,k1(1,i0))
              end do ! l1
            end do ! l

            do l = 1, 2
              p(2,l) = xbd(l,i2)
              p(3,l) = xbd(l,i2+1)
            end do ! l

          else

            call displac (i0,x2,y2)

            do l = 1, 2
              do l1 = 1, 3
                p(2,l) = p(2,l) + t(l,l1)*b(l1,k1(1,i0))
              end do ! l1
            end do ! l

            call displac (i0,x3,y3)

            do  l = 1, 2
              do l1 = 1, 3
                p(3,l) = p(3,l) + t(l,l1)*b(l1,k1(1,i0))
              end do ! l1
            end do ! l

            do l = 1, 2
              p(1,l) = xbd(l,-i2)
            end do ! l

          endif

        endif


C**********************************************************************************************		  
		  	
        call SEARCH()
	
		
        a1 = sqrt((x3 - x2)**2 + (y3 - y2)**2)
        d1 = ((x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1))
        b1 = p(1,1)*(y2 - y3) + p(1,2)*(x3 - x2)
     &     + p(2,1)*(y3 - y1) + p(2,2)*(x1 - x3)
     &     + p(3,1)*(y1 - y2) + p(3,2)*(x2 - x1)
        p(1,3) = (d1 + b1)/a1

c       Assign m0(,2)&o(,1) after iteration according to p(,3)

        o(1,i) = pen(1,1)*snc
		
	if(cn.eq.0) then
        o(1,i) = pen(1,1)*snc*0.5d0
	endif !if(cn
		
	write(*,*) 'o(1,i) JUDGEA y m0(1',o(1,i),m0(1,i)


c       Opening judge: closed contact before iteration

        m0(2,i) = 1
        if (m0(1,i).ne.0) then
          g00 = g0
          if (m0(1,i).ne.m0(2,i)) then
            g00 = 100.d0
          endif
          if (o(1,i).ge.s0*(w(4)-w(2))/g00) then
            m0(2,i) = 0
          endif

c       Opening judge: open contact before iteration

        elseif (o(1,i).ge.-10.0d0*s0*(w(4)-w(2))/g0) then

          m0(2,i) = 0
		
        endif

        if (m0(2,i).eq.0) goto 2228

c       For cases of pre-/post-close cases, keep old lock position (prev. it)

        if (m0(1,i).eq.0) then

c         Find f0() for pre-opened -> closed contact

          if (d1.lt.0.d0.or.abs(b1).lt..000000001d0) then
            t1  =  0.0d0
          else
            t1 = -d1/b1
            if (t1.ge.0.d0) then
              t1 = min(t1,1.d0)
            else
              t1 = 0.d0
            endif
          endif

          x01 = x1 + t1*p(1,1)
          y01 = y1 + t1*p(1,2)
          x02 = x2 + t1*p(2,1)
          y02 = y2 + t1*p(2,2)
          x03 = x3 + t1*p(3,1)
          y03 = y3 + t1*p(3,2)
          b2 = (x03 - x02)**2 + (y03 - y02)**2


c         Find new lock position for pre-opened -> closed contact

          f0(i) = ((x01 - x02)*(x03 - x02) + (y01 - y02)*(y03 - y02))/b2
		  
c         f0(i) = 0.5d0

        endif

c       Find moving distance of lock point along reference line 23 (sliding)

        x0  = (1.d0 - f0(i))*x2 + f0(i)*x3
        y0  = (1.d0 - f0(i))*y2 + f0(i)*y3
        x01 = (1.d0 - f0(i))*p(2,1) + f0(i)*p(3,1)
        y01 = (1.d0 - f0(i))*p(2,2) + f0(i)*p(3,2)
        s1  = (x1 - x0)*(x3 - x2) + (y1 - y0)*(y3 - y2)
     &      + (x3 - x2)*(p(1,1) - x01) + (y3 - y2)*(p(1,2) - y01)
     &      + (x1 - x0)*(p(3,1)-p(2,1)) + (y1 - y0)*(p(3,2)-p(2,2))
        q1  = sgn(o(2,i))
        o(2,i) = s1/a1
		
	write(*,*) 'TANGENCIAL y bandera',o(2,i),m0(1,i)
		
c	pause

c       Judge it sliding or locked in the e-a contact
c         if pre-close-sliding contact moves reversely now, locked!

       if (m0(1,i).eq.1.and.o(2,i)*q1.le.-1.d-7.and.dm(5,i).gt.0.9d0)
     $     goto 2226

c       If penetration distance is numerically too small, sliding!

        if (abs(o(1,i)).lt.1.d-14) goto 2228

c        Set sliding criterion
c        Find phi again for boundary vertex -> ball contact

        if (i.gt.n2.and.i2.lt.0) then

          if (o(2,i).ge.0.d0) then
            dm(5,i) = min(d(6,i1),bd(5,-i2))
            dm(6,i) = min(d(7,i1),bd(6,-i2))
          else
            dm(5,i) = min(d(6,i1),bd(5,-i2+1))
            dm(6,i) = min(d(7,i1),bd(6,-i2+1))
          endif
        endif
		
c************************************MODIFICACION************************************

        a2 =  -o(1,i)*g0*tan(dd*dm(5,i)) + dm(6,i)
        a3 =  o(2,i)*g0


c        a2 = fn(i)*tan(dd*dm(5,i))
c        a3 = fr(i)



        if (m0(1,i).eq.1.and.dm(5,i).le.9.5d0) goto 2228

c       If friction resistance (phi()) <= driving force, sliding!

        if (abs(a3).ge.abs(a2)) then
	write(*,*) 'o(1 o(2',o(1,i),o(2,i)
	write(*,*) 'F spring',a3,'rozamiento',a2,'SLIP'
c	pause
	go to 2228
	else
	write(*,*) 'F spring',a3,'rozamiento',a2,'STICK'
c       pause
	endif

2226      m0(2,i) = 2

2228    continue

	write(*,*) 'JUDGE A',o(1,i),o(2,i),f0(i),m0(1,i),m0(2,i)
		
c		pause

      end do ! i

c     Convergence check

      i7 = 0
      i8 = 0
      do i = 1, n2b

c       open->close case?

        if (m0(1,i).eq.0 .and. m0(2,i).gt.0 ) then
          i7 = i7 + 1

c       close->open case?

        elseif (m0(2,i).eq.0 .and. m0(1,i).gt.0) then
          i8 = i8 + 1
c       e-a contact with sliding-locked change?

        elseif (m0(1,i).ne.m0(2,i)) then
          i7 = i7 + 1
        endif
      end do ! i
	  
      if (i7 + i8.le.0) then
c       m9 = -1
        m9 = min(-1,-abs(m9))
      endif
	
C Solo para verificacion	    
	  do i=1,n2b
      if (i7 + i8.le.0) then
c	  write(*,*) 'CONTACTO CONVERGENTE'
      endif	  
      end do !i

      call timer('judgea',13,2)
      end
