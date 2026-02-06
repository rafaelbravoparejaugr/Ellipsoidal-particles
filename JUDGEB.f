
c     rem ******************** 
c     rem invasion judge before iteration
c     rem ******************** called by main

      subroutine judgeb (i9,n2,n2b,d0)

c     inputs -- i9,n2,m(),d(),d0
c     create -- s(),m0(,2)(when i9=1),f0()
c     output -- m0(,2)(when i9=1),f0()

      implicit     none

      integer      i9,n2,n2b
      real*8       d0

      integer      i,i1,i2,ic, j1
      real*8       a1,a2,x1,y1,x2,y2,x3,y3, bdr,psi
      real*8	   xm,ym,rm1,rm2,thetam,phim,xs,ys
      real*8       rs1,rs2,thetas,phis,xcs,ycs

      include     'param.h'
      include     'zbd.h'
      include     'zd.h'
      include     'zdm.h'
      include     'zf0.h'
      include     'zk1.h'
      include     'zm0.h'
      include     'zo.h'
      include     'zrd.h'
      include     'zs.h'
      include     'zxyt23.h'

c     md: maximum number of disks
c     mi: maximum number of invasions

c     i9: current time step number
c     n2: number of d-d contacts at start of current step
c     n2b: n2 + number of d-b contacts at start of current step
c     d0: criterion of possible contact distance of every two disks
c     m0(,2): contact index for i9=1 - 0:open;2:e-a close;1:a-a 1st close
c             3:a-a 2nd close
c     dm(,1-5): contact ball number/boundary no.,dips,contact friction angle
c               at end and start of previous step
c     d(): coordinates of disk center, their radia/friction angle/accumulated
c          rotation/radia*
c     rd(,3,2): coordinates of 3 contact-involved points
c     s(1-6): temporary use to compute s(1),s(4) - contact distances for
c             1st,2nd contacts; s(2),s(5) - lock positions for 1st,2nd contacts
c     f0(): saved as s(2), lock position, for 1st contact

c     Set 1 as initial values of m0(,2) = 1 only for i9 = 1

      call timer('judgeb',14,1)
      if (i9.eq.1) then
        do i = 1, n2b
          m0(2,i) = 1
        end do ! i
      endif

      do i = 1, n2b
        i1 = dm(1,i)
        j1 = dm(2,i)
		
        if (i.le.n2) then

! Contacto Elipse-Elipse

!rd son las coordenadas que se lleva la rutina SEARCH
!rdd son las coordenadas guardadas que se utilizan más adelante para cada contacto
		
          rd(1,1,1) = d(1,i1) 
          rd(1,2,1) = d(2,i1) 
          rd(1,3,1) = d(3,i1) 
          rd(1,4,1) = d(4,i1) 		  
          rd(1,5,1) = d(5,i1) 		  

          rdd(1,1,i) = rd(1,1,1) 
          rdd(1,2,i) = rd(1,2,1) 
          rdd(1,3,i) = rd(1,3,1)
          rdd(1,4,i) = rd(1,4,1) 		  
          rdd(1,5,i) = rd(1,5,1) 	

		  
          rd(2,1,1) = d(1,j1) 
          rd(2,2,1) = d(2,j1) 
          rd(2,3,1) = d(3,j1) 
          rd(2,4,1) = d(4,j1) 		  
          rd(2,5,1) = d(5,j1) 

          rdd(2,1,i) = rd(2,1,1) 
          rdd(2,2,i) = rd(2,2,1) 
          rdd(2,3,i) = rd(2,3,1)
          rdd(2,4,i) = rd(2,4,1) 		  
          rdd(2,5,i) = rd(2,5,1) 
		  
      	  xm     = rdd(1,1,i)
	  ym     = rdd(1,2,i)
          rm1    = rdd(1,3,i)
	  rm2    = rdd(1,4,i)
	  thetam = rdd(1,5,i)
	  phim   = 0.0d0 

          xs     = rdd(2,1,i)
	  ys     = rdd(2,2,i)
          rs1    = rdd(2,3,i)
	  rs2    = rdd(2,4,i)
	  thetas = rdd(2,5,i)
	  phis   = 0.0d0

		  
          cn = 1
		  
		  
        else
          if (j1.gt.0) then
	
! Contacto Elipse-Frontera		  
			  	  
            rd(1,1,1) = d(1,i1) 
            rd(1,2,1) = d(2,i1) 
            rd(1,3,1) = d(3,i1) 
            rd(1,4,1) = d(4,i1) 			
            rd(1,5,1) = d(5,i1) 			
			
            rdd(1,1,i) = rd(1,1,1) 
            rdd(1,2,i) = rd(1,2,1) 
            rdd(1,3,i) = rd(1,3,1)
            rdd(1,4,i) = rd(1,4,1) 		  
            rdd(1,5,i) = rd(1,5,1) 			
			
			
            rd(2,1,1) = bd(1,j1)
            rd(2,2,1) = bd(2,j1)
            rd(2,3,1) = bd(1,j1+1)
            rd(2,4,1) = bd(2,j1+1)
            rd(2,5,1) = 0.0d0	
			
            rdd(2,1,i) = rd(2,1,1) 
            rdd(2,2,i) = rd(2,2,1) 
            rdd(2,3,i) = rd(2,3,1)
            rdd(2,4,i) = rd(2,4,1) 		  
            rdd(2,5,i) = 0.0d0	
			
	    xm     = rdd(1,1,i)
	    ym     = rdd(1,2,i)
	    rm1    = rdd(1,3,i)
	    rm2    = rdd(1,4,i)
	    thetam = rdd(1,5,i)
	    phim   = 0.0d0
																	
            cn = 0			

          else
          endif
        endif

        call SEARCH()

        s(1) = snc*pen(1,1)*0.5d0

	rdd(1,6,i)=pen(3,1)	
	
        phim=pen(3,1)
			
	if(cn.eq.1) then

        s(1) = snc*pen(1,1)		
	
        rdd(2,6,i)=pen(4,1)
	
        phis=pen(4,1)		
	
        endif
		
		
C****************************************************************************************************
        if (i.le.n2) then

! Contacto Elipse-Elipse

	x1 = (rm1*cos(phim)*cos(thetam) - rm2*sin(phim)*sin(thetam)) + xm
	y1 = (rm1*cos(phim)*sin(thetam) + rm2*sin(phim)*cos(thetam)) + ym

	xcs = (rs1*cos(phis)*cos(thetas) - rs2*sin(phis)*sin(thetas)) + xs
	ycs = (rs1*cos(phis)*sin(thetas) + rs2*sin(phis)*cos(thetas)) + ys
	
	
	x2 = xcs - rs1*sin(phis)*cos(thetas) - rs2*cos(phis)*sin(thetas) 
	y2 = ycs - rs1*sin(phis)*sin(thetas) + rs2*cos(phis)*cos(thetas) 
	x3 = xcs + rs1*sin(phis)*cos(thetas) + rs2*cos(phis)*sin(thetas) 
	y3 = ycs + rs1*sin(phis)*sin(thetas) - rs2*cos(phis)*cos(thetas) 	


        i1 = dm(1,i)
        i2 = dm(2,i)

        if (((i1.eq.2).and.(i2.eq.3)).or.((i1.eq.2).and.(i2.eq.3))) then	
        xt23 = xcs 
        yt23 = ycs
        endif	  
		  
		  
        cn = 1
		  
		  
        else

        if (j1.gt.0) then
	
! Contacto Elipse-Frontera		  
			  	  			
	xm     = rdd(1,1,i)
	ym     = rdd(1,2,i)
	rm1    = rdd(1,3,i)
	rm2    = rdd(1,4,i)
	thetam = rdd(1,5,i)

	x2     = rdd(2,1,i)
	y2     = rdd(2,2,i)
	x3     = rdd(2,3,i)
	y3     = rdd(2,4,i)
	

	x1 = (rm1*cos(phim)*cos(thetam) - rm2*sin(phim)*sin(thetam)) + xm
	y1 = (rm1*cos(phim)*sin(thetam) + rm2*sin(phim)*cos(thetam)) + ym																			
        cn = 0			

        else
        endif
        endif		
		
C****************************************************************************************************


		
c	write(*,*) 'Judgeb.f',rd(1,1,1),rd(1,2,1),rd(2,1,1),rd(2,2,1)
		 
		
        a2 = (x3 - x2)**2.0d0 + (y3 - y2)**2.0d0
        a1 = sqrt(a2)
		
        s(2) = ((x1 - x2)*(x3 - x2) + (y1 - y2)*(y3 - y2))/a2


c       write(*,*) 'JUDGEB',s(1),s(2),x1,y1,x2,y2,x3,y3	
c	pause

        if (s(1).lt.0.000001*d0.or.i9.gt.1) goto 210

c       Set this m0(,2) = 0 only if i9 = 1 and s(1)>.000001*d0

        m0(2,i) = 0

c       Save lock position parameter into f0() at start of current step
c       for i9>1 and m0(,2)previous = 2(locked), keep previous f0()

c       Corrections for clusters (4/8/96)

210     continue

c       Special concern - clusters with 2 disks contacting side

        if (i.gt.n2.and.k1(2,i1).ne.1) then
           if (k1(1,i1).eq.k1(1,int(dm(1,i-1))).or.
     &        k1(1,i1).eq.k1(1,int(dm(1,i+1)))) then
              if (j1.eq.dm(2,i-1).or.j1.eq.dm(2,i+1)) then
                 ic = 1
              endif
           endif
        endif

        if (.not.(abs(m0(2,i)).eq.2 .and. i9.gt.1)) then
        f0(i) = s(2)
        endif

c       Set other values of m0(,2) = 2 only if i9 = 1

        if (m0(2,i).ne.0.and.i9.eq.1) then
          m0(2,i) = 2
        endif
		
	write(*,*) 'JUDGE B',s(1),s(2),m0(1,i),m0(2,i)
c		pause
		
      end do ! i

c     Check friction angle for boundary vertex -> disk contact for i9 > 1

      if (i9.ne.1) then
        do i = n2+1, n2b
          i2 = dm(2,i)
          if (i2.le.0) then
            i1 = dm(1,i)
            if (o(2,i).ge.0.d0) then
              dm(5,i) = min(d(6,i1),bd(5,-i2))
              dm(6,i) = min(d(7,i1),bd(6,-i2))
            else
              dm(5,i) = min(d(6,i1),bd(5,-i2+1))
              dm(6,i) = min(d(7,i1),bd(6,-i2+1))
            endif
          endif
        end do ! i
      endif
	  

      call timer('judgeb',14,2)
      end
