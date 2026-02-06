
c     rem ******************** 26000
c     rem equation of punishment
c     rem ******************** called by punish

      subroutine pun1(i,n2,n2p,i1,i2,s1,s2,s5,s6)

c     inputs -- i,i1,i2,l1,l2,l3,d(),f0()
c     create -- b1,s1,s2,t(),s()
c     output -- s1,s2,s()

      implicit  none

      include  'param.h'
      include  'zf0.h'
      include  'znp0.h'
      include  'zrd.h'
      include  'zs.h'
      include  'zt.h'
      include  'zxbd.h'

      integer   i,n2,n2p,i1,i2
      real*8    s1,s2,s5,s6

      integer   i0,ix, j
      real*8    x,y, x1,x2,x3, y1,y2,y3, x21,x32,x13 
      real*8    y12,y23,y31
      real*8    b1,b2, s3,s4, oms3,om2s3
      real*8    xm,ym,rm1,rm2,thetam,phim,xs,ys,rs1,rs2
      real*8    thetas,phis,xcs,ycs


c     mi: maximum number of invasions

c     i       : Invasion number
c     i1,i2   : Number of contacted balls/boundary
c     rd(,3,2): Coordinates of 3 contact-involved points
c     f0(i)   : Lock position of invasion i at end of previous iteration
c     t()     : Displacement function matrix of (x,y) in disk i0
c     s()     : Punish geometry terms

c     call timer('pun1',22,1)


c Cálculo de x1,y1,x2,y2,x3,y3 según ecuación de la elipse

      if(i.le.n2c) then

! Contacto Elipse-Elipse

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
	
	cn=1


	else

! Contacto Elipse-Frontera

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
	
	cn=0
		
	endif !if (i.le.n2c)
	
	
	rd(1,1,1)    = rdd(1,1,i)
	rd(1,2,1)    = rdd(1,2,i)
	rd(1,3,1)    = rdd(1,3,i)
	rd(1,4,1)    = rdd(1,4,i)
	rd(1,5,1)    = rdd(1,5,i)


	rd(2,1,1)    = rdd(2,1,i)
	rd(2,2,1)    = rdd(2,2,i)
	rd(2,3,1)    = rdd(2,3,i)
	rd(2,4,1)    = rdd(2,4,i)	
	rd(2,5,1)    = rdd(2,5,i)

      y23 = y2 - y3
      y12 = y1 - y2
      y31 = y3 - y1

      x32 = x3 - x2
      x21 = x2 - x1
      x13 = x1 - x3

      b2 = 1.0d0/(x32**2 + y23**2)

      b1 = sqrt(b2)


c     Compute penetration distance at start of current step, s1 & s2

c      s1 = (x21*y31 - y12*x13)*b1
	 
	  write(*,*) 'Pun1.f',xm-x1,ym-y1,xm,ym
	  write(*,*) 'Pun1.f',xs-0.5*(x2+x3),ys-0.5*(y2+y3),xs,ys


c		pause
	
c	  write(*,*) 'Pun1.f',x1,y1,x2,y2,x3,y3

c	  pause
	  	   
	  call SEARCH()
		  				  		
	  if (cn.eq.0) then							    
	  s1 = pen(1,1)*snc*0.5d0  	  
	  else
	  s1 = pen(1,1)*snc	  
      endif

	 write(*,*) 's1',s1
	  
      s3 =  f0(i)
      s4 =  0.0d0

      oms3  = 1.0d0 - s3
      om2s3 = 1.0d0 - 2.0d0*s3

      s2 = ((x1 - oms3*x2 - s3*x3)*x32 - (y1 - oms3*y2 - s3*y3)*y23)*b1

      if (s2.ne.0.0d0) then
        s4 = (-x21*x32 - y12*y23)*b2
        if (s3.eq.s4) s2 = 0.0d0
      endif

c     Compute punishment geometry terms, s(), according to points l1,l2,l3
c         for disk i1 with contact vertex l1, and s5 & s6

      s5 = 0.0d0
      s6 = 0.0d0

      if (i.le.n2  .or.(i.gt.n2 .and.i2.gt.0)) i0 = i1
      if (i.gt.n2 .and. i.le.n2p.and.i2.lt.0)  i0 = np0(-i2)
      if (i.gt.n2p.and.i2.lt.0) goto 211

      x = x1
      y = y1

      call displac (i0,x,y)

      do j = 1, 3
        s(j)   = (y23*t(1,j) + x32*t(2,j))*b1
        s(j+6) = (x32*t(1,j) - y23*t(2,j))*b1
      end do ! j

c     For disk i2 with edge vertices l2,l3

211   if (i.gt.n2p.and.i2.gt.0) then

        s5 = (y31*xbd(1,i2  ) + x13*xbd(2,i2)
     &      + y12*xbd(1,i2+1) + x21*xbd(2,i2+1))*b1
	 
        s6 = ((-x1 + 2.d0*oms3*x2 - om2s3*x3)*xbd(1,i2) + 
     &        (-y1 + 2.d0*oms3*y2 - om2s3*y3)*xbd(2,i2) + 
     &        (x1 - om2s3*x2 - 2.d0*s3*x3)*xbd(1,i2+1) + 
     &        (y1 - om2s3*y2 - 2.d0*s3*y3)*xbd(2,i2+1))*b1

        goto 511

      endif

      if (i.le.n2) then
        i0 = i2
      else
        if (i2.lt.0) i0 = i1
        if (i.le.n2p.and.i2.gt.0) i0 = np0(i2)
      endif

      x = x2
      y = y2

      ix = i0

      call displac (i0,x,y)

      do j = 1, 3
        s(j+3) = y31*t(1,j) + x13*t(2,j)
        s(j+9) = (-x1 + 2.d0*oms3*x2 - om2s3*x3)*t(1,j)
     &         + (-y1 + 2.d0*oms3*y2 - om2s3*y3)*t(2,j)
        s(j+12) = oms3*(x32*t(1,j) - y23*t(2,j))
      end do ! j
	  

      x = x3
      y = y3

      call displac (i0,x,y)

      do j = 1, 3
        s(j+3)  = (s(j+3)  +  y12*t(1,j) + x21*t(2,j))*b1
        s(j+9)  = (s(j+9)  + (x1 - om2s3*x2 - 2.d0*s3*x3)*t(1,j)
     &                     + (y1 - om2s3*y2 - 2.d0*s3*y3)*t(2,j))*b1
        s(j+12) = (s(j+12) + s3*(x32*t(1,j) - y23*t(2,j)))*b1
      end do ! j
	 

      if (i.gt.n2p.and.i2.lt.0) then

        s5 = (y23*xbd(1,-i2) + x32*xbd(2,-i2))*b1
        s6 = (x32*xbd(1,-i2) - y23*xbd(2,-i2))*b1

      endif

511   continue


c	  write(*,*) 'pxy1',s(1),s(2),s(3)
c	  write(*,*) 'pxy2',s(4),s(5),s(6)
	  	

c     call timer('pun1',22,2)
      end
