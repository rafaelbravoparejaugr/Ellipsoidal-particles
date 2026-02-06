	MODULE data_storage
! The data is stored in this module.
! This is equivalent to a COMMON area in old Fortran.

	IMPLICIT NONE
	INTEGER, PARAMETER, PRIVATE  :: dp = SELECTED_REAL_KIND(12, 60)

	REAL (dp), SAVE  :: x(100), y(100)

	END MODULE data_storage



	SUBROUTINE fcn(m, n, p, fvec, fjac, iflag)
	! Calculate either residuals or the Jacobian matrix.

	USE data_storage
	IMPLICIT NONE
	INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
	INTEGER, INTENT(IN)        :: m, n
	REAL (dp), INTENT(IN)      :: p(:)
	REAL (dp), INTENT(IN OUT)  :: fvec(:)
	REAL (dp), INTENT(OUT)     :: fjac(:,:)
	INTEGER, INTENT(IN OUT)    :: iflag

	! Local variables

	REAL (dp), PARAMETER  :: one = 1.0_dp
	INTEGER               :: i
	REAL (dp)             :: expntl, temp


	real *8 e1x,e1y,e2x,e2y,xs,ys,xm,ym,xcm,ycm,xcs,ycs,cx1,cy1,cx2
	real *8 cy2
	real *8 thetam,thetas,rs1,rs2,rm1,rm2,dargs,dargm,dphis,dphim,dmod
	real *8 ddp,phis,phim,cp1,cp2,argm,args,dpd,dxcm,dycm,dxcs,dycs
	real *8 xs2,ys2,dcx1,dcy1,dcx2,dcy2,xctm,yctm,xcnm,ycnm

      include     'param.h'	  
      include     'zrd.h'	


	 if (cn.eq.1) then
	 
!Contacto Elipse-Elipse
	 
			xm=     rd(1,1,1)
			ym=     rd(1,2,1)
			rm1=    rd(1,3,1)
			rm2=    rd(1,4,1)
			thetam= rd(1,5,1)
			xs=     rd(2,1,1)
			ys=     rd(2,2,1)
			rs1=    rd(2,3,1)
			rs2=    rd(2,4,1)
			thetas= rd(2,5,1)

		else

!Contacto Elipse-Frontera

			xm =     rd(1,1,1)
			ym =     rd(1,2,1)
			rm1 =    rd(1,3,1)
			rm2 =    rd(1,4,1)
			thetam = rd(1,5,1)
			xs  =    rd(2,1,1)
			ys  =    rd(2,2,1)
			xs2 =    rd(2,3,1)
			ys2 =    rd(2,4,1)
			rs1 =    rm1
			rs2 =    rm2
			thetas = 2.0d0*atan2((ys2-ys),(xs2-xs))-thetam
			  			  
! Cálculo de la normal y tangente a la frontera

		e1x = -(ys - ys2)/sqrt((xs - xs2)**2 + (ys - ys2)**2)
		e1y =  (xs - xs2)/sqrt((xs - xs2)**2 + (ys - ys2)**2)	
		e2x =  (xs - xs2)/sqrt((xs - xs2)**2 + (ys - ys2)**2)
		e2y =  (ys - ys2)/sqrt((xs - xs2)**2 + (ys - ys2)**2)
		xcs = (e1y*e2x*xm - e1x*(e2y*xs + e2x*(-ys + ym)))
     &          /(e1y*e2x - e1x*e2y)
		ycs = (e1y*(e2y*(-xs + xm) + e2x*ys) - e1x*e2y*ym)
     &          /(e1y*e2x - e1x*e2y)

		xs = 2.0d0*(xcs-xm)+xm
		ys = 2.0d0*(ycs-ym)+ym		
		
		
		endif !if (cn.le.n2c) then
			  

	IF (iflag == 1) THEN

	
	e1x =  (xs - xm)/sqrt((xs - xm)**2 + (ys - ym)**2)
	e1y =  (ys - ym)/sqrt((xs - xm)**2 + (ys - ym)**2)
	e2x = -(ys - ym)/sqrt((xs - xm)**2 + (ys - ym)**2)
	e2y =  (xs - xm)/sqrt((xs - xm)**2 + (ys - ym)**2)
	
	cx1 = cos(p(1)-thetam)*e1x + sin(p(1)-thetam)*e2x
	cy1 = cos(p(1)-thetam)*e1y + sin(p(1)-thetam)*e2y

	cx2 = -cos(p(1)-thetas)*e1x - sin(p(1)-thetas)*e2x
	cy2 = -cos(p(1)-thetas)*e1y - sin(p(1)-thetas)*e2y

	phim = atan2(rm2*cy1,rm1*cx1) 
	phis = atan2(rs2*cy2,rs1*cx2) 
	
	xcm = (rm1*cos(phim)*cos(thetam) - rm2*sin(phim)*sin(thetam)) + xm
	ycm = (rm1*cos(phim)*sin(thetam) + rm2*sin(phim)*cos(thetam)) + ym

	xcs = (rs1*cos(phis)*cos(thetas) - rs2*sin(phis)*sin(thetas)) + xs
	ycs = (rs1*cos(phis)*sin(thetas) + rs2*sin(phis)*cos(thetas)) + ys

	fvec =  sqrt((xcs - xcm)**2 + (ycs - ycm)**2)

	 pen(3,1)=phim
	 pen(4,1)=phis
	 
	xctm=-rm1*sin(phim)*cos(thetam) - rm2*cos(phim)*sin(thetam)
	yctm=-rm1*sin(phim)*sin(thetam) + rm2*cos(phim)*cos(thetam) 
	
	xcnm = -yctm/sqrt(xctm**2+yctm**2)
	ycnm =  xctm/sqrt(xctm**2+yctm**2) 

	 

	if(((xcs-xcm)*xcnm+(ycs-ycm)*ycnm).lt.0.0d0) then
	snc =  1.0d0
	else
	snc = -1.0d0	
	endif
	
	
	ELSE IF (iflag == 2) THEN

  	e1x =  (xs - xm)/sqrt((xs - xm)**2 + (ys - ym)**2)
	e1y =  (ys - ym)/sqrt((xs - xm)**2 + (ys - ym)**2)
	e2x = -(ys - ym)/sqrt((xs - xm)**2 + (ys - ym)**2)
	e2y =  (xs - xm)/sqrt((xs - xm)**2 + (ys - ym)**2)

	cx1 = cos(p(1)-thetam)*e1x+sin(p(1)-thetam)*e2x
	cy1 = cos(p(1)-thetam)*e1y+sin(p(1)-thetam)*e2y	

	cx2 = -cos(p(1)-thetas)*e1x-sin(p(1)-thetas)*e2x
	cy2 = -cos(p(1)-thetas)*e1y-sin(p(1)-thetas)*e2y	

	
	dcx1= -sin(p(1)-thetam)*e1x+cos(p(1)-thetam)*e2x
	dcy1= -sin(p(1)-thetam)*e1y+cos(p(1)-thetam)*e2y

	dcx2=  sin(p(1)-thetam)*e1x-cos(p(1)-thetam)*e2x
	dcy2=  sin(p(1)-thetam)*e1y-cos(p(1)-thetam)*e2y
	

	phim = atan2(rm2*cy1,rm1*cx1) 
	phis = atan2(rs2*cy2,rs1*cx2) 

	xcm = (rm1*cos(phim)*cos(thetam) - rm2*sin(phim)*sin(thetam)) + xm
	ycm = (rm1*cos(phim)*sin(thetam) + rm2*sin(phim)*cos(thetam)) + ym

	xcs = (rs1*cos(phis)*cos(thetas) - rs2*sin(phis)*sin(thetas)) + xs
	ycs = (rs1*cos(phis)*sin(thetas) + rs2*sin(phis)*cos(thetas)) + ys	
	
	argm=rm2*cy1/(rm1*cx1)	
	args=rs2*cy2/(rs1*cx2)

	
	dargm=rm2/rm1*(dcy1*cx1-cy1*dcx1)/cx1**2	
	dargs=rs2/rs1*(dcy2*cx2-cx2*dcx2)/cx2**2
	
	dphim=1.0d0/(1.0d0+argm**2)*dargm
	dphis=1.0d0/(1.0d0+args**2)*dargs
	

	dxcm = (-rm1*sin(phim)*dphim*cos(thetam) 
     &  - rm2*cos(phim)*dphim*sin(thetam))
	dycm = (-rm1*sin(phim)*dphim*sin(thetam)
     &  + rm2*cos(phim)*dphim*cos(thetam)) 

	dxcs = (-rs1*sin(phis)*dphis*cos(thetas)
     &  - rs2*cos(phis)*dphis*sin(thetas))
	dycs = (-rs1*sin(phis)*dphis*sin(thetas)
     &  + rs2*cos(phis)*dphis*cos(thetas)) 	


	dmod = 2.0d0*(xcs-xcm)*(dxcs-dxcm)+2.0d0*(ycs-ycm)*(dycs-dycm)
	
	ddp = 1.0d0/(2.0d0*sqrt((xcs-xcm)**2+(ycs-ycm)**2))*dmod
	
	
	
	fjac(1:m,1) = ddp 
	
	pen(3,1)=phim
	pen(4,1)=phis
 

	END IF
	

	RETURN
	END SUBROUTINE fcn

















