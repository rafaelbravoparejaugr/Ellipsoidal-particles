

	subroutine SEARCH()
	


	USE Levenberg_Marquardt
	USE data_storage
	IMPLICIT NONE
			
	INTERFACE
	  SUBROUTINE fcn(m, n, x, fvec, fjac, iflag)
		IMPLICIT NONE
		INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
		INTEGER, INTENT(IN)        :: m, n
		REAL (dp), INTENT(IN)      :: x(:)
		REAL (dp), INTENT(IN OUT)  :: fvec(:)
		REAL (dp), INTENT(OUT)     :: fjac(:,:)
		INTEGER, INTENT(IN OUT)    :: iflag
	  END SUBROUTINE fcn
	END INTERFACE


	REAL (8), ALLOCATABLE  :: fvec(:), fjac(:,:)
	INTEGER, PARAMETER      :: n = 1                 ! The number of parameters
	INTEGER                 :: info, iostatus, ipvt(n), m
	REAL (8)               :: p(n), tol = 0.0000000000000000000001d0

	integer i

	include     'param.h'	  
	include     'zrd.h'	

	m = 1

	ALLOCATE( fvec(m), fjac(m,n) )

	! Set starting values for parameters.

	p(1) = 0.00000000000001d0

	CALL lmder1(fcn, m, n, p, fvec, fjac, tol, info, ipvt)
		 
		 pen(1,1)=fvec(1)		 		 

	END SUBROUTINE 

















