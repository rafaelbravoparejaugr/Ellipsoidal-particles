
c     rem ******************** 19000 ***
c     rem pcg equation solver
c     rem ******************** called by main

      subroutine tridc(i00)

c     input  -- n3,n0,a(),f(),c0(),k4(,1-2),k()
c     create -- b(),r(),f(),q(),e(),a(),q(),e(),v()
c     output -- b(),r(),a(),v()

      implicit     none

      include     'param.h'
      include     'za.h'
      include     'zb.h'
      include     'zc0.h'
      include     'ze.h'
      include     'zf.h'
      include     'zk.h'
      include     'zk4.h'
      include     'zn.h'
      include     'zq.h'
      include     'zr.h'
      include     'zt.h'
      include     'zv.h'

      logical      noconv, notdiag
      integer      i,i1,itr, j, l,l1, i00
      real*8       rr0,rri, alp,bet,gam, temp1
      real*8       p(3,md), z(3,md)

      real*8       dot

c     md       : maximum number of disks
c     mi       : maximum number of invasions
c     ma       : maximum number of elements in a[]
c     mk       : maximum number of elements in k[]

c     n3       : number of nonzero elements in a()
c     n0       : number of disks
c     a()      : (1) stiffness matrix at current iteration; (2) resultant 
c                lower triangle matrix in a=ldl
c     b()      : stiffness matrix at current iteration saved for next
c                iteration
c     f()      : free terms at current iteration
c     r()      : free terms at current iteration, excluding c0()
c     c0()     : friction force at current iteration
c     k4(i,1-2): starting number of k(),number of k() for i-th disk
c     k()      : disk contact matrix
c     q(),e()  : 6*6 matrix for temporary use
c     v()      : diagonal matrix after triangular decomposition: a=ldl

c     itr      : iterations of solution;
c     e1       : solution error;
c     e2       : sum of solutions

c     Convergence tolerance

      real*8     tol
      data       tol /1.d-6/

c     Inverse 3*3 diagonal submatrices

      call timer('tridc',31,1)
      notdiag = i00.ne.2

c      write(2,3000) (i,(f(j,i),j=1,3),i=1,n0)
c3000  format(' Forces '/(i6,1p,3e12.4))
      do i = 1, n0
        j  = k4(1,i)
        do l = 1, 3
          do l1 = 1, 3
            q(l1,l) = a(l1,l,j)
          end do ! l
        end do ! l

c       Invert q-array: Store in e-array

        e(1,1) =  q(2,2)*q(3,3) - q(2,3)*q(3,2)
        e(1,2) =  q(3,2)*q(1,3) - q(3,3)*q(1,2)
        e(1,3) =  q(1,2)*q(2,3) - q(1,3)*q(2,2)

        temp1  =  (q(1,1)*e(1,1) + q(2,1)*e(1,2) + q(3,1)*e(1,3))

        if(temp1.eq.0.0d0) then
          i1 = k4(1,n0) + k4(2,n0) - 1
          write(*,2000) i,temp1,n0,j,i1,((q(l,l1),l1=1,3),l=1,3)
          write(2,2000) i,temp1,n0,j,i1,((q(l,l1),l1=1,3),l=1,3)
          stop
        endif

        temp1  =  1.d0/temp1
        e(1,1) =  e(1,1)*temp1
        e(1,2) =  e(1,2)*temp1
        e(1,3) =  e(1,3)*temp1
        e(2,1) = (q(2,3)*q(3,1) - q(2,1)*q(3,3))*temp1
        e(2,2) = (q(3,3)*q(1,1) - q(3,1)*q(1,3))*temp1
        e(2,3) = (q(1,3)*q(2,1) - q(1,1)*q(2,3))*temp1
        e(3,1) = (q(2,1)*q(3,2) - q(2,2)*q(3,1))*temp1
        e(3,2) = (q(3,1)*q(1,2) - q(3,2)*q(1,1))*temp1
        e(3,3) = (q(1,1)*q(2,2) - q(1,2)*q(2,1))*temp1

        do l = 1, 3
          do l1 = 1, 3
            v(l1,l,i) = e(l1,l)
          end do ! l1
        end do ! l
      end do ! i

c     Set initial values of free terms and diagonal term

      do i = 1, n0
        j = k4(1,i)
        i1 = k(j)
        do l = 1, 3
          r(l,i) = f(l,i) + c0(l,i)
          do l1 = 1, 3
            r(l,i) = r(l,i) - a(l1,l,j)*b(l1,i1)
          end do ! l1
        end do ! l
      end do ! i

c     Compute initial residual

      if(notdiag) then
      do i = 1, n0
c       Operation on upper and lower triangle

        do j = k4(1,i) + 1, k4(1,i) + k4(2,i) - 1
          i1 = k(j)
          do l1 = 1, 3
            do l = 1, 3
              r(l ,i ) = r(l ,i ) - a(l1,l,j)*b(l1,i1)
              r(l1,i1) = r(l1,i1) - a(l1,l,j)*b(l ,i )
            end do ! l1
          end do ! l
        end do ! j
      end do ! i
      endif ! notdiag

      do i = 1,n0
        do l = 1, 3
          z(l,i) = 0.d0
          do l1 = 1, 3
            z(l,i) = z(l,i) + v(l1,l,i)*r(l1,i)
          end do ! l1
          p(l,i) = z(l,i)
        end do ! l
      end do ! i

      gam = dot(r,z,3*n0)
      rr0 = sqrt(dot(r,r,3*n0))*tol
      if(rr0.le.0.0d0) then
        do i = 1,n0
          do l = 1,3
            b(l,i) = p(l,i)
          end do
        end do
        noconv = .false.
      else
        noconv = .true.
      endif

c     Iteration starts

      itr     =  0
      do while( noconv )

c       Compute: z = A*p

        do i = 1, n0
          do l = 1, 3
            z(l,i) = 0.0d0
          end do ! l
          j  = k4(1,i)
          i1 = k(j)
          do l = 1, 3
            do l1 = 1, 3
              z(l,i) = z(l,i) + a(l1,l,j)*p(l1,i1)
            end do ! l1
          end do ! l
        end do ! i

        if(notdiag) then
        do i = 1, n0

c         Operation on upper and lower triangle

          do j = k4(1,i) + 1, k4(1,i) + k4(2,i) - 1
            i1 = k(j)
            do l1 = 1, 3
              do l = 1, 3
                z(l ,i ) = z(l ,i ) + a(l1,l,j)*p(l1,i1)
                z(l1,i1) = z(l1,i1) + a(l1,l,j)*p(l ,i )
              end do ! l1
            end do ! l
          end do ! j
        end do ! i
        endif ! notdiag

        alp = gam/dot(p,z,3*n0)

        do i = 1, n0
          do l = 1, 3
            b(l,i) = b(l,i) + alp*p(l,i)
            r(l,i) = r(l,i) - alp*z(l,i)
          end do ! l
        end do ! i

        rri = sqrt(dot(r,r,3*n0))

        do i = 1,n0
          do l = 1, 3
            z(l,i) = 0.d0
            do l1 = 1, 3
              z(l,i) = z(l,i) + v(l1,l,i)*r(l1,i)
            end do ! l1
          end do ! l
        end do ! i

        bet = gam
        gam = dot(r,z,3*n0)
        bet = gam/bet

        do i = 1, n0
          do l = 1, 3
            p(l,i) = z(l,i) + bet*p(l,i)
          end do ! l
        end do ! i

c       Error and branch of iteration

        itr = itr + 1
        if (rri.lt.rr0 .or. itr.gt.500) noconv = .false.

      end do ! while noconv

	write(*,*) 'b',b(1,1),b(2,1),b(3,1)   

c      write(2,3001) itr,(i,(b(j,i),j=1,3),i=1,n0)
c3001  format(' solution, iterations = ',i8/(i6,1p,3e12.4))

c     write(99,*) 'Iters=',itr,' rr0,rri=',rr0,rri
      if (itr.gt.500) then
        write(*,*) 'No convergence in TRIDC:'
        stop
      endif
      call timer('tridc',31,2)

2000  format(' ERROR: Tridc: determinant(',i6,') =',1p,e13.5/
     &       'Equation Blocks =',i6,' Block =',i6,' Total Blocks =',i6//
     &      (5x,1p,3e13.5/))

      end

      function dot(a,b,nn)

      implicit none

      integer  n,nn
      real*8   dot, a(*),b(*)

      dot = 0.0d0
      do n = 1,nn
        dot = dot + a(n)*b(n)
      end do ! n

      end
