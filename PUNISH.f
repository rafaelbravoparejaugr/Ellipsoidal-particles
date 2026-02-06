c     rem ******************** 18000
c     rem add and subtract punishment matrices
c     rem ******************** called by main

      subroutine punish (idp,damp,n2,n2p,n2b,g0,k00,dd)

c     inputs -- n0,n2,m0(,1-2),g0,dm(),k1(),k4(,1-2),a(),f(),k(),
c               k00,o(,1-2)
c     create -- c0(),q(2,4),s1,s2,s(),a(),f(),c0()
c     output -- a(),f(),c0()

      implicit      none

      include      'param.h'
      include      'za.h'
      include      'zc0.h'
      include      'zdm.h'
      include      'zf.h'
      include      'zh.h'
      include      'zk.h'
      include      'zk1.h'
      include      'zk4.h'
      include      'zm0.h'
      include      'zn.h'
      include      'znp0.h'
      include      'zo.h'
      include      'zq.h'
      include      'zs.h'

      integer       idp,n2,n2p,n2b,k00
      real*8        damp,g0,dd

      logical       lflag, gflag, iflag
      integer       i,i1,i2,i3, j,j1,j2,j3, l
      real*8        temp1,temp2, vnr,vsh, s1,s2,s4,s5,s6

      integer       indf
      real*8        sgn

c     md: maximum number of disks
c     mi: maximum number of invasions
c     ma: maximum number of elements in a[]
c     mk: maximum number of elements in k[]

c     n0        : number of disks
c     n2        : number of d-d contacts at start of current step
c     n2b       : n2 + number of d-b contacts at start of current step
c     g0        : penalty (contact stiffness) (e10 - e12)
c     k00       : 0 for m9=0, then m0(,1)=0; 1 for m9>1, then m(,1)=m(,2)
c     dd        : conversion factor between radian and degree
c     c0()      : friction force terms of each iteration
c     m0(,1-2)  : contact indices of at end of previous and at start of current
c                 iterations
c     q()       : temporary use for saving punish index of each invasion
c     dm(,1-6)  : ball #,boundary #,dips,phi,coh of a contact
c                 at end and start of previous step
c     k1()      : permutation matrix -  k1(old dk no)=new dk no
c     s(30)     : punish geometry terms
c     k4(i,1-2) : starting number of k(),number of k() for i-th disk
c     a()       : stiffness matrix at current iteration
c     f()       : free terms at current iteration
c     k()       : disk contact matrix
c     o(,1-2)   : distance of penetration, sliding distance w.r.t reference line
c              of each invasion at start of current iteration

c     Set zero initial values of friction force term at start of each iteration

      call timer('punish',23,1)
      do j = 1, 3
        do i = 1, n0
          c0(j,i) = 0.d0
        end do ! i
      end do ! j

c     Add/substract punishment matrix according to contact conditions

      do 2000 i = 1, n2b

c       Set punish index according to change phase between successive it.s

        do j = 1, 2
          if (m0(j,i).eq.0) then
            q(j,1) = 0.d0
            q(j,2) = 0.d0
          elseif (m0(j,i).eq.1) then
            q(j,1) = 1.d0
            q(j,2) = 0.d0
          else
            q(j,1) = 1.d0
            q(j,2) = 1.d0
          endif
        end do ! j

        do j = 1, 2
          q(1,j) = q(2,j) - q(1,j)
        end do ! j

c       Set punish stiffness

        q(1,1) = q(1,1)*g0

c       Change shear stiffness to equal normal stiffness

c       q(1,2) = 0.1d0*q(1,2)*g0
        q(1,2) = q(1,2)*g0

        if (q(1,1).ne.0.d0 .or. q(1,2).ne.0.d0) goto 2304
        if (m0(2,i).eq.1.and.
     &       (m0(1,i).ne.0.or.k00.eq.0)) goto 2304

        goto 2000

2304    continue

        i1 = dm(1,i)
        i2 = dm(2,i)

c       Call subroutine to find s1,s2,s5,s6,s()

        call pun1 (i,n2,n2p,i1,i2,s1,s2,s5,s6)

        lflag = .not.(i.gt.n2p .and. i2.lt.0)
        gflag = .not.(i.gt.n2p .and. i2.gt.0)
        iflag = cex.gt.0.0d0

        i1    = k1(1,i1)
        j1    = i1

        if (i.le.n2) then
          i2 = k1(1,i2)
          j2 = i2
        elseif (i.le.n2p) then
          if (i2.gt.0) then
            j2 = np0( i2) - nd + nd0
          else
            j1 = np0(-i2) - nd + nd0
            j2 = i1
          endif
        else
          if (i2.lt.0) j2 = i1
        endif

        if (q(1,1).eq.0.d0 .and. q(1,2).eq.0.d0) goto 2389

c       Add punishment to term(j1,j1)

        if (lflag .and. iflag) then
          j3 = k4(1,j1)
          do j = 1, 3
            temp1 = q(1,1)*s(j)
            temp2 = q(1,2)*s(j+6)
            do l = 1, 3
              a(l,j,j3) = a(l,j,j3) + temp1*s(l) + temp2*s(l+6)
            end do
          end do ! j
        endif

c       Add punishment to term(j2,j2)

        if (gflag .and. iflag) then
          j3 = k4(1,j2)
          do j = 1, 3
            temp1 = q(1,1)*s(j+3)
            temp2 = q(1,2)*s(j+9)
            do l = 1, 3
              a(l,j,j3) = a(l,j,j3) + temp1*s(l+3) + temp2*s(l+9)
            end do
          end do ! j
        endif

        if (i.gt.n2p) goto 2376

c         locate punishment only to upper triangle
c         add punishment to term (j1,j2) if j2>=j1

        if (j1.lt.j2 .and. iflag) then

          i3 = indf(j1,j2)
          do j = 1, 3
            temp1 = q(1,1)*s(j)
            temp2 = q(1,2)*s(j+6)
            do l = 1, 3
              a(l,j,i3) = a(l,j,i3) + temp1*s(l+3) + temp2*s(l+9)
            end do
          end do ! j
          goto 2376

c         Add punishment to term (j2,j1) if j1>j2

        endif

        if(iflag) then
          i3 = indf(j2,j1)
          do j = 1, 3
            temp1 = q(1,1)*s(j+3)
            temp2 = q(1,2)*s(j+9)
            do l = 1, 3
              a(l,j,i3) = a(l,j,i3) + temp1*s(l) + temp2*s(l+6)
            end do
          end do ! j
        endif ! iflag

c       Add punishment to free term j1 and j2

2376    continue

        vnr = 0.d0
        vsh = 0.d0
        if (idp.eq.1) then
          if (lflag) then
            do j = 1,3
              vnr = vnr + h(j,j1)*s(j)
              vsh = vsh + h(j,j1)*s(j+6)
            end do
          end if
          if (gflag) then
            do j = 1,3
              vnr = vnr + h(j,j2)*s(j+3)
              vsh = vsh + h(j,j2)*s(j+9)
            end do
          end if

        endif

        if(idp.eq.1) then
          temp1 = q(1,1)*(s1 + s5 + damp*vnr)
          temp2 = q(1,2)*(s2 + s6 + damp*vsh)
        else
          temp1 = q(1,1)*(s1 + s5)
          temp2 = q(1,2)*(s2 + s6)
        endif
        if (lflag) then
          do j = 1, 3
            f(j,j1) = f(j,j1) - temp1*s(j) - temp2*s(j+6)
          end do
        endif
        if (gflag) then
          do j = 1, 3
            f(j,j2) = f(j,j2) - temp1*s(j+3) - temp2*s(j+9)
          end do
        endif

c       Compute sliding friction force

2389    continue

        if (m0(2,i).ne.1.or.
     &     (m0(1,i).eq.0 .and. k00.eq.1)) goto 2000

        s4 = -o(1,i)
        if (s4.ge.0.d0) then
          s4 = (s4*g0*tan(dd*dm(5,i)) + dm(6,i))*sgn(o(2,i))
        else
          goto 2000
        endif

        if (i.le.n2p) then
          do j = 1, 3
            c0(j,j1) = c0(j,j1) - s4*s(j+6)
            c0(j,j2) = c0(j,j2) + s4*s(j+12)
          end do ! j
        else
          if (i2.gt.0) then
            do j = 1, 3
              c0(j,j1) = c0(j,j1) - s4*s(j+6)
            end do ! j
          else
            do j = 1, 3
              c0(j,j1) = c0(j,j1) + s4*s(j+12)
            end do ! j
          endif
        endif

2000  continue

      call timer('punish',23,2)
      end
