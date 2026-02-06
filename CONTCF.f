
      subroutine contcf (n2,n2b,g0,dd,mc,ncf,q02,i9,n2p)

      implicit         none

      integer          n2,n2b,mc,ncf,i9,n2p
      real*8           g0,dd,q02

      real*8           sgn

      integer          i,i1,i2,i1c, j,j1,j1c, lmc,nmc,nc, nn
      real*8           e1,e2, psi, bdr, x2,x3, y2,y3, dx,dy, ss,sn

      include         'param.h'
      include         'zb.h'
      include         'zbd.h'
      include         'zbd0.h'
      include         'zcf.h'
      include         'zd.h'
      include         'zdm.h'
      include         'zdz.h'
      include         'zenergy.h'
      include         'zh0.h'
      include         'zib.h'
      include         'zicltyp.h'
      include         'zk0b.h'
      include         'zk1.h'
      include         'zm0.h'
      include         'zmem3.h'
      include         'zn.h'
      include         'zo.h'
      include         'zrd.h'
      include         'zstr.h'

c     mp: maximum number of input points/boundary vertices
c     mi: maximum number of invasions

c     n0: number of disks
c     n2b: n2 + number of d-b contacts at start of current step
c     mc: number of contact forces at end of current step
c     m0(,2): contact indices at end of current step
c     o(,2): sliding distance w.r.t reference line
c              of each invasion at start of current iteration
c     dm(,1-5): contact ball number/boundary no.,dips,contact friction angle
c               at end and start of previous step
c     mm(,1-4): mm(,1-3) are vertex numbers of real contact and mm(,4) contact
c               type at end of current step
c     cf(): contact forces - normal and shear forces
c     topnf: normal force on 2nd brdy requested (topcap)
c     topsf: shear force on 2nd brdy requested (topcap)
c     botnf: normal force on 1st brdy requested (base)
c     botsf: shear force on 1st brdy requested (base)
c     nmc: number of contacts in n2b that are not open
c    ***8/13/96 correction: add i1c and j1c for use in b()
c        b() expecting cluster number not disk number
c     i1c=cluster number for disk 1 (i1)
c     j1c=cluster number for disk2 (j1)

c     Compute contact forces--need n2b,m0(),o(),g0,dd,dm()

      call timer('contcf',4,1)
      mc  = 0
      psi = 0.0001
      bdr = 1.d0/cos(psi)

      do i = 1, n2b

        if (m0(2,i).ne.0) then
          i1 = dm(1,i)
          j1 = dm(2,i)
          i1c = k1(1,i1)
          j1c = k1(1,j1)
          mc = mc + 1
          if (i.le.n2) then
            lmc = mc
            rd(1,1,mc) = d(1,i1) + d(3,i1)*cos(dm(3,i) + b(3,i1c))
            rd(1,2,mc) = d(2,i1) + d(3,i1)*sin(dm(3,i) + b(3,i1c))
            rd(2,1,mc) = d(1,j1) + d(3,j1)*bdr*cos(dm(4,i) + b(3,j1c)
     &                           + psi)
            rd(2,2,mc) = d(2,j1) + d(3,j1)*bdr*sin(dm(4,i) + b(3,j1c)
     &                           + psi)
            rd(3,1,mc) = d(1,j1) + d(3,j1)*bdr*cos(dm(4,i) + b(3,j1c)
     &                           - psi)
            rd(3,2,mc) = d(2,j1) + d(3,j1)*bdr*sin(dm(4,i) + b(3,j1c)
     &                           - psi)
          else
            if (j1.gt.0) then
              rd(1,1,mc) = d(1,i1) + d(3,i1)*cos(dm(3,i) + b(3,i1c))
              rd(1,2,mc) = d(2,i1) + d(3,i1)*sin(dm(3,i) + b(3,i1c))
              rd(2,1,mc) = bd(1,j1)
              rd(2,2,mc) = bd(2,j1)
              rd(3,1,mc) = bd(1,j1+1)
              rd(3,2,mc) = bd(2,j1+1)
            else
              rd(1,1,mc) = bd(1,-j1)
              rd(1,2,mc) = bd(2,-j1)
              rd(2,1,mc) = d(1,i1) + d(3,i1)*bdr*cos(dm(3,i) + b(3,j1c)
     &                             + psi)
              rd(2,2,mc) = d(2,i1) + d(3,i1)*bdr*sin(dm(3,i) + b(3,j1c)
     &                             + psi)
              rd(3,1,mc) = d(1,i1) + d(3,i1)*bdr*cos(dm(3,i) + b(3,j1c)
     &                             - psi)
              rd(3,2,mc) = d(2,i1) + d(3,i1)*bdr*sin(dm(3,i) + b(3,j1c)
     &                             - psi)
            endif
          endif

          if (m0(2,i).eq.1) then

            cf(2,mc) = -o(1,i)*g0

            if (o(1,i).gt.0.d0) cf(2,mc) = 0.d0
            cf(3,mc) = (cf(2,mc)*tan(dd*dm(5,i)) + dm(6,i))*sgn(o(2,i))
            if (cf(3,mc).ne.0.d0) cf(1,mc) = 1.d0
            if (cf(3,mc).eq.0.d0) cf(1,mc) = sgn(o(2,i))

c         e-a locked case: m0(,2) = 2

          else
            cf(1,mc) = 2.d0
            cf(2,mc) =-o(1,i)*g0
            cf(3,mc) = o(2,i)*g0
c           cf(3,mc) = o(2,i)*g0*1.d-1
          endif
        endif
      end do ! i

      nmc = mc

c     Compute total contact forces along selected boundary segments

      if (ncf.eq.0) goto 999

      do i = 1, ncf
        i1 = k0b(1,ib(i,1))
        i2 = i1 + ib(i,2)
        x2 = bd(1,i2)
        y2 = bd(2,i2)
        x3 = bd(1,i2+1)
        y3 = bd(2,i2+1)
        nn = 0
        sn = 0.d0
        ss = 0.d0
        do j = lmc+1, mc
          if (rd(2,1,j).eq.x2.and.rd(2,2,j).eq.y2  .and.
     &        rd(3,1,j).eq.x3.and.rd(3,2,j).eq.y3) then
            nn = nn + 1
            sn = sn + cf(2,j)
            ss = ss + cf(3,j)
          endif
        end do ! j
        dx = bd(1,i2) - bd0(1,i2)
        dy = bd(2,i2) - bd0(2,i2)

c       Save normal force on "endcaps" for output in membrane
c       this assumes first bdry requested is bottom cap and
c       second bdry is topcap (poly or bdry)

        if (i.eq.1) then
          botnf = sn
          botsf = ss
        else
          topnf = sn
          topsf = ss
        endif

      end do ! i

c     Compute energy stored in spring strain

c     energy = 1/2*k*x^2 = 1/2*f^2/k

      if (ieng0.eq.1) then
        if (i9.eq.1.or.mod(i9,ieng).eq.0) then
          e1    = 0.d0
          e2    = 0.d0
          esprn = 0.d0
          esprs = 0.d0
          mc    = 0.d0
          do i = 1,n2b
            if (m0(2,i).ne.0) then
              mc    = mc + 1
              e1    = (cf(2,mc)*cf(2,mc))/(2.d0*g0)
              esprn = esprn + e1
              e2    = (cf(3,mc)*cf(3,mc))/(g0*0.2d0)
              esprs = esprs + e2
            endif
          end do ! i
        endif
      endif

c     Output data for stress calculations (into file 18 'strdat')

      if (istr0.eq.1) then
        if (i9.eq.1.or.mod(i9,istr).eq.0) then
          write (18,*) nd,np,nd0
          nc = 0

          do i = 1,nd
            write (18,*) (d(j,i),j=1,3),icltyp(i),(k1(j,i),j=1,2)
          end do ! i

          do i = 1,nd0
            write (18,402) dz(1,i),dz(2,i),h0(1,i)
          end do ! i

          write (18,*) q02
          write (18,*) n2,n2p,n2b,nmc,lmc
          do i = 1,n2b
            if (m0(2,i).ne.0) then
              nc = nc + 1
              write (18,401) (cf(j,nc),j=1,3),dm(1,i),dm(2,i),
     &                   dm(4,i),rd(1,1,nc),rd(1,2,nc)
            endif
          end do ! i
        endif
      endif

999   return
      call timer('contcf',4,2)

c     Formats

401   format (1x,f5.0,2e20.8,2x,2f7.0,3e20.8)
402   format (1x,3e20.8)

      end
