      program dda

      implicit         none

      include         'param.h'
      include         'files.h'
      include         'zb.h'
      include         'zc.h'
      include         'zd.h'
      include         'zdm.h'
      include         'ze.h'
      include         'zh.h'
      include         'zh0.h'
      include         'zib.h'
      include         'zk.h'
      include         'zm0.h'
      include         'zn.h'
      include         'znc.h'
      include         'znam.h'
      include         'znp0.h'
      include         'zo.h'
      include         'zp.h'
      include         'zq.h'
      include         'zscl.h'
      include         'zt.h'
      include         'zw.h'
      include         'zenergy.h'
      include         'zicltyp.h'
      include         'zmemb.h'
      include         'zmem2.h'
      include         'zmem4.h'
      include         'zmem6.h'
      include         'zstr.h'
      include         'zstatic.h'

      real*4           etime,tt(2),tim0
      character        title*11, date*24

      logical          exsti,exstt,exstp
      integer          i,i0d,i00,i0c,i9,j,k00
      integer          istep,icstep,irstep,idp,ik,innb,innbd0
      integer          iro,irs,ist
c     integer          j1,j2
      integer          m00,m3,m9,mc,mcode,m3b, nreset,i9last
      integer          ng,nst,ntd,ncf,n2,n2b,n2p,n3,n5,n9
      real*8           d0,damp,dd,ddt
      real*8           e1,e2,e3,ekint,ekinr,epot,est,ev2
      real*8           g0,g1,g2,g3, pi, g1max
      real*8           q01,q02,q03
      real*8           rmax, s0, NBeta, NGamma
      real*8           ji, rhof, rug, nu, us

      real*8           tm(2000)

      real*8           sgn

c     mp : maximum number of input points/boundary vertices
c     md : maximum number of disks+polygons
c     mi : maximum number of invasions
c     ma : maximum number of elements in a[]
c     mk : maximum number of elements in k[]
c     mbb: maximum number of polygons/boundaries

c     n0     : nd0 + np
c     nd     : total number of disks
c     nd0    : number of disk groups
c     np     : number of rigid polygons
c     nb     : number of rigid boundaries
c     ip(1,1): n1, number of fixed points
c     ip(2,1): n1a, number of 1-d constraint points
c     ip(3,1): n1b, number of spring points
cxxx  n1c    : number of one-direction moving constraint points xxxx
c     ip(4,1): n1c, number of pin constraint points
c     ip(5,1): nbt, number of bolting points
c     ip(6,1): n1x, number of no rotation points
c     n2     : number of d-d contacts at start of current step
c     n2p    : n2 + number of d-p contacts at start of current step
c     n2b    : n2p + number of d-b contacts at start of current step
c     n3     : number of non-zero 6*6 matrices of a
c     n5     : maximum number of time steps allowed
c     ip(7,1): n6, number of load points
c     ip(8,1): n7, number of measured points
c     n9     : total number of iterations
c     m3     : number of d-d contacts at end of previous step
c     m3b    : m3 + number of d-p + d-b contacts at end of previous step
c     m9     : number of step iterations
c     i00    : 0 for static; 1 for implicit dynamic; 2 for explicit dynamic
c     i9     : current time step number
c     k00    : 0 for before step iterations, with m0(,1)=0; 1 for during step
c              iterations, with m(,1)=m(,2)
c     g0     : penalty (contact stiffness) (e10 - e12)
c     g1     : time step (in sec)
c     g2     : ratio of maximum displacement per step (.0001 - .01)
c     g3     : current time step used
c     dd     : conversion factor between radian and degree
c     d0     : criterion of possible contact distance of every two disks
c     s0     = .0001 as criteria of open-close in sub judgea
c     innb   :  interval of # timesteps when nearest neighbor search is done
c     innbdo : multiple of search dist. for search for neighbors
c     icstep :  interval of timesteps to output bdry contact forces

c     Open files for outputs


      open ( 2,file='result')
      open ( 3,file='grf1')
      open ( 4,file='out1')
      open ( 8,file='kinetic')
      open (11,file='out')
      open (15,file='strsstrn')
      open (16,file='force')
      open (17,file='rotate')
      open (18,file='strdat')
      open (19,file='xdebug')

c     Start timing

      call timer('dda',0,0)

      write(2,2004) date

c     Input data file name

      write (*,3000)
c     read  (*,'(a)') fname  
      fname='cont'

      finte = fname
      fext  = 'int'
      call addext(finte,fext,20,4)

      ftemp = fname
      fext  = 'tem'
      call addext(ftemp,fext,20,4)

      fpara = fname
      fext  = 'par'
      call addext(fpara,fext,20,4)

      inquire(file = finte, exist = exsti )
      inquire(file = ftemp, exist = exstt )
      inquire(file = fpara, exist = exstp )

      if(.not.exsti .or. .not.exstt .or. .not.exstp ) then
        stop ' Data file does not exist'
      endif

      open (7,file = finte)

c     Enter number of output at specified time, ntd

      read (7,*) ntd

c     Enter tm(j),j=1,ntd

      read (7,*) (tm(j),j=1,ntd)

c     Statics control; ist = 1 for statics, 0 for normal

      read (7,*) ist

c     Statics criterion, est (e-8 - e-16)? '

      if (ist.eq.1) then
        read (7,*) est
      endif

c     Option to output contact forces along selected boundary segments
c     i0c = 1 for contact forces, 0 for none

      read (7,*) i0c

      if (i0c.eq.1) then
        open (12,file='contact')
        write (12,2000)

c       Interval of timesteps to ouput forces:   '

        read (7,*) icstep

c       Number of boundary forces to output, ncf:
c         ib(,1): boundary number
c         ib(,2): vertex number

        read (7,*) ncf
        do i = 1, ncf
          read (7,*) (ib(i,j),j = 1,2)
        end do ! i
      else
        ncf = 0
      endif

c     Damping control code, idp = 0:n;1-3:y

      read (7,*) idp
      if (idp.ne.0) then

c       Damping rati0/g0
        read (7,*) damp
      endif

c     Name for restart files and # char. (max=4)'

      read (7,*) rsname,nc

c     Interval for nearest neighbor search:

      read (7,*) innb

c     Multiple of search distance for nearest neighbor&search

      read (7,*) innbd0

c     Particle rotations output, iro = 0:no ;1 yes

      read (7,*) iro
      if (iro.eq.1) then
c        nterval of timestep to output rotation
        read (7,*) irstep
      else
        irstep = 0
      endif

c     Enter info for stress membrane boundaries,m00 =1,yes;=0,no

      read (7,*) m00
      if (m00.eq.1) then

c       Title of test(max.=11 char):  '

        read (7,*) title

c       Increment for output of strs-strn data:

        read (7,*) istep

c       N.B. Only have vertical boundaries so far

        read (7,*) m01
        if (m01.eq.1) then

c         boundary number of top cap:

          read (7,*) nbtopc

c         Vertex number for vertex on bottom of top cap:

          read (7,*) nbvtopc
        else

c         Control point # of constraint point on bottom of top cap:

          read (7,*) mcp

c         Final   sigma 1 (atm):
c         Initial sigma 1 (atm):

          read (7,*) sig1app
          read (7,*) sig1int
          if (sig1app.ne.sig1int) then

c           # steps to turn on sigma 1:

            read (7,*) itosig1
          else
            itosig1 = 1
          endif
        endif

c       Final   sigma 3 (atm):
c       Initial sigma 3 (atm):

        read (7,*) sig3app
        read (7,*) sig3int
        if (sig3app.ne.sig3int) then
c         # steps to turn on force:
          read (7,*) itosig3
        else
          itosig3 = 1
        endif

c       Convert atmospheres to pascals

        sig1app = sig1app*101330.0d0
        sig3app = sig3app*101330.0d0
        sig1int = sig1int*101330.0d0
        sig3int = sig3int*101330.0d0

      endif

c     Contact data for average stress data

      read (7,*) istr0
      if (istr0.eq.1) then
c       interval of timestep to output data:
        read (7,*) istr
      endif

c     Calculate energy of particles, ieng0 = 1

      read (7,*) ieng0
      if (ieng0.eq.1) then
c       Interval of timestep to output data:
        read (7,*) ieng
      endif

      i9 = 0
      n9 = 0
      pi = acos(-1.d0)
      dd = pi/180.d0
      s0 = 0.0001d0

      call datain (pi,q02,i00,i0d,n5,g0,g1,g2,ng,ntd,rmax,nst,NBeta,
     &NGamma,ji,rhof,rug,nu,us)

      d0 = 2.001d0/0.8d0*(w(4) - w(2))*g2

      write(19,*) 'd0 :',d0
      write(19,*) '1/2 window height :',(w(4) - w(2))

      ik  = 0
      irs = 0
      q01 = 0.d0
      q02 = 0.d0

      call grap  (ik,i9,q01,q02,q03,n9,ntd,ng,mc)

      call integ (i9,pi)

c     Running time

      tim0 = etime(tt)

c     Output for verification purposes

      write(19,*) 'nd, nd0, np, nb, n0,nbv',nd,nd0,np,nb,n0,nbv

      nreset = 0
      i9last = 0
      g1max  = g1
      g3old  = 0.0d0
      do ik = 1, ntd

c       Update state for explicit calculations

555     if(i00.eq.2) then
          call result (g2,g3,i9,q01,q02,q03,irstep,i00,NBeta,NGamma)
c         do i = 1,n0
c           do j = 1,3
c             b(j,i) = 0.0d0
c           end do ! j
c         end do ! i
        endif

        if (mod(i9,100).eq.0.0) then
          write(*,*) 'Time step = ',i9
          write(*,*) 'Real time = ',q02,' CPU time =',etime(tt) - tim0
        endif
        i9 = i9 + 1
        if (i9.gt.n5) then
          write(2,*) 'i9 = ',i9,' > n5 = ',n5
          stop
        endif
        g3  = g1
        ddt = tm(ik) - q02
        if (ddt.le.1.d-8) then
          q02 = tm(ik)
          i9  = i9 - 1

cTAYLOR          call drawpt (n2b,i9,0,0,.false.)
          write(*,*) 'Time step = ',i9
          write(*,*) 'Real time = ',q02,' CPU time =',etime(tt) - tim0

          call grap   (ik,i9,q01,q02,q03,n9,ntd,ng,mc)

          call restr  (q02,i00,i0d,n5,g0,g1,g2,nst,irs)

          call fabric (1,i9,q02,n2b,dd)

c         Output displacements of measured points at specified time

          write (11,2001) i9,q02,etime(tt) - tim0
          do i = ip(7,2) + 1,ip(8,2)
            write (11,2002) i-ip(7,2),-c(1,i),-c(2,i),
     $                      sqrt(c(1,i)**2 + c(2,i)**2)
          end do ! i

          goto 1000
        endif

        if (ddt.gt.0.and.ddt.le.g3) g3 = ddt

c       Check over penetration of particles - changed by RLT 10/97

        do i = 1, n2b
          if (o(1,i).lt.-0.10d0*d0 .and. i00.le.1) then
            g0 = g0*1.25d0
            write(2,*) 'Penalty increased to:',g0
            go to 1111
c         if (o(1,i).lt.-0.25d0*d0) then
c           j1 = dm(1,i)
c           j2 = dm(2,i)
c           call drawpt (n2b,i9,1,i,.true.)
c           write(*,*) 'Time step = ',i9
c           write(*,*) 'Real time = ',q02,' CPU time =',etime(tt)-tim0
c           call grap (ntd,i9,q01,q02,q03,n9,ntd,ng,mc)
c           write(*,*) ' Over penetration of particle',i
c           write(2,*) ' Over penetration of particle',i
c           goto 1009
          endif
        end do ! i
1111    continue
        if (i9.eq.1.or.mod(i9,10).eq.0) then
           mcode = 1
        else
           mcode = 0
        endif
        if (m00.eq.1) then
          call membrane (pi,mcode,i9,title,istep,q02)
        endif

        call contact (pi,d0,i9,n2,n2p,n2b,m3,m3b,innb,innbd0)

        call fabric (0,i9,q02,n2b,dd)

        call transf (n2,n2b,i9,m3,m3b)

        call sparse (n2,n2p)

        call point (q02)

c       Reset stiffness/free terms before new iteration upon changing g3

566     call judgeb (i9,n2,n2b,d0)

        call moving (q02,d0,g3)

        call inertia (n3,n2b,i00,g3,NBeta,NGamma)

        call drags (ji, rhof, rug, nu, us)

        if(ip(1,2).gt.0) then
          call fixed (g0)
        endif

        if(ip(3,1).ne.0) then
          call spring (g0)
        endif

        if(ip(4,1).ne.0) then
          call pin (g0)
        endif

        if(ip(5,1).ne.0) then
          call bolting (g0)
        endif

        if(ip(6,1).ne.0) then
          call xrotation (g0)
        endif

        if(ip(7,1).ne.0) then
          call ptload
        endif

        if(o1.ne.0.0d0 .or. o2.ne.0.0d0) then
          call volumn
        endif

        if(o3.ne.0.0d0 .or. o4.ne.0.0d0) then
          call mass
        endif

        k00 = 0
        m9  = 0

c       Iteration starts

575     m9 = m9 + 1
        n9 = n9 + 1

        call punish (idp,damp,n2,n2p,n2b,g0,k00,dd)

        call tridc(i00)

        call judgea (n2,n2p,n2b,g0,s0,dd,m9)

        if(i00.le.1) then
          call ratio (g2,q01,q03)
        else
c         write(2,*) 'B +',(( b(j,i),j=1,3),i=1,n0)
          do i = 1,n0
            do j = 1,3
              ha(j,i) = ha(j,i) + b(j,i)
c             ha(j,i) = b(j,i)
            end do ! j
          end do ! i
c         write(2,*) 'HA+',((ha(j,i),j=1,3),i=1,n0)
          if(m9.ge.2) m9 = -1
        endif

c       Adjust time step to keep convergence in 6 or less iterations: RLT 10/97

        k00 = 1
        if (m9.ge.6) then
          g3 = 0.5d0*g3
          nreset = nreset + 1
          if(nreset.ge.10) then
            if(float(nreset).gt.0.5*float(i9-i9last)) then
              nreset = 0
              i9last = i9
              g1 = g1*0.5d0
              write(2,*) ' dt(g1) reset to ',g1
            else
              nreset = 0
              i9last = i9
            endif
          endif
          write (2,2005) i9,m9,g3
          m9 = 0
          do i = 1, n2b
            m0(2,i) =  m0(3,i)
            o(2,i)  =  sgn(o(3,i))
            o(1,i)  = -abs(o(3,i))
          end do ! i
          goto 566
        endif

        if (m9.gt.0) goto 575

c       Check if q01 exceeds given g2 limit too much?

        if (q01.ge.3.d0 .and. i00.le.1) then
          write (2,*) 'Moved too much at time step = ',i9,' dmax:',
     &                 q01*g2*(w(4)-w(2))
c         pause

          g3 = g3/sqrt(q01)
          do i = 1, n2b
            m0(2,i) =  m0(3,i)
            o(2,i)  =  sgn(o(3,i))
            o(1,i)  = -abs(o(3,i))
          end do ! i
          goto 566
        endif

c       Try increasing time increment - iterations 3 or less - RLT 10/97

        if(abs(m9).le.3 .and. i00.le.1) then
          g1 = min(g1max,g1*1.1d0)
        endif

c       Update state for implicit calculations

        if(i00.le.1) then     
          call result (g2,g3,i9,q01,q02,q03,irstep,i00,NBeta,NGamma)
        endif

        call integ (i9,pi)

c       Calculate kinetic energy and potential energy

        if (ieng0.eq.1) then
          if (i9.eq.1.or.mod(i9,ieng).eq.0) then
            e1    = 0.d0
            ekint = 0.d0
            e2    = 0.d0
            ekinr = 0.d0
            e3    = 0.d0
            epot  = 0.d0
            do i = 1,n0
              e1    = 0.5d0*h0(2,i)*(h(1,i)*h(1,i) + h(2,i)*h(2,i))
              ekint = ekint + e1
              e2    = 0.5d0*h0(2,i)*h(3,i)*h(3,i)
              ekinr = ekinr + e2
              e3    = h0(2,i)*o4*d(2,i)
              epot  = epot + e3
            end do ! i
          endif
        endif

c       Save accumulated time into q02

        q02 = q02 + g3
        call static (ist,est,ev2,i9,rmax,q02)
        call contcf (n2,n2b,g0,dd,mc,ncf,q02,i9,n2p)

c       Output energies

        if (ieng0.eq.1) then
          if (i9.eq.1.or.mod(i9,ieng).eq.0.0) then
         if (i9.eq.1) then
          write(19,*) 'time   sprn    sprs    kint   kinr   pot'
         endif
        if(i9.eq.1.or.mod(i9,icstep).eq.0.0) then
          write(19,2003) q02,esprn,esprs,ekint,ekinr,epot
          endif
        endif
      endif

      if (ist.eq.-1) then
cTAYLOR        call drawpt (n2b,i9,0,0,.true.)
        write(*,*) 'Time step = ',i9
        write(*,*) 'Real time = ',q02,' CPU time =',etime(tt)-tim0
        call grap   (ntd,i9,q01,q02,q03,n9,ntd,ng,mc)
        call restr  (q02,i00,i0d,n5,g0,g1,g2,nst,irs)
        call fabric (1,i9,q02,n2p,dd)
        write(*,*) 'ist  =  -1; running stop at ik:',ik
        write(2,*) 'ist  =  -1; running stop at ik:',ik
        goto 1009
      endif


      g3old = g3
      goto 555

1000  continue

      end do ! ik

1009  continue

      write(*,*) 'End of Problem'
      write(*,*) 'Time step = ',i9
      write(*,*) 'Real time = ',q02,' CPU time =',etime(tt)-tim0

      write(2,*) 'End of Problem'
      write(2,*) 'Time step = ',i9
      write(2,*) 'Real time = ',q02,' CPU time =',etime(tt)-tim0
cTAYLOR      call drawpt (n2b,i9,0,-1,.true.)
cTAYLOR      call drawpt (n2b,i9,0,-2,.false.)

c     Running time

      call restr (q02,i00,i0d,n5,g0,g1,g2,nst,irs)

c     Stop timing

      call timer('dda',0,3)

c     I/O formats

2000  format (8x,'q02  no. cont',7x,
     &       'Normal F.',8x,'Shear F.',7x,'Ev2')

2001  format (' < Step >',i4,4x,'Accum. Time (sec)',f12.6,' CPU',f14.4)

2002  format (' Point Disp:',i2,3e15.5)

2003  format (1x,6e15.5)

2004  format ('  DDA - D i s c r e t e    D e f o r m a t i o n',
     &        '    A n a l y s i s'//15x,'Analysis Date: ',a//)

2005  format('  Time Step =',i7,' # Iterations =',i2,
     &       ' New dt =',1p,e11.4)

3000  format(' Input data file name:',$)

      end

      subroutine timer(sname,num,isw)

      implicit none

      character*10 name(50), sname*(*)
      real*8       ttim(50),t0(50),total
      real*4       etime,tt(2)
      integer      num,isw,i,maxn,calls(50)

c     return
      if(isw.eq.0) then
        do i = 1,50
          ttim(i)  = 0.0d0
          name(i)  = ' '
          calls(i) = 0
        end do
        maxn = 0
      elseif(isw.eq.1) then
        t0(num)    = etime(tt)
        name(num)  = sname
        calls(num) =  calls(num) + 1
      elseif(isw.eq.2) then
       ttim(num)  = ttim(num) + (dble(etime(tt)) - t0(num))
        maxn       = max(maxn,num)
      else
        write(2,2000) (i,name(i),ttim(i),calls(i),i=1,maxn)
        total = 0.0d0
        do i = 1,maxn
          total = total + ttim(i)
        end do
        write(2,2001) total
      endif

2000  format('    Number   Routine      CPU-time     No. Calls'/
     &      (i10,3x,a10,f11.4,2x,i12))

2001  format(/13x,'Total CPU:',f11.4)

      end
