

c     rem ******************** 2000
c     rem read data files and find window limits
c     rem ******************** called by main

      subroutine datain (pi,q02,i00,i0d,n5,g0,g1,g2,ng,ntd,rmax,nst,
     & NBeta,NGamma,ji,rhof,rug,nu,us)

c     inputs -- pi
c     create -- n0,n1,n1a,n1b,n1c,n7,n6,nv,g(),c(),k0(),d(),nb,k0b(),bd()
c               i00,n5,g0,g1,g2,sp(),o0,o1,o2,o3,o4,ng,e(),w(),w5,h(),
c               fact,r0,z1,z2
c     output -- n0,n1,n1a,n1b,n7,n6,nv,g(),c(),k0(),d(),nb,k0b(),bd(),i00,n5,
c               g0,g1,g2,sp(),o0,o1,o2,o3,o4,ng,w(),h(),fact,r0,z1,z2

      implicit         none

      include         'param.h'
      include         'files.h'
      include         'zabd.h'
      include         'zb.h'
      include         'zbd.h'
      include         'zbd0.h'
      include         'zc.h'
      include         'zd.h'
      include         'zdz.h'
      include         'ze.h'
      include         'zg.h'
      include         'zh.h'
      include         'zh0.h'
      include         'zicltyp.h'
      include         'zjd.h'
      include         'zk0b.h'
      include         'zk1.h'
      include         'zk2.h'
      include         'zlp.h'
      include         'zn.h'
      include         'znp0.h'	  
      include         'zp.h'
      include         'zpl.h'
      include         'zs.h'
      include         'zscl.h'
      include         'zsp.h'
      include         'zsty.h'
      include         'zvbd.h'
      include         'zw.h'

      real*8           pi,q02,g0,g1,g2,rmax,NBeta,NGamma
      integer          i00,i0d,n5,ng,ntd,nst

      integer          i,i01,i1,i2,ib, j,j1,jl,ju, k, l0
      real*8           dist2,dll,dxx,dyy, ss,sx,sy, vv,w5
      real*8           ji,rhof,rug,nu,us

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     ms      : Maximum number of surface types

c     n0      : Number of disks + polygons
c     nd      : Number of disks
c     np      : Number of rigid polygons
c     nb      : Number of rigid boundaries
c     n1      : Number of fixed points
c     n1a     : Number of 1-d constraint points
c     n1b     : Number of spring points
c     n1c     : Number of one-direction moving constraint points
c     nbt     : Number of bolting points
c     n1x     : Number of no rotation points
c     n5      : Maximum number of time steps allowed
c     n6      : Number of load points
c     n7      : Number of measured points
c     nbv     : Number of polygon/boundary vertices in bd()
c     ng      : Control code for graphic output of contact forces - 0: once
c                 at end; others: every specified time
c     i00     : 0 for static;
c               1 for implicit dynamic;
c               2 for explicit dynamic
c     g0      : Penalty (contact stiffness) (e10 - e12)
c     g1      : Time step (in sec)
c     g2      : Ratio of maximum displacement per step (.0001 - .01)
c     o0      : Unit disk area mass
c     o3      : Unit mass force in x dirction
c     o4      : Unit mass force in y dirction
c     w5      : Average disk aera
c     g(,1-3) : Coordinates/belonged disk no. of given points
c     d()     : x,y,radius,phi,coh of disk
c     icltyp(): Cluster type that disk belongs to
c     k0b()   : Beginning registered vertex no. of boundary i in bd()
c     bd()    : x,y,length,angle,phi,coh of boundary vertices
c     bd0()   : x,y,surface type of boundary vertices at t=0
c     sp(,1-6): lx(t0),ly(c),del(delx),dely,k:ratio to g0, t0 of control points
c     e(4,2)  : Temporary use to compute window limits
c     w()     : Window limits
c     c()     : Displacement errors for fixed/spring/measured points
c     h()     : Velocity vector of each disk
c     ha()    : Acceleration vector of each disk (explicit only)
c     jd()    : Surface type of each disk
c     sty()   : Surface type #, phi, coh of each type
c     lp()    : Registered # of time-depenedent point "loads"
c     pl(,1-3): Time, Delx, Dely  of time-depenedent point "loads"
c-----[--+---------+---------+---------+---------+---------+---------+-]

c     Read data from 'ftemp'

      call timer('datain',5,1)

      open (7,file = ftemp)

      read (7,*) nd,(ip(j,1),j=1,8),np,nd0

      ip(1,2) = ip(1,1)
      do i = 2,8
        ip(i,2) = ip(i-1,2)+ip(i,1)
      end do ! i
      nt = ip(8,2)

      call exceed (nd+np,md,'md','in datain......')
      call exceed (nt   ,mp,'mp','in datain......')

      do i = 1, nt
        read (7,*) (g(j,i),j=1,3),c(1,i),c(2,i)
      end do ! i

      do i = 1, nd
        read (7,*) (d(j,i),j=1,5),jd(i),icltyp(i),k1(1,i),k1(2,i)
      end do ! i

      read (7,*) nb, nbv
      call exceed (np+nb,mp,'mp','in datain......')
      do i = 1, np+nb
        read (7,*) k0b(1,i)
      end do ! i

      do i = 1, nbv
        read (7,*) (bd0(j,i),j=1,3)
        do j = 1, 2
          bd(j,i) = bd0(j,i)
        end do ! j
      end do ! i

      read (7,*) q02
      close (7)

c     Create nd0,k2()

      nd0 = 0
      do i = 1, nd
        if (k1(2,i).eq.1) then
          nd0 = nd0+1
          k2(1,nd0) = i
        endif
        if (k1(2,i+1).eq.1.or.i.eq.nd) k2(2,nd0) = k1(2,i)
      end do ! i

      n0 = nd0+np

      call exceed (n0,md,'md','in datain......')

c     Compute dz(,2),h0(,1-3)

      do i = 1, nd0
        i1 = k2(1,i)
        if (k2(2,i).eq.1) then
          dz(1,i) = d(1,i1)
          dz(2,i) = d(2,i1)
          dz(5,i) = d(5,i1)		  
          s(1)    = d(3,i1)*d(3,i1)
c area e inercia de la Elipse	
          h0(1,i) = 4.0d0/3.0d0*pi*d(3,i1)*d(4,i1)*d(4,i1)
          h0(2,i) = 4.0d0/3.0d0*pi*d(3,i1)*d(4,i1)*d(4,i1)
          h0(3,i) = h0(2,i)/5*(d(3,i1)**2+d(4,i1)**2)
          hm(1,i) = 4.0d0/3.0d0*pi*d(3,i1)*d(4,i1)*d(4,i1)
          hm(2,i) = 4.0d0/3.0d0*pi*d(3,i1)*d(4,i1)*d(4,i1) 
          hm(3,i) = hm(2,i)/5*(d(3,i1)**2+d(4,i1)**2)
        else
        endif
      end do ! i

c     Set np0() for polygons

      do i = 1, np
        i1 = k0b(1,i)
        i2 = bd(1,i1)
        jl = i1+1
        ju = i1+i2+1
        do j = jl, ju
          np0(j) = i+nd
        end do ! j
      end do ! i

c     Read data from 'fpara'

      open (7,file = fpara)

c     i00: 0: statics; 1: implicit dynamic; 2: explicit dynamics

      read (7,*) i00

c     Newmark Parameters NBeta and NGamma
        
      read (7,*) NBeta
 
      read (7,*) NGamma

      read (7,*) ji !factor de van rijn

      read (7,*) rhof !densidad del fluido

      read (7,*) rug !rugosidad en tanto por uno

      read (7,*) nu !viscosidad

      read (7,*) us !velocidad u estrella


c     n5:  maximum number of time steps (5 - 1000)'
      read (7,*) n5

c     Penalty g0: (e10 - e12)'
      read (7,*) g0

c     Time step g1
      read (7,*) g1

c     g2*(w(4)-w(2)): maximum step vertex displacement

c     Ratio of max. displacement g2: (.0001 - .01)'
      read (7,*) g2

c     Enter time-sequence/constant terms for control points (excluding
c       measured points)

      write(19,*) 'ip(7,2)',ip(7,2)

c     Time-sequence for outputs

      i1 = 0
      do i = 1, ip(7,2)
        if (i.gt.ip(3,2).and.i.le.ip(5,2)) then
          i2 = mod(i-ip(3,2),2)
          if (i2.eq.0) goto 125
        endif
        read (7,*) lp(1,i)
        i1      = i1+lp(1,i)
        lp(2,i) = i1
        lp(1,i) = lp(2,i)-lp(1,i)+1

125     continue

      end do ! i

       write(19,*) 'lp(,1):',(lp(1,k),k=1,ip(7,2))
       write(19,*) 'lp(,2):',(lp(2,k),k=1,ip(7,2))

c     Time-sequence terms for control points

      do i = 1, ip(7,2)
         write(19,*) 'control point:',i
        if (i.gt.ip(3,2).and.i.le.ip(5,2)) then
          i2 = mod(i-ip(3,2),2)
          if (i2.eq.0) goto 130
        endif
        do j = lp(1,i), lp(2,i)
          if (i.le.ip(1,2).or.(i.gt.ip(3,2).and.i.le.ip(4,2)).or.
     $        i.gt.ip(6,2)) then

c           Fixed points/1st: time, delx, dely, k
c             time, delx, dely, k (=100)

            if (i.le.ip(1,2).and.j.eq.lp(1,i)) then
              read (7,*) (pl(k,j),k=1,3),sp(5,i)

c           Pin points/1st: time, delx, dely, t0, coh, k
c             time, delx, dely, t0, coh, k (=100) -'

            elseif (i.gt.ip(3,2).and.i.le.ip(4,2).and.j.eq.lp(1,i)) then
              read (7,*) (pl(k,j),k=1,3),sp(1,i),sp(2,i),sp(5,i)

c           Fixed-pin/others;load points: time, delx;px, dely;py
c             time, delx;px, dely;py -'

            else
              read (7,*) (pl(k,j),k=1,3)
            endif
          else

c           1-d cos.;spring points/1st: time, del, lx, ly, k
c             time, del, lx, ly, k -'

            if (i.gt.ip(1,2).and.i.le.ip(3,2).and.j.eq.lp(1,i)) then
              read (7,*) (pl(k,j),k=1,2),sp(1,i),sp(2,i),sp(5,i)
               write(19,*) 'pl(1-2),sp(1,2,5):',(pl(k,j),k=1,2),sp(1,i),
     $                    sp(2,i),sp(5,i)
c           Bolting;x-rot. points/1st: time, del;theta, t0, k;k_theta
c             time, del;theta, t0, k;k_theta(=100) -'

            elseif (i.gt.ip(4,2).and.i.le.ip(6,2).and.j.eq.lp(1,i)) then
              read (7,*) (pl(k,j),k=1,2),sp(6,i),sp(5,i)
               write(19,*) 'pl(1-2),sp(6,5):',(pl(k,j),k=1,2),sp(6,i),
     $          sp(5,i)
c           Points/others: time, del;theta
c             time, del;theta -'
            else
              read (7,*) (pl(k,j),k=1,2)
               write(19,*) 'pl(1-2):',(pl(k,j),k=1,2)
            endif
          endif
        end do ! j

130     continue
      end do ! i

c     Compute initial (lx,ly) of bolting points

      jl = ip(4,2)+1
      do i = jl, ip(5,2), 2
        dxx = g(1,i+1)-g(1,i)
        dyy = g(2,i+1)-g(2,i)
        dll = sqrt(dxx*dxx+dyy*dyy)
        if (dll.eq.0.d0) then
          sp(1,i) = 1.d0
          sp(2,i) = 0.d0
        else
          sp(1,i) = dxx/dll
          sp(2,i) = dyy/dll
        endif
      end do ! i

c     Unit mass force      x       y - o1, o2'
c     Unit mass force      x       y - o3, o4'

      read (7,*) o1,o2
      read (7,*) o3,o4

c     Control code of o0, i0d: 0--all; 1--each disk
c       All disks/polygons, enter o0(i)

      read (7,*) i0d
      if (i0d.eq.0) then
        read (7,*) o0(1)
        do i = 1, n0
          o0(i) = o0(1)
        end do ! i
      else
        do i = 1, n0
          read (7,*) o0(i)
        end do ! i
      endif
      do i = 1, nd0
        hm(1,i) = o0(i)*hm(1,i)
        hm(2,i) = o0(i)*hm(2,i)
        hm(3,i) = o0(i)*hm(3,i)
        h0(1,i) = o0(i)*h0(1,i)
        h0(2,i) = o0(i)*h0(2,i)
        h0(3,i) = o0(i)*h0(3,i)
      end do ! i

c     Input properties of surface types
c       number of surface types given < ms
c       surface type no., phi, coh of surface type

      read (7,*) nst
      do i = 1, nst
        read (7,*) (sty(j,i),j=1,3)
      end do ! i

      do i = 1, nd
        do j = 1, nst
          j1 = sty(1,j)
          if (jd(i).eq.j1) then
            d(6,i) = sty(2,j1)
            d(7,i) = sty(3,j1)
          endif
        end do ! j
      end do ! i

      write (4,*) '  Friction angle is:  ',d(6,1)
      write (*,*) '  Friction angle is:  ',d(6,1)

      do i = 1, nbv
        do j = 1, nst
          j1 = sty(1,j)
          if (bd0(3,i).eq.j1) then
            bd(5,i) = sty(2,j1)
            bd(6,i) = sty(3,j1)
          endif
        end do ! j
      end do ! i

      if (i00.gt.0) then

c       Zero explicit acceleration

        if(i00.eq.2) then
          do i = 1, n0
            do j = 1, 3
               ha(j,i) = 0.0d0
            end do ! j
          end do ! i
        endif ! i00 .eq. 2

c       Input control code of h(): -1, 0 or +number: '
c         i01: -1 for all h()=0; 0 for all same h(); +number for input h()
c               of selected disks/polygons

        read (7,*) i01
        if (i01.ne.0) then
          do i = 1, n0
            do j = 1,3
               h(j,i) = 0.d0
            end do ! j
          end do ! i

c         Input ib,h() for selected disk/polygons

          if (i01.gt.0) then
            do i = 1, i01
              read (7,*) ib,(h(j,ib),j = 1,3)
		do j = 	1,3
		h(j,i) = h(j,i)		  
		end do ! j
            end do ! i
          endif
        else
          read (7,*) (h(j,1),j = 1,3)
          do i = 2, n0
            do j = 1, 3
              h(j,i) = h(j,1)		  
            end do ! i
          end do ! i
        endif
      endif

c     Input types & inform. of moving boundaries
c       control code of moving bound.: -1, 0 or +number: '
c       i01    : -1 all bd=fixed; 0 all same moving; +number input
c                 moving type of selected disks
c       k0b(,2): 0 - fixed; 1 - translation; 2 - rotation w.r.t (x0,y0);
c                3 - translation/rotation
c       k0b(,3): number of velocity (rotation) or accelration (r. a.)
c                variations in t domain (+: vol(rot); -:acc(r. a.))

      read (7,*) i01
      if (i01.ne.0) then
        do i = np+1, np+nb
          k0b(2,i) = 0
        end do ! i
        if (i01.gt.0) then
          do i = 1, i01
            read (7,*) ib,k0b(2,ib+np)
            if (k0b(2,ib+np).eq.1) then
              read (7,*) k0b(3,ib+np)
              if (k0b(3,ib+np).gt.0) then

c               vbd(,,1-3): vx, vy, ttermined
                do j = 1, k0b(3,ib+np)
                  read (7,*) (vbd(j,k,ib+np),k=1,3)
                end do ! j
              else

c               vbd(,1,1-2): init. vx, vy
c               abd(,,1-3) : ax, ay, ttermined

                read (7,*) vbd(1,1,ib+np),vbd(1,2,ib+np)
                do j = 1, -k0b(3,ib+np)
                  read (7,*) (abd(j,k,ib+np),k=1,3)
                end do ! j
              endif
            else if (k0b(2,ib+np).eq.2) then

c             vbd(,0,1-2): rotation (x,y)

              read (7,*) k0b(3,ib+np)
              read (7,*) vbd(0,1,ib+np),vbd(0,2,ib+np)
              if (k0b(3,ib+np).gt.0) then

c               vbd(,,1-3): w, ttermined

                do j = 1, k0b(3,ib+np)
                  read (7,*) vbd(j,1,ib+np),vbd(j,3,ib+np)
                end do ! j
              else

c               vbd(,1,1-2): init. vx, vy
c               abd(,,1-3): alp, ttermined

                read (7,*) vbd(1,1,ib+np)
                do j = 1, -k0b(3,ib+np)
                  read (7,*) abd(j,1,ib+np),abd(j,3,ib+np)
                end do ! j
              endif
            endif
          end do ! i
        endif
      else
c       Enter k0b(,3): no. of vel/acc. vari. - '
        read (7,*) k0b(2,np+1)
        if (k0b(2,np+1).eq.1) then

c         velocity profile for all bound:'

          read (7,*) k0b(3,np+1)
          if (k0b(3,np+1).gt.0) then

c           vbd(,,1-3): vx, vy, ttermined

            do j = 1, k0b(3,np+1)
              read (7,*) (vbd(j,k,np+1),k=1,3)
            end do ! j
          else

c           vbd(,1,1-2): init. vx, vy
c           abd(,,1-3): ax, ay, ttermined

            read (7,*) vbd(1,1,np+1),vbd(1,2,np+1)
            do j = 1, -k0b(3,np+1)
              read (7,*) (abd(j,k,np+1),k=1,3)
            end do ! j
          endif

c       Enter k0b(,3): no. of rot/r.a. vari. - '

        else if (k0b(2,np+1).eq.2) then

c         vbd(,0,1-2): rotation (x,y)

          read (7,*) k0b(3,np+1)
          read (7,*) vbd(0,1,np+1),vbd(0,2,np+1)
          if (k0b(3,np+1).gt.0) then

c           vbd(,,1-3): w, ttermined

            do j = 1, k0b(3,np+1)
              read (7,*) vbd(j,1,np+1),vbd(j,2,np+1)
            end do ! j
          else

c           vbd(,1,1-2): init. vx, vy
c           abd(,,1-3): alp, ttermined

            read (7,*) vbd(1,1,np+1),vbd(1,2,np+1)
            do j = 1, -k0b(3,np+1)
              read (7,*) abd(j,1,np+1),abd(j,2,np+1)
            end do ! j
          endif
        endif
        do i = np+1, np+nb
          k0b(2,i+np) = k0b(2,np+1)
          k0b(3,i+np) = k0b(3,np+1)
          do j = 1, abs(k0b(3,np+1))
            do k = 0,3
              vbd(j,k,i+np) = vbd(j,k,np+1)
              if (k.ne.0) then
                abd(j,k,i+np) = abd(j,k,np+1)
              endif
            end do ! k
          end do ! j
        end do ! i
      endif
      close (7)

c     Input graphic output control codes

      ng = 1

c     Compute window limits

      e(1,1) = d(1,1)-d(3,1)
      e(1,2) = d(1,1)+d(3,1)
      e(2,1) = d(2,1)-d(3,1)
      e(2,2) = d(2,1)+d(3,1)
      do i = 2, nd
        e(1,1) = min(e(1,1),d(1,i)-d(3,i))
        e(1,2) = max(e(1,2),d(1,i)+d(3,i))
        e(2,1) = min(e(2,1),d(2,i)-d(3,i))
        e(2,2) = max(e(2,2),d(2,i)+d(3,i))
      end do ! i
      do i = 1, np+nb
        i1 = k0b(1,i)
        j1 = bd(1,i1)
        do j = i1+1, i1+j1
          e(1,1) = min(e(1,1),bd(1,j))
          e(1,2) = max(e(1,2),bd(1,j))
          e(2,1) = min(e(2,1),bd(2,j))
          e(2,2) = max(e(2,2),bd(2,j))
        end do ! j
      end do ! i

      e(3,1) = (e(1,2)-e(1,1))*0.5d0
      e(3,2) = (e(2,2)-e(2,1))*0.5d0
      e(4,1) = (e(1,1)+e(1,2))*0.5d0
      e(4,2) = (e(2,1)+e(2,2))*0.5d0

c     vv : 1.0--1.3; ratio of window

      vv = 1.03d0+10.d0*g2/.8d0
      if (e(3,1).gt.1.3d0*e(3,2)) then
        w(1) = e(4,1)-e(3,1)*vv
        w(3) = e(4,1)+e(3,1)*vv
        w(2) = e(4,2)-e(3,1)/1.3d0*vv
        w(4) = e(4,2)+e(3,1)/1.3d0*vv
      else
        w(1) = e(4,1)-e(3,2)*1.3d0*vv
        w(3) = e(4,1)+e(3,2)*1.3d0*vv
        w(2) = e(4,2)-e(3,2)*vv
        w(4) = e(4,2)+e(3,2)*vv
      endif

      w(5) = (w(1)+w(3))*0.5d0
      w(6) = (w(2)+w(4))*0.5d0
      w(7) = (w(4)-w(2))*0.5d0

      w(1) = w(1) - 0.65*w(7)
      w(2) = w(2) - 0.20*w(7)
      w(4) = w(4) + 0.30*w(7)

c     Compute largest disk radius

      rmax = 0.d0
      do i = 1, nd
        rmax = max(rmax,d(3,i))
      end do ! i
      if (nd0.eq.0) rmax = 1.d0

c     Compute average group disk area

      w5 = 0.d0
      do i = 1, n0
        w5 = w5+h0(1,i)
      end do ! i
      w5 = w5/n0

      write(*,*) '  Average area =', w5

c     Output into 'result'

      write (2,*) 'Number of rigid disks                   :',nd
      write (2,*) 'Number of rigid group disks             :',nd0
      write (2,*) 'Number of rigid polygons                :',np
      write (2,*) 'Number of rigid group disks and polygons:',n0
      write (2,*) 'Number of rigid boundaries              :',nb
      write (2,*) 'Beta parameter	                   :',NBeta
      write (2,*) 'Gamma parameter	                   :',NGamma
      write (2,*) 'Maximum number of steps allowed         :',n5
      write (2,*) 'Contact stiffness,         g0           :',g0
      write (2,*) 'Time step assigned,        g1           :',g1
      write (2,*) 'Max. disp. ratio per step, g2           :',g2
      write (2,*) 'Number of fixed points                  :',ip(1,1)
      write (2,*) 'Number of 1-d constraint points         :',ip(2,1)
      write (2,*) 'Number of spring points                 :',ip(3,1)
      write (2,*) 'Number of pin constraint pair points    :',ip(4,1)
      write (2,*) 'Number of bolting points                :',ip(5,1)
      write (2,*) 'Number of no-rotation constraint points :',ip(6,1)
      write (2,*) 'Number of loading points                :',ip(7,1)
      write (2,*) 'Number of measured points               :',ip(8,1)
      write (2,*) 'Uniform unit disk mass of 1st disk      :',o0(1)
      write (2,*) 'x   y   of volume force                 :',o1,o2
      write (2,*) 'x   y   of mass force                   :',o3,o4

c     Compute other scaled factors

      fact = abs(3.d0/w(7))
      r0   = 0.1d0/fact
      z1   = 0.8d0/fact
      l0   = int(log10(z1))
      if (z1.ge.1.d0) then
        z1 = int(z1/10.d0**l0)
        z2 = float(l0)
      else
        z1 = int(z1/10.d0**(l0-1))
        z2 = float(l0-1)
      endif

c     Output initial data into 'grf1'

      write (3,*) nd+np,nd,np,nd0,ip(1,1),ip(2,1)+ip(3,1)+ip(4,1),
     &            ip(5,1),ip(6,1),ip(7,1),ip(8,1),nb,nbv,ntd,ng
      write (3,*) w(5),w(6),w(7),fact,r0,z1,z2

      do i = 1, np+nb
        write (3,*) k0b(1,i)
      end do ! i

c     Initialize solution at i9=0

      do i = 1, n0
        do j = 1, 3
          b(j,i) = 0.d0
        end do ! j
      end do ! i
      call timer('datain',5,2)

      end
