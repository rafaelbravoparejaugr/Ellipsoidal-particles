
c     rem stress controlled membrane

      subroutine membrane (pi,mcode,i9,title,istep,q02)

      implicit       none

      include       'param.h'

      character      title*11
      integer        mcode,i9,istep
      real*8         pi, q02

      integer        idmr(mr),iddcon(10,mr)
      real*8         dmi(5,mr),dmlt(mr),smn(2),dinte(2,2,mr)
      real*8         dmn(2),dpl(2,mr*2),spl(2,mr/2)
      real*8         smp(2),dpli(2,mr)

      integer        i,ii,i1, j,jj, k,kk,k1, m,mm,mi1, nmemb,nmrange
      real*8         conmem1,conmem2,conmem4,conmem5
      real*8         botl,botxl,botr,botxr,areabot
      real*8         topl,topxl,topr,topxr,areatop, trapvol
      real*8         total,dangle,sangle, halfpi
      real*8         dd,dml,h,tvol,tvoll,tvolr,tvol0
      real*8         dist1,dist2,dist3,dist4,dist5,dist6,dist7,distmin
      real*8         edge,edgel,edger,topcap0,topcap2,topcap3,topcapp
      real*8         pts1in,pts1f,q02p
      real*8         pts3inx,pts3iny,pts3fx,pts3fy,ptsg3x,ptsg3y
      real*8         rmin,rmax,r0,r1,r2
      real*8         sml,smltl,smltr,ylenl,ylenr,xavgl,xavgr, xl,xr
      real*8         x0,y0,xx,yy,x1,y1,x2,y2,xd,yd, xmemr
      real*8         height,area,darea,dvol,tdvol,voidr
      real*8         astrn,vstrn,tvstrn
      real*8         norf1,sig1,sig3,sig3x,sig3l,sig3lx,sig3r,sig3rx
      real*8         xforcl,yforcl,xforcr,yforcr,partsig1,sig1adj

      include       'zbd.h'
      include       'zd.h'
      include       'zg.h'
      include       'zk0b.h'
      include       'zlp.h'
      include       'zmemb.h'
      include       'zmem2.h'
      include       'zmem3.h'
      include       'zmem4.h'
      include       'zmem5.h'
      include       'zmem6.h'
      include       'zn.h'
      include       'zpl.h'

c     md: maximum number of disks
c     mr: maximum number of disks in membrane range

c     mcode:  1=update entire membrane,0=update point loads only
c     sig1app: stress to be applied on top cap through point load
c     sig3app: stress on membranes
c     sig1int,sig3int:  initial stresses on specimen
c     itosig1,itosig3: # increments to apply additional stress on specimen
c     idmr(mr): global disk #s of disks in membrane range
c     dmi(5,mr): range #,x-coor.,y-coor.,rad.,global# of disks in membrane
c     dpli(mr,2): x and y point loads in initial membranes
c     dpl(mr*2,2): dpli arrays combined for each set of membranes
c     sangle:  angle of segment bet. disk centers in membrane
c     dangle:  angle of membrane in disk
c     iddcon(*,mr): # disks in range "in contact" with this one, global
c                   disk #s of disks "in contact"
c     xmem: x-coord. of membrane boundary - for vertical membranes
c     xmemr: x-coord. of edge of membrane range - for vertical membranes
c     ybase: y-coord. of start of membrane boundary
c     conmem1: width of membrane range=conmem1*rmax
c     conmem2:  search distance for "contacting" disks
c     conmem3: constant for distance to memb. bdry disk can be
c     conmem4: clearance between top disk in membrane & topcap
c     conmem5: min dist. to base for bottom disk in membrane=conmem5*rmin
c     rmax: largest disk radius
c     rmin: smallest disk radius
c     nmrange: # of disks in membrane search range
c     nmemb: # disks in membrane
c     ndml: # disks in left membrane
c     ndmr: # disks in right membrane
c     dist1: distance from center of disks to start of membrane
c     dist2: distance between edge of disks in membrane range
c     dist3:
c     dist4:  distance from disk in membrane to membrane bdry
c     dist5:  length used for membrane compl. check -dist. bet.
c              every other disk in membrane
c     dist6:
c     dist7: distance from top of disk to top cap
c     ip71o:  ip(7,1) before first membrane made
c     iporg(8,1): orginal point info before membrane point loads
c     topl: y-coord of top of left membrane
c     topr: y-coord of top of right membrane
c     botl: y-coord of bottom of left membrane
c     botr: y-coord of bottom of right membrane
c     mcp: ip# for control point on top cap
c     topcap: y coord. of bottom of topcap
c     topcap0: y coord. of bottom of topcap at start of test
c     topcapp:  y-coord. of topcap last timestep
c     topcap2:  disp. of topcap since last timestep
c     topcap3:  disp. of topcap since start of test
c     height0:  initial sample height
c     height: current sample height
c     astrn: axial strain of sample
c     vstrn: volumetric strain (uses aver.width)
c     vol0: initial sample volume(area)
c     vol:  current sample volume(area)
c     dvol: change in volume(area)
c     area0:  initial sample area (width)
c     area: current sample area (width)-average width
c     topnf: normal force on 2nd brdy requested (topcap)
c     topsf: shear force on 2nd brdy requested (topcap)
c     botnf: normal force on 1st brdy requested (base)
c     botsf: shear force on 1st brdy requested (base)
c     norf1: average of normal force on topcap and base
c     sig1: average vertical stress
c     sig3: average horizontal stress
c     xforcl: total of horz. point loads on disks in left memb.
c     xforcr: total of horz. point loads on disks in right memb.
c     yforcl: total of vert. point loads on disks in left memb.
c     yforcr: total of vert. point loads on disks in right memb.
c     xavgl: average x-coor of left membrane
c     xavgr: average x-ccor of right membrane
c     darea: total area of all disks
c     voidr: void ratio of specimen (uses average width of sample)
c     m01:  type of test (1:strain controlled,2:force controlled)
c     nbtopc:  bdry # of topcap for strain controlled tests
c     nbvtopc: vertex # of pt on bottom of topcap for strain cont. test
c     areatop: area at the topcap
c     areabot: area at the bottom cap

      call timer('membrane',16,1)
      conmem1 = 10.0d0
      conmem4 = 0.25d0
      conmem5 = 0.10d0
      dd      = pi/180.0d0
      halfpi  = 0.5d0*pi

c     Save ip(7,1-2) & orig. loc. of top cap,g(2,mcp)

      if (i9.eq.1) then
        do i = 1,8
          iporg(i,1) = ip(i,1)
          iporg(i,2) = ip(i,2)
        end do ! i
        if (m01.eq.1) then
          topcap0 = bd(2,k0b(1,nbtopc) + nbvtopc)
          topcap  = topcap0
        else
          topcap0 = g(2,mcp)
          topcap  = topcap0
        endif

        xcenter = (edger + edgel)/2

c     Locate topcap

      else

        topcapp = topcap
        if (m01.eq.1) then
          topcap = bd(2,k0b(1,nbtopc) + nbvtopc)
        else
          topcap = g(2,mcp)
        endif
      endif

c     Locate edge of sample
c       if mcode = 0 only update point loads

      if (mcode.ne.0) then
        edgel = d(1,1) - d(3,1)
        edger = d(1,1) + d(3,1)
        rmax  = d(3,1)
        rmin  = d(3,1)
        do i = 2,nd
          xl = d(1,i) - d(3,i)
          xr = d(1,i) + d(3,i)
          edgel = min(edgel,xl)
          edger = max(edger,xr)
          rmax  = max(rmax,d(3,i))
          rmin  = min(rmin,d(3,i))
        end do ! i
        
        conmem2 = 1.5d0*rmin
      end if
c     Zero variables

      smltl = 0.d0
      smltr = 0.d0
      ylenl = 0.d0
      ylenr = 0.d0
      xavgl = 0.d0
      xavgr = 0.d0

c     Locate disks in membrane search range

      do mm = 1,2
        if (mcode.ne.0) then
          if (mm.eq.1) then
            edge  = edgel
            xmemr = edgel + conmem1*rmax
          endif
          if (mm.eq.2) then
            edge  = edger - conmem1*rmax
            xmemr = edger
          endif
          j = 0
          do i = 1,nd
            if (d(1,i).lt.xmemr.and.d(1,i).gt.edge) then
              j       = j + 1
              idmr(j) = i
            endif
          end do ! i
          nmrange = j

          jj = 0
          dist1   = 0
          dist2   = 0
          distmin = 0

c         Zero iddcon and dmi arrays

          do i = 1,nmrange
            do j = 1,10
              iddcon(j,i) = 0
            end do ! j
            do j  =  1,5
              dmi(j,i) = 0.0d0
            end do ! j
          end do ! i

          ii = 0
          do i = 1,nmrange

c           Find initial disk in membrane

            x1 = d(1,idmr(i))
            y1 = d(2,idmr(i))
            r1 = d(3,idmr(i))

c           1st disk must be "on the base"

            if ((y1-r1-ybase).lt.(conmem5*rmin)) then
              ii = ii + 1

c             Find disk closest to edge

              if (mm.eq.1) then
                dist1 = x1 - edgel
              else
                dist1 = edger - x1
              endif
              if (ii.eq.1.or.dist1.lt.distmin) then
                distmin  = dist1
                dmi(1,1) = i
                dmi(2,1) = x1
                dmi(3,1) = y1
                dmi(4,1) = r1
                dmi(5,1) = idmr(i)
              endif
            endif

c           Find contacting disks in membrane (build iddcon(,))

            do j = i+1,nmrange
              x2    = d(1,idmr(j))
              y2    = d(2,idmr(j))
              r2    = d(3,idmr(j))
              dist2 = sqrt((x1 - x2)*(x1 - x2)
     &                   + (y1 - y2)*(y1 - y2)) - r1 - r2
              if (dist2.le.conmem2) then
                iddcon(1,i)    = iddcon(1,i) + 1
                jj             = iddcon(1,i)
                iddcon(jj+1,i) = j
                iddcon(1,j)    = iddcon(1,j) + 1
                jj             = iddcon(1,j)
                iddcon(jj+1,j) = i
              endif
            end do ! j
          end do ! i

c         Locate remaining disks in membrane
c           i = membrane disk number

          do i = 1,nmrange

c           ii = range number of ith disk in membrane

            ii = dmi(1,i)
            jj = 0
            if (iddcon(1,ii).ge.1) then
              do j = 2,iddcon(1,ii)+1
                m = idmr(iddcon(j,ii))
                if (d(2,m).gt.dmi(3,i)) then
                  jj = jj + 1
                  if (mm.eq.1) then
                    dist6 = dmi(2,i+1) - d(1,m)
                  else
                    dist6 = d(1,m) - dmi(2,i+1)
                  endif
                  if (jj.eq.1.or.dist6.gt.0.0) then
                    dmi(1,i+1) = iddcon(j,ii)
                    dmi(2,i+1) = d(1,m)
                    dmi(3,i+1) = d(2,m)
                    dmi(4,i+1) = d(3,m)
                    dmi(5,i+1) = m
                  endif
                endif
              end do ! j
            endif

            if (jj.eq.0) then

c             If there is no disk contact that is above last one
c             in membrane then find next one up

              kk      = 0
              distmin = 0.0d0
              x0      = dmi(2,i)
              y0      = dmi(3,i)
              do k = 1,nmrange
                if (d(2,idmr(k)).gt.dmi(3,i)) then
                  kk    = kk + 1
                  xx    = d(1,idmr(k))
                  yy    = d(2,idmr(k))
                  dist3 = sqrt((xx - x0)*(xx - x0)
     &                       + (yy - y0)*(yy - y0))
                  if (kk.eq.1.or.dist3.lt.distmin) then
                    distmin    = dist3
                    dmi(1,i+1) = k
                    dmi(2,i+1) = xx
                    dmi(3,i+1) = yy
                    dmi(4,i+1) = d(3,idmr(k))
                    dmi(5,i+1) = idmr(k)
                  endif
                endif
              end do ! k
              if (kk.eq.0) then

c               If didn't assign next disk (i+1) then #memb = i

                nmemb = i
                goto 799
              endif
            endif

c           Stop membrane if top cap has been reached

            dist7 = topcap - dmi(3,i+1) - dmi(4,i+1)
            if (dist7.lt.conmem4*rmin) then
              nmemb = i + 1
              goto 799
            endif
          end do ! i
799       continue

c         Remove disks from membrane that are too far from edge

          i1 = 2
          do i = 2,nmemb-1
            x0 = dmi(2,i1)
            y0 = dmi(3,i1)
            r0 = dmi(4,i1)
            x1 = dmi(2,i1-1)
            x2 = dmi(2,i1+1)
            y1 = dmi(3,i1-1)
            y2 = dmi(3,i1+1)
            r1 = dmi(4,i1-1)
            r2 = dmi(4,i1+1)
            if (mm.eq.1) then
              dist4 = (x0 - r0)
     &              - ((y0 - y1)/(y2 - y1)*((x2 - r2) - (x1 - r1))
     &              + (x1 - r1))
            else
              dist4 = (y0 - y1)/(y2 - y1)*((x2 + r2)
     &              - (x1 + r1))+(x1 + r1) - (x0 + r0)
            endif

c           Maximum distance from membrane will depend opening bet.disks

            dist5 = rmax*0.5d0
            if (dist4.gt.4*dist5) then
              do k = i1,nmemb-1
                dmi(1,k) = dmi(1,k+1)
                dmi(2,k) = dmi(2,k+1)
                dmi(3,k) = dmi(3,k+1)
                dmi(4,k) = dmi(4,k+1)
                dmi(5,k) = dmi(5,k+1)
              end do ! k
              nmemb = nmemb - 1
            else
              i1 = i1 + 1
            endif
          end do ! i

        endif

c       Calculate segments, normals, and point loads
c       Start here if updating point loads only


c       update dmi array (update disk loc.) from last update of membrane

        if (mcode.eq.0) then
          if     (mm.eq.1) then
            nmemb = ndml
          elseif (mm.eq.2) then
            nmemb = ndmr
          endif
          if (mm.eq.1) then
            do i = 1,nmemb
              dmi(1,i) = dmem(1,i)
              dmi(2,i) = d(1,int(dmem(5,i)))
              dmi(3,i) = d(2,int(dmem(5,i)))
              dmi(4,i) = dmem(4,i)
              dmi(5,i) = dmem(5,i)
            end do ! i
          else
            do i = 1,nmemb
              dmi(1,i) = dmem(1,i+ndml)
              dmi(2,i) = d(1,int(dmem(5,i+ndml)))
              dmi(3,i) = d(2,int(dmem(5,i+ndml)))
              dmi(4,i) = dmem(4,i+ndml)
              dmi(5,i) = dmem(5,i+ndml)
            end do ! i
          endif
        endif

c       Compute average sample area (width)

        do i = 1,nmemb
          if (mm.eq.1) then
            xavgl = xavgl + (dmi(2,i) - dmi(4,i))
          else
            xavgr = xavgr + (dmi(2,i) + dmi(4,i))
          endif
        end do ! i
        if (mm.eq.1) xavgl = xavgl/nmemb
        if (mm.eq.2) xavgr = xavgr/nmemb

c       Compute volume of sample using trapazoid method

        trapvol = 0.d0
        do i = 1,nmemb
          if (i.eq.1) then
            if (mm.eq.1) then
              trapvol = (xcenter - dmi(2,i) + dmi(4,i))*(dmi(3,i)
     &                 - ybase)
            else
              trapvol = (dmi(2,i) + dmi(4,i) - xcenter)*(dmi(3,i)
     &                 - ybase)
            endif
          endif
          if (mm.eq.1) then
            x1      = dmi(2,i)   - dmi(4,i)
            x2      = dmi(2,i+1) - dmi(4,i+1)
            h       = dmi(3,i+1) - dmi(3,i)
            trapvol = trapvol + (xcenter - (x1 + x2)*0.5d0)*h
          else
            x1      = dmi(2,i)   + dmi(4,i)
            x2      = dmi(2,i+1) + dmi(4,i+1)
            h       = dmi(3,i+1) - dmi(3,i)
            trapvol = trapvol + ((x1 + x2)*0.5d0 - xcenter)*h
          endif
        end do ! i
        if (mm.eq.1) then
          tvoll = trapvol
        else
          tvolr = trapvol
          tvol  = tvoll + tvolr
        endif

c       Zero arrays

        do i = 1,nmemb
          dpli(1,i)    = 0.d0
          dpli(2,i)    = 0.d0
          dmlt(i)      = 0.d0
          dinte(1,1,i) = 0.d0
          dinte(1,2,i) = 0.d0
          dinte(2,1,i) = 0.d0
          dinte(2,2,i) = 0.d0
        end do ! i
        smp(1) = 0.d0
        smp(2) = 0.d0
        dml    = 0.d0
        dmn(1) = 0.d0
        dmn(2) = 0.d0
        sml    = 0.d0

        do i = 1,nmemb

c         Calc length of segment between centers of disks i & i + 1

          if (i.ne.nmemb) then
            xx      = dmi(2,i+1) - dmi(2,i)
            yy      = dmi(3,i+1) - dmi(3,i)
            dmlt(i) = sqrt(xx*xx + yy*yy)

c           Calc angle of segment

            if (xx.ne.0.0) then
              if (xx.gt.0.0) then
                sangle = atan(yy/xx)
              else
                sangle = halfpi + atan(abs(xx)/yy)
              endif
            else
              sangle = halfpi
            endif

c           Find normal

            if (sangle.le.halfpi) then
              smp(1) = cos(sangle)
              smp(2) = sin(sangle)
            else
              smp(1) = sin(sangle - halfpi)
              smp(2) = cos(sangle - halfpi)
            endif

c           Find intersection points

            if (sangle.le.halfpi) then
              dinte(2,1,i)   = dmi(2,i) + dmi(4,i)*smp(1)
              dinte(2,2,i)   = dmi(3,i) + dmi(4,i)*smp(2)
              dinte(1,1,i+1) = dmi(2,i) + (dmlt(i) - dmi(4,i+1))*smp(1)
              dinte(1,2,i+1) = dmi(3,i) + (dmlt(i) - dmi(4,i+1))*smp(2)
            else
              dinte(2,1,i)   = dmi(2,i) - dmi(4,i)*smp(1)
              dinte(2,2,i)   = dmi(3,i) + dmi(4,i)*smp(2)
              dinte(1,1,i+1) = dmi(2,i) - (dmlt(i) - dmi(4,i+1))*smp(1)
              dinte(1,2,i+1) = dmi(3,i) + (dmlt(i) - dmi(4,i+1))*smp(2)
            endif

c           If 1st disk in membrane find bottom point

            if (i.eq.1) then
              dinte(1,1,1) = dmi(2,i)
              dinte(1,2,1) = dmi(3,i) - dmi(4,i)

c             Save location of bottom of membrane

              if (mm.eq.1) then
                botl    = dinte(1,2,1)
                botxl   = dmi(2,i) - dmi(4,i)
              else
                botr    = dinte(1,2,1)
                botxr   = dmi(2,i) + dmi(4,i)
                areabot = botxr - botxl
              endif
            endif
          else

c           If last disk in membrane find top point

            dinte(2,1,i) = dmi(2,i)
            dinte(2,2,i) = dmi(3,i) + dmi(4,i)

c           Save vertical location of top of membrane

            if (mm.eq.1) then
              topl  = dinte(2,2,i)
              dist7 = topcap - topl
              if (dist7.gt.0.0) then
                topl      = topcap
                dpli(1,i) = dpli(1,i) + sig3app*dist7
                ylenl     = ylenl + dist7
              endif

c             Save xloc of membrane for calc of area at topcap

              topxl = dmi(2,i) - dmi(4,i)
            else
              topr  = dinte(2,2,i)
              dist7 = topcap - topr
              if (dist7.gt.0.0) then
                topr      = topcap
                dpli(1,i) = dpli(1,i) - 1*sig3app*dist7
                ylenr     = ylenr + dist7
              endif
              topxr   = dmi(2,i) + dmi(4,i)
              areatop = topxr - topxl

            endif
          endif

c         Calc length of segment in disk i

          xd  = dinte(2,1,i) - dinte(1,1,i)
          yd  = dinte(2,2,i) - dinte(1,2,i)
          dml = sqrt(xd*xd + yd*yd)

c         Calc angle and normal of segment in disk i

          if (xd.ne.0.0) then
            if (xd.gt.0.0) then
              dangle = atan(yd/xd)
            else
              dangle = halfpi + atan(abs(xd)/yd)
            endif
          else
            dangle = halfpi
          endif

          if (mm.eq.1) then
            if (dangle.le.halfpi) then
              dmn(1) =  sin(dangle)
              dmn(2) = -cos(dangle)
            else
              dmn(1) =  sin(pi - dangle)
              dmn(2) =  cos(pi - dangle)
            endif
          else
            if (dangle.le.halfpi) then
              dmn(1) = -sin(dangle)
              dmn(2) =  cos(dangle)
            else
              dmn(1) = -sin(pi - dangle)
              dmn(2) = -cos(pi - dangle)
            endif
          endif
          if (mm.eq.1) then
            ylenl = ylenl + abs(dmn(1))*dml
          else
            ylenr = ylenr + abs(dmn(1))*dml
          endif

c         Calc point load for disk i

          if (xd.ge.0.0) then
            dpli(1,i) = dpli(1,i) + sig3app*dmn(1)*dml
            dpli(2,i) = dpli(2,i) + sig3app*dmn(2)*dml
          else
            dpli(1,i) = dpli(1,i) + sig3app*dmn(1)*dml
            dpli(2,i) = dpli(2,i) + sig3app*dmn(2)*dml
          endif

c         Calc open segment between disks i & i + 1

          if (i.ne.nmemb) then

c           Calc length of segment

            sml = dmlt(i) - dmi(4,i) - dmi(4,i+1)
            if (sml.le.0.0) then
              if (mm.eq.1) then
                smltl = smltl + sml
              else
                smltr = smltr + sml
              endif
            endif

            if (sml.gt.0.0) then

              if (mm.eq.1) then
                if (sangle.le.halfpi) then
                  smn(1) =  sin(sangle)
                  smn(2) = -cos(sangle)
                else
                  smn(1) =  sin(pi - sangle)
                  smn(2) =  cos(pi - sangle)
                endif
              else
                if (sangle.le.halfpi) then
                  smn(1) = -sin(sangle)
                  smn(2) =  cos(sangle)
                else
                  smn(1) = -sin(pi - sangle)
                  smn(2) = -cos(pi - sangle)
                endif
              endif
              if (mm.eq.1) then
                ylenl = ylenl + abs(smn(1))*sml
              else
                ylenr = ylenr + abs(smn(1))*sml
              endif

c             Calc point load from open segment

              spl(1,i) = sig3app*smn(1)*sml
              spl(2,i) = sig3app*smn(2)*sml

c             Apply 1/2 point load to disks above and below open segment

              dpli(1,i)   = dpli(1,i)   + spl(1,i)*0.5d0
              dpli(2,i)   = dpli(2,i)   + spl(2,i)*0.5d0
              dpli(1,i+1) = dpli(1,i+1) + spl(1,i)*0.5d0
              dpli(2,i+1) = dpli(2,i+1) + spl(2,i)*0.5d0
            endif
          endif
        end do ! i

        if (i9.eq.2500.or.i9.eq.5000) then
          write(19,*) 'disk point loads are: '
          do i = 1,nmemb
            write(19,*) dmi(5,i),dpli(1,i),dpli(2,i)
          end do ! i
        endif

        if (mm.eq.1) then
          ndml = nmemb
          do i = 1,nmemb
            do j = 1,5
              dmem(j,i) = dmi(j,i)
            end do ! j
            dpl(1,i) = dpli(1,i)
            dpl(2,i) = dpli(2,i)
          end do ! i
        else
          ndmr = nmemb
          do i = 1,nmemb
            do j = 1,5
              dmem(j,i+ndml) = dmi(j,i)
            end do ! j
            dpl(1,i+ndml) = dpli(1,i)
            dpl(2,i+ndml) = dpli(2,i)
          end do ! i
        endif
      end do ! mm

c     Update ip(7-8,1-2),nt

      ip(7,1) = iporg(7,1) + ndml + ndmr
      do i = 7,8
        ip(i,2) = ip(i-1,2) + ip(i,1)
      end do ! i
      nt = ip(8,2)

c     Update lp(,1-2) array
c     k1: beginning number of new load points

      k1  = iporg(7,2) + 1
      mi1 = lp(2,iporg(7,2))
      do i = k1,ip(7,2)
        lp(1,i) = 2
        mi1     = mi1 + lp(1,i)
        lp(2,i) = mi1
        lp(1,i) = lp(2,i) - lp(1,i) + 1
      end do ! i

c     Update pl(,1-3) array

      do i = k1,ip(7,2)
        do j = lp(1,i),lp(2,i)
          if (j.eq.lp(1,i)) then
            pl(1,j) = 0.0d0
          else
            pl(1,j) = 100.0d0
          endif
          pl(2,j) = dpl(1,i-k1+1)
          total   = pl(2,j)
          pl(3,j) = dpl(2,i-k1+1)

c*****    for testing purposes remove vertical comp. of membrane

c         pl(3,j) = 0.d0

c         Compute portion of ptloads to maintain initial stresses

          pts3inx = (sig3int/sig3app)*pl(2,j)
          pts3iny = (sig3int/sig3app)*pl(3,j)

c         Compute portion of ptloads to turn on incrementally to obtain
c         final stress state

          pts3fx = pl(2,j) - pts3inx
          pts3fy = pl(3,j) - pts3iny

c         Increase stresses incrementally from initial to final state

          if (i9.lt.itosig3) then
            ptsg3x  = pts3fx/itosig3
            ptsg3y  = pts3fy/itosig3
            pl(2,j) = pts3inx + ptsg3x*i9
            pl(3,j) = pts3iny + ptsg3y*i9
          endif

        end do ! j
      end do ! i

c     Update g(,1-3) array

      do i = k1,ip(7,2)
        g(1,i) = dmem(2,i-k1+1)
        g(2,i) = dmem(3,i-k1+1)
        g(3,i) = dmem(5,i-k1+1)
      end do ! i

c     Compute total & incr. disp. of top cap

      topcap2 = topcapp - topcap
      topcap3 = topcap0 - topcap

c     Compute height,area,vol,axial strain,sig3,sig1

      if (i9.eq.1) then
        height0 = topcap0 - ybase
        area0   = xavgr - xavgl
        vol0    = area0*height0
        tvol0   = tvol
      endif
      height = topcap - ybase
      area   = xavgr - xavgl
      vol    = area*height
      dvol   = vol0 - vol
      tdvol  = tvol0 - tvol
      astrn  = topcap3/height0
      vstrn  = dvol/vol0
      tvstrn = tdvol/tvol0
      norf1  = (topnf + botnf)/2
      sig1   = norf1/area
      xforcl = 0.0
      yforcl = 0.0
      do i = lp(1,k1),lp(2,k1+ndml-1)
        xforcl = xforcl + pl(2,i)
        yforcl = yforcl + pl(3,i)
      end do ! i
      xforcr = 0.0
      yforcr = 0.0
      do i = lp(1,k1+ndml),lp(2,k1+ndml+ndmr-1)
        xforcr = xforcr + pl(2,i)
        yforcr = yforcr + pl(3,i)
      end do ! i

      xforcl = xforcl*0.5d0
      xforcr = xforcr*0.5d0
      yforcl = yforcl*0.5d0
      yforcr = yforcr*0.5d0
      sig3l  = xforcl/(topl - botl)
      sig3r  = xforcr/(topr - botr)
      sig3lx = xforcl/ylenl
      sig3rx = xforcr/ylenr
      sig3   = (sig3l  - sig3r)*0.5d0
      sig3x  = (sig3lx - sig3rx)*0.5d0

c     Adjust load on top cap to apply sig1app
c        * for force controlled problems only
c        * assumes contol pt for point load is last one before
c          pts for membrane

      if (m01.eq.2) then
        i = iporg(7,2)
        pl(3,i) = -sig1app*area

c       Compute portion of ptloads to maintain initial stresses

        pts1in = (sig1int/sig1app)*pl(3,i)

c       Compute portion of ptloads to turn on incrementally to obtain
c       final stress state

        pts1f = pl(3,i) - pts1in

c       Increase stress incrementally from initial to final state

        if (i9.lt.itosig1) then
          partsig1 = pts1f/itosig1
          pl(3,i)  = pts1in + partsig1*i9
        endif

        sig1adj = -pl(3,i)/(area*101330.0d0)

      endif

c     Compute void ratio of specimen
c     Compute total area of all disks

      if (i9.eq.1) then
        darea = 0.d0
        do i = 1,nd
          darea = darea + pi*d(3,i)*d(3,i)
        end do ! i
      endif

c     Compute void ratio

      voidr = (vol - darea)/darea

c     Output stress strain data

      if (i9.eq.1) then
        write (15,*)  title
        write (15,*) 'Initial void ratio:   ',voidr
        write (15,*) 'Stress on membrane,sig3(atm):',sig3app/101330
        write (15,*) 'Initial sample height       :',height0
        write (15,*) 'Initial sample area(width)  :',area0
        write (15,*) 'Initial sample volume(area) :',vol0
        write (15,*) 'Initial sample volume(area) :',vol0
        write (15,2002)
      endif

      if (i9.eq.1.or.mod(i9,istep).eq.0) then
       write (15,2003) q02,astrn,vstrn,sig3/101330.,sig1/101330.,
     &                 voidr,sig1adj
      endif

c     Output info for verification purposes

      if (i9.eq.1.or.mod(i9,istep).eq.0) then
        if (i9.eq.1) then
          write (16,2000)
        endif
        write (16,2001) q02,areatop,topnf,areabot,botnf,yforcl,yforcr
      endif

      q02p = q02

      call timer('membrane',16,2)
c     Formats

2000  format('      time    top area  top force    bot area',
     &       ' bot force      mforcl     mforcr')

2001  format (1p,7e12.6)

2002  format(7x,'time',11x,'ea',11x,'ev',6x,'sigma-3',6x,
     &          'sigma-1   void ratio  sigma-1 appl')

2003  format (1p,7e13.5)

      end
