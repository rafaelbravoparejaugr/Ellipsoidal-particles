
c     rem ******************** xxxx
c     rem update moving boundaries and compute sugg. g3
c     rem ******************** called by main

      subroutine moving (q02,d0,g3)

c     inputs -- pi
c     create -- n0,n1,n1a,n1b,n1c,n7,n6,nv,g(),c(),k0(),d(),nb,k0b(),bd()
c               i00,n5,g0,g1,g2,sp(),o0,o1,o2,o3,o4,ng,xs,e(),w(),w5,h(),
c               fact,r0,z1,z2
c     output -- n0,n1,n1a,n1b,n7,n6,nv,g(),c(),k0(),d(),nb,k0b(),bd(),i00,n5,
c               g0,g1,g2,sp(),o0,o1,o2,o3,o4,ng,xs,w(),h(),fact,r0,z1,z2
c     mp: maximum number of input points/boundary vertices

      implicit     none

      include     'param.h'

      real*8       q02,d0,g3

      integer      i,i1,i2, j
      real*8       dt,dx,dy, cc,ss, dmax,dtot, g3p1,g3p2

      include     'zabd.h'
      include     'zbd.h'
      include     'zk0b.h'
      include     'zn.h'
      include     'zvbd.h'
      include     'zxbd.h'

      call timer('moving',17,1)

101   dmax = 0.d0
      do i = np+1, np+nb
        i1 = k0b(1,i)
        i2 = bd(1,i1)
        if (k0b(2,i).eq.0) goto 100
        if (k0b(2,i).eq.1) then

c         For constant vel.

          if (k0b(3,i).gt.0) then
            do j = 1, k0b(3,i)
              if (q02.le.vbd(j,3,i).and.q02+g3.le.vbd(j,3,i)) then
                dx = vbd(j,1,i)*g3
                dy = vbd(j,2,i)*g3
                goto 111
              endif
              if (q02.le.vbd(j,3,i).and.q02+g3.gt.vbd(j,3,i)) then
                g3p1 = vbd(j,3,i) - q02
                g3p2 = q02 + g3 - vbd(j,3,i)
                dx = vbd(j,1,i)*g3p1 + vbd(j+1,1,i)*g3p2
                dy = vbd(j,2,i)*g3p1 + vbd(j+1,2,i)*g3p2
                goto 111
              endif
            end do ! j

c         For constant acc.

          else
            do j = 1, -k0b(3,i)
              if (q02.le.abd(j,3,i).and.q02+g3.le.abd(j,3,i)) then
                dx = (vbd(1,1,i) + abd(j,1,i)*g3*0.5d0)*g3
                dy = (vbd(1,2,i) + abd(j,2,i)*g3*0.5d0)*g3
                vbd(1,1,i) = vbd(1,1,i) + abd(j,1,i)*g3
                vbd(1,2,i) = vbd(1,2,i) + abd(j,2,i)*g3
                goto 111
              endif
              if (q02.le.abd(j,3,i).and.q02+g3.gt.abd(j,3,i)) then
                g3p1 = abd(j,3,i) - q02
                g3p2 = q02 + g3 - abd(j,3,i)
                dx = (vbd(1,1,i) + abd(j,1,i)*g3p1*0.5d0)*g3p1
                dy = (vbd(1,2,i) + abd(j,2,i)*g3p1*0.5d0)*g3p1
                vbd(1,1,i) = vbd(1,1,i) + abd(j,1,i)*g3p1
                vbd(1,2,i) = vbd(1,2,i) + abd(j,2,i)*g3p1
                dx = dx + (vbd(1,1,i) + abd(j+1,1,i)*g3p2*0.5d0)*g3p2
                dy = dy + (vbd(1,2,i) + abd(j+1,2,i)*g3p2*0.5d0)*g3p2
                vbd(1,1,i) = vbd(1,1,i) + abd(j+1,1,i)*g3p2
                vbd(1,2,i) = vbd(1,2,i) + abd(j+1,2,i)*g3p2
                goto 111
              endif
            end do ! j
          endif

111       dtot = dx*dx + dy*dy
          dmax = max(dmax,dtot)

          do j = i1+1, i1+i2
            xbd(1,j) = dx
            xbd(2,j) = dy
          end do ! j

        elseif (k0b(2,i).eq.2) then

c         For constant angular vel.

          if (k0b(3,i).gt.0) then
            do j = 1, k0b(3,i)
              if (q02.le.vbd(j,3,i).and.q02+g3.le.vbd(j,3,i)) then
                dt = vbd(j,1,i)*g3
                goto 141
              endif
              if (q02.le.vbd(j,3,i).and.q02+g3.gt.vbd(j,3,i)) then
                g3p1 = vbd(j,3,i) - q02
                g3p2 = q02 + g3 - vbd(j,3,i)
                dt = vbd(j,1,i)*g3p1 + vbd(j+1,1,i)*g3p2
                goto 141
              endif
            end do ! j

c         For constant angular acc.

          else
            do j = 1, -k0b(3,i)
              if (q02.le.abd(j,3,i).and.q02+g3.le.abd(j,3,i)) then
                dt = vbd(1,1,i)*g3 + abd(j,1,i)*g3*g3*0.5d0
                vbd(1,1,i) = vbd(1,1,i) + abd(j,1,i)*g3
                goto 141
              endif
              if (q02.le.abd(j,3,i).and.q02+g3.gt.abd(j,3,i)) then
                g3p1 = abd(j,3,i) - q02
                g3p2 = q02 + g3 - abd(j,3,i)
                dt = (vbd(1,1,i) + abd(j,1,i)*g3p1*0.5d0)*g3p1
                vbd(1,1,i) = vbd(1,1,i) + abd(j,1,i)*g3p1
                dt = dt + (vbd(1,1,i) + abd(j+1,1,i)*g3p2*0.5d0)*g3p2
                vbd(1,1,i) = vbd(1,1,i) + abd(j+1,1,i)*g3p2
                goto 141
              endif
            end do ! j
          endif

141       cc = cos(dt)
          ss = sin(dt)
          do j = i1+1, i1+i2
            dx = bd(1,j) - vbd(0,1,i)
            dy = bd(2,j) - vbd(0,2,i)
            xbd(1,j) = vbd(0,1,i) + dx*cc - dy*ss - bd(1,j)
            xbd(2,j) = vbd(0,2,i) + dx*ss + dy*cc - bd(2,j)
            dtot = xbd(1,j)*xbd(1,j) + xbd(2,j)*xbd(2,j)
            dmax = max(dmax,dtot)
          end do ! j

        endif

100     continue

      end do ! i

      if (dmax.gt.0.25d0*d0*d0) then
        g3 = g3*0.5d0
        goto 101
      endif

      call timer('moving',17,2)

      end
