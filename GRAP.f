
c     rem ******************** 2200a
c     rem and save d() into #2 #3 and compute contact forces
c     rem ******************** called by main

      subroutine grap (ik,i9,q01,q02,q03,n9,ntd,ng,mc)

c     inputs -- n0,b(),d(),dz(),g(),c(),
c               i9,q01-2,ni,ng,n5,g0,dd,dm()
c     create --
c     output --

      implicit         none

      integer          ik,i9,ntd,ng,mc
      real*8           q01,q02,q03

      integer          i,i1, j, n9
      real*8           uu,vv

      include         'param.h'
      include         'zbd.h'
      include         'zc.h'
      include         'zcf.h'
      include         'zd.h'
      include         'zdz.h'
      include         'zg.h'
      include         'zh.h'
      include         'zicltyp.h'
      include         'zk1.h'
      include         'zn.h'
      include         'zrd.h'
      include         'zrotat.h'
      include         'zw.h'

c     mp: maximum number of input points/boundary vertices
c     md: maximum number of disks
c     mi: maximum number of invasions

c     n0: number of disks
c     ng: control code for graphic output of contact forces - 0: once in the
c         end; others: every specified time
c     n5: maximum number of time steps allowed
c     mc: number of contact forces at end of current step
c     d(): coordinates of disk center, their radia/friction angle/accumulated
c          rotation/radia*
c     dz(): coordinates of cluster centers
c     w(): window limits
c     g(): information of given points
c     c(): displacement errors for fixed/measured points
c     q01-2: max. disp. ratio and actual time step/accum. time at current step
c     mm(,1-4): mm(,1-3) are vertex numbers of real contact and mm(,4) contact
c               type at end of current step
c     cf(): contact forces - normal and shear forces
c     accrot(): accumulated rotations of particles
c     save d(),g(),q01-2,mc,mm(),cf() into 'grf1' at epecified time

      call timer('grap',10,1)
      do i = 1, nbv
        write (3,*) (bd(j,i),j=1,2)
      end do ! i

      if (i9.eq.0) then
        do i = 1, nd
          write (3,*) (d(j,i),j=1,3),icltyp(i)
        end do ! i
      else
        do i = 1, nd
          write (3,*) (d(j,i),j=1,3),icltyp(i)
        end do ! i
      endif

      write (3,*) nt
      do i = 1, nt
        write (3,*) g(1,i),g(2,i),g(3,i),c(1,i),c(2,i)
      end do ! i

      do i = 1,nd0
        write (3,*) dz(1,i),dz(2,i)
      end do ! i

      if (i9.ne.0) then

        write (3,*) i9,q01,q02,q03

c       Output results

        write (2,*) '++++++++++ Step ++++++++++',i9
        write (2,*) 'Total iterations            =',n9
        write (2,*) 'Accumulated time (sec)      =',q02
        write (2,*) 'Ratio of max. disp. to given=',q01
        write (2,*) 'Max. step rotation          =',q03

c       Output velocity vectors of nd disks + np polygons

        if (ng.ne.0 .or. ik.eq.ntd) then

          do i = 1, nd
            i1 = k1(1,i)
            uu = h(1,i1) - (d(2,i) - dz(2,i1))*h(3,i1)
            vv = h(2,i1) + (d(1,i) - dz(1,i1))*h(3,i1)
            write (3,*) uu,vv,h(3,i1)
          end do ! i

          do i = nd0+1, n0
            write (3,*) h(1,i),h(2,i),h(3,i)
          end do ! i

          write (3,*) mc
          do i = 1, mc
            write (3,*) (rd(j,1,i),rd(j,2,i),j=1,3)
            write (3,*) (cf(j,i),j=1,3)
          end do ! i

c         Output accumulated rotations of particles

          do i = 1,nd0
             write (3,*) accrot(i)
          end do ! i

        endif

      endif

      call timer('grap',10,2)
      end
