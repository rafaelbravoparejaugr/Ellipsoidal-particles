
c     rem ********************
c     rem output "temp" and input file for restart analysis
c     rem ********************

      subroutine restr (q02,i00,i0d,n5,g0,g1,g2,nst,irs)
	  
c     inputs -- n0,n1,n1a,n1b,n1c,n7,n6,nv,g(),c(),k0(),d(),i00,n5,g0,g1,g2,
c               sp(),o0,o1,o2,o3,o4,h()
c     create -- none
c     output -- n0,n1,n1a,n1b,n1c,n7,n6,nv,g(),c(),k0(),d(),i00,n5,g0,g1,g2,
c               sp(),o0,o1,o2,o3,o4,h()

      implicit         none

      include         'param.h'
      include         'zbd.h'
      include         'zbd0.h'
      include         'zc.h'
      include         'zchl.h'
      include         'zd.h'
      include         'zg.h'
      include         'zh.h'
      include         'zicltyp.h'
      include         'zjd.h'
      include         'zk0b.h'
      include         'zk1.h'
      include         'zlp.h'
      include         'zmem2.h'  
      include         'zn.h'
      include         'znam.h'
      include         'znc.h'
      include         'zncl.h'
      include         'zp.h'
      include         'zpl.h'
      include         'zsp.h'
      include         'zsty.h'

      integer          i00,i0d,n5,nst,irs
      real*8           q02,g0,g1,g2

      integer          i,i2, j, k, ntorg

      character        rsfname*14, dfrs*14

c     mp: maximum number of input points/boundary vertices
c     md: maximum number of disks

c     7/19/95 removed interactive parts, and do not reset c(,1-2),q02, h()
c     *8/30/95 added output of nt and ip(1-*8,1) minus membrane point loads
c     'temp': input geometry file for next restart analysis

c     Output data to 'temp'

      call timer('restr',25,1)
      irs=irs+1

c     Obtain restart file names

      if (irs.eq.1) then
         call pschar
      endif
      rsfname=rsname(1:nc)//chl(irs)(1:ncl(irs))//'.rs'
      open (7,file=rsfname)
      write (7,*) nd,(iporg(j,1),j=1,8),np,nd0

c     Calculate nt from control points minus membrane point loads
c        ntorg=total # control pts minus membrane pt.loads

      ntorg=0
      do i=1,8
        ntorg=ntorg+iporg(i,1)
      end do ! i

      do i=1,ntorg
        c(1,i)=0.d0
        c(2,i)=0.d0
        write (7,*) (g(j,i),j=1,3),c(1,i),c(2,i)
      end do ! i

      do i=1, nd
        write (7,*) (d(j,i),j=1,5),jd(i),icltyp(i),k1(1,i),k1(2,i)
      end do ! i
      write (7,*) nb,nbv

      do i=1, nb+np
        write (7,*) k0b(1,i)
      end do ! i

      do i=1, nbv
        write (7,*) (bd(j,i),j=1,2),bd0(3,i)
      end do ! i

      write (7,*) q02
      close (7)

      dfrs = rsname(1:nc)//chl(irs)(1:ncl(irs))//'.df'

c     comment out making of parameter files for big run w/lots output
c     commenting everything between here and where file is closed

      open (7,file=dfrs)
      write (7,*) i00
      write (7,*) n5
      write (7,*) g0
      write (7,*) g1
      write (7,*) g2

      do i=1, ip(7,2)
        if (i.gt.ip(3,2).and.i.le.ip(5,2)) then
          i2=mod(i-ip(3,2),2)
          if (i2.eq.0) goto 125
        endif
        write (7,*) lp(2,i)-lp(1,i)+1
125     continue
      end do ! i

      do i=1, ip(7,2)

        if (i.gt.ip(3,2).and.i.le.ip(5,2)) then
          i2=mod(i-ip(3,2),2)
          if (i2.eq.0) goto 130
        endif

        do j=lp(1,i), lp(2,i)

          if (i.le.ip(1,2).or.(i.gt.ip(3,2).and.i.le.ip(4,2)).or.
     &        i.gt.ip(6,2)) then
            if (i.le.ip(1,2).and.j.eq.lp(1,i)) then
              write (7,*) (pl(k,j),k=1,3),sp(5,i)
            elseif (i.gt.ip(3,2).and.i.le.ip(4,2).and.j.eq.lp(1,i)) then
              write (7,*) (pl(k,j),k=1,3),sp(1,i),sp(2,i),sp(5,i)
            else
              write (7,*) (pl(k,j),k=1,3)
            endif
          else
            if (i.gt.ip(1,2).and.i.le.ip(3,2).and.j.eq.lp(1,i)) then
              write (7,*) (pl(k,j),k=1,2),sp(1,i),sp(2,i),sp(5,i)
            elseif (i.gt.ip(4,2).and.i.le.ip(6,2).and.j.eq.lp(1,i)) then
              write (7,*) (pl(k,j),k=1,2),sp(6,i),sp(5,i)
            else
              write (7,*) (pl(k,j),k=1,2)
            endif
          endif
        end do ! j

130     continue
      end do ! i

      write (7,*) o1,o2
      write (7,*) o3,o4
      write (7,*) i0d
      if (i0d.eq.0) then
        write (7,*) o0(1)
      else
        do i=1, n0
          write (7,*) o0(i)
        end do ! i
      endif
      write (7,*) nst
      do i=1, nst
        write (7,*) (sty(j,i),j=1,3)
      end do ! i
      write (7,*) n0
      do i=1, n0
        write (7,*) i,(h(j,i),j=1,3)
      end do ! i

      write (7,*) i00
      close (7)

      call timer('restr',25,2)
      end
