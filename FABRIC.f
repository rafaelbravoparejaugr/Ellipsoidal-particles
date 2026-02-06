
c     rem compute fabric data: average coordinate no.&contact frequency

      subroutine fabric (i0,i9,q02,n2b,dd)

      implicit     none

      integer      i0,i9,n2b
      real*8       q02,dd

      integer      i,j, ni,nn,nnn
      real*8       acn,a0,dag

      include     'param.h'
      include     'zdm.h'
      include     'zn.h'
      include     'zo.h'

      real*8       agi(0:90)

c     Initiate angle interval [0,180) at i9 = 1

      call timer('fabric',8,1)

      if (i9.eq.1) then
        dag = 5.d0
        nn  = int(180.d0/dag)
        do i = 0, nn
          agi(i) = dble(i)*dag
        end do ! i
      endif

c     Compute average coordinate number: effective

      if (i0.eq.1 .or. i9.eq.1) then
        nnn = 0
        do i = 1, n2b
          if (o(1,i).le.0.d0) then
            nnn = nnn + 1
          endif
        end do ! i
        acn = dble(nnn)/n0
        if (i9.eq.1) then
          write (2,*) 'Initial geometry: average coord. number:',acn
        else
          write (2,*) 'q02:',q02,',  average coord. number:',acn
        endif
        write (2,*) 'Interval','     median','   frequency'

c       Compute contact vector frequency

        do i = 1, nn
          ni = 0

c         Set a0 within [0,180)

          do j = 1, n2b
            a0 = dm(3,j)/dd
            if (a0.lt.0.0d0.and.a0.ge.-180.d0) a0 = a0 + 180.d0
            if (a0.lt.360.d0.and.a0.ge.180.d0) a0 = a0 - 180.d0
            if (a0.lt.540.d0.and.a0.ge.360.d0) a0 = a0 - 360.d0

            if (o(1,j).gt.0.d0) then
              if (a0.ge.agi(i-1).and.a0.lt.agi(i)) then
                ni = ni + 1
              endif
            endif
          end do ! j
          if (nnn.ne.0) then
            write (2,*) i,(agi(i-1)+agi(i))*0.5d0,dble(ni)/nnn
          else
            write (2,*) i,(agi(i-1)+agi(i))*0.5d0,dble(nnn)
          endif
        end do ! i

      end if

      call timer('fabric',8,2)
      end
