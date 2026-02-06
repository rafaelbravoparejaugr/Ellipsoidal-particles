
c     rem ******************** totally rewritten
c     rem possible contacts of disk-disk and disk-boundary
c     rem ******************** called by main

      subroutine contact (pi,d0,i9,n2,n2p,n2b,m3,m3b,innb,innbd0)

c     inputs -- pi,n0,k0(),d(),d0,i9,n2,m0(),m(),o(),f0()
c     create -- a(),b(),m2,m3,a()
c     output -- b(),m2,m3,a()

      implicit      none

      integer       i9,n2,n2p,n2b,m3,m3b,innb,innbd0
      real*8        pi,twopi,halfpi,d0

      integer       i,i0,i1,j,j1, k,k0
      real*8        a1,a1p, a2,a2p

      include      'param.h'
      include      'za.h'
      include      'zbd.h'
      include      'zm0.h'
      include      'zd.h'
      include      'zdm.h'
      include      'ze.h'
      include      'zf0.h'
      include      'zk1.h'
      include      'zk0b.h'
      include      'zn.h'
      include      'znnb.h'
      include      'zo.h'
      include      'zrd.h'

c     md: maximum number of disks
c     mi: maximum number of invasions

c     n0      : number of disks
c     n2      : number of possible contacts between balls at start of
c               current step
c     n2p     : n2 + number of d-p contacts at start of current step
c     n2b     : n2b + number of d-b contacts at start of current step
c     i9      : current time step number
c     m3      : number of d-d contacts at end of previous step
c     m3b     : m3 + number of d-p + d-b contacts at end of previous step
c     d0      : criterion of possible contact distance of every two disks
c     d()     : coordinates of disk center, their radia/friction
c               angle/accumulated rotation/radia*
c     k0b()   : beginning registered vertex no. of boundary i in bd()
c     bd()    : coordinates/friction angle/length of boundary vertices/dip
c     a()     : (1) inclination of line; (2) contact information (dm(),m0(),
c               o(),f0()) at end of previous step
c     m0(,2)  : contact indices at end of previous current step
c     dm(,1-6): ball #,boundary #,dips,phi,coh of a contact
c               at end and start of previous step
c     o(,1-2) : penetration & parallel moving distance at end of previous step
c     f0()    : lock position at end of previous step
c     innb    : interval of # timesteps nearest neighboe search is done
c     innbd0  : multiple of search dist. for nnb search
c

c     Save contact information of previous step into a()
      call timer('contact',3,1)

      twopi  = 2.0d0*pi
      halfpi = 0.5d0*pi

      if(i9.gt.1) then
        m3 = 0
        do i = 1, n2
          if (m0(2,i).ne.0) then
            m3 = m3+1
            do j = 1, 2
              a(j,1,m3) = dm(j,i)
            end do ! j
            a(3,2,m3) = m0(2,i)
            a(1,3,m3) = o(1,i)
            a(2,3,m3) = o(2,i)
            a(3,3,m3) = f0(i)
          endif
        end do ! i

        m3b = m3

        do i = n2+1, n2b
          if (m0(2,i).ne.0) then
            m3b = m3b+1
            do j = 1, 2
              a(j,1,m3b) = dm(j,i)
            end do ! j
            a(3,2,m3b) = m0(2,i)
            a(1,3,m3b) = o(1,i)
            a(2,3,m3b) = o(2,i)
            a(3,3,m3b) = f0(i)
          endif
        end do ! i

      endif ! i9 > 1

c     Determine possible disk contacts, exclude those with same group

      n2 = 0

c     Determine nearest neighbor array nnb(,)

      if (i9.eq.1 .or. mod(i9,innb).eq.0) then

        do i = 1,nd-1
          nnb(1,i) = 0
          do j = i+1,nd
            if(k1(1,i).ne.k1(1,j)) then
              a2 = sqrt((d(1,i)-d(1,j))**2
     &   +(d(2,i)-d(2,j))**2)-max(d(3,i),d(4,i))
     &  -max(d(3,j),d(4,j))
		
              if (a2.le.(innbd0*d0)) then
                nnb(1,i) = nnb(1,i)+1
                if (nnb(1,i).eq.nnnb) then
                  write(*,*) 'nnb array not large enough'
                  pause
                endif
                nnb(nnb(1,i)+1,i) = j
              endif
            endif
          end do ! j
        end do ! i

        if (i9.eq.1) then
          write(19,*) 'nearest neighbor arrays'
          write(19,*) 'timestep: ',i9
          do i = 1,nd-1
            write(19,*) i,nnb(1,i)
          end do ! i
        endif
      endif

      do i = 1, nd-1
        do j = 2, nnb(1,i)+1
          k = nnb(j,i)

          rd(1,1,1) = d(1,i) 
          rd(1,2,1) = d(2,i) 
          rd(1,3,1) = d(3,i) 
          rd(1,4,1) = d(4,i) 		  
          rd(1,5,1) = d(5,i) 		  
		
          rd(2,1,1) = d(1,k) 
          rd(2,2,1) = d(2,k) 
          rd(2,3,1) = d(3,k) 
          rd(2,4,1) = d(4,k) 		  
          rd(2,5,1) = d(5,k) 
		  
          cn=1
		  
          write(*,*) 'Contact.f ELEL 1',d(1,i),
     &    d(2,i),d(3,i),d(4,i),d(5,i)
          write(*,*) 'Contact.f ELEL 2',d(1,k),
     &    d(2,k),d(3,k),d(4,k),d(5,k)

		 		  
          call SEARCH()
		  
          a2  = pen(1,1)*snc
		  
          write(*,*) 'DISTANCIA refinada',a2

          if (a2.le.d0) then

          n2 = n2+1
			
          write(*,*) 'Contact found',n2
	  
c	     pause				
			
            if (d(3,i).le.d(3,k)) then
              dm(1,n2) = i
              dm(2,n2) = k
            else
              dm(1,n2) = i
              dm(2,n2) = k
            endif
            i0 = dm(1,n2)
            k0 = dm(2,n2)

            e(1,1) = d(1,i0)-d(1,k0)
            e(1,2) = d(2,i0)-d(2,k0)

            if (e(1,1).eq.0.d0) then
              if(e(1,2).gt.0.d0) then
                e(1,4)  =  halfpi
              elseif(e(1,2).lt.0.d0) then
                e(1,4)  = -halfpi
              endif
            else
              e(1,4) = atan(e(1,2)/abs(e(1,1)))
              if (e(1,1).lt.0.d0) e(1,4) = pi - e(1,4)
            endif
            if (e(1,4).lt.0.d0) e(1,4) = twopi + e(1,4)
            e(1,5) = e(1,4) + pi
            if (e(1,5).gt.twopi) e(1,5) = e(1,5) - twopi
            dm(3,n2) = e(1,5)
            dm(4,n2) = e(1,4)
            dm(5,n2) = min(d(6,i),d(6,k))
            dm(6,n2) = min(d(7,i),d(7,k))
            call exceed (n2,mi,'mi','in contact.....')
          endif
        end do ! j
      end do ! i

c     Determine possible contacts between disks and rigid polygons

      n2p = n2

c     Determine possible contacts between disks and rigid boundaries

      n2b = n2p
	  
      n2c = n2b
	  
      do i = 1, nd
        do j = np+1, np+nb
          i1 = k0b(1,j)
          j1 = bd(1,i1)
          do k = i1+1, i1+j1-1
            if (k.gt.i1+1) a1p = a1

c           a2: distance between disk and boundary
c           a1: disk center projection on boundary

          rd(1,1,1) = d(1,i) 
          rd(1,2,1) = d(2,i) 
          rd(1,3,1) = d(3,i) 
          rd(1,4,1) = d(4,i) 		  
          rd(1,5,1) = d(5,i) 		  
		  

          rd(2,1,1) = bd(1,k)
          rd(2,2,1) = bd(2,k)
          rd(2,3,1) = bd(1,k+1) 
          rd(2,4,1) = bd(2,k+1)		  
          rd(2,5,1) = 0.0d0
		  
          cn = 0

          write(*,*) 'Contact.f',rd(1,1,1),rd(1,2,1)
		  
          call SEARCH()

          a2 = pen(1,1)*snc*0.5d0         

          a1 = ((d(1,i)-bd(1,k))*(bd(1,k+1)-bd(1,k))+(d(2,i)-bd(2,k))*
     &         (bd(2,k+1)-bd(2,k)))/(bd(3,k)*bd(3,k))
	 
			
           write(*,*) 'DISTrefinadaFRONT',a2,'a1',a1
		 
c          pause
	 
           if (a2.gt.d0.or.a2.lt.-0.3d0*d0.or.a1.lt.0.d0.or.a1.gt.1.d0)
     &     goto 430

           write(*,*) 'd0',d0	  

c          Disk corner to boundary side contacts

           n2b       = n2b + 1
			
           write(*,*) 'Contact found',n2b
	  
c	   pause		
	  		
            dm(1,n2b) = i
            dm(2,n2b) = k
            dm(3,n2b) = bd(4,k)-halfpi
            dm(4,n2b) = bd(4,k)-halfpi
            dm(5,n2b) = min(d(6,i),bd(5,k))
            dm(6,n2b) = min(d(7,i),bd(6,k))
            goto 439

430         if (k.eq.i1+1) goto 420

            if (a1p.le.1.d0.or.a1.ge.0.d0) goto 420

            a2 = sqrt((d(1,i)-bd(1,k))**2
     &               +(d(2,i)-bd(2,k))**2)-max(d(3,i),d(4,i))
            if (a2.gt.d0) goto 420

c           Boundary vertex to disk perimeter contacts

            n2b       =  n2b + 1		
            dm(1,n2b) =  i
            dm(2,n2b) = -k
            e(1,1) = bd(1,k)-d(1,i)
            e(1,2) = bd(2,k)-d(2,i)
            if (e(1,1).eq.0.d0)then
              if(e(1,2).gt.0.d0) then
                e(1,4)  =  halfpi
              elseif(e(1,2).lt.0.d0) then
                e(1,4)  = -halfpi
              endif
            else
              e(1,4) = atan(e(1,2)/abs(e(1,1)))
              if (e(1,1).lt.0.d0) e(1,4) = pi-e(1,4)
            endif
            if (e(1,4).lt.0.d0) e(1,4) = twopi+e(1,4)
            dm(3,n2b) = e(1,4)
            dm(4,n2b) = 0.d0
            dm(5,n2b) = min(d(6,i),bd(5,k))
            dm(6,n2b) = min(d(7,i),bd(6,k))

439         call exceed (n2b,mi,'mi','in contact.....')

420         continue
          end do ! k
        end do ! j
      end do ! i
    

c     Zero m0(2,); o(1-2,); f0()



      do i = 1, n2b
        m0(2,i) = 0
        o(1,i)  = 0.0d0
        o(2,i)  = 0.0d0
        f0(i)   = 0.0d0
      end do ! i

      call timer('contact',3,2)
      end
