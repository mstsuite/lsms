c     =============================================================
      subroutine setup_boundary_cluster(i_seed,num_seeds,
     >                          seedx,seedy,seedz,
     >                          vp,ipvp,num_neighbors,rad)
c     =============================================================
c
      implicit   none
c

      integer    max_planes
      parameter (max_planes = 100)
c
      integer    i_seed
      integer    num_seeds
      integer    num_neighbors
      integer    id_neighbor(max_planes)
      integer    tempid(max_planes)
      integer    j_seed
      integer    nm1
      integer    isigma

      integer    n1
      integer    n2
      integer    n3
      integer    i0
      integer    j0
      integer    k0
      integer    i
      integer    j
      integer    k
      integer    n
      integer    id
      integer    nb
      integer    nbt
      integer    checked(max_planes)
c
      real*8     bravais_1(3)
      real*8     bravais_2(3)
      real*8     bravais_3(3)
      real*8     seedx(num_seeds)
      real*8     seedy(num_seeds)
      real*8     seedz(num_seeds)
      real*8     shift_x
      real*8     shift_y
      real*8     shift_z
      integer ipvp
      real*8     vp(3,ipvp)
      real*8     temp(3,max_planes)
      real*8     tmpr2(max_planes)
      real*8     x0
      real*8     y0
      real*8     z0
      real*8     x
      real*8     y
      real*8     z
      real*8     r2
      real*8     xt
      real*8     yt
      real*8     zt
      real*8     r2t
      real*8     half
      real*8     tol
      real*8     rad(num_seeds)
      real*8     rad2my
      real*8     rlam
      real*8     one
c
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (tol=1.0d-06)
c
c     ***************************************************************
c     Calculate boundary planes of the Voronoi polyhedron............
c     ***************************************************************
c
      rad2my=rad(i_seed)**2
      x0=seedx(i_seed)
      y0=seedy(i_seed)
      z0=seedz(i_seed)

c
      j=0

         do 110 j_seed=1,num_seeds
	    x=half*(seedx(j_seed) - x0)
	    y=half*(seedy(j_seed) - y0)
	    z=half*(seedz(j_seed) - y0)
	    r2=x*x+y*y+z*z
	    nb=j_seed
	    if(r2.lt.tol*tol) goto 110
c           --------------------------------------------------------------------
c           code for radical plane construction
c           --------------------------------------------------------------------
            rlam=one+.25d0*(rad2my-rad(j_seed)**2)/r2
            x=rlam*x
            y=rlam*y
            z=rlam*z
            r2=rlam**2*r2
c           --------------------------------------------------------------------
	    id=1
	    do k=1,j
	       id=k
	       if(abs(x-temp(1,k)).lt.tol .and. 
     >            abs(y-temp(2,k)).lt.tol .and.
     >            abs(z-temp(3,k)).lt.tol) then
	          goto 110
	       endif
	       if(r2.lt.tmpr2(k)) then
	          goto 100
	       endif
	    enddo
	    id=j+1
 100        continue
	    do k=id,j
	       xt=temp(1,k)
	       yt=temp(2,k)
	       zt=temp(3,k)
	       r2t=tmpr2(k)
	       nbt=tempid(k)
               temp(1,k)=x
               temp(2,k)=y
               temp(3,k)=z
               tmpr2(k)=r2
	       tempid(k)=nb
               x=xt
               y=yt
               z=zt
               r2=r2t
	       nb=nbt
	    enddo
	    if(j.lt.max_planes) then
	       j=j+1
               temp(1,j)=x
               temp(2,j)=y
               temp(3,j)=z
               tmpr2(j)=r2
	       tempid(j)=nb
	    endif
 110     continue

      nm1=j
c
c     ================================================================
c     reduce nm1 to speed up the process of looking for boundary planes
c     Warning: it may cause problems in some special situtaions.......
c     ================================================================
      nm1=min(nm1,max_planes)
c
      num_neighbors=0
      k=nm1
      j=1
      id=1
      n=0
      do while(id.le.nm1)
c        ============================================================
c        check if [temp(x,j),temp(y,j),temp(z,j)] is an possible boundary 
c        plane. It is a preliminary check before chkbnd..............
c        ============================================================
	 isigma=0
         do i=1,n
            if(id.eq.checked(i)) then
               isigma=1
            endif
         enddo
         if(isigma.eq.0) then
c           ---------------------------------------------------------
	    call filter_edge(temp(1,j),temp(2,j),temp(3,j),tmpr2(j),j,
     >                       temp,tmpr2,k,isigma,i)
c           ---------------------------------------------------------
            if(isigma.eq.1 .and. i.gt.j) then
               n=n+1
	    if(n.gt.max_planes) then
	      write(6,'(''setup_boundary: n > max_planes'')')
	      stop 'setup_boundary'
	    endif
               checked(n)=i+nm1-k
            endif
         endif
         if(isigma.eq.1) then
            num_neighbors=num_neighbors+1
	    if(num_neighbors.gt.ipvp) then
	      write(6,'(''setup_boundary: num_neighbors > ipvp'')')
	      stop 'setup_boundary'
	    endif
            vp(1,num_neighbors)=temp(1,j)
            vp(2,num_neighbors)=temp(2,j)
            vp(3,num_neighbors)=temp(3,j)
	    id_neighbor(num_neighbors)=tempid(j)
	    j=j+1
	 else
	    k=k-1
	    do i=j,k
	       temp(1,i)=temp(1,i+1)
	       temp(2,i)=temp(2,i+1)
	       temp(3,i)=temp(3,i+1)
	    enddo
	    do i=j,k
	       tmpr2(i)=tmpr2(i+1)
	    enddo
	    do i=j,k
	       tempid(i)=tempid(i+1)
	    enddo
         endif
	 id=id+1
      enddo
c
      do j=1,num_neighbors
	 temp(1,j)=vp(1,j)
	 temp(2,j)=vp(2,j)
	 temp(3,j)=vp(3,j)
	 tmpr2(j)=vp(1,j)*vp(1,j)+vp(2,j)*vp(2,j)+vp(3,j)*vp(3,j)
      enddo
      do j=1,num_neighbors
	 tempid(j)=id_neighbor(j)
      enddo
      nm1=num_neighbors
c
      num_neighbors=0
      k=nm1
      j=1
      id=1
      do while(id.le.nm1)
c        ============================================================
c        check if [temp(x,j),temp(y,j),temp(z,j)] is an actual boundary 
c        plane
c        ------------------------------------------------------------
         call chkbnd(temp(1,j),temp(2,j),temp(3,j),tmpr2(j),j,
     >               temp,tmpr2,k,isigma)
c        ------------------------------------------------------------
         if(isigma.eq.1) then
            num_neighbors=num_neighbors+1
	    if(num_neighbors.gt.ipvp) then
	      write(6,'(''setup_boundary: num_neighbors > ipvp'')')
	      stop 'setup_boundary'
	    endif
            vp(1,num_neighbors)=temp(1,j)
            vp(2,num_neighbors)=temp(2,j)
            vp(3,num_neighbors)=temp(3,j)
            id_neighbor(num_neighbors)=tempid(j)
	    j=j+1
	 else
	    k=k-1
	    do i=j,k
	       temp(1,i)=temp(1,i+1)
	       temp(2,i)=temp(2,i+1)
	       temp(3,i)=temp(3,i+1)
	    enddo
	    do i=j,k
	       tmpr2(i)=tmpr2(i+1)
	    enddo
	    do i=j,k
	       tempid(i)=tempid(i+1)
	    enddo
         endif
	 id=id+1
      enddo
c
      if(num_neighbors.gt.max_planes) then
         write(6,'('' The number of neighbors exceeds the limit:'',
     >             '' num_neighbors,max_planes ='',2i5)')
     &         num_neighbors,max_planes
	 stop'setup_boundary'
      endif
c
c     ===============================================================
c     write(6,'(/,6x,''BOUNDARY PLANES:'')')
c     write(6,'(  6x,2i5,3f10.5)')
c    >(j,id_neighbor(j),vp(1,j),vp(2,j),vp(3,j),j=1,num_neighbors)
c     ===============================================================
c
      return
      end
