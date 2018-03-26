      subroutine srcore(rel,ecore,rho,nc,lcore,v,nr,r,nsrch
     *     ,thresh,nrp,eunit)
c
c----- rcore performs an energy search for a core state with estimated
c      energy ecore and angular momentum lcore.  thresh determines the
c      accuracy to which the energy is calculated, which is returned in
c      ecore.
c
c----- calculates charge density from the core state's wavefunctn
c      including contribution of small component of the scalar
c      relativistic wavefunction.
c
c
      implicit real*8 (a-h,o-z)
c
c---- calculation mode:
c
      logical rel               ! scalar relativistic if true
c
c---- dimensions
c
      integer nrp
c
c---- output:
c
      real*8 ecore              ! core energy
      real*8 rho(nrp)           ! core charge density
c
c---- input:
c
      integer lcore             ! core angular momentum number
      integer nr
      real*8 v(nr,nrp)             ! potential
      real*8 r(nrp)             ! radial grid
      integer nsrch             ! maximum number of searches
      real*8 thresh             ! tolerance for core energy search
      real*8 eunit
c
c---- work space:
c
      integer isrch             ! actual number of searches performed
c
c rho is used to store the large wavefunction in se
c
c n0 is the number of nodes in the radial wavefuction
      if(lcore.lt.0) then
      n0=nc+lcore
      else
      n0=nc-lcore-1
      endif
c
c----- initialize core search parameters
c
      isrch=0
c We don't yet have an element with z=107
      elow=-107.d0*107.d0/(eunit*nc*nc)
      ehigh=0.d0
    5 isrch=isrch+1
      if (isrch.gt.nsrch) go to 60
      call se (ecore,lcore,nodes,d0,v,nr,r,rho,rel,nrp,eunit)
      if(nodes.gt.n0) then
      ehigh=min(ecore,ehigh)
      ecore=0.5d0*(ehigh+elow)
      goto 5
      else if(nodes.lt.n0) then
      elow=max(ecore,elow)
      ecore=0.5d0*(ehigh+elow)
      goto 5
      endif
      if (d0.eq.0.0d0) return
      emn=ecore
      epl=ecore
      em0=ecore
      ep0=ecore
      dp0=d0
      dm0=d0
c
c Assume ecore was reasonably close
	dp=0.1d0*(ehigh-epl)
	dm=0.1d0*(emn-elow)
      isrch=0
c
c----- locate sign change in the determinant
c
   10 isrch=isrch+1
      if (isrch.gt.nsrch) go to 60
      epl=epl+dp
      call se (epl,lcore,nodes,dpl,v,nr,r,rho,rel,nrp,eunit)
      if(nodes.eq.n0) then
      ecore=epl
      if (dpl.eq.0.0d0) return
      if (dp0*dpl.lt.0.0d0) then
      emn=ep0
      dmn=dp0
      goto 20
      endif
      ep0=epl
      dp0=dpl
	dp=0.5d0*(ehigh-epl)
      else
	ehigh=epl
      if(dp.gt.0.d0) then
	dp=-0.5d0*dp
      else
	dp=0.5d0*dp
      endif
      endif
      emn=emn-dm
      call se (emn,lcore,nodes,dmn,v,nr,r,rho,rel,nrp,eunit)
      if(nodes.eq.n0) then
      ecore=emn
      if (dmn.eq.0.0d0) return
      if (dm0*dmn.lt.0.0d0) then
      epl=em0
      dpl=dm0
      goto 20
      endif
      em0=emn
      dm0=dmn
      else
	elow=emn
      if(dm.gt.0.d0) then
	dm=-0.5d0*dm
      else
	dm=0.5d0*dm
      endif
      endif
c     write(6,'(''nc,epl,emn='',2i3,1p2e15.6)') nc,isrch,epl,emn
      go to 10
c
c----- determine ecore to within thresh by successively halving the
c      interval
c
   20 isrch=isrch+1
      e0=ecore
      ecore=(epl*dmn-emn*dpl)/(dmn-dpl)
      call se (ecore,lcore,nodes,d0,v,nr,r,rho,rel,nrp,eunit)
      if (d0.eq.0.0d0) return
      if (d0*dpl.lt.0.0d0) then
      emn=ecore
      dmn=d0
      elseif (d0*dmn.lt.0.0d0) then
      epl=ecore
      dpl=d0
      endif
c     write(6,'(''SRCORE: nc,ecore,d0='',2i3,1p2e15.6)')
c    &     nc,isrch,ecore,d0
      if (abs(e0-ecore).ge.thresh) go to 20
      return
c
c----- if a sign change in the determinant is not found within nsrch
c      tries, the calculation is stopped
c
   60 write (6,70) ecore,epl,emn
      stop
c
   70 format('No zero found for core state near energy ',3f11.5)
c
      end
