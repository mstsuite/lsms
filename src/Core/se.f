c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine se (e,l,nodes,det,v,nr,r,p,rel,nrp,eunit)
c     =================================================================
c
c----- routine integrates full relativistic se equation for a core state
c      at energy e.
c      the two coupled first order differential equations are
c
c      p                   p'=(e-v+2*c*c)*q/c-(l+1)*p/r
c      q                   q'=-(e-v)*p/c+(l-1)*p/r
c l in this case is really kappa which is -l-1 for sp dn and l for sp up.
c
c if rel is false then solve the Schrodinger equation
c
      implicit real*8 (a-h,o-z)
c     implicit automatic (a-h,o-z)
c
      parameter (big2=250.0d0)
c
      logical rel
c
      complex*16 ecomp
c
c
      integer nr
      integer nv
      parameter (nv=2)
      real*8 r(nrp),v(nr,nrp),p(nrp),y(nv),yp(nv),dy(nv)
c     real*8 r(*),v(nr,nrp),p(*),y(nv),yp(nv),dy(nv)
      real*8 eunit
c     =====================================================================
c     dimensioned to stop array bounds problems............................
c     =====================================================================
      integer isymm(1)
      real*8  vso(1,1)
c
c----- note watol->1/watol to avoid divisions in bsstep
      parameter (watol=1.0d+8)
c
      external dfv
c
c----- initialisation
c
      ecomp=e
c      ereal=e
      nodes=0
c
      escale=0.5d0*eunit
      allp1=abs(l*(l+1))
c
      do k=1,nrp
        p(k)=0.d0
      enddo
c
      if(r(1).eq.0.d0) then
	ist=2
      else
	ist=1
      endif
c
c     -------------------------------------------------------------
      call initwave(rel,r,v(1,1),nr,nv,ist,l,ecomp,y,escale,.true.,
     &              b,allp1_dfv)
c     -------------------------------------------------------------
c
      call dfv(r(ist),y,dy,nv,nr,r(ist),ecomp,v(1,ist),
     &     vso,1,isymm,eunit,b,allp1_dfv)
      p(ist)=y(1)*y(1)
      if(rel) p(ist)=p(ist)+y(2)*y(2)
c
      vk=v(2,ist)
      if(nr.eq.4) vk=vk+v(nr,ist)
      gkp=allp1/(r(ist)*r(ist))-eunit*(e-vk/r(ist))
c
c----- the outwards integration
c
c     write(6,*) ' outwards integration'
c     write (6,'(i5,f10.5,4e15.7)') 1,r(ist),y(1),y(2),dy(1),dy(2)
      do 40 k=ist+1,nrp
      yp(1)=y(1)
      yp(2)=y(2)
      kp=k-1
      rfrom=r(k-1)
      rto=r(k)
      x=r(k)
      htry=rto-rfrom
c
c----- advance integration using bs integrator
c
   45 continue
c     -------------------------------------------------------------
        call bsstep (y,dy,nv,nr,rfrom,htry,watol,yp,hdid,hnext,idid,
     &    r(kp),ecomp,v(1,kp),vso,1,isymm,eunit,dfv,b,allp1_dfv)
c     -------------------------------------------------------------
      if (abs(rfrom-rto).lt.1.0d-8*r(nrp)) then
      p(k)=y(1)*y(1)
      if(rel) p(k)=p(k)+y(2)*y(2)
      if(y(1)*yp(1).lt.0.d0.or.y(1).eq.0.d0) nodes=nodes+1
c     write (6,'(i5,f10.5,4e15.7)') k,x,y(1),y(2),dy(1),dy(2)
      else
      htry=min(hnext,rto-rfrom)
      goto 45
      endif
!
!       call bulirsch_stoer(y,dy,nv,rfrom,htry,
!    &                      r(kp),ecomp,v(1,kp),eunit,b,allp1,nv)
!       p(k)=y(1)*y(1)
!       if(rel) p(k)=p(k)+y(2)*y(2)

c
      vk=v(2,k)
      if(nr.eq.4) vk=vk+v(nr,k)
      gk=allp1/(x*x)-eunit*(e-vk/r(k))
      if (gk.gt.0.0d0.and.gkp.lt.0.0d0) then
c
c----- at the classical turning point find the value and slope of the
c      wavefunction
c
      kstop=k
      pout=y(1)
      dpout=y(2)
c     write (6,140) e,l,kstop
      go to 50
      endif
      gkp=gk
   40 continue
c
c----- if we have reached here then state not a core
c
        nodes=-1
	return
c
   50 continue
c
c----- the inwards integration
c
c
c----- initialise the inwards integration
c
c first find out where we should start
      y(1)=0.d0
      y(2)=0.d0
      x=r(nrp)
      call dfv(x,y,dy,nv,nr,x,ecomp,v(1,nrp),vso,1,isymm,eunit,
     &         b,allp1_dfv)
      frel=dy(1)
      phi=dy(2)
	if(rel) then
	  rto=sqrt((l*l-big2)/(l*l/(x*x)-phi))
	else
	  rto=sqrt((allp1-big2)/(allp1/(x*x)-phi))
	endif
	if(rto.gt.x) then
	  kn=nrp+1
	  x=rto
	  phi=big2/(x*x)
          rtphi=sqrt(phi)
          y(1)=exp(-rtphi*x)
          if(rel) then
          y(2)=(l/x-rtphi)*y(1)/frel
          else
          y(2)=-(1.0d0/x+rtphi)*y(1)/frel
          endif
	  goto 70
	endif
      do k=nrp,kstop,-1
c
c
c------ start integration with wkb approximation - the wavefunction is
c       in its exponential tail beyond the classical turning point,
c       start integration when wavefunction is bigger than exp(-big)
c
c dfv returns frel and phi in dy if both y(1) and y(2) are set to zero.
      y(1)=0.d0
      y(2)=0.d0
      x=r(k)
      call dfv(x,y,dy,nv,nr,r(k-1),ecomp,v(1,k-1),vso,1,isymm,eunit,
     &         b,allp1_dfv)
      frel=dy(1)
      phi=dy(2)
c
c----- check to see if we really are beyond the classical turning point
c
      if (phi.le.0.0) then
      write (6,160) e,l,k
      stop
      endif
c
      rtphi=sqrt(phi)
c************************************************************
c On CRAY the argument of exp() cannot be less than -708.396
c************************************************************
      if(rtphi*x.gt.708.396d0) then
        y(1)=0.d0
      else
        y(1)=exp(-rtphi*x)
      endif
c     y(1)=exp(-rtphi*x)
      if(rel) then
      y(2)=(l/x-rtphi)*y(1)/frel
      else
      y(2)=-(1.0d0/x+rtphi)*y(1)/frel
      endif
      p(k)=y(1)*y(1)
      if(rel) p(k)=p(k)+y(2)*y(2)
      if(phi*x*x.le.big2) then
	kn=k
	rto=x
	goto 70
      endif
      enddo
c
c----- if we have reached here then state not a core
c
        nodes=-1
	return
   70 continue
c     write (6,*) 'SE: inwards integration'
c     write (6,'(i5,f10.5,4e15.7)') kn,x,y(1),y(2),dy(1),dy(2)
c
c----- integrate towards the classical turning point
c
      do 120 k=kn-1,kstop,-1
      yp(1)=y(1)
      yp(2)=y(2)
      rfrom=rto
      rto=r(k)
      x=r(k)
      htry=rto-rfrom
c
c----- advance integration using bs integrator
c
  145 continue
        call bsstep (y,dy,nv,nr,rfrom,htry,watol,yp,hdid,hnext,idid,
     &    r(k),ecomp,v(1,k),vso,1,isymm,eunit,dfv,b,allp1_dfv)
      if (abs(rfrom-rto).lt.1.0d-8*r(nrp)) then
      p(k)=y(1)*y(1)
      if(rel) p(k)=p(k)+y(2)*y(2)
      if(y(1)*yp(1).lt.0.d0.or.y(1).eq.0.d0.and.k.gt.kstop)
     &  nodes=nodes+1
c     write (6,'(''SE:'',i5,f10.5,4e15.7)') k,x,y(1),y(2),dy(1),dy(2)
      else
      htry=-min(abs(hnext),abs(rto-rfrom))
      goto 145
      endif
  120 continue
c
c----- calculate slope from inwards integration
c
      pin=y(1)
      dpin=y(2)
c
c----- match up the wavefunction for kstop to nrp
c
      scale=pout/pin
c     write(6,*) 'SE: scale:',scale
      dpin=dpin*scale
      scale=scale*scale
c
c     do 130 k=kstop,min(kn,nrp)  # malcolm had this
      do 130 k=kstop,nrp
      p(k)=scale*p(k)
  130 continue
c
      det=dpin-dpout
c     write (6,170) e,l,dpin,dpout,det
c
      return
c
  140 format(' classical turning point reached in outwards integration'/
     1 ' core e = ',f10.4,', k = ',i2,', grid point = ',i5)
  160 format(' unable to start inwards integration'/
     1 ' core e = ',f10.4,', k = ',i2,', grid point = ',i5)
  170 format('SE: core e =',f10.4,' k = ',i2/' slope (in)  = ',e12.5/
     1 ' slope (out) = ',e12.5/' difference in slopes (det) = ',e12.5)
      end
c
c
c**********************************************************************
c
