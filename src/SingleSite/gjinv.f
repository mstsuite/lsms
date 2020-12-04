      subroutine gjinv(a,n,nmax,detl)
c=====================
c
      implicit real*8 (a-h,o-z)
      complex*16 a(nmax,n),q,detl,cpi
      integer ind(10000)
      data cpi/(0.d0,3.1415926535897932d0)/
c
c gauss-jordan inversion with partial pivoting
c
      detl=(0.d0,0.d0)
      do 1 i=1,n
      k=i
      amaxi=cdabs(a(i,i))
c
      do 2 j=i+1,n
      atest=cdabs(a(i,j))
      if(atest.le.amaxi) goto 2
      k=j
      amaxi=atest
    2 continue
c
      ind(i)=k
      if(k.eq.i) goto 4
      detl=detl-cpi
      do j=1,n
        q=a(j,i)
        a(j,i)=a(j,k)
        a(j,k)=qa
      end do
c
    4 q=a(i,i)
      detl=detl+cdlog(q)
c
      q=(1.d0,0.d0)/q
      a(i,i)=(1.d0,0.d0)
      do j=1,n
        a(j,i)=a(j,i)*q
      end do
c
      do 6 j=1,n
      if(j.eq.i) goto 6
      q=a(i,j)
      a(i,j)=(0.d0,0.d0)
      do k=1,n
        a(k,j)=a(k,j)-a(k,i)*q
      end do
    6 continue
    1 continue
c
      do 8 i=n,1,-1
      j=ind(i)
      if(i.eq.j) goto 8
      do k=1,n
        q=a(i,k)
        a(i,k)=a(j,k)
        a(j,k)=q
      end do
    8 continue
c
      return
      end
