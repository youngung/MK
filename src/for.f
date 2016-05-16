c----------------------------------------------------------------------c
c     Von Mises
      subroutine vm(s,phi,dphi,d2phi)
      implicit none
      integer i
      real*8 s(6),h,phi,dff,dphi(6),d2phi(6,6),d2h(6,6),d2ff
Cf2py intent(in,out) s
Cf2py intent(out) phi,dphi,d2phi
      h        = 5d-1*((s(1)-s(2))**2+s(1)**2+s(2)**2+6*s(6)**2)
      phi      = h**5d-1        ! yield surface
      s(:)     = s(:)/phi       ! stress on the yield locus
!     analytic solutions of first and second derivatives of the yield surface (plastic potential)
!     dff = 1./(2*(h**0.5))
      dff      = 1d0/(2*phi)
      dphi(1)  = dff*(2*s(1)-s(2))
      dphi(2)  = dff*(2*s(2)-s(1))
      dphi(6)  = dff*6*s(6)
      d2h(1,1) = 2d0
      d2h(1,2) =-1d0
      d2h(2,2) = 2d0
      d2h(6,6) = 6d0
      d2ff     = -(phi**(-3d0))/4
      d2phi(1,1) = d2ff*dphi(1)*dphi(1) + dff*d2h(1,1)
      d2phi(1,2) = d2ff*dphi(1)*dphi(2) + dff*d2h(1,2)
      d2phi(1,6) = 0d0
      d2phi(2,1) = d2phi(1,2)
      d2phi(2,2) = d2ff*dphi(2)*dphi(2) + dff*d2h(2,2)
      d2phi(2,6) = 0d0
      d2phi(6,1) = 0d0
      d2phi(6,2) = 0d0
      d2phi(6,6) = d2ff*dphi(6)*dphi(6) + dff*d2h(6,6)
      dphi(6) = dphi(6)/2
      return                    !! returns phi, dphi, d2phi
      end subroutine vm
c----------------------------------------------------------------------c
      subroutine hqe(s,r0,r90,phi,dphi,d2phi)
      implicit none
      real*8 s(3), r0, r90
      real*8 h,phi,dphi(3),d2phi, a, b, dff
cf2py intent(in,out) s
cf2py intent(in)     r0, r90
cf2py intent(out)    phi, dphi, d2phi
      a    =  ((r0 * (1. + r90)) /   (r90*(1.+r0)))
      b    =      (2.*r0)        /       (1.+r0)

      h    = s(1)**2 + a * (s(2)**2)
     $     - b *s(1)*s(2)
      phi  = h**0.5             ! yield surface
      s(:) = s(:) / phi         ! stress on the yield locus

!     dff     = 1.  / ( 2  * (h**0.5))
      dff     = 1.  / ( 2. * phi)
      dphi(:) = 0.
      dphi(1) = dff * ( 2.     *s(1) -  b  * s(2) )
      dphi(2) = dff * ( 2. * a *s(2) -  b  * s(1) )
      return
      end subroutine hqe
c----------------------------------------------------------------------c
      subroutine swift(e,ks,n,e0,m,sig,dsig,dm,qq)
      implicit none
      real*8 e,sig,dsig,dm
      real*8 ks,e0,n,qq,m
Cf2py intent(out) sig,dsig,n,m,dm,qq
Cf2py intent(in) e,ks,n,e0,m,qq

! swift
      sig  = ks*(e+e0)**n
      dsig = n*ks*(e+e0)**(n-1)
      dm = 0.0
c     m: strain rate sensitivty
c     qq: ??

c$$$      write(*,*)'sig:',sig
c$$$      write(*,*)'dsig:',dsig
c$$$      write(*,*)'n:',n
c$$$      write(*,*)'dm:',dm
c$$$      write(*,*)'qq:',qq
      return
      end subroutine swift
c----------------------------------------------------------------------c
c     a: jacobian
c     x: res
c     b: f
c     x = a:b
      subroutine gauss(ndim,b,a,x)
      implicit none
      integer ndim,i,k,l,m
      real*8 a(ndim,ndim),b(ndim),x(ndim)
      real*8 c(ndim),d,coef
Cf2py intent(in) ndim
cf2py intent(out) x
cf2py intent(in,out) b,a
cf2py depend(ndim) a,b,x,c

! triangularization of matrix a
      do k = 1 , ndim - 1
         i = k
         do while (a(i,k).eq.0.0.and.i.le.ndim)
            i = i + 1
         end do
         if(i.gt.ndim) then
            write(6,*)'###################################'
            write(6,*)'i m p o r t a n t    m e s s a g e '
            write(6,*)'###################################'
            write(6,*)
            write(6,*)'     determinant = 0 in gauss'
            write(6,*)'      does not converge with'
            write(6,*)'        current input data'
            write(6,*)'     try small chnage in input'
            write(6,*)
            write(6,*)'###################################'
            stop -1
         else
            if(i.eq.k) then
               continue
            else
               do l = k , ndim
                  c(l) = a(i,l)
                  a(i,l) = a(k,l)
                  a(k,l) = c(l)
               end do
               d = b(i)
               b(i) = b(k)
               b(k) = d
            end if
         end if
         do m = k+1 , ndim
            coef = a(m,k)/a(k,k)
            do l = k , ndim
               a(m,l) = a(m,l) - a(k,l)*coef
            end do
            b(m) = b(m) - b(k)*coef
         end do
      end do

! solution of linear system a*x = b
      x(ndim) = b(ndim)/a(ndim,ndim)
      do k = ndim-1 , 1 , -1
         x(k) = b(k)/a(k,k)
         do l = k+1 , ndim
            x(k) = x(k) - a(k,l)*x(l)/a(k,k)
         end do
      end do

      return
      end subroutine gauss

c----------------------------------------------------------------------c
      subroutine norme(ndim,vec,norm)
      implicit none
      integer ndim,i
      real*8 vec(ndim),norm
cf2py intent(in) ndim, vec
cf2py intent(out) norm
cf2py depend(ndim) vec

      norm = 0.0
      do i = 1 , ndim
         norm = norm + vec(i)*vec(i)
      end do
      norm = norm**0.5
      return
      end subroutine norme
