      subroutine gauss(ndim,a,b,x)
        implicit none
        integer ndim
        real*8 a(ndim,ndim),b(ndim),c(ndim)
        
!       integer, parameter:: dp=kind(0.d0)

!       integer,  intent(in)  :: ndim
!       real(dp), intent(in)  :: a(ndim,ndim), b(ndim)
!       real(dp), intent(out) :: x(ndim)
!       real(dp)              :: c(ndim)
!       integer               :: i,k,l,m
!       real(dp)              :: d,coef
! cf2py depend(ndim) a, b, x, c

!     Triangularization of matrix a
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
      end subroutine
