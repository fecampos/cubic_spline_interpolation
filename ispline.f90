      subroutine ispline(u,newy,nu,x,y,b,c,d,n)

      implicit none

      integer, intent(in) :: n, nu

      real, intent(in) :: u(nu), x(n), y(n), b(n), c(n), d(n)

      real, intent(out) :: newy(nu)

      integer :: i, j, k

      real :: dx, uu(n)

      !$OMP PARALLEL DO
      do j = 1,nu
        uu = x-u(j)                
        where(uu.gt.0)
          uu = 0
        else where
          uu = 1
        end where
        k = sum(uu)
        if (k.eq.0) then
          newy(j) = y(1)
        elseif (k.eq.n) then
          newy(j) = y(n)
        else
          dx = u(j)-x(k)       
          newy(j) = y(k)+dx*(b(k)+dx*(c(k) + dx*d(k))) 
        end if
      end do
      !$OMP END PARALLEL DO

      end subroutine
