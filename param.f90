      module param

      implicit none

      character(len=*),parameter :: file_in = "in.nc"
      character(len=*),parameter :: file_clim = "in_clim.nc"
     
      character(len=*),parameter :: t_NAME = "time"
      character(len=*),parameter :: x_NAME = "longitude"
      character(len=*),parameter :: y_NAME = "latitude"
      character(len=*),parameter :: z_NAME = "depth"

      character(len=*),parameter :: xi_NAME = "thetao"
      character(len=*),parameter :: mxi_NAME = "temp"

      integer, parameter :: nx = 361, ny = 6, nz = 43, nt = 1

      integer :: i, j, k

      real, parameter :: pi = 3.1415927

      real, parameter :: missing_val = -32767, sf_sla = 0.000732444226741791, af_sla = 21

      real, parameter :: zo = 5, zf = 3000, dz = 5

      real :: T(nt), X(nx), Y(ny), Z(nz), xi(nx,ny,nz,nt), mxi(nx,ny,nz,nt), yi(nx,ny,nz)

      integer :: msk3d(nx,ny,nz), msk2d(nx,ny), nn, newnz

      real :: b(nx,ny,nz), c(nx,ny,nz), d(nx,ny,nz)

      integer :: ncid, retval, tvarid, xvarid, yvarid, zvarid, xivarid, mxivarid

      real, allocatable :: newxi(:,:,:), newZ(:)

      end module
