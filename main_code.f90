      program main_code

      use netcdf

      use param

      implicit none

      retval = nf90_open(file_in, NF90_NOWRITE, ncid)

      retval = nf90_inq_varid(ncid, t_NAME, tvarid)
      retval = nf90_inq_varid(ncid, x_NAME, xvarid)
      retval = nf90_inq_varid(ncid, y_NAME, yvarid)
      retval = nf90_inq_varid(ncid, z_NAME, zvarid)

      retval = nf90_get_var(ncid, tvarid, T)
      retval = nf90_get_var(ncid, xvarid, X)
      retval = nf90_get_var(ncid, yvarid, Y)
      retval = nf90_get_var(ncid, zvarid, Z)

      retval = nf90_inq_varid(ncid, xi_NAME, xivarid)
      retval = nf90_get_var(ncid, xivarid, xi)

      retval = nf90_close(ncid)
     
      msk3d = 1

      where(xi.ne.missing_val) 
        xi = xi*sf_sla+af_sla
      end where
      
      where(xi(:,:,:,1).eq.missing_val)
        msk3d = 0
      end where

      msk2d = sum(msk3d,3)    

      yi = xi(:,:,:,1)

      !$OMP PARALLEL DO
      do j = 1,ny
        do i =1,nx
          if (msk2d(i,j).ne.0) then
            yi(i,j,msk2d(i,j)+1:nz) = xi(i,j,msk2d(i,j),1)
          end if
        end do
      end do
      !$OMP END PARALLEL DO 

      newnz = (zf-zo)/dz+1

      allocate(newZ(newnz))

      allocate(newxi(nx,ny,newnz))

      !$OMP PARALLEL DO
      do i = 1,newnz
        newZ(i) = zo+(i-1)*dz
      end do
      !$OMP END PARALLEL DO     

      !$OMP PARALLEL DO
      do  j = 1,ny
        do i = 1,nx
          call spline(Z,yi(i,j,:),b(i,j,:),c(i,j,:),d(i,j,:),nz) 
        end do
      end do     
      !$OMP END PARALLEL DO
   
      !$OMP PARALLEL DO
      do  j = 1,ny
        do i = 1,nx
          call ispline(newZ,newxi(i,j,:),newnz,Z,yi(i,j,:),b(i,j,:),&
                 & c(i,j,:),d(i,j,:),nz)
        end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO
      do k = 1,newnz
        do j = 1,ny
          do i = 1,nx
            if (Z(msk2d(i,j)).le.newZ(k)) then
              newxi(i,j,k) = missing_val
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      call write_iso_vert(nx,ny,newnz,nt,X,Y,newZ,T,missing_val,newxi)

      deallocate(newxi)
      deallocate(newZ)

      end program
