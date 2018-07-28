!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2017, Quan Phung                                       *
!***********************************************************************
! Subroutine to translate general 2RDM to pseudocanonical 2RDM
! Written by Quan Phung, Leuven, July 2017

subroutine chemps2_tran2pdm( NACT, LXMAT, CHEMROOT )

USE HDF5
USE ISO_C_BINDING

IMPLICIT NONE

CHARACTER(LEN=30) :: file_2rdm
CHARACTER(LEN=30) :: file_2rdm_tran

INTEGER, INTENT(IN) :: chemroot
INTEGER, INTENT(IN) :: nact
INTEGER             :: nact3, nact4

REAL(8), INTENT(IN) :: lxmat(nact*nact)

REAL*8, DIMENSION( 1 : nact*nact*nact*nact ), TARGET :: twordm

INTEGER( HID_T )   :: file_h5, group_h5, space_h5, dset_h5 ! Handles
INTEGER( HSIZE_T ) :: dimension
INTEGER(4)         :: hdferr
TYPE( C_PTR )      :: f_ptr

integer :: i, j, k, l, ind, ITER, iErr
logical :: irdm
integer :: ii, jj, kk, ll, p
real(8), allocatable :: outrdm(:)
real(8), allocatable :: tmprdm(:), tmprdm2(:)
character(len=10) :: rootindex

#ifdef _MOLCAS_MPP_
  EXTERNAL Is_Real_Par, KING
  Logical KING
  Logical Is_Real_Par
#include "mpif.h"
#endif

#ifdef _MOLCAS_MPP_
  if ( KING() ) then
#endif

nact3 = nact  * nact * nact
nact4 = nact3 * nact

! Check if files exist
write(rootindex,"(I2)") chemroot-1
file_2rdm="molcas_2rdm.h5.r"//trim(adjustl(rootindex))
file_2rdm=trim(adjustl(file_2rdm))
call f_inquire(file_2rdm, irdm)
if (.NOT. irdm) then
  write(6,'(1X,A15,I3,A16)') 'CHEMPS2> Root: ',CHEMROOT,' :: No 2RDM file'
  call abend()
endif

! Open and read HDF5 files
CALL h5open_f( hdferr )
CALL h5fopen_f( file_2rdm, H5F_ACC_RDONLY_F, file_h5, hdferr )
CALL h5gopen_f( file_h5, "2-RDM", group_h5, hdferr )
CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
CALL h5dget_space_f( dset_h5, space_h5, hdferr )
f_ptr = C_LOC( twordm( 1 ) )
CALL h5dread_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
CALL h5dclose_f( dset_h5 , hdferr )
CALL h5sclose_f( space_h5, hdferr )
CALL h5gclose_f( group_h5, hdferr )
CALL h5fclose_f( file_h5,  hdferr )


allocate(outrdm (nact4))
allocate(tmprdm (nact ))
allocate(tmprdm2(nact ))

ind = 1

! Combine RDM with LMAT
DO ITER=1,4

  p=0
  ind = 1

  do i=1,nact
    do j=1,nact
      do k=1,nact
        p = p + 1

        do l=1,nact
          tmprdm(l) = twordm(ind+l-1)
        enddo
        ! Transform twordm with lxmat
        call dgemv_('T', nact, nact, 1.0d0, lxmat, nact, tmprdm, 1, 0.0d0, tmprdm2, 1)

        ! Reorder elements
        do l=1,nact
          outrdm(p+(l-1)*nact3)=tmprdm2(l)
        enddo
        ind = ind + nact

      enddo
    enddo
  enddo

  call dcopy_(nact4,outrdm,1,twordm,1)

ENDDO

! Create a tran file and overwrite the old rdm

file_2rdm_tran="molcas_2rdm.h5.r"//trim(adjustl(rootindex))//".tran"
file_2rdm_tran=trim(adjustl(file_2rdm_tran))
call fcopy(file_2rdm,file_2rdm_tran,iErr)

CALL h5fopen_f( file_2rdm_tran, H5F_ACC_RDWR_F, file_h5, hdferr )
CALL h5gopen_f( file_h5, "2-RDM", group_h5, hdferr )
CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
CALL h5dget_space_f( dset_h5, space_h5, hdferr )
f_ptr = C_LOC( twordm( 1 ) )
CALL h5dwrite_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
CALL h5dclose_f( dset_h5 , hdferr )
CALL h5sclose_f( space_h5, hdferr )
CALL h5gclose_f( group_h5, hdferr )
CALL h5fclose_f( file_h5,  hdferr )

!call fcopy(file_2rdm_tran,file_2rdm,iErr)

deallocate(outrdm )
deallocate(tmprdm )
deallocate(tmprdm2)

#ifdef _MOLCAS_MPP_
  end if
#endif

end subroutine
