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
! Subroutine to translate general 3RDM and F.4RDM to corresponding pseudocanonical RDM
! Written by Quan Phung, Leuven, July 2017

subroutine chemps2_tran3pdm( NACT, LXMAT, CHEMROOT, TRAN3RDM )

USE HDF5
USE ISO_C_BINDING

IMPLICIT NONE

CHARACTER(LEN=30) :: file_3rdm
CHARACTER(LEN=30) :: file_3rdm_tran

CHARACTER(LEN=30) :: file_f4rdm
CHARACTER(LEN=30) :: file_f4rdm_tran

INTEGER, INTENT(IN) :: chemroot
INTEGER, INTENT(IN) :: nact
LOGICAL, INTENT(IN) :: tran3rdm

INTEGER             :: nact2, nact3, nact4, nact5, nact6

REAL(8), INTENT(IN) :: lxmat(nact*nact)

REAL*8, DIMENSION( 1 : nact*nact*nact*nact*nact*nact ), TARGET :: threerdm

INTEGER( HID_T )   :: file_h5, group_h5, space_h5, dset_h5 ! Handles
INTEGER( HSIZE_T ) :: dimension
INTEGER(4)         :: hdferr
TYPE( C_PTR )      :: f_ptr

integer :: i, j, k, l, m, n, iter, ind, iErr, jErr
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

nact2 = nact  * nact
nact3 = nact2 * nact
nact4 = nact3 * nact
nact5 = nact4 * nact
nact6 = nact5 * nact

! Check if files exist
write(rootindex,"(I2)") chemroot-1
file_3rdm  = "molcas_3rdm.h5.r"//trim(adjustl(rootindex))
file_f4rdm ="molcas_f4rdm.h5.r"//trim(adjustl(rootindex))

file_3rdm  = trim(adjustl( file_3rdm))
file_f4rdm = trim(adjustl(file_f4rdm))

call f_inquire( file_3rdm, iErr)
call f_inquire(file_f4rdm, jErr)

if ((.NOT. iErr) .OR. (.NOT. jErr)) then
  write(6,'(1X,A15,I3,A16)') 'CHEMPS2> Root: ',CHEMROOT,' :: No 3RDM or F.4RDM file'
  call abend()
endif

! Open and read HDF5 files
CALL h5open_f( hdferr )

if (tran3rdm) then
  CALL h5fopen_f(  file_3rdm, H5F_ACC_RDONLY_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5,   "3-RDM", group_h5, hdferr )
else
  CALL h5fopen_f( file_f4rdm, H5F_ACC_RDONLY_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5, "F.4-RDM", group_h5, hdferr )
endif

CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
CALL h5dget_space_f( dset_h5, space_h5, hdferr )
f_ptr = C_LOC( threerdm( 1 ) )
CALL h5dread_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
CALL h5dclose_f( dset_h5 , hdferr )
CALL h5sclose_f( space_h5, hdferr )
CALL h5gclose_f( group_h5, hdferr )
CALL h5fclose_f( file_h5,  hdferr )

allocate(outrdm (nact6))
allocate(tmprdm (nact ))
allocate(tmprdm2(nact ))

ind = 1

! Combine RDM with LMAT
DO ITER=1,6

  p=0
  ind = 1

  do i=1,nact
    do j=1,nact
      do k=1,nact
        do l=1,nact
          do m=1,nact
            p = p + 1

            do n=1,nact
              tmprdm(n) = threerdm(ind+n-1)
            enddo
            call dgemv_('T', nact, nact, 1.0d0, lxmat, nact, tmprdm, 1, 0.0d0, tmprdm2, 1)
            do n=1,nact
              outrdm(p+(n-1)*nact5)=tmprdm2(n)
            enddo
            ind = ind + nact
          enddo
        enddo
      enddo
    enddo
  enddo

  call dcopy_(nact6,outrdm,1,threerdm,1)

ENDDO


! Creat a tran file and overwrite the old rdm

file_3rdm_tran  =  "molcas_3rdm.h5.r"//trim(adjustl(rootindex))//".tran"
file_f4rdm_tran = "molcas_f4rdm.h5.r"//trim(adjustl(rootindex))//".tran"

file_3rdm_tran  = trim(adjustl( file_3rdm_tran))
file_f4rdm_tran = trim(adjustl(file_f4rdm_tran))

if (tran3rdm) then
  call fcopy( file_3rdm, file_3rdm_tran,iErr)
else
  call fcopy(file_f4rdm,file_f4rdm_tran,iErr)
endif

if (tran3rdm) then
  CALL h5fopen_f(  file_3rdm_tran, H5F_ACC_RDWR_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5,   "3-RDM", group_h5, hdferr )
else
  CALL h5fopen_f( file_f4rdm_tran, H5F_ACC_RDWR_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5, "F.4-RDM", group_h5, hdferr )
endif

CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
CALL h5dget_space_f( dset_h5, space_h5, hdferr )
f_ptr = C_LOC( threerdm( 1 ) )
CALL h5dwrite_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
CALL h5dclose_f( dset_h5 , hdferr )
CALL h5sclose_f( space_h5, hdferr )
CALL h5gclose_f( group_h5, hdferr )
CALL h5fclose_f( file_h5,  hdferr )

if (tran3rdm) then
  call fcopy( file_3rdm_tran, file_3rdm,iErr)
else
  call fcopy(file_f4rdm_tran,file_f4rdm,iErr)
endif

deallocate(outrdm )
deallocate(tmprdm )
deallocate(tmprdm2)

#ifdef _MOLCAS_MPP_
  end if
#endif

end subroutine
