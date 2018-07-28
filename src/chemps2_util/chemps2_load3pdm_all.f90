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
! Copyright (C) 2018, Quan Phung                                       *
!***********************************************************************
! Subroutine to load all elements of 3-RDM

subroutine chemps2_load3pdm_all( NAC, G3T, chemroot, trans )

  USE HDF5
  USE ISO_C_BINDING
#ifdef _MOLCAS_MPP_
  USE MPI
#endif

  IMPLICIT NONE
  INTEGER, INTENT(IN)   :: NAC, chemroot
  LOGICAL, INTENT(IN)   :: trans
  REAL*8, INTENT(OUT)   :: G3T(NAC, NAC, NAC, NAC, NAC, NAC)

  CHARACTER(LEN=30) :: file_3rdm
  LOGICAL           :: irdm, jrdm
  CHARACTER(LEN=50) :: imp

  INTEGER( HID_T )   :: file_h5, group_h5, space_h5, dset_h5 ! Handles
  INTEGER(4)         :: hdferr
  TYPE( C_PTR )      :: f_ptr
  character(len=10) :: rootindex

#ifdef _MOLCAS_MPP_
  EXTERNAL Is_Real_Par, KING
  Logical KING
  Logical Is_Real_Par
#endif

  INTEGER :: ip1, ip2, ip3, iq1, iq2, iq3, idx, iG3, ierr

  REAL*8, DIMENSION( 1 : NAC * NAC * NAC * NAC * NAC * NAC ), TARGET :: buffer

  write(rootindex,"(I2)") chemroot-1

  if (trans) then
    file_3rdm="molcas_3rdm.h5.r"//trim(adjustl(rootindex))//".tran"
  else
    file_3rdm="molcas_3rdm.h5.r"//trim(adjustl(rootindex))
  endif

  file_3rdm=trim(adjustl(file_3rdm))
  call f_inquire(file_3rdm, irdm)
  if ((.NOT. irdm)) then
#ifdef _MOLCAS_MPP_
     if ( KING() ) then
       write(6,'(1X,A15,I3,A16)') 'CHEMPS2> Root: ',CHEMROOT,' :: No 3RDM file'
       call abend()
     endif
#endif
     imp="ln -sf ../" // file_3rdm // " ."
     imp=trim(adjustl(imp))
     call systemf(imp,iErr)
     write(6,'(1X,A46)') 'CHEMPS2> Automatically symbolic link 3RDM file'

  endif

  CALL h5open_f( hdferr )
  CALL h5fopen_f( file_3rdm, H5F_ACC_RDONLY_F, file_h5, hdferr )
  CALL h5gopen_f( file_h5, "3-RDM", group_h5, hdferr )
  CALL h5dopen_f( group_h5, "elements", dset_h5, hdferr )
  CALL h5dget_space_f( dset_h5, space_h5, hdferr )
  f_ptr = C_LOC( buffer( 1 ) )
  CALL h5dread_f( dset_h5, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
  CALL h5dclose_f( dset_h5 , hdferr )
  CALL h5sclose_f( space_h5, hdferr )
  CALL h5gclose_f( group_h5, hdferr )
  CALL h5fclose_f( file_h5,  hdferr )

  do ip1=0,NAC-1
    do ip2=0,NAC-1
      do ip3=0,NAC-1
        do iq1=0,NAC-1
          do iq2=0,NAC-1
            do iq3=0,NAC-1
              idx = ip1 + NAC * ( ip2 + NAC * ( ip3 + NAC * ( iq1 + NAC * ( iq2 + NAC * iq3 ) ) ) )
              G3T(ip1+1, iq1+1, ip2+1, iq2+1, ip3+1, iq3+1 ) = buffer( 1 + idx )
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine
