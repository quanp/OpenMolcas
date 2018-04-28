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
! Subroutine to translate general 2RDM to pseudocanonical 2RDM
! Written by Quan Phung, Leuven, 2018

subroutine block_tran2pdm_txt( NACT, LXMAT, CHEMROOT )

IMPLICIT NONE

CHARACTER(LEN=50) :: file_2rdm

INTEGER, INTENT(IN) :: chemroot
INTEGER, INTENT(IN) :: nact
INTEGER             :: nact2, nact3, nact4

REAL(8), INTENT(IN) :: lxmat(nact*nact)

REAL*8, DIMENSION( 1 : nact*nact*nact*nact ) :: twordm

integer :: i, j, k, l, iter, ind, LU
integer :: p
real(8), allocatable :: outrdm(:)
real(8), allocatable :: tmprdm(:), tmprdm2(:)
character(len=10) :: rootindex

nact2 = nact  * nact
nact3 = nact2 * nact
nact4 = nact3 * nact



write(6,*) 'BLOCK> Translate 2RDM into canonical basis'
call block_load2pdm_txt( NACT, twordm, CHEMROOT, .FALSE. )

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


! Make a new 2RDM file
write(rootindex,"(I2)") chemroot-1
file_2rdm="./node0/spatial_twopdm."//trim(adjustl(rootindex))//"."//trim(adjustl(rootindex))//".txt.trans"
file_2rdm=trim(adjustl(file_2rdm))

LU=41
call molcas_open(LU,file_2rdm)

write(LU,*) nact
do i=1,nact
  do j=1,nact
    do k=1,nact
      do l=1,nact
        write(LU,'(4I3,E20.12)') i-1, j-1, k-1, l-1, twordm( (i-1)*nact**3 + (l-1)*nact**2 + (j-1)*nact + k )/2.0d0
      enddo
    enddo
  enddo
enddo
close(LU)

deallocate(outrdm )
deallocate(tmprdm )
deallocate(tmprdm2)

end subroutine
