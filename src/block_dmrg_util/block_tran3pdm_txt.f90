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
! Subroutine to translate general 3RDM to corresponding pseudocanonical RDM
! Written by Quan Phung, Leuven, 2018

subroutine block_tran3pdm_txt( NACT, LXMAT, CHEMROOT )

IMPLICIT NONE

CHARACTER(LEN=50) :: file_3rdm

INTEGER, INTENT(IN) :: chemroot
INTEGER, INTENT(IN) :: nact
INTEGER             :: nact2, nact3, nact4, nact5, nact6

REAL(8), INTENT(IN) :: lxmat(nact*nact)

REAL*8, DIMENSION( 1 : nact*nact*nact*nact*nact*nact ) :: threerdm

integer :: i, j, k, l, m, n, iter, ind, LU
integer :: p
real(8), allocatable :: outrdm(:)
real(8), allocatable :: tmprdm(:), tmprdm2(:)
character(len=10) :: rootindex

nact2 = nact  * nact
nact3 = nact2 * nact
nact4 = nact3 * nact
nact5 = nact4 * nact
nact6 = nact5 * nact

write(6,*) 'BLOCK> Translate 3RDM into canonical basis'
call block_load3pdm_txt( NACT, threerdm, CHEMROOT, .FALSE. )

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


! Make a new 3RDM file
write(rootindex,"(I2)") chemroot-1
file_3rdm="./node0/spatial_threepdm."//trim(adjustl(rootindex))//"."//trim(adjustl(rootindex))//".txt.trans"
file_3rdm=trim(adjustl(file_3rdm))

LU=41
call molcas_open(LU,file_3rdm)

write(LU,*) nact
do i=1,nact
  do j=1,nact
    do k=1,nact
      do l=1,nact
        do m=1,nact
          do n=1,nact
            ind = (i-1)*nact5 + (n-1)*nact4 + (j-1)*nact3 + (m-1)*nact2 + (k-1)*nact + l
            write(LU,'(6I3,E20.12)') i-1, j-1, k-1, l-1, m-1, n-1, threerdm(ind)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
close(LU)

deallocate(outrdm )
deallocate(tmprdm )
deallocate(tmprdm2)

end subroutine
