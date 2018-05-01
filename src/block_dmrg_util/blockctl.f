************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2014, Naoki Nakatani                                   *
************************************************************************
      Subroutine BlockCtl(LW1,TUVX,IFINAL,IRST)
************************************************************************
*                                                                      *
*     Description                                                      *
*                                                                      *
*     This is an alternative of DavCtl for DMRG-CASSCF in DMRGCtl      *
*                                                                      *
************************************************************************
*                                                                      *
*     DMRG control section                                             *
*                                                                      *
*     calling arguments:                                               *
*     LW1     : Memory pointer to active Fock matrix                   *
*               array of real*8                                        *
*     TUVX    : array of real*8                                        *
*               two-electron integrals (tu!vx)                         *
*     IFINAL  : integer                                                *
*               termination flag                                       *
*     IRST    : integer                                                *
*               restart flag of DMRG                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     N. Nakatani, Hokkaido University, Japan, 2014                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension LW1(*), TUVX(*)

      Integer iChMolpro(8)
      Character*3 Label
      character(len=150) :: imp1
      character(len=3) :: block_nprocs, char_lroots

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='BLOCKCTL')
      Call qEnter(ROUTINE)

C Local print level (if any)
      IPRLEV=IPRLOC(3)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

#ifndef _NEW_BLOCK_

* Default setting for the first iteraction
      JRST=IRST
      ThDMRG=THRE*1.0d1
      ThNoise=1.0d-4

      If (JRST.NE.0) Then
        If (IFINAL.EQ.0) Then
          If (ABS(ROTMAX).GE.1.0d-1) Then
* No restart
            JRST=0
            ThDMRG=THRE*1.0d1
            ThNoise=1.0d-4
          Else If (ABS(ROTMAX).GE.0.5d-1) Then
* Full-restart with noise
            JRST=2
            ThDMRG=THRE*5.0d0
            ThNoise=1.0d-5
          Else If (ABS(ROTMAX).GE.0.1d-1) Then
* Full-restart with smaller noise
            JRST=2
            ThDMRG=THRE
            ThNoise=1.0d-6
          Else
* Restart without noise
            JRST=1
            ThDMRG=THRE
            ThNoise=0.0d0
          End If
        Else If (IFINAL.EQ.1) Then
* Full-restart for CI-only
          JRST=2
          ThDMRG=THRE
          ThNoise=1.0d-6
        Else
* Full-restart for the final wavefunction
          If (iOrbTyp.EQ.2) Then
* with OutOrb = Canonical, this will be problematic for LMO-DMRG calc.
            JRST=2
            ThDMRG=THRE
            ThNoise=1.0d-6
          Else
* by default, just restarting
            JRST=1
            ThDMRG=THRE
            ThNoise=0.0d0
          End If
        End If
      End If

* Load symmetry info from RunFile
      Call Get_iScalar('NSYM',nIrrep)
      Call Get_iArray('Symmetry operations',iOper,nIrrep)
      Call Get_iScalar('Rotational Symmetry Number',iSigma)

* Get character table to convert MOLPRO symmetry format
      Call MOLPRO_ChTab(nSym,Label,iChMolpro)

      NRDM_ORDER=2
      If (NACTEL.EQ.1) NRDM_ORDER=1

* Convert orbital symmetry into MOLPRO format
      Call Getmem('OrbSym','Allo','Inte',lOrbSym,NAC)
      iOrb=1
      Do iSym=1,nSym
        Do jOrb=1,NASH(iSym)
          iWork(lOrbSym+iOrb-1)=iChMolpro(iSym)
          iOrb=iOrb+1
        End Do
      End Do
      lSymMolpro=iChMolpro(lSym)

* Compute DMRG
      Call block_calldmrg(JRST,lRoots,NAC,NACTEL,ISPIN-1,
     &                    Label,lSymMolpro,iWork(lOrbSym),
     &                    0.0d0,LW1,TUVX,MxDMRG,NRDM_ORDER,
     &                    ThDMRG,ThNoise,ENER(1,ITER),HFOCC,
     &                    NRS2T)

      If (IFINAL.EQ.2 .AND. Do3RDM .AND. NACTEL.GT.2) Then
* Compute 3RDM for DMRG-cu4-CASPT2
        Call block_calldmrg(1,lRoots,NAC,NACTEL,ISPIN-1,
     &                      Label,lSymMolpro,iWork(lOrbSym),
     &                      0.0d0,LW1,TUVX,MxDMRG,3,
     &                      THRE,0.0d0,ENER(1,ITER),HFOCC,
     &                      NRS2T)
      End If

      Call Getmem('OrbSym','Free','Inte',lOrbSym,NAC)

#endif

#ifdef _NEW_BLOCK_

* Load symmetry info from RunFile
      iOper = 0
      Call Get_iScalar('NSYM',nIrrep)
      Call Get_iArray('Symmetry operations',iOper,nIrrep)
      Call Get_iScalar('Rotational Symmetry Number',iSigma)

* Get character table to convert MOLPRO symmetry format
      Call MOLPRO_ChTab(nSym,Label,iChMolpro)

* Convert orbital symmetry into MOLPRO format
      Call Getmem('OrbSym','Allo','Inte',lOrbSym,NAC)
      iOrb=1
      Do iSym=1,nSym
        Do jOrb=1,NASH(iSym)
          iWork(lOrbSym+iOrb-1)=iChMolpro(iSym)
          iOrb=iOrb+1
        End Do
      End Do
      lSymMolpro=iChMolpro(lSym)

      NRDM_ORDER=2
      If (NACTEL.EQ.1) NRDM_ORDER=1


**********************
*  WRITEOUT FCIDUMP  *
**********************

      LINSIZE = ( NAC * ( NAC + 1 ) ) / 2
      NUM_TEI = ( LINSIZE * ( LINSIZE + 1 ) ) / 2
      Call FCIDUMP_OUTPUT( NAC, NACTEL, ISPIN-1,
     &                     lSymMolpro, iWork(lOrbSym),
     &                     0.0d0, LW1, TUVX,
     &                     LINSIZE, NUM_TEI )


      Call Getmem('OrbSym','Free','Inte',lOrbSym,NAC)

*************************
*  WRITEOUT INPUT FILE  *
*************************
      LUTOTE = isFreeUnit(30)
      call molcas_open(LUTOTE,'dmrg.conf')

      write(LUTOTE,'(A4,A4)')  'sym ', Label
      write(LUTOTE,'(A16)') 'orbitals FCIDUMP'
      write(LUTOTE,*)

      write(LUTOTE,'(A5,I3)') 'nelec', NACTEL
      write(LUTOTE,'(A6,I3)') 'nroots', lroots
      write(LUTOTE,'(A4,I3)') 'spin', ISPIN-1
      write(LUTOTE,'(A5,I3)') 'irrep', lSymMolpro
      write(LUTOTE,*)
      write(LUTOTE,'(A13)') 'gaopt default'
      write(LUTOTE,'(A18)') 'warmup local_4site'
      write(LUTOTE,'(A7,I6)') 'maxiter', 50
      write(LUTOTE,'(A6)') 'twopdm'
      write(LUTOTE,'(A11)') 'memory 10 g'

! Always use restart option from 2nd iteration
      if (IRST>0 .or. blockrestart.EQV..TRUE.) then
        write(LUTOTE,'(A11)') 'fullrestart'
        write(LUTOTE,'(A8)') 'schedule'
        write(LUTOTE,'(2I6,2E12.5)') 0, MxDMRG, 1.0e-7, 0.0
        write(LUTOTE,'(A3)') 'end'
        write(LUTOTE,'(A17,I6)') 'twodot_to_onedot ', two2one
      else
        write(LUTOTE,'(A4,I6)') 'maxM', MxDMRG
        write(LUTOTE,'(A16)') 'schedule default'
      endif

      if (sum(hfocc) .EQ. NACTEL) then
        write(6,*)  'BLOCK> Using user-specified ROHF guess'
        write(LUTOTE,'(A7)',ADVANCE='NO') 'hf_occ '
        do ihfocc=1,NAC-1
          write(LUTOTE,'(I3)', ADVANCE='NO') HFOCC(ihfocc)
        enddo
        write(LUTOTE,'(I3)') HFOCC(NAC)
      else
        write(6,*)  'BLOCK> Using default guess'
        write(LUTOTE,'(A15)') 'hf_occ integral'
      endif
      write(LUTOTE,*)

      If (IFINAL.EQ.2 .AND. Do3RDM .AND. NACTEL.GT.2) Then
        write(LUTOTE,'(A8)') 'threepdm'
      endif


      close(LUTOTE)

      call get_environment_variable("MOLCAS_BLOCK",
     &                                block_nprocs, status=ierr)
        if (ierr.NE.0) then
          imp1="block.spin_adapted dmrg.conf >dmrg.out
     &         2>dmrg.err"
        else
          imp1 = "mpirun -np "//trim(adjustl(block_nprocs))//
     &   " block.spin_adapted dmrg.conf >dmrg.out 2>dmrg.err"
        endif

      call systemf(imp1,iErr)

      write(char_lroots,'(I3)') lroots
      imp1='grep "Sweep Energy" dmrg.out | tail -'//
     & trim(adjustl(char_lroots)) //' |
     &      cut -c 88- > block.energy'
      call systemf(imp1,iErr)

      LUTOTE = isFreeUnit(30)
      call molcas_open(LUTOTE,'block.energy')

      do chemroot=1,lroots
        read(LUTOTE,*) ENER(chemroot,ITER)
      enddo
      close(LUTOTE)

#endif

      Call qExit(ROUTINE)

      Return
      End
