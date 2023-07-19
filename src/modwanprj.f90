
! Copyright (C) 2015 Jon Lafuente and Manh Duc Le; 2017-18 Arsenii Gerasimov,
! Yaroslav Kvashnin and Lars Nordstrom. This file is distributed under the terms
! of the GNU General Public License. See the file COPYING for license details.

module modwanprj

!-------------------------------------------!
!   Wannier projector interface variables   !
!-------------------------------------------!
!band indices for wan projectors
integer, allocatable :: idxwanprj(:,:)
!number of state for the ik wan projector 
integer, allocatable :: projstwanprj(:)
! maximum number of nst wan prj
integer maxnstwanprj
! lm size for wanprj
integer ldwanprj
! 
integer nprojwanprj
! global array for subulm
complex(8), allocatable :: subulmwanprj(:,:,:,:)
!global array for sublm
integer, allocatable :: sublmwanprj(:,:,:)
!adnj edit - global dmft variables
!--------------------------------!
!     DMFT (TRIQS) variables     !
!--------------------------------!
!correlated energy from DFMT to be used for adjusted total energy
real(8) :: engydmft=0.d0
!dmft density matrix array
complex, allocatable :: dmatkdmft(:,:,:)

!--------------------------------!
!  Wannier Projector variables   !
!--------------------------------!
!number of desired projectors 
integer norb
!integer array containing the input species, atom and angular l number
!for each desired projector
integer, allocatable :: orb(:,:)
!reduced lm indices used in the wannier porjectors
integer, allocatable :: rorblm(:,:)
!the min and max correlated energy window limits
real(8) emincor, emaxcor
!logical to determine if projector is outputted in lm or irreducible lm basis
logical cubic
!logical to determine if band indices are going to be used for correlated window
logical wanind
!
! index to (l,m) pairs! note this var has been deleated in 8.6.4.
! for convenient, I add it here so consistent with 6.2.8. This var has been generated in 
! initwanprj.f90
integer, allocatable :: idxlm(:,:)
!integer whchtagn
!!!crpa
integer potmode
!
complex(8), allocatable :: symlmwanprj(:,:,:,:)


end module

