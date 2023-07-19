
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getwanprjk(fname,ikp,vpl,ld,nst,nproj,wanprj)
use modmain
use modwanprj
implicit none
! arguments
character(*), intent(in) :: fname
real(8),intent(in) :: vpl(3)
integer, intent(in) :: ikp,ld,nst,nproj
complex(8), intent(out) :: wanprj(ld,nst,nspinor,norb,nproj)
! local variables
integer recl,ld_,nst_,nspinor_,norb_,nproj_
integer isym
integer ik
integer i
real(8) vkl_(3),t1
complex(8) z1
! automatic arrays

! allocatable arrays
integer, allocatable :: map(:)
real(8), allocatable :: vgpl(:,:)
complex(8), allocatable :: cf_(:,:,:),x(:)
! find the equivalent reduced q-point and symmetry which rotates vql to vpl
if (ikp.gt.0) then
  ik=ikp
else
! find the equivalent k-point number and crystal symmetry element
  call findkpt(vpl,isym,ik)
end if
! find the record length
inquire(iolength=recl) vkl(:,ik),ld,nst,nspinor,norb,nproj,wanprj
!$OMP CRITICAL(u180)
do i=1,2
  open(800,file=trim(fname),form='UNFORMATTED',access='DIRECT',recl=recl,err=10)
  read(800,rec=ik,err=10) vkl_,ld_,nst_,nspinor_,norb_,nproj_,wanprj
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(getwanrpjk): unable to read from ",A)') trim(fname)
    write(*,*)
    stop
  end if
  close(800)
end do
!$OMP END CRITICAL(u180)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getwanrpjk): differing vectors for k-point ",I8)') ik
  write(*,'(" current : ",3G18.10)') vkl(:,ik)
  write(*,'(" file    : ",3G18.10)') vkl_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (ld.ne.ld_) then
  write(*,*)
  write(*,'("Error(getwanrpjk): differing ld for k-point ",I8)') ik
  write(*,'(" current : ",I8)') ld
  write(*,'(" file    : ",I8)') ld_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (nst.ne.nst_) then
  write(*,*)
  write(*,'("Error(getwanrpjk): differing nst for k-point ",I8)') ik
  write(*,'(" current : ",I8)') nst
  write(*,'(" file    : ",I8)') nst_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (nspinor.ne.nspinor_) then
  write(*,*)
  write(*,'("Error(getwanrpjk): differing nspinor for k-point ",I8)') ik
  write(*,'(" current : ",I8)') nspinor
  write(*,'(" file    : ",I8)') nspinor_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (norb.ne.norb_) then
  write(*,*)
  write(*,'("Error(getwanrpjk): differing norb for k-point ",I8)') ik
  write(*,'(" current : ",I8)') norb
  write(*,'(" file    : ",I8)') norb_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
if (nproj.ne.nproj_) then
  write(*,*)
  write(*,'("Error(getwanrpjk): differing nproj for k-point ",I8)') ik
  write(*,'(" current : ",I8)') nproj
  write(*,'(" file    : ",I8)') nproj_
  write(*,'(" in file ",A)') trim(fname)
  write(*,*)
  stop
end if
return
end subroutine

