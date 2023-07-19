
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putwanprjk(fext,ik,ld,nst,nproj,wanprj)
use modmain
use modwanprj
implicit none
! arguments
character(*), intent(in) :: fext
integer, intent(in) :: ik,ld,nst,nproj
complex(8), intent(in) :: wanprj(ld,nst,nspinor,norb,nproj)
! local variables
integer recl,i
character(256) fname
! find the record length
inquire(iolength=recl) vkl(:,ik),ld,nst,nspinor,norb,nproj,wanprj
fname=trim('WANPRJ_ELK')//trim(fext)
!$OMP CRITICAL(u122)
do i=1,2
  open(805,file=trim(fname),form='UNFORMATTED',access='DIRECT',recl=recl,err=10)
  write(805,rec=ik,err=10) vkl(:,ik),ld,nst,nspinor,norb,nproj,wanprj
  close(805)
  exit
10 continue
  if (i.eq.2) then
    write(*,*)
    write(*,'("Error(putwanprjk): unable to write to ",A)') trim(fname)
    write(*,*)
    stop
  end if
  close(805)
end do
!$OMP END CRITICAL(u122)
return
end subroutine

