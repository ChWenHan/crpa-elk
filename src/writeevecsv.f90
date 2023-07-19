
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevecsv(ik,evecsv)
use modmain
implicit none
! arguments
!character(*), intent(in) :: fext
integer, intent(in) :: ik
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer recl,i,ist,jst
character(256) fname,fik
! find the record length
!$OMP CRITICAL(u122)
whchcount=whchcount+1
write(fik,'("_IK",I6.6)')ik
fname=trim(scrpath)//'EVECSV_readable'//trim(fik)//'.OUT'
open(122,file=trim(fname),form='FORMATTED')
do ist=1,nstsv
  do jst=1,nstsv
      write(122,'(3G18.10,2I4.2,2G18.10)') vkl(:,ik),ist,jst,evecsv(ist,jst)
  end do
end do
close(122)
write(6,*)'inwriteevecsv',whchcount
write(6,*)'iscl',iscl
!$OMP END CRITICAL(u122)
return
end subroutine

