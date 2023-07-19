
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeoccsv
use modmain
implicit none
! arguments
! local variables
integer ik,ist1
character(256) fname
! find the record length
!$OMP CRITICAL(u130)
open(130,file='occsv_readable.out',form='FORMATTED',access='DIRECT', &
  action='WRITE')
do ik=1,nkpt
  do ist1=1,nstsv
    write(130,'(I8.4,3X,I8.4,G18.10)') ik,ist1,occsv(ist1,ik)
  end do
end do
close(130)
!$OMP END CRITICAL(u130)
return
end subroutine

