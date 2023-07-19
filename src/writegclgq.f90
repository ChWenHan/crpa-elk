! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
subroutine writegclgq(iq,gclgq)
! !USES:
use modmain
implicit none
! local variables
integer, intent(in) :: iq
real(8), intent(in) :: gclgq(ngrf)
integer ig
character(256) filename
! write out the valence eigenvalues
!$OMP CRITICAL(u180)
write(filename,'("GCLGQ_Q",I6.6,".OUT")') iq
open(50,file=trim(filename),form='FORMATTED')
!write(50,'(I6," : nqpt . iq : ",I6.6)') nqpt, iq
do ig=1,ngrf
  write(50,'(I8.6,2X,G18.10)') ig,gclgq(ig)
end do
close(50)
!$OMP END CRITICAL(u180)
return
end subroutine