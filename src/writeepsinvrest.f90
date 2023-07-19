! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
subroutine writeepsinvrest(iq,fpre,epsir)
! !USES:
use modmain
use modomp
implicit none
! local variables
integer(4), intent(in) :: iq
character(256),intent(in) :: fpre
complex(8), intent(in) :: epsir(ngrf,ngrf,nwrf)
integer ig,jg,iw
character(256) filename
! write out the valence eigenvalues
!$OMP CRITICAL(u180)
!write(*,*)iq
write(*,'("Info(writeepsinvr) : write Eps/Pol for q :",I4.2)')iq
write(filename,'("_Q",I6.6,".OUT")') iq
filename=trim(fpre)//trim(filename)
open(50,file=trim(filename),form='FORMATTED')
!write(50,'("ig & ngrf:",I6.3," nqpt :",I6," iq : ",I6.4," wrf, epsir")') ngrf,nqpt, iq
do ig=1,ngrf
  do jg=1,ngrf
    do iw=1,nwrf
      write(50,'(3I8.6,2X,2G18.10)') ig,jg,iw,epsir(ig,jg,iw)
    end do
  end do
end do
close(50)
!$OMP END CRITICAL(u180)
return
end subroutine