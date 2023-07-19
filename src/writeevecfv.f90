
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevecfv(ik,evecfv)
use modmain
implicit none
! arguments
!character(*), intent(in) :: fext
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
! local variables
integer recl,i,imat,istfv,ispnfv
character(256) fname,fik
! find the record length
!$OMP CRITICAL(u122)
write(fik,'("_IK",I6.6)')ik
fname=trim(scrpath)//'EVECFV_readable'//trim(fik)//'.OUT'
whchcount=whchcount+1
  open(122,file=trim(fname),form='FORMATTED')
  do imat=1,nmatmax
    do istfv=1,nstfv
      do ispnfv=1,nspnfv
        write(122,'(3G18.10,I6.4,2I4.2,2G18.10)') vkl(:,ik),imat,istfv,ispnfv,evecfv(imat,istfv,ispnfv)
      end do
    end do
  end do
  close(122)
  write(6,*)'inwriteevecfv',whchcount
  write(6,*)'iscl',iscl
!$OMP END CRITICAL(u122)
return
end subroutine

