!write kmesh
subroutine writekmesh
! !USES:
use modmain
!
implicit none
! local variables
integer i1,i2,i3,ik
integer indnr, indir
real(8) v1,v2,v3
open(50,file='NRKPTS.OUT',form='FORMATTED')
!write(50,'("vkl;  index_ir below")') 
! do i3=0,ngridk(3)-1
!   v3=dble(i3)/dble(ngridk(3))
!   do i2=0,ngridk(2)-1
!     v2=dble(i2)/dble(ngridk(2))
!     do i1=0,ngridk(1)-1
!       v1=dble(i1)/dble(ngridk(1))
!       indnr=ivkiknr(i1,i2,i3)
!       indir=ivkik(i1,i2,i3)
!       write(50,'(3G18.10,2I8)') vkl(:,indnr),indnr,indir
!     end do
!   end do
! end do
do ik=1,nkptnr
  i1=ivk(1,ik);i2=ivk(2,ik);i3=ivk(3,ik)
  indir=ivkik(i1,i2,i3)
  write(50,'(3G18.10,I8)') vkl(:,ik),indir
end do
close(50)
return
end subroutine