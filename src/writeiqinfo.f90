!whch: write iqinfo

subroutine writeiqinfo
!
use modmain
use modomp
implicit none
!local vars
integer(8) i,j,k,iq
character(256) filename1, filename2,filename3
!
!$OMP CRITICAL(u180)
write(*,*)
write(*,'("Info(writeiqinfo) : map from non-reduced grid to reduced index")' )
!
write(filename1,'("IQMAP.OUT")') 
open(50,file=trim(filename1),form='FORMATTED')
do i=0,ngridq(1)-1
  do j=0,ngridq(2)-1
    do k=0,ngridq(3)-1
      write(50,'(3I4.2,I8.6)') i,j,k,ivqiq(i,j,k)
    end do
  end do
end do
close(50)
!
write(*,*)
write(*,'("Info(writeiqinfo) : map from non-reduced grid to non-reduced index")' )
!
write(filename2,'("IQMAPNR.OUT")') 
open(51,file=trim(filename2),form='FORMATTED')
do i=0,ngridq(1)-1
  do j=0,ngridq(2)-1
    do k=0,ngridq(3)-1
      write(51,'(3I4.2,I8.6)') i,j,k,ivqiqnr(i,j,k)
    end do
  end do
end do
close(51)
!
! write(*,*)
! write(*,'("Info(writeivkinfo) : locations of k-points on integer grid")' )
! !
! write(filename3,'("IVK.OUT")') 
! open(52,file=trim(filename1),form='FORMATTED')
! do ik=1,nkptnr
!   write(52,'(I8.6,I4.4)') ivk(:,ik)
! end do
! close(52)
!$OMP END CRITICAL(u180)
return
end subroutine