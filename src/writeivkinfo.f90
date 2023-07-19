!whch: write ivkik

subroutine writeivkinfo
!
use modmain
use modomp
implicit none
!local vars
integer(8) i,j,k,ik
character(256) filename1, filename2,filename3
!
!$OMP CRITICAL(u180)
write(*,*)
write(*,'("Info(writeivkinfo) : map from (i1,i2,i3) to (non)reduced k-point index")' )
!
write(filename1,'("IVKIK.OUT")') 
open(50,file=trim(filename1),form='FORMATTED')
do i=0,ngridk(1)-1
  do j=0,ngridk(2)-1
    do k=0,ngridk(3)-1
      write(50,'(3I4.2,I8.6)') i,j,k,ivkik(i,j,k)
    end do
  end do
end do
close(50)
!
write(filename2,'("IVKIKNR.OUT")') 
open(51,file=trim(filename2),form='FORMATTED')
do i=0,ngridk(1)-1
  do j=0,ngridk(2)-1
    do k=0,ngridk(3)-1
      write(51,'(3I4.2,I8.6)') i,j,k,ivkiknr(i,j,k)
    end do
  end do
end do
close(51)
!
write(*,*)
write(*,'("Info(writeivkinfo) : locations of k-points on integer grid")' )
!
write(filename3,'("IVK.OUT")') 
open(52,file=trim(filename3),form='FORMATTED')
do ik=1,nkptnr
  write(52,'(3I4.2)') ivk(:,ik)
end do
close(52)
!$OMP END CRITICAL(u180)
return
end subroutine