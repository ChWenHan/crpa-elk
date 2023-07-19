subroutine writesubusymlm
!
use modmain
use modomp
use modwanprj
!
character(256) fname
integer ld1,iorb,iproj,ld2
!
!$OMP CRITICAL(writemnn_1)
fname='SUBLM.OUT'
open(888,file=fname,form='FORMATTED')
do ld1=1,ldwanprj
  do iorb=1,norb
    do iproj=1,nprojwanprj
      write(888,'(3I3.2,I4.3)')ld1,iorb,iproj,sublmwanprj(ld1,iorb,iproj)
    end do
  end do
end do
close(888)
!
fname='SUBULM.OUT'
open(888,file=fname,form='FORMATTED')
do ld1=1,ldwanprj
  do ld2=1,ldwanprj
    do iorb=1,norb
      do iproj=1,nprojwanprj
        write(888,'(4I3.2,2G18.10)')ld1,ld2,iorb,iproj,subulmwanprj(ld1,ld2,iorb,iproj)
      end do
    end do
  end do
end do
close(888)
!
fname='SYMLM.OUT'
open(888,file=fname,form='FORMATTED')
do ld1=1,ldwanprj
  do ld2=1,ldwanprj
    do iorb=1,norb
      do iproj=1,nprojwanprj
        write(888,'(4I3.2,2G18.10)')ld1,ld2,iorb,iproj,symlmwanprj(ld1,ld2,iorb,iproj)
      end do
    end do
  end do
end do
close(888)
!$OMP END CRITICAL(writemnn_1)
return
end subroutine