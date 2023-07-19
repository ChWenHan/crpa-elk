subroutine writemmmq(ik,iq,mmmq)
use modmain
use modomp
use modwanprj
!
integer(4), intent(in) :: ik,iq
complex(8), intent(in) :: mmmq(ldwanprj,ldwanprj,nspinor,nspinor,ngrf)
!
character(256) fname
integer ld1,ld2,ispn,jspn,ig
!
!$OMP CRITICAL(writemnn_50)
write(fname,'("MMMQ_K",I6.6,"_Q",I6.6,".OUT")')ik,iq
write(6,*)fname
open(50,file=trim(fname),form='FORMATTED')
do ld1=1,ldwanprj
  do ld2=1,ldwanprj
    do ispn=1,nspinor
      do jspn=1,nspinor
        do ig=1,ngrf
          write(50,'(2I3.2,2I2.1,I6.5,2G18.10)')ld1,ld2,ispn,jspn,ig,mmmq(ld1,ld2,ispn,jspn,ig)
        end do
      end do
    end do
  end do
end do
close(50)
!$OMP END CRITICAL(writemnn_50)
return
end subroutine