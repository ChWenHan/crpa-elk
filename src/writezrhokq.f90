subroutine writezrhokq(ikp,iqp,vkpl,vqpl,zrhokq)
! calculate the matrix of <ψ_{k+q,n,σ}|exp(-i·(q+G)·r))|ψ_{k,n',σ'}>
use modmain
use modomp
implicit none
!
integer,intent(in) :: ikp,iqp
real(8),intent(in) :: vkpl(3),vqpl(3)
complex(8), intent(in) :: zrhokq(nstsv,nstsv,nspinor,nspinor,ngrf)
!local vars
integer ik,jk,iq,ist1,ist2,ig
integer ikq,jkq,ispn,jspn,isym
character(256) fname,fpre
real(8) vkql(3)
!
call findqpt(vqpl,isym,iq)
call findkpt(vkpl,isym,jk)
vkql(:)=vkpl(:)+vqpl(:)
call findkpt(vkql,isym,jkq)
! equivalent reduced k-points for k and k+q
ik=ivkik(ivk(1,jk),ivk(2,jk),ivk(3,jk))
ikq=ivkik(ivk(1,jkq),ivk(2,jkq),ivk(3,jkq))
!$OMP CRITICAL(writezrho)
write(fname,'("ZRHO_IK",I6.6,"_IQ",I6.6,".OUT")')ikp,iqp
!write(*,*)'in writezrhokq',ikp,iqp,fname
open(888,file=trim(fname),form='FORMATTED')
do ist1=1,nstsv
  do ist2=1,nstsv
    do ispn=1,nspinor
      do jspn=1,nspinor
        do ig=1,ngrf
          write(888,'(2I4.3,2I2.1,I6.5,2G18.10)')ist1,ist2,ispn,jspn,ig,zrhokq(ist1,ist2,ispn,jspn,ig)
        end do
      end do
    end do
  end do
end do
close(888)
!$OMP END CRITICAL(writezrho)
return
end subroutine