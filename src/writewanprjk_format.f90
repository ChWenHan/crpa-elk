subroutine writewanprjk_format(ikp,vkpl,fpre,wanprj)
!
use modmain
use modomp
use modwanprj 
!
implicit none
!
integer, intent(in) :: ikp
real(8), intent(in) :: vkpl(3)
character(256), intent(in) :: fpre
complex(8), intent(in) :: wanprj(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj)
!
character(256) fname
integer iorb,is,l,lmmax,i,j,ia,ias,ld1,ist1,ispn,ik,isym
!$OMP CRITICAL(writewanprj_1)
call findkpt(vkpl,isym,ik)
write(fname,'("WANPRJ_K",I6.6)')ikp
open(888,file=trim(fname)//trim(fpre)//'.OUT',form='FORMATTED',STATUS='REPLACE')
do iorb=1,norb
    is=orb(iorb,1)
    l=orb(iorb,2)
    lmmax=2*l+1
    i=l**2+1
    j=(l+1)**2
    do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ld1=1,ldwanprj
            do ist1=1,projstwanprj(ik)
                do ispn=1,nspinor
                    write(888,'(I3.2,I4.3,I2.1,I2.1,I3.2,2G18.10)')ld1,ist1,ispn,iorb,ia,&
                    &wanprj(ld1,ist1,ispn,iorb,ia)
                end do
            end do
        end do
    end do
end do
close(888)
!$OMP END CRITICAL(writewanprj_1)
return
end subroutine