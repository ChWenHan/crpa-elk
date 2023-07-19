subroutine wantest
use modmain
use modwanprj
!
implicit none
!
integer iq,ik,nkqst,nkst,ikq,isym,idxkq(nstsv),idxk(nstsv),jk
integer ikd,ikqd
real(8) vkql(3)
character(256) fpre
!
complex(8) ,allocatable :: wanprjkq(:,:,:,:,:),wanprjk(:,:,:,:,:)
!

call initwanprj_globalvars
allocate(wanprjkq(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj))
allocate(wanprjk(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj))
!
write(6,*)'wantest init finish'
!
open(50,file='WANTEST.OUT',form='formatted')
do iq=1,nqpt
  do ik=1,nkptnr
    !write(6,'(2I6.4)')ik,iq
    vkql=vkl(:,ik)+vql(:,iq)
    call findkpt(vkql,isym,ikq)
    call findkpt(vkl(:,ik),isym,jk)
    nkst=projstwanprj(jk)
    nkqst=projstwanprj(ikq)
    idxkq(1:nkqst)=idxwanprj(1:nkqst,ikq)
    idxk(1:nkst)=idxwanprj(1:nkst,jk)
    call wanprojkpl(vkql,idxkq(1:nkqst),ldwanprj,nprojwanprj,maxnstwanprj,nkqst,&
       subulmwanprj,sublmwanprj,symlmwanprj,wanprjkq)
    !write(6,*)wanprjkq(1,1,1,1,1)
    write(fpre,'("_IQ",I6.6)')iq
    call writewanprjk_format(ik,vkql,fpre,wanprjkq)
    ikd=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
    ikqd=ivkik(ivk(1,ikq),ivk(2,ikq),ivk(3,ikq))
    write(50,'(6I4.3,3G18.10,2I3.2)')iq,ik,ikq,jk,ikd,ikqd,vkql
    call wanprojkpl(vkl(:,ik),idxk(1:nkst),ldwanprj,nprojwanprj,maxnstwanprj,nkst,&
    subulmwanprj,sublmwanprj,symlmwanprj,wanprjk)
    write(fpre,'("_IQLOOP",I6.6)')iq
    call writewanprjk_format(ik,vkl(:,ik),fpre,wanprjk)
  end do
end do
close(50)
!
deallocate(wanprjkq)
return
end subroutine