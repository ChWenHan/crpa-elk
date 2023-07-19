subroutine writewanprjcrpa
!
use modmain
use modomp
use modmpi
use modwanprj
!
implicit none
!
integer iorb,ik,is,l,lmmax,i,j,ia,ias,ld1,ist1,ispn
integer ikq,jk,jkq,isym,idxk(nstsv),idxkq(nstsv),iq,nkst,nkqst
real(8) vkql(3)
character(256) fname,fpre
complex(8), allocatable :: wanprj(:,:,:,:,:,:),wanprj_testk(:,:,:,:,:)
complex(8), allocatable :: wanprj_testkq(:,:,:,:,:)
call initwanprj_globalvars
write(6,*)'wan proj used for crpa initialized'
allocate(wanprj(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj,nkpt))
allocate(wanprj_testk(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj))
allocate(wanprj_testkq(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj))
wanprj(:,:,:,:,:,:)=0.d0
! call wannier projector routines 
! wanprj been put to binary file.
call wanproj(ldwanprj,nprojwanprj,maxnstwanprj,idxwanprj,projstwanprj,&
    &sublmwanprj,subulmwanprj,wanprj,nkpt)
if (mp_mpi) then
    do ik=1,nkpt
        write(fname,'("WANPRJ_K",I6.6,".OUT")')ik
        open(888,file=trim(fname),form='FORMATTED',STATUS='REPLACE')
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
                            &wanprj(ld1,ist1,ispn,iorb,ia,ik)
                        end do
                    end do
                end do
            end do
        end do
        close(888)
    end do
end if
! do iq=1,nqpt
!     do jk=1,nkptnr
!         vkql(:)=vkl(:,jk)+vql(:,iq)
!         call findkpt(vkql,isym,jkq)
!         ik=ivkik(ivk(1,jk),ivk(2,jk),ivk(3,jk))
!         ikq=ivkik(ivk(1,jkq),ivk(2,jkq),ivk(3,jkq))
!         nkst=projstwanprj(ik)
!         nkqst=projstwanprj(ikq)
!         idxk(1:nkst)=idxwanprj(1:nkst,ik)
!         idxkq(1:nkqst)=idxwanprj(1:nkqst,ikq)       
!         call wanprojkpl(vkql,idxkq(1:nkqst),ldwanprj,nprojwanprj,maxnstwanprj,nkst,&
!         subulmwanprj,sublmwanprj,symlmwanprj,wanprj_testkq)
!         write(fpre,'("_",I4.4,"_IK",I4.4,"_IQ")')jk,iq
!         call writewanprjk_format(iq,fpre,wanprj_testkq)
!     end do
! end do
! do jk=1,nkptnr
!     call findkpt(vkql,isym,jk)
!     ik=ivkik(ivk(1,jk),ivk(2,jk),ivk(3,jk))
!     nkst=projstwanprj(ik)
!     idxk(1:nkst)=idxwanprj(1:nkst,ik)    
!     call wanprojkpl(vkl(:,jk),idxk(1:nkst),ldwanprj,nprojwanprj,maxnstwanprj,nkst,&
!     subulmwanprj,sublmwanprj,symlmwanprj,wanprj_testk)
!     write(fpre,'("_",I4.4,"_IK")')jk
!     call writewanprjk_format(jk,fpre,wanprj_testk)
! end do



deallocate(wanprj)
return
end subroutine