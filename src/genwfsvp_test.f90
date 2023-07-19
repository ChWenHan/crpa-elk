subroutine genwfsvp_test
! calculate the matrix of <ψ_{k+q,n,σ}|exp(-i·(q+G)·r))|ψ_{k,n',σ'}>
use modmain
!use modgw
use modomp
!use modwanprj
implicit none
! arguments
!integer, intent(in) :: ikp
! real(8),intent(in) :: vkpl(3),vqpl(3)
! complex(8), intent(out) :: zrhokq(nstsv,nstsv,nspinor,nspinor,ngrf)
! local variables
integer ik,jk,iq,ist1,ist2
integer ikq,jkq,ispn,jspn
integer isym
integer nst,nthd
integer ngpk(nspnfv),ngpkq(nspnfv)
real(8) vkql(3),vc(3),t1,t2
real(8) vkpl(3),vqpl(3)
complex(8) zrhokq(nstsv,nstsv,nspinor,nspinor,ngrf)
real(8) coeff
complex(8) z1,z2,z
character(256) fpre,fname1,fname2,fname,fnamewan
! automatic arrays
integer idx(nstsv)
! allocatable arrays
!integer(8), allocatable :: lock(:)
integer, allocatable :: igpigk(:,:),igpigkq(:,:)
real(8), allocatable :: vgqc(:,:),gqc(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: wfmtk(:,:,:,:),wfirk(:,:,:)
complex(8), allocatable :: wfmtkq(:,:,:,:),wfirkq(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:,:)
complex(8), allocatable :: zfgq(:),zrho(:,:,:,:)

!complex(8), allocatable :: gs(:,:),wc(:,:),zv(:)
! external functions
complex(8) zfinp
external zfinp
call init0
call init1
call init2
call init3
! allocate local arrays
allocate(jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(vgqc(3,ngrf),gqc(ngrf))
allocate(wfmtkq(npcmtmax,natmtot,nspinor,nstsv),wfirkq(ngtc,nspinor,nstsv))
allocate(wfmtk(npcmtmax,natmtot,nspinor,nstsv),wfirk(ngtc,nspinor,nstsv))
allocate(igpigk(ngkmax,nspnfv),igpigkq(ngkmax,nspnfv))
vqpl=0.d0
vkpl=0.d0
whchtagn=0
tevecsv=.true.
write(6,*)'in genwfsvp_test,'
write(6,*)'input set vqpl',vqpl
write(6,*)'input set vvkpl',vkpl
!find out q-point index 
call findqpt(vqpl,isym,iq)
write(6,*)'findqpt vql maping to iq',iq
!generate q related quantities
call gengqf(ngrf,vqc(:,iq),vgqc,gqc,jlgqr,ylmgq,sfacgq)
! k+q-vector in lattice coordinates
call findkpt(vkpl,isym,jk)
vkql(:)=vkpl(:)+vqpl(:)
call findkpt(vkql,isym,jkq)
!$OMP CRITICAL(u1)
write(6,'("input kvec map to kind :",I6.5," k+q map to kind :",I6.5,".OUT")')jk,jkq
!$OMP END CRITICAL(u1)
! equivalent reduced k-points for k and k+q
ik=ivkik(ivk(1,jk),ivk(2,jk),ivk(3,jk))
ikq=ivkik(ivk(1,jkq),ivk(2,jkq),ivk(3,jkq))
! index to all states
nst=nstsv
do ist1=1,nstsv
  idx(ist1)=ist1
end do
!$OMP CRITICAL(u1)
write(6,'("Gen wavef for eq irk :",I6.5,"and eq ir-k+q :",I6.5,".OUT")')ik,ikq
write(6,*)' ' 
!$OMP END CRITICAL(u1)
! calculate the wavefunctions for all states of the input k-point
call genwfsvp(.false.,.false.,nst,idx,ngdgc,igfc,vkl(:,ik),ngpk,igpigk,wfmtk,ngtc, &
  wfirk)
!$OMP CRITICAL(u1)
write(6,'("value wfmtk111 irk :",I6.5,"and eq ir-k+q :",I6.5,".OUT")')ik,ikq
write(6,*)wfmtk(1,1,1,1) 
write(6,*)'vkl for gen wfmtk:',vkl(:,ik)
!$OMP END CRITICAL(u1)
! calculate the wavefunctions for all states of the input k+q-point
call genwfsvp(.false.,.false.,nst,idx,ngdgc,igfc,vkl(:,ikq),ngpkq,igpigkq,wfmtkq,ngtc, &
  wfirkq)
!$OMP CRITICAL(u1)
write(6,'("value wfmtkq 111 irk :",I6.5,"and eq ir-k+q :",I6.5,".OUT")')ik,ikq
write(6,*)wfmtkq(1,1,1,1) 
write(6,*)'vkql for gen wfmtkq:',vkl(:,iq)
!$OMP END CRITICAL(u1)
! !write(*,*)'in genzrhokq flag1'
! ! loop over non-reduced k-point set
! zrhokq=complex(0.d0,0.d0)
! do ist1=1,nstsv
! ! calculate the matrix of <ψ_{k+q,n,σ}|exp(-i·(q+G)·r))|ψ_{k,n',σ'}>
!   call holdthd(nstsv,nthd)
! !$OMP PARALLEL DEFAULT(SHARED) &
! !$OMP PRIVATE(zrhomt,zrhoir,zfgq) &
! !$OMP PRIVATE(ist2,ispn,jspn) &
! !$OMP NUM_THREADS(nthd)
!   allocate(zrhomt(npcmtmax,natmtot,nstsv),zrhoir(ngtc,nstsv))
!   allocate(zfgq(ngrf))
! !$OMP DO
!   do ist2=1,nstsv
!     do ispn=1,nspinor
!       do jspn=1,nspinor
!         call genzrho(.true.,.false.,ngtc,wfmtkq(:,:,ispn,ist1),wfirkq(:,ispn,ist1), &
!         wfmtk(:,:,jspn,ist2),wfirk(:,jspn,ist2),zrhomt(:,:,ist2),zrhoir(:,ist2))
!         call zftzf(ngrf,jlgqr,ylmgq,ngvc,sfacgq,zrhomt(:,:,ist2),zrhoir(:,ist2), &
!         zfgq)
!         zrhokq(ist1,ist2,ispn,jspn,:)=conjg(zfgq(:))
!       end do
!     end do
!   end do
! !$OMP END DO
!   deallocate(zrhomt,zrhoir,zfgq)
! !$OMP END PARALLEL
!   call freethd(nthd)
! end do
! !$OMP CRITICAL(u1)
! write(6,'("value zrhokq11111 irk :",I6.5,"and eq ir-k+q :",I6.5,".OUT")')ik,ikq
! write(6,*)zrhokq(1,1,1,1,1) 
! !$OMP END CRITICAL(u1)
deallocate(wfmtk,wfirk,wfmtkq,wfirkq,jlgqr,ylmgq,sfacgq)
!deallocate(zfgq,zrhoir,zrhomt)
deallocate(vgqc,gqc)
deallocate(igpigk,igpigkq)
return
end subroutine