subroutine genzrhokq(iqp,vkpl,vqpl,zzrhokq)
! calculate the matrix of <ψ_{k+q,n,σ}|exp(-i·(q+G)·r))|ψ_{k,n',σ'}>
use modmain
!use modgw
use modomp
!use modwanprj
implicit none
! arguments
integer, intent(in) :: iqp
real(8),intent(in) :: vkpl(3),vqpl(3)
complex(8), intent(out) :: zzrhokq(nstsv,nstsv,nspinor,nspinor,ngrf)
! local variables
integer ik,jk,iq,ist1,ist2
integer ikq,jkq,ispn,jspn
integer isym
integer nst,nthd
integer ngpk(nspnfv),ngpkq(nspnfv)
real(8) vkql(3),vc(3),t1,t2
real(8) coeff
complex(8) z1,z2,z
character(256) fpre,fname1,fname2,fname,fnamewan
integer(4) ERR_MESSAGE
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
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(jlgqr(njcmax,nspecies,ngrf))
if (allocated(ylmgq)) then 
!  write(6,*) 'ylmgq alloc'
  deallocate(ylmgq)
else
!  write(6,*) 'ylmgq not alloc'
end if
ALLOCATE(ylmgq(1:lmmaxo,1:ngrf),STAT=ERR_MESSAGE)
IF(ERR_MESSAGE.NE.0) PRINT *,'ALLOCATION ERROR'
allocate(sfacgq(ngrf,natmtot))
allocate(vgqc(3,ngrf),gqc(ngrf))
allocate(wfmtkq(npcmtmax,natmtot,nspinor,nstsv),wfirkq(ngtc,nspinor,nstsv))
allocate(wfmtk(npcmtmax,natmtot,nspinor,nstsv),wfirk(ngtc,nspinor,nstsv))
allocate(igpigk(ngkmax,nspnfv),igpigkq(ngkmax,nspnfv))
!find out q-point index, this is for dbug. Sometime the default int type is dif thus may cause bug
call findqpt(vqpl,isym,iq)
if (iq.ne.iqp) then
  write(6,'("Error(genzrhokq) Q point wrong for IQ: ",I6.6)')iqp
end if

!generate G+q related quantities
! gqc,jlgqr,ylmgq,sfacgq
call gengqf(ngrf,vqc(:,iq),vgqc,gqc,jlgqr,ylmgq,sfacgq)
! k+q-vector in lattice coordinates
call findkpt(vkpl,isym,jk)
vkql(:)=vkpl(:)+vqpl(:)
call findkpt(vkql,isym,jkq)
! equivalent reduced k-points for k and k+q/not used, just for check
ik=ivkik(ivk(1,jk),ivk(2,jk),ivk(3,jk))
ikq=ivkik(ivk(1,jkq),ivk(2,jkq),ivk(3,jkq))
!write(6,'("Info(genzrhokq) iqp,ikp: ",2I6.4)')iqp,ik
! index to all states/will set it to corrlated states 
nst=nstsv
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! calculate the wavefunctions for all states of the input k-point
call genwfsvp(.false.,.false.,nst,idx,ngdgc,igfc,vkpl,ngpk,igpigk,wfmtk,ngtc, &
 wfirk)
! calculate the wavefunctions for all states of the input k+q-point
call genwfsvp(.false.,.false.,nst,idx,ngdgc,igfc,vkql,ngpkq,igpigkq,wfmtkq,ngtc, &
 wfirkq)
! loop over non-reduced k-point set
zzrhokq(:,:,:,:,:)=zzero
do ist1=1,nstsv
! calculate the matrix of <ψ_{k+q,n,σ}|exp(-i·(q+G)·r))|ψ_{k,n',σ'}>
  call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zfgq) &
!$OMP PRIVATE(ist2,ispn,jspn) &
!$OMP NUM_THREADS(nthd)
  allocate(zrhomt(npcmtmax,natmtot,nstsv),zrhoir(ngtc,nstsv))
  allocate(zfgq(ngrf))
!$OMP DO
  do ist2=1,nstsv
    do ispn=1,nspinor
      do jspn=1,nspinor
        call genzrho(.true.,.false.,ngtc,wfmtkq(:,:,ispn,ist1),wfirkq(:,ispn,ist1), &
        wfmtk(:,:,jspn,ist2),wfirk(:,jspn,ist2),zrhomt(:,:,ist2),zrhoir(:,ist2))
        ! call genzrho(.true.,.false.,ngtc,wfmtkq(:,:,:,ist1),wfirkq(:,:,ist1), &
        ! wfmtk(:,:,:,ist2),wfirk(:,:,ist2),zrhomt(:,:,ist2),zrhoir(:,ist2))
        !matrix of <ψ_{k+q,n,σ}|exp(i·(q+G)·Rα))|ψ_{k,n',σ'}>, so take the conjg
        ! equation for M_G,n,n'
        call zftzf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,zrhomt(:,:,ist2),zrhoir(:,ist2), &
        zfgq)
        !zrhokq(istk,istkq,ispnk,ispnkq)
        zzrhokq(ist1,ist2,jspn,ispn,:)=conjg(zfgq(:))
      end do
    end do
  end do
!$OMP END DO
  deallocate(zrhomt,zrhoir,zfgq)
!$OMP END PARALLEL
  call freethd(nthd)
end do
if (allocated(ylmgq)) deallocate(ylmgq)
deallocate(wfmtk,wfirk,wfmtkq,wfirkq,jlgqr,sfacgq)
!deallocate(zfgq,zrhoir,zrhomt)
deallocate(vgqc,gqc)
deallocate(igpigk,igpigkq)
return
end subroutine