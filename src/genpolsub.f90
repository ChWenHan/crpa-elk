
! Copyright (C) 2010 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genpolsub(poldr,t3hw,ik,lock,ssr,vqpl,gclgq,jlgqr,ylmgq,sfacgq,nm,vchi0)
use modmain
use modomp
use modwanprj
implicit none
! local variables
logical, intent(in) :: t3hw
integer, intent(in) :: ik,poldr
integer(8), intent(in) :: lock(nwrf)
real(8), intent(in) :: ssr,vqpl(3),gclgq(ngrf)
real(8), intent(in) :: jlgqr(njcmax,nspecies,ngrf)
complex(8), intent(in) :: ylmgq(lmmaxo,ngrf)
complex(8), intent(in) :: sfacgq(ngrf,natmtot)
integer, intent(in) :: nm
complex(8), intent(inout) :: vchi0(nm,nm,nwrf)
! local variables
logical tq0,tssr
integer isym,jk,jkq,iw
integer nst,nstq,ist,jst,kst,lst,nstsub,nstqsub
integer nstall, nstqall
integer nm2,ig,jg,i,j,nthd
real(8) vkql(3),ei,ej,eij,t1,t2
complex(8) a(3,3),z1
! automatic arrays
integer idx(nstsv),idxq(nstsv)
integer idxall(nstsv),idxqall(nstsv)
integer idxsub(nstsv),idxqsub(nstsv)
integer ngp(nspnfv),ngpq(nspnfv)
! allocatable arrays
integer, allocatable :: igpig(:,:),igpqig(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:),wfirq(:,:,:)
complex(8), allocatable :: zrhomt(:,:),zrhoir(:),zrhoig(:)
complex(8), allocatable :: pmat(:,:,:),zw(:),b(:,:)
! check if q=0
tq0=.false.
if (sum(abs(vqpl(:))).lt.epslat) tq0=.true.
! k+q-vector in lattice coordinates
vkql(:)=vkl(:,ik)+vqpl(:)
! equivalent reduced k-points for k and k+q
jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
call findkpt(vkql,isym,jkq)
! count and index states at k and k+q in energy window
nst=0
nstsub=0
nstall=0
do ist=1,nstsv
  if (abs(evalsv(ist,jk)-efermi).lt.emaxrf) then
    nstall=nstall+1
    idxall(nstall)=ist
    if (poldr.eq.0) then ! pick in the window/P_d
      if(((evalsv(ist,jk)-efermi).ge.emincor) .and. & 
        ((evalsv(ist,jk)-efermi).le.emaxcor)) then
          nstsub=nstsub+1
          nst=nstsub
          idxsub(nstsub)=ist
          idx(nst)=ist
        end if
    elseif (poldr.eq.1) then ! pick not in the window P(but not contain r-d transitions)
      if(((evalsv(ist,jk)-efermi).lt.emincor) .or. & 
        ((evalsv(ist,jk)-efermi).gt.emaxcor)) then
          nstsub=nstsub+1
          nst=nstsub
          idxsub(nstsub)=ist
          idx(nst)=ist
      end if
    elseif (poldr.eq.2) then
      nst=nstall
      idx(nst)=ist
    end if
  end if
end do
nstqall=0
nstqsub=0
do jst=1,nstsv
  if (abs(evalsv(jst,jkq)-efermi).lt.emaxrf) then
    nstqall=nstqall+1
    idxqall(nstqall)=jst
    if (poldr.eq.0) then ! pick in the window/P_d
      if(((evalsv(jst,jkq)-efermi).ge.emincor) .and. & 
        ((evalsv(jst,jkq)-efermi).le.emaxcor)) then
          nstqsub=nstqsub+1
          nstq=nstqsub
          idxqsub(nstqsub)=jst
          idxq(nstq)=jst
        end if
    elseif (poldr.eq.1) then ! pick not in the window P(but not contain r-d transitions)
      if(((evalsv(jst,jkq)-efermi).lt.emincor) .or. & 
        ((evalsv(jst,jkq)-efermi).gt.emaxcor)) then
          nstqsub=nstqsub+1
          nstq=nstqsub
          idxqsub(nstqsub)=jst
          idxq(nstq)=jst
      end if
    elseif(poldr.eq.2) then
      nstq=nstqall
      idxq(nstqall)=jst
    end if
  end if
end do

! write(6,*)'in genpolsub, idx is'
! write(6,*)idx(1:nst)



! do ist=1,nstsv
!   if (abs(evalsv(ist,jk)-efermi).lt.emaxrf) then
!     nst=nst+1
!     idx(nst)=ist
!   end if
! end do
! nstq=0
! do jst=1,nstsv
!   if (abs(evalsv(jst,jkq)-efermi).lt.emaxrf) then
!     nstq=nstq+1
!     idxq(nstq)=jst
!   end if
! end do
! set flag for scissor operator
if (abs(ssr).gt.1.d-8) then
  tssr=.true.
else
  tssr=.false.
end if
! generate the wavefunctions for all states at k and k+q in energy window
allocate(igpig(ngkmax,nspnfv))
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfir(ngtc,nspinor,nst))
call genwfsvp(.false.,.false.,nst,idx,ngdgc,igfc,vkl(:,ik),ngp,igpig,wfmt,ngtc,&
  wfir)
deallocate(igpig)
allocate(igpqig(ngkmax,nspnfv))
allocate(wfmtq(npcmtmax,natmtot,nspinor,nstq),wfirq(ngtc,nspinor,nstq))
call genwfsvp(.false.,.false.,nstq,idxq,ngdgc,igfc,vkql,ngpq,igpqig,wfmtq,ngtc,&
  wfirq)
deallocate(igpqig)
! read the momentum matrix elements from file for q=0
if (tq0) then
  allocate(pmat(nstsv,nstsv,3))
  call getpmat(vkl(:,ik),pmat)
! divide by unit cell volume
  t1=1.d0/omega
  pmat(:,:,:)=t1*pmat(:,:,:)
end if
nm2=nm**2
call holdthd(nst,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhomt,zrhoir,zrhoig,zw,b) &
!$OMP PRIVATE(jst,kst,lst,ei,ej,eij,t1,t2) &
!$OMP PRIVATE(iw,ig,jg,z1,i,j,a) &
!$OMP NUM_THREADS(nthd)
allocate(zrhomt(npcmtmax,natmtot),zrhoir(ngtc))
allocate(zrhoig(ngrf),zw(nwrf))
if (tq0.and.t3hw) then
  allocate(b(-1:ngrf,-1:ngrf))
else
  allocate(b(ngrf,ngrf))
end if
!$OMP DO
do ist=1,nst
  kst=idx(ist)
  ei=evalsv(kst,jk)
  do jst=1,nstq
    lst=idxq(jst)
    t1=wkptnr*omega*(occsv(kst,jk)-occsv(lst,jkq))
    if (abs(t1).lt.1.d-8) cycle
    ej=evalsv(lst,jkq)
    eij=ei-ej
! scissor operator
    if (tssr) then
      t2=eij
      if ((ei.gt.efermi).and.(ej.le.efermi)) then
        eij=eij+ssr
      else if ((ei.le.efermi).and.(ej.gt.efermi)) then
        eij=eij-ssr
      end if
      t2=eij/t2
! scale the momentum matrix elements for q=0
      if (tq0) pmat(kst,lst,1:3)=t2*pmat(kst,lst,1:3)
    end if
! frequency-dependent part in response function formula for all frequencies
    do iw=1,nwrf
      zw(iw)=t1/(eij+wrf(iw))
    end do
! compute the complex density in G+q-space
    call genzrho(.true.,.true.,ngtc,wfmt(:,:,:,ist),wfir(:,:,ist), &
      wfmtq(:,:,:,jst),wfirq(:,:,jst),zrhomt,zrhoir)
    call zftzf(ngrf,jlgqr,ylmgq,ngrf,sfacgq,zrhomt,zrhoir,zrhoig)
! Hermitian part of body
    do jg=1,ngrf
      do ig=1,jg-1
        b(ig,jg)=conjg(b(jg,ig))
      end do
      z1=gclgq(jg)*conjg(zrhoig(jg))
      do ig=jg,ngrf
        b(ig,jg)=gclgq(ig)*zrhoig(ig)*z1
      end do
    end do
! case of q=0
    if (tq0) then
      if (t3hw) then
        b(-1:1,-1:1)=0.d0
! calculate 3 x ngrf wings of matrix
        t1=-sqrt(fourpi)/eij
        do i=-1,1
          z1=t1*pmat(kst,lst,i+2)
          b(i,2:)=z1*conjg(zrhoig(2:))*gclgq(2:)
          b(2:,i)=conjg(b(i,2:))
        end do
      else
! use trace of 3 x 3 head of matrix
        t1=sum(dble(pmat(kst,lst,1:3))**2+aimag(pmat(kst,lst,1:3))**2)/3.d0
        b(1,1)=(fourpi/eij**2)*t1
! wings of matrix
        t1=-sqrt(fourpi)/eij
        z1=(t1/3.d0)*(pmat(kst,lst,1)+pmat(kst,lst,2)+pmat(kst,lst,3))
        b(1,2:)=z1*conjg(zrhoig(2:))*gclgq(2:)
        b(2:,1)=conjg(b(1,2:))
      end if
    end if
! add to body and wings of the response function
    do iw=1,nwrf
      call omp_set_lock(lock(iw))
      call zaxpy(nm2,zw(iw),b,1,vchi0(1,1,iw),1)
      call omp_unset_lock(lock(iw))
    end do
! calculate 3 x 3 head
    if (tq0.and.t3hw) then
      t1=-fourpi/eij
      zw(1:nwrf)=zw(1:nwrf)/wrf(1:nwrf)
      do j=1,3
        do i=1,3
          a(i,j)=t1*pmat(kst,lst,i)*conjg(pmat(kst,lst,j))
        end do
      end do
! add to the head of the response function
      do iw=1,nwrf
        call omp_set_lock(lock(iw))
        vchi0(1:3,1:3,iw)=vchi0(1:3,1:3,iw)+a(1:3,1:3)*zw(iw)
        call omp_unset_lock(lock(iw))
      end do
    end if
! end loop over jst
  end do
! end loop over ist
end do
!$OMP END DO
deallocate(zrhomt,zrhoir,zrhoig,zw,b)
!$OMP END PARALLEL
call freethd(nthd)
deallocate(wfmt,wfir,wfmtq,wfirq)
if (tq0) deallocate(pmat)
end subroutine

