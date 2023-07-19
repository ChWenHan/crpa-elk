
! Copyright (C) 2011 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwfsvp(tsh,tgp,nst,idx,ngdg,igf,vpl,ngp,igpig,wfmt,ld,wfir)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh,tgp
integer, intent(in) :: nst,idx(*),ngdg(3),igf(*)
real(8), intent(in) :: vpl(3)
integer, intent(out) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(out) :: wfmt(npcmtmax,natmtot,nspinor,nst)
integer, intent(in) :: ld
complex(8), intent(out) :: wfir(ld,nspinor,nst)
! local variables
integer ispn
integer igk,i,imat
real(8) vl(3),vc(3)
! automatic arrays
real(8) vgpl(3,ngkmax,nspnfv),vgpc(3,ngkmax),gpc(ngkmax)
! allocatable arrays
complex(8), allocatable :: sfacgp(:,:),apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
allocate(sfacgp(ngkmax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
! loop over first-variational spins
do ispn=1,nspnfv
  vl(:)=vpl(:)
  vc(:)=bvec(:,1)*vpl(1)+bvec(:,2)*vpl(2)+bvec(:,3)*vpl(3)
! spin-spiral case
  if (spinsprl) then
    if (ispn.eq.1) then
      vl(:)=vl(:)+0.5d0*vqlss(:)
      vc(:)=vc(:)+0.5d0*vqcss(:)
    else
      vl(:)=vl(:)-0.5d0*vqlss(:)
      vc(:)=vc(:)-0.5d0*vqcss(:)
    end if
  end if
! generate the G+p-vectors
  call gengkvec(ngvc,ivg,vgc,vl,vc,gkmax,ngkmax,ngp(ispn),igpig(:,ispn), &
   vgpl(:,:,ispn),vgpc,gpc)
! generate structure factors for G+p-vectors
  call gensfacgp(ngp(ispn),vgpc,ngkmax,sfacgp)
! find the matching coefficients
  call match(ngp(ispn),vgpc,gpc,sfacgp,apwalm(:,:,:,:,ispn))
end do
! get the first- and second-variational eigenvectors from file
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
call getevecfv(filext,0,vpl,vgpl,evecfv)
call getevecsv(filext,0,vpl,evecsv)
! calculate the second-variational wavefunctions
call genwfsv(tsh,tgp,nst,idx,ngdg,igf,ngp,igpig,apwalm,evecfv,evecsv,wfmt,ld, &
 wfir)

! write(ftag,'("vgpc_",I6.6,".OUT")')whchtagn
! whchtagn=whchtagn+1
! open(150,file=trim(ftag),form='FORMATTED')
! do igk=1,ngkmax
!   write(150,'(3G18.10)')vgpc(:,igk)
! end do
! close(150)

! !
! write(ftag,'("gpc_",I6.6,".OUT")')whchtagn
! whchtagn=whchtagn+1
! open(150,file=trim(ftag),form='FORMATTED')
! do igk=1,ngkmax
!   write(150,'(G18.10)')gpc(igk)
! end do
! close(150)
! !
! write(ftag,'("sfacgp_",I6.6,".OUT")')whchtagn
! whchtagn=whchtagn+1
! open(150,file=trim(ftag),form='FORMATTED')
! do igk=1,ngkmax
!   do i=1,natmtot
!     write(150,'(2G18.10)')sfacgp(igk,i)
!   end do
! end do
! close(150)
! write(ftag,'("evecfv_",I6.6,".OUT")')whchtagn
! whchtagn=whchtagn+1
! open(150,file=trim(ftag),form='FORMATTED')
! do imat=1,nmatmax
!   do i=1,nstfv
!     do ispn=1,nspnfv
!       write(150,'(I8.6,2I4.3,2G18.10)')imat,i,ispn,evecfv(imat,i,ispn)
!     end do
!   end do
! end do
! close(150)
! write(ftag,'("evecsv_",I6.6,".OUT")')whchtagn
! whchtagn=whchtagn+1
! open(150,file=trim(ftag),form='FORMATTED')
! do imat=1,nstsv
!   do i=1,nstsv
!     write(150,'(2I4.3,2G18.10)')imat,i,evecsv(imat,i)
!   end do
! end do
! close(150)
! write(6,'("value wfmtk111 om gemwfsvp :",I6.5,"and eq ir-k+q :",I6.5,".OUT")')
! write(6,*)'wftmt1',wfmt(1,1,1,1) 
! write(6,*)'nst1',nst
! write(6,*)'idx1',idx(1)
! write(6,*)'ngdg',ngdg
! write(6,*)'igf',igf(1)
! write(6,*)'ngp',ngp(1)
! write(6,*)'igpig',igpig(1,1)
! write(6,*)'evecfv',evecfv(1,1,1)
! write(6,*)'evecsv',evecsv(1,1)
! write(6,*)'apwalm',apwalm(1,1,1,1,1)

deallocate(sfacgp)
deallocate(apwalm,evecfv,evecsv)
end subroutine

