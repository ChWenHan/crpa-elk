
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genepsinvrest
use modmain
use modmpi
use modomp
use modwanprj
implicit none
! local variables
integer iq,ik,ig,iw
integer n,nthd
character(256) fpre
! allocatable arrays
integer(8), allocatable :: lock(:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:),epsi(:,:,:),epsid(:,:,:)
complex(8), allocatable :: epsi_standard(:,:,:)
! allocate local arrays
allocate(vgqc(3,ngrf),gqc(ngrf),gclgq(ngrf))
allocate(jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(epsi(ngrf,ngrf,nwrf),epsid(ngrf,ngrf,nwrf),epsi_standard(ngrf,ngrf,nwrf))
! initialise the OpenMP locks
allocate(lock(nwrf))
do iw=1,nwrf
  call omp_init_lock(lock(iw))
end do
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! loop over q-points
do iq=1,nqpt
  if (mp_mpi) write(*,'("Info(genepsinvrest): ",I6," of ",I6," q-points")') iq,nqpt
! generate the G+q-vectors and related functions
  call gengqf(ngrf,vqc(:,iq),vgqc,gqc,jlgqr,ylmgq,sfacgq)
  !write(6,*)'in genepsinvrest, after gengqf'
  !write(6,*)vqc(:,iq)
! generate the regularised Coulomb Green's function in G+q-space
  call gengclgq(.true.,iq,ngrf,gqc,gclgq)
  !write(6,*)'in genepsinvrest, after gengclgq'
! use the symmetric form
  gclgq(:)=sqrt(gclgq(:))
! zero the response function (stored in epsi)
  epsi(:,:,:)=0.d0
  epsid(:,:,:)=0.d0
  epsi_standard(:,:,:)=0.d0
  !write(6,*)'in genepsinvrest, after init'
  call holdthd(nkptnr/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do ik=1,nkptnr
! distribute among MPI processes
    if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
! compute v^1/2 chi0 v^1/2
    !for Pol_d
    call genpolsub(0,.false.,ik,lock,0.d0,vql(:,iq),gclgq,jlgqr,ylmgq,sfacgq, &
      ngrf,epsid)
    !for Pol_all
    call genpolsub(2,.false.,ik,lock,0.d0,vql(:,iq),gclgq,jlgqr,ylmgq,sfacgq, &
      ngrf,epsi)
    ! call genvchi0(.false.,ik,lock,0.d0,vql(:,iq),gclgq,jlgqr,ylmgq,sfacgq, &
    !   ngrf,epsi_standard)
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! add epsi from each process and redistribute
  if (np_mpi.gt.1) then
    n=ngrf*ngrf*nwrf
    call mpi_allreduce(mpi_in_place,epsi,n,mpi_double_complex,mpi_sum,mpicom, &
      ierror)
  end if
  fpre='POL'
  call writeepsinvrest(iq,fpre,epsi)
  ! fpre='POL_STANDARD'
  ! call writeepsinvrest(iq,fpre,epsi_standard)
  fpre='POLD'
  call writeepsinvrest(iq,fpre,epsid)
  epsi=epsi(:,:,:)-epsid(:,:,:)
  fpre='POLR'
  call writeepsinvrest(iq,fpre,epsi)
  ! negate and add delta(G,G')
  epsi(:,:,:)=-epsi(:,:,:)
  do ig=1,ngrf
    epsi(ig,ig,:)=epsi(ig,ig,:)+1.d0
  end do
!-------------------------------------!
!     invert epsilon over G-space     !
!-------------------------------------!
  call holdthd(nwrf,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
  do iw=1,nwrf
    call zminv(ngrf,epsi(:,:,iw))
  end do
!$OMP END PARALLEL DO
  call freethd(nthd)
! write inverse RPA epsilon to EPSINV.OUT
  if (mp_mpi) call putepsinvr(iq,epsi)
  fpre='EPSIR'
  call writeepsinvrest(iq,fpre,epsi)
! end loop over q-points
end do
! destroy the OpenMP locks
do iw=1,nwrf
  call omp_destroy_lock(lock(iw))
end do
deallocate(lock)
deallocate(vgqc,gqc,gclgq,jlgqr)
deallocate(ylmgq,sfacgq,epsi,epsid,epsi_standard)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

