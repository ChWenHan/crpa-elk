subroutine writewanrho
use modmain
!use modgw
use modomp
use modwanprj
!use modmpi
implicit none
!
integer iq,ik,iw
integer nthd
logical twzrkq
!
complex(8),allocatable :: mmm(:,:,:,:,:,:),mmmq(:,:,:,:,:)
complex(8), allocatable :: zzrhokq(:,:,:,:,:)

!
call initwanprj_globalvars
twzrkq=.false.
whchtagn=10
allocate(mmm(ldwanprj,ldwanprj,nspinor,nspinor,ngrf,nqpt))
mmm=0.d0
do iq=1,nqpt
  !write(*,*)'writemnn',iq
!$OMP CRITICAL(writemnn_1)
  write(*,'("Info(writewanrho): ",I6," of ",I6," q-points")') iq,nqpt
!$OMP END CRITICAL(writemnn_1)
  call holdthd(nkptnr,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zzrhokq,mmmq) &
!$OMP NUM_THREADS(nthd)
  allocate(zzrhokq(nstsv,nstsv,nspinor,nspinor,ngrf))
  allocate(mmmq(ldwanprj,ldwanprj,nspinor,nspinor,ngrf))
!$OMP DO
  do ik=1,nkptnr
!$OMP CRITICAL(writemnn_1)
!    write(*,'("Info(writewanrho): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(writemnn_1)
    call genzrhokq(iq,vkl(:,ik),vql(:,iq),zzrhokq)
    if (twzrkq) then
      call writezrhokq(ik,iq,vkl(:,ik),vql(:,iq),zzrhokq)
    end if
    call genwanrhokq_matop(vkl(:,ik),vql(:,iq),zzrhokq,mmmq)
!$OMP CRITICAL(writemnn_1)
    mmm(:,:,:,:,:,iq)=mmmq*wkptnr+mmm(:,:,:,:,:,iq)
!$OMP END CRITICAL(writemnn_1)
  end do
!$OMP END DO
  deallocate(zzrhokq,mmmq)
!$OMP END PARALLEL
  call freethd(nthd)
!$OMP CRITICAL(writemnn_1)
  call writemmmq(0,iq,mmm(:,:,:,:,:,iq))
!$OMP END CRITICAL(writemnn_1)
end do
deallocate(mmm)
return
end subroutine