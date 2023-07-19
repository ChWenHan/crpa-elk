subroutine writecrpaumat
use modmain
!use modgw
use modomp
use modwanprj
use modmpi
implicit none
! local variables
integer iq,ik,iw,ild1,ild2,ild3,ild4
integer nthd
integer ld1,ld2,ispn,jspn
logical twzrkq
character(256) fpre
complex(8), allocatable :: zrhokq(:,:,:,:,:),zrhoq(:,:,:,:,:,:),zrho(:,:,:,:,:)
complex(8),allocatable :: mmm(:,:,:,:,:,:),mmmq(:,:,:,:,:)
complex(8),allocatable :: umatq(:,:,:,:,:)
complex(8), allocatable :: wanprj(:,:,:,:,:,:)
complex(8), allocatable :: umat(:,:,:,:,:)
complex(8),allocatable :: uu(:),u(:),jj(:),j(:),u2ind(:,:,:),j2ind(:,:,:)

! call init0
! call init1
! call init2
! call init3
call initwanprj_globalvars
twzrkq=.false.
!
! call initwanprj
! allocate(wanprj(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj,nkpt))
! wanprj(:,:,:,:,:,:)=0.d0
! ! call wannier projector routines 
! ! wanprj been put to binary file.
! call wanproj(ldwanprj,nprojwanprj,maxnstwanprj,idxwanprj,projstwanprj,&
!     &sublmwanprj,subulmwanprj,wanprj,nkpt)
! deallocate(wanprj)
! call genwrest
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
end do
allocate(mmm(ldwanprj,ldwanprj,nspinor,nspinor,ngrf,nqpt))
allocate(umat(ldwanprj*nspinor,ldwanprj*nspinor,ldwanprj*nspinor,ldwanprj*nspinor,nwrf))
allocate(zrhoq(nstsv,nstsv,nspinor,nspinor,ngrf,nqpt))
allocate(zrho(nstsv,nstsv,nspinor,nspinor,ngrf))
mmm=zzero
umat=zzero
zrhoq=zzero
do iq=1,nqpt
  allocate(umatq(ldwanprj*nspinor,ldwanprj*nspinor,ldwanprj*nspinor,ldwanprj*nspinor,nwrf))
  !write(*,*)'writemnn',iq
  call holdthd(nkptnr,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zrhokq,mmmq) &
!$OMP NUM_THREADS(nthd)
  allocate(zrhokq(nstsv,nstsv,nspinor,nspinor,ngrf))
  allocate(mmmq(ldwanprj,ldwanprj,nspinor,nspinor,ngrf))
!$OMP DO
  do ik=1,nkptnr
!$OMP CRITICAL(writemnn_1)
    !write(*,'("Info(writecrpaumat): ",I6," of ",I6," k-points")') ik,nkptnr
!$OMP END CRITICAL(writemnn_1)
    call genzrhokq(iq,vkl(:,ik),vql(:,iq),zrhokq)
!$OMP CRITICAL(writemnn_1)
    if (twzrkq) then
      call writezrhokq(ik,iq,vkl(:,ik),vql(:,iq),zrhokq)
    end if
    !write(6,*)'run genwanrhokq'
!$OMP END CRITICAL(writemnn_1)
    call genwanrhokq(ik,iq,vkl(:,ik),vql(:,iq),zrhokq,mmmq)
   !call genwanrhokq_matop(vkl(:,ik),vql(:,iq),zrhokq,mmmq)
! !$OMP CRITICAL(writemnn_1)
!     write(*,'("Info(writecrpaumat): wanrhokq generated and wrote",I6," of ",I6," k-points")') ik,nkptnr
!     call writemmmq(ik,iq,mmmq)
! !$OMP END CRITICAL(writemnn_1)
!$OMP CRITICAL(writemnn_1)
    mmm(:,:,:,:,:,iq)=mmmq*wkptnr+mmm(:,:,:,:,:,iq)
    zrhoq(:,:,:,:,:,iq)=zrhokq*wkptnr+zrhoq(:,:,:,:,:,iq)
!$OMP END CRITICAL(writemnn_1)
  end do
!$OMP END DO
  deallocate(zrhokq,mmmq)
!$OMP END PARALLEL
  call freethd(nthd)
  ! if (np_mpi.gt.1) then
  !   n=nwrf*(ldwanprj*nspinor**4)
  !   call mpi_allreduce(mpi_in_place,umat,n,mpi_double_complex,mpi_sum,mpicom, &
  !    ierror)
  ! end if
!$OMP CRITICAL(writemnn_1)
  call writemmmq(0,iq,mmm(:,:,:,:,:,iq))
  call writezrhokq(0,iq,vkl(:,1),vql(:,iq),zrhoq(:,:,:,:,:,iq))
  zrho(:,:,:,:,:)=zrhoq(:,:,:,:,:,iq)*wqpt(iq)+zrho(:,:,:,:,:)
  !write(*,'("Info(writecrpaumat): run genumatq ",I6," of ",I6," k-points")') iq,nqpt
  call genumatq(iq,vql(:,iq),mmm,umatq)
  !write(*,*) umatq(1,1,1,1,1,1,1)
  !write(*,'("Info(writecrpaumat): writeumatq ",I6," of ",I6," k-points")') iq,nqpt
  umatq=umatq
  umat=umat+umatq*wqpt(iq)
  !write(*,*) umatq(1,1,1,1,1,1,1)
  call writeumatq(iq,umatq*omega)
!$OMP END CRITICAL(writemnn_1)
  deallocate(umatq)
end do
call writeumatq(0,umat*omega)
call writezrhokq(0,0,vkl(:,1),vql(:,1),zrho)

!--------------------below is for output---------------------
allocate(uu(nwrf),u(nwrf),jj(nwrf),j(nwrf))
allocate(u2ind(ldwanprj*nspinor,ldwanprj*nspinor,nwrf))
allocate(j2ind(ldwanprj*nspinor,ldwanprj*nspinor,nwrf))
uu=zzero
jj=zzero
j2ind=zzero
u2ind=zzero
do ld1=1,ldwanprj
  do ld2=1,ldwanprj
    do ispn=1,nspinor
      ild1=ld1+(ispn-1)*nspinor
      do jspn=1,nspinor
        ild2=ld2+(jspn-1)*nspinor
        u2ind(ild1,ild2,1:nwrf)=umat(ild1,ild2,ild1,ild2,1:nwrf)
        if (ld1.ne.ld2) j2ind(ild1,ild2,1:nwrf)=umat(ild1,ild2,ild2,ild1,1:nwrf)
        uu=uu+umat(ild1,ild2,ild1,ild2,1:nwrf)
        if(ld1.ne.ld2) jj=jj+umat(ild1,ild2,ild2,ild1,1:nwrf)
      end do
    end do
  end do
end do
open(50,file='CRPAU.OUT',form='FORMATTED')
do iw=1,nwrf
  write(50,'(2G18.10,2G18.10)')wrf(iw),uu(iw)*(1.d0/dble(nspinor**2))*(1.d0/(dble(ldwanprj))**2)*omega
end do 
close(50)
open(50,file='CRPAU4IND.OUT',form='FORMATTED')
do ld1=1,ldwanprj
  do ld2=1,ldwanprj
    do ispn=1,nspinor
      ild1=ld1+(ispn-1)*nspinor
      do jspn=1,nspinor
        ild2=ld2+(jspn-1)*nspinor
        do iw=1,nwrf
          write(50,'(2I3.2,2I2.1,2G18.10,2G18.10)')ld1,ld2,ispn,jspn,&
          wrf(iw),&
          umat(ild1,ild2,ild1,ild2,iw)*omega
        end do
      end do
    end do
  end do
end do
close(50)
!
! do ispn=1,nspinor
!   do jspn=1,nspinor
!   open(50,file='CRPA4indU.OUT',form='FORMATTED')
!   do ld1=1,ldwanprj
!     do ld2=1,ldwanprj
!       write(50,'(2G18.10,2G18.10)')wrf(iw)*ha_ev,uu(iw)*(1.d0/dble(nspinor**2))*(1.d0/(dble(ldwanprj))**2)*ha_ev
!     end do
!   end do
!   end do
! end do
! close(50)
open(50,file='CRPAU2IND_W0.OUT',form='FORMATTED')
open(51,file='CRPAJ2IND_W0.OUT',form='FORMATTED')
do ld1=1,ldwanprj
  do ld2=1,ldwanprj
    do ispn=1,nspinor
      ild1=ld1+(ispn-1)*nspinor
      do jspn=1,nspinor
        ild2=ld2+(jspn-1)*nspinor
          write(50,'(2I3.2,2I2.1,2G18.10)')ld1,ld2,ispn,jspn,&
          umat(ild1,ild2,ild1,ild2,1)*omega*Ha_ev
          if (ld1.ne.ld2) then
            write(51,'(2I3.2,2I2.1,2G18.10)')ld1,ld2,ispn,jspn,&
            umat(ild1,ild2,ild2,ild1,1)*omega*Ha_ev
          else
            write(51,'(2I3.2,2I2.1,2G18.10)')ld1,ld2,ispn,jspn,&
            zzero
          end if
      end do
    end do
  end do
end do
close(50)
close(51)



u=(1.d0/dble(nspinor**2))*(1.d0/(dble(ldwanprj))**2)*uu*omega*Ha_ev
j=(1.d0/dble(nspinor**2))*(1.d0/(dble(ldwanprj)*dble(ldwanprj-1)))*jj*omega*Ha_ev
write(6,'("U in eV is :",G18.10)') real(u(1))
write(6,'("J in eV is :",G18.10)') real(j(1))
deallocate(mmm,umat)
deallocate(u,j,uu,jj,u2ind,j2ind,zrho)
return
end subroutine