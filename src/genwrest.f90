subroutine genwrest
use modmain
use modwanprj
use modmpi
!use modvars
implicit none
! local variables
integer ik,iw
! initialise global variables
call init0
call init1
call init2
call init3

! pmpd=.true.
!whch write some info
call writeivkinfo
call writeiqinfo
call writekmesh
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
if (potmode.ne.0) then
  call genpmat
  ! generate the inverse dielectric function and write to file
  call genepsinvrest
end if
!write wrf
if (mp_mpi) then
  open(50,file=trim('WRF.OUT'),form='FORMATTED')
  do iw=1,nwrf
    write(50,'(2G18.10)') wrf(iw)
  end do
  close(50)
end if
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(writeepsinvrest):")')
  write(*,'(" inverse RPA dielectric function, eps^(-1)_r(G,G'',q,w), written to &
    &EPSINVR_Q....OUT")')
end if
!
return
end subroutine