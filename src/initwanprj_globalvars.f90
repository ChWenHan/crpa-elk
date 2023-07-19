subroutine initwanprj_globalvars
!This routine generates the Wannier projectors of all user input orbitals
use modmain
use moddftu
use modmpi
use modomp
use modwanprj

implicit none
!local variables
integer iorb,lm,ia,is,l,ias,lmmax,rlm,m
integer ld,ik,ist,nproj,nst,ispn
integer imin,imax,iscp
integer lm1,nea,iea,ja,ka,i,isym,lspl,la
logical, allocatable :: done(:)
complex(8), allocatable :: subulm(:,:,:,:)
integer, allocatable :: sublm(:,:,:)
integer, allocatable :: idx(:,:), projst(:)
integer, allocatable :: idxiea(:)
complex(8), allocatable :: a(:,:)
complex(8), allocatable :: z(:),zz(:),z1(:,:)
complex(8), allocatable :: symlm(:,:,:,:)
if(norb.le.0) then
  write(*,*) 'No projected orbital specified. Stopping.'
  stop
endif
!ensure correct the correlated window bounds
if(emincor.gt.emaxcor) then
  write(*,*) 'Lower correlated window bound greater than upper bound. Stopping.'
  stop
endif
!initalise the universal variables 
call init0
call init1
call init2
call init3
if (allocated(idxlm)) deallocate(idxlm)
allocate(idxlm(0:lmaxapw,-lmaxapw:lmaxapw))
lm=0
do l=0,lmaxapw
  do m=-l,l
    lm=lm+1
    idxlm(l,m)=lm
  end do
end do
! read density and potentials from file
call readstate
! Fourier transform Kohn-Sham potential to G-space
call genvsig
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
!generate radial wavefunctions
call genapwfr
!generate LO radial wavefuntions
call genlofr
!calculate the radial overlap integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! generate the spin-orbit coupling radial functions
call gensocfr
!check if input indices are correct
if(wanind) then
  !round input "indices" values to nearest integer
  imin=nint(emincor)
  imax=nint(emaxcor)
  iscp=0
  nst=nstsv
  if(spcpl) iscp=1
  !if spins are decoupled in Hamiltonian - spins are in two blocks
  if((spinpol).and.(.not.spcpl)) nst=nstfv
  if((imin.lt.1).or.(imin.gt.nst).or. &
      (imax.lt.1).or.(imax.gt.nst)) then
    write(*,*) 'Wrong band indices for correlated energy window used.'
    write(*,*) 'Please changes these indices.'
    stop
  endif 
endif
!if outputting projectors into irreducible basis
if(cubic) lmirep=.true.
ld=0
nproj=0
! rotation matrix to convert projector basis (initalised to zero)
do iorb=1,norb
  is = orb(iorb,1)
  l = orb(iorb,2)
  if(is.le.0) then
    write(*,*) '(Error(initwanprj)) incorrect species for orbital. Stopping.'
    stop
  endif
  if(l.lt.0) then
    write(*,*) '(Error(initwanprj)) incorrect l for orbital. Stopping.'
    stop
  endif
  ld=max(ld,2*l+1)
  nproj=max(nproj,natoms(is))
enddo
! set global value
ldwanprj=ld
nprojwanprj=nproj
!allocate global array
if (allocated(subulmwanprj)) deallocate(subulmwanprj)
allocate(subulmwanprj(ldwanprj,ldwanprj,norb,nprojwanprj))
if (allocated(sublmwanprj)) deallocate(sublmwanprj)
allocate(sublmwanprj(ldwanprj,norb,nprojwanprj))
!finish allocate global array
allocate(subulm(ld,ld,norb,nproj))
allocate(sublm(ld,norb,nproj))
subulm(:,:,:,:)=0.d0
sublm(:,:,:)=0.d0
lmmax=0
do iorb=1,norb
  is = orb(iorb,1)
  l = orb(iorb,2)
  lmmax=2*l+1
  do ia=1,natoms(is)
    rlm=orb(iorb,3)
    ias=idxas(ia,is)
    sublm(1:rlm,iorb,ia)=rorblm(iorb,1:rlm)
    call su2lm(lmmax,l,rlm,ias,subulm(:,:,iorb,ia),sublm(:,iorb,ia))
!check that the correct lm values have been inputted
    do lm=1,rlm
      if((sublm(lm,iorb,ia).le.0).or.(sublm(lm,iorb,ia).gt.lmmax)) then
        write(*,*) '(Error(initwanprj)) incorrect lm values input. Stopping'
        write(*,*) rlm,sublm(lm,iorb,ia)
        stop
      endif
    enddo
  enddo
enddo
subulmwanprj=subulm
sublmwanprj=sublm
!allocate the band indices in the correlated window
allocate(idx(nstsv,nkpt))
idx(:,:)=0
!allocate the number of band indices in the correlated window for each ik
allocate(projst(nkpt))
projst(:)=0
!allocate global array
if (allocated(idxwanprj)) deallocate(idxwanprj)
allocate(idxwanprj(nstsv,nkpt))
if (allocated(projstwanprj)) deallocate(projstwanprj)
allocate(projstwanprj(nkpt))
do ik=1,nkpt
!read the second variational energy eigenvalues from EVALSV.OUT
  call getevalsv(filext,ik,vkl(:,ik),evalsv(:,ik))
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
!put occupations for denser k-mesh projector calc
! count and index states at k in energy window
  nst=0
!k index for reading in eigenvectors in genwfsvpwan
! for band index inputs
  if(wanind) then
    !get indices for both spins (iscp is for spin coupled systems)
    do ispn=1,nspinor-iscp
      do ist=imin,imax
        nst=nst+1
        idx(nst,ik)=ist+(ispn-1)*nstfv
      end do  
    end do   
! for energy window bounds input
  else
    do ist=1,nstsv
      if(((evalsv(ist,ik)-efermi).ge.emincor) .and. & 
          ((evalsv(ist,ik)-efermi).le.emaxcor)) then
        nst=nst+1
        idx(nst,ik)=ist
      end if
    end do
  end if
  projst(ik)=nst
enddo
!set value for gloabl array
idxwanprj=idx
projstwanprj=projst
!Output energy window
if (mp_mpi) then
  if(wanind) then
    write(*,'("Energy window boundary indices:")') 
    write(*,*) imin, imax
  else 
    write(*,'("Energy window:")') 
    write(*,*) emincor, emaxcor
  endif
  write(*,'("")') 
endif
!nst is now the maximum number of band indices
nst=maxval(projst(:))
!set global array
maxnstwanprj=nst
!
allocate(symlm(ldwanprj,ldwanprj,norb,nprojwanprj))
do iorb=1,norb
! general projector info
  allocate(done(natmtot),idxiea(natmtot))
  done(:)=.false.
  is=orb(iorb,1)
  l=orb(iorb,2)
  lmmax=2*l+1
  lm1=(l+1)**2
  allocate(a(lm1,lm1),z1(lm1,lm1))
!set up identity matrix
  a(:,:)=0.d0
  do lm=1,lm1
    a(lm,lm)=1.d0
  end do
!loop over atoms
  do ia=1,natoms(is)
    if(done(ia)) cycle
    ias=idxas(ia,is)
!make an array of equivalent atom indices
    nea=0
    do ja=1,natoms(is)
      if(.not.eqatoms(ia,ja,is)) cycle
      nea=nea+1
      idxiea(nea)=ja
    enddo
    ka=idxiea(1)
!determine the symmetry operation required to transform
!the equivalent atom ja to idxiea(1). 
    do i=1,nea
      ja=idxiea(i)
! find first symmetry matrix which transforms ia to ja.
      do isym=1,nsymcrys
! index to spatial rotation lattice symmetry
        lspl=lsplsymc(isym)
! check that the crystal symmetry is the same as the site
! symmetry
        la=ieqatom(ja,is,isym)
! want to find a proper symmetric matrix (i.e. does not 
! include inversion symmetry) which transfroms atom la
! to ka. This is for local <-> global coordinate transformation. 
        if((ka.eq.la).and.(symlatd(lspl).gt.0.d0)) then
! calculate the symmerty matrix in lm basis
          call rotzflm(symlatc(:,:,lspl),l,l,lm1,lm1,lm1,a,z1)
! Keep desired array subset 
          symlm(1:lmmax,1:lmmax,iorb,ja)=z1(1:lmmax,1:lmmax)
! found symmetry, skip for this index
          done(ja)=.true.
! Exiting the loop
          exit
        endif
      enddo
    enddo
  enddo
  deallocate(done,idxiea,a,z1)
enddo
if (allocated(symlmwanprj)) deallocate(symlmwanprj)
allocate(symlmwanprj(ldwanprj,ldwanprj,norb,nprojwanprj))
symlmwanprj=symlm
call writewanvars
call writesubusymlm
deallocate(subulm,sublm,symlm)
return
end subroutine
