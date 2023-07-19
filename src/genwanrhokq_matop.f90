subroutine genwanrhokq_matop(vkpl,vqpl,zrhokq,mmm)
! calculate the M_{G,m,m'}(k,q)=
!             Σ_{ist1=1,n,ist2=1,n'}
!              {conjg(wanprj_kq(ist1,m))*zrhokq_G(ist1,ist2)
!             *wanprj_k(kst2,m')}
use modmain
!use modgw
use modomp
use modwanprj
implicit none
! arguments
!integer, intent(in) :: ikp
real(8),intent(in) :: vkpl(3),vqpl(3)
complex(8), intent(in) :: zrhokq(nstsv,nstsv,nspinor,nspinor,ngrf)
complex(8),intent(out) :: mmm(ldwanprj,ldwanprj,nspinor,nspinor,ngrf)
! local variables
integer ik,jk,iq,ist1,ist2
integer ikq,jkq,ispn,jspn
integer ld1,ld2,ia,iorb,is,l,lmmax
integer i,j,ias,ig,istk,istkq,irk
integer nkst,nkqst
integer iii(3)
integer idxk(nstsv),idxkq(nstsv)
integer isym
integer nthd
real(8) vkql(3),vc(3),t1,t2
real(8) coeff
complex(8) z1,z2,z
character(256) fnamewan,fpre
! automatic arrays
integer idx(nstsv)
! allocatable arrays
!integer(8), allocatable :: lock(:)
complex(8), allocatable :: zrho(:,:,:,:)
complex(8), allocatable :: wanprjkq(:,:,:,:,:),wanprjk(:,:,:,:,:)

!complex(8), allocatable :: gs(:,:),wc(:,:),zv(:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(wanprjk(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj))
allocate(wanprjkq(ldwanprj,maxnstwanprj,nspinor,norb,nprojwanprj))
!find out q-point index 
fnamewan=trim('WANPRJ_ELK')//trim('.OUT')
call findqpt(vqpl,isym,iq)
!generate q related quantities
! k+q-vector in lattice coordinates
call findkpt(vkpl,isym,jk)
!find the index of non-reduced kset
iii=vkpl*ngridk
irk=ivkiknr(iii(1),iii(2),iii(3))
!
vkql(:)=vkpl(:)+vqpl(:)
call findkpt(vkql,isym,jkq)
! equivalent reduced k-points for k and k+q
ik=ivkik(ivk(1,jk),ivk(2,jk),ivk(3,jk))
!map q point to reduced-kpoint set.
ikq=ivkik(ivk(1,jkq),ivk(2,jkq),ivk(3,jkq))
! get wanproj for k and k+q(from equivalent reduced k and k+q)
nkst=projstwanprj(ik)
nkqst=projstwanprj(ikq)
idxk(1:nkst)=idxwanprj(1:nkst,ik)
idxkq(1:nkqst)=idxwanprj(1:nkqst,ikq)
!generate the wanproj on the fly
! notice that the rotation of evecsv(k+Q) has been applied when
! reading the evecsv(see getevecsv and wanprojkpl)
call wanprojkpl(vkpl,idxk(1:nkst),ldwanprj,nprojwanprj,maxnstwanprj,nkst,&
  subulmwanprj,sublmwanprj,symlmwanprj,wanprjk)
call wanprojkpl(vkql,idxkq(1:nkqst),ldwanprj,nprojwanprj,maxnstwanprj,nkqst,&
  subulmwanprj,sublmwanprj,symlmwanprj,wanprjkq)
! write(fpre,'("_IQ",I6.6)')iq
! call writewanprjk_format(irk,vkql,fpre,wanprjkq)
! write(fpre,'("_IQ",I6.6)')0
! call writewanprjk_format(irk,vkpl,fpre,wanprjk)
!!!!!!!!!
! rotation matrix only apply to the evecf/sv
! not sure why the original wanprojk not useing getevecf/sv with vector vkl,
! but using the index of vkl.(getevecsv(ikp=0))
! not sure if k+Q vec can be always mapped to the reduced-set withthe  
! wavefunc keep same.
! the eigen value would keep same so the band indx for eigen value would
! be same as the reduced k set

! if generate the wanproj on the fly, use below code
! nstk=projstwanprj(ik)
! nstkq=projstwanprj(ikq)
!

! call wanprojkpl(vkpl,idxk,ldwanprj,nprojwanprj,maxnstwanprj,nkst,&
!   subulmwanprj,sublmwanprj,symlmwanprj,wanprjk)
! fpre='K'
! call writewanprjk_format(ik,fpre,wanprjk)
! call wanprojkpl(vkql,idxkq,ldwanprj,nprojwanprj,maxnstwanprj,nkqst,&
!   subulmwanprj,sublmwanprj,symlmwanprj,wanprjkq)
! fpre='KQ'
! call writewanprjk_format(ikq,fpre,wanprjkq)

! loop over ig,jg
! calculate the matrix of  <ψ_{k+q,n,σ}|exp(-i·(q+G)·r))|ψ_{k,n',σ'}> * conjg(Cm,n,kq)* Cm',n',k
! try matrix format.
mmm=0.d0
call holdthd(ngrf,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(iorb,is,l,lmmax,i,j,ia,ias,ispn,jspn) &
!$OMP PRIVATE(ld1,ld2,ist1,ist2,istkq,istk) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do ig=1,ngrf
  do iorb=1,norb
    is=orb(iorb,1)
    l=orb(iorb,2)
    lmmax=2*l+1
    i=l**2+1
    j=(l+1)**2
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ispn=1,nspinor
        do jspn=1,nspinor
          do ld1=1,lmmax
            do ld2=1,lmmax
              do ist1=1,nkqst
                istkq=idxkq(ist1)
                do ist2=1,nkst
                  istk=idxk(ist2)
                  mmm(ld1,ld2,ispn,jspn,ig)=mmm(ld1,ld2,ispn,jspn,ig)+&
                  zrhokq(istkq,istk,ispn,jspn,ig)*&
                  conjg(wanprjkq(ld1,ist1,ispn,iorb,ia))*&
                  wanprjk(ld2,ist2,jspn,iorb,ia)
          !conjg(wanprjkq(1:lmmax,1:nkqst))*zrho(1:nkqst,1:nkst)*wanprjkq(1:lmmax,1:nkqst,jspn,iorb,ia)
          ! call zgemm('N','N',lmmax,nkst,nkst,zone, &
          !   zrhokq(:,idxk,ispn,jspn,ig),lmmax,z1,nkst,zzero,z,lmmax)
          ! call zgemm('C','N',lmmax,lmmax,mst,zone, &
          !   z,lmmax,z2,lmmax,zzero,z3,lmmax)
                end do!over ist2/nkst
              end do!over ist1/nkqst
            end do!over ld2/ldk
          end do!over ld1/ldkq
        end do!over jspn/spink
      end do!over ispn/spinkq
    end do!over ia  
  end do!over iorb
end do!over ig
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
!call freethd(nthd)
!deallocate(zfgq,zrhoir,zrhomt)
deallocate(wanprjk,wanprjkq)
return


end subroutine