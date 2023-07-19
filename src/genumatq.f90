subroutine genumatq(vqpl,mmmq,umatq)
!calculate the cRPA 4 index umat(lm,lm,lm,lm,q,k) for iq∈nqpt & ik∈nkptnr
!For input: 
!vkpl: k vector in lat coord
!vqpl: qvec in lat coord
!mmm:  M_G^{m,m'}=conjg(wanprj_kq)*zrho_G*wanprj_k
!ukq: 
use modmain
!use modgw
use modomp
use modwanprj
implicit none
! arguments
!integer, intent(in) :: ikp
real(8),intent(in) :: vqpl(3)
complex(8), intent(in) :: mmmq(ldwanprj,ldwanprj,nspinor,nspinor,ngrf)
complex(8),intent(out) :: umatq(ldwanprj,ldwanprj,ldwanprj,ldwanprj,nspinor,nspinor,nspinor,nspinor,nwrf)
!local vars
integer ig,jg,ld1,ld2,ld3,ld4
integer ispn1,ispn2,ispn3,ispn4
integer iq,isym
integer nthd
complex(8) z1,z2,z(nwrf)
character(256)fpre
!allocatable
complex(8),allocatable :: epsi(:,:,:)
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqr(:,:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
!
allocate(vgqc(3,ngrf),gqc(ngrf),gclgq(ngrf))
allocate(jlgqr(njcmax,nspecies,ngrf))
allocate(ylmgq(lmmaxo,ngrf),sfacgq(ngrf,natmtot))
allocate(epsi(ngrf,ngrf,nwrf))
!get epsr
if (potmode.eq.0) then !use bare col
  epsi=complex(0.d0,0.d0)
  do ig=1,ngrf
      epsi(ig,ig,:)=complex(1.d0,0.d0)
  end do
elseif (potmode.eq.1) then ! use epsir/ U
  call getcfgq('EPSINV_R.OUT',vqpl,ngrf,nwrf,epsi)
elseif (potmode.eq.2) then ! use full epsiv./W
  call getcfgq('EPSINV.OUT',vqpl,ngrf,nwrf,epsi)
end if


!find q indx
call findqpt(vqpl,isym,iq)
!fpre='EPSINV_R_READBY_GETCFGQ'
!call writeepsinvrest(iq,fpre,epsi)
!
call gengqf(ngrf,vqpl,vgqc,gqc,jlgqr,ylmgq,sfacgq)

call gengclgq(.true.,iq,ngrf,gqc,gclgq)
call writegclgq(iq,gclgq)
!
!write(*,'("Info(genumatq): after generate dep quant")') 

umatq(:,:,:,:,:,:,:,:,:)=complex(0.d0,0.d0)
! call holdthd(ngrf,nthd)
! !$OMP PARALLEL DEFAULT(SHARED) &
! !$OMP PRIVATE(jg,ld1,ld2,ld3,ld4,ispn,jspn) &
! !$OMP PRIVATE(z1,z2) &
! !$OMP NUM_THREADS(nthd)
! !$OMP DO


do ig=1,ngrf
! !$OMP CRITICAL(writemnn_1)
!   write(*,'("Info(genumatq): in ig loop of ig:", I6.6)')ig 
! !$OMP END CRITICAL(writemnn_1)
  z1=sqrt(gclgq(ig))
  do jg=1,ngrf
    z2=sqrt(gclgq(jg))*z1
    do ld1=1,ldwanprj
      do ld2=1,ldwanprj
        do ld3=1,ldwanprj
          do ld4=1,ldwanprj
            do ispn1=1,nspinor
              do ispn2=1,nspinor
                do ispn3=1,nspinor
                  do ispn4=1,nspinor
                    !Σ conjg(M _m1m3,σσ',G (q)) * (M _m2m4,σσ',G' (q)) * 1/|q+G||q+G'| *
                    ! epsiv_r(G,G',q,w=1:nw) = u_m1m2m3m4,σσ'(q,w=1:nw)
                    z=&
                    conjg(mmmq(ld3,ld1,ispn3,ispn1,ig))*&
                    mmmq(ld4,ld2,ispn2,ispn4,jg)*&
                    epsi(ig,jg,1:nwrf)*z2
                    !
                    umatq(ld1,ld2,ld3,ld4,ispn1,ispn2,ispn3,ispn4,1:nwrf)=&
                    z+umatq(ld1,ld2,ld3,ld4,ispn1,ispn2,ispn3,ispn4,1:nwrf)
                  end do
                end do
              end do
            end do!over ispn
          end do!overld4
        end do!over ld3
      end do!over ld2
    end do!overld1
  end do!over jg
end do!over ig
! !$OMP END DO
! !$OMP END PARALLEL
! call freethd(nthd)
!
!write(*,*)'genumatq:after umatq1111', umatq(1,1,1,1,1,1,1,1,1)
deallocate(epsi,ylmgq,sfacgq,jlgqr,vgqc,gqc,gclgq)
return
end subroutine
