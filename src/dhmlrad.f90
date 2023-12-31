
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dhmlrad
use modmain
use modphonon
implicit none
! local variables
integer is,ias
integer nr,nri,iro
integer ir,npi,i
integer l1,l2,l3,m2,lm2
integer io,jo,ilo,jlo
real(8) t1
complex(8) zsm
! begin loops over atoms and species
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  npi=npmti(is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
  do l1=0,lmaxapw
    do io=1,apword(l1,is)
      do l3=0,lmaxapw
        do jo=1,apword(l3,is)
          lm2=0
          do l2=0,lmaxi
            do m2=-l2,l2
              lm2=lm2+1
              zsm=0.d0
              i=lm2
              do ir=1,nri
                t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*wrmt(ir,is)
                zsm=zsm+t1*dvsmt(i,ias)
                i=i+lmmaxi
              end do
              do ir=iro,nr
                t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*wrmt(ir,is)
                zsm=zsm+t1*dvsmt(i,ias)
                i=i+lmmaxo
              end do
              dhaa(lm2,jo,l3,io,l1,ias)=zsm
            end do
          end do
          do l2=lmaxi+1,lmaxo
            do m2=-l2,l2
              lm2=lm2+1
              zsm=0.d0
              i=npi+lm2
              do ir=iro,nr
                t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*wrmt(ir,is)
                zsm=zsm+t1*dvsmt(i,ias)
                i=i+lmmaxo
              end do
              dhaa(lm2,jo,l3,io,l1,ias)=zsm
            end do
          end do
        end do
      end do
    end do
  end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do l3=0,lmaxapw
      do io=1,apword(l3,is)
        lm2=0
        do l2=0,lmaxi
          do m2=-l2,l2
            lm2=lm2+1
            zsm=0.d0
            i=lm2
            do ir=1,nri
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*wrmt(ir,is)
              zsm=zsm+t1*dvsmt(i,ias)
              i=i+lmmaxi
            end do
            do ir=iro,nr
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*wrmt(ir,is)
              zsm=zsm+t1*dvsmt(i,ias)
              i=i+lmmaxo
            end do
            dhloa(lm2,io,l3,ilo,ias)=zsm
          end do
        end do
        do l2=lmaxi+1,lmaxo
          do m2=-l2,l2
            lm2=lm2+1
            zsm=0.d0
            i=npi+lm2
            do ir=iro,nr
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*wrmt(ir,is)
              zsm=zsm+t1*dvsmt(i,ias)
              i=i+lmmaxo
            end do
            dhloa(lm2,io,l3,ilo,ias)=zsm
          end do
        end do
      end do
    end do
  end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do jlo=1,nlorb(is)
      l3=lorbl(jlo,is)
      lm2=0
      do l2=0,lmaxi
        do m2=-l2,l2
          lm2=lm2+1
          zsm=0.d0
          i=lm2
          do ir=1,nri
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*wrmt(ir,is)
            zsm=zsm+t1*dvsmt(i,ias)
            i=i+lmmaxi
          end do
          do ir=iro,nr
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*wrmt(ir,is)
            zsm=zsm+t1*dvsmt(i,ias)
            i=i+lmmaxo
          end do
          dhlolo(lm2,jlo,ilo,ias)=zsm
        end do
      end do
      do l2=lmaxi+1,lmaxo
        do m2=-l2,l2
          lm2=lm2+1
          zsm=0.d0
          i=npi+lm2
          do ir=iro,nr
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*wrmt(ir,is)
            zsm=zsm+t1*dvsmt(i,ias)
            i=i+lmmaxo
          end do
          dhlolo(lm2,jlo,ilo,ias)=zsm
        end do
      end do
    end do
  end do
! end loops over atoms and species
end do
end subroutine

