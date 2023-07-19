subroutine writeumatq(iq,umatq)
use modmain
use modwanprj
implicit none
!
integer(4), intent(in) :: iq
complex(8), intent(in) :: umatq(ldwanprj,ldwanprj,ldwanprj,ldwanprj,nspinor,nspinor,nspinor,nspinor,nwrf)
!
! local variables
integer ld1,ld2,ld3,ld4,ispn1,ispn2,ispn3,ispn4,iw
character(256) filename
!$OMP CRITICAL(u180)
write(filename,'("UMAT_Q",I6.6,".OUT")') iq
open(50,file=trim(filename),form='FORMATTED')
!write(50,'("ig & ngrf:",I6.3," nqpt :",I6," iq : ",I6.4," wrf, epsir")') ngrf,nqpt, iq
do ld1=1,ldwanprj
  do ld2=1,ldwanprj
    do ld3=1,ldwanprj
      do ld4=1,ldwanprj
        do ispn1=1,nspinor
          do ispn2=1,nspinor
            do ispn3=1,nspinor
              do ispn4=1,nspinor
                do iw=1,nwrf
                  write(50,'(4I3.2,4I2.1,I6.5,2G18.10)') ld1,ld2,ld3,ld4,ispn1,&
                  ispn2,ispn3,ispn4,iw,&
                  umatq(ld1,ld2,ld3,ld4,ispn1,ispn2,ispn3,ispn4,iw)
                end do
              end do
            end do 
          end do
        end do
      end do
    end do
  end do
end do
close(50)
!$OMP END CRITICAL(u180)
return
end subroutine