subroutine writewanvars
use modmain
use modwanprj
implicit none
! arguments
! local variables
integer it1,it2,ik,ist1
!real(8) dw,w,e
character(256) fname
  write(fname,'("WANPRJ_VAR1.OUT")') 
  open(50,file=trim(fname),form='FORMATTED')
  write(50,'(A20,5X,G18.10)') 'maxnstwanprj', maxnstwanprj
  write(50,'(A20,5X,G18.10)') 'nprojwanprj(nproj)', nprojwanprj
  write(50,'(A20,5X,G18.10)') 'ldwanprj', ldwanprj
  write(50,'(A20,5X,G18.10)') 'norb:number of projs',norb
  write(50,'(A20,5X,G18.10)') 'ngrf', ngrf
  write(50,'(A20,5X,G18.10)') 'nwrf', nwrf
  close(50)
  !whch---- modified
  write(fname,'("WAN_IDX.OUT")') 
  open(51,file=trim(fname),form='FORMATTED')
  write(51,'("ik,ist1,idxwanprj")') 
  do ik=1,nkpt
    do ist1=1,projstwanprj(ik)
      write(51,'(3I10.8)') ik,ist1,idxwanprj(ist1,ik)
    end do
  end do
  close(51)
  write(fname,'("WAN_PRJST.OUT")') 
  open(51,file=trim(fname),form='FORMATTED')
  write(51,'("ik,projstwanprj")') 
  do ik=1,nkpt
    write(51,'(2I6.4)') ik,projstwanprj(ik)
  end do
  close(51)

return
end subroutine