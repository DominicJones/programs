module sparse_m
use const_m
implicit none
contains



function coo2csr(i,j,ja,ia) result(kij)
  integer,intent(in)::i,j
  integer,dimension(:),intent(in)::ja,ia
  integer::kij

  do kij=ia(i),ia(i+1)-1
    if(ja(kij)==j)return
  end do
  kij=-1
end function



function csr_diag(a_ij,ja,ia) result(diag)
  integer,dimension(:)::ja,ia
  real(rp),dimension(:)::a_ij

  real(rp),dimension(size(ia)-1)::diag
  integer::i,kii,n
  n=size(ia)-1

  do i=1,n
    kii=coo2csr(i,i,ja,ia)
    diag(i)=a_ij(kii)
  end do
end function



subroutine csr_transpose(a_ij,ja,ia,at_ij)
  integer,dimension(:)::ja,ia
  real(rp),dimension(:)::a_ij,at_ij
  integer::i,kij,j,kji,n
  n=size(ia)-1

  do i=1,n
    do kij=ia(i),ia(i+1)-1
      j=ja(kij)
      kji=coo2csr(j,i,ja,ia)
      at_ij(kji)=a_ij(kij)
    end do
  end do
end subroutine



subroutine csr_to_csc(ja,ia,jao,iao,a,ao)
  integer,dimension(:)::ia,iao,ja,jao
  real(rp),dimension(:),optional::a,ao

  integer i,j,k,n,next,ipos

  n=size(ia)-1
  iao = 0

  do i=1, n
    do k=ia(i), ia(i+1)-1
      j = ja(k)+1
      iao(j) = iao(j)+1
    end do
  end do

  iao(1) = 1
  do i=1,n
    iao(i+1) = iao(i) + iao(i+1)
  end do

  do i=1,n
    do k=ia(i),ia(i+1)-1
      j = ja(k)
      next = iao(j)
      if(present(a).and.present(ao))then
        ao(next) = a(k)
      end if
      jao(next) = i
      iao(j) = next+1
    end do
  end do

  do i=n,1,-1
    iao(i+1) = iao(i)
  end do
  iao(1) = 1
end subroutine



subroutine matvec_prod(a_ij,phi,ja,ia,res)
  integer,dimension(:)::ja,ia
  real(rp),dimension(:)::a_ij
  real(rp),dimension(:)::phi

  real(rp),dimension(:),intent(inout)::res
  integer::i,k,n
  n=size(ia)-1

  do i=1,n
    res(i)=0.d0
    do k=ia(i),ia(i+1)-1
      res(i)=res(i)+a_ij(k)*phi(ja(k))
    end do
  end do
end subroutine



subroutine print_csr_matrix(ios,fnm,ja,ia,a_ij)
  integer,intent(in)::ios
  character(*)::fnm
  integer,dimension(:)::ja,ia
  real(rp),dimension(:),optional::a_ij

  integer::i,j,kij
  integer,dimension(size(ia))::ib
  integer,dimension(size(ja))::jb
  real(rp),dimension(size(ja))::b_ij

  if(present(a_ij))then
    call csr_to_csc(ja,ia,jb,ib,a_ij,b_ij)
  else
    b_ij=1
    call csr_to_csc(ja,ia,jb,ib)
  end if

  print "(1x,a,a)","creating: ",trim(fnm)
  open(ios,file=trim(fnm))

  write(ios,*) "# name: A"
  write(ios,*) "# type: sparse matrix"
  write(ios,*) "# nnz:",size(jb)
  write(ios,*) "# rows:",size(ib)-1
  write(ios,*) "# columns:",size(ib)-1

  do i=1,size(ib)-1
    do kij=ib(i),ib(i+1)-1
      j=jb(kij)
      write(ios,*)j,i,b_ij(kij)
    end do
  end do
  close(ios)
  print "(1x,a,a)","closed:   ",trim(fnm)
  print "(1x,a)","$> octave; load "//trim(fnm)//"; spy(A)"
end subroutine

end module
