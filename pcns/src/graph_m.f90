module graph_m
use constants_m

implicit none

type graph_t
  integer::mblk=0 ! number of blocks
  integer::mrow=0 ! number of rows

  integer,dimension(:),pointer::nblk=>null() ! number of block entries
  integer,dimension(:),pointer::blk=>null() ! block start index

  integer,dimension(:),pointer::nrow=>null() ! number of row entries
  integer,dimension(:),pointer::row=>null() ! row start index

  integer,dimension(:),pointer::col=>null()

  integer,dimension(:),pointer::ival=>null()
  real(rp),dimension(:),pointer::rval=>null()
end type
contains


! C
subroutine graph_to_dual ( &
  m_grph, grph_ia, grph_ja, &
  m_dual, dual_ia, dual_ja)

  integer,intent(in)::m_grph
  integer,intent(out)::m_dual
  integer,dimension(:),intent(in)::grph_ia, grph_ja
  integer,dimension(:),pointer::dual_ia, dual_ja

  integer::i, jj1, jjn, jj, j, jjd;
  integer::nnz_grph, nnz_dual;
  integer,dimension(:),allocatable::cnt;


  nnz_grph = grph_ia(m_grph + 1) - 1;


  m_dual = 0;
  do i = 1, m_grph;
    jj1 = grph_ia(i);
    jjn = grph_ia(i + 1) - 1;
    do jj = jj1, jjn;
      j = grph_ja(jj);
      m_dual = max (j, m_dual);
    end do
  end do


  allocate (cnt(m_dual));
  cnt = 0;


  do i = 1, m_grph;
    jj1 = grph_ia(i);
    jjn = grph_ia(i + 1) - 1;
    do jj = jj1, jjn;
      j = grph_ja(jj);
      if (j > 0) cnt(j) = cnt(j) + 1;
    end do
  end do


  allocate (dual_ia(m_dual+1));
  allocate (dual_ja(nnz_grph));


  dual_ia(1) = 1;
  do i = 1, m_dual;
    dual_ia(i + 1) = dual_ia(i) + cnt(i);
    cnt(i) = 0;
  end do


  do i = 1, m_grph;
    jj1 = grph_ia(i);
    jjn = grph_ia(i + 1) - 1;
    do jj = jj1, jjn;
      j = grph_ja(jj);
      if (j > 0) then
        jjd = dual_ia(j) + cnt(j);
        dual_ja(jjd) = i;
        cnt(j) = cnt(j) + 1;
      end if
    end do
  end do
end subroutine


! C
subroutine print_range (name, mrow, row)
  character(*),intent(in)::name
  integer,intent(in)::mrow
  integer,dimension(:),intent(in)::row

  integer::i,jj1,jjn,jj,j

  write (*,"(/a)",advance="yes") trim(name)//" graph:"
  do i = 1, mrow
    jj1 = row(i)
    jjn = row(i + 1) - 1
    write (*,"(a,i0,1x,i0,a,i0)",advance="no") &
         "row [", jj1, jjn, "] ",i
    write (*,"(a)",advance="yes") ""
  end do

  write (*,"(a)",advance="yes") ""
end subroutine


! C
subroutine print_graph ( &
  name, mrow, row, col, idx, val)

  character(*),intent(in)::name
  integer,intent(in)::mrow
  integer,dimension(:),intent(in)::row, col
  integer,dimension(:),intent(in),optional::idx
  real(rp),dimension(:),intent(in),optional::val

  integer::i,jj1,jjn,jj,j

  write (*,"(/a)",advance="yes") trim(name)//" graph:"
  do i = 1, mrow
    jj1 = row(i)
    jjn = row(i + 1) - 1
    if (present(idx)) then
      write (*,"(a,i0,1x,i0,a,i0,a,i0,a)",advance="no") &
        "row [", jj1, jjn, "] ",i," -> ",idx(i),", col("
    else
      write (*,"(a,i0,1x,i0,a,i0,a)",advance="no") &
        "row [", jj1, jjn, "] ",i,", col("
    end if
    do jj = jj1, jjn
      j = col(jj)
      if (present(val)) then
        write (*,"(1x,g12.4)",advance="no") val(jj)
      else
        write (*,"(1x,i0)",advance="no") j
      end if
    end do
    write (*,"(a)",advance="yes") ")"
  end do

!   write (*,"(a/)",advance="yes") "end of "//trim(name)//" graph"
  write (*,"(a)",advance="yes") ""
end subroutine




subroutine print_blockgraph ( &
  name, mblk, blk, row, col)

  character(*),intent(in)::name
  integer,intent(in)::mblk
  integer,dimension(:),intent(in)::blk,row, col

  integer::i,jj1,jjn,jj,j,kk1,kkn,kk,k

  write (*,"(/a)",advance="yes") trim(name)//" block graph:"

  do i = 1, mblk
!     ty = mesh%fmt_ty(i)
!     dms = mesh%fmt_dim(i)
    jj1 = blk(i)
    jjn = blk(i + 1) - 1
    write (*,"(/a,1x,i0)",advance="yes") "block",i
    j = 0
    do jj = jj1, jjn
      j = j + 1
      kk1 = row(jj)
      kkn = row(jj + 1) - 1
      write (*,"(a,1x,i0,a)",advance="no") "row",j,", col("
      do kk = kk1, kkn
        k = col(kk)
        write (*,"(1x,i0)",advance="no") k
      end do
      write (*,"(a)",advance="yes") ")"
    end do
  end do

!   write (*,"(a/)",advance="yes") "end of "//trim(name)//" block graph"
  write (*,"(a)",advance="yes") ""
end subroutine
end module
