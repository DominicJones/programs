module connectivity_m
implicit none
contains



!==============================================================================!
subroutine connectivity(msh_fmt,msh)
!------------------------------------------------------------------------------!
  use mesh_data_m
  use mesh_format_m
  use list_m
  use alloc_m
  use sort_m
  use sparse_m
!------------------------------------------------------------------------------!
  type(msh_fmt_t),intent(in)::msh_fmt
  type(msh_t),intent(inout)::msh
!------------------------------------------------------------------------------!
  type(elm_fmt_t),dimension(:),pointer::elm_fmt
  type(fce_fmt_t),dimension(:),pointer::fce_fmt

  type(lst_t),dimension(:),allocatable::fce_vrt,fce_elm!,fce_vrt2
  integer,dimension(:),allocatable::fce_tag

  integer::k,k2,i,ib,j,l,m1,m2,n,n_dx,n_fce_vrt,n_elm_vrt,ty
!------------------------------------------------------------------------------!

  ! create temporary face lists
  n=msh%n_fce*2
  call alloc_lst(fce_vrt,n)
  call alloc_lst(fce_elm,n,2)
  call alloc(fce_tag,i=n)

  ! - add boundary faces
  k=0
  do i=1,msh%n_bfce
    k=k+1
    n_fce_vrt=size(msh%fce_vrt(i)%lst)
    fce_tag(k)=n_fce_vrt
    allocate(fce_vrt(k)%lst(n_fce_vrt))
    fce_vrt(k)%lst=msh%fce_vrt(i)%lst
  end do

  ! - add element faces
  elm_fmt=>msh_fmt%elm_fmt
  do i=1,msh%n_elm
    n_dx=msh%elm_tag(i)%lst(2)
    n_elm_vrt=size(msh%elm_vrt(i)%lst)

    ty=msh%elm_tag(i)%lst(3)

    fce_fmt=>elm_fmt(ty)%fce_fmt
    do j=1,size(fce_fmt)
      n_fce_vrt=fce_fmt(j)%n_vrt
      k=k+1
      fce_tag(k)=n_fce_vrt
      fce_elm(k)%lst(1)=i
      allocate(fce_vrt(k)%lst(n_fce_vrt))
      fce_vrt(k)%lst=msh%elm_vrt(i)%lst(fce_fmt(j)%vrt)
    end do
  end do


  ! find repeated faces
  call generate_graph(fce_vrt,fce_elm,fce_tag)


  ! complete the permenent face lists
  call alloc_lst(msh%fce_elm,msh%n_fce,2)

  ! - boundary faces: fce_elm list
  ib=msh%n_elm
  do k=1,msh%n_bfce
    k2=fce_tag(k); ib=ib+1
    msh%fce_elm(k)%lst(1)=fce_elm(k2)%lst(1)
    msh%fce_elm(k)%lst(2)=-ib
    msh%fce_vrt(k)%lst=fce_vrt(k2)%lst
  end do

  ! - element faces: fce_elm & fce_vrt list
  k=msh%n_bfce
  do k2=msh%n_bfce+1,size(fce_tag)
    if(fce_tag(k2)>0.and.fce_elm(k2)%lst(1)>0.and.fce_elm(k2)%lst(2)>0)then
      k=k+1
      msh%fce_elm(k)%lst=fce_elm(k2)%lst
      allocate(msh%fce_vrt(k)%lst(size(fce_vrt(k2)%lst)))
      msh%fce_vrt(k)%lst=fce_vrt(k2)%lst
    end if
  end do


  ! matrix representation
  call list_to_matrix(msh%fce_elm(msh%n_bfce+1:),msh%elm_ja,msh%elm_ia)
!   call print_csr_matrix(10,"elm_mat.dat",msh%elm_ja,msh%elm_ia)


  ! create temporary boundary egde lists
  n=msh%n_bedg*2
  call alloc_lst(fce_vrt,n)
  call alloc_lst(fce_elm,n,2)
  call alloc(fce_tag,i=n)

  ! - add egdes
  k=0
  do i=1,msh%n_bfce
    n_fce_vrt=size(msh%fce_vrt(i)%lst)
    do j=1,n_fce_vrt
      k=k+1
      fce_tag(k)=msh%n_dx-1
      fce_elm(k)%lst(1)=i

      if(msh%n_dx==3)then
        allocate(fce_vrt(k)%lst(2))
        m1=j; m2=j+1; if(j==n_fce_vrt)m2=1
        fce_vrt(k)%lst(1)=msh%fce_vrt(i)%lst(m1)
        fce_vrt(k)%lst(2)=msh%fce_vrt(i)%lst(m2)
      end if

      if(msh%n_dx==2)then
        allocate(fce_vrt(k)%lst(1))
        fce_vrt(k)%lst(1)=msh%fce_vrt(i)%lst(j)
      end if
    end do
  end do


  ! find repeated edges
  call generate_graph(fce_vrt,fce_elm,fce_tag)


  ! complete the permenent face lists
  call alloc_lst(msh%edg_fce,msh%n_bedg,2)
  call alloc_lst(msh%edg_vrt,msh%n_bedg)

  ! - edg_fce & edg_vrt list
  k=0
  do k2=1,size(fce_tag)
    if(fce_tag(k2)>0.and.fce_elm(k2)%lst(1)>0.and.fce_elm(k2)%lst(2)>0)then
      k=k+1
      msh%edg_fce(k)%lst=fce_elm(k2)%lst
      allocate(msh%edg_vrt(k)%lst(size(fce_vrt(k2)%lst)))
      msh%edg_vrt(k)%lst=fce_vrt(k2)%lst
    end if
  end do


  ! matrix representation
  call list_to_matrix(msh%edg_fce,msh%fce_ja,msh%fce_ia)
!   call print_csr_matrix(10,"fce_mat.dat",msh%fce_ja,msh%fce_ia)


  call invert_list(msh%elm_vrt,msh%vrt_elm)

  call arc_csr_indices(msh%elm_ja,msh%elm_ia,msh%fce_elm,msh%fce_tag,4)
  call node_csr_indices(msh%elm_ja,msh%elm_ia,msh%elm_tag,4)
end subroutine



!==============================================================================!
subroutine generate_graph(fce_vrt,fce_elm,fce_tag)
!------------------------------------------------------------------------------!
  use alloc_m
  use list_m
  use sort_m
!------------------------------------------------------------------------------!
  type(lst_t),dimension(:),intent(in)::fce_vrt
  type(lst_t),dimension(:),intent(inout)::fce_elm
  integer,dimension(:),intent(inout)::fce_tag
!------------------------------------------------------------------------------!
  type(lst_t),dimension(:),allocatable::vrt_fce,fce_vrt2
  integer,dimension(size(fce_vrt))::arc_exists
  integer::ifce,ivrt_loc,ivrt,n_vrt_fce,ivrt_fce,jfce_loc,jfce
  integer::n_mtch,i,iarc,iflst,flst(299),summ,m_mtch,n
!------------------------------------------------------------------------------!

  n = size(fce_vrt)
  call alloc_lst(fce_vrt2,n)

  do i = 1, size(fce_vrt2)
    n = size(fce_vrt(i)%lst)
    allocate(fce_vrt2(i)%lst(n))
    fce_vrt2(i)%lst = fce_vrt(i)%lst
!     print *, i, n, fce_vrt2(i)%lst
    call hpsorti (n, fce_vrt2(i)%lst)
  end do


  call invert_list(fce_vrt2,vrt_fce)

  m_mtch = 0
  arc_exists = 0

  do ifce = 1, size(fce_vrt2)
    ivrt_loc = shortest_list(fce_vrt2(ifce),vrt_fce)
    ivrt = fce_vrt2(ifce)%lst(ivrt_loc)
    n_vrt_fce = size(vrt_fce(ivrt)%lst)

    do ivrt_fce = 1, n_vrt_fce
      jfce = vrt_fce(ivrt)%lst(ivrt_fce)


      if (jfce == ifce) cycle
      if (fce_elm(jfce)%lst(2) > 0) cycle
      if (size(fce_vrt2(jfce)%lst) /= size(fce_vrt2(ifce)%lst)) cycle

      summ = 0
      do i = 1, size(fce_vrt2(jfce)%lst)
        summ = summ + abs(fce_vrt2(jfce)%lst(i)-fce_vrt2(ifce)%lst(i))
        if (summ > 0) exit
      end do

      if (summ == 0) then
        fce_elm(ifce)%lst(2) = fce_elm(jfce)%lst(1)
        arc_exists(ifce) = jfce
        m_mtch = m_mtch + 1
        exit
      end if

    end do
  end do

  fce_tag=arc_exists
  print *, "match: ", size(fce_vrt2), m_mtch, size(fce_vrt2) - m_mtch*2
end subroutine



!==============================================================================!
subroutine arc_csr_indices(elm_ja,elm_ia,fce_elm,fce_tag,offset)
!------------------------------------------------------------------------------!
  use list_m
  use sort_m
  use sparse_m
!------------------------------------------------------------------------------!
  integer,dimension(:),intent(in)::elm_ja, elm_ia
  type(lst_t),dimension(:),intent(in)::fce_elm
  type(lst_t),dimension(:),intent(inout)::fce_tag
  integer,intent(in)::offset
!------------------------------------------------------------------------------!
  integer::k,i,j,kii,kjj,kij,kji
!------------------------------------------------------------------------------!

  do k=1,size(fce_elm)
    i=fce_elm(k)%lst(1)
    j=fce_elm(k)%lst(2)

    kii=coo2csr(i,i,elm_ja,elm_ia)

    if(j>0)then
      kjj=coo2csr(j,j,elm_ja,elm_ia)
      kij=coo2csr(i,j,elm_ja,elm_ia)
      kji=coo2csr(j,i,elm_ja,elm_ia)
    else
      kjj=0; kij=0; kji=0
    end if

    fce_tag(k)%lst(offset+1)=kii
    fce_tag(k)%lst(offset+2)=kjj
    fce_tag(k)%lst(offset+3)=kij
    fce_tag(k)%lst(offset+4)=kji
  end do
end subroutine



!==============================================================================!
subroutine node_csr_indices(elm_ja,elm_ia,elm_tag,offset)
!------------------------------------------------------------------------------!
  use list_m
  use sort_m
  use sparse_m
!------------------------------------------------------------------------------!
  integer,dimension(:),intent(in)::elm_ja, elm_ia
  type(lst_t),dimension(:),intent(inout)::elm_tag
  integer,intent(in)::offset
!------------------------------------------------------------------------------!
  integer::i,kii
!------------------------------------------------------------------------------!

  do i=1,size(elm_tag)
    kii=coo2csr(i,i,elm_ja,elm_ia)
    elm_tag(i)%lst(offset+1)=kii
  end do
end subroutine

end module
