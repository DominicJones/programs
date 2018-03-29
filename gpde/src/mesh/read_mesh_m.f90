module read_mesh_m
implicit none
contains



subroutine read_vertices(ios,fnm,n_vrt,key,x,dx)
  use const_m
  use mesh_data_m

  integer, intent(in)::ios
  character(*), intent(in)::fnm
  integer, intent(inout)::n_vrt
  integer,dimension(:), intent(in)::key
  real(rp),dimension(:,:), intent(inout)::x
  real(rp),dimension(:,:), intent(inout), optional::dx

  real(rp),dimension(3)::x3
  integer::i,j,n_dim

  print "(1x,a,a)","opening: ",trim(fnm)
  open(ios,file=trim(fnm))

  read(ios,*) n_vrt
  n_dim = size(x,1)

  if(present(dx))then
    do i=1,n_vrt
      read(ios,*) j, x3(:n_dim)
      dx(:,j) = x3(:n_dim) - x(:,j)
    end do
  else
    do i=1,n_vrt
      read(ios,*) j, x3(:n_dim)
      x(:,j) = x3(:n_dim)
    end do
  end if

  close(ios)
  print "(1x,a,a)","closed:   ",trim(fnm)
end subroutine



subroutine read_mesh(ios,int_fmt,msh_fmt,msh,x0)
  use mesh_format_m
  use mesh_data_m
  use sort_m
  use hash_m
  use alloc_m

  integer,intent(in)::ios
  type(msh_fmt_t),intent(in)::int_fmt,msh_fmt
  type(msh_t),intent(inout)::msh
  real(rp),dimension(:,:),allocatable::x0

  character(50)::fnm,str
  integer,dimension(3)::n_elm_dx,n_elm_fce
  integer,dimension(10)::tag
  integer::ty,i,k,m,mm,n,idx,idx_elm,idx_fce,ifce,iset,key,iv
  integer::n_vol_fce,n_elm,n_vol,n_fce,n_bfce,n_vrt,n_dx,n_bset
  integer::n_tag,n_bvrt,n_evrt,n_dx_max
  type(hash_data_t),dimension(:),allocatable::vrt_tbl,elm_tbl


  n_dx_max=0
  n_elm_dx=0
  n_elm_fce=0


  fnm=trim(msh%path)//trim(msh%fnm)
  print "(1x,a,a)","opening: ",trim(fnm)
  open(ios,file=trim(fnm),status="old")

  ! set-up:
  !   vrt_tag, x_vrt, fce_tag, fce_vrt, elm_tag, elm_vrt
  if(trim(msh_fmt%ty)=="gmsh")call read_gmsh()
  if(trim(msh_fmt%ty)=="gambit")call read_gambit()

  close(ios)
  print "(1x,a,a)","closed:  ",trim(fnm)


  ! number of boundary vertices and edges
  n_bvrt=0
  do k=1,msh%n_bfce
    n_bvrt=n_bvrt+size(msh%fce_vrt(k)%lst)
  end do
  msh%n_bvrt=n_bvrt
  msh%n_bedg=n_bvrt/2


  print "(1x,a,3i7)", "  No. dimensions:    ", msh%n_dx
  print "(1x,a,3i7)", "  No. vertices:      ", msh%n_vrt
  print "(1x,a,3i7)", "  No. boundary egdes:", msh%n_bedg
  print "(1x,a,3i7)", "  No. faces:         ", msh%n_fce
  print "(1x,a,3i7)", "  No. boundary faces:", msh%n_bfce
  print "(1x,a,3i7)", "  No. elements:      ", msh%n_elm
  print "(1x,a,3i7)", "  No. elm. + b.faces:", msh%n_ebf
contains



subroutine read_gmsh()
  do
    read(ios,*)str
    if(trim(str)=="$Elements")then
      read(ios,*)n_elm

      do i=1,n_elm
        read(ios,*)key,ty,n_tag,tag(:n_tag)

        n_dx=msh_fmt%elm_fmt(ty)%n_dx
        n_fce=msh_fmt%elm_fmt(ty)%n_fce
        n_elm_dx(n_dx)=n_elm_dx(n_dx)+1
        n_elm_fce(n_dx)=n_elm_fce(n_dx)+n_fce

        n_dx_max=max(n_dx_max,n_dx)
      end do
      exit
    end if
  end do


  rewind(ios)


  n_dx=n_dx_max
  n_vol=n_elm_dx(n_dx)
  n_vol_fce=n_elm_fce(n_dx)
  n_bfce=n_elm_dx(n_dx-1)
  n_fce=(n_vol_fce+n_bfce)/2

  msh%n_dx=n_dx
  msh%n_elm=n_vol
  msh%n_fce=n_fce
  msh%n_bfce=n_bfce
  msh%n_ebf=n_vol+n_bfce

  call alloc_lst(msh%fce_tag,n_fce,8)
  call alloc_lst(msh%fce_vrt,n_fce)

  call alloc_lst(msh%elm_tag,n_vol,5)
  call alloc_lst(msh%elm_vrt,n_vol)


  do
    read(ios,*)str
    if(trim(str)=="$Nodes")then
      read(ios,*)n_vrt

      msh%n_vrt=n_vrt
      call alloc_lst(msh%vrt_tag,n_vrt,4)
      call alloc(x0,i=n_dx,j=n_vrt)

      call initialise_hash(vrt_tbl,n_vrt*15/10)

      do idx=1,n_vrt
        read(ios,*)key

        msh%vrt_tag(idx)%lst(1)=key
        call insert_hash_value(key,idx,vrt_tbl)

        backspace(ios)
        read(ios,*)key,x0(:,idx)
      end do
      exit
    end if
  end do


  do
    read(ios,*)str
    if(trim(str)=="$Elements")then
      read(ios,*)n_elm

      idx_elm=0; idx_fce=0
      do i=1,n_elm
        read(ios,*)key,ty,n_tag,tag(:n_tag)

        n_dx=msh_fmt%elm_fmt(ty)%n_dx
        n_vrt=msh_fmt%elm_fmt(ty)%n_vrt

        if(n_dx==n_dx_max)then
          idx_elm=idx_elm+1

          msh%elm_tag(idx_elm)%lst(1)=key
          msh%elm_tag(idx_elm)%lst(2)=n_dx
          msh%elm_tag(idx_elm)%lst(3)=ty
          msh%elm_tag(idx_elm)%lst(4)=tag(1)
! print *, n_dx,key,ty,tag(1)

          allocate(msh%elm_vrt(idx_elm)%lst(n_vrt))

          backspace(ios)
          read(ios,*)key,ty,n_tag,tag(:n_tag),msh%elm_vrt(idx_elm)%lst

          do iv=1,n_vrt
            key=msh%elm_vrt(idx_elm)%lst(iv)
            call get_hash_value(key,msh%elm_vrt(idx_elm)%lst(iv),vrt_tbl)
          end do


        else if(n_dx==(n_dx_max-1))then
          idx_fce=idx_fce+1

          msh%fce_tag(idx_fce)%lst(1)=key
          msh%fce_tag(idx_fce)%lst(2)=n_dx
          msh%fce_tag(idx_fce)%lst(3)=ty
          msh%fce_tag(idx_fce)%lst(4)=tag(1)
! print *, n_dx,key,ty,tag(1)

          allocate(msh%fce_vrt(idx_fce)%lst(n_vrt))

          backspace(ios)
          read(ios,*)key,ty,n_tag,tag(:n_tag),msh%fce_vrt(idx_fce)%lst

          do iv=1,n_vrt
            key=msh%fce_vrt(idx_fce)%lst(iv)
            call get_hash_value(key,msh%fce_vrt(idx_fce)%lst(iv),vrt_tbl)
          end do
        end if
      end do
      exit
    end if
  end do

  call finalise_hash(vrt_tbl)
end subroutine



subroutine read_gambit()
!   do
!     read(ios,*)str
!     if(trim(str)=="NUMNP")then
!       read(ios,*)n_vrt,n_elm,n,n_bset
!       exit
!     end if
!   end do
!
!
!   do
!     read(ios,*)str
!     if(trim(str)=="ELEMENTS")then
!       do i=1,n_elm
!         read(ios,*)key,ty,n
!
!         n_dx=msh_fmt%elm_fmt(ty)%n_dx
!         n_fce=msh_fmt%elm_fmt(ty)%n_fce
!         n_elm_dx(n_dx)=n_elm_dx(n_dx)+1
!         n_elm_fce(n_dx)=n_elm_fce(n_dx)+n_fce
!         msh%n_dx=max(msh%n_dx,n_dx)
!       end do
!       exit
!     end if
!   end do
!
!
!   do iset=1,n_bset
!   do
!     read(ios,*)str
!     if(trim(str)=="BOUNDARY")then
!       read(ios,*)tag(1),n,n_bfce
!
!       do i=1,n_bfce
!         read(ios,*)key,ty,ifce
!
!         n_dx=msh_fmt%elm_fmt(ty)%n_dx-1
!         n_elm_dx(n_dx)=n_elm_dx(n_dx)+1
!         n_elm_fce(n_dx)=n_elm_fce(n_dx)+1
!         msh%n_dx=max(msh%n_dx,n_dx)
!       end do
!       exit
!     end if
!   end do
!   end do
!
!
!   rewind(ios)
!
!
!   n_dx=msh%n_dx
!   n_vol=n_elm_dx(msh%n_dx)
!   n_vol_fce=n_elm_fce(msh%n_dx)
!   msh%n_bfce=n_elm_dx(msh%n_dx-1)
!   n_fce=(n_vol_fce+msh%n_bfce)/2
!   msh%n_elm=n_vol
!   msh%n_fce=n_fce
!   msh%n_ebf=n_vol+msh%n_bfce
!
!   allocate(msh%vrt(n_vrt))
!   allocate(msh%fce(n_fce))
!   allocate(msh%elm(n_vol))
!
!   call initialise_hash(vrt_tbl,n_vrt*15/10)
!   call initialise_hash(elm_tbl,n_vol*15/10)
!
!
!   call alloc(x0,i=n_dx,j=n_vrt)
!   do
!     read(ios,*)str
!     if(trim(str)=="NODAL")then
!       do idx=1,n_vrt
!         read(ios,*)key
!         vrt=>msh%vrt(idx)
!         vrt%key=key
!         call insert_hash_value(key,idx,vrt_tbl)
!
!         backspace(ios)
!         read(ios,*)key,x0(:,idx)
!       end do
!       exit
!     end if
!   end do
!
!
!   do
!     read(ios,*)str
!     if(trim(str)=="ELEMENTS")then
!
!       idx_elm=0
!       do i=1,size(msh%elm)
!         read(ios,*)key,ty,n_vrt
!
!         n_dx=msh_fmt%elm_fmt(ty)%n_dx
!         n_vrt=msh_fmt%elm_fmt(ty)%n_vrt
!
!         idx_elm=idx_elm+1
!         elm=>msh%elm(idx_elm)
!         elm%key=key
!         elm%ty=ty
!         elm%n_dx=n_dx
!         call insert_hash_value(key,idx_elm,elm_tbl)
!
!         allocate(elm%vrt(n_vrt))
!         backspace(ios)
!         read(ios,*)key,ty,n_vrt,elm%vrt(:)
!         do iv=1,n_vrt
!           key=elm%vrt(iv)
!           call get_hash_value(key,elm%vrt(iv),vrt_tbl)
!         end do
!       end do
!       exit
!     end if
!   end do
!
!
!   idx_fce=0
!   do iset=1,n_bset
!   do
!     read(ios,*)str
!     if(trim(str)=="BOUNDARY")then
!       read(ios,*)tag(1),n,n_bfce
!
!       do i=1,n_bfce
!         read(ios,*)key,ty,ifce
!
!         n_dx=msh_fmt%elm_fmt(ty)%n_dx-1
!         n_vrt=msh_fmt%elm_fmt(ty)%fce_fmt(ifce)%n_vrt
!
!         idx_fce=idx_fce+1
!         fce=>msh%fce(idx_fce)
!         fce%key=n_vol+idx_fce
!         fce%ty=elm_type(msh_fmt,elm_name(n_dx,n_vrt))
!         fce%ph_ty=tag(1)
!         fce%n_dx=n_dx
!
!         allocate(fce%vrt(n_vrt))
!
!         call get_hash_value(key,idx_elm,elm_tbl)
!         fce%vrt=msh%elm(idx_elm)%vrt( &
!           msh_fmt%elm_fmt(ty)%fce_fmt(ifce)%vrt )
!
!         do iv=1,n_vrt
!           key=fce%vrt(iv)
!           call get_hash_value(key,fce%vrt(iv),vrt_tbl)
!         end do
!       end do
!       exit
!     end if
!   end do
!   end do
!
!   call finalise_hash(vrt_tbl)
!   call finalise_hash(elm_tbl)
end subroutine
end subroutine

end module
