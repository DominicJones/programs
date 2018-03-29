module mesh_format_m
implicit none


type::fce_fmt_t
  integer::n_dx,fce_ty
  integer::n_vrt
  integer,dimension(:),pointer::vrt
end type


type::elm_fmt_t
  character(20)::name
  integer::n_dx,elm_ty
  integer::n_vrt,n_fce
  type(fce_fmt_t),dimension(:),pointer::fce_fmt
end type


type::msh_fmt_t
  character(50)::path,fnm,ty
  integer::n_ty,is
  type(elm_fmt_t),dimension(:),pointer::elm_fmt
end type


! boundary physical types
integer,dimension(2),parameter::inlet=(/10,19/)
integer,dimension(2),parameter::outlt=(/20,29/)
integer,dimension(2),parameter::wall=(/30,39/)
integer,dimension(2),parameter::symm=(/40,49/)


contains



subroutine read_mesh_format(ios,msh_fmt)
  integer,intent(in)::ios
  type(msh_fmt_t),intent(inout)::msh_fmt
  type(elm_fmt_t),pointer::elm_fmt
  type(fce_fmt_t),pointer::fce_fmt
  character(50)::fnm,name
  integer,parameter::max_ty=20
  integer::ity,ifce
  integer::n_dx,elm_ty

  allocate(msh_fmt%elm_fmt(max_ty))


!   msh_fmt%fnm=trim(msh_fmt%ty)//".fmt"
  fnm=trim(msh_fmt%path)//trim(msh_fmt%fnm)
  print "(1x,a,a)","opening: ",trim(fnm)
  open(ios,file=trim(fnm),status="old")


  read(ios,*) msh_fmt%n_ty
!   allocate(msh_fmt%elm_fmt(msh_fmt%n_ty))


  read(ios,*) msh_fmt%is
  if (msh_fmt%is==0) msh_fmt%is=1


  do ity=1,msh_fmt%n_ty
!     elm_fmt=>msh_fmt%elm_fmt(ity)
!     read(ios,*) elm_fmt%name,elm_fmt%n_dx,elm_fmt%elm_ty
!     read(ios,*) elm_fmt%n_vrt,elm_fmt%n_fce

    read(ios,*) name,n_dx,elm_ty
    elm_fmt=>msh_fmt%elm_fmt(elm_ty)
    elm_fmt%name=name
    elm_fmt%n_dx=n_dx
    elm_fmt%elm_ty=elm_ty

    read(ios,*) elm_fmt%n_vrt,elm_fmt%n_fce


    allocate(elm_fmt%fce_fmt(elm_fmt%n_fce))
    do ifce=1,elm_fmt%n_fce
      fce_fmt=>elm_fmt%fce_fmt(ifce)
      read(ios,*) fce_fmt%n_vrt
      backspace(ios)
      allocate(fce_fmt%vrt(fce_fmt%n_vrt))
      read(ios,*) fce_fmt%n_vrt,fce_fmt%vrt(:)
      fce_fmt%vrt=fce_fmt%vrt+msh_fmt%is
    end do
  end do


  close(ios)
  print "(1x,a,a)","closed:  ",trim(fnm)
end subroutine



! function elm_name(n_dx,n_vrt)
!   integer,intent(in)::n_dx,n_vrt
!   character(len=20)::elm_name
!
!   if(n_dx==1)then
!     select case(n_vrt)
!     case(2)
!       elm_name="line"
!     end select
!   else if(n_dx==2)then
!     select case(n_vrt)
!     case(3)
!       elm_name="triangle"
!     case(4)
!       elm_name="quadrangle"
!     end select
!   else if(n_dx==3)then
!     select case(n_vrt)
!     case(4)
!       elm_name="tetrahedron"
!     case(5)
!       elm_name="pyramid"
!     case(6)
!       elm_name="prism"
!     case(8)
!       elm_name="hexahedron"
!     end select
!   end if
! end function



! ! should be a hash table look-up
! function elm_type(msh_fmt,elm_name)
!   type(msh_fmt_t)::msh_fmt
!   character(len=*)::elm_name
!   integer::elm_type
!   integer::i
!
!   elm_type=-1
!   do i=1,size(msh_fmt%elm_fmt)
!     if(msh_fmt%elm_fmt(i)%name(1:2)==elm_name(1:2))then
!       elm_type=msh_fmt%elm_fmt(i)%elm_ty
!       exit
!     end if
!   end do
! end function

end module
