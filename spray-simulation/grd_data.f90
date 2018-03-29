module grd_data_m
  use precision_m
  use cfd_data_m
  implicit none


  ! constants
  integer::n_dx=1,n_tl=1
  real(wp)::l_sc=1.
  logical::cyl=.false.


  ! verticies
  type::vrt_t
    integer::n_fce=0
    real(wp),dimension(:),allocatable::x
    integer,dimension(:),allocatable::fce_lst
  end type


  ! faces
  type::fce_t
    integer::ip=0,in=0,ib=0
    integer::n_vrt=0,bnd_v=0,bnd_ty=0
    real(wp)::s=0.,l_pn=0.,w_fp=0.
    real(wp),dimension(:),allocatable::x
    real(wp),dimension(:),allocatable::n
    real(wp),dimension(:),allocatable::dx_p
    real(wp),dimension(:),allocatable::dx_n
    integer,dimension(:),allocatable::vrt_lst

    real(wp)::spd=0.,frc_dns=0.
    type(phs_t),dimension(:),pointer::phs=>null() !! req. cfd_data_m
  end type


  ! volumes
  type::vol_t
    integer::v_ty=0,n_vrt=0,n_fce=0,n_bnd=0
    real(wp)::v=0.,v_r=0.,s_f=0.
    logical::bnd_v=.false.
    real(wp),dimension(:),allocatable::x
    integer,dimension(:),allocatable::vrt_lst
    integer,dimension(:),allocatable::bnd_lst
    integer,dimension(:),pointer::fce_lst=>null()
    integer,dimension(:),pointer::ng_lst=>null()

    real(wp)::frc_vol=0.,co_num=0.,cs=0.
    real(wp),dimension(:),pointer::p=>null()
!     real(wp)::p=0.
    real(wp),dimension(:),allocatable::dpdx
    type(phs_t),dimension(:),pointer::phs=>null() !! req. cfd_data_m
  end type


  ! grids
  type::grd_t
    integer::n_vrt=0,n_fce=0,n_vol=0,n_ph=0
    integer::n_bf=0,n_vbf=0,n_bs=1,n_vc=0,n_grp=0
    type(vrt_t),dimension(:),allocatable::vrt
    type(fce_t),dimension(:),allocatable::fce
    type(vol_t),dimension(:),allocatable::vol

    type(iph_t),dimension(:),pointer::iph=>null() !! req. cfd_data_m
  end type

  integer::n_grd=1
  type(grd_t),dimension(:),pointer::grd=>null()


  ! volume types
  type::fce_lst_t
    integer::n_vrt=0
    integer,dimension(:),allocatable::vrt_lst
  end type

  type::vol_ty_t
    integer::n_dx=0,n_vrt=0,n_fce=0,v_ty=0
    integer,dimension(:),allocatable::vrt_lst
    type(fce_lst_t),dimension(:),allocatable::fce_lst
  end type

  integer::n_vol_ty=0
  type(vol_ty_t),dimension(:),pointer::vol_ty=>null()
end module
