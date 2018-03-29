module geom_data_m
use const_m
implicit none

type::geom_t
  integer::id=0,n_gd=2
  real(rp)::lc_min=1,lc_max=1
  real(rp)::vol_min=1,vol_max=1

  ! shape(msh%n_dx, msh%n_vrt)
  real(rp),dimension(:,:),allocatable::dx
  real(rp),dimension(:,:),allocatable::x

  ! shape(msh%n_dx, msh%n_fce)
  real(rp),dimension(:,:),allocatable::norm
  real(rp),dimension(:,:),allocatable::x_fc
  real(rp),dimension(:,:),allocatable::x_pc
  real(rp),dimension(:,:),allocatable::x_nc

  ! shape(msh%n_fce)
  real(rp),dimension(:),allocatable::area
  real(rp),dimension(:),allocatable::l_pn
  real(rp),dimension(:),allocatable::w_fp

  ! shape(msh%n_dx, msh%n_elm)
  real(rp),dimension(:,:),allocatable::x_vc

  ! shape(msh%n_elm)
  real(rp),dimension(:),allocatable::vol
  real(rp),dimension(:),allocatable::volr

  ! shape(msh%n_dx, msh%n_dx, msh%n_elm)
  real(rp),dimension(:,:,:),allocatable::dx_inv

  ! shape(msh%n_ebf)
  real(rp),dimension(:),allocatable::d_wall
end type

end module
