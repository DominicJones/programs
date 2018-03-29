module pde_data_m
use const_m
implicit none

type::pde_t
  integer::id=0
  integer::run=1
  integer::n_cmp=1
  integer::ldwall=0,lspd=1,lcc=1,lcd=1,lsbc=1,llim=0
  integer::lcls=1,lc_sch=0
  integer::k_iter=10,n_iter=1000

  real(rp)::cbl=0.2d0,dbl=1.d0
  real(rp)::urf=0.5d0,clim=5.d0
  real(rp)::reduc=1.d-6

  integer::nls=25,lls=2
  real(rp)::rls=1.d-1

  real(rp)::ref_mag=1

  ! shape(n_cmp)
  real(rp),dimension(:),allocatable::ref

  ! shape(msh%n_ebf)
  real(rp),dimension(:),allocatable::cc

  ! shape(msh%n_ebf)
  real(rp),dimension(:),allocatable::cd

  ! shape(n_cmp, msh%n_elm)
  real(rp),dimension(:,:),allocatable::rhs

  ! shape(n_cmp, msh%n_ebf)
  real(rp),dimension(:,:),allocatable::phi

  ! shape(msh%n_dx, n_cmp, msh%n_elm)
  real(rp),dimension(:,:,:),allocatable::grad
end type


type::bnd_t
  integer::id=31
  character(50)::qnty="wall force"
  real(rp),dimension(:),allocatable::ref
  real(rp),dimension(:),allocatable::val
end type

logical,save::last_fp_iter=.false.

end module
