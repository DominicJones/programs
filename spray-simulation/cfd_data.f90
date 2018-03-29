module cfd_data_m
  use precision_m
  implicit none


  ! constants
  real(wp)::m_sc=0.,t_sc=0.
  real(wp)::c_gas=8.314472
  real(wp)::c_grv=-9.81
  real(wp)::c_sph=4./3.*pi

  logical::dbg=.false.
  logical::turb=.false.,spr=.false.
  real(wp)::dt=0.,time=0.
  real(wp)::c_in=0.,c_out=0.
  real(wp)::mx_res_n=0.
  integer::n_grad=2


  ! pproperty ID data
  type::ipr_t
    logical::s_tl=.false.
    integer::n_tl=1
    real(wp),dimension(:),pointer::rho=>null()
  end type


  ! equation ID data
  type::ieq_t
    integer::n_sc=1,n_vc=1,n_tl=1
    integer::n_st=1,n_cp=1
    integer::c_sch=0,n_iter=2
    real(wp)::lb_cut=1.e-8
    real(wp)::c_bl=0.,t_bl=1.
    real(wp)::urf=0.7,reduc=0.1
    logical::solv=.true.
    logical::rhs=.false.,ddx=.false.
    integer,dimension(:),pointer::st_id=>null()
    real(wp),dimension(:),pointer::c_in=>null()
    real(wp),dimension(:),pointer::c_out=>null()
    real(wp),dimension(:),pointer::res_n=>null()
    real(wp),dimension(:),pointer::mod_cst=>null()
    real(wp),dimension(:),pointer::pr_num=>null()
    real(wp),dimension(:),pointer::phi_max=>null()
  end type


  ! phase ID data
  type::iph_t
    integer::n_pr=0,n_eq=0
    type(ipr_t),dimension(:),pointer::ipr=>null()
    type(ieq_t),dimension(:),pointer::ieq=>null()
  end type


  ! property data
  type::prp_t
    real(wp),dimension(:),pointer::rho=>null()
  end type


  ! equation data
  type eqn_cof_t
    real(wp),dimension(:),pointer::frc=>null()
    real(wp),dimension(:),pointer::dns=>null()
    real(wp),dimension(:),pointer::vis=>null()
  end type
  type::eqn_t
    type(eqn_cof_t),dimension(:),pointer::cof=>null()
    real(wp),dimension(:,:,:),pointer::phi=>null()
    real(wp),dimension(:,:,:),pointer::ddx=>null()
    real(wp),dimension(:,:),pointer::rhs=>null()
    logical,dimension(:),pointer::f_st=>null()
  end type


  ! phase data
  type::phs_t
    type(prp_t),dimension(:),pointer::prp=>null()
    type(eqn_t),dimension(:),pointer::eqn=>null()
  end type
end module
