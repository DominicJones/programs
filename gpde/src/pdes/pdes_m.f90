module pdes_m
implicit none
contains



!==============================================================================!
subroutine nearest_wall(msh,geom,dwall,spd,cc,cd,a_ij,pnorm)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use pde_data_m
  use transport_m
  use gradient_m
  use matrix_utils_m
  use limiters_m
  use lin_eqns_solvers_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  type(pde_t),intent(inout)::dwall

  real(rp),dimension(:),intent(inout)::spd
  real(rp),dimension(:),intent(inout)::cc
  real(rp),dimension(:),intent(inout)::cd
  real(rp),dimension(:),intent(inout)::a_ij
  real(rp),dimension(:),intent(inout)::pnorm
!------------------------------------------------------------------------------!
  integer::i,ic
  real(rp),dimension(1,msh%n_ebf)::phi_cp
  real(rp),dimension(msh%n_dx,1,msh%n_elm)::grad_cp
!------------------------------------------------------------------------------!
  pnorm=0
  a_ij=0; dwall%rhs=0; spd=0

!   call gradient(msh,geom,dwall%phi,dwall%grad,geom%n_gd,.false.)

  phi_cp = dwall%phi
  grad_cp = dwall%grad
  call gradient(msh,geom,phi_cp,grad_cp,geom%n_gd,.false.)
  dwall%phi = phi_cp
  dwall%grad = grad_cp


  do i=1,msh%n_elm
    dwall%rhs(:,i)=dwall%rhs(:,i) + geom%vol(i)
  end do

  call spd_conv_diff(msh,geom,dwall,spd,cc,cd,a_ij,dwall%rhs)

  ic=1
  call lin_eqns_solver( &
    a_ij,dwall%phi(ic,1:msh%n_elm),dwall%rhs(ic,:), &
    pnorm(ic),msh%elm_ja,msh%elm_ia, &
    dwall%nls,dwall%rls,dwall%lls,dwall%id)
end subroutine



!==============================================================================!
subroutine mesh_perturbation(msh,geom,delx,spd,cc,cd,a_ij,pnorm)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use pde_data_m
  use transport_m
  use gradient_m
  use matrix_utils_m
  use limiters_m
  use lin_eqns_solvers_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  type(pde_t),intent(inout)::delx

  real(rp),dimension(:),intent(inout)::spd
  real(rp),dimension(:),intent(inout)::cc
  real(rp),dimension(:),intent(inout)::cd
  real(rp),dimension(:),intent(inout)::a_ij
  real(rp),dimension(:),intent(inout)::pnorm
!------------------------------------------------------------------------------!
  integer::i,ic
!------------------------------------------------------------------------------!
  pnorm=0
  a_ij=0; delx%rhs=0; spd=0

  call gradient(msh,geom,delx%phi,delx%grad,geom%n_gd,.false.)

  call spd_conv_diff(msh,geom,delx,spd,cc,cd,a_ij,delx%rhs)

  do ic=1,msh%n_dx
    call lin_eqns_solver( &
      a_ij,delx%phi(ic,1:msh%n_elm),delx%rhs(ic,:), &
      pnorm(ic),msh%elm_ja,msh%elm_ia, &
      delx%nls,delx%rls,delx%lls,delx%id)
  end do
end subroutine



!==============================================================================!
subroutine momentum_equation(msh,geom,vel,pres,spd,cc,cd,apr,a_ij,pnorm)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use pde_data_m
  use transport_m
  use gradient_m
  use matrix_utils_m
  use limiters_m
  use lin_eqns_solvers_m
  use sparse_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  type(pde_t),intent(inout)::vel
  type(pde_t),intent(inout)::pres

  real(rp),dimension(:),intent(inout)::spd
  real(rp),dimension(:),intent(inout)::cc
  real(rp),dimension(:),intent(inout)::cd
  real(rp),dimension(:),intent(inout)::apr
  real(rp),dimension(:),intent(inout)::a_ij
  real(rp),dimension(:),intent(inout)::pnorm
!------------------------------------------------------------------------------!
  integer::i,ic
!   real(rp),dimension(msh%n_dx,msh%n_elm)::umin,umax,lim
!------------------------------------------------------------------------------!
  pnorm=0
  a_ij=0; vel%rhs=0

  call gradient(msh,geom,pres%phi,pres%grad,geom%n_gd,.true.)
  call gradient(msh,geom,vel%phi,vel%grad,geom%n_gd,.false.)

!   if(pres%llim>0)then
!     call local_extrema(msh%fce_elm,pres%phi, &
!       umin(:1,:),umax(:1,:))
!     call grad_limiter(msh%fce_elm,geom,pres%phi,pres%grad, &
!       umin(:1,:),umax(:1,:),lim(:1,:),pres%clim,pres%llim)
!   end if

  do i=1,msh%n_elm
    vel%rhs(:,i)=vel%rhs(:,i)-pres%grad(:,1,i)*geom%vol(i)
  end do

  call spd_conv_diff(msh,geom,vel,spd,cc,cd,a_ij,vel%rhs)

  do ic=1,msh%n_dx
    call lin_eqns_solver( &
      a_ij,vel%phi(ic,1:msh%n_elm),vel%rhs(ic,:), &
      pnorm(ic),msh%elm_ja,msh%elm_ia, &
      vel%nls,vel%rls,vel%lls,vel%id)
  end do

  apr=csr_diag(a_ij,msh%elm_ja,msh%elm_ia)
  apr=1.d0/(apr+small)
end subroutine



!==============================================================================!
subroutine continuity_equation(msh,geom,vel,pres,pp,spd,cc,apr,a_ij,pnorm)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use pde_data_m
  use transport_m
  use gradient_m
  use matrix_utils_m
  use lin_eqns_solvers_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  type(pde_t),intent(inout)::vel
  type(pde_t),intent(inout)::pres
  type(pde_t),intent(inout)::pp

  real(rp),dimension(:),intent(inout)::spd
  real(rp),dimension(:),intent(inout)::cc
  real(rp),dimension(:),intent(inout)::apr
  real(rp),dimension(:),intent(inout)::a_ij
  real(rp),intent(inout)::pnorm
!------------------------------------------------------------------------------!
  integer::i,j,k,ic
  real(rp)::c_lim,c_rat,dpres
  real(rp),dimension(msh%n_dx)::dvel
!------------------------------------------------------------------------------!
  pnorm=0;
  a_ij=0; pp%rhs=0

  call pp_setup(msh,geom,vel,pres,spd,cc,apr,a_ij,pp%rhs)

  pp%phi=0
  call lin_eqns_solver( &
    a_ij,pp%phi(1,:msh%n_elm),pp%rhs(1,:), &
    pnorm,msh%elm_ja,msh%elm_ia, &
    pp%nls,pp%rls,pp%lls,pp%id)

  call gradient(msh,geom,pp%phi,pp%grad,geom%n_gd,.true.)

  do i=1,msh%n_elm
    dvel=-pp%grad(:,1,i)
    dvel=dvel*apr(i)*geom%vol(i)
    vel%phi(:,i)=vel%phi(:,i)+dvel

    dpres=pres%urf*(pp%phi(1,i)-pp%phi(1,1))
    pres%phi(1,i)=pres%phi(1,i)+dpres
  end do
end subroutine



!==============================================================================!
subroutine scalar_equation(msh,geom,scal,spd,cc,cd,a_ij,pnorm)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use pde_data_m
  use transport_m
  use gradient_m
  use matrix_utils_m
  use limiters_m
  use lin_eqns_solvers_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(in)::geom
  type(pde_t),intent(inout)::scal

  real(rp),dimension(:),intent(inout)::spd
  real(rp),dimension(:),intent(inout)::cc
  real(rp),dimension(:),intent(inout)::cd
  real(rp),dimension(:),intent(inout)::a_ij
  real(rp),intent(inout)::pnorm
!------------------------------------------------------------------------------!
  integer::i
!------------------------------------------------------------------------------!
  pnorm=0;
  a_ij=0; scal%rhs=0;

  call gradient(msh,geom,scal%phi,scal%grad,geom%n_gd,.false.)
  call spd_conv_diff(msh,geom,scal,spd,cc,cd,a_ij,scal%rhs)

  call lin_eqns_solver( &
    a_ij,scal%phi(1,1:msh%n_elm),scal%rhs(1,:), &
    pnorm,msh%elm_ja,msh%elm_ia,scal%nls,scal%rls,scal%lls,scal%id)
end subroutine

end module
