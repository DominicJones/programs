module algorithms_m
implicit none
contains



!==============================================================================!
subroutine state_objective(msh,geom,geom1,cc,cd, &
  dwall,delx,vel,pp,pres,spal,bnd,state,obj)
!------------------------------------------------------------------------------!
  use const_m
  use mesh_data_m
  use geom_data_m
  use geometry_m
  use pde_data_m
  use pdes_m
  use pde_utils_m
  use vector_utils_m
  use spalart_allmaras_m
  use misc_utils_m
  use lin_eqns_solvers_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  type(msh_t),intent(in)::msh
  type(geom_t),intent(inout)::geom,geom1
  real(rp),dimension(:),intent(inout)::cc
  real(rp),dimension(:),intent(inout)::cd
  type(pde_t),intent(inout)::dwall,delx,vel,pp,pres,spal
  type(bnd_t),intent(inout)::bnd
  real(rp),dimension(:,:),intent(in)::state
  real(rp),intent(out)::obj
!------------------------------------------------------------------------------!
  real(rp),dimension(size(msh%elm_ja))::a_ij
  real(rp),dimension(msh%n_ebf)::cdmp
  real(rp),dimension(msh%n_ebf)::cdsa
  real(rp),dimension(msh%n_ebf)::cdns
  real(rp),dimension(msh%n_ebf)::d_wall
  real(rp),dimension(msh%n_ebf)::vis_t
  real(rp),dimension(msh%n_fce)::spd
  real(rp),dimension(msh%n_elm)::apr
  real(rp),dimension(msh%n_dx)::xtmp,dx,norm
  real(rp)::pnorm(10),norm_gl,bl
  integer::fp_iter,n_fp_iter(2),fp_exit,fp_count
  integer::i
  character(50)::fnm
!------------------------------------------------------------------------------!
  spd=0.d0

  call nearest_wall(msh,geom,dwall,spd,dwall%cc,dwall%cd,a_ij,pnorm)
  call approximate_dwall(msh,geom,dwall,geom%d_wall)


  geom1%x=geom%x
  call geometry(msh,geom1)


  ! solve flow field
  if(vel%run/=0)then
  vis_t=0.d0

  if(spal%run/=0)then
    print *,"nearest_wall_dist"
    call nearest_wall(msh,geom1,dwall,spd,dwall%cc,dwall%cd,a_ij,pnorm)
    call approximate_dwall(msh,geom1,dwall,geom1%d_wall)
  end if


  print *,"navier_stokes:"
  do fp_iter=1,vel%n_iter
    cdns=cd+vis_t

    call momentum_equation(msh,geom1,vel,pres,spd,cc,cdns,apr,a_ij,pnorm)

    if(pp%run/=0)then
      call continuity_equation(msh,geom1,vel,pres,pp,spd,cc,apr,a_ij,pnorm(4))
    end if

    if(spal%run/=0)then
      call sa_turb_visc(msh,cc,cd,spal%phi(1,:),vis_t)

      if(spal%run==10)then
        call nr_wall_turb_visc(msh,geom1,cc,cd,vel%phi,vis_t)
      end if

      call sa_turb_diff(msh,cc,cd,spal%phi(1,:),cdsa)

      spal%rhs=0.d0
      call sa_turb_source(msh,d_wall,geom%vol,cc,cd, &
        vis_t,spal%phi(1,:),vel%grad,spal%rhs(1,:))

      call scalar_equation(msh,geom1,spal,spd,cc,cdsa,a_ij,pnorm(5))

      do i=1,msh%n_elm
        spal%phi(1,i)=max(zero,spal%phi(1,i))
      end do
    end if

    call fp_iter_ctrl('F',vel%id,fp_iter,vel%n_iter,vel%k_iter,fp_exit,vel%reduc,pnorm,norm_gl)
    if(fp_exit > 1)exit
  end do
  end if


  call objective(msh,geom1,cc,cdns,vel,pres,bnd,obj)
end subroutine

end module
