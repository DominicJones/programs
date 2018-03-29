!==============================================================================!
subroutine gpde_main(fl1, meshfile, fl2, mshffile, fl3, casedir)
!------------------------------------------------------------------------------!
  use const_m
  use alloc_m
  use list_m

  use mesh_format_m
  use mesh_data_m
  use read_mesh_m
  use connectivity_m
  use write_mesh_m
  use write_field_m

  use vector_utils_m
  use matrix_utils_m
  use misc_utils_m

  use geom_data_m
  use geometry_m
  use gradient_m

  use pde_data_m
  use transport_m

  use lin_eqns_solvers_m
  use pdes_m
  use pde_utils_m

  use algorithms_m
  use spalart_allmaras_m
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  integer::fl1, fl2, fl3
  character(fl1)::meshfile
  character(fl2)::mshffile
  character(fl3)::casedir
!------------------------------------------------------------------------------!
  character(250)::arg,fnm,run
  type(msh_fmt_t)::int_fmt,msh_fmt
  type(msh_t)::msh
  type(geom_t)::geom,geom1
  type(pde_t)::dwall,delx,vel,pp,pres,spal
  type(bnd_t)::bnd
  real(rp),dimension(:),allocatable::cc,cd,vis_t,rtmp,rtmp2,rtmp3
  real(rp)::del,state(2),obj,time(10),spd(2),area(2)
  real(rp),dimension(:,:),allocatable::init_x,init_vel,init_spal
  integer::i,j,k,l,ib,ic,id,ip,ph_ty,it=0,n_bvrt,ios1,ios2,runi
  integer,dimension(:),allocatable::bvrt_idx
!==============================================================================!
  call split_fnm(meshfile,msh%path,msh%fnm)
  call split_fnm(mshffile,msh_fmt%path,msh_fmt%fnm)
  int_fmt%path=msh_fmt%path; int_fmt%fnm="gmsh.fmt"

  if(index(msh%fnm,".msh",.true.)>0)msh_fmt%ty="gmsh"
  if(index(msh%fnm,".neu",.true.)>0)msh_fmt%ty="gambit"

  call read_mesh_format(10,int_fmt)
  call read_mesh_format(10,msh_fmt)

  call read_mesh(10,int_fmt,msh_fmt,msh,geom%x)

  call connectivity(msh_fmt,msh)
!------------------------------------------------------------------------------!
  call alloc(geom%dx, ref=geom%x)

  call alloc(geom%x_fc, i=msh%n_dx,j=msh%n_fce)
  call alloc(geom%x_pc, i=msh%n_dx,j=msh%n_fce)
  call alloc(geom%x_nc, i=msh%n_dx,j=msh%n_fce)
  call alloc(geom%norm, i=msh%n_dx,j=msh%n_fce)

  call alloc(geom%l_pn, i=msh%n_fce)
  call alloc(geom%w_fp, i=msh%n_fce)
  call alloc(geom%area, i=msh%n_fce)

  call alloc(geom%x_vc, i=msh%n_dx,j=msh%n_elm)
  call alloc(geom%vol,  i=msh%n_elm)
  call alloc(geom%volr, i=msh%n_elm)

  call alloc(geom%dx_inv, i=msh%n_dx,j=msh%n_dx,k=msh%n_elm)
  call alloc(geom%d_wall, i=msh%n_ebf)


  call alloc(geom1%x, ref=geom%x)

  call alloc(geom1%x_fc, ref=geom%x_fc)
  call alloc(geom1%x_pc, ref=geom%x_pc)
  call alloc(geom1%x_nc, ref=geom%x_nc)
  call alloc(geom1%norm, ref=geom%norm)

  call alloc(geom1%l_pn, ref=geom%l_pn)
  call alloc(geom1%w_fp, ref=geom%w_fp)
  call alloc(geom1%area, ref=geom%area)

  call alloc(geom1%x_vc, ref=geom%x_vc)
  call alloc(geom1%vol,  ref=geom%vol)
  call alloc(geom1%volr, ref=geom%volr)

  call alloc(geom1%dx_inv, ref=geom%dx_inv)
  call alloc(geom1%d_wall, ref=geom%d_wall)
!------------------------------------------------------------------------------!

  call geometry(msh,geom)

  geom%vol_min=minval(geom%vol)
  geom%vol_max=maxval(geom%vol)
  print "(1x,a,3es12.4)","geom[sum,min,max]:", &
    sum(geom%vol),geom%vol_min,geom%vol_max

  call alloc(grad_wk,i=msh%n_dx,j=msh%n_dx,k=msh%n_elm)

  call ls_grad_matrix(msh,geom)
  call test_gradient(msh,geom,1)
  geom%n_gd=1


  ! Re = (cc / cd) * U_{ref} * L_{ref}
  call alloc(cc, i=msh%n_ebf); cc=1
  call alloc(cd, i=msh%n_ebf); cd=one/1000;


  ! distance to nearest wall equation
  dwall%id=5; dwall%run=0
  dwall%n_cmp=1; dwall%k_iter=1
  dwall%reduc=1.d-5; dwall%n_iter=10
  dwall%lls=2; dwall%rls=1.d-5; dwall%nls=25
  dwall%urf=1; dwall%lspd=0; dwall%lcc=0; dwall%lsbc=0; dwall%ldwall=1
  call alloc(dwall%cc,  i=msh%n_ebf)
  call alloc(dwall%cd,  i=msh%n_ebf)
  call alloc(dwall%phi,  i=dwall%n_cmp,j=msh%n_ebf)
  call alloc(dwall%grad, i=msh%n_dx,j=dwall%n_cmp,k=msh%n_elm)
  call alloc(dwall%rhs,  i=dwall%n_cmp,j=msh%n_elm)
  dwall%cc=0.d0; dwall%cd=1.d0


  ! mesh deflection equation
  delx%id=1; delx%run=0
  delx%n_cmp=msh%n_dx; delx%k_iter=1
  delx%reduc=1.d-5; delx%n_iter=10
  delx%lls=2; delx%rls=1.d-1; delx%nls=25
  delx%urf=0.95; delx%lspd=0; delx%lcc=0; delx%lsbc=0
  call alloc(delx%phi,  i=delx%n_cmp,j=msh%n_ebf)
  call alloc(delx%grad, i=msh%n_dx,j=delx%n_cmp,k=msh%n_elm)
  call alloc(delx%rhs,  i=delx%n_cmp,j=msh%n_elm)


  ! momentum equation
  vel%id=2;
  vel%n_cmp=msh%n_dx; vel%k_iter=1
  vel%reduc=1.d-4; vel%n_iter=50
  vel%lls=2; vel%rls=1.d-1; vel%nls=5
  vel%urf=0.55; vel%cbl=0.; vel%dbl=1
  vel%lc_sch=1
  call alloc(vel%ref,  i=vel%n_cmp)
  call alloc(vel%phi,  i=vel%n_cmp,j=msh%n_ebf)
  call alloc(vel%grad, i=msh%n_dx,j=vel%n_cmp,k=msh%n_elm)
  call alloc(vel%rhs,  i=vel%n_cmp,j=msh%n_elm)
  vel%ref(1)=1


  ! pressure correction equation
  pp%id=3; pp%run=1
  pp%n_cmp=1; pp%lspd=0
  pp%lls=1; pp%rls=1.d-1; pp%nls=10
  call alloc(pp%phi,  i=pp%n_cmp,j=msh%n_ebf)
  call alloc(pp%grad, i=msh%n_dx,j=pp%n_cmp,k=msh%n_elm)
  call alloc(pp%rhs,  i=pp%n_cmp,j=msh%n_elm)


  ! pressure field (not a PDE)
  pres%n_cmp=1; pres%urf=0.1
  call alloc(pres%phi,  i=pres%n_cmp,j=msh%n_ebf)
  call alloc(pres%grad, i=msh%n_dx,j=pres%n_cmp,k=msh%n_elm)


  ! Spalart-Allmaras equation (0: off, 1: on, 10: on with wall func)
  spal%id=4; spal%run=0
  spal%n_cmp=1; spal%lspd=0
  spal%lls=2; spal%rls=1.d-1; spal%nls=25
  spal%urf=0.5; spal%cbl=0; spal%dbl=1
  call alloc(spal%phi,  i=spal%n_cmp,j=msh%n_ebf)
  call alloc(spal%grad, i=msh%n_dx,j=spal%n_cmp,k=msh%n_elm)
  call alloc(spal%rhs,  i=spal%n_cmp,j=msh%n_elm)
  spal%ref_mag=cd(1)/cc(1)
  spal%phi=spal%ref_mag * 100


  ! objective function
!   bnd%qnty="pressure loss"; bnd%id=0
!   bnd%qnty="power loss"; bnd%id=0
  bnd%qnty="wall force"; bnd%id=31
!   bnd%qnty="flow uniformity"; bnd%id=21


  ! boundary conditions
  area = zero; spd = zero
  call alloc(bnd%ref, i=msh%n_dx)
  call alloc(bnd%val, i=msh%n_dx)

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    ! wall boundary conditions
    if(ph_ty==31)then
      vel%phi(1,ib)=1.00
    end if

    ! inlet boundary conditions
    if(ph_ty==11)then
!       vel%phi(1,ib)=1.00
!       vel%phi(:,ib)=-geom%norm(:,k) * 1.00
    end if

    select case(ph_ty)
    case(inlet(1):inlet(2))
      area(1) = area(1) + geom%area(k)
      spd(1) = spd(1) + sum(vel%phi(:,ib) * geom%norm(:,k))
      spal%phi(1,ib)=0

    case(outlt(1):outlt(2))
      area(2) = area(2) + geom%area(k)

    case(wall(1):wall(2))
      dwall%phi(1,ib)=0
      spal%phi(1,ib)=0
    end select
  end do

  print *, "USING MOVING WALL BOUNDARY CONDITION"
!   print *, "USING INLET BOUNDARY CONDITION"


  spd(2) = (-1) * spd(1) * (area(1) / area(2))

  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    ph_ty=msh%fce_tag(k)%lst(4)

    select case(ph_ty)
    case(outlt(1):outlt(2))
      vel%phi(:,ib) = geom%norm(:,k) * (area(1) / area(2))
      vel%phi(:,i) = vel%phi(:,ib)
    end select
  end do


  call state_objective(msh,geom,geom1,cc,cd, &
       dwall,delx,vel,pp,pres,spal,bnd,geom%dx,obj)


  fnm=trim(msh%path)//trim(casedir)//"vel.msh"
  call write_field(10,trim(fnm),msh,vel%phi,msh%n_dx)

  fnm=trim(msh%path)//trim(casedir)//"pres.msh"
  call write_field(10,trim(fnm),msh,pres%phi,1)
end subroutine
