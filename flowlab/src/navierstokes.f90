subroutine navierstokes ( &
  fn_len, filename, &
  mesh_blk, mesh_n_blk, mesh_bnd_blk, mesh_dmn_blk, &

  mesh_n_dim, mesh_n_vrt, mesh_vrt_idx, geom_x_vrt, &

  mshf_n_elm, mshf_elm_idx, mshf_elm_ty, mshf_elm_dim, &
  mshf_elm_fce_ptr, mshf_fce_vrt_ptr, mshf_fce_vrt_lst, &

  mesh_n_bnd, mesh_bnd_idx, mesh_bnd_tag_ptr, &
  mesh_bnd_tag_lst, mesh_bnd_vrt_ptr, mesh_bnd_vrt_lst, &

  mesh_n_gbnd, mesh_gbnd_idx, mesh_gbnd_tag_ptr, &
  mesh_gbnd_tag_lst, mesh_gbnd_vrt_ptr, mesh_gbnd_vrt_lst, &

  mesh_n_dmn, mesh_dmn_idx, mesh_dmn_tag_ptr, &
  mesh_dmn_tag_lst, mesh_dmn_vrt_ptr, mesh_dmn_vrt_lst, &

  mesh_n_gdmn, mesh_gdmn_idx, mesh_gdmn_tag_ptr, &
  mesh_gdmn_tag_lst, mesh_gdmn_vrt_ptr, mesh_gdmn_vrt_lst, &

  mesh_n_xbnd, mesh_xbnd_vrt_ptr, mesh_xbnd_vrt_lst, &
  mesh_xbnd_elm_ptr, mesh_xbnd_elm_lst, &

  mesh_n_xdmn, mesh_xdmn_elm_ptr, mesh_xdmn_elm_lst, &
  mesh_xdmn_vrt_ptr, mesh_xdmn_vrt_lst, &

  mesh_bnd_neig_ptr, mesh_bnd_neig_lst, &
  mesh_dmn_neig_ptr, mesh_dmn_neig_lst &
)

!   use mpi
  use constants_m
  use alloc_m
  use solv_m
  use hash_m
  use geometry_m
  use transport_m
  use fson

  implicit none

  integer, intent(in) :: fn_len
  character(len=fn_len), intent(in) :: filename
  integer, intent(in) :: mesh_blk
  integer, intent(in) :: mesh_n_blk
  integer, dimension(mesh_n_blk + 1), intent(in), target :: mesh_bnd_blk
  integer, dimension(mesh_n_blk + 1), intent(in), target :: mesh_dmn_blk


  integer, intent(in) :: mesh_n_dim, mesh_n_vrt
  integer, dimension(mesh_n_vrt), intent(in), target :: mesh_vrt_idx
  real(dp), dimension(3, mesh_n_vrt), intent(in), target :: geom_x_vrt


  integer, intent(in) :: mshf_n_elm
  integer, dimension(20), intent(in), target :: mshf_elm_idx
  integer, dimension(mshf_n_elm), intent(in), target :: mshf_elm_ty
  integer, dimension(mshf_n_elm), intent(in), target :: mshf_elm_dim
  integer, dimension(mshf_n_elm+1), intent(in), target :: mshf_elm_fce_ptr
  integer, dimension(mshf_elm_fce_ptr(mshf_n_elm + 1) + 1), &
    intent(in), target :: mshf_fce_vrt_ptr
  integer, dimension(mshf_fce_vrt_ptr(mshf_elm_fce_ptr(mshf_n_elm + 1) + 1)), &
    intent(in), target :: mshf_fce_vrt_lst


  integer, intent(in) :: mesh_n_bnd
  integer, dimension(mesh_n_bnd), intent(in), target :: mesh_bnd_idx
  integer, dimension(mesh_n_bnd + 1), intent(in), target :: mesh_bnd_tag_ptr
  integer, dimension(mesh_bnd_tag_ptr(mesh_n_bnd + 1)), &
    intent(in), target :: mesh_bnd_tag_lst
  integer, dimension(mesh_n_bnd + 1), intent(in), target :: mesh_bnd_vrt_ptr
  integer, dimension(mesh_bnd_vrt_ptr(mesh_n_bnd + 1)), &
    intent(in), target :: mesh_bnd_vrt_lst

  integer, intent(in) :: mesh_n_gbnd
  integer, dimension(mesh_n_gbnd), intent(in), target :: mesh_gbnd_idx
  integer, dimension(mesh_n_gbnd + 1), intent(in), target :: mesh_gbnd_tag_ptr
  integer, dimension(mesh_gbnd_tag_ptr(mesh_n_gbnd + 1)), &
    intent(in), target :: mesh_gbnd_tag_lst
  integer, dimension(mesh_n_gbnd + 1), intent(in), target :: mesh_gbnd_vrt_ptr
  integer, dimension(mesh_gbnd_vrt_ptr(mesh_n_gbnd + 1)), &
    intent(in), target :: mesh_gbnd_vrt_lst

  integer, intent(in) :: mesh_n_dmn
  integer, dimension(mesh_n_dmn), intent(in), target :: mesh_dmn_idx
  integer, dimension(mesh_n_dmn + 1), intent(in), target :: mesh_dmn_tag_ptr
  integer, dimension(mesh_dmn_tag_ptr(mesh_n_dmn + 1)), &
    intent(in), target :: mesh_dmn_tag_lst
  integer, dimension(mesh_n_dmn + 1), intent(in), target :: mesh_dmn_vrt_ptr
  integer, dimension(mesh_dmn_vrt_ptr(mesh_n_dmn + 1)), &
    intent(in), target :: mesh_dmn_vrt_lst

  integer, intent(in) :: mesh_n_gdmn
  integer, dimension(mesh_n_gdmn), intent(in), target :: mesh_gdmn_idx
  integer, dimension(mesh_n_gdmn + 1), intent(in), target :: mesh_gdmn_tag_ptr
  integer, dimension(mesh_gdmn_tag_ptr(mesh_n_gdmn + 1)), &
    intent(in), target :: mesh_gdmn_tag_lst
  integer, dimension(mesh_n_gdmn + 1), intent(in), target :: mesh_gdmn_vrt_ptr
  integer, dimension(mesh_gdmn_vrt_ptr(mesh_n_gdmn + 1)), &
    intent(in), target :: mesh_gdmn_vrt_lst


  integer, intent(in) :: mesh_n_xbnd
  integer, dimension(mesh_n_xbnd + 1), intent(in), target :: mesh_xbnd_vrt_ptr
  integer, dimension(mesh_xbnd_vrt_ptr(mesh_n_xbnd + 1)), &
    intent(in), target :: mesh_xbnd_vrt_lst
  integer, dimension(mesh_n_xbnd + 1), intent(in), target :: mesh_xbnd_elm_ptr
  integer, dimension(mesh_xbnd_elm_ptr(mesh_n_xbnd + 1)), &
    intent(in), target :: mesh_xbnd_elm_lst

  integer, intent(in) :: mesh_n_xdmn
  integer, dimension(mesh_n_xdmn + 1), intent(in), target :: mesh_xdmn_vrt_ptr
  integer, dimension(mesh_xdmn_vrt_ptr(mesh_n_xdmn + 1)), &
    intent(in), target :: mesh_xdmn_vrt_lst
  integer, dimension(mesh_n_xdmn + 1), intent(in), target :: mesh_xdmn_elm_ptr
  integer, dimension(mesh_xdmn_elm_ptr(mesh_n_xdmn + 1)), &
    intent(in), target :: mesh_xdmn_elm_lst


  integer, dimension(mesh_n_bnd + 1), intent(in), target :: mesh_bnd_neig_ptr
  integer, dimension(mesh_bnd_neig_ptr(mesh_n_bnd + 1)), &
    intent(in), target :: mesh_bnd_neig_lst
  integer, dimension(mesh_n_dmn + 1), intent(in), target :: mesh_dmn_neig_ptr
  integer, dimension(mesh_dmn_neig_ptr(mesh_n_dmn + 1)), &
    intent(in), target :: mesh_dmn_neig_lst


  integer :: err = 0
  integer, dimension(:), allocatable :: mpi_rqst
  integer, dimension(:,:), allocatable :: mpi_stat
  integer :: blk, n_blk, oblk
  type (mshf_t) :: mshf
  type (mesh_t) :: bnd, dmn
  type (geom_t) :: xbnd, xdmn

!   integer :: i, ii, j, jj1, jjn, k, m, n
  integer, dimension(1) :: lst1
  integer, dimension(:), allocatable :: itmp
  real(dp) :: dnum
  real(rp) :: cc, cd, summ, res0(10), gl_res0, time(2), res_drop
  real(rp), dimension(:), allocatable :: stmp
  real(rp), dimension(:,:), allocatable :: vtmp
  type(pde_t) :: vel, pp, pres
  integer :: iter, n_iter, fp_exit, n_gg
  integer :: ii, m, n, nnz, ipp, n_pp
  integer :: i,j,k,kii,kij,kji,kjj,ic,ib,ph_ty
  integer :: j1, j11, j1l, jb, jb1, jbl, jn, jnl, jj1, jjn, jj
  integer :: bnd1,bndl,bndu,dmn1,dmnl,dmnu
  integer :: f_len
  character(sl) :: out_file
  type(sfld_t) :: area
  integer,dimension(:),allocatable::bvrt_idx
  type(fson_value), pointer :: json

  json => fson_parse(trim(filename))
  call fson_get(json, "Material.Density", cc)
  call fson_get(json, "Material.Viscosity", cd)
  call fson_get(json, "Solver.Iteration-Limit", n_iter)
  call fson_get(json, "Solver.Residual-Drop", res_drop)
  call fson_get(json, "Equation.Velocity.Relaxation", vel%urf)
  call fson_get(json, "Equation.Pressure.Relaxation", pp%urf)
  call fson_get(json, "Model.Navier-Stokes.Convection-Blending", vel%bl_c)
  call fson_destroy(json)

  !! initialize communications
  blk = 0; n_blk = 1;
#ifdef USE_MPI
  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)
#endif

  n_proc = n_blk; proc = blk
#ifdef USE_MPI
  allocate (mpi_rqst(n_blk)); mpi_rqst = 0
  allocate (mpi_stat(mpi_status_size, n_blk)); mpi_stat = 0
#endif

  !! mesh format
  mshf%n_elm = mshf_n_elm
  mshf%elm_idx => mshf_elm_idx
  mshf%elm_ty => mshf_elm_ty
  mshf%elm_dim => mshf_elm_dim
  mshf%fce_vrt%blk => mshf_elm_fce_ptr
  mshf%fce_vrt%row => mshf_fce_vrt_ptr
  mshf%fce_vrt%col => mshf_fce_vrt_lst


  !! adjust mesh format to Fortran indexing
  mshf%elm_idx = eoshift (mshf%elm_idx, shift = 1)
  mshf%elm_idx = mshf%elm_idx + 1
  mshf%fce_vrt%blk = mshf%fce_vrt%blk + 1
  mshf%fce_vrt%row = mshf%fce_vrt%row + 1
  mshf%fce_vrt%col = mshf%fce_vrt%col + 1


  !! basic geometry
  dmn%n_dim = mesh_n_dim
  dmn%n_vrt = mesh_n_vrt
  dmn%vrt_idx => mesh_vrt_idx

!   xdmn%x_vrt => geom_x_vrt
  allocate (xdmn%x_vrt(dmn%n_dim, dmn%n_vrt))
  xdmn%x_vrt = geom_x_vrt(:dmn%n_dim,:)


  !! mesh boundary
  bnd%blk = mesh_blk
  bnd%n_blk = mesh_n_blk
  bnd%elm_blk => mesh_bnd_blk
  bnd%n_elm = mesh_n_bnd
  bnd%elm_idx => mesh_bnd_idx
  bnd%elm_tag%row => mesh_bnd_tag_ptr
  bnd%elm_tag%col => mesh_bnd_tag_lst
  bnd%elm_vrt%row => mesh_bnd_vrt_ptr
  bnd%elm_vrt%col => mesh_bnd_vrt_lst
  bnd%elm_neig%row => mesh_bnd_neig_ptr
  bnd%elm_neig%col => mesh_bnd_neig_lst

  bnd%n_ghst = mesh_n_gbnd
  bnd%ghst_idx => mesh_gbnd_idx
  bnd%ghst_tag%row => mesh_gbnd_tag_ptr
  bnd%ghst_tag%col => mesh_gbnd_tag_lst
  bnd%ghst_vrt%row => mesh_gbnd_vrt_ptr
  bnd%ghst_vrt%col => mesh_gbnd_vrt_lst

  bnd%n_fce = mesh_n_xbnd
  bnd%fce_elm%row => mesh_xbnd_elm_ptr
  bnd%fce_elm%col => mesh_xbnd_elm_lst
  bnd%fce_vrt%row => mesh_xbnd_vrt_ptr
  bnd%fce_vrt%col => mesh_xbnd_vrt_lst


  !! adjust mesh boundary to Fortran indexing
  bnd%blk = bnd%blk + 1
  bnd%elm_blk = bnd%elm_blk + 1

  bnd%elm_tag%row = bnd%elm_tag%row + 1
  bnd%elm_vrt%row = bnd%elm_vrt%row + 1
  bnd%elm_neig%row = bnd%elm_neig%row + 1
  bnd%elm_neig%col = bnd%elm_neig%col + 1

  bnd%ghst_tag%row = bnd%ghst_tag%row + 1
  bnd%ghst_vrt%row = bnd%ghst_vrt%row + 1

  bnd%fce_elm%row = bnd%fce_elm%row + 1
  bnd%fce_vrt%row = bnd%fce_vrt%row + 1


  !! mesh domain
  dmn%blk = mesh_blk
  dmn%n_blk = mesh_n_blk
  dmn%elm_blk => mesh_dmn_blk
  dmn%n_elm = mesh_n_dmn
  dmn%elm_idx => mesh_dmn_idx
  dmn%elm_tag%row => mesh_dmn_tag_ptr
  dmn%elm_tag%col => mesh_dmn_tag_lst
  dmn%elm_vrt%row => mesh_dmn_vrt_ptr
  dmn%elm_vrt%col => mesh_dmn_vrt_lst
  dmn%elm_neig%row => mesh_dmn_neig_ptr
  dmn%elm_neig%col => mesh_dmn_neig_lst

  dmn%n_ghst = mesh_n_gdmn
  dmn%ghst_idx => mesh_gdmn_idx
  dmn%ghst_tag%row => mesh_gdmn_tag_ptr
  dmn%ghst_tag%col => mesh_gdmn_tag_lst
  dmn%ghst_vrt%row => mesh_gdmn_vrt_ptr
  dmn%ghst_vrt%col => mesh_gdmn_vrt_lst

  dmn%n_fce = mesh_n_xdmn
  dmn%fce_elm%row => mesh_xdmn_elm_ptr
  dmn%fce_elm%col => mesh_xdmn_elm_lst
  dmn%fce_vrt%row => mesh_xdmn_vrt_ptr
  dmn%fce_vrt%col => mesh_xdmn_vrt_lst


  !! adjust mesh domain to Fortran indexing
  dmn%blk = dmn%blk + 1
  dmn%elm_blk = dmn%elm_blk + bnd%elm_blk(bnd%n_blk + 1)

  dmn%elm_tag%row = dmn%elm_tag%row + 1
  dmn%elm_vrt%row = dmn%elm_vrt%row + 1
  dmn%elm_neig%row = dmn%elm_neig%row + 1
  dmn%elm_neig%col = dmn%elm_neig%col + 1

  dmn%ghst_tag%row = dmn%ghst_tag%row + 1
  dmn%ghst_vrt%row = dmn%ghst_vrt%row + 1

  dmn%fce_elm%row = dmn%fce_elm%row + 1
  dmn%fce_vrt%row = dmn%fce_vrt%row + 1


  !! setup vertex index mapping
  call hash_initialise (dmn%vi2i, dmn%n_vrt)
  do i = 1, dmn%n_vrt
    call hash_insert (dmn%vi2i, dmn%vrt_idx(i), i)
  end do


  !! setup ghost index and block mappings
  call hash_initialise (bnd%gi2i, bnd%n_ghst)
  call hash_initialise (bnd%gi2b, bnd%n_ghst)
  do i=1,bnd%n_ghst
    jjn = bnd%ghst_tag%row(i + 1) - 1
    oblk = (-1) * bnd%ghst_tag%col(jjn)
    call hash_insert (bnd%gi2b, bnd%ghst_idx(i), oblk)
    call hash_insert (bnd%gi2i, bnd%ghst_idx(i), i)
  end do

  call hash_initialise (dmn%gi2i, dmn%n_ghst)
  call hash_initialise (dmn%gi2b, dmn%n_ghst)
  do i=1,dmn%n_ghst
    jjn = dmn%ghst_tag%row(i + 1) - 1
    oblk = (-1) * dmn%ghst_tag%col(jjn)
    call hash_insert (dmn%gi2b, dmn%ghst_idx(i), oblk)
    call hash_insert (dmn%gi2i, dmn%ghst_idx(i), i)
  end do


#ifdef USE_MPI
  call init_ghost_update (bnd, dmn)
#endif

  !! initialize domain geometry
  allocate (xdmn%x_vc%dmn(dmn%n_dim, dmn%n_elm))
  allocate (xdmn%x_vc%ldmn(dmn%n_dim, dmn%n_locl))
  allocate (xdmn%x_vc%gdmn(dmn%n_dim, dmn%n_ghst))

  allocate (xdmn%vol%dmn(dmn%n_elm))
  allocate (xdmn%vol%ldmn(dmn%n_locl))
  allocate (xdmn%vol%gdmn(dmn%n_ghst))

  allocate (xdmn%gl_min%dmn(dmn%n_dim, dmn%n_elm))
  allocate (xdmn%gl_min%ldmn(dmn%n_dim, dmn%n_locl))
  allocate (xdmn%gl_min%gdmn(dmn%n_dim, dmn%n_ghst))
  allocate (xdmn%gl_max%dmn(dmn%n_dim, dmn%n_elm))
  allocate (xdmn%gl_max%ldmn(dmn%n_dim, dmn%n_locl))
  allocate (xdmn%gl_max%gdmn(dmn%n_dim, dmn%n_ghst))
  allocate (xdmn%gl%dmn(dmn%n_dim, dmn%n_elm))
  allocate (xdmn%gl%ldmn(dmn%n_dim, dmn%n_locl))
  allocate (xdmn%gl%gdmn(dmn%n_dim, dmn%n_ghst))

  allocate (xdmn%ls%dmn(dmn%n_dim, dmn%n_dim, dmn%n_elm))
  allocate (xdmn%ls%ldmn(dmn%n_dim, dmn%n_dim, dmn%n_locl))
  allocate (xdmn%ls%gdmn(dmn%n_dim, dmn%n_dim, dmn%n_ghst))
  allocate (xdmn%gg%dmn(dmn%n_dim, dmn%n_dim, dmn%n_elm))
  allocate (xdmn%gg%ldmn(dmn%n_dim, dmn%n_dim, dmn%n_locl))
  allocate (xdmn%gg%gdmn(dmn%n_dim, dmn%n_dim, dmn%n_ghst))

  allocate (xdmn%x_fc(dmn%n_dim, dmn%n_fce))
  allocate (xdmn%x_pc(dmn%n_dim, dmn%n_fce))
  allocate (xdmn%x_nc(dmn%n_dim, dmn%n_fce))
  allocate (xdmn%norm(dmn%n_dim, dmn%n_fce))
  allocate (xdmn%volr(dmn%n_elm))
  allocate (xdmn%area(dmn%n_fce))
  allocate (xdmn%l_pn(dmn%n_fce))
  allocate (xdmn%w_fp(dmn%n_fce))


  !! calculate domain geometry
  call geometry (bnd, dmn, xbnd, xdmn)


  !! check total volume
#ifdef USE_MPI
  if (n_blk > 1) then
    summ = xdmn%vol_sum
    call MPI_Reduce (summ, xdmn%vol_sum, 1, mpi_rp, &
      mpi_sum, 0, mpi_comm_world, err)
  end if
#endif

  if (blk == 0) print "(a,1x,g12.4)","Volume", xdmn%vol_sum


  !! compute LS gradient matrix
  call ls_grad_matrix (bnd, dmn, xbnd, xdmn)


  !! initialize linear solver
  call parms_init (dmn%n_elm)


  !! initialize partial differential equations
  vel%n_cmp = dmn%n_dim

  allocate (vel%bnd_c(vel%n_cmp, 49))
  allocate (vel%bnd_d(vel%n_cmp, 49))

  allocate (vel%grad%dmn(dmn%n_dim, vel%n_cmp, dmn%n_elm))
  allocate (vel%grad%ldmn(dmn%n_dim, vel%n_cmp, dmn%n_locl))
  allocate (vel%grad%gdmn(dmn%n_dim, vel%n_cmp, dmn%n_ghst))

  allocate (vel%phi%dmn(vel%n_cmp, dmn%n_elm))
  allocate (vel%phi%ldmn(vel%n_cmp, dmn%n_locl))
  allocate (vel%phi%gdmn(vel%n_cmp, dmn%n_ghst))
  allocate (vel%phi%bnd(vel%n_cmp, bnd%n_elm))

  allocate (vel%rhs(vel%n_cmp, dmn%n_elm))
  allocate (vel%spd(dmn%n_fce))

  allocate (vel%cc%dmn(dmn%n_elm))
  allocate (vel%cc%ldmn(dmn%n_locl))
  allocate (vel%cc%gdmn(dmn%n_ghst))
  allocate (vel%cc%bnd(bnd%n_elm))

  allocate (vel%cd%dmn(dmn%n_elm))
  allocate (vel%cd%ldmn(dmn%n_locl))
  allocate (vel%cd%gdmn(dmn%n_ghst))
  allocate (vel%cd%bnd(bnd%n_elm))

  allocate (vel%diag%dmn(dmn%n_elm))
  allocate (vel%diag%ldmn(dmn%n_locl))
  allocate (vel%diag%gdmn(dmn%n_ghst))

  m = dmn%n_elm
  nnz = dmn%elm_neig%row(m + 1) - 1
  allocate (vel%mat%idx(m))
  allocate (vel%mat%val(nnz))
  vel%mat%idx = dmn%elm_idx - dmn%elm_blk(1) + 1
  vel%mat%row => dmn%elm_neig%row
  vel%mat%col => dmn%elm_neig%col

  vel%calc_spd = 1
  vel%bl_d = 1
  vel%spd = zero
  vel%cc%dmn = cc; vel%cc%gdmn = cc; vel%cc%bnd = cc
  vel%cd%dmn = cd; vel%cd%gdmn = cd; vel%cd%bnd = cd
  vel%phi%dmn = zero; vel%phi%gdmn = zero; vel%phi%bnd = zero

!   call random_number (vel%phi%dmn)
#ifdef USE_MPI
  if (n_blk > 1) call ghost_update_v (bnd, dmn, vel%phi, 0, 1)
#endif

  vel%rhs = zero
  vel%mat%val = zero


  !! initialize pressure correction equation
  pp%n_cmp = 1

  allocate (pp%grad%dmn(dmn%n_dim, pp%n_cmp, dmn%n_elm))
  allocate (pp%grad%ldmn(dmn%n_dim, pp%n_cmp, dmn%n_locl))
  allocate (pp%grad%gdmn(dmn%n_dim, pp%n_cmp, dmn%n_ghst))

  allocate (pp%phi%dmn(pp%n_cmp, dmn%n_elm))
  allocate (pp%phi%ldmn(pp%n_cmp, dmn%n_locl))
  allocate (pp%phi%gdmn(pp%n_cmp, dmn%n_ghst))
  allocate (pp%phi%bnd(pp%n_cmp, bnd%n_elm))

  allocate (pp%ref(pp%n_cmp))

  pp%spd => vel%spd
  pp%cc%dmn => vel%cc%dmn
  pp%cc%ldmn => vel%cc%ldmn
  pp%cc%gdmn => vel%cc%gdmn
  pp%cc%bnd => vel%cc%bnd
  pp%rhs => vel%rhs
  pp%diag%dmn => vel%diag%dmn
  pp%diag%ldmn => vel%diag%ldmn
  pp%diag%gdmn => vel%diag%gdmn

  pp%mat%idx => vel%mat%idx
  pp%mat%val => vel%mat%val
  pp%mat%row => vel%mat%row
  pp%mat%col => vel%mat%col

  pv_scheme = 1
!   pv_scheme = 2
  pp%ref = zero
  pp%phi%dmn = zero; pp%phi%gdmn = zero; pp%phi%bnd = zero


  !! initialize pressure field
  pres%n_cmp = 1

  allocate (pres%grad%dmn(dmn%n_dim, pres%n_cmp, dmn%n_elm))
  allocate (pres%grad%ldmn(dmn%n_dim, pres%n_cmp, dmn%n_locl))
  allocate (pres%grad%gdmn(dmn%n_dim, pres%n_cmp, dmn%n_ghst))

  allocate (pres%phi%dmn(pres%n_cmp, dmn%n_elm))
  allocate (pres%phi%ldmn(pres%n_cmp, dmn%n_locl))
  allocate (pres%phi%gdmn(pres%n_cmp, dmn%n_ghst))
  allocate (pres%phi%bnd(pres%n_cmp, bnd%n_elm))

  pres%phi%dmn = zero; pres%phi%gdmn = zero; pres%phi%bnd = zero


  call set_boundary_conditions()


  ! workspace
  allocate (wrk(3))
  allocate (wrk(1)%d(size(vel%mat%val)))
  allocate (wrk(2)%d(dmn%n_elm))
  allocate (wrk(3)%d(dmn%n_elm))


  ! iteratively solve PDEs
  n_gg = 1
  n_pp = 1

  do iter = 1, n_iter


    ! initialize workspace
    vel%rhs = zero
    vel%mat%val = zero


    ! pressure gradient
    call ls_gradient (bnd, dmn, xbnd, xdmn, pres%phi, pres%grad, 1)
!     call gg_gradient (bnd, dmn, xbnd, xdmn, pres%phi, pres%grad, 1, n_gg)
#ifdef USE_MPI
    if (n_blk > 1) call ghost_update_t (bnd, dmn, pres%grad, 0, 1)
#endif
!     call gl_local_extrema (bnd, dmn, xbnd, xdmn, pres%phi)
!     call grad_limiter (bnd, dmn, xbnd, xdmn, pres%phi, pres%grad, 2)
!     if (n_blk > 1) call ghost_update_t (bnd, dmn, pres%grad, 0, 1)


    ! velocity gradient
    call ls_gradient (bnd, dmn, xbnd, xdmn, vel%phi, vel%grad, 0)
!     call gg_gradient (bnd, dmn, xbnd, xdmn, vel%phi, vel%grad, 0, n_gg)
#ifdef USE_MPI
    if (n_blk > 1) call ghost_update_t (bnd, dmn, vel%grad, 0, 1)
#endif
!     call gl_local_extrema (bnd, dmn, xbnd, xdmn, vel%phi)
!     call grad_limiter (bnd, dmn, xbnd, xdmn, vel%phi, vel%grad, 2)
!     if (n_blk > 1) call ghost_update_t (bnd, dmn, vel%grad, 0, 1)


    ! momentum equation
    call ns_source (bnd, dmn, xbnd, xdmn, pres, vel)
    call construct_pde (bnd, dmn, xbnd, xdmn, vel)
#ifdef USE_MPI
    if (n_blk > 1) call ghost_update_s (bnd, dmn, vel%diag, 0, 1)
#endif
    do i = 1, vel%n_cmp
      if (i == 1) wrk(1)%d = real (vel%mat%val, dp)
      wrk(2)%d = real (vel%rhs(i,:), dp)
      wrk(3)%d = real (vel%phi%dmn(i,:), dp)

      call parms_solv ( &
        dmn%n_elm, vel%mat%idx, vel%mat%row, vel%mat%col, &
        wrk(1)%d, wrk(2)%d, wrk(3)%d, dnum)

      vel%phi%dmn(i,:) = real (wrk(3)%d, rp)
      res0(i) = dnum

#ifndef USE_MPI
      call lin_eqns_solver ( &
        vel%mat%val,vel%phi%dmn(i,:),vel%rhs(i,:), &
        res0(i),vel%mat%col,vel%mat%row,25,0.01_rp,2,i)
#endif
    end do

#ifdef USE_MPI
    if (n_blk > 1) call ghost_update_v (bnd, dmn, vel%phi, 0, 1)
#endif

    ! initialize workspace
    pp%rhs = zero
    pp%mat%val = zero


    ! pressure correction equation
    call construct_pp (bnd, dmn, xbnd, xdmn, vel, pres, pp)


!     do ipp = 1, n_pp
      i = 4
      pp%phi%dmn(1,:) = zero

      wrk(1)%d = real (pp%mat%val, dp)
      wrk(2)%d = real (pp%rhs(1,:), dp)
      wrk(3)%d = real (pp%phi%dmn(1,:), dp)

      call parms_solv ( &
        dmn%n_elm, pp%mat%idx, pp%mat%row, pp%mat%col, &
        wrk(1)%d, wrk(2)%d, wrk(3)%d, dnum)

      pp%phi%dmn(1,:) = real (wrk(3)%d, rp)
      res0(i) = dnum

#ifndef USE_MPI
      call lin_eqns_solver ( &
        pp%mat%val,pp%phi%dmn(1,:),pp%rhs(1,:), &
        res0(i),vel%mat%col,vel%mat%row,25,0.01_rp,1,i)
#endif

#ifdef USE_MPI
      if (n_blk > 1) call ghost_update_v (bnd, dmn, pp%phi, 0, 1)
#endif

      ! pressure correction gradient
      call ls_gradient (bnd, dmn, xbnd, xdmn, pp%phi, pp%grad, 1)
!       call gg_gradient (bnd, dmn, xbnd, xdmn, pp%phi, pp%grad, 1, n_gg)
#ifdef USE_MPI
      if (n_blk > 1) call ghost_update_t (bnd, dmn, pp%grad, 0, 1)
#endif
!       call gl_local_extrema (bnd, dmn, xbnd, xdmn, pp%phi)
!       call grad_limiter (bnd, dmn, xbnd, xdmn, pp%phi, pp%grad, 2)
!       if (n_blk > 1) call ghost_update_t (bnd, dmn, pp%grad, 0, 1)


      ! correct velocity and pressure
      call ns_correction (bnd, dmn, xbnd, xdmn, pp, vel, pres)
#ifdef USE_MPI
      if (n_blk > 1) call ghost_update_v (bnd, dmn, vel%phi, 0, 1)
      if (n_blk > 1) call ghost_update_v (bnd, dmn, pres%phi, 0, 1)
#endif

!       pp%rhs = zero
!       call correct_skew_pp (bnd, dmn, xbnd, xdmn, vel, pres, pp)
!     end do


    if (blk == 0) then
      call fp_iter_ctrl('F', 1, iter, n_iter, 1, fp_exit, res_drop, res0, gl_res0)
    end if

!     if (mod (iter, n_iter / 10) == 0) call write_fields ()

#ifdef USE_MPI
    call MPI_Bcast (fp_exit, 1, mpi_integer, 0, mpi_comm_world, err)
#endif
    if (fp_exit > 1) exit
  end do


  call write_fields ()


contains
subroutine write_fields ()
  out_file = dir (filename)
  out_file = trim(out_file)//"vel.msh"
  if (n_blk > 1) then
    out_file = trim(out_file)//"_p"//trim(itoa(blk))
  end if
  f_len = len_trim(out_file)

  call write_gmsh_vector_field (f_len, out_file, dmn%n_elm, vel%n_cmp, dmn%elm_idx, vel%phi%dmn)


  out_file = dir (filename)
  out_file = trim(out_file)//"spd.msh"
  if (n_blk > 1) then
    out_file = trim(out_file)//"_p"//trim(itoa(blk))
  end if
  f_len = len_trim(out_file)

  wrk(3)%d = real (vvec_mag(vel%phi%dmn), dp)
  call write_gmsh_field (f_len, out_file, dmn%n_elm, dmn%elm_idx, wrk(3)%d)


  out_file = dir (filename)
  out_file = trim(out_file)//"pres.msh"
  if (n_blk > 1) then
    out_file = trim(out_file)//"_p"//trim(itoa(blk))
  end if
  f_len = len_trim(out_file)

  wrk(3)%d = real (pres%phi%dmn(1,:), dp)
  call write_gmsh_field (f_len, out_file, dmn%n_elm, dmn%elm_idx, wrk(3)%d)
end subroutine


subroutine set_boundary_conditions()
  integer :: bnd1,bndl,bndu,dmn1,dmnl,dmnu
  integer :: i,j,k,kii,kij,kji,kjj,ic,ib,ph_ty
  integer :: j1, j11, j1l, jb, jb1, jbl, jn, jnl, jj1, jjn, jj
  character(len=256):: bnd_type
  integer :: bnd_os, bnd_id
  real(rp),dimension(4) :: bnd_value
  type(fson_value), pointer :: json

  json => fson_parse(trim(filename))
  call fson_get(json, "Boundary.Type", bnd_type)
  call fson_get(json, "Boundary.Index", bnd_id)
  call fson_get(json, "Boundary.Velocity.X-Direction", bnd_value(1))
  call fson_get(json, "Boundary.Velocity.Y-Direction", bnd_value(2))
  call fson_get(json, "Boundary.Velocity.Z-Direction", bnd_value(3))
  call fson_get(json, "Boundary.Velocity.Magnitude", bnd_value(4))
  call fson_destroy(json)

  print "(a,i2,4g12.4)", trim(bnd_type), bnd_id, bnd_value

  call boundary_offset(bnd_type, bnd_os)
  bnd_id = bnd_id + bnd_os

  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  k = 0
  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j11 = (j1 - dmn1) + 1
    j1l = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)

    if (jn < 0) then
      jb1 = ((-1)*jn - bnd1) + 1
      jbl = ((-1)*jn - bndl) + 1
      ph_ty = bnd%elm_tag%row(jbl)
      ph_ty = bnd%elm_tag%col(ph_ty+1)

      if (ph_ty == bnd_id) then
         if (vec_mag(bnd_value(1:3)) < 1.0e-20) then
            vel%phi%bnd(:,jbl) = -xdmn%norm(:,i) * bnd_value(4)
         else
            vel%phi%bnd(1,jbl) = bnd_value(1) * bnd_value(4)
            vel%phi%bnd(2,jbl) = bnd_value(2) * bnd_value(4)
            vel%phi%bnd(3,jbl) = bnd_value(3) * bnd_value(4)
         end if
      end if
    end if
  end do
end subroutine

end subroutine
