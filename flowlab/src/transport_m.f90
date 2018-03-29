module transport_m
use constants_m!, rp => rsp, mpi_rp => mpi_rsp
use hash_m
use graph_m
use geometry_m
use vector_utils_m
use matrix_utils_m

implicit none

integer, dimension(2), parameter :: inlet = (/10, 19/)
integer, dimension(2), parameter :: outlt = (/20, 29/)
integer, dimension(2), parameter ::  wall = (/30, 39/)
integer, dimension(2), parameter ::  symm = (/40, 49/)

type mat_t
  integer :: id = 0
  integer, dimension(:), pointer :: idx => null()
  integer, dimension(:), pointer :: row => null()
  integer, dimension(:), pointer :: col => null()
  real(rp), dimension(:), pointer :: val => null()
end type

type pde_t
  integer :: id = 0, n_cmp = 1, std_bc = 1
  integer :: calc_spd = 0, calc_c = 1, calc_d = 1
  real(rp) :: urf = 0.4, bl_c = 0.7, bl_d = one
  real(rp), dimension(:,:), pointer :: bnd_c => null()
  real(rp), dimension(:,:), pointer :: bnd_d => null()
  real(rp), dimension(:), pointer :: spd => null()
  type(sfld_t) :: cc, cd
  real(rp), dimension(:,:), pointer :: rhs => null()
  type(sfld_t) :: diag
  type(mat_t) :: mat
  type(vfld_t) :: phi
  real(rp), dimension(:), pointer :: ref => null()
  type(tfld_t) :: grad
end type

integer :: pv_scheme = 1
logical, private :: last_fp_iter

contains



subroutine boundary_offset(bnd_type, bnd_id)
  character(len=256) :: bnd_type
  integer :: bnd_id
  bnd_id = 0;
  if (trim(bnd_type) == "Inlet")    bnd_id = inlet(1)
  if (trim(bnd_type) == "Outlet")   bnd_id = outlt(1)
  if (trim(bnd_type) == "Wall")     bnd_id = wall(1)
  if (trim(bnd_type) == "Symmetry") bnd_id = symm(1)
end subroutine


subroutine construct_pde (bnd, dmn, xbnd, xdmn, pde)
  type(mesh_t), intent(in) :: bnd, dmn
  type(geom_t), intent(in) :: xbnd, xdmn
  type(pde_t), intent(inout) :: pde

  real(rp), dimension(dmn%n_dim, dmn%n_dim) :: fgrad,tr
  real(rp), dimension(dmn%n_dim) :: dx,tang,norm1,delta
  real(rp), dimension(pde%n_cmp) :: c_imp,c_exp,c_cor
  real(rp), dimension(pde%n_cmp) :: d_imp,d_exp,d_cor
  real(rp), dimension(pde%n_cmp) :: vphi,nphi,fphi
  real(rp) :: w_fp1, w_fn1, area1, reduc, summ
  real(rp) :: c,c_ij,c_ji,fcc,d,d_ij,fcd,urfr,c_sum(3)
  integer :: n_blk, blk, err
  integer :: i,j,k,kii,kij,kji,kjj,ic,ib,ph_ty
  integer :: j1, j11, j1l, jb, jb1, jbl, jn, jnl, jj1, jjn
  integer::bnd1,bndl,bndu,dmn1,dmnl,dmnu

  blk = 1; n_blk = 1;
#ifdef USE_MPI
  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)
#endif

  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  pde%diag%dmn = zero
  pde%bnd_c = zero
  pde%bnd_d = zero


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j11 = (j1 - dmn1) + 1
    j1l = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)

    kii = csr_index (j1l, j11, pde%mat%col, pde%mat%row)

    area1 = xdmn%area(i)
    norm1 = xdmn%norm(:,i)


    ! internal and internal-ghost faces:
    select case (jn)
    case (1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb1 = (jn - dmn1) + 1
        jbl = (jn - dmnl) + 1

        kjj = csr_index (jbl, jb1, pde%mat%col, pde%mat%row)
        kij = csr_index (j1l, jb1, pde%mat%col, pde%mat%row)
        kji = csr_index (jbl, j11, pde%mat%col, pde%mat%row)


        w_fp1 = xdmn%w_fp(i)
        w_fn1 = one - w_fp1

        delta = xdmn%x_vc%dmn(:,jbl) - xdmn%x_vc%dmn(:,j1l)
        delta = unit_vec (delta)
        delta = norm1 - delta

        fgrad = w_fp1 * pde%grad%dmn(:,:,jbl) + w_fn1 * pde%grad%dmn(:,:,j1l)
        dx = w_fp1 * xdmn%x_vc%dmn(:,jbl) + w_fn1 * xdmn%x_vc%dmn(:,j1l)
        dx = xdmn%x_fc(:,i) - dx
        do ic = 1, pde%n_cmp
          fphi(ic) = dot_product (fgrad(:,ic), dx)
        end do
        fphi = w_fp1 * pde%phi%dmn(:,jbl) + w_fn1 * pde%phi%dmn(:,j1l) + fphi

        fcc = w_fp1 * pde%cc%dmn(jbl) + w_fn1 * pde%cc%dmn(j1l)
        fcd = w_fp1 * pde%cd%dmn(jbl) + w_fn1 * pde%cd%dmn(j1l)

        c = fcc * area1 * pde%spd(i)
        d = fcd * area1 / xdmn%l_pn(i)

        c_ij = min (c, zero)
        c_ji = -max (c, zero)
        d_ij = -d

        if (c > zero) then
          c_imp = pde%phi%dmn(:,j1l)
        else
          c_imp = pde%phi%dmn(:,jbl)
        end if
        c_exp = fphi
        c_cor = pde%bl_c * c * (c_imp - c_exp)

        do ic = 1, pde%n_cmp
          d_cor(ic) = dot_product (fgrad(:,ic), delta)
        end do
        d_cor = d_cor * fcd * area1

        if(pde%calc_spd /= 0 .and. pde%calc_d /= 0)then
          do ic = 1, pde%n_cmp
            d_exp(ic) = dot_product (fgrad(ic,:), norm1)
          end do
          d_exp = d_exp * fcd * area1
          d_cor = d_cor + d_exp
        end if


        ! matrix off-diagonal contribution
        pde%mat%val(kij) = pde%mat%val(kij) + c_ij + d_ij
        pde%mat%val(kji) = pde%mat%val(kji) + c_ji + d_ij

        ! matrix diagonal contribution
        pde%mat%val(kii) = pde%mat%val(kii) - c_ij - d_ij
        pde%mat%val(kjj) = pde%mat%val(kjj) - c_ji - d_ij

        ! right-hand side contribution
        pde%rhs(:,j1l) = pde%rhs(:,j1l) + c_cor + d_cor
        pde%rhs(:,jbl) = pde%rhs(:,jbl) - c_cor - d_cor

        pde%diag%dmn(j1l) = pde%diag%dmn(j1l) + pde%mat%val(kij) !!
        pde%diag%dmn(jbl) = pde%diag%dmn(jbl) + pde%mat%val(kji) !!


      ! internal-ghost face:
      else
        jb1 = (jn - dmn1) + 1
        call hash_get_value (dmn%gi2i, jn, jbl)

        kij = csr_index (j1l, jb1, pde%mat%col, pde%mat%row)

        w_fp1 = half
        w_fn1 = one - w_fp1

        delta = xdmn%x_vc%gdmn(:,jbl) - xdmn%x_vc%dmn(:,j1l)
        delta = unit_vec (delta)
        delta = norm1 - delta

        fgrad = w_fp1 * pde%grad%gdmn(:,:,jbl) + w_fn1 * pde%grad%dmn(:,:,j1l)
        dx = w_fp1 * xdmn%x_vc%gdmn(:,jbl) + w_fn1 * xdmn%x_vc%dmn(:,j1l)
        dx = xdmn%x_fc(:,i) - dx
        do ic = 1, pde%n_cmp
          fphi(ic) = dot_product (fgrad(:,ic), dx)
        end do
        fphi = w_fp1 * pde%phi%gdmn(:,jbl) + w_fn1 * pde%phi%dmn(:,j1l) + fphi

        fcc = w_fp1 * pde%cc%gdmn(jbl) + w_fn1 * pde%cc%dmn(j1l)
        fcd = w_fp1 * pde%cd%gdmn(jbl) + w_fn1 * pde%cd%dmn(j1l)

        c = fcc * area1 * pde%spd(i)
        d = fcd * area1 / xdmn%l_pn(i)

        c_ij = min (c, zero)
        c_ji = -max (c, zero)
        d_ij = -d

        if (c > zero) then
          c_imp = pde%phi%dmn(:,j1l)
        else
          c_imp = pde%phi%gdmn(:,jbl)
        end if
        c_exp = fphi
        c_cor = pde%bl_c * c * (c_imp - c_exp)

        do ic = 1, pde%n_cmp
          d_cor(ic) = dot_product (fgrad(:,ic), delta)
        end do
        d_cor = d_cor * fcd * area1

        if(pde%calc_spd /= 0 .and. pde%calc_d /= 0)then
          do ic = 1, pde%n_cmp
            d_exp(ic) = dot_product (fgrad(ic,:), norm1)
          end do
          d_exp = d_exp * fcd * area1
          d_cor = d_cor + d_exp
        end if


        ! matrix off-diagonal contribution
        pde%mat%val(kij) = pde%mat%val(kij) + c_ij + d_ij

        ! matrix diagonal contribution
        pde%mat%val(kii) = pde%mat%val(kii) - c_ij - d_ij

        ! right-hand side contribution
        pde%rhs(:,j1l) = pde%rhs(:,j1l) + c_cor + d_cor

        pde%diag%dmn(j1l) = pde%diag%dmn(j1l) + pde%mat%val(kij)
      end if


    case default
      jb1 = ((-1)*jn - bnd1) + 1
      jbl = ((-1)*jn - bndl) + 1

      ph_ty = bnd%elm_tag%row(jbl)
      ph_ty = bnd%elm_tag%col(ph_ty+1)

      fgrad = pde%grad%dmn(:,:,j1l)
      fphi = pde%phi%bnd(:,jbl)
      vphi = pde%phi%dmn(:,j1l) !!

      fcc = pde%cc%bnd(jbl)
      fcd = pde%cd%bnd(jbl)

      c = fcc * area1 * pde%spd(i)
      d = fcd * area1 / xdmn%l_pn(i)

      c_ij = min (c, zero)
      d_ij = -d

      d_imp = fphi - pde%phi%dmn(:,j1l)

      if (pde%calc_spd /= 0 .and. pde%calc_d /= 0) then
        select case (ph_ty)
        case (wall(1):wall(2))
          d_exp = fphi - (vphi - norm1 * dot_product (vphi, norm1))
          d_exp = d_exp / xdmn%l_pn(i)
        case (symm(1):symm(2))
          d_exp = 2 * (zero - norm1 * dot_product (vphi, norm1))
          d_exp = d_exp / xdmn%l_pn(i)
        case default
          tr = transpose (fgrad)
          d_exp = matmul (fgrad + tr, norm1)
        end select
      else
        do ic = 1, pde%n_cmp
          d_exp(ic) = sum (fgrad(:,ic) * norm1)
        end do
      end if

      d_cor = fcd * area1 * d_exp - d * d_imp

      if (pde%std_bc /= 0) then
        select case (ph_ty)
        case (outlt(1):outlt(2))
          fphi = vphi
          c_ij = zero
          d_ij = zero
        case (symm(1):symm(2))
          if (pde%calc_spd == 0) then
            fphi = vphi
            c_ij = zero
            d_ij = zero
          end if
          if (pde%calc_spd /= 0 .and. pde%calc_d /= 0) then
            fphi = vphi - norm1 * dot_product (vphi, norm1)
          end if
        end select
      end if

      ! matrix diagonal contribution
      pde%mat%val(kii) = pde%mat%val(kii) - c_ij - d_ij

      ! right-hand side contribution
      pde%rhs(:,j1l) = pde%rhs(:,j1l) - (c_ij + d_ij) * fphi + d_cor

      ! update boundary value
      pde%phi%bnd(:,jbl) = fphi

      ! store boundary fluxes
      pde%bnd_c(:,ph_ty) = pde%bnd_c(:,ph_ty) + c * fphi
      pde%bnd_d(:,ph_ty) = pde%bnd_d(:,ph_ty) + d * d_exp
    end select
  end do

  urfr = one / pde%urf

  do i = 1, dmn%n_elm
    j1l = i
    j11 = i + (dmnl - dmn1)
    kii = csr_index (j1l, j11, pde%mat%col, pde%mat%row)

    vphi = pde%phi%dmn(:,j1l)
    pde%mat%val(kii) = pde%mat%val(kii) * urfr
    pde%rhs(:,j1l) = pde%rhs(:,j1l) &
      + (one - pde%urf) * pde%mat%val(kii) * vphi

    if (pv_scheme == 1) then
      pde%diag%dmn(j1l) = one / pde%mat%val(kii)
    else if (pv_scheme == 2) then
      pde%diag%dmn(j1l) = pde%mat%val(kii) + pde%diag%dmn(j1l)
      pde%diag%dmn(j1l) = one / (pde%diag%dmn(j1l) + small)
    end if
  end do
end subroutine




subroutine flux_out_corr (bnd, dmn, xbnd, xdmn, vel, c_rat)
  type(mesh_t), intent(in) :: bnd, dmn
  type(geom_t), intent(in) :: xbnd, xdmn
  type(pde_t), intent(in) :: vel
  real(rp), intent(out) :: c_rat

  real(rp), parameter :: c_lim = ten**(3)
  real(rp) :: c, c_in, c_out, c_io0(2), c_io1(2)
  integer :: n_blk, err
  integer :: i, jj1, jjn, j1, j11, j1l, jn, jb1, jbl, ph_ty
  integer :: bnd1, bndl, bndu, dmn1, dmnl, dmnu

  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  c_in = zero; c_out = zero

  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j11 = (j1 - dmn1) + 1
    j1l = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)

    ! boundary faces:
    select case(jn)
    case(:-1)
      jb1 = ((-1)*jn - bnd1) + 1
      jbl = ((-1)*jn - bndl) + 1

      ph_ty = bnd%elm_tag%row(jbl)
      ph_ty = bnd%elm_tag%col(ph_ty+1)

      select case (ph_ty)
      case (inlet(1):inlet(2))
        c = vel%cc%bnd(jbl) * xdmn%area(i) * vel%spd(i)
        c_in = c_in + c
      case (outlt(1):outlt(2))
        c = vel%cc%bnd(jbl) * xdmn%area(i) * vel%spd(i)
        c_out = c_out + c
      end select
    end select
  end do

  n_blk = 1;
#ifdef USE_MPI
  call MPI_Comm_size (mpi_comm_world, n_blk, err)
#endif

  if (n_blk > 1) then
    c_io0(1) = c_in
    c_io0(2) = c_out

#ifdef USE_MPI
    call MPI_Allreduce (c_io0, c_io1, 2, mpi_rp, &
      mpi_sum, mpi_comm_world, err)
#endif

    c_in = c_io1(1)
    c_out = c_io1(2)
  end if

  c_rat = c_in / (c_out + small)
  c_rat = abs (c_rat)
  c_rat = min (c_lim, c_rat)
end subroutine




subroutine construct_pp (bnd, dmn, xbnd, xdmn, vel, pres, pp)
  type(mesh_t), intent(in) :: bnd, dmn
  type(geom_t), intent(in) :: xbnd, xdmn
  type(pde_t), intent(inout) :: vel
  type(pde_t), intent(in) :: pres
  type(pde_t), intent(inout) :: pp

  real(rp), dimension(dmn%n_dim, dmn%n_dim) :: fgrad
  real(rp), dimension(dmn%n_dim) :: dx, fvel, norm1
  real(rp) :: phi_pc, phi_nc, c_sum(3), summ, reduc
  real(rp) :: w_fp1, w_fn1, fcc, fapr, fvol, c, d_ij, d_cor, c_rat
  integer :: n_blk, blk, err
  integer :: kii, kij, kji, kjj
  integer :: ic, i, jj1, jjn, j1, j11, j1l, jn, jb1, jbl, ph_ty
  integer :: bnd1, bndl, bndu, dmn1, dmnl, dmnu

  blk = 1; n_blk = 1;
#ifdef USE_MPI
  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)
#endif

  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  c_sum = zero


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j11 = (j1 - dmn1) + 1
    j1l = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)

    kii = csr_index (j1l, j11, pp%mat%col, pp%mat%row)

    norm1 = xdmn%norm(:,i)


    ! internal and internal-ghost faces:
    select case(jn)
    case(1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb1 = (jn - dmn1) + 1
        jbl = (jn - dmnl) + 1

        kjj = csr_index (jbl, jb1, pp%mat%col, pp%mat%row)
        kij = csr_index (j1l, jb1, pp%mat%col, pp%mat%row)
        kji = csr_index (jbl, j11, pp%mat%col, pp%mat%row)

        w_fp1 = xdmn%w_fp(i)
        w_fn1 = one - w_fp1

        fgrad = w_fp1 * vel%grad%dmn(:,:,jbl) + w_fn1 * vel%grad%dmn(:,:,j1l)
        dx = w_fp1 * xdmn%x_vc%dmn(:,jbl) + w_fn1 * xdmn%x_vc%dmn(:,j1l)
        dx = xdmn%x_fc(:,i) - dx
        do ic = 1, dmn%n_dim
          fvel(ic) = dot_product (fgrad(:,ic), dx)
        end do
        fvel = w_fp1 * vel%phi%dmn(:,jbl) + w_fn1 * vel%phi%dmn(:,j1l) + fvel
        vel%spd(i) = dot_product (fvel, norm1)

        fcc = w_fp1 * vel%cc%dmn(jbl) + w_fn1 * vel%cc%dmn(j1l)
        fapr = w_fp1 * vel%diag%dmn(jbl) + w_fn1 * vel%diag%dmn(j1l)
        d_ij = -fcc * fapr * xdmn%area(i)**2

        c = -fcc * xdmn%area(i) * vel%spd(i)

        fgrad(:,1) = w_fp1 * pres%grad%dmn(:,1,jbl) + w_fn1 * pres%grad%dmn(:,1,j1l)
        dx = xdmn%x_vc%dmn(:,jbl) - xdmn%x_vc%dmn(:,j1l)
        phi_pc = pres%phi%dmn(1,j1l)
        phi_nc = pres%phi%dmn(1,jbl)
        d_cor = phi_nc - phi_pc
        d_cor = d_cor - dot_product (fgrad(:,1), dx)
        d_cor = -d_cor * d_ij


        ! matrix off-diagonal contribution
        pp%mat%val(kij) = pp%mat%val(kij) + d_ij
        pp%mat%val(kji) = pp%mat%val(kji) + d_ij

        ! matrix diagonal contribution
        pp%mat%val(kii) = pp%mat%val(kii) - d_ij
        pp%mat%val(kjj) = pp%mat%val(kjj) - d_ij

        ! right-hand side contribution
        pp%rhs(1,j1l) = pp%rhs(1,j1l) + c + d_cor
        pp%rhs(1,jbl) = pp%rhs(1,jbl) - c - d_cor


      ! internal-ghost face:
      else
        jb1 = (jn - dmn1) + 1
        call hash_get_value (dmn%gi2i, jn, jbl)

        kij = csr_index (j1l, jb1, pp%mat%col, pp%mat%row)

        w_fp1 = half
        w_fn1 = one - w_fp1

        fgrad = w_fp1 * vel%grad%gdmn(:,:,jbl) + w_fn1 * vel%grad%dmn(:,:,j1l)
        dx = w_fp1 * xdmn%x_vc%gdmn(:,jbl) + w_fn1 * xdmn%x_vc%dmn(:,j1l)
        dx = xdmn%x_fc(:,i) - dx
        do ic = 1, dmn%n_dim
          fvel(ic) = dot_product (fgrad(:,ic), dx)
        end do
        fvel = w_fp1 * vel%phi%gdmn(:,jbl) + w_fn1 * vel%phi%dmn(:,j1l) + fvel
        vel%spd(i) = dot_product (fvel, norm1)

        fcc = w_fp1 * vel%cc%gdmn(jbl) + w_fn1 * vel%cc%dmn(j1l)
        fapr = w_fp1 * vel%diag%gdmn(jbl) + w_fn1 * vel%diag%dmn(j1l)
        d_ij = -fcc * fapr * xdmn%area(i)**2

        c = -fcc * xdmn%area(i) * vel%spd(i)

        fgrad(:,1) = w_fp1 * pres%grad%gdmn(:,1,jbl) + w_fn1 * pres%grad%dmn(:,1,j1l)
        dx = xdmn%x_vc%gdmn(:,jbl) - xdmn%x_vc%dmn(:,j1l)
        phi_pc = pres%phi%dmn(1,j1l)
        phi_nc = pres%phi%gdmn(1,jbl)
        d_cor = phi_nc - phi_pc
        d_cor = d_cor - dot_product (fgrad(:,1), dx)
        d_cor = -d_cor * d_ij

        ! matrix off-diagonal contribution
        pp%mat%val(kij) = pp%mat%val(kij) + d_ij

        ! matrix diagonal contribution
        pp%mat%val(kii) = pp%mat%val(kii) - d_ij

        ! right-hand side contribution
        pp%rhs(1,j1l) = pp%rhs(1,j1l) + c + d_cor
      end if


    case default
      jb1 = ((-1)*jn - bnd1) + 1
      jbl = ((-1)*jn - bndl) + 1

      ph_ty = bnd%elm_tag%row(jbl)
      ph_ty = bnd%elm_tag%col(ph_ty+1)

      fvel = vel%phi%bnd(:,jbl)
      vel%spd(i) = dot_product (fvel, norm1)

      c = -vel%cc%bnd(jbl) * xdmn%area(i) * vel%spd(i)

      select case (ph_ty)
      case (outlt(1):outlt(2))
        if (vel%spd(i) < zero) vel%spd(i) = zero
        c = zero
      end select

      pp%rhs(1,j1l) = pp%rhs(1,j1l) + c

    end select
  end do


  call flux_out_corr (bnd, dmn, xbnd, xdmn, vel, c_rat)


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j11 = (j1 - dmn1) + 1
    j1l = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)

    kii = csr_index (j1l, j11, pp%mat%col, pp%mat%row)

    select case(jn)
    case(:-1)
      jb1 = ((-1)*jn - bnd1) + 1
      jbl = ((-1)*jn - bndl) + 1

      ph_ty = bnd%elm_tag%row(jbl)
      ph_ty = bnd%elm_tag%col(ph_ty+1)

      select case (ph_ty)
      case (outlt(1):outlt(2))
        vel%spd(i) = vel%spd(i) * c_rat
        c = -vel%cc%bnd(jbl) * xdmn%area(i) * vel%spd(i)
        pp%rhs(1,j1l) = pp%rhs(1,j1l) + c
      end select
    end select
  end do
end subroutine




subroutine correct_skew_pp (bnd, dmn, xbnd, xdmn, vel, pres, pp)
  type(mesh_t), intent(in) :: bnd, dmn
  type(geom_t), intent(in) :: xbnd, xdmn
  type(pde_t), intent(in) :: vel
  type(pde_t), intent(in) :: pres
  type(pde_t), intent(inout) :: pp

  real(rp), dimension(dmn%n_dim) :: fgrad, dx
  real(rp) :: phi_pc, phi_nc, c_sum(3), summ, reduc
  real(rp) :: w_fp1, fcc, fapr, fvol, c, d_ij, d_cor, c_cor
  integer :: n_blk, blk, err
  integer :: kii, kij, kji, kjj
  integer :: i, jj1, jjn, j1, j11, j1l, jn, jb1, jbl, ph_ty
  integer :: bnd1, bndl, bndu, dmn1, dmnl, dmnu

  blk = 1; n_blk = 1;
#ifdef USE_MPI
  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)
#endif

  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j11 = (j1 - dmn1) + 1
    j1l = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)

    kii = csr_index (j1l, j11, pp%mat%col, pp%mat%row)


    ! internal and internal-ghost faces:
    select case(jn)
    case(1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb1 = (jn - dmn1) + 1
        jbl = (jn - dmnl) + 1

        kjj = csr_index (jbl, jb1, pp%mat%col, pp%mat%row)
        kij = csr_index (j1l, jb1, pp%mat%col, pp%mat%row)
        kji = csr_index (jbl, j11, pp%mat%col, pp%mat%row)

        dx = xdmn%x_pc(:,i) - xdmn%x_vc%dmn(:,j1l)
        phi_pc = dot_product (pp%grad%dmn(:,1,j1l), dx)

        dx = xdmn%x_nc(:,i) - xdmn%x_vc%dmn(:,jbl)
        phi_nc = dot_product (pp%grad%dmn(:,1,jbl), dx)

        d_ij = pp%mat%val(kij)
        d_cor = -d_ij * (phi_nc - phi_pc)

        pp%rhs(1,j1l) = pp%rhs(1,j1l) + d_cor
        pp%rhs(1,jbl) = pp%rhs(1,jbl) - d_cor

      ! internal-ghost face:
      else
        jb1 = (jn - dmn1) + 1
        call hash_get_value (dmn%gi2i, jn, jbl)

        kij = csr_index (j1l, jb1, pp%mat%col, pp%mat%row)

        dx = xdmn%x_pc(:,i) - xdmn%x_vc%dmn(:,j1l)
        phi_pc = dot_product (pp%grad%dmn(:,1,j1l), dx)

        dx = xdmn%x_nc(:,i) - xdmn%x_vc%gdmn(:,jbl)
        phi_nc = dot_product (pp%grad%gdmn(:,1,jbl), dx)

        d_ij = pp%mat%val(kij)
        d_cor = -d_ij * (phi_nc - phi_pc)

        pp%rhs(1,j1l) = pp%rhs(1,j1l) + d_cor
      end if
    end select
  end do
end subroutine




subroutine ns_source (bnd, dmn, xbnd, xdmn, pres, vel)
  type(mesh_t), intent(in) :: bnd, dmn
  type(geom_t), intent(in) :: xbnd, xdmn
  type(pde_t), intent(inout) :: pres
  type(pde_t), intent(inout) :: vel

  real(rp), dimension(dmn%n_dim) :: dpres

  integer :: n_blk, blk, err
  integer :: i, jj1, jjn, j1, j11, j1l, jn, jb1, jbl, ph_ty
  integer :: bnd1, bndl, bndu, dmn1, dmnl, dmnu

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

!$AD II-LOOP
  do i = 1, dmn%n_elm
    j1l = i
    j11 = i + (dmnl - dmn1)

    dpres = -pres%grad%dmn(:,1,j1l)
    dpres = dpres * xdmn%vol%dmn(j1l)

    vel%rhs(:,j1l) = vel%rhs(:,j1l) + dpres
  end do
end subroutine




subroutine ns_correction (bnd, dmn, xbnd, xdmn, pp, vel, pres)
  type(mesh_t), intent(in) :: bnd, dmn
  type(geom_t), intent(in) :: xbnd, xdmn
  type(pde_t), intent(inout) :: pp
  type(pde_t), intent(inout) :: vel
  type(pde_t), intent(inout) :: pres

  real(rp), dimension(dmn%n_dim) :: dvel
  real(rp) :: dpres, urf

  integer :: n_blk, blk, err
  integer :: i, jj1, jjn, j1, j11, j1l, jn, jb1, jbl, ph_ty
  integer :: bnd1, bndl, bndu, dmn1, dmnl, dmnu

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  blk = 1; n_blk = 1;
#ifdef USE_MPI
  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)
#endif

  if (blk == 0) then
    pp%ref(1) = pp%phi%dmn(1,1)
  end if

#ifdef USE_MPI
  if (n_blk > 1) then
    call MPI_Bcast (pp%ref(1), 1, mpi_rp, &
      0, mpi_comm_world, err)
  end if
#endif

  urf = pp%urf

  do i = 1, dmn%n_elm
    j1l = i
    j11 = i + (dmnl - dmn1)

    dvel = -pp%grad%dmn(:,1,j1l)
    dvel = dvel * xdmn%vol%dmn(j1l) * vel%diag%dmn(j1l)

    dpres = pp%phi%dmn(1,j1l) - pp%ref(1)
    dpres = dpres * urf

    vel%phi%dmn(:,j1l) = vel%phi%dmn(:,j1l) + dvel
    pres%phi%dmn(1,j1l) = pres%phi%dmn(1,j1l) + dpres
  end do
end subroutine



subroutine fp_iter_ctrl (ctrl, id, iter, n_iter, k_iter, fp_exit, cutoff, norm, norm_gl)
  integer, intent(in)::id
  integer, intent(in)::iter,n_iter,k_iter
  real(rp),dimension(:), intent(inout)::norm
  real(rp), intent(in)::cutoff
  integer, intent(inout)::fp_exit
  real(rp), intent(out)::norm_gl
  character, intent(in)::ctrl

  real(rp),save::norm_mx
  integer::i


  last_fp_iter=.false.

  if((ctrl=='F'.and.iter==1) .or. (ctrl=='R'.and.iter==n_iter))then
    fp_exit=0
    norm_mx=zero
  end if


  norm_gl=maxval(norm)
  norm=zero
  norm_mx=max(norm_mx,norm_gl)
  norm_gl=norm_gl/(norm_mx+small)

  if(norm_gl<cutoff .or. (ctrl=='F'.and.iter>n_iter-2))then
    fp_exit=fp_exit+1
    if(fp_exit==1)last_fp_iter=.true.
  end if


  if(k_iter>=0 .and. mod(iter,k_iter)==0)then
    print "(a,i1,a,i6,a,es12.4)",ctrl,id," iteration:",iter,",  Norm:",norm_gl
  end if
end subroutine
end module
