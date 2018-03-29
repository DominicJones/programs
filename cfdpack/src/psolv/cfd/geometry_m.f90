module geometry_m
use constants_m
use hash_m
use graph_m
use vector_utils_m
implicit none

type mshf_t
  integer :: n_elm = 0
  integer,dimension(:),pointer::elm_idx => null()
  integer,dimension(:),pointer::elm_ty => null()
  integer,dimension(:),pointer::elm_dim => null()
  type(graph_t)::fce_vrt
end type

type mesh_t
  integer::blk = 0
  integer::n_blk = 0
  integer::n_dim = 0
  integer::n_vrt = 0
  integer::n_elm = 0
  integer::n_ghst = 0
  integer::n_locl = 0
  integer::n_fce = 0
  integer,dimension(:),pointer::vrt_idx => null()
  integer,dimension(:),pointer::elm_blk => null()
  integer,dimension(:),pointer::elm_idx => null()
  integer,dimension(:),pointer::ghst_blk => null()
  integer,dimension(:),pointer::ghst_idx => null()
  integer,dimension(:),pointer::locl_blk => null()
  integer,dimension(:),pointer::locl_idx => null()
  type(graph_t)::elm_vrt, elm_tag, elm_neig
  type(graph_t)::ghst_vrt, ghst_tag
  type(graph_t)::fce_vrt, fce_elm
  type(hash_t),dimension(:),allocatable::vi2i
  type(hash_t),dimension(:),allocatable::gi2i, gi2b
end type

type sfld_t
  integer :: id = 0
  real(rp),dimension(:),pointer::bnd => null()
  real(rp),dimension(:),pointer::gbnd => null()
  real(rp),dimension(:),pointer::lbnd => null()
  real(rp),dimension(:),pointer::dmn => null()
  real(rp),dimension(:),pointer::gdmn => null()
  real(rp),dimension(:),pointer::ldmn => null()
end type

type vfld_t
  integer :: id = 0
  real(rp),dimension(:,:),pointer::bnd => null()
  real(rp),dimension(:,:),pointer::gbnd => null()
  real(rp),dimension(:,:),pointer::lbnd => null()
  real(rp),dimension(:,:),pointer::dmn => null()
  real(rp),dimension(:,:),pointer::gdmn => null()
  real(rp),dimension(:,:),pointer::ldmn => null()
end type

type tfld_t
  integer :: id = 0
  real(rp),dimension(:,:,:),pointer::bnd => null()
  real(rp),dimension(:,:,:),pointer::gbnd => null()
  real(rp),dimension(:,:,:),pointer::lbnd => null()
  real(rp),dimension(:,:,:),pointer::dmn => null()
  real(rp),dimension(:,:,:),pointer::gdmn => null()
  real(rp),dimension(:,:,:),pointer::ldmn => null()
end type

type geom_t
  real(rp)::vol_min = 0
  real(rp)::vol_max = 0
  real(rp)::vol_sum = 0
  type(sfld_t)::vol
  type(vfld_t)::x_vc, gl_min, gl_max, gl
  type(tfld_t)::ls, gg
  real(rp),dimension(:,:),pointer::x_vrt => null()
  real(rp),dimension(:,:),pointer::x_fc => null()
  real(rp),dimension(:,:),pointer::x_pc => null()
  real(rp),dimension(:,:),pointer::x_nc => null()
  real(rp),dimension(:,:),pointer::norm => null()
  real(rp),dimension(:),pointer::l_pn => null()
  real(rp),dimension(:),pointer::w_fp => null()
  real(rp),dimension(:),pointer::area => null()
  real(rp),dimension(:),pointer::volr => null()
end type

contains

subroutine geometry (bnd, dmn, xbnd, xdmn)
  type(mesh_t),intent(in)::bnd, dmn
  type(geom_t),intent(inout)::xbnd, xdmn

  integer::n_dim,bc,bnd1,bndl,bndu,dmn1,dmnl,dmnu
  integer::i,jj1,jjn,jj,j,j1,jn,jb,jv,kk1,kkn,skk,kk,k,k1,k2,kv,k1v,l
  integer::j11,j1l,jb1,jbl
  real(rp)::fvol,area1,l_pn1,area_sum,vol_sum,bc_area(100)
  real(rp),dimension(3)::edge1,edge,cprod,unit_z
  real(rp),dimension(dmn%n_dim)::dx,x1,x2,x_fc1,x_pc1,x_nc1,norm1


  n_dim = dmn%n_dim
  unit_z = (/zero, zero, one/)

  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)


  do i = 1, dmn%n_elm
    jj1 = dmn%elm_vrt%row(i)
    jjn = dmn%elm_vrt%row(i+1) - 1

    x1 = zero
    do jj = jj1, jjn
      j = dmn%elm_vrt%col(jj)
      call hash_get_value (dmn%vi2i, j, jv)
      x1 = x1 + xdmn%x_vrt(:,jv)
    end do

    x1 = x1 / (jjn - jj1 + 1)
    xdmn%x_vc%dmn(:,i) = x1
    xdmn%vol%dmn(i) = zero
  end do


  call ghost_update_v (bnd, dmn, xdmn%x_vc, 0, 1)


  bc_area = zero
  area_sum = zero
  vol_sum = zero

  l = 0
  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j11 = (j1 - dmn1) + 1
    j1l = (j1 - dmnl) + 1
    j1 = (j1 - dmnl) + 1

    jn = dmn%fce_elm%col(jjn)

    kk1 = dmn%fce_vrt%row(i)
    kkn = dmn%fce_vrt%row(i+1) - 1


    ! face centre
    x1 = zero
    do kk = kk1, kkn
      k = dmn%fce_vrt%col(kk)
      call hash_get_value (dmn%vi2i, k, kv)
      x1 = x1 + xdmn%x_vrt(:,kv)
    end do

    x1 = x1 / (kkn - kk1 + 1)
    x_fc1 = x1
    xdmn%x_fc(:,i) = x1


    ! face normal and area
    k1 = dmn%fce_vrt%col(kk1)
    call hash_get_value (dmn%vi2i, k1, k1v)
    x1 = zero

    do kk = kk1 + 1, kkn
      k = dmn%fce_vrt%col(kk)
      call hash_get_value (dmn%vi2i, k, kv)
      edge = zero
      edge(:n_dim) = xdmn%x_vrt(:,kv) - xdmn%x_vrt(:,k1v)

      select case(kkn - kk1 + 1)
      case(2)
        cprod = cross_prod (edge, unit_z)
        x1 = cprod (:n_dim)
      case(3:)
        select case(kk - (kk1 + 1))
        case(1:)
          cprod = cross_prod (edge1, edge)
          x1 = x1 + half * cprod
        end select
        edge1 = edge
      end select
    end do

    ! select case (n_dim)
    ! case (2)
    !   k1 = dmn%fce_vrt%col(kk1)
    !   call hash_get_value (dmn%vi2i, k1, k1v)

    !   k = dmn%fce_vrt%col(kk1 + 1)
    !   call hash_get_value (dmn%vi2i, k, kv)

    !   edge = zero
    !   edge(:n_dim) = xdmn%x_vrt(:,kv) - xdmn%x_vrt(:,k1v)

    !   cprod = cross_prod (edge, unit_z)
    !   x1 = cprod(:n_dim)

    ! case (3)
    !   k1 = dmn%fce_vrt%col(kk1)
    !   call hash_get_value (dmn%vi2i, k1, k1v)
    !   edge1 = xdmn%x_vrt(:,k1v) - x_fc1

    !   x1 = zero
    !   do kk = kk1, kkn
    !     k = k1
    !     if (kk < kkn) k = dmn%fce_vrt%col(kk + 1)
    !     call hash_get_value (dmn%vi2i, k, kv)

    !     edge = xdmn%x_vrt(:,kv) - x_fc1
    !     cprod = cross_prod (edge1, edge)
    !     x1 = x1 + half * cprod

    !     edge1 = edge
    !   end do
    ! end select

    area1 = vec_mag (x1)
    norm1 = x1 / area1
    xdmn%area(i) = area1
    xdmn%norm(:,i) = norm1


    ! volume
    fvol = x_fc1(1) * norm1(1) * area1
    xdmn%vol%dmn(j1) = xdmn%vol%dmn(j1) + fvol


    ! aux. pole position
    dx = x_fc1 - xdmn%x_vc%dmn(:,j1)
    x_pc1 = x_fc1 - norm1 * dot_product (dx,norm1)
    xdmn%x_pc(:,i) = x_pc1


    ! internal and internal-ghost faces:
    select case(jn)
    case(1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb = (jn - dmnl) + 1

        ! volume
        xdmn%vol%dmn(jb) = xdmn%vol%dmn(jb) - fvol


        ! aux. neig. position
        dx = x_fc1 - xdmn%x_vc%dmn(:,jb)
        x_nc1 = x_fc1 - norm1 * dot_product (dx,norm1)
        xdmn%x_nc(:,i) = x_nc1


        ! pole to neig. length
        dx = x_nc1 - x_pc1
        l_pn1 = vec_mag (dx)


        ! interpolation factor
        dx = x_fc1 - x_pc1
        xdmn%w_fp(i) = vec_mag(dx) / l_pn1


      ! internal-ghost face:
      else
        jb1 = (jn - dmn1) + 1
        call hash_get_value (dmn%gi2i, jn, jbl)

        call hash_get_value (dmn%gi2i, jn, jb)


        ! aux. neig. position
        dx = x_fc1 - xdmn%x_vc%gdmn(:,jb)
        x_nc1 = x_fc1 - norm1 * dot_product (dx,norm1)
        xdmn%x_nc(:,i) = x_nc1


        ! pole to neig. length
        dx = x_nc1 - x_pc1
        l_pn1 = vec_mag (dx)


        ! interpolation factor
!         dx = x_fc1 - x_pc1
!         xdmn%w_fp(i) = vec_mag(dx) / l_pn1
        xdmn%w_fp(i) = half
      end if


    ! boundary faces:
    case default
      l = l + 1
      jb = ((-1)*jn - bndl) + 1


      ! checks
      bc = bnd%elm_tag%row(jb)
      bc = bnd%elm_tag%col(bc+1)
      bc_area(bc) = bc_area(bc) + area1

      area_sum = area_sum + area1
      vol_sum = vol_sum + fvol


      ! pole to face length
      dx = x_fc1 - x_pc1
      l_pn1 = vec_mag (dx)
    end select


    xdmn%l_pn(i) = l_pn1
  end do


  call ghost_update_s (bnd, dmn, xdmn%vol, 0, 1)


  xdmn%vol_min = large
  xdmn%vol_max = zero
  xdmn%vol_sum = zero

  do i = 1, dmn%n_elm
    xdmn%volr(i) = 1 / (xdmn%vol%dmn(i) + small)
    xdmn%vol_min = min (xdmn%vol_min, xdmn%vol%dmn(i))
    xdmn%vol_max = max (xdmn%vol_max, xdmn%vol%dmn(i))
    xdmn%vol_sum = xdmn%vol_sum + xdmn%vol%dmn(i)
  end do

  if(xdmn%vol_min < small)then
    print "(a)","invalid volume(s) found"
  end if
end subroutine




subroutine gl_local_extrema (bnd, dmn, xbnd, xdmn, phi)
  use vector_utils_m
  use matrix_utils_m

  type(mesh_t),intent(in)::bnd, dmn
  type(geom_t),intent(inout)::xbnd, xdmn
  type(vfld_t),intent(inout)::phi

  integer::bnd1,bndl,bndu,dmn1,dmnl,dmnu
  integer::i,jj1,jjn,j1,jn,ic,jb,jc, n_cmp
  real(rp),dimension(dmn%n_dim)::dx
  real(rp),dimension(dmn%n_dim,dmn%n_dim)::mat

  n_cmp = size (phi%dmn, 1)

  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  do i = 1, dmn%n_elm
    do ic = 1, n_cmp
      xdmn%gl_min%dmn(ic,i) = phi%dmn(ic,i)
      xdmn%gl_max%dmn(ic,i) = phi%dmn(ic,i)
    end do
  end do


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j1 = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)


    ! internal and internal-ghost faces:
    select case(jn)
    case(1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb = (jn - dmnl) + 1

        do ic = 1, n_cmp
          xdmn%gl_min%dmn(ic,j1) = min (xdmn%gl_min%dmn(ic,j1), phi%dmn(ic,jb))
          xdmn%gl_min%dmn(ic,jb) = min (xdmn%gl_min%dmn(ic,jb), phi%dmn(ic,j1))

          xdmn%gl_max%dmn(ic,j1) = max (xdmn%gl_max%dmn(ic,j1), phi%dmn(ic,jb))
          xdmn%gl_max%dmn(ic,jb) = max (xdmn%gl_max%dmn(ic,jb), phi%dmn(ic,j1))
        end do


      ! internal-ghost face:
      else
        call hash_get_value (dmn%gi2i, jn, jb)

        do ic = 1, n_cmp
          xdmn%gl_min%dmn(ic,j1) = min (xdmn%gl_min%dmn(ic,j1), phi%gdmn(ic,jb))
          xdmn%gl_max%dmn(ic,j1) = max (xdmn%gl_max%dmn(ic,j1), phi%gdmn(ic,jb))
        end do
      end if

    end select
  end do
end subroutine




subroutine grad_limiter (bnd, dmn, xbnd, xdmn, phi, grad, ctrl)
  use vector_utils_m
  use matrix_utils_m

  type(mesh_t),intent(in)::bnd, dmn
  type(geom_t),intent(inout)::xbnd, xdmn
  type(vfld_t),intent(inout)::phi
  type(tfld_t),intent(inout)::grad
  integer, intent(in) :: ctrl

  integer::bnd1,bndl,bndu,dmn1,dmnl,dmnu
  integer::i,jj1,jjn,j1,jn,ic,jb,jc, n_cmp
  real(rp),dimension(dmn%n_dim)::dx
  real(rp),dimension(dmn%n_dim,dmn%n_dim)::mat

  real(rp)::fac, fac3, pwr
  real(rp),dimension(2)::eps2
  real(rp),dimension(dmn%n_dim,2)::d1min, d1max, d2, limval

  fac = 5.
  fac3 = fac**3
  pwr = 3./dmn%n_dim

  n_cmp = size (phi%dmn, 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  xdmn%gl%dmn = one


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j1 = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)


    ! internal and internal-ghost faces:
    select case(jn)
    case(1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb = (jn - dmnl) + 1

        eps2(1) = fac3*xdmn%vol%dmn(j1)**pwr
        eps2(2) = fac3*xdmn%vol%dmn(jb)**pwr

        do ic = 1, n_cmp
          d1min(ic,1) = xdmn%gl_min%dmn(ic,j1) - phi%dmn(ic,j1)
          d1min(ic,2) = xdmn%gl_min%dmn(ic,jb) - phi%dmn(ic,jb)

          d1max(ic,1) = xdmn%gl_max%dmn(ic,j1) - phi%dmn(ic,j1)
          d1max(ic,2) = xdmn%gl_max%dmn(ic,jb) - phi%dmn(ic,jb)
        end do

        dx = xdmn%x_fc(:,i) - xdmn%x_vc%dmn(:,j1)
        do ic = 1, n_cmp
          d2(ic,1) = sum(grad%dmn(:,ic,j1)*dx)
          d2(ic,2) = sum(grad%dmn(:,ic,jb)*(-dx))
        end do

        select case(ctrl)
        case (1)
          do ic = 1, n_cmp
            limval(ic,1) = barth (d2(ic,1),d1min(ic,1),d1max(ic,1))
            limval(ic,2) = barth (d2(ic,2),d1min(ic,2),d1max(ic,2))
          end do
        case (2)
          do ic = 1, n_cmp
            limval(ic,1) = venkat (d2(ic,1),d1min(ic,1),d1max(ic,1),eps2(1))
            limval(ic,2) = venkat (d2(ic,2),d1min(ic,2),d1max(ic,2),eps2(2))
          end do
        end select

        do ic = 1, n_cmp
          xdmn%gl%dmn(ic,j1) = min(xdmn%gl%dmn(ic,j1), limval(ic,1))
          xdmn%gl%dmn(ic,jb) = min(xdmn%gl%dmn(ic,jb), limval(ic,2))
        end do


      ! internal-ghost face:
      else
        call hash_get_value (dmn%gi2i, jn, jb)


        eps2(1) = fac3*xdmn%vol%dmn(j1)**pwr
        eps2(2) = fac3*xdmn%vol%dmn(jb)**pwr

        do ic = 1, n_cmp
          d1min(ic,1) = xdmn%gl_min%dmn(ic,j1) - phi%dmn(ic,j1)
          d1max(ic,1) = xdmn%gl_max%dmn(ic,j1) - phi%dmn(ic,j1)
        end do

        dx = xdmn%x_fc(:,i) - xdmn%x_vc%dmn(:,j1)
        do ic = 1, n_cmp
          d2(ic,1) = sum(grad%dmn(:,ic,j1)*dx)
        end do

        select case(ctrl)
        case (1)
          do ic = 1, n_cmp
            limval(ic,1) = barth (d2(ic,1),d1min(ic,1),d1max(ic,1))
          end do
        case (2)
          do ic = 1, n_cmp
            limval(ic,1) = venkat (d2(ic,1),d1min(ic,1),d1max(ic,1),eps2(1))
          end do
        end select

        do ic = 1, n_cmp
          xdmn%gl%dmn(ic,j1) = min(xdmn%gl%dmn(ic,j1), limval(ic,1))
        end do

      end if
    end select
  end do


  do i = 1, dmn%n_elm
    do ic = 1, n_cmp
      grad%dmn(:,ic,i) = grad%dmn(:,ic,i) * xdmn%gl%dmn(ic,i)
    end do
  end do
contains


function barth(d2,d1min,d1max) result(res)
  real(rp):: d2, d1min, d1max

  real(rp):: num, den, res

  if (d2 > small) then
    num    = d1max
    den    = d2
    res = min(one,num/(den+small))
  else if (d2 < -small) then
    num    = d1min
    den    = d2
    res = min(one,num/(den+small))
  else
    res = one
  endif
end function


function venkat(d2,d1min,d1max,eps2) result(res)
  real(rp):: d2, d1min, d1max, eps2

  real(rp):: num, den, res

  if (d2 > small) then
    num    = (d1max*d1max+eps2)*d2 + 2*d2*d2*d1max
    den    = d2*(d1max*d1max+2*d2*d2+d1max*d2+eps2)
    res = num/(den+small)
  else if (d2 < -small) then
    num    = (d1min*d1min+eps2)*d2 + 2*d2*d2*d1min
    den    = d2*(d1min*d1min+2*d2*d2+d1min*d2+eps2)
    res = num/(den+small)
  else
    res = one
  endif
end function
end subroutine




subroutine ls_grad_matrix (bnd, dmn, xbnd, xdmn)
  use vector_utils_m
  use matrix_utils_m

  type(mesh_t),intent(in)::bnd, dmn
  type(geom_t),intent(inout)::xbnd, xdmn

  integer::bnd1,bndl,bndu,dmn1,dmnl,dmnu
  integer::i,jj1,jjn,j1,jn,ic,jb,jc
  real(rp),dimension(dmn%n_dim)::dx
  real(rp),dimension(dmn%n_dim,dmn%n_dim)::mat

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  xdmn%ls%dmn = zero


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j1 = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)


    ! internal and internal-ghost faces:
    select case(jn)
    case(1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb = (jn - dmnl) + 1

        dx = xdmn%x_vc%dmn(:,jb) - xdmn%x_vc%dmn(:,j1)

        do jc = 1, dmn%n_dim
          do ic = 1, dmn%n_dim
            xdmn%ls%dmn(ic,jc,j1) = xdmn%ls%dmn(ic,jc,j1) + dx(ic) * dx(jc)
            xdmn%ls%dmn(ic,jc,jb) = xdmn%ls%dmn(ic,jc,jb) + dx(ic) * dx(jc)
          end do
        end do


      ! internal-ghost face:
      else
        call hash_get_value (dmn%gi2i, jn, jb)

        dx = xdmn%x_vc%gdmn(:,jb) - xdmn%x_vc%dmn(:,j1)

        do jc = 1, dmn%n_dim
          do ic = 1, dmn%n_dim
            xdmn%ls%dmn(ic,jc,j1) = xdmn%ls%dmn(ic,jc,j1) + dx(ic) * dx(jc)
          end do
        end do
      end if


    ! boundary faces:
    case default
      dx = xdmn%x_fc(:,i) - xdmn%x_vc%dmn(:,j1)

      do jc = 1, dmn%n_dim
        do ic = 1, dmn%n_dim
          xdmn%ls%dmn(ic,jc,j1) = xdmn%ls%dmn(ic,jc,j1) + dx(ic) * dx(jc)
        end do
      end do
    end select
  end do


! UPDATE FROM GHOST VALUES


  do i = 1, dmn%n_elm
    mat = xdmn%ls%dmn(:,:,i)
    call matrix_inv (mat, xdmn%ls%dmn(:,:,i))
  end do
end subroutine




subroutine ls_gradient (bnd, dmn, xbnd, xdmn, phi, grad, neu_bnd)
  use vector_utils_m
  use matrix_utils_m

  type(mesh_t),intent(in)::bnd, dmn
  type(geom_t),intent(inout)::xbnd, xdmn

  type(vfld_t),intent(inout)::phi
  type(tfld_t),intent(inout)::grad
  integer,intent(in)::neu_bnd

  real(rp),dimension(size(phi%dmn,1))::dphi
  real(rp),dimension(dmn%n_dim)::dx
  integer::i,jj1,jjn,j1,jn,jb,ic,jc,n_cmp,n_dmn
  integer::bnd1,bndl,bndu,dmn1,dmnl,dmnu


  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  n_cmp = size(phi%dmn,1)
  grad%dmn = zero


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j1 = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)


    ! internal and internal-ghost faces:
    select case(jn)
    case(1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb = (jn - dmnl) + 1

        dx = xdmn%x_vc%dmn(:,jb) - xdmn%x_vc%dmn(:,j1)
        dphi = phi%dmn(:,jb) - phi%dmn(:,j1)

        do ic=1, n_cmp
          grad%dmn(:,ic,j1) = grad%dmn(:,ic,j1) + dx * dphi(ic)
          grad%dmn(:,ic,jb) = grad%dmn(:,ic,jb) + dx * dphi(ic)
        end do


      ! internal-ghost face:
      else
        call hash_get_value (dmn%gi2i, jn, jb)

        dx = xdmn%x_vc%gdmn(:,jb) - xdmn%x_vc%dmn(:,j1)
        dphi = phi%gdmn(:,jb) - phi%dmn(:,j1)

        do ic=1, n_cmp
          grad%dmn(:,ic,j1) = grad%dmn(:,ic,j1) + dx * dphi(ic)
        end do
      end if


    case default
      jb = ((-1)*jn - bndl) + 1

      dx = xdmn%x_fc(:,i) - xdmn%x_vc%dmn(:,j1)
      dphi = phi%bnd(:,jb) - phi%dmn(:,j1)
      if (neu_bnd /= 0) dphi = zero

      do ic=1, n_cmp
        grad%dmn(:,ic,j1) = grad%dmn(:,ic,j1) + dx * dphi(ic)
      end do
    end select
  end do


  do i=1,dmn%n_elm
    do ic=1, n_cmp
      grad%dmn(:,ic,i) = matmul (xdmn%ls%dmn(:,:,i), grad%dmn(:,ic,i))
    end do
  end do


  if(neu_bnd /= 0)then
    do i = 1, dmn%n_fce
      jj1 = 2 * (i - 1) + 1
      jjn = 2 * (i)

      j1 = dmn%fce_elm%col(jj1)
      j1 = (j1 - dmnl) + 1
      jn = dmn%fce_elm%col(jjn)


      ! internal and internal-ghost faces:
      select case(jn)
      case(:-1)
        jb = ((-1)*jn - bndl) + 1

        dx = xdmn%x_fc(:,i) - xdmn%x_vc%dmn(:,j1)

        do ic=1, n_cmp
          phi%bnd(ic,jb) = phi%dmn(ic,j1) + dot_product (grad%dmn(:,ic,j1), dx)
        end do
      end select
    end do
  end if
end subroutine



subroutine gg_gradient (bnd, dmn, xbnd, xdmn, phi, grad, neu_bnd, n_iter)
  use mpi
  use vector_utils_m
  use matrix_utils_m

  type(mesh_t),intent(in)::bnd, dmn
  type(geom_t),intent(inout)::xbnd, xdmn

  type(vfld_t),intent(inout)::phi
  type(tfld_t),intent(inout)::grad
  integer,intent(in)::neu_bnd
  integer,intent(in)::n_iter

  real(rp),dimension(dmn%n_dim, size(phi%dmn,1))::fgrad
  real(rp),dimension(size(phi%dmn,1))::vphi, nphi, fphi
  real(rp),dimension(dmn%n_dim)::vdx, ndx, dx
  real(rp)::w_fp1,w_fn1
  integer::n_blk, err
  integer::i,j,jj1,jjn,j1,jn,jb,ic,jc,n_cmp,n_dmn
  integer::bnd1,bndl,bndu,dmn1,dmnl,dmnu
  integer::iter, m_iter

  call MPI_Comm_size (mpi_comm_world, n_blk, err)

  bnd1 = bnd%elm_blk(1)
  bndl = bnd%elm_blk(bnd%blk)
  bndu = bnd%elm_blk(bnd%blk + 1)

  dmn1 = dmn%elm_blk(1)
  dmnl = dmn%elm_blk(dmn%blk)
  dmnu = dmn%elm_blk(dmn%blk + 1)

  n_cmp = size(phi%dmn,1)
  xdmn%gg%dmn = zero
  xdmn%gg%gdmn = zero

  m_iter = n_iter
  m_iter = max (m_iter, 1)
  m_iter = min (m_iter, 3)


  do iter = 1, m_iter
  grad%dmn = zero


  if(neu_bnd /= 0)then
    do i = 1, dmn%n_fce
      jj1 = 2 * (i - 1) + 1
      jjn = 2 * (i)

      j1 = dmn%fce_elm%col(jj1)
      j1 = (j1 - dmnl) + 1
      jn = dmn%fce_elm%col(jjn)

      ! boundary faces:
      select case(jn)
      case(:-1)
        jb = ((-1)*jn - bndl) + 1
        vdx = xdmn%x_fc(:,i) - xdmn%x_vc%dmn(:,j1)
        do ic = 1, n_cmp
          phi%bnd(ic,jb) = phi%dmn(ic,j1) + dot_product (xdmn%gg%dmn(:,ic,j1), vdx)
        end do
      end select
    end do
  end if


  do i = 1, dmn%n_fce
    jj1 = 2 * (i - 1) + 1
    jjn = 2 * (i)

    j1 = dmn%fce_elm%col(jj1)
    j1 = (j1 - dmnl) + 1
    jn = dmn%fce_elm%col(jjn)


    ! internal and internal-ghost faces:
    select case(jn)
    case(1:)


      ! internal face:
      if (jn >= dmnl .and. jn < dmnu) then
        jb = (jn - dmnl) + 1

        w_fp1 = xdmn%w_fp(i)
        w_fn1 = one - w_fp1
        do ic = 1, n_cmp
          fgrad(:,ic) = w_fp1 * xdmn%gg%dmn(:,ic,jb) + w_fn1 * xdmn%gg%dmn(:,ic,j1)
        end do
        dx = w_fp1 * xdmn%x_vc%dmn(:,jb) + w_fn1 * xdmn%x_vc%dmn(:,j1)
        dx = xdmn%x_fc(:,i) - dx
        do ic = 1, n_cmp
          fphi(ic) = dot_product (fgrad(:,ic), dx)
        end do
        fphi = w_fp1 * phi%dmn(:,jb) + w_fn1 * phi%dmn(:,j1) + fphi

        vdx = xdmn%area(i) * xdmn%norm(:,i) * xdmn%volr(j1)
        ndx = xdmn%area(i) * xdmn%norm(:,i) * xdmn%volr(jb)

        do ic = 1, n_cmp
          grad%dmn(:,ic,j1) = grad%dmn(:,ic,j1) + vdx * fphi(ic)
          grad%dmn(:,ic,jb) = grad%dmn(:,ic,jb) - ndx * fphi(ic)
        end do

      ! internal-ghost face:
      else
        call hash_get_value (dmn%gi2i, jn, jb)

        w_fp1 = half
        w_fn1 = one - w_fp1
        do ic = 1, n_cmp
          fgrad(:,ic) = w_fp1 * xdmn%gg%gdmn(:,ic,jb) + w_fn1 * xdmn%gg%dmn(:,ic,j1)
        end do
        dx = w_fp1 * xdmn%x_vc%gdmn(:,jb) + w_fn1 * xdmn%x_vc%dmn(:,j1)
        dx = xdmn%x_fc(:,i) - dx
        do ic = 1, n_cmp
          fphi(ic) = dot_product (fgrad(:,ic), dx)
        end do
        fphi = w_fp1 * phi%gdmn(:,jb) + w_fn1 * phi%dmn(:,j1) + fphi




        vdx = xdmn%area(i) * xdmn%norm(:,i) * xdmn%volr(j1)

        do ic = 1, n_cmp
          grad%dmn(:,ic,j1) = grad%dmn(:,ic,j1) + vdx * fphi(ic)
        end do
      end if


    case default
      jb = ((-1)*jn - bndl) + 1

      vdx = xdmn%area(i) * xdmn%norm(:,i) * xdmn%volr(j1)
      fphi = phi%bnd(:,jb)

      do ic = 1, n_cmp
        grad%dmn(:,ic,j1) = grad%dmn(:,ic,j1) + vdx * fphi(ic)
      end do
    end select
  end do


  if(iter /= n_iter)then
    xdmn%gg%dmn(:,:n_cmp,:) = grad%dmn
    if (n_blk > 1) call ghost_update_t (bnd, dmn, xdmn%gg, 0, 1)
  end if
  end do
end subroutine




subroutine init_ghost_update (bnd, dmn)
  use mpi
  use hash_m

  type(mesh_t),intent(inout)::bnd, dmn

  integer :: err = 0
  integer :: i, ii, j, jj1, jjn, k, m, n
  integer :: blk, n_blk, oblk
  integer, dimension(:), allocatable :: mpi_rqst
  integer, dimension(:,:), allocatable :: mpi_stat


  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)

  allocate (mpi_rqst(n_blk)); mpi_rqst = 0
  allocate (mpi_stat(mpi_status_size, n_blk)); mpi_stat = 0


  !! setup ghost element index block offsets (ascending ghst_idx assumed)
  allocate (bnd%ghst_blk(n_blk + 1))
  bnd%ghst_blk = 0

  do i = 1, bnd%n_ghst
    call hash_get_value (bnd%gi2b, bnd%ghst_idx(i), oblk)
    bnd%ghst_blk(oblk + 1) = bnd%ghst_blk(oblk + 1) + 1
  end do

  bnd%ghst_blk(1) = 1
  do i = 2, n_blk + 1
    bnd%ghst_blk(i) = bnd%ghst_blk(i) + bnd%ghst_blk(i - 1)
  end do


  allocate (dmn%ghst_blk(n_blk + 1))
  dmn%ghst_blk = 0

  do i = 1, dmn%n_ghst
    call hash_get_value (dmn%gi2b, dmn%ghst_idx(i), oblk)
    dmn%ghst_blk(oblk + 1) = dmn%ghst_blk(oblk + 1) + 1
  end do

  dmn%ghst_blk(1) = 1
  do i = 2, n_blk + 1
    dmn%ghst_blk(i) = dmn%ghst_blk(i) + dmn%ghst_blk(i - 1)
  end do


  !! determine the ghosts possessed by this [partition]
  allocate (bnd%locl_blk(n_blk + 1))
  bnd%locl_blk = 0

  do i = 0, n_blk - 1
    if (i == blk) cycle
    call MPI_Irecv (bnd%locl_blk(i+2), 1, MPI_INTEGER, i, &
      mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
  end do

  do i = 0, n_blk - 1
    if (i == blk) cycle
    n = bnd%ghst_blk(i + 2) - bnd%ghst_blk(i + 1)
    call MPI_Send (n, 1, MPI_INTEGER, i, &
      100, mpi_comm_world, err)
  end do

  call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)

  bnd%locl_blk(1) = 1
  do i = 2, n_blk + 1
    bnd%locl_blk(i) = bnd%locl_blk(i) + bnd%locl_blk(i - 1)
  end do


  allocate (dmn%locl_blk(n_blk + 1))
  dmn%locl_blk = 0

  do i = 0, n_blk - 1
    if (i == blk) cycle
    call MPI_Irecv (dmn%locl_blk(i+2), 1, MPI_INTEGER, i, &
      mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
  end do

  do i = 0, n_blk - 1
    if (i == blk) cycle
    n = dmn%ghst_blk(i + 2) - dmn%ghst_blk(i + 1)
    call MPI_Send (n, 1, MPI_INTEGER, i, &
      100, mpi_comm_world, err)
  end do

  call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)

  dmn%locl_blk(1) = 1
  do i = 2, n_blk + 1
    dmn%locl_blk(i) = dmn%locl_blk(i) + dmn%locl_blk(i - 1)
  end do


  !! which elements are they (what are their indices)?
  n = bnd%locl_blk(n_blk + 1) - 1
  bnd%n_locl = n

  allocate (bnd%locl_idx(n))
  bnd%locl_idx = 0

  do i = 0, n_blk - 1
    jj1 = bnd%locl_blk(i + 1)
    jjn = bnd%locl_blk(i + 2) - 1
    n = jjn - jj1 + 1
    if (i == blk .or. n == 0) cycle
    call MPI_Irecv (bnd%locl_idx(jj1:jjn), n, MPI_INTEGER, i, &
      mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
  end do

  do i = 0, n_blk - 1
    jj1 = bnd%ghst_blk(i + 1)
    jjn = bnd%ghst_blk(i + 2) - 1
    n = jjn - jj1 + 1
    if (i == blk .or. n == 0) cycle
    call MPI_Send (bnd%ghst_idx(jj1:jjn), n, MPI_INTEGER, i, &
      100, mpi_comm_world, err)
  end do

  call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)

  n = dmn%locl_blk(n_blk + 1) - 1
  dmn%n_locl = n

  allocate (dmn%locl_idx(n))
  dmn%locl_idx = 0

  do i = 0, n_blk - 1
    jj1 = dmn%locl_blk(i + 1)
    jjn = dmn%locl_blk(i + 2) - 1
    n = jjn - jj1 + 1
    if (i == blk .or. n == 0) cycle
    call MPI_Irecv (dmn%locl_idx(jj1:jjn), n, MPI_INTEGER, i, &
      mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
  end do

  do i = 0, n_blk - 1
    jj1 = dmn%ghst_blk(i + 1)
    jjn = dmn%ghst_blk(i + 2) - 1
    n = jjn - jj1 + 1
    if (i == blk .or. n == 0) cycle
    call MPI_Send (dmn%ghst_idx(jj1:jjn), n, MPI_INTEGER, i, &
      100, mpi_comm_world, err)
  end do

  call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)

!   call print_graph ("boundary ghost indices", &
!      n_blk, bnd%ghst_blk, bnd%ghst_idx)
!
!   call print_graph ("boundary local indices", &
!      n_blk, bnd%locl_blk, bnd%locl_idx)
end subroutine




subroutine ghost_update_s (bnd, dmn, phi, ubnd, udmn)
  use mpi
  use hash_m

  type(mesh_t),intent(in)::bnd, dmn
  type(sfld_t),intent(inout)::phi
  integer,intent(in)::ubnd, udmn

  integer :: err = 0
  integer :: i, ii, j, jj1, jjn, k, m, n
  integer :: blk, n_blk, oblk
  integer, dimension(:), allocatable :: mpi_rqst
  integer, dimension(:,:), allocatable :: mpi_stat


  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)

  allocate (mpi_rqst(n_blk)); mpi_rqst = 0
  allocate (mpi_stat(mpi_status_size, n_blk)); mpi_stat = 0


  !! send the requested elmement values
  if (ubnd /= 0) then
    n = bnd%locl_blk(n_blk + 1) - 1
    do i = 1, n
      j = bnd%locl_idx(i)
      k = j - bnd%elm_blk(blk+1) + 1
      phi%lbnd(i) = phi%bnd(k)
    end do


    do i = 0, n_blk - 1
      jj1 = bnd%ghst_blk(i + 1)
      jjn = bnd%ghst_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Irecv (phi%gbnd(jj1:jjn), n, mpi_rp, i, &
        mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
    end do

    do i = 0, n_blk - 1
      jj1 = bnd%locl_blk(i + 1)
      jjn = bnd%locl_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Send (phi%lbnd(jj1:jjn), n, mpi_rp, i, &
        100, mpi_comm_world, err)
    end do

    call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)
  end if


  if (udmn /= 0) then
    n = dmn%locl_blk(n_blk + 1) - 1
    do i = 1, n
      j = dmn%locl_idx(i)
      k = j - dmn%elm_blk(blk+1) + 1
      phi%ldmn(i) = phi%dmn(k)
    end do


    do i = 0, n_blk - 1
      jj1 = dmn%ghst_blk(i + 1)
      jjn = dmn%ghst_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Irecv (phi%gdmn(jj1:jjn), n, mpi_rp, i, &
        mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
    end do

    do i = 0, n_blk - 1
      jj1 = dmn%locl_blk(i + 1)
      jjn = dmn%locl_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Send (phi%ldmn(jj1:jjn), n, mpi_rp, i, &
        100, mpi_comm_world, err)
    end do

    call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)
  end if


!   do i = 1, size(phi%gbnd)
!     j = bnd%ghst_idx(i)
!     call hash_get_value (bnd%gi2b, j, k)
!     print *, blk+1, "phi_ghst",i,j,k,phi%gbnd(i)
!   end do
end subroutine




subroutine ghost_update_v (bnd, dmn, phi, ubnd, udmn)
  use mpi
  use hash_m

  type(mesh_t),intent(in)::bnd, dmn
  type(vfld_t),intent(inout)::phi
  integer,intent(in)::ubnd, udmn

  integer :: err = 0
  integer :: i, ii, j, jj1, jjn, k, m, n
  integer :: n1, n2
  integer :: blk, n_blk, oblk
  integer, dimension(:), allocatable :: mpi_rqst
  integer, dimension(:,:), allocatable :: mpi_stat


  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)

  allocate (mpi_rqst(n_blk)); mpi_rqst = 0
  allocate (mpi_stat(mpi_status_size, n_blk)); mpi_stat = 0

  n1 = size (phi%dmn, 1)


  !! send the requested elmement values
  if (ubnd /= 0) then
    n = bnd%locl_blk(n_blk + 1) - 1
    do i = 1, n
      j = bnd%locl_idx(i)
      k = j - bnd%elm_blk(blk+1) + 1
      phi%lbnd(:,i) = phi%bnd(:,k)
    end do


    do i = 0, n_blk - 1
      jj1 = bnd%ghst_blk(i + 1)
      jjn = bnd%ghst_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Irecv (phi%gbnd(1,jj1), n1*n, mpi_rp, i, &
        mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
    end do

    do i = 0, n_blk - 1
      jj1 = bnd%locl_blk(i + 1)
      jjn = bnd%locl_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Send (phi%lbnd(1,jj1), n1*n, mpi_rp, i, &
        100, mpi_comm_world, err)
    end do

    call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)
  end if


  if (udmn /= 0) then
    n = dmn%locl_blk(n_blk + 1) - 1
    do i = 1, n
      j = dmn%locl_idx(i)
      k = j - dmn%elm_blk(blk+1) + 1
      phi%ldmn(:,i) = phi%dmn(:,k)
    end do

    do i = 0, n_blk - 1
      jj1 = dmn%ghst_blk(i + 1)
      jjn = dmn%ghst_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Irecv (phi%gdmn(1,jj1), n1*n, mpi_rp, i, &
        mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
    end do

    do i = 0, n_blk - 1
      jj1 = dmn%locl_blk(i + 1)
      jjn = dmn%locl_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Send (phi%ldmn(1,jj1), n1*n, mpi_rp, i, &
        100, mpi_comm_world, err)
    end do

    call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)
  end if


!   do i = 1, size(phi%gbnd)
!     j = bnd%ghst_idx(i)
!     call hash_get_value (bnd%gi2b, j, k)
!     print *, blk+1, "phi_ghst",i,j,k,phi%gbnd(i)
!   end do
end subroutine




subroutine ghost_update_t (bnd, dmn, phi, ubnd, udmn)
  use mpi
  use hash_m

  type(mesh_t),intent(in)::bnd, dmn
  type(tfld_t),intent(inout)::phi
  integer,intent(in)::ubnd, udmn

  integer :: err = 0
  integer :: i, ii, j, jj1, jjn, k, m, n
  integer :: n1, n2
  integer :: blk, n_blk, oblk
  integer, dimension(:), allocatable :: mpi_rqst
  integer, dimension(:,:), allocatable :: mpi_stat


  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)

  allocate (mpi_rqst(n_blk)); mpi_rqst = 0
  allocate (mpi_stat(mpi_status_size, n_blk)); mpi_stat = 0

  n1 = size (phi%dmn, 1)
  n2 = size (phi%dmn, 2)


  !! send the requested elmement values
  if (ubnd /= 0) then
    n = bnd%locl_blk(n_blk + 1) - 1
    do i = 1, n
      j = bnd%locl_idx(i)
      k = j - bnd%elm_blk(blk+1) + 1
      phi%lbnd(:,:,i) = phi%bnd(:,:,k)
    end do


    do i = 0, n_blk - 1
      jj1 = bnd%ghst_blk(i + 1)
      jjn = bnd%ghst_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Irecv (phi%gbnd(1,1,jj1), n1*n2*n, mpi_rp, i, &
        mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
    end do

    do i = 0, n_blk - 1
      jj1 = bnd%locl_blk(i + 1)
      jjn = bnd%locl_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Send (phi%lbnd(1,1,jj1), n1*n2*n, mpi_rp, i, &
        100, mpi_comm_world, err)
    end do

    call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)
  end if


  if (udmn /= 0) then
    n = dmn%locl_blk(n_blk + 1) - 1
    do i = 1, n
      j = dmn%locl_idx(i)
      k = j - dmn%elm_blk(blk+1) + 1
      phi%ldmn(:,:,i) = phi%dmn(:,:,k)
    end do

    do i = 0, n_blk - 1
      jj1 = dmn%ghst_blk(i + 1)
      jjn = dmn%ghst_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Irecv (phi%gdmn(1,1,jj1), n1*n2*n, mpi_rp, i, &
        mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
    end do

    do i = 0, n_blk - 1
      jj1 = dmn%locl_blk(i + 1)
      jjn = dmn%locl_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Send (phi%ldmn(1,1,jj1), n1*n2*n, mpi_rp, i, &
        100, mpi_comm_world, err)
    end do

    call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)
  end if


!   do i = 1, size(phi%gbnd)
!     j = bnd%ghst_idx(i)
!     call hash_get_value (bnd%gi2b, j, k)
!     print *, blk+1, "phi_ghst",i,j,k,phi%gbnd(i)
!   end do
end subroutine
end module
