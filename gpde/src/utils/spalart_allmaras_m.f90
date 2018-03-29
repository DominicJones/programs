module spalart_allmaras_m
use const_m
use turbulence_m
implicit none

! model parameters
real(rp),parameter,private::sigma=2*third
real(rp),parameter,private::c_prod=2
real(rp),parameter,private::c_v=7.1

real(rp),parameter,dimension(2),private::&
  c_b=(/0.1355_rp,0.622_rp/)

real(rp),parameter,dimension(3),private::&
  c_w=(/c_b(1)/kappa**2+(one+c_b(2))/sigma,0.3_rp,2._rp/)

real(rp),parameter,dimension(4),private::&
  c_t=(/1._rp,2._rp,1.2_rp,0.5_rp/)

contains



subroutine sa_turb_visc(msh,dens,visc,vis_sa,vis_t)
  use mesh_data_m
  type(msh_t),intent(in)::msh
  real(rp),dimension(:),intent(inout)::dens,visc,vis_sa,vis_t

  real(rp)::cof(2),term(3)
  integer::i,k,ib,ph_ty

  cof(1)=c_v**3

!$AD II-LOOP
  do i=1,msh%n_elm
    cof(2)=dens(i)/visc(i)
    term(1)=vis_sa(i)*cof(2)
!     term(1)=vis_sa(i)/(visc(i)/dens(i))
    term(2)=term(1)**3
    term(3)=term(2)/(term(2)+cof(1))
    vis_t(i)=term(3)*dens(i)*vis_sa(i)
  end do
end subroutine



subroutine sa_turb_diff(msh,dens,visc,vis_sa,cd)
  use mesh_data_m
  type(msh_t),intent(in)::msh
  real(rp),dimension(:),intent(inout)::dens,visc,vis_sa,cd

  real(rp)::cof(2),term(2)
  integer::i,k,ib

  cof(1)=one/sigma
  cof(2)=cof(1)*c_b(2)

!$AD II-LOOP
  do i=1,msh%n_elm
    term(1)=cof(1)*(visc(i)+dens(i)*vis_sa(i))
    term(2)=cof(2)*dens(i)
    cd(i)=term(1)+term(2)
  end do

!$AD II-LOOP
  do k=1,msh%n_bfce
    i=msh%fce_elm(k)%lst(1)
    ib=abs(msh%fce_elm(k)%lst(2))
    term(1)=cof(1)*(visc(ib)+dens(ib)*vis_sa(ib))
    term(2)=cof(2)*dens(ib)
    cd(ib)=term(1)+term(2)
  end do
end subroutine



subroutine sa_turb_source(msh,d_wall,vol,dens,visc,vis_t,vis_sa,gradu,rhs)
  use mesh_data_m
  type(msh_t),intent(in)::msh
  real(rp),dimension(:),intent(inout)::d_wall,vol
  real(rp),dimension(:),intent(inout)::dens,visc,vis_t,vis_sa
  real(rp),dimension(:,:,:),intent(inout)::gradu
  real(rp),dimension(:),intent(inout)::rhs

  real(rp)::dudy,dvdx,stress,rotate,production,destruction
  real(rp)::cof(3),term(13)
  integer::i,i1,i2

  cof(1)=c_v**3
  cof(2)=kappa**2
  cof(3)=c_w(3)**6


!$AD II-LOOP
  do i=1,msh%n_elm

    ! \Chi
    term(1)=vis_sa(i)/(visc(i)/dens(i))

    ! f_{v1}
    term(2)=term(1)**3
    term(3)=term(2)/(term(2)+cof(1))

    ! f_{v2}
    term(4)=term(1)/(1.d0+term(1)*term(3))
    term(5)=one-term(4)

    ! production
    stress=zero
    rotate=zero
!$AD NO-II-LOOP # leave as NO-
    do i1=1,msh%n_dx
!$AD NO-II-LOOP # leave as NO-
      do i2=1,msh%n_dx
        dudy=gradu(i2,i1,i)
        dvdx=gradu(i1,i2,i)
        stress=stress+2*(half*(dudy+dvdx))**2
        rotate=rotate+2*(half*(dudy-dvdx))**2
      end do
    end do
    stress=sqrt(stress)
    rotate=sqrt(rotate)
    term(6)=rotate
    term(6)=term(6)+c_prod*min(zero,stress-rotate)

    term(7)=vis_sa(i)/(cof(2)*d_wall(i)**2)*term(5)
    term(8)=term(6)+term(7)
    term(9)=c_b(1)*dens(i)*term(8)*vis_sa(i)
    production=term(9)

    ! destruction
    term(10)=vis_sa(i)/(term(8)*cof(2)*d_wall(i)**2)
    term(11)=term(10)+c_w(2)*(term(10)**6-term(10))
    term(12)=term(11)*((one+cof(3))/(term(11)**6+cof(3)))**(one/6)
    term(13)=c_w(1)*dens(i)*term(12)*(vis_sa(i)/d_wall(i))**2
    destruction=term(13)

    ! source term
    rhs(i)=rhs(i)+(production-destruction)*vol(i)
  end do
end subroutine
end module
