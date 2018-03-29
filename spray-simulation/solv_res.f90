module solv_res_m
  use solv_data_m
  implicit none
  contains
  subroutine solv_res(phi,res0,res0_n,not_r0)
    type(phi_t),dimension(:),intent(inout)::phi
    real(wp),dimension(:,:),intent(inout)::res0
    real(wp),dimension(:),intent(inout)::res0_n
    logical,intent(in)::not_r0
    integer::ip,in,ifce

    if(.not.not_r0)then
    do ip=1,size(pl)
      res0(ip,:)=res0(ip,:)-pl(ip)*phi(ip)%vc
    end do
    do ifce=1,size(ng)
      if(ng(ifce)%in>0)then
        ip=ng(ifce)%ip
        in=ng(ifce)%in
        res0(ip,:)=res0(ip,:)-ng(ifce)%p*phi(in)%vc
        res0(in,:)=res0(in,:)-ng(ifce)%n*phi(ip)%vc
      end if
    end do
    end if

    res0_n=zero
    do ip=1,size(pl)
      res0_n=res0_n+abs(res0(ip,:))
    end do
  end subroutine
end module
