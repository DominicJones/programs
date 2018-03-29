module solv_precond_m
  use solv_data_m
  implicit none
  contains
  subroutine solv_precond(imat,c,c_r)
    integer,intent(in)::imat
    real(wp),dimension(:),intent(inout)::c,c_r
    integer::ip,in,ivf,ifce

    do ip=1,size(pl)
      c(ip)=pl(ip)

      if(imat>0)then
        do ivf=1,size(v_lst(ip)%fce_lst)
          in=v_lst(ip)%ng_lst(ivf)
          if(in>0.and.in<ip)then
            ifce=v_lst(ip)%fce_lst(ivf)
            if(imat==1)then
              c(ip)=c(ip)-ng(ifce)%n**2*c(in)
            else if(imat==2)then
              c(ip)=c(ip)-ng(ifce)%p*c(in)*ng(ifce)%n
            end if
          end if
        end do
      end if

      c_r(ip)=c(ip)
      c(ip)=one/(c(ip)+small)
    end do
  end subroutine
end module
