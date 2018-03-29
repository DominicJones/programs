module solv_direct_m
  use solv_data_m
  implicit none
  contains
  subroutine solv_direct(imat,c,c_r,x,b)
    integer,intent(in)::imat
    real(wp),dimension(:),intent(in)::c,c_r,b
    real(wp),dimension(:),intent(inout)::x
    integer::ip,in,ivf,ifce

    if(imat>0)then
      do ip=1,size(pl)
        x(ip)=b(ip)

        do ivf=1,size(v_lst(ip)%fce_lst)
          in=v_lst(ip)%ng_lst(ivf)
          if(in>0.and.in<ip)then
            ifce=v_lst(ip)%fce_lst(ivf)
            x(ip)=x(ip)-ng(ifce)%n*x(in)
          end if
        end do

        x(ip)=x(ip)*c(ip)
      end do

      x=x*c_r

      do ip=size(pl),1,-1
        do ivf=1,size(v_lst(ip)%fce_lst)
          in=v_lst(ip)%ng_lst(ivf)
          if(in>0.and.in>ip)then
            ifce=v_lst(ip)%fce_lst(ivf)
            x(ip)=x(ip)-ng(ifce)%p*x(in)
          end if
        end do

        x(ip)=x(ip)*c(ip)
      end do

    else
      x=b*c
    end if
  end subroutine
end module
