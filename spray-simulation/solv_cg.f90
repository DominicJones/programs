module solv_cg_m
  use solv_res_m
  use solv_precond_m
  use solv_direct_m
  implicit none
  contains
  subroutine solv_cg(imat,reduc,n_iter,phi,res0,res0_n,not_r0)
    integer,intent(in)::imat,n_iter
    real(wp),intent(in)::reduc
    logical,intent(in)::not_r0
    type(phi_t),dimension(:),intent(inout)::phi
    real(wp),dimension(:,:),intent(inout)::res0
    real(wp),dimension(:),intent(inout)::res0_n
    integer::ivc,iter,ip,in,ifce,num
    real(kind=wp),dimension(4)::w_sc
    real(kind=wp),dimension(size(pl))::c,c_r
    real(kind=wp),dimension(size(pl),3)::w_vec

    ! calculate residual vector and 1-norm
    call solv_res(phi,res0,res0_n,not_r0)
    if(maxval(res0_n)<small)return

    ! calculate preconditioning vector
    call solv_precond(imat,c,c_r)

    do ivc=1,size(res0_n)

      ! w_sc: res_n=1 alpha=2 beta=3 s=4
      w_sc=zero; w_sc(4)=large

      ! w_vec: res=1 p=2 z=3
      w_vec(:,1)=res0(:,ivc)
      w_vec(:,2:3)=zero

      ! iterate:
      num=0
      do iter=1,n_iter
        num=num+1

        ! solve: M z = res
        call solv_direct(imat,c,c_r,w_vec(:,3),w_vec(:,1))

        ! beta = 1 / s
        w_sc(3)=1./w_sc(4)

        ! s = p.dot.z
        w_sc(4)=dot_product(w_vec(:,1),w_vec(:,3))

        ! beta = beta * s
        w_sc(3)=w_sc(3)*w_sc(4)

        ! p = z + beta.p
        w_vec(:,2)=w_vec(:,3)+w_sc(3)*w_vec(:,2)

        ! calculate: z = A.p
        w_vec(:,3)=pl*w_vec(:,2)
        do ifce=1,size(ng)
          if(ng(ifce)%in>0)then
            ip=ng(ifce)%ip
            in=ng(ifce)%in
            w_vec(ip,3)=w_vec(ip,3)+ng(ifce)%p*w_vec(in,2)
            w_vec(in,3)=w_vec(in,3)+ng(ifce)%n*w_vec(ip,2)
          end if
        end do

        ! alpha = s / p.dot.z
        w_sc(2)=dot_product(w_vec(:,2),w_vec(:,3))
        w_sc(2)=w_sc(4)/(w_sc(2)+small)

        ! phi = phi + alpha.p
        ! res = res - alpha.z
        ! res_n = sum |res|
        do ip=1,size(pl)
          phi(ip)%vc(ivc)=phi(ip)%vc(ivc) &
            +w_sc(2)*w_vec(ip,2)
        end do
        w_vec(:,1)=w_vec(:,1)-w_sc(2)*w_vec(:,3)
        w_sc(1)=sum(abs(w_vec(:,1)))

        ! res_n = res_n / res0_n
        w_sc(1)=w_sc(1)/(res0_n(ivc)+small)
        if(w_sc(1)<reduc)exit
        w_sc(1)=zero
      end do
!       print*,'solv_cg:',ivc,num
    end do
  end subroutine
end module
