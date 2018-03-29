module solv_bicg_m
  use solv_res_m
  use solv_precond_m
  use solv_direct_m
  implicit none
  contains
  subroutine solv_bicg(imat,reduc,n_iter,phi,res0,res0_n,not_r0)
    integer,intent(in)::imat,n_iter
    real(wp),intent(in)::reduc
    logical,intent(in)::not_r0
    type(phi_t),dimension(:),intent(inout)::phi
    real(wp),dimension(:,:),intent(inout)::res0
    real(wp),dimension(:),intent(inout)::res0_n
    integer::ivc,iter,ip,in,ifce,num
    real(wp),dimension(6)::w_sc
    real(wp),dimension(size(pl),5)::w_vec
    real(wp),dimension(size(pl))::c,c_r

    ! calculate residual vector and 1-norm
    call solv_res(phi,res0,res0_n,not_r0)
    if(maxval(res0_n)<small)return

    ! calculate preconditioning vector
    call solv_precond(imat,c,c_r)

    ! loop over each component
    do ivc=1,size(res0_n)

      ! initialise workspace scalars:
      ! res_n=1 alpha=2 beta0=3 beta=4 gam=5 omeg=6
      w_sc=zero
      w_sc(2:3)=one
      w_sc(5)=one

      ! initialise workspace vectors:
      ! res=1 p=2 u=3 v=4 z=5
      w_vec(:,1)=res0(:,ivc)
      w_vec(:,2:5)=zero

      ! iterate:
      num=0
      do iter=1,n_iter
        num=num+1

        ! beta = res.dot.res0
        ! omega = beta.gamma / alpha.beta0
        ! beta0 = beta
        w_sc(4)=dot_product(w_vec(:,1),res0(:,ivc))
        w_sc(6)=w_sc(4)*w_sc(5)/(w_sc(2)*w_sc(3)+small)
        w_sc(3)=w_sc(4)
!         print*,w_sc(4),w_sc(6)

        ! p = res + omega.(p - alpha.u)
        w_vec(:,2)=w_vec(:,1)+w_sc(6) &
          *(w_vec(:,2)-w_sc(2)*w_vec(:,3))

        ! solve: M z = p
        call solv_direct(imat,c,c_r,w_vec(:,5),w_vec(:,2))!; exit

        ! calculate: u = A.z
        w_vec(:,3)=pl*w_vec(:,5)
        do ifce=1,size(ng)
          if(ng(ifce)%in>0)then
            ip=ng(ifce)%ip
            in=ng(ifce)%in
            w_vec(ip,3)=w_vec(ip,3)+ng(ifce)%p*w_vec(in,5)
            w_vec(in,3)=w_vec(in,3)+ng(ifce)%n*w_vec(ip,5)
          end if
        end do

        ! gamma = beta / u.dot.res0
        w_sc(5)=dot_product(w_vec(:,3),res0(:,ivc))
        w_sc(5)=w_sc(4)/(w_sc(5)+small)

        ! phi = phi + gamma.z
        ! res = res - gamma.u
        do ip=1,size(pl)
          phi(ip)%vc(ivc)=phi(ip)%vc(ivc) &
            +w_sc(5)*w_vec(ip,5)
        end do
        w_vec(:,1)=w_vec(:,1)-w_sc(5)*w_vec(:,3)

        ! solve: M z = res
        call solv_direct(imat,c,c_r,w_vec(:,5),w_vec(:,1))

        ! calculate:  v = A.z
        w_vec(:,4)=pl*w_vec(:,5)
        do ifce=1,size(ng)
          if(ng(ifce)%in>0)then
            ip=ng(ifce)%ip
            in=ng(ifce)%in
            w_vec(ip,4)=w_vec(ip,4)+ng(ifce)%p*w_vec(in,5)
            w_vec(in,4)=w_vec(in,4)+ng(ifce)%n*w_vec(ip,5)
          end if
        end do

        ! alpha = v.dot.res / v.dot.v
        w_sc(1)=dot_product(w_vec(:,4),w_vec(:,1))
        w_sc(2)=dot_product(w_vec(:,4),w_vec(:,4))
        w_sc(2)=w_sc(1)/(w_sc(2)+small)

        ! phi = phi + alpha.z
        ! res = res - alpha.v
        ! res_n = sum |res|
        do ip=1,size(pl)
          phi(ip)%vc(ivc)=phi(ip)%vc(ivc) &
            +w_sc(2)*w_vec(ip,5)
        end do
        w_vec(:,1)=w_vec(:,1)-w_sc(2)*w_vec(:,4)
        w_sc(1)=sum(abs(w_vec(:,1)))

!         print*,w_sc(5),w_sc(2),w_sc(3),w_sc(1)

        ! res_n = res_n / res0_n
        w_sc(1)=w_sc(1)/(res0_n(ivc)+small)
        if(w_sc(1)<reduc)exit
        w_sc(1)=zero
      end do
!       print*,'solv_bicg:',ivc,num
    end do
  end subroutine
end module
