module grd_geom_m
  use grd_conn_m
  implicit none
  contains



  ! calculate grid geometric quantities
  subroutine grd_geom(igrd)
    integer,intent(in)::igrd
    real(kind=wp)::r_p,v_p,t_vol=0.,t_s_f=0.
    real(kind=wp),dimension(:,:),allocatable::vec
    real(kind=wp),dimension(:),allocatable::sc
    real(wp),dimension(n_dx)::n_sum
    integer::i,ivol,ifce,ifce1,ivrt,ivrt1,ivrt2,ibnd,ip,in
    integer::ivf,ifv,ivv,n_fv,ivf2
    integer,dimension(:),allocatable::vrt_lst
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p,vol_n
    type(fce_t),pointer::fce_p
    type(vrt_t),pointer::vrt_p
    real(wp),dimension(n_dx)::x_p_aux,x_n_aux


    print*,'calculating geometry...'

    grd_p=>grd(igrd)
    allocate(vec(1,n_dx))

    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)


      ! volume centre
      vec(1,:)=0.
      do ivrt=1,vol_p%n_vrt
        ivrt2=vol_p%vrt_lst(ivrt)
        vec(1,:)=vec(1,:)+grd_p%vrt(ivrt2)%x
      end do

      allocate(vol_p%x(n_dx))
      vol_p%x=vec(1,:)/vol_p%n_vrt


      ! set boundary types
      if(vol_p%bnd_v)then
      do ivf=1,vol_p%n_fce
        if(vol_p%bnd_lst(ivf)>0)then
          n_fv=vol_ty(vol_p%v_ty)%fce_lst(ivf)%n_vrt
          allocate(vrt_lst(n_fv))

          do ifv=1,n_fv
            ivv=vol_ty(vol_p%v_ty)%fce_lst(ivf)%vrt_lst(ifv)
            vrt_lst(ifv)=vol_p%vrt_lst(ivv)
          end do

          do ivf2=1,vol_p%n_fce
            ifce=vol_p%fce_lst(ivf2)
            fce_p=>grd_p%fce(ifce)
            if(match_vrt(fce_p%vrt_lst,vrt_lst))then
              fce_p%bnd_ty=vol_p%bnd_lst(ivf)
              grd_p%n_bf=grd_p%n_bf+1
              fce_p%ib=grd_p%n_vol+grd_p%n_bf
            end if
          end do
          deallocate(vrt_lst)
        end if
      end do
      end if
    end do

    grd_p%n_vbf=grd_p%n_vol+grd_p%n_bf
    print*,'no. b. faces:',grd_p%n_bf
    print*,'no. i. faces:',grd_p%n_fce-grd_p%n_bf
    print*,'no. vol+b.f.:',grd_p%n_vbf


    do ifce=1,grd_p%n_fce
      fce_p=>grd_p%fce(ifce)
      ip=fce_p%ip
      vol_p=>grd_p%vol(ip)

      allocate(fce_p%x(n_dx)); fce_p%x=0.
      allocate(fce_p%n(n_dx)); fce_p%n=0.
      allocate(fce_p%dx_p(n_dx)); fce_p%dx_p=0.


      ! 2d: normal, face centre & area
      r_p=0.
      if(n_dx==2)then
        if(allocated(vec))deallocate(vec)
        allocate(vec(2,3))
        if(allocated(sc))deallocate(sc)
        allocate(sc(0:1))

        ivrt1=fce_p%vrt_lst(1)
        ivrt2=fce_p%vrt_lst(2)
        vec(1,1:2)=grd_p%vrt(ivrt2)%x-grd_p%vrt(ivrt1)%x
        vec(1,3)=0.
        sc(0)=vec_mag(vec(1,:))

        vec(2,:)=(/0.,0.,1./)
        vec(1,:)=cross_prod(vec(1,:),vec(2,:))
        fce_p%n=vec(1,1:2)/vec_mag(vec(1,1:2))
        fce_p%s=sc(0)

        if(cyl)then
          r_p=0.5*(grd_p%vrt(ivrt1)%x(2) &
            +grd_p%vrt(ivrt2)%x(2))
          fce_p%s=fce_p%s*r_p
        end if

        vec(1,1:2)=grd_p%vrt(ivrt1)%x+grd_p%vrt(ivrt2)%x
        fce_p%x=0.5*vec(1,1:2)
      end if


      ! 3d: normal, face centre & area
      if(n_dx==3)then
        if(allocated(vec))deallocate(vec)
        allocate(vec(0:fce_p%n_vrt,3)); vec=0.
        if(allocated(sc))deallocate(sc)
        allocate(sc(0:1)); sc=0.

        ivrt1=fce_p%vrt_lst(1)
        do ivrt=2,fce_p%n_vrt
          ivrt2=fce_p%vrt_lst(ivrt)
          vec(ivrt,:)=grd_p%vrt(ivrt2)%x-grd_p%vrt(ivrt1)%x

          if(ivrt>2)then
            vec(1,:)=0.5*cross_prod(vec(ivrt-1,:),vec(ivrt,:))
            sc(1)=vec_mag(vec(1,:))
            vec(0,:)=vec(0,:)+vec(1,:)/sc(1)
            sc(0)=sc(0)+sc(1)
          end if
        end do

        vec(0,:)=vec(0,:)/(fce_p%n_vrt-2)
        fce_p%n=vec(0,:)
        fce_p%s=sc(0)

        vec(0,:)=0.
        do ivrt=1,fce_p%n_vrt
          ivrt1=fce_p%vrt_lst(ivrt)
          vec(0,:)=vec(0,:)+grd_p%vrt(ivrt1)%x
        end do
        fce_p%x=vec(0,:)/fce_p%n_vrt
      end if


      ! volume
      v_p=fce_p%x(1)*fce_p%s*fce_p%n(1)
      vol_p%v=vol_p%v+v_p

      ! front face area
      if(n_dx==2.and.cyl)then
        vol_p%s_f=vol_p%s_f+v_p/(r_p+small)
      end if

      ! aux. pole location
      x_p_aux=fce_p%x-fce_p%n &
        *dot_product(fce_p%x-vol_p%x,fce_p%n)

      ! pole to aux. pole displacement
      fce_p%dx_p=x_p_aux-vol_p%x


      ! internal faces
      if(fce_p%in>0)then
        in=fce_p%in
        vol_n=>grd_p%vol(in)

        allocate(fce_p%dx_n(n_dx)); fce_p%dx_n=0.


        ! volume
        v_p=fce_p%x(1)*fce_p%s*fce_p%n(1)
        vol_n%v=vol_n%v-v_p

        ! front face area
        if(n_dx==2.and.cyl)then
          vol_n%s_f=vol_n%s_f-v_p/(r_p+small)
        end if

        ! aux. pole location
        x_n_aux=fce_p%x-fce_p%n &
          *dot_product(fce_p%x-vol_n%x,fce_p%n)

        ! neig. to aux. neig. displacement
        fce_p%dx_n=x_n_aux-vol_n%x


        ! pole to neig. length
        fce_p%l_pn=vec_mag(x_n_aux-x_p_aux)

        ! interpolation factor
        fce_p%w_fp=vec_mag(fce_p%x-x_p_aux) &
          /fce_p%l_pn


      ! boundary faces
      else

        ! pole to face length
        fce_p%l_pn=vec_mag(fce_p%x-x_p_aux)

        ! link face to volume
        fce_p%bnd_v=ip
      end if
    end do


    ! volume reciprocal
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      vol_p%v_r=1./vol_p%v
      t_vol=t_vol+vol_p%v
      t_s_f=t_s_f+vol_p%s_f
    end do

    print*,'volume:',t_vol
    if(n_dx==2.and.cyl)then
      print*,'frontal area:',t_s_f
    end if
  end subroutine
end module
