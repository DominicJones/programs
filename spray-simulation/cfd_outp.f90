module cfd_outp_m
  use grd_constr_m
  implicit none
  contains
!
! write equation output file for Tecplot
!
  subroutine eqn_outp(igrd,iph,ieq,ict,dir)
    integer,intent(in)::igrd,iph,ieq,ict
    character(len=*),intent(in)::dir
    integer::iunit=100,l0,l
    character(len=60)::dx_var,et_kind,cv_vars,title
    character(len=4)::ict_ch,ict_ch0='0000'
    character(len=8),dimension(2)::int_ch
    character(len=500)::phi_var
    character(len=1)::igrd_ch,iph_ch,ieq_ch
    integer,dimension(:),allocatable::vrt_lst
    integer::idx,ivrt,ivol,isc,ivc,n_vrt,l1,l2
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(ieq_t),pointer::ieq_p
    type(eqn_t),pointer::eqn_p

    ! title
    write(igrd_ch,'(i1)')igrd
    write(iph_ch,'(i1)')iph
    write(ieq_ch,'(i1)')ieq
    title='grd'//igrd_ch// &
      '%'//'phs'//iph_ch// &
      '%'//'eqn'//ieq_ch
    if(n_tl>1.and.ict>0)then
      write(ict_ch,'(i4)')ict
      ict_ch=adjustl(ict_ch)
      l0=len(trim(ict_ch0))
      l=len(trim(ict_ch))
      ict_ch=ict_ch0(1:l0-l)//ict_ch(1:l)
      title=trim(title)//'@'//trim(ict_ch)
    end if
    open(iunit,file=trim(dir)//trim(title)//'.dat')

    ! header
    if(n_dx==2)then
      dx_var='"x_[1]","x_[2]"'
      et_kind='quadrilateral'; n_vrt=4
    end if
    if(n_dx==3)then
      dx_var='"x_[1]","x_[2]","x_[3]"'
      et_kind='brick'; n_vrt=8
    end if
    allocate(vrt_lst(n_vrt)); vrt_lst=0

    grd_p=>grd(igrd)
    ieq_p=>grd_p%iph(iph)%ieq(ieq)

    phi_var=''
    do isc=1,ieq_p%n_sc
      write(int_ch(1),'(i2)')isc
      do ivc=1,ieq_p%n_vc
        write(int_ch(2),'(i2)')ivc
        phi_var=trim(phi_var)//',"phi_['// &
          trim(adjustl(int_ch(1)))//','//trim(adjustl(int_ch(2)))//']"'
      end do
    end do

    write(iunit,*)'title = "',trim(title),'"'
    write(iunit,*)'variables = '//trim(dx_var)//trim(phi_var)
    write(iunit,*)'zone n = ',grd_p%n_vrt,', e = ',grd_p%n_vol,',&
      & f = feblock, et = '//trim(et_kind)
    write(int_ch(1),'(i2)')n_dx+1
    cv_vars=trim(adjustl(int_ch(1)))//'-'
    write(int_ch(1),'(i2)')99
    cv_vars=trim(cv_vars)//trim(adjustl(int_ch(1)))
    write(iunit,*)'varlocation = (['//trim(cv_vars)//']=cellcentered)'

    ! coordinates of vertices
    do idx=1,n_dx
      do ivrt=1,grd_p%n_vrt
        write(iunit,*)grd(igrd)%vrt(ivrt)%x(idx)
      end do
    end do

    ! equation variable
    do isc=1,ieq_p%n_sc
      do ivc=1,ieq_p%n_vc
        do ivol=1,grd(igrd)%n_vol
          vol_p=>grd_p%vol(ivol)
          eqn_p=>vol_p%phs(iph)%eqn(ieq)
          write(iunit,*)eqn_p%phi(isc,ivc,1)
        end do
      end do
    end do

    ! volume vertex arrangement
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      call grd_postpr(vrt_lst1=vol_p%vrt_lst,vrt_lst2=vrt_lst)
      write(iunit,*)vrt_lst
    end do
    close(iunit)
    print*,'written: '//trim(dir)//trim(title)//'.dat'
  end subroutine
!
! write property output file for Tecplot
!
  subroutine prp_outp(igrd,iph,ipr,ict,dir,p)
    integer,intent(in)::igrd,iph,ipr,ict
    character(len=*),intent(in)::dir
    logical,intent(in),optional::p
    integer::iunit=100,l0,l
    character(len=60)::dx_var,et_kind,cv_vars,title
    character(len=4)::ict_ch,ict_ch0='0000'
    character(len=8),dimension(2)::int_ch
    character(len=200)::phi_var,phi_res
    character(len=1)::igrd_ch,iph_ch,ipr_ch
    integer,dimension(:),allocatable::vrt_lst
    integer::idx,ivrt,ivol,isc,ivc,n_vrt,l1,l2
    type(grd_t),pointer::grd_p
    type(vol_t),pointer::vol_p
    type(ipr_t),pointer::ipr_p
    type(prp_t),pointer::prp_p

    ! title
    write(igrd_ch,'(i1)')igrd
    write(iph_ch,'(i1)')iph
    write(ipr_ch,'(i1)')ipr
    title='grd'//igrd_ch// &
      '%'//'phs'//iph_ch// &
      '%'//'prp'//ipr_ch
    if(n_tl>1.and.ict>0)then
      write(ict_ch,'(i4)')ict
      ict_ch=adjustl(ict_ch)
      l0=len(trim(ict_ch0))
      l=len(trim(ict_ch))
      ict_ch=ict_ch0(1:l0-l)//ict_ch(1:l)
      title=trim(title)//'@'//trim(ict_ch)
    end if
    open(iunit,file=trim(dir)//trim(title)//'.dat')

    ! header
    if(n_dx==2)then
      dx_var='"x_[1]","x_[2]"'
      et_kind='quadrilateral'; n_vrt=4
    end if
    if(n_dx==3)then
      dx_var='"x_[1]","x_[2]","x_[3]"'
      et_kind='brick'; n_vrt=8
    end if
    allocate(vrt_lst(n_vrt)); vrt_lst=0

    grd_p=>grd(igrd)
    phi_var=',"rho_[1]"'

    write(iunit,*)'title = "',trim(title),'"'
    write(iunit,*)'variables = '//trim(dx_var)//trim(phi_var)
    write(iunit,*)'zone n = ',grd_p%n_vrt,', e = ',grd_p%n_vol,',&
      & f = feblock, et = '//trim(et_kind)
    write(int_ch(1),'(i2)')n_dx+1
    cv_vars=trim(adjustl(int_ch(1)))//'-'
    write(int_ch(1),'(i2)')99
    cv_vars=trim(cv_vars)//trim(adjustl(int_ch(1)))
    write(iunit,*)'varlocation = (['//trim(cv_vars)//']=cellcentered)'

    ! coordinates of vertices
    do idx=1,n_dx
      do ivrt=1,grd_p%n_vrt
        write(iunit,*)grd(igrd)%vrt(ivrt)%x(idx)
      end do
    end do

    ! property variable
    do ivol=1,grd(igrd)%n_vol
      vol_p=>grd_p%vol(ivol)
      if(present(p).and.p)then
        write(iunit,*)vol_p%p
      else
        prp_p=>vol_p%phs(iph)%prp(ipr)
        write(iunit,*)prp_p%rho(1)
      end if
    end do

    ! volume vertex arrangement
    do ivol=1,grd_p%n_vol
      vol_p=>grd_p%vol(ivol)
      call grd_postpr(vrt_lst1=vol_p%vrt_lst,vrt_lst2=vrt_lst)
      write(iunit,*)vrt_lst
    end do
    close(iunit)
    print*,'written: '//trim(dir)//trim(title)//'.dat'
  end subroutine
end module
