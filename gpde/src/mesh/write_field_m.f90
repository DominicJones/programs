module write_field_m
implicit none
contains



subroutine write_field(ios,fnm,msh,phi,n_cmp,time,tstep)
  use const_m
  use sort_m
  use mesh_data_m

  character(*),intent(in)::fnm
  integer,intent(in)::ios,n_cmp
  type(msh_t),intent(in)::msh
  real(rp),dimension(:,:),intent(in)::phi
  real(rp),intent(in),optional::time
  integer,intent(in),optional::tstep

  character(50),dimension(1)::str_tag
  real(rp),dimension(1)::real_tag
  integer,dimension(4)::int_tag
  logical::vrtdata
  integer::i,ib,n_vrt,n_elm,n_ebf,key


  ! posible field data sizes
  n_vrt=msh%n_vrt
  n_elm=msh%n_elm
  n_ebf=msh%n_ebf

  ! tags
  str_tag(1)=trim(fnm) ! name of view

  real_tag=0
  if(present(time))real_tag(1)=time ! time

  int_tag=0
  if(present(tstep))int_tag(1)=tstep ! time step index
  int_tag(2)=1 ! number of field components (1|3|^9)
  if(n_cmp>1)int_tag(2)=3
  int_tag(3)=size(phi)/n_cmp ! number of elements
  int_tag(4)=0 ! partition index

  vrtdata=.false.
  if(int_tag(3)==n_vrt)vrtdata=.true.


  ! write file
  print "(1x,a,a)","creating: ",trim(fnm)
  open(ios,file=trim(fnm))


  write(ios,"(a)")"$MeshFormat"
  write(ios,"(3i7)")2,0,8
  write(ios,"(a)")"$EndMeshFormat"


  if(vrtdata)then
    write(ios,"(a)")"$NodeData"
    print "(1x,a)","writing vertex field data"
  else
    write(ios,"(a)")"$ElementData"
    if(size(phi)/n_cmp==n_ebf)then
      print "(1x,a)","writing element+boundary field data"
    else
      print "(1x,a)","writing element field data"
    end if
  end if


  write(ios,"(i7)")size(str_tag)
  do i=1,size(str_tag)
    write(ios,"(a)")str_tag(i)
  end do

  write(ios,"(i7)")size(real_tag)
  do i=1,size(real_tag)
    write(ios,"(es14.4)")real_tag(i)
  end do

  write(ios,"(i7)")size(int_tag)
  do i=1,size(int_tag)
    write(ios,"(i7)")int_tag(i)
  end do


  if(vrtdata)then
    do i=1,msh%n_vrt
      key=msh%vrt_tag(i)%lst(1)
      if(n_cmp==1)then
        write(ios,"(i7,3es14.4)")key,phi(1,i)
      else if(n_cmp==2)then
        write(ios,"(i7,3es14.4)")key,phi(:,i),0.d0
      else if(n_cmp==3)then
        write(ios,"(i7,3es14.4)")key,phi(:,i)
      end if
    end do

  else

    do i=1,msh%n_elm
      key=msh%elm_tag(i)%lst(1)
      if(n_cmp==1)then
        write(ios,"(i7,3es14.4)")key,phi(1,i)
      else if(n_cmp==2)then
        write(ios,"(i7,3es14.4)")key,phi(:,i),0.d0
      else if(n_cmp==3)then
        write(ios,"(i7,3es14.4)")key,phi(:,i)
      end if
    end do

    if(size(phi)/n_cmp==n_ebf)then
      ib=msh%n_elm
      do i=1,msh%n_bfce
        key=msh%fce_tag(i)%lst(1)
        ib=ib+1
        if(n_cmp==1)then
          write(ios,"(i7,3es14.4)")key,phi(1,ib)
        else if(n_cmp==2)then
          write(ios,"(i7,3es14.4)")key,phi(:,ib),0.d0
        else if(n_cmp==3)then
          write(ios,"(i7,3es14.4)")key,phi(:,ib)
        end if
      end do
    end if

  end if


  if(vrtdata)then
    write(ios,"(a)")"$EndNodeData"
  else
    write(ios,"(a)")"$EndElementData"
  end if


  close(ios)
  print "(1x,a,a)","closed:   ",trim(fnm)
end subroutine
end module
