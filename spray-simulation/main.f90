program main
  use grd_geom_m
  use grd_solv_m
  use cfd_solve_m
  implicit none
  integer::igrd
  real(wp)::t1,t2
  character(len=50)::src_dir,cas_dir,outp_dir

  call cpu_time(t1)


  ! define vertex ordering
  src_dir='../src/'
!   call grd_prepr(trim(src_dir)//'gambit.vrt')
!   call grd_postpr(trim(src_dir)//'tecplot.vrt')
  call grd_prepr('../etc/gambit.vrt')
  call grd_postpr('../etc/tecplot.vrt')


  ! set-up grid
  cas_dir='../cas/2d_spray/'
  outp_dir=trim(cas_dir)//'output/'

  n_tl=3
  n_grd=1
  l_sc=ten**(3)
  cyl=.true.
  allocate(grd(n_grd))

  igrd=1
  call grd_read(igrd,trim(cas_dir)//'closed.neu')
  call grd_conn(igrd)
  call grd_geom(igrd)
  call grd_solv(igrd)


  ! solve flow
  call c_std(igrd)
  call cfd_solve(igrd,trim(outp_dir))


  call cpu_time(t2)
  print*,'cpu time (min):',(t2-t1)/60._wp
end program
