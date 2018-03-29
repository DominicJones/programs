!$AD USE stack_m
module prec_m
use mpi

implicit none

integer, parameter :: rsp = kind (1.e0)
integer, parameter :: rdp = kind (1.d0)

integer, parameter :: mpi_rsp = mpi_real
integer, parameter :: mpi_rdp = mpi_double_precision
end module
