include "abaqus_fortran_cpp_interface.f90"

subroutine umat(stress, statev, ddsdde, sse, spd, scd, &
    rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, & 
    dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props_array, &
    nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, &
    layer, kspt, jstep, kinc)

  use abaqus_fortran_cpp_interface
  use , intrinsic :: ISO_C_BINDING
  implicit none

  real(c_double), intent(inout) :: stress(ntens), statev(nstatv), ddsdde(ntens,ntens)
  real(c_double), intent(inout) :: sse, spd, scd, rpl, ddsddt(ntens), drplde(ntens)
  real(c_double), intent(inout) :: drpldt
  real(c_double), intent(in) :: stran(ntens), dstran(ntens), time(2), dtime, temp
  real(c_double), intent(in) :: dtemp, predef(1), dpred(1)
  character(c_char), intent(in) :: cmname(80)
  integer(c_int), intent(in) :: ndi, nshr, ntens, nstatv, nprops
  real(c_double), intent(in) :: props_array(nprops), coords(3), drot(3, 3), pnewdt
  real(c_double), intent(in) :: celent, dfgrd0(3,3), dfgrd1(3,3)
  integer(c_int), intent(in) :: noel, npt, layer, kspt, jstep(4), kinc

  call umat_cpp_interface(stress, statev, ddsdde, sse, spd, scd, &
    rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, & 
    dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props_array, &
    nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, &
    layer, kspt, jstep, kinc)

end subroutine umat