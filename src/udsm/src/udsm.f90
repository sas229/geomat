include "plaxis_fortran_cpp_interface.f90"

subroutine User_Mod(IDTask, iMod, IsUndr, iStep, iTer, iEl, Int, &
    X, Y, Z, Time0, dTime, Props, Sig0, Swp0, StVar0, dEps, D, &
    Bulk_W, Sig, Swp, StVar, ipl, nStat, NonSym, iStrsDep, &
    iTimeDep, iTang, iPrjDir, iPrjLen, iAbort)

  use plaxis_fortran_cpp_interface
  use, intrinsic :: ISO_C_BINDING
  implicit none

  integer(c_int), intent(inout) :: IDTask, iMod, IsUndr, iStep, iTer
  integer(c_int), intent(in) :: iEl, Int
  real(c_double), intent(in) :: X, Y, Z, Time0, dTime
  real(c_double), intent(in) :: Props(50), Sig0(20), Swp0
  real(c_double), intent(inout) :: StVar0(*)
  real(c_double), intent(in) :: dEps(12)
  real(c_double), intent(inout) :: D(6,6), Bulk_W
  real(c_double), intent(inout) :: Sig(6), Swp, StVar(*)
  integer(c_int), intent(inout) :: ipl, nStat, NonSym
  integer(c_int), intent(inout) :: iStrsDep, iTimeDep, iTang
  character(c_char), intent(in) :: iPrjDir(*)
  integer(c_int), intent(in) :: iPrjLen
  integer(c_int), intent(inout) :: iAbort

  call udsm_cpp_interface(IDTask, iMod, IsUndr, iStep, iTer, &
    iEl, Int, X, Y, Z, Time0, dTime, Props, Sig0, Swp0, StVar0, &
    dEps, D, Bulk_W, Sig, Swp, StVar, ipl, nStat, NonSym, &
    iStrsDep, iTimeDep, iTang, iPrjDir, iPrjLen, iAbort)

end subroutine User_Mod
