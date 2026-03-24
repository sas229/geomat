module plaxis_fortran_cpp_interface
    interface
      subroutine udsm_cpp_interface(IDTask, iMod, IsUndr, iStep, iTer, &
          iEl, Int, X, Y, Z, Time0, dTime, Props, Sig0, Swp0, StVar0, &
          dEps, D, Bulk_W, Sig, Swp, StVar, ipl, nStat, NonSym, &
          iStrsDep, iTimeDep, iTang, iPrjDir, iPrjLen, iAbort) &
      bind(C, name="udsm_cpp_interface")
      use, intrinsic :: ISO_C_BINDING
      implicit none
      integer(c_int) :: IDTask, iMod, IsUndr, iStep, iTer, iEl, Int
      real(c_double) :: X, Y, Z, Time0, dTime
      real(c_double) :: Props(50), Sig0(20), Swp0, StVar0(*)
      real(c_double) :: dEps(12)
      real(c_double) :: D(6,6), Bulk_W
      real(c_double) :: Sig(6), Swp, StVar(*)
      integer(c_int) :: ipl, nStat, NonSym, iStrsDep, iTimeDep, iTang
      character(c_char) :: iPrjDir(*)
      integer(c_int) :: iPrjLen, iAbort
      end subroutine udsm_cpp_interface
    end interface
  end module plaxis_fortran_cpp_interface
