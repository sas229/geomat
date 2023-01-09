module abaqus_fortran_cpp_interface
    interface
      subroutine umat_cpp_interface(stress, statev, ddsdde, sse, spd, scd, &
          rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, & 
          dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props_array, &
          nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, &
          layer, kspt, jstep, kinc) & 
      bind(C, name="umat_cpp_interface")  
      use , intrinsic :: ISO_C_BINDING
      implicit none
      real(c_double) :: stress(ntens), statev(nstatv), ddsdde(ntens,ntens), &
          sse, spd, scd, rpl, ddsddt(ntens), drplde(ntens), drpldt, &
          stran(ntens), dstran(ntens), time(2), dtime, temp, dtemp, & 
          predef(1), dpred(1), props_array(nprops), coords(3), drot(3, 3), &
          pnewdt, celent, dfgrd0(3,3), dfgrd1(3,3)
      integer(c_int) :: ndi, nshr, ntens, nstatv, nprops, &
          noel, npt, layer, kspt, jstep(4), kinc
      character(c_char) :: cmname(*)
      end subroutine umat_cpp_interface
    end interface
  end module abaqus_fortran_cpp_interface