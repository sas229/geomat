include "umat.f90"

program main
    use, intrinsic :: ISO_C_BINDING
    
    implicit none
    integer(c_int) :: ndi, nshr, ntens, nstatv, nprops, &
        noel, npt, layer, kspt, jstep(4), kinc, i
    real(c_double) :: stress(6), statev(2), ddsdde(6,6), &
        sse, spd, scd, rpl, ddsddt(6), drplde(6), drpldt, &
        stran(6), dstran(6), time(2), dtime, temp, dtemp, & 
        predef(1), dpred(1), props_array(7), coords(3), drot(3, 3), &
        pnewdt, celent, dfgrd0(3,3), dfgrd1(3,3)
    character(kind=c_char, len=80) :: cmname
    
    ! Set the parameters, initial stress and strain increment.
    ntens = 6
    nstatv = 0
    nprops = 7
    props_array = (/500.0, 0.2, 0.0, 33.0, 5.0, 29.0, 1.1/)
    dstran = (/-0.005, 0.0025, 0.0025, 0.0, 0.0, 0.0/)
    stress = (/-50, -50, -50, 0, 0, 0/)
    cmname = C_CHAR_"C2MC"//C_NULL_CHAR

    ! Run some increments.
    print *, "Model name from Fortran: ", cmname
    write(*,*) stress(1), stress(2), stress(3), stress(4), stress(5), stress(6)
    do while (i < 100)
        call umat(stress, statev, ddsdde, sse, spd, scd, &
        rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, & 
        dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props_array, &
        nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, &
        layer, kspt, jstep, kinc)
        write(*,*) stress(1), stress(2), stress(3), stress(4), stress(5), stress(6)
        i = i + 1
    end do

end program main