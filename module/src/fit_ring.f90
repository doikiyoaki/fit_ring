module fit_ring ! only small char, with underscore
    implicit none

    real(8), parameter :: PI = atan(1.0_8)*4.0_8
    real(8), parameter :: G_grav = 6.67259e-8    ! gravitational constant
    real(8), parameter :: c_light = 2.99792458e10 ! light speed
    real(8), parameter :: h_planck = 6.6260755e-27 ! planck constant
    real(8), parameter :: k_boltz = 1.380658e-16 ! boltzman constant

contains

subroutine raytrace_ring( &
    & npixx, npixy, & ! optional
    & ring_r0, ring_tau0, ring_w, ring_h, & 
    & incl, posang, &
    & Temperature_T0, Temperature_r0, Temperature_pow, &
    & wavelength, x_grid, y_grid, Intensity)
    ! CGS unit except for the grid (any unit is okay)
    ! incl, posang in radians

    ! coord type contents
    integer, intent(in) :: npixx
    integer, intent(in) :: npixy

    ! ring type contents
    real(8), intent(in) :: ring_r0
    real(8), intent(in) :: ring_tau0
    real(8), intent(in) :: ring_w
    real(8), intent(in) :: ring_h

    ! angle type contents
    real(8), intent(in) :: incl
    real(8), intent(in) :: posang

    !Temperature type contents
    real(8), intent(in) :: Temperature_T0
    real(8), intent(in) :: Temperature_r0
    real(8), intent(in) :: Temperature_pow

    ! other arguments
    real(8), intent(in) :: x_grid(npixx, npixy), y_grid(npixx, npixy)
    real(8), intent(out) :: Intensity(npixx, npixy)

    real(8) :: truncation, zmax, dz
    integer :: n_integral
    real(8) :: dtau
    real(8) :: y_fixed, x_integ, y_integ_min, y_integ_mid, y_integ_max, z_integ_min, z_integ_mid, z_integ_max!x_fixed
    real(8) :: r_integ_min, r_integ_mid, r_integ_max
    integer :: i, j, k
    real(8) :: wavelength
    real(8) :: HC_p_KbLambda
    real(8) :: posang90
    
    posang90 = posang + PI/2.0_8

    n_integral = 10 ! 5 (~ n_truncation) is enough
    truncation = 4.0

    HC_p_KbLambda = h_planck * c_light / k_boltz / wavelength ! calculate often used constant
    zmax = ring_h * truncation ! upper limit for integration
    dz = zmax / n_integral

    do i = 1, npixx
        do j = 1, npixy
            x_integ =  x_grid(i,j)*cos(posang90) + y_grid(i,j)*sin(posang90)
            y_fixed = -x_grid(i,j)*sin(posang90) + y_grid(i,j)*cos(posang90)

            do k = -n_integral, n_integral-1
                z_integ_min =  k     * dz
                z_integ_mid = (k+0.5)* dz
                z_integ_max = (k+1.0)* dz
                !x_integ = x_fixed
                y_integ_min = y_fixed + z_integ_min*tan(incl)
                y_integ_mid = y_fixed + z_integ_mid*tan(incl)
                y_integ_max = y_fixed + z_integ_max*tan(incl)
                r_integ_min = sqrt(x_integ**2 + y_integ_min**2)
                r_integ_mid = sqrt(x_integ**2 + y_integ_mid**2)
                r_integ_max = sqrt(x_integ**2 + y_integ_max**2)
                dtau = ring_tau0 /(sqrt(2.0 * PI)* ring_h) &                                          ! mutual factor for integration
                *(   exp( -(r_integ_min -ring_r0)**2/(2*ring_w**2) -z_integ_min**2/(2*ring_h**2)) &   ! simpson fitst term
                +4.0*exp( -(r_integ_mid -ring_r0)**2/(2*ring_w**2) -z_integ_mid**2/(2*ring_h**2)) &   ! simpson midterm
                +    exp( -(r_integ_max -ring_r0)**2/(2*ring_w**2) -z_integ_max**2/(2*ring_h**2))) *dz/6.0 ! simpson integration for tau
                Intensity(i,j) = Intensity(i,j) + (1.0-exp(-dtau)) * &
                (Planck_K_short(T_r(Temperature_T0, Temperature_r0, Temperature_pow, r_integ_mid), HC_p_KbLambda) - Intensity(i,j))!& !(1.0_8- exp(-dtau)) * &
                !((Planck_K_short(T_model_T0*((r_integ_mid/T_model_r0)**T_model_pow), HC_p_KbLambda))  - Intensity(i,j))     ! simple integration
                                                                                                        ! it is Okay because T dependence on r is very very very shallow
            end do
        end do
    end do
end subroutine

real(8) function Planck_K_short(Temperature, HC_p_KbLambda)
    real(8) :: Temperature ! [cm]
    real(8) :: HC_p_KbLambda
    Planck_K_short =  HC_p_KbLambda/(exp(HC_p_KbLambda/Temperature)-1.0_8)
end function

real(8) function Planck_K(Temperature, lambda)
    real(8) :: Temperature ! [cm]
    real(8) :: lambda
    real(8) :: HC_p_KbLambda
    HC_p_KbLambda = h_planck*c_light/k_boltz/lambda
    Planck_K =  HC_p_KbLambda/(exp(HC_p_KbLambda/Temperature)-1.0_8)
end function

real(8) function T_r(Temperature_T0, Temperature_r0, Temperature_pow, r)
    real(8), intent(in) :: Temperature_T0
    real(8), intent(in) :: Temperature_r0
    real(8), intent(in) :: Temperature_pow
    real(8), intent(in) :: r
    T_r = Temperature_T0 * ((r/Temperature_r0) ** Temperature_pow)
    if(T_r > 1000.0)then
        T_r = 1000.0
    end if
end function


subroutine convolution(image_in, kernel, image_out, nx, ny, nkx, nky)
    real(8), intent(in) :: image_in(nx, ny)
    real(8), intent(in) :: kernel(nkx, nky)
    real(8), intent(out) :: image_out(nx, ny)
    integer, intent(in) :: nx, ny, nkx, nky
    real(8) :: image_in_extended(nx+nkx-1, ny+nky-1)
    integer :: i, j, k, l
    image_in_extended = 0.0_8
    image_out = 0.0_8
    image_in_extended((nkx+1)/2:nx+(nkx-1)/2,(nky+1)/2:ny+(nky-1)/2) = image_in

    do i = 1, nx
        do j = 1, ny
            do k = 1, nkx
                do l = 1, nky
                    image_out(i,j) = image_out(i,j) + image_in_extended(i+k-1,j+l-1)*kernel(k,l)
                end do
            end do
        end do
    end do
end subroutine


subroutine interpolate_data(x1d, y1d, image, x_interp, y_interp, out, nx, ny, n_interp)
    ! grid: image[x,y]
    ! use image.T from python to pass image
    integer, intent(in) :: nx, ny, n_interp
    real(8), intent(in) :: x1d(nx), y1d(ny)
    real(8), intent(in) :: image(nx, ny)
    real(8), intent(in) :: x_interp(n_interp), y_interp(n_interp)
    real(8), intent(out) :: out(n_interp)

    integer :: i, j, k
    real(8) :: dx, dy
    real(8) :: r_x, r_y

    dx = x1d(2) - x1d(1)
    dy = y1d(2) - y1d(1)

    do k = 1, n_interp
        i = int((x_interp(k) - x1d(1))/dx) + 1
        j = int((y_interp(k) - y1d(1))/dy) + 1
        r_x = (x_interp(k) - x1d(i))/dx
        r_y = (y_interp(k) - y1d(j))/dy
        out(k) = (1.0_8-r_x)*(1.0_8-r_y)*image(i,j) &
               + (1.0_8-r_x)*r_y        *image(i,j+1) &
               + r_x        *(1.0_8-r_y)*image(i+1,j) &
               + r_x       *r_y         *image(i+1,j+1)
    end do
end subroutine

end module fit_ring