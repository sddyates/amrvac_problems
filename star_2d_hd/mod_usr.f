module mod_usr
  use mod_hd

  implicit none

  ! Custom variables can be defined here
  ! ...

  integer, parameter :: dp = selected_real_kind(15, 300)

  ! Assigning constants.
  real(kind=dp), parameter :: pi = 3.1415926536_dp
  real(kind=dp), parameter :: kb = 1.38064852e-16_dp ![erg K-1]
  real(kind=dp), parameter :: mp = 1.67262171e-24_dp ![g]
  real(kind=dp), parameter :: c = 3.0e+5_dp          ![km s-1]
  real(kind=dp), parameter :: Big_G = 6.67259e-8_dp  ![cm3 g-1 s-2]
  real(kind=dp), parameter :: Msun = 1.989e+33_dp    ![g]
  real(kind=dp), parameter :: Rsun = 6.955e+10_dp    ![cm]
  real(kind=dp), parameter :: L_sun = 3.846e+33_dp
  ! Stellar parameters
  real(kind=dp), parameter :: Cs_p = 4.0_dp
  real(kind=dp), parameter :: Mratio = 26.6_dp
  real(kind=dp), parameter :: Lratio = 1.15e+5_dp
  real(kind=dp), parameter :: Bcgs = 300.0_dp
  real(kind=dp), parameter :: T = 36.3e+3_dp
  real(kind=dp), parameter :: mu = 1.09_dp
  real(kind=dp), parameter :: a = 0.6_dp
  real(kind=dp), parameter :: b = 0.8_dp
  real(kind=dp), parameter :: Qfac = 700.0_dp
  real(kind=dp), parameter :: a_eff = 0.55_dp
  real(kind=dp), parameter :: beta = 0.0_dp

  real(kind=dp), parameter :: UNIT_DENSITY = 1.0e-12_dp
  real(kind=dp), parameter :: UNIT_LENGTH = 6.955e+10_dp*9.0_dp
  real(kind=dp), parameter :: UNIT_VEL = 1.0e+5_dp
  real(kind=dp), parameter :: UNIT_TIME = UNIT_LENGTH/UNIT_VEL
  real(kind=dp), parameter :: UNIT_MASS = UNIT_DENSITY*UNIT_LENGTH**3
  real(kind=dp), parameter :: UNIT_G = Big_G*UNIT_DENSITY*UNIT_TIME**2
  real(kind=dp), parameter :: UNIT_B = sqrt(4.0_dp*pi*UNIT_DENSITY*UNIT_VEL*UNIT_VEL)
  real(kind=dp), parameter :: UNIT_L = UNIT_TIME**(-3)*UNIT_DENSITY*UNIT_LENGTH**5

  real(kind=dp), parameter :: M_star = (Mratio*Msun/UNIT_MASS)
  real(kind=dp), parameter :: Edd = (2.6e-5_dp*(Lratio)*(1.0_dp/Mratio))
  real(kind=dp), parameter :: L = (Lratio*L_sun/UNIT_L)
  real(kind=dp), parameter :: M_dot = (1.0_dp+a_eff)**(-1.0_dp/a_eff) * a_eff * (1.0-a_eff)**(-1)*&
                              (L/(c*c)) * ( Qfac * Edd * (1.0_dp - Edd)**(-1) )**(a_eff**(-1) - 1.0_dp)
  real(kind=dp), parameter :: cs = sqrt(2.0_dp*kb*T/mp)/UNIT_VEL
  real(kind=dp), parameter :: Bq = Bcgs/UNIT_B
  real(kind=dp), parameter :: v_esc = sqrt(2.0_dp*UNIT_G*M_star*(1.0_dp - Edd))
  real(kind=dp), parameter :: v_inf = v_esc * sqrt((a/(1.0_dp - a)))

  real(kind=dp) :: vv
  real(kind=dp) :: rho
  real(kind=dp) :: vx
  real(kind=dp) :: vy
  real(kind=dp) :: prs
  real(kind=dp) :: eng

contains

  subroutine usr_init()
    use mod_usr_methods

    call set_coordinate_system("Cartesian")

    usr_init_one_grid => initial_conditions
    usr_internal_bc => intern_bc
    usr_gravity => gravity
    usr_source => source_terms
    !usr_get_dt => get_usr_dt

    hd_gamma = 1.05_dp

    call hd_activate()

  end subroutine usr_init

  subroutine initial_conditions(ixGmin1, ixGmin2, ixGmax1, ixGmax2,& 
                                ixmin1, ixmin2,  ixmax1, ixmax2, w, x)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1, ixGmin2, ixGmax1, ixGmax2, &
                           ixmin1, ixmin2, ixmax1, ixmax2
    real(kind=dp), intent(in) :: x(ixGmin1:ixGmax1, ixGmin2:ixGmax2, ndim)
    real(kind=dp), intent(inout) :: w(ixGmin1:ixGmax1, ixGmin2:ixGmax2, nw)
    integer :: i, j
    real(kind=dp) :: r

    do j=ixmin2,ixmax2
      do i=ixmin1,ixmax1

        r = sqrt(x(i,j,1)**2 + x(i,j,2)**2)

        vv = v_inf*(1.0 - (1.0/r))**b

        if (r > 1.0) then
          rho = (M_dot/(4.0*pi*vv*r*r))
          vx = vv*x(i,j,1)/r
          vy = vv*x(i,j,2)/r
          prs = cs*cs*rho/hd_gamma
          eng = prs/(rho*(hd_gamma - 1.0_dp)) + 0.5*rho*(vx**2 + vy**2)
        else if (r > 0.5 .and. r <= 1.0) then
          rho = M_dot/(4.0*pi*(cs/Cs_p))
          vx = 0.0
          vy = 0.0
          prs = cs*cs*rho/hd_gamma
          eng = prs/(rho*(hd_gamma - 1.0_dp)) + 0.5*rho*(vx**2 + vy**2)
        else if (r <= 0.5) then
          rho = M_dot/(4.0*pi*(cs/Cs_p))
          vx = 0.0
          vy = 0.0
          prs = cs*cs*rho/hd_gamma
          eng = prs/(rho*(hd_gamma - 1.0_dp)) + 0.5*rho*(vx**2 + vy**2)
        end if

        w(i,j,rho_) = rho
        w(i,j,mom(1)) = rho*vx
        w(i,j,mom(2)) = rho*vy
        w(i,j,e_) = eng

      end do
    end do

  end subroutine initial_conditions

  subroutine intern_bc(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
        ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,level
    real(kind=dp), intent(in) :: qt
    real(kind=dp), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp), intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    integer :: i, j
    real(kind=dp) :: r

    do j=ixImin2,ixImax2
      do i=ixImin1,ixImax1

        r = sqrt(x(i,j,1)**2 + x(i,j,2)**2)

        if (r <= 0.5_dp) then
          rho = M_dot/(4.0_dp*pi*(cs/Cs_p))
          vx = 0.0_dp
          vy = 0.0_dp
          prs = cs*cs*rho/hd_gamma
          eng = prs/(rho*(hd_gamma - 1.0_dp)) + 0.5_dp*rho*(vx**2 + vy**2)
          w(i,j,rho_) = rho
          w(i,j,mom(1)) = rho*vx
          w(i,j,mom(2)) = rho*vy
          w(i,j,e_) = eng
        else if (r <= 1.0_dp .and. r > 0.5_dp) then
          rho = M_dot/(4.0_dp*pi*(cs/Cs_p))
          vx = (cs/Cs_p)*x(i,j,1)/r
          vy = (cs/Cs_p)*x(i,j,2)/r
          prs = cs*cs*rho/hd_gamma
          eng = prs/(rho*(hd_gamma - 1.0_dp)) + 0.5_dp*rho*(vx**2 + vy**2)
          !print*, sqrt(vx**2 + vy**2), rho, r
          w(i,j,rho_) = rho
          w(i,j,mom(1)) = rho*vx
          w(i,j,mom(2)) = rho*vy
          w(i,j,e_) = eng
          ! This needs to be done for 4 quadrents,
          !if x(i+1,j,1) outside and x(i,j,1) is inside, then
          !  w(i,j,mom(1)) = w(i+1,j,mom(1))
          !  w(i-1,j,mom(1)) = w(i+1,j,mom(1))
          !end if

          !if x(i,j+1,2) outside and x(i,j,2) is inside, then
          !  w(i,j,mom(2)) = w(i,j+1,mom(2))
          !  w(i,j-1,mom(2)) = w(i,j+1,mom(2))
          !end if
        end if

      end do
    end do

  end subroutine intern_bc

  subroutine gravity(ixImin1,ixImin2,ixImax1,ixImax2,&
                     ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmax1,ixOmax2
    real(kind=dp), intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp), intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp), intent(out) :: gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,ndim)
    integer :: i, j
    real(kind=dp) :: r, GM

    do j=ixImin2,ixImax2
      do i=ixImin1,ixImax1

        r = sqrt(x(i,j,1)**2 + x(i,j,2)**2)
        GM = -UNIT_G*(M_star - Edd)/r/r

        if (r > 0.5_dp) then
          gravity_field(i, j, 1) = GM*x(i,j,1)/r
          gravity_field(i, j, 2) = GM*x(i,j,2)/r
        else if (r <= 0.5_dp) then
          gravity_field(i, j, 1) = 0.0!GM*x(i,j,1)/r
          gravity_field(i, j, 2) = 0.0!GM*x(i,j,2)/r
        end if

      end do
    end do

  end subroutine gravity

  subroutine source_terms(qdt, ixImin1, ixImin2, ixImax1, ixImax2, ixOmin1,&
                          ixOmin2, ixOmax1, ixOmax2, iwmin, iwmax, &
                          qtC, wCT, qt, w, x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1, ixImin2, ixImax1, ixImax2
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmax1, ixOmax2
    integer, intent(in) :: iwmin, iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    call source_cak(ixImin1, ixImin2, ixImax1, ixImax2,&
                    ixOmin1, ixOmin2, ixOmax1, ixOmax2,&
                    w, wCT, x, qdt)

  end subroutine source_terms

  subroutine get_new_dt(w, ixGmin1, ixGmin2, ixGmax1, ixGmax2, ixmin1, ixmin2, &
                        ixmax1, ixmax2, dtnew, dx1, dx2, x)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1, ixGmin2, ixGmax1, ixGmax2
    integer, intent(in) :: ixmin1, ixmin2, ixmax1, ixmax2
    double precision, intent(in) :: dx1,dx2 
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    double precision, intent(in) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision, intent(inout) :: dtnew
 
    !dtnew = ????????
 
  end subroutine get_new_dt

  subroutine source_cak(ixImin1, ixImin2, ixImax1, ixImax2,&
                          ixOmin1, ixOmin2, ixOmax1, ixOmax2,&
                          w, wCT, x, qdt)
    use mod_global_parameters

    integer :: i, j, p
    integer, intent(in) :: ixImin1, ixImin2, ixImax1, ixImax2
    integer, intent(in) :: ixOmin1, ixOmin2, ixOmax1, ixOmax2
    real(kind=dp), intent(in) :: qdt
    real(kind=dp), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp), intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    real(kind=dp), intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp) :: q(ixImin1:ixImax1,ixImin2:ixImax2,2)
    real(kind=dp) :: r, ddr, vr, P00, P11, P22, P01, P12, dx1, dx2, dr 
    real(kind=dp) :: vrI(2), dvdr, nu2_c, B, sigma, f, gL, gLx, gLy
    real(kind=dp) :: ke, Afac

    q(ixImin1:ixImax1,ixImin2:ixImax2,1) = wCT(ixImin1:ixImax1,ixImin2:ixImax2,mom(1))&
                                           /wCT(ixImin1:ixImax1,ixImin2:ixImax2,1)
    q(ixImin1:ixImax1,ixImin2:ixImax2,2) = wCT(ixImin1:ixImax1,ixImin2:ixImax2,mom(2))&
                                           /wCT(ixImin1:ixImax1,ixImin2:ixImax2,1)

    ke = (4.0*pi*UNIT_G*M_star*c*Edd)/L
    Afac = (1.0/(1.0-a))*((ke*L*Qfac)/(4.0*pi*c))

    do j=ixImin2+1,ixImax2-1
      do i=ixImin1+1,ixImax1-1

        r = sqrt(x(i,j,1)**2 + x(i,j,2)**2)

        if (r < 1.0) then

          w(i,j,rho_) = wCT(i,j,rho_)
          w(i,j,mom(1)) = wCT(i,j,mom(1))
          w(i,j,mom(2)) = wCT(i,j,mom(2))
          w(i,j,e_) = wCT(i,j,e_)

        else if (r > 1.0) then

          do p = 1,2

            dx1 = x(i+1,j,1) - x(i,j,1)
            dx2 = x(i,j+1,2) - x(i,j,2)
            dr = sqrt(dx1**2 + dx2**2)

            if (p == 1) then
              ddr = -0.9*dr
            else
              ddr = 0.9*dr
            end if

            call Interpolate2D(ixImin1, ixImin2, ixImax1, ixImax2,&
                               i, j, q, x, r, ddr, vrI(p))
          end do

          vr = q(i,j,1)*x(i,j,1)/r + q(i,j,2)*x(i,j,2)/r
          P00 = vrI(1)
          P11 = vr
          P22 = vrI(2)
          P01 = (0.5_dp*abs(ddr)*P00 + 0.5_dp*abs(ddr)*P11)/abs(ddr)
          P12 = (0.5_dp*abs(ddr)*P11 + 0.5_dp*abs(ddr)*P22)/abs(ddr)
          dvdr = abs(((1.0_dp/12.0_dp)*P00) - ((2.0_dp/3.0_dp)*P01) + ((2.0_dp/3.0_dp)*P12) &
                 - ((1.0_dp/12.0_dp)*P22))/abs(0.5_dp*ddr)
          !dvdr = abs((vrI[1] - vrI[0])/(2.0*ddr))

          nu2_c = 1.0_dp - 1.0_dp/r**2
          B = wCT(i,j,rho_)*Qfac*c*ke
          sigma = (r/abs(vr))*dvdr-1.0_dp
          f = ( ((1.0_dp + sigma)**(1.0_dp + a) - (1.0_dp + sigma*nu2_c)**(1.0_dp + a))&
              /((1.0_dp + a)*(1.0_dp - nu2_c)*sigma*(1.0_dp + sigma)**a))
          gL = f*Afac*r**(-2)*(dvdr/B)**a

          gLx = gL*x(i,j,1)/r
          gLy = gL*x(i,j,2)/r

          vx = q(i,j,1) + gLx*qdt 
          vy = q(i,j,2) + gLy*qdt
          prs = cs**2 * wCT(i,j,rho_)

          eng = prs/(hd_gamma - 1.0_dp) + 0.5_dp*wCT(i,j,rho_)*(vx**2 + vy**2)

          !deng = 0.5_dp*wCT(i,j,rho_)*(vx**2 + vy**2)

          !print*, "change in vel=", q(i,j,1) + gLx*qdt, q(i,j,2) + gLy*qdt
          !print*, "dvdr=", dvdr, "f=", f, "sigma=", sigma, "B=", B, "nu2_c=", nu2_c, "qdt=", qdt
          !print*, "vrI(1)=", vrI(1), "vr=", vr, "vrI(2)=", vrI(2), "gL=", gL
          !print*, ""

          w(i,j,rho_) = wCT(i,j,rho_)
          w(i,j,mom(1)) = wCT(i,j,rho_)*vx
          w(i,j,mom(2)) = wCT(i,j,rho_)*vy
          !w(i,j,p_) = wCT(i,j,p_)
          w(i,j,e_) = eng

          !print*, "vx=", vx, "vy=", vy, "eng=", eng, "rho=",w(i,j,rho_), "prs=",wCT(i,j,p_), "old_e=",wCT(i,j,e_) 
          

        end if

      end do
    end do

  end subroutine source_cak

  subroutine Interpolate2D(ixImin1, ixImin2, ixImax1, ixImax2,&
                           i, j, q, x, r, ddr, vrI)

    !       _____________________________
    !  j+1 |              |              | 
    !      |              |              |
    !      |              |     *        |
    !      |              |              |
    !      |              |              |
    !    j |______________|______________|
    !      |              |              | 
    !      |              |              |
    !      |              |              |
    !      |              |              |
    !      |              |              |
    !  j-1 |______________|______________|
    !  
    !    i - 1          i            i + 1
    !  
    !  
    !  yb 3               4
    !      |```|`````````|
    !      |   |         |
    !   yI |___|_________|
    !      |   |         |
    !      |   |         |
    !      |___|_________|
    !  ya 1    xI         2
    !      xa            xb
    !
    ! The interpolation points are always between xa, xb and ya, yb.
    !       
    ! ddr is the distance forwards (backwards) from the point 
    ! that the velocity gradient is needed to give the place 
    ! to perform the interpolation.
    !
    !  r -> r±ddr
    !
    ! |'''''''''''''''''''''''''|
    ! |                         |
    ! |                         |
    ! |      * r+ddr            |
    ! |     /|                  |
    ! |    / |                  |
    ! |   /  | r+ddr*cos(theta) |
    ! |  /   |                  |
    ! | /    |                  |
    ! |/_____|__________________|
    ! r     r+ddr*sin(theta)
    !
    ! xI and yI are the interpolation componets. 
    ! 

    use mod_global_parameters

    real(kind=dp), intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    real(kind=dp), intent(in) :: q(ixImin1:ixImax1,ixImin2:ixImax2,2)
    real(kind=dp), intent(in) :: r, ddr
    integer, intent(in) :: ixImin1, ixImin2, ixImax1, ixImax2
    real(kind=dp), intent(inout) :: vrI ! final interpolated radial-velocity. 

    integer :: u, s, i, j
    real(kind=dp) :: Ntot ! total volume of interpolation "space".
    real(kind=dp) :: N(4) ! Areas used to waight nearest neighbours (NN).
    real(kind=dp) :: V(2, 4) ! Array to hold velocity componants at NN.
    real(kind=dp) :: xa=0.0_dp, xb=0.0_dp ! bracketing x values.
    real(kind=dp) :: ya=0.0_dp, yb=0.0_dp ! bracketing y values.
    real(kind=dp) :: theta, xI, yI

    integer :: tag = 0

    N(:) = 0.0_dp
    V(:, :) = 0.0_dp
    theta = atan2(x(i,j,1), x(i,j,2)) ! convert to polar from Cartesian.
    xI = (r+ddr)*sin(theta)
    yI = (r+ddr)*cos(theta)

    !
    ! Eath of the following if statments checks which quadrent 
    ! the interpolation point is in and gets the components 
    ! of the velocity and the bracketing x and y values.
    !
    if (xI > x(i,j,1) .and. yI > x(i,j,2)) then
      tag = 1
      xa = x(i,j,1)
      xb = x(i+1,j,1)
      ya = x(i,j,2)
      yb = x(i,j+1,2)
      do u = 1, 2
        V(u, 1) = q(i,j,u)
        V(u, 2) = q(i+1,j,u)
        V(u, 3) = q(i,j+1,u)
        V(u, 4) = q(i+1,j+1,u)
      end do
    else if (xI < x(i,j,1) .and. yI > x(i,j,2)) then
      tag = 2
      xa = x(i-1,j,1)
      xb = x(i,j,1)
      ya = x(i,j,2)
      yb = x(i,j+1,2)
      do u = 1, 2
        V(u, 1) = q(i-1,j,u)
        V(u, 2) = q(i,j,u)
        V(u, 3) = q(i-1,j+1,u)
        V(u, 4) = q(i,j+1,u)
      end do
    else if (xI < x(i,j,1) .and. yI < x(i,j,2)) then
      tag = 3
      xa = x(i-1,j,1)
      xb = x(i,j,1)
      ya = x(i,j-1,2)
      yb = x(i,j,2)
      do u = 1, 2
        V(u, 1) = q(i-1,j-1,u)
        V(u, 2) = q(i,j-1,u)
        V(u, 3) = q(i-1,j,u)
        V(u, 4) = q(i,j,u)
      end do
    else if (xI > x(i,j,1) .and. yI < x(i,j,2)) then
      tag = 4
      xa = x(i,j,1)
      xb = x(i+1,j,1)
      ya = x(i,j-1,2)
      yb = x(i,j,2)
      do u = 1, 2
        V(u, 1) = q(i,j-1,u)
        V(u, 2) = q(i+1,j-1,u)
        V(u, 3) = q(i,j,u)
        V(u, 4) = q(i,j+1,u)
      end do
    end if

    ! Find total volume.
    N(1) = (xb - xI)*(yb - yI)
    N(2) = (xI - xa)*(yb - yI)
    N(3) = (xb - xI)*(yI - ya)
    N(4) = (xI - xa)*(yI - ya)
    Ntot = sum(N)
    ! Normalise volumes by total.
    N(:) = N(:)/Ntot
    vrI = (sum(V(1, :)*N(:))*x(i,j,1) + sum(V(2, :)*N(:))*x(i,j,2))/r

    !if (isnan(vrI)) then
    !  print*, "ddr=", ddr
    !  print*, "x(i,j,1)=", x(i,j,1), "x(i,j,2)=", x(i,j,2)
    !  print*, "xa=", xa, "xI=", xI, "xb=", xb
    !  print*, "ya=", ya, "yI=", yI, "yb=", yb
    !  print*, "V(1, :)=", V(1, :)
    !  print*, "V(2, :)=", V(2, :)
    !  print*, "tag =", tag, "vrI =", vrI, "Ntot=", Ntot, "N(:)=", N(:)
    !end if

  end subroutine Interpolate2D

end module mod_usr
