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
  real(kind=dp), parameter :: Q = 700.0_dp
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
                              (L/(c*c)) * ( Q * Edd * (1.0_dp - Edd)**(-1) )**(a_eff**(-1) - 1.0_dp)
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
    usr_set_B0 => background_field

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
        else if (r > 0.5 .and. r <= 1.0) then
          rho = M_dot/(4.0*pi*(cs/Cs_p))
          vx = 0.0
          vy = 0.0
          prs = cs*cs*rho/hd_gamma
        else if (r <= 0.5) then
          rho = M_dot/(4.0*pi*(cs/Cs_p))
          vx = 0.0
          vy = 0.0
          prs = cs*cs*rho/hd_gamma
        end if

        eng = prs/(rho*(hd_gamma - 1.0_dp)) + 0.5*rho*(vx**2 + vy**2)

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
          vx = 0.0
          vy = 0.0
          prs = cs*cs*rho/hd_gamma
          eng = prs/(rho*(hd_gamma - 1.0_dp)) + 0.5*rho*(vx**2 + vy**2)

          w(i,j,rho_) = rho
          w(i,j,mom(1)) = rho*vx
          w(i,j,mom(2)) = rho*vy
          w(i,j,e_) = eng

        else if (r <= 1.0 .and. r > 0.5) then

          rho = M_dot/(4.0*pi*(cs/Cs_p))
          vx = (cs/Cs_p)*x(i,j,1)/r
          vy = (cs/Cs_p)*x(i,j,2)/r
          vx = w(i,j,mom(1))/w(i,j,rho_)
          vy = w(i,j,mom(2))/w(i,j,rho_)
          prs = cs*cs*rho/hd_gamma
          eng = prs/(rho*(hd_gamma - 1.0_dp)) + 0.5*rho*(vx**2 + vy**2)

          w(i,j,rho_) = rho
          !w(i,j,mom(1)) = rho*vx
          !w(i,j,mom(2)) = rho*vy
          w(i,j,e_) = eng

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
          gravity_field(i, j, 1) = GM*x(i,j,1)/r
          gravity_field(i, j, 2) = GM*x(i,j,2)/r
        end if

      end do
    end do

  end subroutine gravity

  subroutine background_field(ixImin1,ixImin2,ixImax1,ixImax2,&
                              ixOmin1,ixOmin2,ixOmax1,ixOmax2,x,wB0)
    use mod_global_parameters

    integer, intent(in) :: ixImin1, ixImin2, ixImax1, ixImax2,&
                           ixOmin1, ixOmin2, ixOmax1, ixOmax2
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

  end subroutine specialset_B0

end module mod_usr
