!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module ocean_rough_mod

!-----------------------------------------------------------------------

use          mpp_mod, only: input_nml_file

use       fms_mod, only: error_mesg, FATAL,  mpp_error, &
                         check_nml_error, mpp_pe, mpp_root_pe, &
                         write_version_number, stdlog
use constants_mod, only: grav, vonkarm

implicit none
private

public :: compute_ocean_roughness, fixed_ocean_roughness, &
          cal_z0_hwrf17, cal_zt_hwrf17, read_ocean_rough_scheme

!-----------------------------------------------------------------------
character(len=256) :: version = '$Id$'
character(len=256) :: tagname = '$Name$'
!-----------------------------------------------------------------------
!----- namelist -----

  real    :: roughness_init = 0.00044   ! not used in this version
  real    :: roughness_min  = 1.e-6
  real    :: charnock       = 0.032
  
  real    :: roughness_mom   = 5.8e-5
  real    :: roughness_heat  = 5.8e-5   ! was 4.00e-4
  real    :: roughness_moist = 5.8e-5
!  real, parameter :: zcoh1 = 0.0       ! Beljaars 1994 values
!  real, parameter :: zcoq1 = 0.0
! real, parameter :: zcoh1 = 1.4e-5
! real, parameter :: zcoq1 = 1.3e-4
  real            :: zcoh1 = 0.0 !miz
  real            :: zcoq1 = 0.0 !miz
  real            :: v10m  = 32.5 !jhc
  real            :: v10n  = 17.5 !jhc
  logical :: do_highwind     = .false.
  logical :: do_cap40        = .false.

  character(len=32) :: rough_scheme = 'fixed'   ! possible values:
                                                !   'fixed'
                                                !   'charnock'
                                                !   'beljaars'

namelist /ocean_rough_nml/ roughness_init, roughness_heat,  &
                           roughness_mom,  roughness_moist, &
                           roughness_min,                   &
                           charnock,                        &
                           rough_scheme, do_highwind,       &!miz
                           do_cap40, zcoh1, zcoq1,          &!sjl
                           v10m, v10n                        !jhc


!-----------------------------------------------------------------------

  logical :: do_init = .true.

!-----------------------------------------------------------------------
! ---- constants ----

! ..... high wind speed - rough sea
  real, parameter :: zcom1 = 1.8e-2    ! Charnock's constant
! ..... low wind speed - smooth sea
  real, parameter :: gnu   = 1.5e-5
  real, parameter :: zcom2 = 0.11
  real, parameter :: zcoh2 = 0.40
  real, parameter :: zcoq2 = 0.62


contains

!#######################################################################

 subroutine compute_ocean_roughness ( ocean, u_star,  &
                                      rough_mom, rough_heat, rough_moist )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(in)  :: u_star(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:), rough_moist(:,:)

!-----------------------------------------------------------------------
!  computes ocean roughness for momentum using wind stress
!  and sets roughness for heat/moisture using namelist value
!-----------------------------------------------------------------------

   real, dimension(size(ocean,1),size(ocean,2)) :: ustar2, xx1, xx2, w10 !miz
   real, dimension(size(ocean,1),size(ocean,2)) :: ustar, xx3, u10n, z0, zt, z1,  &
                                                   alpha_v, reynolds_rough
   real:: zt1
   integer :: i, j, n, iter
   real :: ustar_min, m, b, u_max, rough_mom_init

   if (do_init) call ocean_rough_init


   if (trim(rough_scheme) == 'fixed') then

!  --- set roughness for momentum and heat/moisture ---

      call fixed_ocean_roughness ( ocean, rough_mom, rough_heat, &
                                          rough_moist )


!  --- compute roughness for momentum, heat, moisture ---

   else if (trim(rough_scheme) == 'beljaars' .or. &
            trim(rough_scheme) == 'charnock') then

      where (ocean)
          ustar2(:,:) = max(gnu*gnu,u_star(:,:)*u_star(:,:))          
          xx1(:,:) = gnu / sqrt(ustar2(:,:))
          xx2(:,:) = ustar2(:,:) / grav
      elsewhere
          rough_mom   = 0.0
          rough_heat  = 0.0
          rough_moist = 0.0
      endwhere

      if (trim(rough_scheme) == 'charnock') then
          where (ocean)
              rough_mom  (:,:) = charnock * xx2(:,:)
              rough_mom  (:,:) = max( rough_mom(:,:), roughness_min )
              rough_heat (:,:) = rough_mom  (:,:)
              rough_moist(:,:) = rough_mom  (:,:)
          endwhere
      else if (trim(rough_scheme) == 'beljaars') then
! --- SJL ---- High Wind correction following Moon et al 2007 ------
          if (do_highwind) then       !  Moon et al. formular
              do j=1,size(ocean,2)
                 do i=1,size(ocean,1)
                    if ( ocean(i,j) ) then
                      w10(i,j) = 2.458 + u_star(i,j)*(20.255-0.56*u_star(i,j))  ! Eq(7) Moon et al.
                      if ( w10(i,j) > 12.5 ) then
                           rough_mom(i,j) = 0.001*(0.085*w10(i,j) - 0.58)    ! Eq(8b) Moon et al.
! SJL mods: cap the growth of z0 with w10 up to 40 m/s
! z0 (w10=40) = 2.82E-3
                           if(do_cap40) rough_mom(i,j) = min( rough_mom(i,j), 2.82E-3)
                      else
                           rough_mom(i,j) = 0.0185/grav*u_star(i,j)**2  ! (8a) Moon et al.
                      endif
                           zt1 = min( 1., max(0., (w10(i,j)-v10n)/(v10m-v10n)) )
                           rough_moist(i,j) = zcoq1*zt1*xx2(i,j) + zcoq2 * xx1(i,j)
                           rough_heat (i,j) = zcoh1*zt1*xx2(i,j) + zcoh2 * xx1(i,j)
!                 --- lower limit on roughness? ---
                      rough_mom  (i,j) = max( rough_mom  (i,j), roughness_min )
                      rough_heat (i,j) = max( rough_heat (i,j), roughness_min )
                      rough_moist(i,j) = max( rough_moist(i,j), roughness_min )
                    endif
                 enddo
              enddo
! SJL -----------------------------------------------------------------------------------
          else
!     --- Beljaars scheme ---
          where (ocean)
              rough_mom  (:,:) = zcom1 * xx2(:,:) + zcom2 * xx1(:,:)
              rough_heat (:,:) = zcoh1 * xx2(:,:) + zcoh2 * xx1(:,:)
              rough_moist(:,:) = zcoq1 * xx2(:,:) + zcoq2 * xx1(:,:)
!             --- lower limit on roughness? ---
              rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )
              rough_heat (:,:) = max( rough_heat (:,:), roughness_min )
              rough_moist(:,:) = max( rough_moist(:,:), roughness_min )
          endwhere
          endif
      endif

   else if (trim(rough_scheme) == 'coare3.5') then    !IH -ref Edson et al, JPO, Aug 2013 

      !if rough_mom is available in input (rough_mon_in) then set rough_mom1_init = rough_mom_in
      rough_mom_init = 1.e-03
      ustar_min = 1.e-05 ! added by IH just in case ustar is unphysically small 
      m = 0.0017         !Edson 2013
      b = -0.005         !Edson 2013
      u_max = 18.0       !Edson 2013
      ustar(:,:)=u_star(:,:)
      where (ocean)
          ustar(:,:)  = max(ustar(:,:), ustar_min)  ! IH
	  ustar2(:,:) = ustar(:,:)*ustar(:,:)
          xx1(:,:)    = gnu/ustar(:,:)
          xx2(:,:)    = ustar2(:,:)/grav
	  xx3(:,:)    = ustar(:,:)/vonkarm 
      endwhere	  
      !if rough_mom is available in input then set z0(:,:) to this input value
      z0(:,:) = rough_mom_init
      iter = 5      ! IH - overkill
      do j=1,size(ocean,2)
     	do i=1,size(ocean,1)
           if ( ocean(i,j) ) then
              do n = 1, iter
	      	 u10n(i,j) = xx3(i,j)*log(10/z0(i,j))             ! "neutral" 10m wind
	    	 alpha_v(i,j) = m*min(u10n(i,j),u_max) + b;     ! Charnock coefficient, Edson 2013
	    	 z1(i,j) = zcom2*xx1(i,j) + alpha_v(i,j)*xx2(i,j) ! Edson 2013
	    	 z0(i,j) = z1(i,j);                               !IH -- iteration
	      enddo
           endif
        enddo
      enddo

      where (ocean)
      	  rough_mom  (:,:) = z0(:,:)
      	  rough_mom  (:,:) = max( rough_mom  (:,:), roughness_min )

	  reynolds_rough(:,:) = ustar(:,:)*rough_mom(:,:)/gnu
	  rough_heat (:,:)    = 5.5e-05*(reynolds_rough(:,:)**(-0.6))
	  rough_heat (:,:)    = min(1.1e-04, rough_heat(:,:))
	  rough_moist(:,:)    = rough_heat (:,:)
      elsewhere
          rough_mom   = 0.0
          rough_heat  = 0.0
          rough_moist = 0.0
      endwhere

   else if (trim(rough_scheme) == 'hwrf17') then

      rough_mom_init = 1.e-03
      ustar_min = 1.e-05
      ustar(:,:) = u_star(:,:)
      where (ocean)
          ustar(:,:)  = max(ustar(:,:), ustar_min)
          ustar2(:,:) = ustar(:,:)*ustar(:,:)
          xx1(:,:)    = gnu/ustar(:,:)
          xx2(:,:)    = ustar2(:,:)/grav
          xx3(:,:)    = ustar(:,:)/vonkarm
      endwhere

      z0(:,:) = rough_mom_init
      iter = 5

      do j = 1, size(ocean, 2)
        do i = 1, size(ocean, 1)
           if ( ocean(i,j) ) then
              do n = 1, iter
               u10n(i,j) = xx3(i,j) * log(10 / z0(i,j))
               call cal_z0_hwrf17(u10n(i,j), z0(i,j))
              enddo
              call cal_zt_hwrf17(u10n(i,j), zt(i,j))
              rough_mom  (i,j) = z0(i,j)
              rough_heat (i,j) = zt(i,j)
              rough_moist(i,j) = rough_heat(i,j)
           endif
        enddo
      enddo

      where (.not. ocean)
          rough_mom   = 0.0
          rough_heat  = 0.0
          rough_moist = 0.0
      endwhere

   else
      call mpp_error(FATAL, '==>Error from ocean_rough_mod(compute_ocean_roughness): '//&
            'Unknown roughness scheme (case sensitive): ' //trim(rough_scheme))
   endif

!-----------------------------------------------------------------------

 end subroutine compute_ocean_roughness

!#######################################################################

 subroutine fixed_ocean_roughness ( ocean, rough_mom, rough_heat, rough_moist )

 logical, intent(in)  :: ocean(:,:)
 real,    intent(out) :: rough_mom(:,:), rough_heat(:,:), rough_moist(:,:)

   if (do_init) call ocean_rough_init

    where (ocean)
       rough_mom   = roughness_mom
       rough_heat  = roughness_heat
       rough_moist = roughness_moist
    endwhere

 end subroutine fixed_ocean_roughness

!#######################################################################

 subroutine ocean_rough_init

   integer :: unit, ierr, io

!   ----- read and write namelist -----

  read (input_nml_file, nml=ocean_rough_nml, iostat=io)
  ierr = check_nml_error(io, 'ocean_rough_nml')

!------- write version number and namelist ---------

    if ( mpp_pe() == mpp_root_pe() ) then
         call write_version_number(version, tagname)
         unit = stdlog()
         write (unit,nml=ocean_rough_nml)
         write (unit,11)
    endif

!------ constants -----

    roughness_moist = max (roughness_moist, roughness_min)
    roughness_heat  = max (roughness_heat , roughness_min)
    roughness_mom   = max (roughness_mom  , roughness_min)

    do_init = .false.

11 format (/,'namelist option USE_FIXED_ROUGH is no longer supported', &
           /,'use variable ROUGH_SCHEME instead')

 end subroutine ocean_rough_init


!#######################################################################

 subroutine read_ocean_rough_scheme (ocean_rough_scheme)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! to read in the ocean roughness scheme used
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  character(len=32), intent(inout) :: ocean_rough_scheme
  integer :: ierr, io

  read (input_nml_file, nml=ocean_rough_nml, iostat=io)
  ierr = check_nml_error(io, 'ocean_rough_nml')
  ocean_rough_scheme = rough_scheme

 end subroutine read_ocean_rough_scheme

!#######################################################################

 subroutine cal_z0_hwrf17 (uref,z0)

! This subroutine is orginally from HWRF model (2017 version, znot_m_v6)
! Implemented here by Kun Gao & Baoqiang Xiang

! Calculate areodynamical roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Cd-U10 relationship from COARE V3.5 (Edson et al. 2013)
! For high winds, try to fit available observational data
!
! uref(m/s)   :   wind speed at 10-m height
! z0 (meter)  :   areodynamical roughness scale over water
!

    real, intent(in) :: uref
    real, intent(out):: z0
    real  :: p13, p12, p11, p10
    real  :: p25, p24, p23, p22, p21, p20
    real  :: p35, p34, p33, p32, p31, p30
    real  :: p40

    p13 = -1.296521881682694e-02
    p12 =  2.855780863283819e-01
    p11 = -1.597898515251717e+00
    p10 = -8.396975715683501e+00

    p25 =  3.790846746036765e-10
    p24 =  3.281964357650687e-09
    p23 =  1.962282433562894e-07
    p22 = -1.240239171056262e-06
    p21 =  1.739759082358234e-07
    p20 =  2.147264020369413e-05

    p35 =  1.840430200185075e-07
    p34 = -2.793849676757154e-05
    p33 =  1.735308193700643e-03
    p32 = -6.139315534216305e-02
    p31 =  1.255457892775006e+00
    p30 = -1.663993561652530e+01

    p40 =  4.579369142033410e-04

    if (uref >= 0.0 .and.  uref <= 6.5 ) then
      z0 = exp( p10 + p11*uref + p12*uref**2 + p13*uref**3)
    elseif (uref > 6.5 .and. uref <= 15.7) then
      z0 = p25*uref**5 + p24*uref**4 + p23*uref**3 + p22*uref**2 + p21*uref + p20
    elseif (uref > 15.7 .and. uref <= 53.0) then
      z0 = exp( p35*uref**5 + p34*uref**4 + p33*uref**3 + p32*uref**2 + p31*uref + p30 )
    elseif ( uref > 53.0) then
      z0 = p40
    else
      print*, 'Wrong input uref value:',uref
    endif
!
 end subroutine cal_z0_hwrf17

!#######################################################################

 subroutine cal_zt_hwrf17 (uref,zt)

! This subroutine is orginally from HWRF model (2017 version, znot_t_v6)
! Implemented here by Kun Gao & Baoqiang Xiang

! Calculate scalar roughness over water with input 10-m wind
! For low-to-moderate winds, try to match the Ck-U10 relationship from COARE algorithm
! For high winds, try to retain the Ck-U10 relationship of FY2015 HWRF
!
! uref(m/s)   :   wind speed at 10-m height
! zt(meter)   :   scalar roughness scale over water
!

    real, intent(in) :: uref
    real, intent(out):: zt

    real  :: p00
    real  :: p15, p14, p13, p12, p11, p10
    real  :: p25, p24, p23, p22, p21, p20
    real  :: p35, p34, p33, p32, p31, p30
    real  :: p45, p44, p43, p42, p41, p40
    real  :: p56, p55, p54, p53, p52, p51, p50
    real  :: p60

    p00 =  1.100000000000000e-04

    p15 = -9.144581627678278e-10
    p14 =  7.020346616456421e-08
    p13 = -2.155602086883837e-06
    p12 =  3.333848806567684e-05
    p11 = -2.628501274963990e-04
    p10 =  8.634221567969181e-04

    p25 = -8.654513012535990e-12
    p24 =  1.232380050058077e-09
    p23 = -6.837922749505057e-08
    p22 =  1.871407733439947e-06
    p21 = -2.552246987137160e-05
    p20 =  1.428968311457630e-04

    p35 =  3.207515102100162e-12
    p34 = -2.945761895342535e-10
    p33 =  8.788972147364181e-09
    p32 = -3.814457439412957e-08
    p31 = -2.448983648874671e-06
    p30 =  3.436721779020359e-05


    p45 = -3.530687797132211e-11
    p44 =  3.939867958963747e-09
    p43 = -1.227668406985956e-08
    p42 = -1.367469811838390e-05
    p41 =  5.988240863928883e-04
    p40 = -7.746288511324971e-03

    p56 = -1.187982453329086e-13
    p55 =  4.801984186231693e-11
    p54 = -8.049200462388188e-09
    p53 =  7.169872601310186e-07
    p52 = -3.581694433758150e-05
    p51 =  9.503919224192534e-04
    p50 = -1.036679430885215e-02

    p60 =  4.751256171799112e-05

    if (uref >= 0.0 .and. uref < 5.9 ) then
      zt = p00
    elseif (uref >= 5.9 .and. uref <= 15.4) then
      zt = p15*uref**5 + p14*uref**4 + p13*uref**3 + p12*uref**2 + p11*uref + p10
    elseif (uref > 15.4 .and. uref <= 21.6) then
      zt = p25*uref**5 + p24*uref**4 + p23*uref**3 + p22*uref**2 + p21*uref + p20
    elseif (uref > 21.6 .and. uref <= 42.2) then
      zt = p35*uref**5 + p34*uref**4 + p33*uref**3 + p32*uref**2 + p31*uref + p30
    elseif ( uref > 42.2 .and. uref <= 53.3) then
      zt = p45*uref**5 + p44*uref**4 + p43*uref**3 + p42*uref**2 + p41*uref + p40
    elseif ( uref > 53.3 .and. uref <= 80.0) then
      zt = p56*uref**6 + p55*uref**5 + p54*uref**4 + p53*uref**3 + p52*uref**2 + p51*uref + p50
    elseif ( uref > 80.0) then
      zt = p60
    else
      print*, 'Wrong input uref value:',uref
    endif
!
  end subroutine cal_zt_hwrf17

end module ocean_rough_mod

