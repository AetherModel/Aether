
! These codes are used to call MSIS from C++
! Developed by A. Ridley at University of Michigan, April, 2023

! Initialize MSIS:
! It would be nice to feed the path into this, but I am not going
! to worry about character arrays at this point.

subroutine init_msis() bind(C, name = 'init_msis')

  use msis_init, only : msisinit
  
  implicit none

  call msisinit(parmpath='',parmfile='msis21.parm')
  
end subroutine init_msis

! Call MSIS for one time / location:
! The variables that we pass in could be one type if we want (like
! doubles or whatever), while those used to call MSIS are the specific
! variable type that MSIS needs. It just makes sure that there is
! nothing lost in translation.

subroutine call_msis_f(iYear, iDay, second, gLonDeg, gLatDeg, altKm, &
     f107in, f107ain, apin, &
     density, temperature) bind(C, name = 'call_msis_f')

  use msis_init, only : msisinit
  
  implicit none

  integer :: iYear, iDay, i
  real :: second, gLonDeg, gLatDeg, altKm, f107in, f107ain, apin
  real :: density(10), temperature(2)
  
  integer :: iyd, mass
  real(4) :: sec, alt, glat, glong, stl, f107a, f107, ap(7), apd
  real(4) :: d(10),t(2)

  iyd = iYear * 1000 + iDay
  sec = second
  glat = gLatDeg
  glong = gLonDeg
  alt = altKm
  stl = mod(glong + second/(3600.0*15.0), 24.0)
  f107 = f107in
  f107a = f107ain
  ap(1) = apin
  
  call gtd8d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)

  do i = 1, 10
     density(i) = d(i)
  enddo
  temperature(1) = t(1)
  temperature(2) = t(2)
  
end subroutine call_msis_f
