!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModEIE_Interface

  use ModKind
  use ModCharSize
  
  real, allocatable, dimension(:,:,:) :: EIEr3_HaveLats, EIEr3_HaveMLTs
  real, allocatable, dimension(:,:,:) :: EIEr3_HavePotential
  real, allocatable, dimension(:,:,:) :: EIEr3_HaveEFlux
  real, allocatable, dimension(:,:,:) :: EIEr3_HaveAveE
  real, allocatable, dimension(:,:,:) :: EIEr3_HaveIonEFlux
  real, allocatable, dimension(:,:,:) :: EIEr3_HaveIonAveE

  real (Real8_)       :: EIEd_CurrentTime
  integer             :: EIEi_HavenLats
  integer             :: EIEi_HavenMLTs
  integer             :: EIEi_HavenBLKs
  integer             :: EIEi_HavenTimes

  real (Real8_)               :: IOd_NeedTime = -1.0e32
  real, allocatable, dimension(:,:) :: IOr2_NeedLats, IOr2_NeedMLTs
  real, allocatable, dimension(:,:) :: IOr2_NeedPotential
  real, allocatable, dimension(:,:) :: IOr2_NeedEFlux
  real, allocatable, dimension(:,:) :: IOr2_NeedAveE
  real, allocatable, dimension(:,:) :: IOr2_NeedIonEFlux
  real, allocatable, dimension(:,:) :: IOr2_NeedIonAveE
  integer                           :: IOi_NeednLats = 0
  integer                           :: IOi_NeednMLTs = 0
  integer                           :: IOi_NeednTimes
  logical                           :: IOi_ResetGrid = .true.
  integer, allocatable, dimension(:,:,:) :: IOi3_InterpolationIndices
  real, allocatable, dimension(:,:,:)    :: IOr3_InterpolationRatios
  real :: IOr_NeedIMFBz   = -1.0e32
  real :: IOr_NeedIMFBy   = -1.0e32 
  real :: IOr_NeedSWV     = -1.0e32 
  real :: IOr_NeedSWN     = -1.0e32 
  real :: IOr_NeedHPI     = -1.0e32 
  real :: IOr_NeedHPINorm = -1.0e32 
  real :: IOr_NeedKp      = -1.0e32 
  logical :: IOl_IsNorth  = .true.
  real :: IOr_NeedAE     = -1.0e32 
  real :: IOr_NeedAU     = -1.0e32 
  real :: IOr_NeedAL     = -1.0e32 

!  real (Real8_)               :: UAd_NeedTime = -1.0e32
  real, allocatable, dimension(:,:) :: UAr2_NeedLats, UAr2_NeedMLTs
  real, allocatable, dimension(:,:) :: UAr2_NeedPotential
  real, allocatable, dimension(:,:) :: UAr2_NeedEFlux
  real, allocatable, dimension(:,:) :: UAr2_NeedAveE
  real, allocatable, dimension(:,:) :: UAr2_NeedIonEFlux
  real, allocatable, dimension(:,:) :: UAr2_NeedIonAveE
  integer                           :: UAi_NeednLats
  integer                           :: UAi_NeednMLTs
  integer                           :: UAi_NeednTimes
  integer, allocatable, dimension(:,:,:) :: UAi3_InterpolationIndices
  real, allocatable, dimension(:,:,:)    :: UAr3_InterpolationRatios
  logical :: UAl_IsNorth  = .true.

  integer                           :: iDebugLevel = 0
  integer                           :: iProc = 0

  integer, parameter                :: EIE_Closest_     = 1
  integer, parameter                :: EIE_After_       = 2
  integer, parameter                :: EIE_Interpolate_ = 3

  character (len=iCharLenIE_) :: EIE_NameOfEFieldModel
  character (len=iCharLenIE_) :: EIE_NameOfAuroralModel
  character (len=iCharLenIE_) :: EIE_NameOfSolarModel
  character (len=iCharLenIE_) :: EIE_NameOfModelDir

  logical :: UAl_UseGridBasedEIE
  logical :: UseGridBasedEIE

  logical :: IsFixedTilt = .false.

end module ModEIE_Interface
