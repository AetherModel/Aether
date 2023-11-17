!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module ModAMIE_Interface

  use ModCharSize
  character (len=iCharLenIE_) :: AMIE_FileName
  integer :: AMIE_nLats, AMIE_nMlts, AMIE_nTimes

  ! For a single file
  real*4, allocatable,dimension(:) :: AMIE_Lats, AMIE_MLTs
  real*8, allocatable,dimension(:,:) :: AMIE_Time
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_Potential
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_PotentialY
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_Value
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_EFlux
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_AveE
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_IonEFlux
  real*4, allocatable,dimension(:,:,:,:) :: AMIE_IonAveE

  integer, parameter :: AMIE_Closest_     = 1
  integer, parameter :: AMIE_After_       = 2
  integer, parameter :: AMIE_Interpolate_ = 3

  integer :: AMIE_iDebugLevel = 0

  integer :: AMIE_South_ = 1
  integer :: AMIE_North_ = 2

  integer, parameter :: potential_ = 1
  integer, parameter :: eflux_     = 2
  integer, parameter :: avee_      = 3
  integer, parameter :: potentialy_ = 4
  integer, parameter :: ioneflux_ = 5
  integer, parameter :: ionavee_ = 6

end Module ModAMIE_Interface
