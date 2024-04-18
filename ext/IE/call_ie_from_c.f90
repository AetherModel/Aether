
! These codes are used to call the ionospheric electrodynamics 
! library developed for GITM from C++
! Developed by A. Ridley at University of Michigan, October, 2023

! Initialize the Ionospheric Electrodynamics Library:
! It would be nice to feed the path into this, but I am not going
! to worry about character arrays at this point.

subroutine convert_int_to_string(inArray, outString, length)
  implicit none
  integer, intent(in) :: length
  integer, dimension(length), intent(in) :: inArray
  character(length), intent(out) :: outString

  integer :: i
  do i = 1, length
    if (inArray(i) > 1) then
      outString(i:i) = achar(inArray(i))
    else
      outString(i:i) = ' '
    endif
  enddo

end subroutine convert_int_to_string


subroutine ie_init_library(ieDir_intArray, &
                          eField_intArray, &
                          aurora_intArray, &
                          iError) bind(C, name = 'ie_init_library')

    use ModCharSize
    use ModErrors
    implicit none
    integer, dimension(iCharLenIE_), intent(in) :: ieDir_intArray
    integer, dimension(iCharLenIE_), intent(in) :: eField_intArray
    integer, dimension(iCharLenIE_), intent(in) :: aurora_intArray
    integer, intent(out) :: iError

    character (len=iCharLenIE_) :: ieDir
    character (len=iCharLenIE_) :: eField
    character (len=iCharLenIE_) :: aurora
    character (len=iCharLenIE_), dimension(100) :: StringInputLines

    iError = 0;

    call convert_int_to_string(ieDir_intArray, ieDir, iCharLenIE_)
    call convert_int_to_string(eField_intArray, eField, iCharLenIE_)
    call convert_int_to_string(aurora_intArray, aurora, iCharLenIE_)

    StringInputLines(1) = "#BACKGROUND"
    StringInputLines(2) = ieDir
    StringInputLines(3) = eField
    StringInputLines(4) = aurora
    StringInputLines(5) = "unknown"  
    StringInputLines(6) = ""
    StringInputLines(7) = "#DEBUG"
    StringInputLines(8) = "0"
    StringInputLines(9) = "0"
    StringInputLines(10) = ""
    StringInputLines(11) = "#FIXTILT"
    StringInputLines(12) = "F"
    StringInputLines(13) = ""
    StringInputLines(14) = "#END"
    StringInputLines(15) = ""

    call EIE_set_inputs(StringInputLines)  
    call EIE_Initialize(iError)

    if (iError > 0) then
      write(*,*) 'Error in ie_init_library : ', iError, " -> ", &
        cErrorCodes(iError)
    endif

end subroutine ie_init_library


!----------------------------------------------------------------------

subroutine ie_set_nxs(nXsIn) bind(C, name = 'ie_set_nxs')
  use ModEIE_Interface
  implicit none
  integer, intent(in) :: nXsIn
  if (nXsIn /= IOi_NeednMLTs) IOi_ResetGrid = .true.
  IOi_NeednMLTs = nXsIn
end subroutine ie_set_nxs

!----------------------------------------------------------------------

subroutine ie_set_nys(nYsIn) bind(C, name = 'ie_set_nys')
  use ModEIE_Interface
  implicit none
  integer, intent(in) :: nYsIn
  if (nYsIn /= IOi_NeednLats) IOi_ResetGrid = .true.
  IOi_NeednLats = nYsIn
end subroutine ie_set_nys

!----------------------------------------------------------------------

subroutine ie_set_time(TimeIn) bind(C, name = 'ie_set_time')
  use ModEIE_Interface
  implicit none
!  Integer, parameter :: Real8_ = selected_real_kind(14,200)
  real (kind=Real8_), intent(in) :: TimeIn
  IOd_NeedTime = TimeIn
end subroutine ie_set_time

!----------------------------------------------------------------------

subroutine ie_set_imfbz(IMFBzIn) bind(C, name = 'ie_set_imfbz')
  use ModEIE_Interface
  implicit none
  real, intent(in) :: IMFBzIn
  IOr_NeedIMFBz = IMFBzIn
end subroutine ie_set_imfbz

!----------------------------------------------------------------------

subroutine ie_set_imfby(IMFByIn) bind(C, name = 'ie_set_imfby')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: IMFByIn
  IOr_NeedIMFBy = IMFByIn
end subroutine ie_set_imfby

!----------------------------------------------------------------------

subroutine ie_set_swv(SWVIn) bind(C, name = 'ie_set_swv')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: SWVIn
  IOr_NeedSWV = abs(SWVIn)
end subroutine 

!----------------------------------------------------------------------

subroutine ie_set_swn(SWNIn) bind(C, name = 'ie_set_swn')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: SWNIn
  IOr_NeedSWN = SWNIn
end subroutine ie_set_swn

!----------------------------------------------------------------------

subroutine ie_set_hpi(HPIIn) bind(C, name = 'ie_set_hpi')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: HPIIn
  IOr_NeedHPI = HPIIn
  IOr_NeedHPINorm = 2.09 * ALOG(HPIIn) * 1.0475
end subroutine ie_set_hpi

!----------------------------------------------------------------------

subroutine ie_set_hp_from_ae(AeIn) bind(C, name = 'ie_set_hp_from_ae')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: AeIn
  IOr_NeedHPI =  0.102 * AeIn + 8.953
  IOr_NeedHPINorm = 2.09 * ALOG(IOr_NeedHPI) * 1.0475
end subroutine ie_set_hp_from_ae

!----------------------------------------------------------------------

subroutine ie_set_ae(AeIn) bind(C, name = 'ie_set_ae')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: AeIn
  IOr_NeedAe = AeIn
end subroutine ie_set_ae

!----------------------------------------------------------------------

subroutine ie_set_au(AuIn) bind(C, name = 'ie_set_au')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: AuIn
  IOr_NeedAu = AuIn
end subroutine ie_set_au

!----------------------------------------------------------------------

subroutine ie_set_al(AlIn) bind(C, name = 'ie_set_al')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: AlIn
  IOr_NeedAL = AlIn
end subroutine ie_set_al

!----------------------------------------------------------------------

subroutine ie_set_kp(KpIn) bind(C, name = 'ie_set_kp')
!  use ModKind
  use ModEIE_Interface
  implicit none
  real, intent(in) :: KpIn
  IOr_NeedKp = KpIn
end subroutine ie_set_kp

!----------------------------------------------------------------------

subroutine ie_set_north() bind(C, name = 'ie_set_north')
!  use ModKind
  use ModEIE_Interface
  implicit none
  IOl_IsNorth = .true.
end subroutine ie_set_north

!----------------------------------------------------------------------

subroutine ie_set_south() bind(C, name = 'ie_set_south')
!  use ModKind
  use ModEIE_Interface
  implicit none
  IOl_IsNorth = .false.
end subroutine ie_set_south

!----------------------------------------------------------------------

subroutine ie_set_mlts(MLTsIn, iError) bind(C, name = 'ie_set_mlts')

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(IOi_NeednMLTs, IOi_NeednLats), intent(in) :: MLTsIn

  integer :: i, j

  iError = 0
  if (IOi_ResetGrid .and. allocated(IOr2_NeedMLTs)) deallocate(IOr2_NeedMLTs)
  if (.not. allocated(IOr2_NeedMLTs)) then
    allocate(IOr2_NeedMLTs(IOi_NeednMLTs, IOi_NeednLats), stat=iError)
    if (iError /= 0) then
        iError = ecAllocationError_
        return
    endif
  endif
  do i=1,IOi_NeednMLTs
    do j=1,IOi_NeednLats
      IOr2_NeedMLTs(i, j) = mod((MLTsIn(i, j)+24.0),24.0)
    enddo
  enddo

end subroutine ie_set_mlts

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine ie_set_lats(LatsIn, iError) bind(C, name = 'ie_set_lats')

  use ModErrors
  use ModEIE_Interface

  implicit none

  integer, intent(out) :: iError
  real, dimension(IOi_NeednMLTs, IOi_NeednLats), intent(in) :: LatsIn

  integer :: i, j

  iError = 0
  if (IOi_ResetGrid .and. allocated(IOr2_NeedLats)) deallocate(IOr2_NeedLats)
  if (.not. allocated(IOr2_NeedLats)) then
    allocate(IOr2_NeedLats(IOi_NeednMLTs, IOi_NeednLats), stat=iError)
    if (iError /= 0) then
        iError = ecAllocationError_
        return
    endif
  endif
  do i=1,IOi_NeednMLTs
    do j=1,IOi_NeednLats
      IOr2_NeedLats(i, j) = LatsIn(i, j) * 180.0 / 3.1415927
    enddo
  enddo

end subroutine ie_set_lats

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine ie_update_grid(iError) bind(C, name = 'ie_update_grid')

    implicit none

    integer, intent(out) :: iError

    iError = 0

    call IO_Calc_Interpolation(iError)

end subroutine ie_update_grid

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine ie_get_potential(PotentialOut, iError) bind(C, name = 'ie_get_potential')
  
    use ModEIE_Interface

    implicit none

    integer, intent(out) :: iError
    real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: PotentialOut

    iError = 0

    call IO_GetPotential(PotentialOut, iError)

end subroutine ie_get_potential

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine ie_get_eflux(EfluxOut, iError) bind(C, name = 'ie_get_eflux')
  
    use ModEIE_Interface

    implicit none

    integer, intent(out) :: iError
    real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: EfluxOut

    iError = 0

    call IO_GetEFlux(EFluxOut, iError)

end subroutine ie_get_eflux


!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine ie_get_avee(AveEOut, iError) bind(C, name = 'ie_get_avee')
  
    use ModEIE_Interface

    implicit none

    integer, intent(out) :: iError
    real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: AveEOut

    iError = 0

    call IO_GetAveE(AveEOut, iError)

end subroutine ie_get_avee


!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine ie_get_electron_diffuse_aurora(EFluxOut, AveEOut, iError) &
            bind(C, name = 'ie_get_electron_diffuse_aurora')
  
    use ModEIE_Interface

    implicit none

    integer, intent(out) :: iError
    real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: EFluxOut
    real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: AveEOut

    iError = 0

    call IO_GetElectronDiffuseAurora(EFluxOut, AveEOut, iError)

end subroutine ie_get_electron_diffuse_aurora




