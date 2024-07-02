!-- This module defines an FFT class that in this particular case provide 
!   a wrapper for FFTW library.

module FFT_FFTW__Base

  use Basics
  
  implicit none
  private
  
  include 'fftw3.f'
  
  public :: &
    Initialize, &
    Compute, &
    Finalize
  
  type, public :: FFT_FFTW_Base
    integer(KIB) :: &
      FORWARD   = FFTW_FORWARD, &
      BACKWARD  = FFTW_BACKWARD, &
      Direction = 0, &
      Handle    = 0
    integer(KI), dimension(:), allocatable :: &
      nData
    real(KR) :: &
      Normalization = 1.0_KR
    complex(KC), dimension(:), pointer :: &
      Data_1D => null()
    complex(KC), dimension(:,:), pointer :: &
      Data_2D => null()
    complex(KC), dimension(:,:,:), pointer :: &
      Data_3D => null()
    logical(KL) :: &
      LocalDataStorage = .false.
  end type FFT_FFTW_Base
  
  interface Initialize
    module procedure Initialize_FFT_FFTW
  end interface Initialize
  
  interface Compute
    module procedure Compute_FFT_FFTW
  end interface Compute
  
  interface Finalize
    module procedure Finalize_FFT_FFTW
  end interface Finalize
  
    private :: &
      Initialize_FFT_FFTW_Complex_1D, &
      Initialize_FFT_FFTW_Complex_Many_1D

contains
  

  subroutine Initialize_FFT_FFTW( &
               FFT_FFTW, Direction, Data_3D_Option, Data_2D_Option, &
               Data_1D_Option, DataShapeOption, nDataOption)
    
    !-- Initialize and setup FFTW "plan" for transform operation.

    type(FFT_FFTW_Base), intent(inout) :: &
      FFT_FFTW
    integer(KIB) :: &
      Direction
    complex(KC), dimension(:,:,:), intent(in), target, optional :: &
      Data_3D_Option
    complex(KC), dimension(:,:), intent(in), target, optional :: &
      Data_2D_Option
    complex(KC), dimension(:), intent(in), target, optional :: &
      Data_1D_Option
    integer(KI), dimension(:), intent(in), optional :: &
      DataShapeOption
    integer, intent(in), optional :: &
      nDataOption
    
    if(present(Data_1D_Option) .or. present(nDataOption))then
      call Initialize_FFT_FFTW_Complex_1D( &
             FFT_FFTW, Direction, Data_1D_Option, nDataOption)
    else if(present(Data_3D_Option) .or. present(Data_2D_Option) &
            .or. present(DataShapeOption)) &
    then
      call Initialize_FFT_FFTW_Complex_Many_1D( &
               FFT_FFTW, Direction, Data_3D_Option, Data_2D_Option, &
               DataShapeOption)
    end if
    
  end subroutine Initialize_FFT_FFTW


  subroutine Compute_FFT_FFTW(FFT_FFTW)
  
    type(FFT_FFTW_Base), intent(inout) :: &
      FFT_FFTW
    
    if(kind((1.0_KC,1.0_KC)) == kind((1.0,1.0)))then
      call SFFTW_EXECUTE(FFT_FFTW%Handle)
    elseif(kind((1.0_KC,1.0_KC)) == kind((1.0d0,1.0d0)))then
      call DFFTW_EXECUTE(FFT_FFTW%Handle)
    else
      !-- Long double not implemented
    end if

    FFT_FFTW%Normalization = sqrt(product(1.0_KR/FFT_FFTW%nData))

  end subroutine Compute_FFT_FFTW
  
  
  subroutine Finalize_FFT_FFTW(FFT_FFTW)
    
    !-- Free-up FFTW plan
  
    type(FFT_FFTW_Base), intent(inout) :: &
      FFT_FFTW
    
    if(KC == kind((1.0,1.0)))then
      call SFFTW_DESTROY_PLAN(FFT_FFTW%Handle)
    else if(KC == kind((1.0d0,1.0d0)))then
      call DFFTW_DESTROY_PLAN(FFT_FFTW%Handle)
    else
      !-- Long double not implemented
    end if

    if(associated(FFT_FFTW%Data_3D))then
      if(FFT_FFTW%LocalDataStorage)then
        deallocate(FFT_FFTW%Data_3D)
      else
        FFT_FFTW%Data_3D => null()
      end if
    end if
    if(associated(FFT_FFTW%Data_2D))then
      if(FFT_FFTW%LocalDataStorage)then
        deallocate(FFT_FFTW%Data_2D)
      else
        FFT_FFTW%Data_2D => null()
      end if
    end if
    if(associated(FFT_FFTW%Data_1D))then
      if(FFT_FFTW%LocalDataStorage)then
        deallocate(FFT_FFTW%Data_1D)
      else
        FFT_FFTW%Data_1D => null()
      end if
    end if
    FFT_FFTW%LocalDataStorage = .false.

    deallocate(FFT_FFTW%nData)

    FFT_FFTW%Handle    = 0
    FFT_FFTW%Direction = 0

  end subroutine Finalize_FFT_FFTW


  subroutine Initialize_FFT_FFTW_Complex_1D( &
               FFT_FFTW, Direction, Data_1D_Option, nDataOption)
    
    type(FFT_FFTW_Base), intent(inout) :: &
      FFT_FFTW
    integer(KIB) :: &
      Direction
    complex(KC), dimension(:), intent(in), target, optional :: &
      Data_1D_Option
    integer, intent(in), optional :: &
      nDataOption
    
    complex(KC), dimension(:), pointer :: &
      DataHead

    FFT_FFTW%Direction = Direction
    
    allocate(FFT_FFTW%nData(1))
    if(present(Data_1D_Option))then
      FFT_FFTW%nData = size(Data_1D_Option)
    else if(present(nDataOption))then
      FFT_FFTW%nData = nDataOption
    end if

    if(present(Data_1D_Option))then
      FFT_FFTW%Data_1D => Data_1D_Option
      FFT_FFTW%LocalDataStorage = .false.
    else if(present(nDataOption))then
      allocate(FFT_FFTW%Data_1D(nDataOption))
      FFT_FFTW%LocalDataStorage = .true.
    end if

    DataHead => FFT_FFTW%Data_1D(1:1)

    if(kind(DataHead(1)) == kind((1.0,1.0)))then
      call SFFTW_PLAN_DFT_1D( &
             FFT_FFTW%Handle, FFT_FFTW%nData, DataHead, DataHead, Direction, &
             FFTW_ESTIMATE)
    else if(kind(DataHead(1)) == kind((1.0d0,1.0d0)))then
      call DFFTW_PLAN_DFT_1D( &
             FFT_FFTW%Handle, FFT_FFTW%nData, DataHead, DataHead, Direction, &
             FFTW_ESTIMATE)
    else
      !-- Long double not implemented
    end if
    
  end subroutine Initialize_FFT_FFTW_Complex_1D
  
  
  subroutine Initialize_FFT_FFTW_Complex_Many_1D( &
               FFT_FFTW, Direction, Data_3D_Option, Data_2D_Option, &
               DataShapeOption)
    
    !-- nElementsPerTransform = DataShape(1)

    type(FFT_FFTW_Base), intent(inout) :: &
      FFT_FFTW
    integer(KIB), intent(in) :: &
      Direction
    complex(KC), dimension(:,:,:), intent(in), target, optional :: &
      Data_3D_Option
    complex(KC), dimension(:,:), intent(in), target, optional :: &
      Data_2D_Option
    integer(KI), dimension(:), intent(in), optional :: &
      DataShapeOption
    
    integer(KI), parameter :: &
      N_TRANSFORM_DIMENSIONS = 1, &
      STRIDE = 1
    integer(KI) :: &
      nTransforms
    integer(KI), dimension(N_TRANSFORM_DIMENSIONS) :: &
      nElementsPerTransform, &
      nElements
    integer(KI), dimension(:), allocatable :: &
      DataShape
    complex(KC), dimension(:), pointer :: &
      DataHead
    
    FFT_FFTW%Direction = Direction
    
    if(present(Data_3D_Option))then
      allocate(DataShape(3))
      DataShape = shape(Data_3D_Option)
    else if(present(Data_2D_Option))then
      allocate(DataShape(2))
      DataShape = shape(Data_2D_Option)
    else if(present(DataShapeOption))then
      allocate(DataShape(size(DataShapeOption)))
      DataShape = DataShapeOption
    end if

    if(size(DataShape) > 1)then
      nTransforms = product(DataShape(2:))
    else
      nTransforms = 1
    end if
    nElementsPerTransform = DataShape(1)
    nElements = product(DataShape)
    
    allocate(FFT_FFTW%nData(N_TRANSFORM_DIMENSIONS))
    FFT_FFTW%nData = nElementsPerTransform
    
    if(present(Data_3D_Option))then
      FFT_FFTW%Data_3D => Data_3D_Option
      FFT_FFTW%LocalDataStorage = .false.
    else if(present(Data_2D_Option))then
      FFT_FFTW%Data_2D => Data_2D_Option
      FFT_FFTW%LocalDataStorage = .false.
    else if(size(DataShape) == 3)then
      allocate(FFT_FFTW%Data_3D(DataShape(1),DataShape(2),DataShape(3)))
      FFT_FFTW%LocalDataStorage = .true.
    else if(size(DataShape) == 2)then
      allocate(FFT_FFTW%Data_2D(DataShape(1),DataShape(2)))
      FFT_FFTW%LocalDataStorage = .true.
    end if

    select case(size(DataShape))
      case (2)
        DataHead => FFT_FFTW%Data_2D(1:1,1)
      case (3)
        DataHead => FFT_FFTW%Data_3D(1:1,1,1)
    end select

    if(kind(DataHead(1)) == kind((1.0,1.0)))then
      call SFFTW_PLAN_MANY_DFT( &
             FFT_FFTW%Handle, &
             N_TRANSFORM_DIMENSIONS, nElementsPerTransform, nTransforms, &
             DataHead, nElements, STRIDE, product(nElementsPerTransform), &
             DataHead, nElements, STRIDE, product(nElementsPerTransform), &
             Direction, FFTW_ESTIMATE)
    else if(kind(DataHead(1)) == kind((1.0d0,1.0d0)))then
      call DFFTW_PLAN_MANY_DFT( &
             FFT_FFTW%Handle, &
             N_TRANSFORM_DIMENSIONS, nElementsPerTransform, nTransforms, &
             DataHead, nElements, STRIDE, product(nElementsPerTransform), &
             DataHead, nElements, STRIDE, product(nElementsPerTransform), &
             Direction, FFTW_ESTIMATE)
    else
      !-- Long double not implemented
    end if

    deallocate(DataShape)

  end subroutine Initialize_FFT_FFTW_Complex_Many_1D
  
  
end module FFT_FFTW__Base
