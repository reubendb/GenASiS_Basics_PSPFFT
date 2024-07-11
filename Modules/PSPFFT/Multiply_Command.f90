!-- This module provides abstraction to multiply the contant arrays 
!   by overloading a generic "Multiply()" subroutine.

module Multiply_Command

  use Basics
    
  implicit none
  private
  
  public :: &
    Multiply
    
  interface Multiply
    module procedure MultiplyComplex_3D_InPlace
    module procedure MultiplyReal_3D_ScalarInPlace
  end interface Multiply
  
contains


  subroutine MultiplyComplex_3D_InPlace(A, B)
    
    complex(KDC), dimension(:,:,:), intent(inout) :: &
      A
    complex(KDC), dimension(:,:,:), intent(in) :: &
      B
    
    A = A * B
  
  end subroutine MultiplyComplex_3D_InPlace
  
  
  subroutine MultiplyReal_3D_ScalarInPlace(A, B)
    
    real(KDR), dimension(:,:,:), intent(inout) :: &
      A
    real(KDR), intent(in) :: &
      B
      
    A = A * B
  
  end subroutine MultiplyReal_3D_ScalarInPlace
  
  
end module Multiply_Command
