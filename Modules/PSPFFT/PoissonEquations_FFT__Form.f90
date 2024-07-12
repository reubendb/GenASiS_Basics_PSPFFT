!-- This module defines a class that provides an abstraction of a Poisson
!   equation and compute its solution.

module PoissonEquations_FFT__Form

  use Basics
  use Multiply_Command
  use Transpose_Command
  use FFT_FFTW__Base, &
        FFT_Base => FFT_FFTW_Base
  use LaplacianIsolated_FFT__Form, &
        LaplacianForm => LaplacianIsolated_FFT_Form
  
  implicit none
  private
  
  public :: &
    Create, &
    Solve, &
    Destroy
  
  type, public :: PoissonEquations_FFT_Form
    type(CommunicatorForm), pointer :: &
      Communicator
    type(LaplacianForm), pointer :: &
      Laplacian
  end type PoissonEquations_FFT_Form
  
  interface Create 
   
    module subroutine Create_L ( Laplacian, Communicator, CellWidth, nCells )
     implicit none 
      type ( LaplacianForm ), pointer :: &
        Laplacian
      type ( CommunicatorForm ), pointer :: &
        Communicator
      real ( KDR ), dimension ( 3 ), intent ( in ) :: &
        CellWidth
      integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nCells
    end subroutine
  
    module subroutine Create_C ( Communicator, MPI_COMM, CommunicatorName )
     implicit none 
      type ( CommunicatorForm ), pointer :: &
        Communicator
      integer ( KDI ), intent ( in ) :: &
        MPI_COMM
      character ( * ), intent ( in ) :: &
        CommunicatorName
    end subroutine
  
  end interface Create
  
  interface Solve
    module procedure Solve_PE
  end interface Solve
  
  interface Destroy
   
   module subroutine Destroy_L ( Laplacian, DestroyHandleOption )
     implicit none 
      type ( LaplacianForm ), pointer :: &
        Laplacian
      logical ( KDL ), intent ( in ), optional :: &
        DestroyHandleOption
    end subroutine
  
   end interface Destroy

 !interface
 
  !module subroutine Create ( ) 
 
   !type ( PoissonEquations_FFT_Form ), Po

  !end subroutine
 
 !end interface 
 
contains


  subroutine Create_PE &
               ( PE, CellWidth, nCells, MPI_COMM, VerbosityOption )
    
    type ( PoissonEquations_FFT_Form ), pointer :: &
      PE
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      CellWidth
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nCells
    integer ( KDI ), intent ( in ) :: &
      MPI_COMM
    character ( * ), intent ( in ), optional :: &
      VerbosityOption
  
    character ( LDL ) :: &
      Verbosity
    
    allocate ( PE )
    call Create ( PE % Communicator, MPI_COMM, 'PSPFFT_Communicator' )
    
    Verbosity = 'INFO_1'
    if ( present ( VerbosityOption ) ) Verbosity = VerbosityOption
    call CONSOLE % SetDisplayRank ( 0 )
    call CONSOLE % SetVerbosity ( Verbosity )
    call CONSOLE % Initialize ( ProcessRank = PE % Communicator % Rank )
    
    call Create ( PE % Laplacian, PE % Communicator, CellWidth, nCells )
    
  end subroutine Create_PE
  
  
  subroutine Solve_PE ( PE, SourceSolution )
     
    !-- Solve the Poisson equation by multiplyin the Green's function with
    !   the Source in the Fourier space, then transform the resulting 
    !   potential back to the real space. 
    
    type( PoissonEquations_FFT_Form ), intent( inout ) :: &
      PE
    type( Real_3D_Form ), dimension( : ), intent( inout ), target :: &
      SourceSolution
      
    integer ( KDI ) :: &
      iDim, &
      iSolution
    real ( KDR ) :: &
      Normalization
    logical ( KDL ) :: &
      Forward, &
      Backward 
    type ( LaplacianForm ), pointer :: &
      L
      
    L => PE % Laplacian
    L % InputOutput => SourceSolution
    
    Forward  = .true.
    Backward = .not.Forward
    
    call Load ( L )
      
    call Compute ( L % FFT_Forward ( 1 ) )
    
    L % FFT_Forward ( 2 ) % Data_3D &
      = Transpose ( L % Communicator_XY, Forward, L % FFT_Forward(1) % Data_3D )
    call Compute ( L % FFT_Forward(2) )
    
    L % FFT_Forward(3) % Data_3D &
      = Transpose ( L % Communicator_YZ, Forward, L % FFT_Forward(2) % Data_3D )
    call Compute ( L % FFT_Forward(3) )
    
    call Multiply ( L % FFT_Forward(3) % Data_3D, L % GreensFunction_Z )
    
    call Compute ( L % FFT_Backward(3) )
    L % FFT_Backward(2) % Data_3D &
      = Transpose ( L % Communicator_YZ, Backward, L % FFT_Backward(3) % Data_3D )

    call Compute ( L % FFT_Backward(2) )
    L % FFT_Backward(1) % Data_3D &
      = Transpose ( L % Communicator_XY, Backward, L % FFT_Backward(2) % Data_3D )

    call Compute ( L % FFT_Backward(1) )
    
    call Store ( L )

    Normalization = 1.0_KDR
    do iDim = 1, 3
      Normalization = Normalization * L % FFT_Forward(iDim) % Normalization
      Normalization = Normalization * L % FFT_Backward(iDim) % Normalization
    end do

    do iSolution = 1, size ( L % InputOutput )
      call Multiply ( L % InputOutput ( iSolution ) % Value, Normalization )
    end do
  
  end subroutine Solve_PE
  
  
  subroutine Destroy_PE ( PE )
    
    type ( PoissonEquations_FFT_Form ), pointer :: &
      PE
    
    call Destroy ( PE % Laplacian )
    call Destroy ( PE % Communicator, DestroyHandleOption = .false. )
    deallocate ( PE )
  
  end subroutine Destroy_PE
  

end module PoissonEquations_FFT__Form
