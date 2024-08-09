program HomogeneousSpheroid_3D_Program
 
  use Basics
  use PSPFFT
  
  implicit none

  integer ( KDI ) :: &
    iVrbl, &
    iSource, &
    iNS, &   !-- iNumericalSolution
    iAS, &   !-- iAnalyticalSolution
    Error, &
    nProcs, &
    nProcsRoot
  integer ( KDI ), dimension ( 3 ) :: &
    nTotalCells
  real ( KDR ) :: &
    Pi, &
    GravitationalConstant = 1.0_KDR, &
    MyAnalyticalSum, MyDifferenceSum, &
    AnalyticalSum, DifferenceSum, &
    L1_Norm
  real ( KDR ), dimension(3) :: &
    CellWidth
  real ( KDR ), dimension ( :, :, : ), pointer :: &
    SS_3D   !-- SourceSolution_3D
  character ( 2 ) :: &
    nProcsString
  character ( LDF ) :: &
    DataDirectory
  type ( GridImageStreamForm ) :: &
    GIS
  type ( StructuredGridImageForm ) :: & 
    SGI    
  type ( Real_3D_Form ), dimension ( 3 ) :: &
    R_3D
  type ( PSPFFT_Form ), pointer :: &
    PS
  
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'HomogeneousSpheroid' )
  
  associate ( PH => PROGRAM_HEADER )
  
  nProcsRoot = PH % Communicator % Size ** ( 1.0_KDR / 3 ) + 0.5_KDR
  write(nProcsString, fmt='(i2)') PH % Communicator % Size
  
  DataDirectory = 'Data/Data_'//trim(adjustl(nProcsString)) // 'proc/'
  
  call Show ( DataDirectory, 'DataDirectory' )
  
  call GIS % Initialize &
         ( PH % Name, CommunicatorOption = PH % Communicator, &
           WorkingDirectoryOption = DataDirectory )
  call GIS % Open ( GIS % ACCESS_READ, NumberOption = 0 )
  
  call SGI % Initialize ( GIS )
  call SGI % Read ( )
  
  call GIS % Close ( )
  
  associate ( S => SGI % Storage ( 1 ) )
  do iVrbl = 1, S % nVariables
    if ( trim ( S % Variable ( iVrbl ) ) == 'Source' ) &
      iSource = iVrbl
    if ( trim ( S % Variable ( iVrbl ) ) == 'NumericalSolution' ) &
      iNS = iVrbl
    if ( trim ( S % Variable ( iVrbl ) ) == 'AnalyticalSolution' ) &
      iAS = iVrbl
  end do
  
  S % Value ( :, iNS ) &
    = S % Value ( :, iSource ) &
        * 4.0_KDR * CONSTANT % PI * CONSTANT % GRAVITATIONAL
  
  CellWidth ( 1 ) = SGI % NodeCoordinate_1 ( 2 ) - SGI % NodeCoordinate_1 ( 1 )
  CellWidth ( 2 ) = SGI % NodeCoordinate_2 ( 2 ) - SGI % NodeCoordinate_2 ( 1 )
  CellWidth ( 3 ) = SGI % NodeCoordinate_3 ( 2 ) - SGI % NodeCoordinate_3 ( 1 )
  
  !-- Rank-remaaped S % Value to 3D, then copy it to R_3D % Value
  SS_3D ( 1 : SGI % nCells ( 1 ), &
          1 : SGI % nCells ( 2 ), &
          1 : SGI % nCells ( 3 ) ) => S % Value ( :, iNS )
          
  call R_3D ( 3 ) % Initialize ( SGI % nCells )
  call Copy ( SS_3D, R_3D ( 3 ) % Value )

  nTotalCells = SGI % nCells * nProcsRoot
  
  !-- FIXME: This needs the ported version to work
  
  call Create ( PS, CellWidth, nTotalCells, PH % Communicator % Handle, &
                'INFO_5' )
  
  call Show ( R_3D ( 3 ) % Value ( 1 : 10, 5 : 10, 13 : 14 ), '<<<< Source ' )
  call Solve ( PS, R_3D ( 3 : 3 ) )
  call Show ( R_3D ( 3 ) % Value ( 1 : 10, 5 : 10,  13 : 14  ), '<<<< Solution ' )
  
  !-- Copy back the solution
  call Copy ( R_3D ( 3 ) % Value, SS_3D )
  
  !-- Calculate L1 norm of relative error between Analytical and Numerical
  !-- Solution
  MyAnalyticalSum = sum ( abs ( S % Value ( :, iAS ) ) )
  MyDifferenceSum = sum ( abs ( S % Value ( :, iNS ) &
                                -  S % Value ( :, iAS ) ) )
  !-- FIXME: NEED REDUCTION for MultiProc
  AnalyticalSum = MyAnalyticalSum
  DifferenceSum = MyDifferenceSum
  L1_Norm = DifferenceSum / AnalyticalSum
  
  call Show ( L1_Norm, 'L1_Norm Error' )
  
  end associate   !-- S
  
  end associate   !-- PH

  deallocate ( PROGRAM_HEADER )
  
end program HomogeneousSpheroid_3D_Program
