program HomogeneousSpheroid_3D_Program

  use PSPFFT
  use IO
  
  implicit none
  
  include 'mpif.h'

  integer :: &
    Error, &
    nProcs, &
    nProcsRoot
  integer, dimension(3) :: &
    nTotalCells
  real(KR) :: &
    Pi, &
    GravitationalConstant = 1.0_KR, &
    MyAnalyticalSum, MyDifferenceSum, &
    AnalyticalSum, DifferenceSum, &
    L1_Norm
  real(KR), dimension(3) :: &
    CellWidth
  character(2) :: &
    nProcsString
  character(LL), dimension(3) :: &
    VariableName
  character(LF) :: &
    DataDirectory
  type(IO_Form), pointer :: &
    FileRead, &
    FileWrite
  type(QuadMeshForm), pointer :: &
    Mesh
  type(ArrayReal_3D_Base), dimension(3) :: &
    AR
  type(VariableArrayReal_3D_GroupBase), dimension(1) :: &
    VAR
  type(PSPFFT_Form), pointer :: &
    PS
  
  call MPI_INIT(Error)
  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, Error)
  write(nProcsString, fmt='(i2)') nProcs
  nProcsRoot = nProcs**(1.0_KR/3) + 0.5_KR
  DataDirectory = 'Data/Data_'//trim(adjustl(nProcsString)) // 'proc/'
  
  VariableName = (/'Source                        ', &
                   'AnalyticalSolution            ', & 
                   'NumericalSolution             '/)
  
  !-- VAR(1) uses AR as its storage, which at this point is unallocated 
  call Initialize(VAR(1), AR, VariableOption = VariableName)
  
  !-- Create IO objects read data file from, and write the solution to
  call Create( &
         FileRead, 'HomogeneousSpheroid_3D', MPI_COMM_Option = MPI_COMM_WORLD)
  call Open( &
         FileRead, FileRead%READ_ACCESS, &
         WorkingDirectoryOption = trim(DataDirectory), NumberOption = 0)
  call Create( &
         FileWrite, 'HomogeneousSpheroid_3D', MPI_COMM_Option = MPI_COMM_WORLD)
  call Open( &
         FileWrite, FileWrite%WRITE_ACCESS, &
         WorkingDirectoryOption = trim(DataDirectory), &
         SeriesOption = .true., NumberOption = 1)
  
  !-- Create Mesh object then read mesh & variable. The mesh is uses VAR
  !   to store its value (which in turns points to AR as the real storage) 
  !   Storage in AR is allocated during the Read() method as requested by its 
  !   last argument.
  call Create(Mesh, VAR)
  call Modify(Mesh, Directory = '/', Name = '')
  call Read(Mesh, FileRead, AllocateStorageOption = .true.)
  
  !-- Copy the "Source" the "NumericalSolution" because our solver solves
  !   the Poisson equation in-place:
  !-- VAR(1)%Data(1)%Data is a rank 3 array with the index i,j,k
  !   gives the value of the source on mesh grid i,j,k
  
  Pi = acos(-1.0_KR)
  VAR(1)%Data(3)%Data &
    = VAR(1)%Data(1)%Data * 4.0_KR * Pi * GravitationalConstant
  
  !-- Create Poisson Solver Object
  CellWidth(1) = Mesh%NodeCoordinate_1(2) - Mesh%NodeCoordinate_1(1)
  CellWidth(2) = Mesh%NodeCoordinate_2(2) - Mesh%NodeCoordinate_2(1)
  CellWidth(3) = Mesh%NodeCoordinate_3(2) - Mesh%NodeCoordinate_3(1)
  nTotalCells = Mesh%VariableArrayGroup_3D(1)%nData(1)%Data * nProcsRoot
  call Create(PS, CellWidth,  nTotalCells, MPI_COMM_WORLD)
  
  !-- Solve the Poisson Equations with Solve() subroutine, passing an array
  !   of sources as argument. Although currently only array of size 1 is 
  !   supported, this is done for future compatibility where multiple sources 
  !   may be solved at once. 
  !   Source is replaced by Solution upon subroutine exit 
  
  call Solve(PS, VAR(1)%Data(3:3))
  
  !-- Calculate L1 norm of relative error between Analytical Solution
  !-- and Numerical solution
  
  MyAnalyticalSum = sum(abs(VAR(1)%Data(2)%Data))
  MyDifferenceSum = sum(abs(VAR(1)%Data(3)%Data - VAR(1)%Data(2)%Data))
  call Reduce( &
         MyAnalyticalSum, FileRead%Communicator, REDUCTION_SUM, 0, AnalyticalSum)
  call Reduce( &
         MyDifferenceSum, FileRead%Communicator, REDUCTION_SUM, 0, DifferenceSum)
  L1_Norm = DifferenceSum / AnalyticalSum
  call Show(L1_Norm, 'L1_Norm Error')
  
  !-- Write the mesh with the Solution (i.e. potential) for the source
  
  call Write(Mesh, FileWrite)
  
  !-- Cleanup
  call Destroy(PS)
  
  call Destroy(Mesh)
  
  call Close(FileWrite)
  call Destroy(FileWrite)
  
  call Close(FileRead)
  call Destroy(FileRead)  
  
  call Finalize(VAR(1))
  
  call MPI_FINALIZE(Error)
  
end program HomogeneousSpheroid_3D_Program
