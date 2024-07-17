module QuadMeshSilo_Form

  use pspfft
  use InputOutputStreamSilo_Form, &
        InputOutputStreamForm => InputOutputStreamSiloForm
  
  implicit none
  private
  
  include 'silo_f9x.inc'
  
  public :: &
    Create, &
    Modify, &
    Write, &
    Read, &
    Destroy
  
  type, public :: QuadMeshSiloForm
    integer(KDI) :: &
      lName       = 0, &
      lDirectory  = 0, &
      nDimensions = 0
    integer(KDI), dimension(:), allocatable :: &
      nNodes
    real(KDR), dimension(:), allocatable :: &
      NodeCoordinate_1, &
      NodeCoordinate_2, &
      NodeCoordinate_3
    character(LDL) :: &
      Name              = '', &
      Directory         = '', &
      CoordinateLabel_1 = 'X                             ', &
      CoordinateLabel_2 = 'Y                             ', &
      CoordinateLabel_3 = 'Z                             '
    type(LabeledValueBase) :: &
      CoordinateUnit_1, &
      CoordinateUnit_2, &
      CoordinateUnit_3
    type(VariableArrayReal_3D_GroupBase), dimension(:), allocatable :: &
      VariableArrayGroup_3D
    !-- If we ever need to write 2D QuadMesh, we will need to define
    !   VariableArrayReal_2D_GroupBase type for VariableArrayGroup_2D
    !   and overload the Create() and Write() subroutines accordingly
  end type QuadMeshSiloForm
  
  interface Create
    module procedure Create_QM
  end interface Create
  
  interface Modify
    module procedure Modify_QM_3D_ForWriting
    module procedure Modify_QM_3D_ForReading
  end interface Modify
  
  interface Write
    module procedure Write_QM 
  end interface Write
  
  interface Read
    module procedure Read_QM
  end interface Read
  
  interface Destroy
    module procedure Destroy_QM
  end interface Destroy
  
    private :: &
      WriteMesh, &
      WriteVariableArrayGroup_3D, &
      ReadMesh, &
      ReadVariableArrayGroup_3D, &
      WriteMeshMultiMesh, &
      WriteVariableArrayGroupMultiMesh
  
contains

  
  subroutine Create_QM(QM, VariableArrayGroup)
    
    type(QuadMeshSiloForm), pointer :: &
      QM
    type(VariableArrayReal_3D_GroupBase), dimension(:), intent(in) :: &
      VariableArrayGroup
      
    integer(KDI) :: &
      iGroup
    
    allocate(QM)
    
    QM%nDimensions = 3
    
    allocate(QM%VariableArrayGroup_3D(size(VariableArrayGroup)))
    do iGroup = 1, size(VariableArrayGroup)
      call Initialize( &
             QM%VariableArrayGroup_3D(iGroup), VariableArrayGroup(iGroup))
    end do
    
  end subroutine Create_QM
  
  
  subroutine Modify_QM_3D_ForWriting( &
               QM, CoordinateUnit_1, CoordinateUnit_2, CoordinateUnit_3, &
               Directory, Name, NodeCoordinate_1, NodeCoordinate_2, &
               NodeCoordinate_3, CoordinateLabel_1_Option, &
               CoordinateLabel_2_Option, CoordinateLabel_3_Option)
               
    type(QuadMeshSiloForm), intent(inout) :: &
      QM
    type(LabeledValueBase), intent(in) :: &
      CoordinateUnit_1, &
      CoordinateUnit_2, &
      CoordinateUnit_3
    character(*), intent(in) :: &
      Directory, &
      Name
    real(KDR), dimension(:), intent(in) :: &
      NodeCoordinate_1, &
      NodeCoordinate_2, &
      NodeCoordinate_3
    character(*), intent(in), optional :: &
      CoordinateLabel_1_Option, &
      CoordinateLabel_2_Option, &
      CoordinateLabel_3_Option
      
    QM%lName      = len_trim(Name)
    QM%lDirectory = len_trim(Directory)
      
    allocate(QM%nNodes(QM%nDimensions))
    QM%nNodes &
      = (/ size(NodeCoordinate_1), size(NodeCoordinate_2), &
           size(NodeCoordinate_3) /)
    
    allocate(QM%NodeCoordinate_1(QM%nNodes(1)))
    allocate(QM%NodeCoordinate_2(QM%nNodes(2)))
    allocate(QM%NodeCoordinate_3(QM%nNodes(3)))
    QM%NodeCoordinate_1 = NodeCoordinate_1
    QM%NodeCoordinate_2 = NodeCoordinate_2
    QM%NodeCoordinate_3 = NodeCoordinate_3
    
    QM%Name      = Name
    QM%Directory = Directory
    if(trim(QM%Directory) == '/') QM%Directory = ''
    
    if(present(CoordinateLabel_1_Option)) &
      QM%CoordinateLabel_1 = CoordinateLabel_1_Option
    if(present(CoordinateLabel_2_Option)) &
      QM%CoordinateLabel_2 = CoordinateLabel_2_Option
    if(present(CoordinateLabel_1_Option)) &
      QM%CoordinateLabel_3 = CoordinateLabel_3_Option
      
    QM%CoordinateUnit_1 = CoordinateUnit_1
    QM%CoordinateUnit_2 = CoordinateUnit_2
    QM%CoordinateUnit_3 = CoordinateUnit_3
    
  end subroutine Modify_QM_3D_ForWriting
  
  
  subroutine Modify_QM_3D_ForReading( &
               QM, Directory, Name, CoordinateUnit_1_Option, &
               CoordinateUnit_2_Option, CoordinateUnit_3_Option)
               
    type(QuadMeshSiloForm), intent(inout) :: &
      QM
    character(*), intent(in) :: &
      Directory, &
      Name
    type(LabeledValueBase), intent(in), optional :: &
      CoordinateUnit_1_Option, &
      CoordinateUnit_2_Option, &
      CoordinateUnit_3_Option
    
    
    QM%lName      = len_trim(Name)
    QM%lDirectory = len_trim(Directory)
      
    allocate(QM%nNodes(QM%nDimensions))
    
    QM%Name      = Name
    QM%Directory = Directory
    if(trim(QM%Directory) == '/') QM%Directory = ''
    
    if(present(CoordinateUnit_1_Option) &
       .and. present(CoordinateUnit_2_Option) &
       .and. present(CoordinateUnit_3_Option))then
      QM%CoordinateUnit_1 = CoordinateUnit_1_Option
      QM%CoordinateUnit_2 = CoordinateUnit_2_Option
      QM%CoordinateUnit_3 = CoordinateUnit_3_Option
    end if
    
  end subroutine Modify_QM_3D_ForReading
  
  
  subroutine Write_QM(QM, IOS, TimeUnitOption, TimeOption, CycleNumberOption)
  
    type(QuadMeshSiloForm), intent(inout) :: &
      QM
    type(InputOutputStreamForm), intent(inout) :: &
      IOS
    type(LabeledValueBase), intent(in), optional :: &
      TimeUnitOption
    real(KDR), intent(in), optional :: &
      TimeOption
    integer(KDI), intent(in), optional :: &
      CycleNumberOption

    integer(KDI) :: &
      iGroup
    character(LDF) :: &
      WorkingDirectory
      
    if(IOS%MeshBlockHandle == -1) return
      
    WorkingDirectory = IOS%CurrentDirectory
    
    if(QM%lDirectory > 0) call MakeDirectory(IOS, QM%Directory)
    if(QM%lName > 0) call MakeDirectory(IOS, QM%Name)
    
    call WriteMesh(QM, IOS, TimeUnitOption, TimeOption, CycleNumberOption)
    
    if(QM%nDimensions == 3)then
      do iGroup = 1, size(QM%VariableArrayGroup_3D)
        call WriteVariableArrayGroup_3D( &
               QM, IOS, QM%VariableArrayGroup_3D(iGroup))
      end do
    end if

    call ChangeDirectory(IOS, WorkingDirectory)
    
  end subroutine Write_QM
  
  
  subroutine Read_QM(QM, IOS, AllocateStorageOption)
  
    type(QuadMeshSiloForm), intent(inout) :: &
      QM
    type(InputOutputStreamForm), intent(inout) :: &
      IOS
    logical(KDL), intent(in), optional :: &
      AllocateStorageOption

    integer(KDI) :: &
      iGroup
    character(LDF) :: &
      WorkingDirectory
      
    if(IOS%MeshBlockHandle == -1) return
      
    WorkingDirectory = IOS%CurrentDirectory
    
    call ChangeDirectory(IOS, QM%Directory)
    call ChangeDirectory(IOS, QM%Name)
    
    call ReadMesh(QM, IOS)
    
    if(QM%nDimensions == 3)then
      do iGroup = 1, size(QM%VariableArrayGroup_3D)
        call ReadVariableArrayGroup_3D( &
               QM, IOS, QM%VariableArrayGroup_3D(iGroup), &
               AllocateStorageOption)
      end do
    end if
    
    call ChangeDirectory(IOS, WorkingDirectory)
    
  end subroutine Read_QM
  
  
  subroutine Destroy_QM(QM)
    
    type(QuadMeshSiloForm), pointer :: &
      QM
    
    integer(KDI) :: &
      iGroup
      
    if(.not.associated(QM)) return
    
    if(allocated(QM%NodeCoordinate_3)) deallocate(QM%NodeCoordinate_3)
    if(allocated(QM%NodeCoordinate_2)) deallocate(QM%NodeCoordinate_2)
    if(allocated(QM%NodeCoordinate_1)) deallocate(QM%NodeCoordinate_1)
    
    if(allocated(QM%nNodes)) deallocate(QM%nNodes)
    
    do iGroup = 1, size(QM%VariableArrayGroup_3D)
      call Finalize(QM%VariableArrayGroup_3D(iGroup))
    end do
    deallocate(QM%VariableArrayGroup_3D)
    
    deallocate(QM)
  
  end subroutine Destroy_QM
  
  
  subroutine WriteMesh(QM, IOS, TimeUnitOption, TimeOption, CycleNumberOption)
  
    type(QuadMeshSiloForm), intent(inout) :: &
      QM
    type(InputOutputStreamForm), intent(inout) :: &
      IOS
    type(LabeledValueBase), intent(in), optional :: &
      TimeUnitOption
    real(KDR), intent(in), optional :: &
      TimeOption
    integer(KDI), intent(in), optional :: &
      CycleNumberOption

    integer(KDI) :: &
      nSiloOptions, &
      SiloOptionList, &
      Error
    real(KDR) :: &
      Time
      
    nSiloOptions = 3
    if(len_trim(QM%CoordinateUnit_1%Label) > 0) &
      nSiloOptions = nSiloOptions + 1
    if(len_trim(QM%CoordinateUnit_2%Label) > 0) &
      nSiloOptions = nSiloOptions + 1
    if(len_trim(QM%CoordinateUnit_3%Label) > 0) &
      nSiloOptions = nSiloOptions + 1
    if(present(TimeOption)) &
      nSiloOptions = nSiloOptions + 1
    if(present(CycleNumberOption)) & 
      nSiloOptions = nSiloOptions + 1
    
    Error = DBMKOPTLIST(nSiloOptions, SiloOptionList)
    
    Error = DBADDCOPT( &
              SiloOptionList, DBOPT_XLABEL, &
              trim(QM%CoordinateLabel_1), len_trim(QM%CoordinateLabel_1))
    Error = DBADDCOPT( &
              SiloOptionList, DBOPT_YLABEL, &
              trim(QM%CoordinateLabel_2), len_trim(QM%CoordinateLabel_2))
    Error = DBADDCOPT( &
              SiloOptionList, DBOPT_ZLABEL, &
              trim(QM%CoordinateLabel_3), len_trim(QM%CoordinateLabel_3))
    
    if(len_trim(QM%CoordinateUnit_1%Label) > 0) &
      Error = DBADDCOPT( &
                SiloOptionList, DBOPT_XUNITS, &
                trim(QM%CoordinateUnit_1%Label), &
                len_trim(QM%CoordinateUnit_1%Label))
    if(len_trim(QM%CoordinateUnit_2%Label) > 0) &   
      Error = DBADDCOPT( &
                SiloOptionList, DBOPT_YUNITS, &
                trim(QM%CoordinateUnit_2%Label), &
                len_trim(QM%CoordinateUnit_2%Label))
    if(len_trim(QM%CoordinateUnit_3%Label) > 0) &   
      Error = DBADDCOPT( &
                SiloOptionList, DBOPT_ZUNITS, &
                trim(QM%CoordinateUnit_3%Label), &
                len_trim(QM%CoordinateUnit_3%Label))
    if(present(CycleNumberOption)) &
      Error = DBADDIOPT(SiloOptionList, DBOPT_CYCLE, CycleNumberOption)
    if(present(TimeOption))then
      Time = TimeOption
      if(present(TimeUnitOption)) Time = Time / TimeUnitOption
      Error = DBADDIOPT(SiloOptionList, DBOPT_DTIME, Time)
    end if

    Error = DBPUTQM( &
              IOS%MeshBlockHandle, 'Mesh', 4, &
              QM%CoordinateLabel_1, len_trim(QM%CoordinateLabel_1), &
              QM%CoordinateLabel_2, len_trim(QM%CoordinateLabel_2), &
              QM%CoordinateLabel_3, len_trim(QM%CoordinateLabel_3), &
              QM%NodeCoordinate_1 / QM%CoordinateUnit_1%Value, &
              QM%NodeCoordinate_2 / QM%CoordinateUnit_2%Value, &
              QM%NodeCoordinate_3 / QM%CoordinateUnit_3%Value, &
              QM%nNodes, QM%nDimensions, DB_DOUBLE, DB_COLLINEAR, &
              SiloOptionList, Error)
    
    call WriteMeshMultiMesh(QM, IOS, SiloOptionList)
    
    Error = DBFREEOPTLIST(SiloOptionList)
    
  end subroutine WriteMesh
  
  
  subroutine WriteVariableArrayGroup_3D(QM, IOS, VAG)
    
    type(QuadMeshSiloForm), intent(inout) :: &
      QM
    type(InputOutputStreamForm), intent(inout) :: &
      IOS
    type(VariableArrayReal_3D_GroupBase), intent(in) :: &
      VAG
      
    integer(KDI) :: &
      iVariable, &
      iSelected, &
      SiloOptionList, &
      Centering, &
      Error
    character(LDF) :: &
      Mesh
      
    call MakeDirectory(IOS, VAG%Name)
    
    Mesh = ''
    if(QM%lDirectory > 0) Mesh = trim(Mesh) // trim(QM%Directory) // '/'
    if(QM%lName > 0) Mesh = trim(Mesh) // trim(QM%Name) // '/'
    Mesh = trim(Mesh) // 'Mesh'
    
    do iSelected = 1, size(VAG%Selected)
      
      iVariable = VAG%Selected(iSelected)
      
      if(len_trim(VAG%Unit(iVariable)%Label) > 0)then
        Error = DBMKOPTLIST(1, SiloOptionList)
        Error = DBADDCOPT( &
                  SiloOptionList, DBOPT_UNITS, &
                  trim(VAG%Unit(iVariable)%Label), &
                  len_trim(VAG%Unit(iVariable)%Label))
      else
        SiloOptionList = DB_F77NULL
      end if
      
      if(all(QM%nNodes == VAG%nData(iVariable)%Data))then
        Centering = DB_NODECENT
      else
        Centering = DB_ZONECENT
      end if 
      
      Error = DBPUTQV1( &
                IOS%MeshBlockHandle, trim(VAG%Variable(iVariable)), &
                VAG%lVariable(iVariable), Mesh, len_trim(Mesh), &
                VAG%Data(iVariable)%Data/VAG%Unit(iVariable)%Value, &
                VAG%nData(iVariable)%Data, QM%nDimensions, DB_F77NULL, 0, &
                DB_DOUBLE, Centering, SiloOptionList, Error)
      
      if(SiloOptionList /= DB_F77NULL) Error = DBFREEOPTLIST(SiloOptionList)
    
    end do
    
    call WriteVariableArrayGroupMultiMesh(QM, VAG, IOS)
    
    call ChangeDirectory(IOS, '../')
    
  end subroutine WriteVariableArrayGroup_3D
  
  
  subroutine ReadMesh(QM, IOS)
  
    type(QuadMeshSiloForm), intent(inout) :: &
      QM
    type(InputOutputStreamForm), intent(inout) :: &
      IOS
      
    call SiloWrappers_DBGetQuadMeshMtdt(IOS%MeshBlockHandle, QM%nNodes)
    
    allocate(QM%NodeCoordinate_1(QM%nNodes(1)))
    allocate(QM%NodeCoordinate_2(QM%nNodes(2)))
    allocate(QM%NodeCoordinate_3(QM%nNodes(3)))
    
    call SiloWrappers_DBGetQuadMeshCrdnt( &
           IOS%MeshBlockHandle, QM%NodeCoordinate_1, QM%NodeCoordinate_2, &
           QM%NodeCoordinate_3)
    
    QM%NodeCoordinate_1 = QM%NodeCoordinate_1 * QM%CoordinateUnit_1%Value
    QM%NodeCoordinate_2 = QM%NodeCoordinate_2 * QM%CoordinateUnit_2%Value
    QM%NodeCoordinate_3 = QM%NodeCoordinate_3 * QM%CoordinateUnit_3%Value
    
  end subroutine ReadMesh
  
  
  subroutine ReadVariableArrayGroup_3D(QM, IOS, VAG, AllocateStorageOption)
    
    type(QuadMeshSiloForm), intent(inout) :: &
      QM
    type(InputOutputStreamForm), intent(inout) :: &
      IOS
    type(VariableArrayReal_3D_GroupBase), intent(inout) :: &
      VAG
    logical(KDL), intent(in), optional :: &
      AllocateStorageOption
    
    integer(KDI) :: &
      iSelected, &
      iVariable
    logical(KDL) :: &
      AllocateStorage
    
    AllocateStorage = .false.
    if(present(AllocateStorageOption)) &
      AllocateStorage = AllocateStorageOption
    
    call ChangeDirectory(IOS, VAG%Name)
    
    do iSelected = 1, size(VAG%Selected)
      iVariable = VAG%Selected(iSelected)
      if(AllocateStorage)then
        call Initialize(VAG%Data(iVariable), QM%nNodes-1)
        VAG%nData(iVariable)%Data = shape(VAG%Data(iVariable)%Data)
      end if
      call SiloWrappers_DBGetQuadMeshVrbl( &
             IOS%MeshBlockHandle, trim(VAG%Variable(iVariable))//char(0), &
             VAG%Data(iVariable)%Data)
      VAG%Data(iVariable)%Data &
        = VAG%Data(iVariable)%Data * VAG%Unit(iVariable)%Value
    end do
    
    call ChangeDirectory(IOS, '../')
    
  end subroutine ReadVariableArrayGroup_3D
  
  
  subroutine WriteMeshMultiMesh(QM, IOS, SiloOptionList)

    type(QuadMeshSiloForm), intent(in) :: &
      QM
    type(InputOutputStreamForm), intent(in) :: &
      IOS
    integer(KDI), intent(in) :: &
      SiloOptionList

    integer(KDI) :: &
      iBlock, &
      nBlocks, &
      Error
    integer(KDI), dimension(:), allocatable :: &
      lBlockName, &
      SiloMeshType
    character(LDF), dimension(:), allocatable :: &
      BlockName
    
    if(IOS%Communicator%Rank /= 0) return
    if(QM%nDimensions == 1) return
    
    nBlocks = IOS%Communicator%Size
    allocate(BlockName(nBlocks))
    allocate(lBlockName(nBlocks))
    allocate(SiloMeshType(nBlocks))
      
    do iBlock = 1, nBlocks
      BlockName(iBlock) &
        = trim(IOS%MeshBlocksPrefix(iBlock)) &
          // trim(IOS%CurrentDirectory) // 'Mesh'
      lBlockName(iBlock)  = len_trim(BlockName(iBlock))
      SiloMeshType(iBlock) = DB_QUAD_RECT
    end do
      
    Error = DBSET2DSTRLEN(LF)
    Error = DBPUTMMESH( &
              IOS%MultiMeshHandle, 'Mesh', 4, nBlocks, &
              BlockName, lBlockName, SiloMeshType, SiloOptionList, Error)

    deallocate(BlockName)
    deallocate(lBlockName)
    deallocate(SiloMeshType)

  end subroutine WriteMeshMultiMesh


  subroutine WriteVariableArrayGroupMultiMesh(QM, VAG, IOS)

    type(QuadMeshSiloForm), intent(in) :: &
      QM
    type(VariableArrayReal_3D_GroupBase), intent(in) :: &
      VAG
    type(InputOutputStreamForm), intent(in) :: &
      IOS

    integer(KDI) :: &
      iVariable, &
      iBlock, &
      nBlocks, &
      Error
    integer(KDI), dimension(:), allocatable :: &
      lBlockVariable, &
      SiloVariableType
    character(LDF), dimension(:), allocatable :: &
      BlockVariable
    
    if(QM%nDimensions == 1) return
    if(IOS%Communicator%Rank /= 0) return
    
    nBlocks = IOS%Communicator%Size
    allocate(BlockVariable(nBlocks))
    allocate(lBlockVariable(nBlocks))
    allocate(SiloVariableType(nBlocks))
      
    do iVariable = 1, size(VAG%Selected)

      do iBlock = 1, nBlocks
        BlockVariable(iBlock) & 
          = trim(IOS%MeshBlocksPrefix(iBlock)) &
            // trim(IOS%CurrentDirectory) &
            // trim(VAG%Variable(VAG%Selected(iVariable)))
        lBlockVariable(iBlock) = len_trim(BlockVariable(iBlock))
        SiloVariableType(iBlock) = DB_QUADVAR
      end do
      
      Error = DBSET2DSTRLEN(LF)
      Error = DBPUTMVAR( &
                IOS%MultiMeshHandle, &
                trim(VAG%Variable(VAG%Selected(iVariable))), &
                VAG%lVariable(VAG%Selected(iVariable)), &
                nBlocks, BlockVariable, lBlockVariable, SiloVariableType, &
                DB_F77NULL, Error)
    
    end do
      
    !-- FIXME: vector is not supported in VariableArrayGroup yet
    
    deallocate(BlockVariable)
    deallocate(lBlockVariable)
    deallocate(SiloVariableType)

  end subroutine WriteVariableArrayGroupMultiMesh
                 
               
end module QuadMeshSilo_Form
