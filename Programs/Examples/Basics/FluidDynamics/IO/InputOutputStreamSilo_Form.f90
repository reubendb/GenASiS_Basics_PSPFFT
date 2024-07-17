module InputOutputStreamSilo_Form

  use pspfft

  implicit none
  private 
  
  include 'silo_f9x.inc'
  
  public :: &
    Create, &
    Open, &
    MakeDirectory, &
    ChangeDirectory, &
    Close, &
    Destroy
  
  type, public :: InputOutputStreamSiloForm
    integer(KI) :: &
      WRITE_ACCESS, &
      APPEND_ACCESS, &
      READ_ACCESS
    integer(KI) :: &
      nBlocks, &
      lName, &
      lDescription, &
      Number, &
      AccessMode, &
      VariableType, &
      MeshBlockHandle = -1, &
      MultiMeshHandle = -1
    character(LL), dimension(:), allocatable :: &
      DirectoryList, &
      VariableList
    character(LF) :: &
      Name, &
      Description, &
      CurrentDirectory
    character(LF), dimension(:), allocatable :: &
      MeshBlocksPrefix
    logical(KL) :: &
      Parallel
    type(CommunicatorForm), pointer :: &
      Communicator => null()
  end type InputOutputStreamSiloForm
  
  interface Create
    module procedure Create_IOSS
  end interface Create
  
  interface Open
    module procedure Open_IOSS
  end interface Open
  
  interface Close
    module procedure Close_IOSS
  end interface Close
  
  interface Destroy
    module procedure Destroy_IOSS
  end interface Destroy
  
    private :: &
      Open_IOSS_ForWriting, &
      Append_IOSS, &
      Open_IOSS_ForReading
    
      private :: &
        TrimDoubleSlash
    
  integer(KI), private, parameter :: &
    ACCESS_MODE_UNOPENED = -1, &
    ACCESS_MODE_WRITE  = 1, &
    ACCESS_MODE_APPEND = 2, &
    ACCESS_MODE_READ   = 3
  integer(KI), private, parameter :: &
    FIELD_TYPE_CURVE  = 1, &
    FIELD_TYPE_UCDVAR = 2
  
contains

  
  subroutine Create_IOSS(IOSS, Name, DescriptionOption, MPI_COMM_Option)
    
    type(InputOutputStreamSiloForm), pointer :: &
      IOSS
    character(*), intent(in) :: &
      Name
    character(*), intent(in), optional :: &
      DescriptionOption
    integer(KI), intent(in), optional :: &
      MPI_COMM_Option
    
    allocate(IOSS)
    
    IOSS%WRITE_ACCESS  = ACCESS_MODE_WRITE
    IOSS%APPEND_ACCESS = ACCESS_MODE_APPEND
    IOSS%READ_ACCESS   = ACCESS_MODE_READ

    IOSS%lName = len_trim(Name)

    IOSS%lDescription = IOSS%lName
    if(present(DescriptionOption)) &
      IOSS%lDescription = len_trim(DescriptionOption)

    IOSS%Number       = -1
    IOSS%AccessMode   = ACCESS_MODE_UNOPENED
    IOSS%VariableType = -1

    IOSS%Name = trim(Name)
    
    IOSS%Description = IOSS%Name
    if(present(DescriptionOption)) IOSS%Description = DescriptionOption
    
    if(present(MPI_COMM_Option))then
      IOSS%Parallel = .true.
      call Create(IOSS%Communicator, MPI_COMM_Option, 'IOSS_Communicator')
    else
      IOSS%Parallel = .false.
      IOSS%Communicator => null()
    end if
    
  end subroutine Create_IOSS
  
  
  subroutine Open_IOSS( &
               IOSS, AccessMode, WorkingDirectoryOption, SeriesOption, &
               NumberOption, BlockNumberOption)
               
    type(InputOutputStreamSiloForm), intent(inout) :: &
      IOSS
    integer(KDI), intent(in) :: &
      AccessMode
    character(*), intent(in), optional :: &
      WorkingDirectoryOption
    logical(KDL), intent(in), optional :: &
      SeriesOption
    integer(KDI), intent(in), optional :: &
      NumberOption, &
      BlockNumberOption
      
    if(AccessMode == ACCESS_MODE_WRITE) &
      call Open_IOSS_ForWriting( &
             IOSS, WorkingDirectoryOption,  SeriesOption, NumberOption)
      
    if(AccessMode == ACCESS_MODE_APPEND) &
      call Append_IOSS(IOSS, SeriesOption)
    
    if(AccessMode == ACCESS_MODE_READ) &
      call Open_IOSS_ForReading( &
             IOSS, WorkingDirectoryOption, NumberOption, BlockNumberOption)
    
  end subroutine Open_IOSS
  
  
  subroutine MakeDirectory(IOSS, Name)
    
    type(InputOutputStreamSiloForm), intent(inout) :: &
      IOSS
    character(*), intent(in) :: &
      Name
    
    integer(KDI) :: &
      Error
      
    if(IOSS%MeshBlockHandle == -1) return;
    
    if(IOSS%AccessMode /= ACCESS_MODE_WRITE &
       .and. IOSS%AccessMode /= ACCESS_MODE_APPEND) &
      return
    
    if(len_trim(Name) == 0) return
    if(trim(Name) == '/') return
    if(trim(Name) == '.' .or. trim(Name) == './') return
    if(trim(Name) == '..' .or. trim(Name) == '../') return
    
    Error = DBMKDIR(IOSS%MeshBlockHandle, trim(Name), len_trim(Name), Error)
    Error = DBSETDIR(IOSS%MeshBlockHandle, trim(Name), len_trim(Name))
    
    IOSS%CurrentDirectory = trim(IOSS%CurrentDirectory) // trim(Name) // '/'
    call TrimDoubleSlash(IOSS%CurrentDirectory)
    
    if(.not.IOSS%Parallel) return
    
    if(IOSS%Communicator%Rank == 0)then
      if(IOSS%MultiMeshHandle == -1) return;
      Error = DBMKDIR(IOSS%MultiMeshHandle, trim(Name), len_trim(Name), Error)
      Error = DBSETDIR(IOSS%MultiMeshHandle, trim(Name), len_trim(Name))
    end if
    
  end subroutine MakeDirectory
  
  
  subroutine ChangeDirectory(IOSS, PathName, SuccessOption)
    
    type(InputOutputStreamSiloForm), intent(inout) :: &
      IOSS
    character(*), intent(in) :: &
      PathName
    logical(KDL), intent(out), optional :: &
      SuccessOption
    
    integer(KDI) :: &
      Error, &
      LastSlash
    character(LDF) :: &
      Buffer
    
    if(len_trim(PathName) == 0) return
    
    !-- Early return means change directory fails
    if(present(SuccessOption)) SuccessOption = .false.
    
    if(IOSS%MeshBlockHandle == -1) return
    
    if(IOSS%AccessMode == ACCESS_MODE_UNOPENED) return

    Error = DBSETDIR(IOSS%MeshBlockHandle, trim(PathName), len_trim(PathName))
    if(Error /= 0) return
    
    !-- Set IOSS%CurrentDirectory to reflect the new directory we are in
    if(Pathname(:1) == '/')then
      IOSS%CurrentDirectory = trim(Pathname) // '/'
      call TrimDoubleSlash(IOSS%CurrentDirectory)
    else if(PathName(:2) == '..')then
      Buffer = IOSS%CurrentDirectory
      LastSlash = index(Buffer(:len_trim(Buffer)-1), '/', Back=.true.)
      IOSS%CurrentDirectory = IOSS%CurrentDirectory(:LastSlash)
      if(len_trim(PathName) > 3) &
        IOSS%CurrentDirectory &
          = trim(IOSS%CurrentDirectory) // trim(PathName(4:)) // '/'
      call TrimDoubleSlash(IOSS%CurrentDirectory)
    else
      IOSS%CurrentDirectory &
        = trim(IOSS%CurrentDirectory) // trim(PathName) // '/'
      call TrimDoubleSlash(IOSS%CurrentDirectory)
    end if
    
    if(IOSS%Parallel)then
      if(IOSS%Communicator%Rank == 0)then
        if(IOSS%MultiMeshHandle == -1) return
        Error = DBSETDIR( &
                  IOSS%MultiMeshHandle, trim(PathName), len_trim(PathName))
        if(Error /= 0) return
      end if
    end if
      
    if(present(SuccessOption)) SuccessOption = .true.
    
  end subroutine ChangeDirectory
  
  
  subroutine Close_IOSS(IOSS)
    
    type(InputOutputStreamSiloForm), intent(inout) :: &
      IOSS
    
    integer(KDI) :: &
      Error
    
    if(IOSS%MeshBlockHandle == -1) return;       
    
    Error = DBCLOSE(IOSS%MeshBlockHandle)
    if(Error == 0) IOSS%MeshBlockHandle = -1
      
    if(IOSS%MultiMeshHandle /= -1)then
      Error = DBCLOSE(IOSS%MultiMeshHandle)
      if(Error == 0) IOSS%MultiMeshHandle = -1
    end if
    
    IOSS%AccessMode = ACCESS_MODE_UNOPENED
    
    if(allocated(IOSS%MeshBlocksPrefix)) &
      deallocate(IOSS%MeshBlocksPrefix)
    
  end subroutine Close_IOSS
  
  
  subroutine Destroy_IOSS(IOSS)
  
    type(InputOutputStreamSiloForm), pointer :: &
      IOSS
    
    if(.not.associated(IOSS)) return
    
    if(allocated(IOSS%MeshBlocksPrefix)) deallocate(IOSS%MeshBlocksPrefix)
    if(allocated(IOSS%VariableList)) deallocate(IOSS%VariableList)
    if(allocated(IOSS%DirectoryList)) deallocate(IOSS%DirectoryList)
    
    if(associated(IOSS%Communicator)) &
      call Destroy(IOSS%Communicator, DestroyHandleOption = .false.)
    
    deallocate(IOSS)
      
  end subroutine Destroy_IOSS
  
  
  subroutine Open_IOSS_ForWriting( &
               IOSS, WorkingDirectoryOption, SeriesOption, NumberOption)
    
    type(InputOutputStreamSiloForm), intent(inout) :: &
      IOSS
    character(*), intent(in), optional :: &
      WorkingDirectoryOption
    logical(KDL), intent(in), optional :: &
      SeriesOption
    integer(KDI), intent(in), optional :: &
      NumberOption
      
    integer(KDI) :: &
      iBlock, &
      lPathName, &
      Error
    logical(KDL) :: &
      Series
    character(LDF) :: &
      WorkingDirectory, &
      MeshBlockDirectory, &
      PathSuffix, &
      PathName
    character(LN+1) :: &
      FileNumberString, &
      BlockNumberString
    
    IOSS%Number     = IOSS%Number + 1
    IOSS%AccessMode = ACCESS_MODE_WRITE
    
    WorkingDirectory = './'
    if(present(WorkingDirectoryOption)) &
      WorkingDirectory =  WorkingDirectoryOption
    
    Series = .true.
    if(present(SeriesOption)) Series = SeriesOption
    
    if(present(NumberOption)) IOSS%Number = NumberOption
    
    if(Series)then
      write(FileNumberString, fmt='(a1,i7.7)') '_', IOSS%Number
    else
      FileNumberString = ''
    end if
    
    if(IOSS%Parallel)then
      write(BlockNumberString, fmt='(a1,i7.7)') '_', IOSS%Communicator%Rank
    else
      BlockNumberString = ''
    end if
    
    if(IOSS%Parallel)then
    
      IOSS%nBlocks = IOSS%Communicator%Size
      
      MeshBlockDirectory &
        = trim(WorkingDirectory) //  trim(adjustl(IOSS%Name)) &
          // trim(FileNumberString) // '_MeshBlocks/'
      if(IOSS%Communicator%Rank == 0)then
        call SiloWrappers_MakeDirectory( &
               trim(MeshBlockDirectory)//char(0), Error)
      end if
      call Block(IOSS%Communicator)
      
      !-- Create file for a mesh block and thereby obtain 
      !   IOSS%MeshBlockHandle
      
      PathSuffix &
        = trim(FileNumberString) // '_MeshBlock' // trim(BlockNumberString) &
          // '.silo'
      PathName  &
        = trim(MeshBlockDirectory) // trim(adjustl(IOSS%Name)) &
          // trim(PathSuffix)
      lPathName = len_trim(PathName)
        
      Error = DBCREATE( &
                PathName, lPathName, DB_CLOBBER, DB_LOCAL, &
                IOSS%Description, IOSS%lDescription, DB_PDB, &
                IOSS%MeshBlockHandle)
      
      if(IOSS%Communicator%Rank == 0)then
        
        !-- Create multi-block mesh file and thereby obtain 
        !   IOSS%MultiMeshHandle

        PathSuffix = '_MultiMesh' // trim(FileNumberString) // '.silo'
        PathName  &
          = trim(WorkingDirectory) // trim(adjustl(IOSS%Name)) &
            // trim(PathSuffix)
        lPathName = len_trim(PathName)
        
        Error = DBCREATE( &
                  PathName, lPathName, DB_CLOBBER, DB_LOCAL, &
                  IOSS%Description, IOSS%lDescription, DB_PDB, &
                  IOSS%MultiMeshHandle)
                  
      end if
      
      if(IOSS%Communicator%Rank == 0)then
        allocate(IOSS%MeshBlocksPrefix(IOSS%nBlocks))
        do iBlock = 1, IOSS%nBlocks
          write(BlockNumberString, fmt='(a1,i7.7)') '_', iBlock - 1
          PathSuffix &
            = trim(FileNumberString) // '_MeshBlock' &
              // trim(BlockNumberString) // '.silo'
          IOSS%MeshBlocksPrefix(iBlock) &
            = trim(adjustl(IOSS%Name)) // trim(FileNumberString) &
              // '_MeshBlocks/' // trim(adjustl(IOSS%Name)) &
              // trim(PathSuffix) // ':'
        end do
      end if
    
    else
    
      IOSS%nBlocks = 0
      PathSuffix = trim(FileNumberString) // '.silo'
      PathName  = trim(adjustl(IOSS%Name)) // trim(PathSuffix)
      lPathName = len_trim(PathName)
      
      !-- Create file for a mesh block and thereby obtain IOSS%MeshBlockHandle
      
      Error = DBCREATE( &
                PathName, lPathName, DB_CLOBBER, DB_LOCAL, &
                IOSS%Description, IOSS%lDescription, DB_PDB, &
                IOSS%MeshBlockHandle)
    
    end if
    
    IOSS%CurrentDirectory = '/'
    
  end subroutine Open_IOSS_ForWriting
  
  
  subroutine Append_IOSS(IOSS, SeriesOption)
  
    type(InputOutputStreamSiloForm), intent(inout) :: &
      IOSS
    logical(KDL), intent(in), optional :: &
      SeriesOption
      
    integer(KI) :: &
      iBlock, &
      lPathName, &
      Error
    logical(KDL) :: &
      Series
    character(LN+1) :: &
      FileNumberString, &
      BlockNumberString
    character(LDF):: &
      PathName, &
      PathSuffix
    
    if(IOSS%Parallel)then
      IOSS%nBlocks = IOSS%Communicator%Size
    else
      IOSS%nBlocks = 0
    end if
    
    IOSS%AccessMode = ACCESS_MODE_APPEND
    
    Series = .true.
    if(present(SeriesOption)) Series = SeriesOption
    
    if(Series)then
      write(FileNumberString, fmt='(a1,i7.7)') '_', IOSS%Number
    else
      FileNumberString = ''
    end if
    
    if(IOSS%Parallel)then
      write(BlockNumberString, fmt='(a1,i7.7)') '_', IOSS%Communicator%Rank
    else
      BlockNumberString = ''
    end if
    
    !-- Open file for a mesh block and thereby obtain IOSS%MeshBlockHandle
    
    if(IOSS%Parallel)then
      PathSuffix &
        = trim(FileNumberString) // '_MeshBlock' // trim(BlockNumberString) &
          // '.silo'
      PathName = trim(adjustl(IOSS%Name)) // trim(PathSuffix)
      lPathName = len_trim(PathName)
      
      Error = DBOPEN( &
                PathName, lPathName, DB_PDB, DB_APPEND, IOSS%MeshBlockHandle)
      
      if(IOSS%Communicator%Rank == 0)then

        !-- Open multi-block mesh file and thereby obtain &
        !   IOSS%MultiMeshHandle

        PathSuffix = '_MultiMesh' // trim(FileNumberString) // '.silo'
        PathName = trim(adjustl(IOSS%Name)) // trim(PathSuffix)
        lPathName = len_trim(PathName)
        
        Error = DBOPEN( &
                  PathName, lPathName, DB_PDB, DB_APPEND, &
                  IOSS%MultiMeshHandle)
        
      end if

      if(IOSS%Communicator%Rank == 0)then
        allocate(IOSS%MeshBlocksPrefix(IOSS%nBlocks))    
        do iBlock = 1, IOSS%nBlocks
          write(BlockNumberString, fmt='(a1,i7.7)') '_', iBlock - 1
          PathSuffix &
            = trim(FileNumberString) // '_MeshBlock' &
              // trim(BlockNumberString) // '.silo'
          IOSS%MeshBlocksPrefix(iBlock) &
            = trim(adjustl(IOSS%Name)) // trim(FileNumberString) &
              // '_MeshBlocks/' // trim(adjustl(IOSS%Name)) &
              // trim(PathSuffix) // ':'
        end do
      end if
      
    else
    
      IOSS%nBlocks = 0
      PathSuffix = trim(FileNumberString) // '.silo'
      PathName  = trim(adjustl(IOSS%Name)) // trim(PathSuffix)
      lPathName = len_trim(PathName)
      
      !-- Create file for a mesh block and thereby obtain IOSS%MeshBlockHandle
      
      Error = DBOPEN( &
                PathName, lPathName, DB_PDB, DB_APPEND, IOSS%MeshBlockHandle)
      
    end if
    
    IOSS%CurrentDirectory = '/'
      
  end subroutine Append_IOSS
  
  
  subroutine Open_IOSS_ForReading( &
               IOSS, WorkingDirectoryOption, NumberOption, BlockNumberOption)
    
    type(InputOutputStreamSiloForm), intent(inout) :: &
      IOSS
    character(*), intent(in), optional :: &
      WorkingDirectoryOption
    integer(KDI), intent(in), optional :: &
      NumberOption, &
      BlockNumberOption
         
    integer(KDI) :: &
      lPathName, &
      lWorkingDirectory, &
      BlockNumber, &
      Error
    logical(KDL) :: &
      MeshBlocksExist, &
      MultiMeshExists
    character(LDF) :: &
      MeshBlocksPrefix, &
      MultiMeshPrefix, &
      PathName, &
      WorkingDirectory, &
      PathSuffix
    character(LN+1) :: &
      FileNumberString, &
      BlockNumberString
    
      
    if(IOSS%AccessMode == ACCESS_MODE_WRITE &
       .or. IOSS%AccessMode == ACCESS_MODE_APPEND) &
      return
      
    if(IOSS%AccessMode == ACCESS_MODE_READ) call Close(IOSS)
    
    WorkingDirectory = ''
    if(present(WorkingDirectoryOption)) &
      WorkingDirectory = WorkingDirectoryOption 
    lWorkingDirectory = len_trim(WorkingDirectory)
    
    FileNumberString = ''
    if(present(NumberOption))then
      IOSS%Number = NumberOption
      write(FileNumberString, fmt='(a1,i7.7)') '_', IOSS%Number
    end if
    
    if(IOSS%Parallel)then
      BlockNumber = IOSS%Communicator%Rank
      if(present(BlockNumberOption)) BlockNumber = BlockNumberOption
      write(BlockNumberString, fmt='(a1,i7.7)') '_', BlockNumber
      !-- Setting of IOSS%nBlocks for Parallel case is deferred until 
      !   MultiMesh handle is acquired
    else
      BlockNumberString = ''
      IOSS%nBlocks = 0
    end if
        
    IOSS%AccessMode = ACCESS_MODE_READ

    if(IOSS%Parallel)then
      
      !-- Try to open mesh block files
      MultiMeshPrefix = ''
      PathSuffix &
       = trim(FileNumberString) // '_MeshBlock' // trim(BlockNumberString) &
         // '.silo'
      MeshBlocksPrefix &
        = trim(adjustl(IOSS%Name)) // trim(FileNumberString) // '_MeshBlocks/'
      PathName  &
        = trim(WorkingDirectory) // trim(MeshBlocksPrefix) &
          // trim(adjustl(IOSS%Name)) // trim(PathSuffix)
      lPathName = len_trim(PathName)
      
      inquire(file = trim(PathName), exist = MeshBlocksExist)
      
      if(.not.MeshBlocksExist)then
        call Show('Trying to open non-existing input file', CONSOLE_ERROR)
        call Show(trim(PathName), 'File', CONSOLE_ERROR)
        call Abort(IOSS%Communicator)
      end if
    
    else
      
      !-- Try to open mesh block files
      
      MeshBlocksPrefix = ''
      MultiMeshPrefix = ''
      PathSuffix = trim(FileNumberString) // '.silo'
      PathName  &
        = trim(WorkingDirectory) // trim(MeshBlocksPrefix) &
          // trim(adjustl(IOSS%Name)) // trim(PathSuffix)
      lPathName = len_trim(PathName)
      
      inquire(file = trim(PathName), exist = MeshBlocksExist)
      
      if(.not.MeshBlocksExist)then
        call Show('Trying to open non-existing input file', CONSOLE_ERROR)
        call Show(trim(PathName), 'File', CONSOLE_ERROR)
      end if
    
    end if
  
    Error = DBOPEN( &
              trim(PathName), lPathName, DB_PDB, DB_READ, &
              IOSS%MeshBlockHandle)
    
    !-- Open MultiMesh file
    
    if(IOSS%Parallel)then
      
      PathSuffix = '_MultiMesh' // trim(FileNumberString) // '.silo'
      PathName &
        = trim(WorkingDirectory) // trim(MultiMeshPrefix) &
          // trim(adjustl(IOSS%Name)) // trim(PathSuffix)
      lPathName = len_trim(PathName)
      
      inquire(file = trim(PathName), exist = MultiMeshExists)
      if(.not.MultiMeshExists)then
        call Show('Trying to open non-existing input file', CONSOLE_ERROR)
        call Show(trim(PathName), 'File', CONSOLE_ERROR)
        call Abort(IOSS%Communicator)
      end if
      
      Error = DBOPEN( &
                trim(PathName), lPathName, DB_PDB, DB_READ, &
                IOSS%MultiMeshHandle)
      
      !-- Assume number of blocks is consistent with the job size
      IOSS%nBlocks = IOSS%Communicator%Size
      
    end if
      
    IOSS%CurrentDirectory = '/'
    
  end subroutine Open_IOSS_ForReading
  
  
  subroutine TrimDoubleSlash(Path)
  
    character(*), intent(inout) :: &
      Path
      
    integer(KDI) :: &
      TrimmedLength
      
    TrimmedLength = len_trim(Path)
    
    if(TrimmedLength <= 1) return
    
    if(Path(TrimmedLength-1:TrimmedLength) == '//')then
      Path = Path(:TrimmedLength-1)
    end if
  
  end subroutine TrimDoubleSlash
  
    
end module InputOutputStreamSilo_Form
