!-- This module defines a class to provide an abstraction of the isolated
!   laplacian operator with the Green's function necessary for to solve the 
!   the resulting Poisson equation. The object also provides the data storage 
!   of the source (and solution) in form of pillars.

module LaplacianIsolated_FFT__Form
  
  use Basics
  use Transpose_Command
  use FFT_FFTW__Base, &
        FFT_Base => FFT_FFTW_Base

  implicit none
  private
  
  public :: &
    Create, &
    Load, &
    Store, &
    Destroy
  
  type, public :: LaplacianIsolated_FFT_Form
    integer(KDI), dimension(3) :: &
      nCellsBrick, &
      PillarWidth, PillarHeight
    complex(KDC), dimension(:,:,:), allocatable :: &
      GreensFunction_Z
    type(Real_3D_Form), dimension(:), pointer :: &
      InputOutput
    type(CommunicatorForm), pointer :: &
      Communicator_X, Communicator_XY, Communicator_YZ
    type(FFT_Base), dimension(3) :: &
      FFT_Forward, FFT_Backward
  end type LaplacianIsolated_FFT_Form
  
  interface Create
    module procedure Create_L
  end interface Create
  
  interface Load
    module procedure Load_L
  end interface Load
  
  interface Store
    module procedure Store_L
  end interface Store
  
  interface Destroy
    module procedure Destroy_L
  end interface Destroy
  
    private :: &
      CreateCommunicators, &
      ComputeGreensFunction

      private :: &
        ComputeGreensFunction_X, &
        ComputeGreensFunction_Y, &
        ComputeGreensFunction_Z
  
contains

  
  subroutine Create_L(L, C, CellWidth, nCells)
  
    !-- Setup and allocate the "pillars" data storage (see manuscript about 
    !   pillars ).
  
    type(LaplacianIsolated_FFT_Form), pointer :: &
      L
    type(CommunicatorForm), intent(in) :: &
      C
    real(KDR), dimension(3), intent(in) :: &
      CellWidth
    integer(KDI), dimension(3), intent(in) :: &
      nCells
      
    integer(KDI) :: &
      iDim, jDim, kDim, &
      nRanksRoot
    integer(KDI), dimension(3) :: &
      PW, PH, &
      Factor
    
    allocate(L)
    
    nRanksRoot = C%Size**(1.0_KDR/3) + 0.5_KDR
    L%nCellsBrick = nCells / nRanksRoot
    
    if(any(mod(L%nCellsBrick, nRanksRoot) /= 0))then
      call Show( &
             'FFT requires that brick widths be divisible by nRanksRoot', &
             IgnorabilityOption = CONSOLE % ERROR)
      call Show( & 
             L%nCellsBrick(1), 'brick width 1', &
             IgnorabilityOption = CONSOLE % ERROR)
      call Show( & 
             L%nCellsBrick(2), 'brick width 2', &
             IgnorabilityOption = CONSOLE % ERROR)
      call Show( & 
             L%nCellsBrick(3), 'brick width 3', &
             IgnorabilityOption = CONSOLE % ERROR)
      call Show( &
             nRanksRoot, 'nRanksRoot', &
             IgnorabilityOption = CONSOLE % ERROR)
      call Abort()
    end if
    
    do iDim = 1, 3
      jDim = mod(iDim, 3) + 1
      kDim = mod(jDim, 3) + 1
      if(mod(iDim, 2) == 1)then
        L%PillarWidth(iDim) = L%nCellsBrick(jDim) / nRanksRoot
        L%PillarHeight(iDim) = L%nCellsBrick(kDim)
      else
        L%PillarWidth(iDim) = L%nCellsBrick(jDim) 
        L%PillarHeight(iDim) = L%nCellsBrick(kDim) / nRanksRoot
      end if
    end do
    PW = L%PillarWidth
    PH = L%PillarHeight
    
    allocate(L%GreensFunction_Z(2*nCells(3),2*PW(3),2*PH(3)))

    call CreateCommunicators(L, C, nCells, nRanksRoot)

    Factor = 1
    do iDim = 1, 3
      Factor(2) = 2
      Factor = cshift(Factor, 1)
      call Initialize( &
             L%FFT_Forward(iDim), L%FFT_Forward(iDim)%FORWARD, &
             DataShapeOption &
               = (/Factor(1)*nCells(iDim), Factor(2)*PW(iDim), &
                   Factor(3)*PH(iDim)/))
      call Initialize( &
             L%FFT_Backward(iDim), L%FFT_Backward(iDim)%BACKWARD, &
             Data_3D_Option = L%FFT_Forward(iDim)%Data_3D)
    end do
    
    call ComputeGreensFunction(L, CellWidth, nCells, C%Rank)

  end subroutine Create_L
  
  
  subroutine Load_L(L)

    !-- Distribute the source, in bricks, to x pillars

    type(LaplacianIsolated_FFT_Form), intent(inout) :: &
      L
    
    integer(KDI) :: &
      iRank, &
      iSource, &
      oData, &
      oBuffer, &
      nSources, &
      nData, &
      ChunkSize
    real(KDR), dimension(:), allocatable :: &
      SendBuffer, &
      ReceiveBuffer
    type(CommunicatorForm), pointer :: &
      C
    
    C => L%Communicator_X  
    
    nSources = size(L%InputOutput)

    if(nSources > 1) stop    
    
    nData = 0
    do iSource = 1, nSources
      nData = nData + size(L%InputOutput(iSource)%Data)
    end do
    allocate(SendBuffer(nData))
    allocate(ReceiveBuffer(nData))
    
    call Clear(L%FFT_Forward(1)%Data_3D)
    call Clear(L%FFT_Forward(2)%Data_3D)
    call Clear(L%FFT_Forward(3)%Data_3D)
    
    ChunkSize = product(L%nCellsBrick) * nSources / C%Size
    
    oBuffer = 0
    do iRank = 0, C%Size-1
      oData = iRank * L%nCellsBrick(2) / C%Size 
      do iSource = 1, nSources
        SendBuffer(oBuffer+1:oBuffer+ChunkSize/nSources) &
          = reshape( &
              L%InputOutput(iSource) &
                %Data(:,oData+1:oData+L%nCellsBrick(2)/C%Size,:), &
              (/ ChunkSize / nSources /))
        oBuffer = oBuffer + ChunkSize / nSources
      end do
    end do
    
    call AllToAll(SendBuffer, C, ChunkSize, ReceiveBuffer)
    
    oBuffer = 0
    do iRank = 0, C%Size-1
      oData = iRank * L%nCellsBrick(1) 
      do iSource = 1, nSources
        L%FFT_Forward(1)%Data_3D(oData+1:oData+L%nCellsBrick(1),:,:) &
          = reshape( &
              ReceiveBuffer(oBuffer+1:oBuffer+ChunkSize/nSources), &
              (/L%nCellsBrick(1), L%PillarWidth(1), L%PillarHeight(1)/))
        oBuffer = oBuffer + ChunkSize / nSources
      end do
    end do
    
    deallocate(SendBuffer)
    deallocate(ReceiveBuffer)
  
  end subroutine Load_L
  
  
  subroutine Store_L(L)
  
    !-- Distribute the solution, in x pillars, to bricks

    type(LaplacianIsolated_FFT_Form), intent(inout) :: &
      L
    
    integer(KDI) :: &
      iRank, &
      iSolution, &
      oData, &
      oBuffer, &
      nSolutions, &
      nData, &
      ChunkSize
    real(KDR), dimension(:), allocatable :: &
      SendBuffer, &
      ReceiveBuffer
    type(CommunicatorForm), pointer :: &
      C
    
    C => L%Communicator_X  
    
    nSolutions = size(L%InputOutput)

    if(nSolutions > 1) stop
    
    nData = 0
    do iSolution = 1, nSolutions
      nData = nData + size(L%InputOutput(iSolution)%Data)
    end do
    allocate(SendBuffer(nData))
    allocate(ReceiveBuffer(nData))
    
    ChunkSize = product(L%nCellsBrick) * nSolutions / C%Size
    
    oBuffer = 0 
    do iRank = 0, C%Size-1
      oData = iRank * L%nCellsBrick(1) 
      do iSolution = 1, nSolutions
        SendBuffer(oBuffer+1:oBuffer+ChunkSize/nSolutions) &
          = reshape( &
              real(L%FFT_Backward(1) &
                    %Data_3D(oData+1:oData+L%nCellsBrick(1),:,:), KDR), &
              (/ChunkSize / nSolutions/))
        oBuffer = oBuffer + ChunkSize / nSolutions
      end do
    end do
    
    call AllToAll(SendBuffer, C, ChunkSize, ReceiveBuffer)  
    
    oBuffer = 0
    do iRank = 0, C%Size-1
      oData = iRank * L%nCellsBrick(2) / C%Size 
      do iSolution = 1, nSolutions
        L%InputOutput(iSolution) &
          %Data(:,oData+1:oData+L%nCellsBrick(2)/C%Size,:) &
            = reshape( &
                ReceiveBuffer(oBuffer+1:oBuffer+ChunkSize/nSolutions), &
                (/L%nCellsBrick(1), L%PillarWidth(1), L%nCellsBrick(3)/))
        oBuffer = oBuffer + ChunkSize / nSolutions
      end do
    end do
    
    deallocate(SendBuffer)
    deallocate(ReceiveBuffer)
  
  end subroutine Store_L
  
  
  subroutine Destroy_L(L)
  
    type(LaplacianIsolated_FFT_Form), pointer :: &
      L
    
    integer(KDI) :: &
      iDim
    
    do iDim = 1, 3
      call Finalize(L%FFT_Backward(iDim))
      call Finalize(L%FFT_Forward(iDim))
    end do
    
    call Destroy(L%Communicator_YZ)
    call Destroy(L%Communicator_XY)
    call Destroy(L%Communicator_X)
    
    deallocate(L%GreensFunction_Z)

    deallocate(L)
    
  end subroutine Destroy_L
  

  subroutine CreateCommunicators(L, C, nCells, nRanksRoot)
  
    !-- Create sub-communicators for processors in pillars and slabs necessary
    !   for the brick to pillar distribution (and vice-versa) and pillar 
    !   transpose in a slab.

    type(LaplacianIsolated_FFT_Form), intent(inout) :: &
      L
    type(CommunicatorForm), intent(in) :: &
      C
    integer(KDI), dimension(3), intent(in) :: &
      nCells
    integer(KDI), intent(in) :: &
      nRanksRoot

    integer(KDI) :: &
      iRank, &
      nRanksPerSlab_XY, nRanksPerSlab_YZ, &
      BrickPosition_Y, BrickPosition_Z, &
      Pillar_X, &
      Slab_XY, Slab_YZ
    integer(KDI), dimension(3) :: &
      PW, PH
    integer(KDI), dimension(:), allocatable :: &
      Pillar_X_Rank, &
      Slab_XY_Rank, Slab_YZ_Rank

    PW = L%PillarWidth
    PH = L%PillarHeight
    
    BrickPosition_Y = mod(C%Rank/nRanksRoot, nRanksRoot)
    BrickPosition_Z = C%Rank / (nRanksRoot**2)
    allocate(Pillar_X_Rank(nRanksRoot))
    Pillar_X = BrickPosition_Y + BrickPosition_Z*nRanksRoot + 1
    Pillar_X_Rank &
      = (/(iRank, iRank = (Pillar_X-1)*nRanksRoot, Pillar_X*nRanksRoot - 1)/)
    call Create( &
           L%Communicator_X, NameOption='Pillar_X', C_Option=C, &
           RanksOption=Pillar_X_Rank)
    deallocate(Pillar_X_Rank)
    
    nRanksPerSlab_XY = nCells(2) / PW(1)
    nRanksPerSlab_YZ = nCells(3) / PH(1)
    allocate(Slab_XY_Rank(nRanksPerSlab_XY))
    allocate(Slab_YZ_Rank(nRanksPerSlab_YZ))
    Slab_XY = C%Rank / nRanksPerSlab_XY + 1
    Slab_YZ = mod(C%Rank, nRanksPerSlab_XY) + 1
    Slab_XY_Rank &
      = (/(iRank, & 
           iRank = (Slab_XY-1)*nRanksPerSlab_XY,Slab_XY*nRanksPerSlab_XY - 1)/)
    Slab_YZ_Rank = (/(iRank, iRank = Slab_YZ-1, C%Size-1, nRanksPerSlab_XY)/)
    call Create( &
           L%Communicator_XY, NameOption='Slab_XY', C_Option=C, &
           RanksOption=Slab_XY_Rank)
    call Create( &
           L%Communicator_YZ, NameOption='Slab_YZ', C_Option=C, &
           RanksOption=Slab_YZ_Rank)
    deallocate(Slab_YZ_Rank)
    deallocate(Slab_XY_Rank)
    
  end subroutine CreateCommunicators


  subroutine ComputeGreensFunction(L, CellWidth, nCells, ProcessRank)
    
    type(LaplacianIsolated_FFT_Form), intent(inout) :: &
      L
    real(KDR), dimension(3), intent(in) :: &
      CellWidth
    integer(KDI), dimension(3), intent(in) :: &
      nCells
    integer(KDI), intent(in) :: &
      ProcessRank
      
    integer(KDI) :: &
      nRanksPerSlab_XY, &
      RankPosition_Y, RankPosition_Z
    integer(KDI), dimension(3) :: &
      PW, PH
    type(FFT_Base) :: &
      FFT_XIA_X, FFT_XAI_X, FFT_XII_X, FFT_XYI_Y
      
    !-- The triplet in the comments and/or variables below indicate whether 
    !   the data pillar is an Active (A) or Inactive (I) in each of the 
    !   three-dimensions. The last letter indicates the direction of the 
    !   pillar (i.e. X, Y, or Z).
    
    PW = L%PillarWidth
    PH = L%PillarHeight
    
    call Initialize( &
           FFT_XIA_X, FFT_XIA_X%FORWARD, &
           DataShapeOption = (/2*nCells(1),1,PH(1)/) )
    call Initialize( &
           FFT_XAI_X, FFT_XAI_X%FORWARD, &
           DataShapeOption = (/2*nCells(1),PW(1),1/) )
    call Initialize( &
           FFT_XII_X, FFT_XII_X%FORWARD, nDataOption = 2*nCells(1))
    call Initialize( &
           FFT_XYI_Y, FFT_XYI_Y%FORWARD, &
           DataShapeOption = (/2*nCells(2),1,2*PH(2)/) )
    
    nRanksPerSlab_XY = nCells(2) / PW(1)
    RankPosition_Y = mod(ProcessRank, nRanksPerSlab_XY)
    RankPosition_Z = ProcessRank / nRanksPerSlab_XY

    call ComputeGreensFunction_X( &
           L, FFT_XIA_X, FFT_XAI_X, FFT_XII_X, CellWidth, nCells, &
           RankPosition_Y, RankPosition_Z)
    call ComputeGreensFunction_Y( &
           L, FFT_XIA_X, FFT_XAI_X, FFT_XII_X, FFT_XYI_Y, nCells, &
           RankPosition_Y)
    call ComputeGreensFunction_Z(L, FFT_XYI_Y, nCells, RankPosition_Z)

    call Copy(L%FFT_Forward(3)%Data_3D, L%GreensFunction_Z)
    
    call Finalize(FFT_XYI_Y)
    call Finalize(FFT_XII_X)
    call Finalize(FFT_XAI_X)
    call Finalize(FFT_XIA_X)
    
  end subroutine ComputeGreensFunction
  

  subroutine ComputeGreensFunction_X( &
               L, FFT_XIA_X, FFT_XAI_X, FFT_XII_X, H, nCells, &
               RankPosition_Y, RankPosition_Z)

    type(LaplacianIsolated_FFT_Form), intent(inout) :: &
      L
    type(FFT_Base), intent(inout) :: &
      FFT_XIA_X, FFT_XAI_X, FFT_XII_X
    real(KDR), dimension(3), intent(in) :: &
      H
    integer(KDI), dimension(3), intent(in) :: &
      nCells
    integer(KDI), intent(in) :: &
      RankPosition_Y, RankPosition_Z

    integer(KDI) :: &
      iCell, jCell, kCell, &
      P, Q, R
    integer(KDI), dimension(3) :: &
      PW, PH
    real(KDR) :: &
      Pi

    !-- The triplet in the comments and/or variables below indicate whether 
    !   the data pillar is an Active (A) or Inactive (I) in each of the 
    !   three-dimensions. The last letter indicates the direction of the 
    !   pillar (i.e. X, Y, or Z).

    PW = L%PillarWidth
    PH = L%PillarHeight
    
    Pi = acos(-1.0_KDR)
    
    !-- fill in AAA_X
    do kCell = 1, PH(1)
      do jCell = 1, PW(1)
        do iCell = 1, nCells(1)

          P = iCell - 1
          Q = (jCell - 1) + RankPosition_Y * PW(1)
          R = (kCell - 1) + RankPosition_Z * PH(1)
          
          if(P + Q + R == 0)then
            L%FFT_Forward(1)%Data_3D(iCell,jCell,kCell) &
              = - product(H) / (4.0_KDR * Pi * minval(H))
          else
            L%FFT_Forward(1)%Data_3D(iCell,jCell,kCell) &
              = - product(H) &
                  / (4.0_KDR * Pi &
                     * sqrt(H(1)**2*P**2 + H(2)**2*Q**2 + H(3)**2*R**2))
          end if

        end do
      end do
    end do
    
    !-- fill in IAA_X
    do kCell = 1, PH(1)
      do jCell = 1, PW(1)

        P = nCells(1)
        Q = (jCell - 1) + RankPosition_Y * PW(1)
        R = (kCell - 1) + RankPosition_Z * PH(1)
        
        L%FFT_Forward(1)%Data_3D(nCells(1)+1,jCell,kCell) &
          = - product(H) &
              / (4.0_KDR * Pi &
                 * sqrt(H(1)**2*P**2 + H(2)**2*Q**2 + H(3)**2*R**2))
          
        L%FFT_Forward(1)%Data_3D(nCells(1)+2:2*nCells(1), jCell, kCell) &
          = L%FFT_Forward(1)%Data_3D(nCells(1):2:-1, jCell, kCell)

      end do
    end do
    
    !-- fill in AIA_X and IIA_X
    do kCell = 1, PH(1)
      do iCell = 1, nCells(1)

        P = (iCell - 1)
        Q = nCells(2)
        R = (kCell - 1) + RankPosition_Z * PH(1)
        
        FFT_XIA_X%Data_3D(iCell, 1, kCell) &
          = - product(H) &
              / (4.0_KDR * Pi &
                 * sqrt(H(1)**2*P**2 + H(2)**2*Q**2 + H(3)**2*R**2))
          
        if(iCell == 1) then
          P = nCells(1)
          FFT_XIA_X%Data_3D(nCells(1)+1, 1, kCell) &
            = - product(H) &
                / (4.0_KDR * Pi &
                   * sqrt(H(1)**2*P**2 + H(2)**2*Q**2 + H(3)**2*R**2))
        end if

      end do

      FFT_XIA_X%Data_3D(nCells(1)+2:,1,kCell) &
        = FFT_XIA_X%Data_3D(nCells(1):2:-1,1,kCell)

    end do
    
    !-- fill in AAI_X, IAI_X, AII_X, III_X
    do jCell = 1, PW(1)
      do iCell = 1, nCells(1)

        P = (iCell - 1)
        Q = (jCell - 1) + RankPosition_Y * PW(1)
        R = nCells(3)
        
        FFT_XAI_X%Data_3D(iCell, jCell, 1) &
          = - product(H) &
              / (4.0_KDR * Pi &
                 * sqrt(H(1)**2*P**2 + H(2)**2*Q**2 + H(3)**2*R**2))
        
        if(iCell == 1)then
          P = nCells(1)
          FFT_XAI_X%Data_3D(nCells(1)+1, jCell, 1) &
            = - product(H) &
                / (4.0_KDR * Pi &
                   * sqrt(H(1)**2*P**2 + H(2)**2*Q**2 + H(3)**2*R**2))
        end if
        
        if(jCell == 1)then
          Q = nCells(2)
          FFT_XII_X%Data_1D(iCell) &
            = - product(H) &
                / (4.0_KDR * Pi &
                   * sqrt(H(1)**2*P**2 + H(2)**2*Q**2 + H(3)**2*R**2))
        end if
        
        if(jCell == 1 .and. iCell == 1)then
          P = nCells(1)
          Q = nCells(2)
          FFT_XII_X%Data_1D(nCells(1)+1) &
            = - product(H) &
                / (4.0_KDR * Pi &
                   * sqrt(H(1)**2*P**2 + H(2)**2*Q**2 + H(3)**2*R**2))
        end if

      end do

      FFT_XAI_X%Data_3D(nCells(1)+2:,jCell,1) &
        = FFT_XAI_X%Data_3D(nCells(1):2:-1,jCell,1)

    end do

    FFT_XII_X%Data_1D(nCells(1)+2:) = FFT_XII_X%Data_1D(nCells(1):2:-1)

    !-- Transform
    call Compute(L%FFT_Forward(1))
    call Compute(FFT_XIA_X)
    call Compute(FFT_XAI_X)
    call Compute(FFT_XII_X)
    
  end subroutine ComputeGreensFunction_X


  subroutine ComputeGreensFunction_Y( &
               L, FFT_XIA_X, FFT_XAI_X, FFT_XII_X, FFT_XYI_Y, nCells, &
               RankPosition_Y)

    type(LaplacianIsolated_FFT_Form), intent(inout) :: &
      L
    type(FFT_Base), intent(inout) :: &
      FFT_XIA_X, FFT_XAI_X, FFT_XII_X, FFT_XYI_Y
    integer(KDI), dimension(3), intent(in) :: &
      nCells
    integer(KDI), intent(in) :: &
      RankPosition_Y

    integer(KDI) :: &
      jCell, kCell, &
      oCell
    integer(KDI), dimension(3) :: &
      PW, PH
    logical(KDL) :: &
      Forward

    !-- The triplet in the comments and/or variables below indicate whether 
    !   the data pillar is an Active (A) or Inactive (I) in each of the 
    !   three-dimensions. The last letter indicates the direction of the 
    !   pillar (i.e. X, Y, or Z).

    PW = L%PillarWidth
    PH = L%PillarHeight
    
    Forward = .true.
    
    !-- Transpose to Y pillar
    L%FFT_Forward(2)%Data_3D &
      = Transpose(L%Communicator_XY, Forward, L%FFT_Forward(1)%Data_3D)
    do kCell = 1, PH(2)
      do jCell = 1, PW(2)
        oCell = RankPosition_Y * PW(1)
        L%FFT_Forward(2)%Data_3D(nCells(2)+1,jCell,kCell) &
          = FFT_XIA_X%Data_3D(oCell+kCell,1,jCell)
        oCell = oCell + nCells(1)
        L%FFT_Forward(2)%Data_3D(nCells(2)+1,jCell,PH(2)+kCell) &
          = FFT_XIA_X%Data_3D(oCell+kCell,1,jCell)
      end do
    end do
    L%FFT_Forward(2)%Data_3D(nCells(2)+2:,:,:) &
      = L%FFT_Forward(2)%Data_3D(nCells(2):2:-1,:,:)
    
    FFT_XYI_Y%Data_3D &
      = Transpose(L%Communicator_XY, Forward, FFT_XAI_X%Data_3D)
    do kCell = 1, PH(2)
      oCell = RankPosition_Y * PW(1)
      FFT_XYI_Y%Data_3D(nCells(2)+1,1,kCell) &
        = FFT_XII_X%Data_1D(oCell+kCell)
      oCell = oCell + nCells(1)
      FFT_XYI_Y%Data_3D(nCells(2)+1,1,PH(2)+kCell) &
        = FFT_XII_X%Data_1D(oCell+kCell)
    end do
    FFT_XYI_Y%Data_3D(nCells(2)+2:,:,:) = FFT_XYI_Y%Data_3D(nCells(2):2:-1,:,:)

    !-- Transform
    call Compute(L%FFT_Forward(2))
    call Compute(FFT_XYI_Y)
    
  end subroutine ComputeGreensFunction_Y


  subroutine ComputeGreensFunction_Z(L, FFT_XYI_Y, nCells, RankPosition_Z)

    type(LaplacianIsolated_FFT_Form), intent(inout) :: &
      L
    type(FFT_Base), intent(inout) :: &
      FFT_XYI_Y
    integer(KDI), dimension(3), intent(in) :: &
      nCells
    integer(KDI), intent(in) :: &
      RankPosition_Z

    integer(KDI) :: &
      jCell, kCell, &
      oCell
    integer(KDI), dimension(3) :: &
      PW, PH
    logical(KDL) :: &
      Forward 

    !-- The triplet in the comments and/or variables below indicate whether 
    !   the data pillar is an Active (A) or Inactive (I) in each of the 
    !   three-dimensions. The last letter indicates the direction of the 
    !   pillar (i.e. X, Y, or Z).

    PW = L%PillarWidth
    PH = L%PillarHeight
    
    Forward = .true.
    
    !-- Transpose to Z Pillar
    L%FFT_Forward(3)%Data_3D &
      = Transpose(L%Communicator_YZ, Forward, L%FFT_Forward(2)%Data_3D)
    do kCell = 1, PH(3)
      do jCell = 1, PW(3)
        oCell = RankPosition_Z * PW(2)
        L%FFT_Forward(3)%Data_3D(nCells(3)+1,jCell,kCell) &
          = FFT_XYI_Y%Data_3D(oCell+kCell,1,jCell)
        L%FFT_Forward(3)%Data_3D(nCells(3)+1,PW(3)+jCell,kCell) &
          = FFT_XYI_Y%Data_3D(oCell+kCell,1,PH(2)+jCell)
        oCell = oCell + nCells(2)
        L%FFT_Forward(3)%Data_3D(nCells(3)+1,jCell,PH(3)+kCell) &
          = FFT_XYI_Y%Data_3D(oCell+kCell,1,jCell)
        L%FFT_Forward(3)%Data_3D(nCells(3)+1,PW(3)+jCell,PH(3)+kCell) &
          = FFT_XYI_Y%Data_3D(oCell+kCell,1,PH(2)+jCell)
      end do
    end do
    L%FFT_Forward(3)%Data_3D(nCells(3)+2:,:,:) &
      = L%FFT_Forward(3)%Data_3D(nCells(3):2:-1,:,:)

    !-- Transform
    call Compute(L%FFT_Forward(3))
    
  end subroutine ComputeGreensFunction_Z


end module LaplacianIsolated_FFT__Form
