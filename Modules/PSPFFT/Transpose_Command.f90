!-- This module provides abstraction to transpose the pillar storage in 
!   in parallel over the whole computational domain.  

module Transpose_Command

  use Basics

  implicit none
  
  public :: &
    Transpose
  
  interface Transpose
    module procedure TransposeComplexDoublePillar
  end interface Transpose

contains

  
  function TransposeComplexDoublePillar( &
             C, Forward, SourcePillar) result(TargetPillar)
               
    type(CommunicatorForm), intent(in) :: &
      C
    logical(KL), intent(in) :: &
      Forward
    complex(KC), dimension(:,:,:), intent(in) :: &
      SourcePillar
    complex(KC), dimension(:,:,:), allocatable :: &
      TargetPillar
      
    integer(KI) :: &
      iRank, &
      iDatum, jData, kData, &
      iBuffer, &
      iSubpillar, &
      oData1, oData3, &
      oBuffer, &
      oPlane, &
      oLine, &
      ChunkSize, &
      HeightFactor
    integer(KI), dimension(3) :: &
      SourceShape, TargetShape
    complex(KC), dimension(:), allocatable :: &
      SendBuffer, &
      ReceiveBuffer
    
    !-- Pillars have their full dimension as the first index so transpose
    !-- becomes a permutation
    
    SourceShape = shape(SourcePillar)
    
    allocate(SendBuffer(size(SourcePillar)))
    allocate(ReceiveBuffer(size(SourcePillar)))
    
    if(Forward)then
      TargetShape(1) = SourceShape(2) * C%Size * 2
      TargetShape(2) = SourceShape(3)
      TargetShape(3) = SourceShape(1) / C%Size
    else !-- Reverse
      TargetShape(1) = SourceShape(3) * C%Size
      TargetShape(2) = SourceShape(1) / 2 / C%Size
      TargetShape(3) = SourceShape(2)
    end if
    
    allocate(TargetPillar(TargetShape(1),TargetShape(2),TargetShape(3)))
    call Clear(TargetPillar)
    
    if(Forward)then
      ChunkSize = product(SourceShape) / C%Size 
      HeightFactor = 1
    else
      ChunkSize = product(TargetShape) / C%Size 
      HeightFactor = 2
    end if
    
    do iRank = 0, C%Size-1

      !-- portion from "active" half of pillar
      oBuffer = iRank * ChunkSize
      oData1 = iRank * SourceShape(1) / (C%Size * 2)
      oData3 = 0
      SendBuffer(oBuffer+1:oBuffer+ChunkSize/2) &
        = reshape( &
            SourcePillar(oData1+1:oData1+SourceShape(1)/(C%Size*2),:, &
                         oData3+1:oData3+SourceShape(3)/HeightFactor), &
            (/ ChunkSize / 2 /))
      
      !-- portion from "inactive" half of pillar
      oBuffer = oBuffer + ChunkSize / 2
      if(Forward)then
        oData1 = oData1 + SourceShape(1) / 2
      else
        oData3 = SourceShape(3) / 2
      end if
      SendBuffer(oBuffer+1:oBuffer+ChunkSize/2) &
        = reshape( &
            SourcePillar(oData1+1:oData1+SourceShape(1)/(C%Size*2),:, &
                         oData3+1:oData3+SourceShape(3)/HeightFactor), &
            (/ ChunkSize / 2 /))

    end do
  
    call AllToAll(SendBuffer, C, ChunkSize, ReceiveBuffer)
    
    if(Forward)then
      iBuffer = 0
      do iRank = 0, C%Size-1
        oData3 = 0
        oData1 = iRank * TargetShape(1) / (C%Size * 2)
        do iSubpillar = 1, 2
          oPlane = TargetShape(3)/2 * TargetShape(1)/(C%Size*2)
          oLine = TargetShape(3)/2
          oBuffer = (iSubpillar - 1) * oPlane * TargetShape(2) &
                    + iRank * product(TargetShape) / (C%Size*2)
          !-- experiment shows that accessing TargetPillar contiguously gives
          !   better performance on Jaguar
          do kData = 1, TargetShape(3)/2
            do jData = 1, TargetShape(2)
              do iDatum = 1, TargetShape(1) / (C%Size * 2)
                iBuffer = oBuffer + (jData-1)*oPlane + (iDatum-1)*oLine + kData
                TargetPillar(oData1+iDatum,jData,oData3+kData) &
                  = ReceiveBuffer(iBuffer)
              end do
            end do
          end do
          oData3 = TargetShape(3) / 2
        end do
      end do
    else !-- Reverse
      iBuffer = 0
      do iRank = 0, C%Size-1
        oData1 = iRank * TargetShape(1) / (C%Size * 2)
        do iSubpillar = 1, 2
          !-- experiment shows that accessing ReceiveBuffer contiguously gives
          !   better performance on Jaguar
          do iDatum = 1, TargetShape(1) / (C%Size * 2)
            do kData = 1, TargetShape(3) 
              do jData = 1, TargetShape(2)
                iBuffer = iBuffer + 1
                TargetPillar(oData1+iDatum,jData,kData) &
                  = ReceiveBuffer(iBuffer)
              end do
            end do
          end do
          oData1 = oData1 + TargetShape(1) / 2  
        end do
      end do
    end if
    
    deallocate(SendBuffer)
    deallocate(ReceiveBuffer)
    
    !-- TargetPillar as result array is automatically deallocated 
    
  end function TransposeComplexDoublePillar


end module Transpose_Command
