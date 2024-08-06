program Transpose_Command_Test

  use Basics
  use PSPFFT
  use Transpose_Command
  
  implicit none 

  type ( CommunicatorForm ), allocatable :: &
    C
  logical ( KDL ) :: & 
    Forward
  complex ( KDC ), dimension ( 2, 2, 2 ) :: & 
    SourcePillar 
  complex ( KDC ), dimension ( 2, 2, 2 ) :: & 
    TargetPillar 
  
  print *
  print *, 'Transpose Test'
  print * 

  allocate ( C )
  call C % Initialize ( ) 
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetDisplayRank ( 0 )


  deallocate ( C )

end program Transpose_Command_Test
