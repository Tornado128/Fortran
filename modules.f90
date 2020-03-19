
module rutina
  
  INTERFACE
     SUBROUTINE ggem_regularized(delSIGMAX,delSIGMAY,SIZEB, DELTAX, DELTAY, GRIDX, GRIDY, UXF, UYF)
       implicit none
integer, parameter :: dp  = kind(1.0D0)
  integer, parameter :: dpc = kind( (1.d0,1.d0) )
       INTEGER, INTENT(IN) :: GRIDX, GRIDY
       REAL(dp), INTENT(IN) :: SIZEB, DELTAX, DELTAY 
	REAL(dp), DIMENSION(GRIDX,GRIDY), INTENT(IN):: delSIGMAX, delSIGMAY
  REAL(dp), DIMENSION(GRIDX,GRIDY), INTENT(OUT):: UXF, UYF
     end subroutine ggem_regularized
  END INTERFACE
  


end module rutina
