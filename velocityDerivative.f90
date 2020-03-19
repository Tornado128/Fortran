  SUBROUTINE velocityDerivative(UXF, UYF, DELTAX, DELTAY, GRIDX, GRIDY, delUXX, delUYY, delUYX, delUXY)
  IMPLICIT NONE
 ! this subroutine calculates the first and second dervative of U_x and U_y (Div(V) and Grad(V))
integer, parameter :: dp  = kind(1.0D0)
INTEGER*8, INTENT(IN) :: GRIDX, GRIDY
REAL*8,  INTENT(IN) :: DELTAX, DELTAY
REAL(dp), DIMENSION(GRIDX,GRIDY), INTENT(IN):: UXF, UYF
REAL*8, DIMENSION(GRIDX,GRIDY), INTENT(OUT) :: delUXY, delUYX, delUXX, delUYY
Integer I, J

DO I=1,GRIDX
	DO J=1,GRIDY
		IF ((I .GT. 1) .AND. (I .LT. GRIDX)) THEN
			delUXX(I,J)=(UXF(I+1,J)-UXF(I-1,J))/(DELTAX+DELTAX)		! d(U_x)/dx
	 	END IF
		delUXX(GRIDX,J) = (UXF(1,J)-UXF(GRIDX-1,J))/(DELTAX+DELTAX)		
		delUXX(1,J)     = (UXF(2,J)-UXF(GRIDX,J))/(DELTAX+DELTAX)		
		delUYY(I,J)     = - delUXX(I,J)
		
		IF ((J .GT. 1) .AND. (J .LT. GRIDY)) THEN!
		        delUXY(I,J)=(UXF(I,J+1)-UXF(I,J-1))/(DELTAY+DELTAY)
		END IF
		IF ((I .GT. 1) .AND. (I .LT. GRIDX)) THEN
		        delUYX(I,J)=(UYF(I+1,J)-UYF(I-1,J))/(DELTAX+DELTAX)
		END IF

		delUXY(I,GRIDY) = (UXF(I,1)-UXF(I,GRIDY-1))/(DELTAY+DELTAY)
		delUXY(I,1)     = (UXF(I,2)-UXF(I,GRIDY))/(DELTAY+DELTAY)
		delUYX(GRIDX,J) = (UYF(1,J)-UYF(GRIDX-1,J))/(DELTAX+DELTAX)
		delUYX(1,J)     = (UYF(2,J)-UYF(GRIDX,J))/(DELTAX+DELTAX)
	END DO
END DO



!DO I=1,GRIDX
!	DO J=1,GRIDY
!
!		IF (UXF(I,J) > 0.0) THEN
!			IF (I .GT. 2) THEN
!				delUXX(I,J)=( 3.0*UXF(I,J)-4.0*UXF(I-1,J)+UXF(I-2,J) )/(2.0*DELTAX)
!			END IF
!			delUXX(2,J)=( 3.0*UXF(2,J)-4.0*UXF(1,J)+UXF(GRIDX,J) )/(2.0*DELTAX)
!			delUXX(1,J)=( 3.0*UXF(1,J)-4.0*UXF(GRIDX,J)+UXF(GRIDX-1,J) )/(2.0*DELTAX)
!			IF (J .GT. 2) THEN
!				delUXY(I,J)=( 3.0*UXF(I,J)-4.0*UXF(I,J-1)+UXF(I,J-2) )/(2.0*DELTAY)
!			END IF
!			delUXY(I,2)=( 3.0*UXF(I,2)-4.0*UXF(I,1)+UXF(I,GRIDY) )/(2.0*DELTAY)
!			delUXY(I,1)=( 3.0*UXF(I,1)-4.0*UXF(I,GRIDY)+UXF(I,GRIDY-1) )/(2.0*DELTAY)
!		END IF
!		IF (UXF(I,J) < 0.0) THEN
!			IF (I .LT. GRIDX-1) THEN
!				delUXX(I,J)=( -UXF(I+2,J)+4.0*UXF(I+1,J)-3.0*UXF(I,J) )/(2.0*DELTAX)
!			END IF
!			delUXX(GRIDX,J)   = ( -UXF(2,J)+4.0*UXF(1,J)-3.0*UXF(GRIDX,J) )/(2.0*DELTAX)
!			delUXX(GRIDX-1,J) = ( -UXF(1,J)+4.0*UXF(GRIDX,J)-3.0*UXF(GRIDX-1,J) )/(2.0*DELTAX)
!			IF (J .LT. GRIDY-1) THEN
!				delUXY(I,J)=( -UXF(I,J+2)+4.0*UXF(I,J+1)-3.0*UXF(I,J) )/(2.0*DELTAY)
!			END IF
!			delUXY(I,GRIDY)=( -UXF(I,2)+4.0*UXF(I,1)-3.0*UXF(I,GRIDY) )/(2.0*DELTAY)
!			delUXY(I,GRIDY-1)=( -UXF(I,1)+4.0*UXF(I,GRIDY)-3.0*UXF(I,GRIDY-1) )/(2.0*DELTAY)
!		END IF 
!		
!
!		IF (UYF(I,J) > 0.0) THEN
!			IF (J .GT. 2) THEN
!				delUYY(I,J)=( 3.0*UYF(I,J)-4.0*UYF(I,J-1)+UYF(I,J-2) )/(2.0*DELTAY)
!			END IF
!			delUYY(I,2)=( 3.0*UYF(I,2)-4.0*UYF(I,1)+UYF(I,GRIDY) )/(2.0*DELTAY)
!			delUYY(I,1)=( 3.0*UYF(I,1)-4.0*UYF(I,GRIDY)+UYF(I,GRIDY-1) )/(2.0*DELTAY)
!			IF (I .GT. 2) THEN
!				delUYX(I,J)=( 3.0*UYF(I,J)-4.0*UYF(I-1,J)+UYF(I-2,J) )/(2.0*DELTAX)
!			END IF
!			delUYX(2,J)=( 3.0*UYF(2,J)-4.0*UYF(1,J)+UYF(GRIDX,J) )/(2.0*DELTAX)
!			delUYX(1,J)=( 3.0*UYF(1,J)-4.0*UYF(GRIDX,J)+UYF(GRIDX-1,J) )/(2.0*DELTAX)
!		END IF
!		IF (UYF(I,J) < 0.0) THEN
!			IF (J .LT. GRIDY-1) THEN
!				delUYY(I,J)=( -UYF(I,J+2)+4.0*UYF(I,J+1)-3.0*UYF(I,J) )/(2.0*DELTAY)
!			END IF
!			delUYY(I,GRIDY)  = ( -UYF(I,2)+4.0*UYF(I,1)-3.0*UYF(I,GRIDY) )/(2.0*DELTAY)
!			delUYY(I,GRIDY-1)= ( -UYF(I,1)+4.0*UYF(I,GRIDY)-3.0*UYF(I,GRIDY-1) )/(2.0*DELTAY)
!			IF (I .LT. GRIDX-1) THEN
!				delUYX(I,J)=( -UYF(I+2,J)+4.0*UYF(I+1,J)-3.0*UYF(I,J) )/(2.0*DELTAX)
!			END IF
!			delUYX(GRIDX,J)=( -UYF(2,J)+4.0*UYF(1,J)-3.0*UYF(GRIDX,J) )/(2.0*DELTAX)
!			delUYX(GRIDX-1,J)=( -UYF(1,J)+4.0*UYF(GRIDX,J)-3.0*UYF(GRIDX-1,J) )/(2.0*DELTAX)
!		END IF
!
!	END DO
 !END DO


END SUBROUTINE velocityDerivative
