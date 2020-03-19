SUBROUTINE  fourier(SIZEB,conc,sigmaxy,gridx,gridy,deltax,deltay,concF,sigmaxyF,MODES, WAVEXX, WAVEYY)
IMPLICIT NONE 

INTEGER,   PARAMETER  :: dp  = kind(1.0D0)
INTEGER*8, INTENT(IN) :: GRIDX, GRIDY, MODES
REAL*8,    INTENT(IN) :: DELTAX, DELTAY, SIZEB
REAL(dp),  PARAMETER  :: TWOPI = 6.283185307179586476925286766559005768394
REAL(dp),  PARAMETER  :: PI    = 3.141592653589793238462643383279502884197
REAL*8,    DIMENSION(GRIDX,GRIDY), INTENT(IN)  :: conc, sigmaxy
REAL*8,    DIMENSION(MODES,modes), INTENT(OUT) :: concF, sigmaxyF
REAL*8     concFSIN, concFCOS, sigmaxyFSIN, sigmaxyFCOS
REAL(dp),  DIMENSION(MODES), INTENT(IN) :: WAVEXX, WAVEYY
REAL*8 X, Y
INTEGER I, J, P, L, W

concF=0.0
sigmaxyF=0.0
concFSIN=0.0
sigmaxyFSIN=0.0
concFCOS=0.0
sigmaxyFCOS=0.0

DO P=1, MODES
	DO L=1, MODES
		concFCOS    = 0.0
		sigmaxyFCOS = 0.0
		concFSIN    = 0.0
		sigmaxyFSIN = 0.0
		DO I=1, GRIDX
			DO J=1, GRIDY
				concFCOS   = concFCOS    + CONC(I,J)*COS( WAVEXX(P)*(SIZEB/GRIDX)*(I-1)+WAVEYY(L)*(SIZEB/GRIDY)*(J-1) )
				concFSIN   = concFSIN    - CONC(I,J)*SIN( WAVEXX(P)*(SIZEB/GRIDX)*(I-1)+WAVEYY(L)*(SIZEB/GRIDY)*(J-1) )
				sigmaxyFCOS= sigmaxyFCOS + SIGMAXY(I,J)*COS( WAVEXX(P)*(SIZEB/GRIDX)*(I-1)+WAVEYY(L)*(SIZEB/GRIDY)*(J-1) )
				sigmaxyFSIN= sigmaxyFSIN - SIGMAXY(I,J)*SIN( WAVEXX(P)*(SIZEB/GRIDX)*(I-1)+WAVEYY(L)*(SIZEB/GRIDY)*(J-1) )
			END DO
		END DO
		concF   (P,L) = ((concFCOS**2.0 + concFSIN**2.0) **0.5)/(gridx*gridy)
		sigmaxyF(P,L) = ((sigmaxyFCOS**2.0  + sigmaxyFSIN**2.0 ) **0.5)/(gridx*gridy)
	END DO
END DO


END SUBROUTINE fourier
 

