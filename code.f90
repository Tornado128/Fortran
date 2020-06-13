PROGRAM code
IMPLICIT NONE
			!constant zero shear viscosity
 
INTEGER, 	PARAMETER 				:: DP    = kind(1.0D0)
INTEGER, 	PARAMETER 				:: DPC   = kind( (1.d0,1.d0) )
REAL*8,  	PARAMETER 				:: SIZEB =50
REAL*8,	  	PARAMETER 				:: ALPHA =-1.0
INTEGER*8, 	PARAMETER 				:: GRIDX =64
INTEGER*8, 	PARAMETER 				:: GRIDY =64
INTEGER*8, 	PARAMETER 				:: stepN =45000
INTEGER*8,	PARAMETER				:: TRACER=400
INTEGER*8,	PARAMETER				:: sqrtTRACER=20
REAL*8,  	PARAMETER 				:: DELTAt=0.01
INTEGER*8, 	PARAMETER 				:: PHID  =60
INTEGER*8,	PARAMETER 				:: MODES =15
REAL(DP), 	PARAMETER 				:: TWOPI = 6.283185307179586476925286766559005768394
REAL(DP), 	PARAMETER 				:: PI    = 3.1415926535897932384626433832795028841977
REAL*8,  	PARAMETER  				:: Dr=0.005, Dt=0.45
REAL*8,  	PARAMETER 				:: GAMM=1.0

REAL*8,		PARAMETER				:: ZETA	 =0.8
REAL*8,		PARAMETER				:: De	 =0.8

INTEGER*8	, PARAMETER 				:: PW=1, PWnx=2, PWny=3, PWux=4, PWuy=5, PW2=6, PW3=7, PW4=8, pw5=9, pw6=10, pw7=11
INTEGER*8	, PARAMETER				:: pw8=12, i8=13, psp=14, PSP3=15, PSP4=16, PSP5=17, PSP6=18
INTEGER*8	, PARAMETER 				:: H1=19, H2=20, H3=21, H4=22, H5=23, H6=24, I1=25, I2=26,I3=27,I4=28,I5=29,I6=30
INTEGER*8	, PARAMETER				:: I7=31, i9=32, PSP1=33, PSP20=34, pwee1=35, pwee2=36, PSPPW=37
INTEGER*8	, PARAMETER 				:: g1=38, g2=39, g3=40, g4=41, g5=42, g6=43, ghg=44, FGF=45, FGF1=46, FGF2=47
INTEGER*8	, PARAMETER				:: FGF3=48, FGF4=49, FGF5=50, PSP2=51, kjk1=52, kjk2=53, maxim=54, W99=55
INTEGER*8	, PARAMETER				:: W1=56, W2=57, W3=58, W4=59, W5=60, PSP7=61, W12=62, W13=63, W14=64, W9=65
INTEGER*8	, PARAMETER				:: W32=66, W33=67, W34=68, W52=69, W53=70, W54=71, h7=72, W7=73, L1=74, L2=75
INTEGER*8	, PARAMETER				:: L3=76, L4=77, L5=78, sty=79, STY1=80, STY2=81, STY3=82, STY4=83, STT=84,psp8=85
INTEGER*8	, PARAMETER				:: kjk3=86, kjk4=87, kjk5=88, kjk6=89, kjk7=90, psp9=91, psp10=92, psp11=93, psp12=94
INTEGER*8	, PARAMETER				:: psp13=95, psp14=96, psp15=97

REAL*8		, PARAMETER 				:: DELTAX=SIZEB/(GRIDX)
REAL*8		, PARAMETER 				:: DELTAY=SIZEB/(GRIDY)
REAL*8		, PARAMETER 				:: TOLERANCE=0.01
INTEGER*8	, PARAMETER 				:: ITERATION=100

INTEGER*8 						   I, J, A, P, k, TIME, NUMER, kkk,za, L, W, auxi, time1, delTsampling, t
REAL(DP)	, DIMENSION(GRIDX,GRIDY,PHID)  		:: ERROR, delXPSI, delX2PSI, delYPSI, delY2PSI, delO2PSI, delOPSI
REAL(DP)	, DIMENSION(GRIDX,GRIDY) 		:: SIGMAXY, SIGMAXX, SIGMAYY, CONC, ORIX, ORIY, ANG, delSIGMAY
REAL(DP)	, DIMENSION(GRIDX,GRIDY) 		:: delSIGMAX, UXF, UYF , DELdotNCX, DELdotNCY, sigmaXYm, sigmaXXm, TAUpxyM, TAUpxxM
REAL(DP)	, DIMENSION(GRIDX,GRIDY) 		:: TAUpXX, TAUpYY, TAUpXY, delTAUpXXx, delTAUpXYx, delTAUpYYx
REAL(DP)	, DIMENSION(GRIDX,GRIDY) 		:: delTAUpXXy, delTAUpXYy, delTAUpYYy, ERRORtauP, UelongationX, UshearX
REAL(DP)	, DIMENSION(MODES,MODES) 		:: sigmaxyF, concF
REAL(DP)	, DIMENSION(INT( GRIDX/2 + 1),GRIDY) 	:: concFF, stressXXff, stressYYff, stressXYff, sigmaXYff, ORIXff, ORIYff
REAL(DP)	, DIMENSION(INT( GRIDX/2 + 1),GRIDY) 	:: sigmaXXff, sigmaYYff, TauPxyff, TauPxxff, TauPyyff, FNSCff, TRACEff, &
							   CONVECxxFF, CONVECxyFF, LToldrXXff, LToldrXYff, orientXXff, orientXYff, &
							   elongationXff, shearXff
REAL(DP)	, DIMENSION(int( GRIDX/2 + 1),GRIDY) 	:: UxFOUR, UyFOUR, delSIGMAxFFF, delSIGMAyFFF, contild
REAL(dp)	, DIMENSION(GRIDX,GRIDY) 		:: delUXY, delUYX, delUXX, delUYY, CONTINUITY, KOM, CONVECxx, CONVECxy, &
							   LToldrXX, LToldrXY, orientXX, orientXY
REAL(dp)	, DIMENSION(GRIDX,GRIDY) 		:: PPYY, PPXY, PPXX, TAUpXX1, TAUpYY1, TAUpXY1, TAUpXXi, TAUpYYi, TAUpXYi
REAL(dp)	, DIMENSION(GRIDX,GRIDY,PHID) 		:: PSI1, PSI, PSIi, PP, PPP
REAL(dp)	, DIMENSION(MODES) 			:: WAVEXX, WAVEYY
REAL(DP)	, DIMENSION(TRACER,stepN+1)		:: Xtracer, Ytracer, UXtracer, UYtracer
REAL(DP)	, DIMENSION(TRACER)		  	:: Xinitial, Yinitial
REAL(DP)	, DIMENSION(TRACER,stepN+1)		:: XtracerNOp, YtracerNOp
REAL(DP)	, DIMENSION(stepN/2)			:: MSD, Nsample


REAL(DP)  						D1, D2, D3, D4, D5, D6, D7, D8, D9, FXY, ANGLE, Pphi, D11, D12, uDOTn, CONAV, &
							delUxyV, delUv, delUdif, addPSI3
REAL*8  						DA, CONTIN, wat, VEL, ENTROPY, ERRORM1, ERRORM2, nxyN
REAL*8  						RAD, AAA, IPI, mixing, uSHEARdotELONG
REAL(DP) 						X, Y, ERRORM, addPSI, addPSI1, addPSI2, addCON, addCON1, addCON2, POWER
REAL*8 							time_begin, time_end, time_begin1, time_end1, time_begin2, time_end2
REAL*8							time_begin3, time_end3, elapsed, elapsed1, elapsed2, elapsed3
REAL*8							dx, dy, posiX, posiY
integer*8						nume, Xbelow, Ybelow, Xabove, Yabove, XXbelow, YYbelow, XXabove, YYabove

REAL, DIMENSION(2) 					:: tarray
REAL 							:: result
INTEGER*8						:: CLOCK
real(kind=8)						:: random
integer(kind=4)						:: n
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
 CALL RANDOM_SEED(SIZE=n)
 ALLOCATE(seed(n))

 CALL 							   ETIME(tarray, result)
! print *, result
! print *, tarray(1)
! print *, tarray(2)
print*, 'H=', ZETA,'     ', 'De=', De,'     ', 'Dr=', Dr,'     ', 'Dt=', Dt

 TAUpXX = 0
 TAUpYY = 0
 TAUpXY = 0
time1=0
W=0
DO P=1,MODES
	DO L=1,MODES
		W = W + 1

		WAVEXX=0!MODES*TWOPI/SIZEB
		WAVEYY=0!MODES*TWOPI/SIZEB
		IF (W .NE. 1) THEN
			WAVEXX(P)=(P-1)*TWOPI/SIZEB
			WAVEYY(L)=(L-1)*TWOPI/SIZEB
		END IF
		D1=2*RAND()-1
		D2=2*RAND()-1
		D3=2*RAND()-1
		D4=2*RAND()-1
		D5=2*RAND()-1
		D6=2*RAND()-1
		D7=2*RAND()-1
		D9=0.01*(2*RAND()-1)
		D8=(TWOPI/4.0)*RAND()
		DO I=1,GRIDX
			X=(SIZEB/(GRIDX))*(I-1)

			DO J=1,GRIDY
				Y  = (SIZEB/(GRIDY))*(J-1)

				FXY= COS( WAVEXX(P)*X+WAVEYY(L)*Y+D8 )
				DO A=1,PHID
					ANGLE = TWOPI*(A-1)/(PHID-1)
					Pphi  = D1*SIN(ANGLE)+D2*COS(ANGLE) + D3*SIN(2.0*ANGLE)+ D4*COS(2.0*ANGLE) &
					 + D5*SIN(3.0*ANGLE) + D6*COS(3.0*ANGLE) + D7
					PSI(I,J,A) = PSI(I,J,A) + (D9*FXY*Pphi)
				END DO
			END DO
		END DO
	END DO
END DO

 !addPSI1=SUM(PSI)
 !PSI = PSI - addPSI1/(gridx*gridy*phid)
 !addPSI2=SUM(PSI)

 PSI = PSI/TWOPI + 1.0/TWOPI

 CALL   	estimation(ALPHA, GRIDX, GRIDY, DELTAX, DELTAY, PHID, PSI, SIGMAXX, SIGMAYY, SIGMAXY, &
	    	CONC, ORIX, ORIY, delXPSI, delX2PSI, delYPSI, delY2PSI, delO2PSI, delOPSI, delSIGMAX, &
		delSIGMAY, DELdotNCX, DELdotNCY, ENTROPY )

 CALL 		ggm(MODES,SIZEB, DELTAX, DELTAY, GRIDX, GRIDY, UXF, UYF, UxFOUR, UyFOUR, SIGMAXY, SIGMAXX, &
		SIGMAYY, TAUpXY, TAUpXX, TAUpYY, contild, elongationXff, shearXff, UelongationX, UshearX, &
		sigmaXYm, sigmaXXm, TAUpxyM, TAUpxxM, ZETA)

 CALL 		velocityDerivative(UXF, UYF, DELTAX, DELTAY, GRIDX, GRIDY, delUXX, delUYY, delUYX, delUXY)

 CALL 		tauPderivative(ZETA,De,UXF,UYF,GRIDX,GRIDY,DELTAX,DELTAY,delUXX, delUYY, delUYX, delUXY, &
		TAUpXX,TAUpYY,TAUpXY, delTAUpXXx, delTAUpXYx, delTAUpYYx, delTAUpXXy, delTAUpXYy, delTAUpYYy)

 		!addCON1=SUM(CONC)
		!CONC=(CONC - (addCON1-gridx*gridy)/(gridx*gridy) )
		CONTINUITY = ( delUXX+delUYY ) /( (UXF**2.0 + UYF**2.0)**0.5 )

 CALL 		fftwCheck(conc,sigmaxy,TAUpXY,SIZEB, DELTAX, DELTAY, GRIDX, GRIDY, concFF &
	     	,sigmaxx,TAUpXX, sigmayy, TAUpYY, stressXXff, stressYYff, stressXYff, mixing &
	     	,sigmaxxFF, sigmaxyFF, sigmayyFF, TAUpXYff, TAUpYYff, TAUpXXff, ORIX, ORIY, ORIXff &
	     	,ORIYff, TRACEff, FNSCff, LToldrXX, LToldrXY, LToldrXXff, LToldrXYff  &
		,CONVECxx, CONVECxy, CONVECxxFF, CONVECxyFF, orientXX, orientXY, orientXXff, orientXYff)

 CALL CHDIR("/home/yaser/Desktop/RESEARCH/NonLinearDynamics/OldroydB2constant/matlabForProcessingData")

 OPEN(h1    ,	FILE='nx.txt')
 OPEN(h2    ,	FILE='ny.txt')
 OPEN(h3    ,	FILE='CONCENTRATION.txt')
 OPEN(h4    ,	FILE='SIGMAXY1.txt')
 OPEN(h5    ,	FILE='ux.txt')
 OPEN(h6    ,	FILE='uy.txt')
 OPEN(h7    ,   FILE='deldotnc0.txt')
 OPEN(STY1  ,	FILE='delU1.txt')

 OPEN(kjk1  ,	FILE='fourierStress.txt')
 OPEN(PWEE2 ,	FILE='fourierTauP.txt')
 OPEN(PWEE1 ,	FILE='fourier.txt')

 OPEN(PSP   ,	FILE='totalCONC.txt'  )
 OPEN(PSP1  ,	FILE='totalUX2.txt'   )
 OPEN(PSP2  ,	FILE='totalUY2.txt'   )
 OPEN(PSP3  ,	FILE='totalNXX.txt'   )
 OPEN(PSP4  ,	FILE='totalNYY.txt'   )
 OPEN(PSP5  ,	FILE='totalSIGMA.txt' )
 OPEN(PSP11  ,	FILE='totalSIGMAxx.txt' )
 OPEN(PSP6  ,	FILE='totalTRACE.txt' )
 OPEN(PSP9  ,	FILE='totalTAUXX.txt' )
 OPEN(PSP10  ,	FILE='totalTAUYY.txt' )
 OPEN(PSP8  ,	FILE='totalTAUpXY.txt')
 OPEN(kjk3  ,	FILE='tracer.txt'     )
 OPEN(kjk4  ,	FILE='MSD.txt'	      )



 OPEN(MAXIM ,	FILE='CHECKING.txt')
 OPEN(W13   ,   FILE='SIGMAXX3.txt')
 OPEN(W14   ,   FILE='SIGMAXX4.txt')
 OPEN(W5    ,	FILE='TRACE1.txt')
 OPEN(W52   ,	FILE='TRACE2.txt')
 OPEN(W53   ,	FILE='TRACE3.txt')
 OPEN(W54   ,	FILE='TRACE4.txt')

 OPEN(PSP7  ,	FILE='totalENTROPY.txt')

 OPEN(W9    ,	FILE='COMPARISON.txt')

 open(STT   ,   FILE='orientationFF.txt')
 open(kjk5  ,	FILE='UshearX.txt')
 open(kjk6  ,	FILE='UelongationX.txt')
 open(kjk7  ,	FILE='shearDOTelongation.txt')
 open(psp12 ,	FILE='sigmaXYm.txt')
 open(psp13 ,	FILE='sigmaXXm.txt')
 open(psp14 ,	FILE='taupXYm.txt')
 open(psp15 ,	FILE='taupXXm.txt')
!sigmaXYm(I,J), sigmaXXm(I,J), TAUpxyM(I,J), TAUpxxM(I,J)
 DO I=1, GRIDX
 	DO J=1, GRIDY

		LToldrXX(I,J)=-TAUpXX(I,J)/(4.0*De) - delUXX(I,J)*ZETA/(2.0*De)

		convecXX(I,J)=-UXF(I,J)*delTAUpXXx(I,J)-UYF(I,J)*delTAUpXXy(I,J) + 2.0*TAUpXX(I,J)*delUXX(I,J)+ 2.0*TAUpXY(I,J)*delUYX(I,J)

		LToldrXY(I,J)=-TAUpXY(I,J)/(4.0*De) - ( delUXY(I,J)+delUYX(I,J) )*ZETA/(4.0*De)

		convecXY(I,J)=-delTAUpXYx(I,J)*UXF(I,J)-delTAUpXYy(I,J)*UYF(I,J) +TAUpXX(I,J)*delUXY(I,J) + delUYX(I,J)*TAUpYY(I,J)

		orientXX(I,J)=2.0*TAUpXX(I,J)*delUXX(I,J)+ 2.0*TAUpXY(I,J)*delUYX(I,J)
		orientXY(I,J)=TAUpXX(I,J)*delUXY(I,J) + delUYX(I,J)*TAUpYY(I,J)

	 	write (h1,	*)   	ORIX(I,J)
	 	write (h2,	*)   	ORIY(I,J)
	 	write (h3,	*)   	CONC(I,J)
		write (h4,	*)   	SIGMAXY	(I,J), SIGMAXX(I,J), SIGMAYY(I,J), SIGMAXX(I,J)+SIGMAYY(I,J), SIGMAXY(I,J)-TAUpXY(I,J), &
			 		SIGMAXX(I,J)-TAUpXX(I,J)
	 	write (h5,	*)   	UXF(I,J)
		write (h6,	*)   	UYF(I,J)
		WRITE (h7,	*)	-DELdotNCX(I,J)-DELdotNCY(I,J)
		WRITE (STY1,	*)  	SQRT( delUXX(I,J)*delUXX(I,J)+delUYY(I,J)*delUYY(I,J) ), &
					SQRT( delUXY(I,J)*delUXY(I,J)+delUYX(I,J)*delUYX(I,J) ), &
					delUXX(I,J)-delUYY(I,J)
		WRITE (W5,      *)	TAUpXX(I,J)+TAUpYY(I,J), TAUpXX(I,J)-TAUpYY(I,J),LToldrXXff(I,J),LToldrXYff(I,J), &
					CONVECxx(I,J),CONVECxy(I,J), TAUpXX(I,J), TAUpXY(I,J)

		WRITE (PSP,	*)  	CONC(I,J)
		WRITE (PSP1,	*) 	UXF(I,J)
		WRITE (PSP2,	*) 	UYF(I,J)
		WRITE (PSP3,	*) 	ORIX(I,J)
		WRITE (PSP4,	*) 	ORIY(I,J)
		WRITE (PSP5,	*) 	SIGMAXY(I,J)
		WRITE (PSP6,    *)	TAUpXX(I,J)+TAUpYY(I,J)
		WRITE (PSP8,	*)	TAUpXY(I,J)
		WRITE (PSP9,    *)	TAUpXX(I,J)
		WRITE (PSP10,	*)	TAUpYY(I,J)
		WRITE (kjk5,	*)	UshearX(I,J)
		WRITE (kjk6,	*)	UelongationX(I,J)
		WRITE (PSP11,	*)	SIGMAXX(I,J)
		WRITE (PSP12,	*)	sigmaXYm(I,J)
		WRITE (PSP13,	*)	sigmaXXm(I,J)
		WRITE (PSP14,	*)	TAUpxyM(I,J)
		WRITE (PSP15,	*)	TAUpxxM(I,J)


 	END DO
 END DO

 WRITE (PSP7,*) ENTROPY

WRITE (KJK1,*) stressXXff(1,1)  ,stressYYff(1,1),  stressXYff(1,1)
WRITE (KJK1,*) stressXXff(2,1)  ,stressYYff(2,1),  stressXYff(2,1)
WRITE (KJK1,*) stressXXff(1,2)  ,stressYYff(1,2),  stressXYff(1,2)
WRITE (KJK1,*) stressXXff(2,2)  ,stressYYff(2,2),  stressXYff(2,2)
WRITE (KJK1,*) stressXXff(3,1)  ,stressYYff(3,1),  stressXYff(3,1)
WRITE (KJK1,*) stressXXff(1,3)  ,stressYYff(1,3),  stressXYff(1,3)
WRITE (KJK1,*) stressXXff(3,2)  ,stressYYff(3,2),  stressXYff(3,2)
WRITE (KJK1,*) stressXXff(2,3)  ,stressYYff(2,3),  stressXYff(2,3)
WRITE (KJK1,*) stressXXff(3,3)  ,stressYYff(3,3),  stressXYff(3,3)
WRITE (KJK1,*) stressXXff(4,1)  ,stressYYff(4,1),  stressXYff(4,1)
WRITE (KJK1,*) stressXXff(1,4)  ,stressYYff(1,4),  stressXYff(1,4)
WRITE (KJK1,*) stressXXff(4,2)  ,stressYYff(4,2),  stressXYff(4,2)
WRITE (KJK1,*) stressXXff(2,4)  ,stressYYff(2,4),  stressXYff(2,4)
WRITE (KJK1,*) stressXXff(3,4)  ,stressYYff(3,4),  stressXYff(3,4)
WRITE (KJK1,*) stressXXff(4,4)  ,stressYYff(4,4),  stressXYff(4,4)
WRITE (KJK1,*) stressXXff(5,1)  ,stressYYff(5,1),  stressXYff(5,1)
WRITE (KJK1,*) stressXXff(5,3)  ,stressYYff(5,3),  stressXYff(5,3)
WRITE (KJK1,*) stressXXff(1,6)  ,stressYYff(1,6),  stressXYff(1,6)
WRITE (KJK1,*) stressXXff(7,1)  ,stressYYff(7,1),  stressXYff(7,1)
WRITE (KJK1,*) stressXXff(1,8)  ,stressYYff(1,8),  stressXYff(1,8)
WRITE (KJK1,*) stressXXff(5,5)  ,stressYYff(5,5),  stressXYff(5,5)
WRITE (KJK1,*) stressXXff(6,6)  ,stressYYff(6,6),  stressXYff(6,6)
WRITE (KJK1,*) stressXXff(10,5)  ,stressYYff(10,5),  stressXYff(10,5)
WRITE (KJK1,*) stressXXff(13,13)  ,stressYYff(13,13),  stressXYff(13,13)

WRITE (STT,*) orientXXff(1,1)  ,  orientXYff(1,1), UxFOUR(1,1), elongationXff(1,1), shearXff(1,1)
WRITE (STT,*) orientXXff(2,1)  ,  orientXYff(2,1), UxFOUR(2,1), elongationXff(2,1), shearXff(2,1)
WRITE (STT,*) orientXXff(1,2)  ,  orientXYff(1,2), UxFOUR(1,2),elongationXff(1,2), shearXff(1,2)
WRITE (STT,*) orientXXff(2,2)  ,  orientXYff(2,2), UxFOUR(2,2),elongationXff(2,2), shearXff(2,2)
WRITE (STT,*) orientXXff(3,1)  ,  orientXYff(3,1), UxFOUR(3,1),elongationXff(3,1), shearXff(3,1)
WRITE (STT,*) orientXXff(1,3)  ,  orientXYff(1,3), UxFOUR(1,3),elongationXff(1,3), shearXff(1,3)
WRITE (STT,*) orientXXff(3,2)  ,  orientXYff(3,2), UxFOUR(3,2),elongationXff(3,2), shearXff(3,2)
WRITE (STT,*) orientXXff(2,3)  ,  orientXYff(2,3), UxFOUR(2,3),elongationXff(2,3), shearXff(2,3)
WRITE (STT,*) orientXXff(3,3)  ,  orientXYff(3,3), UxFOUR(3,3),elongationXff(3,3), shearXff(3,3)
WRITE (STT,*) orientXXff(4,1)  ,  orientXYff(4,1), UxFOUR(4,1),elongationXff(4,1), shearXff(4,1)
WRITE (STT,*) orientXXff(1,4)  ,  orientXYff(1,4), UxFOUR(1,4),elongationXff(1,4), shearXff(1,4)
WRITE (STT,*) orientXXff(4,2)  ,  orientXYff(4,2), UxFOUR(4,2),elongationXff(4,2), shearXff(4,2)
WRITE (STT,*) orientXXff(2,4)  ,  orientXYff(2,4), UxFOUR(2,4),elongationXff(2,4), shearXff(2,4)
WRITE (STT,*) orientXXff(3,4)  ,  orientXYff(3,4), UxFOUR(3,4),elongationXff(3,4), shearXff(3,4)
WRITE (STT,*) orientXXff(4,4)  ,  orientXYff(4,4), UxFOUR(4,4),elongationXff(4,4), shearXff(4,4)
WRITE (STT,*) orientXXff(5,1)  ,  orientXYff(5,1), UxFOUR(5,1),elongationXff(5,1), shearXff(5,1)
WRITE (STT,*) orientXXff(5,3)  ,  orientXYff(5,3), UxFOUR(5,3),elongationXff(5,3), shearXff(5,3)
WRITE (STT,*) orientXXff(1,6)  ,  orientXYff(1,6), UxFOUR(1,6),elongationXff(1,6), shearXff(1,6)
WRITE (STT,*) orientXXff(7,1)  ,  orientXYff(7,1), UxFOUR(7,1),elongationXff(7,1), shearXff(7,1)
WRITE (STT,*) orientXXff(1,8)  ,  orientXYff(1,8), UxFOUR(1,8),elongationXff(1,8), shearXff(1,8)
WRITE (STT,*) orientXXff(5,5)  ,  orientXYff(5,5), UxFOUR(5,5),elongationXff(5,5), shearXff(5,5)
WRITE (STT,*) orientXXff(6,6)  ,  orientXYff(6,6), UxFOUR(6,6),elongationXff(6,6), shearXff(6,6)
WRITE (STT,*) orientXXff(10,5) ,  orientXYff(10,5), UxFOUR(10,5),elongationXff(10,5), shearXff(10,5)
WRITE (STT,*) orientXXff(13,13),  orientXYff(13,13), UxFOUR(13,13),elongationXff(13,13), shearXff(13,13)

WRITE (PWEE2,*) TauPxxff(1,1),TauPyyff(1,1),TauPXYff(1,1),convecxxFF(1,1),convecxyFF(1,1),LToldrXXff(1,1),LToldrXYff(1,1)
WRITE (PWEE2,*) TauPxxff(2,1),TauPyyff(2,1),TauPXYff(2,1),convecxxFF(2,1),convecxyFF(2,1),LToldrXXff(2,1),LToldrXYff(2,1)
WRITE (PWEE2,*) TauPxxff(1,2),TauPyyff(1,2),TauPXYff(1,2),convecxxFF(1,2),convecxyFF(1,2),LToldrXXff(1,2),LToldrXYff(1,2)
WRITE (PWEE2,*) TauPxxff(2,2),TauPyyff(2,2),TauPXYff(2,2),convecxxFF(2,2),convecxyFF(2,2),LToldrXXff(2,2),LToldrXYff(2,2)
WRITE (PWEE2,*) TauPxxff(3,1),TauPyyff(3,1),TauPXYff(3,1),convecxxFF(3,1),convecxyFF(3,1),LToldrXXff(3,1),LToldrXYff(3,1)
WRITE (PWEE2,*) TauPxxff(1,3),TauPyyff(1,3),TauPXYff(1,3),convecxxFF(1,3),convecxyFF(1,3),LToldrXXff(1,3),LToldrXYff(1,3)
WRITE (PWEE2,*) TauPxxff(3,2),TauPyyff(3,2),TauPXYff(3,2),convecxxFF(3,2),convecxyFF(3,2),LToldrXXff(3,2),LToldrXYff(3,2)
WRITE (PWEE2,*) TauPxxff(2,3),TauPyyff(2,3),TauPXYff(2,3),convecxxFF(2,3),convecxyFF(2,3),LToldrXXff(2,3),LToldrXYff(2,3)
WRITE (PWEE2,*) TauPxxff(3,3),TauPyyff(3,3),TauPXYff(3,3),convecxxFF(3,3),convecxyFF(3,3),LToldrXXff(3,3),LToldrXYff(3,3)
WRITE (PWEE2,*) TauPxxff(4,1),TauPyyff(4,1),TauPXYff(4,1),convecxxFF(4,1),convecxyFF(4,1),LToldrXXff(4,1),LToldrXYff(4,1)
WRITE (PWEE2,*) TauPxxff(1,4),TauPyyff(1,4),TauPXYff(1,4),convecxxFF(1,4),convecxyFF(1,4),LToldrXXff(1,4),LToldrXYff(1,4)
WRITE (PWEE2,*) TauPxxff(4,2),TauPyyff(4,2),TauPXYff(4,2),convecxxFF(4,2),convecxyFF(4,2),LToldrXXff(4,2),LToldrXYff(4,2)
WRITE (PWEE2,*) TauPxxff(2,4),TauPyyff(2,4),TauPXYff(2,4),convecxxFF(2,4),convecxyFF(2,4),LToldrXXff(2,4),LToldrXYff(2,4)
WRITE (PWEE2,*) TauPxxff(3,4),TauPyyff(3,4),TauPXYff(3,4),convecxxFF(3,4),convecxyFF(3,4),LToldrXXff(3,4),LToldrXYff(3,4)
WRITE (PWEE2,*) TauPxxff(4,4),TauPyyff(4,4),TauPXYff(4,4),convecxxFF(4,4),convecxyFF(4,4),LToldrXXff(4,4),LToldrXYff(4,4)
WRITE (PWEE2,*) TauPxxff(5,1),TauPyyff(5,1),TauPXYff(5,1),convecxxFF(4,1),convecxyFF(5,1),LToldrXXff(5,1),LToldrXYff(5,1)
WRITE (PWEE2,*) TauPxxff(5,3),TauPyyff(5,3),TauPXYff(5,3),convecxxFF(5,3),convecxyFF(5,3),LToldrXXff(5,3),LToldrXYff(5,3)
WRITE (PWEE2,*) TauPxxff(1,6),TauPyyff(1,6),TauPXYff(1,6),convecxxFF(1,6),convecxyFF(1,6),LToldrXXff(1,6),LToldrXYff(1,6)
WRITE (PWEE2,*) TauPxxff(7,1),TauPyyff(7,1),TauPXYff(7,1),convecxxFF(7,1),convecxyFF(7,1),LToldrXXff(7,1),LToldrXYff(7,1)
WRITE (PWEE2,*) TauPxxff(1,8),TauPyyff(1,8),TauPXYff(1,8),convecxxFF(1,8),convecxyFF(1,8),LToldrXXff(1,8),LToldrXYff(1,8)
WRITE (PWEE2,*) TauPxxff(5,5),TauPyyff(5,5),TauPXYff(5,5),convecxxFF(5,5),convecxyFF(5,5),LToldrXXff(5,5),LToldrXYff(5,5)
WRITE (PWEE2,*) TauPxxff(6,6),TauPyyff(6,6),TauPXYff(6,6),convecxxFF(6,6),convecxyFF(6,6),LToldrXXff(6,6),LToldrXYff(6,6)
WRITE (PWEE2,*) TauPxxff(10,5),TauPyyff(10,5),TauPXYff(10,5),convecxxFF(10,5),convecxyFF(10,5), 			&
		LToldrXXff(10,5),LToldrXYff(10,5)
WRITE (PWEE2,*) TauPxxff(13,13),TauPyyff(13,13),TauPXYff(13,13),convecxxFF(13,13),convecxyFF(13,13), 			&
		LToldrXXff(13,13),LToldrXYff(13,13)

WRITE (PWEE1,*) sigmaXXff(1,1)  ,sigmaYYff(1,1),  sigmaXYff(1,1), concFF(1,1), SQRT(UxFOUR(1,1)**2.0+UyFOUR(1,1)**2.0), &
		SQRT(ORIXff(1,1)**2.0+ORIYff(1,1)**2.0), TRACEff(1,1), FNSCff(1,1)
WRITE (PWEE1,*) sigmaXXff(2,1)  ,sigmaYYff(2,1),  sigmaXYff(2,1), concFF(2,1), SQRT(UxFOUR(2,1)**2.0+UyFOUR(2,1)**2.0), &
		SQRT(ORIXff(2,1)**2.0+ORIYff(2,1)**2.0), TRACEff(2,1), FNSCff(2,1)
WRITE (PWEE1,*) sigmaXXff(1,2)  ,sigmaYYff(1,2),  sigmaXYff(1,2), concFF(1,2), SQRT(UxFOUR(1,2)**2.0+UyFOUR(1,2)**2.0),	&
		SQRT(ORIXff(1,2)**2.0+ORIYff(1,2)**2.0), TRACEff(1,2), FNSCff(1,2)
WRITE (PWEE1,*) sigmaXXff(2,2)  ,sigmaYYff(2,2),  sigmaXYff(2,2), concFF(2,2), SQRT(UxFOUR(2,2)**2.0+UyFOUR(2,2)**2.0),	&
		SQRT(ORIXff(2,2)**2.0+ORIYff(2,2)**2.0), TRACEff(2,2), FNSCff(2,2)
WRITE (PWEE1,*) sigmaXXff(3,1)  ,sigmaYYff(3,1),  sigmaXYff(3,1), concFF(3,1), SQRT(UxFOUR(3,1)**2.0+UyFOUR(3,1)**2.0),	&
		SQRT(ORIXff(3,1)**2.0+ORIYff(3,1)**2.0), TRACEff(3,1), FNSCff(3,1)
WRITE (PWEE1,*) sigmaXXff(1,3)  ,sigmaYYff(1,3),  sigmaXYff(1,3), concFF(1,3), SQRT(UxFOUR(1,3)**2.0+UyFOUR(1,3)**2.0),	&
		SQRT(ORIXff(1,3)**2.0+ORIYff(1,3)**2.0), TRACEff(1,3), FNSCff(1,3)
WRITE (PWEE1,*) sigmaXXff(3,2)  ,sigmaYYff(3,2),  sigmaXYff(3,2), concFF(3,2), SQRT(UxFOUR(3,2)**2.0+UyFOUR(3,2)**2.0),	&
		SQRT(ORIXff(3,2)**2.0+ORIYff(3,2)**2.0), TRACEff(3,2), FNSCff(3,2)
WRITE (PWEE1,*) sigmaXXff(2,3)  ,sigmaYYff(2,3),  sigmaXYff(2,3), concFF(2,3), SQRT(UxFOUR(2,3)**2.0+UyFOUR(2,3)**2.0),	&
		SQRT(ORIXff(2,3)**2.0+ORIYff(2,3)**2.0), TRACEff(2,3), FNSCff(2,3)
WRITE (PWEE1,*) sigmaXXff(3,3)  ,sigmaYYff(3,3),  sigmaXYff(3,3), concFF(3,3), SQRT(UxFOUR(3,3)**2.0+UyFOUR(3,3)**2.0), &
		SQRT(ORIXff(3,3)**2.0+ORIYff(3,3)**2.0), TRACEff(3,3), FNSCff(3,3)
WRITE (PWEE1,*) sigmaXXff(4,1)  ,sigmaYYff(4,1),  sigmaXYff(4,1), concFF(4,1), SQRT(UxFOUR(4,1)**2.0+UyFOUR(4,1)**2.0), &
		SQRT(ORIXff(4,1)**2.0+ORIYff(4,1)**2.0), TRACEff(4,1), FNSCff(4,1)
WRITE (PWEE1,*) sigmaXXff(1,4)  ,sigmaYYff(1,4),  sigmaXYff(1,4), concFF(1,4), SQRT(UxFOUR(1,4)**2.0+UyFOUR(1,4)**2.0), &
		SQRT(ORIXff(1,4)**2.0+ORIYff(1,4)**2.0), TRACEff(1,4), FNSCff(1,4)
WRITE (PWEE1,*) sigmaXXff(4,2)  ,sigmaYYff(4,2),  sigmaXYff(4,2), concFF(4,2), SQRT(UxFOUR(4,2)**2.0+UyFOUR(4,2)**2.0), &
		SQRT(ORIXff(4,2)**2.0+ORIYff(4,2)**2.0), TRACEff(4,2), FNSCff(4,2)
WRITE (PWEE1,*) sigmaXXff(2,4)  ,sigmaYYff(2,4),  sigmaXYff(2,4), concFF(2,4), SQRT(UxFOUR(2,4)**2.0+UyFOUR(2,4)**2.0), &
		SQRT(ORIXff(2,4)**2.0+ORIYff(2,4)**2.0), TRACEff(2,4), FNSCff(2,4)
WRITE (PWEE1,*) sigmaXXff(3,4)  ,sigmaYYff(3,4),  sigmaXYff(3,4), concFF(3,4), SQRT(UxFOUR(3,4)**2.0+UyFOUR(3,4)**2.0), &
		SQRT(ORIXff(3,4)**2.0+ORIYff(3,4)**2.0), TRACEff(3,4), FNSCff(3,4)
WRITE (PWEE1,*) sigmaXXff(4,4)  ,sigmaYYff(4,4),  sigmaXYff(4,4), concFF(4,4), SQRT(UxFOUR(4,4)**2.0+UyFOUR(4,4)**2.0), &
		SQRT(ORIXff(4,4)**2.0+ORIYff(4,4)**2.0), TRACEff(4,4), FNSCff(4,4)
WRITE (PWEE1,*) sigmaXXff(5,1)  ,sigmaYYff(5,1),  sigmaXYff(5,1), concFF(5,1), SQRT(UxFOUR(5,1)**2.0+UyFOUR(5,1)**2.0), &
		SQRT(ORIXff(5,1)**2.0+ORIYff(5,1)**2.0), TRACEff(5,1), FNSCff(5,1)
WRITE (PWEE1,*) sigmaXXff(5,3)  ,sigmaYYff(5,3),  sigmaXYff(5,3), concFF(5,3), SQRT(UxFOUR(5,3)**2.0+UyFOUR(5,3)**2.0), &
		SQRT(ORIXff(5,3)**2.0+ORIYff(5,3)**2.0), TRACEff(5,3), FNSCff(5,3)
WRITE (PWEE1,*) sigmaXXff(1,6)  ,sigmaYYff(1,6),  sigmaXYff(1,6), concFF(1,6), SQRT(UxFOUR(1,6)**2.0+UyFOUR(1,6)**2.0), &
		SQRT(ORIXff(1,6)**2.0+ORIYff(1,6)**2.0), TRACEff(1,6), FNSCff(1,6)
WRITE (PWEE1,*) sigmaXXff(7,1)  ,sigmaYYff(7,1),  sigmaXYff(7,1), concFF(7,1), SQRT(UxFOUR(7,1)**2.0+UyFOUR(7,1)**2.0), &
		SQRT(ORIXff(7,1)**2.0+ORIYff(7,1)**2.0), TRACEff(7,1), FNSCff(7,1)
WRITE (PWEE1,*) sigmaXXff(1,8)  ,sigmaYYff(1,8),  sigmaXYff(1,8), concFF(1,8), SQRT(UxFOUR(1,8)**2.0+UyFOUR(1,8)**2.0), &
		SQRT(ORIXff(1,8)**2.0+ORIYff(1,8)**2.0), TRACEff(1,8), FNSCff(1,8)
WRITE (PWEE1,*) sigmaXXff(5,5)  ,sigmaYYff(5,5),  sigmaXYff(5,5), concFF(5,5), SQRT(UxFOUR(5,5)**2.0+UyFOUR(5,5)**2.0), &
		SQRT(ORIXff(5,5)**2.0+ORIYff(5,5)**2.0), TRACEff(5,5), FNSCff(5,5)
WRITE (PWEE1,*) sigmaXXff(6,6)  ,sigmaYYff(6,6),  sigmaXYff(6,6), concFF(6,6), SQRT(UxFOUR(6,6)**2.0+UyFOUR(6,6)**2.0), &
		SQRT(ORIXff(6,6)**2.0+ORIYff(6,6)**2.0), TRACEff(6,6), FNSCff(6,6)
WRITE (PWEE1,*) sigmaXXff(10,5)  ,sigmaYYff(10,5),  sigmaXYff(10,5), concFF(10,5), 					&
		SQRT(UxFOUR(10,5)**2.0+UyFOUR(10,5)**2.0), SQRT(ORIXff(10,5)**2.0+ORIYff(10,5)**2.0), TRACEff(10,5), FNSCff(10,5)
WRITE (PWEE1,*) sigmaXXff(13,13)  ,sigmaYYff(13,13),  sigmaXYff(13,13), concFF(13,13),					&
		SQRT(UxFOUR(13,13)**2.0+UyFOUR(13,13)**2.0), SQRT(ORIXff(13,13)**2.0+ORIYff(13,13)**2.0), TRACEff(13,13), FNSCff(13,13)


NUME=0
 DO I=1,sqrtTRACER
 	 DO J=1,sqrtTRACER
 		NUME=NUME+1						! index for tracer
 		Xtracer(NUME,1)=(I-.50001)*(SIZEB-DELTAx)/sqrtTRACER	! initialization of the x-location of tracers
		Ytracer(NUME,1)=(J-.50001)*(SIZEB-DELTAy)/sqrtTRACER	! initialization of the y-location of tracers

 		posiX=Xtracer(NUME,1)/  DELTAX				! the ranges of the index corresponds to their x-location
 		posiY=Ytracer(NUME,1)/  DELTAY				! the ranges of the index corresponds to their y-location

 		Xbelow=floor(posiX)+1
 		Ybelow=floor(posiY)+1
 		Xabove=ceiling(posiX)+1
 		Yabove=ceiling(posiY)+1

		!!!! bilinear interpolation
 		UXtracer(NUME,1)=(Xabove-posiX)*(Yabove-posiY)*UXF(Xbelow,Ybelow)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
 				+(posiX-Xbelow)*(Yabove-posiY)*UXF(Xabove,Ybelow)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
 				+(Xabove-posiX)*(posiY-Ybelow)*UXF(Xbelow,Yabove)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
 				+(posiX-Xbelow)*(posiY-Ybelow)*UXF(Xabove,Yabove)/((Xabove-Xbelow)*(Yabove-Ybelow))

 		UYtracer(NUME,1)=(Xabove-posiX)*(Yabove-posiY)*UYF(Xbelow,Ybelow)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
 				+(posiX-Xbelow)*(Yabove-posiY)*UYF(Xabove,Ybelow)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
 				+(Xabove-posiX)*(posiY-Ybelow)*UYF(Xbelow,Yabove)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
 				+(posiX-Xbelow)*(posiY-Ybelow)*UYF(Xabove,Yabove)/((Xabove-Xbelow)*(Yabove-Ybelow))

 		!WRITE (kjk3,*)	0, nume, Xtracer(NUME,1), Ytracer(NUME,1), UXtracer(NUME,1), UYtracer(NUME,1)


 		!XtracerNOp(NUME)=Xtracer0(NUME)+UXtracer(NUME,1)*DELTAt
 		!YtracerNOp(NUME)=Ytracer0(NUME)+UYtracer(NUME,1)*DELTAt

		XtracerNOp(NUME,1)=Xtracer(NUME,1)
		YtracerNOp(NUME,1)=Ytracer(NUME,1)

 		!if (Xtracer(NUME,1)>sizeb) then
 		!	dx=Xtracer(NUME,1)-sizeb
		!	Xtracer(NUME,1)=dx
 		!end if
 		!if (Ytracer(NUME,1)>sizeb) then
 		!	dy=Ytracer(NUME,1)-sizeb
 		!	Ytracer(NUME,1)=dy
 		!end if
		!if (Xtracer(NUME,1)<0) then
 		!	dx=Xtracer(NUME,1)
 		!	Xtracer(NUME,1)=sizeb+dx
 		!end if
 		!if (Ytracer(NUME,1)<0) then
 		!	dy=Ytracer(NUME,1)
 		!	Ytracer(NUME,1)=sizeb+dy
 		!end if

 		!Xtracer0(NUME)=Xtracer(NUME,1)					! X position of tracers at t=1 will be used at t=2
		!Ytracer0(NUME)=Ytracer(NUME,1)					! Y position of tracers at t=1 will be used at t=2

 	 END DO
 END DO

  OPEN(PW8	,FILE='D007.txt')
  open(GHG	,FILE='conce.txt')
  OPEN(FGF	,FILE='uDOTn.txt')
  OPEN(FGF1	,FILE='POWER.txt')
  OPEN(FGF2	,FILE='JAHAT.txt')
  OPEN(FGF3	,FILE='JAHAT1.txt')
  OPEN(FGF4	,FILE='JAHAT2.txt')
  OPEN(FGF5	,FILE='JAHAT3.txt')
  OPEN(PSP20	,FILE='CONTINUITY.txt')
  OPEN(W7	,FILE='LNvsNLtaup.txt')
  OPEN(W99	,FILE='STRESS.txt')

 CALL CPU_TIME ( time_begin )
 DO TIME=1, stepN
	NUMER	=0
	ERRORM	=1.0
	PSIi	=PSI
	TAUpXXi	=TAUpXX
	TAUpYYi	=TAUpYY
	TAUpXYi	=TAUpXY
	DO WHILE ( (ERRORM .GT. TOLERANCE ) .and. (NUMER < ITERATION) )
		NUMER	=NUMER+1
		IF (NUMER .EQ. 1) THEN
			DO I=1, GRIDX
				DO J=1, GRIDY
					PPYY(I,J) =-TAUpYYi(I,J)/(4.0*De) - delUYY(I,J)*ZETA/(2.0*De) &
					- UYF(I,J)*delTAUpYYy(I,J) -UXF(I,J)*delTAUpYYx(I,J) + 2.0*TAUpYYi(I,J)*delUYY(I,J) &
					+ 2.0*TAUpXYi(I,J)*delUXY(I,J)

					PPXX(I,J) =-TAUpXXi(I,J)/(4.0*De) - delUXX(I,J)*ZETA/(2.0*De) &
					- UXF(I,J)*delTAUpXXx(I,J)- UYF(I,J)*delTAUpXXy(I,J) + 2.0*TAUpXXi(I,J)*delUXX(I,J) &
					+ 2.0*TAUpXYi(I,J)*delUYX(I,J)

					PPXY(I,J) =-TAUpXYi(I,J)/(4.0*De) - ( delUXY(I,J)+delUYX(I,J) )*ZETA/(4.0*De) &
					-delTAUpXYx(I,J)*UXF(I,J) -delTAUpXYy(I,J)*UYF(I,J) +TAUpXXi(I,J)*delUXY(I,J) &
					+ delUYX(I,J)*TAUpYYi(I,J)

					TAUpXX1(I,J)=DELTAt*PPXX(I,J) + TAUpXXi(I,J)
					TAUpXX(I,J) =TAUpXX1(I,J)

					TAUpYY1(I,J)=DELTAt*PPYY(I,J) + TAUpYYi(I,J)
					TAUpYY(I,J) =TAUpYY1(I,J)

					TAUpXY1(I,J)=DELTAt*PPXY(I,J) + TAUpXYi(I,J)
					TAUpXY(I,J) =TAUpXY1(I,J)

					DO A=1, PHID

						PP(I,J,A) = 2.0*GAMM*( PSIi(I,J,A) ) * (     delUXX(I,J) * COS(4.0*PI*(A-1)/(PHID-1)) &
						+ 0.5*( delUXY(I,J)+delUYX(I,J) ) * SIN(4.0*PI*(A-1)/(PHID-1))   )&
						-( COS(2.0*PI*(A-1)/(PHID-1)) )*delXPSI(I,J,A) - ( SIN(2.0*PI*(A-1)/(PHID-1)) ) * delYPSI(I,J,A) &
			   			+ Dt* delY2PSI(I,J,A) + Dt* delX2PSI(I,J,A) + Dr * delO2PSI(I,J,A) &
					   	- UXF(I,J)*delXPSI(I,J,A) - UYF(I,J)*delYPSI(I,J,A) &
						+ ( (SIN(2.0*PI*(A-1)/(PHID-1))) **2.0 ) * (delOPSI(I,J,A)*delUXY(I,J)) &
						- ( (COS(2.0*PI*(A-1)/(PHID-1))) **2.0 ) * (delOPSI(I,J,A)*delUYX(I,J)) &
						+ SIN(4.0*PI*(A-1)/(PHID-1))*delOPSI(I,J,A)*delUXX(I,J)

						PSI1(I,J,A)= DELTAt*PP(I,J,A) + PSIi(I,J,A)

						PSI(I,J,A) = PSI1(I,J,A)

					END DO
				END DO
			END DO
			CALL CPU_TIME ( time_end1 )
			elapsed1 =elapsed1 + time_end1 - time_begin1
		END IF
		CALL CPU_TIME ( time_begin1 )
		CALL   estimation(ALPHA, GRIDX, GRIDY, DELTAX, DELTAY, PHID, PSI, SIGMAXX, SIGMAYY, SIGMAXY, &
	    	       CONC, ORIX, ORIY, delXPSI, delX2PSI, delYPSI, delY2PSI, delO2PSI, delOPSI, delSIGMAX, &
		       delSIGMAY, DELdotNCX, DELdotNCY, ENTROPY )
		CALL CPU_TIME ( time_end1 )
		elapsed1 =elapsed1 + time_end1 - time_begin1

		CALL CPU_TIME ( time_begin2 )
		call	ggm(MODES,SIZEB, DELTAX, DELTAY, GRIDX, GRIDY, UXF, UYF, UxFOUR, UyFOUR, SIGMAXY, SIGMAXX, &
			SIGMAYY, TAUpXY, TAUpXX, TAUpYY, contild, elongationXff, shearXff, UelongationX, UshearX, &
			sigmaXYm, sigmaXXm, TAUpxyM, TAUpxxM, ZETA)

		CALL CPU_TIME ( time_end2 )
		elapsed2 =elapsed2+ time_end2 - time_begin2

		CALL CPU_TIME ( time_begin1 )
		CALL 	velocityDerivative(UXF, UYF, DELTAX, DELTAY, GRIDX, GRIDY, delUXX, delUYY, delUYX, delUXY)
		CALL 	tauPderivative(ZETA,De,UXF,UYF,GRIDX,GRIDY,DELTAX,DELTAY,delUXX, delUYY, delUYX, delUXY, &
			TAUpXX,TAUpYY,TAUpXY, delTAUpXXx, delTAUpXYx, delTAUpYYx, delTAUpXXy, delTAUpXYy, delTAUpYYy)
		!addCON1=SUM(CONC)
		!CONC=(CONC - (addCON1-gridx*gridy)/(gridx*gridy) )
		CALL CPU_TIME ( time_end1 )
		elapsed1 =elapsed1 + time_end1 - time_begin1

		CALL CPU_TIME ( time_begin1 )
	 	DO I=1, GRIDX
			DO J=1, GRIDY
				kom(i,j)=0.0
				LToldrXX(I,J)=  -TAUpXX(I,J)/(4.0*De) - delUXX(I,J)*ZETA/(2.0*De)

				convecXX(I,J)=  -UXF(I,J)*delTAUpXXx(I,J)-UYF(I,J)*delTAUpXXy(I,J) + &
						2.0*TAUpXX(I,J)*delUXX(I,J)+ 2.0*TAUpXY(I,J)*delUYX(I,J)

				LToldrXY(I,J)=  -TAUpXY(I,J)/(4.0*De) - ( delUXY(I,J)+delUYX(I,J) )*ZETA/(4.0*De)

				convecXY(I,J)=  -delTAUpXYx(I,J)*UXF(I,J)-delTAUpXYy(I,J)*UYF(I,J)  + &
					        TAUpXX(I,J)*delUXY(I,J) + delUYX(I,J)*TAUpYY(I,J)

				orientXX(I,J)=2.0*TAUpXX(I,J)*delUXX(I,J)+ 2.0*TAUpXY(I,J)*delUYX(I,J)
				orientXY(I,J)=TAUpXX(I,J)*delUXY(I,J) + delUYX(I,J)*TAUpYY(I,J)

				PPYY(I,J) =-TAUpYY(I,J)/(4.0*De) - delUYY(I,J)*ZETA/(2.0*De) 					&
				- UYF(I,J)*delTAUpYYy(I,J)- UXF(I,J)*delTAUpYYx(I,J) + 2.0*TAUpYY(I,J)*delUYY(I,J) 		&
				+ 2.0*TAUpXY(I,J)*delUXY(I,J)

				PPXX(I,J) =-TAUpXX(I,J)/(4.0*De) - delUXX(I,J)*ZETA/(2.0*De) 					&
				- UXF(I,J)*delTAUpXXx(I,J)- UYF(I,J)*delTAUpXXy(I,J) + 2.0*TAUpXX(I,J)*delUXX(I,J) 		&
				+ 2.0*TAUpXY(I,J)*delUYX(I,J)

				PPXY(I,J) =-TAUpXY(I,J)/(4.0*De) - ( delUXY(I,J)+delUYX(I,J) )*ZETA/(4.0*De) 			&
				-delTAUpXYx(I,J)*UXF(I,J) -delTAUpXYy(I,J)*UYF(I,J) +TAUpXX(I,J)*delUXY(I,J) 			&
				+ delUYX(I,J)*TAUpYY(I,J)

				TAUpXX1(I,J)=DELTAt*PPXX(I,J) + TAUpXXi(I,J)
				TAUpYY1(I,J)=DELTAt*PPYY(I,J) + TAUpYYi(I,J)
				TAUpXY1(I,J)=DELTAt*PPXY(I,J) + TAUpXYi(I,J)

				ERRORtauP(I,J) = 1.00000 * (  ABS( (TAUpXX1(I,J)-TAUpXX(I,J))/ TAUpXX1(I,J) ) + &
					  	     	      ABS( (TAUpYY1(I,J)-TAUpYY(I,J))/ TAUpYY1(I,J) ) + &
					       	              ABS( (TAUpXY1(I,J)-TAUpXY(I,J))/ TAUpXY1(I,J) ) )

				TAUpXX(I,J)=TAUpXX1(I,J)
				TAUpYY(I,J)=TAUpYY1(I,J)
				TAUpXY(I,J)=TAUpXY1(I,J)

				DO A=1, PHID

					PPP(I,J,A) = 2.0*GAMM*( PSI(I,J,A) ) * (     delUXX(I,J) * COS(4.0*PI*(A-1)/(PHID-1)) &
					+ 0.5*( delUXY(I,J)+delUYX(I,J) ) * SIN(4.0*PI*(A-1)/(PHID-1))   )&
					-( COS(2.0*PI*(A-1)/(PHID-1)) )*delXPSI(I,J,A) - ( SIN(2.0*PI*(A-1)/(PHID-1)) ) * delYPSI(I,J,A) &
		   			+ Dt* delY2PSI(I,J,A) + Dt* delX2PSI(I,J,A) + Dr * delO2PSI(I,J,A) &
				   	- UXF(I,J)*delXPSI(I,J,A) - UYF(I,J)*delYPSI(I,J,A) &
					+ ( SIN(2.0*PI*(A-1)/(PHID-1)) **2.0 ) * (delOPSI(I,J,A)*delUXY(I,J)) &
					- ( COS(2.0*PI*(A-1)/(PHID-1)) **2.0 ) * (delOPSI(I,J,A)*delUYX(I,J)) &
					+ SIN(4.0*PI*(A-1)/(PHID-1)) * delOPSI(I,J,A) * delUXX(I,J)


					PSI1(I,J,A)=DELTAt*PPP(I,J,A) + PSIi(I,J,A)

					ERROR(I,J,A) = ABS( (PSI1(I,J,A)-PSI(I,J,A))/( PSI(I,J,A)-1.0/TWOPI ) )

					PSI(I,J,A) = PSI1(I,J,A)
					IF (A .LT. PHID) THEN
						DA=(2.0*PI)/(PHID-1)
						kom(I,J) = kom(I,J) + 										&
							   0.5*ALPHA*( delUXX(I,J)*COS(4.0*PI*(A-1)/(PHID-1)) +			 		&
					   		   0.5*(delUXY(I,J)+delUYX(I,J))*SIN(4.0*PI*(A-1)/(PHID-1)) )* PSI(I,J,A) *DA + 	&
					 	  	   0.5*ALPHA*( delUXX(I,J)*COS(4.0*PI*A/(PHID-1)) + 					&
							   0.5*(delUXY(I,J)+delUYX(I,J))*SIN(4.0*PI*A/(PHID-1)) ) *PSI(I,J,A+1) * DA
					END IF
				END DO
				CONTINUITY(I,J) = ( delUXX(I,J)+delUYY(I,J) ) /( (UXF(I,J)**2.0 + UYF(I,J)**2.0)**0.5 )
			END DO
		END DO
		ERRORM1=maxval(error)
		ERRORM2=maxval(ERRORtauP)
		IF (ERRORM1>ERRORM2) THEN
			ERRORM=ERRORM1
		END IF
		IF (ERRORM2>ERRORM1) THEN
			ERRORM=ERRORM2
		END IF
		CALL CPU_TIME ( time_end1 )
		elapsed1 =elapsed1 + time_end1 - time_begin1
	END DO
	!addPSI3=SUM(TAUpXY)/(gridx*gridy)
	!TAUpXY=TAUpXY - addPSI3
	!addPSI3=SUM(TAUpXX)/(gridx*gridy)
	!TAUpXX=TAUpXX - addPSI3
	!addPSI3=SUM(TAUpYY)/(gridx*gridy)
	!TAUpYY=TAUpYY - addPSI3
	NUME=0
	DO I=1,sqrtTRACER
		DO J=1,sqrtTRACER
			NUME=NUME+1

			Xtracer(NUME,TIME+1)=Xtracer(NUME,TIME)+UXtracer(NUME,TIME)*DELTAt
			Ytracer(NUME,TIME+1)=Ytracer(NUME,TIME)+UYtracer(NUME,TIME)*DELTAt
			XtracerNOp(NUME,TIME+1)=XtracerNOp(NUME,TIME)+UXtracer(NUME,TIME)*DELTAt
			YtracerNOp(NUME,TIME+1)=YtracerNOp(NUME,TIME)+UYtracer(NUME,TIME)*DELTAt

			if (Xtracer(NUME,TIME+1)>sizeb) then
				dx=Xtracer(NUME,TIME+1)-sizeb
				Xtracer(NUME,TIME+1)=dx
			end if
			if (Ytracer(NUME,TIME+1)>sizeb) then
				dy=Ytracer(NUME,TIME+1)-sizeb
				Ytracer(NUME,TIME+1)=dy
			end if
			if (Xtracer(NUME,TIME+1)<0) then
				dx=Xtracer(NUME,TIME+1)
				Xtracer(NUME,TIME+1)=sizeb+dx
			end if
			if (Ytracer(NUME,TIME+1)<0) then
				dy=Ytracer(NUME,TIME+1)
				Ytracer(NUME,TIME+1)=sizeb+dy
			end if

			posiX=Xtracer(NUME,TIME+1)/DELTAX	!!!!the ranges of the index
			posiY=Ytracer(NUME,TIME+1)/DELTAY	!!!!the ranges of the index

			Xbelow=floor(posiX)+1
			Ybelow=floor(posiY)+1
			Xabove=ceiling(posiX)+1
			Yabove=ceiling(posiY)+1

			XXbelow=Xbelow
			YYbelow=Ybelow
			XXabove=Xabove
			YYabove=Yabove

			if (XXabove > gridx) then
				XXabove=1
			end if
			if (YYabove > gridy) then
				YYabove=1
			end if

			!!!! bilinear interpolation
			UXtracer(NUME,time+1)=(Xabove-posiX)*(Yabove-posiY)*UXF(XXbelow,YYbelow)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
					     +(posiX-Xbelow)*(Yabove-posiY)*UXF(XXabove,YYbelow)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
					     +(Xabove-posiX)*(posiY-Ybelow)*UXF(XXbelow,YYabove)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
					     +(posiX-Xbelow)*(posiY-Ybelow)*UXF(XXabove,YYabove)/((Xabove-Xbelow)*(Yabove-Ybelow))

			UYtracer(NUME,time+1)=(Xabove-posiX)*(Yabove-posiY)*UYF(XXbelow,YYbelow)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
					     +(posiX-Xbelow)*(Yabove-posiY)*UYF(XXabove,YYbelow)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
					     +(Xabove-posiX)*(posiY-Ybelow)*UYF(XXbelow,YYabove)/((Xabove-Xbelow)*(Yabove-Ybelow))	&
					     +(posiX-Xbelow)*(posiY-Ybelow)*UYF(XXabove,YYabove)/((Xabove-Xbelow)*(Yabove-Ybelow))

			!!!WRITE (kjk3,*)	nume, Xtracer(NUME,TIME), Ytracer(NUME,TIME), UXtracer(NUME,TIME), UYtracer(NUME,TIME)



			!if (time > stepN/2) then
			!	auxMSD=auxMSD &
			!	+( (XtracerNOp(NUME,TIME+1)-XtracerNOp(NUME,TIME))**2.0+(YtracerNOp(NUME,TIME+1)-YtracerNOp(NUME,TIME))**2.0)
			!end if

		END DO
	END DO
	CALL CPU_TIME ( time_begin3 )
	CALL fftwCheck(conc,sigmaxy,TAUpXY,SIZEB, DELTAX, DELTAY, GRIDX, GRIDY, concFF &
	     	,sigmaxx,TAUpXX, sigmayy, TAUpYY, stressXXff, stressYYff, stressXYff, mixing &
	     	,sigmaxxFF, sigmaxyFF, sigmayyFF, TAUpXYff, TAUpYYff, TAUpXXff, ORIX, ORIY, ORIXff &
	     	,ORIYff, TRACEff, FNSCff, LToldrXX, LToldrXY, LToldrXXff, LToldrXYff  &
		,CONVECxx, CONVECxy, CONVECxxFF, CONVECxyFF, orientXX, orientXY, orientXXff, orientXYff)


	IF (TIME .EQ. 1) THEN
		DO A=1, PHID
			WRITE (FGF2,*) PSI(GRIDX/4, GRIDY/2, A)
 		END DO
	END IF
	IF (TIME .EQ. stepN/3) THEN
		DO A=1, PHID
			WRITE (FGF3,*) PSI(GRIDX/4, GRIDY/2, A)
 		END DO
	END IF
	IF (TIME .EQ. 2*stepN/3) THEN
		DO A=1, PHID
			WRITE (FGF4,*) PSI(GRIDX/4, GRIDY/2, A)
 		END DO
	END IF
	IF (TIME .EQ. stepN) THEN
		DO A=1, PHID
			WRITE (FGF5,*) PSI(GRIDX/4, GRIDY/2, A)
 		END DO
	END IF

	AUXI     = MOD(time,int(1/DELTAt)+1)

	IF ( AUXI .EQ. 0 ) THEN
		time1=time1+1
		uDOTn=0.0
		CONAV=0.0
		POWER=0.0
		WAT  =0.0
		VEL  =0.0
		nxyN =0.0
		delUv=0.0
		delUxyV=0.0
		delUdif=0.0
		uSHEARdotELONG=0.0
		DO I=1, GRIDX
			DO J=1, GRIDY
				uDOTn=uDOTn + ( CONC(I,J)*UXF(I,J)*ORIX(I,J)+CONC(I,J)*UYF(I,J)*ORIY(I,J) ) *DELTAX*DELTAY
				VEL=VEL     + SQRT(UXF(I,J)**2.0+UYF(I,J)**2.0)/(GRIDX*GRIDY)
				CONAV=CONAV + CONC(I,J)*DELTAX*DELTAY
				nxyN =nxyN  + SQRT( ORIX(I,J)*ORIX(I,J)+ORIY(I,J)*ORIY(I,J) )/(GRIDX*GRIDY)
				WAT = WAT   - KOM(I,J)*DELTAX*DELTAY
				POWER=POWER + 2.0 * ( 2.0*delUXX(I,J)**2.0 + 0.5*(delUXY(I,J)+delUYX(I,J))**2.0 )  *DELTAX*DELTAY
				delUv=delUv     + SQRT( delUXX(I,J)*delUXX(I,J)+delUYY(I,J)*delUYY(I,J) )/(GRIDX*GRIDY)
				delUxyV=delUxyV + SQRT( delUXY(I,J)*delUXY(I,J)+delUYX(I,J)*delUYX(I,J) )/(GRIDX*GRIDY)
				delUdif=delUdif + abs ( delUXX(I,J)-delUYY(I,J) )/(GRIDX*GRIDY)
				uSHEARdotELONG=uSHEARdotELONG+(UelongationX(I,J)*UshearX(I,J))/abs( UelongationX(I,J)*UshearX(I,J) )

			END DO
		END DO
		write (FGF,*) uDOTn/CONAV
		write (FGF1,*) WAT, power
		WRITE (PSP7,*) ENTROPY
		WRITE (PSP20,*)  MAXVAL(ABS(contild)), MAXVAL(ABS(CONTINUITY))
		write (sty, *) nxyN, vel, delUxyV, delUv, delUdif
		write (kjk7,*) uSHEARdotELONG

		NUME=0
		DO I=1,sqrtTRACER
			DO J=1,sqrtTRACER
				NUME=NUME+1
				WRITE (kjk3,*)	TIME1, nume, Xtracer(NUME,TIME), Ytracer(NUME,TIME), UXtracer(NUME,TIME), &
					 		UYtracer(NUME,TIME)
			END DO
		END DO
		DO I=1, GRIDX
			 DO J=1, GRIDY
				WRITE (PSP,	*)  	CONC(I,J)
				WRITE (PSP1,	*) 	UXF(I,J)
				WRITE (PSP2,	*) 	UYF(I,J)
				WRITE (PSP3,	*) 	ORIX(I,J)
				WRITE (PSP4,	*) 	ORIY(I,J)
				WRITE (PSP5,	*) 	SIGMAXY(I,J)
				WRITE (PSP11,	*)	SIGMAXX(I,J)
				WRITE (PSP6,    *)	TAUpXX(I,J)+TAUpYY(I,J)
				WRITE (PSP8,	*)	TAUpXY(I,J)
				WRITE (kjk5,	*)	UshearX(I,J)
				WRITE (kjk6,	*)	UelongationX(I,J)
				WRITE (PSP9,    *)	TAUpXX(I,J)
				WRITE (PSP10,	*)	TAUpYY(I,J)
				WRITE (PSP12,	*)	sigmaXYm(I,J)
				WRITE (PSP13,	*)	sigmaXXm(I,J)
				WRITE (PSP14,	*)	TAUpxyM(I,J)
				WRITE (PSP15,	*)	TAUpxxM(I,J)
			END DO
		END DO

		WRITE (KJK1,*) stressXXff(1,1)  ,stressYYff(1,1),  stressXYff(1,1)
		WRITE (KJK1,*) stressXXff(2,1)  ,stressYYff(2,1),  stressXYff(2,1)
		WRITE (KJK1,*) stressXXff(1,2)  ,stressYYff(1,2),  stressXYff(1,2)
		WRITE (KJK1,*) stressXXff(2,2)  ,stressYYff(2,2),  stressXYff(2,2)
		WRITE (KJK1,*) stressXXff(3,1)  ,stressYYff(3,1),  stressXYff(3,1)
		WRITE (KJK1,*) stressXXff(1,3)  ,stressYYff(1,3),  stressXYff(1,3)
		WRITE (KJK1,*) stressXXff(3,2)  ,stressYYff(3,2),  stressXYff(3,2)
		WRITE (KJK1,*) stressXXff(2,3)  ,stressYYff(2,3),  stressXYff(2,3)
		WRITE (KJK1,*) stressXXff(3,3)  ,stressYYff(3,3),  stressXYff(3,3)
		WRITE (KJK1,*) stressXXff(4,1)  ,stressYYff(4,1),  stressXYff(4,1)
		WRITE (KJK1,*) stressXXff(1,4)  ,stressYYff(1,4),  stressXYff(1,4)
		WRITE (KJK1,*) stressXXff(4,2)  ,stressYYff(4,2),  stressXYff(4,2)
		WRITE (KJK1,*) stressXXff(2,4)  ,stressYYff(2,4),  stressXYff(2,4)
		WRITE (KJK1,*) stressXXff(3,4)  ,stressYYff(3,4),  stressXYff(3,4)
		WRITE (KJK1,*) stressXXff(4,4)  ,stressYYff(4,4),  stressXYff(4,4)
		WRITE (KJK1,*) stressXXff(5,1)  ,stressYYff(5,1),  stressXYff(5,1)
		WRITE (KJK1,*) stressXXff(5,3)  ,stressYYff(5,3),  stressXYff(5,3)
		WRITE (KJK1,*) stressXXff(1,6)  ,stressYYff(1,6),  stressXYff(1,6)
		WRITE (KJK1,*) stressXXff(7,1)  ,stressYYff(7,1),  stressXYff(7,1)
		WRITE (KJK1,*) stressXXff(1,8)  ,stressYYff(1,8),  stressXYff(1,8)
		WRITE (KJK1,*) stressXXff(5,5)  ,stressYYff(5,5),  stressXYff(5,5)
		WRITE (KJK1,*) stressXXff(6,6)  ,stressYYff(6,6),  stressXYff(6,6)
		WRITE (KJK1,*) stressXXff(10,5)  ,stressYYff(10,5),  stressXYff(10,5)
		WRITE (KJK1,*) stressXXff(13,13)  ,stressYYff(13,13),  stressXYff(13,13)

		WRITE (STT,*) orientXXff(1,1)  ,  orientXYff(1,1), UxFOUR(1,1), elongationXff(1,1), shearXff(1,1)
		WRITE (STT,*) orientXXff(2,1)  ,  orientXYff(2,1), UxFOUR(2,1), elongationXff(2,1), shearXff(2,1)
		WRITE (STT,*) orientXXff(1,2)  ,  orientXYff(1,2), UxFOUR(1,2),elongationXff(1,2), shearXff(1,2)
		WRITE (STT,*) orientXXff(2,2)  ,  orientXYff(2,2), UxFOUR(2,2),elongationXff(2,2), shearXff(2,2)
		WRITE (STT,*) orientXXff(3,1)  ,  orientXYff(3,1), UxFOUR(3,1),elongationXff(3,1), shearXff(3,1)
		WRITE (STT,*) orientXXff(1,3)  ,  orientXYff(1,3), UxFOUR(1,3),elongationXff(1,3), shearXff(1,3)
		WRITE (STT,*) orientXXff(3,2)  ,  orientXYff(3,2), UxFOUR(3,2),elongationXff(3,2), shearXff(3,2)
		WRITE (STT,*) orientXXff(2,3)  ,  orientXYff(2,3), UxFOUR(2,3),elongationXff(2,3), shearXff(2,3)
		WRITE (STT,*) orientXXff(3,3)  ,  orientXYff(3,3), UxFOUR(3,3),elongationXff(3,3), shearXff(3,3)
		WRITE (STT,*) orientXXff(4,1)  ,  orientXYff(4,1), UxFOUR(4,1),elongationXff(4,1), shearXff(4,1)
		WRITE (STT,*) orientXXff(1,4)  ,  orientXYff(1,4), UxFOUR(1,4),elongationXff(1,4), shearXff(1,4)
		WRITE (STT,*) orientXXff(4,2)  ,  orientXYff(4,2), UxFOUR(4,2),elongationXff(4,2), shearXff(4,2)
		WRITE (STT,*) orientXXff(2,4)  ,  orientXYff(2,4), UxFOUR(2,4),elongationXff(2,4), shearXff(2,4)
		WRITE (STT,*) orientXXff(3,4)  ,  orientXYff(3,4), UxFOUR(3,4),elongationXff(3,4), shearXff(3,4)
		WRITE (STT,*) orientXXff(4,4)  ,  orientXYff(4,4), UxFOUR(4,4),elongationXff(4,4), shearXff(4,4)
		WRITE (STT,*) orientXXff(5,1)  ,  orientXYff(5,1), UxFOUR(5,1),elongationXff(5,1), shearXff(5,1)
		WRITE (STT,*) orientXXff(5,3)  ,  orientXYff(5,3), UxFOUR(5,3),elongationXff(5,3), shearXff(5,3)
		WRITE (STT,*) orientXXff(1,6)  ,  orientXYff(1,6), UxFOUR(1,6),elongationXff(1,6), shearXff(1,6)
		WRITE (STT,*) orientXXff(7,1)  ,  orientXYff(7,1), UxFOUR(7,1),elongationXff(7,1), shearXff(7,1)
		WRITE (STT,*) orientXXff(1,8)  ,  orientXYff(1,8), UxFOUR(1,8),elongationXff(1,8), shearXff(1,8)
		WRITE (STT,*) orientXXff(5,5)  ,  orientXYff(5,5), UxFOUR(5,5),elongationXff(5,5), shearXff(5,5)
		WRITE (STT,*) orientXXff(6,6)  ,  orientXYff(6,6), UxFOUR(6,6),elongationXff(6,6), shearXff(6,6)
		WRITE (STT,*) orientXXff(10,5) ,  orientXYff(10,5), UxFOUR(10,5),elongationXff(10,5), shearXff(10,5)
		WRITE (STT,*) orientXXff(13,13),  orientXYff(13,13), UxFOUR(13,13),elongationXff(13,13), shearXff(13,13)

		WRITE (PWEE2,*) TauPxxff(1,1),TauPyyff(1,1),TauPXYff(1,1),convecxxFF(1,1),convecxyFF(1,1),LToldrXXff(1,1),LToldrXYff(1,1)
		WRITE (PWEE2,*) TauPxxff(2,1),TauPyyff(2,1),TauPXYff(2,1),convecxxFF(2,1),convecxyFF(2,1),LToldrXXff(2,1),LToldrXYff(2,1)
		WRITE (PWEE2,*) TauPxxff(1,2),TauPyyff(1,2),TauPXYff(1,2),convecxxFF(1,2),convecxyFF(1,2),LToldrXXff(1,2),LToldrXYff(1,2)
		WRITE (PWEE2,*) TauPxxff(2,2),TauPyyff(2,2),TauPXYff(2,2),convecxxFF(2,2),convecxyFF(2,2),LToldrXXff(2,2),LToldrXYff(2,2)
		WRITE (PWEE2,*) TauPxxff(3,1),TauPyyff(3,1),TauPXYff(3,1),convecxxFF(3,1),convecxyFF(3,1),LToldrXXff(3,1),LToldrXYff(3,1)
		WRITE (PWEE2,*) TauPxxff(1,3),TauPyyff(1,3),TauPXYff(1,3),convecxxFF(1,3),convecxyFF(1,3),LToldrXXff(1,3),LToldrXYff(1,3)
		WRITE (PWEE2,*) TauPxxff(3,2),TauPyyff(3,2),TauPXYff(3,2),convecxxFF(3,2),convecxyFF(3,2),LToldrXXff(3,2),LToldrXYff(3,2)
		WRITE (PWEE2,*) TauPxxff(2,3),TauPyyff(2,3),TauPXYff(2,3),convecxxFF(2,3),convecxyFF(2,3),LToldrXXff(2,3),LToldrXYff(2,3)
		WRITE (PWEE2,*) TauPxxff(3,3),TauPyyff(3,3),TauPXYff(3,3),convecxxFF(3,3),convecxyFF(3,3),LToldrXXff(3,3),LToldrXYff(3,3)
		WRITE (PWEE2,*) TauPxxff(4,1),TauPyyff(4,1),TauPXYff(4,1),convecxxFF(4,1),convecxyFF(4,1),LToldrXXff(4,1),LToldrXYff(4,1)
		WRITE (PWEE2,*) TauPxxff(1,4),TauPyyff(1,4),TauPXYff(1,4),convecxxFF(1,4),convecxyFF(1,4),LToldrXXff(1,4),LToldrXYff(1,4)
		WRITE (PWEE2,*) TauPxxff(4,2),TauPyyff(4,2),TauPXYff(4,2),convecxxFF(4,2),convecxyFF(4,2),LToldrXXff(4,2),LToldrXYff(4,2)
		WRITE (PWEE2,*) TauPxxff(2,4),TauPyyff(2,4),TauPXYff(2,4),convecxxFF(2,4),convecxyFF(2,4),LToldrXXff(2,4),LToldrXYff(2,4)
		WRITE (PWEE2,*) TauPxxff(3,4),TauPyyff(3,4),TauPXYff(3,4),convecxxFF(3,4),convecxyFF(3,4),LToldrXXff(3,4),LToldrXYff(3,4)
		WRITE (PWEE2,*) TauPxxff(4,4),TauPyyff(4,4),TauPXYff(4,4),convecxxFF(4,4),convecxyFF(4,4),LToldrXXff(4,4),LToldrXYff(4,4)
		WRITE (PWEE2,*) TauPxxff(5,1),TauPyyff(5,1),TauPXYff(5,1),convecxxFF(4,1),convecxyFF(5,1),LToldrXXff(5,1),LToldrXYff(5,1)
		WRITE (PWEE2,*) TauPxxff(5,3),TauPyyff(5,3),TauPXYff(5,3),convecxxFF(5,3),convecxyFF(5,3),LToldrXXff(5,3),LToldrXYff(5,3)
		WRITE (PWEE2,*) TauPxxff(1,6),TauPyyff(1,6),TauPXYff(1,6),convecxxFF(1,6),convecxyFF(1,6),LToldrXXff(1,6),LToldrXYff(1,6)
		WRITE (PWEE2,*) TauPxxff(7,1),TauPyyff(7,1),TauPXYff(7,1),convecxxFF(7,1),convecxyFF(7,1),LToldrXXff(7,1),LToldrXYff(7,1)
		WRITE (PWEE2,*) TauPxxff(1,8),TauPyyff(1,8),TauPXYff(1,8),convecxxFF(1,8),convecxyFF(1,8),LToldrXXff(1,8),LToldrXYff(1,8)
		WRITE (PWEE2,*) TauPxxff(5,5),TauPyyff(5,5),TauPXYff(5,5),convecxxFF(5,5),convecxyFF(5,5),LToldrXXff(5,5),LToldrXYff(5,5)
		WRITE (PWEE2,*) TauPxxff(6,6),TauPyyff(6,6),TauPXYff(6,6),convecxxFF(6,6),convecxyFF(6,6),LToldrXXff(6,6),LToldrXYff(6,6)
		WRITE (PWEE2,*) TauPxxff(10,5),TauPyyff(10,5),TauPXYff(10,5),convecxxFF(10,5),convecxyFF(10,5), 			&
				LToldrXXff(10,5),LToldrXYff(10,5)
		WRITE (PWEE2,*) TauPxxff(13,13),TauPyyff(13,13),TauPXYff(13,13),convecxxFF(13,13),convecxyFF(13,13), 			&
				LToldrXXff(13,13),LToldrXYff(13,13)

		WRITE (PWEE1,*) sigmaXXff(1,1)  ,sigmaYYff(1,1),  sigmaXYff(1,1), concFF(1,1), SQRT(UxFOUR(1,1)**2.0+UyFOUR(1,1)**2.0), &
				SQRT(ORIXff(1,1)**2.0+ORIYff(1,1)**2.0), TRACEff(1,1), FNSCff(1,1)
		WRITE (PWEE1,*) sigmaXXff(2,1)  ,sigmaYYff(2,1),  sigmaXYff(2,1), concFF(2,1), SQRT(UxFOUR(2,1)**2.0+UyFOUR(2,1)**2.0), &
				SQRT(ORIXff(2,1)**2.0+ORIYff(2,1)**2.0), TRACEff(2,1), FNSCff(2,1)
		WRITE (PWEE1,*) sigmaXXff(1,2)  ,sigmaYYff(1,2),  sigmaXYff(1,2), concFF(1,2), SQRT(UxFOUR(1,2)**2.0+UyFOUR(1,2)**2.0),	&
				SQRT(ORIXff(1,2)**2.0+ORIYff(1,2)**2.0), TRACEff(1,2), FNSCff(1,2)
		WRITE (PWEE1,*) sigmaXXff(2,2)  ,sigmaYYff(2,2),  sigmaXYff(2,2), concFF(2,2), SQRT(UxFOUR(2,2)**2.0+UyFOUR(2,2)**2.0),	&
				SQRT(ORIXff(2,2)**2.0+ORIYff(2,2)**2.0), TRACEff(2,2), FNSCff(2,2)
		WRITE (PWEE1,*) sigmaXXff(3,1)  ,sigmaYYff(3,1),  sigmaXYff(3,1), concFF(3,1), SQRT(UxFOUR(3,1)**2.0+UyFOUR(3,1)**2.0),	&
				SQRT(ORIXff(3,1)**2.0+ORIYff(3,1)**2.0), TRACEff(3,1), FNSCff(3,1)
		WRITE (PWEE1,*) sigmaXXff(1,3)  ,sigmaYYff(1,3),  sigmaXYff(1,3), concFF(1,3), SQRT(UxFOUR(1,3)**2.0+UyFOUR(1,3)**2.0),	&
				SQRT(ORIXff(1,3)**2.0+ORIYff(1,3)**2.0), TRACEff(1,3), FNSCff(1,3)
		WRITE (PWEE1,*) sigmaXXff(3,2)  ,sigmaYYff(3,2),  sigmaXYff(3,2), concFF(3,2), SQRT(UxFOUR(3,2)**2.0+UyFOUR(3,2)**2.0),	&
				SQRT(ORIXff(3,2)**2.0+ORIYff(3,2)**2.0), TRACEff(3,2), FNSCff(3,2)
		WRITE (PWEE1,*) sigmaXXff(2,3)  ,sigmaYYff(2,3),  sigmaXYff(2,3), concFF(2,3), SQRT(UxFOUR(2,3)**2.0+UyFOUR(2,3)**2.0),	&
				SQRT(ORIXff(2,3)**2.0+ORIYff(2,3)**2.0), TRACEff(2,3), FNSCff(2,3)
		WRITE (PWEE1,*) sigmaXXff(3,3)  ,sigmaYYff(3,3),  sigmaXYff(3,3), concFF(3,3), SQRT(UxFOUR(3,3)**2.0+UyFOUR(3,3)**2.0), &
				SQRT(ORIXff(3,3)**2.0+ORIYff(3,3)**2.0), TRACEff(3,3), FNSCff(3,3)
		WRITE (PWEE1,*) sigmaXXff(4,1)  ,sigmaYYff(4,1),  sigmaXYff(4,1), concFF(4,1), SQRT(UxFOUR(4,1)**2.0+UyFOUR(4,1)**2.0), &
				SQRT(ORIXff(4,1)**2.0+ORIYff(4,1)**2.0), TRACEff(4,1), FNSCff(4,1)
		WRITE (PWEE1,*) sigmaXXff(1,4)  ,sigmaYYff(1,4),  sigmaXYff(1,4), concFF(1,4), SQRT(UxFOUR(1,4)**2.0+UyFOUR(1,4)**2.0), &
				SQRT(ORIXff(1,4)**2.0+ORIYff(1,4)**2.0), TRACEff(1,4), FNSCff(1,4)
		WRITE (PWEE1,*) sigmaXXff(4,2)  ,sigmaYYff(4,2),  sigmaXYff(4,2), concFF(4,2), SQRT(UxFOUR(4,2)**2.0+UyFOUR(4,2)**2.0), &
				SQRT(ORIXff(4,2)**2.0+ORIYff(4,2)**2.0), TRACEff(4,2), FNSCff(4,2)
		WRITE (PWEE1,*) sigmaXXff(2,4)  ,sigmaYYff(2,4),  sigmaXYff(2,4), concFF(2,4), SQRT(UxFOUR(2,4)**2.0+UyFOUR(2,4)**2.0), &
				SQRT(ORIXff(2,4)**2.0+ORIYff(2,4)**2.0), TRACEff(2,4), FNSCff(2,4)
		WRITE (PWEE1,*) sigmaXXff(3,4)  ,sigmaYYff(3,4),  sigmaXYff(3,4), concFF(3,4), SQRT(UxFOUR(3,4)**2.0+UyFOUR(3,4)**2.0), &
				SQRT(ORIXff(3,4)**2.0+ORIYff(3,4)**2.0), TRACEff(3,4), FNSCff(3,4)
		WRITE (PWEE1,*) sigmaXXff(4,4)  ,sigmaYYff(4,4),  sigmaXYff(4,4), concFF(4,4), SQRT(UxFOUR(4,4)**2.0+UyFOUR(4,4)**2.0), &
				SQRT(ORIXff(4,4)**2.0+ORIYff(4,4)**2.0), TRACEff(4,4), FNSCff(4,4)
		WRITE (PWEE1,*) sigmaXXff(5,1)  ,sigmaYYff(5,1),  sigmaXYff(5,1), concFF(5,1), SQRT(UxFOUR(5,1)**2.0+UyFOUR(5,1)**2.0), &
				SQRT(ORIXff(5,1)**2.0+ORIYff(5,1)**2.0), TRACEff(5,1), FNSCff(5,1)
		WRITE (PWEE1,*) sigmaXXff(5,3)  ,sigmaYYff(5,3),  sigmaXYff(5,3), concFF(5,3), SQRT(UxFOUR(5,3)**2.0+UyFOUR(5,3)**2.0), &
				SQRT(ORIXff(5,3)**2.0+ORIYff(5,3)**2.0), TRACEff(5,3), FNSCff(5,3)
		WRITE (PWEE1,*) sigmaXXff(1,6)  ,sigmaYYff(1,6),  sigmaXYff(1,6), concFF(1,6), SQRT(UxFOUR(1,6)**2.0+UyFOUR(1,6)**2.0), &
				SQRT(ORIXff(1,6)**2.0+ORIYff(1,6)**2.0), TRACEff(1,6), FNSCff(1,6)
		WRITE (PWEE1,*) sigmaXXff(7,1)  ,sigmaYYff(7,1),  sigmaXYff(7,1), concFF(7,1), SQRT(UxFOUR(7,1)**2.0+UyFOUR(7,1)**2.0), &
				SQRT(ORIXff(7,1)**2.0+ORIYff(7,1)**2.0), TRACEff(7,1), FNSCff(7,1)
		WRITE (PWEE1,*) sigmaXXff(1,8)  ,sigmaYYff(1,8),  sigmaXYff(1,8), concFF(1,8), SQRT(UxFOUR(1,8)**2.0+UyFOUR(1,8)**2.0), &
				SQRT(ORIXff(1,8)**2.0+ORIYff(1,8)**2.0), TRACEff(1,8), FNSCff(1,8)
		WRITE (PWEE1,*) sigmaXXff(5,5)  ,sigmaYYff(5,5),  sigmaXYff(5,5), concFF(5,5), SQRT(UxFOUR(5,5)**2.0+UyFOUR(5,5)**2.0), &
				SQRT(ORIXff(5,5)**2.0+ORIYff(5,5)**2.0), TRACEff(5,5), FNSCff(5,5)
		WRITE (PWEE1,*) sigmaXXff(6,6)  ,sigmaYYff(6,6),  sigmaXYff(6,6), concFF(6,6), SQRT(UxFOUR(6,6)**2.0+UyFOUR(6,6)**2.0), &
				SQRT(ORIXff(6,6)**2.0+ORIYff(6,6)**2.0), TRACEff(6,6), FNSCff(6,6)
		WRITE (PWEE1,*) sigmaXXff(10,5)  ,sigmaYYff(10,5),  sigmaXYff(10,5), concFF(10,5), 					&
		SQRT(UxFOUR(10,5)**2.0+UyFOUR(10,5)**2.0), SQRT(ORIXff(10,5)**2.0+ORIYff(10,5)**2.0), TRACEff(10,5), FNSCff(10,5)
		WRITE (PWEE1,*) sigmaXXff(13,13)  ,sigmaYYff(13,13),  sigmaXYff(13,13), concFF(13,13),					&
		SQRT(UxFOUR(13,13)**2.0+UyFOUR(13,13)**2.0), SQRT(ORIXff(13,13)**2.0+ORIYff(13,13)**2.0), TRACEff(13,13), FNSCff(13,13)

		WRITE (PW8,*) 	MAXVAL(TAUpXX+TAUpYY) ,MAXVAL(TAUpXY), MAXVAL(SIGMAXY), MAXVAL(TAUpXX),  MAXVAL(SIGMAXX), 			&
				MAXVAL(TAUpYY), MAXVAL(SIGMAYY)

		WRITE (GHG,*) 	TIME1, MAXVAL(CONC), MINVAL(conc), VEL, MINVAL(TAUpXX+TAUpYY), MINVAL(TAUpXY), MINVAL(SIGMAXY),			&
				MINVAL(TAUpXX), MINVAL(SIGMAXX), MINVAL(TAUpYY),  MINVAL(SIGMAYY)

		WRITE (W99,*)	MAXVAL(SIGMAXX-TAUpXX), MAXVAL(SIGMAYY-TAUpYY), MAXVAL(SIGMAXY-TAUpXY)					, 	&
				MINVAL(SIGMAXX-TAUpXX), MINVAL(SIGMAYY-TAUpYY), MINVAL(SIGMAXY-TAUpXY)

		WRITE (W7, *) 	SUM(ABS(TAUpXY))/(GRIDX*GRIDY), SUM(ABS(SIGMAXY))/(GRIDX*GRIDY), 						&
				SUM(ABS(TAUpXX))/(GRIDX*GRIDY), SUM(ABS(SIGMAXX))/(GRIDX*GRIDY),						&
				SUM(ABS(TAUpYY))/(GRIDX*GRIDY), SUM(ABS(SIGMAYY))/(GRIDX*GRIDY)

		WRITE (W9, *)	TIME1, VEL, MAXVAL(SIGMAXY-TAUpXY), MINVAL(SIGMAXY-TAUpXY), MINVAL(CONC), MAXVAL(CONC)

		!PRINT      *, 	TIME1, VEL, nxyN ,MAXVAL(SIGMAXY), MINVAL(SIGMAXY), MINVAL(CONC), MAXVAL(CONC), delUxyV/delUv, delUdif/delUv, &
		!		delUdif, delUv

		PRINT      *, 	TIME1, numer, VEL ,SUM(ABS(SIGMAXY))/(GRIDX*GRIDY), SUM(ABS(SIGMAXX))/(GRIDX*GRIDY),&
				SUM(ABS(TAUpXY))/(GRIDX*GRIDY),	SUM(ABS(TAUpXX))/(GRIDX*GRIDY),MAXVAL(CONC), MINVAL(CONC)

	END IF
	VEL  =0.0
	DO I=1, GRIDX
		DO J=1, GRIDY
			VEL=VEL+SQRT(UXF(I,J)**2.0+UYF(I,J)**2.0)/(GRIDX*GRIDY)
		END DO
	END DO

	!PRINT      *, 	TIME1, VEL, nxyN ,MAXVAL(SIGMAXY), MINVAL(SIGMAXY), MAXVAL(TAUpXX+TAUpYY), MINVAL(TAUpXX+TAUpYY), 		&
	!	         MAXVAL(SIGMAXY-TAUpXY), MINVAL(SIGMAXY-TAUpXY), MAXVAL(SIGMAXX-TAUpXX), MINVAL(SIGMAXX-TAUpXX)
	!PRINT*, 	TIME1, VEL, MINVAL(TAUpXX+TAUpYY), MAXVAL(TAUpXX+TAUpYY), MAXVAL(SIGMAXY-TAUpXY), MINVAL(SIGMAXY-TAUpXY), 	&
	!		SUM(ABS(SIGMAXY))/(GRIDX*GRIDY), SUM(ABS(TAUpXY))/(GRIDX*GRIDY) ,MINVAL(CONC), MAXVAL(CONC)


	IF (TIME .EQ. stepN/3) THEN
		 OPEN(i1,FILE='nx2.txt')
		 OPEN(i2,FILE='ny2.txt')
		 OPEN(i3,FILE='CONCENTRATION2.txt')
		 OPEN(i4,FILE='SIGMAXY2.txt')
		 OPEN(i5,FILE='ux2.txt')
		 OPEN(i6,FILE='uy2.txt')
		 OPEN(I8,FILE='deldotnc1.txt')
		 OPEN(STY2,FILE='delU2.txt')
		 DO I=1, GRIDX
		 	DO J=1, GRIDY

			 	write (i1,*)  ORIX	(I,J)
			 	write (i2,*)  ORIY	(I,J)
			 	write (i3,*)  CONC	(I,J)
			 	write (i4,*)  SIGMAXY(I,J), SIGMAXX(I,J), SIGMAYY(I,J), SIGMAXX(I,J)+SIGMAYY(I,J), SIGMAXY(I,J)-TAUpXY(I,J), &
			 		      SIGMAXX(I,J)-TAUpXX(I,J)
				write (i5,*)  UXF	(I,J)
				write (i6,*)  UYF	(I,J)
				write (i8,*)  -DELdotNCX(I,J)	-DELdotNCY	(I,J)
				WRITE (STY2,*)  SQRT( delUXX(I,J)*delUXX(I,J)+delUYY(I,J)*delUYY(I,J) ), &
						SQRT( delUXY(I,J)*delUXY(I,J)+delUYX(I,J)*delUYX(I,J) ), &
						delUXX(I,J)-delUYY(I,J)
				WRITE (W52,*) TAUpXX(I,J)+TAUpYY(I,J), TAUpXX(I,J)-TAUpYY(I,J),LToldrXXff(I,J),LToldrXYff(I,J), &
					      CONVECxx(I,J),CONVECxy(I,J), TAUpXX(I,J), TAUpXY(I,J)
			END DO
		END DO

	END IF
	IF (TIME .EQ. 2*stepN/3) THEN
	         OPEN(g1,FILE='nx3.txt')
		 OPEN(g2,FILE='ny3.txt')
		 OPEN(g3,FILE='CONCENTRATION3.txt')
		 OPEN(g4,FILE='SIGMAXY3.txt')
		 OPEN(g5,FILE='ux3.txt')
		 OPEN(g6,FILE='uy3.txt')
		 OPEN(I9,FILE='deldotnc2.txt')
 		 OPEN(STY3,FILE='delU3.txt')
		 DO I=1, GRIDX
		 	DO J=1, GRIDY

			 	write (g1,*)  ORIX	(I,J)
			 	write (g2,*)  ORIY	(I,J)
			 	write (g3,*)  CONC	(I,J)
			 	write (g4,*)  SIGMAXY	(I,J), SIGMAXX(I,J), SIGMAYY(I,J), SIGMAXX(I,J)+SIGMAYY(I,J), SIGMAXY(I,J)-TAUpXY(I,J), &
			 		      SIGMAXX(I,J)-TAUpXX(I,J)
				write (g5,*)  UXF	(I,J)
				write (g6,*)  UYF	(I,J)
				write (I9,*)  -DELdotNCX(I,J)	-DELdotNCY	(I,J)
				WRITE (STY3,*)  SQRT( delUXX(I,J)*delUXX(I,J)+delUYY(I,J)*delUYY(I,J) ), &
						SQRT( delUXY(I,J)*delUXY(I,J)+delUYX(I,J)*delUYX(I,J) ), &
						delUXX(I,J)-delUYY(I,J)
				WRITE (W53,*) TAUpXX(I,J)+TAUpYY(I,J), TAUpXX(I,J)-TAUpYY(I,J),LToldrXXff(I,J),LToldrXYff(I,J), &
					      CONVECxx(I,J),CONVECxy(I,J), TAUpXX(I,J), TAUpXY(I,J)
			END DO
		END DO

	END IF

	addPSI1=SUM(PSI)/(gridx*gridy*phid)
	PSI = PSI - (addPSI1-1.0/TWOPI)

	CALL CPU_TIME ( time_end3 )
	elapsed3 =elapsed3 + time_end3 - time_begin3
END DO

do delTsampling=1,stepN/2,1
	Nsample(delTsampling)=0
	do t=stepN/2+1,stepN+1,delTsampling
		do i=1,tracer
			Nsample(delTsampling)=Nsample(delTsampling)+1
			if ((delTsampling+t)<stepN+2) then
				msd(delTsampling)=msd(delTsampling)+&
				( (XtracerNOp(i,t+delTsampling)-XtracerNOp(i,t))**2.0+(YtracerNOp(i,t+delTsampling)-YtracerNOp(i,t))**2.0)
			end if
		end do
	end do
end do
do k=1,stepN/2
	write(kjk4,*) Nsample(k), msd(k) ,msd(k)/(Nsample(k))
end do

 CALL CPU_TIME ( time_end )
 print*,  elapsed1, elapsed2, elapsed3, time_end - time_begin
 call ETIME(tarray, result)
 print *, result
 print *, tarray(1)
 print *, tarray(2)
 CLOSE (FGF2)
 CLOSE (FGF3)
 CLOSE (FGF4)
 CLOSE (FGF5)
 CLOSE (FGF)
 CLOSE (FGF1)
 CLOSE (PW8)
 CLOSE (GHG)
 CLOSE (KJK1)

 OPEN(PWnx,FILE='nx4.txt')
 OPEN(PWny,FILE='ny4.txt')
 OPEN(PW2 ,FILE='CONCENTRATION4.txt')
 OPEN(PW3 ,FILE='SIGMAXY4.txt')
 OPEN(PWux,FILE='ux4.txt')
 OPEN(PWuy,FILE='uy4.txt')
 OPEN(i7,FILE='deldotnc3.txt')
 OPEN(STY4,FILE='delU4.txt')
 DO I=1, GRIDX
 	DO J=1, GRIDY
		write (PWnx,*) ORIX	(I,J)
		write (PWny,*) ORIY	(I,J)
		write (PW2,*)  CONC	(I,J)
		write (PW3,*)  SIGMAXY	(I,J), SIGMAXX(I,J), SIGMAYY(I,J), SIGMAXX(I,J)+SIGMAYY(I,J), SIGMAXY(I,J)-TAUpXY(I,J), &
			       SIGMAXX(I,J)-TAUpXX(I,J)
		write (PWux,*) UXF	(I,J)
		write (PWuy,*) UYF	(I,J)
		write (I7,*)   -DELdotNCX(I,J)	-DELdotNCY	(I,J)
		WRITE (W54,*)  TAUpXX(I,J)+TAUpYY(I,J), TAUpXX(I,J)-TAUpYY(I,J),LToldrXXff(I,J),LToldrXYff(I,J), &
			       CONVECxx(I,J),CONVECxy(I,J), TAUpXX(I,J), TAUpXY(I,J)
		WRITE (STY4,*) SQRT( delUXX(I,J)*delUXX(I,J)+delUYY(I,J)*delUYY(I,J) ), &
			       SQRT( delUXY(I,J)*delUXY(I,J)+delUYX(I,J)*delUYX(I,J) ), &
			       delUXX(I,J)-delUYY(I,J)
 	END DO
END DO
 CLOSE (PWnx)
 CLOSE (PWny)
 CLOSE (PW2)
 CLOSE (PW3)
 CLOSE (PWux)
 CLOSE (PWuy)
 CLOSE (i7)
 CLOSE (i8)
 CLOSE (i9)

 CLOSE (h1)
 CLOSE (h2)
 CLOSE (h3)
 CLOSE (h4)
 CLOSE (h5)
 CLOSE (h6)
 CLOSE (i1)
 CLOSE (i2)
 CLOSE (i3)
 CLOSE (i4)
 CLOSE (i5)
 CLOSE (i6)
 CLOSE (g1)
 CLOSE (g2)
 CLOSE (g3)
 CLOSE (g4)
 CLOSE (g5)
 CLOSE (g6)
 close (psp)
 CLOSE (PSP1)
 CLOSE (PSP2)
 CLOSE (PSP3)
 CLOSE (PSP4)
 CLOSE (PSP5)
 CLOSE (PSP6)
 CLOSE (MAXIM)
 close (psp20)
 close (PWEE2)
 close (PWEE1)
 CLOSE (W13)
 CLOSE (W14)
 CLOSE (W4)
 CLOSE (W5)
 CLOSE (W52)
 CLOSE (W53)
 CLOSE (W54)
 CLOSE (PSP7)
 CLOSE (h7)
 CLOSE (W7)
 CLOSE (W99)
 CLOSE (W9)
 close (sty)
 close (sty1)
 close (sty2)
 close (sty3)
 close (sty4)
 close (stt)
 CLOSE (PSP8)
 close (kjk3)
 close (kjk4)
 close (kjk5)
 close (kjk6)
 close (kjk7)
 close (psp9)
 close (psp10)
 close (psp11)
 close (psp12)
 close (psp13)
 close (psp14)
 close (psp15)
 CALL CHDIR("/home/yaser/Desktop/RESEARCH/NonLinearDynamics/OldroydB2constant")

END PROGRAM CODE
