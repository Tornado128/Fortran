  SUBROUTINE    fftwCheck(conc,sigmaxy,TAUpXY,SIZEB, DELTAX, DELTAY, GRIDX, GRIDY, concFF &
	     	,sigmaxx,TAUpXX, sigmayy, TAUpYY, stressXXff, stressYYff, stressXYff, mixing &
	     	,sigmaxxFF, sigmaxyFF, sigmayyFF, TAUpXYff, TAUpYYff, TAUpXXff, ORIX, ORIY, ORIXff &
	     	,ORIYff, TRACEff, FNSCff, LToldrXX, LToldrXY, LToldrXXff, LToldrXYff  &
		,CONVECxx, CONVECxy, CONVECxxFF, CONVECxyFF, orientXX, orientXY, orientXXff, orientXYff)

  IMPLICIT NONE

  INTEGER*8, INTENT(IN) :: GRIDX, GRIDY
  integer, parameter :: dp  = kind(1.0D0)
  integer, parameter :: dpc = kind( (1.d0,1.d0) )					     
  REAL(dp), INTENT(IN) :: SIZEB, DELTAX, DELTAY 
  REAL(dp), DIMENSION(GRIDX,GRIDY), INTENT(IN) 			:: conc, sigmaxy, TAUpXY, sigmaxx,TAUpXX, sigmayy,TAUpYY, ORIX, ORIY, &
								   LToldrXX, LToldrXY, CONVECxx, CONVECxy, orientXX, orientXY
  COMPLEX(dpc), DIMENSION(int( GRIDX/2 + 1),GRIDY) 		:: concFFF, stressXXfff, stressYYfff, stressXYfff, sigmaXXfff, 	&
						   		   sigmaYYfff, sigmaXYfff, TauPxyfff, TauPxxfff, TauPyyfff, 	&
								   ORIXfff, ORIYfff, TRACEfff, FNSCfff, LToldrXXfff, LToldrXYfff,  &
								   CONVECxxFFF, CONVECxyFFF, orientXXfff, orientXYfff
  REAL(dpc), DIMENSION(int( GRIDX/2 + 1),GRIDY), INTENT(OUT)	:: concFF, stressXXff, stressYYff, stressXYff, ORIXff, ORIYff, 	&
						   		   sigmaXXff, sigmaYYff, sigmaXYff, TauPxyff, TauPxxff, TauPyyff,  &
							           TRACEff, FNSCff, CONVECxxFF, CONVECxyFF, LToldrXXff, LToldrXYff, &
								   orientXXff, orientXYff
  REAL(DP), PARAMETER :: TWOPI = 6.283185307179586476925286766559005768394
  REAL(DP), PARAMETER :: PI    = 3.141592653589793238462643383279502884197
  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64) 
  integer i, j
  integer*8 plan
  REAL*8 R, K1, K2, KK, MC1, MC2, MS1, MS2, mixing   				            ! SEE MODULES.F90  ! optimaaly spped up FFT
  REAL*8, DIMENSION(:,:),ALLOCATABLE:: TMP1
  complex(dpc) RHX, RHY
  COMPLEX(dpc), DIMENSION(:,:),ALLOCATABLE:: TMP2  			     !Fourier Space 
 
  integer, parameter :: pw5=1
  INTEGER*8, DIMENSION(6) :: ICC          ! It was 9. I made it 6!
     icc(1) = GRIDX
     icc(2) = GRIDY
     icc(3) = 1
     icc(4) = int( GRIDX/2 + 1);
     icc(5) = int( GRIDY/2 + 1);                          
     icc(6) = 1      
 
                                   
  allocate(   TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2))   )	

  TMP1 = sigmaxy - TAUpXY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE) ! Check 2
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  stressXYfff = TMP2/(GRIDx*GRIDy)

  TMP1 = conc
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  concFFF = TMP2/(GRIDx*GRIDy)

  TMP1 = ORIX
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  ORIXfff = TMP2/(GRIDx*GRIDy)

  TMP1 = ORIY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  ORIYfff = TMP2/(GRIDx*GRIDy)

  TMP1 = sigmaxx - TAUpXX
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  stressXXfff = TMP2/(GRIDx*GRIDy)
  
  TMP1 = sigmayy - TAUpYY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  stressYYfff = TMP2/(GRIDx*GRIDy)



  TMP1 = sigmaxy
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  sigmaxyFFF = TMP2/(GRIDx*GRIDy)

  TMP1 = sigmaxx
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  sigmaXXfff = TMP2/(GRIDx*GRIDy)
  
  TMP1 = sigmayy
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  sigmaYYfff = TMP2/(GRIDx*GRIDy)




  TMP1 = TAUpXY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  TAUpXYfff = TMP2/(GRIDx*GRIDy)

  TMP1 = TAUpXX
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  TAUpXXfff = TMP2/(GRIDx*GRIDy)
  
  TMP1 = TAUpYY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  TAUpYYfff = TMP2/(GRIDx*GRIDy)

  TMP1 = TAUpXX+TAUpYY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  TRACEfff = TMP2/(GRIDx*GRIDy)
  
  TMP1 = TAUpXX-TAUpYY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  FNSCfff = TMP2/(GRIDx*GRIDy)

  TMP1 = LToldrXX
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  LToldrXXfff = TMP2/(GRIDx*GRIDy)
  
  TMP1 = LToldrXY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  LToldrXYfff = TMP2/(GRIDx*GRIDy)

  TMP1 = convecXX
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  convecXXfff = TMP2/(GRIDx*GRIDy)
  
  TMP1 = convecXY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  convecXYfff = TMP2/(GRIDx*GRIDy)

  TMP1 = orientXX
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  orientXXfff = TMP2/(GRIDx*GRIDy)
  
  TMP1 = orientXY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE)
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  orientXYfff = TMP2/(GRIDx*GRIDy)

 mixing =0.0
 DO I=1, ICC(4)
	DO J=1, ICC(2)
		MC1=AIMAG( concFFF(I,J)    )
		MC2=REAL ( concFFF(I,J)    )
		concFF(I,J)=( MC1**2.0+MC2**2.0 )**0.5
		MS1=AIMAG( stressYYfff(I,J) )
		MS2=REAL ( stressYYfff(I,J) )
		stressYYff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( stressXXfff(I,J) )
		MS2=REAL ( stressXXfff(I,J) )
		stressXXff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( stressXYfff(I,J) )
		MS2=REAL ( stressXYfff(I,J) )
		stressXYff(I,J)=( MS1**2.0+MS2**2.0 )**0.5


		MS1=AIMAG( sigmaxyFFF(I,J) )
		MS2=REAL ( sigmaxyFFF(I,J) )
		sigmaxyFF(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( sigmayyFFF(I,J) )
		MS2=REAL ( sigmayyFFF(I,J) )
		sigmayyFF(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( sigmaxxFFF(I,J) )
		MS2=REAL ( sigmaxxFFF(I,J) )
		sigmaxxFF(I,J)=( MS1**2.0+MS2**2.0 )**0.5


		MS1=AIMAG( TAUpXYfff(I,J) )
		MS2=REAL ( TAUpXYfff(I,J) )
		TAUpXYff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( TAUpXXfff(I,J) )
		MS2=REAL ( TAUpXXfff(I,J) )
		TAUpXXff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( TAUpYYfff(I,J) )
		MS2=REAL ( TAUpYYfff(I,J) )
		TAUpYYff(I,J)=( MS1**2.0+MS2**2.0 )**0.5

		MS1=AIMAG( ORIXfff(I,J) )
		MS2=REAL ( ORIXfff(I,J) )
		ORIXff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( ORIYfff(I,J) )
		MS2=REAL ( ORIYfff(I,J) )
		ORIYff(I,J)=( MS1**2.0+MS2**2.0 )**0.5

		MS1=AIMAG( TRACEfff(I,J) )
		MS2=REAL ( TRACEfff(I,J) )
		TRACEff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( FNSCfff(I,J) )
		MS2=REAL ( FNSCfff(I,J) )
		FNSCff (I,J)=( MS1**2.0+MS2**2.0 )**0.5


		MS1=AIMAG( convecXXfff(I,J) )
		MS2=REAL ( convecXXfff(I,J) )
		convecXXff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( convecXYfff(I,J) )
		MS2=REAL ( convecXYfff(I,J) )
		convecXYff (I,J)=( MS1**2.0+MS2**2.0 )**0.5

		MS1=AIMAG( LToldrXXfff(I,J) )
		MS2=REAL ( LToldrXXfff(I,J) )
		LToldrXXff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( LToldrXYfff(I,J) )
		MS2=REAL ( LToldrXYfff(I,J) )
		LToldrXYff(I,J)=( MS1**2.0+MS2**2.0 )**0.5

		MS1=AIMAG( orientXXfff(I,J) )
		MS2=REAL ( orientXXfff(I,J) )
		orientXXff(I,J)=( MS1**2.0+MS2**2.0 )**0.5
		MS1=AIMAG( orientXYfff(I,J) )
		MS2=REAL ( orientXYfff(I,J) )
		orientXYff(I,J)=( MS1**2.0+MS2**2.0 )**0.5

		mixing       = mixing + SQRT( concFF(I,J)*concFF(I,J)/SQRT( 1.0+(I*I+J*J)*4.0*PI*PI/(SIZEB*SIZEB) ))
	END DO
 END DO 

deallocate( tmp1, tmp2)

end subroutine fftwCheck
