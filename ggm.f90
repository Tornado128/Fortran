  SUBROUTINE ggm(MODES,SIZEB, DELTAX, DELTAY, GRIDX, GRIDY, UXF, UYF, UxFOUR, UyFOUR, SIGMAXY, SIGMAXX, &
		SIGMAYY, TAUpXY, TAUpXX, TAUpYY, contild, elongationXff, shearXff, UelongationX, UshearX, &
		sigmaXYm, sigmaXXm, TAUpxyM, TAUpxxM, ZETA)
  IMPLICIT NONE

  INTEGER*8, INTENT(IN) :: GRIDX, GRIDY, MODES
  integer, parameter :: dp  = kind(1.0D0)
  integer, parameter :: dpc = kind( (1.d0,1.d0) )					     
  REAL(dp), INTENT(IN) :: SIZEB, DELTAX, DELTAY, ZETA
  REAL(dp), DIMENSION(GRIDX,GRIDY), INTENT(IN):: SIGMAXY, SIGMAXX, SIGMAYY
  REAL(DP), PARAMETER :: TWOPI = 6.283185307179586476925286766559005768394
  REAL(DP), PARAMETER :: PI    = 3.141592653589793238462643383279502884197
  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64) 
  integer i, j, numer1
  integer*8 plan
  REAL*8 R, K1, K2, KK , RRR, MC1, MC2, DUMXX, DUMYY, DENOM
  REAL*8 SIZEBX, SIZEBY        				            ! SEE MODULES.F90  ! optimaaly spped up FFT
  REAL*8, DIMENSION(:,:),ALLOCATABLE:: TMP1
  COMPLEX(dpc) PRESSURE
  COMPLEX(dpc), DIMENSION(:,:),ALLOCATABLE:: TMP2  			     !Fourier Space 
  COMPLEX(dpc), DIMENSION(:,:),ALLOCATABLE:: RHOX, RHOY, RHOXX, RHOYY		     !Fourier Space
  COMPLEX(dpc), DIMENSION(int( GRIDX/2 + 1),GRIDY) :: Utildex, Utildey
  COMPLEX(dpc), DIMENSION(int( GRIDX/2 + 1),GRIDY), intent(out) :: contild 
  REAL(dp), DIMENSION(GRIDX,GRIDY), INTENT(OUT) :: UXF, UYF, UshearX, UelongationX, sigmaXYm, sigmaXXm, TAUpxyM, TAUpxxM
  REAL(dp), DIMENSION(GRIDX,GRIDY), INTENT(IN)  :: TAUpXY, TAUpXX, TAUpYY
  REAL(dp), DIMENSION(int( GRIDX/2 + 1),GRIDY), INTENT(OUT) :: UxFOUR, UyFOUR, shearXff, elongationXff
  COMPLEX(dpc), DIMENSION(int( GRIDX/2 + 1),GRIDY) :: sigmaxyFFFF, sigmaxxFFFF, TAUpXYFFFF, TAUpXXFFFF, TAUpYYFFFF, &
						      SIGMAYYFFFF, shearXFFFF, elongationXFFFF
  integer, parameter :: pw5=1
  INTEGER*8, DIMENSION(6) :: ICC          ! It was 9. I made it 6!
     icc(1) = GRIDX
     icc(2) = GRIDY
     icc(3) = 1
     icc(4) = int( GRIDX/2 + 1);
     icc(5) = int( GRIDY/2 + 1);                          
     icc(6) = 1     
  
 
                                   

 allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)),       &     ! How should I change the dimension of tmp1?
       RHOX(ICC(4),ICC(2)), RHOY(ICC(4),ICC(2)), RHOXX(ICC(4),ICC(2)), RHOYY(ICC(4),ICC(2)) )

  TMP1 = SIGMAXY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE) ! Check 2
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  sigmaxyFFFF = TMP2/(GRIDx*GRIDy)

  TMP1 = SIGMAXX
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE) ! Check 2
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  sigmaxxFFFF = TMP2/(GRIDx*GRIDy)

  TMP1 = SIGMAYY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE) ! Check 2
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  sigmayyFFFF = TMP2/(GRIDx*GRIDy)

  TMP1 = TAUpXY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE) ! Check 2
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  TAUpXYFFFF = TMP2/(GRIDx*GRIDy)

  TMP1 = TAUpYY
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE) ! Check 2
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  TAUpYYFFFF = TMP2/(GRIDx*GRIDy)

  TMP1 = TAUpXX
  CALL dfftw_plan_dft_r2c_2d( plan, ICC(1), ICC(2), TMP1, TMP2, FFTW_ESTIMATE) ! Check 2
  CALL dfftw_execute(plan)
  CALL dfftw_destroy_plan(plan)
  TAUpXXFFFF = TMP2/(GRIDx*GRIDy)


RHOXx(1,1)=0.0
RHOYy(1,1)=0.0
Utildex(1,1)=0
Utildey(1,1)=0
 DO j = 1, icc(2)
        DO i = 1, icc(4)

           	k1 = dble(i - 1)*(twopi/sizeb)
           	if (j <= icc(5) ) k2 = dble(j - 1)*(twopi/sizeb)
           	if (j >  icc(5) ) k2 = -dble(2*icc(5) - j - 1)*(twopi/sizeb)

		!IF ((K1 .EQ. 0) .AND. (K2 .EQ. 0)) THEN
		!	K1 = MODES *(twopi/sizeb)
		!	K2 = MODES *(twopi/sizeb)
		!END IF
		
		!if (i == icc(4)) k1 = 0.d0 ! correction made by PTU, 09/06/07
		!if (j == icc(5)) k2 = 0.d0

		kk = k1*k1 + k2*k2   !k^2

	   	DUMXX = k1/( (kk**0.5))
           	DUMYY = k2/( (kk**0.5))
	   	DENOM = ( KK )**0.5
		IF (KK .NE. 0) THEN
			RHOXx(I,J) = ( CMPLX(0,1)/DENOM ) * ( (sigmaxyFFFF(I,J)-TAUpXYFFFF(I,J))*( DUMYY**3.0-DUMYY*DUMXX*DUMXX ) &
			+DUMXX*(DUMYY**2.0)*( sigmaxxFFFF(I,J)-TAUpXXFFFF(I,J) ) - DUMXX*(DUMYY**2.0)*( sigmayyFFFF(I,J)-TAUpYYFFFF(I,J) ) )/(1-ZETA)

			shearXffff(i,j)=( CMPLX(0,1)/DENOM )*( (sigmaxyFFFF(I,J)-TAUpXYFFFF(I,J))*( DUMYY**3.0-DUMYY*DUMXX*DUMXX ))
			elongationXffff(i,j)=( CMPLX(0,1)/DENOM )*(DUMXX*(DUMYY**2.0)*( sigmaxxFFFF(I,J)-TAUpXXFFFF(I,J) ) 	&
			- DUMXX*(DUMYY**2.0)*( sigmayyFFFF(I,J)-TAUpYYFFFF(I,J) ) )
			

			RHOYy(I,J) = ( CMPLX(0,1)/DENOM ) * ( (sigmaxyFFFF(I,J)-TAUpXYFFFF(I,J))*( DUMXX**3.0-DUMXX*DUMYY*DUMYY ) &
			+DUMYY*(DUMXX**2.0)*( sigmayyFFFF(I,J)-TAUpYYFFFF(I,J) ) - DUMYY*(DUMXX**2.0)*( sigmaxxFFFF(I,J)-TAUpXXFFFF(I,J) ) )/(1-ZETA)

			Utildex(I,J)=RHOXx(I,J)
			Utildey(I,J)=RHOYy(I,J)

			contild(I,J)=CMPLX(0,1)*RHOXx(I,J)*DUMXX+CMPLX(0,1)*RHOYy(I,J)*DUMYY
		END IF
	END DO
END DO

DO I=1, ICC(4)
	DO J=1, ICC(2)
		MC1=AIMAG( Utildex(I,J)    )
		MC2=REAL ( Utildex(I,J)    )
		UxFOUR(I,J)=( MC1**2.0+MC2**2.0 )**0.5
		MC1=AIMAG( Utildey(I,J)    )
		MC2=REAL ( Utildey(I,J)    )
		UyFOUR(I,J)=( MC1**2.0+MC2**2.0 )**0.5

		MC1=AIMAG( shearXffff(I,J)    )
		MC2=REAL ( shearXffff(I,J)    )
		shearXff(I,J)=( MC1**2.0+MC2**2.0 )**0.5
		MC1=AIMAG( elongationXffff(I,J)    )
		MC2=REAL ( elongationXffff(I,J)    )
		elongationXff(I,J)=( MC1**2.0+MC2**2.0 )**0.5
	END DO
END DO 

DO I=1, ICC(4)
	DO J=1, ICC(2)
		if (I .eq. J) then
			sigmaxyFFFF(I,J)=0.0
			 TAUpxyFFFF(I,J)=0.0
		end if
		if ((I .eq. 1) .or. (J .eq. 1)) then
			 TAUpYYFFFF(I,J)=0.0
			 TAUpXXFFFF(I,J)=0.0
			sigmaXXFFFF(I,J)=0.0
			sigmaYYFFFF(I,J)=0.0
		end if
	end do
end do

  deallocate(tmp1, tmp2)
  allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)) )
  tmp2 = sigmaxyFFFF
  call dfftw_plan_dft_c2r_2d( plan, icc(1), icc(2), tmp2, tmp1, FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  sigmaXYm(1:icc(1),1:icc(2)) = tmp1 

  !uxf(icc(1),:)=uxf(1,:)
  !uxf(:,icc(2))=uxf(:,1)

  deallocate(tmp1,  tmp2)
  allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)) )
  !uy-far REAL space
  tmp2 = TAUpxyFFFF
  call dfftw_plan_dft_c2r_2d( plan, icc(1), icc(2), tmp2, tmp1, FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  TAUpxyM(1:icc(1),1:icc(2)) = tmp1			




  deallocate(tmp1, tmp2)
  allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)) )
  tmp2 = sigmaxxFFFF
  call dfftw_plan_dft_c2r_2d( plan, icc(1), icc(2), tmp2, tmp1, FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  sigmaXXm(1:icc(1),1:icc(2)) = tmp1 

  !uxf(icc(1),:)=uxf(1,:)
  !uxf(:,icc(2))=uxf(:,1)

  deallocate(tmp1,  tmp2)
  allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)) )
  !uy-far REAL space
  tmp2 = TAUpxxFFFF
  call dfftw_plan_dft_c2r_2d( plan, icc(1), icc(2), tmp2, tmp1, FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  TAUpxxM(1:icc(1),1:icc(2)) = tmp1





  deallocate(tmp1, tmp2)
  allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)) )
  tmp2 = RHOXx
  call dfftw_plan_dft_c2r_2d( plan, icc(1), icc(2), tmp2, tmp1, FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  uxf(1:icc(1),1:icc(2)) = tmp1 

  !uxf(icc(1),:)=uxf(1,:)
  !uxf(:,icc(2))=uxf(:,1)

  deallocate(tmp1,  tmp2)
  allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)) )
  !uy-far REAL space
  tmp2 = RHOYy
  call dfftw_plan_dft_c2r_2d( plan, icc(1), icc(2), tmp2, tmp1, FFTW_ESTIMATE)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)
  uyf(1:icc(1),1:icc(2)) = tmp1

  !!uyf(icc(1),:)=uyf(1,:)
  !!uyf(:,icc(2))=uyf(:,1)

  !deallocate(tmp1, tmp2)
  !allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)) )
  !tmp2 = shearXffff
  !call dfftw_plan_dft_c2r_2d( plan, icc(1), icc(2), tmp2, tmp1, FFTW_ESTIMATE)
  !call dfftw_execute(plan)
  !call dfftw_destroy_plan(plan)
  !UshearX(1:icc(1),1:icc(2)) = tmp1 

  !uxf(icc(1),:)=uxf(1,:)
  !uxf(:,icc(2))=uxf(:,1)

  !deallocate(tmp1,  tmp2)
  !allocate( TMP1(ICC(1),ICC(2)), TMP2(ICC(4),ICC(2)) )
  !!uy-far REAL space
  !tmp2 = elongationXffff
  !call dfftw_plan_dft_c2r_2d( plan, icc(1), icc(2), tmp2, tmp1, FFTW_ESTIMATE)
  !call dfftw_execute(plan)
  !call dfftw_destroy_plan(plan)
  !UelongationX(1:icc(1),1:icc(2)) = tmp1

  !!uyf(icc(1),:)=uyf(1,:)
  !!uyf(:,icc(2))=uyf(:,1)

deallocate( tmp1, tmp2, rhox, rhoy, RHOXX, RHOYY)

end subroutine ggm
