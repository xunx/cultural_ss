PROGRAM LAND
IMPLICIT NONE
!$DEBUG

! THIS PROGRAM COMPUTES A STEADY-STATE EQUILIBRIUM OF THE `LAND' PAPER,
! ENTITLED `SOCIAL SECURITY IN AN OVERLAPPING GENERATIONS ECONOMY WITH LAND'
! FORTHCOMING IN REVIEW OF ECONOMIC DYNAMICS. IN PARTICULAR, THIS CODE 
! GENERATES THE ROW `THETA=0.4' IN TABLE 3.
! SIMPLY BY CHANGING A STATEMENT BELOW INVOLVING THETA, ONE CAN GENERATE 
! ALL OF THE STEADY-STATES IN TABLE 3.

!******************
!
!   Initial Guess Values
!
!******************

REAL, PARAMETER :: K0   = 2.30  !  Capital stock
REAL, PARAMETER :: BEQ0 = 0.0420         !  Accidental bequests

!******************
!
!   Set Parameters
!
!******************

REAL, PARAMETER :: TOLK = 0.001       !  Convergence tolerance for capital stock
REAL, PARAMETER :: TOLB = 0.001       !  Convergence tolerance for bequests
REAL, PARAMETER :: GRADK = 0.2        !  Convergence gradient for capital stock
REAL, PARAMETER :: GRADB = 0.6        !  Convergence gradient for bequests
INTEGER, PARAMETER :: MAXITER = 50       !  Maximum number of iterations for convergence

REAL, PARAMETER :: ALPHA  = 0.690     !  Labor exponent in production function
REAL, PARAMETER :: TFP    = 0.9     !  Multiplicative constant in production function
REAL, PARAMETER :: GROWTH = 0.0165    !  Growth rate of per capita output
REAL, PARAMETER :: DEP    = 0.069     !  Depreciation rate

REAL, PARAMETER :: BETA   = 0.978     !  Subjective discount factor
REAL, PARAMETER :: GAMMA  = 2.0       !  Risk aversion parameter

!PARAMETER (THETA  = 0.40)      !  Social security replacement rate

INTEGER, PARAMETER :: MAXAGE = 65        !  Maximum age allowed
INTEGER, PARAMETER :: RETAGE = 45        !  Retirement age  
REAL, PARAMETER :: HBAR   = 1.00      !  Exogenous hours of work
REAL, PARAMETER :: RHO    = 0.012     !  Population growth rate

REAL, PARAMETER :: AMAX   = 40.0      !  Maximum permissible asset
INTEGER, PARAMETER :: NGRID  = 4097      !  Number of points on asset grid

INTEGER, PARAMETER :: KLOW   = 10        !  10 times smallest initial capital in fixed-point diagram
INTEGER, PARAMETER :: KHIGH  = 60        !  10 times largest initial capital in fixed-point diagram

REAL, PARAMETER :: CMIN   = 0.00005   !  Minimum permissible consumption
!PARAMETER (CMAX   = 12.0)      !  Maximum permissible consumption
! below, we use JMAX = 1.5 * AMAX
REAL, PARAMETER :: CINCR  = 0.001     !  Consumption increment for utility tabulation
INTEGER, PARAMETER :: JMAX = 60000	! JMAX = 1.5 * AMAX / CINCR

!********************
!
!   Data Type Declarations and Dimension Statements
!
!********************

INTEGER AGE, IL, IS, ISKIP, IU, JAMAX, MODE, IA, CON2, I, J

REAL  ATR, BEQ1, BEQEND, EXDEM, EXPCTDUT, INT, INV, K, K1, CON1, CON3, CUM
REAL  KDEV, L, OUTPUT, SS, STAX, VMAX, WAGE, WAGEZ, X3
REAL BEQ, CONS, THETA, AGROWTH, AEND, BEQDEV, ABEG, ACONS, UTIL, VTEMP
INTEGER ITHETA, ITER, ITERINC, JA


INTEGER IDCR(RETAGE:MAXAGE,NGRID)      !  Asset decision rules for retirees
INTEGER IDCW(RETAGE-1,NGRID)    !  Asset decision rules for working-age agents

REAL A(NGRID)              !  Asset levels
REAL ACROSS(MAXAGE)         !  Cross-sectional age-assets profile
REAL ALONG(MAXAGE)          !  Longitudinal age-assets profile
REAL CCROSS(MAXAGE)         !  Cross-sectional age-consumption profile
REAL CLONG(MAXAGE)          !  Longitudinal age-consumption profile
REAL CUMS(MAXAGE)           !  Unconditional survival probabilities, age 1 to age j
REAL EFFCROSS(RETAGE-1)       !  Cross-sectional age-earnings profile
REAL EFFLONG(RETAGE-1)        !  Longitudinal age-earnings profile for given cohort
REAL ICROSS(MAXAGE)         !  Cross-sectional age-income profile
REAL ILONG(MAXAGE)          !  Longitudinal age-income profile
REAL MU(MAXAGE)             !  Age distribution of population
REAL S(MAXAGE)              !  Conditional survival probabilities, age j-1 to age j
REAL UT(JMAX)             !Values in utility function lookup table
REAL VR(RETAGE:MAXAGE,NGRID)           !  Value function for retirees
REAL VW(RETAGE-1,NGRID)         !  Value function for working-age agents

REAL DLONG(MAXAGE)			! position of optimal Kt+1 for each age. dlong(1) is end of period capital at age 1 


!*********************
!0. Open Files


MODE = 2

IF (MODE==1) THEN   ! This is HP linux
	OPEN(UNIT=7,FILE='/home/sean/Desktop/Dropbox/cultural_ss/comboeff.txt')
	OPEN(UNIT=8,FILE='/home/sean/Desktop/Dropbox/cultural_ss/psurv.txt')
	OPEN(UNIT=10,FILE='/home/sean/Desktop/Dropbox/cultural_ss/result.txt')
	OPEN(UNIT=18,FILE='/home/sean/Desktop/Dropbox/cultural_ss/profile.txt')
ELSE IF (MODE==2) THEN  ! This is my bootcamp Windows
	OPEN(UNIT=7,FILE='C:\Users\Sean\Desktop\Dropbox\cultural_ss\comboeff.txt')
	OPEN(UNIT=8,FILE='C:\Users\Sean\Desktop\Dropbox\cultural_ss\psurv.txt')
	OPEN(UNIT=10,FILE='C:\Users\Sean\Desktop\Dropbox\cultural_ss\result.txt')
	OPEN(UNIT=18,FILE='C:\Users\Sean\Desktop\Dropbox\cultural_ss\profile.txt')
END IF





!***********************
! 1.1 tabulate spaces & misc variables

	READ(7,*) ( EFFCROSS(AGE), AGE=1,RETAGE-1 )
	READ(8,*) ( S(AGE), AGE=1,MAXAGE )

!   Tabulate asset levels

     A = (/ ( (FLOAT(IA-1)*AMAX)/FLOAT(NGRID-1), IA=1,NGRID ) /)

!   Tabulate utility function

     DO I=1,JMAX
        CONS = CMIN + (I-1)*CINCR
        UT(I) = (CONS**(1.0-GAMMA))/(1.0-GAMMA)
     END DO

!   Longitudinal age-earnings profile for given cohort

     EFFLONG = (/ ( ((1+GROWTH)**(AGE-1))*EFFCROSS(AGE), AGE=1,RETAGE-1 ) /)

!   Growth rate of aggregate output

     AGROWTH = (1.0+GROWTH)*(1.0+RHO) - 1.0


!***********************
! 1.2 values for K, L, beq



! initial guess for K, beq

     K = K0
     BEQ = BEQ0



! compute pop to find L
!  Unconditional survival probabilities


     CUMS(1) = 1.0
     DO J=2,MAXAGE
        CUMS(J) = CUMS(J-1)*S(J)
     END DO
!  Age distribution of population

     CUM = 0.0
     DO AGE=1,MAXAGE
        CUM = CUM + CUMS(AGE)/((1.0+RHO)**(AGE-1))
     END DO

     MU(1) = 1.0/CUM
     DO AGE=2,MAXAGE
        MU(AGE) = S(AGE)*MU(AGE-1)/(1.0+RHO)
     END DO
!   Labor input

     L = 0.0
     DO AGE=1,RETAGE-1
        L = L + HBAR*EFFCROSS(AGE)*MU(AGE)
     END DO
 




!**************************************
! starting parameter loop

! 2.0 ss loop
 DO ITHETA = 4,4,1                            ! 1025-point grid
     THETA = FLOAT(ITHETA)/10                     ! 1025-point grid
     PRINT *, 'THETA =', THETA                      ! 1025-point grid


     ITER = 1


 200 ITERINC = 0



! 2.1 Factor prices

     WAGE = ALPHA*TFP*(L**(ALPHA-1.0))*(K**(1.0-ALPHA))
	 PRINT*, 'WAGE=', WAGE
     INT = (1.0-ALPHA)*TFP*(L**ALPHA)*(K**(-ALPHA)) - DEP
	 PRINT*, 'INT=', INT
     IF (INT<AGROWTH) THEN
        PRINT *, "Interest rate is less than growth rate."
        PRINT *, "Program terminates."
          WRITE(10,*) "WARNING:  Interest rate is less than growth rate."
 !       GO TO 999
     END IF

     ATR = 1.0 + INT    



! 2.2 ss benefits and ss tax rate                                                                 
!   Social security benefits for a given cohort, constant across their
!      retirement years

     CUM = 0.0
     DO AGE=1,RETAGE-1
        CUM = CUM + EFFLONG(AGE)
     END DO
     SS = THETA*CUM*WAGE*HBAR/(FLOAT(RETAGE-1))

!   Social security tax rate

     CUM = 0.0                                    !  Downweight aggregate
     DO AGE=RETAGE,MAXAGE                         !  benefits because benefits
        CUM = CUM + MU(AGE)*(1+GROWTH)**(1-AGE)   !  decline cross-sectionally
     END DO                                       !  with age
     STAX = SS*CUM
     CUM = 0.0                                    !  Divide aggregate benefits
     DO AGE=1,RETAGE-1                            !  by total revenues, which
        CUM = CUM + MU(AGE)*EFFCROSS(AGE)         !  depend on cross-sectional
     END DO                                       !  age-earnings profile of
     STAX = STAX/(WAGE*HBAR*CUM)             !  this period's workers
            





!2.3 Find decision rules


     CALL DECRULE01                     

!   Calculate expected lifetime utility

     EXPCTDUT = VW(1,1)



!      DO AGE=1,RETAGE-1
! 		 DO IA=1,NGRID
!            IF (IDCW(AGE,IA)>=NGRID) THEN	! note it only checks the points on the decision rule
!               PRINT *, 'Error:  Maximum asset limit is binding.'	! if it's off the decision rule, it doesn't matter
!               PRINT *, 'AGE =', AGE
!               PRINT *, 'IA =', IA
!               GO TO 999
!            END IF
! 	   	 END DO
!      END DO



!2.4 Compute age profiles and average lifetime utility


     CALL PROFILE01

	 DO AGE=1,MAXAGE
		 IF (DLONG(AGE)>=NGRID) THEN
             PRINT *, 'Error:  Maximum asset limit is binding.'	! if it's off the decision rule, it doesn't matter
             PRINT *, 'AGE =', AGE
			 GO TO 999
		 END IF
	 END DO

!2.5 sum up assets and bequests

     AEND = 0.0
     BEQEND = 0.0

     DO AGE=1,MAXAGE-1
        AEND = AEND + ACROSS(AGE)*MU(AGE)
        BEQEND = BEQEND + ACROSS(AGE)*MU(AGE)*(1.0-S(AGE+1))
     END DO

     OUTPUT = TFP*(L**ALPHA)*(K**(1.0-ALPHA))
     K1 = AEND/(1.0+AGROWTH)
     BEQ1 = BEQEND/(1.0+AGROWTH)

!   BEQ is the bequest received per effective labor unit that enters the
!       budget constraint.  The land component is still valued at its price
!       from the previous period.  Thus, a single rate of return applies to
!       the land and capital components of both bequests and beginning-of-
!       period assets.
!   BEQEND is bequests per period-t effective labor unit left at the end
!       of period t.
!   BEQ1 is bequests per period-(t+1) effective labor unit received at the
!       beginning of period t+1, with the land component still valued at 
!       its period-t price.  This is the appropriate quantity to compare 
!       with BEQ in judging convergence.

     KDEV = ABS(K-K1)/K
     BEQDEV = ABS(BEQ-BEQ1)/BEQ

     PRINT *, 'Iteration', ITER
     PRINT *, 'Initial capital =', K
     PRINT *, 'Ending capital =', K1
     PRINT *, 'Relative change in capital =', KDEV
     PRINT *, 'Initial bequests =', BEQ
     PRINT *, 'Ending bequests =', BEQ1
     PRINT *, 'Relative change in bequests =', BEQDEV

     IF ( (KDEV>TOLK) .OR. (BEQDEV>TOLB) ) THEN
        K = (1 - GRADK)*K + GRADK*K1
        BEQ = (1 - GRADB)*BEQ + GRADB*BEQ1
        ITERINC = 1
     END IF

     IF (ITERINC>0) THEN

        ITER = ITER + ITERINC

        IF (ITER>MAXITER) THEN
           PRINT *, 'Maximum number of iterations exceeded.'
           PRINT *, 'Program terminates.'
           GO TO 799
        END IF

        GO TO 200

     END IF

!*********************************
! 4. final summary to show results

!   Compute average beginning-of-period assets

!   ABEG is assets per period-t effective labor unit with the land component
!      still valued at its price from period t-1.  The capital gain on land
!      occurs after ABEG is measured, and is included in the rate of return
!      in the agent's budget constraint.

 799 ABEG = 0.0

!      DO AGE=2,RETAGE-1
!         DO IA=1,NGRID
!               ABEG = ABEG + A(IA)*MU(AGE)*((1+GROWTH)**(1-AGE))
!         END DO
!      END DO
! 
!      DO AGE=RETAGE,MAXAGE
!         DO IA = 1,NGRID
!            ABEG = ABEG + A(IA)*MU(AGE)*((1+GROWTH)**(1-AGE))
!         END DO
!      END DO

!   Compute average consumption

     ACONS = 0.0
     DO AGE=1,MAXAGE
          ACONS = ACONS + CCROSS(AGE)*MU(AGE)
     END DO

!   Check adding-up constraints

     INV = (AGROWTH+DEP)*K
     EXDEM = ACONS + INV - OUTPUT

     WRITE(10,*) 'FINAL RESULTS'
     WRITE(10,*)

     WRITE(10,*) 'Labor input =', L
     WRITE(10,*)

     WRITE(10,*) 'Wage rate =', WAGE
     WRITE(10,*) 'Interest rate =', INT
     WRITE(10,*)

     WRITE(10,*) 'Social security benefit =', SS
     WRITE(10,*) 'Social security tax rate =', STAX
     WRITE(10,*)

     WRITE(10,*) 'Output =', OUTPUT
     print*, 'Output =', OUTPUT
     WRITE(10,*) 'Consumption =', ACONS
     WRITE(10,*) 'Investment =', INV
     WRITE(10,*) 'Excess demand =', EXDEM
     WRITE(10,*)

     WRITE(10,*) 'Capital (from ending assets) =', K1
!      WRITE(10,*) 'Capital (from beginning assets) =', ABEG+BEQ1/(1.0+AGROWTH)
     WRITE(10,*)

!   Because ABEG and BEQ1 have land valued at last period's price, the value
!      of land deducted to get the capital stock must also be valued at last
!      period's price.

     WRITE(10,*) 'Capital-Output Ratio', K1/OUTPUT
     WRITE(10,*)
     WRITE(10,*) 'Bequests (received) =', BEQ1
     WRITE(10,*)
     
     WRITE(10,*) 'Expected lifetime utility (value function) =', EXPCTDUT
!      WRITE(10,*) 'Average lifetime utility (invariant distribution) =', AVGUTIL

	 !	Print out the life cycle profiles

	DO AGE=1,MAXAGE
	WRITE(18,88) AGE, ALONG(AGE), CLONG(AGE), ILONG(AGE)
88	FORMAT (1X,I2,1X,3(F10.7,X))
	END DO
     




 END DO	! theta loop



!****************
!3 subroutines
CONTAINS


! 3.3 search five
SUBROUTINE SRCHFIVE01
                                          

DO JA=IL,IU,ISKIP
    CON2 = 0
   CONS = X3 - A(JA)
   EVAL:  IF (CONS>=CMIN) THEN
      CON1 = (CONS - CMIN)/CINCR + 1.0
      CON2 = FLOOR(CON1)                             !  Integer
      CON3 = CON1 - FLOAT(CON2)                        !  Real
      UTIL = (1.0-CON3)*UT(CON2) + CON3*UT(CON2+1)

      IF (AGE<RETAGE-1) THEN
         VTEMP = UTIL + BETA*S(AGE+1)*VW(AGE+1,JA)
      ELSE
         VTEMP = UTIL + BETA*S(AGE+1)*VR(AGE+1,JA)
      END IF
      
      UPDATE:  IF (VTEMP>=VMAX) THEN
         VMAX = VTEMP
         JAMAX = JA
      END IF UPDATE

   END IF EVAL

END DO

END SUBROUTINE

! 3.2 bracket

SUBROUTINE BRACKET01



     X3 = ATR*( A(IA) + BEQ*(1.0+GROWTH)**(AGE-1) ) + WAGEZ
     VMAX = -1.E6
     JAMAX = 1
     IL = 1
     IU = NGRID
     ISKIP = (NGRID - 1)/4
101  CALL SRCHFIVE01                   !  Updates VMAX and JAMAX
     NARROW01:  IF (ISKIP>1) THEN

        NARROW02:  IF ((JAMAX>1).AND.(JAMAX<NGRID)) THEN
           IL = JAMAX - ISKIP
           IU = JAMAX + ISKIP
           ISKIP = ISKIP/2
        ELSE IF (JAMAX==1) THEN
           IF (ISKIP>=4) THEN
              ISKIP = ISKIP/4
              IU = IL + 4*ISKIP
           ELSE IF (ISKIP==2) THEN
              ISKIP = 1
              IU = 2
           ELSE
              PRINT *, 'Error in Subroutine BRACKET at NARROW02'
           END IF
        ELSE
           IF (ISKIP>=4) THEN
              ISKIP = ISKIP/4
              IL = IU - 4*ISKIP
           ELSE IF (ISKIP==2) THEN
              ISKIP = 1
              IL = NGRID - 1
           ELSE
              PRINT *, 'Error in Subroutine BRACKET at NARROW02'
           END IF
        END IF NARROW02

          IF (IL<1) PRINT *, 'Error:  IL<1 in Subroutine BRACKET'
          IF (IL>NGRID) PRINT *, 'Error:  IL>NGRID in Subroutine BRACKET'

        GO TO 101

     END IF NARROW01

     IF (AGE<RETAGE) THEN
        VW(AGE,IA) = VMAX
        IDCW(AGE,IA) = JAMAX
     ELSE
        VR(AGE,IA) = VMAX
        IDCR(AGE,IA) = JAMAX
     END IF

END SUBROUTINE


! 3.1 decrule
SUBROUTINE DECRULE01



!   Initialize value function and decision rules

     DO AGE=1,RETAGE-1
        DO IA=1,NGRID
              VW(AGE,IA) = -10000.
              IDCW(AGE,IA) = -1
        END DO
     END DO
     DO AGE=RETAGE,MAXAGE
        DO IA=1,NGRID
           VR(AGE,IA) = -10000.
           IDCR(AGE,IA) = -1
        END DO
     END DO
!   Value function and asset choice for last period

     DO IA=1,NGRID
        CONS = ATR*( A(IA) + BEQ*(1.0+GROWTH)**(MAXAGE-1) ) + SS
        VR(MAXAGE,IA) = (CONS**(1.0-GAMMA))/(1.0-GAMMA)
        IDCR(MAXAGE,IA) = 1
     END DO
!   Remaining retirees
     DO AGE=MAXAGE-1,RETAGE,-1

        WAGEZ = SS      
        DO IA=1,NGRID
           CALL BRACKET01  ! Finds optimal asset choice for this state
        END DO

     END DO

!   Working-age agents

     DO AGE=RETAGE-1,1,-1

        WAGEZ = (1.0-STAX)*WAGE*HBAR*EFFLONG(AGE)
        DO IA=1,NGRID   
           CALL BRACKET01
        END DO

     END DO

     END SUBROUTINE



! 3.4 profile
SUBROUTINE PROFILE01



!   Compute longitudinal profiles for a given cohort,
!      and average lifetime utility

	 ALONG(:) = 0.0
	 CLONG(:) = 0.0
	 ILONG(:) = 0.0
	 
	 DLONG(1) = IDCW(1,1)
	 DO AGE=2,MAXAGE
		 IF (AGE<RETAGE) DLONG(AGE) = IDCW(AGE,DLONG(AGE-1))
		 IF (AGE>=RETAGE) DLONG(AGE) = IDCR(AGE,DLONG(AGE-1))
	 END DO

	 
	 DO AGE=1,MAXAGE
		 ALONG(AGE) = A( DLONG(AGE) )
		 IF (AGE<RETAGE) ILONG(AGE) = (1.0-STAX)*WAGE*HBAR*EFFLONG(AGE)
		 IF (AGE>=RETAGE) ILONG(AGE) = SS 
     END DO
     CLONG(1) = ATR*BEQ + ILONG(1) - ALONG(1)
	 DO AGE=2,MAXAGE
         CLONG(AGE) = ATR*( ALONG(AGE-1) + BEQ*(1.0+GROWTH)**(AGE-1) ) + ILONG(AGE) - ALONG(AGE)
     END DO


!   Compute cross-sectional profiles for a given time period

     DO AGE=1,MAXAGE

        ACROSS(AGE) = ALONG(AGE)*(1+GROWTH)**(1-AGE)
        CCROSS(AGE) = CLONG(AGE)*(1+GROWTH)**(1-AGE)
        ICROSS(AGE) = ILONG(AGE)*(1+GROWTH)**(1-AGE)
 
     END DO
 
     END SUBROUTINE






 999 END PROGRAM
