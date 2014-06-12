PROGRAM LAND
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
PARAMETER (BEQ0 = 0.0420)         !  Accidental bequests

!******************
!
!   Set Parameters
!
!******************

PARAMETER (TOLK = 0.001)       !  Convergence tolerance for capital stock
PARAMETER (TOLB = 0.001)       !  Convergence tolerance for bequests
PARAMETER (GRADK = 0.2)        !  Convergence gradient for capital stock
PARAMETER (GRADB = 0.6)        !  Convergence gradient for bequests
PARAMETER (MAXITER = 50)       !  Maximum number of iterations for convergence

PARAMETER (ALPHA  = 0.690)     !  Labor exponent in production function
PARAMETER (TFP    = 1.005)     !  Multiplicative constant in production function
PARAMETER (GROWTH = 0.0165)    !  Growth rate of per capita output
PARAMETER (DEP    = 0.069)     !  Depreciation rate

PARAMETER (BETA   = 0.978)      !  Subjective discount factor
PARAMETER (GAMMA  = 2.0)       !  Risk aversion parameter

!PARAMETER (THETA  = 0.40)      !  Social security replacement rate
PARAMETER (PHI    = 0.30)      !  Unemployment insurance replacement ratio

PARAMETER (MAXAGE = 65)        !  Maximum age allowed
INTEGER, PARAMETER :: RETAGE = 45        !  Retirement age  
PARAMETER (HBAR   = 1.00)      !  Exogenous hours of work
PARAMETER (RHO    = 0.012)     !  Population growth rate

PARAMETER (AMAX   = 40.0)      !  Maximum permissible asset
PARAMETER (NGRID  = 4097)      !  Number of points on asset grid

PARAMETER (KLOW   = 10)        !  10 times smallest initial capital in fixed-point diagram
PARAMETER (KHIGH  = 60)        !  10 times largest initial capital in fixed-point diagram

PARAMETER (CMIN   = 0.00005)   !  Minimum permissible consumption
!PARAMETER (CMAX   = 12.0)      !  Maximum permissible consumption
! below, we use JMAX = 1.5 * AMAX
PARAMETER (CINCR  = 0.001)     !  Consumption increment for utility tabulation

!********************
!
!   Data Type Declarations and Dimension Statements
!
!********************

INTEGER AGE, IL, IS, ISKIP, IU, JAMAX

REAL  ATR, BEQ1, BEQEND, EXDEM, EXPCTDUT, INT, INV, K, K1
REAL  KDEV, L, OUTPUT, SS, STAX, VMAX, WAGE, WAGEZ, X3

INTEGER IDCR(:,:)      !  Asset decision rules for retirees
INTEGER IDCW(:,:,:)    !  Asset decision rules for working-age agents

REAL A(:)              !  Asset levels
REAL ACROSS(:)         !  Cross-sectional age-assets profile
REAL ALONG(:)          !  Longitudinal age-assets profile
REAL CCROSS(:)         !  Cross-sectional age-consumption profile
REAL CLONG(:)          !  Longitudinal age-consumption profile
REAL CUMS(:)           !  Unconditional survival probabilities, age 1 to age j
REAL EFFCROSS(:)       !  Cross-sectional age-earnings profile
REAL EFFLONG(:)        !  Longitudinal age-earnings profile for given cohort
REAL ICROSS(:)         !  Cross-sectional age-income profile
REAL ILONG(:)          !  Longitudinal age-income profile
REAL MU(:)             !  Age distribution of population
REAL P(2,2)            !  Employment state transition probabilities
REAL S(:)              !  Conditional survival probabilities, age j-1 to age j
REAL UT(:)             !
   !Values in utility function lookup table
REAL VR(:,:)           !  Value function for retirees
REAL VW(:,:,:)         !  Value function for working-age agents
REAL YR(:,:)           !  Age-dependent distribution of retirees across asset
                       !     states
REAL YW(:,:,:)         !  Age-dependent distribution of working-age agents
                       !     across asset and employment states

ALLOCATABLE  A, ACROSS, ALONG, CCROSS, CLONG, CUMS, EFFCROSS, EFFLONG, ICROSS
ALLOCATABLE  IDCR, IDCW, ILONG, MU, S, UT, VR, VW, YR, YW

!*********************
!
!   Open Files
!
!*********************

OPEN(UNIT=7,FILE='/home/sean/Desktop/Dropbox/cultural_ss/comboeff.txt')
OPEN(UNIT=8,FILE='/home/sean/Desktop/Dropbox/cultural_ss/psurv.txt')
OPEN(UNIT=10,FILE='/home/sean/Desktop/Dropbox/cultural_ss/result.txt')
OPEN(UNIT=18,FILE='/home/sean/Desktop/Dropbox/cultural_ss/profile.txt')

!*********************
!
!   Read Data
!
!*********************

     ALLOCATE ( EFFCROSS(RETAGE-1), S(MAXAGE) )
     READ(7,*) ( EFFCROSS(AGE), AGE=1,RETAGE-1 )
     READ(8,*) ( S(AGE), AGE=1,MAXAGE )
     P = RESHAPE((/0.94, 0.94, 0.06, 0.06/),(/2,2/))

!*********************
!
!   Preliminary Calculations
!
!*********************

     K = K0
     BEQ = BEQ0

!   Tabulate asset levels

     ALLOCATE ( A(NGRID) )
     A = (/ ( (FLOAT(IA-1)*AMAX)/FLOAT(NGRID-1), IA=1,NGRID ) /)

!   Tabulate utility function

     JMAX = 1.5*AMAX/CINCR
     ALLOCATE ( UT(JMAX) )
     DO I=1,JMAX
        CONS = CMIN + (I-1)*CINCR
        UT(I) = (CONS**(1.0-GAMMA))/(1.0-GAMMA)
     END DO

!  Unconditional survival probabilities

     ALLOCATE ( CUMS(MAXAGE), MU(MAXAGE) )

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
        L = L + 0.94*HBAR*EFFCROSS(AGE)*MU(AGE)
     END DO
 
!  INSERT BLOCK C HERE TO COMPUTE A TABLE

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$
!$  Beginning of BLOCK C
!$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!   Supplemental code for creating a table by varying a single parameter value

!   Also requires BLOCK D

!   To use 1025-point grid,
!      (1)  comment out the relevant parameter statement above
!      (2)  set up DO-loop in terms of the parameter value
!      (3)  implicitly let the ending values of K and BEQ for one row of the
!           table be the initial values for the next row
!      (4)  write parameter value, K, and BEQ to unit 11
!      (5)  possibly write other output to unit 12

!   To use 4097-point grid,
!      (1)  comment out the relevant parameter statement above
!      (2)  set up DO-loop in terms of NROW index
!      (3)  read parameter value, K, and BEQ for each row of the table
!           from data in unit 11
!      (4)  write output to unit 12

     OPEN(UNIT=11,FILE='k_beq.out')
     OPEN(UNIT=12,FILE='table.out')

 DO ITHETA = 4,4,1                            ! 1025-point grid
     THETA = FLOAT(ITHETA)/10                     ! 1025-point grid
     PRINT *, 'THETA =', THETA                      ! 1025-point grid

! DO IROW = 1,39                                    ! 4097-point grid
!     READ(11,*) BETA, K, BEQ                      ! 4097-point grid
!     PRINT *, 'ROW', IROW, '      BETA =', BETA   ! 4097-point grid

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$
!$  End of BLOCK C
!$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!   Unemployment insurance tax rate

     UTAX = (0.06/0.94)*PHI

!   Growth rate of aggregate output

     AGROWTH = (1.0+GROWTH)*(1.0+RHO) - 1.0

!   Longitudinal age-earnings profile for given cohort

     ALLOCATE ( EFFLONG(RETAGE-1) )
     EFFLONG = (/ ( ((1+GROWTH)**(AGE-1))*EFFCROSS(AGE), AGE=1,RETAGE-1 ) /)

     WRITE(10,*) 'PARAMETER VALUES'
     WRITE(10,*)

     WRITE(10,*) '  Convergence tolerance for capital stock:', TOLK
     WRITE(10,*) '  Convergence tolerance for bequests:', TOLB
     WRITE(10,*) '  Convergence gradient for capital stock:', GRADK
     WRITE(10,*) '  Convergence gradient for bequests:', GRADB
     WRITE(10,*) '  Maximum number of iterations for convergence:', MAXITER
     WRITE(10,*)

     WRITE(10,*) '  Labor exponent in production function:', ALPHA
     WRITE(10,*) '  Multiplicative constant in production function:', TFP
     WRITE(10,*) '  Growth rate of per capita output:', GROWTH
     WRITE(10,*) '  Depreciation rate:', DEP
     WRITE(10,*)

     WRITE(10,*) '  Subjective discount factor:', BETA
     WRITE(10,*) '  Risk aversion parameter:', GAMMA
     WRITE(10,*)

     WRITE(10,*) '  Social security replacement ratio:', THETA
     WRITE(10,*) '  Unemployment insurance replacement rate:', PHI
     WRITE(10,*)

     WRITE(10,*) '  Maximum age allowed:', MAXAGE
     WRITE(10,*) '  Retirement age:', RETAGE
     WRITE(10,*) '  Exogenous hours of work:', HBAR
     WRITE(10,*) '  Population growth rate:', RHO
     WRITE(10,*)

     WRITE(10,*) '  Maximum permissible asset:', AMAX
     WRITE(10,*) '  Number of points on asset grid:', NGRID
     WRITE(10,*)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$
!$  Beginning of BLOCK A
!$
!$  Choose appropriate version from SOLVE.SUP or FIXEDPT.SUP
!$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!*********************
!
!   Iterate to Convergence
!
!*********************

     ITER = 1

     WRITE(10,*) 'ITERATION RESULTS'
     WRITE(10,*)
 200 ITERINC = 0

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$
!$  End of BLOCK A
!$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!
!   More Preliminary Calculations
!

!   Factor prices

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
     STAX = STAX/(0.94*WAGE*HBAR*CUM)             !  this period's workers
            
!
!   Main Calculations
!

!   Find decision rules for all ages and states

     ALLOCATE ( IDCR(RETAGE:MAXAGE,NGRID), IDCW(RETAGE-1,NGRID,2))
     ALLOCATE ( VR(RETAGE:MAXAGE,NGRID), VW(RETAGE-1,NGRID,2) )
     CALL DECRULE01                     

!   Calculate expected lifetime utility

     EXPCTDUT = 0.94*VW(1,1,1) + 0.06*VW(1,1,2)
     DEALLOCATE ( VR, VW )

!   Find invariant distribution

     ALLOCATE ( YR(RETAGE:MAXAGE,NGRID), YW(RETAGE-1,NGRID,2) )
     CALL INVAR01

     DO AGE=1,RETAGE-1
        DO IA=1,NGRID
           IF ((IDCW(AGE,IA,1)>=NGRID).AND.(YW(AGE,IA,1)>0)) THEN
              PRINT *, 'Error:  Maximum asset limit is binding.'
              PRINT *, 'AGE =', AGE
              PRINT *, 'IA =', IA
              WRITE(10,*) 'AGE =', AGE
              WRITE(10,*) 'IA =', IA
              GO TO 999
           END IF
        END DO
     END DO

!   Compute age profiles and average lifetime utility

     ALLOCATE ( ACROSS(MAXAGE), CCROSS(MAXAGE), ICROSS(MAXAGE) )
     ALLOCATE ( ALONG(MAXAGE), CLONG(MAXAGE), ILONG(MAXAGE) )
     CALL PROFILE01

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$
!$  Beginning of BLOCK B
!$
!$  Choose appropriate version from SOLVE.SUP or FIXEDPT.SUP
!$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!   Compute average end-of-period assets and bequests

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

     WRITE(10,*) 'Iteration', ITER
     PRINT *, 'Iteration', ITER
     WRITE(10,*) 'Initial capital =', K
     PRINT *, 'Initial capital =', K
     WRITE(10,*) 'Ending capital =', K1
     PRINT *, 'Ending capital =', K1
     WRITE(10,*) 'Relative change in capital =', KDEV
     PRINT *, 'Relative change in capital =', KDEV
     WRITE(10,*) 'Initial bequests =', BEQ
     PRINT *, 'Initial bequests =', BEQ
     WRITE(10,*) 'Ending bequest =', BEQ1
     PRINT *, 'Ending bequests =', BEQ1
     WRITE(10,*) 'Relative change in bequests =', BEQDEV
     PRINT *, 'Relative change in bequests =', BEQDEV
     WRITE(10,*)
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
           WRITE(10,*) 'Maximum number of iterations exceeded.'
           WRITE(10,*) 'Program terminates.'
           GO TO 799
        END IF
        DEALLOCATE (IDCR, IDCW, YR, YW)
        DEALLOCATE (ACROSS, ALONG, CCROSS, CLONG, ICROSS, ILONG)
        GO TO 200

     END IF

!*********************
!
!   Summary Calculations
!
!*********************

!   Compute average beginning-of-period assets

!   ABEG is assets per period-t effective labor unit with the land component
!      still valued at its price from period t-1.  The capital gain on land
!      occurs after ABEG is measured, and is included in the rate of return
!      in the agent's budget constraint.

 799 ABEG = 0.0

     DO AGE=2,RETAGE-1
        DO IA=1,NGRID
           DO IS=1,2
              ABEG = ABEG + A(IA)*YW(AGE,IA,IS)*MU(AGE)*((1+GROWTH)**(1-AGE))
           END DO
        END DO
     END DO

     DO AGE=RETAGE,MAXAGE
        DO IA = 1,NGRID
           ABEG = ABEG + A(IA)*YR(AGE,IA)*MU(AGE)*((1+GROWTH)**(1-AGE))
        END DO
     END DO

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
     WRITE(10,*) 'Unemployment insurance tax rate =', UTAX
     WRITE(10,*)

     WRITE(10,*) 'Wage rate =', WAGE
     WRITE(10,*) 'Interest rate =', INT
     WRITE(10,*)

     WRITE(10,*) 'Social security benefit =', SS
     WRITE(10,*) 'Social security tax rate =', STAX
     WRITE(10,*)

     WRITE(10,*) 'Output =', OUTPUT
     WRITE(10,*) 'Consumption =', ACONS
     WRITE(10,*) 'Investment =', INV
     WRITE(10,*) 'Excess demand =', EXDEM
     WRITE(10,*)

     WRITE(10,*) 'Capital (from ending assets) =', K1
     WRITE(10,*) 'Capital (from beginning assets) =', ABEG+BEQ1/(1.0+AGROWTH)
     WRITE(10,*)

!   Because ABEG and BEQ1 have land valued at last period's price, the value
!      of land deducted to get the capital stock must also be valued at last
!      period's price.

     WRITE(10,*) 'Capital-Output Ratio', K1/OUTPUT
     WRITE(10,*)
     WRITE(10,*) 'Bequests (received) =', BEQ1
     WRITE(10,*)
     
     WRITE(10,*) 'Expected lifetime utility (value function) =', EXPCTDUT
     WRITE(10,*) 'Average lifetime utility (invariant distribution) =', AVGUTIL

	 !	Print out the life cycle profiles

	DO AGE=1,MAXAGE
	WRITE(18,88) AGE, ALONG(AGE), CLONG(AGE), ILONG(AGE)
88	FORMAT (1X,I2,1X,3(F10.7,X))
	END DO
     
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$
!$  End of BLOCK B
!$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!  INSERT BLOCK D HERE TO COMPUTE A TABLE

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$
!$  Beginning of BLOCK D
!$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!   Supplemental code for creating a table by varying a single parameter value

!   Also requires BLOCK C

     K = K1                                       ! 1025-point grid
     BEQ = BEQ1                                   ! 1025-point grid
     WRITE(11,*) THETA, K1, BEQ1                   ! 1025-point grid
     WRITE(12,*) THETA, OUTPUT, K1, ACONS, INT, WAGE, EXPCTDUT
     WRITE(10,*)
     WRITE(10,*)

     DEALLOCATE ( EFFLONG )
     DEALLOCATE ( IDCR, IDCW, YR, YW )
     DEALLOCATE ( ACROSS, ALONG, CCROSS, CLONG, ICROSS, ILONG )

 END DO

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$
!$  End of BLOCK D
!$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!*********************
!
!   Internal Functions and Subroutines
!
!*********************
     CONTAINS
!******************************************************************

SUBROUTINE SRCHFIVE01

!******************
!
!   Finds optimum asset choice over several (generally five)
!         possibly non-contiguous grid points
!
!******************
!
!   Uses:    Parameters:   BETA, CINCR, CMIN, RETAGE
!            Variables:    AGE, IL, ISKIP, IU, VMAX, X3
!            Arrays:       A(:), P(:,:), S(:), UT(:), VR(:,:), VW(:,:,:)
!
!   Returns:   JAMAX, VMAX
!
!   Local:     CONS, DC, JA, JC, UTIL, VTEMP, XC
!
!   Calls:     None
!
!******************
                                                                  

DO JA=IL,IU,ISKIP

   CONS = X3 - A(JA)
   EVAL:  IF (CONS>=CMIN) THEN
      XC = (CONS - CMIN)/CINCR + 1.0
      JC = XC                             !  Integer
      DC = XC - JC                        !  Real
      UTIL = (1.0-DC)*UT(JC) + DC*UT(JC+1)

      IF (AGE<RETAGE-1) THEN
         VTEMP = UTIL + BETA*S(AGE+1)*( P(IS,1)*VW(AGE+1,JA,1) + &
                 P(IS,2)*VW(AGE+1,JA,2) )
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

!******************************************************************

SUBROUTINE BRACKET01

!******************
!
!   Finds global optimum asset choice for a single state
!
!******************
!
!   Uses:    Parameters:   NGRID, RETAGE
!            Variables:    AGE, ATR, BEQ, IA, IS, WAGEZ
!            Arrays:       None
!
!   Returns:   IDCR(:,:), IDCW(:,:,:), VR(:,:), VW(:,:,:)
!
!   Local:     IL, ISKIP, IU, JAMAX, VMAX, X3
!
!   Calls:     SRCHFIVE01
!
!******************

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
        VW(AGE,IA,IS) = VMAX
        IDCW(AGE,IA,IS) = JAMAX
     ELSE
        VR(AGE,IA) = VMAX
        IDCR(AGE,IA) = JAMAX
     END IF

END SUBROUTINE

!******************************************************************

SUBROUTINE DECRULE01

!******************
!
!   Finds optimal asset choice for all states
!
!******************
!
!   Uses:    Parameters:   GAMMA, HBAR, MAXAGE, NGRID, PHI, RETAGE,
!            Variables:    ATR, BEQ, SS, STAX, UTAX, WAGE
!            Arrays:       A(:), EFFLONG(:)
!
!   Returns:   IDCR(:,:), IDCW(:,:,:), VR(:,:), VW(:,:,:)
!
!   Local:     CONS, WAGEZ
!
!   Calls:     BRACKET01
!
!******************

!   Initialize value function and decision rules

     DO AGE=1,RETAGE-1
        DO IA=1,NGRID
           DO IS=1,2
              VW(AGE,IA,IS) = -10000.
              IDCW(AGE,IA,IS) = -1
           END DO
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

        IS = 1                           !  Employed
        WAGEZ = (1.0-STAX-UTAX)*WAGE*HBAR*EFFLONG(AGE)
        DO IA=1,NGRID   
           CALL BRACKET01
        END DO

        IS = 2                                  !  Unemployed
        WAGEZ = PHI*WAGE*HBAR*EFFLONG(AGE)  
        DO IA=1,NGRID
           CALL BRACKET01  ! Finds optimal asset choice for this state
        END DO

     END DO

     END SUBROUTINE

!******************************************************************

SUBROUTINE INVAR01

!******************
!
!   Finds invariant distribution
!
!******************
!
!   Uses:    Parameters:   MAXAGE, NGRID, RETAGE
!            Variables:    None
!            Arrays:       IDCR(:,:), IDCW(:,:,:), P(:,:)
!
!   Returns:   YR(:,:), YW(:,:,:)
!
!   Local:     JA
!
!   Calls:     None
!
!******************

!   Initialize age-dependent distributions

     DO AGE=1,RETAGE-1
        DO IA=1,NGRID
           DO IS=1,2
              YW(AGE,IA,IS)=0
           END DO
        END DO
     END DO
     DO AGE=RETAGE,MAXAGE
        DO IA=1,NGRID
           YR(AGE,IA)=0
        END DO
     END DO
     YW(1,1,1) = 0.94
     YW(1,1,2) = 0.06

!   Recursively compute Y(AGE,:,:) from Y(AGE-1,:,:) for
!      working-age agents
     DO AGE=2,RETAGE-1
        DO IA=1,NGRID

           DO IS=1,2
              JA = IDCW(AGE-1,IA,IS)
              YW(AGE,JA,1) = YW(AGE,JA,1) + YW(AGE-1,IA,IS)*P(IS,1)
              YW(AGE,JA,2) = YW(AGE,JA,2) + YW(AGE-1,IA,IS)*P(IS,2)

           END DO

        END DO

     END DO

!   New retirees

     DO IA=1,NGRID

        DO IS=1,2
           JA = IDCW(RETAGE-1,IA,IS)
           YR(RETAGE,JA) = YR(RETAGE,JA) + YW(RETAGE-1,IA,IS)
        END DO
     END DO

!   Previous retirees

     DO AGE=RETAGE+1,MAXAGE

        DO IA=1,NGRID
           JA = IDCR(AGE-1,IA)
           YR(AGE,JA) = YR(AGE,JA) + YR(AGE-1,IA)
        END DO

     END DO

 END SUBROUTINE

!******************************************************************

SUBROUTINE PROFILE01

!******************
!
!   Finds age profiles for income, consumption, and assets, and
!      computes average lifetime utility
!
!******************
!
!   Uses:    Parameters:   CINCR, CMIIN, HBAR, NGRID, PHI, RETAGE
!            Variables:    ATR, BEQ, SS, STAX, UTAX, WAGE
!            Arrays:       A(:), EFFLONG(:), IDCR(:,:), IDCW(:,:,:),
!                          UT(:), YR(:,:), YW(:,:,:)
!
!   Returns:   ACROSS(:), ALONG(:), AVGUTIL, CCROSS(:), CLONG(:),
!              ICROSS(:), ILONG(:)
!
!   Local:     CONS, DC, JA, JC, UTIL, WAGEZ, XC
!
!   Calls:     None
!
!******************

!   Compute longitudinal profiles for a given cohort,
!      and average lifetime utility

     AVGUTIL = 0.0
     DO AGE=1,RETAGE-1                 !  Working-age agents
        ALONG(AGE) = 0.0
        CLONG(AGE) = 0.0
        ILONG(AGE) = 0.0

        DO IA=1,NGRID
           DO IS=1,2

              !  Assets
              JA = IDCW(AGE,IA,IS)
              ALONG(AGE) = ALONG(AGE) + A(JA)*YW(AGE,IA,IS)

              !  Income
              IF(IS==1) THEN
                 WAGEZ = (1.0-STAX-UTAX)*WAGE*HBAR*EFFLONG(AGE)
              ELSE
                 WAGEZ = PHI*WAGE*HBAR*EFFLONG(AGE)
              END IF
              ILONG(AGE) = ILONG(AGE) + WAGEZ*YW(AGE,IA,IS)

              !  Consumption
              CONS = ATR*( A(IA) + BEQ*(1.0+GROWTH)**(AGE-1) ) + WAGEZ - A(JA)
              CLONG(AGE) = CLONG(AGE) + CONS*YW(AGE,IA,IS)
              !  Average lifetime utility
              XC = (CONS - CMIN)/CINCR + 1.0
              JC = XC
              DC = XC - JC
              UTIL = (1.0 - DC)*UT(JC) + DC*UT(JC+1)
              AVGUTIL = AVGUTIL + UTIL*YW(AGE,IA,IS)*CUMS(AGE)*(BETA**(AGE-1))

           END DO

        END DO

     END DO

     DO AGE=RETAGE,MAXAGE              !  Retirees

        ALONG(AGE) = 0.0
        CLONG(AGE) = 0.0
        ILONG(AGE) = SS

        DO IA=1,NGRID

           !  Assets
           JA = IDCR(AGE,IA)
           ALONG(AGE) = ALONG(AGE) + A(JA)*YR(AGE,IA)

           !  Consumption
           CONS = ATR*( A(IA) + BEQ*(1.0+GROWTH)**(AGE-1) ) + SS - A(JA)
           CLONG(AGE) = CLONG(AGE) + CONS*YR(AGE,IA)

           !  Average lifetime utility
           XC = (CONS - CMIN)/CINCR + 1.0
           JC = XC
           DC = XC - JC
           UTIL = (1.0 - DC)*UT(JC) + DC*UT(JC+1)
           AVGUTIL = AVGUTIL + UTIL*YR(AGE,IA)*CUMS(AGE)*(BETA**(AGE-1))

        END DO

     END DO

!   Compute cross-sectional profiles for a given time period

     DO AGE=1,MAXAGE

        ACROSS(AGE) = ALONG(AGE)*(1+GROWTH)**(1-AGE)
        CCROSS(AGE) = CLONG(AGE)*(1+GROWTH)**(1-AGE)
        ICROSS(AGE) = ILONG(AGE)*(1+GROWTH)**(1-AGE)
 
     END DO
 
     END SUBROUTINE

!*****************************************

 999 END PROGRAM
