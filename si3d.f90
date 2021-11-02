 
!*************************************************************************
                               PROGRAM si3d
!*************************************************************************
!
!  Purpose:  Main program for the semi-implicit 3-d hydrodynamic model.
!            The model uses a leapfrog-trapezoidal finite-difference
!            numerical scheme.  
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!-------------------------------------------------------------------------
!                         U. S. Geological Survey
!                             Sacramento, CA
!-------------------------------------------------------------------------

   USE si3d_procedures
   IMPLICIT NONE

   !.....External function declaration.....
   REAL, EXTERNAL :: TIMER             ! Use timing routine from NSPCG

   !.....Local variables.....
   REAL :: utime, stime, utime1, stime1, utime2, stime2, utime3, stime3,  &
         & ttime1, ttime2, atime, ltime1, ltime2
   REAL :: TimeStart, TimeEnd
   INTEGER :: maxcount = 1E6
   INTEGER :: iter, itemp, i, j, niter1
 
   !.....Retrieve begin time of run.....
   stime = TIMER(0.0)
 
   !.....Read input parameters.....
   CALL input

   ! ... Output run parameters for checking purposes 
   IF (idbg == 1) PRINT *, " Before entry to SUB outr"
   CALL outr

   !.....Allocate space for arrays.....
   IF (idbg == 1) PRINT *, " Before entry to SUB allocate_space"
   CALL AllocateSpace  

   !.....Read bathymetry file and setup logical mask arrays.....
   IF (idbg == 1) PRINT *, " Before entry to SUB bathy"
   CALL bathy

   !.....Define initial conditions.....
   IF (idbg == 1) PRINT *, " Before entry to SUB init"
   CALL init
   
   !.....Open files with lateral BC data and assign values at t=0..... 
   IF (idbg == 1) PRINT *, " Before entry to SUB openbc0"
   CALL openbc0

   ! ... Open files with surface BC data and assign values at t=0.....
   IF (idbg == 1) PRINT *, " Before entry to SUB surfbc0"
   CALL surfbc0 

   !.....Output initial conditions.....
   IF (idbg == 1) PRINT *, " Before entry to SUB outt"
   IF (nnodes> 0) CALL outt   ! Point profiles
   IF (idbg == 1) PRINT *, " Before entry to SUB outv"
   IF (iox   > 0) CALL outv   ! Cross-sections
   IF (idbg == 1) PRINT *, " Before entry to SUB outh"
   IF (iop   > 0) CALL outh   ! Horiz-planes
   IF (idbg == 1) PRINT *, " Before entry to SUB outw"
   IF (iop   > 0) CALL outw   ! Wind field
   IF (idbg == 1) PRINT *, " Before entry to SUB outz"
   IF (iotr  > 0) CALL outz   ! Tracers in 3D domain
   IF (idbg == 1) PRINT *, " Before entry to SUB outs"
   IF (ipxml > 0) CALL outs   ! Three-dimensional xml outputs
   IF (idbg == 1) PRINT *, " Before entry to SUB outp"
   IF (ipxml < 0) CALL outp   ! Three-dimensional outputs for ptrack
   IF (idbg == 1) PRINT *, " Before entry to SUB outnb"
   IF (ioNBO > 0) CALL outNB
   IF (idbg == 1) PRINT *, " Before entry to SUB outscalarbalance"
   IF (iobal > 0) CALL OutScalarEnergyBalance

   ltime1 = TIMER(0.0)                ! Retrieve begin time of loop

   !.....Integrate over time.....
   IF (idbg == 1) PRINT *, "Before entering loop over time"
   DO n = 1, nts
      TimeStart = TIMER(0.0)
      its = its + idt
      thrs = its/3600.
      CALL compute_date (idt)
      niter1 = niter
      lastiter = 0
      IF ( n == 1 ) THEN
         ! Use at least two iterations to start computations
         IF ( niter <  2 ) niter1 = 2
         IF ( niter >= 2 ) niter1 = niter
         GO TO 1
      END IF

      !.....Solve leapfrog step.....
      istep = 1
      IF (idbg == 1) PRINT *, " Before entry into SUB fd in leapfrog step"
      IF((itrap == 0) .OR. (MOD(n,MAX(itrap,1)) /= 0)) lastiter = 1 
      CALL fd
      IF((itrap == 0) .OR. (MOD(n,MAX(itrap,1)) /= 0)) GO TO 2 
 
      !.....Solve trapezoidal step.....
      CALL settrap
    1 istep = 2
      DO iter = 1, niter1
         IF ( iter == niter1 ) lastiter = 1; 
         IF (idbg == 1) PRINT *, " Before entry into SUB fd in trap step"
         CALL fd
         IF (idbg == 1) PRINT *, " After exit from SUB fd in trapezoidal step"
         IF(iter < niter1) CALL settrap2
         IF(iter > 20) THEN
            PRINT *, " ERROR--Too many iterations requested"
            EXIT
         END IF
      END DO

      !.....Save information for next time step.....
    2 CALL save

      !.....Output results.....
    3 IF((nnodes > 0) .AND. (MOD(n,MAX(ipt,  1)) == 0)) CALL outt
      IF (idbg == 1) PRINT *, " After entry to SUB outt"
      IF((iox    > 0) .AND. (MOD(n,MAX(iox,  1)) == 0)) CALL outv
      IF (idbg == 1) PRINT *, " After entry to SUB outv"
      IF((iop    > 0) .AND. (MOD(n,MAX(iop,  1)) == 0)) CALL outh
      IF (idbg == 1) PRINT *, " After entry to SUB outh"
      IF((iop    > 0) .AND. (MOD(n,MAX(iop,  1)) == 0)) CALL outw
      IF (idbg == 1) PRINT *, " After entry to SUB outw"
      IF((iotr   > 0) .AND. (MOD(n,MAX(iotr, 1)) == 0)) CALL outz
      IF (idbg == 1) PRINT *, " After entry to SUB outz"
      IF((ipxml  > 0) .AND. (MOD(n,MAX(ipxml,1)) == 0)) CALL outs
      IF (idbg == 1) PRINT *, " After entry to SUB outs"
      IF((ipxml  < 0) .AND. (MOD(n,MAX(apxml,1)) == 0)) CALL outp
      !IF((ioNBO  > 0) .AND. (MOD(n,MAX(ioNBO,1)) == 0)) CALL outNB ! MassBalanceCheck
      IF (idbg == 1) PRINT *, " After entry to SUB outp"
	  CALL outNB
      IF (idbg == 1) PRINT *, " After entry to SUB outNB"
      IF(iobal   > 0  .AND. (MOD(n,MAX(ipt  ,1)) == 0)) CALL OutScalarEnergyBalance
      IF (idbg == 1) PRINT *, " After entry to SUB outenergy"

      !.....Write to log file and check if job should be stopped.....
      CALL outlog
      CALL check_stopfile

      !.....End loop over time.....
      IF(n > maxcount) THEN
         PRINT *, " ERROR--A maximum of 1 million time steps is allowed"
         EXIT
      END IF
      TimeEnd = TIMER(0.0)

      PRINT *, 'Time spent in step ', n, ' = ', TimeEnd - TimeStart, ' seconds'

   END DO
   ltime2 = TIMER(0.0)
   PRINT '(A,F10.3,A)', " Loop time =", ltime2-ltime1, " seconds"

   !.....Output CPU times.....
   CALL cputimes
   utime = TIMER(0.0)
   
   !.....Print total time of run.....
   PRINT '(A, F10.3, A)',  &
   & " Total CPU Time =", utime-stime,         " seconds"
   
   !.....End program.....
   PRINT '(A)', " "
   PRINT *, " *****PROGRAM TERMINATED NORMALLY"

END PROGRAM si3d


