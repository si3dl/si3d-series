!************************************************************************
  MODULE si3d_BoundaryConditions
!************************************************************************
!
!  Purpose: Define Boundary Conditions to use in the solution of NS eqns.
!           Includes openbc routines, originally included in si3d_proc and
!           surfbc_heat routines created to represent heat fluxes through
!           the free surface.
!
!-------------------------------------------------------------------------

  USE si3d_Types 

  IMPLICIT NONE
  SAVE

  CONTAINS

!************************************************************************
SUBROUTINE openbc0
!************************************************************************
!
!  Purpose: This routine is called at the beginning of the program
!           to open any files with boundary condition data, to
!           read the boundary condition data, and to 
!           assign the initial boundary values at time t=0.  
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL :: areatot, uavg, vavg, sbc, qbc, hbc
   INTEGER :: i, j, k, l, nn, is, ie, js, je, ks, ke, & 
              kmx, kmy, nwlayers, ios, istat
   INTEGER :: nbiid, noNGB, maxiptsNB, niNGB, icl, m1, m2
   INTEGER :: npts_wse, npts_flw, npts_sal, nptsOpenBC
   CHARACTER(LEN=14) :: openbcfmt, openbcfile
   REAL, ALLOCATABLE, DIMENSION(:,:) :: inputvar 

   ! ... Return if no open boundaries exist
   IF (nopen <= 0) RETURN

   !               -----Read files with bc data-----

   ! ---- Initialize No. of nested grid boundaries to zero
   noNGB = 0; maxiptsNB = -1;

   ! ----Loop over nopen to Open & read files with bc data-----  
   DO nn = 1, nopen

      SELECT CASE (itype(nn))

      CASE(1:3) ! Observed values at boundaries provided at txt files .....

      ! ... Construct file name (beware that nopen cannot be > 99)
      openbcfile = "openbc0 .txt"
      IF ( nn <  10 ) WRITE ( openbcfile(8:8), FMT='(I1)' ) nn
      IF ( nn >= 10 ) WRITE ( openbcfile(7:8), FMT='(I2)' ) nn

      ! ... Open IO unit
      OPEN (UNIT=i52, FILE=openbcfile, STATUS="old", IOSTAT=ios)
      IF (ios /= 0) CALL open_error ( "Error opening "//openbcfile, ios )

      ! Skip over first five header records in open boundary condition file
      READ (UNIT=i52, FMT='(////)', IOSTAT=ios)
      IF (ios /= 0) CALL input_error ( ios, 20 )

      ! Read number of points in file from seventh header record
      READ (UNIT=i52, FMT='(10X,I7)', IOSTAT=ios) nptsOpenBC
      IF (ios /= 0) CALL input_error ( ios, 21 )

      ! Allocate space for the array of data first time in loop
      IF ( nn == 1) THEN
        ALLOCATE ( varsOpenBC(nopen, ntr+2, nptsOpenBC), STAT=istat )
        IF (istat /= 0) CALL allocate_error ( istat, 16 )
      ENDIF

      ! Write the format of the data records into an internal file
      WRITE (UNIT=openbcfmt, FMT='("(10X,",I3,"G11.2)")') ntr+2

 
      ! Read data array
      DO j = 1, nptsOpenBC
         READ (UNIT=i52, FMT=openbcfmt, IOSTAT=ios) &
              (varsOpenBC(nn,i,j), i=1,ntr+2)
         IF (ios /= 0) CALL input_error ( ios, 22 )
      END DO

      ! ... Close IO unit
      CLOSE (i52)

      CASE(4:) ! Nested boundary conditions constructed in previous runs on coarser grid

        ! ... Open file
        openbcfile = "nbifilex0    "
        IF ( nn <  10) WRITE ( openbcfile(10:12), FMT='(I1,"  ")' ) nn
        IF ( nn >= 10) WRITE ( openbcfile( 9:12), FMT='(I2,"  ")' ) nn
        nbiid = nbiid0 + nn;
        OPEN(UNIT=nbiid,file=openbcfile,FORM='UNFORMATTED',IOSTAT=ios)
        IF (ios /= 0) CALL open_error ( "Error opening "//openbcfile, ios )

        !... Read type of boundary (needs to agree with input file)
        READ (nbiid) isdNBI(nn)
        !... Read time information (no. of frames) 
        READ (nbiid) nfrNBI(nn)
        !... Read spatial information (no. of cells)  
        READ (nbiid) iptNBI(nn), ntrNBI(nn)

        ! ... Check whether the length of simulations in the fine & coarse
        !     grids
        IF ( nfrNBI(nn) * dtsecOpenBC < tl) THEN
          PRINT * , 'Length of Coarse & Fine grid runs DISAGREE'
          PRINT * , 'Nested Boundary File no. = ', nn
          PRINT * , 'nframes in NB file = ', nfrNBI(nn)
          PRINT * , 'Time (s) between frames  = ', dtsecOpenBC
          PRINT * , 'Length of time simulated = ', tl 
          STOP
        ENDIF

        ! ... Increase the no. of nested grid boundaries by 1, update 
        !     max. no. of points within nested grid boundaries, and
        !     check for consistency between information from the coarse 
        !     and information required by the fine grids
        noNGB = noNGB + 1
        maxiptsNB = MAX(maxiptsNB, iptNBI(nn))
        IF ( ntrNBI(nn) .NE. ntr) THEN
          PRINT *, '******************************************************'
          PRINT *, 'STOP - No. of tracers different in coarse & fine grids'
          PRINT *, 'Please CHECK Coarse No. Tracers(',nn,')=', ntrNBI(nn)
          PRINT *, '             Fine   No. Tracers(',nn,')=', ntr
          PRINT *, '******************************************************'
          STOP
        ENDIF
        IF ( isdNBI(nn) .NE. iside(nn)) THEN
          PRINT *, '******************************************************'
          PRINT *, 'STOP - Boundary sides in coarse & fine grids dis-agree'
          PRINT *, 'Please CHECK Coarse Grid SIDE(',nn,')=', isdNBI(nn)
          PRINT *, '             Fine   Grid SIDE(',nn,')=', iside (nn)
          PRINT *, '******************************************************'
          STOP
        ENDIF

      END SELECT

   ENDDO

   ! -------------- Nested Grid Boundaries READ & CHECK -------------------
   IF ( noNGB > 0) THEN

     ! ... Allocate space for arrays holding nested grid boundaries 
     ALLOCATE ( uhNGB(maxiptsNB,noNGB), uhNGBp(maxiptsNB,noNGB), &
                vhNGB(maxiptsNB,noNGB), vhNGBp(maxiptsNB,noNGB), &
                scNGB(maxiptsNB,noNGB), scNGBp(maxiptsNB,noNGB), &
                STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 31 )
     ALLOCATE ( kNGB (maxiptsNB,noNGB), & 
                iNGB (maxiptsNB,noNGB), &
                jNGB (maxiptsNB,noNGB), STAT=istat)
     IF (istat /= 0) CALL allocate_error ( istat, 32 )

     ! ... Initialize arrays to zero
     uhNGB  = 0.0E0
     vhNGB  = 0.0E0
     scNGB  = 0.0E0
     uhNGBp = 0.0E0
     vhNGBp = 0.0E0
     scNGBp = 0.0E0
     kNGB   = 0

     ! ... Allocate space & initialize arrays if tracers are modelled
     IF (ntr > 0) THEN
       ALLOCATE ( trNGB (maxiptsNB,noNGB,ntr), &
                  trNGBp(maxiptsNB,noNGB,ntr), STAT=istat)
       IF (istat /= 0) CALL allocate_error ( istat, 33 )
       trNGB  = 0.0E0
       trNGBp = 0.0E0
     ENDIF

     ! Read frames 1 & 2 for nested grid boundaries
     niNGB = 0
     DO nn = 1, nopen

       IF ( itype(nn) < 4) CYCLE
       nbiid = nbiid0 + nn;
       niNGB = niNGB  + 1 ;

       ! ... Allocate space for temporary input variable array 
       ALLOCATE( inputvar ( iptNBI(nn), 5+ntr ), STAT=istat )
       IF (istat /= 0) CALL allocate_error (istat,34)

       ! ... Read FIRST FRAME & store variables ...............................
       READ(nbiid) thrsNGBp, &
           ((inputvar(m1,m2),m2=1,5+ntr),m1=1,iptNBI(nn))

       SELECT CASE (iside(nn))
       CASE(1,3)
         DO icl = 1, iptNBI(nn)
           iNGB  (icl,niNGB) = inputvar(icl,1)
           jNGB  (icl,niNGB) = inputvar(icl,2)
           kNGB  (icl,niNGB) = inputvar(icl,3)
           uhNGBp(icl,niNGB) = inputvar(icl,4)
           scNGBp(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGBp(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       CASE(2,4)
         DO icl = 1, iptNBI(nn)
           iNGB  (icl,niNGB) = inputvar(icl,1)
           jNGB  (icl,niNGB) = inputvar(icl,2)
           kNGB  (icl,niNGB) = inputvar(icl,3)
           vhNGBp(icl,niNGB) = inputvar(icl,4)
           scNGBp(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGBp(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       END SELECT

       ! ... Check grid size for consistency between coarse & fine grid .........
       SELECT CASE (iside(nn))
       CASE(1,3)
         i   = isbc(nn); 
         js  = jsbc(nn); 
         je  = jebc(nn);
         icl = 0
         DO j = js, je
           DO k = k1, kmz(i,j)
             icl = icl + 1      
             IF(kNGB(icl,niNGB) .NE. k) THEN
                PRINT *, '******************************************************'
                PRINT *, 'STOP - Coarse & fine grids numbering NOT CONSISTENT'
                PRINT *, 'Boundary No. ', nn
                PRINT *, 'Please CHECK Coarse Grid (    k) =', iNGB(icl,niNGB), &
                &                                              jNGB(icl,niNGB), &
                &                                              kNGB(icl,niNGB)
                PRINT *, '             Fine   Grid (i,j,k) =', i,j,k
                PRINT *, '******************************************************'
                STOP
             ENDIF
           ENDDO
         ENDDO
         IF (icl .NE. iptNBI(nn) ) THEN
           PRINT *, '******************************************************'
           PRINT *, 'STOP - Coarse & fine grids numbering NOT CONSISTENT'
           PRINT *, 'Boundary No. ', nn
           PRINT *, 'Cells in output nested boundary = ', iptNBI(nn)
           PRINT *, 'Cells in fine grid at  boundary = ', icl
           PRINT *, '******************************************************'
           STOP
         ENDIF

       CASE(2,4)
         j   = jsbc(nn); 
         is  = isbc(nn); 
         ie  = iebc(nn);
         icl = 0
         DO i = is, ie
           DO k = k1, kmz(i,j)
             icl = icl + 1      
             IF(kNGB(icl,niNGB) .NE. k) THEN
                PRINT *, '******************************************************'
                PRINT *, 'STOP - Coarse & fine grids numbering NOT CONSISTENT'
                PRINT *, 'Boundary No. ', nn
                PRINT *, 'Please CHECK Coarse Grid (i,j,k) =', iNGB(icl,niNGB), &
                &                                              jNGB(icl,niNGB), &
                &                                              kNGB(icl,niNGB)
                PRINT *, '             Fine   Grid (i,j,k) =', i,j,k
                PRINT *, '******************************************************'
                STOP
             ENDIF
           ENDDO
         ENDDO
         IF (icl .NE. iptNBI(nn) ) THEN
           PRINT *, '******************************************************'
           PRINT *, 'STOP - Coarse & fine grids numbering NOT CONSISTENT'
           PRINT *, 'Boundary No. ', nn
           PRINT *, 'Cells in output nested boundary = ', iptNBI(nn)
           PRINT *, 'Cells in fine grid at  boundary = ', icl
           PRINT *, '******************************************************'
           STOP
         ENDIF

       END SELECT      
 
       ! ... Read SECOND FRAME & store variables ...............................
       READ(nbiid) thrsNGB, & 
                 ((inputvar(m1,m2),m2=4,5+ntr),m1=1,iptNBI(nn))
       SELECT CASE (iside(nn))
       CASE(1,3)
         DO icl = 1, iptNBI(nn)
           uhNGB(icl,niNGB) = inputvar(icl,4)
           scNGB(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGB(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       CASE(2,4)
         DO icl = 1, iptNBI(nn)
           vhNGB(icl,niNGB) = inputvar(icl,4)
           scNGB(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGB(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       END SELECT

       DEALLOCATE (inputvar)

     ENDDO

   ENDIF         

   !           -----Assign boundary values at time t=0.0-----

   !.....Initialize to zero boundary flows & thickness at time n+1
   uhEB   = 0.0  ; huEB   = 0.0  ;
   uhWB   = 0.0  ; huWB   = 0.0  ;
   vhNB   = 0.0  ; hvNB   = 0.0  ;
   vhSB   = 0.0  ; hvSB   = 0.0  ;
   niNGB  = 0    ;

   !.....Loop over open boundaries..........................................
   DO nn = 1, nopen

      SELECT CASE ( itype(nn) )
 
      ! ..... CASE 1 -- wse specified ......................................
      CASE (1)

         ! Get first boundary value of zeta & make sure the i,j locations
         ! are wett boundary cells (Assume first value applies for t=0.0. 
         ! It should be consistent with the initial condition for zeta 
         ! defined in SUBROUTINE init)
         sbc = varsOpenBC(nn,1,1);
        
         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)

            ! ... Retrieve (i,j) index where BC is specified
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            ! ... Assign boundary condition
            sp(i,js:je) = sbc
            ! ... Make sure i,j locations are wett boundary cells
            DO j = js, je   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i-1,j) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - West bdry. not a bdry. point"
                STOP 
              END IF

              nwlayers = (kmz(i,j) - k1z(i,j)) + 1
              IF(nwlayers <= 0) THEN  ! dry point on boundary is not allowed
                PRINT *, " ERROR--dry point on west boundary" 
                PRINT *, "  "
                PRINT *, "  "
                PRINT *, " ****STOPPING si3d for boundary condition error"
                STOP 
              END IF
            END DO

         ! North boundary
         CASE (2)

            ! ... Retrieve (i,j) index where BC is specified
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! ... Assign boundary condition
            sp(is:ie,j) = sbc
            ! ... Make sure i,j locations are wett boundary cells
            DO i = is, ie   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j+1) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - North bdry. not a bdry. point"
                STOP 
              ENDIF
              nwlayers = (kmz(i,j) - k1z(i,j)) + 1
              IF (nwlayers < 1) THEN  ! dry point on boundary is not allowed
                PRINT *, " ERROR--dry point on north boundary" 
                PRINT *, "  "
                PRINT *, "  "
                PRINT *, " ****STOPPING si3d for boundary condition error"
                STOP 
              ENDIF
            END DO

         ! East boundary
         CASE (3)

            ! ... Retrieve (i,j) index where BC is specified
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            ! ... Assign boundary condition
            sp(i,js:je) = sbc
            ! ... Make sure i,j locations are wett boundary cells
            DO j = js, je   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i+1,j) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - East bdry. not a bdry. point"
                STOP 
              ENDIF
              nwlayers = (kmz(i,j) - k1z(i,j)) + 1
              IF ( nwlayers < 1) THEN ! dry point on boundary is not allowed
                PRINT *, " ERROR--dry point on east boundary" 
                PRINT *, "  "
                PRINT *, "  "
                PRINT *, " ****STOPPING si3d for boundary condition error"
                STOP 
              ENDIF
            ENDDO

         ! South boundary
         CASE (4)

            ! ... Retrieve (i,j) index where BC is specified
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! ... Assign boundary condition
            sp(is:ie,j) = sbc
            ! ... Make sure i,j locations are wett boundary cells
            DO i = is, ie 
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j-1) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - South bdry. not a bdry. point"
                STOP 
              ENDIF
              nwlayers = (kmz(i,j) - k1z(i,j)) + 1
              IF (nwlayers < 1) THEN ! dry point on boundary is not allowed
                PRINT *, " ERROR--dry point on south boundary" 
                PRINT *, "  "
                PRINT *, "  "
                PRINT *, " ****STOPPING si3d for boundary condition error"
                STOP 
              ENDIF
            END DO 

         END SELECT

      !..... CASE 2 -- Free surface flow specified .........................
      CASE (2)
      
         ! Get first boundary value of flow (in units of m**3/sec)
         ! (Assume first value applies for t=0.0. It should be consistent
         !  with the initial condition for uh, vh, u,and v)
         qbc = varsOpenBC(nn,1,1);

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
 
            ! Get i-, j- indexes for bdry. point
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)

            ! Make sure i,j locations are wett boundary cells
            DO j = js, je   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i-1,j) ) THEN
                 PRINT *, "  "
                 PRINT *, " ****STOPPING - West bdry. not a bdry. point"
                 STOP
              END IF
            ENDDO

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO j = js,je
              IF(hhs(i,j)>hbc) hbc = hhs(i,j)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(i,j), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; huWB(:,js:je) = ZERO
            DO j = js,je
              kmx = kmz(i,j)
              DO k = k1, kmx
                huWB (k,j)=AMIN1(zlevel(k+1),hhs(i,j)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(huWB(k,j) <= HMIN) huWB(k,j) = ZERO;
                areatot = areatot + huWB(k,j) * dy
              ENDDO
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; uhWB(:,js:je) = 0.0E0
            DO j = js, je
              kmx = kmz(i,j)
              DO k = k1, kmx
                uhWB(k,j) = uavg * huWB(k,j)
              END DO
            ENDDO

         ! East boundary
         CASE (3)
 
            ! ... Get i-, j- indexes for bdry. point
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            ! ... Make sure i,j locations are wett boundary cells
            DO j = js, je   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i+1,j) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - East bdry. not a bdry. point"
                STOP
              END IF
            ENDDO
            ! ... Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO j = js,je
              IF(hhs(i,j)>hbc) hbc = hhs(i,j)
            ENDDO
            ! ... Define free surface location at bdry.
            sbc = MAX(s(i,j), -hbc+dzmin)
            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; huEB(:,js:je) = ZERO
            DO j = js,je
              kmx = kmz(i,j)
              DO k = k1, kmx
                huEB (k,j)=AMIN1(zlevel(k+1),hhs(i,j)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(huEB(k,j) <= HMIN) huEB(k,j) = ZERO;
                areatot = areatot + huEB(k,j) * dy
              ENDDO
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; uhEB(:,js:je) = 0.0E0
            DO j = js, je
              kmx = kmz(i,j)
              DO k = k1, kmx
                 uhEB(k,j) = uavg * huEB(k,j)
              END DO
            ENDDO
            

         ! North boundary
         CASE (2)
 
            ! ... Get i-, j- indexes for bdry. point
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! ... Make sure i,j locations are wett boundary cells
            DO i = is, ie   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j+1) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - North bdry. not a bdry. point"
                STOP 
              END IF
            ENDDO
            ! ... Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO i = is,ie
              IF(hhs(i,j)>hbc) hbc = hhs(i,j)
            ENDDO
            ! ... Define free surface location at bdry.
            sbc = MAX(s(i,j), -hbc+dzmin)
            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; hvNB(:,is:ie) = ZERO
            DO i = is,ie
              kmy = kmz(i,j)
              DO k = k1, kmy
                hvNB (k,i)=AMIN1(zlevel(k+1),hhs(i,j)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(hvNB(k,i) <= HMIN) hvNB(k,i) = ZERO;
                areatot = areatot + hvNB(k,i) * dy
              ENDDO
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; vhNB(:,is:ie) = 0.0E0
            DO i = is, ie
              kmy = kmz(i,j)
              DO k = k1z(i,j), kmy
                 vhNB(k,i) = vavg * hvNB(k,i)
              END DO
            ENDDO

         ! South boundary
         CASE (4)
 
            ! ... Get i-, j- indexes for bdry. point
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! ... Make sure i,j locations are wett boundary cells
            DO i = is, ie   
              IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j-1) ) THEN
                PRINT *, "  "
                PRINT *, " ****STOPPING - South bdry. not a bdry. point"
                STOP 
              END IF
            ENDDO
            ! ... Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO i = is,ie
              IF(hhs(i,j)>hbc) hbc = hhs(i,j)
            ENDDO
            ! ... Define free surface location at bdry.
            sbc = MAX(s(i,j), -hbc+dzmin)
            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; hvSB(:,is:ie) = ZERO
            DO i = is,ie
              kmy = kmz(i,j)
              DO k = k1, kmy
                hvSB (k,i)=AMIN1(zlevel(k+1),hhs(i,j)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(hvSB(k,i) <= HMIN) hvSB(k,i) = ZERO;
                areatot = areatot + hvSB(k,i) * dy
              ENDDO
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; vhSB(:,is:ie) = 0.0E0
            DO i = is, ie
              kmy = kmz(i,j)
              DO k = k1, kmy
                 vhSB(k,i) = vavg * hvSB(k,i)
              END DO
            ENDDO

         END SELECT         

      !..... CASE 3 -- Submerged flow specified.............................
      CASE (3)
      
         ! Get first boundary value of flow (in units of m**3/sec)
         ! (Assume first value applies for t=0.0. It should be consistent
         !  with the initial condition for uh, vh, u,and v)
         qbc = varsOpenBC(nn,1,1);

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
 
            ! ... Get i-, j- indexes for bdry. point
            i  = isbc(nn); 
            j  = jsbc(nn); 
            l  = ij2l(i,j)
            ks = iebc(nn); 
            ke = jebc(nn); 
            ! ... Make sure i,j locations are wett boundary cells
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i-1,j) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - West bdry. not a bdry. point"
              STOP
            END IF
            ! ... Make sure k- location is wett, define thickness of bdry. cell
            !     & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              huWB (k,j)= hp(k,l) 
              IF ( huWB(k,j) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Spec. West bdry. below free surface"
                STOP
              ENDIF
              areatot = areatot + huWB(k,j) * dy
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; 
            DO k = ks, ke
              uhWB(k,j) = uavg * huWB(k,j)
            ENDDO

         ! East boundary
         CASE (3)
 
            ! ... Get i-, j- indexes for bdry. point
            i  = isbc(nn); 
            j  = jsbc(nn); 
            l  = ij2l(i,j)
            ks = iebc(nn); 
            ke = jebc(nn); 

            ! ... Make sure i,j locations are wett boundary cells
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i+1,j) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - East bdry. not a bdry. point"
              STOP
            END IF
            ! ... Make sure k- location is wett, define thickness of bdry. cell
            !     & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              huEB (k,j)= hp(k,l) 
              IF ( huEB(k,j) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - East bdry. below free surface"
                STOP
              ENDIF
              areatot = areatot + huEB(k,j) * dy
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; 
            DO k = ks, ke
              uhEB(k,j) = uavg * huEB(k,j)
            ENDDO

         ! North boundary
         CASE (2)
 
            ! ... Get i-, j- indexes for bdry. point
            j  = jsbc(nn); 
            i  = isbc(nn);
            l  = ij2l(i,j) 
            ks = iebc(nn);
            ke = jebc(nn);
            ! ... Make sure i,j locations are wett boundary cells
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j+1) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - North bdry. not a bdry. point"
              STOP 
            END IF
            ! ... Make sure k- location is wett, define thickness of bdry. cell
            !     & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              hvNB (k,i)= hp(k,l) 
              IF ( hvNB(k,i) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged North bdry. on a DRY CELL"
                STOP
              ENDIF
              areatot = areatot + hvNB(k,i) * dx
            ENDDO
            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; 
            DO k = ks, ke
              vhNB(k,i) = vavg * hvNB(k,i)
            ENDDO

         ! South boundary
         CASE (4)
 
            ! ... Get i-, j- indexes for bdry. point
            j  = jsbc(nn); 
            i  = isbc(nn);
            l  = ij2l(i,j) 
            ks = iebc(nn);
            ke = jebc(nn);
            ! ... Make sure i,j locations are wett boundary cells
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j-1) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - South bdry. not a bdry. point"
              STOP 
            END IF
            ! ... Make sure k- location is wett, define thickness of bdry. cell
            !     & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              hvSB (k,i)= hp(k,l) 
              IF ( hvSB(k,i) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged South bdry. on a DRY CELL"
                STOP
              ENDIF
              areatot = areatot + hvSB(k,i) * dx
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; 
            DO k = ks, ke
              vhSB(k,i) = vavg * hvSB(k,i)
            ENDDO

         END SELECT

      !..... CASE 4 -- Nested grid boundaries ..............................
      CASE (4)

        niNGB = niNGB + 1;
        SELECT CASE (iside(nn))
        CASE(1)
          i   = isbc(nn); 
          js  = jsbc(nn); 
          je  = jebc(nn);
          icl = 0
          DO j = js, je; 
            l   = ij2l(i,j)
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i-1,j) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - West bdry. not a bdry. point"
              STOP
            END IF
            DO k = k1, kmz(i,j)
              icl = icl + 1      
              uhWB(k,j) = uhNGBp(icl,niNGB)
              huWB(k,j) = hp(k,l)
            ENDDO 
          ENDDO 
        CASE(3)
          i   = isbc(nn); 
          js  = jsbc(nn); 
          je  = jebc(nn);
          icl = 0
          DO j = js, je; 
            l  = ij2l(i,j)
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i+1,j) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - East bdry. not a bdry. point"
              STOP
            END IF
            DO k = k1, kmz(i,j)
              icl = icl + 1      
              uhEB(k,j) = uhNGBp(icl,niNGB)
              huEB(k,j) = hp(k,l)
            ENDDO
          ENDDO
        CASE(2)
          j   = jsbc(nn); 
          is  = isbc(nn); 
          ie  = iebc(nn);
          icl = 0
          DO i = is, ie; 
            l  = ij2l(i,j)
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j+1) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - North bdry. not a bdry. point"
              STOP
            END IF
            DO k = k1, kmz(i,j)
              icl = icl + 1      
              vhNB(k,i) = vhNGBp(icl,niNGB)
              hvNB(k,i) = hp(k,l)
            ENDDO
          ENDDO
        CASE(4)
          j   = jsbc(nn); 
          is  = isbc(nn); 
          ie  = iebc(nn);
          icl = 0
          DO i = is, ie; 
            l  = ij2l(i,j)
            IF ( (.NOT. mask2d(i,j)) .OR. mask2d(i,j-1) ) THEN
              PRINT *, "  "
              PRINT *, " ****STOPPING - South bdry. not a bdry. point"
              STOP
            END IF
            DO k = k1, kmz(i,j)
              icl = icl + 1       
              vhSB(k,i) = vhNGBp(icl,niNGB)
              hvSB(k,i) = hp(k,l)
            ENDDO
          ENDDO
        END SELECT   
          
      END SELECT
   END DO

   ! ... Initialize flow variables at time n-1 & n
   uhEBp  = uhEB ; huEBp  = huEB ;
   uhWBp  = uhWB ; huWBp  = huWB ;
   vhNBp  = vhNB ; hvNBp  = hvNB ;
   vhSBp  = vhSB ; hvSBp  = hvSB ;
   uhEBpp = uhEBp; huEBpp = huEBp;
   uhWBpp = uhWBp; huWBpp = huWBp;
   vhNBpp = vhNBp; hvNBpp = hvNBp;
   vhSBpp = vhSBp; hvSBpp = hvSBp;

END SUBROUTINE openbc0

!************************************************************************
SUBROUTINE openbcUVH
!************************************************************************
!
!  Purpose: To assign values of water surface elevation or velocity
!           along open boundaries at the new (n+1) time level. 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL    :: areatot, uavg, vavg, sbc, qbc, hbc, dthrs_wse, dthrs_flw
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kmx, kmy, ks, ke
   INTEGER :: icl, m1, m2, niNGB, nbid, istat
   REAL    :: weight

   ! ... Return if no open boundaries exist
   IF (nopen <= 0) RETURN
   
   ! ... Initialize to zero counter for nested grid boundaries
   niNGB = 0;

   !.....Loop over open boundaries.....
   DO nn = 1, nopen

      SELECT CASE ( itype(nn) )
 
      !.....Case 1 -- wse specified.......................................
      CASE (1)

         ! Get new boundary value of zeta (NewOpenBC)
         dthrs_wse = dtsecOpenBC/3600.
         sbc = parab(0.,thrs,varsOpenBC(nn,1,:),dthrs_wse)

         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )
 
         ! West boundary
         CASE (1)
            i  = isbc(nn); 
            js = jsbc(nn); 
            je = jebc(nn)
            s(i,js:je) = sbc

         ! North boundary
         CASE (2)
            j  = jsbc(nn); 
            is = isbc(nn); 
            ie = iebc(nn)
            s(is:ie,j) = sbc

         ! East boundary
         CASE (3)
            i  = isbc(nn); 
            js = jsbc(nn); 
            je = jebc(nn)
            s(i,js:je) = sbc
 
         ! South boundary              
         CASE (4)
            j  = jsbc(nn); 
            is = isbc(nn);  
            ie = iebc(nn)          
            s(is:ie,j) = sbc

         END SELECT

      !.....Case 2 -- Free surface flow specified........................
      CASE (2)
      
         ! Get new boundary value of flow (in units of m**3/sec) (NewOpenBC)
         dthrs_flw = dtsecOpenBC/3600.
         qbc = parab(0.,thrs,varsOpenBC(nn,1,:),dthrs_flw)

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
 
            ! Get i-, j- indexes for bdry. point
            i  = isbc(nn); 
            js = jsbc(nn); 
            je = jebc(nn)

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO j = js,je
              IF(hhs(i,j)>hbc) hbc = hhs(i,j)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(i,j), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; huWB(:,js:je) = ZERO
            DO j = js,je
              kmx = kmz(i,j)
              DO k = k1, kmx
                huWB (k,j)=AMIN1(zlevel(k+1),hhs(i,j)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(huWB(k,j) <= HMIN) huWB(k,j) = ZERO;
                areatot = areatot + huWB(k,j) * dy
              ENDDO
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; uhWB(:,js:je) = 0.0E0
            DO j = js, je
              kmx = kmz(i,j)
              DO k = k1, kmx
                 uhWB(k,j) = uavg * huWB(k,j)
                 !PRINT *, j,k, uhWB(k,j), huWB(k,j)
              END DO
            ENDDO

         ! East boundary
         CASE (3)
 
            ! Get i-, j- indexes for bdry. point
            i  = isbc(nn); 
            js = jsbc(nn); 
            je = jebc(nn);

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO j = js,je
              IF(hhs(i,j)>hbc) hbc = hhs(i,j)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(i,j), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; huEB(:,js:je) = ZERO
            DO j = js,je
              kmx = kmz(i,j)
              DO k = k1, kmx
                huEB (k,j)=AMIN1(zlevel(k+1),hhs(i,j)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(huEB(k,j) <= HMIN) huEB(k,j) = ZERO;
                areatot = areatot + huEB(k,j) * dy
              ENDDO
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; uhEB(:,js:je) = 0.0E0
            DO j = js, je
              kmx = kmz(i,j)
              DO k = k1, kmx
                 uhEB(k,j) = uavg * huEB(k,j)
              END DO
            ENDDO

         ! North boundary
         CASE (2)
 
            ! Get i-, j- indexes for bdry. point
            j  = jsbc(nn); 
            is = isbc(nn); 
            ie = iebc(nn)

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO i = is,ie
              IF(hhs(i,j)>hbc) hbc = hhs(i,j)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(i,j), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; hvNB(:,is:ie) = ZERO
            DO i = is,ie
              kmy = kmz(i,j)
              DO k = k1, kmy
                hvNB (k,i)=AMIN1(zlevel(k+1),hhs(i,j)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(hvNB(k,i) <= HMIN) hvNB(k,i) = ZERO;
                areatot = areatot + hvNB(k,i) * dy
              ENDDO
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; vhNB(:,is:ie) = 0.0E0
            DO i = is, ie
              kmy = kmz(i,j)
              DO k = k1, kmy
                 vhNB(k,i) = vavg * hvNB(k,i)
               END DO
            ENDDO

         ! South boundary
         CASE (4)
 
            ! Get i-, j- indexes for bdry. point
            j  = jsbc(nn); 
            is = isbc(nn); 
            ie = iebc(nn)

            ! Get max. depth in section with bdry. values specified
            hbc = -1.e6; 
            DO i = is,ie
              IF(hhs(i,j)>hbc) hbc = hhs(i,j)
            ENDDO

            ! ... Define free surface location at bdry.
            sbc = MAX(s(i,j), -hbc+dzmin)

            ! ... Define thickness of bdry. wet cells & total area
            areatot = 0.0; hvSB(:,is:ie) = ZERO
            DO i = is,ie
              kmy = kmz(i,j)
              DO k = k1, kmy
                hvSB (k,i)=AMIN1(zlevel(k+1),hhs(i,j)) -        &
                &          AMAX1(zlevel(  k),-sbc)
                IF(hvSB(k,i) <= HMIN) hvSB(k,i) = ZERO;
                areatot = areatot + hvSB(k,i) * dy
              ENDDO
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; vhSB(:,is:ie) = 0.0E0
            DO i = is, ie
              kmy = kmz(i,j)
              DO k = k1, kmy
                 vhSB(k,i) = vavg * hvSB(k,i)
              END DO
            ENDDO
         END SELECT

      !.....Case 3 -- Submerged flow specified...........................
      CASE (3)

         ! Get first boundary value of flow (in units of m**3/sec)
         ! (Assume first value applies for t=0.0. It should be consistent
         !  with the initial condition for uh, vh, u,and v)
         ! Get new boundary value of flow (in units of m**3/sec) (NewOpenBC)
         dthrs_flw = dtsecOpenBC/3600.
         qbc = parab(0.,thrs,varsOpenBC(nn,1,:),dthrs_flw)

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
 
            ! Get i-, j- indexes for bdry. point
            i  = isbc(nn); 
            j  = jsbc(nn); 
            l  = ij2l(i,j)
            ks = iebc(nn); 
            ke = jebc(nn); 

            ! Make sure k- location is wett, define thickness of bdry. cell
            ! & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              huWB (k,j)= hp(k,l) 
              IF ( huWB(k,j) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged West bdry. on DRY CELL"
                STOP
              ENDIF
              areatot = areatot + huWB(k,j) * dy
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; 
            DO k = ks, ke
              uhWB(k,j) = uavg * huWB(k,j)       
            ENDDO

         ! East boundary
         CASE (3)
 
            ! Get i-, j- indexes for bdry. point
            i  = isbc(nn); 
            j  = jsbc(nn); 
            l  = ij2l(i,j)
            ks = iebc(nn); 
            ke = jebc(nn); 

            ! Make sure k- location is wett, define thickness of bdry. cell
            ! & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              huEB (k,j)= hp(k,l) 
              IF ( huEB(k,j) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged East bdry. on DRY CELL"
                STOP
              ENDIF
              areatot = areatot + huEB(k,j) * dy
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            uavg = qbc/areatot; 
            DO k = ks, ke
              uhEB(k,j) = uavg * huEB(k,j)
            ENDDO

         ! North boundary
         CASE (2)
 
            ! Get i-, j- indexes for bdry. point
            j  = jsbc(nn); 
            i  = isbc(nn);
            l  = ij2l(i,j) 
            ks = iebc(nn);
            ke = jebc(nn);

            ! Make sure k- location is wett, define thickness of bdry. cell
            ! & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              hvNB (k,i)= hp(k,l) 
              IF ( hvNB(k,i) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged North bdry. on a DRY CELL"
                STOP
              ENDIF
              areatot = areatot + hvNB(k,i) * dx
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; 
            DO k = ks, ke
              vhNB(k,i) = vavg * hvNB(k,i)
            ENDDO

         ! South boundary
         CASE (4)
 
            ! Get i-, j- indexes for bdry. point
            j  = jsbc(nn); 
            i  = isbc(nn);
            l  = ij2l(i,j) 
            ks = iebc(nn);
            ke = jebc(nn);

            ! Make sure k- location is wett, define thickness of bdry. cell
            ! & total area for outflow-inflow section
            areatot = 0.0
            DO k = ks, ke
              hvSB (k,i)= hp(k,l) 
              IF ( hvSB(k,i) <= ZERO ) THEN 
                PRINT *, "  "
                PRINT *, " ****STOPPING - Submerged South bdry. on a DRY CELL"
                STOP
              ENDIF
              areatot = areatot + hvSB(k,i) * dx
            ENDDO

            ! ... Define xs average velocity and estimate uh from there
            !     assuming velocity is uniform over xs
            vavg = qbc/areatot; 
            DO k = ks, ke
              vhSB(k,i) = vavg * hvSB(k,i)
            ENDDO

         END SELECT    

      !.....Case 4 -- Nested grid boundaries specified ..................
      CASE (4)
      
         ! ... Update counter of nested grid boundaries
         niNGB = niNGB + 1;

         ! ... Define weighting coefficients for records 
         weight  = (thrs - thrsNGBp)/(thrsNGB-thrsNGBp)

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE (iside(nn))
         CASE(1)
           i   = isbc(nn); 
           js  = jsbc(nn); 
           je  = jebc(nn);
           icl = 0
           DO j = js, je;
             l = ij2l(i,j) 
             DO k = k1, kmz(i,j)
               icl = icl + 1      
               uhWB(k,j) = uhNGB (icl,niNGB)*    weight + &
                           uhNGBp(icl,niNGB)*(1.-weight)
               !huWB(k,i)= h(k,l)
               huWB(k,j)= h(k,l)
             ENDDO 
           ENDDO 
         CASE(3)
           i   = isbc(nn); 
           js  = jsbc(nn); 
           je  = jebc(nn);
           icl = 0
           DO j = js, je; 
             l = ij2l(i,j) 
             DO k = k1, kmz(i,j)
               icl = icl + 1      
               uhEB(k,j) = uhNGB (icl,niNGB)*    weight + &
                           uhNGBp(icl,niNGB)*(1.-weight)
               huEB(k,j) = h(k,l)
             ENDDO
           ENDDO
         CASE(2)
           j   = jsbc(nn); 
           is  = isbc(nn); 
           ie  = iebc(nn);
           icl = 0
           DO i = is, ie; 
             l = ij2l(i,j) 
             DO k = k1, kmz(i,j)
               icl = icl + 1      
               vhNB(k,i) = vhNGB (icl,niNGB)*    weight + &
                           vhNGBp(icl,niNGB)*(1.-weight)
               !hvNB(k,j) = h(k,l)
               hvNB(k,i) = h(k,l)
             ENDDO
           ENDDO
         CASE(4)
           j  = jsbc(nn); 
           is = isbc(nn); 
           ie = iebc(nn);
           icl = 0
           DO i = is, ie; 
             l = ij2l(i,j) 
             DO k = k1, kmz(i,j)
               icl = icl + 1      
               vhSB(k,i) = vhNGB (icl,niNGB)*    weight + &
                           vhNGBp(icl,niNGB)*(1.-weight)
               !hvSB(k,j) = h(k,l)
               hvSB(k,i) = h(k,l)
             ENDDO
           ENDDO
         END SELECT   
 
      END SELECT
   END DO

END SUBROUTINE openbcUVH

!************************************************************************
SUBROUTINE readbcNGB
!************************************************************************
!
!  Purpose: To read in a new frame with nested grid boundary conditions
!           if needed
!
!------------------------------------------------------------------------

   REAL, ALLOCATABLE, DIMENSION (:,:) :: inputvar
   INTEGER :: icl, m1, m2, niNGB, nbiid, istat, nn

   IF (nopen <= 0 .OR. ioNBTOGGLE <= 0) RETURN

   ! ... Initialize counter for nested grid boundaries
   niNGB = 0; 

   ! ... Loop over open boundaries
   DO nn = 1, nopen 

     ! ... Cycle if not an embedded boundary
     IF ( itype(nn) < 4 ) CYCLE

     ! ... Update counter for nested grid boundaries ...    
     niNGB = niNGB + 1; 

     ! Save and read new frame if thrs > thrsNGB .......
     IF (thrs > thrsNGB) THEN

       ! ... Save variables from previous time .........
       thrsNGBp = thrsNGB; 
       uhNGBp = uhNGB; 
       vhNGBp = vhNGB; 
       scNGBp = scNGB; 
       trNGBp = trNGB;

       ! ... Set file ID ...............................            
       nbiid = nbiid0 + nn

       ! ... Allocate space for temporary input variable array 
       ALLOCATE( inputvar ( iptNBI(nn), 5+ntr ), STAT=istat )
       IF (istat /= 0) CALL allocate_error (istat,34)

       ! ... Read variables for NEXT FRAME variables ...
       READ(nbiid) thrsNGB, & 
           ((inputvar(m1,m2),m2=4,5+ntr),m1=1,iptNBI(nn))
  
       ! ... Assign variables  
       SELECT CASE (iside(nn))
       CASE(1,3)
         DO icl = 1, iptNBI(nn)
           uhNGB(icl,niNGB) = inputvar(icl,4)
           scNGB(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGB(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       CASE(2,4)
         DO icl = 1, iptNBI(nn)
           vhNGB(icl,niNGB) = inputvar(icl,4)
           scNGB(icl,niNGB) = inputvar(icl,5) 
           IF (ntr > 0) THEN
             trNGB(icl,1:ntr,niNGB) = inputvar(icl,6:5+ntr) 
           ENDIF
         ENDDO
       END SELECT
          
       DEALLOCATE (inputvar)

     ENDIF
   ENDDO

END SUBROUTINE readbcNGB

!************************************************************************
SUBROUTINE MODqqddrr4openBC
!************************************************************************
!
!  Purpose: To adjust the matrix coefficients used in the soln of the
!           continuity equation to account for open boundary conditions,
!           either wse or flow bdries.
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: i, j, nn, is, ie, js, je, ks, ke
   REAL    :: dt1, dtdx1, dtdy1

   !.....Constants.....
   dtdx1 = dtdx*tz; 
   dtdy1 = dtdy*tz

   !.....Loop over open boundaries.....
   DO nn = 1, nopen
      SELECT CASE ( itype(nn) )
 
      !.....Case 1 -- wse specified.....
      CASE (1)
      
         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )
         ! West boundary
         CASE (1)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            ! Adjust [qq] array for column of nodes just inside boundary
            qq(i+1,js:je) = qq(i+1,js:je) + sx(i,js:je)*s(i,js:je)
            ! Set sx(i,js:je) to zero in matrix (not really necessary)
            sx(i,js:je) = 0.0
            
         ! North boundary
         CASE (2)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! Adjust [qq] array for row of nodes just inside boundary
            qq(is:ie,j-1) = qq(is:ie,j-1) + sy(is:ie,j-1)*s(is:ie,j)
            ! Set sy(is:ie,j-1) to zero in matrix (necessary)
            sy(is:ie,j-1) = 0.0

         ! East boundary
         CASE (3)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            ! Adjust [qq] array for column of nodes just inside boundary
            qq(i-1,js:je) = qq(i-1,js:je) + sx(i-1,js:je)*s(i,js:je)
            ! Set sx(i-1,js:je) to zero in matrix (necessary)
            sx(i-1,js:je) = 0.0

         ! South boundary
         CASE (4)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! Adjust [qq] array for row of nodes just inside boundary
            qq(is:ie,j+1) = qq(is:ie,j+1) + sy(is:ie,j)*s(is:ie,j)
            ! Set sy(is:ie,j) to zero in matrix (not really necessary)
            sy(is:ie,j) = 0.0
            
         END SELECT

      !.....Case 2,4  -- Free surface flow specified.....
      CASE (2,4)

         ! Identify flow boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            ! Adjust [qq] array for flow rate into the
            ! column of nodes just inside the boundary 
            DO j = js, je
               qq(i,j) = qq(i,j) + dtdx1*SUM(uhWB  (k1:kmz(i,j),j)) &
                                 + dtdx1*SUM(uhWBpp(k1:kmz(i,j),j))
            END DO

         ! North boundary
         CASE (2)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            ! Adjust [qq] array for flow rate into the
            ! row of nodes just inside the boundary
            DO i = is, ie
               qq(i,j) = qq(i,j) - dtdy1*SUM(vhNB  (k1:kmz(i,j),i)) &
                                 - dtdy1*SUM(vhNBpp(k1:kmz(i,j),i))
            END DO

         ! East boundary
         CASE (3)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            ! Adjust [qq] array for flow rate into the
            ! column of nodes just inside the boundary
            DO j = js, je
               qq(i,j) = qq(i,j) - dtdx1*SUM(uhEB  (k1:kmz(i,j),j)) &
                                 - dtdx1*SUM(uhEBpp(k1:kmz(i,j),j)) 
            END DO

         ! South boundary
         CASE (4)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)          
            ! Adjust [qq] array for flow rate into the
            ! row of nodes just inside the boundary
            DO i = is, ie
               qq(i,j) = qq(i,j) + dtdy1*SUM(vhSB  (k1:kmz(i,j),i)) &
                                 + dtdy1*SUM(vhSBpp(k1:kmz(i,j),i))
            END DO
         END SELECT

      !.....Case 3 -- Submerged flow specified.....
      CASE (3)

         ! Identify flow boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
            i  = isbc(nn); 
            j  = jsbc(nn);
            ks = iebc(nn);
            ke = jebc(nn);
            ! Adjust [qq] array for flow rate into the
            ! column of nodes just inside the boundary 
            qq(i,j) = qq(i,j) + dtdx1*SUM(uhWB  (ks:ke,j)) &
                              + dtdx1*SUM(uhWBpp(ks:ke,j))

         ! North boundary
         CASE (2)
            j  = jsbc(nn); 
            i  = isbc(nn); 
            ks = iebc(nn);
            ke = jebc(nn);
            ! Adjust [qq] array for flow rate into the
            ! row of nodes just inside the boundary
            qq(i,j) = qq(i,j) - dtdy1*SUM(vhNB  (ks:ke,i)) &
                              - dtdy1*SUM(vhNBpp(ks:ke,i))

         ! East boundary
         CASE (3)
            i  = isbc(nn); 
            j  = jsbc(nn);
            ks = iebc(nn);
            ke = jebc(nn);
            ! Adjust [qq] array for flow rate into the
            ! column of nodes just inside the boundary
            qq(i,j) = qq(i,j) - dtdx1*SUM(uhEB  (ks:ke,j)) &
                              - dtdx1*SUM(uhEBpp(ks:ke,j)) 

         ! South boundary
         CASE (4)
            j  = jsbc(nn); 
            i  = isbc(nn); 
            ks = iebc(nn);
            ke = jebc(nn);
            ! Adjust [qq] array for flow rate into the
            ! row of nodes just inside the boundary
            qq(i,j) = qq(i,j) + dtdy1*SUM(vhSB  (ks:ke,i)) &
                              + dtdy1*SUM(vhSBpp(ks:ke,i))
         END SELECT

      END SELECT
   END DO

END SUBROUTINE MODqqddrr4openBC

!************************************************************************
SUBROUTINE MODcoef4openBC 
!************************************************************************
!
!  Purpose: To adjust the matrix coefficients used in the soln of the
!           continuity equation to account for open boundary conditions.
!
!------------------------------------------------------------------------
 
   !.....Local variables.....
   INTEGER :: i, j, nn, is, ie, js, je, m

   !.....Loop over open boundaries.....
   DO nn = 1, nopen
      SELECT CASE ( itype(nn) )
 
      !.....Case 1 -- wse specified.....
      CASE (1)

         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            DO j = js, je
              m = (i-ifirst)*ibdwd + (j-jfirst)+1
              coef(m,1) = 1.E1
              coef(m,2) = 0.
              coef(m,3) = 0.
              rhs(m)    = 1.E1 * s(i,j) 
              zeta(m)   = s(i,j)
            ENDDO           

         ! North boundary
         CASE (2)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            DO i = is, ie
              m = (i-ifirst)*ibdwd + (j-jfirst)+1
              coef(m,1) = 1.E1
              coef(m,2) = 0.
              coef(m,3) = 0.
              rhs(m)    = 1.E1 * s(i,j) 
              zeta(m)   = s(i,j)
            ENDDO        
   
         ! East boundary
         CASE (3)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            DO j = js, je
              m = (i-ifirst)*ibdwd + (j-jfirst)+1
              coef(m,1) = 1.E1
              coef(m,2) = 0.
              coef(m,3) = 0.
              rhs(m)    = 1.E1 * s(i,j)
              zeta(m)   = s(i,j) 
            ENDDO           
 
         ! South boundary
         CASE (4)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            DO i = is, ie
              m = (i-ifirst)*ibdwd + (j-jfirst)+1
              coef(m,1) = 1.E1
              coef(m,2) = 0.
              coef(m,3) = 0.
              rhs(m)    = 1.E1*s(i,j) 
              zeta(m)   = s(i,j)
            ENDDO           

         END SELECT

      !.....Case 2&3 -- Free surface & submerged flow specified.....
      CASE (2:)

      END SELECT
   END DO

END SUBROUTINE MODcoef4openBC


!************************************************************************
SUBROUTINE MODcoefA4openBC 
!************************************************************************
!
!  Purpose: To adjust the matrix coefficients used in the soln of the
!           continuity equation to account for open boundary conditions.
!
!------------------------------------------------------------------------
 
   !.....Local variables.....
   INTEGER :: i, j, nn, is, ie, js, je, m

   !.....Loop over open boundaries.....
   DO nn = 1, nopen
      SELECT CASE ( itype(nn) )
 
      !.....Case 1 -- wse specified.....
      CASE (1)

         ! Identify wse boundary as on the west, north, east, or south
         SELECT CASE ( iside(nn) )

         ! West boundary
         CASE (1)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            DO j = js, je
              m = ij2l(i,j)
              coeffA(m,1) = 1.E1
              coeffA(m,2) = 0.
              coeffA(m,3) = 0.
              coeffA(m,4) = 0.
              coeffA(m,5) = 0.
              rhs(m)    = 1.E1 * s(i,j) 
              zeta(m)   = s(i,j)
            ENDDO           

         ! North boundary
         CASE (2)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            DO i = is, ie
              m = ij2l(i,j)
              coeffA(m,1) = 1.E1
              coeffA(m,2) = 0.
              coeffA(m,3) = 0.
              coeffA(m,4) = 0.
              coeffA(m,5) = 0.
              rhs(m)    = 1.E1 * s(i,j) 
              zeta(m)   = s(i,j)
            ENDDO        
   
         ! East boundary
         CASE (3)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            DO j = js, je
              m = ij2l(i,j)
              coeffA(m,1) = 1.E1
              coeffA(m,2) = 0.
              coeffA(m,3) = 0.
              coeffA(m,4) = 0.
              coeffA(m,5) = 0.
              rhs(m)    = 1.E1 * s(i,j)
              zeta(m)   = s(i,j) 
            ENDDO           
 
         ! South boundary
         CASE (4)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            DO i = is, ie
              m = ij2l(i,j)
              coeffA(m,1) = 1.E1
              coeffA(m,2) = 0.
              coeffA(m,3) = 0.
              coeffA(m,4) = 0.
              coeffA(m,5) = 0.
              rhs(m)    = 1.E1*s(i,j) 
              zeta(m)   = s(i,j)
            ENDDO           

         END SELECT

      !.....Case 2&3 -- Free surface & submerged flow specified.....
      CASE (2:)

      END SELECT

   END DO

END SUBROUTINE MODcoefA4openBC

!************************************************************************
SUBROUTINE MODexmom4openBCX 
!************************************************************************
!
!  Purpose: To recompute ex matrix for velocity columns located along or 
!           next to open boundaries 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL :: twodt1
   REAL :: advx, advy, uE, uW, vN, vS, wU, wD, &
           ubdry, vbdry, uhbdry, vhbdry
   INTEGER :: i, j, k, l, istat, kmx, kmy, k1x, k1y, k1ne
   INTEGER :: nn, is, ie, js, je, ks, ke, nwlayers

   !.....Constant.....
   twodt1 = twodt*tz
 
   !.....Loop over open boundaries.....
   DO nn = 1, nopen
 
     SELECT CASE ( itype(nn) )
 
     !....... Case 1 -- wse specified .....................................
     CASE (1)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)

         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute the layer number for the bottom wet u-pt
           kmx = MIN(kmz(i+1,j), kmz(i,j))
           k1x =                 k1u(i,j)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmx-k1x) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term (only using advection)
           ex(:,l) = 0.0
           DO k = k1x,kmx
                                                                   
            ! Horizontal advection - Upwind differencing  
            uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
            uW = uhp(k,        l  ) +uhp(k  ,    l )
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* upp(k,lEC(l)) -           &
                     (uW+ABS(uW))* upp(k,    l ) -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO
         END DO

       ! ..... North boundary ......
       CASE (2)

       ! ..... East boundary ......
       CASE (3)
         i  = isbc(nn)-1; 
         js = jsbc(nn)  ; 
         je = jebc(nn)  ;
         DO j = js, je

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute the layer number for the bottom wet u-pt
           kmx = MIN(kmz(i+1,j), kmz(i,j))
           k1x =                 k1u(i,j)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmx-k1x) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term (only using advection)
           ex(:,l) = 0.0
           DO k = k1x,kmx
                                                                   
            ! Horizontal advection - Upwind differencing  
            uE = uhp(k,        l  ) +uhp(k  ,    l )
            uW = uhp(k,        l  ) +uhp(k  ,lWC(l)) 
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* upp(k,    l ) -           &
                     (uW+ABS(uW))* upp(k,lWC(l)) -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv) 
  
           END DO
         ENDDO           

       ! ..... South boundary ......
       CASE (4)

       END SELECT

     !....... Case 2,4 -- Free surface flow specified ......................
     CASE (2,4)

      ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)

         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn)
         DO j = js, je

           ! ... Map (i,j) into l-index
           l = ij2l(i,j); 

           ! ... Cycle if W-column is dry
           IF (.NOT. mask2d(i+1,j)) CYCLE

           ! Compute the layer number for the bottom wet u-pt
           kmx = MIN(kmz(i+1,j), kmz(i,j))
           k1x =                 k1u(i,j)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmx-k1x) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term (only using advection)
           ex(:,l) = 0.0
           DO k = k1x,kmx
                                                                   
            ! Horizontal advection - Upwind differencing  
            uhbdry = (uhWB(k,j)+uhWBpp(k,j))/2.
            IF (huWBpp(k,j)>ZERO) THEN
               ubdry  = uhWBpp(k,j)/huWBpp(k,j)
            ELSE
               ubdry  = 0.0
            ENDIF
            uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
            uW = uhp(k,        l  ) +uhbdry
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* upp(k,lEC(l)) -           &
                     (uW+ABS(uW))* ubdry         -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO
         END DO

       ! ..... North boundary ......
       CASE (2)

       ! ..... East boundary ......
       CASE (3)

         i  = isbc(nn)-1; 
         js = jsbc(nn)  ; 
         je = jebc(nn)  ; 
         DO j = js, je

           ! ... Map (i,j) into l-index
           l = ij2l(i,j); 

           ! ... Cycle if W-column is dry
           IF (.NOT. mask2d(i+1,j)) CYCLE

           ! Compute the layer number for the bottom wet u-pt
           kmx = MIN(kmz(i+1,j), kmz(i,j))
           k1x =                 k1u(i,j)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmx-k1x) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term (only using advection)
           ex(:,l) = 0.0
           DO k = k1x,kmx
                                                                   
            ! Horizontal advection - Upwind differencing  
            uhbdry = (uhEB(k,j)+uhEBpp(k,j))/2.
            IF (huEBpp(k,j)>ZERO) THEN
               ubdry  = uhEBpp(k,j)/huEBpp(k,j)
            ELSE
               ubdry  = 0.0
            ENDIF
            uE = uhp(k,        l  ) +uhbdry
            uW = uhp(k,    lWC(l) ) +uhp(k  ,    l ) 
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* ubdry         -           &
                     (uW+ABS(uW))* upp(k,lWC(l)) -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
		      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv) 
  
           END DO
         ENDDO           

       ! ..... South boundary ......
       CASE (4)

       END SELECT

     !....... Case 3 -- Submerged Flow specified ..........................
     CASE (3)

      ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)

         i  = isbc(nn); 
         j  = jsbc(nn); 
         ks = iebc(nn);
         ke = jebc(nn);
         l = ij2l(i,j);

         ! Compute explicit term (use only advection)
         ex(ks:ke,l) = 0.0

         DO k = ks,ke
                                                                   
            ! Horizontal advection - Upwind differencing  
            uhbdry = (uhWBpp(k,j)+uhWB  (k,j))/2.
            ubdry  =  uhWBpp(k,j)/huWBpp(k,j)
            uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
            uW = uhp(k,        l  ) +uhbdry
            vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
            vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
            wU = wp (k,    lEC(l) ) +wp (k  ,    l )
            wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
            advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                     (uE-ABS(uE))* upp(k,lEC(l)) -           &
                     (uW+ABS(uW))* ubdry         -           &
                     (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                  +( (vN+ABS(vN))* upp(k,    l ) +           &
                     (vN-ABS(vN))* upp(k,lNC(l)) -           &
                     (vS+ABS(vS))* upp(k,lSC(l)) -           &
                     (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                       (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * upp(k+1,l) +           &
                       (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
            !.....Final explicit term.....
            ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
         END DO

       ! ..... North boundary ......
       CASE (2)

         i  = isbc(nn); 
         j  = jsbc(nn); 
         ks = iebc(nn);
         ke = jebc(nn);

         ! Compute explicit term (use only advection) for u-points
		 ! to the East of a water column with North flow specified
         IF (mask2d(i+1,j)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);
           ex(ks:ke,l) = 0.0
           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             vhbdry = (vhNBpp(k,i)+vhNB  (k,i))/2.
             ! vbdry  =  vhNBpp(k,i)/hvNBpp(k,i)
             uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhp(k  ,    l )
             vN = vhp(k,    lEC(l) ) +vhbdry
             vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
             wU = wp (k,    lEC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
             advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                      (uE-ABS(uE))* upp(k,lEC(l)) -           &
                      (uW+ABS(uW))* upp(k,lWC(l)) -           & ! 30062008
                      (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* upp(k,    l ) +           &
                      (vN-ABS(vN))* upp(k,lNC(l)) -           &
                      (vS+ABS(vS))* upp(k,lSC(l)) -           &
                      (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

             ! Vertical advection - Upwind for near bdry. cells
             advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                        (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * upp(k+1,l) +           &
                        (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO 

         ENDIF

         ! Compute explicit term (use only advection) for u-points
         ! to the West of a water column with North flow specified
         IF (mask2d(i-1,j)) THEN 

           l = ij2l(i-1,j);

           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0
           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             vhbdry = (vhNBpp(k,i)+vhNB  (k,i))/2.
             ! vbdry  =  vhNBpp(k,i)/hvNBpp(k,i)
             uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhp(k  ,    l )
             vN = vhbdry             +vhp(k  ,    l )
             vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
             wU = wp (k,    lEC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
             advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                      (uE-ABS(uE))* upp(k,lEC(l)) -           &
                      (uW+ABS(uW))* upp(k,lWC(l)) -           &  ! 30062008
                      (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* upp(k,    l ) +           &
                      (vN-ABS(vN))* upp(k,lNC(l)) -           &
                      (vS+ABS(vS))* upp(k,lSC(l)) -           &
                      (vS-ABS(vS))* upp(k,    l ) ) / fourdy  

             ! Vertical advection - Upwind for near bdry. cells
             advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                        (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * upp(k+1,l) +           &
                        (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO

         ENDIF

       ! ..... East boundary ......
       CASE (3)

         i  = isbc(nn)-1; 
         j  = jsbc(nn)  ; 
         ks = iebc(nn)  ;
         ke = jebc(nn)  ;

         ! ... Map (i,j) into l-index
         l = ij2l(i,j);

         ! Compute explicit term (only using advection)
         ex(ks:ke,l) = 0.0
         DO k = ks,ke
                                                                   
           ! Horizontal advection - Upwind differencing  
           uhbdry = (uhEB(k,j)+uhEBpp(k,j))/2.
           ubdry  = uhEBpp(k,j)/huEBpp(k,j)
           uE = uhp(k,        l  ) +uhbdry
           uW = uhp(k,        l  ) +uhp(k  ,    l )
           vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
           vS = vhp(k,lSC(lEC(l))) +vhp(k  ,lSC(l))
           wU = wp (k,    lEC(l) ) +wp (k  ,    l )
           wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
           advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                    (uE-ABS(uE))* ubdry         -           &
                    (uW+ABS(uW))* upp(k,lWC(l)) -           &
                    (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                 +( (vN+ABS(vN))* upp(k,    l ) +           &
                    (vN-ABS(vN))* upp(k,lNC(l)) -           &
                    (vS+ABS(vS))* upp(k,lSC(l)) -           &
                    (vS-ABS(vS))* upp(k,    l ) ) / fourdy  
			      		                                  
           ! Vertical advection - Upwind for near bdry. cells
           advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                      (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                    -((wD+ABS(wD)) * upp(k+1,l) +           &
                      (wD-ABS(wD)) * upp(k  ,l)) / 4.

           !.....Final explicit term.....
           ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv) 
  
         END DO

       ! ..... South boundary ......
       CASE (4)

         i  = isbc(nn); 
         j  = jsbc(nn); 
         ks = iebc(nn);
         ke = jebc(nn);

         ! Compute explicit term (use only advection) for u-points
		 ! to the East of a water column with North flow specified
         IF (mask2d(i+1,j)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             vhbdry = (vhSBpp(k,i)+vhSB  (k,i))/2.
             ! vbdry  =  vhSBpp(k,i)/hvSBpp(k,i)
             uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhp(k  ,    l )
             vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
             vS = vhp(k,lSC(lEC(l))) +vhbdry
             wU = wp (k,    lEC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
             advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                      (uE-ABS(uE))* upp(k,lEC(l)) -           &
                      (uW+ABS(uW))* upp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* upp(k,    l ) +           &
                      (vN-ABS(vN))* upp(k,lNC(l)) -           &
                      (vS+ABS(vS))* upp(k,lSC(l)) -           &
                      (vS-ABS(vS))* upp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                        (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * upp(k+1,l) +           &
                        (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
               
           END DO

         ENDIF


         ! Compute explicit term (use only advection) for u-points
         ! to the West of a water column with North flow specified
         IF (mask2d(i-1,j)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i-1,j);
 
           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             vhbdry = (vhSBpp(k,i)+vhSB  (k,i))/2.
             ! vbdry  =  vhSBpp(k,i)/hvSBpp(k,i)
             uE = uhp(k,    lEC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhp(k  ,    l )
             vN = vhp(k,    lEC(l) ) +vhp(k  ,    l )
             vS = vhbdry             +vhp(k  ,lSC(l))
             wU = wp (k,    lEC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lEC(l) ) +wp (k+1,    l )        
             advx = ( (uE+ABS(uE))* upp(k,    l ) +           &
                      (uE-ABS(uE))* upp(k,lEC(l)) -           &
                      (uW+ABS(uW))* upp(k,lWC(l)) -           &
                      (uW-ABS(uW))* upp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* upp(k,    l ) +           &
                      (vN-ABS(vN))* upp(k,lNC(l)) -           &
                      (vS+ABS(vS))* upp(k,lSC(l)) -           &
                      (vS-ABS(vS))* upp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advx=advx+((wU+ABS(wU)) * upp(k  ,l) +           &
                        (wU-ABS(wU)) * upp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * upp(k+1,l) +           &
                        (wD-ABS(wD)) * upp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = uhpp(k,l) - twodt1*(advx*iadv)
                
           END DO

         ENDIF

       END SELECT ! End select from boundary locations (EWNS)

     END SELECT ! End select from boundary types (wse/flow/submerged)

   END DO ! End loop over open boundaries

END SUBROUTINE MODexmom4openBCX

!************************************************************************
SUBROUTINE MODexmom4openBCY 
!************************************************************************
!
!  Purpose: To recompute ex matrix for cells located next to cells with 
!           open boundary conditions specified 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL :: twodt1
   REAL :: advx, advy, uE, uW, vN, vS, wU, wD, &
           ubdry, vbdry, uhbdry, vhbdry
   INTEGER :: i, j, k, l, istat, kmx, kmy, k1x, k1y, k1ne
   INTEGER :: nn, is, ie, js, je, ks, ke, nwlayers

   !.....Constant.....
   twodt1 = twodt*tz
 
   !.....Loop over open boundaries.....
   DO nn = 1, nopen
 
     SELECT CASE ( itype(nn) )
 
     !....... Case 1 -- wse specified......................................
     CASE (1)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)

       ! ..... North boundary ......
       CASE (2)
         
         j  = jsbc(nn)-1; 
         is = isbc(nn)  ; 
         ie = iebc(nn)  ;
         DO i = is, ie

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute the layer number for top & bottom wet v-pt
           kmy = MIN(kmz(i,j+1), kmz(i,j))
           k1y =                 k1v(i,j)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmy-k1y) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term
           ex(:,l) = 0.0	     
           DO k = k1y,kmy
                                                                 
            ! Horizontal advection Upwind differencing 
            uE = uhp(k,    lNC(l) ) +uhp(k,    l )
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
            vN = vhp(k,        l  ) +vhp(k,    l )
            vS = vhp(k,        l  ) +vhp(k,lSC(l))
            wU = wp (k,    lNC(l) ) +wp (k,    l )
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
            advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                & +( (vN+ABS(vN))* vpp(k,    l ) +           &
            	&    (vN-ABS(vN))* vpp(k,    l ) -           &
                &    (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * vpp(k+1,l) +           &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO
         END DO

       ! ..... East boundary
       CASE (3)

       ! ..... South boundary ......
       CASE (4)

         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn)
         DO i = is, ie

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute the layer number for top & bottom wet v-pt
           kmy = MIN(kmz(i,j+1), kmz(i,j))
           k1y =                 k1v(i,j)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmy-k1y) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term
           ex(:,l) = 0.0	     
           DO k = k1y,kmy
                                                                 
            ! Horizontal advection Upwind differencing 
            uE = uhp(k,    lNC(l) ) +uhp(k,    l )
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
            vN = vhp(k,    lNC(l) ) +vhp(k,    l )
            vS = vhp(k,        l  ) +vhp(k,    l )
            wU = wp (k,    lNC(l) ) +wp (k,    l )
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
            advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                & +( (vN+ABS(vN))* vpp(k,    l ) +           &
            	&    (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                &    (vS+ABS(vS))* vpp(k,    l ) -           &
                &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind in near bdry. cells
            advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * vpp(k+1,l) +           &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO
         END DO 

       END SELECT

     !....... Case 2,4 -- Free surface flow specified.......................
     CASE (2,4)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ...... *** need to include code here
       CASE (1)

       ! ..... North boundary ......
       CASE (2)
         
         j  = jsbc(nn)-1; 
         is = isbc(nn)  ; 
         ie = iebc(nn)  ;
         DO i = is, ie

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! ... Cycle if N-column is dry
           IF ( .NOT. mask2d(i,j+1) ) CYCLE

           ! Compute the layer number for top & bottom wet v-pt
           kmy = MIN(kmz(i,j+1), kmz(i,j))
           k1y =                 k1v(i,j)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmy-k1y) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term
           ex(:,l) = 0.0	     
           DO k = k1y,kmy
                                                                 
            ! Horizontal advection Upwind differencing 
            vhbdry = (vhNB(k,i)+vhNBpp(k,i))/2.
            IF (hvNBpp(k,i)>ZERO) THEN
               vbdry  = vhNBpp(k,i)/hvNBpp(k,i)
            ELSE
               vbdry  = 0.0
            ENDIF
            uE = uhp(k,    lNC(l) ) +uhp(k,    l )
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
            vN = vhp(k,        l  ) +vhbdry
            vS = vhp(k,        l  ) +vhp(k,lSC(l))
            wU = wp (k,    lNC(l) ) +wp (k,    l )
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
            advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                & +( (vN+ABS(vN))* vpp(k,    l ) +           &
            	&    (vN-ABS(vN))* vbdry         -           &
                &    (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind for near bdry. cells
            advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * vpp(k+1,l) +           &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO
         END DO

       ! ..... East boundary ...... 
       CASE (3)

       ! ..... South boundary ......
       CASE (4)

         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! ... Cycle if N-column is dry
           IF ( .NOT. mask2d(i,j+1) ) CYCLE

           ! Compute the layer number for top & bottom wet v-pt
           kmy = MIN(kmz(i,j+1), kmz(i,j))
           k1y =                 k1v(i,j)

           ! Compute number of wet layers & skip for dry cells
           nwlayers = (kmy-k1y) + 1
           IF(nwlayers < 1) CYCLE

           ! Compute explicit term
           ex(:,l) = 0.0	     
           DO k = k1y,kmy
                                                                 
            ! Horizontal advection Upwind differencing 
            vhbdry = (vhSB(k,i)+vhSBpp(k,i))/2.
            IF (hvSBpp(k,i)>ZERO) THEN
               vbdry  = vhSBpp(k,i)/hvSBpp(k,i)
            ELSE
               vbdry  = 0.0
            ENDIF
            uE = uhp(k,    lNC(l) ) +uhp(k,    l )
            uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
            vN = vhp(k,    lNC(l) ) +vhp(k,    l )
            vS = vhp(k,        l  ) +vhbdry
            wU = wp (k,    lNC(l) ) +wp (k,    l )
            wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
            advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
                &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                & +( (vN+ABS(vN))* vpp(k,    l ) +           &
            	&    (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                &    (vS+ABS(vS))* vbdry         -           &
                &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

            ! Vertical advection - Upwind in near bdry. cells
            advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                       (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                     -((wD+ABS(wD)) * vpp(k+1,l) +           &
                       (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
            !.....Final explicit term.....
            ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO
         END DO 

       END SELECT

     !....... Case 3 -- Submerged flow specified...........................
     CASE (3)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)

         i  = isbc(nn); 
         j  = jsbc(nn); 
         ks = iebc(nn);
         ke = jebc(nn);

         ! Compute explicit term (use only advection) for v-points
		 ! to the North of a water column with EW flow specified
         IF (mask2d(i,j+1)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             uhbdry = (uhWBpp(k,j)+uhWB  (k,j))/2.
             ! ubdry  =  uhWBpp(k,j)/huWBpp(k,j) 
             uE = uhp(k,    lNC(l) ) +uhp(k  ,    l )
             uW = uhp(k,lWC(lNC(l))) +uhbdry
             vN = vhp(k,    lNC(l) ) +vhp(k  ,    l )
             vS = vhp(k,    lSC(l) ) +vhp(k  ,    l )
             wU = wp (k,    lNC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lNC(l) ) +wp (k+1,    l )        
             advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                      (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                      (uW+ABS(uW))* vpp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* vpp(k,    l ) +           &
                      (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                      (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                      (vS-ABS(vS))* vpp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                        (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * vpp(k+1,l) +           &
                        (wD-ABS(wD)) * vpp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
               
           END DO

         ENDIF

         ! Compute explicit term (use only advection) for v-points
		 ! to the South of a water column with EW flow specified
         IF (mask2d(i,j-1)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j-1);
 
           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             uhbdry = (uhWBpp(k,j)+uhWB  (k,j))/2.
             ! ubdry  =  uhWBpp(k,j)/huWBpp(k,j) 
             uE = uhp(k,    lNC(l) ) +uhp(k  ,    l )
             uW = uhp(k,    lWC(l) ) +uhbdry
             vN = vhp(k,    lNC(l) ) +vhp(k  ,    l )
             vS = vhp(k,    lSC(l) ) +vhp(k  ,    l )
             wU = wp (k,    lNC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lNC(l) ) +wp (k+1,    l )        
             advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                      (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                      (uW+ABS(uW))* vpp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* vpp(k,    l ) +           &
                      (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                      (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                      (vS-ABS(vS))* vpp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                        (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * vpp(k+1,l) +           &
                        (wD-ABS(wD)) * vpp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO

         ENDIF 

       ! ..... North boundary ......
       CASE (2)
         
         j  = jsbc(nn)-1; 
         i  = isbc(nn); 
         ks = iebc(nn);
         ke = jebc(nn); 
         l = ij2l(i,j);

         ! Compute explicit term
         ex(ks:ke,l) = 0.0	     
         DO k = ks,ke
                                                                 
           ! Horizontal advection Upwind differencing 
           vhbdry = (vhNB(k,i)+vhNBpp(k,i))/2.
           vbdry  = vhNBpp(k,i)/hvNBpp(k,i)
           uE = uhp(k,    lNC(l) ) +uhp(k,    l )
           uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
           vN = vhp(k,        l  ) +vhbdry
           vS = vhp(k,        l  ) +vhp(k,lSC(l))
           wU = wp (k,    lNC(l) ) +wp (k,    l )
           wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
           advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
               &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
               &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
               &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
               & +( (vN+ABS(vN))* vpp(k,    l ) +           &
               &    (vN-ABS(vN))* vbdry         -           &
               &    (vS+ABS(vS))* vpp(k,lSC(l)) -           &
               &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

           ! Vertical advection - Upwind for near bdry. cells
           advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                      (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                    -((wD+ABS(wD)) * vpp(k+1,l) +           &
                      (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
           !.....Final explicit term.....
           ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
         END DO

       ! ..... East boundary ......
       CASE (3)

         i  = isbc(nn); 
         j  = jsbc(nn); 
         ks = iebc(nn);
         ke = jebc(nn);

         ! Compute explicit term (use only advection) for v-points
		 ! to the North of a water column with EW flow specified
         IF (mask2d(i,j+1)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j);

           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             uhbdry = (uhEBpp(k,j)+uhEB  (k,j))/2.
             ! ubdry  =  uhEBpp(k,j)/huEBpp(k,j) 
             uE = uhp(k,    lNC(l) ) +uhbdry
             uW = uhp(k,lWC(lNC(l))) +uhp(k  ,lWC(l))
             vN = vhp(k,    lNC(l) ) +vhp(k  ,    l )
             vS = vhp(k,    lSC(l) ) +vhp(k  ,    l )
             wU = wp (k,    lNC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lNC(l) ) +wp (k+1,    l )        
             advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                      (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                      (uW+ABS(uW))* vpp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* vpp(k,    l ) +           &
                      (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                      (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                      (vS-ABS(vS))* vpp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                        (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * vpp(k+1,l) +           &
                        (wD-ABS(wD)) * vpp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
               
           END DO

         ENDIF

         ! Compute explicit term (use only advection) for v-points
		 ! to the South of a water column with EW flow specified
         IF (mask2d(i,j-1)) THEN 

           ! ... Map (i,j) into l-index
           l = ij2l(i,j-1);
 
           ! Compute explicit term (use only advection)
           ex(ks:ke,l) = 0.0

           DO k = ks,ke
                                                                   
             ! Horizontal advection - Upwind differencing  
             uhbdry = (uhEBpp(k,j)+uhEB  (k,j))/2.
             ! ubdry  =  uhEBpp(k,j)/huEBpp(k,j) 
             uE = uhp(k,        l  ) +uhbdry
             uW = uhp(k,lWC(lNC(l))) +uhp(k  ,lWC(l))
             vN = vhp(k,    lNC(l) ) +vhp(k  ,    l )
             vS = vhp(k,    lSC(l) ) +vhp(k  ,    l )
             wU = wp (k,    lNC(l) ) +wp (k  ,    l )
             wD = wp (k+1,  lNC(l) ) +wp (k+1,    l )        
             advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
                      (uE-ABS(uE))* vpp(k,lEC(l)) -           &
                      (uW+ABS(uW))* vpp(k,lWC(l)) -           & !30062008
                      (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
                   +( (vN+ABS(vN))* vpp(k,    l ) +           &
                      (vN-ABS(vN))* vpp(k,lNC(l)) -           &
                      (vS+ABS(vS))* vpp(k,lSC(l)) -           &
                      (vS-ABS(vS))* vpp(k,    l ) ) / fourdy   

             ! Vertical advection - Upwind for near bdry. cells
             advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                        (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                      -((wD+ABS(wD)) * vpp(k+1,l) +           &
                        (wD-ABS(wD)) * vpp(k  ,l)) / 4.
			      		                                  
             !.....Final explicit term.....
             ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
           END DO

         ENDIF 

       ! ..... South boundary ......
       CASE (4)

         j  = jsbc(nn); 
         i  = isbc(nn); 
         ks = iebc(nn);
         ke = jebc(nn);
         l = ij2l(i,j);

         ! Compute explicit term
         ex(ks:ke,l) = 0.0	     
         DO k = ks,ke
                                                                 
           ! Horizontal advection Upwind differencing 
           vhbdry = (vhSB(k,i)+vhSBpp(k,i))/2.
           vbdry  =  vhSBpp(k,i)/hvSBpp(k,i)
           uE = uhp(k,    lNC(l) ) +uhp(k,    l )
           uW = uhp(k,lWC(lNC(l))) +uhp(k,lWC(l))
           vN = vhp(k,    lNC(l) ) +vhp(k,    l )
           vS = vhp(k,        l  ) +vhbdry
           wU = wp (k,    lNC(l) ) +wp (k,    l )
           wD = wp (k+1,  lNC(l) ) +wp (k+1  ,l )        
           advy = ( (uE+ABS(uE))* vpp(k,    l ) +           &
               &    (uE-ABS(uE))* vpp(k,lEC(l)) -           &
               &    (uW+ABS(uW))* vpp(k,lWC(l)) -           &
               &    (uW-ABS(uW))* vpp(k,    l ) ) / fourdx  &
               & +( (vN+ABS(vN))* vpp(k,    l ) +           &
               &    (vN-ABS(vN))* vpp(k,lNC(l)) -           &
               &    (vS+ABS(vS))* vbdry         -           &
               &    (vS-ABS(vS))* vpp(k,    l ) ) / fourdy  

           ! Vertical advection - Upwind in near bdry. cells
           advy=advy+((wU+ABS(wU)) * vpp(k  ,l) +           &
                      (wU-ABS(wU)) * vpp(k-1,l)) / 4.       &
                    -((wD+ABS(wD)) * vpp(k+1,l) +           &
                      (wD-ABS(wD)) * vpp(k  ,l)) / 4.
                   
           !.....Final explicit term.....
           ex(k,l) = vhpp(k,l) - twodt1*(advy*iadv)
                
         END DO 

       END SELECT
     END SELECT
   END DO
 

END SUBROUTINE MODexmom4openBCY

!************************************************************************
SUBROUTINE MODvel4openBC 
!************************************************************************
!
!  Purpose: To adjust the velocity estimates provided by the code at time
!           n+1, for the presence of open boundaries
!
!------------------------------------------------------------------------
 
   !.....Local variables.....
   INTEGER :: i, j, k, l, nn, is, ie, js, je, ks, ke, m, kms
   REAL    :: deltaU

   !.....Loop over open boundaries.....
   DO nn = 1, nopen

      SELECT CASE ( itype(nn) )
 
      !....... Case 1 -- wse specified.....................................
      CASE (1)

         SELECT CASE ( iside(nn) )

         ! ... West boundary
         CASE (1)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            DO j = js, je
              l = ij2l(i,j)
              vh (:,l) = 0.0;
              wp (:,l) = wp(:, lEC(l));               
            ENDDO           

         ! ... North boundary
         CASE (2)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            DO i = is, ie
              l = ij2l(i,j)
              uh (:,l) = 0.0
              wp (:,l) = wp(:,lSC(l));
            ENDDO           

         ! ... East boundary
         CASE (3)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            DO j = js, je
              l  = ij2l(i,j)
              vh (:,l) = 0.0;
              wp (:,l) = wp(:, lWC(l)); 
            ENDDO           
 
         ! ... South boundary
         CASE (4)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            DO i = is, ie
              l = ij2l(i,j)
              uh (:,l) = 0.0
              wp (:,l) = wp(:,lNC(l));
            ENDDO           

         END SELECT

      !....... Case 2 -- Free surface flow specified ......................
      CASE (2,4)

         SELECT CASE ( iside(nn) )

         ! ... West boundary 
         CASE (1)
           i = isbc(nn); js = jsbc(nn); je = jebc(nn)
           DO j = js, je
             l   = ij2l(i,j)         
             kms = kmz(i,j); 
             wp(kms+1,l) = 0.0
             DO k = kms,k1,-1
               wp(k,l) = wp  (k+1,l)                            &
                   &   -(uh  (k  ,l)-uhWB  (k  ,    j )+        &
                         uhpp(k  ,l)-uhWBpp(k  ,    j ))/twodx  &
                   &   -(vh  (k  ,l)-vh    (k  ,lSC(l))+        &
                         vhpp(k  ,l)-vhpp  (k  ,lSC(l)))/twody
             END DO           
            ENDDO           

         ! ... North boundary
         CASE (2)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            DO i = is, ie
             l   = ij2l(i,j)         
             kms = kmz(i,j); 
             wp(kms+1,l) = 0.0
             DO k = kms,k1,-1
               wp(k,l) = wp  (k+1,l)                         &
                   &   -(uh    (k  ,l)-uh   (k,lWC(l))+       &
                         uhpp  (k  ,l)-uhpp (k,lWC(l)))/twodx &
                   &   -(vhNB  (k  ,i)-vh   (k  ,lSC(l))+        &
                         vhNBpp(k  ,i)-vhpp (k  ,lSC(l)))/twody
             END DO           
            ENDDO           

         ! ... East boundary
         CASE (3)
            i = isbc(nn); js = jsbc(nn); je = jebc(nn)
            DO j = js, je
             l   = ij2l(i,j)         
             kms = kmz(i,j); 
             wp(kms+1,l) = 0.0
             DO k = kms,k1,-1
               wp(k,l) = wp  (k+1,l)                         &
                   &   -(uhEB  (k  ,j)-uh   (k,lWC(l))+       &
                         uhEBpp(k  ,j)-uhpp (k,lWC(l)))/twodx &
                   &   -(vh    (k  ,l)-vh   (k,lSC(l))+       &
                         vhpp  (k  ,l)-vhpp (k,lSC(l)))/twody
             ENDDO           
            ENDDO

         ! ... South boundary
         CASE (4)
            j = jsbc(nn); is = isbc(nn); ie = iebc(nn)
            DO i = is, ie
             l   = ij2l(i,j)         
             kms = kmz(i,j); 
             wp(kms+1,l) = 0.0
             DO k = kms,k1,-1
               wp(k,l) = wp  (k+1,l)                           &
                   &   -(uh  (k  ,l)-uh    (k  ,lWC(l))+       &
                         uhpp(k  ,l)-uhpp  (k  ,lWC(l)))/twodx &
                   &   -(vh  (k  ,l)-vhSB  (k  ,    i )+       &
                         vhpp(k  ,l)-vhSBpp(k  ,    i ))/twody
             ENDDO 
            ENDDO          

         END SELECT

      !....... Case 3 -- Submerged flow specified .........................
      CASE (3)

         SELECT CASE ( iside(nn) )

         ! ... West boundary 
         CASE (1)

           i  = isbc(nn); 
           j  = jsbc(nn);
           ks = iebc(nn); 
           ke = jebc(nn); 
           l  = ij2l(i,j)         
           DO k = ks, ke
             deltaU = (uhWB(k,j)+uhWBpp(k,j))/twodx
             DO kms = k1, k
               wp(kms,l) = wp(kms,l)+deltaU                         
             END DO
           END DO

         ! ... North boundary
         CASE (2)

           i  = isbc(nn); 
           j  = jsbc(nn); 
           ks = iebc(nn); 
           ke = jebc(nn);
           l  = ij2l(i,j)        
           DO k = ks, ke
             deltaU = (vhNB(k,i)+vhNBpp(k,i))/twody
             DO kms = k1, k
               wp(kms,l) = wp(kms,l)-deltaU
             END DO    
           END DO           

         ! ... East boundary
         CASE (3)

           i  = isbc(nn); 
           j  = jsbc(nn); 
           ks = iebc(nn); 
           ke = jebc(nn); 
           l  = ij2l(i,j) 
           DO k = ks, ke
             deltaU = (uhEB(k,j)+uhEBpp(k,j))/twodx
             DO kms = k1, k
               wp(kms,l) = wp(kms,l)-deltaU   
             END DO
           END DO

         ! ... South boundary
         CASE (4)

           i  = isbc(nn); 
           j  = jsbc(nn); 
           ks = iebc(nn); 
           ke = jebc(nn);
           l  = ij2l(i,j)
           DO k = ks, ke
             deltaU = (vhSB(k,i)+vhSBpp(k,i))/twody
             DO kms = k1, k
               wp(kms,l) = wp(kms,l)+deltaU
             END DO    
           END DO           

         END SELECT

      END SELECT
   END DO

END SUBROUTINE MODvel4openBC

!************************************************************************
SUBROUTINE openbcSCA
!************************************************************************
!
!  Purpose: To assign values of active scalars along open wse boundaries 
!           at new (n+1) time level. 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL :: salbc, ss, delsal1, delsal2, delsal, dthrs_sal, advx
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kms, k1s

   ! ... Return if no open boundaries are specified
   IF (nopen <= 0) RETURN

   !.....Loop over open boundaries.....
   DO nn = 1, nopen

     SELECT CASE ( itype(nn) )      
 
     !.....Case 1 -- scalar specified on a wse BC .........................
     CASE (1)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(nn,2,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)
 
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je

           l = ij2l(i,j);

           ! Compute top & bottom layers
           kms = kmz(i,j)
           k1s = k1z(i,j)

           ! Compute explicit term (only using advection)
           DO k = k1s,kms
                                                                   
            ! Horizontal advection - Upwind differencing  
            IF ( uh(k,l) > 0.0 .AND. salbc > 0) THEN 
              sal(k,l) = salbc
            ELSE
              sal(k,l) = sal(k,lEC(l))
            ENDIF       
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)
 
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie

           l = ij2l(i,j);

           ! Compute top & bottom layers
           kms = kmz(i,j)
           k1s = k1z(i,j)

           ! Compute explicit term (only using advection)
           DO k = k1s,kms
                                                                   
            ! Horizontal advection - Upwind differencing  
            IF ( vh(k,lSC(l)) < 0.0 .AND. salbc > 0) THEN 
              sal(k,l) = salbc 
            ELSE
              sal(k,l) = sal(k,lSC(l))
            ENDIF

           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
 
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je

           l = ij2l(i,j);

           ! Compute top & bottom layers

           kms = kmz(i,j)
           k1s = k1z(i,j)

           ! Compute explicit term (only using advection)
           DO k = k1s,kms
                                                                   
            ! Horizontal advection - Upwind differencing  
            IF ( uh(k,lWC(l)) < 0.0 .AND. salbc > 0) THEN 
              sal(k,l) = salbc  
            ELSE
              sal(k,l) = sal(k,lWC(l))
            ENDIF       

           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)
 
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie

           l = ij2l(i,j);

           ! Compute top & bottom layers
           kms = kmz(i,j)
           k1s = k1z(i,j)

           ! Compute explicit term (only using advection)
           DO k = k1s,kms
                                                                   
            ! Horizontal advection - Upwind differencing  
            IF ( vh(k,l) > 0.0 .AND. salbc > 0) THEN 
              sal(k,l) = salbc  
            ELSE
              sal(k,l) = sal(k,lNC(l))
            ENDIF

           ENDDO
         ENDDO

        ENDSELECT

      ! ... Case 2,3,4,... -- scalar specified on a flow BC ...............
      CASE (2:)

      END SELECT

   ENDDO


END SUBROUTINE openbcSCA

!************************************************************************
SUBROUTINE MODexsal4openbc
!************************************************************************
!
!  Purpose: Modifies the ex arrays in the scalar transport equation to 
!           account for inflows/outflows. 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL    :: salbc, dthrs_sal, uE, uW, vN, vS, scE, scW, scN, scS
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kms, k1s
   INTEGER :: icl, m1, m2, niNGB, nbid
   REAL    :: weight  

   ! ... Return if no open boundaries are specified
   IF (nopen <= 0) RETURN

   ! ... Initialize indexes
   niNGB = 0

   !.....Loop over open boundaries.....
   DO nn = 1, nopen

     SELECT CASE ( itype(nn) )      
 
     !.....Case 1 -- scalar specified on wse BC.....
     CASE (1)

     !.....Case 2 -- scalar specified on free surface flow BC..............
     CASE (2)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(nn,2,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)
 
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je
           l = ij2l(i,j);
           kms = kmz(i,j)
           k1s = k1z(i,j)
           DO k = k1s,kms                                                                 
             ! ... Define velocity at boundary face
             uW = uhWB(k,j) + uhWBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uW >= 0.0 .AND. salbc > 0) THEN 
               scW = salbc  
             ELSE
               scW = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex- array 
             ex(k,l) = ex(k,l) +  uW * scW / twodx
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)
 
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie
           l = ij2l(i,j);
           kms = kmz(i,j)
           k1s = k1z(i,j)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             vN = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vN < 0.0 .AND. salbc > 0) THEN 
               scN = salbc  
             ELSE
               scN = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  vN * scN / twody 
           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
 
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je
           l = ij2l(i,j);
           kms = kmz(i,j)
           k1s = k1z(i,j)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             uE = uhEB(k,j) + uhEBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uE < 0.0 .AND. salbc > 0) THEN 
               scE = salbc  
             ELSE
               scE = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  uE * scE / twodx 
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)
 
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie
           l = ij2l(i,j);
           kms = kmz(i,j)
           k1s = k1z(i,j)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             vS = vhSB(k,i) + vhSBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vS > 0.0 .AND. salbc > 0) THEN 
               scS = salbc  
             ELSE
               scS = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) +  vS * scS / twody
           ENDDO
         ENDDO

        ENDSELECT

     ! ... Case 3 -- scalar specified on submerged flow BC ................
     CASE (3)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(nn,2,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)
 
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jsbc(nn);
         k1s= iebc(nn);
         kms= jebc(nn);
         DO j = js, je
           l = ij2l(i,j);
           DO k = k1s,kms               
             ! ... Define velocity at boundary face
             uW = uhWB(k,j) + uhWBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uW > 0.0 .AND. salbc > 0) THEN 
               scW = salbc  
             ELSE
               scW = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex- array 
             ex(k,l) = ex(k,l) +  uW * scW / twodx 
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)

         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = isbc(nn);
         k1s= iebc(nn);
         kms= jebc(nn);
         DO i = is, ie
           l = ij2l(i,j);
           DO k = k1s,kms              
             ! ... Define velocity at boundary face
             vN = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vN < 0.0 .AND. salbc > 0) THEN 
               scN = salbc  
             ELSE
               scN = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  vN * scN / twody 
           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
 
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jsbc(nn);
         k1s= iebc(nn);
         kms= jebc(nn);
         DO j = js, je
           l = ij2l(i,j);
           DO k = k1s,kms               
             ! ... Define velocity at boundary face
             uE = uhEB(k,j) + uhEBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uE < 0.0 .AND. salbc > 0) THEN 
               scE = salbc  
             ELSE
               scE = salpp(k,l)
             ENDIF
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  uE * scE / twodx 
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)

         i  = isbc(nn); 
         j  = jsbc(nn); 
         k1s= iebc(nn);
         kms= jebc(nn);
         l = ij2l(i,j);
         DO k = k1s, kms              
           ! ... Define velocity at boundary face
           vS = vhSB(k,i) + vhSBpp(k,i);                               
           ! ... Define scalar   at boundary face  
           IF ( vS > 0.0 .AND. salbc > 0) THEN 
             scS = salbc  
           ELSE
             scS = salpp(k,l)
           ENDIF
           ! ... Include boundary flux in ex 
           ex(k,l) = ex(k,l) +  vS * scS / twody 
         ENDDO
 
       END SELECT

       ! ... Case 4 - Active scalars at nested boundaries .................
       CASE (4) 

         ! ... Counter of nested grid boundaries
         niNGB = niNGB + 1;

         ! ... Define weighting coefficients for records 
         weight  = (thrs - thrsNGBp)/(thrsNGB-thrsNGBp)

         ! Identify boundary as on the west, north, east, or south
         SELECT CASE (iside(nn))
         CASE(1)
           i = isbc(nn); js = jsbc(nn); je = jebc(nn);
           icl = 0
           DO j = js, je; 
             l   = ij2l(i,j)
             DO k = k1, kmz(i,j)
             icl = icl + 1      
             ! ... Define scalar at boundary face  
             scW = scNGB (icl,niNGB)*    weight + &
                   scNGBp(icl,niNGB)*(1.-weight)
             ! ... Define velocity at boundary face
             uW  = uhWB(k,j) + uhWBpp(k,j); 
             ! ... Re-define scalar at boundary face if needed                               
             IF ( uW <= 0.0 ) scW = salpp(k,l)
             ! ... Include boundary flux in ex- array 
             ex(k,l) = ex(k,l) +  uW * scW / twodx
           ENDDO; ENDDO 
         CASE(3)
           i = isbc(nn); js = jsbc(nn); je = jebc(nn);
           icl = 0
           DO j = js, je; 
             l   = ij2l(i,j)
             DO k = k1, kmz(i,j)
             icl = icl + 1 
             ! ... Define scalar at boundary face  
             scE = scNGB (icl,niNGB)*    weight + &
                   scNGBp(icl,niNGB)*(1.-weight)
             ! ... Define velocity at boundary face
             uE  = uhEB(k,j) + uhEBpp(k,j);                               
             ! ... Redefine scalar at boundary face if needed 
             IF ( uE >= 0.0 ) scE = salpp(k,l)
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  uE * scE / twodx       
           ENDDO; ENDDO
         CASE(2)
           j = jsbc(nn); is = isbc(nn); ie = iebc(nn);
           icl = 0
           DO i = is, ie; 
             l   = ij2l(i,j)
             DO k = k1, kmz(i,j)
             icl = icl + 1      
             ! ... Define scalar at boundary face  
             scN = scNGB (icl,niNGB)*    weight + &
                   scNGBp(icl,niNGB)*(1.-weight)
             ! ... Define velocity at boundary face
             vN  = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Redefine scalar at boundary face if needed  
             IF ( vN >= 0.0 ) scN = salpp(k,l)
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  vN * scN / twody 
           ENDDO; ENDDO
         CASE(4)
           j = jsbc(nn); is = isbc(nn); ie = iebc(nn);
           icl = 0
           DO i = is, ie; 
             l   = ij2l(i,j)
             DO k = k1, kmz(i,j)
             icl = icl + 1      
             ! ... Define scalar at boundary face  
             scS = scNGB (icl,niNGB)*    weight + &
                   scNGBp(icl,niNGB)*(1.-weight)
             ! ... Define velocity at boundary face
             vS  = vhSB(k,i) + vhSBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vS <= 0.0 ) scS = salpp(k,l)
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) +  vS * scS / twody
           ENDDO; ENDDO

         END SELECT        

     END SELECT

   ENDDO

END SUBROUTINE MODexsal4openbc

!************************************************************************
SUBROUTINE openbcTracer (itr)
!************************************************************************
!
!  Purpose: To assign values of active scalars along open wse boundaries 
!           at new (n+1) time level. 
!
!------------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT (IN) :: itr

   !.....Local variables.....
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kms, k1s
   REAL    :: salbc, ss, delsal1, delsal2, delsal, dthrs_sal, &
              advx, uE, uW, vN, vS, wU, wD, twodt1

   ! ... Return if no open boundaries are specified
   IF (nopen <= 0) RETURN

   ! ... Constants used in solution
   twodt1 = twodt*tz

   !.....Loop over open boundaries.....
   DO nn = 1, nopen

     SELECT CASE ( itype(nn) )      
 
     !.....Case 1 -- scalar specified on a wse BC.....
     CASE (1)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(nn,2+itr,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)
 
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je

           ! ... Map 3D into 2D array index
           l = ij2l(i,j);
           ! Compute top & bottom layers
           kms = kmz(i,j)
           k1s = k1z(i,j)
           ! Compute explicit term (only using advection)
           DO k = k1s,kms                                                                 
             ! Horizontal advection - Upwind differencing  
             IF ( uh(k,l) > 0.0 .AND. salbc > 0.0) THEN 
               tracer(k,l,itr) = salbc  
             ELSE
               tracer(k,l,itr) = tracer(k,lEC(l),itr)
             ENDIF       
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)
 
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie
           ! ... Map 3D into 2D array counter
           l = ij2l(i,j);
           ! Compute top & bottom layers
           kms = kmz(i,j)
           k1s = k1z(i,j)
           ! Compute explicit term (only using advection)
           DO k = k1s,kms                                                                  
             ! Horizontal advection - Upwind differencing  
             IF ( vh(k,lSC(l)) < 0.0 .AND. salbc > 0.0) THEN 
               tracer(k,l,itr) = salbc  
             ELSE
               tracer(k,l,itr) = tracer(k,lSC(l),itr)
             ENDIF
           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
 
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je
           ! ... Map 3D into 2D array index
           l = ij2l(i,j);
           ! Compute top & bottom layers
           kms = kmz(i,j)
           k1s = k1z(i,j)
           ! Compute explicit term (only using advection)
           DO k = k1s,kms                                                        
             ! Horizontal advection - Upwind differencing  
             IF ( uh(k,lWC(l)) < 0.0 .AND. salbc > 0.0) THEN 
               tracer(k,l,itr) = salbc  
             ELSE
               tracer(k,l,itr) = tracer(k,lWC(l),itr)
             ENDIF       
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)
 
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie
           ! ... Map 3D into 2D array index
           l = ij2l(i,j);
           ! Compute top & bottom layers
           kms = kmz(i,j)
           k1s = k1z(i,j)
           ! Compute explicit term (only using advection)
           DO k = k1s,kms                                                               
             ! Horizontal advection - Upwind differencing  
             IF ( vh(k,l) > 0.0 .AND. salbc > 0.0) THEN 
               tracer(k,l,itr) = salbc  
             ELSE
               tracer(k,l,itr) = tracer(k,lNC(l),itr)
             ENDIF
           ENDDO
         ENDDO

        ENDSELECT

      CASE (2:)

        CYCLE

      END SELECT

   ENDDO


END SUBROUTINE openbcTracer

!************************************************************************
SUBROUTINE MODexTracer4openbc (itr)
!************************************************************************
!
!  Purpose: Modifies the ex arrays in the scalar transport equation to 
!           account for inflows/outflows. 
!
!------------------------------------------------------------------------

   ! ... Arguments
   INTEGER, INTENT (IN) :: itr

   !.....Local variables.....
   REAL    :: salbc, dthrs_sal, uE, uW, vN, vS, scE, scW, scN, scS
   INTEGER :: i, j, k, l, ios, nn, is, ie, js, je, kms, k1s
   INTEGER :: icl, niNGB
   REAL    :: weight  

   ! ... Return if no open boundaries are specified
   IF (nopen <= 0) RETURN

   niNGB = 0

   !.....Loop over open boundaries.....
   DO nn = 1, nopen

     SELECT CASE ( itype(nn) )      

     CASE (1) 

       CYCLE
 
     !.....Case 2 -- scalar specified on free surface flow BC............
     CASE (2)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(nn,2+itr,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je
           l = ij2l(i,j);
           kms = kmz(i,j)
           k1s = k1z(i,j)
           DO k = k1s,kms                                                                 
             ! ... Define velocity at boundary face
             uW = uhWB(k,j) + uhWBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uW > 0.0 .AND. salbc > 0.0) THEN 
               scW = salbc  
             ELSE
               scW = tracerpp(k,l,itr)
             ENDIF
             ! ... Include boundary flux in ex- array 
             ex(k,l) = ex(k,l) +  uW * scW / twodx
           ENDDO
         ENDDO

       ! ..... North boundary ......
       CASE (2)
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie
           l = ij2l(i,j);
           kms = kmz(i,j)
           k1s = k1z(i,j)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             vN = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vN < 0.0 .AND. salbc > 0.0) THEN 
               scN = salbc  
             ELSE
               scN = tracerpp(k,l,itr)
             ENDIF
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  vN * scN / twody 
           ENDDO
         ENDDO

       ! ..... East boundary ......
       CASE (3)
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         DO j = js, je
           l = ij2l(i,j);
           kms = kmz(i,j)
           k1s = k1z(i,j)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             uE = uhEB(k,j) + uhEBpp(k,j);                               
             ! ... Define scalar   at boundary face  
             IF ( uE < 0.0 .AND. salbc > 0.0) THEN 
               scE = salbc  
             ELSE
               scE = tracerpp(k,l,itr)
             ENDIF
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  uE * scE / twodx 
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE (4)
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         DO i = is, ie
           l = ij2l(i,j);
           kms = kmz(i,j)
           k1s = k1z(i,j)
           DO k = k1s,kms
             ! ... Define velocity at boundary face
             vS = vhSB(k,i) + vhSBpp(k,i);                               
             ! ... Define scalar   at boundary face  
             IF ( vS > 0.0 .AND. salbc > 0.0) THEN 
               scS = salbc  
             ELSE
               scS = tracerpp(k,l,itr)
             ENDIF
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) +  vS * scS / twody
           ENDDO
         ENDDO

        ENDSELECT

     ! ... Case 3 -- scalar specified on submerged flow BC ..............
     CASE (3)

       !.....Get new boundary value of scalar.....
       ! ... Set value of scalar at OpenBC nn - 
       dthrs_sal = dtsecOpenBC/3600.
       salbc = parab(0.,thrs,varsOpenBC(nn,2+itr,:),dthrs_sal)

       ! Identify wse boundary as on the west, north, east, or south
       SELECT CASE ( iside(nn) )

       ! ..... West boundary ......
       CASE (1)
         i   = isbc(nn); 
         j   = jsbc(nn); 
         k1s = iebc(nn);
         kms = jebc(nn);
         l   = ij2l(i,j);
         DO k = k1s,kms               
           ! ... Define velocity at boundary face
           uW = uhWB(k,j) + uhWBpp(k,j);                               
           ! ... Define scalar   at boundary face  
           IF ( uW > 0.0 .AND. salbc > 0.0) THEN 
             scW = salbc  
           ELSE
             scW = tracerpp(k,l,itr)
           ENDIF
           ! ... Include boundary flux in ex- array 
           ex(k,l) = ex(k,l) +  uW * scW / twodx 
         ENDDO

       ! ..... North boundary ......
       CASE (2)
         j   = jsbc(nn); 
         i   = isbc(nn); 
         k1s = iebc(nn);
         kms = jebc(nn);
         l   = ij2l(i,j);
         DO k = k1s,kms              
           ! ... Define velocity at boundary face
           vN = vhNB(k,i) + vhNBpp(k,i);                               
           ! ... Define scalar   at boundary face  
           IF ( vN < 0.0 .AND. salbc > 0.0) THEN 
             scN = salbc  
           ELSE
             scN = tracerpp(k,l,itr)
           ENDIF
           ! ... Include boundary flux in ex 
           ex(k,l) = ex(k,l) -  vN * scN / twody 
         ENDDO

       ! ..... East boundary ......
       CASE (3)
         i   = isbc(nn); 
         j   = jsbc(nn); 
         k1s = iebc(nn);
         kms = jebc(nn);
         l   = ij2l(i,j);
         DO k = k1s,kms               
           ! ... Define velocity at boundary face
           uE = uhEB(k,j) + uhEBpp(k,j);                               
           ! ... Define scalar   at boundary face  
           IF ( uE < 0.0 .AND. salbc > 0.0) THEN 
             scE = salbc  
           ELSE
             scE = tracerpp(k,l,itr)
           ENDIF
           ! ... Include boundary flux in ex 
           ex(k,l) = ex(k,l) -  uE * scE / twodx 
         ENDDO

       ! ..... South boundary ......
       CASE (4)
         i   = isbc(nn); 
         j   = jsbc(nn); 
         k1s = iebc(nn);
         kms = jebc(nn);
         l   = ij2l(i,j);
         DO k = k1s, kms              
           ! ... Define velocity at boundary face
           vS = vhSB(k,i) + vhSBpp(k,i);                               
           ! ... Define scalar   at boundary face  
           IF ( vS > 0.0 .AND. salbc > 0.0) THEN 
             scS = salbc  
           ELSE
             scS = tracerpp(k,l,itr)
           ENDIF
           ! ... Include boundary flux in ex 
           ex(k,l) = ex(k,l) +  vS * scS / twody 
         ENDDO
 
       END SELECT

     ! ... Case 4 - Nested grid boundaries ...............................
     CASE (4)

       ! ... Counter of nested grid boundaries
       niNGB = niNGB + 1;

       ! ... Define weighting coefficients for records 
       weight  = (thrs - thrsNGBp)/(thrsNGB-thrsNGBp)

       ! Identify boundary as on the west, north, east, or south
       SELECT CASE (iside(nn))

       ! ..... West boundary ......
       CASE (1)
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         icl = 0
         DO j = js, je; 
           l   = ij2l(i,j)
           DO k = k1, kmz(i,j)
             icl = icl + 1      
             ! ... Define scalar at boundary face  
             scW = trNGB (icl,niNGB,itr)*    weight + &
                   trNGBp(icl,niNGB,itr)*(1.-weight)
             ! ... Define velocity at boundary face
             uW  = uhWB(k,j) + uhWBpp(k,j); 
             ! ... Re define scalar at boundary face if needed                               
             IF ( uW < 0.0 ) scW = tracerpp(k,l,itr)
             ! ... Include boundary flux in ex- array 
             ex(k,l) = ex(k,l) +  uW * scW / twodx
           ENDDO; ENDDO 

       ! ..... East boundary ......
       CASE (3)
         i  = isbc(nn); 
         js = jsbc(nn); 
         je = jebc(nn);
         icl = 0
         DO j = js, je; 
           l   = ij2l(i,j)
           DO k = k1, kmz(i,j)
           icl = icl + 1 
           ! ... Define scalar at boundary face  
           scE = trNGB (icl,niNGB,itr)*    weight + &
                 trNGBp(icl,niNGB,itr)*(1.-weight)
           ! ... Define velocity at boundary face
           uE  = uhEB(k,j) + uhEBpp(k,j);                               
           ! ... Redefine scalar at boundary face if needed 
           IF ( uE >= 0.0 ) scE = tracerpp(k,l,itr)
           ! ... Include boundary flux in ex 
           ex(k,l) = ex(k,l) -  uE * scE / twodx       
         ENDDO; ENDDO

       ! ..... North boundary ......
       CASE (2)
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         icl = 0
         DO i = is, ie; 
           l   = ij2l(i,j)
           DO k = k1, kmz(i,j)
             icl = icl + 1      
             ! ... Define scalar at boundary face  
             scN = trNGB (icl,niNGB,itr)*    weight + &
                   trNGBp(icl,niNGB,itr)*(1.-weight)
             ! ... Define velocity at boundary face
             vN  = vhNB(k,i) + vhNBpp(k,i);                               
             ! ... Redefine scalar at boundary face if needed  
             IF ( vN >= 0.0 ) scN = tracerpp(k,l,itr)
             ! ... Include boundary flux in ex 
             ex(k,l) = ex(k,l) -  vN * scN / twody 
           ENDDO
         ENDDO

       ! ..... South boundary ......
       CASE(4)
         j  = jsbc(nn); 
         is = isbc(nn); 
         ie = iebc(nn);
         icl = 0
         DO i = is, ie; 
           l   = ij2l(i,j)
           DO k = k1, kmz(i,j)
           icl = icl + 1      
           ! ... Define scalar at boundary face  
           scS = trNGB (icl,niNGB,itr)*    weight + &
                 trNGBp(icl,niNGB,itr)*(1.-weight)
           ! ... Define velocity at boundary face
           vS  = vhSB(k,i) + vhSBpp(k,i);                               
           ! ... Define scalar at boundary face if needed  
           IF ( vS <= 0.0 ) scS = tracerpp(k,l,itr)
           ! ... Include boundary flux in ex 
           ex(k,l) = ex(k,l) +  vS * scS / twody
         ENDDO; ENDDO
       END SELECT      

     END SELECT

   ENDDO

END SUBROUTINE MODexTracer4openbc

!************************************************************************
SUBROUTINE surfbc0
!************************************************************************
!
!  Purpose: This routine is called at the beginning of the program
!           to open files with heat boundary condition data, to
!           read the boundary condition time series data, to 
!           assign heat_sources at time t=0. It uses
!           the same scheme as used in openbc routines to read bc values. 
!
!------------------------------------------------------------------------

   !.....Local variables.....
   INTEGER :: ios, istat, nn, j, npsurfbc, nvsurfbc, imet
   CHARACTER(LEN=14) :: surfbcfmt, metxyfmt

   SELECT CASE (ifSurfBC) 

   ! ... Surface boundary conditions set to constant values (no heat flux)
   CASE (0) 

     RETURN

   ! .... Surface boundary conditions read from files - PRE-PROCESS mode
   CASE (1)

     !               ----- Open files with heatflux surface bc data-----   
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----
     ! Skip over first six header records in salinity boundary condition file 
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 101 )
   
     ! Read number of points in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npsurfbc
     IF (ios /= 0) CALL input_error ( ios, 102 )
   
     ! Allocate space for the array of data
     ALLOCATE ( surfbc1(nvSurfbcP,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 103 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbcP
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
       READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
            (surfbc1(nn,j), nn = 1, nvSurfbcP)
       IF (ios /= 0) CALL input_error ( ios, 104 )
     END DO
   
     !           ----- Assign heat flux terms at time t=0.0-----
     eta = surfbc1(1,1)
     Qsw = surfbc1(2,1)
     Qn  = surfbc1(3,1) 

     !           ----- Assign momentum flux terms at t=0.0 -----
     cdw = surfbc1(4,1)
     uair= surfbc1(5,1)
     vair= surfbc1(6,1)

     ! ... Set heat sources for each cell
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   ! .... Surface boundary conditions RUN-TIME (I) mode
   CASE (2)

     !               ----- Open files with heatflux surface bc data-----   
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----
     ! Skip over first six header records in salinity boundary condition file 
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 105 )
   
     ! Read number of points in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npSurfbc
     IF (ios /= 0) CALL input_error ( ios, 106 )
   
     ! Allocate space for the array of data
     ALLOCATE ( surfbc1(nvSurfbcR,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 107 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbcR
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
       READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
            (surfbc1(nn,j), nn = 1, nvSurfbcR)
       IF (ios /= 0) CALL input_error ( ios, 108 )
     END DO
   
     !           ----- Assign heat flux terms at time t=0.0-----
     eta = surfbc1(1,1)
     Qsw = surfbc1(2,1)
     Ta  = surfbc1(3,1)        ! Data in oC
     Pa  = surfbc1(4,1)        ! Pascals
     Rh  = surfbc1(5,1)        ! fraction (i.e. < 1)
     Cc  = surfbc1(6,1)        ! fraction (i.e. < 1)

     !           ----- Assign momentum flux terms at t=0.0 -----
     cdw = surfbc1(7,1)
     uair= surfbc1(8,1)
     vair= surfbc1(9,1)

     ! ... Set heat sources for each cell
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   ! .... Surface boundary conditions RUN-TIME (II) mode
   CASE (3)

     !               ----- Open files with heatflux surface bc data-----   
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----
     ! Skip over first six header records in salinity boundary condition file 
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 109 )
   
     ! Read number of points in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npSurfbc
     IF (ios /= 0) CALL input_error ( ios, 110 )
   
     ! Allocate space for the array of data
     ALLOCATE ( surfbc1(nvSurfbcR,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 111 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbcR
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
       READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
            (surfbc1(nn,j), nn = 1, nvSurfbcR)
       IF (ios /= 0) CALL input_error ( ios, 112 )
     END DO
   
     !           ----- Assign heat flux terms at time t=0.0-----
     eta = surfbc1(1,1)        ! m-1
     Qsw = surfbc1(2,1)        ! W/m2
     Ta  = surfbc1(3,1)        ! Data in oC
     Pa  = surfbc1(4,1)        ! Pascals
     Rh  = surfbc1(5,1)        ! fraction (i.e. < 1)
     Qlw = surfbc1(6,1)        ! W/m2

     !           ----- Assign momentum flux terms at t=0.0 -----
     cdw = surfbc1(7,1)
     uair= surfbc1(8,1)
     vair= surfbc1(9,1)

     ! ... Set heat sources for each cell
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   CASE (10) ! Space & Time varying met. variables - Heat budget on run-time (I) mode

     !               ----- Open files with heatflux surface bc data-----
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----

     ! Skip over first six header records in boundary condition file
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 113 )  

     ! Read number of met stations in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) nmetstat
     IF (ios /= 0) CALL input_error ( ios, 114 )

     ! Allocate space for metxy matrix containing x,y locations 
     ! of all met stations & variables that store individual met records
     ! for a given time step for each of the variables
     ALLOCATE ( metxy (nmetstat*2), &
                Qsw2D (nmetstat  ), Ta2D  (nmetstat),           &
                RH2D  (nmetstat  ), Cc2D  (nmetstat),           &
                uair2D(nmetstat  ), vair2D(nmetstat), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 115 )

     ! Allocate space for weightst matrix containing weighting coefficients 
     ! assigned to each of the met stations for each grid point
     ALLOCATE ( weightst(im1,jm1,nmetstat), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 116 )

     ! Write the format of the data records into an internal file
     WRITE (UNIT=metxyfmt, FMT='("(10X,",I3,"G11.2)")') nmetstat*2

     ! Read x,y position (in grid units) for each met station
     READ (UNIT=i53, FMT=metxyfmt, IOSTAT=ios)(metxy(nn), nn = 1, nmetstat*2)
     IF (ios /= 0) CALL input_error ( ios, 117 )
 
     ! ... Initialize Interpolation Schemes (weights for Barnes)
     ! CALL InitializeInterpolationMethods ! Barnes Corrected
     CALL InitializeBarnesInterpolation

     ! Read number of records in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npSurfbc
     IF (ios /= 0) CALL input_error ( ios, 110 )
   
     ! Allocate space for the array of data
     nvSurfbc = 6 * nmetstat + 2 
     ALLOCATE ( surfbc1(nvSurfbc,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 111 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbc
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
         READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
		      (surfbc1(nn,j), nn = 1, nvSurfbc)
         IF (ios /= 0) CALL input_error ( ios, 112 )
     END DO
 
     !           ----- Assign surfBC at time t=0.0-----
     eta           = surfbc1(1             ,1)
     Pa            = surfbc1(2             ,1)
     DO imet = 1, nmetstat
       Qsw2D (imet)  = surfbc1((imet-1)*6 + 3,1)
       Ta2D  (imet)  = surfbc1((imet-1)*6 + 4,1) 
       RH2D  (imet)  = surfbc1((imet-1)*6 + 5,1) 
       Cc2D  (imet)  = surfbc1((imet-1)*6 + 6,1) 
       uair2D(imet)  = surfbc1((imet-1)*6 + 7,1)
       vair2D(imet)  = surfbc1((imet-1)*6 + 8,1)
     ENDDO

     ! ... Distribute heat and momentum sources entering through free surface
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   CASE (11) ! Space & Time varying met. variables - Heat budget on run-time (II) mode

     !               ----- Open files with heatflux surface bc data-----
     OPEN (UNIT=i53, FILE='surfbc.txt', STATUS="OLD", IOSTAT=ios)
     IF (ios /= 0) CALL open_error ( "Error opening surfbc.txt", ios )

     !               -----Read files with heatflux surface bc data-----

     ! Skip over first six header records in boundary condition file
     READ (UNIT=i53, FMT='(/////)', IOSTAT=ios)
     IF (ios /= 0) CALL input_error ( ios, 113 )  

     ! Read number of met stations in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) nmetstat
     IF (ios /= 0) CALL input_error ( ios, 114 )

     ! Allocate space for metxy matrix containing x,y locations 
     ! of all met stations & variables that store individual met records
     ! for a given time step for each of the variables
     ALLOCATE ( metxy (nmetstat*2),                             &
                Qsw2D (nmetstat  ), Ta2D  (nmetstat),           &
                RH2D  (nmetstat  ), Qlw2D (nmetstat),           &
                uair2D(nmetstat  ), vair2D(nmetstat), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 115 )

     ! Allocate space for weightst matrix containing weighting coefficients 
     ! assigned to each of the met stations for each grid point
     ALLOCATE ( weightst(im1,jm1,nmetstat), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 116 )

     ! Write the format of the data records into an internal file
     WRITE (UNIT=metxyfmt, FMT='("(10X,",I3,"G11.2)")') nmetstat*2

     ! Read x,y position (in grid units) for each met station
     READ (UNIT=i53, FMT=metxyfmt, IOSTAT=ios)(metxy(nn), nn = 1, nmetstat*2)
     IF (ios /= 0) CALL input_error ( ios, 117 )
 
     ! ... Initialize Interpolation Schemes (weights for Barnes)
     ! CALL InitializeInterpolationMethods ! Barnes Corrected
     CALL InitializeBarnesInterpolation

     ! Read number of records in file from seventh header record
     READ (UNIT=i53, FMT='(10X,I7)', IOSTAT=ios) npSurfbc
     IF (ios /= 0) CALL input_error ( ios, 110 )
   
     ! Allocate space for the array of data
     nvSurfbc = 6 * nmetstat + 2 
     ALLOCATE ( surfbc1(nvSurfbc,npSurfbc), STAT=istat )
     IF (istat /= 0) CALL allocate_error ( istat, 111 )
   
     ! Write the format of the data records into an internal file
     WRITE (UNIT=surfbcfmt, FMT='("(10X,",I3,"G11.2)")') nvSurfbc
   
     ! Read data array and store it in memory
     DO j = 1, npSurfbc
         READ (UNIT=i53, FMT=surfbcfmt, IOSTAT=ios) &
		      (surfbc1(nn,j), nn = 1, nvSurfbc)
         IF (ios /= 0) CALL input_error ( ios, 112 )
     END DO
 
     !           ----- Assign surfBC at time t=0.0-----
     eta           = surfbc1(1             ,1)
     Pa            = surfbc1(2             ,1)
     DO imet = 1, nmetstat
       Qsw2D (imet)  = surfbc1((imet-1)*6 + 3,1)
       Ta2D  (imet)  = surfbc1((imet-1)*6 + 4,1) 
       RH2D  (imet)  = surfbc1((imet-1)*6 + 5,1) 
       Qlw2D (imet)  = surfbc1((imet-1)*6 + 6,1) 
       uair2D(imet)  = surfbc1((imet-1)*6 + 7,1)
       vair2D(imet)  = surfbc1((imet-1)*6 + 8,1)
     ENDDO

     ! ... Distribute heat and momentum sources entering through free surface
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   END SELECT
   
END SUBROUTINE surfbc0

!************************************************************************
SUBROUTINE surfbc
!************************************************************************
!
!  Purpose: To define heat & momentum fluxes through the free surface
!           at each time  
!
!------------------------------------------------------------------------

   !.....Local variables.....
   REAL    :: dthrs_surfbc
   INTEGER :: i, j, k, ios, nn, is, ie, js, je, kb, isalin, itest, imet

   SELECT CASE (ifSurfBC) 

   !               ----- No surface bc data----------------
   CASE (0)

      uair = -wa * SIN(pi*phi/180.);
      vair = -wa * COS(pi*phi/180.);
      cdw  =  cw

   !               ----- Use surface bc data from file ----
   CASE(1) ! Heat budget on preprocess mode - shortwave radiative and 
           ! net heat fluxes (including longwave & sensible & latent) as input 
           ! Space uniform & time varying

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF

     !.....Interpolate heat & momentum flux vars. to time n .....
     dthrs_surfbc = dtSurfbc/3600.
     eta = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
     Qsw = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
     Qn  = parab(0.,thrs,surfbc1(3,:),dthrs_surfbc) 
     cdw = parab(0.,thrs,surfbc1(4,:),dthrs_surfbc)
     uair= parab(0.,thrs,surfbc1(5,:),dthrs_surfbc)
     vair= parab(0.,thrs,surfbc1(6,:),dthrs_surfbc) 

     ! ... Calculate 3D-spatially variable sources 
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   CASE (2) ! Heat budget on run-time mode (I) - shortwave fluxes as input; 
            ! longwave & latent & sensible heat fluxes calculated.
            ! space uniform & time varying 

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF

     !.....Interpolate heat & momentum flux vars. values to present time step .....
     dthrs_surfbc = dtSurfbc/3600.
     eta = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
     Qsw = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
     Ta  = parab(0.,thrs,surfbc1(3,:),dthrs_surfbc)  
     Pa  = parab(0.,thrs,surfbc1(4,:),dthrs_surfbc)
     Rh  = parab(0.,thrs,surfbc1(5,:),dthrs_surfbc)
     Cc  = parab(0.,thrs,surfbc1(6,:),dthrs_surfbc) 
     cdw = parab(0.,thrs,surfbc1(7,:),dthrs_surfbc)
     uair= parab(0.,thrs,surfbc1(8,:),dthrs_surfbc)
     vair= parab(0.,thrs,surfbc1(9,:),dthrs_surfbc) 

     ! ... Calculate 3D-spatially variable heat sources
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   CASE (3) ! Heat budget on run-time mode (II) - radiative (short & longwave)
            ! fluxes as input; latent & sensible heat fluxes calculated. 
            ! Space uniform & time varying 

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF

     !.....Interpolate heat & momentum flux vars. values to present time step .....
     dthrs_surfbc = dtSurfbc/3600.
     eta = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
     Qsw = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
     Ta  = parab(0.,thrs,surfbc1(3,:),dthrs_surfbc)  
     Pa  = parab(0.,thrs,surfbc1(4,:),dthrs_surfbc)
     Rh  = parab(0.,thrs,surfbc1(5,:),dthrs_surfbc)
     Qlw = parab(0.,thrs,surfbc1(6,:),dthrs_surfbc) 
     cdw = parab(0.,thrs,surfbc1(7,:),dthrs_surfbc)
     uair= parab(0.,thrs,surfbc1(8,:),dthrs_surfbc)
     vair= parab(0.,thrs,surfbc1(9,:),dthrs_surfbc) 

     ! ... Calculate 3D-spatially variable heat sources
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   CASE (10) ! Spatially & time varying surface BC - Heat Budget calculated
              ! on RUN-TIME (I) mode

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF
     
     !.....Interpolate heat & momentum flux vars. values to present time step .....
     DO imet = 1, nmetstat
       dthrs_surfbc  = dtSurfbc/3600.
       eta           = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
       Pa            = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
       Qsw2D (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 3,:),dthrs_surfbc)
       Ta2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 4,:),dthrs_surfbc) 
       RH2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 5,:),dthrs_surfbc) 
       Cc2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 6,:),dthrs_surfbc) 
       uair2D(imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 7,:),dthrs_surfbc)
       vair2D(imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 8,:),dthrs_surfbc) 
     ENDDO

     ! ... Distribute heat and momentum sources entering through free surface
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   CASE (11) ! Spatially & time varying surface BC - Heat Budget calculated
              ! on RUN-TIME (II) mode

     !.....Return from subroutine on trapezoidal steps (except if n=1).....
     IF (n > 1) THEN; IF (istep == 2) RETURN; END IF
     
     !.....Interpolate heat & momentum flux vars. values to present time step .....
     DO imet = 1, nmetstat
       dthrs_surfbc  = dtSurfbc/3600.
       eta           = parab(0.,thrs,surfbc1(1,:),dthrs_surfbc)
       Pa            = parab(0.,thrs,surfbc1(2,:),dthrs_surfbc)
       Qsw2D (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 3,:),dthrs_surfbc)
       Ta2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 4,:),dthrs_surfbc) 
       RH2D  (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 5,:),dthrs_surfbc) 
       Qlw2D (imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 6,:),dthrs_surfbc) 
       uair2D(imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 7,:),dthrs_surfbc)
       vair2D(imet)  = parab(0.,thrs,surfbc1((imet-1)*6 + 8,:),dthrs_surfbc) 
     ENDDO

     ! ... Distribute heat and momentum sources entering through free surface
     CALL DistributeQsw
     CALL DistributeMomentumHeatSources

   END SELECT

END SUBROUTINE surfbc

!************************************************************************
 SUBROUTINE DistributeQsw
!************************************************************************
!
!  Purpose: To apportion the solar irradiance penetrating the 
!           lake through the surface among the layers in each
!           water column
!
!------------------------------------------------------------------------

  ! ... Local variables
  INTEGER:: i, j, k , l, kb, nwlayers, k1s, kms
  REAL   :: remFr, zfromt0
  REAL, DIMENSION (im1,jm1) :: htot
  REAL, DIMENSION (1  :km1) :: zfromt

  !.....Sweep over wet pressure points
  DO l = 1, lm

       ! ... Map 2D-l into 3D(i,j) indexes .....
       i = l2i(l); j = l2j(l);

       ! ... Define top & bottom wet layers.....
       kms = kmz(i,j);
       k1s = k1z(i,j);    
       nwlayers = (kms-k1s) + 1

       SELECT CASE (nwlayers)

       CASE (1)

         QswFr(k1s,l) = 1.

       CASE (2:)

         ! ... Compute the array of vertical distances from the 
         !     free surface to the top of each layer
         zfromt(k1s) = 0.0
         DO k = k1s+1, kms+1          
           zfromt(k)    = zfromt(k-1) + hp(k-1,l)
           QswFr(k-1,l) = SolarFr(zfromt(k-1)) - SolarFr(zfromt(k))
         END DO
         ! ... Redistribute Qsw not absorbed once it reaches the bottom 
         !     Probably it is not needed in deep lakes but could account for
         !     overheating for shallower lakes. 
         remFr = 1. - SUM ( QswFr (k1s:kms,l) ) 
         IF ( remFr > 1.e-6 ) THEN
           QswFr (k1s:kms,l) = QswFr (k1s:kms,l) + remFr / nwlayers
         END IF

       END SELECT

    END DO

END SUBROUTINE DistributeQsw

!************************************************************************
REAL FUNCTION SolarFr ( depth ) 
!************************************************************************
!
!  Purpose: To calculate the attenuation of solar irradiance penetrating  
!           the water column through the free surface. It uses formulation
!           proposed in Henderson-Sellers' Engineering Limnology Eq. 2.25 
!           According to the authors this equation is only valid for 
!           eta (attenuation coefficient) larger than 0.1
!
!------------------------------------------------------------------------

  ! ... Arguments
  REAL, INTENT (IN) :: depth

  ! ... Local variables
  REAL            :: BetaSol, z
  REAL, PARAMETER :: zA = 0.60
  
  z = depth
  IF ( eta .GE. 0.1) THEN
    BetaSol = 0.265*LOG (eta)+0.614; 
    IF ( z<zA ) THEN
      SolarFr = (1.-BetaSol*z/zA)
    ELSE
      SolarFr = (1.-BetaSol)*EXP(-eta*(z-zA)); 
    ENDIF
  ELSE
    SolarFr = EXP(-eta*z)
  ENDIF

END FUNCTION SolarFr

!************************************************************************
 SUBROUTINE DistributeMomentumHeatSources
!************************************************************************
!
!  Purpose: To construct a 2D met field from discrete variables & 
!           to contruct heat sources for each computational cell - Include
!           corrections to Barnes Method (Steve Andrews)
!
!------------------------------------------------------------------------

  ! ... Local variables
  INTEGER         :: i, j, k,l, k1s, kms, irg, imet, ii, jj, nF
  REAL            :: esw, ea, Lv, EmAir, Qbri, Qbra
  REAL            :: coeffE, coeffH, ws, wsi, ws0 
  REAL            :: sumw, uairdum, vairdum              ! SWA
  REAL            :: uresid(nmetstat), vresid(nmetstat)  ! SWA

  ! ... Local arrays containing parameters
  INTEGER, DIMENSION (6) :: range 
  REAL,    DIMENSION (5) :: a_d, b_d, p_d
  REAL,    DIMENSION (5) :: a_h, b_h, p_h, c_h
  REAL,    DIMENSION (5) :: a_e, b_e, p_e, c_e
  
  ! .... Parameter definition
  REAL, PARAMETER :: EmWater = 0.97
  REAL, PARAMETER :: StephanBoltzman = 5.6697e-8
  REAL, PARAMETER :: Al = 0.03
  REAL, PARAMETER :: Cp = 4181.6
  REAL, PARAMETER :: SpecificHeatAir = 1012.0

  ! ... Define parameters on first call
  IF ( n == 1) THEN
    range = (/  0.000, 2.2000, 5.0000, 8.0000, 25.00, 50.0 /)
    a_d   = (/  0.000, 0.7710, 0.8670, 1.2000, 0.000 /)
    b_d   = (/  1.080, 0.0858, 0.0667, 0.0250, 0.073 /)
    p_d   = (/ -0.150, 1.0000, 1.0000, 1.0000, 1.000 /)
    a_h   = (/  0.000, 0.9270, 1.1500, 1.1700, 1.652 /)
    b_h   = (/  1.185, 0.0521, 0.0100, 0.0075,-0.017 /)
    c_h   = (/  0.000, 0.0000, 0.0000,-4.5E-4, 0.000 /)
    p_h   = (/ -0.157, 1.0000, 1.0000, 1.0000, 1.000 /)
    a_e   = (/  0.000, 0.9690, 1.1800, 1.1960, 1.680 /)
    b_e   = (/  1.230, 0.0521, 0.0100, 0.0080,-0.016 /)
    c_e   = (/  0.000, 0.0000, 0.0000,-4.0E-4, 0.000 /)
    p_e   = (/ -0.160, 1.0000, 1.0000, 1.0000, 1.000 /)
  ENDIF

  ! ... Initialize Heat Source
  HeatSource = 0.0E0

  SELECT CASE (ifSurfbc)

  CASE (1) ! Heat budget on PRE-PROCESS mode - spatially uniform conditions

    DO l = 1,lm ;
 
      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define top & bottom wet layers........................
      k1s = k1z(i,j);    
      kms = kmz(i,j);

      ! ... Add non-penetrative components to surface layer ......
      HeatSource(k1s,l)=HeatSource(k1s,l)+ (Qn-Qsw)
 
      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) + & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
       IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (2) ! Heat budget on run-time mode (I) - spatially uniform conditions

    ! ... Define variables used in cals. and equal to all surface cells
    ea    = saturated_vapor_pressure (Ta+273.) * Rh
    EmAir = 0.642 * ( ea / (Ta+273.) )**0.14285714 		
    EmAir = EmAir * ( 1. + 0.17 * Cc**2. ) 
    Qbri  = EmAir * StephanBoltzman * (Ta+273.) ** 4. 
    Qbra  = Qbri * ( 1. - Al ) 

    ! ... Loop over surface cells 
    DO l = 1, lm;

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define bulk aerodynamic coefficients for neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ws = SQRT(uair(i,j)**2. + vair(i,j)**2.)    
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Define top & bottom wet layers........................
      k1s = k1z(i,j);    
      kms = kmz(i,j);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (3) ! Heat budget on run-time mode (II) - spatially uniform conditions

    ! ... Define variables used in cals. and equal to all surface cells
    ea    = saturated_vapor_pressure (Ta+273.) * Rh

    ! ... Loop over surface cells 
    DO l = 1, lm;

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Define bulk aerodynamic coefficients for neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ws = SQRT(uair(i,j)**2. + vair(i,j)**2.)    
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Define top & bottom wet layers........................
      k1s = k1z(i,j);    
      kms = kmz(i,j);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation (incoming LW as measured) 
      HeatSource(k1s,l)  = Qlw * (1.- Al) -                      &
                           EmWater * StephanBoltzman *           &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO

  CASE (10) ! Interpolation

    ! ... Loop over surface cells 
    DO l = 1, lm;

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Interpolate met records to grid points
      Qsw = 0.0; Ta = 0.0; RH = 0.0; Cc = 0.0; uair(i,j) = 0.0; vair(i,j) = 0.0;
      DO imet = 1, nmetstat
        Qsw  = Qsw    + WeightSt(i,j,imet) * Qsw2D (imet)
        Ta   = Ta     + WeightSt(i,j,imet) * Ta2D  (imet)
        RH   = RH     + WeightSt(i,j,imet) * RH2D  (imet)
        Cc   = Cc     + WeightSt(i,j,imet) * Cc2D  (imet)
        uair(i,j) = uair(i,j) + WeightSt(i,j,imet) * uair2D(imet)
        vair(i,j) = vair(i,j) + WeightSt(i,j,imet) * vair2D(imet)
      ENDDO 

      ! ... Wind speed over the water column
      ws  = SQRT (uair(i,j)**2.+vair(i,j)**2.);
        
      ! ... Define bulk aerodynamic coefficients for heat transfer under neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ! ... Define range of wind speed 
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Calculate drag coefficient (Amorocho & DeVries)  
      cdw(i,j) = 0.0015 * 1./(1.0+ EXP((12.5-ws)/1.56))+0.00104      

      ! ... Define variables used in cals. and equal to all surface cells
      ea    = saturated_vapor_pressure (Ta+273.) * Rh
      EmAir = 0.642 * ( ea / (Ta+273.) )**0.14285714 		
      EmAir = EmAir * ( 1. + 0.17 * Cc**2. ) 
      Qbri  = EmAir * StephanBoltzman * (Ta+273.) ** 4. 
      Qbra  = Qbri * ( 1. - Al ) 

      ! ... Define top & bottom wet layers........................
      k1s = k1z(i,j);    
      kms = kmz(i,j);

      ! ... Add non-penetrative components to surface cells
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO
  
  CASE (11) ! Interpolation

    ! Calculate residuals from first Barnes pass at met station points (SWA)
    DO k = 1, nmetstat
      i = NINT(metxy((k-1)*2+1))
      j = NINT(metxy((k-1)*2+2))
      sumw=SUM(Weightst(i,j,:))
      uairdum=0.0; vairdum=0.0
      DO imet = 1, nmetstat
        uairdum=uairdum + WeightSt(i,j,imet) / sumw * uair2D(imet)
        vairdum=vairdum + WeightSt(i,j,imet) / sumw * vair2D(imet)
      ENDDO
      uresid(k)=uair2D(k)-uairdum
      vresid(k)=vair2D(k)-vairdum
    ENDDO

    ! ... Interpolate other variables & compute heat fluxes
    DO l = 1, lm;

      ! ... Map 2D-l into 3D(i,j) indexes ........................
      i = l2i(l); j = l2j(l);

      ! ... Interpolate met records to grid points
      Qsw = 0.0; Ta = 0.0; RH = 0.0; Qlw = 0.0; 
      uair(i,j) = 0.0E0; vair(i,j) = 0.0E0; 
      DO imet = 1, nmetstat
        !Qsw  = Qsw    + WeightSt(i,j,imet) * Qsw2D (imet) ! SWA20101127
        !Ta   = Ta     + WeightSt(i,j,imet) * Ta2D  (imet) ! SWA20101127 
        !RH   = RH     + WeightSt(i,j,imet) * RH2D  (imet) ! SWA20101127
        !Qlw  = Qlw    + WeightSt(i,j,imet) * Qlw2D (imet) ! SWA20101127
        !uair(i,j) = uair(i,j) + WeightSt(i,j,imet) * uair2D(imet) ! SWA20101127
        !vair(i,j) = vair(i,j) + WeightSt(i,j,imet) * vair2D(imet) ! SWA20101127
        Qsw  = Qsw    + WeightSt(i,j,imet) / sumw * Qsw2D (imet) ! SWA20101127 
        Ta   = Ta     + WeightSt(i,j,imet) / sumw * Ta2D  (imet) ! SWA20101127
        RH   = RH     + WeightSt(i,j,imet) / sumw * RH2D  (imet) ! SWA20101127
        Qlw  = Qlw    + WeightSt(i,j,imet) / sumw * Qlw2D (imet) ! SWA20101127
        uair(i,j) = uair(i,j) + WeightSt(i,j,imet) / sumw * uair2D(imet) ! SWA20101127
        vair(i,j) = vair(i,j) + WeightSt(i,j,imet) / sumw * vair2D(imet) ! SWA20101127
      ENDDO

      ! Second pass Barnes interpolation for winds
      sumw=sum(WeightSt(i,j,:)**(1.0/gammaB))
      DO imet = 1, nmetstat
        uair(i,j) = uair(i,j) + WeightSt(i,j,imet)**(1.0/gammaB)/sumw * uresid(imet)
        vair(i,j) = vair(i,j) + WeightSt(i,j,imet)**(1.0/gammaB)/sumw * vresid(imet)
      ENDDO

      ! ... Wind speed over the water column  
      ws  = SQRT (uair(i,j)**2.+vair(i,j)**2.);

      ! ... Define bulk aerodynamic coefficients for heat transfer under neutral conditions 
      !     based on Kondo, 1975. Air-sea bulk transfer coefficients in diabatic
      !     conditions. In Boundary-Layer Meteorology 9 (1975) 91-112
      ! ... Define range of wind speed 
      IF      (ws>range(1).AND.ws<=range(2)) THEN; irg = 1
      ELSE IF (ws>range(2).AND.ws<=range(3)) THEN; irg = 2
      ELSE IF (ws>range(3).AND.ws<=range(4)) THEN; irg = 3
      ELSE IF (ws>range(4).AND.ws<=range(5)) THEN; irg = 4
      ELSE                                       ; irg = 5; END IF
      coeffH = (a_h(irg)+b_h(irg)*ws**p_h(irg)+c_h(irg)*(ws-8.)**2.)*1.e-3
      coeffE = (a_e(irg)+b_e(irg)*ws**p_e(irg)+c_e(irg)*(ws-8.)**2.)*1.e-3

      ! ... Calculate drag coefficient (Amorocho & DeVries)  
      cdw(i,j) = 0.0015 * 1./(1.0+ EXP((12.5-ws)/1.56))+0.00104      

      ! ... Define variables used in cals. and equal to all surface cells
      Qbra  = Qlw * ( 1. - Al ) 

      ! ... Define top & bottom wet layers........................
      k1s = k1z(i,j);    
      kms = kmz(i,j);

      ! ... Add non-penetrative components to surface cells
      ea  = saturated_vapor_pressure(Ta         +273.) * Rh
      esw = saturated_vapor_pressure(salp(k1s,l)+273.)
      Lv  = latent_heat_vaporization(salp(k1s,l)+273.)
      ! a. Long wave radiation 
      HeatSource(k1s,l)  = Qbra - EmWater * StephanBoltzman *    &
                           (salp(k1s,l)+273.) ** 4.
      ! b. Latent & Sensible heat fluxes
      HeatSource(k1s,l) = HeatSource(k1s,l) -                    &
      & rhoair*             Lv * coeffE*(esw-ea)*0.622/Pa * ws - & 
      & rhoair*SpecificHeatAir * coeffH*(salp(k1s,l)-Ta ) * ws

      ! ... Add penetrative components to water column & 
      !     express heat sources in temperature units ............
      DO k = k1s, kms
        HeatSource(k,l) = (HeatSource(k,l) +                     & 
          Qsw*QswFr(k,l))/((rhop(k,l)+1000.)*Cp)
        IF (Qsw<-100.) THEN; HeatSource(k,l) = 0.0E0; ENDIF
      END DO

    ENDDO
  
  END SELECT

END SUBROUTINE DistributeMomentumHeatSources

!************************************************************************
REAL FUNCTION saturated_vapor_pressure (T)
!************************************************************************

  IMPLICIT NONE
  REAL, INTENT (IN) :: T

  saturated_vapor_pressure = 2.1718e10 * EXP( -4157. / (T - 33.91) ) 

  END FUNCTION saturated_vapor_pressure 

!************************************************************************
REAL FUNCTION latent_heat_vaporization (T)
!************************************************************************

  IMPLICIT NONE
  REAL, INTENT (IN) :: T

  latent_heat_vaporization = 1.91846e6 * ( T / (T - 33.91) ) ** 2. 

  END FUNCTION latent_heat_vaporization

!***********************************************************************
SUBROUTINE InitializeInterpolationMethods
!***********************************************************************
!
!  Purpose: Solves for the weighting parameters used in Barnes interpolation sch.
!
!-----------------------------------------------------------------------

  ! ... Local variables
  INTEGER, PARAMETER :: intmethod = 0
  INTEGER :: i, j, imet, istat
  REAL    :: xi, yi, di, da, sumw
  ! .... Variables and parameters used in multiquadric interpolation method
  REAL    :: dsq, dfacto
  REAL, DIMENSION (nmetstat) :: xob, yob, verrsq
  REAL, PARAMETER :: c=0.05,smoo=0.1

  SELECT CASE (intmethod) 

  CASE (0) ! Barnes interpolation 

  DO i = i1, im; DO j = j1, jm

    ! Find average distance from grip point to all met stations
    da = 0.0
    DO imet = 1, nmetstat
      xi = metxy((imet-1)*2+1) - FLOAT(i);  
      yi = metxy((imet-1)*2+2) - FLOAT(j);
      di = SQRT(xi**2.+yi**2.)
      da = da + di
    ENDDO
    da = da / nmetstat
        
    ! Estimate weight based on distance to met station imet & da
    sumw = 0.0E0
    DO imet = 1, nmetstat
      xi = metxy((imet-1)*2+1) - FLOAT(i);  
      yi = metxy((imet-1)*2+2) - FLOAT(j);
      di = SQRT(xi**2.+yi**2.)
      weightst(i,j,imet) = EXP(-4.60517018598809 * (di**2.) / (da**2.))
      sumw = sumw + weightst(i,j,imet)
    ENDDO
        
    ! Normalize the weight
    DO imet = 1, nmetstat
      weightst(i,j,imet)=weightst(i,j,imet)/sumw
    ENDDO

  ENDDO; ENDDO

  CASE (1) ! Inverse squared distance 

  DO i = i1, im; DO j = j1, jm

    sumw = 0.0E0
    DO imet = 1, nmetstat
       xi = metxy((imet-1)*2+1) - FLOAT(i);  
       yi = metxy((imet-1)*2+2) - FLOAT(j);
       di = (xi**2.+yi**2.)
       weightst(i,j,imet) = MIN(1./di,1.E20)
       sumw = sumw + weightst(i,j,imet)
    ENDDO
        
    ! Normalize the weight
    DO imet = 1, nmetstat
      weightst(i,j,imet)=weightst(i,j,imet)/sumw
    ENDDO

  ENDDO; ENDDO

  CASE(2) ! ... Multiquadric interpolation 

    ! ... Allocate space for arrays used only in MQ interpolation
    ALLOCATE ( indice (nmetstat), mqQij(nmetstat, nmetstat), STAT=istat )
    IF (istat /= 0) CALL allocate_error ( istat, 120 )

    ! Allocate space for variables that store individual met records
    ! for a given time step for each of the variables - that will be 
    ! used for interpolation - This allows for smoothing while the
    ! Barnes method does not.
    ALLOCATE ( QswMQ (nmetstat  ), QlwMQ (nmetstat),           &
               RHMQ  (nmetstat  ), TaMQ  (nmetstat),           &
               uairMQ(nmetstat  ), vairMQ(nmetstat), STAT=istat )
    IF (istat /= 0) CALL allocate_error ( istat, 121 )

    ! ... Define verrsq to give more weight to the masurement station at the center
    verrsq = 0.1
    !verrsq(4) = 0.001;

    ! ... factor for normalizing distances in grid units
    dfacto=1.0/(FLOAT(MAX(im,jm))-1.)

    ! ... Initialize values of x,y meteorological observations
    DO imet = 1, nmetstat
       xob(imet) = metxy((imet-1)*2+1) * dfacto;  
       yob(imet) = metxy((imet-1)*2+2) * dfacto;
    END DO

    ! ... Fill the coordinate matrix weightst  
    DO i=1,nmetstat; DO j=i,nmetstat
       dsq=(xob(i)-xob(j))**2.+(yob(i)-yob(j))**2.
       mqQij(i,j)= - SQRT ( 1.0 + dsq / (c ** 2. ) )
       IF (i == j) THEN
          mqQij(i,i)=mqQij(i,j)+FLOAT(nmetstat)*smoo*verrsq(i)
       ELSE
          mqQij(j,i)=mqQij(i,j)
       END IF
    END DO; END DO

    ! ... Interpolate to normalized grid with origin at 0,0
    DO i = 1,im1; DO j = 1, jm1
      DO imet = 1,nmetstat
        xi = metxy((imet-1)*2+1) - float(i); xi = xi * dfacto 
        yi = metxy((imet-1)*2+2) - float(j); yi = yi * dfacto
        di = xi**2.+yi**2.
        weightst(i,j,imet) = - SQRT ( 1.0 + di / (c ** 2. ) )
      ENDDO
  ENDDO; ENDDO
  
  END SELECT

END SUBROUTINE InitializeInterpolationMethods

!***********************************************************************
SUBROUTINE InitializeBarnesInterpolation
!***********************************************************************
!
!  Purpose: Solves for the weighting parameters used in Barnes 
!  interpolation sch. Implementation of Barnes Scheme of Steve Andrews.
!
!-----------------------------------------------------------------------

  ! ... Local variables
  INTEGER :: i, j, imet, istat
  REAL    :: xi, yi, di, da !, sumw
  REAL    :: x1,y1,x2,y2,dist1(nmetstat-1),kappa1,delNfactor
  INTEGER :: counter1

  gammaB=0.3
  delNfactor=1.0;

  ! Find average minimum distance of a met station to all other met stations
  ! Distance is in grid units
  da=0.0
  dist1(:)=1.0E5;
  DO imet = 1, nmetstat
    x1 = metxy((imet-1)*2+1)
    y1 = metxy((imet-1)*2+2)
    counter1=1
    DO i = 1, nmetstat
      IF (imet /= i) THEN
        x2 = metxy((i-1)*2+1)
        y2 = metxy((i-1)*2+2)
        dist1(counter1) = SQRT((x1-x2)**2.0+(y1-y2)**2.0)
        counter1=counter1+1
      ENDIF
    ENDDO
    da = da + MINVAL(dist1)
  ENDDO
  da = da / REAL(nmetstat)

  ! Calculate smoothing scale length (from Koch 1983)
  IF (gammaB==0.2) THEN
     kappa1 = 5.0515*4.0/pi**2.0 * (delNfactor*da)**2.0
  ELSEIF (gammaB==0.3) THEN
     kappa1 = 3.5132*4.0/pi**2.0 * (delNfactor*da)**2.0
  ELSEIF (gammaB==0.5) THEN
     kappa1 = 2.3838*4.0/pi**2.0 * (delNfactor*da)**2.0
  ELSEIF (gammaB==0.6) THEN
     kappa1 = 2.1153*4.0/pi**2.0 * (delNfactor*da)**2.0
  ELSEIF (gammaB==0.8) THEN
     kappa1 = 1.7826*4.0/pi**2.0 * (delNfactor*da)**2.0
  ENDIF

  DO i = i1, im; DO j = j1, jm

    ! Find average distance from grid point to all met stations
    !da = 0.0
    !DO imet = 1, nmetstat
    !  xi = metxy((imet-1)*2+1) - FLOAT(i);  
    !  yi = metxy((imet-1)*2+2) - FLOAT(j);
    !  di = SQRT(xi**2.+yi**2.)
    !  da = da + di
    !ENDDO
    !da = da / nmetstat

    ! Estimate weight based on distance to met station imet & da
    !sumw = 0.0E0
    DO imet = 1, nmetstat
      xi = metxy((imet-1)*2+1) - FLOAT(i);
      yi = metxy((imet-1)*2+2) - FLOAT(j);
      !di = SQRT(xi**2.+yi**2.)
      di = xi**2.+yi**2.
      !weightst(i,j,imet) = EXP(-4.60517018598809 * (di**2.) / (da**2.))
      weightst(i,j,imet) = EXP(-di / kappa1)
      !sumw = sumw + weightst(i,j,imet)
    ENDDO

    ! Normalize the weight
    !DO imet = 1, nmetstat
    !  weightst(i,j,imet)=weightst(i,j,imet)/sumw
    !ENDDO

  ENDDO; ENDDO

END SUBROUTINE InitializeBarnesInterpolation

!***********************************************************************
  SUBROUTINE EstimateMQMetVars 
!***********************************************************************
!
!  The routine implements a multiquadric interpolation (e.g. Nuss &
!  Titley, Mon. Weather Rev., 122, 1611-1631, July, 1994) to
!  get a regular gridded array of interpolated values from
!  irregularly spaced observations.  This version assumes that
!  the same sites are always used and that data are never missing.
!
!-----------------------------------------------------------------------

   ! ... Local variables
   REAL, DIMENSION ( nmetstat) :: alpha
   REAL                        :: dp
   INTEGER                     :: is
    
   ! ... Find out the decomposition of mqQij and store it on first call
   IF ( n == 0 ) THEN; CALL ludcmp(dp); ENDIF

   alpha = Qsw2D ; CALL lubksb(alpha); QswMQ  = alpha;
   alpha = Ta2D  ; CALL lubksb(alpha); TaMQ   = alpha;
   alpha = RH2D  ; CALL lubksb(alpha); RHMQ   = alpha;
   alpha = Qlw2D ; CALL lubksb(alpha); QlwMQ  = alpha;
   alpha = uair2D; CALL lubksb(alpha); uairMQ = alpha;
   alpha = vair2D; CALL lubksb(alpha); vairMQ = alpha;

   END SUBROUTINE EstimateMQmetVars

!***********************************************************************
  SUBROUTINE ludcmp(D)
!***********************************************************************
!
!  LU matrix decomposition solver from Numerical Recipes
!
!-----------------------------------------------------------------------

      INTEGER, PARAMETER :: NMAX=100
      REAL, PARAMETER    :: TINY=1.0e-20
      REAL, DIMENSION (NMAX) :: VV

      ! ... Local variables
      INTEGER :: NN, NP
      INTEGER :: i, j, k, imax
      REAL, INTENT (INOUT) :: D 
      REAL :: AAMAX, SUM, DUM
      REAL, ALLOCATABLE, DIMENSION (:,:) :: A
      INTEGER, ALLOCATABLE, DIMENSION (:) :: INDX

      ALLOCATE ( A ( nmetstat, nmetstat ) )
      ALLOCATE ( INDX ( nmetstat ) )  

      A = mqQij
      INDX = indice
      NN = nmetstat
      NP = nmetstat      

      D=1.
      DO 12 I=1,NN
        AAMAX=0.
        DO 11 J=1,NN
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,NN
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,NN
          SUM=A(I,J)
          IF (J>1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM>=AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J/=IMAX)THEN
          DO 17 K=1,NN
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J/=NN)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,NN
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(NN,NN)==0.)A(NN,NN)=TINY

  mqQij = A
  indice = INDX

  DEALLOCATE (A)
  DEALLOCATE (INDX)

  END SUBROUTINE ludcmp

!***********************************************************************
  SUBROUTINE LUBKSB(B)
!***********************************************************************
!
!  LU equation solver from Numerical Recipes
!
!-----------------------------------------------------------------------
    
      REAL, DIMENSION (:), INTENT  (INOUT) :: B 

      ! ... Local variables
      INTEGER :: NN, NP 
      INTEGER :: ii, i, j, ll 
      REAL :: SUM
      REAL   , ALLOCATABLE, DIMENSION (:,:) :: A
      INTEGER, ALLOCATABLE, DIMENSION (:  ) :: INDX

      ! ... Allocate space
      ALLOCATE ( A ( nmetstat, nmetstat ) )
      ALLOCATE ( INDX ( nmetstat ) )  

      ! ... Initialize local variables
      A = mqQij
      INDX = indice
      NN = nmetstat
      NP = nmetstat      

      II=0
      DO 12 I=1,NN
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II/=0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM/=0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=NN,1,-1
        SUM=B(I)
        IF(I<NN)THEN
          DO 13 J=I+1,NN
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE

  DEALLOCATE ( A )
  DEALLOCATE (INDX)

  END SUBROUTINE lubksb

!***********************************************************************
FUNCTION parab ( frstpt, x, fx, dx )
!***********************************************************************
!
!  Purpose: To interpolate parabolically between the functional values
!           within the array fx. 
!
!-----------------------------------------------------------------------

   REAL, DIMENSION(:), INTENT(IN) :: fx      ! Assumed-shape array
   REAL, INTENT(IN) :: frstpt, x, dx
   REAL :: parab
   REAL :: om, theta
   INTEGER :: m

   m = (x - frstpt)/dx
   om = m
   theta = (x - frstpt - om*dx) / dx
   IF (m == 0) THEN
      m = 2
      theta = theta - 1.0
   ELSE
      m = m + 1
   END IF
   parab=fx(m)+0.5*theta*(fx(m+1)-fx(m-1)+theta*(fx(m+1)+fx(m-1)-2.0*fx(m)))

END FUNCTION parab

!***********************************************************************
FUNCTION linear ( frstpt, x, fx, dx )
!***********************************************************************
!
!  Purpose: To interpolate linearly between the functional values
!           within the array fx. 
!
!-----------------------------------------------------------------------

   REAL, DIMENSION(:), INTENT(IN) :: fx      ! Assumed-shape array
   REAL, INTENT(IN) :: frstpt, x, dx
   REAL :: linear
   REAL :: om, theta
   INTEGER :: m

   m = FLOOR((x - frstpt)/dx)
   theta = (x - frstpt - m*dx) / dx
   linear=(1-theta)*fx(m+1)+theta*fx(m+2)


END FUNCTION linear

!************************************************************************
SUBROUTINE allocate_error ( istat, ierror_code )
!************************************************************************
!
!  Purpose: Prints messages regarding any errors encountered during
!           allocation of space for model arrays. Stops program after
!           messages.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!------------------------------------------------------------------------

   !.....Arguments.....
   INTEGER, INTENT(IN) :: istat, ierror_code

   !.....Print error messages....
   PRINT *, " Program could not allocate space for arrays"
   PRINT '(" istat= ", I5, "  error code=", I3)', istat, ierror_code
   PRINT *, "  "
   PRINT *, "  "
   PRINT *, " ****STOPPING si3d due to allocate error"
   STOP

END SUBROUTINE allocate_error


!************************************************************************
SUBROUTINE open_error ( string, iostat )
!************************************************************************
!
!  Purpose: Prints 'string' regarding any error encountered during the
!           opening of a file. Also prints the iostat error code. The
!           program is stopped after printing messages.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!------------------------------------------------------------------------

   !.....Arguments.....
   CHARACTER(LEN=*), INTENT(IN) :: string
   INTEGER, INTENT(IN) :: iostat

   !.....Print error messages....
   PRINT '(A)', string
   PRINT '(" The iostat error number is ", I5)', iostat
   PRINT *, "  "
   PRINT *, "  "
   PRINT *, " ****STOPPING si3d due to OPEN error"
   STOP

END SUBROUTINE open_error


!************************************************************************
SUBROUTINE input_error ( ios, ierror_code )
!************************************************************************
!
!  Purpose: Prints messages regarding any errors encountered during
!           reading of the input file. Stops program after messages.
!
!  Revisions:
!    Date            Programmer        Description of revision
!    ----            ----------        -----------------------
!
!------------------------------------------------------------------------

   !.....Arguments.....
   INTEGER, INTENT(IN) :: ios, ierror_code

   !.....Determine type of read error.....
   SELECT CASE (ios)

   !.....End of file (ios=-1).....
   CASE (:-1)
      PRINT *, " Unexpected end-of-file encountered while reading file"

   !.....Error during reading (ios>0).....
   CASE (1:)
      PRINT *, " Error during reading data"
   END SELECT

   !.....Print ios and error statement number.....
   PRINT '(" The iostat error number is ", I5)',     ios
   PRINT '(" Error on read statement number ", I3)', ierror_code
   PRINT *, "  "
   PRINT *, "  "
   PRINT *, " ****STOPPING si3d due to read error"
   STOP

END SUBROUTINE input_error


END MODULE si3d_BoundaryConditions
