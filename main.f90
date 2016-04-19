
PROGRAM FlatPlateParallel
    USE MPI
    USE Constants_Model
    IMPLICIT NONE
    INTEGER:: I, J, k
    ! for Cartesian grid information
    REAL(KIND=8), ALLOCATABLE:: MeshX( : ), MeshY( : )
    REAL(KIND=8), ALLOCATABLE:: DeltaX( : ), DeltaY( : )
    ! 
    ! non-uniform background flow field
    REAL(KIND=8), ALLOCATABLE:: u_0( :, : ), v_0( :, : ), rho_0( :, : ), c_0( :, : )
    REAL(KIND=8), ALLOCATABLE:: Irho_0( :, : ), LNrho_0( :, : )
    REAL(KIND=8), ALLOCATABLE:: u_0y( :, : ), Irho_0y( :, : ), LNrho_0y( :, : )
    ! 
    ! acoustic field: nonconservative variables and conservative variables
    ! output to file for main processor
    REAL(KIND=8), ALLOCATABLE:: TotalP( :, : ), TotalROU( :, : )
    REAL(KIND=8), ALLOCATABLE:: TotalU( :, : ), TotalV( :, : )
    REAL(KIND=8), ALLOCATABLE:: TotalW( :, : ), TotalOmiga( :, : )
    ! 
    ! calculation
    REAL(KIND=8), ALLOCATABLE:: U( :, : ), V( :, : )
    REAL(KIND=8), ALLOCATABLE:: Omiga( :, : )
    REAL(KIND=8), ALLOCATABLE:: P( :, : ), ROU( :, : )
    REAL(KIND=8), ALLOCATABLE:: Du( :, : ), Dv( :, : ), Dw( :, : )
    REAL(KIND=8), ALLOCATABLE:: Uold( :, : ), Vold( :, : )
    REAL(KIND=8), ALLOCATABLE:: Pold( :, : )
    REAL(KIND=8), ALLOCATABLE:: Q1( :, : ), Q2( :, : ), Q3( :, : ), Q4( :, : )
    REAL(KIND=8), ALLOCATABLE:: F1( :, : ), F2( :, : ), F3( :, : ), F4( :, : )
    REAL(KIND=8), ALLOCATABLE:: G1( :, : ), G2( :, : ), G3( :, : ), G4( :, : )
    REAL(KIND=8), ALLOCATABLE:: H1( :, : ), H2( :, : ), H3( :, : ), H4( :, : )
    ! 
    ! PML zone 
    REAL(KIND=8):: SigmaX1( -(NPML-1) : 0 )
    REAL(KIND=8):: SigmaX2( 1 : NPML )
    REAL(KIND=8):: SigmaY1( -(NPML-1) : 0 )
    REAL(KIND=8):: SigmaY2( 1 : NPML )
    REAL(KIND=8), ALLOCATABLE:: Au11( :, : ), Au12( :, : )
    REAL(KIND=8), ALLOCATABLE:: Au13( :, : ), Au14( :, : )
    REAL(KIND=8), ALLOCATABLE:: Au21( :, : ), Au22( :, : )
    REAL(KIND=8), ALLOCATABLE:: Au23( :, : ), Au24( :, : )
    REAL(KIND=8), ALLOCATABLE:: Au31( :, : ), Au32( :, : )
    REAL(KIND=8), ALLOCATABLE:: Au33( :, : ), Au34( :, : )
    REAL(KIND=8), ALLOCATABLE:: S1( :, : ), S2( :, : )
    REAL(KIND=8), ALLOCATABLE:: S3( :, : ), S4( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu11( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu12( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu13( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu14( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu21( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu22( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu23( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu24( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu31( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu32( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu33( :, : )
    REAL(KIND=8), ALLOCATABLE:: SAu34( :, : )
    !!
    !MPI variables definition
    INTEGER:: IERR, NUMPROCS
    INTEGER:: MYID, MYROOT
    INTEGER:: MYLEFT, MYRIGHT, MYUPPER, MYLOWER
    INTEGER:: PX, PY
    INTEGER:: HTYPE, VTYPE
    INTEGER:: XTYPE, YTYPE, ZTYPE
    INTEGER, ALLOCATABLE:: STATUS( :, : ), REQ(:)
    INTEGER:: COUNT
    INTEGER,ALLOCATABLE:: BLOCKLENS( : ), INDICES( : )
    INTEGER:: SENDCNT
    INTEGER, ALLOCATABLE:: RECVCNT( : )
    INTEGER, ALLOCATABLE:: DISPLS( : )
    ! 
    ! grid number for each processor( MPI )
    ! NA: ghost point for MPI information exchange
    ! 		Here using DRP scheme, NA = 3.
    INTEGER:: XN1, YN2, ZN3
    ! 
    ! coordinate transformation
    REAL(KIND=8), ALLOCATABLE:: Jacobi( :, : )
    REAL(KIND=8), ALLOCATABLE:: KexiX( : ), EitaY( : ), TaoZ( : )
    ! 
    ! loop control variables
    INTEGER:: Tstep0, K0, MaxTimeStep, Loops, Flag
    REAL(KIND=8):: MaxResidual, MaxResidual0
    ! 
    ! output flag
    INTEGER:: IPX, IPY, IPZ
    INTEGER:: SizeN1, SizeN2, SizeN3
    INTEGER:: OUTPUTFlag, FLAG0, SIZE0
    CHARACTER(LEN=80):: FilenameINPUT, FilenameOUTPUT, FilenameFile
!    INTEGER, ALLOCATABLE:: ProbeFlag( : ), ProbeFlagRoot( : )
!    REAL(KIND=8), ALLOCATABLE:: ProbePressure0( : )
!    REAL(KIND=8), ALLOCATABLE:: ProbePressure1( : )
!    REAL(KIND=8), ALLOCATABLE:: ProbePressure( : )
    REAL(KIND=8):: STARTTIME, ENDTIME
    !!
    !step- 1: initialize the parallel processors
    CALL MPI_INIT( IERR )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS, IERR )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
    ! 
    MYROOT = 0
    ! 
    ! for data exchanger
    ALLOCATE( RECVCNT( NUMPROCS ) )
    ALLOCATE( DISPLS( NUMPROCS ) )
    ! 
    ! for series I/O: output results( MPI_ISEND, MPI_IRECV )
    IF ( MYID == MYROOT ) THEN
        SIZE0 = NUMPROCS - 1
        ALLOCATE( REQ( SIZE0 ) )
        ALLOCATE( STATUS( MPI_STATUS_SIZE, SIZE0 ) )
    ELSE
        SIZE0 = 1
        ALLOCATE( REQ( SIZE0:SIZE0 ) )
        ALLOCATE( STATUS( MPI_STATUS_SIZE, SIZE0 ) )
    END IF
    ! 
    IF ( ( NPX * NPY * NPZ ) .NE. NUMPROCS ) THEN
        WRITE( *, * ) 'The processor grid distribution doesn''t ', &
    &          ' conform with the processor summation! '
        STOP
        CALL MPI_FINALIZE( IERR )
    END IF
    !!
    !get the current processor coordinates( PX, PY, PZ )
    PZ = MYID / ( NPX * NPY )
    PY = ( MYID - PZ * NPX * NPY ) / ( NPX )
    PX = ( MYID - PZ * NPX * NPY ) - PY * NPX
    !determine the topology relationship of all the processors
    !!
    MYLEFT = MYID - 1
    IF ( MOD( MYID, NPX ) .EQ. 0 ) MYLEFT = MPI_PROC_NULL
    MYRIGHT = MYID + 1
    IF ( MOD( MYRIGHT, NPX ) .EQ. 0 ) MYRIGHT = MPI_PROC_NULL
    MYREAR = MYID + NPX
    IF ( MYREAR .GE. ( PZ+1 )*NPX*NPY ) MYREAR = MPI_PROC_NULL
    MYFORWARD = MYID - NPX
    IF ( MYFORWARD .LT. PZ*NPX*NPY ) MYFORWARD = MPI_PROC_NULL
    MYUPPER = MYID + NPX * NPY
    IF ( MYUPPER .GE. NUMPROCS ) MYUPPER = MPI_PROC_NULL
    MYLOWER = MYID - NPX * NPY
    IF( MYLOWER .LT. 0 ) MYLOWER = MPI_PROC_NULL 
    !!
    ! 
    OUTPUTFlag = 0
    FilenameINPUT = './INPUT/'
    FilenameOUTPUT = './OUTPUT/'
    FilenameFile = 'Cylinder.dat'
    DO I = 1, NumProbe
        ProbeAngle( I ) = 2.0d0 * PI * I / NumProbe
    END DO
    ! 
    !step- 2: read the number for surface mesh
    IF ( MYID == MYROOT ) THEN
        OPEN( UNIT = 10, FILE = TRIM( TRIM( FilenameINPUT ) // TRIM( FilenameFile ) ) )
        READ( UNIT = 10, FMT = * ) LN
        CLOSE( UNIT = 10 )
    END IF
    !!
    !bcast LN to every processor
    COUNT = 1
    CALL MPI_BCAST( LN, COUNT, MPI_INTEGER, MYROOT, &
    &       MPI_COMM_WORLD, IERR )
    !!
    !step- 3: allocate array variables for immersed body
    !   allocate variable related to the immersed boundary
    ALLOCATE( EffectRange( LN, 6 ), NormalVector( LN, 3 ) )
    ALLOCATE( Shapes( LN, 3 ), ShapesForP( LN, 3 ), Cell_S( LN ) )
    ALLOCATE( AV( LN, LN ), AT( LN, LN ) )
    ALLOCATE( LagForce( LN, 3 ), SurfUVW( LN, 3 ) )
    PLN = CEILING( 1.0d0 * LN / NUMPROCS )
    ALLOCATE( AProcessor( LN, PLN ) )
    ALLOCATE( ATotal( LN, NUMPROCS * PLN ) )
    ALLOCATE( DuDt0( LN ) )
    ALLOCATE( DuDt( LN ) )
    ALLOCATE( DuDt1( LN ) )
    ALLOCATE( DuDtFlag( LN+1 ) )
    ALLOCATE( DuDtFlagRoot( LN ) )
    ! 
!    ALLOCATE( ProbePressure0( NumProbe ) )
!    ALLOCATE( ProbePressure1( NumProbe ) )
!    ALLOCATE( ProbePressure( NumProbe ) )
!    ALLOCATE( ProbeFlag( NumProbe+1 ) )
!    ALLOCATE( ProbeFlagRoot( NumProbe ) )
    !!
    !read the original body coordinates for the main processor
    IF ( MYID == 0 ) THEN
        OPEN( UNIT = 10, FILE = TRIM( TRIM( FilenameINPUT ) // TRIM( FilenameFile ) ) )
        READ( UNIT = 10, FMT = * ) LN
        IF ( ModelFlag == 0 ) THEN
            DO I = 1, LN
                READ( UNIT = 10, FMT = * ) shapes( I, : )
!                 READ( UNIT=10, FMT=* ) TempShapes(1,1:4)
!                 Shapes(I,:) = TempShapes(1,1:3)
!                 CELL_S(I) = TempShapes(1,4)
            END DO
            Shapes( : , 3 ) = Centre( 3 )
        ELSE
            DO I = 1, LN
                READ( UNIT = 10, FMT = * ) shapes( I, : )
!                 READ( UNIT=10, FMT = * ) TempShapes(1,1:4) 
!                 Shapes(I,:) = TempShapes( 1,1:3 )
!                 CELL_S(I) = TempShapes(1,4)
            END DO
        END IF
        CLOSE( UNIT = 10 )
    END IF
    ! 
    ! bcast the shapes to other processors
    !bcast LN to every processor
    COUNT = LN*3
    CALL MPI_BCAST( Shapes(1,1), COUNT, MPI_DOUBLE_PRECISION, MYROOT, &
    &       MPI_COMM_WORLD, IERR )
    !  f
    !step- 4: input the analysis model
    !IF ( MYID == MYROOT ) THEN
    CALL Models( LRef, R_PML0, Shapes, NormalVector, Cell_S, LN, AOA, &
    &       Centre, ModelFlag )
    !END IF
    !!
    !step-5: get the size for Cartesian mesh
    CALL GetMeshSize( N1, N2, N3, StartEnd, CutPoints, DeltaXYZ, &
    &           R_PML, R_PML0, Shapes, Cell_S, LN, Centre, MaxX, MaxY, &
    &           MaxZ, SolidRelativePosi, q0, LRef, NPML, NPX, NPY, NPZ, &
    &           ModelFlag, MYID, FilenameOUTPUT )
    !!
!    IF ( MYID == MYROOT ) WRITE( *, * ) N1, N2, N3
    ! 
    !get the mesh number for each processor
    XN1 = ( 2 * NPML + N1 ) / NPX
    YN2 = ( 2 * NPML + N2 ) / NPY
    IF ( ModelFlag == 0 ) THEN
        ! 2-D Model
        ZN3 = 1
    ELSE
        ! 3-D Model
        ZN3 = ( 2 * NPML + N3 ) / NPZ
    END IF
    !!
    !-----------------------------------------------------------------------
    ! 					-	M	A	I	N								****
    ! 						-	Z	O	N	E	-						****
    ! 							-	V	A	R	I	A	T	I	O	N-	****
    ! 					-	P	M	L	-								****
    ! 						-	Z	O	N	E	-						****
    !-----------------------------------------------------------------------
    !PML equation: 
    ! : Q_t+(F-F0)_x + ( G-G0 )_y + ( H-H0 )_z + sigma(x) * Au(1) + sigma(y)
    !	* Au(2) + sigma(z) * Au(3) + beita * sigma(x) * (F-F0) = 0  ---- (1)
    ! : Au(1)_t + sigma(x) * Au(1) + (F-F0)_x + beita * sigma(x) * (F-F0) 
    !	= 0 -------------------------------------------------------------(2)
    ! : Au(2)_t + sigma(y) * Au(2) + (G-G0)_y = 0 -----------------------(3)
    ! : Au(3)_t + sigma(z) * Au(3) + (H-H0)_z = 0 -----------------------(4)
    !-----------
    !-----------
    ! Because the governing equation doesn't include viscous terms, the 
    ! equations (1)~(4) only need to be solved.
    ! (Here can refer to the paper< Absorbing boundary conditions for nonl- 
    ! inear Euler and Navier-stokes equations based on the perfectly match-
    ! -ed layer techniques >. Fang Q. Hu, X.D Li, D.K Lin, 2008) for equation
    ! (38)~(44)
    !----------
    !----------
    !step-6: allocate array variables for Cartesian mesh of each processor
    ALLOCATE( MeshX( -(NPML-1) : (N1+NPML) ) )
    ALLOCATE( MeshY( -(NPML-1) : (N2+NPML) ) )
    ALLOCATE( MeshZ( -(NPML-1) : (N3+NPML) ) )
    ALLOCATE( DeltaX( -(NPML-1) : (N1+NPML-1) ) )
    ALLOCATE( DeltaY( -(NPML-1) : (N2+NPML-1) ) )
    ALLOCATE( DeltaZ( -(NPML-1) : (N3+NPML-1) ) )
    ALLOCATE( Jacobi( -(NPML-1):(N1+NPML), -(NPML-1):(N2+NPML), -(NPML-1):(N3+NPML) ) )
    ALLOCATE( KexiX( -(NPML-1) : ( N1 + NPML ) ) )
    ALLOCATE( EitaY( -(NPML-1) : ( N2 + NPML ) ) )
    ALLOCATE( TaoZ( -(NPML-1) : ( N3 + NPML ) ) )
    ALLOCATE( U0( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( V0( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( W0( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( P0( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( ROU0( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( C0( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( U( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( V( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( W( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Omiga( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( P( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( ROU( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Du( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Dv( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Dw( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Q1( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Q2( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Q3( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Q4( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F1( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F2( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F3( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( F4( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G1( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G2( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G3( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( G4( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H1( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H2( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H3( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( H4( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Au11( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au12( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au13( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au14( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au21( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au22( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au23( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au24( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au31( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au32( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au33( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Au34( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( S1( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( S2( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( S3( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( S4( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu11( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu12( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu13( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu14( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu21( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu22( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu23( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu24( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu31( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu32( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu33( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( SAu34( 1 : XN1, 1 : YN2, 1 : ZN3 ) )
    ALLOCATE( Uold( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Vold( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Wold( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ALLOCATE( Pold( -(NA-1) : (XN1+NA), -(NA-1) : (YN2+NA), -(NA-1) : (ZN3+NA) ) )
    ! 
    ! total variables for outputting to file in main processor
    IF ( MYID == MYROOT ) THEN
        SizeN1 = N1 + 2 * NPML
        SizeN2 = N2 + 2 * NPML
        IF ( ModelFlag == 0 ) THEN
            ! 2-D Model
            SizeN3 = 1
        ELSE
            ! 3-D Model
            SizeN3 = N3 + 2 * NPML
        END IF
        ALLOCATE( TotalP( SizeN1, SizeN2, SizeN3 ) )
        ALLOCATE( TotalRou( SizeN1, SizeN2, SizeN3 ) )
        ALLOCATE( TotalU( SizeN1, SizeN2, SizeN3 ) )
        ALLOCATE( TotalV( SizeN1, SizeN2, SizeN3 ) )
        ALLOCATE( TotalW( SizeN1, SizeN2, SizeN3 ) )
        ALLOCATE( TotalOmiga( SizeN1, SizeN2, SizeN3 ) )
    END IF
    ! 
    ! -----------------------------------------------------------------------------------
    ! Define new datatype for output operation and the MPI message transfer
    ! -----------------------------------------------------------------------------------
    ! 
    ! create new datatype to output variables to file from each processor
    ! new datatype for sending process
    COUNT = YN2 * ZN3
    ALLOCATE( BLOCKLENS( COUNT ) )
    ALLOCATE( INDICES( COUNT ) )
    BLOCKLENS = XN1
    DO J = 1, ZN3
        DO I = 1, YN2
            INDICES( YN2*(J-1)+I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO
    END DO
    CALL MPI_TYPE_INDEXED( COUNT, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, HTYPE, IERR )
    CALL MPI_TYPE_COMMIT( HTYPE, IERR )
    ! 
    ! new datattype for receving process
    DO J = 1, ZN3
        DO I = 1, YN2
            INDICES( YN2*(J-1)+I ) = (J-1) * SizeN1 * SizeN2 + (I-1) * SizeN1
        END DO
    END DO
    CALL MPI_TYPE_INDEXED( COUNT, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, VTYPE, IERR )
    CALL MPI_TYPE_COMMIT( VTYPE, IERR )
    !   
    ! create new datatype for data exchange on the boundary for each processor
    ! TYPE- 1: XTYPE ------for forward and rearward side data exchange
    DEALLOCATE( BLOCKLENS, INDICES )
    COUNT = NA * ZN3
    ALLOCATE( BLOCKLENS( COUNT ) )
    ALLOCATE( INDICES( COUNT ) )
    BLOCKLENS = XN1
    DO J = 1, ZN3
        DO I = 1, NA
            INDICES( NA*(J-1)+I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO 
    END DO
    CALL MPI_TYPE_INDEXED( COUNT, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, XTYPE, IERR )
    CALL MPI_TYPE_COMMIT( XTYPE, IERR )
    ! 
    ! TYPE- 2: YTYPE----------for left and right side data exchange
    DEALLOCATE( BLOCKLENS, INDICES )
    COUNT = YN2 * ZN3
    ALLOCATE( BLOCKLENS( COUNT ) )
    ALLOCATE( INDICES( COUNT ) )
    BLOCKLENS = NA
    DO J = 1, ZN3
        DO I = 1, YN2
            INDICES( YN2*(J-1)+I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO 
    END DO
    CALL MPI_TYPE_INDEXED( COUNT, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, YTYPE, IERR )
    CALL MPI_TYPE_COMMIT( YTYPE, IERR )
    ! 
    ! TYPE- 3: ZTYPE ---------for upper and lower side data exchange
    DEALLOCATE( BLOCKLENS, INDICES )
    COUNT = YN2 * NA
    ALLOCATE( BLOCKLENS( COUNT ) )
    ALLOCATE( INDICES( COUNT ) )
    BLOCKLENS = XN1
    DO J = 1, NA
        DO I = 1, YN2
            INDICES( YN2*(J-1)+I ) = (J-1) * (XN1+2*NA) * (YN2+2*NA) + (I-1) * (XN1+2*NA)
        END DO
    END DO
    CALL MPI_TYPE_INDEXED( COUNT, BLOCKLENS, INDICES, MPI_DOUBLE_PRECISION, ZTYPE, IERR )
    CALL MPI_TYPE_COMMIT( ZTYPE, IERR )
    ! -----------------------------------------------------------------------------------
    ! ------------------------------ new datatype for MPI -------------------------------
    !----
    !-----------------------------------------------------------------------
    !                   -   M   A   I   N                               ****
    !                       -   Z   O   N   E   -                       ****
    !                           -   V   A   R   I   A   T   I   O   N-  ****
    !                   -   P   M   L   -                               ****
    !                       -   Z   O   N   E   -                       ****
    !-----------------------------------------------------------------------
    !!
    !step- 7: grid generation
    CALL Mesh( MeshX, MeshY, MeshZ, DeltaX, DeltaY, DeltaZ, N1, N2, N3, &
    &       NPML, StartEnd, CutPoints, DeltaXYZ, q0, ModelFlag, MYID, FilenameOUTPUT )
    ! 
    CALL NumericalProbe( Corner, SideLength, RelativePosi, ProbeCoordinate, NumProbe, &
    &           Rp, ProbeAngle, MeshX, MeshY, MeshZ, N1, N2, N3, NPML, ModelFlag )
    !!
    !step- 7.1: coordinate transformation/ get Jacobi matrix
    CALL TransformCoordinate( Jacobi, KexiX, EitaY, TaoZ, &
    &   MeshX, MeshY, MeshZ, N1, N2, N3, NPML, DeltaXYZ, ModelFlag )
    ! 
    !step- 7.2: get the absorbing coefficient SigmaX, SigmaY, SigmaZ for PML zone
    CALL AbsorbingCoefficient( Beita, SigmaX1, SigmaX2, SigmaY1, SigmaY2, &
    &           SigmaZ1, SigmaZ2, AFA, SIGMA0, NPML, MeshX, MeshY, MeshZ, &
    &                   N1, N2, N3, BeitaFlag, ModelFlag )
    !!
    !step- 8: get the background flow field
    CALL GetMeanFlow( U0, V0, W0, P0, ROU0, C0, NA, XN1, YN2, ZN3, &
    &    MeshX, MeshY, MeshZ, Centre, NPML, N1, N2, N3, PX, PY, PZ, &
    &       NPX, NPY, NPZ, ModelFlag )
    !!
    !step- 9: initialize acoutic field
    CALL InitializeAcousticField( U, V, W, P, ROU, Au11, Au12, Au13, Au14, &
    &       Au21, Au22, Au23, Au24, Au31, Au32, Au33, Au34, NA, XN1, YN2, &
    &           ZN3, MeshX, MeshY, MeshZ, C0, NPML, N1, N2, N3, PX, PY, PZ, &
    &                NPX, NPY, NPZ, ModelFlag )
    !!
    !step- 10: enter loops
    !!----------------------------Algorithm---------------------
    !here we solve conservative form of acoustic propagation
    !	governing equations, the equation group can be got in 
    !	< Acoustic Scattering in non-uniform flow >. Formula (7)
    !	Derived by Cheng Long( 2015 )
    !!----------------------------------------------------------
    DeltaT = DeltaXYZ * CFL
    Tstep0 = floor( 0.1d0 / DeltaT )
    K0 = 0
    MaxTimeStep = CEILING( TotalTime / DeltaT )
    IF ( MYID == MYROOT ) THEN
        OPEN( UNIT = 11, FILE = TRIM( TRIM( FilenameOUTPUT ) &
    &                   // TRIM( 'NumericalProbeAB.dat' ) ) )
    END IF
    Q1 = 0.0d0
    Q2 = 0.0d0
    Q3 = 0.0d0
    Q4 = 0.0d0
    F1 = 0.0d0
    F2 = 0.0d0
    F3 = 0.0d0
    F4 = 0.0d0
    G1 = 0.0d0
    G2 = 0.0d0
    G3 = 0.0d0
    G4 = 0.0d0
    H1 = 0.0d0
    H2 = 0.0d0
    H3 = 0.0d0
    H4 = 0.0d0