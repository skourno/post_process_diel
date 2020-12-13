!----------------------------------------------------------------------------------
!                               input.f90  -  description
!                             ------------------------------
!                  begin    : Fri 18.11.2016
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
!
! This is the standard input module where the input file "input_particles.txt" is read.
! It also contains the appropriate routine to read the trajectory and the topology (if 
! needed)
! The loading of the topology and the trajectory is done in "read_cases.f90", which is
! a general use module.
!
! Reading formats:
!   IRead (1|2) --> (lammps dump.dip | .traj file)
!   ITopo (0|1) --> (no connectivity | .topo file)
!
!   Reading happens like this: 
!     Read all frames until IStartFrame. Then process each frame until NFramesProcess 
!     is the total number of processed frames. You can use NFramesStride to process 
!     the frames with the desired step (default value NFramesStride = 1).
!
!----------------------------------------------------------------------------------
!
!  17.01.17   1.0   Skou    readInput, read_sys      
!  21.02.20   2.0   Skou    updated to be used in dielectric analysis code      
!
!----------------------------------------------------------------------------------
 Module Input
    Use nrtype
    Implicit None

    ! output file path and name
    Character    :: CFileTraj*240          ! Original trajectory.
    Character    :: CFileTopo*240          ! Input Topology File
    Character    :: CFileOutput*240        ! Output file. Depending on the user specified mode of execution 
                                           ! an appropriate suffix is added at the end. (see the write_xxxx
                                           ! modules for more info

    Integer(I4B) :: IRead                  ! Trajectory read format
    Integer(I4B) :: ITopo                  ! Topology read format (or 0 for absence)
    Integer(I4B) :: IStartFrame            ! Start process at this frame
    Integer(I4B) :: ILastFrame             ! Don't process frames with greater index than this one
    Integer(I4B) :: NFramesStride          ! Use this step to process frames  

    Integer(I4B) :: NOptions   = 4         ! Number of histogram options given to the user
    Integer(I4B) :: IOption(4) = 0         ! 0 means that the user does not wish for this histo to be generated
    Integer(I4B) :: NBins(4)   = 0         ! Number of bins for this histogram
    Real(DP)     :: BLeft(4)   = ZERO      ! Left histogram boundary
    Real(DP)     :: BRight(4)  = ZERO      ! Right histogram boundary


 contains

    !------------------------------------------------------------------------------
    !  Read input file. Check input variables for consistency
    !
    !------------------------------------------------------------------------------
    Subroutine ReadInput
        Implicit None
        Integer(I4B)   :: NLinesRead, IUn, i
        Integer(I4B)   :: NModesRequested = 0
        Character      :: CDummy*16

        NLinesRead = 0
        IUn        = NGetFileUnit()
        open(IUn, File='input_dielpp.txt', status='OLD', ERR=67)

        call ReadCommentLines(IUn, 3)
        NLinesRead  = NLinesRead + 3

        read(IUn, *, ERR=66) CDummy, CFileTraj    ;  read(IUn,*, ERR=66)  ;  NLinesRead = NLinesRead + 2
        read(IUn, *, ERR=66) CDummy, CFileTopo    ;  read(IUn,*, ERR=66)  ;  NLinesRead = NLinesRead + 2
        read(IUn, *, ERR=66) CDummy, CFileOutput  ;  read(IUn,*, ERR=66)  ;  NLinesRead = NLinesRead + 2

        !print*, TRIM(CFileTraj)
        !print*, TRIM(CFileTopo)
        !print*, TRIM(CFileOutput)

        call ReadCommentLines(IUn, 14)
        NLinesRead  = NLinesRead + 14

        read(IUn, *, ERR=66) IRead, ITopo, IStartFrame, ILastFrame, NFramesStride
        NLinesRead  = NLinesRead + 1

        if ( (IRead   < 1).OR.(IRead   > 2) ) STOP "ReadInput: ERROR - IRead format not supported"
        if ( (ITopo   < 1).OR.(ITopo   > 1) ) STOP "ReadInput: ERROR - ITopo not supported"
        
        if ( IStartFrame    <= 0)                 STOP "ReadInput: ERROR - IStartFrame loading instructions "
        if ( ILastFrame     <  IStartFrame)       STOP "ReadInput: ERROR - ILastFrame loading instructions"
        if ( NFramesStride  <= 0)                 STOP "ReadInput: ERROR - IFramesStride loading instructions"
        if ( MOD(ILastFrame,NFramesStride).NE.0 ) STOP "ReadInput: ERROR - ILastFrame needs to be a multiplier of NFramesStride"


        call ReadCommentLines(IUn, 6)
        NLinesRead  = NLinesRead + 6

        do i=1,NOptions
            read(IUn,*, ERR=66) IOption(i), NBins(i), BLeft(i), BRight(i)
            NLinesRead = NLinesRead + 1

            Call ReadCommentLines(IUn, 1)
            NLinesRead = NLinesRead + 1
        enddo

        close(IUn)

        return

        66 print*, 'ReadInput: Read Error - LinesRead ', NLinesRead ; close(IUn) ; STOP
        67 print*, 'ReadInput: Open  Problem'                       ; close(IUn) ; STOP
    End Subroutine ReadInput

    !------------------------------------------------------------------------------
    !  Read system. This routine reads 1 frame of the InputTraj file as specified 
    !  in the input file. For the first entry it also reads the topology (if ITopo
    !  .NE. 0).   
    !
    !------------------------------------------------------------------------------
    Subroutine Read_Sys
        Use System,      Only : Atom, AtomType, Cell, Topo 
        Use Read_Cases,  Only : Read_LAMMPS, Read_LAMMPS_dip, Read_Topology
        Implicit None
        Integer(I4B)         :: IUnTraj, IUnTopo
        Integer(I4B)         :: NFramesRead =  0
        Logical              :: LFirstEntry      = .TRUE.
        SAVE                 :: LFirstEntry, NFramesRead, IUnTraj, IUnTopo

        if (LFirstEntry) then
            LFirstEntry = .FALSE.
            IUnTraj     =  NGetFileUnit()
            IUnTopo     =  NGetFileUnit()

            open(IUnTraj, File=TRIM(CFileTraj), status='OLD',access='SEQUENTIAL', ERR=67)
            open(IUnTopo, File=TRIM(CFileTopo), status='OLD',access='SEQUENTIAL', ERR=77)
                
            ! do stuff here
            call Read_Topology(ITopo, IUnTopo, Topo, AtomType)

            if (.NOT. Allocated(Atom    )) Allocate(    Atom(Topo%NAt     ))
            if (.NOT. Allocated(AtomType)) Allocate(AtomType(Topo%NAtTypes))
        endif



        Select Case (IRead)
            Case(1) ; Call Read_LAMMPS(     IRead, IUnTraj, Atom, AtomType, Cell, Topo)
            Case(2) ; Call Read_LAMMPS_dip( IRead, IUnTraj, Atom, AtomType, Cell, Topo)

            Case Default ; print*, "Read_Sys: ERROR - This trajectory format is not supported. IRead =", IRead 
                         ; close(IUnTraj) ; close(IUnTopo) ; STOP
        End Select

        NFramesRead      = NFramesRead + 1

        ! for both of these cases the program will not be reading a frame again.
        if (ILastFrame == NFramesRead) then
            close(IUnTraj)
            close(IUnTopo)
        endif

        return
        
        67 print*, 'Read_Sys: Trajectory file open problem'   ; close(IUnTraj) ; STOP
        77 print*, 'Read_Sys: Topology file open problem'     ; close(IUnTopo) ; STOP
    End Subroutine Read_Sys

    !------------------------------------------------------------------------------
    !  Read a number of comment lines from a given file
    !
    !------------------------------------------------------------------------------
    Subroutine ReadCommentLines(IFileUnit, NLines)
        Implicit None
        Integer(I4B) :: IFileUnit, NLines, i ; Character :: CDummy*1
        Intent(IN)   :: IFileUnit, NLines
        do i=1, NLines  ;  read(IFileUnit,*) CDummy  ;  enddo
    End Subroutine ReadCommentLines

    !------------------------------------------------------------------------------
    ! Give sequential file unit numbers on each call, starting from "NFileUnit"
    !
    !------------------------------------------------------------------------------
    Function NGetFileUnit()
        Implicit None
        Integer(I4B)       :: NGetFileUnit
        Integer(I4B), Save :: NFileUnit = 20
        NGetFileUnit = NFileUnit
        NFileUnit    = NFileUnit + 1
    End Function NGetFileUnit

End Module Input
