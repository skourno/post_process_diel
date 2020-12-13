!----------------------------------------------------------------------------------
!                             read_cases.f90  -  description
!                             ------------------------------
!                  begin    : Tue 22.11.2016
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
!
! This module contains various routines to read different file formats. The loaded
! trajectories are saved into the "system" module using the appropriate data types.
!
!----------------------------------------------------------------------------------
!
! xx.xx.2016           
!
!----------------------------------------------------------------------------------
 Module Read_Cases
 	Use nrtype
 	Implicit None
 
 contains


    !------------------------------------------------------------------------------
    !  Reads LAMMPS dump file   
    !
    !------------------------------------------------------------------------------
    Subroutine Read_LAMMPS (IRead, IUnTraj, Atom, AtomType, Cell, Topo)
        Use system,      Only : Bead, BeadType, Topology, Domain
        Use Utilities,   Only : Fold_Pos_In_Cell
        Implicit None
        Type(Bead),     Allocatable   :: Atom(:)
        Type(BeadType), Allocatable   :: AtomType(:)
        Type(Domain)                  :: Cell
        Type(Topology)                :: Topo
        Integer(I4B)                  :: IRead, IUnTraj
        Integer(I4B)                  :: iAt, NAtomsDump, id, IType, i
        Real(DP)                      :: coords(3)

        call ReadCommentLines(IUnTraj, 3)
        
        read(IUnTraj,*, ERR=66) NAtomsDump

        if (NAtomsDump .NE. Topo%NAt) STOP "Read_LAMMPS: ERROR - Number of atoms mismatch"

        call ReadCommentLines(IUnTraj, 1)

        do i = 1,3 ; read(IUnTraj,*, ERR=66) Cell%R0(i), Cell%Rl(i) ; enddo

        Cell%Edge(:) = Cell%Rl(:) - Cell%R0(:)
        
        call ReadCommentLines(IUnTraj, 1)


        do iAt = 1, Topo%NAt
            read(IUnTraj,*, ERR=66) id, IType, coords(:)

            if (IType > Topo%NAtTypes) STOP "Read_LAMMPS: ERROR - invalid atom type detected in the dump file"

            Atom(id)%IType = IType
            Atom(id)%R(:)  = coords(:) / Cell%Edge(:)

            ! fold for safety. sometimes LAMMPS outputs coordinates outside the box
            Call Fold_Pos_In_Cell(Atom(id)%R, Cell%Edge)
        enddo

        return
        66 print*, 'Read_LAMMPS: Read Error'; close(IUnTraj) ; STOP

    End Subroutine Read_LAMMPS

    !------------------------------------------------------------------------------
    !  Reads LAMMPS dump.dip file   
    !
    !------------------------------------------------------------------------------
    Subroutine Read_LAMMPS_dip (IRead, IUnTraj, Atom, AtomType, Cell, Topo)
        Use system,      Only : Bead, BeadType, Topology, Domain
        Use Utilities,   Only : Fold_Pos_In_Cell
        Implicit None
        Type(Bead),     Allocatable   :: Atom(:)
        Type(BeadType), Allocatable   :: AtomType(:)
        Type(Domain)                  :: Cell
        Type(Topology)                :: Topo
        Integer(I4B)                  :: IRead, IUnTraj
        Integer(I4B)                  :: iAt, NAtomsDump, id, IType, i
        Real(DP)                      :: coords(3), dMom(3)

        call ReadCommentLines(IUnTraj, 3)
        
        read(IUnTraj,*, ERR=66) NAtomsDump

        if (NAtomsDump .NE. Topo%NAt) STOP "Read_LAMMPS: ERROR - Number of atoms mismatch"

        call ReadCommentLines(IUnTraj, 1)

        do i = 1,3 ; read(IUnTraj,*, ERR=66) Cell%R0(i), Cell%Rl(i) ; enddo

        Cell%Edge(:) = Cell%Rl(:) - Cell%R0(:)
        
        call ReadCommentLines(IUnTraj, 1)


        do iAt = 1, Topo%NAt
            read(IUnTraj,*, ERR=66) id, IType, coords(:), dMom(:)

            if (IType > Topo%NAtTypes) STOP "Read_LAMMPS: ERROR - invalid atom type detected in the dump file"

            Atom(id)%IType   = IType
            Atom(id)%R(:)    = coords(:) - Cell%R0(:) ! the box origin is set to (0,0,0)
            Atom(id)%dMom(:) = dMom(:)

            !print*, id, Atom(id)%IType, Atom(id)%R(:), Atom(id)%dMom(:)   

            ! fold for safety. sometimes LAMMPS outputs coordinates outside the box
            Call Fold_Pos_In_Cell(Atom(id)%R, Cell%Edge)
        enddo

        return
        66 print*, 'Read_LAMMPS_dip: Read Error'; close(IUnTraj) ; STOP

    End Subroutine Read_LAMMPS_dip

    !------------------------------------------------------------------------------
    !
    !  Read the topology file. FORMAT:
    !
    !     #   NAt         NBonds         NAng        NDih
    !         ###          ###           ###          ###
    !        
    !     # NAtTypes     NBondTypes    NAngTypes   NDihTypes
    !          ###          ###          ###           ### 
    !           
    !     # epsilon
    !         eps.1.1          eps.1.2   ...    eps.1.NAtTypes
    !         eps.2.1             ...    ...        ...
    !          ...                ...    ...        ...
    !         eps.NAtTypes.1      ...    ...    eps.NAtTypes.NAtTypes
    !
    !     # sigma      mass      charge 
    !        sigma.1    mass.1    charge.1
    !          ...       ...        ...
    !        sigma.NAtTypes   ...      ...
    !
    !  NOTE1: The topo file should only include parameters and quantities that 
    !         remain the same for the whole simulation. If this is not TRUE then 
    !         expect a separate call of the routine accompanied with comments etc
    !  
    !  NOTE2: The topology file should be open. Only it's index is passed (IUnTopo)
    !
    !------------------------------------------------------------------------------
    Subroutine Read_Topology(ITopo, IUnTopo, Topo, AtomType)
        Use System,               Only : Topology, BeadType
        Implicit None
        Integer(I4B)                  :: ITopo, IUnTopo
        Type(Topology)                :: Topo
        Type(BeadType), Allocatable   :: AtomType(:)
        Integer(I4B)                  :: i

        Select Case (ITopo) 
            Case (1) ! simple topo file
                
                call ReadCommentLines(IUnTopo, 1)
                read(IUnTopo, *, ERR=66) Topo%NAt, Topo%NBonds, Topo%NAng, Topo%NDih

                !print*, Topo%NAt, Topo%NBonds, Topo%NAng, Topo%NDih 

                call ReadCommentLines(IUnTopo, 1)
                read(IUnTopo, *, ERR=66) Topo%NAtTypes, Topo%NBondTypes, Topo%NAngTypes, Topo%NDihTypes

                !print*, Topo%NAtTypes, Topo%NBondTypes, Topo%NAngTypes, Topo%NDihTypes

                Allocate( AtomType(Topo%NAtTypes) )

                call ReadCommentLines(IUnTopo, 1)
                do i = 1, Topo%NAtTypes
                    read(IUnTopo, *, ERR=66) AtomType(i)%eps, AtomType(i)%sigma, AtomType(i)%mass, AtomType(i)%charge, AtomType(i)%dipMom
                    !print*, AtomType(i)%eps, AtomType(i)%sigma, AtomType(i)%mass, AtomType(i)%charge 
                enddo

            Case Default ; print*, "Read_Topology: ERROR - Invalid ITopo value"
        End Select 
        
        write(*,*) 
        write(*,*) "--------------------------------------"
        write(*,*) "  > Updated TOPOLOGY:   "
        write(*,*) "    NAt          =", Topo%NAt
        write(*,*) "    NBonds       =", Topo%NBonds
        write(*,*) "    NAng         =", Topo%NAng
        write(*,*) "    NDih         =", Topo%NDih
        write(*,*) "    NAtTypes     =", Topo%NAtTypes
        write(*,*) "    NBondTypes   =", Topo%NBondTypes
        write(*,*) "    NAngTypes    =", Topo%NAngTypes
        write(*,*) "    NDihTypes    =", Topo%NDihTypes 
        write(*,*) "--------------------------------------"
        write(*,*)

        return
        66 print*, 'Read_Topology: Read Error'; close(IUnTopo) ; STOP

    End Subroutine Read_Topology

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

 End Module Read_Cases