!----------------------------------------------------------------------------------
!                               system.f90  -  description
!                             ------------------------------
!                  begin    : Mon 21.02.2018
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
!
! This is a module that is used for saving coordinates and other parameters that 
! define the loaded system. From here you can access the atoms and loop at them.
!
! All the units are dimesionless, following the same reduction rules as . "Real" units are described 
! in the LAMMPS manual.
!
!----------------------------------------------------------------------------------
!
!  16.01.17   1.0   Skou    Bead, BeadType, Topology, Domain      
!  21.02.18   2.0   Skou    add_atoms_to_system
!
!----------------------------------------------------------------------------------
 Module System
 	Use nrtype
 	Implicit None

 	Type Bead
 		Integer(I4B)              :: IType        = 0     ! Type of the bead
 		Integer(I4B)              :: NBonds       = 0     ! number of bonded beads
 		Integer(I4B), Allocatable :: JBead(:)             ! indices of the bonded beads, 1-NBonds

 		Real(DP)                  :: R(3)         = ZERO  ! xyz coordinates
        	Real(DP)                  :: dMom(3)      = ZERO  ! dipole moment vector
 	End Type Bead

 	Type BeadType
 		Real(DP)                  ::  sigma      = ZERO   
 		Real(DP)                  ::  charge     = ZERO
	        Real(DP)                  ::  dipMom     = ZERO
 		Real(DP)                  ::  eps        = ZERO    
 		Real(DP)                  ::  mass       = ZERO
 	EndType BeadType

 	Type Topology
 		Integer(I4B)              :: NAt         = 0      ! number of atoms
 		Integer(I4B)              :: NBonds      = 0      ! number of bonds
 		Integer(I4B)              :: NAng        = 0      ! number of angles
 		Integer(I4B)              :: NDih        = 0      ! number of dihedrals
 		Integer(I4B)              :: NMol        = 0      ! number of molecules

 		Integer(I4B)              :: NAtTypes    = 0
 		Integer(I4B)              :: NBondTypes  = 0
 		Integer(I4B)              :: NAngTypes   = 0
 		Integer(I4B)              :: NDihTypes   = 0
  	End Type Topology

  	Type Domain
  		Real(DP)                  :: R0(3)       = ZERO   ! box origin
  		Real(DP)                  :: Rl(3)       = ZERO   ! box end  (Rl(i) - R0(i) = edge vector in i dimension)
  		Real(DP)                  :: Edge(3)     = ZERO   ! domain dimensions
  	End Type Domain


 	Type(Bead),     Allocatable   :: Atom(:)
 	Type(BeadType), Allocatable   :: AtomType(:)
 	Type(Domain)                  :: Cell
 	Type(Topology)                :: Topo
 	Real(DP)                      :: Dielectric  = ONE    ! initialize to vacuum

 contains

    !------------------------------------------------------------------------------
    !  Adds atoms to the system. 
    !  The coordinates of the new particle are always initialized to ZERO. The goal
    !  of this routine is just to expand the capacity of the arrays and initialize
    !  everything appropriately. Meaningful coordinates are to be generated elsewhere
    !  according to the desired application.
    !
    !  IType     : is the type of the atoms that will be inserted
    !  NNewAtoms : is the number of atoms you desire to be added to the system 
    !
    !  NOTE(!!!) : The topology is not changed! 
    !------------------------------------------------------------------------------
    Subroutine add_atoms_to_system(NNewAtoms, IType)
    	Implicit None
    	Integer(I4B)             :: IType
    	Integer(I4B)             :: NNewAtoms
    	Integer(I4B)             :: NOriginalAtoms, iAt
    	Type(Bead), Allocatable  :: AtomBackup(:)
    	Intent(IN)               :: IType, NNewAtoms

    	if ((IType     > Topo%NAtTypes).OR.(IType     <= 0)) STOP "add_atoms_to_system: ERROR - invalid IType"
    	if ((NNewAtoms > Topo%NAt     ).OR.(NNewAtoms <= 0)) STOP "add_atoms_to_system: ERROR - invalid NNewAtoms"

    	NOriginalAtoms = Topo%NAt

    	! do a backup of the old Atom array
    	Allocate(AtomBackup(NOriginalAtoms))
    	AtomBackup(:)  = Atom(:)

    	! purge the old Atom array and create a new one with the appropriate dimensions
    	Deallocate(Atom)
    	Allocate(Atom(NOriginalAtoms+NNewAtoms))

    	! transfer back the original data
    	Atom(1:NOriginalAtoms)  =  AtomBackup(:)

    	! fill in the added elements
    	do iAt = NOriginalAtoms+1, NOriginalAtoms+NNewAtoms
    		Atom(iAt)%IType  =  IType
    		Atom(iAt)%NBonds =  0      ! just particles for this routine
    		Atom(iAt)%R(:)   =  ZERO   ! initialize to zero. see description 
    	enddo

    End Subroutine add_atoms_to_system

 End Module System 
