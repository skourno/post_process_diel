!----------------------------------------------------------------------------------
!                              compute.f90  -  description
!                             ------------------------------
!                  begin    : Thu 27.02.2020
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
!
! This modules containes a series of computations that were initially used for
! the dielectric constant analysis of Stockmayer fluids. 
!
!----------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------
Module compute
	Use nrtype
 	Implicit None
 
contains


  !------------------------------------------------------------------------------
  ! Compute the dielectric constant of the configuration
  !
  !------------------------------------------------------------------------------
	subroutine compute_dielectric(Atom, AtomType, Cell, Topo, sum_dmom_sq, sum_dmom)
		Use system,  Only :  Bead, BeadType, Domain, Topology

 		Type(Bead),     Allocatable   :: Atom(:)
 		Type(BeadType), Allocatable   :: AtomType(:)
 		Type(Domain)                  :: Cell
 		Type(Topology)                :: Topo

 		Integer(I4B)                  :: iAt, jAt
 		Real(DP)                      :: sum_dmom(3), sum_dmom_sq
 		Intent(IN)                    :: Atom, AtomType, Cell, Topo
 		Intent(OUT)                   :: sum_dmom_sq, sum_dmom

 		sum_dmom(:) = ZERO

 		do iAt = 1, Topo%NAt
 			sum_dmom  = sum_dmom  + Atom(iAt)%dMom(:)
 		enddo

 		sum_dmom_sq = DOT_PRODUCT(sum_dmom,sum_dmom)

 		!print*, sum_dmom_sq, sum_dmom
	end subroutine compute_dielectric

  !------------------------------------------------------------------------------
  ! Compute relative P_1 average
  !
  !------------------------------------------------------------------------------
	subroutine compute_ave_relative_P1(Atom, AtomType, Cell, Topo)
		Use system,    Only :  Bead, BeadType, Domain, Topology
		Use Utilities, Only :  Angle

 		Type(Bead),     Allocatable   :: Atom(:)
 		Type(BeadType), Allocatable   :: AtomType(:)
 		Type(Domain)                  :: Cell
 		Type(Topology)                :: Topo

 		Integer(I4B)                  :: iAt, jAt, npairs
 		Real(DP)                      :: angle_sum, angle_ave


 		do iAt = 1, Topo%NAt
 			do jAt = iAt+1, Topo%NAt
 				angle_sum = angle_sum + COS(Angle(Atom(iAt)%dMom,Atom(jAt)%dMom)) 
 			enddo
 		enddo

 		npairs     = (Topo%NAt * (Topo%NAt - 1))/2
 		angle_ave  = angle_sum / npairs

 		!print*, angle_ave, angle_sum, npairs
	end subroutine compute_ave_relative_P1

End Module compute