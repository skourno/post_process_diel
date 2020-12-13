MODULE nrtype
!	Symbolic names for kind types of 4-, 2-, and 1-byte integers:
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
!	Symbolic names for kind types of single- and double-precision reals:
	INTEGER, PARAMETER :: SP = KIND(1.0)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
!	Symbolic names for kind types of single- and double-precision complex:
	INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
!	Symbolic name for kind type of default logical:
	INTEGER, PARAMETER :: LGT = KIND(.true.)
!	Frequently used mathematical constants (with precision to spare):
	REAL(DP), PARAMETER :: PI     =  3.141592653589793238462643383279502884197_dp
	REAL(DP), PARAMETER :: TWOPI  =  6.283185307179586476925286766559005768394_dp
	REAL(DP), PARAMETER :: ZERO   =  0.0_dp
	REAL(DP), PARAMETER :: ONE    =  1.0_dp
    REAL(DP), PARAMETER :: TWO    =  2.0_dp
    REAL(DP), PARAMETER :: THREE  =  3.0_dp
    REAL(DP), PARAMETER :: FOUR   =  4.0_dp
    REAL(DP), PARAMETER :: FIVE   =  5.0_dp
    REAL(DP), PARAMETER :: SIX    =  6.0_dp
    REAL(DP), PARAMETER :: EIGHT  =  8.0_dp
    REAL(DP), PARAMETER :: TEN    = 10.0_dp
    REAL(DP), PARAMETER :: DOZEN  = 12.0_dp
    REAL(DP), PARAMETER :: TWFOUR = 24.0_dp

    REAL(DP), PARAMETER :: HALF           = 0.5_dp
    REAL(DP), PARAMETER :: ONEHALF        = 1.5_dp
    REAL(DP), PARAMETER :: QUARTER        = 0.25_dp
    REAL(DP), PARAMETER :: RAD_TO_DEGREES = 180.0_dp / PI
    REAL(DP), PARAMETER :: DEGREES_TO_RAD = PI / 180.0_dp
END MODULE nrtype
