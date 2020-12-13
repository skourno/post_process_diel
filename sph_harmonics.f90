!----------------------------------------------------------------------------------
!                               sph_harmonics.f90  -  description
!                             ------------------------------
!                  begin    : Tue 19.05.2020
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
! This module contains all the necessary routines and functions to perform spherical
! harmonics calculations.
!----------------------------------------------------------------------------------
!
!  19.05.2020   1.0   Skou    Legendre polynomials        
!
!----------------------------------------------------------------------------------
Module Sph_harmonics
  Use nrtype

Contains

!-------------------------------------------------------
! Returns the corresponding Legendre polynomial of order 
! l with argument x
!-------------------------------------------------------
  Function LPoly(x, l)
    Implicit None
    Real(DP)      :: LPoly
    Real(DP)      :: x
    Integer(I4B)  :: l
    Intent(IN)    :: x, l
    
    Select Case(l)
    
    Case(0) ; LPoly = 1
    Case(1) ; LPoly = x
    Case(2) ; LPoly = HALF*(THREE*x*x   - ONE)
    Case(3) ; LPoly = HALF*(FIVE *x*x*x - THREE*x)
    
    Case(: -1) ; STOP "ERROR - LPoly : Requested negative Legendre Polynomial order"
    Case(4  :) ; STOP "ERROR - LPoly : Requested unsupported Legendre Polynomial order"
    End Select

  End Function LPoly

!-------------------------------------------------------
! Returns the corresponding Spherical Harmonic of order l,m
! provided th and ph as arguments
!-------------------------------------------------------
  Function SphHarm(th,ph,l,m)
    Use Utilities
    Implicit None
    COMPLEX(DPC) :: SphHarm 
    Real(DP)     :: th, ph
    Integer(I4B) :: l, m
    
    Intent(IN)   :: th, ph, l, m

    if     (l == 0) then
      
      if (m .ne. 0)            STOP 'ERROR - SphHarm : invalid m value'
      SphHarm = CMPLX( HALF * (1 / SQRT(PI)) , ZERO )    
    elseif (l == 1) then

      if ((m < -1).OR.(m > 1)) STOP 'ERROR - SphHarm : invalid m value'
      if (m == -1) SphHarm = + HALF * SQRT(THREE / (TWO*PI)) * SIN(th) * exp(  CMPLX(ZERO,-ph)  )
      if (m ==  0) SphHarm = +CMPLX(HALF * SQRT(THREE /      PI ) * COS(th) , ZERO )
      if (m ==  1) SphHarm = - HALF * SQRT(THREE / (TWO*PI)) * SIN(th) * exp(  CMPLX(ZERO,+ph)  ) 

    endif


  End Function SphHarm
!-------------------------------------------------------
! Returns the corresponding Spherical Harmonic of order l,m
! provided the xyz vector
!-------------------------------------------------------
  Function SphHarm_xyz(xyz,l,m)
    Use Utilities
    Implicit None
    COMPLEX(DPC) :: SphHarm_xyz 
    Real(DP)     :: xyz(3), sphCoords(3)
    Integer(I4B) :: l, m
    
    Intent(IN)   :: xyz, l, m

    sphCoords   = Sph_Coords(xyz)    

    SphHarm_xyz = sphHarm(sphCoords(2), sphCoords(3), l, m)

  End Function SphHarm_xyz


End Module Sph_harmonics
