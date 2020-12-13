!----------------------------------------------------------------------------------
!                               utilities.f90  -  description
!                             ------------------------------
!                  begin    : Tue 12.01.2017
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
! This module contains various utility routines for analyzing molecular trajectories
! It is designed to work within my platform, so it requires the module "system" to 
! compile. However, most of the routines can be used independently. To do that you 
! should comment out all the routines that utilize data from the system module.
!----------------------------------------------------------------------------------
!
!  16.01.17   1.0   Skou    Folding, Minimum Image, Vectors, Rotation        
!
!----------------------------------------------------------------------------------
Module Utilities
  Use nrtype

Contains

!------------------------------------------------------------------------------------------
! Fold 'Diffussive' or 'Unfolded' coordinates into the cell with dimensions 'Cell_Dim'.
! Then check the folded coordinate (commented).
!------------------------------------------------------------------------------------------
  Subroutine Fold_Pos_In_Cell(Rdp, Cell_Dim)
    Implicit None
    Real(DP)      :: Rdp(3), Cell_Dim(3)
    Intent(IN)    :: Cell_Dim
    Intent(INOUT) :: Rdp

    !Call Displace_Pos_on_Cell_Edge(Rdp)

    Where(Rdp >= Cell_Dim) Rdp = Rdp - Cell_Dim *   INT(Rdp/Cell_Dim)
    Where(Rdp <  ZERO)     Rdp = Rdp + Cell_Dim * (-INT(Rdp/Cell_Dim) + 1)
  End Subroutine Fold_Pos_In_Cell

!--------------------------------------------------------------------------------------------------------
! Fold 'R' coordinates into the cell.
!
! NOTES : 'CellDimensions' of 'Input' module is used. Watch for misuse.
!--------------------------------------------------------------------------------------------------------
  Subroutine FoldPos(R)
    Use System, Only : Cell
    Implicit None
    Real(DP)      :: R(3)
    Intent(INOUT) :: R

    !Call Displace_Pos_on_Cell_Edge(Rdp)

    Where(R > Cell%Edge)  R = R - Cell%Edge *   INT(R/Cell%Edge)
    Where(R < ZERO)       R = R + Cell%Edge * (-INT(R/Cell%Edge) + 1)
  End Subroutine FoldPos

!----------------------------------------------------------------------------------
! Return the Image of coordinate R.
! We assume that the simulation Cell extends between (0.,0.,0.) and 'CellDim(:)'
!----------------------------------------------------------------------------------
  Subroutine Get_Image_FromPos(R, CellDim, Img)
    Implicit None
    Integer(I4B)            :: Img(3)
    Real(DP)                :: R(3), CellDim(3)
    Intent(IN)              :: R, CellDim
    Intent(OUT)             :: Img

    Img(:) = 0
    Where(R > CellDim)  Img =  INT(R/CellDim)
    Where(R < ZERO)     Img =  INT(R/CellDim) - 1
    !print*,IMG
  End Subroutine Get_Image_FromPos

!-------------------------------------------------------
! Cross Product "V = Va x Vb" of two vectors "Va, Vb"
!-------------------------------------------------------
  Function CrossPr(Va, Vb)
    Implicit None
    Real(DP)      :: Va(3), Vb(3), CrossPr(3)
    Intent(IN)    :: Va,    Vb
    CrossPr(1)  = Va(2)*Vb(3) - Va(3)*Vb(2)
    CrossPr(2)  = Va(3)*Vb(1) - Va(1)*Vb(3)
    CrossPr(3)  = Va(1)*Vb(2) - Va(2)*Vb(1)
  End Function CrossPr

!----------------------------------------------------------------------------------
! Rotate a given vector a around an axis defined by vector n by an angle theta and
! produce the rotated unit vector b
!----------------------------------------------------------------------------------
    Subroutine Rotate(r, n, theta, b)
       Implicit None
       Real(DP)    :: r(3), n(3), b(3), r_unit(3), n_unit(3), r_ext_n(3)
       Real(DP)    :: theta, n_norm, r_norm, n_dot_r
       Intent(IN)  :: r, n, theta
       Intent(OUT) :: b

       n_norm    = Sqrt(Sum(n**2))
       n_unit(:) = n(:)/n_norm

       r_norm     = Sqrt(Sum(r**2))
       r_unit(:)  = r(:)/r_norm

       n_dot_r    = Dot_Product(n_unit,r_unit)

       r_ext_n(:) = CrossPr(r_unit, n_unit)

       b(:) = dcos(theta)*r_unit(:) + n_dot_r*(ONE-dcos(theta))*n_unit(:) + dsin(theta)*r_ext_n(:)
       b(:) = b(:) / Sqrt(Sum(b**2))
     End Subroutine Rotate

!----------------------------------------------------------------------------------
! Calculate the distance between two 3D Points by applying the Min. Im. Criterion.
! When the Min.Im.Cr. uses images to get the correct distance then DR points from
! A to B as well. The Cell dimensions must be provided
!----------------------------------------------------------------------------------
  Function DistMinIm(PositionB, PositionA, CellDim)
    Implicit None
    Real(DP)   :: PositionA(3), PositionB(3), CellDim(3), DR(3), DistMinIm
    Intent(IN) :: PositionA, PositionB, CellDim

    DR        = PositionB - PositionA
    DR        = DR - DNINT(DR / CellDim) * CellDim  ! Min. Im. Criterion
    DistMinIm = SQRT(SUM(DR**2))

  End Function DistMinIm

!----------------------------------------------------------------------------------
! Calculate the distance between two 3D Points without applying the Min. Im. Criterion.
!----------------------------------------------------------------------------------
  Function Distance(PositionB, PositionA)
    Implicit None
    Real(DP)   :: PositionA(3), PositionB(3), DR(3), Distance
    Intent(IN) :: PositionA, PositionB

    DR = PositionB - PositionA
    Distance = SQRT(SUM(DR**2))
  End Function Distance

!-------------------------------------------------------
! Returns magnitude of vector V
!-------------------------------------------------------
  Function Magnitude(V)
    Implicit None
    Integer    :: i
    Real*8     :: V(3), Magnitude
    Intent(IN) :: V
    Magnitude = 0.D0

    do i = 1,3
        Magnitude = V(i)**2 + Magnitude
    enddo

    Magnitude = SQRT(Magnitude)
  End Function Magnitude

!-------------------------------------------------------
! Returns the angle between vectors v & u
!-------------------------------------------------------
  Function Angle(v,u)
    Implicit None
    Real(DP)   :: Angle
    Real(DP)   :: v(3), u(3), dot, costheta
    Intent(IN) :: v, u

    dot       =  DOT_PRODUCT(v,u)
    costheta  =  dot / ( Magnitude(u)*Magnitude(v) )
    angle     =  DACOS(costheta)

  End Function Angle

!-------------------------------------------------------
! Returns the factorial of the argument n
!-------------------------------------------------------
  Function Factorial(n)
    Implicit None
    Integer(I4B) :: Factorial 
    Integer(I4B) :: n, i
    
    Intent(IN)   :: n
    
    Factorial = 1
    do i = 1,n
      Factorial = Factorial * i
    enddo
  End Function Factorial

!-------------------------------------------------------
! Returns the spherical coordinates from xyz coordinates
! The result is the vector (r,th,ph)
!-------------------------------------------------------
  Function Sph_Coords(xyz) 
    Implicit None
    Real(DP)     :: Sph_Coords(3), xyz(3)
    Real(DP)     :: x, y, z, r, th, ph 
    
    Intent(IN)   :: xyz
    
    x = xyz(1)
    y = xyz(2)
    z = xyz(3)


    r  = SQRT(x*x + y*y + z*z)
    th = ACOS(z/r)
    ph = ATAN(y/x)

    if (ph < 0) ph = 2*PI + ph

    if ((th < ZERO).OR.(th >   PI)) STOP 'Sph_Coords : Out of bounds theta'
    if ((ph < ZERO).OR.(ph > 2*PI)) then
      print*, 'Sph_Coords : Out of bounds phi', xyz, r,th,ph
      STOP
    endif

    
    Sph_Coords = (/ r , th , ph /)

  End Function Sph_Coords

!-------------------------------------------------------
! Returns unitary vectors with the z axis pointing
! along the vector v
!-------------------------------------------------------
  Function Rotate_z_axis_to_align_with(v) Result(unitVecs)
    Implicit None
    Real(DP) :: unitVecs(3,3), v(3)
    Real(DP) :: x(3), y(3), z(3)


    z(:) = v(:) / Magnitude(v) ! z axis unitary

    if (z(1) < 1.D-10) z(1) = z(1) + 1.D-10 

    y(1) = (-z(2)-z(3)) / z(1) ! choose a unitary y at random. there are infinite choices
    y(2) = ONE
    y(3) = ONE
    y(:) = y(:) / Magnitude(y)

    if (y(1) < 1.D-10) y(1) = y(1) + 1.D-10

    x(:) = CrossPr(y,z)

    unitVecs(:,1) = x
    unitVecs(:,2) = y 
    unitVecs(:,3) = z

  End Function Rotate_z_axis_to_align_with

!-----------------------------------------------------------------------------------------------------------------------------
! Given the points 'A', 'B' inside the cell, and the FOLDED coordinates 'RA, RB', 'DR' is the vector connecting 'RA'
! (inside the cell) with the nearest image of 'B', and it points from 'RA' to 'RB'. The nearest image of 'B' can be situated
! inside (folded coordinates) or outside the box. The coordinates of this image are 'RA + DR'.
! When the nearest image of 'B' to 'A' is out of the box 'LNearImgOutOfBox = .TRUE.', otherwise 'FALSE'
!
! ATTENTION : 'RA, RB' should be given as FOLDED coordinates
!-----------------------------------------------------------------------------------------------------------------------------
  Subroutine Get_Min_Img_Conn_Vec(RA, RB, DR, CellDim, LNearImgOutOfBox)
    Implicit None
    Real(DP)     :: RB(3),      & !Folded coordinates
                    RA(3),      & !   //       //
                    DR(3),      & !RB - RA
                    DR_Img(3),  & !Has a > ZERO component when the Min.Im.Crit. connects 'RA' with an image of 'B' out of the simulation box
                    CellDim(3) 
    Integer(I4B) :: Img(3)
    Logical      :: LNearImgOutOfBox
    OPTIONAL     :: LNearImgOutOfBox
    Intent(IN)   :: RA, RB, CellDim
    Intent(OUT)  :: DR, LNearImgOutOfBox

    DR               = RB - RA

    !This vector has (real) NonZero Components (+-1.D0, +-2.D0, ...) where '|DR(i)| > CellDimensions(i)/2'  (Min.Im.Crit.),
    !and zero components (0.D0) where '|DR(i)| < CellDimensions(i)/2'
    DR_Img = DNINT(DR / CellDim)

    !When the nearest image of 'B' to 'A' is out of the box 'LNearImgOutOfBox = .TRUE.' 
    if (PRESENT(LNearImgOutOfBox)) then
      Img(:) = NINT(DR_Img)
      if(ANY(Img /= 0)) LNearImgOutOfBox = .TRUE.
    endif

    !'DR' points from 'RA' to 'RB' and connects 'RA' (by definition inside the cell) with the specific image of
    !'B' that has the minimum distance to 'RA'.
    DR = DR - DR_Img * CellDim
  End Subroutine Get_Min_Img_Conn_Vec



End Module Utilities
