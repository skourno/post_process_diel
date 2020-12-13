!----------------------------------------------------------------------------------
!                             histograms.f90  -  description
!                             ------------------------------
!                  begin    : Tue 15.01.2017
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
! There are two type of histograms that can be allocated and processed by these routines.
! 1) Standard histograms       (only one value is needed by PutInHistogram)
! 2) 'MeanValue' Histograms    (two values are needed by PutInHistogram)
! 
! The latter consist of boxes where a given X value is stored additively to a box with
! index 'IBox' according to a given corresponding 'ValueForBoxIndex' value. 'HistVis(:)' 
! is allocated and used only for this type of histograms. 'HistVis(:)'stores the # of hits 
! (visits) to each box. 'HistArr(:)' stores the sum of the variables added at each box 
! during runtime. This way a Mean Value Histogram (plot) of the average of Y for different 
! values of X can be acquired at the end of the calculations.
!
! ATTENTION : Save Histogram variables in routines which are repeatedly called
!             'HistVis' component is allocated and used only for 'MeanValue' Histograms
!
! NOTES     : 'WriteHistogram' is called for both types. Put 'IMean=0' for the 2nd type.
!             'PutInHistogram(HistName, x, CType, ValueForBoxIndex)' should be called
!              without the last two optional arguments in case of type 1 Histograms.            
!
!             For Type 1 Histograms call all subroutines without their optional 
!             'CType' and 'ValueForBoxIndex' arguments.
!             
!             The optional variable 'CAuxVarLabel, AuxVar' of ''WriteHistogram' 
!             can be used if the printing of additional variables is required, either for 
!             Type 1 (Standard) or Type 2 (MeanValue) histograms.
!
!             'NormalizeHistogram' does different things on the hist arrays according to Type.
! 
!             For Type 2 Histograms you should use 'CType = "MeanVal"' when calling
!             the folowing subroutines,
!                 NewHistogram(NBox, BoundLeft, BoundRight, CType)
!                 PutInHistogram(HistName, x, CType, ValueForBoxIndex)
!                 NormalizeHistogram(HistName, CType)
!                 WriteHistogram(HistName, CFileName, IMean, CAuxVarLabel, AuxVar, CType)
!             For Type 2 Histograms, do not use Ctype when calling
!                 DelHistogram(HistName)
!----------------------------------------------------------------------------------
!
!  16.01.17   1.0   Skou    Standard and MeanValue Histograms   
!
!----------------------------------------------------------------------------------   
Module Histograms
  USE  nrtype
  
  TYPE Histogram
    REAL(DP),     ALLOCATABLE :: HistArr(:)   !HistArr(NBox), Histogram Array Values
    REAL(DP),     ALLOCATABLE :: HistVis(:)   !HistVis(NBox), Box Visits
    INTEGER(I4B)              :: NBox
    REAL(DP)                  :: Visits
    REAL(DP)                  :: OutOfRangeLeft
    REAL(DP)                  :: OutOfRangeRight
    REAL(DP)                  :: BoxLength
    REAL(DP)                  :: BoundLeft
    REAL(DP)                  :: BoundRight   
  END TYPE Histogram

  CONTAINS

!------------------------------------------------------------------------------------
! Constructor for a new Histogram object.
! It should be called without the last *OPTIONAL ARGUMENT* in case of *TYPE 1* 
! Histogram.'CType' is essentially used in order to discriminate between Type 1 and
! 2 histograms and take different action when called.
!------------------------------------------------------------------------------------
     FUNCTION NewHistogram(NBox, BoundLeft, BoundRight, CType) Result(HistName)
        IMPLICIT None
        Type(Histogram)          :: HistName
        Character,     Optional  :: CType*(*)
        INTEGER(I4B), INTENT(IN) :: NBox 
        REAL(DP),     INTENT(IN) :: BoundLeft, BoundRight
        
        !Check for wrong boundaries
        if(BoundLeft >= BoundRight) then
           print*,'Boundaries ', BoundLeft, BoundRight
           STOP 'NewHistogram: Wrong Boundaries'
        endif
        ALLOCATE(HistName%HistArr(NBOX))
        HistName%NBox       = NBox
        HistName%BoundLeft  = BoundLeft
        HistName%BoundRight = BoundRight
        HistName%BoxLength  = (BoundRight - BoundLeft) / DBLE(NBox)
        HistName%HistArr    = ZERO  ! Initialize histogram array to zero       
        HistName%Visits     = ZERO  ! Initialize # of calls to PutInHistogram, (# of values thrown in boxes) 
        HistName%OutOfRangeLeft  = ZERO          
        HistName%OutOfRangeRight = ZERO          

        !Used only for Mean value Histograms
        if(Present(CType)) then
           Allocate(HistName%HistVis(NBOX))
           HistName%HistVis = ZERO
        endif
     END FUNCTION NewHistogram

!-------------------------------------------------------------------------------------
! Put a given value to the  corresponding box of the histogram array.
! It should be called without the last two *OPTIONAL ARGUMENTS* in case of *TYPE 1* 
! Histograms. 'CType' is essentially used in order to discriminate between Type 1 and
! 2 histograms and take different action when called.
!-------------------------------------------------------------------------------------
     SUBROUTINE PutInHistogram(HistName, x, CType, ValueForBoxIndex)
        IMPLICIT None
        TYPE(Histogram), INTENT(INOUT) :: HistName
        CHARACTER,       OPTIONAL      :: CType*(*)
        REAL(DP),        OPTIONAL      :: ValueForBoxIndex
        REAL(DP)                       :: x, ValBoxInd, ValBoxAdd, &
                                          DistInBoxLengths, DistFromLeftBound
        INTEGER(I4B)                   :: IBox
        INTENT(IN)                     :: x, CType, ValueForBoxIndex

        if(Present(CType)) then          !Mean value Histogram Case
           ValBoxInd = ValueForBoxIndex
           ValBoxAdd = x
        else                             !Standard Histogram Case
           ValBoxInd = x
           ValBoxAdd = ONE
        endif

        if( ValBoxInd < HistName%BoundLeft ) then
          HistName%OutOfRangeLeft  = HistName%OutOfRangeLeft  + ONE
          return
        endif
        if( ValBoxInd > HistName%BoundRight) then 
          HistName%OutOfRangeRight = HistName%OutOfRangeRight + ONE
          return
        endif

        DistFromLeftBound = Abs(ValBoxInd - HistName%BoundLeft)
        DistInBoxLengths  = DistFromLeftBound/ HistName%BoxLength
        IBox              = INT(DistInBoxLengths) + 1 
        !Avoid 'array out of bounds' case when x = BoundRight, or BoundLeft (IBox indetermined).
        if(IBox > HistName%NBox)    IBox = HistName%NBox
        !print*, IBox, DistFromLeftBound, HistName%BoundRight , x
        !print*, IBox, HistName%NBox, x, DistFromLeftBound, DistInBoxLengths 

        HistName%HistArr(IBox) = HistName%HistArr(IBox) + ValBoxAdd
        HistName%Visits        = HistName%Visits        + ONE

        if(Present(CType)) HistName%HistVis(IBox) = HistName%HistVis(IBox) + ONE
     END SUBROUTINE PutInHistogram

!------------------------------------------------------------------------
! For Type 1 histograms performs a Normalization of the histogram array
! For Type 2 histograms calculates the mean values for each box
!
! Add checking for box with zero visits (Mask the array). Do not normalize there
! Add option for normalization to one
!------------------------------------------------------------------------
     SUBROUTINE NormalizeHistogram(HistName, CType)
        IMPLICIT None
        TYPE(Histogram), INTENT(INOUT) :: HistName
        CHARACTER,       OPTIONAL      :: CType*(*)
        INTENT(IN)                     :: CType

        if(Present(CType)) then
           Where(HistName%HistArr /= ZERO) HistName%HistArr = HistName%HistArr / HistName%HistVis        
        else
           Where(HistName%HistArr /= ZERO) HistName%HistArr = HistName%HistArr / HistName%Visits 
           if(ANY(HistName%HistArr /= ZERO)) &
              HistName%HistArr = HistName%HistArr / SUM(HistName%HistArr * HistName%BoxLength)
        endif              
     END SUBROUTINE NormalizeHistogram

!------------------------------------------------------------------------
! Deallocation of histogram array, and vanishing of histogram components
!------------------------------------------------------------------------
     Subroutine DelHistogram(HistName)
        Implicit None
        TYPE(Histogram), Intent(INOUT) :: HistName
        
        if(Allocated(HistName%HistArr)) then
           if(Allocated(HistName%HistVis)) Deallocate(HistName%HistVis)
           Deallocate(HistName%HistArr)
           HistName%NBox            = 0
           HistName%Visits          = ZERO     
           HistName%BoundLeft       = ZERO
           HistName%BoundRight      = ZERO
           HistName%BoxLength       = ZERO
           HistName%OutOfRangeLeft  = ZERO 
           HistName%OutOfRangeRight = ZERO          
        else
           STOP 'DelHistogram: Unallocated Histogram'
        endif
     End Subroutine DelHistogram

!-----------------------------------------------------------------------------------------------------
! Histogram array printing. 
!
! NOTES : An *OPTIONAL ARRAY* 'AuxVar' and an *OPTIONAL LABEL* 'CAuxVarLabel'
!         with the corresponding variable names can be given and written on output. 
!         This costruct makes possible the *external* calculation and *passing* of 
!         some other variables (e.g. other moments of the distribution) for printing. 
!        
!         IMean = (0/1)->(No/Yes) : Writes the Area, Avg, Deviation, SqAvg , RMS
!                                   of the distribution. Use 0 for rdf histograms.
!
! ATTENTION : In case of a 'Type 2' histogram, if we want to print out HistName%HistVis(:)
!             as a 3rd column, then we have to call WriteHistogram(......, CType='2').
!             CAuxVarLabel, AuxVar keywords may be ommited, or included when calling
!             the subroutine. However, 'CType= 'any string' has to be 
!             included in the call of WriteHistograms, when it is called without 
!             CAuxVarLabel, AuxVar arguments, because of positional 'matching'
!             of optional arguments. In 'Type_2' histograms no printouts for boxes without visits.
!-----------------------------------------------------------------------------------------------------
     Subroutine WriteHistogram(HistName, CFileName, IMean, CAuxVarLabel, AuxVar, CType)
        Use Input, ONLY : NGetFileUnit
        Implicit None
        Type(Histogram)            :: HistName
        Character                  :: CFileName*(*), CAuxVarLabel*(*), CLabel*27, CType*(*)
        Integer(I4B)               :: IMean, IUn, i
        Real(DP)                   :: AuxVar(:), BoxPositionOffset, AvgVal, Diff, &
                                      Area, BoxPos, BoxVal, AvgValSq, Deviation 
        Optional                   :: CAuxVarLabel, AuxVar, CType
        Intent(IN)                 :: HistName, CFileName, IMean, CAuxVarLabel, AuxVar, CType

        IUn = NGetFileUnit() 
        open(IUn, File=CFileName, Status='UNKNOWN')

        !Initial comment and aux variable writings
        write(IUn, '("# ")')

        if(Present(CAuxVarLabel)) write(IUn, '("# ", A, /, "# ")')  CAuxVarLabel

        if(Present(AuxVar)) then 
           do i=1, Size(AuxVar)  ;  write(IUn, '("# ", E17.8)')  AuxVar(i)  ;  enddo
        endif
        write(IUn, '("# ")')
        write(IUn, '("# ", I5, " Boxes, BoxLength= ", F8.3, " , Range ", E15.5, " - ", E15.5)') &
              HistName%NBox, HistName%BoxLength, HistName%BoundLeft, HistName%BoundRight
       
        CLabel = ' values were out of range '
        if( HistName%OutOfRangeLeft > ZERO) &
          write(IUn, 35) HistName%OutOfRangeLeft, HistName%Visits , CLabel, 'Left'
        if( HistName%OutOfRangeRight > ZERO) &
          write(IUn, 35) HistName%OutOfRangeRight, HistName%Visits , CLabel, 'Right'
        35 format('# ',/,'# ' E15.5, ' from ', E15.5, A, 2x, A)
        
        !Calculate and write mean values
        BoxPositionOffset = Half * HistName%BoxLength
        if(IMean == 1) then
          Area     = Sum(HistName%HistArr) * HistName%BoxLength
          if(Area == ZERO) then
             print*, 'WriteHistogram: ', Trim(CFileName), ' Area = ZERO. Void Histogram.'
             return
          endif
          AvgVal   = ZERO
          AvgValSq = ZERO
          do i=1, HistName%NBox
             BoxPos   = HistName%BoundLeft + (i*HistName%BoxLength)-BoxPositionOffset
             BoxVal   = BoxPos   * HistName%HistArr(i) * HistName%BoxLength   !Prob= Bin area
             AvgVal   = AvgVal   + BoxVal
             AvgValSq = AvgValSq + BoxVal * BoxPos
          enddo
          AvgVal    = AvgVal   / Area
          AvgValSq  = AvgValSq / Area      
          Diff      = AvgValSq - AvgVal**2
          if(Diff > ZERO) then
             Deviation = Sqrt(AvgValSq - AvgVal**2)
          elseif(Diff == ZERO) then
             Deviation = ZERO
          elseif(Diff < ZERO .AND. ABS(Diff) < 1.D-5) then   !Vanishing small but negative
             Deviation = ZERO
          elseif(Diff < ZERO) then
             STOP 'WriteHistogram: Deviation < 0'
          endif
             write(IUn, '("#",/, "# Area, Avg, Deviation, SqAvg, RMS",/, "# ", 5(2x,E15.8))'), &
                          Area, AvgVal, Deviation, AvgValSq, Sqrt(AvgValSq)
        endif

        write(IUn, "('# ')")

        do i=1, HistName%NBox
           BoxPos = HistName%BoundLeft + (i*HistName%BoxLength)-BoxPositionOffset
           if(Present(CType)) then
              !if(HistName%HistVis(i) > ZERO) write(IUn, 101) BoxPos, HistName%HistArr(i), HistName%HistVis(i)
              if(HistName%HistVis(i) > ZERO) write(IUn, 101) BoxPos, HistName%HistArr(i), HistName%HistVis(i)/HistName%Visits
 101       format(x, E15.5, 2x, E15.5, 2x, E15.5)
           else
             write(IUn, 100) BoxPos, HistName%HistArr(i)
           endif
 100       format(x, E15.5, 2x, E15.5)
        enddo
        close(IUn)
        
     End Subroutine WriteHistogram

!--------------------------------------------------------------------------------------
! Given a normalized histogram, return the Avg, SqAvg, StDev
!--------------------------------------------------------------------------------------
     Subroutine Avg_SqAvg_StDv_Hist(HistName, AvgVal, Deviation, AvgValSq)
        Use Input, ONLY : NGetFileUnit
        Implicit None
        Type(Histogram) :: HistName
        Integer(I4B)    :: i
        Real(DP)        :: BoxPositionOffset, AvgVal, Diff, &
                           Area, BoxPos, BoxVal, AvgValSq, Deviation 
        Intent(IN)      :: HistName
        Intent(OUT)     :: AvgVal, Deviation, AvgValSq 

        !Calculate and write mean values
        BoxPositionOffset = Half * HistName%BoxLength
        Area     = Sum(HistName%HistArr) * HistName%BoxLength
        if(Abs(Area - ONE) >= 0.02_dp) then
          print*, 'Area ', Area
          STOP 'Avg_SqAvg_StDv_Hist: Area <= 0.98'  !check
        endif

        AvgVal   = ZERO
        AvgValSq = ZERO
        do i=1, HistName%NBox
           BoxPos   = HistName%BoundLeft + (i*HistName%BoxLength)-BoxPositionOffset
           BoxVal   = BoxPos   * HistName%HistArr(i) * HistName%BoxLength   !Prob= Bin area
           AvgVal   = AvgVal   + BoxVal
           AvgValSq = AvgValSq + BoxVal * BoxPos
        enddo
        AvgVal    = AvgVal   / Area
        AvgValSq  = AvgValSq / Area      
        Diff      = AvgValSq - AvgVal**2
        if(Diff > ZERO) then
           Deviation = Sqrt(AvgValSq - AvgVal**2)
        elseif(Diff == ZERO) then
           Deviation = ZERO
        elseif(Diff < ZERO) then
           STOP 'Avg_SqAvg_StDv_Hist: Deviation < 0'
        endif
        return      
     End Subroutine Avg_SqAvg_StDv_Hist

!-------------------------------------------------------------------------------------
! Return the box 'IBox' where the input variable 'x' is contained. 
!-------------------------------------------------------------------------------------
     FUNCTION GetHistogramIBox(HistName, x) RESULT(IBox)
        IMPLICIT None
        TYPE(Histogram), INTENT(IN) :: HistName
        REAL(DP)                    :: x, ValBoxInd, DistInBoxLengths, DistFromLeftBound
        INTEGER(I4B)                :: IBox
        INTENT(IN)                  :: x

        ValBoxInd = x

        if(ValBoxInd < HistName%BoundLeft)  then ;  print*, ValBoxInd, HistName%BoundLeft  
                                                 ;  STOP  'GetHistogramIBox: x out of range, Left  '  ;  endif
        if(ValBoxInd > HistName%BoundRight) then ;  print*, ValBoxInd, HistName%BoundRight 
                                                 ;  STOP  'GetHistogramIBox: x out of range, Right '  ;  endif

        DistFromLeftBound = Abs(ValBoxInd - HistName%BoundLeft)
        DistInBoxLengths  = DistFromLeftBound/ HistName%BoxLength
        IBox              = INT(DistInBoxLengths) + 1 
        !Avoid 'array out of bounds' case when x = BoundRight, or BoundLeft (IBox indetermined).
        if(IBox > HistName%NBox)    IBox = HistName%NBox
        !print*, IBox, DistFromLeftBound, HistName%BoundRight , x
        !print*, IBox, HistName%NBox, x, DistFromLeftBound, DistInBoxLengths 
     END FUNCTION GetHistogramIBox

END MODULE Histograms