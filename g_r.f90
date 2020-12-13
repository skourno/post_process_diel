!----------------------------------------------------------------------------------
!                                g_r.f90  -  description
!                             ------------------------------
!                  begin    : Tue 16.01.2017
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
! This module contains g(r) distribution routines.
!----------------------------------------------------------------------------------
!
!  17.01.17   1.0   Skou    g(r)    distribution
!  17.06.20   2.0   Skou    g110(r) distribution     
!
!----------------------------------------------------------------------------------
    
Module g_r
  Use nrtype
  Implicit None

Contains

!-------------------------------------------------------------------------------------------------------------------
! Calculates total and inter g(r) radial distribution function (RDF) 
!-------------------------------------------------------------------------------------------------------------------
  Subroutine g_r_Dist(ISwitch, NBox, BoundLeft, BoundRight)
    Use Input,     Only: CFileOutput
    Use Utilities, Only: DistMinIm
    Use System,    Only: Atom, AtomType, Topo, Cell
    Use Histograms
    Implicit None
    Type(Histogram), Save :: g_r 
    Integer(I4B)          :: ISwitch, NBox, NFr, i, iAt, jAt, n, DNCount
    Real(DP)              :: BoundLeft, BoundRight, Dist, AvgCellEdge(3), AvgVol, DNormInter, r
    Optional              :: NBox, BoundLeft, BoundRight
    Intent(IN)            :: ISwitch, NBox, BoundLeft, BoundRight
    Save                  :: NFr, AvgCellEdge, DNCount


    Select Case (ISwitch)    
    
    Case(1)
      if(.NOT. (Present(NBox) .AND. Present(BoundLeft) .AND. Present(BoundRight))) & 
         STOP 'g_r_Distr: Initialization problem'
      
      g_r            = NewHistogram(NBox, BoundLeft, BoundRight) 
      NFr            = 0
      DNCount        = 0
      AvgCellEdge(:) = ZERO

    Case(2)
      NFr            = NFr + 1                   
      
      !loop over all atoms
      do iAt= 1, Topo%NAt
        do jAt = iAt+1, Topo%NAt
          Dist     = DistMinIm(Atom(iAt)%R, Atom(jAt)%R, Cell%Edge)
          Call PutInHistogram(g_r, Dist)
          DNCount  = DNCount + 1
        enddo
      enddo   

      AvgCellEdge(:) = (AvgCellEdge(:)*(NFr-1) + Cell%Edge(:)) / NFr
            
    Case(3)      

      ! we have to normalize with the atoms we expect to be in a shell. However, as we get further away
      ! from a particle the spherical shells become bigger.
      ! Thus, the we take into account the correct shell volume.
      AvgVol      = Product(AvgCellEdge)        !Avg Volume                
    
      DNormInter  = (DNCount / NFr) * (4.D0/3.D0)* PI / AvgVol
      if(ANY(g_r%HistArr /= ZERO)) &
         g_r%HistArr    = (g_r%HistArr/DBLE(NFr))/DNormInter
         
      do i =1, g_r%NBOX
         r              = DBLE(i)*g_r%BoxLength
         g_r%HistArr(i) = g_r%HistArr(i)/(r**3 -((r - g_r%BoxLength)**3))  

      end do 

      Call WriteHistogram(g_r, Trim(CFileOutput)//'_g_r_all.txt', 0)

    End Select        
   
  End Subroutine g_r_Dist


!-------------------------------------------------------------------------------------------------------------------
! Calculates <Ylm* Yl-m*>(r) histograms. The last argument is a string that specifies the type of lm call. At this point
! the following options are supported:
!
! 1,+1 : l=1, m=+1
! 1,+0 : l=1, m=+0
! 1,-1 : l=1, m=-1
!
!-------------------------------------------------------------------------------------------------------------------
  Subroutine Ylm_Dist(ISwitch, NBox, BoundLeft, BoundRight, Clm)
    Use Input,         Only: CFileOutput
    Use Utilities,     Only: DistMinIm, Get_Min_Img_Conn_Vec, Rotate_z_axis_to_align_with
    Use Sph_Harmonics, Only: SphHarm_xyz
    Use System,        Only: Atom, AtomType, Topo, Cell
    Use Histograms
    Implicit None
    Type(Histogram), Allocatable, Save :: Yprod(:,:), ph
    Real(DP),        Allocatable, Save :: Yprod_ave(:,:,:,:) ! two first indices are l,m; last two are bin value and counter
    Integer(I4B)          :: ISwitch, NBox, NFr, i, iAt, jAt, n, DNCount
    Integer(I4B)          :: l, m, IBox
    Real(DP)              :: BoundLeft, BoundRight, Dist, AvgCellEdge(3), AvgVol, DNormInter, r, sphCoords(3)
    Real(DP)              :: DistFromLeftBound, DistInBoxLengths, r12(3), unitVecsRelative(3,3), dMom1_rel(3), dMom2_rel(3)
    Complex(DPC)          :: Ylpm, Ylnm, ComplexAve(-1:1)
    Logical               :: LFirstEntry = .TRUE.
    Optional              :: NBox, BoundLeft, BoundRight
    Character*4           :: Clm
    Intent(IN)            :: ISwitch, NBox, BoundLeft, BoundRight, Clm
    Save                  :: NFr, AvgCellEdge, DNCount, LFirstEntry, ComplexAve


    Select Case (ISwitch)    
    Case(1)

      if (LFirstEntry) then
        if (Allocated(Yprod)) Deallocate(Yprod)
        Allocate(Yprod(1,-1:1))


        Allocate(Yprod_ave(1,-1:1,NBox,2))
        Yprod_ave = ZERO
        LFirstEntry = .FALSE.
      endif

      if(.NOT. (Present(NBox) .AND. Present(BoundLeft) .AND. Present(BoundRight))) & 
         STOP 'Ylm_Distr: Initialization problem'
      
      Select Case(Clm)
      Case("1,+1")
        Yprod(1,+1) = NewHistogram(NBox, BoundLeft, BoundRight) 
      Case("1,+0")
        Yprod(1,+0) = NewHistogram(NBox, BoundLeft, BoundRight) 
      Case("1,-1")
        Yprod(1,-1) = NewHistogram(NBox, BoundLeft, BoundRight) 
      Case Default 
        STOP "Ylm_Distr: Invalid l,m indices" 
      End Select 

      NFr            = 0
      DNCount        = 0
      ComplexAve     = (ZERO,ZERO)

    Case(2)
      NFr            = NFr + 1                   

      Select Case(Clm)
      Case("1,+1") ; l = 1 ; m = +1
      Case("1,+0") ; l = 1 ; m = +0
      Case("1,-1") ; l = 1 ; m = -1
      Case Default 
        STOP "Ylm_Distr: Invalid l,m indices" 
      End Select 

      
      !loop over all atoms
      do iAt= 1, Topo%NAt
        do jAt = iAt+1, Topo%NAt
          DNCount  = DNCount + 1

          Dist     = DistMinIm(Atom(iAt)%R, Atom(jAt)%R, Cell%Edge)

          ! we KNOW that all the Ylm histos will have the same bounds. As we set it this way abovr
          if( dist < Yprod(l,m)%BoundLeft ) CYCLE
          if( dist > Yprod(l,m)%BoundRight) CYCLE 

          DistFromLeftBound = Abs(dist - Yprod(l,m)%BoundLeft)
          DistInBoxLengths  = DistFromLeftBound / Yprod(l,m)%BoxLength
          IBox              = INT(DistInBoxLengths) + 1 

          !Avoid 'array out of bounds' case when x = BoundRight, or BoundLeft (IBox indetermined).
          if(IBox > Yprod(l,m)%NBox)    IBox = Yprod(l,m)%NBox

          ! We need to change the frame of reference to be the axis between the two atoms
          call Get_Min_Img_Conn_Vec(Atom(iAt)%R, Atom(jAt)%R, r12, Cell%Edge)
          unitVecsRelative  = Rotate_z_axis_to_align_with(r12)

          dMom1_rel(1) = DOT_PRODUCT(unitVecsRelative(:,1), Atom(iAt)%dMom)
          dMom1_rel(2) = DOT_PRODUCT(unitVecsRelative(:,2), Atom(iAt)%dMom)
          dMom1_rel(3) = DOT_PRODUCT(unitVecsRelative(:,3), Atom(iAt)%dMom)

          dMom2_rel(1) = DOT_PRODUCT(unitVecsRelative(:,1), Atom(jAt)%dMom)
          dMom2_rel(2) = DOT_PRODUCT(unitVecsRelative(:,2), Atom(jAt)%dMom)
          dMom2_rel(3) = DOT_PRODUCT(unitVecsRelative(:,3), Atom(jAt)%dMom)

          Ylpm              = SphHarm_xyz(dMom1_rel, l, +m)
          Ylnm              = SphHarm_xyz(dMom2_rel, l, -m)

          Ylpm              = CONJG(Ylpm)
          Ylnm              = CONJG(Ylnm)
          
          Yprod_ave(l,m,IBox,1) = Yprod_ave(l,m,IBox,1) + Ylpm*Ylnm
          Yprod_ave(l,m,IBox,2) = Yprod_ave(l,m,IBox,2) + 1
          ComplexAve(m)         = ComplexAve(m) + Ylpm*Ylnm

          !if (IBox == NBox) print*, m, ComplexAve(m) / Yprod_ave(l,m,IBox,2)
        enddo
      enddo   

      AvgCellEdge(:) = (AvgCellEdge(:)*(NFr-1) + Cell%Edge(:)) / NFr
            
    Case(3)  

      Select Case(Clm)
      Case("1,+1") ; l = 1 ; m = +1
      Case("1,+0") ; l = 1 ; m = +0
      Case("1,-1") ; l = 1 ; m = -1
      Case Default 
        STOP "Ylm_Distr: Invalid l,m indices" 
      End Select 

      do i =1, Yprod(l,m)%NBOX

         if (Yprod_ave(l,m,i,2) == 0) then  
           Yprod(l,m)%HistArr(i) = ZERO
           CYCLE
         endif

         Yprod(l,m)%HistArr(i) = Yprod_ave(l,m,i,1) / Yprod_ave(l,m,i,2)   
      enddo 


      Select Case(Clm)
      Case("1,+1") ; Call WriteHistogram(Yprod(l,m), Trim(CFileOutput)//'_Y1p1.txt', 0)
      Case("1,+0") ; Call WriteHistogram(Yprod(l,m), Trim(CFileOutput)//'_Y1p0.txt', 0)
      Case("1,-1") ; Call WriteHistogram(Yprod(l,m), Trim(CFileOutput)//'_Y1n1.txt', 0)
      Case Default 
        STOP "Ylm_Distr: Invalid l,m indices" 
      End Select 

    End Select        
  End Subroutine Ylm_Dist


!-------------------------------------------------------------------------------------------------------------------
! Calculates the cos gamma12 (r) histogram 
!-------------------------------------------------------------------------------------------------------------------
  Subroutine gamma12_r_Dist(ISwitch, NBox, BoundLeft, BoundRight)
    Use Input,     Only: CFileOutput
    Use System,    Only: Atom, AtomType, Topo, Cell
    Use Utilities, Only: Angle, DistMinIm
    Use Histograms
    Implicit None
    Real(DP), Allocatable, Save  :: gamma12(:,:)
    Type(Histogram),       Save  :: gamma12_r 
    Integer(I4B)          :: ISwitch, NBox, NFr, i, iAt, jAt, IBox
    Real(DP)              :: BoundLeft, BoundRight, Dist
    Real(DP)              :: cosgamma12, dip1(3), dip2(3)
    Real(DP)              :: DistFromLeftBound, DistInBoxLengths

    Optional              :: NBox, BoundLeft, BoundRight
    Intent(IN)            :: ISwitch, NBox, BoundLeft, BoundRight
    Save                  :: NFr

    Select Case (ISwitch)    
    
    Case(1)
      if(.NOT. (Present(NBox) .AND. Present(BoundLeft) .AND. Present(BoundRight))) & 
         STOP 'gamma12_r_Dist: Initialization problem'
      
      gamma12_r      = NewHistogram(NBox, BoundLeft, BoundRight) 
      NFr            = 0

      Allocate(gamma12(NBox,2))
      gamma12(:,:)   = ZERO
    Case(2)
      NFr            = NFr + 1                   
      
      !loop over all atoms
      do iAt= 1, Topo%NAt
        do jAt = iAt+1, Topo%NAt

          dip1 = Atom(iAt)%dMom(:)
          dip2 = Atom(jAt)%dMom(:)

          Dist     = DistMinIm(Atom(iAt)%R, Atom(jAt)%R, Cell%Edge)


          if( dist < gamma12_r%BoundLeft ) CYCLE
          if( dist > gamma12_r%BoundRight) CYCLE 

          cosgamma12 = COS( Angle(dip1,dip2) )

          DistFromLeftBound = Abs(dist - gamma12_r%BoundLeft)
          DistInBoxLengths  = DistFromLeftBound/ gamma12_r%BoxLength
          IBox              = INT(DistInBoxLengths) + 1 

          !print*, IBox, Dist, cosgamma12

          !Avoid 'array out of bounds' case when x = BoundRight, or BoundLeft (IBox indetermined).
          if(IBox > gamma12_r%NBox)    IBox = gamma12_r%NBox
          !print*, IBox, DistFromLeftBound, gamma12_r%BoundRight , x
          !print*, IBox, gamma12_r%NBox, x, DistFromLeftBound, DistInBoxLengths 

          gamma12(IBox,1) = gamma12(IBox,1) + cosgamma12
          gamma12(IBox,2) = gamma12(IBox,2) + 1

        enddo
      enddo   
            
    Case(3)      
      ! we have to normalize with the atoms we counted to be in the sphericql shell
    
      do i =1, gamma12_r%NBOX

         if (gamma12(i,2) == 0) then  
           gamma12_r%HistArr(i) = ZERO
           CYCLE
         endif

         gamma12_r%HistArr(i) = gamma12(i,1) / gamma12(i,2)   

      end do 

      Call WriteHistogram(gamma12_r, Trim(CFileOutput)//'_cosgamma_r.txt', 0)

    End Select        
   
  End Subroutine gamma12_r_Dist


!!-------------------------------------------------------------------------------------------------------------------
!! Calculates the g110(r) distribution. This calculation involves the calculation of the g(r), so avoid running both
!! at the same time
!!-------------------------------------------------------------------------------------------------------------------
!  Subroutine g110_r_Dist(ISwitch, NBox, BoundLeft, BoundRight)
!    Use Input,     Only: CFileOutput
!    Use Utilities, Only: DistMinIm
!    Use System,    Only: Atom, AtomType, Topo, Cell
!    Use Histograms
!    Implicit None
!    Type(Histogram), Save :: g110_r
!    Integer(I4B)          :: ISwitch, NBox, NFr, i, iAt, jAt, n, DNCount
!    Real(DP)              :: BoundLeft, BoundRight, Dist, AvgCellEdge(3), AvgVol, DNormInter, r
!    Optional              :: NBox, BoundLeft, BoundRight
!    Intent(IN)            :: ISwitch, NBox, BoundLeft, BoundRight
!    Save                  :: NFr, AvgCellEdge, DNCount
!
!    Select Case (ISwitch)    
!    
!    Case(1)
!      if(.NOT. (Present(NBox) .AND. Present(BoundLeft) .AND. Present(BoundRight))) & 
!         STOP 'g110_r_Dist : Initialization problem'
!
!      Call g_r_Dist(1, NBox, BoundLeft, BoundRight)
!      
!      g110_r         = NewHistogram(NBox, BoundLeft, BoundRight) 
!      NFr            = 0
!      DNCount        = 0
!
!    Case(2)
!      NFr            = NFr + 1   
!
!      Call g_r_Dist(2, NBox, BoundLeft, BoundRight)
!
!      !loop over all atoms
!      do iAt= 1, Topo%NAt
!        do jAt = iAt+1, Topo%NAt
!
!          !do i = 
!
!
!
!
!          !sum      = ...      
!          Call PutInHistogram(g110_r, sum)
!          DNCount  = DNCount + 1
!        enddo
!      enddo   
!
!      AvgCellEdge(:) = (AvgCellEdge(:)*(NFr-1) + Cell%Edge(:)) / NFr
!            
!    Case(3)      
!      Call g_r_Dist(3, NBox, BoundLeft, BoundRight)
!
!      ! we have to normalize with the atoms we expect to be in a shell. However, as we get further away
!      ! from a particle the spherical shells become bigger.
!      ! Thus, the we take into account the correct shell volume.
!      AvgVol      = Product(AvgCellEdge)        !Avg Volume                
!      DNormInter  = (DNCount / NFr) * (4.D0/3.D0)* PI / AvgVol
!    
!      if(ANY(g_r%HistArr /= ZERO)) &
!         g110_r%HistArr    = (g110_r%HistArr/DBLE(NFr))/DNormInter
!         
!      do i =1, g_r%NBOX
!         r                 = DBLE(i)*g_r%BoxLength
!         g110_r%HistArr(i) = g110_r%HistArr(i)/(r**3 -((r - g110_r%BoxLength)**3))    
!      end do 
!
!      Call WriteHistogram(g_r, Trim(CFileOutput)//'_g110_r.txt', 0)
!
!    End Select        
!   
!  End Subroutine g_r_Dist

End Module g_r
