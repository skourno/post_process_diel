!----------------------------------------------------------------------------------
!                               flow.f90  -  description
!                             ------------------------------
!                  begin    : Fri 16.11.2016
!                  author   : Spiros Kournopoulos
!                  email    : sk1916@ic.ac.uk, skournopoulos@msn.com
!----------------------------------------------------------------------------------
!
! This module is designed to provide an interface between the input of an analyze
! program and the calculation routines.
!
!----------------------------------------------------------------------------------
!
!  17.01.17   1.0   Skou    FrameAct, Action 
!
!----------------------------------------------------------------------------------

Module Flow
  Use nrtype
  Use g_r
  Implicit None

Contains

  !------------------------------------------------------------------------------------------------------
  ! Called once for each frame.   ISwitch = (1, 2, 3) : (Initialize, Collect, Write)
  !
  !------------------------------------------------------------------------------------------------------
  Subroutine FrameAct(ISwitch)
    Use Input,    Only : NOptions, IOption
    Implicit None
    Integer(I4B) :: i, ISwitch
    Intent(IN)   :: ISwitch

    do i=1, NOptions
      if(IOption(i).NE.0) then
        Call Action(i, ISwitch)
      endif
    enddo 
End Subroutine FrameAct

!---------------------------------------------------------------------------------
! Leads to the correct action based on the input directives
!---------------------------------------------------------------------------------
Subroutine Action(i,ISwitch)
  Use Input,  Only : NBins, BLeft, BRight, IOption
  Use System, Only : Atom
  Implicit None
  Integer(I4B) :: i, ISwitch
  Intent(IN)   :: i, ISwitch


  if(.NOT. Allocated(Atom)) STOP 'Action: g_r: system not loaded'                

  Select Case (i)


    Case(1)  ; Call g_r_Dist(ISwitch, NBins(i), BLeft(i), BRight(i))
    Case(2)  ; Call g_r_Dist(ISwitch, NBins(i), BLeft(i), BRight(i))   ! this calculation - g110 - requires multiple histograms (g(r) and Ylm)
             ; Call Ylm_Dist(ISwitch, NBins(i), BLeft(i), BRight(i),"1,+1") 
             ; Call Ylm_Dist(ISwitch, NBins(i), BLeft(i), BRight(i),"1,+0") 
             ; Call Ylm_Dist(ISwitch, NBins(i), BLeft(i), BRight(i),"1,-1") 
    Case(3)  ; Call gamma12_r_Dist(ISwitch, NBins(i), BLeft(i), BRight(i))

    Case Default  ;  STOP 'Action: Wrong Call'
  End Select
End Subroutine Action


End Module Flow
