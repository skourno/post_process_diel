program main
	Use nrtype
	Use input,   Only : IStartFrame, ILastFrame, NFramesStride ! variables
	Use input,   Only : ReadInput  , Read_Sys                      !
	Use system,  Only : Atom, AtomType, Cell, Topo
	Use compute, Only : compute_dielectric, compute_ave_relative_P1
	Use Flow,    Only : FrameAct	

	Implicit None 

	Integer(I4B) :: iFr, NAnalyzed
	Real(DP)     :: dmom(3), dmom_ave(3), dmom_ave_sq
	Real(DP)     :: dmom_sq, dmom_sq_ave
	Logical      :: LFirstEntry = .TRUE.


	call ReadInput()

	dmom_ave(:)      = ZERO
	dmom_sq_ave      = ZERO
	NAnalyzed        = 0
    
	print*, '!-------------------------------------------'
	print*, '! Dielectric Post-Process version 2.0 (dev)'
	print*, '!-------------------------------------------'


    !Read the first 'IStartFrames-1' frames which will not be processed
    do iFr=1, IStartFrame-1
       Call Read_Sys()
       print*, 'SKIPPING FRAME', iFr
    enddo


    !Read and process the frames up to 'NFramesInput' with a stride of 'NFramesRelaxInterval'
	do iFr=IStartFrame, ILastFrame
		call Read_Sys()

		if (MOD(iFr, NFramesStride).NE. 0) then
			print*, 'SKIPPING FRAME', iFr
			CYCLE 
		endif

		call compute_dielectric(Atom, AtomType, Cell, Topo, dmom_sq, dmom)

		dmom_ave    = dmom_ave    + dmom
		dmom_sq_ave = dmom_sq_ave + dmom_sq
		NAnalyzed   = NAnalyzed   + 1

		write(*,111) iFr, NAnalyzed, dmom_sq, dmom
        
        !Initialize histograms.
        if(LFirstEntry) then
          Call FrameAct(1)
          LFirstEntry = .FALSE.
        endif
        !Collect Histogram data
        Call FrameAct(2)

        !print*, 'PROCESSED FRAME', iFr
    enddo

    !Write Histograms
    Call FrameAct(3)
    write(*,*) '###'
    write(*,*) '### Histo output written ###'
    write(*,*) '###'


	

	dmom_sq_ave = dmom_sq_ave / NAnalyzed
	dmom_ave    = dmom_ave    / NAnalyzed
	dmom_ave_sq = DOT_PRODUCT(dmom_ave, dmom_ave)

	write(*,*) '##### Final result #####'
	write(*,112) dmom_sq_ave, dmom_ave_sq, dmom_ave_sq / dmom_sq_ave

	111 FORMAT(I3, 2X, I3, 5X, 4(F12.3))
	112 FORMAT(2(F12.3), F14.6) 

end program main