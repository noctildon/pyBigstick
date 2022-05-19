!===================================================================
!
!  subroutine menu
!
!  master routine to set up flags to direct flow of program
!
! SUBROUTINES CALLED:
!	wfn_ropen_file  
!     	call read_wfn_header
!........................................................................
subroutine menu

  use menu_choices
  use io
  use verbosity
  use nodeinfo
  use flagger
  use wfn_mod
  use bmpi_mod
  use jumpNbody
  use onebodypot
  use densities, only : pndensities
  use coupledmatrixelements,only: dens2bflag,diag_den2b
  use configurations,only:configout
  use flags3body,only:threebodycheck
  use TRstuff
  use jumpstart
  use tracy
  implicit none

  integer(4) :: ierr
  character :: ychar

! File to log interactive inputs to so the user can create an input file
  open(unit=autoinputfile,file='autoinput.bigstick',status='unknown')  ! always do this

  auto_input = .false.
  pndensities = .false.
  dens2bflag = .false.
  diag_den2b = .false.
  useTRphase = useTRphase_def
  dotflag = .false.
  
  if ( iproc == 0 ) then
     print*,' '
     print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
     print*,' *                                                                         * '
     print*,' *               OPTIONS (choose one)                                      * '
     print*,' * (i) Input automatically read from "autoinput.bigstick" file             * '
     print*,' *  (note: autoinput.bigstick file created with each nonauto run)          * '
     print*,' * (n) Compute spectrum (default); (ns) to suppress eigenvector write up   * '
     print*,' * (j) Carry out jumpstart run to set MPI timing (MPI only)                * '
     print*,' * (d) Densities: Compute spectrum + all one-body densities                * '
     print*,' * (2) Two-body density from previous wfn (default p-n format)             * '

!     print*,' * (dx[m]) Densities: Compute one-body densities from previous run (.wfn)  * '
!     print*,' *     optional m enables mathematica output                               * '
!     print*,' * (p) Compute spectrum + single-particle occupations; (ps) to suppress wfn* '
!     print*,' * (occ) single-particle occupations (from previous wfn)                   * '
     print*,' * (x) eXpectation value of a scalar Hamiltonian (from previous wfn)       * '
     print*,' * (o) Apply a one-body (transition) operator to previous wfn and write out* '
     print*,' * (s) Strength function (using starting pivot )                           * '
!     print*,' * (r) Restart from a previous run                                         * '
!     print*,' * (a) Apply a scalar Hamiltonian to a previous wfn and write out          * '
     print*,' * (g) Apply the resolvant 1/(E-H) to a previous wfn and write out         * '  ! ADDED 7.8.8; restarted 7.9.10
	 
!     print*,' * (y) Cycle back to this menu after calculation completed (not active)    * '
!     print*,' * (v) Overlap of initial states with final states                         * '
     print*,' * (m) print information for Modeling parallel distribution                * '
!     print*,' * (b) Create file with full basis information (for postprocessing)        * '
!     print*,' * (f) Self-consistent mean-field approximation (prepare pivot)            * '    ! ADDED 7.7.4
     print*,' * (l) print license and copyright information                             * '
	 
     print*,' * (?) Print out all options                                               * '
     print*,' *                                                                         * '
     print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
  end if
1 continue

  if( auto_input)then
     read(autoinputfile,'(a)')menu_char
  else  if ( iproc == 0 ) then
     print*,' '
     print*,' Enter choice '
     read(5,'(a)')menu_char
     if(menu_char /= 'i' .and. menu_char /= 'I')then
         write(autoinputfile,'(a,"    ! menu choice ")')menu_char
     end if
  end if

  call BMPI_BARRIER(icomm,ierr)
  call BMPI_BCAST(menu_char,3,0,icomm,ierr)

  auto_readin = .false.
  ham_readin = .false.
  op_readin  = .false.
  print4modelinfo = .false.
  strengthflag = .false.
  greenflag = .false.
  densityflag = .false.
  trdensout   = .false.
  spoccflag   = .false.
  menu_dx_omathematica = .true.  ! flag enabling mathematica output of density matrices
  get_JT       = .true.
  onebodyonly  = .false.
  printouthamflag = .false.
  printouttrans1flag = .false.
  baseASCII=.false.
  configout = .false.
  modeldensities=.false.
  
  threebodycheck=.false.   ! a beta version in 7.9.2
  blockstrengthflag=.false.
  
  writejumpstart = .false.
  readjumpstart  = .false.
  strengthnorm = .true.
  dotflag = .true.
  
  
  select case (menu_char)

  case ('?')
     if ( iproc == 0 ) then
        print*,' '
        print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
        print*,' *                                                                         * '
        print*,' *               OPTIONS (choose 1)                                        * '
        print*,' * (i) Input automatically read from "autoinput.bigstick" file             * '
        print*,' *  (note: autoinput.bigstick file created with each nonauto run)          * '
        print*,' * (n) Compute spectrum (default); (ns) to suppress eigenvector write up   * '
        print*,' * (ne) Compute energies but NO observable (i.e. J or T)                   * '   !ADDED 7.6.8
        print*,' * (np) Compute spectrum starting from prior pivot (first N vectors)       * '   !ADDED 7.7.7. modified in 7.9.11
!        print*,' * (nc) Compute spectrum starting from prior pivot, choose vectors         * '   !ADDED 7.9.11
        print*,' * (j) Carry out jumpstart run to set MPI timing (MPI only)                * '
        print*,' * (d) Densities: Compute spectrum + all one-body densities                * '
        print*,' * (dx[m]) Densities: Compute one-body densities from previous run (.wfn)  * '
        print*,' *     optional m enables mathematica output                               * '
        print*,' * (dp) Densities in proton-neutron format                                 * '
        print*,' * (dxp) Compute one-body densities from prior run (.wfn) in p-n format.   * '
		
        print*,' * (2) Two-body density from previous wfn (default p-n format)             * '  ! ADDED 7.9.1
        print*,' * (2d) Two-body density from previous wfn, only initial=final, Jt=0       * '  ! added 7.9.5
!        print*,' * (2i) Two-body density from previous wfn (alternate isopin format)       * '  ! NOT YET IMIPLEMENTED
        print*,' * (3) Normal spectrum but using three-body forces (beta version)          * '
		
		
        print*,' * (p) Compute spectrum + single-particle occupations,(ps) to supress wfn  * '
        print*,' * (occ) single-particle occupations (from previous wfn)                   * '
        print*,' * (x) eXpectation value of a scalar Hamiltonian (from previous wfn)       * '
        print*,' * (o) Apply a one-body (transition) operator to previous wfn and write out* '
        print*,' * (oo) Apply a one-body (transition) operator with enforced orthogonality * '
        print*,' *     (that is, forced resulting pivot to be orthogonal to starting pivot * '
        print*,' *      to eliminate transition strength to original state)                * '
        print*,' * (s),(sn) Strength function (using starting pivot ) wfn out normalized   * '
        print*,' * (su) Strength function (using starting pivot ) wfn out unnormalized     * '
        print*,' * (sb) Strength function (using block of starting pivots )                * '
!        print*,' * (r) Restart from a previous run                                         * '  ! NOW OBSOLETE
        print*,' * (a) Apply a scalar Hamiltonian to a previous wfn and write out          * '
        print*,' * (g) Apply resolvent 1/(E-H) to a previous wfn and write out             * '  ! ADDED 7.8.8; restarted 7.9.10
        print*,' * (gv) Apply resolvent 1/(E-H) to a previous wfn, then take dot prod      * '  ! ADDED 7.9.10

!       print*,' * (y) Cycle back to this menu after calculation completed (not active)    * '
        print*,' * (v) Overlap of initial states with final states                         * '
        print*,' * (vs) Relative entropy between initial states and final states                         * '
		
        print*,' * (m) print information for Modeling parallel distribution                * '
        print*,' * (md) Modeling parallel distribution for 1-body densities                * '
        print*,' * (m3) print information for Modeling parallel distribution w/3-body force* '
!        print*,' * (f) Self-consistent mean-field approximation (prepare pivot)            * ' NOW OBSOLETE AND DEPRECATED
        print*,' * (t) create TRDENS-readable file for post processing                     * '
		print*,' * (b) Create binary file with full basis information (for postprocessing) * '
		print*,' * (ba) Create ASCII file with full basis information (for postprocessing) * '  ! ADDED 7.9.1
		
        print*,' * (wh) write out Hamiltonian matrix to a file and stop                    * '  ! ADDED 7.8.5
        print*,' * (wo) write out one-body transition matrix to a file and stop            * '  ! ADDED 7.8.5
!        print*,' * (wm) write out uncoupled 1+2-body operator m.e.s to file and stop       * '  ! ADDED 7.9.12

        print*,' * (co) Compute configurations (partitions)                                * '   ! ADDED 7.9.1

        print*,' * (c) Compute traces                                                      * '
        print*,' * (l) print license and copyright information                             * '
		
!       print*,' * (l) Time Hamiltonian loops  (obsolete?)                                 * '
        print*,' * (?) Print out all options                                               * '
    
        print*,' *                                                                         * '
        print*,' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ' 
     end if
     goto 1
     
  case ('i','I')
     if(iproc == 0)print*,' Reading automatically from "autoinput.bigstick" file '
     auto_input = .true.
     goto 1
  case ('d','D')
     menu_char = 'd'
     ham_readin = .true.
     if ( iproc == 0 )print*,' Compute spectrum + all one-body densities ' 
     densityflag = .true.
     write_wfn = .true.
  case ('dp','DP')
     menu_char = 'd'
     ham_readin = .true.
     if ( iproc == 0 )print*,' Compute spectrum + all one-body densities in proton-neutron format ' 
     densityflag = .true.
     write_wfn = .true.
	 pndensities = .true.
  case ('dx','DX')
     menu_dx_omathematica = .false.
     menu_char = 'dx'
     auto_readin = .true.
     ham_readin = .false.
     if ( iproc == 0)then
		 print*,' Compute one-body densities from previous wfn ' 
		 print*,' Output written to .dres file in isospin format '
		 print*,' NEW: You will need to select initial, final range '
		
	 end if
     densityflag = .true.
	 onebodyonly =.true.
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
	 
  case ('dxp','DXP')
     menu_dx_omathematica = .false.
     menu_char = 'dx'
     auto_readin = .true.
    ham_readin = .false.
   	 pndensities = .true.
		
        if ( iproc == 0)then
   		 print*,' Compute one-body densities from previous wfn ' 
   		 print*,' Output written to .dres file in proton-neuton format '
   		 print*,' NEW: You will need to select initial, final range '
		
   	 end if
     densityflag = .true.
   	 onebodyonly =.true.
        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)	 
	 
  case ('dxm', 'DXM')
        menu_dx_omathematica = .true.
        menu_char = 'dx'
        auto_readin = .true.
        ham_readin = .false.
        if ( iproc == 0 )print*,' Compute one-body densities from previous wfn ' 
        if ( iproc == 0 )print*,' and write out in Mathematica-readable format ' 

        densityflag = .true.
   	 onebodyonly =.true.
        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)	 
  
  case ('2','2i','2d')                  ! ADDED IN 7.9.1  -- TWO BODY DENSITIES FROM EXISTING WFNs
     if(menu_char=='2i')then
	    pndensities=.false.
     else
	    pndensities = .true.
     end if
	 if(menu_char=='2d' .or. menu_char=='2i')then
		 diag_den2b = .true.
	 end if
	 
     menu_char = '2b'
     auto_readin = .true.
     ham_readin = .false.
     if ( iproc == 0 )then
		 print*,' Compute two-body densities from previous wfn ' 
!		 if(nproc> 1)print*,' THIS DOES NOT YET WORK IN  MPI '
		 
		 if(diag_den2b)then
			 print*,' (Only for initial=final, Jt=0)'
			 if(pndensities)then
				 print*,' (Computed in proton-neutron formalism )'
			 else
				 print*,' (Computed in isospin formalism )'
			 end if
		 end if
	 end if
     densityflag = .true.
	 dens2bflag = .true.
	 useTRphase = .false.  ! for now
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)	 

  case ('b','B','ba','BA')         ! ADDED IN 7.7.2
     if ( iproc == 0 )print*,' Create basis (.bas) file for later postprocessing ' 

     if(menu_char =='ba' .or. menu_char=='BA') then
		 if(nproc > 1)then
			 if(iproc==0)then
				 print*,' Sorry, this option will not work in MPI '
			 end if
		     call BMPI_Abort(icomm,101,ierr)
			 stop
		 end if
		 baseASCII=.true.
		 if(iproc==0)print*,' (writing out in human-readable ASCII)'
	 end if
     menu_char = 'b'
	 writeout=.true.

!  case ('r','R')
!     menu_char = 'r'
!     auto_readin = .true.
!    ham_readin = .true.
!     call open_restart_files
!     call read_wfn_header(wfnfile,.true.)
!     if ( iproc == 0 ) then
!        print*,' Restart from previous run ' 
!     end if
     
  case ('a','A')
     menu_char = 'a'
     ham_readin = .true.
     auto_readin =.true.
     write_wfn = .true.
     if ( iproc == 0 ) then
        print*,' Apply a scalar Hamiltonian to a previous wfn and write out ' 
        print*,' (You must read in a previously computed wfn) '
     end if
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
     
  case ('o','O')
     menu_char = 'o'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.true.
     write_wfn = .true.
	 onebodyonly = .true.

     if ( iproc == 0 )print*,' Apply a one-body operator to a previous wfn and write out ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('oo','OO')
     menu_char = 'oo'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.true.
     write_wfn = .true.
	 onebodyonly = .true.

     if ( iproc == 0 )print*,' Apply a one-body operator to a previous wfn and write out ' 
     if ( iproc == 0 )print*,' (Will force resulting wfn to be orthogonal to original) ' 

     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('wo','WO')    ! write out one-body operator matrix to file trans1.dat '
        menu_char = 'wo'
        ham_readin = .false.
        auto_readin = .true.
        op_readin  =.true.
        write_wfn = .true.
   	 onebodyonly = .true.
     printouttrans1flag = .true.
	 

        if ( iproc == 0 )then
			print*,' Writing one-body transition operator to file trans1.dat '
			print*,' Kluge: must read in existing .wfn file to set up '
		end if
		 
        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)
	

  case ('v','V','vs','VS')
     dotflag = .true.
     if(menu_char=='vs' .or. menu_char=='VS')dotflag=.false.
     menu_char = 'v'
     ham_readin = .false.
     auto_readin = .true.
     op_readin  =.false.
     write_wfn = .false.
     writeout = .false.

     if ( iproc == 0 )then
	print*,' Read in a wavefunction file, select one initial state ' 
        print*,' then read in second waverfunction file and compute overlap with all final states '
        print*,' '
        print*,' INITIAL STATE '

     endif
     
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('x','X')
     menu_char = 'x'
     auto_readin = .true.
     ham_readin = .true.
     if ( iproc == 0 ) print*,' Compute expectation value of Hamiltonian operator using previous wfn ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('p','P', 'ps', 'PS')
     ! s suppresses writing wfn - takes a long time
	 if(.not. (menu_char(2:2) == 's' .or. menu_char(2:2) == 'S')) write_wfn = .true.
     menu_char = 'n'
     if ( iproc == 0 )print*,' Compute spectrum + single-particle occupations ' 
     ham_readin = .true.
     spoccflag = .true.

  case ('occ', 'OCC')
     menu_char = 'p'
     auto_readin = .true.
     ham_readin = .false.
     spoccflag = .true.
     if ( iproc == 0 ) print*,' Compute single-particle occupations using previous wfn ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)

  case ('m','M')
     menu_char = 'm'
     ham_readin = .false.
     if ( iproc == 0 )print*,' Print information needed to model parallel distribution ' 
     print4modelinfo = .true.

  case ('md','MD')
     menu_char = 'm'
	 modeldensities=.true.
	 densityflag = .true.
     ham_readin = .false.
     if ( iproc == 0 )print*,' Print information needed to model parallel distribution for 1-body densities' 
     print4modelinfo = .true.   
	 
  case ('m3','M3')
        menu_char = 'm'
        ham_readin = .false.
        if ( iproc == 0 )print*,' Print information needed to model parallel distribution ' 
        print4modelinfo = .true.
	    threebodycheck=.true.
		
  case ('s','S','sn','SN')
     menu_char = 's'
     ham_readin = .true.
     auto_readin = .true.
     strengthflag = .true.
	 strengthnorm = .true.
     if ( iproc == 0 )print*,' Compute strength function distribution using previous wfn ' 
	 

     if( iproc==0)then ! added 7.9.7
		 write(6,*)' Output wfns will  be normalized '
		 write(6,*)' To change this, use option (su) instead'
	 end if
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
     write_wfn = .true.

     case ('su','SU')
        menu_char = 's'
		strengthnorm = .false.
        ham_readin = .true.
        auto_readin = .true.
        strengthflag = .true.
        if ( iproc == 0 )print*,' Compute strength function distribution using previous wfn ' 
	 
        if(iproc==0)then ! added 7.9.7
   		 write(6,*)' Output wfns will not be normalized, but prior normalization x strength included '
   		 write(6,*)' To change this, use option (sn) instead '
   	 end if

        call wfn_ropen_file(oldwfnfile)
        call read_wfn_header(oldwfnfile,.true.)
        write_wfn = .true.	 
  
  case ('sb','SB')
     menu_char = 's'
     ham_readin = .true.
     auto_readin = .true.
     strengthflag = .true.
	 blockstrengthflag=.true.
     if ( iproc == 0 )print*,' Compute strength function distribution using block of previous wfns ' 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
     write_wfn = .true.
  	 
  case ('g','G','gv','GV','Gv')
     if(menu_char=='gv' .or. menu_char=='GV'.or. menu_char=='Gv')then
		dotflag=.true.
	 end if
     menu_char = 'g'
     ham_readin = .true.
     auto_readin = .true.
     greenflag = .true.
	 strengthnorm=.false.
     if ( iproc == 0 )print*,' Apply resolvant/Green function 1/(E-H) using previous wfn ' 
	 if(iproc==0 .and. dotflag)print*,' (You will have a chance to take a dot product with other wave functions at the end) '
	 
     call wfn_ropen_file(oldwfnfile)
     call read_wfn_header(oldwfnfile,.true.)
     write_wfn = .true.
	 
  case ('t','T')
     menu_char = 'n'
     ham_readin = .true.
     if ( iproc == 0 ) print*,' Will create output usable for TRDENS post-processing '
     trdensout = .true.
     write_wfn = .true.

  case ('ns','NS')
     menu_char = 'ns'
     if ( iproc == 0 ) then 
        print*,' Compute spectrum only '
        print*,' Suppress eigenvector write up '
     end if
     ham_readin = .true.
     write_wfn = .false.
	 
  case ('ne','NE')   ! ADDED 7.6.8
        menu_char = 'ne'
        if ( iproc == 0 ) then 
           print*,' Compute spectrum only '
           print*,' Suppress eigenvector write up '
		   print*,' No observables, no J, no T'
        end if
        ham_readin = .true.
        write_wfn = .false.	 
		get_JT     = .false.
		
  case('np','NP')   ! ADDED 7.7.7  -- option to read with a specified pivot 
      menu_char = 'np'
      if ( iproc == 0 )print*,' Compute spectrum starting from a prior pivot ' 
      ham_readin = .true.
      write_wfn = .true. 
      auto_readin = .true.
      call wfn_ropen_file(oldwfnfile)
      call read_wfn_header(oldwfnfile,.true.)

case ('wh','WH')   ! ADDED 7.8.5
      menu_char = 'wh'
      if ( iproc == 0 ) then 
         print*,' writing Hamiltonian to file ham.dat'

      end if
      ham_readin = .true.
      write_wfn = .false.	 
	get_JT     = .false.
    printouthamflag = .true.
	
    case ('co','CO')
       menu_char = 'n'
       if ( iproc == 0 ) then 
          print*,' Compute configurations '
       end if
       ham_readin = .true.
       write_wfn = .true.
	   configout = .true.
	  
  case ('c','C')
     menu_char = 'c'
     if ( iproc == 0 ) then 
        print*,' Compute traces (moments) '
		print*,' '
		print*,' Do you want variance (2nd moment?) (y/n)'
		if(auto_input)then
			read(autoinputfile,'(a)')ychar
		else
		    read(5,'(a)')ychar
			write(autoinputfile,'(a)')ychar
		end if
		if(ychar=='y' .or.ychar=='Y')then
			getvariance=.true.
			write(6,*)' Computing both centroid and variance '
		else
			getvariance=.false.
			write(6,*)' Computing only the centroid '
		end if
     end if
	 call BMPI_Bcast(getvariance,1,0,icomm,ierr)
     ham_readin = .true.
     write_wfn = .false.
	 threebodycheck=.true.

  case ('l','L')

  call printlicense
  goto 1

!  case ('f','F')
!        menu_char = 'f'
!		meanie=.true.
 !       if ( iproc == 0 )print*,' Running self-consistent mean-field ' 
  !      ham_readin = .true.
 !       write_wfn = .true.
!        get_JT    = .false.
		
  case ('3')
	menu_char = 'n'
	if ( iproc == 0 )print*,' Compute spectrum only ' 
	if ( iproc == 0 )print*,' Allowing three-body ' 
    threebodycheck=.true.
	ham_readin = .true.
	write_wfn = .true.
  case ('j','J')
    menu_char='ne'		   		
    if(iproc==0)then
		
		if(nproc==1)then
			print*,' '
			print*,' Warning! Cannot carry out jumpstart on 1 MPI process '
			print*,' You MUST run on the same number of MPI processes as for the full run '
			print*,' This is to obtain timing information '
			print*,' '
			write(logfile,*)' '
			write(logfile,*)' Warning! Cannot carry out jumpstart on 1 MPI process '
			write(logfile,*)' You MUST run on the same number of MPI processes as for the full run '
			write(logfile,*)' This is to obtain timing information '
			write(logfile,*)' '
			close(logfile)
			stop
		end if
			
		print*,'  ............................ '
		print*,'  Carrying out a jumpstart run '
		print*,'  ............................ '
		
		write(logfile,*)' ............................ '
		write(logfile,*)' Carrying out a jumpstart run '
		write(logfile,*)' ............................ '
	end if
	ham_readin = .true.
	
	writejumpstart = .true.
	readjumpstart  = .false.
	get_JT = .false.  ! suppress eigenvectors, getting J^2, T^2
  
  case default
     menu_char = 'n'
     if ( iproc == 0 )print*,' Compute spectrum only ' 
     ham_readin = .true.
     write_wfn = .true.
  end select

  return
end subroutine menu
!=====================================================================

subroutine printlicense
	
	print*,' '

    print*,'  This code is distributed under MIT Open Source License. '
	print*,' Copyright (c) 2017 Lawrence Livermore National Security '
	
	print*,' and San Diego State University Research Foundation  '
	print*,' '
			
	print*,' Permission is hereby granted, free of charge, to any person obtaining '
	print*,' a copy of this software and associated documentation files (the "Software"),'
    print*,' to deal in the Software without restriction, including without limitation'
	print*,' the rights to use, copy, modify, merge, publish, distribute, sublicense, '
	print*,' and/or sell copies of the Software, and to permit persons to whom the Software '
	print*,' is furnished to do so, subject to the following conditions: '
	print*,' '
	print*,' The above copyright notice and this permission notice shall be included '
	print*,' in all copies or substantial portions of the Software. '
	print*,' '
	print*,' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,'
	print*,' INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A '
	print*,' PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT '
	print*,' HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION '
	print*,' OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE '
	print*,' OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'

	print*,' '
	
	return
	
	
end subroutine printlicense
