!========================================================

!=====================================================================
!
!   COUNT UP STATISTICS OF OPERATIONS BETWEEN SECTORS AND BETWEEN "FRAGMENTS"
!   --- must do this before calculating the 
!   distribution of work across MPI compute nodes
!  NOTE: work betwen sectors is equivalent to work in an opbundle
!  as each opbundle is associated with a given "sector jump"
!  i.e. a set of jumps between sectors

! NOTE IN 7.7.1 these are not really needed because all the information is 
! contained in opbundles; therefore these are slated for obsolescence
! However for now drafting opbundles first checks this information;
! this must be replaced
!
! IMPORTANT: the operations counted here are NOT necessarily 
!  the operations in a bundle; bundles are defined by sector jumps
!  while these sum over *all* conjugate neutron sector jumps
!  Hence these routines are less useful.
!
!=====================================================================
!
! master routine for computing statistics of operations between sectors
! 
! subroutine master_op_stat
!
! master subroutine for tracking statistics 
!        on operations between sectors
!
! started 8/2011 by CWJ @ SDSU
!
!  CALLED BY
!     main routine  in bigstick_main.f90
!     reset_jumps4obs in blanczoslib2.f90
!
!  SUBROUTINES CALLED:  (all in this file)
!     count_XXXop_stats
!     count_XXYop_stats
!     count_XXop_stats
!     count_PNop_stats
!     count_SPEop_stats
!     count_total_ops
!     getfragmentstats
!
subroutine master_op_stat

  use operation_stats
  use sectors
  use flags3body
  implicit none

  integer :: is,fs
  integer(8) :: nops
  integer :: aerr

  if(.not.allocated(opstat)) then
     allocate ( opstat( nsectors(1), nsectors(1) ), stat=aerr  )
     if(aerr /= 0) call memerror("master_op_stat")
  end if
  opstat(:,:)%nopPP = 0
  opstat(:,:)%nopNN = 0
  opstat(:,:)%nopPN = 0
  opstat(:,:)%nopPPP= 0
  opstat(:,:)%nopNNN= 0
  opstat(:,:)%nopPPN= 0
  opstat(:,:)%nopPNN= 0
  opstat(:,:)%nopSPE = 0

  do is = 1,nsectors(1)

     do fs = 1,is    ! this ordering is mandatory in bigstick
        if(threebody)then
        call count_XXXop_stats(1, is,fs,nops)
        opstat(is,fs)%nopPPP = nops
        opstat(fs,is)%nopPPP = nops
        call count_XXXop_stats(2, is,fs,nops)
        opstat(is,fs)%nopNNN = nops
        opstat(fs,is)%nopNNN = nops

        call count_XXYop_stats(1,is,fs,nops)
        opstat(is,fs)%nopPPN = nops
        opstat(fs,is)%nopPPN = nops

        call count_XXYop_stats(2,is,fs,nops)
        opstat(is,fs)%nopPNN = nops
        opstat(fs,is)%nopPNN = nops

         else

         call count_XXop_stats(1,is,fs,nops)
         opstat(is,fs)%nopPP = nops
         opstat(fs,is)%nopPP = nops
         call count_XXop_stats(2,is,fs,nops)
         opstat(is,fs)%nopNN =  nops
         opstat(fs,is)%nopNN =  nops

         call count_PNop_stats( is,fs,nops)
         opstat(is,fs)%nopPN =  nops
         opstat(fs,is)%nopPN =  nops
 
         endif
         opstat(is,fs)%nopSPE = 0
         opstat(fs,is)%nopSPE = 0
      end do
      call count_SPEop_stats(is,nops)
      opstat(is,is)%nopSPE = nops
   end do
   call count_total_ops
   call getfragmentstats

   return

end subroutine master_op_stat
!===============================================================
!
! counts SPE operations in a sector; always diagonal
!
! CALLED By: master_op_stat
!
! INPUT
!
subroutine count_SPEop_stats(is,nops)
   use sectors
   use jumpNbody
   use system_parameters

   implicit none
   integer is
   integer(8) :: nops
   integer :: xsj  ! sector jump
   integer :: ncs
   integer csj    ! conjugate sector(jumps)
   integer cs
 
   xsj = is !subsectorlist(is)%parentsector
   cs  = xsd(1)%sector(xsj)%csector(1)
   ncs = xsd(1)%sector(xsj)%ncsectors
   csj  = xsd(1)%sector(xsj)%csector(ncs) 
   nops = int(xsd(2)%sector(cs)%xsdend-xsd(2)%sector(csj)%xsdstart+1,8)* & 
                   int(xsd(1)%sector(is)%xsdend-xsd(1)%sector(is)%xsdstart+1,8)

   return
end subroutine count_SPEop_stats

!===============================================================
!
! subroutine count_XXop_stats
!
! counts operations between two sectors
!
! INPUT:
!   is = initial (proton) sector NOTE: always have is >= fs
!   fs = final (proton) sector
!  itx = species (PP or NN )
!
! OUTPUT:
!   nops  = # of noperations
!
! CALLED BY: master_op_stat
!
subroutine count_XXop_stats(itx, is,fs,nops)

  use sectors
  use jumpNbody
  use system_parameters

  implicit none
  integer :: xsj  ! sector jump
  integer is,fs
  integer itx, ity
  integer(8) :: nops
  integer(8) ncsd    ! # of conjugate SDs
  integer csj    ! conjugate sector(jumps)
  integer ic
  integer cs
  integer(8) njumps

  ity = 3 -itx

  nops = 0
  if(np(itx) < 2)return

  if(itx == 2 .and. is /= fs)return   ! NN can't change proton sector

!..... LOOP OVER SECTOR JUMPS......... 

  do xsj = 1,x2bjump(itx)%nsectjumps
     if ( itx == 1)then   ! pp -sector jumps
         if(    is == x2bjump(itx)%isector(xsj) .and.  & 
                fs == x2bjump(itx)%fsector(xsj) ) then
            njumps = x2bjump(itx)%sjmp(xsj)%njumps
            ncsd = 0

!.............COUNT UP CONJUGATE SPECTATOR SDs

            do ic = 1,x2bjump(itx)%csjmp(xsj)%ncjmps  ! loop over conjugate sectors
                 cs = x2bjump(itx)%csjmp(xsj)%cjump(ic)
                 ncsd = ncsd + xsd(ity)%sector(cs)%nxsd
            end do
            nops = nops + njumps*ncsd
         endif

     else
         do ic = 1,x2bjump(itx)%csjmp(xsj)%ncjmps
            cs = x2bjump(itx)%csjmp(xsj)%cjump(ic)
            njumps = x2bjump(itx)%sjmp(xsj)%njumps

            if(cs == is .and. cs == fs) then
              ncsd =xsd(ity)%sector(cs)%nxsd
              nops = nops + njumps*ncsd
            end if

         end do

     end if

  end do ! xsj

  return

end subroutine count_XXop_stats

!===============================================================
!
! subroutine count_PNop_stats
!
! counts PN operations between two sectors
!
! INPUT:
!   is = initial (proton) sector NOTE: always have is >= fs
!   fs = final (proton) sector
!
! OUTPUT:

!   nops  = # of noperations
!
subroutine count_PNop_stats( is,fs,nops)

  use sectors
  use jumpNbody
  use system_parameters

  implicit none
  integer :: psj,nsj  ! proton, neutron sector jump
  integer is,fs,is0,fs0
  integer(8) :: nops
  integer ic
  integer cs
  integer(8) npjumps,nnjumps  ! # of proton, neutron jumps


  nops = 0
!........ NOTE ASSUME is >= fs ; correct for this
  if(is < fs)then
	  is0 = fs
	  fs0 = is
  else
	  is0 = is
	  fs0 = fs
  end if 
  if(np(1)*np(2) < 1)return

!..... LOOP OVER PROTON SECTOR JUMPS......... 

  do psj = 1,x1bjump(1)%nsectjumps
         if(    is0 == x1bjump(1)%isector(psj) .and.  & 
                fs0 == x1bjump(1)%fsector(psj) ) then
            npjumps = x1bjump(1)%sjmp(psj)%njumps

!.............. LOOP OVER CONJUGATE NEUTRON SECTOR JUMPS....

!.............COUNT UP CONJUGATE JUMPS 
            nnjumps = 0
            do ic = 1,x1bjump(1)%csjmp(psj)%ncjmps  ! loop over conjugate sectors
                 nsj = x1bjump(1)%csjmp(psj)%cjump(ic)
                 nnjumps = nnjumps + x1bjump(2)%sjmp(nsj)%njumps
            end do
            nops = nops + npjumps*nnjumps
         endif
  end do ! psj

  return

end subroutine count_PNop_stats

!===============================================================
!
! subroutine count_XXXop_stats
!
! counts operations between two sectors
!
! INPUT:
!   is = initial (proton) sector NOTE: always have is >= fs
!   fs = final (proton) sector
!  itx = species (PPP or NNN )
!
! OUTPUT:

!   nops  = # of noperations
!
subroutine count_XXXop_stats(itx, is,fs,nops)

  use sectors
  use jump3body
  use system_parameters

  implicit none
  integer :: xsj  ! sector jump
  integer is,fs
  integer itx, ity
  integer(8) :: nops
  integer(8) ncsd    ! # of conjugate SDs
  integer csj    ! conjugate sector(jumps)
  integer ic
  integer cs
  integer(8) njumps

  ity = 3 -itx

  nops = 0
  if(np(itx) < 3)return


  if(itx == 2 .and. is /= fs)return   ! NNN can't change proton sector

!..... LOOP OVER SECTOR JUMPS......... 

  do xsj = 1,x3bjump(itx)%nsectjumps
     if ( itx == 1)then   ! ppp -sector jumps
         if(    is == x3bjump(itx)%isector(xsj) .and.  & 
                fs == x3bjump(itx)%fsector(xsj) ) then
            njumps = x3bjump(itx)%sjmp(xsj)%njumps
            ncsd = 0

!.............COUNT UP CONJUGATE SPECTATOR SDs

            do ic = 1,x3bjump(itx)%csjmp(xsj)%ncjmps  ! loop over conjugate sectors
                 cs = x3bjump(itx)%csjmp(xsj)%cjump(ic)
                 ncsd = ncsd + xsd(ity)%sector(cs)%nxsd
            end do
            nops = nops + njumps*ncsd
         endif

     else
         do ic = 1,x3bjump(itx)%csjmp(xsj)%ncjmps
            cs = x3bjump(itx)%csjmp(xsj)%cjump(ic)
            njumps = x3bjump(itx)%sjmp(xsj)%njumps

            if(cs == is .and. cs == fs) then
              ncsd =xsd(ity)%sector(cs)%nxsd
              nops = nops + njumps*ncsd

            end if

         end do

     end if

  end do ! xsj

  return

end subroutine count_XXXop_stats
!===============================================================
!
! subroutine count_XXYop_stats
!
! counts PPN or PNN operations between two sectors
!
! INPUT:
!   itx = species of X
!   is = initial (proton) sector NOTE: always have is >= fs
!   fs = final (proton) sector
!
! OUTPUT:

!   nops  = # of noperations
!
subroutine count_XXYop_stats(itx, is,fs,nops)

  use sectors
  use jumpNbody
  use jumpdef
  use system_parameters

  implicit none
  integer itx,ity
  integer :: xsj,ysj  ! proton, neutron sector jump
  integer is,fs
  integer(8) :: nops
  integer ic
  integer cs
  integer(8) nxjumps,nyjumps  ! # of proton, neutron jumps

  type (jumpsect), pointer :: pNbjump, nNbjump

  ity = 3-itx
  nops = 0
  if((np(itx)-1)*np(ity) < 1)return

  if(itx ==1)then
      pNbjump => x2bjump(1)
      nNbjump => x1bjump(2)

  else
      pNbjump => x1bjump(1)
      nNbjump => x2bjump(2)

  end if

!..... LOOP OVER X SECTOR JUMPS......... 

  do xsj = 1,pNbjump%nsectjumps

         if(   ( is == pNbjump%isector(xsj) .and.  & 
                fs == pNbjump%fsector(xsj)) .or.  & 
                 ( fs == pNbjump%isector(xsj) .and.  & 
                is == pNbjump%fsector(xsj))  ) then
            nxjumps = pNbjump%sjmp(xsj)%njumps

!.............. LOOP OVER CONJUGATE Y SECTOR JUMPS....

!.............COUNT UP CONJUGATE JUMPS 
            nyjumps = 0
            do ic = 1,pNbjump%csjmp(xsj)%ncjmps  ! loop over conjugate sectors
                 ysj = pNbjump%csjmp(xsj)%cjump(ic)
                 nyjumps = nyjumps + nNbjump%sjmp(ysj)%njumps
            end do
            nops = nops + nxjumps*nyjumps
         endif


  end do ! xsj

  return

end subroutine count_XXYop_stats
!===============================================================
!
!  CALLED BY: master_op_stats
! 
subroutine count_total_ops

   use sectors
   use operation_stats
   implicit none
   integer is,fs
   
   allopPP = 0
   allopPN = 0
   allopNN = 0
   allopPPP = 0
   allopPPN = 0
   allopPNN = 0
   allopNNN = 0

   do is = 1,nsectors(1)
      do fs = 1,is
         allopPP = allopPP + opstat(is,fs)%nopPP
         allopPN = allopPN + opstat(is,fs)%nopPN
         allopNN = allopNN + opstat(is,fs)%nopNN
         allopPPP = allopPPP + opstat(is,fs)%nopPPP
         allopPPN = allopPPN + opstat(is,fs)%nopPPN
         allopPNN = allopPNN + opstat(is,fs)%nopPNN
         allopNNN = allopNNN + opstat(is,fs)%nopNNN

      end do ! fs

   end do  ! is

!   allops = allopPP + allopPN + allopNN
!   write(6,*)' total operations: ',allops
!   write(6,'(e12.6)' )dfloat(allops)
   return
end subroutine count_total_ops

!===============================================================
!
! routine to gather opstat data into opfragstat
!
! CALLED BY:
!   master_op_stat
!
subroutine getfragmentstats
  use fragments
  use operation_stats
  implicit none
  integer :: ifrag,jfrag
  integer(8) :: is, js
  integer :: aerr

  if(allocated(opfragstat)) deallocate(opfragstat)
  allocate( opfragstat(nfragments, nfragments), stat=aerr)
  if(aerr /= 0) call memerror("getfragmentstats 1")
  do ifrag = 1,nfragments
     do jfrag = 1,nfragments
        opfragstat(ifrag,jfrag)%nopSPE = 0
        opfragstat(ifrag,jfrag)%nopPP = 0
        opfragstat(ifrag,jfrag)%nopPN = 0
        opfragstat(ifrag,jfrag)%nopNN = 0
        opfragstat(ifrag,jfrag)%nopPPP = 0
        opfragstat(ifrag,jfrag)%nopPPN = 0
        opfragstat(ifrag,jfrag)%nopPNN = 0
        opfragstat(ifrag,jfrag)%nopNNN = 0

        do is = fragmentlist(ifrag)%ssectorstart,fragmentlist(ifrag)%ssectorend
           do js = fragmentlist(jfrag)%ssectorstart,fragmentlist(jfrag)%ssectorend
                 opfragstat(ifrag,jfrag)%nopSPE = opfragstat(ifrag,jfrag)%nopSPE+ opstat(is,js)%nopSPE
                 opfragstat(ifrag,jfrag)%nopPP = opfragstat(ifrag,jfrag)%nopPP+ opstat(is,js)%nopPP
                 opfragstat(ifrag,jfrag)%nopPN = opfragstat(ifrag,jfrag)%nopPN+ opstat(is,js)%nopPN
                 opfragstat(ifrag,jfrag)%nopNN = opfragstat(ifrag,jfrag)%nopNN+ opstat(is,js)%nopNN
                 opfragstat(ifrag,jfrag)%nopPPP = opfragstat(ifrag,jfrag)%nopPPP+ opstat(is,js)%nopPPP
                 opfragstat(ifrag,jfrag)%nopPPN = opfragstat(ifrag,jfrag)%nopPPN+ opstat(is,js)%nopPPN
                 opfragstat(ifrag,jfrag)%nopPNN = opfragstat(ifrag,jfrag)%nopPNN+ opstat(is,js)%nopPNN
                 opfragstat(ifrag,jfrag)%nopNNN = opfragstat(ifrag,jfrag)%nopNNN+ opstat(is,js)%nopNNN
!.................... NEXT WEIGHT TO GET COMBINED OPERATIONS.......................
                 opfragstat(ifrag,jfrag)%optot = opfragstat(ifrag,jfrag)%nopPP*opwtPP  & 
                                               + opfragstat(ifrag,jfrag)%nopPN*opwtPN  & 
                                               + opfragstat(ifrag,jfrag)%nopNN*opwtNN  & 
                                               + opfragstat(ifrag,jfrag)%nopPPP*opwtPPP  & 
                                               + opfragstat(ifrag,jfrag)%nopPPN*opwtPPN  & 
                                               + opfragstat(ifrag,jfrag)%nopPNN*opwtPNN  & 
                                               + opfragstat(ifrag,jfrag)%nopNNN*opwtNNN  &
											   + opfragstat(ifrag,jfrag)%nopSPE*opwtSPE

           end do ! js
        end do  ! is
     end do  ! jf
  end do  ! if

  return
end subroutine getfragmentstats
!==========================================================================================

