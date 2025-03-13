module phases_set_form_mod
!
! Set/calculate average atomic form factor
!
contains
!
!*******************************************************************************
!
subroutine phases_set_form(iscat, npkt, pha_curr, PHA_MAXPTS, PHA_MAXSCAT, PHA_MAXPHA, &
   pha_form, MAX_ISTL, powder_istl, DIF_MAXSCAT, cfact_pure,          &
   diff_table, RAD_DISC, rten)
!
use discus_config_mod
use discamb_mod
!
use lib_metric_mod
use precision_mod
use spline_mod
!
implicit none
!
integer, intent(in) :: iscat
integer, intent(in) :: npkt
integer, intent(in) :: pha_curr
integer, intent(in) :: PHA_MAXPTS
integer, intent(in) :: PHA_MAXSCAT
integer, intent(in) :: PHA_MAXPHA
real(kind=PREC_DP), dimension(0:PHA_MAXPTS, 0:PHA_MAXSCAT, 1:PHA_MAXPHA), intent(inout) :: pha_form
integer, intent(in) :: MAX_ISTL
integer           , dimension(0:MAX_ISTL), intent(in) :: powder_istl
!integer, intent(in) :: CFPKT
integer, intent(in) :: DIF_MAXSCAT
complex(kind=PREC_DP), dimension(0:CFPKT, 1:DIF_MAXSCAT), intent(in) :: cfact_pure
integer, intent(in) :: diff_table
integer, intent(in) :: RAD_DISC   
real(kind=PREC_DP), dimension(3,3), intent(in) :: rten
!
character(len=48) :: ofile
integer :: hl, kl, ll, hh, kh, lh ! Limits on hkl in discamb table
integer :: h, k, l
integer :: i
integer :: ientry
integer :: ipkt              ! Points in initial unevenly spaced DISCAMB form factors
integer :: npkt_equi         ! Points in           evenly spaced DISCAMB form factors
integer :: ilow, ihigh       ! Min max entries in stl table
real(kind=PREC_DP), dimension(3) :: hkl
real(kind=PREC_DP) :: stl
real(kind=PREC_DP) :: signum
real(kind=PREC_DP) :: fract
real(kind=PREC_DP) :: formf
real(kind=PREC_DP) :: xmin, xmax, xstep    ! min, max, step for evenly palces sine(theta)/lambda
real(kind=PREC_DP), dimension(:), allocatable :: xpl_in
real(kind=PREC_DP), dimension(:), allocatable :: ypl_in
real(kind=PREC_DP), dimension(:), allocatable :: xequi
real(kind=PREC_DP), dimension(:), allocatable :: yequi
real(kind=PREC_DP), dimension(:), allocatable :: form_discamb
real(kind=PREC_DP), dimension(:), allocatable :: table ! Table of form factors from DISCAMB
real(kind=PREC_DP), dimension(:), allocatable :: weight! Table of form factors from DISCAMB
!
if(diff_table==RAD_DISC) then    ! Average form factors from Discamb
   allocate(table(0:CFPKT))
   allocate(weight(0:CFPKT))
   table  = 0.0_PREC_DP
   weight = 0.0_PREC_DP
   hl = lbound(disc_list(iscat)%four_form_tsc,1)
   kl = lbound(disc_list(iscat)%four_form_tsc,2)
   ll = lbound(disc_list(iscat)%four_form_tsc,3)
   hh = ubound(disc_list(iscat)%four_form_tsc,1)
   kh = ubound(disc_list(iscat)%four_form_tsc,2)
   lh = ubound(disc_list(iscat)%four_form_tsc,3)
!write(*,*) ' DISCAMB LIMITS ', lbound(disc_list(iscat)%four_form_tsc), ubound(disc_list(iscat)%four_form_tsc)
   ilow  = CFPKT + 1
   ihigh =-1
hkl(1) = hh
hkl(2) = kh
hkl(3) = lh
!write(*,*) ' HKL ', hkl, 0.5_PREC_DP*lib_blen(rten, hkl)
   do l=ll, lh
      do k=kl, kh
         do h=hl, hh
            if(real(disc_list(iscat)%four_form_tsc(h,k,l,1))>0.0) then
               hkl(1) = h
               hkl(2) = k
               hkl(3) = l
               stl = 0.5_PREC_DP*lib_blen(rten, hkl)
               ientry = int(stl/CFINC)
               ilow  = min(ilow, ientry)
               ihigh = max(ihigh, ientry)
               fract  = stl - real(ientry, kind=PREC_DP)*CFINC
               formf  = sqrt(dble(      disc_list(iscat)%four_form_tsc(h,k,l,1)*   &
                                  conjg(disc_list(iscat)%four_form_tsc(h,k,l,1))))
               table(ientry  )  = table(ientry  )  + (1.0_PREC_DP-fract)*formf
               table(ientry+1)  = table(ientry+1)  + (            fract)*formf
               weight(ientry  ) = weight(ientry  ) +  1.0_PREC_DP-fract
               weight(ientry+1) = weight(ientry+1) +              fract
            endif
         enddo
      enddo
   enddo
!write(*,*) ' DID accumulation ', ilow, ihigh, maxval(table), maxval(weight)
   allocate(xpl_in(1   :ihigh))
   allocate(ypl_in(1   :ihigh))
   xpl_in = 0.0_PREC_DP
   ypl_in = 0.0_PREC_DP
   ipkt = 0
   do i=ilow,ihigh
      if(table(i)>0.0_PREC_DP) then
         ipkt = ipkt + 1
         xpl_in(ipkt) = real(i*CFINC,kind=PREC_DP)
         ypl_in(ipkt) = table(i)/weight(i)
      endif
   enddo
!write(*,*) ' DID xpl ypl       ', maxval(xpl_in), maxval(ypl_in)
!write(ofile,'(a,i3.3)') 'discamb_xpl.',iscat
!open(54,file=ofile, status='unknown')
!do i=1, ipkt
!  write(54, '(2(2x,f10.6))') xpl_in(i), ypl_in(i)
!enddo
!close(54)
   deallocate(table)
   deallocate(weight)
   xmin  = real(ilow *CFINC,kind=PREC_DP)
   xmax  = real(ihigh*CFINC,kind=PREC_DP)
   xstep = CFINC
   npkt_equi = nint((xmax-xmin)/xstep) + 1
   allocate(xequi(1: npkt_equi))
   allocate(yequi(1: npkt_equi))
!write(*,*) ' PRIOR TO SPLINE ', ipkt, npkt_equi
   ilow = 1
   call spline_prep(ilow,ipkt, xpl_in, ypl_in, xmin, xmax, xstep, npkt_equi, xequi, yequi)
!write(*,*) ' AFTER    SPLINE ', maxval(xequi), maxval(yequi)
!  call spline_prep(nlow,npkt, xpl_in, ypl_in, xmin, xmax, xstep, npkt_equi, xequi, yequi)
   deallocate(xpl_in)
   deallocate(ypl_in)
   allocate(form_discamb(0:CFPKT))
   form_discamb = 0.0_PREC_DP
   form_discamb(1:ilow-1) = yequi(ilow)
   form_discamb(ilow: npkt_equi) = yequi(ilow: npkt_equi)
!write(ofile,'(a,i3.3)') 'discamb_form.',iscat
!open(54,file=ofile, status='unknown')
!do i=0, npkt_equi
!  write(54, '(2(2x,f10.6))') i*CFINC, form_discamb(i)
!enddo
!close(54)
!write(*,*) ' WRITE ISCAT ', iscat, maxval(form_discamb)
   deallocate(xequi)
   deallocate(yequi)
   do k=0, npkt
      pha_form(k, iscat, pha_curr) = form_discamb(powder_istl(k))
   enddo
   deallocate(form_discamb)
else
   do k=0, npkt
      signum = 1.0D0
      IF(real(cfact_pure(1, iscat),kind=PREC_DP)< 0.0D0) signum = -1.0D0
      pha_form(k, iscat, pha_curr) = sqrt( dble (       cfact_pure (powder_istl (k), iscat)   *   &
                                                 conjg (cfact_pure (powder_istl (k), iscat) )  )  &
                                         )                                                        &
                                     * signum
   enddo
endif
!
end subroutine phases_set_form
!
end module phases_set_form_mod
