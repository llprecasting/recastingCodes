C ************************************************************
C Source for the library implementing a bias function that 
C populates the large pt tale of the leading jet. 
C
C The two options of this subroutine, that can be set in
C the run card are:
C    > (double precision) pT_bias_target : target chargino pt value
C    > (double precision) pT_bias_enhancement_power : exponent
C    > (double precision) pT_bias_min : minimum pT value
C    > (double precision) pT_bias_enhancement_power : exponent
C
C Schematically, the functional form of the enhancement is
C    bias_wgt = [max pT(evt)/pT_bias_target]^enhancement_power
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the 
C relevant parameters for this bias module.
C
C  1000 = pT_bias_target
C  4 = pT_bias_enhancement_power
C  50 = pT_bias_min

      subroutine bias_wgt_custom(p, original_weight, bias_weight)
      implicit none
C
C Parameters
C
          include '../../maxparticles.inc'          
          include 'nexternal.inc'
          include 'run.inc' ! include defition from the run_card (via common-block).

C
C Arguments
C
          double precision p(0:3,nexternal)
          double precision pTmax
          double precision original_weight, bias_weight
C
C local variables
C
          integer i
c
c local variables defined in the run_card
c
c          double precision pT_bias_target
c          double precision pT_bias_enhancement_power
C
C Global variables
C
C
C Mandatory common block to be defined in bias modules
C
          double precision stored_bias_weight
          data stored_bias_weight/1.0d0/          
          logical impact_xsec, requires_full_event_info
C         We only want to bias distributions, but not impact the xsec. 
          data impact_xsec/.False./
C         Of course this module does not require the full event
C         information (color, resonances, helicities, etc..)
          data requires_full_event_info/.False./ 
c          common/bias/stored_bias_weight,impact_xsec,
c     &                requires_full_event_info
          logical is_a_j(nexternal),is_a_l(nexternal),
     &            is_a_b(nexternal),is_a_a(nexternal),
     &            is_a_onium(nexternal),is_a_nu(nexternal),
     &            is_heavy(nexternal),do_cuts(nexternal)
          common/to_specisa/is_a_j,is_a_a,is_a_l,is_a_b,is_a_nu,
     &                      is_heavy,is_a_onium,do_cuts


C
C    Setup the value of the parameters from the run_card    
C
c      include '../bias.inc'

C --------------------
C BEGIN IMPLEMENTATION
C --------------------
          
          bias_weight = 1.0d0
          pTmax = pT_bias_min
          do i=1,nexternal
            if (is_heavy(i)) then
              pTmax = max(pTmax,dsqrt(p(1,i)**2 + p(2,i)**2))
            endif
          enddo
 
          if (pTmax.gt.0d0) then
            bias_weight = (pTmax/pT_bias_target)**pT_bias_enhancement_power
          endif

       return

      end subroutine bias_wgt_custom
