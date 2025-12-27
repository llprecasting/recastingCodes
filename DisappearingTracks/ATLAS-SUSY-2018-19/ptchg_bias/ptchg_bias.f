C ************************************************************
C Source for the library implementing a bias function that 
C populates the large pt tale of the leading jet. 
C
C The two options of this subroutine, that can be set in
C the run card are:
C    > (double precision) pT_bias_target : target chargino pt value
C    > (double precision) pT_bias_enhancement_power : exponent
C
C Schematically, the functional form of the enhancement is
C    bias_wgt = [average pT(evt)/pT_bias_target]^enhancement_power
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the 
C relevant parameters for this bias module.
C
C  parameters = {'pT_bias_target': 1000.0,
C               'pT_bias_enhancement_power': 4.0}
C

      subroutine bias_wgt_custom(p, original_weight, bias_weight)
      implicit none
C
C Parameters
C
          include '../../maxparticles.inc'          
          include '../../nexternal.inc'
          include 'run.inc' ! include defition from the run_card (via common-block).

C
C Arguments
C
          double precision p(0:3,nexternal)
          double precision pTlist(0:1)
          double precision pTavg
          double precision original_weight, bias_weight
C
C local variables
C
          integer i,j
          integer ipdg(nexternal)
          double precision pt(nexternal)
c
c local variables defined in the run_card
c
c          double precision pT_bias_target_chg
c          double precision pT_bias_enhancement_power_chg
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


C
C    Setup the value of the parameters from the run_card    
C
c      include '../bias.inc'

C --------------------
C BEGIN IMPLEMENTATION
C --------------------
          
          bias_weight = 1.0d0
          pTlist = (/ 0d0, 0d0 /)
          do i=1,nexternal
            if (ipdg(i).eq.1000024) then
              pTlist(0) = dsqrt(p(1,i)**2 + p(2,i)**2)
            else if (ipdg(i).eq.-1000024) then
              pTlist(1) = dsqrt(p(1,i)**2 + p(2,i)**2)
            endif
          enddo

          pTavg = (pTlist(0)+pTlist(1))/2.0

          if (pTavg.gt.0.0d0) then
            bias_weight = (pTavg/pT_bias_target_chg)**pT_bias_enhancement_power_chg
          endif

       return

      end subroutine bias_wgt_custom
