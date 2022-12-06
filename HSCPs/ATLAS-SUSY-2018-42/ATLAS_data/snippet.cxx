// Assumes some basic ROOT classes usable
#include <TRandom3.h>
#include <TLorentzVector.h>

// Skeleton classes -- precise definition not given within this snippet.
class TruthEvent; // collection of truth-level particles of the event
class Particle;   // contains 4-vector, PDG particle ID, charge, lifetime, production and decay vertices


//----------------------------------------------------------------------------------------------------
bool acceptParticle( Particle& particle ) {

    // the particle should be a BSM particle
    if( ! particle.isBSMLLP() )                    return false;
    
    // the particle should be charged ( |q| = 1 needs to be assumed so that Bethe-Bloch relation to work )
    if( 1 != std::abs( particle.charge() ) )       return false;
    
    // the particle's acceptance
    if( particle.Pt() < 120 * Unit::GeV )          return false;
    if( std::abs( particle.Eta() ) > 1.8 )         return false;
    if( particle.decay.Perp() < 500. * Unit::mm )  return false;

    return true;
}



//----------------------------------------------------------------------------------------------------
bool isDecayBeforeMuonSpectrometer( Particle& particle ) {
    
    //check radii and z bounding of END of tile calo
    // |z| < 6.0 m, R < 3.9m
    
    if ( std::abs(particle.Z() ) < 6000. * Unit::mm &&
         particle.Perp() < 3900. * Unit::mm ) {
        return true;
    }
    
    return false;
}


//----------------------------------------------------------------------------------------------------
bool isCountedForCaloMet( Particle& particle ) {
    
    bool caloMetFlag { false };
          
    if( particle.isPrompt( particle ) ) {
              
        if( particle.isMuon() || 
            particle.isNeutrino() ||
            particle.isNeutralStableBSM() ||
            particle.isChargedBSMLLP() && isDecayBeforeMuonSpectrometer( particle ) )
        {
            caloMetFlag = true;
        }
    }

    return caloMetFlag;
}



//----------------------------------------------------------------------------------------------------
std::map<std::string, double> processEvent( TruthEvent *event, double target_mass ) {

    auto& particles = event->particles();

    std::vector<Particle> muons;
    std::vector<Particle> signalParticles;
    std::vector<Particle> truthCaloMetParticles;

    // Extract SM muons
    std::copy_if( particles.begin(), particles.end(), std::back_inserter(muons), isMuon );

    // Extract BSM particles satisfying the acceptance
    std::copy_if( particles.begin(), particles.end(), std::back_inserter(signalParticles), acceptParticle );


    
    //===================================================================
    // Calculate truthCaloMET
    
    std::copy_if( particles.begin(), particles.end(), std::back_inserter(truthCaloMetParticles), isCountedForCaloMet );
    
    auto truthCaloMissingFourVec = std::accumulate( truthCaloMetParticles.begin(), truthCaloMetParticles.end(), TLorentzVector(),
                                                [](auto& p1, auto& p2){ return p1.p4() - p2.p4(); } );

    auto truthCaloMet = truthCaloMissingFourVec.Perp();


    //===================================================================
    // Calculate truthMET
    
    auto truthMissingFourVec = truthCaloMetFourVec;

    // 1) subtract muons momentum
    for( auto& muon : muons ) {
        truthMissingFourVec -= muon.p4();
    }

    // 2) stochastic subtraction of BSM charged LLP momentum if it is reconstructed as a muon
    for( auto& signal : signalParticles ) {
        auto expDecayLength = truth.betaGamma() * spped_of_light * ( signal.lifetime() * Unit::ns); // Unit: [m]

        // detector-stable -- R = 10 m
        auto stableLength   = (10.0 * Unit::m) * cosh( signal.Eta() );                              // Unit: [m]
        auto prob_stable    = exp( -expDecayLength / stableLength );
        
        // Estimate the reconstruction efficiency from the HpData record
        auto recoEff = muonRecoEffMap->Interpolate( truth.bEta(), std::abs( truth.Eta() ) );
        
        if( gRandom->Uniform() < prob_stable * recoEff ) {
            
            truthMissingFourVec -= signal.p4();
            
        }
    }
    
    auto truthOfflineMet = truthMissingFourVec.Perp();

    
    //===================================================================
    // Retrieve trigger efficiency using the Figure 3a turn-on curve
    auto triggerPassEfficiency = getTriggerEfficiency( truthCaloMet );
    
    // Retrieve trigger efficiency using the Figure 3b turn-on curve
    auto eventEfficiency       = getEventEfficiency( truthOfflineMet );

    
    // Loop over signal particles
    std::vector<double> passWeight_High;
    std::vector<double> passWeight_Low;

    // mass window efficiency

    
    for( auto& particle : signalParticles ) {
        
        auto betaGamma = particle.beta() * particle.gamma();
        
        // Mutually exclusive categorization
        auto selectionEfficiency_High = getTrackSelectionEfficiency_High( betaGamma );
        auto selectionEfficiency_Low  = getTrackSelectionEfficiency_Low( betaGamma );

        // Here we implicitly assume the target_mass is identical to the signal particle's mass
        auto massWindowEff_High = 0.75 - 0.052 * ( target_mass / Unit::TeV );
        auto massWindowEff_Low  = 0.60;
        
        passWeight_High.emplace_back( selectionEfficiency_High * massWindowEff_High );
        passWeight_Low .emplace_back( selectionEfficiency_Low  * massWindowEff_Low  );
        
    }

    
    auto accumFunc = []( double init, double eff ) { return init * (1.0 - eff); };
    
    auto accum_effHigh  = 1.0 - std::accumulate( passWeight_High.begin(),  passWeight_High.end(),  1.0, accumFunc );
    auto accum_effLow   = 1.0 - std::accumulate( passWeight_Low.begin(),   passWeight_Low.end(),   1.0, accumFunc );

    std::map<std::string, double> weights {
        "SR-Inclusive_High" : event->weight() * triggerPassEfficiency * eventEfficiency * accum_effHigh,
        "SR-Inclusive_Low"  : event->weight() * triggerPassEfficiency * eventEfficiency * accum_effLow
    };

    return weights;
}


