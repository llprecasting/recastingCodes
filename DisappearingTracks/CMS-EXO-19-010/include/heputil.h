
// Following code contains some code from Event.h in HEPUtils:
// -*- C++ -*-
// 
// This file is part of HEPUtils -- https://bitbucket.org/andybuckley/heputils
// Copyright (C) 2013-2018 Andy Buckley <andy.buckley@cern.ch>
//
// Embedding of HEPUtils code in other projects is permitted provided this
// notice is retained and the HEPUtils namespace and include path are changed.
//




#pragma once

#include "ExHEPUtils/BinnedFn.h"
#include "ExHEPUtils/Jet.h"
#include "ExHEPUtils/MathUtils.h"
#include "ExHEPUtils/Vectors.h"
#include "ExHEPUtils/Particle.h"
#include "pdg.h"

//#include "useYoda.h"






namespace HEP {

  // Need a function to calculate the ET from a vector (defined as E sin theta). Honestly I don't know why that is not there
  // can therefore define as E * pT/p = E*pT/sqrt(pT^2 + pz^2) = E/sqrt(1 + pz^2/pT^2)
  // in fastjet this is inline double Et() const {return (_kt2==0) ? 0.0 : _E/sqrt(1.0+_pz*_pz/_kt2);}
  // in HEPUtils we store the mass of the four-vector instead of the energy (go figure) so we use instead
  // pT * sqrt((1 + m^2/p^2)) provided that p^2 is not zero

  


  inline double ET2(P4 p4)
  {
    double p2 = p4.p2();
    if(p2 == 0)
    {
      return 0.0;
    };
       
    return p4.pT2()*(1.0+p4.m2()/p2);
  };

  inline double ET(P4 p4)
  {
    return(sqrt(ET2(p4)));
  };

  inline bool compare_pT(const HEP::Particle* p1, const HEP::Particle* p2)
	{
		return (p1->pT() > p2->pT());  // to sort with highest pt first
	};

  inline bool vertex_match(const HEP::P4 &v1, const HEP::P4 &v2)
  {
    HEP::P4 vdiff = v1-v2;
    if ((v1-v2).p() > 1.0) return false; // vertices closer than a mm are combined
    if ((v1-v2).T() > 1E-3) return false;
    return true;

  }


  // glorified container for final-state particles + jet constituents, with a method to translate production/decay vertices
  class PileupEvent {
    private:
      std::vector<P4*> _jet_constituents;
      std::vector<Particle*> _particles;
      std::vector<P4*> _neutralhadrons;
      std::vector<Particle*> _chargedhadrons;
    public:
   
    void clear() {
      
      #define DELCLEAR(v) do { if (!v.empty()) for (size_t i = 0; i < v.size(); ++i) delete v[i]; v.clear(); } while (0)
      DELCLEAR(_jet_constituents);
      DELCLEAR(_particles);
      DELCLEAR(_neutralhadrons);
      DELCLEAR(_chargedhadrons);
      #undef DELCLEAR
    }
       /// Default constructor
    PileupEvent() { clear(); }

    /// Destructor (cleans up all passed Particles and calculated Jets)
    ~PileupEvent() {
      clear();
    }


      void add_neutral_hadron(double px, double py, double pz, double e) {
        P4* thad= new P4(px, py,pz, e);
        _neutralhadrons.push_back(thad);
    }

    void add_neutral_hadron(const P4 &p4)
    {
      P4* thad = new P4(p4);
      _neutralhadrons.push_back(thad);

    }


    void add_charged_hadron(double px, double py, double pz, double e, int pid) {
        Particle* thad= new Particle(px, py, pz, e, pid);
        _chargedhadrons.push_back(thad);
    }


    void add_charged_hadron(Particle* p)
    {
       
      _chargedhadrons.push_back(p);

    }
    /*
    void add_charged_hadron(const Particle &p4)
    {
       Particle* thad = new Particle(p4);
      _chargedhadrons.push_back(thad);

    }
    */

    void add_jet_constituent(double px, double py, double pz, double e) {
        P4* thad= new P4(px, py,pz, e);
        _jet_constituents.push_back(thad);
    }

    void add_jet_constituent(const P4 &p4)
    {
      P4* thad = new P4(p4);
      _jet_constituents.push_back(thad);

    }

    void add_particle(Particle* p) {
       _particles.push_back(p);

    }

    /// Get jet constituents
    const std::vector<P4*>& get_jet_constituents() const {
      return _jet_constituents;
    }

    /// Get charged hadrons
    const std::vector<Particle*>& get_charged_hadrons() const {
      return _chargedhadrons;
    }

    /// Get neutral hadrons
    const std::vector<P4*>& get_neutral_hadrons() const {
      return _neutralhadrons;
    }

    /// Get particles
    const std::vector<Particle*>& get_particles() const {
      return _particles;
    }

    std::vector<Particle*> translate_particles(double dz, double dt) {
      std::vector<Particle*> out_particles;
      for(auto part : _particles)
      {
        Particle* new_part = new Particle(part);
        new_part->translate_zt(dz,dt);
        out_particles.push_back(new_part);
      }
      return out_particles;
    }

    std::vector<Particle*> translate_charged_hadrons(double dz, double dt) {
      std::vector<Particle*> out_particles;
      for(auto part : _chargedhadrons)
      {
        Particle* new_part = new Particle(part);
        new_part->translate_zt(dz,dt);
        out_particles.push_back(new_part);
      }
      return out_particles;
    }



  };


  /// Simple event class, separating particles into classes
  // I've had to replace the HEPUtils event because it was insufficient



  class Event {

    private:

    /// @name Internal particle / vector containers
    //@{

    /// Event weights
    std::vector<double> _weights;

    /// @name Separate particle collections
    //@{
    std::vector<Particle*> _photons, _electrons, _muons, _taus, _invisibles, _HSCPs, _chargedhadrons;
    //@}

    /// Jets collection (mutable to allow sorting)
    mutable std::vector<Jet*> _jets;

    /// Missing momentum vector
    P4 _pmiss;

    std::vector<P4*> _neutralhadrons;
    //@}
    /// Hide copy assignment, since shallow copies of Particle & Jet pointers create ownership/deletion problems
    /// @todo Reinstate as a deep copy uing cloneTo?
    
    void operator = (const Event& e) {
      clear(); //< Delete current particles
      _weights = e._weights;
      _photons = e._photons;
      _electrons = e._electrons;
      _muons = e._muons;
      _taus = e._taus;
      _invisibles = e._invisibles;
      _HSCPs=e._HSCPs;
      _jets = e._jets;
      _pmiss = e._pmiss;
      _neutralhadrons = e._neutralhadrons;
      _chargedhadrons = e._chargedhadrons;
    }


  public:
    double Average_rho;
    /// Empty the event's particle, jet and MET collections
    void clear() {
      _weights.clear();
      // TODO: indexed loop -> for (Particle* p : particles()) delete p;
      #define DELCLEAR(v) do { if (!v.empty()) for (size_t i = 0; i < v.size(); ++i) delete v[i]; v.clear(); } while (0)
      DELCLEAR(_photons);
      DELCLEAR(_electrons);
      DELCLEAR(_muons);
      DELCLEAR(_taus);
      DELCLEAR(_invisibles);
      DELCLEAR(_HSCPs);
      DELCLEAR(_jets);
      DELCLEAR(_neutralhadrons);
      DELCLEAR(_chargedhadrons);
      #undef DELCLEAR

      _pmiss.clear();
    }
    /// Default constructor
    Event() { clear(); }

    /// Constructor from a list of Particles
    Event(const std::vector<Particle*>& ps, const std::vector<double>& weights=std::vector<double>()) {
      clear();
      _weights = weights;
      add_particles(ps);
    }

    /// Destructor (cleans up all passed Particles and calculated Jets)
    ~Event() {
      clear();
    }

    void add_invisible(Particle* p) {
        _invisibles.push_back(p);
    }

    void add_HSCP(Particle* p) {
        _HSCPs.push_back(p);
    }


    void add_neutral_hadron(double px, double py, double pz, double e) {
        P4* thad= new P4(px, py,pz, e);
        _neutralhadrons.push_back(thad);
    }
    
    void add_neutral_hadron(const P4 &p4)
    {
      P4* thad = new P4(p4);
      _neutralhadrons.push_back(thad);

    }
    void add_neutral_hadrons(const std::vector<P4*> &hadrs)
    {
      for (size_t i = 0; i < hadrs.size(); ++i) add_neutral_hadron(*hadrs[i]);
    }

    

    void add_charged_hadron(double px, double py, double pz, double e, int pid) {
        Particle* thad= new Particle(px, py, pz, e, pid);
        _chargedhadrons.push_back(thad);
    }

    void add_charged_hadron(const P4 &p4, int pid)
    {
       Particle* thad = new Particle(p4, pid);
      _chargedhadrons.push_back(thad);

    }

    void add_charged_hadron(Particle* p)
    {
       
      _chargedhadrons.push_back(p);

    }

    void add_charged_hadrons(const std::vector<Particle*> &hadrs)
    {
      for (size_t i = 0; i < hadrs.size(); ++i) add_charged_hadron(hadrs[i]);
    }
    /*
    void add_charged_hadron(const Particle &p4)
    {
       Particle* thad = new Particle(p4);
      _chargedhadrons.push_back(thad);

    }
    

    void add_charged_hadrons(const std::vector<Particle*> &hadrs)
    {
      for (size_t i = 0; i < hadrs.size(); ++i) add_charged_hadron(*hadrs[i]);
    }
    */
   
    const std::vector<Particle*>& get_charged_hadrons() const {
      return _chargedhadrons;
    }

    /// Get prompt electrons
    const std::vector<P4*>& get_neutral_hadrons() const {
      return _neutralhadrons;
    }

    /// Clone a copy on the heap
    Event* clone() const {
      Event* rtn = new Event();
      cloneTo(rtn);
      return rtn;
    }

    /// Clone a deep copy (new Particles and Jets allocated) into the provided event pointer
    void cloneTo(Event* e) const {
      assert(e != NULL);
      cloneTo(*e);
    }

    /// Clone a deep copy (new Particles and Jets allocated) into the provided event object
    void cloneTo(Event& e) const {
      e.set_weights(_weights);
      const std::vector<Particle*> ps = particles();
      for (size_t i = 0; i < ps.size(); ++i) {
        e.add_particle(new Particle(*ps[i]));
      }
      const std::vector<Jet*> js = jets();
      for (size_t i = 0; i < js.size(); ++i) {
        e.add_jet(new Jet(*js[i]));
      }
      e._pmiss = _pmiss;
    }

    //@}


    ///////////////////////


    /// Set the event weights (also possible directly via non-const reference)
    void set_weights(const std::vector<double>& ws) {
      _weights = ws;
    }

    /// Set the event weights to the single given weight
    void set_weight(double w) {
      _weights.clear();
      _weights.push_back(w);
    }

    /// Get the event weights (const)
    const std::vector<double>& weights() const {
      return _weights;
    }

    /// Get the event weights (non-const)
    std::vector<double>& weights() {
      return _weights;
    }

    /// Get a single event weight -- the nominal, by default
    double weight(size_t i=0) const {
      if (_weights.empty()) {
        if (i == 0) return 1;
        throw std::runtime_error("Trying to access non-default weight from empty weight vector");
      }
      return _weights[i];
    }


    /////////////////


    /// Add a particle to the event
    ///
    /// Supplied particle should be new'd, and Event will take ownership.
    ///
    /// @warning The event takes ownership of all supplied Particles -- even
    /// those it chooses not to add to its collections, which will be
    /// immediately deleted. Accordingly, the pointer passed by user code
    /// must be considered potentially invalid from the moment this function is called.
    ///
    /// @todo "Lock" at some point so that jet finding etc. only get done once
    void add_particle(Particle* p) {
      if (!p->is_prompt())
        {
          if(p->abspid() == 1000024)
          {
            _HSCPs.push_back(p);
          }
          else if(p->pid() == 1000022)
          {
            _invisibles.push_back(p);
          }
          else
          {
            delete p;
          }
          /*
          if( (p->abspid()) < 1000000 || (p->apspid() > 9000000))
          {
            delete p;
          }
          else if(p->abspid() == 1000024)
          {
            _HSCPs.push_back(p);
          }
          else if(p->pid() == 1000022)
          {

          }
          else
          {
            delete p;
          }
          */
          
        }
      else if (p->pid() == 22)
        _photons.push_back(p);
      else if (p->abspid() == 11)
        _electrons.push_back(p);
      else if (p->abspid() == 13)
        _muons.push_back(p);
      else if (p->abspid() == 15)
        _taus.push_back(p);
      else if (p->abspid() == 12 || p->abspid() == 14 || p->abspid() == 16 ||
               p->pid() == 1000022 || p->pid() == 1000039 ||
               in_range(p->pid(), 50, 60)) //< invert definition to specify all *visibles*?
        _invisibles.push_back(p);
      else if(p->abspid() == 1000024)
        _HSCPs.push_back(p);
      else
        delete p;
    }


    /// Add a collection of final state particles to the event
    ///
    /// Supplied particles should be new'd, and Event will take ownership.
    void add_particles(const std::vector<Particle*>& ps) {
      for (size_t i = 0; i < ps.size(); ++i) add_particle(ps[i]);
    }


    /// Get all final state particles
    /// @todo Note the return by value: it's not efficient yet!
    /// @note Overlap of taus and e/mu
    std::vector<Particle*> particles() const {
      // Add together all the vectors of the different particle types
      std::vector<Particle*> rtn;
      // rtn.reserve(visible_particles().size() + _invisibles.size());
      rtn.reserve(_photons.size() + _electrons.size() + _muons.size() + _taus.size() + _invisibles.size());
      #define APPEND_VEC(vec) rtn.insert(rtn.end(), vec.begin(), vec.end())
      // APPEND_VEC(visible_particles());
      APPEND_VEC(_photons);
      APPEND_VEC(_electrons);
      APPEND_VEC(_muons);
      APPEND_VEC(_taus);
      APPEND_VEC(_invisibles);
      #undef APPEND_VEC
      return rtn;
      /// @todo Or use Boost range join to iterate over separate containers transparently... I like this
      /// @todo Cache, or otherwise think about efficiency since this gets called by the destructor
    }


    /// Get visible state particles
    /// @todo Note the return by value: it's not efficient yet!
    /// @note Overlap of taus and e/mu
    std::vector<Particle*> visible_particles() const {
      // Add together all the vectors of the different particle types
      std::vector<Particle*> rtn;
      rtn.reserve(_photons.size() + _electrons.size() + _muons.size() + _taus.size());
      #define APPEND_VEC(vec) rtn.insert(rtn.end(), vec.begin(), vec.end() )
      APPEND_VEC(_photons);
      APPEND_VEC(_electrons);
      APPEND_VEC(_muons);
      APPEND_VEC(_taus);
      #undef APPEND_VEC
      return rtn;
      /// @todo Add together all the vectors of the different particle types
      /// @todo Or use Boost range join to iterate over separate containers transparently... I like this
    }


    /// Get invisible final state particles
    const std::vector<Particle*>& HSCPs() const {
      return _HSCPs;
    }
    /// Get invisible final state particles (non-const)
    std::vector<Particle*>&  HSCPs() {
      return _HSCPs;
    }


    /// Get invisible final state particles
    const std::vector<Particle*>& invisible_particles() const {
      return _invisibles;
    }
    /// Get invisible final state particles (non-const)
    std::vector<Particle*>& invisible_particles() {
      return _invisibles;
    }


    /// Get prompt electrons
    const std::vector<Particle*>& electrons() const {
      return _electrons;
    }
    /// Get prompt electrons (non-const)
    std::vector<Particle*>& electrons() {
      return _electrons;
    }


    /// Get prompt muons
    const std::vector<Particle*>& muons() const {
      return _muons;
    }
    /// Get prompt muons (non-const)
    std::vector<Particle*>& muons() {
      return _muons;
    }


    /// Get prompt (hadronic) taus
    const std::vector<Particle*>& taus() const {
      return _taus;
    }
    /// Get prompt (hadronic) taus (non-const)
    std::vector<Particle*>& taus() {
      return _taus;
    }


    /// Get prompt photons
    const std::vector<Particle*>& photons() const {
      return _photons;
    }
    /// Get prompt photons (non-const)
    std::vector<Particle*>& photons() {
      return _photons;
    }


    /// @name Jets
    /// @todo Why the new'ing? Can we use references more?
    //@{

    /// @brief Get anti-kT 0.4 jets (not including charged leptons or photons)
    const std::vector<Jet*>& jets() const {
      std::sort(_jets.begin(), _jets.end(), _cmpPtDesc);
      return _jets;
    }

    /// @brief Get anti-kT 0.4 jets (not including charged leptons or photons) (non-const)
    std::vector<Jet*>& jets()  {
      std::sort(_jets.begin(), _jets.end(), _cmpPtDesc);
      return _jets;
    }

    /// @brief Set the jets collection
    ///
    /// The Jets should be new'd; Event will take ownership.
    void set_jets(const std::vector<Jet*>& jets) {
      _jets = jets;
    }

    /// @brief Add a jet to the jets collection
    ///
    /// The Jet should be new'd; Event will take ownership.
    void add_jet(Jet* j) {
      _jets.push_back(j);
    }

    //@}


    /// @name Missing energy
    //@{

    /// @brief Get the missing momentum vector
    ///
    /// Not _necessarily_ the sum over momenta of final state invisibles
    const P4& missingmom() const {
      return _pmiss;
    }

    /// @brief Set the missing momentum vector
    ///
    /// Not _necessarily_ the sum over momenta of final state invisibles
    void set_missingmom(const P4& pmiss) {
      _pmiss = pmiss;
    }

    /// Get the missing transverse momentum in GeV
    double met() const {
      return missingmom().pT();
    }

    //@}


  };


}




