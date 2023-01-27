#ifndef ICHiggsTauTau_Photon_hh
#define ICHiggsTauTau_Photon_hh
#include <vector>
#include "UserCode/ICHiggsTauTau/interface/Candidate.hh"
#include "Rtypes.h"

namespace ic {

/**
 * @brief This class stores a subset of the reco::Photon
 * properties which are most commonly used in analysis.
 */
class Photon : public Candidate {
 public:
  Photon();
  virtual ~Photon();
  virtual void Print() const;

  /// @name Properties
  /**@{*/
  /// PF isolation, using all charged particles in a cone with                                                                                               
  /// \f$ \Delta R = 0.3 \f$                                                                                                                                 
  inline float dr03_pfiso_charged_all() const {
    return dr03_pfiso_charged_all_;
  }

  /// PF isolation, using charged hadrons in a cone with \f$ \Delta R = 0.3 \f$                                                                              
  inline float dr03_pfiso_charged() const { return dr03_pfiso_charged_; }

  /// PF isolation, using neutral hadrons in a cone with \f$ \Delta R = 0.3 \f$                                                                              
  inline float dr03_pfiso_neutral() const { return dr03_pfiso_neutral_; }

  /// PF isolation, using photons in a cone with \f$ \Delta R = 0.3 \f$                                                                                      
  inline float dr03_pfiso_gamma() const { return dr03_pfiso_gamma_; }

  /// PF isolation, using charged pileup in a cone with \f$ \Delta R = 0.3 \f$                                                                               
  inline float dr03_pfiso_pu() const { return dr03_pfiso_pu_; }

  /// PF isolation, using all charged particles in a cone with                                                                                               
  /// \f$ \Delta R = 0.4 \f$                                                                                                                                 
  inline float dr04_pfiso_charged_all() const {
    return dr04_pfiso_charged_all_;
  }

  /// PF isolation, using charged hadrons in a cone with \f$ \Delta R = 0.4 \f$                                                                              
  inline float dr04_pfiso_charged() const { return dr04_pfiso_charged_; }

  /// PF isolation, using neutral hadrons in a cone with \f$ \Delta R = 0.4 \f$                                                                              
  inline float dr04_pfiso_neutral() const { return dr04_pfiso_neutral_; }

  /// PF isolation, using photons in a cone with \f$ \Delta R = 0.4 \f$                                                                                      
  inline float dr04_pfiso_gamma() const { return dr04_pfiso_gamma_; }

  /// PF isolation, using charged pileup in a cone with \f$ \Delta R = 0.4 \f$                                                                               
  inline float dr04_pfiso_pu() const { return dr04_pfiso_pu_; }

  /// Hadronic over electromagnetic energy fraction
  inline float had_tower_over_em() const { return had_tower_over_em_; }

  /// Weighted cluster RMS in the \f$\eta\f$-direction, measured in units
  /// of crystals in a 5x5 block centred on the seed crystal
  inline float sigma_IetaIeta() const { return sigma_IetaIeta_; }

  /// True if photon passes the electron veto
  inline bool pass_electron_veto() const { return pass_electron_veto_; }

  inline std::size_t supercluster() const { return supercluster_; }

  /// The ratio of the energy in the 3x3 crystal array centred on the seed
  /// over the energy of the complete supercluster
  inline float r9() const { return r9_; }

  /**@}*/

  /// @name Setters
  /**@{*/
  /// @copybrief dr03_pfiso_charged_all()                                                                                                                    
  inline void set_dr03_pfiso_charged_all(float const& dr03_pfiso_charged_all) {
    dr03_pfiso_charged_all_ = dr03_pfiso_charged_all;
  }

  /// @copybrief dr03_pfiso_charged()                                                                                                                        
  inline void set_dr03_pfiso_charged(float const& dr03_pfiso_charged) {
    dr03_pfiso_charged_ = dr03_pfiso_charged;
  }

  /// @copybrief dr03_pfiso_neutral()                                                                                                                        
  inline void set_dr03_pfiso_neutral(float const& dr03_pfiso_neutral) {
    dr03_pfiso_neutral_ = dr03_pfiso_neutral;
  }

  /// @copybrief dr03_pfiso_gamma()                                                                                                                          
  inline void set_dr03_pfiso_gamma(float const& dr03_pfiso_gamma) {
    dr03_pfiso_gamma_ = dr03_pfiso_gamma;
  }

  /// @copybrief dr03_pfiso_pu()                                                                                                                             
  inline void set_dr03_pfiso_pu(float const& dr03_pfiso_pu) {
    dr03_pfiso_pu_ = dr03_pfiso_pu;
  }

  /// @copybrief dr04_pfiso_charged_all()                                                                                                                    
  inline void set_dr04_pfiso_charged_all(float const& dr04_pfiso_charged_all) {
    dr04_pfiso_charged_all_ = dr04_pfiso_charged_all;
  }
  
  /// @copybrief dr04_pfiso_charged()                                                                                                                        
  inline void set_dr04_pfiso_charged(float const& dr04_pfiso_charged) {
    dr04_pfiso_charged_ = dr04_pfiso_charged;
  }
  
  /// @copybrief dr04_pfiso_neutral()                                                                                                                        
  inline void set_dr04_pfiso_neutral(float const& dr04_pfiso_neutral) {
    dr04_pfiso_neutral_ = dr04_pfiso_neutral;
  }
  
  /// @copybrief dr04_pfiso_gamma()                                                                                                                          
  inline void set_dr04_pfiso_gamma(float const& dr04_pfiso_gamma) {
    dr04_pfiso_gamma_ = dr04_pfiso_gamma;
  }
  
  /// @copybrief dr04_pfiso_pu()                                                                                                                             
  inline void set_dr04_pfiso_pu(float const& dr04_pfiso_pu) {
    dr04_pfiso_pu_ = dr04_pfiso_pu;
  }

  /// @copybrief had_tower_over_em()
  inline void set_had_tower_over_em(float const& had_tower_over_em) {
    had_tower_over_em_ = had_tower_over_em;
  }

  /// @copybrief sigma_IetaIeta()
  inline void set_sigma_IetaIeta(float const& sigma_IetaIeta) {
    sigma_IetaIeta_ = sigma_IetaIeta;
  }

  /// @copybrief pass_electron_veto()
  inline void set_pass_electron_veto(bool const& pass_electron_veto) {
    pass_electron_veto_ = pass_electron_veto;
  }

  inline void set_supercluster(std::size_t const& supercluster) {
    supercluster_ = supercluster;
  }

  /// @copybrief r9()
  inline void set_r9(float const& r9) { r9_ = r9; }


//supercluster_
  /**@}*/

 private:
  float dr03_pfiso_charged_all_;
  float dr03_pfiso_charged_;
  float dr03_pfiso_neutral_;
  float dr03_pfiso_gamma_;
  float dr03_pfiso_pu_;
  float dr04_pfiso_charged_all_;
  float dr04_pfiso_charged_;
  float dr04_pfiso_neutral_;
  float dr04_pfiso_gamma_;
  float dr04_pfiso_pu_;

  float had_tower_over_em_;
  float sigma_IetaIeta_;
  float r9_;

  bool pass_electron_veto_;
  std::size_t supercluster_;

 #ifndef SKIP_CINT_DICT
 public:
  ClassDef(Photon, 3);
 #endif
};

typedef std::vector<ic::Photon> PhotonCollection;
}
/** \example plugins/ICPhotonProducer.cc */
#endif
