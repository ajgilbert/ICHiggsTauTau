#ifndef UserCode_ICHiggsTauTau_ICSingleMetProducer_h
#define UserCode_ICHiggsTauTau_ICSingleMetProducer_h

#include <memory>
#include <vector>
#include <string>
#include "boost/functional/hash.hpp"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/METReco/interface/MET.h"
#include "UserCode/ICHiggsTauTau/interface/Met.hh"

/**
 * @brief See documentation [here](\ref objs-met)
 */
class ICSingleMetProducer : public edm::stream::EDProducer<> {
 public:
  explicit ICSingleMetProducer(const edm::ParameterSet&);
  ~ICSingleMetProducer();

 private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  ic::Met* met_;
  edm::InputTag input_;
  std::string branch_;
  boost::hash<reco::MET const*> met_hasher_;
};

#endif
