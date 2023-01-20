#ifndef UserCode_ICHiggsTauTau_ICElectronHcalDepthCalculator_h
#define UserCode_ICHiggsTauTau_ICElectronHcalDepthCalculator_h

#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"

/**
 * @brief Produces an edm::ValueMap<float> for the electron HCAL energy depth
 *sum
 *
 * **Example usage**
 * @snippet python/default_producers_cfi.py ElectronHcalDepth
 */
class ICElectronHcalDepthCalculator : public edm::one::EDProducer<> {
 public:
  explicit ICElectronHcalDepthCalculator(const edm::ParameterSet &);
  ~ICElectronHcalDepthCalculator();

 private:
  virtual void beginJob();
  virtual void produce(edm::Event &, const edm::EventSetup &);
  virtual void endJob();

  edm::InputTag input_;
};

#endif
