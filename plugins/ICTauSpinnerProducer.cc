#include "UserCode/ICHiggsTauTau/plugins/ICTauSpinnerProducer.hh"
#include <memory>
#include <string>
#include <vector>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "UserCode/ICHiggsTauTau/interface/StaticTree.hh"
#include "UserCode/ICHiggsTauTau/interface/city.h"
#include "UserCode/ICHiggsTauTau/plugins/Consumes.h"
#include "UserCode/ICHiggsTauTau/plugins/PrintConfigTools.h"


ICTauSpinnerProducer::ICTauSpinnerProducer(const edm::ParameterSet& config)
    : input_(config.getParameter<edm::InputTag>("input")),
      branch_(config.getParameter<std::string>("branch")),
      theta_(config.getParameter<std::string>("theta")){
  consumes<edm::View<reco::GenParticle>>(input_);

  PrintHeaderWithProduces(config, input_, branch_);
  TauSpinnerSettingsPDF="NNPDF30_nlo_as_0118";
  Ipp=true;
  Ipol=1;
  nonSM2=0;
  nonSMN=0;
  CMSENE=13000.0;
  bosonPdgId_=25;
  info_ = new ic::EventInfo();
}

ICTauSpinnerProducer::~ICTauSpinnerProducer() {}

std::vector<std::pair<std::string,double>> ICTauSpinnerProducer::SplitString(std::string instring){
  std::vector<std::pair<std::string,double>> out;
  std::stringstream ss(instring);   
  std::string splitstring;  
  while(std::getline(ss, splitstring, ',')){
    double val = std::stod(splitstring);
    if(splitstring.find(".") != std::string::npos) splitstring.replace(splitstring.find("."),1,"p");
    std::string weight_name = "wt_cp_"+splitstring;    
    out.push_back(std::make_pair(weight_name,val)); 
  }
  return out;
}

reco::GenParticle ICTauSpinnerProducer::getBoson(edm::Handle<edm::View<reco::GenParticle> > parts_handle){
  reco::GenParticle boson;
  bool foundBoson = false;
  for (unsigned i = 0; i < parts_handle->size(); ++i) {
    reco::GenParticle const& part = parts_handle->at(i);
    bool isLastCopy = part.isLastCopy();
    int part_pdgid = fabs(part.pdgId());
    if (part_pdgid==bosonPdgId_&&isLastCopy){
        boson = part;
        foundBoson = true;
        break;
    }
  }
  if(!foundBoson) {std::cout << "ICTauSpinnerProducer: Gen boson not found. Throwing exception." << std::endl; throw;}
  return boson;
}

std::vector<reco::GenParticle> ICTauSpinnerProducer::getTaus(reco::GenParticle boson, edm::Handle<edm::View<reco::GenParticle> > parts_handle){
  std::vector<reco::GenParticle> taus;
  unsigned nTaus=0;
  for (unsigned i = 0; i < boson.daughterRefVector().size(); ++i) {
    int daughter_index = static_cast<int>(boson.daughterRefVector().at(i).key());
    reco::GenParticle daughter = parts_handle->at(daughter_index);
    int daughter_pdgid = fabs(daughter.pdgId());
    if (daughter_pdgid != 15) continue;
    taus.push_back(daughter); 
    nTaus++;
  }
  if(nTaus!=2) {std::cout << "ICTauSpinnerProducer: Found " << nTaus << " taus, expected 2. Throwing exception." << std::endl; throw;}
  return taus;
}

void ICTauSpinnerProducer::getTauDaughters(std::vector<reco::GenParticle> &tau_daughters, unsigned &type, reco::GenParticle tau, edm::Handle<edm::View<reco::GenParticle> > parts_handle){
  for (unsigned i = 0; i < tau.daughterRefVector().size(); ++i) {
    int daughter_index = static_cast<int>(tau.daughterRefVector().at(i).key());
    reco::GenParticle daughter = parts_handle->at(daughter_index);
    int daughter_pdgid = fabs(daughter.pdgId());
    if(daughter_pdgid == 22 || daughter_pdgid == 111 || daughter_pdgid == 211 || daughter_pdgid == 321 || daughter_pdgid == 130 || daughter_pdgid == 310 || daughter_pdgid == 11 || daughter_pdgid == 12 || daughter_pdgid == 13 || daughter_pdgid == 14 || daughter_pdgid == 16){
      if(daughter_pdgid == 11) type = 1;
      if(daughter_pdgid == 13) type = 2;
      tau_daughters.push_back(daughter);
    }
    else getTauDaughters(tau_daughters, type, daughter, parts_handle);
  }
}

void ICTauSpinnerProducer::removeGammas(std::vector<TauSpinner::SimpleParticle> &tau_daughters){
  for(int i=tau_daughters.size()-1; i>=0; --i){
    int pdgid = fabs(tau_daughters[i].pdgid());
    if(pdgid == 22) tau_daughters.erase(tau_daughters.begin()+i);   
  }
}

void ICTauSpinnerProducer::removeSecondGamma(std::vector<TauSpinner::SimpleParticle> &tau_daughters){
  int gamma_count=0;
  int gamma_index=0;
  double min_gamma_pt=0;
  for(unsigned i=0; i<tau_daughters.size(); ++i){
    int pdgid = fabs(tau_daughters[i].pdgid());
    if(pdgid == 22){
      double gamma_pt=sqrt(pow(tau_daughters[i].px(),2) + pow(tau_daughters[i].py(),2));
      if(gamma_count==0){ min_gamma_pt = gamma_pt; gamma_index=i;}
      else if (gamma_pt<min_gamma_pt){ min_gamma_pt = gamma_pt; gamma_index=i;}
      gamma_count++;
    }
  }
  if(gamma_count==2) tau_daughters.erase(tau_daughters.begin()+gamma_index);   
}

bool ICTauSpinnerProducer::channelMatch(std::vector<reco::GenParticle> parts, std::vector<int> matches, bool signed_){
  std::vector<int> list;
  for(unsigned i=0; i<parts.size();++i) list.push_back(parts[i].pdgId());
                
  for(unsigned i=0;i<matches.size();i++) {
    bool found = false;  
    int match = matches[i];
    if(!signed_) match = fabs(match);
    
    for(unsigned j=0;j<list.size(); j++){
      int pdgid = list[j];
      if(!signed_) pdgid = fabs(pdgid);
      // if pdgid is found - erese it from the list and search for the next one
      if(pdgid==match) {
        found = true;
        list.erase(list.begin()+j);
        break;
      }
    }
    
    if(!found) return false;
  }
  if(list.size()!=0) return false;
  return true;
    
}

TauSpinner::SimpleParticle ICTauSpinnerProducer::ConvertToSimplePart(reco::GenParticle input_part){
  return TauSpinner::SimpleParticle(input_part.px(), input_part.py(), input_part.pz(), input_part.energy(), input_part.pdgId());
}

void ICTauSpinnerProducer::initialize(){
  Tauolapp::Tauola::initialize();
  LHAPDF::initPDFSetByName(TauSpinnerSettingsPDF);
  TauSpinner::initialize_spinner(Ipp, Ipol, nonSM2, nonSMN,  CMSENE);
}

void ICTauSpinnerProducer::produce(edm::Event& event,
                                  const edm::EventSetup& setup) {
  edm::Handle<edm::View<reco::GenParticle> > parts_handle;
  event.getByLabel(input_, parts_handle);
  
  reco::GenParticle boson = getBoson(parts_handle);
  std::vector<reco::GenParticle> taus = getTaus(boson,parts_handle);
  
  std::vector<reco::GenParticle> tau1_daughters;
  std::vector<reco::GenParticle> tau2_daughters;
  unsigned type1 = 0; 
  unsigned type2 = 0;
  getTauDaughters(tau1_daughters,type1,taus[0],parts_handle);  
  getTauDaughters(tau2_daughters,type2,taus[1],parts_handle); 
  
  TauSpinner::SimpleParticle simple_boson = ConvertToSimplePart(boson);
  TauSpinner::SimpleParticle simple_tau1 = ConvertToSimplePart(taus[0]);
  TauSpinner::SimpleParticle simple_tau2 = ConvertToSimplePart(taus[1]);
  std::vector<TauSpinner::SimpleParticle> simple_tau1_daughters;
  std::vector<TauSpinner::SimpleParticle> simple_tau2_daughters;
  
  for(unsigned i=0; i<tau1_daughters.size(); ++i) simple_tau1_daughters.push_back(ConvertToSimplePart(tau1_daughters[i]));
  for(unsigned i=0; i<tau2_daughters.size(); ++i) simple_tau2_daughters.push_back(ConvertToSimplePart(tau2_daughters[i]));
  // remove the photons for hadronic decays as tauspinner doesn't concider these decay modes
  if(!(type1==1||type1==2)) removeGammas(simple_tau1_daughters);
  else removeSecondGamma(simple_tau1_daughters);
  // for leptonic decays remove the second photon if it exists (remove lowest pT photon)
  if(!(type2==1||type2==2)) removeGammas(simple_tau2_daughters);
  else removeSecondGamma(simple_tau2_daughters);
  
  for(unsigned i=0; i<theta_vec_.size(); ++i){
    double theta_val_ = theta_vec_[i].second;
    std::string weight_name_ = theta_vec_[i].first;
    // Can make this more general by having boson pdgid as input or have option for set boson type
    TauSpinner::setHiggsParametersTR(-cos(2*M_PI*theta_val_),cos(2*M_PI*theta_val_),-sin(2*M_PI*theta_val_),-sin(2*M_PI*theta_val_));
    double weight_ = TauSpinner::calculateWeightFromParticlesH(simple_boson,simple_tau1,simple_tau2,simple_tau1_daughters,simple_tau2_daughters); 
    info_->set_weight(weight_name_,weight_,false);
  }
  
  //if(weight_==1){  
  //  std::cout << "****************************************************************************" << std::endl;
  //  std::cout << "boson = " << simple_boson.pdgid() << "  " << simple_boson.px() << "  " << simple_boson.py() << "  " << simple_boson.pz() << "  " << simple_boson.e() << std::endl;
  //  std::cout << "tau 1 = " << simple_tau1.pdgid() << "  " << simple_tau1.px() << "  " << simple_tau1.py() << "  " << simple_tau1.pz() << "  " << simple_tau1.e() << std::endl;
  //  std::cout << "type = " << type1 << std::endl;
  //  std::cout << "tau 1 daughters:" << std::endl;
  //  for(unsigned i=0; i< tau1_daughters.size(); ++i){
  //    std::cout << tau1_daughters[i].pdgId() << "  " << tau1_daughters[i].px() << "  " << tau1_daughters[i].py() << "  " << tau1_daughters[i].pz() << "  " << tau1_daughters[i].energy() << std::endl;
  //    //std::cout << "status = " << tau1_daughters[i].status() << std::endl;
  //    //std::cout << "flags: " <<std::endl;
  //    //std::cout << " IsPrompt                             =  " << tau1_daughters[i].statusFlags().isPrompt                          () << std::endl;  
  //    //std::cout << " IsDecayedLeptonHadron                =  " << tau1_daughters[i].statusFlags().isDecayedLeptonHadron             () << std::endl;
  //    //std::cout << " IsTauDecayProduct                    =  " << tau1_daughters[i].statusFlags().isTauDecayProduct                 () << std::endl;
  //    //std::cout << " IsPromptTauDecayProduct              =  " << tau1_daughters[i].statusFlags().isPromptTauDecayProduct           () << std::endl;
  //    //std::cout << " IsDirectTauDecayProduct              =  " << tau1_daughters[i].statusFlags().isDirectTauDecayProduct           () << std::endl;
  //    //std::cout << " IsDirectPromptTauDecayProduct        =  " << tau1_daughters[i].statusFlags().isDirectPromptTauDecayProduct     () << std::endl;
  //    //std::cout << " IsDirectHadronDecayProduct           =  " << tau1_daughters[i].statusFlags().isDirectHadronDecayProduct        () << std::endl;
  //    //std::cout << " IsHardProcess                        =  " << tau1_daughters[i].statusFlags().isHardProcess                     () << std::endl;
  //    //std::cout << " FromHardProcess                      =  " << tau1_daughters[i].statusFlags().fromHardProcess                   () << std::endl;
  //    //std::cout << " IsHardProcessTauDecayProduct         =  " << tau1_daughters[i].statusFlags().isHardProcessTauDecayProduct      () << std::endl;
  //    //std::cout << " IsDirectHardProcessTauDecayProduct   =  " << tau1_daughters[i].statusFlags().isDirectHardProcessTauDecayProduct() << std::endl;
  //    //std::cout << " FromHardProcessBeforeFSR             =  " << tau1_daughters[i].statusFlags().fromHardProcessBeforeFSR          () << std::endl;
  //    //std::cout << " IsFirstCopy                          =  " << tau1_daughters[i].statusFlags().isFirstCopy                       () << std::endl;
  //    //std::cout << " IsLastCopy                           =  " << tau1_daughters[i].statusFlags().isLastCopy                        () << std::endl;
  //    //std::cout << " IsLastCopyBeforeFSR                  =  " << tau1_daughters[i].statusFlags().isLastCopyBeforeFSR               () << std::endl;
  //  }
  //  std::cout << "tau 2 = " << simple_tau2.pdgid() << "  " << simple_tau2.px() << "  " << simple_tau2.py() << "  " << simple_tau2.pz() << "  " << simple_tau2.e() << std::endl;
  //  std::cout << "type = " << type2 << std::endl;
  //  std::cout << "tau 2 daughters:" << std::endl;
  //  for(unsigned i=0; i< tau2_daughters.size(); ++i){
  //    std::cout << tau2_daughters[i].pdgId() << "  " << tau2_daughters[i].px() << "  " << tau2_daughters[i].py() << "  " << tau2_daughters[i].pz() << "  " << tau2_daughters[i].energy() << std::endl;
  //    //std::cout << "status = " << tau2_daughters[i].status() << std::endl;
  //    //std::cout << "flags: " <<std::endl;
  //    //std::cout << " IsPrompt                             =  " << tau2_daughters[i].statusFlags().isPrompt                          () << std::endl;  
  //    //std::cout << " IsDecayedLeptonHadron                =  " << tau2_daughters[i].statusFlags().isDecayedLeptonHadron             () << std::endl;
  //    //std::cout << " IsTauDecayProduct                    =  " << tau2_daughters[i].statusFlags().isTauDecayProduct                 () << std::endl;
  //    //std::cout << " IsPromptTauDecayProduct              =  " << tau2_daughters[i].statusFlags().isPromptTauDecayProduct           () << std::endl;
  //    //std::cout << " IsDirectTauDecayProduct              =  " << tau2_daughters[i].statusFlags().isDirectTauDecayProduct           () << std::endl;
  //    //std::cout << " IsDirectPromptTauDecayProduct        =  " << tau2_daughters[i].statusFlags().isDirectPromptTauDecayProduct     () << std::endl;
  //    //std::cout << " IsDirectHadronDecayProduct           =  " << tau2_daughters[i].statusFlags().isDirectHadronDecayProduct        () << std::endl;
  //    //std::cout << " IsHardProcess                        =  " << tau2_daughters[i].statusFlags().isHardProcess                     () << std::endl;
  //    //std::cout << " FromHardProcess                      =  " << tau2_daughters[i].statusFlags().fromHardProcess                   () << std::endl;
  //    //std::cout << " IsHardProcessTauDecayProduct         =  " << tau2_daughters[i].statusFlags().isHardProcessTauDecayProduct      () << std::endl;
  //    //std::cout << " IsDirectHardProcessTauDecayProduct   =  " << tau2_daughters[i].statusFlags().isDirectHardProcessTauDecayProduct() << std::endl;
  //    //std::cout << " FromHardProcessBeforeFSR             =  " << tau2_daughters[i].statusFlags().fromHardProcessBeforeFSR          () << std::endl;
  //    //std::cout << " IsFirstCopy                          =  " << tau2_daughters[i].statusFlags().isFirstCopy                       () << std::endl;
  //    //std::cout << " IsLastCopy                           =  " << tau2_daughters[i].statusFlags().isLastCopy                        () << std::endl;
  //    //std::cout << " IsLastCopyBeforeFSR                  =  " << tau2_daughters[i].statusFlags().isLastCopyBeforeFSR               () << std::endl;
  //        
  //  }
  //  std::cout << "weight = " << weight_ << std::endl;
  //  
  //  
  //}
  //bool tau1_dm = channelMatch(tau1_daughters, {211, -211, -211, 111, 111,16},false);
  //bool tau2_dm = channelMatch(tau2_daughters, {211, -211, -211, 111, 111,16},false);
  //if(tau1_dm) count++;
  //if(tau2_dm) count++;
  //total++; total++;
  //std::cout << count << "/" << total << std::endl;

}

void ICTauSpinnerProducer::beginJob() {
  ic::StaticTree::tree_->Branch(branch_.c_str(), &info_);
  theta_vec_ = SplitString(theta_);  
  initialize();  
}

void ICTauSpinnerProducer::endJob() {}

DEFINE_FWK_MODULE(ICTauSpinnerProducer);