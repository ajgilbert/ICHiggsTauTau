#include <iostream>
#include <fstream>
#include <sstream>
#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/interface/HinvPrint.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/ElectronEffectiveArea.h"
#include "TMVA/Reader.h"
#include "TVector3.h"


namespace ic {

  HinvPrint::HinvPrint(std::string const& name, bool is_data, bool filter, bool runLumiEvt) : ModuleBase(name) {
    runLumiEvt_ = runLumiEvt;
    filter_ = filter;
    is_data_ = is_data;
  }

  HinvPrint::~HinvPrint() {
    ;
  }

  int HinvPrint::PreAnalysis() {
    if (runLumiEvt_) {
      std::ostringstream outName;
      outName << "output_synch/EventList_";
      if (is_data_) outName << "Data";
      else outName << "MC";
      outName << ".txt" ;
      foutList_.open(outName.str());
      if (!foutList_.is_open()){
        std::cerr << " -- Failed to open file " << outName.str() << " for printing event list." << std::endl;
        if (is_data_){
          foutList_ << "Run lumiBlock Event" << std::endl;
        }
        else {
          foutList_ << "Event TotalWeight PUWeight TrigWeight IdIsoWeight" << std::endl;
        }
        throw;
      }
    }
    return 0;
  }

  void HinvPrint::PrintEvent(unsigned run, unsigned lumi, unsigned evt) {
    RunLumiEvent lRLE;
    lRLE.run = run;
    lRLE.lumi = lumi;
    lRLE.evt = evt;
    events_.insert(lRLE);
  }

  int HinvPrint::Execute(TreeEvent *event) {

    EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");

    if (runLumiEvt_) {
      if (is_data_){
        foutList_ << eventInfo->run() << " " << eventInfo->lumi_block() << " " << eventInfo->event() << std::endl;
      }
      else {
        foutList_ << eventInfo->event() << " " << eventInfo->total_weight() << " " << eventInfo->weight("pileup") << " " << eventInfo->weight("trigger")<< " " << eventInfo->weight("idisoTight") << std::endl;
      }
      //std::cout << "-- Event passed: " << eventInfo->run() << " " << eventInfo->lumi_block() << " " << eventInfo->event() << std::endl;
      return 0; 
    }

 
    RunLumiEvent lRLE;
    lRLE.run = eventInfo->run();
    lRLE.lumi = eventInfo->lumi_block();
    if (!is_data_) {
      lRLE.run = 0;
      lRLE.lumi = 0;
    }
    lRLE.evt = eventInfo->event();
    std::cout<<"run: "<<lRLE.run<<", lumi: "<<lRLE.lumi<<", event: "<<lRLE.evt<<std::endl;
    if (events_.find(lRLE) != events_.end() || events_.size()==0) {

      if (filter_) return 0;

      std::ostringstream outName;
      outName << "output_synch/Run" << eventInfo->run() << "_Lumi" << eventInfo->lumi_block() << "_Event" << eventInfo->event() << ".txt" ;
      std::ofstream fout(outName.str());
      if (!fout.is_open()){
        std::cerr << " -- Failed to open file " << outName.str() << " for printing event " << eventInfo->event() << std::endl;
        throw;
      }
      std::vector<Muon*> const& muons = event->GetPtrVec<Muon>("muonsPFlow");
      fout << "-----------------------------------------" << std::endl;
      fout << "event: " <<  eventInfo->event() << " lumi: " << eventInfo->lumi_block() << " run: " << eventInfo->run() << std::endl;
      fout << "-----------------------------------------" << std::endl;
      fout << "nGoodVertices: " << eventInfo->good_vertices() << std::endl;
      fout << "Jet Rho: " << eventInfo->jet_rho() << std::endl;
      fout << "Lepton Rho: " << eventInfo->lepton_rho() << std::endl;
      Met *lMet = event->GetPtr<Met>("pfMetType1");
      fout << "PF MET type0+1 : " << lMet->pt() << std::endl;
      Met *lMetNoMu = event->GetPtr<Met>("metNoMuons");
      fout << "PF MET no muons : " << lMetNoMu->pt() << std::endl;

      fout << "-----------------------------------------" << std::endl;
      fout << "   Printing collection muonsPFlow" << std::endl;
      fout << "-----------------------------------------" << std::endl;

      for (unsigned i = 0; i < muons.size(); ++i) {
        fout << "Muon " << i << " " << muons[i]->vector() << std::endl;
        fout << "-dxyVertex: " << muons[i]->dxy_vertex() << std::endl;
        fout << "-dzVertex: " << muons[i]->dz_vertex() << std::endl;
        fout << "-isGlobalMuon: " << muons[i]->is_global() << std::endl;
        fout << "-isTrackerMuon: " << muons[i]->is_tracker() << std::endl;
        fout << "-numberOfValidPixelHits: " << muons[i]->it_pixel_hits() << std::endl;
        fout << "-numberOfValidMuonHits: " << muons[i]->gt_valid_muon_hits() << std::endl;
        fout << "-trackerLayersWithMeasurement: " << muons[i]->it_layers_with_measurement() << std::endl;
        fout << "-normalizedChi2: " << muons[i]->gt_normalized_chi2() << std::endl;
        fout << "-numberOfMatchedStations: " << muons[i]->matched_stations() << std::endl;
        fout << "-Muon is tight: " << MuonTight(muons[i]) << std::endl;
        fout << "-dr04_pfiso_charged_all " << muons[i]->dr04_pfiso_charged_all() << std::endl;
        fout << "-dr04_pfiso_charged " << muons[i]->dr04_pfiso_charged() << std::endl;
        fout << "-dr04_pfiso_neutral " << muons[i]->dr04_pfiso_neutral() << std::endl;
        fout << "-dr04_pfiso_gamma " << muons[i]->dr04_pfiso_gamma() << std::endl;
        fout << "-dr04_pfiso_pu " << muons[i]->dr04_pfiso_pu() << std::endl;

        double iso =  muons[i]->dr04_pfiso_charged_all() + std::max(muons[i]->dr04_pfiso_neutral() + muons[i]->dr04_pfiso_gamma() - 0.5 * muons[i]->dr04_pfiso_pu(), 0.0);
        iso = iso / muons[i]->pt();

        fout << "-isolation " << iso << std::endl;

      }

      std::vector<Electron*> const& elecs = event->GetPtrVec<Electron>("electrons");
      fout << "-----------------------------------------" << std::endl;
      fout << "   Printing collection electrons" << std::endl;
      fout << "-----------------------------------------" << std::endl;

      for (unsigned i = 0; i < elecs.size(); ++i) {
        fout << "Elec electrons" << i << " " << elecs[i]->vector() << std::endl;
        elecs[i]->Print();
        fout << "-dxyVertex: " << elecs[i]->dxy_vertex() << std::endl;
        fout << "-dzVertex: " << elecs[i]->dz_vertex() << std::endl;
        fout << "-VetoElectronID: " << VetoElectronID(elecs[i]) << std::endl;
        fout << "-TightElectronID=Electron2011WP70ID: " << Electron2011WP70ID(elecs[i]) << std::endl;
        fout << "-dr03_pfiso_charged " << elecs[i]->dr03_pfiso_charged() << std::endl;
        fout << "-dr03_pfiso_neutral " << elecs[i]->dr03_pfiso_neutral() << std::endl;
        fout << "-dr03_pfiso_gamma " << elecs[i]->dr03_pfiso_gamma() << std::endl;	

        float lEffArea = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , elecs[i]->sc_eta() , ElectronEffectiveArea::kEleEAData2012);

        double lIso = elecs[i]->dr03_pfiso_charged() + std::max(elecs[i]->dr03_pfiso_gamma() + elecs[i]->dr03_pfiso_neutral() - eventInfo->lepton_rho() * lEffArea, 0.);
        lIso = lIso / elecs[i]->pt();
        fout << "-isolation " << lIso << std::endl;

      }

    for (unsigned i = 0; i < elecs.size(); ++i) {
      for (unsigned j = 0; j < muons.size(); ++j) {
        fout << "Elec " << i << ", Muon " << j << " DR: " <<
          ROOT::Math::VectorUtil::DeltaR(elecs[i]->vector(), muons[j]->vector()) << std::endl;
      }
    }

    std::vector<Tau*> const& taus = event->GetPtrVec<Tau>("taus");
    fout << "-----------------------------------------" << std::endl;
    fout << "   Printing collection taus" << std::endl;
    fout << "-----------------------------------------" << std::endl;

    for(unsigned itau=0;itau<taus.size();itau++){
      fout << "Tau "<<itau<<" "<<taus[itau]->vector()<<std::endl;
    }

    std::vector<PFJet*> const& jets = event->GetPtrVec<PFJet>("pfJetsPFlow");
//     std::vector<PFJet*> matched_jets;
//     std::vector< std::pair<PFJet*, GenJet*> > matches;
//     if (!eventInfo->is_data()) {
//       std::vector<GenJet*> const& gen_jets = event->GetPtrVec<GenJet>("genJets");
//       matches = MatchByDR(jets, gen_jets, 0.5, true, true);
//       matched_jets = ExtractFirst(matches);
//     }

      fout << "-----------------------------------------" << std::endl;
      fout << "   Printing collection pfJetsPFlow" << std::endl;
      fout << "-----------------------------------------" << std::endl;

    for (unsigned i = 0; i < jets.size(); ++i) {
      fout << "Jet " << i << ": " << jets[i]->vector() << std::endl;
      std::cout << " - Event " << eventInfo->event() << " jet " << i << " pT = " << jets[i]->pt() << " E_T = " << jets[i]->vector().Et() << " mass " << jets[i]->M() << std::endl;
      fout << "-pileup id loose: " << jets[i]->pu_id_mva_loose() << std::endl;
      fout << "-beta: " << jets[i]->beta() << std::endl;
      fout << "-beta_max: " << jets[i]->beta_max() << std::endl;
      //double charged_frac = (jets[i]->charged_em_energy() + jets[i]->charged_had_energy()) / jets[i]->uncorrected_energy();
      //fout << "-charged_frac: " << charged_frac << std::endl;
      fout << "-Uncorrected: " << jets[i]->GetJecFactor("Uncorrected") << std::endl;
      fout << "-L1FastJet: " << jets[i]->GetJecFactor("L1FastJet") << std::endl;
      fout << "-L2Relative: " << jets[i]->GetJecFactor("L2Relative") << std::endl;
      fout << "-L3Absolute: " << jets[i]->GetJecFactor("L3Absolute") << std::endl;
      fout << "-L2L3Residual: " << jets[i]->GetJecFactor("L2L3Residual") << std::endl;
      fout << "-PFJetID: " << PFJetID(jets[i]) << std::endl;
      int n_pu = jets[i]->charged_multiplicity() + jets[i]->neutral_multiplicity() + jets[i]->HF_had_multiplicity() + jets[i]->HF_em_multiplicity();
      double neutralFrac = ( jets[i]->neutral_had_energy() + jets[i]->HF_had_energy() ) / jets[i]->uncorrected_energy();

      fout << "-neutralFrac: " << neutralFrac  << std::endl;
      fout << "-neutral_em_energy_frac: " << jets[i]->neutral_em_energy_frac() << std::endl;
      fout << "-n_pu: " << n_pu  << std::endl;
      fout << "-charged_had_energy_frac: " << jets[i]->charged_had_energy_frac() << std::endl;
      fout << "-charged_multiplicity: " << jets[i]->charged_multiplicity() << std::endl;
      fout << "-charged_em_energy_frac: " << jets[i]->charged_em_energy_frac() << std::endl;

//       std::vector<PFJet*>::const_iterator it = std::find(matched_jets.begin(),matched_jets.end(), jets[i]);
//       if (it != matched_jets.end()) {
//         GenJet const* matched_genjet = NULL;
//         for (unsigned j = 0; j < matches.size(); ++j) {
//           if (matches[j].first == (*it)) matched_genjet = matches[j].second;
//         }
//         fout << "-has matched genjet: " << matched_genjet->vector() << std::endl;
//       }

      for (unsigned j = 0; j < elecs.size(); ++j) {
        fout << "Elec " << j << " DR: " <<
          ROOT::Math::VectorUtil::DeltaR(jets[i]->vector(), elecs[j]->vector()) << std::endl;
      }
      for (unsigned j = 0; j < muons.size(); ++j) {
        fout << "Muon " << j << " DR: " <<
          ROOT::Math::VectorUtil::DeltaR(jets[i]->vector(), muons[j]->vector()) << std::endl;
      }

    }

    std::vector<CompositeCandidate *> const& dijet_vec = event->GetPtrVec<CompositeCandidate>("jjLeadingCandidates");
    if (dijet_vec.size() != 0) {
      fout << "-----------------------------------------" << std::endl;
      fout << "   Printing leading jet pair" << std::endl;
      fout << "-----------------------------------------" << std::endl;

      CompositeCandidate const* dijet = dijet_vec.at(0);
      Candidate const* jet1 = dijet->GetCandidate("jet1");
      Candidate const* jet2 = dijet->GetCandidate("jet2");
      fout << "-Jet " << 1 << ": " << jet1->vector() << std::endl;
      fout << "-Jet " << 2 << ": " << jet2->vector() << std::endl;
      fout << "-Dijet mass " << dijet->M() << std::endl;
      fout << "-Deta " << fabs(jet1->eta() - jet2->eta()) << std::endl;
      fout << "etaprod " << jet1->eta() * jet2->eta() << std::endl;
      fout << "DR " << ROOT::Math::VectorUtil::DeltaR(jet1->vector(),jet2->vector()) << std::endl;
      fout << "Dphi " << fabs(ROOT::Math::VectorUtil::DeltaPhi(jet1->vector(),jet2->vector())) << std::endl;
    }

    /*
    fout << "---Triggers" << std::endl;
    std::vector<std::string> paths = { "triggerObjectsIsoMu17LooseTau20", "triggerObjectsEle22WP90RhoLooseTau20" };
    for (auto path : paths) {
      if (!eventInfo->is_data()) {
        fout << "--" << path << std::endl;
        auto objs = event->GetPtrVec<TriggerObject>(path);
        for (auto cand : objs) {
          cand->Print();
          for (auto label : cand->filters()) {
            fout << "-" << label << std::endl;
          }
        }
      }
    }*/

    fout.close();
    return 0;

    }
    if (filter_) return 1;
    return 0;

  }
  int HinvPrint::PostAnalysis() {
    foutList_.close();
    return 0;
  }

  void HinvPrint::PrintInfo() {
    ;
  }
}
