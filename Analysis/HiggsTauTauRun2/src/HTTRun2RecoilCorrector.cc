#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/interface/HTTRun2RecoilCorrector.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "boost/format.hpp"
#include "TMVA/Reader.h"
#include "TVector3.h"

namespace ic {

  HTTRun2RecoilCorrector::HTTRun2RecoilCorrector(std::string const& name) : ModuleBase(name), 
    channel_(channel::et),
    strategy_(strategy::paper2013),
    mc_(mc::summer12_53X),
    era_(era::data_2012_rereco) {
    met_label_ = "pfMVAMet";
    jets_label_ = "pfJetsPFlow";
    sample_ = "";
    use_quantile_map_ = false;
    store_boson_pt_ = false;
    disable_recoil_corrs = true;
    disable_met_sys = true;
    is_wjets = false;
    met_scale_mode_ = 0;
    met_res_mode_ = 0;
    njets_mode_ = 0;  
    do_recoil_=true;
  }

  HTTRun2RecoilCorrector::~HTTRun2RecoilCorrector() {
    ;
  }

  int HTTRun2RecoilCorrector::PreAnalysis() {
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "HTTRun2RecoilCorrector" << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    std::cout << boost::format(param_fmt()) % "channel"         % Channel2String(channel_);
    std::cout << boost::format(param_fmt()) % "strategy"        % Strategy2String(strategy_);
    std::cout << boost::format(param_fmt()) % "era"             % Era2String(era_);
    std::cout << boost::format(param_fmt()) % "mc"              % MC2String(mc_);
    std::cout << boost::format(param_fmt()) % "met_label"       % met_label_;
    std::cout << boost::format(param_fmt()) % "jets_label"      % jets_label_;

    std::string process_file;
    std::string syst_file;
    if (strategy_ ==strategy::fall15){
      process_file = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/recoilMvaMEt_76X_newTraining_MG5.root";
      syst_file = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/MEtSys.root";
    } else if (strategy_ == strategy::mssmspring16 || strategy_ == strategy::smspring16 || strategy_==strategy::mssmsummer16 || strategy_ == strategy::smsummer16 || strategy_ == strategy::cpsummer16 ||  strategy_ == strategy::legacy16 || strategy_ == strategy::cpdecays16){
      process_file = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/MvaMET_MG_2016BCD_RooT_5.2.root";
      if(met_label_ == "pfMET"){
          process_file = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/TypeI-PFMet_Run2016BtoH.root";
      }
      syst_file    = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/PFMEtSys_2016.root";
    } else if (strategy_ == strategy::cpsummer17 || strategy_ == strategy::cpdecays17) {
      process_file = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/Type1_PFMET_2017.root";
      syst_file    = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/PFMEtSys_2017.root";
    } else if (strategy_ == strategy::cpdecays18) {
      process_file = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/TypeI-PFMet_Run2018.root";
      syst_file    = "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/recoilfits/PFMEtSys_2017.root"; // update when available
    } else{
      std::cerr << "Strategy: " << Strategy2String(strategy_) << " not recognised, an exception will be thrown." << std::endl;
      throw;
    }

   if(met_scale_mode_ > 0 || met_res_mode_ > 0) disable_met_sys =false;
   
    if ( (sample_.find("WJetsToLNu") != sample_.npos) || (sample_.find("W1JetsToLNu") != sample_.npos) || (sample_.find("W2JetsToLNu")!=sample_.npos) || (sample_.find("W3JetsToLNu")!=sample_.npos) || (sample_.find("W4JetsToLNu")!=sample_.npos) || (sample_.find("WG")!=sample_.npos) || (sample_.find("EWKW")!=sample_.npos)){
      disable_recoil_corrs = false;
      is_wjets = true;
    }

    if ( (sample_.find("DY")!=sample_.npos && sample_.find("JetsToLL")!=sample_.npos) || (sample_.find("EWKZ")!=sample_.npos) ) {
      disable_recoil_corrs = false;
    }

    if (sample_.find("HToTauTau")!=sample_.npos || sample_.find("VBFH")!=sample_.npos || sample_.find("GluGluH")!=sample_.npos || sample_.find("WHiggs")!=sample_.npos || sample_.find("ZHiggs")!=sample_.npos){ 
     disable_recoil_corrs = false;
    }
    
    disable_recoil_corrs = disable_recoil_corrs || !do_recoil_;
    std::cout<<"disable recoil ? "<< disable_recoil_corrs<<std::endl;
    if (disable_recoil_corrs && disable_met_sys) {
      std::cout << boost::format(param_fmt()) % "Recoil corrs enabled"      % false;
      std::cout << boost::format(param_fmt()) % "Met systs enabled"         % false;
      return 0;
    } else {
     if(!disable_recoil_corrs){
       corrector_ = new RecoilCorrectorRun2(process_file);
    }
     if(!disable_met_sys){
       metSys_ = new MEtSys(syst_file);
     }
   }

    return 0;
  }

  int HTTRun2RecoilCorrector::Execute(TreeEvent *event) {


    // Get the stuff we need from the event
    Met * met = event->GetPtr<Met>(met_label_);
    std::vector<GenParticle *> parts;
    //bool is_mc = event->ExistsInTree("genParticles");
    if(event->ExistsInTree("genParticles")) parts = event->GetPtrVec<GenParticle>("genParticles");
    std::vector<GenParticle *> sel_gen_parts;
    std::vector<GenParticle *> sel_vis_parts;
    double genpX=0;
    double genpY=0;
    double genM=0;
    double genpT=0;
    double genEta=0;
    double genPhi=0;
    double vispX=0;
    double vispY=0;
    
    // for data events genpX and genpY should always equal zero, which we assume to be correct for Z->mumu events

    for(unsigned i = 0; i < parts.size(); ++i){
      std::vector<bool> status_flags = parts[i]->statusFlags();
      unsigned id = abs(parts[i]->pdgid());
      unsigned status = abs(parts[i]->status());
      if ( (id >= 11 && id <= 16 && status_flags[FromHardProcess] && status==1) || status_flags[IsDirectHardProcessTauDecayProduct]) sel_gen_parts.push_back(parts[i]);
      if ( ( (id == 11 || id == 13 || id == 15) && status_flags[FromHardProcess] && status==1) || (status_flags[IsDirectHardProcessTauDecayProduct] && !(id==12||id==14||id==16))) sel_vis_parts.push_back(parts[i]);
   }

   
  ROOT::Math::PtEtaPhiEVector gen_boson;
   for( unsigned i = 0; i < sel_gen_parts.size() ; ++i){
     genpX+= sel_gen_parts[i]->vector().px();
     genpY+= sel_gen_parts[i]->vector().py();
     gen_boson += sel_gen_parts[i]->vector();
   }
    genpT = gen_boson.pt();
    genM = gen_boson.M();
    genEta = gen_boson.Rapidity();
    genPhi = gen_boson.Phi();


   for( unsigned i = 0; i < sel_vis_parts.size() ; ++i){
     vispX+= sel_vis_parts[i]->vector().px();
     vispY+= sel_vis_parts[i]->vector().py();
   }
  
  if(store_boson_pt_){
    if(!event->ExistsInEvent("genpX")) event->Add("genpX", genpX);
    if(!event->ExistsInEvent("genpY")) event->Add("genpY", genpY);
    event->ForceAdd("vispX", vispX);
    event->ForceAdd("vispY", vispY);
  }

  if(!event->ExistsInEvent("genpT")) event->Add("genpT", genpT);
  if(!event->ExistsInEvent("genM")) event->Add("genM", genM);
  if(!event->ExistsInEvent("genEta")) event->Add("genEta", genEta);
  if(!event->ExistsInEvent("genPhi")) event->Add("genPhi", genPhi);

  float U1=0., U2=0., metU1=0., metU2=0.;


  double Metx = met->vector().px();
  double Mety = met->vector().py();

  corrector_->CalculateU1U2FromMet(Metx,
                       Mety,
                       genpX,
                       genpY,
                       vispX,
                       vispY,
                       U1,
                       U2,
                       metU1,
                       metU2);


  if(!event->ExistsInEvent("U1")) event->Add("U1", (double)U1);
  if(!event->ExistsInEvent("U2")) event->Add("U2", (double)U2);


  //if (disable_recoil_corrs && disable_met_sys) return 0;
  return 0; // do not actually apply recoil corrections for ML project 
 
  std::vector<PFJet*> jets = event->GetPtrVec<PFJet>(jets_label_); // Make a copy of the jet collection
  ic::erase_if(jets,!boost::bind(MinPtMaxEta, _1, 30.0, 4.7));
  unsigned njets = jets.size();
  if(is_wjets) njets+=1;

  //double Metx = met->vector().px();
  //double Mety = met->vector().py();
  double Mete = met->energy();
  double MetPt = met->vector().Pt();
  double MetPhi = met->vector().Phi();
  // add uncorrected MET and MET phi to event
  if(!event->ExistsInEvent("met_norecoil")) event->Add("met_norecoil", MetPt);
  if(!event->ExistsInEvent("met_phi_norecoil")) event->Add("met_phi_norecoil", MetPhi);
  double met_res_e = sqrt(Mete*Mete-Metx*Metx-Mety*Mety);//As we can only set pt/eta/phi/energy, not px or py
  float correctedMetx, correctedMety;
  if(!disable_recoil_corrs){
    if(!use_quantile_map_){
      corrector_->CorrectByMeanResolution(Metx,Mety,genpX,genpY,vispX,vispY,njets,correctedMetx,correctedMety); 
    } else {
      corrector_->Correct(Metx,Mety,genpX,genpY,vispX,vispY,njets,correctedMetx,correctedMety); 
    }
    //Now stick this back into our met object:
    met->set_pt(sqrt(correctedMetx*correctedMetx+correctedMety*correctedMety));
    met->set_phi(atan2(correctedMety,correctedMetx));
    met->set_energy(sqrt(met_res_e*met_res_e+correctedMetx*correctedMetx+correctedMety*correctedMety));
  } else {
    correctedMetx = Metx;
    correctedMety = Mety;
  }
    
  //Apply systematic shifts to MET if requested 
  if(met_scale_mode_ > 0 || met_res_mode_ >0) {
    if((njets_mode_==1 && njets!=0) || (njets_mode_==2 && njets!=1) || (njets_mode_==3 && njets<2)) return 0;
    float met_Shift_x, met_Shift_y;
    MEtSys::SysType sysType;
    MEtSys::SysShift sysShift;
    // MEtSys::ProcessType processType;
    // if(is_wjets || (sample_.find("DY")!=sample_.npos && sample_.find("JetsToLL")!=sample_.npos) || sample_.find("EWKZ")!=sample_.npos  || sample_.find("HToTauTau")!=sample_.npos || sample_.find("VBFH")!=sample_.npos || sample_.find("GluGluH")!=sample_.npos || sample_.find("WHiggs")!=sample_.npos || sample_.find("ZHiggs")!=sample_.npos ){
    //     processType = MEtSys::ProcessType::BOSON;
    // } else if (sample_.find("TT")!=sample_.npos){
    //    processType = MEtSys::ProcessType::TOP;
    // } else {
    //    processType = MEtSys::ProcessType::EWK;
    // }
    
    if(met_scale_mode_ == 1) {
        sysType = MEtSys::SysType::Resolution;
        sysShift = MEtSys::SysShift::Down;
    }
    if(met_scale_mode_ == 2) {
        sysType = MEtSys::SysType::Resolution;
        sysShift = MEtSys::SysShift::Up;
    }
    if(met_res_mode_ == 1) {
        sysType = MEtSys::SysType::Response;
        sysShift = MEtSys::SysShift::Down;
    }
    if(met_res_mode_ == 2) {
        sysType = MEtSys::SysType::Response;
        sysShift = MEtSys::SysShift::Up;
    }
    
    metSys_->ApplyMEtSys(
        correctedMetx,correctedMety, // (float) mva met, use RECOIL CORRECTED value for the Higgs / DY / W+Jets MC
        float(genpX),float(genpY), // (float) transverse momentum of the full leptonic system
        float(vispX),float(vispY), // (float) transverse momentum of the visible leptonic system
        njets, // (int) number of jets : pT > 30 GeV, eta<4.7, loose PF JetID
        // processType, // (int) type of process 
        sysType, // (int) type of systematic uncertainty
        sysShift, // (int) direction of systematic shift
        met_Shift_x,met_Shift_y // (float) shifted value of the mva met
    );
    met->set_pt(sqrt(met_Shift_x*met_Shift_x+met_Shift_y*met_Shift_y));
    met->set_phi(atan2(met_Shift_y,met_Shift_x));
    met->set_energy(sqrt(met_res_e*met_res_e+met_Shift_x*met_Shift_x+met_Shift_y*met_Shift_y));
  
  }


    return 0;
  }
  int HTTRun2RecoilCorrector::PostAnalysis() {

    return 0;
  }

  void HTTRun2RecoilCorrector::PrintInfo() {
    ;
  }
}
