#include <vector>
#include <set>
#include <map>
#include <string>
#include <iostream>
#include <memory>
// ROOT
#include "TH1.h"
// Objects
#include "UserCode/ICHiggsTauTau/interface/Electron.hh"
#include "UserCode/ICHiggsTauTau/interface/Muon.hh"
#include "UserCode/ICHiggsTauTau/interface/Tau.hh"
#include "UserCode/ICHiggsTauTau/interface/Met.hh"
#include "UserCode/ICHiggsTauTau/interface/CompositeCandidate.hh"
//boost
#include <boost/format.hpp>
#include "boost/lexical_cast.hpp"
// Utilities
#include "Utilities/interface/FnRootTools.h"
// HTT-specific modules
#include "HiggsTauTauRun2/interface/HTTSequence.h"
#include "HiggsTauTauRun2/interface/HTTTriggerFilter.h"
#include "HiggsTauTauRun2/interface/HTTEnergyScale.h"
#include "HiggsTauTauRun2/interface/HTTCategories.h"
#include "HiggsTauTauRun2/interface/HTTPairSelector.h"
#include "HiggsTauTauRun2/interface/HTTPairGenInfo.h"
#include "HiggsTauTauRun2/interface/BTagCheck.h"
#include "HiggsTauTauRun2/interface/SVFitTest.h"
#include "HiggsTauTauRun2/interface/HTTRun2RecoilCorrector.h"
#include "HiggsTauTauRun2/interface/HTTLegacyRun2RecoilCorrector.h"
#include "HiggsTauTauRun2/interface/HTTWeights.h"
#include "HiggsTauTauRun2/interface/HTTStitching.h"
#include "HiggsTauTauRun2/interface/JetEnergyUncertainty.h"
#include "HiggsTauTauRun2/interface/JetEnergyResolution.h"
#include "HiggsTauTauRun2/interface/HTTPrint.h"
#include "HiggsTauTauRun2/interface/HTTFilter.h"
#include "HiggsTauTauRun2/interface/EffectiveEvents.h"
#include "HiggsTauTauRun2/interface/BTagWeightLegacyRun2.h"
#include "HiggsTauTauRun2/interface/BTagWeightRun2.h"
#include "HiggsTauTauRun2/interface/HTTGenMatchSelector.h"
#include "HiggsTauTauRun2/interface/HTTFakeFactorWeights.h"
#include "HiggsTauTauRun2/interface/HTTGenAnalysis.h"
#include "HiggsTauTauRun2/interface/TagAndProbe.h"
#include "HiggsTauTauRun2/interface/HTTSmearScale.h"
#include "HiggsTauTauRun2/interface/HTTPreFireWeight.h"
#include "HiggsTauTauRun2/interface/HTTMuonEnergyScale.h"
#include "HiggsTauTauRun2/interface/HTTEventClassifier.h"
#include "HiggsTauTauRun2/interface/HTTEventClassifier.h"
#include "HiggsTauTauRun2/interface/PreselectionFilter.h"

// Generic modules
#include "Modules/interface/LumiMask.h"
#include "Modules/interface/SimpleFilter.h"
#include "Modules/interface/CompositeProducer.h"
#include "Modules/interface/CopyCollection.h"
#include "Modules/interface/OneCollCompositeProducer.h"
#include "Modules/interface/OverlapFilter.h"
#include "Modules/interface/EnergyShifter.h"
#include "Modules/interface/EnergyResShifter.h"
#include "Modules/interface/PileupWeight.h"
#include "Modules/interface/CheckEvents.h"
#include "Modules/interface/GenericModule.h"


namespace ic {

HTTSequence::HTTSequence(std::string& chan, std::string postf, Json::Value const& json, Json::Value const& unmod_json) {
  if(json["output_name"].asString()!=""){output_name=json["output_name"].asString();} else{std::cout<<"ERROR: output_name not set"<<std::endl; exit(1);};
  do_recoil = json["do_recoil"].asBool() && ((output_name.find("DY")!=output_name.npos && output_name.find("JetsToLL")!=output_name.npos) || output_name.find("HToTauTau")!=output_name.npos || output_name.find("WJetsToLNu") != output_name.npos || output_name.find("W1JetsToLNu") != output_name.npos || output_name.find("W2JetsToLNu")!=output_name.npos || output_name.find("W3JetsToLNu")!=output_name.npos || output_name.find("W4JetsToLNu")!=output_name.npos || output_name.find("WG")!=output_name.npos || output_name.find("EWKW")!=output_name.npos || output_name.find("EWKZ")!=output_name.npos || output_name.find("VBFH")!=output_name.npos || output_name.find("GluGluH")!=output_name.npos || output_name.find("ZHiggs")!=output_name.npos || output_name.find("WHiggs")!=output_name.npos || output_name.find("JJH")!=output_name.npos);
  addit_output_folder=json["baseline"]["addit_output_folder"].asString();
  new_svfit_mode = json["new_svfit_mode"].asUInt();
  if(new_svfit_mode > 0){
    if(json["svfit_folder"].asString()!="") {svfit_folder = json["svfit_folder"].asString();} else {std::cout<<"ERROR: svfit_folder not set"<<std::endl; exit(1);};
    if(addit_output_folder.find("BFAKE_") != output_name.npos || addit_output_folder.find("BTAG_") != output_name.npos) svfit_folder=svfit_folder+"/";
    else if(json["baseline"]["jer_mode"].asUInt() > 0 && do_recoil) svfit_folder=svfit_folder+"/";
    else if(json["baseline"]["jes_mode"].asUInt() > 0 && do_recoil) svfit_folder=svfit_folder+"/";
    else svfit_folder=svfit_folder+"/"+addit_output_folder+"/";
  }
  svfit_override = json["svfit_override"].asString();
  if(json["output_folder"].asString()!=""){output_folder=json["output_folder"].asString();} else{std::cout<<"ERROR: output_folder not set"<<std::endl; exit(1);};
  addit_output_folder=json["baseline"]["addit_output_folder"].asString();
  output_folder=output_folder+"/"+addit_output_folder+"/";
  lumimask_output_name=output_name + "_" +chan + "_" + postf; 
  output_name=output_name + "_" + chan + "_" + postf + ".root";
  special_mode=json["special_mode"].asUInt();
  if (special_mode > 0) output_name = "Special_"+boost::lexical_cast<std::string>(special_mode)+"_" + output_name;
  if(!json["make_sync_ntuple"].asBool() && !json["make_mva_ntuple"].asBool()) {
      fs = std::make_shared<fwlite::TFileService>(
       (output_folder+output_name).c_str());
  } else {
      fs = NULL;
  }


  js = json;
  unmod_js = unmod_json;
  channel_str = chan;
  jes_mode=json["baseline"]["jes_mode"].asUInt();
  jer_mode=json["baseline"]["jer_mode"].asUInt();
  jes_corr_mode=json["baseline"]["jes_corr_mode"].asUInt();
  metscale_mode=json["baseline"]["metscale_mode"].asUInt();
  metres_mode=json["baseline"]["metres_mode"].asUInt();
  metcl_mode=json["baseline"]["metcl_mode"].asUInt();
  metuncl_mode=json["baseline"]["metuncl_mode"].asUInt();
  do_reshape=json["baseline"]["do_reshape"].asBool(); 
  use_deep_csv=json["baseline"]["use_deep_csv"].asBool();
  use_deep_jet=json["baseline"]["use_deep_jet"].asBool(); 
  wp_to_check=json["baseline"]["wp_to_check"].asString(); 
  btag_mode=json["baseline"]["btag_mode"].asUInt();
  bfake_mode=json["baseline"]["bfake_mode"].asUInt();
  if(js["mc"].asString()=="" && era_str != "data_2017"){std::cout<<"ERROR: MC type not set"<<std::endl; exit(1);}
  mc_str = js["mc"].asString();
  if(js["era"].asString()=="") {std::cout<<"ERROR: era not set"<<std::endl; exit(1);}
  era_str = js["era"].asString();
  if(js["strategy"].asString()==""){std::cout<<"ERROR: strategy not set"<<std::endl; exit(1);}
  strategy_str = js["strategy"].asString();
  ic::strategy strategy_type  = String2Strategy(strategy_str);
  ic::era  era_type  = String2Era(era_str);
  veto_elec_pt = 10;
  veto_elec_eta = 2.5;
  veto_elec_dxy = 0.045;
  veto_elec_dz = 0.2;
  veto_muon_pt = 10;
  veto_muon_eta = 2.4;
  veto_muon_dxy = 0.045;
  veto_muon_dz = 0.2;
  veto_dielec_pt = 15;
  veto_dielec_eta = 2.5;
  veto_dielec_dxy = 0.045;
  veto_dielec_dz = 0.2;
  veto_dimuon_pt = 15;
  veto_dimuon_eta = 2.4;
  veto_dimuon_dxy = 0.045;
  veto_dimuon_dz = 0.2;
  lead_tau_pt=0;
  lead_min_taus=0;

  elec_dz = 0.2;
  elec_dxy = 0.045;
  muon_dxy = 0.045;
  muon_dz = 0.2;
  tau_dz = 0.2;
  pair_dr = 0.5;

  if(channel_str == "em"){
    pair_dr = 0.3;
    elec_pt = 15.;
    elec_eta = 2.5;
    muon_pt = 15;
    muon_eta = 2.4;
  }
  if (channel_str == "et"){
    min_taus = 1;
    tau_pt  = 20;
    tau_eta = 2.3;
    if(era_type == era::data_2017 || era_type == era::data_2017UL) elec_pt = 28.;
    if(era_type == era::data_2018 || era_type == era::data_2018UL) elec_pt = 33.;
    if(era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) elec_pt = 33.;
    else elec_pt = 25.;
    elec_eta = 2.5;
  } 
  if (channel_str == "mt"){
   //muon_eta = 2.1; 
   tau_pt = 20;
   tau_eta = 2.3;
   min_taus = 1;
   //if(era_type == era::data_2016 || era_type == era::data_2016UL_preVFP || era_type == era::data_2016UL_postVFP || json["make_sync_ntuple"].asBool()) muon_pt = 20.0; 
   //else  muon_pt = 21.0;
    if(era_type == era::data_2017 || era_type == era::data_2017UL  || era_type == era::data_2018 || era_type == era::data_2018UL) muon_pt = 25.;
    if(era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) muon_pt = 25.;
    else muon_pt = 23.;
    muon_eta = 2.4;
  }
  if (channel_str == "tt"){
   tau_pt=40;
   lead_tau_pt=40;
   if(js["do_vbftrg"].asBool()) {
     tau_pt = 20.;
     lead_tau_pt = 20;
   }
   if(era_type == era::data_2022_preEE || era_type == era::data_2022_preEE){
    tau_pt=35;
    lead_tau_pt=35;
   }
   tau_eta = 2.1;
   min_taus = 2;
   lead_min_taus = 1;
  }
  if(channel_str == "zee"){
    if(era_type == era::data_2017 || era_type == era::data_2017UL) elec_pt = 28.;
    if(era_type == era::data_2018 || era_type == era::data_2018UL) elec_pt = 33.;
    if(era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) elec_pt = 33.;
    else elec_pt = 25.;
    elec_eta = 2.5;
  }
  if(channel_str == "zmm"){
    if(era_type == era::data_2017 || era_type == era::data_2017UL  || era_type == era::data_2018 || era_type == era::data_2018UL) muon_pt = 25.;
    if(era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) muon_pt = 25.;
    else muon_pt = 23.;
    muon_eta = 2.4; 
  }
  if(channel_str == "tpzmm"){
    muon_pt = 10;
    muon_eta = 2.4;
  }
  if(channel_str == "tpzee"){
    elec_pt = 10;
    elec_eta = 2.5;
  }
  if(channel_str == "tpmt"){
    muon_pt = 23;
    muon_eta = 2.1;
    min_taus = 1;
    tau_pt = 10; 
    tau_eta = 2.3;
  }
 do_qcd_scale_wts_=false;
 do_qcd_scale_wts_ = js["do_qcd_scale_wts"].asBool();
 if (output_name.find("SUSYGluGluToBBHToTauTau_M") != output_name.npos && output_name.find("NLO") != output_name.npos) do_qcd_scale_wts_=true;
  
 is_data      = json["is_data"].asBool();
 is_embedded  = json["is_embedded"].asBool();
 jets_label   = json["jets"].asString();
 met_label    = json["met"].asString();
 if(json["baseline"]["mass_scale_mode"].asBool()==true){
   mass_shift = json["baseline"]["mass_shift"].asDouble();
 } else mass_shift=1.00;
 do_met_filters = json["do_met_filters"].asBool();
 tau_scale_mode = json["baseline"]["tau_scale_mode"].asBool();
 e_scale_mode = json["baseline"]["e_scale_mode"].asBool();
 e_unc_mode = js["baseline"]["e_unc_mode"].asUInt();
 mu_scale_mode = json["baseline"]["mu_scale_mode"].asBool();
 //Need this to correctly set tau /elec ES
 muon_shift = json["baseline"]["muon_es_shift"].asDouble();
 muon_shift_barrel = json["baseline"]["muon_es_shift_barrel"].asDouble();
 muon_shift_nearendcap = json["baseline"]["muon_es_shift_nearendcap"].asDouble();
 muon_shift_farendcap = json["baseline"]["muon_es_shift_farendcap"].asDouble();
 if(channel_str!="em"){
 tau_shift = json["baseline"]["tau_es_shift"].asDouble();
 } else {
 tau_shift = 1.0;//Just to set it to something
 }
 elec_shift_barrel = json["baseline"]["elec_es_shift_barrel"].asDouble();
 elec_shift_endcap = json["baseline"]["elec_es_shift_endcap"].asDouble();
 tau_shift_1prong0pi0 = 1.0;
 tau_shift_1prong1pi0 = 1.0;
 tau_shift_3prong0pi0 = 1.0;
 tau_shift_3prong1pi0 = 1.0;
 fakeE_tau_shift_0pi = 1.0;
 fakeE_tau_shift_1pi = 1.0;
 fakeE_tau_shift_0pi_endcap = 1.0;
 fakeE_tau_shift_1pi_endcap = 1.0;
 fakeMu_tau_shift_0pi = 1.0;
 fakeMu_tau_shift_1pi = 1.0;
 fakeJ_tau_shift =1.;

 muon_res_shift = 0.;
 elec_res_shift = 0.;
 tau_res_shift = 0.;

 if(strategy_type==strategy::mssmsummer16 || strategy_type == strategy::smsummer16 || strategy_type == strategy::cpsummer16 ||strategy_type == strategy::legacy16 || strategy_type == strategy::cpdecays16 || strategy_type == strategy::cpsummer17 || strategy_type == strategy::cpdecays17 || strategy_type == strategy::cpdecays18){
   fakeE_tau_shift_0pi = json["baseline"]["efaketau_0pi_es_shift"].asDouble();
   fakeE_tau_shift_1pi = json["baseline"]["efaketau_1pi_es_shift"].asDouble();
   fakeE_tau_shift_0pi_endcap = json["baseline"]["efaketau_0pi_es_shift_endcap"].asDouble();
   fakeE_tau_shift_1pi_endcap = json["baseline"]["efaketau_1pi_es_shift_endcap"].asDouble();
   if(!is_embedded){
     tau_shift_1prong0pi0 = json["baseline"]["tau_1prong0pi0_es_shift"].asDouble();
     tau_shift_1prong1pi0 = json["baseline"]["tau_1prong1pi0_es_shift"].asDouble();
     tau_shift_3prong0pi0 = json["baseline"]["tau_3prong0pi0_es_shift"].asDouble();
     tau_shift_3prong1pi0 = json["baseline"]["tau_3prong1pi0_es_shift"].asDouble();
   } else {
     tau_shift_1prong0pi0 = json["baseline"]["embedtau_1prong0pi0_es_shift"].asDouble();
     tau_shift_1prong1pi0 = json["baseline"]["embedtau_1prong1pi0_es_shift"].asDouble();
     tau_shift_3prong0pi0 = json["baseline"]["embedtau_3prong0pi0_es_shift"].asDouble();
     tau_shift_3prong1pi0 = json["baseline"]["embedtau_3prong1pi0_es_shift"].asDouble();
   }
   muon_res_shift = json["baseline"]["muon_res_shift"].asDouble();
   elec_res_shift = json["baseline"]["elec_res_shift"].asDouble();
   tau_res_shift = json["baseline"]["tau_res_shift"].asDouble();
 }
 if(strategy_type == strategy::smsummer16 || strategy_type == strategy::cpsummer16 || strategy_type == strategy::legacy16 || strategy_type == strategy::cpdecays16 || strategy_type == strategy::cpsummer17 || strategy_type == strategy::cpdecays17 || strategy_type == strategy::cpdecays18){
   fakeMu_tau_shift_0pi = json["baseline"]["mufaketau_0pi_es_shift"].asDouble();
   fakeMu_tau_shift_1pi = json["baseline"]["mufaketau_1pi_es_shift"].asDouble();

   if(json["baseline"]["fakeJ_tau_shift"].asString()!="") fakeJ_tau_shift = json["baseline"]["fakeJ_tau_shift"].asDouble();
 }

tau_shift_func_1prong0pi0 = "";
tau_shift_func_1prong1pi0 = "";
tau_shift_func_3prong0pi0 = "";
tau_shift_func_3prong1pi0 = "";
if(!is_embedded){
  if(json["baseline"]["tau_shift_func_1prong0pi0"].asString()!="") tau_shift_func_1prong0pi0=json["baseline"]["tau_shift_func_1prong0pi0"].asString();
  if(json["baseline"]["tau_shift_func_1prong1pi0"].asString()!="") tau_shift_func_1prong1pi0=json["baseline"]["tau_shift_func_1prong1pi0"].asString();
  if(json["baseline"]["tau_shift_func_3prong0pi0"].asString()!="") tau_shift_func_3prong0pi0=json["baseline"]["tau_shift_func_3prong0pi0"].asString();
  if(json["baseline"]["tau_shift_func_3prong1pi0"].asString()!="") tau_shift_func_3prong1pi0=json["baseline"]["tau_shift_func_3prong1pi0"].asString();
}

alt_jes_input_set = json["baseline"]["jes_input_set"].asString();

}

HTTSequence::~HTTSequence() {}


void HTTSequence::BuildSequence(){
  using ROOT::Math::VectorUtil::DeltaR;

  // Set global parameters that get used in multiple places
  ic::channel channel         = String2Channel(channel_str);
  ic::mc mc_type              = String2MC(mc_str);
  ic::era era_type            = String2Era(era_str);
  ic::strategy strategy_type  = String2Strategy(strategy_str);

  std::cout << "-------------------------------------" << std::endl;
  std::cout << "HiggsToTauTau Analysis" << std::endl;
  std::cout << "-------------------------------------" << std::endl;      std::string param_fmt = "%-25s %-40s\n";
  std::cout << boost::format(param_fmt) % "output" % (output_folder+output_name);
  std::cout << boost::format(param_fmt) % "strategy" % Strategy2String(strategy_type);
  std::cout << boost::format(param_fmt) % "era" % Era2String(era_type);
  std::cout << boost::format(param_fmt) % "mc" % MC2String(mc_type);
  std::cout << boost::format(param_fmt) % "channel" % channel_str;
  std::cout << boost::format(param_fmt) % "is_data" % is_data;
  std::cout << boost::format(param_fmt) % "is_embedded" % is_embedded;
  std::cout << boost::format(param_fmt) % "tau_scale_mode" % tau_scale_mode;
  std::cout << boost::format(param_fmt) % "new_svfit_mode" % new_svfit_mode;
  if (new_svfit_mode > 0) {
    std::cout << boost::format(param_fmt) % "svfit_folder" % svfit_folder;
    std::cout << boost::format(param_fmt) % "svfit_override" % svfit_override;
  }
  std::cout << boost::format(param_fmt) % "make_sync_ntuple" % js["make_sync_ntuple"].asBool();
  std::cout << boost::format(param_fmt) % "save_output_jsons" % js["save_output_jsons"].asBool();
  std::cout << boost::format(param_fmt) % "use_deep_csv" % use_deep_csv;
  std::cout << "** Kinematics **" << std::endl;
  std::cout << boost::format(param_fmt) % "elec_pt" % elec_pt;
  std::cout << boost::format(param_fmt) % "elec_eta" % elec_eta;
  std::cout << boost::format(param_fmt) % "elec_dxy" % elec_dxy;
  std::cout << boost::format(param_fmt) % "elec_dz" % elec_dz;
  std::cout << boost::format(param_fmt) % "muon_pt" % muon_pt;
  std::cout << boost::format(param_fmt) % "muon_eta" % muon_eta;
  std::cout << boost::format(param_fmt) % "muon_dxy" % muon_dxy;
  std::cout << boost::format(param_fmt) % "muon_dz" % muon_dz;
  std::cout << boost::format(param_fmt) % "tau_pt" % tau_pt;
  std::cout << boost::format(param_fmt) % "lead_tau_pt" % lead_tau_pt;
  std::cout << boost::format(param_fmt) % "tau_eta" % tau_eta;
  std::cout << boost::format(param_fmt) % "tau_dz" % tau_dz;
  std::cout << boost::format(param_fmt) % "veto_elec_pt" % veto_elec_pt;
  std::cout << boost::format(param_fmt) % "veto_elec_eta" % veto_elec_eta;
  std::cout << boost::format(param_fmt) % "veto_elec_dxy" % veto_elec_dxy;
  std::cout << boost::format(param_fmt) % "veto_elec_dz" % veto_elec_dz;
  std::cout << boost::format(param_fmt) % "veto_dielec_pt" % veto_dielec_pt;
  std::cout << boost::format(param_fmt) % "veto_dielec_eta" % veto_dielec_eta;
  std::cout << boost::format(param_fmt) % "veto_dielec_dxy" % veto_dielec_dxy;
  std::cout << boost::format(param_fmt) % "veto_dielec_dz" % veto_dielec_dz;
  std::cout << boost::format(param_fmt) % "veto_muon_pt" % veto_muon_pt;
  std::cout << boost::format(param_fmt) % "veto_muon_eta" % veto_muon_eta;
  std::cout << boost::format(param_fmt) % "veto_muon_dxy" % veto_muon_dxy;
  std::cout << boost::format(param_fmt) % "veto_muon_dz" % veto_muon_dz;
  std::cout << boost::format(param_fmt) % "veto_dimuon_pt" % veto_dimuon_pt;
  std::cout << boost::format(param_fmt) % "veto_dimuon_eta" % veto_dimuon_eta;
  std::cout << boost::format(param_fmt) % "veto_dimuon_dxy" % veto_dimuon_dxy;
  std::cout << boost::format(param_fmt) % "veto_dimuon_dz" % veto_dimuon_dz;

  auto eventChecker = CheckEvents("EventChecker").set_skip_events(true);
  std::vector<int> to_check =
  {
  };   
  if(strcmp((js["event_check_file"].asString()).c_str(),"")!=0){
    std::ifstream file((js["event_check_file"].asString()).c_str());
    if (!file.is_open()) {
      std::cerr << "Warning: File " << js["event_check_file"].asString() << " cannot be opened." << std::endl;
    } 
    int nums;
    std::string line;
    while(std::getline(file, line)){
      nums = std::atoi(line.c_str());
      to_check.push_back(nums);
      std::cout << nums << std::endl;
    }

    file.close();
  }

  // Defining good-lumi jsons
  std::string data_json = "";
  if ((era_type == era::data_2016 || era_type == era::data_2016UL_preVFP || era_type == era::data_2016UL_postVFP))
      data_json = "input/json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt";
  if ((era_type == era::data_2017 || era_type == era::data_2017UL))
      data_json = "input/json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt";
  if ((era_type == era::data_2018 || era_type == era::data_2018UL))
      data_json = "input/json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt";
  if ((era_type == era::data_2022_preEE || era_type == era::data_2022_postEE))
      //data_json = "input/json/Cert_Collisions2022_355100_362760_Golden_lumicorr.json"; //359571 and 359661 removed https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3  
      data_json = "input/json/Cert_Collisions2022_355100_362760_Golden.json"; 

 if(js["get_effective"].asBool() && js["make_sync_ntuple"].asBool()){
   std::cerr<< "Error: cannot run effective number of event module in make_syncntuple mode"<<std::endl;
   throw;
 }
 if(js["get_effective"].asBool()){
  BuildModule(EffectiveEvents("EffectiveEvents")
    .set_fs(fs.get())
    .set_do_qcd_scale_wts(do_qcd_scale_wts_));
  }else if(is_data && js["lumi_mask_only"].asBool()){
   //LumiMask module  just for the mode in which we only write out the 
   //output jsons
   LumiMask lumiMask = LumiMask("LumiMask")
     .set_produce_output_jsons(lumimask_output_name.c_str())
     .set_input_file(data_json);
 
    BuildModule(lumiMask);
  
  } else{

  HTTPrint httPrint("HTTPrint");
  if(era_type==era::data_2015 || era_type==era::data_2016  
          || era_type == era::data_2017 || era_type == era::data_2018 || era_type == era::data_2016UL_preVFP || era_type == era::data_2016UL_postVFP || era_type == era::data_2017UL || era_type == era::data_2018UL || era_type == era::data_2022_preEE || era_type == era::data_2022_postEE){
    httPrint.set_muon_label(js["muons"].asString());
    httPrint.set_jet_label(jets_label);
  }
  for (auto ch : to_check) {
   eventChecker.CheckEvent(ch);
   httPrint.PrintEvent(ch);
  }
  httPrint.set_skip_events(false);
  if (to_check.size() > 0){
  BuildModule(eventChecker);
  BuildModule(httPrint);  
}

if(!is_data && js["do_gen_analysis"].asBool()){

  // do pu weights for MC samples

  TH1D d_pu = GetFromTFile<TH1D>(js["data_pu_file"].asString(), "/", "pileup");
  TH1D m_pu = GetFromTFile<TH1D>(js["mc_pu_file"].asString(), "/", "pileup");
  //if (js["do_pu_wt"].asBool()&&!is_data&&!is_embedded) {
  //  BuildModule( PileupWeight("PileupWeight")
  //      .set_data(new TH1D(d_pu)).set_mc(new TH1D(m_pu)));
  //}

  // do stitching for dy samples:
  HTTStitching httStitching = HTTStitching("HTTStitching")
      .set_era(era_type)
      .set_fs(fs.get());
  if ((output_name.find("DY") != output_name.npos && output_name.find("JetsToLL-LO") != output_name.npos && !(output_name.find("JetsToLL-LO-10-50") != output_name.npos))){
    httStitching.set_do_dy_soup(true);

    if(era_type == era::data_2016) {
      httStitching.SetDYInputCrossSections(4954, 1012.5, 332.8, 101.8,54.8);
      httStitching.SetDYInputYields(49748967+90972768, 63730337, 19879279, 5857441, 4197868);
    }
    if(era_type == era::data_2016UL_preVFP){
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118);
      httStitching.SetDYInputYields(95170542, 31663215, 12696373, 10285720, 4655680);
    }
    if(era_type == era::data_2016UL_postVFP){
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118);
      httStitching.SetDYInputYields(82444586, 31570459, 14573456, 9553830, 4450325);
    }
    if(era_type == era::data_2017) {
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118); //Target fractions are xs_n-jet/xs_inclusive
      httStitching.SetDYInputYields(48590164+49031214,42205667+33563494,88795+9871382,1147725+5740168,4317756);
    }
    if(era_type == era::data_2017UL) {
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118); //Target fractions are xs_n-jet/xs_inclusive
      httStitching.SetDYInputYields(97237260+101893848, 66063783, 27099640, 20165687, 11267712);
    }
    if(era_type == era::data_2018) {
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118); //Target fractions are xs_n-jet/xs_inclusive
      httStitching.SetDYInputYields(99536185, 66994632, 19918141, 5554767, 2812482);
    }
    if(era_type == era::data_2018UL) {
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118); //Target fractions are xs_n-jet/xs_inclusive
      httStitching.SetDYInputYields(96055748+101373493, 60368985, 27494377, 20466041, 8885353);
    }
    if(era_type == era::data_2022_preEE ) {
      httStitching.SetDYInputCrossSections(5455.0,978.3,315.1,93.7,45.4);
      httStitching.SetDYInputYields(74978096,0,14777357,8699442,3267376);
    }
    if(era_type == era::data_2022_postEE) {
      httStitching.SetDYInputCrossSections(5455.0,978.3,315.1,93.7,45.4);
      httStitching.SetDYInputYields(244244756,50224518,50586120,28833220,9405785);
    }  


  }
  BuildModule(httStitching);

  std::string mass_str = output_name;
  mass_str.erase(0, mass_str.find("_M-")+3);
  mass_str.erase(mass_str.find("_"),mass_str.length()-mass_str.find("_"));
 
  TH2D deepjet_bbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2018_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_b");
 
  BuildModule(HTTGenAnalysis("HTTGenAnalysis")
    .set_fs(fs.get())
    .set_channel_str(channel_str)
    .set_min_jet_pt(30.)
    .set_max_jet_eta(4.7)
    .set_min_e_pt(0)
    .set_min_mu_pt(0)
    .set_min_tau1_pt(0)
    .set_min_tau2_pt(0)
    .set_max_e_eta(1000)
    .set_max_mu_eta(1000)
    .set_max_tau_eta(1000)
    .set_do_theory_uncert(true)
    .set_mssm_mass(mass_str)
    .set_bbtag_eff((TH2F*)(new TH2D(deepjet_bbtag_eff)))
    .set_era(era_type)
  );
  return;  
}

 //Build LumiMask module here - json file definition eralier in the code to be able to 
 //run *just* this module
 if((is_data|| is_embedded) && strategy_type!=strategy::phys14){
   LumiMask lumiMask = LumiMask("LumiMask")
     .set_produce_output_jsons("")
     .set_input_file(data_json);
 
 if(js["save_output_jsons"].asBool()){
  lumiMask.set_produce_output_jsons(lumimask_output_name.c_str());
   }
  BuildModule(lumiMask);
 }

if((strategy_type == strategy::fall15 || strategy_type ==strategy::mssmsummer16 || strategy_type == strategy::smsummer16 || strategy_type == strategy::cpsummer16 || strategy_type == strategy::legacy16 || strategy_type == strategy::cpdecays16 || strategy_type == strategy::cpsummer17 || strategy_type == strategy::cpdecays17 || strategy_type == strategy::cpdecays18) && output_name.find("WGToLNuG")!=output_name.npos){
  SimpleCounter<GenParticle> wgammaStarFilter = SimpleCounter<GenParticle>("WgammaStarSelector")
    .set_input_label("genParticles")
    .set_predicate(
      (bind(&GenParticle::status,_1) == 44) &&
      (bind(&GenParticle::pdgid, _1) == 22))
    .set_min(0).set_max(0);
 
  BuildModule(wgammaStarFilter);
}

if (channel == channel::et) BuildETPairs();
if (channel == channel::mt) BuildMTPairs();
if (channel == channel::em) BuildEMPairs();
if (channel == channel::tt) BuildTTPairs();
if (channel == channel::zee) BuildZEEPairs();
if (channel == channel::zmm) BuildZMMPairs();
if (channel == channel::tpzee) BuildTPZEEPairs();
if (channel == channel::tpzmm) BuildTPZMMPairs();

// Pair DeltaR filtering
BuildModule(SimpleFilter<CompositeCandidate>("PairFilter")
      .set_input_label("ditau").set_min(1)
      .set_predicate([=](CompositeCandidate const* c) {
        return DeltaR(c->at(0)->vector(), c->at(1)->vector())
            > pair_dr;
      }));

bool usePFMET = false;
usePFMET = js["usePFMET"].asBool();

HTTPairSelector httPairSelector = HTTPairSelector("HTTPairSelector")
  .set_channel(channel)
  .set_fs(fs.get())
  .set_pair_label("ditau")
  .set_met_label(met_label)
  .set_era(era_type)
  .set_strategy(strategy_type)
  .set_mva_met_from_vector(mva_met_mode==1)
  .set_faked_tau_selector(faked_tau_selector)
  .set_hadronic_tau_selector(hadronic_tau_selector)
  .set_ztt_mode(ztautau_mode)
  .set_mc(mc_type)
  .set_gen_taus_label(is_embedded ? "genParticlesEmbedded" : "genParticlesTaus")
  .set_scale_met_for_tau((tau_scale_mode > 0 || (moriond_tau_scale && (is_embedded || !is_data) )   ))
  .set_tau_scale(tau_shift)
  .set_use_most_isolated((strategy_type != strategy::paper2013) && (!(channel == channel::zee || channel == channel::zmm || channel == channel::tpzmm || channel == channel::tpzee)))
  .set_use_os_preference((strategy_type == strategy::paper2013) || (channel == channel::zee || channel == channel::zmm || channel == channel::tpzmm || channel == channel::tpzee))
  .set_allowed_tau_modes(allowed_tau_modes)
  .set_metuncl_mode(metuncl_mode)
  .set_metcl_mode(metcl_mode)
  .set_shift_jes(!do_recoil)
  .set_usePFMET(usePFMET);
    
httPairSelector.set_gen_taus_label("genParticles");

BuildModule(httPairSelector);
 


// Lepton Vetoes
if (js["baseline"]["di_elec_veto"].asBool()) BuildDiElecVeto();
if (js["baseline"]["di_muon_veto"].asBool()) BuildDiMuonVeto();
if (js["baseline"]["extra_elec_veto"].asBool()) BuildExtraElecVeto();
if (js["baseline"]["extra_muon_veto"].asBool()) BuildExtraMuonVeto();


if(js["do_preselection"].asBool() && channel != channel::tpzee && channel != channel::tpzmm){
  BuildModule(PreselectionFilter("PreselectionFilter")
   .set_era(era_type)
   .set_channel(channel)
   .set_do_preselection(true)
   .set_dilepton_label("ditau"));
}


// Pileup Weighting
TH1D d_pu = GetFromTFile<TH1D>(js["data_pu_file"].asString(), "/", "pileup");
TH1D m_pu = GetFromTFile<TH1D>(js["mc_pu_file"].asString(), "/", "pileup");
if (js["do_pu_wt"].asBool()&&!is_data&&!is_embedded) {
  BuildModule( PileupWeight("PileupWeight")
      .set_data(new TH1D(d_pu)).set_mc(new TH1D(m_pu)));
}


if(do_met_filters){
  BuildModule(GenericModule("MetFiltersRecoEffect")
    .set_function([=](ic::TreeEvent *event){
      EventInfo *eventInfo = event->GetPtr<EventInfo>("eventInfo");
      std::vector<std::string> met_filters = {"Flag_eeBadScFilter"};
      if (era_type == era::data_2016)
      {
        met_filters.push_back("Flag_HBHENoiseFilter");
	met_filters.push_back("Flag_HBHENoiseIsoFilter");
	met_filters.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
	met_filters.push_back("Flag_goodVertices");
	met_filters.push_back("badMuonFilter");
	met_filters.push_back("Flag_globalTightHalo2016Filter");
      } 
      if (era_type == era::data_2017 || era_type == era::data_2018)
      {
        met_filters.push_back("Flag_HBHENoiseFilter");
        met_filters.push_back("Flag_HBHENoiseIsoFilter");
        met_filters.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
        met_filters.push_back("Flag_goodVertices");
        met_filters.push_back("badMuonFilter");
        met_filters.push_back("Flag_globalSuperTightHalo2016Filter");
	met_filters.push_back("ecalBadCalibReducedMINIAODFilter");
      }
      if (era_type == era::data_2016UL_preVFP || era_type == era::data_2016UL_postVFP) 
      {
        met_filters.push_back("Flag_goodVertices");
	met_filters.push_back("Flag_globalSuperTightHalo2016Filter");
	met_filters.push_back("Flag_HBHENoiseFilter");
	met_filters.push_back("Flag_HBHENoiseIsoFilter");
	met_filters.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
	met_filters.push_back("Flag_BadPFMuonFilter");
	//disabled (missing filter) 
	//met_filters.push_back("Flag_BadPFMuonDzFilter");
      }
      if (era_type == era::data_2017UL || era_type == era::data_2018UL) 
      {
        met_filters.push_back("Flag_goodVertices");
        met_filters.push_back("Flag_globalSuperTightHalo2016Filter");
        met_filters.push_back("Flag_HBHENoiseFilter");
        met_filters.push_back("Flag_HBHENoiseIsoFilter");
        met_filters.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
        met_filters.push_back("Flag_BadPFMuonFilter");
        //disabled (missing filter)
	//met_filters.push_back("Flag_BadPFMuonDzFilter");
	met_filters.push_back("Flag_ecalBadCalibFilter");
      }
      if (era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) 
      {
        met_filters.push_back("Flag_goodVertices");
        met_filters.push_back("Flag_globalSuperTightHalo2016Filter");
        met_filters.push_back("Flag_HBHENoiseFilter");
        met_filters.push_back("Flag_HBHENoiseIsoFilter");
        met_filters.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
        met_filters.push_back("Flag_BadPFMuonFilter");
        //disabled (missing filter)
	//met_filters.push_back("Flag_BadPFMuonDzFilter");
	met_filters.push_back("Flag_ecalBadCalibFilter");
      }      
      bool pass_filters = true;
      for(unsigned i=0;i<met_filters.size();++i){
       pass_filters = pass_filters&& eventInfo->filter_result(met_filters.at(i));
      }
      if(js["make_sync_ntuple"].asBool() and false) {
        event->Add("flagMETFilter", pass_filters);  
        return false;    
      }
      return !pass_filters;
    }));
}

//if((era_type == era::data_2016 || era_type == era::data_2017) && !is_embedded && !is_data) { // we now use the value stored on the miniAOD
//  TH2F prefire_hist; 
//  if(era_type == era::data_2016) prefire_hist = GetFromTFile<TH2F>("input/prefire/L1prefiring_jetpt_2016BtoH.root","/","L1prefiring_jetpt_2016BtoH");
//  else prefire_hist = GetFromTFile<TH2F>("input/prefire/L1prefiring_jetpt_2017BtoF.root","/","L1prefiring_jetpt_2017BtoF");
//  BuildModule(HTTPreFireWeight<PFJet>("HTTPreFireWeight")
//    .set_prefire_hist(new TH2F(prefire_hist)));
//}

BuildModule(CopyCollection<PFJet>("CopyFilteredJets",jets_label,jets_label+"UnFiltered"));

SimpleFilter<PFJet> jetIDFilter = SimpleFilter<PFJet>("JetIDFilter")
.set_input_label(jets_label);
if (era_type == era::data_2016 || era_type == era::data_2016UL_preVFP || era_type == era::data_2016UL_postVFP) {
  jetIDFilter.set_predicate(bind(PFJetID2016, _1));
} else if (era_type == era::data_2017 || era_type == era::data_2017UL) {
  jetIDFilter.set_predicate(bind(PFJetID2017, _1));
} else if (era_type == era::data_2018 || era_type == era::data_2018UL) {
  jetIDFilter.set_predicate(bind(PFJetID2018, _1));
}
else if (era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) {
  jetIDFilter.set_predicate(bind(PFJetID2018, _1));
}
BuildModule(jetIDFilter);

// Apply loose PUJID universally
//BuildModule(SimpleFilter<PFJet>("JetPUIDFilter")
//  .set_input_label(jets_label)
//  .set_predicate([=](PFJet const* jet) {
//    return  PileupJetID(jet, pu_id_training, false, true);
//  })
//);

if (era_type == era::data_2017 || era_type == era::data_2017UL) {
  BuildModule(SimpleFilter<PFJet>("JetEENoiseVetoFilter")
    .set_input_label(jets_label)
    .set_predicate([=](PFJet const* jet) {
      return  (jet->pt() * (jet->uncorrected_energy()/jet->energy())) > 50 ||
        fabs(jet->eta()) > 3.139 ||
        fabs(jet->eta()) < 2.65 ;
    })
  );
}

// nominal JER
if (!is_data && !is_embedded) {
   BuildModule(JetEnergyResolution<PFJet>("JetEnergyResolution")
     .set_input_label(jets_label)
     .set_jer_shift_mode(jer_mode)
     .set_EENoiseFix(era_type == era::data_2017 || era_type == era::data_2017UL)
     .set_shift_met(false)
   );
}

BuildModule(OverlapFilter<PFJet, CompositeCandidate>("JetLeptonOverlapFilter")
  .set_input_label(jets_label)
  .set_reference_label("ditau")
  .set_min_dr(0.5));

if(!is_data) {
  bool do_ngenjets = (output_name.find("GluGluHToTauTauUncorrelatedDecay") != output_name.npos || output_name.find("JJH")!=output_name.npos);
  BuildModule(HTTPairGenInfo("HTTPairGenInfo")
    .set_fs(fs.get())
    .set_write_plots(false)
    .set_ditau_label("ditau")
    .set_channel(channel)
    .set_era(era_type)
    .set_is_embedded(is_embedded)
    .set_ngenjets((do_ngenjets||true) && !is_data && !is_embedded));
}

std::string scalefactor_file;
std::string scalefactor_file_ggh;
std::string scalefactor_file_UL;
if(era_type == era::data_2016) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_legacy_2016.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2016_MGggh.root";
}
if(era_type == era::data_2016UL_preVFP) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_legacy_2016.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2016_MGggh.root";
   scalefactor_file_UL = "input/scale_factors/htt_scalefactors_UL_2016preVFP.root";
}
if(era_type == era::data_2016UL_postVFP) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_legacy_2016.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2016_MGggh.root";
   scalefactor_file_UL = "input/scale_factors/htt_scalefactors_UL_2016postVFP.root";
}
if(era_type == era::data_2017) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_legacy_2017.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2017_MGggh.root";
}
if(era_type == era::data_2017UL) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_legacy_2017.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2017_MGggh.root";
   scalefactor_file_UL = "input/scale_factors/htt_scalefactors_UL_2017.root";   
}
if(era_type == era::data_2018) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_legacy_2018.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2017_MGggh.root";
}
if(era_type == era::data_2018UL) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_legacy_2018.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2017_MGggh.root";
   scalefactor_file_UL = "input/scale_factors/htt_scalefactors_UL_2018.root";
}
if(era_type == era::data_2022_preEE) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_2022_preEE.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2017_MGggh.root";
   scalefactor_file_UL = "";
}
if(era_type == era::data_2022_postEE) {
   scalefactor_file = "input/scale_factors/htt_scalefactors_2022_postEE.root";
   scalefactor_file_ggh = "input/ggh_weights/htt_scalefactors_2017_MGggh.root";
   scalefactor_file_UL = "";
}


TH2F embed_pt_weights_ic;

if(era_type == era::data_2016 || era_type == era::data_2016UL_preVFP || era_type == era::data_2016UL_postVFP) embed_pt_weights_ic = GetFromTFile<TH2F>("input/scale_factors/embed_zmm_shifts_v2.root","/","shifts_2016");
if(era_type == era::data_2017 || era_type == era::data_2017UL) embed_pt_weights_ic = GetFromTFile<TH2F>("input/scale_factors/embed_zmm_shifts_v2.root","/","shifts_2017");
if(era_type == era::data_2018 || era_type == era::data_2018UL) embed_pt_weights_ic = GetFromTFile<TH2F>("input/scale_factors/embed_zmm_shifts_v2.root","/","shifts_2018");
if(era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) embed_pt_weights_ic = GetFromTFile<TH2F>("input/scale_factors/embed_zmm_shifts_v2.root","/","shifts_2018");

HTTWeights httWeights = HTTWeights("HTTWeights")   
 .set_channel(channel)
 .set_era(era_type)
 .set_mc(mc_type)
 .set_do_tau_id_weights(!is_data && !(era_type == era::data_2022_preEE || era_type == era::data_2022_postEE)) // tau ID corrections not derived yet for Run-3
 .set_do_em_qcd_weights(true)
 .set_ditau_label("ditau")
 .set_jets_label(jets_label)
 .set_do_quarkmass_higgspt(output_name.find("JJH")!=output_name.npos)
 .set_do_ps_weights(output_name.find("JJH")!=output_name.npos)
 .set_embed_pt_weights_ic_hist(new TH2F(embed_pt_weights_ic))
 .set_do_nnlops_weights(output_name.find("JJH")!=output_name.npos||output_name.find("GluGluHToTauTauUncorrelatedDecay") != output_name.npos||(output_name.find("GluGlu") != output_name.npos && output_name.find("125") != output_name.npos && output_name.find("SUSY") == output_name.npos));
httWeights.set_strategy(strategy_type);
httWeights.set_scalefactor_file_ggh(scalefactor_file_ggh);
httWeights.set_scalefactor_file(scalefactor_file);
httWeights.set_scalefactor_file_UL(scalefactor_file_UL);
httWeights.set_is_embedded(is_embedded);
if (!is_data ) {
  httWeights.set_do_trg_weights(true).set_trg_applied_in_mc(js["trg_in_mc"].asBool()).set_do_idiso_weights(true);
  if(era_type != era::data_2022_preEE && era_type != era::data_2022_postEE){ //these weights don't exist yet for Run-3
    httWeights.set_do_tau_id_sf(true);
    if(channel == channel::et || channel == channel::mt || channel==channel::tt) {
      httWeights.set_do_etau_fakerate(true);
      httWeights.set_do_mtau_fakerate(true);
    }
    if(channel == channel::et || channel==channel::em || channel==channel::mt || channel==channel::zmm || channel==channel::zee) httWeights.set_do_tracking_eff(true);
  }
}
if ((output_name.find("DY") != output_name.npos && output_name.find("JetsToLL-LO") != output_name.npos && !(output_name.find("JetsToLL-LO-10-50") != output_name.npos)) || output_name.find("DYto2L-4Jets") != output_name.npos){
  httWeights.set_do_zpt_weight(true&&channel!=channel::tpzee&&channel!=channel::tpzmm&&channel!=channel::tpmt&&channel != channel::tpem);
  httWeights.set_do_zpt_weight_NLO(true&&channel!=channel::tpzee&&channel!=channel::tpzmm&&channel!=channel::tpmt&&channel != channel::tpem);
}

if ((output_name.find("DYJetsToLL-NLO") != output_name.npos || output_name.find("DYJetstoLL-NLO") != output_name.npos || output_name.find("DYJetsToLL_0J-NLO") != output_name.npos|| output_name.find("DYJetsToLL_1J-NLO") != output_name.npos || output_name.find("DYJetsToLL_2J-NLO") != output_name.npos )){
  httWeights.set_do_zpt_weight_NLO(true&&channel!=channel::tpzee&&channel!=channel::tpzmm&&channel!=channel::tpmt&&channel != channel::tpem);
  httWeights.set_do_zpt_weight(true&&channel!=channel::tpzee&&channel!=channel::tpzmm&&channel!=channel::tpmt&&channel != channel::tpem);
}


if (output_name.find("TT") != output_name.npos) httWeights.set_do_topquark_weights(true);
if (output_name.find("WJetsToLNu-LO") != output_name.npos || output_name.find("W1JetsToLNu-LO") != output_name.npos || output_name.find("W2JetsToLNu-LO") != output_name.npos ||
     output_name.find("W3JetsToLNu-LO") != output_name.npos || output_name.find("W4JetsToLNu-LO") != output_name.npos){
  httWeights.set_do_tau_fake_weights(true);
}
if(output_name.find("SUSYGluGluToHToTauTau_M") != output_name.npos){
  httWeights.set_do_mssm_higgspt(true); 
  if(era_type == era::data_2016 || era_type == era::data_2016UL_preVFP || era_type == era::data_2016UL_postVFP) httWeights.set_mssm_higgspt_file("input/mssm_higgspt/higgs_pt_2016_v2.root");
  else httWeights.set_mssm_higgspt_file("input/mssm_higgspt/higgs_pt_v2.root");
  std::string mass_str = output_name;
  mass_str.erase(0, mass_str.find("_M-")+3);
  mass_str.erase(mass_str.find("_"),mass_str.length()-output_name.find("_"));
  httWeights.set_mssm_mass(mass_str);
}

if(channel!=channel::tpzee&&channel!=channel::tpzmm) BuildModule(httWeights);

if(channel!=channel::tpzee&&channel!=channel::tpzmm){
  HTTStitching httStitching = HTTStitching("HTTStitching")  
      .set_era(era_type)
      .set_fs(fs.get());
  if (output_name.find("WJetsToLNu-LO") != output_name.npos || output_name.find("W1JetsToLNu-LO") != output_name.npos || output_name.find("W2JetsToLNu-LO") != output_name.npos ||
      output_name.find("W3JetsToLNu-LO") != output_name.npos || output_name.find("W4JetsToLNu-LO") != output_name.npos || output_name.find("WtoLNu-4Jets") != output_name.npos || output_name.find("WJetsToLNu-4Jets") != output_name.npos){
    httStitching.set_do_w_soup(true);
    if(era_type == era::data_2016) {
     httStitching.SetWInputCrossSections(50380,9644.5,3144.5,954.8,485.6); 
     httStitching.SetWInputYields(29514020+57402435, 45283121, 30064264+30374504, 39356879+19798117, 18751462+2073275+9116657);
    }
    if(era_type == era::data_2016UL_preVFP){
     httStitching.SetWInputCrossSections(1.0,0.1522,0.0515,0.0184,0.0103);
     httStitching.SetWInputYields(74676349, 47938929, 26899459, 16517929, 4575554);
    }
    if(era_type == era::data_2016UL_postVFP){
     httStitching.SetWInputCrossSections(1.0,0.1522,0.0515,0.0184,0.0103);
     httStitching.SetWInputYields(82742220, 42019790, 24473569, 17830876, 4816485);
    } 
    if(era_type == era::data_2017) {
      httStitching.SetWInputCrossSections(1.0,0.1522,0.0515,0.0184,0.0103);
      httStitching.SetWInputYields(33043732+44579681,54070624,6511577,19644624,10874617);
    }
    if(era_type == era::data_2017UL){
      httStitching.SetWInputCrossSections(1.0,0.1522,0.0515,0.0184,0.0103);
      httStitching.SetWInputYields(78981243, 46923030, 26941858, 18045982, 9606448);
    }
    if(era_type == era::data_2018) {
      httStitching.SetWInputCrossSections(1.0,0.1522,0.0515,0.0184,0.0103);
      httStitching.SetWInputYields(63389404, 31521085, 16766140, 13815650, 9874750); 
    }
    if(era_type == era::data_2018UL){
      httStitching.SetWInputCrossSections(1.0,0.1522,0.0515,0.0184,0.0103);
      httStitching.SetWInputYields(82360950, 47850048, 27411802, 18297679, 9130068);
    }
    if(era_type == era::data_2022_preEE){
      httStitching.SetWInputCrossSections(55300.,9128.,2922.,861.3,415.4);
      httStitching.SetWInputYields(87890414, 0, 0, 8214152, 1460716);
    }
    if(era_type == era::data_2022_postEE){
      httStitching.SetWInputCrossSections(55300.,9128.,2922.,861.3,415.4);
      httStitching.SetWInputYields(348372911, 42838500, 36384403, 0, 4860871);
    }      
  }
  if ((output_name.find("DY") != output_name.npos && output_name.find("JetsToLL-LO") != output_name.npos && !(output_name.find("JetsToLL-LO-10-50") != output_name.npos))||output_name.find("DYto2L-4Jets_MLL-50") != output_name.npos){
    httStitching.set_do_dy_soup(true);
    if(era_type == era::data_2016) {
      httStitching.SetDYInputCrossSections(4954, 1012.5, 332.8, 101.8,54.8);
      httStitching.SetDYInputYields(49748967+90972768, 63730337, 19879279, 5857441, 4197868);
    }
    if(era_type == era::data_2016UL_preVFP){
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118);
      httStitching.SetDYInputYields(95170542, 31663215, 12696373, 10285720, 4655680);
    }
    if(era_type == era::data_2016UL_postVFP){
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118);
      httStitching.SetDYInputYields(82444586, 31570459, 14573456, 9553830, 4450325);
    }
    if(era_type == era::data_2017) {
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118); //Target fractions are xs_n-jet/xs_inclusive
      httStitching.SetDYInputYields(48590164+49031214,42205667+33563494,88795+9871382,1147725+5740168,4317756);
    }
    if(era_type == era::data_2017UL) {
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118); //Target fractions are xs_n-jet/xs_inclusive
      httStitching.SetDYInputYields(97237260+101893848, 66063783, 27099640, 20165687, 11267712);
    }
    if(era_type == era::data_2018) {
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118); //Target fractions are xs_n-jet/xs_inclusive
      httStitching.SetDYInputYields(99536185, 66994632, 19918141, 5554767, 2812482);
    }
    if(era_type == era::data_2018UL) {
      httStitching.SetDYInputCrossSections(1.0, 0.1641, 0.0571, 0.0208, 0.0118); //Target fractions are xs_n-jet/xs_inclusive
      httStitching.SetDYInputYields(96055748+101373493, 60368985, 27494377, 20466041, 8885353);
    }
    if(era_type == era::data_2022_preEE ) {
      httStitching.SetDYInputCrossSections(5455.0,978.3,315.1,93.7,45.4);
      httStitching.SetDYInputYields(74978096,0,14777357,8699442,3267376);
    }
    if(era_type == era::data_2022_postEE) {
      httStitching.SetDYInputCrossSections(5455.0,978.3,315.1,93.7,45.4);
      httStitching.SetDYInputYields(244244756,50224518,50586120,28833220,9405785);
    } 
  }

  if ((output_name.find("DYJetsToLL-NLO") != output_name.npos || output_name.find("DYJetstoLL-NLO") != output_name.npos || output_name.find("DYJetsToLL_0J-NLO") != output_name.npos|| output_name.find("DYJetsToLL_1J-NLO") != output_name.npos || output_name.find("DYJetsToLL_2J-NLO") != output_name.npos )){
    httStitching.set_do_dy_soup_NLO(true);
    if(era_type == era::data_2016UL_preVFP){
      httStitching.SetDYInputCrossSections_NLO(6404.0,5129.0,951.5,361.4);
      httStitching.SetDYInputYields_NLO(61134271, 64114200, 41519367, 13244344);
    }
    if(era_type == era::data_2016UL_postVFP){
      httStitching.SetDYInputCrossSections_NLO(6404.0,5129.0,951.5,361.4);
      httStitching.SetDYInputYields_NLO(49693903, 61458133, 38644052, 12915771);
    }
    if(era_type == era::data_2017UL) {
      httStitching.SetDYInputCrossSections_NLO(6404.0,5129.0,951.5,361.4);
      httStitching.SetDYInputYields_NLO(132089877, 63067016, 39621895, 14506834);
    }
    if(era_type == era::data_2018UL) {
      httStitching.SetDYInputCrossSections_NLO(6404.0,5129.0,951.5,361.4);
      httStitching.SetDYInputYields_NLO(132321936, 69493624, 43156627, 13793541);
    }
    if(era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) {
      httStitching.SetDYInputCrossSections_NLO(6404.0,5129.0,951.5,361.4);
      httStitching.SetDYInputYields_NLO(132321936, 69493624, 43156627, 13793541);
    }    
  }
  BuildModule(httStitching);   
}

if(js["baseline"]["do_ff_weights"].asBool() && (addit_output_folder=="" || addit_output_folder.find("TSCALE")!=std::string::npos || addit_output_folder.find("ESCALE")!=std::string::npos)){
  BuildModule(HTTFakeFactorWeights("HTTFakeFactorWeights")
      .set_channel(channel)
      .set_era(era_type)
      .set_ditau_label("ditau")
      .set_met_label(met_label)
      .set_jets_label(jets_label)
      .set_strategy(strategy_type)
      .set_categories(js["baseline"]["ff_categories"].asString())
      .set_do_systematics(js["baseline"]["do_ff_systematics"].asBool() && addit_output_folder=="")
      .set_ff_file(js["baseline"]["ff_file"].asString())
      .set_fracs_file(js["baseline"]["ff_fracs_file"].asString())
      .set_is_embedded(is_embedded)
      );
}   

std::vector<std::string> jet_met_uncerts = {""};
if(js["do_jetmet_uncerts"].asBool() && addit_output_folder=="" && !is_embedded) { // only loop for nominal case
  jet_met_uncerts = {"","scale_met_lo","scale_met_hi","res_met_lo","res_met_hi","scale_j_relbal_lo","scale_j_relbal_hi","scale_j_abs_lo","scale_j_abs_hi","scale_j_abs_year_lo","scale_j_abs_year_hi","scale_j_flav_lo","scale_j_flav_hi","scale_j_bbec1_lo","scale_j_bbec1_hi","scale_j_bbec1_year_lo","scale_j_bbec1_year_hi","scale_j_ec2_lo","scale_j_ec2_hi","scale_j_ec2_year_lo","scale_j_ec2_year_hi","scale_j_hf_lo","scale_j_hf_hi","scale_j_hf_year_lo","scale_j_hf_year_hi","scale_j_relsamp_year_lo","scale_j_relsamp_year_hi","res_j_lo","res_j_hi","met_uncl_lo","met_uncl_hi"};
}
if(js["do_jetmet_uncerts"].asBool() && addit_output_folder=="" && is_embedded) { // only loop for nominal case
  jet_met_uncerts = {"","scale_embed_met_hi","scale_embed_met_lo"};
}

for (unsigned i=0; i<jet_met_uncerts.size(); ++i) {
  std::string uncert = jet_met_uncerts[i];
  // from here all jet and met recoil uncertainties should be run in a loop

  int jer_mode_ = jer_mode;
  int jes_mode_ = jes_mode;
  int metscale_mode_ = metscale_mode;
  int metres_mode_ = metres_mode;
  int metuncl_mode_ = 0;
  int embed_metscale_mode_ = 0;

  std::string svfit_folder_ = svfit_folder; 
  std::string add_dir = ""; 
  if (uncert != "") {
    // if we are doing a JET/MET shift then reset the modes
    if(uncert.find("scale_j")!=std::string::npos) jes_mode_ = unmod_js["sequences"][uncert]["baseline"]["jes_mode"].asUInt();
    if(uncert.find("scale_j")!=std::string::npos) alt_jes_input_set = unmod_js["sequences"][uncert]["baseline"]["jes_input_set"].asString();
    if(uncert.find("res_j")!=std::string::npos) jer_mode_ = unmod_js["sequences"][uncert]["baseline"]["jer_mode"].asUInt();
    if(uncert.find("res_met")!=std::string::npos) metres_mode_ = unmod_js["sequences"][uncert]["baseline"]["metres_mode"].asUInt();
    if(uncert.find("scale_met")!=std::string::npos) metscale_mode_ = unmod_js["sequences"][uncert]["baseline"]["metscale_mode"].asUInt();
    if(uncert.find("met_uncl")!=std::string::npos) metuncl_mode_ = unmod_js["sequences"][uncert]["baseline"]["metuncl_mode"].asUInt();
    if(uncert.find("scale_embed_met")!=std::string::npos) embed_metscale_mode_ = unmod_js["sequences"][uncert]["baseline"]["embed_metscale_mode"].asUInt();

    std::string add_dir = unmod_js["sequences"][uncert]["baseline"]["addit_output_folder"].asString(); 
 
    fs_vec.push_back(std::make_shared<fwlite::TFileService>(
       (output_folder+"/"+add_dir+"/"+output_name).c_str()));
 
    if(new_svfit_mode > 0){
      if(!((jes_mode_ && do_recoil) || (jer_mode_ > 0 && do_recoil))) svfit_folder_=svfit_folder+add_dir+"/";
    }

  } else {
      fs_vec.push_back(fs);
  }

  //if processing a sample with no recoil correction uncertainties then skip the recoil uncertainties
  if(!do_recoil && (metscale_mode_>0||metres_mode_>0)) continue;
  //if processing a sample with recoil correction uncertainties then skip the unclustered energy uncertainty
  if(do_recoil && metuncl_mode_>0) continue;


  // first make a deep copy of the jets and the met - we will need to delete these after we are finished with them  

  std::string shift_jets_label = jets_label+"_uncert_"+uncert; 
  std::string shift_met_label = met_label+"_uncert_"+uncert; 
  BuildModule(CopyCollection<PFJet>("CopyJetsForUncert_"+uncert,jets_label,shift_jets_label,true));

  std::string ext = "";
  if(is_embedded && embed_metscale_mode_==1) ext = "_down";
  if(is_embedded && embed_metscale_mode_==2) ext = "_up";

  BuildModule(GenericModule("CopyMetForUncert_"+uncert)
    .set_function([=](ic::TreeEvent *event){
       Met* met = event->GetPtr<Met>(met_label+ext);
       Met *met_copy = new Met(*met);
       event->Add(shift_met_label,met_copy);
       return 0;
    }));

  if(metuncl_mode_>0) {
    std::string shift_type = "NoShift";
    if(metuncl_mode_==1) shift_type = "UnclusteredEnDown";
    if(metuncl_mode_==2) shift_type = "UnclusteredEnUp";
    BuildModule(GenericModule("MetUnclusteredEnergyShift")
      .set_function([=](ic::TreeEvent *event){
         Met* met = event->GetPtr<Met>(shift_met_label);
         double pt = met->GetShiftedMet(shift_type).pt();
         double phi = met->GetShiftedMet(shift_type).phi();
         met->set_vector(ROOT::Math::PtEtaPhiEVector(pt,0,phi,pt));
         return 0;
      }));
  }

  if (jer_mode_>0 && !is_data && !is_embedded) {
     BuildModule(JetEnergyResolution<PFJet>("JetEnergyResolution")
       .set_input_label(shift_jets_label)
       .set_met_label(shift_met_label)
       .set_jer_shift_mode(jer_mode_)
       .set_EENoiseFix(era_type == era::data_2017 || era_type == era::data_2017UL)
       .set_shift_met(!do_recoil)
     );
  }
  
  if (jes_mode_ > 0 && !is_data ){
    std::string jes_input_file  = "";
    std::string jes_input_set  = "";
    if (era_type == era::data_2016) {
      jes_input_file = "input/jec/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt"; 
      jes_input_set  = "Total";
    }
    if (era_type == era::data_2017) {
      jes_input_file = "input/jec/RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt";
      jes_input_set  = "Total";
    }
    if (era_type == era::data_2018) {
      jes_input_file = "input/jec/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt";
      jes_input_set  = "Total";
    }
    if (era_type == era::data_2016UL_preVFP) {
      jes_input_file = "input/jec/RegroupedV2_Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt";
      jes_input_set  = "Total";
    }
    if (era_type == era::data_2016UL_postVFP) {
      jes_input_file = "input/jec/RegroupedV2_Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt";
      jes_input_set  = "Total";
    }
    if (era_type == era::data_2017UL) {
      jes_input_file = "input/jec/RegroupedV2_Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt";
      jes_input_set  = "Total";
    }
    if (era_type == era::data_2018UL) {
      jes_input_file = "input/jec/RegroupedV2_Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt";
      jes_input_set  = "Total";
    }
    if (era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) {
      jes_input_file = "input/jec/RegroupedV2_Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt";
      jes_input_set  = "Total";
    }    

    
    if(alt_jes_input_set!="") jes_input_set = alt_jes_input_set;
   
    std::string jes_input_set_ = jes_input_set;

    BuildModule(JetEnergyUncertainty<PFJet>("JetEnergyUncertainty")
      .set_input_label(shift_jets_label)
      .set_met_label(shift_met_label)
      .set_jes_shift_mode(jes_mode_)
      .set_uncert_file(jes_input_file)
      .set_uncert_set(jes_input_set_)
      .set_EENoiseFix(era_type == era::data_2017 || era_type == era::data_2017)
      .set_shift_met(!do_recoil)
     );
   
  }

  bool hem_fix=false; // disable this "fix" for HEM issue (this is really only used to check impact of HEM problem rather than actually fix it)
  if (hem_fix && !is_data && !is_embedded &&  (era_type == era::data_2018 || era_type == era::data_2018UL)) {
    BuildModule(JetEnergyUncertainty<PFJet>("HEMFix")
      .set_input_label(shift_jets_label)
      .set_met_label(shift_met_label)
      .set_jes_shift_mode(0)
      .set_uncert_file("")
      .set_uncert_set("")
      .set_EENoiseFix(false)
      .set_HEMFix(true)
      .set_shift_met(!do_recoil)
     );
  } 
  
  if(js["do_btag_eff"].asBool()){
     BuildModule(BTagCheck("BTagCheck")
      .set_fs(fs.get())
      .set_channel(channel)
      .set_era(era_type)
      .set_strategy(strategy_type)
      .set_do_legacy(false)
      .set_use_deep_csv(use_deep_csv)
      .set_use_deep_jet(use_deep_jet)
      .set_wp_to_check(wp_to_check)
      .set_jet_label(shift_jets_label));
  }
  
  if(!is_data && !is_embedded){
    TH2F bbtag_eff;
    TH2F cbtag_eff;
    TH2F othbtag_eff;
    TH2F bbtag_eff_alt;
    TH2F cbtag_eff_alt;
    TH2F othbtag_eff_alt;
 

    TH2D deepjet_bbtag_eff;
    TH2D deepjet_cbtag_eff;
    TH2D deepjet_othbtag_eff;
	
 
    if (era_type == era::data_2016){
      bbtag_eff = GetFromTFile<TH2F>("input/btag_sf/deepCSV_efficiencies_leg2016.root","/","btag_eff_b");
      cbtag_eff = GetFromTFile<TH2F>("input/btag_sf/deepCSV_efficiencies_leg2016.root","/","btag_eff_c");
      othbtag_eff = GetFromTFile<TH2F>("input/btag_sf/deepCSV_efficiencies_leg2016.root","/","btag_eff_oth");
  
      bbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/deepCSV_efficiencies_leg2016_loose.root","/","btag_eff_b");
      cbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/deepCSV_efficiencies_leg2016_loose.root","/","btag_eff_c");
      othbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/deepCSV_efficiencies_leg2016_loose.root","/","btag_eff_oth");

      deepjet_bbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2016_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_b");
      deepjet_cbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2016_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_c");
      deepjet_othbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2016_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_oth");
    } else if (era_type == era::data_2016UL_preVFP) {
      bbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpreVFP_medium.root","/","btag_eff_b");
      cbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpreVFP_medium.root","/","btag_eff_c");
      othbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpreVFP_medium.root","/","btag_eff_oth");
  
      bbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpreVFP_loose.root","/","btag_eff_b");
      cbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpreVFP_loose.root","/","btag_eff_c");
      othbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpreVFP_loose.root","/","btag_eff_oth");
	  
    } else if (era_type == era::data_2016UL_postVFP) {
      bbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpostVFP_medium.root","/","btag_eff_b");
      cbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpostVFP_medium.root","/","btag_eff_c");
      othbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpostVFP_medium.root","/","btag_eff_oth");
  
      bbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpostVFP_loose.root","/","btag_eff_b");
      cbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpostVFP_loose.root","/","btag_eff_c");
      othbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2016ULpostVFP_loose.root","/","btag_eff_oth");
	  
    } else if (era_type == era::data_2017) {
      bbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2017.root","/","btag_eff_b");
      cbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2017.root","/","btag_eff_c");
      othbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2017.root","/","btag_eff_oth");
  
      bbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2017_loose.root","/","btag_eff_b");
      cbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2017_loose.root","/","btag_eff_c");
      othbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2017_loose.root","/","btag_eff_oth");

      deepjet_bbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2017_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_b");
      deepjet_cbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2017_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_c");
      deepjet_othbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2017_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_oth");
    } else if (era_type == era::data_2017UL) {
      bbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2017UL_medium.root","/","btag_eff_b");
      cbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2017UL_medium.root","/","btag_eff_c");
      othbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2017UL_medium.root","/","btag_eff_oth");
  
      bbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2017UL_loose.root","/","btag_eff_b");
      cbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2017UL_loose.root","/","btag_eff_c");
      othbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2017UL_loose.root","/","btag_eff_oth");
	  
    } else if (era_type == era::data_2018) {
      bbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2018.root","/","btag_eff_b");
      cbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2018.root","/","btag_eff_c");
      othbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2018.root","/","btag_eff_oth");
  
      bbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2018_loose.root","/","btag_eff_b");
      cbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2018_loose.root","/","btag_eff_c");
      othbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2018_loose.root","/","btag_eff_oth");

      deepjet_bbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2018_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_b");
      deepjet_cbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2018_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_c");
      deepjet_othbtag_eff = GetFromTFile<TH2D>("input/btag_sf/Eff_DeepFlavour_2018_medium_all_proc_DeepFlavour_medium_inclusive_inclusive.root","/","btag_eff_oth");
    } else if (era_type == era::data_2018UL) {
      bbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_medium.root","/","btag_eff_b");
      cbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_medium.root","/","btag_eff_c");
      othbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_medium.root","/","btag_eff_oth");
  
      bbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_loose.root","/","btag_eff_b");
      cbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_loose.root","/","btag_eff_c");
      othbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_loose.root","/","btag_eff_oth");
	  
    } else if (era_type == era::data_2022_preEE || era_type == era::data_2022_postEE) {
      bbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_medium.root","/","btag_eff_b");
      cbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_medium.root","/","btag_eff_c");
      othbtag_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_medium.root","/","btag_eff_oth");
  
      bbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_loose.root","/","btag_eff_b");
      cbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_loose.root","/","btag_eff_c");
      othbtag_eff_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepJet_2018UL_loose.root","/","btag_eff_oth");
	  
    }

    BTagWeightLegacyRun2 BTagWeight = BTagWeightLegacyRun2("BTagWeightLegacyRun2")
     .set_channel(channel)
     .set_era(era_type)
     .set_strategy(strategy_type)
     .set_jet_label(shift_jets_label)
 	 .set_bbtag_eff(new TH2F(bbtag_eff))
     .set_cbtag_eff(new TH2F(cbtag_eff))
     .set_othbtag_eff(new TH2F(othbtag_eff))
     .set_bbtag_eff_alt(new TH2F(bbtag_eff_alt))
     .set_cbtag_eff_alt(new TH2F(cbtag_eff_alt))
     .set_othbtag_eff_alt(new TH2F(othbtag_eff_alt))   
     .set_use_deep_jet(use_deep_jet)
     .set_use_deep_csv(use_deep_csv);

   // for MSSM we use deepjet and use promote demote method
   // BuildModule(BTagWeightRun2("DeepJetBTagWeightRun2")
    // .set_channel(channel)
    // .set_era(era_type)
     //.set_strategy(strategy_type)
    // .set_jet_label(shift_jets_label)
     //.set_bbtag_eff((TH2F*)(new TH2D(deepjet_bbtag_eff)))
    // .set_cbtag_eff((TH2F*)(new TH2D(deepjet_cbtag_eff)))
     //.set_othbtag_eff((TH2F*)(new TH2D(deepjet_othbtag_eff)))
     //.set_do_reshape(do_reshape)
     //.set_use_deep_jet(true)
     //.set_use_deep_csv(false)
     //.set_btag_mode(btag_mode)
     //.set_bfake_mode(bfake_mode));

    if((era_type == era::data_2016 || era_type == era::data_2016UL_preVFP || era_type == era::data_2016UL_postVFP) && (output_name.find("TTTo2L2Nu") != output_name.npos || output_name.find("TTToHadronic") != output_name.npos || output_name.find("TTToSemiLeptonic") != output_name.npos)) {
      BTagWeight.set_do_cp5_2016(true);
      TH2F bbtag_cp5_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2016_CP5.root","/","btag_eff_b");
      TH2F cbtag_cp5_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2016_CP5.root","/","btag_eff_c");
      TH2F othbtag_cp5_eff = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2016_CP5.root","/","btag_eff_oth");
      TH2F bbtag_eff_cp5_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2016_CP5_loose.root","/","btag_eff_b");
      TH2F cbtag_eff_cp5_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2016_CP5_loose.root","/","btag_eff_c");
      TH2F othbtag_eff_cp5_alt = GetFromTFile<TH2F>("input/btag_sf/tagging_efficiencies_deepCSV_2016_CP5_loose.root","/","btag_eff_oth");
      BTagWeight.set_bbtag_eff_cp5(new TH2F(bbtag_cp5_eff));
      BTagWeight.set_cbtag_eff_cp5(new TH2F(cbtag_cp5_eff));
      BTagWeight.set_othbtag_eff_cp5(new TH2F(othbtag_cp5_eff));
      BTagWeight.set_bbtag_eff_cp5_alt(new TH2F(bbtag_eff_cp5_alt));
      BTagWeight.set_cbtag_eff_cp5_alt(new TH2F(cbtag_eff_cp5_alt));
      BTagWeight.set_othbtag_eff_cp5_alt(new TH2F(othbtag_eff_cp5_alt));
   }
 
   BuildModule(BTagWeight);
  }
  
  //if(usePFMET && do_recoil ){ 
  if(true ){ // for ML project we always call this so that we calculate U1 and U2, the corrections do not actually get applied!!
    BuildModule(HTTRun2RecoilCorrector("HTTRun2RecoilCorrector")
     .set_sample(output_name)
     .set_channel(channel)
     .set_mc(mc_type)
     .set_met_label(shift_met_label)
     .set_jets_label(shift_jets_label)
     .set_strategy(strategy_type)
     .set_use_quantile_map(false)
     .set_met_scale_mode(metscale_mode_)
     .set_met_res_mode(metres_mode_)
     .set_store_boson_pt(js["make_sync_ntuple"].asBool()));
  }
  
  if(!usePFMET){
    unsigned njets_mode = js["njets_mode"].asUInt();
  
    BuildModule(HTTLegacyRun2RecoilCorrector("HTTLegacyRun2RecoilCorrector")
     .set_sample(output_name)
     .set_channel(channel)
     .set_mc(mc_type)
     .set_met_label(shift_met_label)
     .set_jets_label(shift_jets_label)
     .set_strategy(strategy_type)
     .set_use_quantile_map(true) // use quantile mapping now
     .set_use_puppimet(true) 
     .set_met_scale_mode(metscale_mode_)
     .set_met_res_mode(metres_mode_)
     .set_store_boson_pt(js["make_sync_ntuple"].asBool())
     .set_njets_mode(njets_mode)
     .set_do_recoil(do_recoil)
     );

  }
  
  int mode = new_svfit_mode==1 && (jes_mode_ > 0 || jer_mode_ > 0) && do_recoil ? 0 : new_svfit_mode;
  
  SVFitTest svFitTest  = SVFitTest("SVFitTest")
    .set_channel(channel)
    .set_strategy(strategy_type)
    .set_mc(mc_type)
    .set_outname(svfit_override == "" ? output_name : svfit_override)
    .set_run_mode(mode)
    .set_fail_mode(0)
    .set_require_inputs_match(false)
    .set_split(40000)
    .set_dilepton_label("ditau")
    .set_met_label(shift_met_label)
    .set_fullpath(svfit_folder_)
    .set_legacy_svfit(true)
    .set_do_preselection(false)
    .set_MC(true)
    .set_do_light(true)
    .set_do_vloose_preselection(js["baseline"]["do_ff_weights"].asBool())
    .set_era(era_type);
   svFitTest.set_legacy_svfit(false);
   //svFitTest.set_do_preselection(!js["make_sync_ntuple"].asBool() && !js["baseline"]["do_faketaus"].asBool());
   svFitTest.set_do_preselection(false); // this is done elsewhere
   svFitTest.set_read_svfit_mt(true);
   svFitTest.set_tau_optimisation(js["tau_id_study"].asBool());
   svFitTest.set_read_all(js["read_all_svfit_files"].asBool());
   svFitTest.set_from_grid(js["svfit_from_grid"].asBool());
  
  BuildModule(svFitTest);
  
  if(channel != channel::tpzmm &&channel !=channel::tpzee && (is_data || js["trg_in_mc"].asBool())){

    bool tau_dataset = output_name.find("Tau") != output_name.npos && is_data; 
    BuildModule(HTTTriggerFilter("HTTTriggerFilter")
        .set_channel(channel)
        .set_mc(mc_type)
        .set_era(era_type)
        .set_strategy(strategy_type)
        .set_is_data(is_data)
        .set_tau_dataset(tau_dataset)
        .set_is_embedded(is_embedded)
        .set_do_singletau(js["do_singletau"].asBool())
        .set_do_filter(false)
        .set_pair_label("ditau"));
  }
  
  
  bool do_mssm_higgspt = output_name.find("SUSYGluGluToHToTauTau_M") != output_name.npos;
  BuildModule(HTTCategories("HTTCategories")
      .set_fs(fs_vec[i].get())
      .set_channel(channel)
      .set_era(era_type)
      .set_strategy(strategy_type)
      .set_ditau_label("ditau")
      .set_met_label(shift_met_label)
      .set_jets_label(shift_jets_label)
      .set_make_sync_ntuple(js["make_sync_ntuple"].asBool())
      .set_make_mva_ntuple(js["make_mva_ntuple"].asBool())
      .set_sync_output_name(js["output_folder"].asString()+"/SYNCFILE_"+output_name)
      .set_mva_output_name(js["output_folder"].asString()+"/MVAFILE_"+output_name)
      .set_is_embedded(is_embedded)
      .set_is_data(is_data)
      .set_systematic_shift(addit_output_folder!="" && i==0)
      //Good to avoid accidentally overwriting existing output files when syncing
      .set_write_tree(!js["make_sync_ntuple"].asBool())
      .set_do_ff_weights(js["baseline"]["do_ff_weights"].asBool())
      .set_do_ff_systematics(i==0 && js["baseline"]["do_ff_systematics"].asBool()&& (addit_output_folder=="" || addit_output_folder.find("TSCALE")!=std::string::npos || addit_output_folder.find("ESCALE")!=std::string::npos || addit_output_folder.find("MUSCALE")!=std::string::npos))
      .set_do_qcd_scale_wts(do_qcd_scale_wts_)
      .set_do_mssm_higgspt(do_mssm_higgspt)
      .set_do_faketaus(js["baseline"]["do_faketaus"].asBool())
      .set_trg_applied_in_mc(js["trg_in_mc"].asBool()));

  // now we need to delete the deep copies of the jet and met collections
  BuildModule(GenericModule("DeleteCopiedJetMET")
    .set_function([=](ic::TreeEvent *event){
       Met* met = event->GetPtr<Met>(shift_met_label);
       delete met;
       event->Remove(shift_met_label);
       std::vector<PFJet*> & jets = event->GetPtrVec<PFJet>(shift_jets_label);
       for (auto j : jets) delete j;
       jets.clear();
       event->Remove(shift_jets_label);
       return 0;
    }));

  
} // end of loop over jet met uncertainties

  if(channel == channel::tpzmm){   
    std::function<bool(Muon const*)> muon_probe_id;
    
    if(era_type == era::data_2017 || era_type == era::data_2017UL || era_type == era::data_2018 || era_type == era::data_2018UL){

        muon_probe_id = [](Muon const* m) {return MuonMedium(m); };
        std::function<bool(Muon const*)> MuonLooseID = [](Muon const* m) { return MuonLoose(m) && m->is_global(); }; 
        // low pT mu leg of e+mu trigger
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_EMLow")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu27")
          .set_tag_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id) 
          .set_probe_trg_objects("triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)
          .set_extra_l1_probe_pt(5.)
          .set_add_name("_em_low")
        );
        // high pT mu leg of e+mu trigger
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_EMHigh")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu27")
          .set_tag_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(23.)
          .set_extra_l1_probe_pt(20.)
          .set_add_name("_em_high")
        );
        // mu leg of mu+tau cross trigger for runs C-F
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_MT_CtoF")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu27")
          .set_tag_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsDoubleMu20")
          .set_probe_trg_filters("hltL3crIsoL1sDoubleMu18erL1f0L2f10QL3f20QL3trkIsoFiltered0p07")
          .set_add_name("_mt_c_to_f")
        );
        // mu leg of mu+tau cross trigger for run B
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_MT_B")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu27")
          .set_tag_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsIsoMu20Tau27")
          .set_probe_trg_filters("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07")
          .set_do_extra(true)
          .set_add_name("_mt_b")
        );
        // single muons trigger measurment
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_single")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu27")
          .set_tag_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsIsoMu27,triggerObjectsIsoMu24")
          .set_probe_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07,hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07")  
        );
        // low pT leg of dimuon trigger used for embedded selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuLow")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu27")
          .set_tag_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")
          .set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsIsoMu27,triggerObjectsIsoMu24")
          .set_probe_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07,hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07")
          .set_add_name("_dimu_low")
          .set_probe_trg_objects("triggerObjectsMu17Mu8DZ,triggerObjectsMu17Mu8DZmass8,triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)

        );
        // high pT leg of dimuon trigger used for embedded selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuhigh")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu27")
          .set_tag_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")
          .set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsIsoMu27,triggerObjectsIsoMu24")
          .set_probe_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07,hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07")
          .set_add_name("_dimu_high")
          .set_probe_trg_objects("triggerObjectsMu17Mu8DZ,triggerObjectsMu17Mu8DZmass8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(17.)
        );
        // use to measure DZ filter efficiency for embedded selection selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuDZ")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu27")
          .set_tag_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")
          .set_probe_id(MuonLooseID)
          .set_tag_id(MuonLooseID)
          .set_probe_trg_objects("triggerObjectsIsoMu27,triggerObjectsIsoMu24")
          .set_probe_trg_filters("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07,hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07")
          .set_add_name("_dimu_dz")
          .set_probe_trg_objects("triggerObjectsMu17Mu8DZ,triggerObjectsMu17Mu8DZmass8,triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)
          .set_do_dzmass(true)
        );
         
    }else if(era_type == era::data_2022_postEE || era_type == era::data_2022_preEE){
        std::string isoMu24FilterName = (era_type == era::data_2022_postEE) ?
         "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered" : "hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p08";
        muon_probe_id = [](Muon const* m) {return MuonMedium(m); };
        std::function<bool(Muon const*)> MuonLooseID = [](Muon const* m) { return MuonLoose(m) && m->is_global(); }; 
        // low pT mu leg of e+mu trigger
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_EMLow")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters(isoMu24FilterName)
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id) 
          .set_probe_trg_objects("triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)
          .set_extra_l1_probe_pt(5.)
          .set_add_name("_em_low")
        );
        // high pT mu leg of e+mu trigger
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_EMHigh")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters(isoMu24FilterName)
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(23.)
          .set_extra_l1_probe_pt(20.)
          .set_add_name("_em_high")
        );
        // mu leg of mu+tau cross trigger for runs C-F
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_MT_CtoF")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters(isoMu24FilterName)
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsDoubleMu20")
          .set_probe_trg_filters("hltL3crIsoL1sDoubleMu18erL1f0L2f10QL3f20QL3trkIsoFiltered0p07")
          .set_add_name("_mt_c_to_f")
        );
        // mu leg of mu+tau cross trigger for run B
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_MT_B")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters(isoMu24FilterName)
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsIsoMu20Tau27")
          .set_probe_trg_filters("hltL3crIsoL1sMu18erTau24erIorMu20erTau24erL1f0L2f10QL3f20QL3trkIsoFiltered0p07")
          .set_do_extra(true)
          .set_add_name("_mt_b")
        );
        // single muons trigger measurment
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_single")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")                
          .set_tag_trg_filters(isoMu24FilterName)
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsIsoMu24")
          .set_probe_trg_filters(isoMu24FilterName)  
        );
        // low pT leg of dimuon trigger used for embedded selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuLow")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters(isoMu24FilterName)
          .set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsIsoMu24")
          .set_probe_trg_filters(isoMu24FilterName)
          .set_add_name("_dimu_low")
          .set_probe_trg_objects("triggerObjectsMu17Mu8DZ,triggerObjectsMu17Mu8DZmass8,triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)

        );
        // high pT leg of dimuon trigger used for embedded selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuhigh")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters(isoMu24FilterName)
          .set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_probe_trg_objects("triggerObjectsIsoMu24")
          .set_probe_trg_filters(isoMu24FilterName)
          .set_add_name("_dimu_high")
          .set_probe_trg_objects("triggerObjectsMu17Mu8DZ,triggerObjectsMu17Mu8DZmass8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(17.)
        );
        // use to measure DZ filter efficiency for embedded selection selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuDZ")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters(isoMu24FilterName)
          .set_probe_id(MuonLooseID)
          .set_tag_id(MuonLooseID)
          .set_probe_trg_objects("triggerObjectsIsoMu24")
          .set_probe_trg_filters(isoMu24FilterName)
          .set_add_name("_dimu_dz")
          .set_probe_trg_objects("triggerObjectsMu17Mu8DZ,triggerObjectsMu17Mu8DZmass8,triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4,hltDiMuon178RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)
          .set_do_dzmass(true)
        );
         
    } else {
      muon_probe_id = [](Muon const* m) {return MuonMedium(m); };
      std::function<bool(Muon const*)> MuonLooseID = [](Muon const* m) { return MuonLoose(m) && m->is_global(); };
      // for mu8 leg of EMu cross-trigger
      BuildModule(TagAndProbe<Muon const*>("TagAndProbe_emLow")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")
          .set_probe_id(muon_probe_id)
          //.set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_add_name("_emLow")          
          .set_probe_trg_objects("triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)
          .set_extra_l1_probe_pt(5.)
      );          
          
          // for mu23 leg of EMu cross-trigger
      BuildModule(TagAndProbe<Muon const*>("TagAndProbe_emHigh")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")
          .set_probe_id(muon_probe_id)
          //.set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_add_name("_emHigh")          
          .set_probe_trg_objects("triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(23.)
          .set_extra_l1_probe_pt(20.)
      );
          
      // for single muon trigger:
      BuildModule(TagAndProbe<Muon const*>("TagAndProbe_Single")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")
          .set_probe_id(muon_probe_id)
          //.set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_add_name("_Single")                    
          .set_probe_trg_objects("triggerObjectsIsoMu22,triggerObjectsIsoTkMu22,triggerObjectsIsoMu22Eta2p1,triggerObjectsIsoTkMu22Eta2p1")
          .set_probe_trg_filters("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09,hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09,hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09,hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09")
      );

          // for muon leg of mu+tau cross trigger:
      BuildModule(TagAndProbe<Muon const*>("TagAndProbe_mt")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")
          .set_probe_id(muon_probe_id)
          .set_tag_id(muon_probe_id)
          .set_add_name("_mt")          
          .set_probe_trg_objects("triggerObjectsIsoMu19LooseTau20SingleL1")
          .set_probe_trg_filters("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09")
      
      );

        // low pT leg of dimuon trigger used for embedded selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuLow")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")
          .set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_add_name("_dimu_low")
          .set_probe_trg_objects("triggerObjectsMu17TkMu8,triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4,hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)
        );

        // high pT leg of dimuon trigger used for embedded selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuhigh")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")
          .set_probe_id(MuonLooseID)
          .set_tag_id(muon_probe_id)
          .set_add_name("_dimu_high")
          .set_probe_trg_objects("triggerObjectsMu17TkMu8,triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4,hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(17.)
        );

        // use to measure DZ filter efficiency for embedded selection selection
        BuildModule(TagAndProbe<Muon const*>("TagAndProbe_DiMuDZ")
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsIsoMu24")
          .set_tag_trg_filters("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")
          .set_probe_id(MuonLooseID)
          .set_tag_id(MuonLooseID)
          .set_add_name("_dimu_dz")
          .set_probe_trg_objects("triggerObjectsMu17TkMu8,triggerObjectsMu17Mu8")
          .set_probe_trg_filters("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4,hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4")
          .set_extra_hlt_probe_pt(8.)
          .set_do_dzmass(true)
        );

    }
  } 

  if(channel == channel::tpzee){
    if(era_type == era::data_2022_postEE || era_type == era::data_2022_preEE){
      std::function<bool(Electron const*)> elec_probe_id = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); };
      std::function<bool(Electron const*)> elec_tag_id = [](Electron const* e) { return ElectronHTTIdFall17V2(e, false); };

      BuildModule(TagAndProbe<Electron const*>("TagAndProbe_single")
          .set_is_embedded(is_embedded)
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsEle30")
          .set_tag_trg_filters("hltEle30WPTightGsfTrackIsoFilter")
          .set_probe_id(elec_probe_id)
          .set_tag_id(elec_tag_id)
          .set_probe_trg_objects("triggerObjectsEle30")
          .set_probe_trg_filters("hltEle30WPTightGsfTrackIsoFilter")
      );


    }
    if(era_type == era::data_2017 || era_type == era::data_2017UL || era_type == era::data_2018 || era_type == era::data_2018UL
    || era_type == era::data_2022_preEE || era_type == era::data_2022_postEE){
      std::function<bool(Electron const*)> elec_probe_id = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); };
      std::function<bool(Electron const*)> elec_tag_id = [](Electron const* e) { return ElectronHTTIdFall17V2(e, false); };

      // for el12 leg of EMu cross-trigger
      BuildModule(TagAndProbe<Electron const*>("TagAndProbe_emLow")
          .set_is_embedded(is_embedded)
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsEle35")
          .set_tag_trg_filters("hltEle35noerWPTightGsfTrackIsoFilter")
          .set_extra_l1_tag_pt(32.) // ensure L1 was not prescaled during data-taking
          .set_probe_id(elec_probe_id)
          .set_tag_id(elec_tag_id)
          .set_probe_trg_objects("triggerObjectsEle24Ele12") 
          .set_probe_trg_filters("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter")
          .set_extra_l1_probe_pt(10.)
          .set_add_name("_emLow")
      );

      // for el23 leg of EMu cross-trigger
      BuildModule(TagAndProbe<Electron const*>("TagAndProbe_emHigh")
          .set_is_embedded(is_embedded)
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsEle35")
          .set_tag_trg_filters("hltEle35noerWPTightGsfTrackIsoFilter")
          .set_extra_l1_tag_pt(32.) // ensure L1 was not prescaled during data-taking
          .set_probe_id(elec_probe_id)
          .set_tag_id(elec_tag_id)
          .set_probe_trg_objects("triggerObjectsEle24Ele12") //Ele23 actually-> 
          .set_probe_trg_filters("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter")
          .set_extra_l1_probe_pt(23.)
          .set_extra_l1_iso_probe_pt(20.)
          .set_add_name("_emHigh")
      );

      if(strategy_type == strategy::cpdecays18) {
        BuildModule(TagAndProbe<Electron const*>("TagAndProbe_ET")
            .set_is_embedded(is_embedded)
            .set_fs(fs.get())
            .set_channel(channel)
            .set_era(era_type)
            .set_strategy(strategy_type)
            .set_ditau_label("ditau")
            .set_tag_trg_objects("triggerObjectsEle35")
            .set_tag_trg_filters("hltEle35noerWPTightGsfTrackIsoFilter")
            //.set_extra_l1_tag_pt(32.) // ensure L1 was not prescaled during data-taking
            .set_probe_id(elec_probe_id)
            .set_tag_id(elec_tag_id)
            .set_probe_trg_objects("triggerObjectsEle24Tau30,triggerObjectsEle24TauHPS30") // for 2018
            .set_probe_trg_filters("hltEle24erWPTightClusterShapeFilterForTau,hltEle24erWPTightClusterShapeFilterForTau") // for 2018
            .set_do_extra(true)
            .set_add_name("_ET")
        );

        BuildModule(TagAndProbe<Electron const*>("TagAndProbe_ET_v2")
            .set_is_embedded(is_embedded)
            .set_fs(fs.get())
            .set_channel(channel)
            .set_era(era_type)
            .set_strategy(strategy_type)
            .set_ditau_label("ditau")
            .set_tag_trg_objects("triggerObjectsEle35")
            .set_tag_trg_filters("hltEle35noerWPTightGsfTrackIsoFilter")
            .set_probe_id(elec_probe_id)
            .set_tag_id(elec_tag_id)
            .set_probe_trg_objects("triggerObjectsDoubleEl24") // use double electron monitoring trigger for 2018
            .set_probe_trg_filters("hltDoubleEle24erWPTightGsfTrackIsoFilterForTau") // for 2018
            .set_add_name("_ET_v2")
        );

        BuildModule(TagAndProbe<Electron const*>("TagAndProbe_single")
            .set_is_embedded(is_embedded)
            .set_fs(fs.get())
            .set_channel(channel)
            .set_era(era_type)
            .set_strategy(strategy_type)
            .set_ditau_label("ditau")
            .set_tag_trg_objects("triggerObjectsEle35")
            .set_tag_trg_filters("hltEle35noerWPTightGsfTrackIsoFilter")
            //.set_extra_l1_tag_pt(32.) // ensure L1 was not prescaled during data-taking
            .set_probe_id(elec_probe_id)
            .set_tag_id(elec_tag_id)
            .set_probe_trg_objects("triggerObjectsEle32,triggerObjectsEle35")
            .set_probe_trg_filters("hltEle32WPTightGsfTrackIsoFilter,hltEle35noerWPTightGsfTrackIsoFilter")
        );

        BuildModule(TagAndProbe<Electron const*>("TagAndProbe_single32")
            .set_is_embedded(is_embedded)
            .set_fs(fs.get())
            .set_channel(channel)
            .set_era(era_type)
            .set_strategy(strategy_type)
            .set_ditau_label("ditau")
            .set_tag_trg_objects("triggerObjectsEle35")
            .set_tag_trg_filters("hltEle35noerWPTightGsfTrackIsoFilter")
            //.set_extra_l1_tag_pt(32.) // ensure L1 was not prescaled during data-taking
            .set_probe_id(elec_probe_id)
            .set_tag_id(elec_tag_id)
            .set_probe_trg_objects("triggerObjectsEle32")
            .set_probe_trg_filters("hltEle32WPTightGsfTrackIsoFilter")
            .set_add_name("_single32")
        );

      } else {

        // for el23 leg of EMu cross-trigger
        BuildModule(TagAndProbe<Electron const*>("TagAndProbe_ET")
            .set_is_embedded(is_embedded)
            .set_fs(fs.get())
            .set_channel(channel)
            .set_era(era_type)
            .set_strategy(strategy_type)
            .set_ditau_label("ditau")
            .set_tag_trg_objects("triggerObjectsEle35")
            .set_tag_trg_filters("hltEle35noerWPTightGsfTrackIsoFilter")
            .set_extra_l1_tag_pt(32.) // ensure L1 was not prescaled during data-taking
            .set_probe_id(elec_probe_id)
            .set_tag_id(elec_tag_id)
            .set_probe_trg_objects("triggerObjectsEle24Tau30") // for 2017
            .set_probe_trg_filters("hltEle24erWPTightGsfTrackIsoFilterForTau") // for 2017
            .set_do_extra(true)
            .set_add_name("_ET")
        );

        BuildModule(TagAndProbe<Electron const*>("TagAndProbe_single")
            .set_is_embedded(is_embedded)
            .set_fs(fs.get())
            .set_channel(channel)
            .set_era(era_type)
            .set_strategy(strategy_type)
            .set_ditau_label("ditau")
            .set_tag_trg_objects("triggerObjectsEle35")
            .set_tag_trg_filters("hltEle35noerWPTightGsfTrackIsoFilter")
            .set_extra_l1_tag_pt(32.) // ensure L1 was not prescaled during data-taking
            .set_probe_id(elec_probe_id)
            .set_tag_id(elec_tag_id)
            .set_probe_trg_objects("triggerObjectsEle27,triggerObjectsEle32L1DoubleEG")
            .set_probe_trg_filters("hltEle27WPTightGsfTrackIsoFilter,hltEle32L1DoubleEGWPTightGsfTrackIsoFilter")
        );
      }

    } else {
        
      std::function<bool(Electron const*)> elec_probe_id;  
      std::function<bool(Electron const*)> elec_tag_id;  
      elec_probe_id = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); };
      elec_tag_id = [](Electron const* e) { return ElectronHTTIdFall17V2(e, false); };
      
          // for Ele12 leg of EMu cross-trigger
      BuildModule(TagAndProbe<Electron const*>("TagAndProbe_emLow")
          .set_is_embedded(is_embedded)
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsEle25GsfTightEta2p1")
          .set_tag_trg_filters("hltEle25erWPTightGsfTrackIsoFilter")
          .set_add_name("_emLow")          
          .set_probe_trg_objects("triggerObjectsEle23Ele12")
          .set_probe_trg_filters("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter")
          .set_extra_l1_probe_pt(10.)
          .set_probe_id(elec_probe_id)
          .set_tag_id(elec_tag_id)
          // em filters hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter -> electron 
      );   
          
          //// for Ele23 leg of EMu cross-trigger
      BuildModule(TagAndProbe<Electron const*>("TagAndProbe_emHigh")
          .set_is_embedded(is_embedded)
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_tag_trg_objects("triggerObjectsEle25GsfTightEta2p1")
          .set_tag_trg_filters("hltEle25erWPTightGsfTrackIsoFilter")
          .set_add_name("_emHigh")          
          .set_probe_trg_objects("triggerObjectsEle23Ele12")
          .set_probe_trg_filters("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter")
          .set_extra_l1_probe_pt(20.)
          .set_extra_l1_iso_probe_pt(18.)
          .set_probe_id(elec_probe_id)
          .set_tag_id(elec_tag_id)
      );
          // for single electron trigger
      BuildModule(TagAndProbe<Electron const*>("TagAndProbe_Single")
          .set_is_embedded(is_embedded)
          .set_fs(fs.get())
          .set_channel(channel)
          .set_era(era_type)
          .set_strategy(strategy_type)
          .set_ditau_label("ditau")
          .set_probe_id(elec_probe_id)
          .set_tag_id(elec_tag_id)
          .set_tag_trg_objects("triggerObjectsEle25GsfTightEta2p1")
          .set_tag_trg_filters("hltEle25erWPTightGsfTrackIsoFilter")          
          .set_add_name("_Single")          
          .set_probe_trg_objects("triggerObjectsEle25GsfTightEta2p1")
          .set_probe_trg_filters("hltEle25erWPTightGsfTrackIsoFilter")
          
      );  
    }
  } 

}
}

// --------------------------------------------------------------------------
// TT Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildTTPairs(){
 
 BuildTauSelection();

 BuildModule(CompositeProducer<Tau, Tau>("TTPairProducer")
      .set_input_label_first(js["taus"].asString())
      .set_input_label_second(js["taus"].asString())
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));
}

  
// --------------------------------------------------------------------------
// ET Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildETPairs() {
  ic::strategy strategy_type  = String2Strategy(strategy_str);

  BuildModule(CopyCollection<Electron>("CopyToSelectedElectrons",
      js["electrons"].asString(), "sel_electrons"));

  std::function<bool(Electron const*)> ElecID = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); };

  BuildModule(SimpleFilter<Electron>("ElectronFilter")
      .set_input_label("sel_electrons").set_min(1)
      .set_predicate([=](Electron const* e) {
        return ElecID(e);
      }));

  if (!is_embedded) {
    BuildModule(HTTSmearScale("ElectronSmearScaleCorrection")
        .set_input_label(js["electrons"].asString())
        .set_e_unc_mode(e_unc_mode)
    );
  }
  if(!is_data && is_embedded){
    BuildModule(HTTEnergyScale("ElectronEnergyScaleCorrection")
        .set_input_label(js["electrons"].asString())
        .set_shift(elec_shift_barrel)
        .set_shift_endcap(elec_shift_endcap)
        .set_strategy(strategy_type)
        .set_channel(channel::em)
        .set_moriond_corrections(moriond_tau_scale));
  }

  if (!is_data){
    if(elec_res_shift!=0.) {

     BuildModule(HTTGenMatchSelector<Electron>("HTTElectronGenMatchSelector")
      .set_input_vec_label(js["electrons"].asString())
      .set_output_vec_label("genmatched_elecs")
      .set_add_gen(true)
      .set_gen_match(mcorigin::tauE));

      BuildModule(EnergyResShifter<Electron>("ElecEnergyResShifter")
      .set_input_label("genmatched_elecs")
      .set_input_label_gen("genmatched_elecs_gen")
      .set_save_shifts(true)
      .set_shift_label("res_elecs")
      .set_shift(elec_res_shift));
    }
  }

  BuildModule(SimpleFilter<Electron>("ElectronFilter")
      .set_input_label("sel_electrons").set_min(1)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > elec_pt    &&
                fabs(e->eta())          < elec_eta   &&
                fabs(e->dxy_vertex())   < elec_dxy   &&
                fabs(e->dz_vertex())    < elec_dz ;

      }));

  BuildTauSelection();

  BuildModule(CompositeProducer<Electron, Tau>("ETPairProducer")
      .set_input_label_first("sel_electrons")
      .set_input_label_second(js["taus"].asString())
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));
}

// --------------------------------------------------------------------------
// MT Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildMTPairs() {

  if (mu_scale_mode > 0){
    BuildModule(HTTMuonEnergyScale("MuonEnergyScaleCorrection")
       .set_input_label("muons")
       .set_far_endcap(muon_shift_farendcap)
       .set_near_endcap(muon_shift_nearendcap)
       .set_barrel(muon_shift_barrel)
       );
  }

  BuildModule(CopyCollection<Muon>("CopyToSelectedMuons",
       js["muons"].asString(), "sel_muons"));

  if (!is_data){
    if(muon_res_shift!=0.) {

     BuildModule(HTTGenMatchSelector<Muon>("HTTMuonGenMatchSelector")
      .set_input_vec_label("sel_muons")
      .set_output_vec_label("genmatched_muons")
      .set_add_gen(true)
      .set_gen_match(mcorigin::tauMu));

      BuildModule(EnergyResShifter<Muon>("MuonEnergyResShifter")
      .set_input_label("genmatched_muons")
      .set_input_label_gen("genmatched_muons_gen")
      .set_save_shifts(true)
      .set_shift_label("res_muons")
      .set_shift(muon_res_shift));
    }
  }

  std::function<bool(Muon const*)> MuonID = [](Muon const* m) { return MuonMedium(m); };

  BuildModule(SimpleFilter<Muon>("MuonFilter")
    .set_input_label("sel_muons").set_min(1)
    .set_predicate([=](Muon const* m) {
      return  m->pt()                 > muon_pt    &&
              fabs(m->eta())          < muon_eta   &&
              fabs(m->dxy_vertex())   < muon_dxy   &&
              fabs(m->dz_vertex())    < muon_dz   &&
              MuonID(m);
    }));

  BuildTauSelection();

  if(js["do_mt_tagandprobe"].asBool()){
    BuildModule(OverlapFilter<Tau,Muon>("TauMuonOverlapFilter")
      .set_input_label(js["taus"].asString())
      .set_reference_label("loose_muons")
      .set_min_dr(0.5));
  }

  BuildModule(CompositeProducer<Muon, Tau>("MTPairProducer")
      .set_input_label_first("sel_muons")
      .set_input_label_second(js["taus"].asString())
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));
}


// --------------------------------------------------------------------------
// EM Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildEMPairs() {
  ic::strategy strategy_type  = String2Strategy(strategy_str);

  if (mu_scale_mode > 0){
    BuildModule(HTTMuonEnergyScale("MuonEnergyScaleCorrection")
       .set_input_label("muons")
       .set_far_endcap(muon_shift_farendcap)
       .set_near_endcap(muon_shift_nearendcap)
       .set_barrel(muon_shift_barrel)
       );
  }

  BuildModule(CopyCollection<Muon>("CopyToSelectedMuons",
      js["muons"].asString(), "sel_muons"));

  if (!is_data){
    if(muon_res_shift!=0.) {

     BuildModule(HTTGenMatchSelector<Muon>("HTTMuonGenMatchSelector")
      .set_input_vec_label("sel_muons")
      .set_output_vec_label("genmatched_muons")
      .set_add_gen(true)
      .set_gen_match(mcorigin::tauMu));

      BuildModule(EnergyResShifter<Muon>("MuonEnergyResShifter")
      .set_input_label("genmatched_muons")
      .set_input_label_gen("genmatched_muons_gen")
      .set_save_shifts(true)
      .set_shift_label("res_muons")
      .set_shift(muon_res_shift));
    }
  }

  std::function<bool(Muon const*)> MuonID = [](Muon const* m) { return MuonMedium(m); };

  BuildModule(SimpleFilter<Muon>("MuonFilter")
      .set_input_label("sel_muons").set_min(1)
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > muon_pt    &&
                fabs(m->eta())          < muon_eta   &&
                fabs(m->dxy_vertex())   < muon_dxy   &&
                fabs(m->dz_vertex())    < muon_dz   &&
                PF04IsolationVal(m, 0.5,0) < 0.5  &&
                MuonID(m);

      }));

  BuildModule(CopyCollection<Electron>("CopyToSelectedElectrons",
      js["electrons"].asString(), "sel_electrons"));


  std::function<bool(Electron const*)> ElecID = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); }; 

  BuildModule(SimpleFilter<Electron>("ElectronFilter")
      .set_input_label("sel_electrons").set_min(1) 
      .set_predicate([=](Electron const* e) {
        return ElecID(e);
      }));

  if (!is_embedded) {
    BuildModule(HTTSmearScale("ElectronSmearScaleCorrection")
        .set_input_label("sel_electrons")
        .set_e_unc_mode(e_unc_mode)
    );
  }
  if(!is_data && is_embedded){
    BuildModule(HTTEnergyScale("ElectronEnergyScaleCorrection")
        .set_input_label("sel_electrons")
        .set_shift(elec_shift_barrel)
        .set_shift_endcap(elec_shift_endcap)
        .set_strategy(strategy_type)
        .set_channel(channel::em)
        .set_moriond_corrections(moriond_tau_scale));
  }

  if (!is_data){
    if(elec_res_shift!=0.) {

     BuildModule(HTTGenMatchSelector<Electron>("HTTElectronGenMatchSelector")
      .set_input_vec_label("sel_electrons")
      .set_output_vec_label("genmatched_elecs")
      .set_add_gen(true)
      .set_gen_match(mcorigin::tauE));

      BuildModule(EnergyResShifter<Electron>("ElecEnergyResShifter")
      .set_input_label("genmatched_elecs")
      .set_input_label_gen("genmatched_elecs_gen")
      .set_save_shifts(true)
      .set_shift_label("res_elecs")
      .set_shift(elec_res_shift));
    }
  }

  BuildModule(SimpleFilter<Electron>("ElectronFilter")
      .set_input_label("sel_electrons").set_min(1)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > elec_pt    &&
                fabs(e->eta())          < elec_eta   &&
                fabs(e->dxy_vertex())   < elec_dxy   &&
                fabs(e->dz_vertex())    < elec_dz ;
      }));
  
  BuildModule(CompositeProducer<Electron, Muon>("EMPairProducer")
      .set_input_label_first("sel_electrons")
      .set_input_label_second("sel_muons")
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));

  BuildModule(SimpleFilter<CompositeCandidate>("EMPairFilter")
      .set_input_label("ditau").set_min(1)
      .set_predicate([=](CompositeCandidate const* c) {
        return PairOneWithPt(c, 24.0);
      }));
}

// --------------------------------------------------------------------------
// ZEE Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildZEEPairs() {
  ic::strategy strategy_type  = String2Strategy(strategy_str);


  BuildModule(CopyCollection<Electron>("CopyToSelectedElectrons",
      js["electrons"].asString(), "sel_electrons"));
  BuildModule(CopyCollection<Electron>("CopyToSelectedElectrons",
      js["electrons"].asString(), "sel_lead_electrons"));

  std::function<bool(Electron const*)> ElecID = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); }; 

  BuildModule(SimpleFilter<Electron>("ElectronIDFilter")
      .set_input_label("sel_electrons").set_min(2)
      .set_predicate([=](Electron const* e) {
        return ElecID(e);
      }));
  BuildModule(SimpleFilter<Electron>("LeadElectronIDFilter")
      .set_input_label("sel_lead_electrons").set_min(1)
      .set_predicate([=](Electron const* e) {
        return ElecID(e) ;
      }));

  if (!is_embedded) {
    BuildModule(HTTSmearScale("ElectronSmearScaleCorrection")
        .set_input_label(js["electrons"].asString())
        .set_e_unc_mode(e_unc_mode)
    );
  }
  if(!is_data && is_embedded){
    BuildModule(HTTEnergyScale("ElectronEnergyScaleCorrection")
        .set_input_label(js["electrons"].asString())
        .set_shift(elec_shift_barrel)
        .set_shift_endcap(elec_shift_endcap)
        .set_strategy(strategy_type)
        .set_channel(channel::em)
        .set_moriond_corrections(moriond_tau_scale));
  }

  BuildModule(SimpleFilter<Electron>("ElectronFilter")
      .set_input_label("sel_electrons").set_min(2)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > 20.        &&
                fabs(e->eta())          < 2.4        &&
                fabs(e->dxy_vertex())   < elec_dxy   &&
                fabs(e->dz_vertex())    < elec_dz ;

      }));
  
  BuildModule(SimpleFilter<Electron>("LeadElectronFilter")
      .set_input_label("sel_lead_electrons").set_min(1)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > elec_pt    &&
                fabs(e->eta())          < elec_eta   &&
                fabs(e->dxy_vertex())   < elec_dxy   &&
                fabs(e->dz_vertex())    < elec_dz ;

      }));

  BuildModule(CopyCollection<Electron>("CopyToFSRElectrons",
      js["electrons"].asString(), "fsr_electrons"));

  BuildModule(SimpleFilter<Electron>("FSRElectronIDFilter")
      .set_input_label("fsr_electrons").set_min(0)
      .set_predicate([=](Electron const* e) {
        return ElecID(e);
      }));

  BuildModule(SimpleFilter<Electron>("FSRElectronFilter")
      .set_input_label("fsr_electrons").set_min(0)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > 10.        &&
                fabs(e->eta())          < 2.4        &&
                fabs(e->dxy_vertex())   < elec_dxy   &&
                fabs(e->dz_vertex())    < elec_dz ;

      }));

  BuildModule(CopyCollection<Muon>("CopyToFSRMuons",
      js["muons"].asString(), "fsr_muons"));

  std::function<bool(Muon const*)> MuonID = [](Muon const* m) { return MuonMedium(m); };

  BuildModule(SimpleFilter<Muon>("FSRMuonFilter")
      .set_input_label("fsr_muons").set_min(0)
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > 5.        &&
                fabs(m->eta())          < 2.4   &&
                fabs(m->dxy_vertex())   < muon_dxy   &&
                fabs(m->dz_vertex())    < muon_dz   &&
                MuonID(m);

      }));



     
  BuildModule(CompositeProducer<Electron, Electron>("ZEEPairProducer")
      .set_input_label_first("sel_lead_electrons")
      .set_input_label_second("sel_electrons")
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));

}

// --------------------------------------------------------------------------
// Tag and probe ZEE Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildTPZEEPairs() {

  BuildModule(CopyCollection<Electron>("CopyToSelectedElectrons",
      js["electrons"].asString(), "sel_electrons"));

  std::function<bool(Electron const*)> ElecID = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); }; //set to false for 80wp //ElectronHTTIdFall17V2 //ElectronHTTIsoIdFall17V2

  BuildModule(SimpleFilter<Electron>("ProbeFilter")
      .set_input_label("sel_electrons").set_min(2)
      .set_max(2)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > elec_pt    &&
                fabs(e->eta())          < elec_eta   &&
                fabs(e->dxy_vertex())   < elec_dxy   &&
                fabs(e->dz_vertex())    < elec_dz    ;

      }));

  BuildModule(CompositeProducer<Electron, Electron>("ZEEPairProducer")
      .set_input_label_first("sel_electrons")
      .set_input_label_second("sel_electrons")
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));

}

// --------------------------------------------------------------------------
// ZMM Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildZMMPairs() {
 
 if (mu_scale_mode > 0){
   BuildModule(HTTMuonEnergyScale("MuonEnergyScaleCorrection")
      .set_input_label("muons")
      .set_far_endcap(muon_shift_farendcap)
      .set_near_endcap(muon_shift_nearendcap)
      .set_barrel(muon_shift_barrel)
      );
 }
 
 BuildModule(CopyCollection<Muon>("CopyToSelectedMuons",
      js["muons"].asString(), "sel_muons"));

  std::function<bool(Muon const*)> MuonID = [](Muon const* m) { return MuonMedium(m); };
   
  BuildModule(SimpleFilter<Muon>("MuonFilter")
      .set_input_label("sel_muons").set_min(2)
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > 20.        &&
                fabs(m->eta())          < 2.4   &&
                fabs(m->dxy_vertex())   < muon_dxy   &&
                fabs(m->dz_vertex())    < muon_dz   &&
                MuonID(m);

      }));
 
 
  BuildModule(CopyCollection<Muon>("CopyToSelectedMuons",
      "sel_muons", "sel_lead_muons"));
  
  BuildModule(SimpleFilter<Muon>("LeadMuonFilter")
      .set_input_label("sel_lead_muons").set_min(1)
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > muon_pt    &&
                fabs(m->eta())          < muon_eta   &&
                fabs(m->dxy_vertex())   < muon_dxy   &&
                fabs(m->dz_vertex())    < muon_dz   &&
                MuonID(m);
      }));

  BuildModule(CopyCollection<Electron>("CopyToFSRElectrons",
      js["electrons"].asString(), "fsr_electrons"));
  
  std::function<bool(Electron const*)> ElecID = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); };

  BuildModule(SimpleFilter<Electron>("FSRElectronIDFilter")
      .set_input_label("fsr_electrons").set_min(0)
      .set_predicate([=](Electron const* e) {
        return ElecID(e);
      }));

  BuildModule(SimpleFilter<Electron>("FSRElectronFilter")
      .set_input_label("fsr_electrons").set_min(0)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > 10.        &&
                fabs(e->eta())          < 2.4        &&
                fabs(e->dxy_vertex())   < elec_dxy   &&
                fabs(e->dz_vertex())    < elec_dz ;

      }));

  BuildModule(CopyCollection<Muon>("CopyToFSRMuons",
      js["muons"].asString(), "fsr_muons"));


  BuildModule(SimpleFilter<Muon>("FSRMuonFilter")
      .set_input_label("fsr_muons").set_min(0)
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > 5.        &&
                fabs(m->eta())          < 2.4   &&
                fabs(m->dxy_vertex())   < muon_dxy   &&
                fabs(m->dz_vertex())    < muon_dz   &&
                MuonID(m);

      }));
  
  BuildModule(CompositeProducer<Muon, Muon>("ZMMPairProducer")
      .set_input_label_first("sel_lead_muons")
      .set_input_label_second("sel_muons")
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));

  // Filter cases where muon 1 is softer than muon 2 which can happen due to different eta cuts for sel_muons and sel_leadmuons 
  BuildModule(SimpleFilter<CompositeCandidate>("ZMMPairFirstLeadingSecondFilter")
      .set_input_label("ditau").set_min(1)
      .set_predicate([=](CompositeCandidate const* c) {
        return PairFirstLeadingSecond(c);
      }));
 
  
  if(js["baseline"]["do_faketaus"].asBool()){
    BuildModule(SimpleFilter<Tau>("TauFilter")
        .set_input_label(js["taus"].asString()).set_min(0)
        .set_predicate([=](Tau const* t) {
          return  t->pt()                     > 15          &&
                  fabs(t->eta())              <  2.3        &&
                  fabs(t->lead_dz_vertex())   <  0.2        &&
                  fabs(t->charge())           == 1          &&
                  t->GetTauID("decayModeFindingNewDMs") > 0.5 && (t->decay_mode()<3 || t->decay_mode()>9);
    
        }));  
    BuildModule(OverlapFilter<Tau,CompositeCandidate>("TauMuonOverlapFilter")
      .set_input_label(js["taus"].asString())
      .set_reference_label("ditau")
      .set_min_dr(0.5));
  }
  
}

// --------------------------------------------------------------------------
// Tag and probe ZMM Pair Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildTPZMMPairs() {
 
 BuildModule(CopyCollection<Muon>("CopyToSelectedMuons",
      js["muons"].asString(), "sel_muons"));

  std::function<bool(Muon const*)> MuonID = [](Muon const* m) { return MuonMedium(m); };

  BuildModule(SimpleFilter<Muon>("ProbeFilter")
      .set_input_label("sel_muons").set_min(2).set_max(2)
      .set_predicate([=](Muon const* m) {
        return  m->pt()                 > muon_pt    &&
                fabs(m->eta())          < muon_eta   &&
                fabs(m->dxy_vertex())   < muon_dxy   &&
                fabs(m->dz_vertex())    < muon_dz  ;

      }));
  
  BuildModule(CompositeProducer<Muon, Muon>("ZMMPairProducer")
      .set_input_label_first("sel_muons")
      .set_input_label_second("sel_muons")
      .set_candidate_name_first("lepton1")
      .set_candidate_name_second("lepton2")
      .set_output_label("ditau"));

  std::function<bool(Electron const*)> ElecID = [](Electron const* e) { return ElectronHTTIdFall17V2(e, true); };

  BuildModule(CopyCollection<Electron>("CopyToFSRElectrons",
      js["electrons"].asString(), "fsr_electrons"));

  BuildModule(SimpleFilter<Electron>("FSRElectronIDFilter")
      .set_input_label("fsr_electrons").set_min(0)
      .set_predicate([=](Electron const* e) {
        return ElecID(e);
      }));

  BuildModule(SimpleFilter<Electron>("FSRElectronFilter")
      .set_input_label("fsr_electrons").set_min(0)
      .set_predicate([=](Electron const* e) {
        return  e->pt()                 > 10.        &&
                fabs(e->eta())          < 2.5        &&
                fabs(e->dxy_vertex())   < elec_dxy   &&
                fabs(e->dz_vertex())    < elec_dz ;

      }));
}


// --------------------------------------------------------------------------
// Tau Selection Sequence
// --------------------------------------------------------------------------
void HTTSequence::BuildTauSelection(){
  
  // filter taus first with loose pT cut - this avoids running more time consuming parts of the code for events with taus that just wont pass the offline cuts anyway 
  double loose_tau_pt = tau_pt*0.8;
  BuildModule(SimpleFilter<Tau>("TauFilterNewDMLoosePT")
     .set_input_label(js["taus"].asString()).set_min(min_taus)
     .set_predicate([=](Tau const* t) {
       return  t->pt()                     >  loose_tau_pt &&
               fabs(t->eta())              <  tau_eta    &&
               fabs(t->lead_dz_vertex())   <  tau_dz     &&
               fabs(t->charge())           == 1          &&
               t->GetTauID("decayModeFindingNewDMs") > 0.5 && (t->decay_mode()<3 || t->decay_mode()>9) &&
               //((t->HasTauID("byVVVLooseDeepTau2017v2p1VSjet") && t->GetTauID("byVVVLooseDeepTau2017v2p1VSjet") > 0.5 && t->GetTauID("byVVVLooseDeepTau2017v2p1VSe") > 0.5 && t->GetTauID("byVLooseDeepTau2017v2p1VSmu") > 0.5) || 
               //(t->HasTauID("byVVVLooseDeepTau2018v2p5VSjet") && t->GetTauID("byVVVLooseDeepTau2018v2p5VSjet") > 0.5 && t->GetTauID("byVVVLooseDeepTau2018v2p5VSe") > 0.5 && t->GetTauID("byVLooseDeepTau2018v2p5VSmu") > 0.5));
               (t->HasTauID("byMediumDeepTau2018v2p5VSjet") && t->GetTauID("byMediumDeepTau2018v2p5VSjet") > 0.5 && t->GetTauID("byVVLooseDeepTau2018v2p5VSe") > 0.5 && t->GetTauID("byTightDeepTau2018v2p5VSmu") > 0.5);
     }));
 
 if (tau_scale_mode > 0 && !is_data){
    // Tau energy scale is applied to genuine taus only - this works by selecting taus matched to generator level hadronic taus and saving these as a collection of pointers. The the shift is then applied to this collection which in turn shifts the tau energy for all corresponding taus in the origional tau collection (this will only work if both collections of taus are stored as pointers!)
    BuildModule(HTTGenMatchSelector<Tau>("HTTGenMatchSelector")
      .set_input_vec_label(js["taus"].asString())
      .set_output_vec_label("genmatched_taus")
      .set_add_gen(true)
      .set_gen_match(mcorigin::tauHad));

    BuildModule(CopyCollection<Tau>("CopyTauHadTo1Prong0Pi",
      "genmatched_taus", "genmatched_taus_1prong0pi0"));
    
    BuildModule(CopyCollection<Tau>("CopyTauHadTo1Prong1Pi",
      "genmatched_taus", "genmatched_taus_1prong1pi0"));

    BuildModule(CopyCollection<Tau>("CopyTauHadTo3Prong0Pi",
      "genmatched_taus", "genmatched_taus_3prong0pi0"));

    BuildModule(CopyCollection<Tau>("CopyTauHadTo3Prong1Pi",
      "genmatched_taus", "genmatched_taus_3prong1pi0"));

    BuildModule(SimpleFilter<Tau>("1Prong0PiTauHadFilter")
      .set_input_label("genmatched_taus_1prong0pi0")
      .set_predicate([=](Tau const* t) {
        return  t->decay_mode() == 0;
      }));
    
    BuildModule(SimpleFilter<Tau>("1Prong1PiTauHadFilter")
      .set_input_label("genmatched_taus_1prong1pi0")
      .set_predicate([=](Tau const* t) {
        return  t->decay_mode() == 1;
      }));

    BuildModule(SimpleFilter<Tau>("3Prong0PiTauHadFilter")
      .set_input_label("genmatched_taus_3prong0pi0")
      .set_predicate([=](Tau const* t) {
        return  t->decay_mode() == 10;
      }));
    
    BuildModule(SimpleFilter<Tau>("3Prong1PiTauHadFilter")
      .set_input_label("genmatched_taus_3prong1pi0")
      .set_predicate([=](Tau const* t) {
        return  t->decay_mode() == 11;
      }));
 
    BuildModule(EnergyShifter<Tau>("TauEnergyShifter1prong0pi0")
    .set_input_label("genmatched_taus_1prong0pi0")
    .set_save_shifts(true) 
    .set_shift_label("scales_taues_1prong0pi0") 
    .set_shift(tau_shift_1prong0pi0)
    .set_shift_func(tau_shift_func_1prong0pi0));

    BuildModule(EnergyShifter<Tau>("TauEnergyShifter1prong1pi0")
    .set_input_label("genmatched_taus_1prong1pi0")
    .set_save_shifts(true) 
    .set_shift_label("scales_taues_1prong1pi0") 
    .set_shift(tau_shift_1prong1pi0)
    .set_shift_func(tau_shift_func_1prong1pi0));

    BuildModule(EnergyShifter<Tau>("TauEnergyShifter3prong0pi0")
    .set_input_label("genmatched_taus_3prong0pi0")
    .set_save_shifts(true) 
    .set_shift_label("scales_taues_3prong0pi0") 
    .set_shift(tau_shift_3prong0pi0)
    .set_shift_func(tau_shift_func_3prong0pi0));

    BuildModule(EnergyShifter<Tau>("TauEnergyShifter3prong1pi0")
    .set_input_label("genmatched_taus_3prong1pi0")
    .set_save_shifts(true)
    .set_shift_label("scales_taues_3prong1pi0")
    .set_shift(tau_shift_3prong1pi0)
    .set_shift_func(tau_shift_func_3prong1pi0));
  }
  // energy resolutions shifting
  if (!is_data){
    if(tau_res_shift!=0.) {
      BuildModule(EnergyResShifter<Tau>("TauEnergyResShifter")
      .set_input_label("genmatched_taus")
      .set_input_label_gen("genmatched_taus_gen")
      .set_save_shifts(true)
      .set_shift_label("res_taus")
      .set_shift(tau_res_shift));
    }
 }
 if (!is_data){
   BuildModule(HTTGenMatchSelector<Tau>("FakeEGenMatchSelector")
     .set_input_vec_label(js["taus"].asString())
     .set_output_vec_label("fakeE_genmatched_taus")
     .set_gen_match(mcorigin::promptE));
   
   BuildModule(CopyCollection<Tau>("CopyTo1Prong0PiFakeE",
     "fakeE_genmatched_taus", "fakeE_genmatched_taus_0pi"));
   
   BuildModule(CopyCollection<Tau>("CopyTo1Prong1PiFakeE",
     "fakeE_genmatched_taus", "fakeE_genmatched_taus_1pi"));

   BuildModule(CopyCollection<Tau>("CopyTo1Prong0PiFakeE_EC",
     "fakeE_genmatched_taus", "fakeE_genmatched_taus_0pi_endcap"));

   BuildModule(CopyCollection<Tau>("CopyTo1Prong1PiFakeE_EC",
     "fakeE_genmatched_taus", "fakeE_genmatched_taus_1pi_endcap"));
   
   BuildModule(SimpleFilter<Tau>("1Prong0PiTauFilterFakeE")
     .set_input_label("fakeE_genmatched_taus_0pi")
     .set_predicate([=](Tau const* t) {
       return  t->decay_mode() == 0 && fabs(t->eta()) < 1.5;
     }));
   
   BuildModule(SimpleFilter<Tau>("1Prong1PiTauFilterFakeE")
     .set_input_label("fakeE_genmatched_taus_1pi")
     .set_predicate([=](Tau const* t) {
       return  t->decay_mode() == 1 && fabs(t->eta()) < 1.5;
     }));
    
   BuildModule(SimpleFilter<Tau>("1Prong0PiEndCapTauFilterFakeE")
     .set_input_label("fakeE_genmatched_taus_0pi_endcap")
     .set_predicate([=](Tau const* t) {
       return  t->decay_mode() == 0 && fabs(t->eta()) >= 1.5;
     }));

   BuildModule(SimpleFilter<Tau>("1Prong1PiEndCapTauFilterFakeE")
     .set_input_label("fakeE_genmatched_taus_1pi_endcap")
     .set_predicate([=](Tau const* t) {
       return  t->decay_mode() == 1 && fabs(t->eta()) >= 1.5;
     }));

   BuildModule(EnergyShifter<Tau>("FakeE1Prong0PiEnergyShifter")
   .set_input_label("fakeE_genmatched_taus_0pi")
   .set_save_shifts(true)
   .set_shift_label("scales_efaketaues_1prong0pi0")
   .set_shift(fakeE_tau_shift_0pi));
   
   BuildModule(EnergyShifter<Tau>("FakeE1Prong1PiEnergyShifter")
   .set_input_label("fakeE_genmatched_taus_1pi")
   .set_save_shifts(true)
   .set_shift_label("scales_efaketaues_1prong1pi0")
   .set_shift(fakeE_tau_shift_1pi));

   BuildModule(EnergyShifter<Tau>("FakeE1Prong0PiEndCapEnergyShifter")
   .set_input_label("fakeE_genmatched_taus_0pi_endcap")
   .set_save_shifts(true)
   .set_shift_label("scales_efaketaues_1prong0pi0_endcap")
   .set_shift(fakeE_tau_shift_0pi_endcap));

   BuildModule(EnergyShifter<Tau>("FakeE1Prong1PiEndCapEnergyShifter")
   .set_input_label("fakeE_genmatched_taus_1pi_endcap")
   .set_save_shifts(true)
   .set_shift_label("scales_efaketaues_1prong1pi0_endcap")
   .set_shift(fakeE_tau_shift_1pi_endcap));

   BuildModule(HTTGenMatchSelector<Tau>("FakeMuGenMatchSelector")
     .set_input_vec_label(js["taus"].asString())
     .set_output_vec_label("fakeMu_genmatched_taus")
     .set_gen_match(mcorigin::promptMu));
  
   BuildModule(HTTGenMatchSelector<Tau>("FakeJGenMatchSelector")
     .set_input_vec_label(js["taus"].asString())
     .set_output_vec_label("fakeJ_genmatched_taus")
     .set_gen_match(mcorigin::fake));
 
   BuildModule(CopyCollection<Tau>("CopyTo1Prong0PiFakeMu",
     "fakeMu_genmatched_taus", "fakeMu_genmatched_taus_0pi"));
   
   BuildModule(CopyCollection<Tau>("CopyTo1Prong1PiFakeMu",
     "fakeMu_genmatched_taus", "fakeMu_genmatched_taus_1pi"));
   
   BuildModule(SimpleFilter<Tau>("1Prong0PiTauFilterFakeMu")
     .set_input_label("fakeMu_genmatched_taus_0pi")
     .set_predicate([=](Tau const* t) {
       return  t->decay_mode() == 0;
     }));
   
   BuildModule(SimpleFilter<Tau>("1Prong1PiTauFilterFakeMu")
     .set_input_label("fakeMu_genmatched_taus_1pi")
     .set_predicate([=](Tau const* t) {
       return  t->decay_mode() == 1;
     }));
    
   BuildModule(EnergyShifter<Tau>("FakeMu1Prong0PiEnergyShifter")
   .set_input_label("fakeMu_genmatched_taus_0pi")
   .set_save_shifts(true)
   .set_shift_label("scales_mufaketaues_1prong0pi0")
   .set_shift(fakeMu_tau_shift_0pi));
   
   BuildModule(EnergyShifter<Tau>("FakeMu1Prong1PiEnergyShifter")
   .set_input_label("fakeMu_genmatched_taus_1pi")
   .set_save_shifts(true)
   .set_shift_label("scales_mufaketaues_1prong1pi0")
   .set_shift(fakeMu_tau_shift_1pi));

 }

 if(fakeJ_tau_shift!=1.) {
   //fakeJ_genmatched_taus
   BuildModule(EnergyShifter<Tau>("FakeJEnergyShifter")
   .set_input_label("fakeJ_genmatched_taus")
   .set_save_shifts(true)
   .set_shift_label("scales_jfaketaues")
   .set_shift(fakeJ_tau_shift));
 } 

 BuildModule(SimpleFilter<Tau>("TauFilterNewDM")
    .set_input_label(js["taus"].asString()).set_min(min_taus)
    .set_predicate([=](Tau const* t) {
      return  t->pt()                     >  tau_pt     &&
              fabs(t->eta())              <  tau_eta    &&
              fabs(t->lead_dz_vertex())   <  tau_dz     &&
              fabs(t->charge())           == 1          &&
              t->GetTauID("decayModeFindingNewDMs") > 0.5 && (t->decay_mode()<3 || t->decay_mode()>9) &&
               ((t->HasTauID("byVVVLooseDeepTau2017v2p1VSjet") && t->GetTauID("byVVVLooseDeepTau2017v2p1VSjet") > 0.5 && t->GetTauID("byVVVLooseDeepTau2017v2p1VSe") > 0.5 && t->GetTauID("byVLooseDeepTau2017v2p1VSmu") > 0.5) || 
               (t->HasTauID("byVVVLooseDeepTau2018v2p5VSjet") && t->GetTauID("byVVVLooseDeepTau2018v2p5VSjet") > 0.5 && t->GetTauID("byVVVLooseDeepTau2018v2p5VSe") > 0.5 && t->GetTauID("byVLooseDeepTau2018v2p5VSmu") > 0.5));

    }));

 }

 void HTTSequence::BuildDiElecVeto() {

   BuildModule(CopyCollection<Electron>("CopyToVetoElecs",
       js["electrons"].asString(), "veto_elecs"));

   BuildModule(GenericModule("VetoElecIDIsoFilter")
      .set_function([=](ic::TreeEvent *event){
      EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
      std::vector<Electron*> & vec = event->GetPtrVec<Electron>("veto_elecs");
      ic::erase_if(vec,!boost::bind(VetoElectronIDFall17,_1, eventInfo->jet_rho())); //lepton_rho
      ic::erase_if(vec,!boost::bind(PF03EAElecIsolation, _1, eventInfo->jet_rho(), 0.3));
      return 0;
   }));    

   SimpleFilter<Electron> vetoElecFilter = SimpleFilter<Electron>("VetoElecFilter")
      .set_input_label("veto_elecs");
   vetoElecFilter.set_predicate([=](Electron const* e) {
     return  e->pt()                 > veto_dielec_pt    &&
            fabs(e->eta())          < veto_dielec_eta   &&
            fabs(e->dxy_vertex())   < veto_dielec_dxy   &&
            fabs(e->dz_vertex())    < veto_dielec_dz;
   });

   BuildModule(vetoElecFilter);

   BuildModule(OneCollCompositeProducer<Electron>("VetoElecPairProducer")
      .set_input_label("veto_elecs").set_output_label("elec_veto_pairs")
      .set_candidate_name_first("elec1").set_candidate_name_second("elec2"));

   HTTFilter<CompositeCandidate> vetoElecPairFilter = HTTFilter<CompositeCandidate>("VetoElecPairFilter")
      .set_veto_name("dielec_veto")
      .set_no_filter(true)
      .set_input_label("elec_veto_pairs").set_min(0).set_max(0)
      .set_predicate([=](CompositeCandidate const* c) {
        return  c->DeltaR("elec1", "elec2") > 0.15 &&
                c->charge() == 0;
      });
   BuildModule(vetoElecPairFilter);
 }

 void HTTSequence::BuildDiMuonVeto() {

  BuildModule(CopyCollection<Muon>("CopyToVetoMuons",
      js["muons"].asString(), "veto_muons"));

  SimpleFilter<Muon> vetoMuonFilter = SimpleFilter<Muon>("VetoMuonFilter")
      .set_input_label("veto_muons");
  vetoMuonFilter.set_predicate([=](Muon const* m) {
    return  m->pt()                 > veto_dimuon_pt    &&
            fabs(m->eta())          < veto_dimuon_eta   &&
            fabs(m->dxy_vertex())   < veto_dimuon_dxy   &&
            fabs(m->dz_vertex())    < veto_dimuon_dz    &&
            m->is_global()                    &&
            m->is_tracker()                   &&
            m->is_pf()                        &&
            PF04IsolationVal(m, 0.5,0) < 0.3;
  });

  BuildModule(vetoMuonFilter);

  BuildModule(OneCollCompositeProducer<Muon>("VetoMuonPairProducer")
      .set_input_label("veto_muons").set_output_label("muon_veto_pairs")
      .set_candidate_name_first("muon1").set_candidate_name_second("muon2"));

  HTTFilter<CompositeCandidate> vetoMuonPairFilter = HTTFilter<CompositeCandidate>("VetoMuonPairFilter")
      .set_input_label("muon_veto_pairs").set_min(0).set_max(0)
      .set_veto_name("dimuon_veto")
      .set_no_filter(true)
      .set_predicate([=](CompositeCandidate const* c){
        return c->DeltaR("muon1", "muon2") > 0.15 &&
	       c->charge() == 0;
      });	
	
  BuildModule(vetoMuonPairFilter);

 }

 void HTTSequence::BuildExtraElecVeto(){
 
   BuildModule(CopyCollection<Electron>("CopyToExtraElecs",
       js["electrons"].asString(), "extra_elecs"));
 
   BuildModule(GenericModule("ExtraElecIsoFilter")
     .set_function([=](ic::TreeEvent *event){
        EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
        std::vector<Electron*> & vec = event->GetPtrVec<Electron>("extra_elecs");
        ic::erase_if(vec,!boost::bind(PF03EAElecIsolation, _1, eventInfo->jet_rho(), 0.3)); //lepton_rho
        return 0;
   }));
 
   HTTFilter<Electron> extraElecFilter = HTTFilter<Electron>("ExtraElecFilter")
       .set_input_label("extra_elecs")
       .set_veto_name("extra_elec_veto")
       .set_min(0).set_max(js["baseline"]["max_extra_elecs"].asUInt());
   extraElecFilter.set_no_filter(true);
   extraElecFilter.set_predicate([=](Electron const* e) {
      return  e->pt()                 > veto_elec_pt    &&
              fabs(e->eta())          < veto_elec_eta   &&
              fabs(e->dxy_vertex())   < veto_elec_dxy   &&
              fabs(e->dz_vertex())    < veto_elec_dz    &&
              ElectronHTTIdFall17V2(e, true);
    });
   BuildModule(extraElecFilter);
 }
 
 
 void HTTSequence::BuildExtraMuonVeto(){
 
   BuildModule(CopyCollection<Muon>("CopyToExtraMuons",
       js["muons"].asString(), "extra_muons"));
 
   HTTFilter<Muon> extraMuonFilter = HTTFilter<Muon>("ExtraMuonFilter")
       .set_input_label("extra_muons")
       .set_veto_name("extra_muon_veto")
       .set_min(0).set_max(js["baseline"]["max_extra_muons"].asUInt());
   std::function<bool(Muon const*)> ExtraMuonID = [](Muon const* m) {return MuonMedium(m); };
   extraMuonFilter.set_no_filter(true);
   extraMuonFilter.set_predicate([=](Muon const* m) {
     return  m->pt()                 > veto_muon_pt    &&
             fabs(m->eta())          < veto_muon_eta   &&
             fabs(m->dxy_vertex())   < veto_muon_dxy   &&
             fabs(m->dz_vertex())    < veto_muon_dz    &&
             ExtraMuonID(m)                     &&
             PF04IsolationVal(m, 0.5,0) < 0.3;
   });
   BuildModule(extraMuonFilter);
 
 
 }
}
