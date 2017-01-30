#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTRun2AnalysisTools.h"
#include <iostream>
#include <vector>
#include <map>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/program_options.hpp"
#include "boost/range/algorithm.hpp"
#include "boost/range/algorithm_ext.hpp"
#include "boost/filesystem.hpp"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/SimpleParamParser.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnRootTools.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/th1fmorph.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTau/interface/HTTConfig.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TEfficiency.h"
#include "TEntryList.h"
#include "TMath.h"
#include "TH1.h"
#include "TLegend.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooPlot.h"

namespace ic {

  HTTRun2Analysis::HTTRun2Analysis(ic::channel ch, std::string year, int verbosity, bool is_sm) : ch_(ch), year_(year), verbosity_(verbosity), is_sm_(is_sm) {
    lumi_ = 1.;
    do_ss_ = false;
    qcd_os_ss_factor_ = 1.06;
    /*if(ch_ == channel::et){
      w_os_ss_factor_ = 4.09;
    } else if(ch_ == channel::mt){
      w_os_ss_factor_ = 3.76;
    } else w_os_ss_factor_=1.0;*/
    using boost::range::push_back;

    //Sample splitting
    if (ch_ == channel::et){
      alias_map_["ztt_sel"] = "(gen_match_2==5)";
      alias_map_["zl_sel"] = "(gen_match_2<5)";
      alias_map_["zj_sel"] = "(gen_match_2==6)";
    }
    if (ch_ == channel::mt){
      alias_map_["ztt_sel"] = "(gen_match_2==5)";
      alias_map_["zl_sel"] = "(gen_match_2<5)";
      alias_map_["zj_sel"] = "(gen_match_2==6)";
    }
    if (ch_ == channel::tt){
      alias_map_["ztt_sel"] = "(gen_match_1==5&&gen_match_2==5)";
      alias_map_["zl_sel"] = "(gen_match_2<6&&gen_match_1<6&&!(gen_match_1==5&&gen_match_2==5))";
      alias_map_["zj_sel"] = "(gen_match_2==6||gen_match_1==6)";
    }
    if (ch_ == channel::em ){
      alias_map_["ztt_sel"] = "(gen_match_1>2 && gen_match_2>3)";
      alias_map_["zll_sel"] = "(gen_match_1<3 || gen_match_2<4)";    
    }
    if (ch_ == channel::zee||ch_==channel::zmm||ch_ == channel::tpzee||ch_==channel::tpzmm||ch_ ==channel::wmnu){
      alias_map_["ztt_sel"] = "";
      alias_map_["zll_sel"] = "1.";
    }
    
    alias_map_["variable_cat"]      = "1";
    alias_map_["inclusive"]                 = "1";
    
    if (ch_ == channel::et || ch_ == channel::mt){
      alias_map_["vbf_run1"] = "(n_jets>=2 && jdeta>3.5 && mjj>500 && n_bjets==0 && n_jetsingap==0)";
      alias_map_["vbf_tight_run1"] = "(n_jets>=2 && jdeta>4. && mjj>700 && n_bjets==0 && n_jetsingap==0 && pt_tt>100)";
      alias_map_["vbf_loose_run1"] = "(!"+alias_map_["vbf_tight_run1"]+" && n_jets>=2 && jdeta>3.5 && mjj>500 && n_bjets==0 && n_jetsingap==0)";    
      if(ch_ == channel::et){
        alias_map_["1jet_loose_run1_et"] = "(!"+alias_map_["vbf_run1"]+" && n_jets>=1 && n_bjets==0 && met>30 && pt_2<=45)";
        alias_map_["1jet_tight_run1_et"] = "(!"+alias_map_["vbf_run1"]+" && n_jets>=1 && n_bjets==0 && met>30 && pt_tt>100 && pt_2>45)";
      } else if(ch_ == channel::mt){
        alias_map_["1jet_loose_run1"] = "(!"+alias_map_["vbf_run1"]+" && n_jets>=1 && n_bjets==0 && pt_2<=45)";
        alias_map_["1jet_medium_run1"] = "(!"+alias_map_["vbf_run1"]+" && n_jets>=1 && n_bjets==0 && pt_tt<=100 && pt_2>45)";
        alias_map_["1jet_tight_run1"] = "(!"+alias_map_["vbf_run1"]+" && n_jets>=1 && n_bjets==0 && pt_tt>100 && pt_2>45)";
      }
      alias_map_["0jet_loose_run1"] = "(n_jets==0 && pt_2<=45)";
      alias_map_["0jet_tight_run1"] = "(n_jets==0 && pt_2>45)";
    }
    else if (ch_ == channel::tt){
      alias_map_["vbf_run1_tt"] = "(n_jets>=2 && jdeta>3.5 && mjj>500 && n_bjets==0 && pt_tt>100 && n_jetsingap==0)";
      alias_map_["1jet_run1_tt"] = "(!"+alias_map_["vbf_run1_tt"]+" && n_jets>=1 && n_bjets==0 && pt_tt>100)";
      alias_map_["1jet_loose_run1_tt"] = "(!"+alias_map_["vbf_run1_tt"]+" && n_jets>=1 && n_bjets==0 && pt_tt>100 && pt_tt<=170)";
      alias_map_["1jet_tight_run1_tt"] = "(!"+alias_map_["vbf_run1_tt"]+" && n_jets>=1 && n_bjets==0 && pt_tt>170 )";
      alias_map_["0jet_run1_tt"] = "(!"+alias_map_["vbf_run1_tt"]+" && !"+alias_map_["1jet_run1_tt"]+" && n_jets>=0 && n_bjets==0)";
    }


    if (ch_ == channel::et || ch_ == channel::mt) {
      // SM Categories
      //SM Categories
        
      alias_map_["inclusive"]         = "1";
      alias_map_["nojet"]           ="(n_jets==0)";
      alias_map_["njet"]           ="(n_jets>0)";
      alias_map_["notwoprong"]       ="(tau_decay_mode_2!=6&&tau_decay_mode_2!=5)";
      alias_map_["baseline"] = "(iso_1<0.1 && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto)";
      
      if(is_sm){
        if(ch_ == channel::mt && year_.find("6")!=year_.npos) alias_map_["baseline"] = "(iso_1<0.15 && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto)";
        if(ch_ == channel::et && year_.find("6")!=year_.npos) alias_map_["baseline"] = "(iso_1<0.1  && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto)";    
      } else{
        if(ch_ == channel::mt && year_.find("6")!=year_.npos) alias_map_["baseline"] = "(iso_1<0.15 && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto)";
        if(ch_ == channel::et && year_.find("6")!=year_.npos) alias_map_["baseline"] = "(iso_1<0.1  && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto)";
      }
//      alias_map_["baseline"]          = "1";
      alias_map_["incvlelm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_vloose_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["incvletm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_vloose_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["inclelm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_loose_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["incletm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_loose_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["incmelm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_medium_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["incmetm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_medium_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["inctelm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_tight_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["inctetm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_tight_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["incvtelm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_vtight_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["incvtetm"]         = "(iso_1<0.1&&iso_2<1.5 && antie_vtight_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["dbriso1p5"]        = "(iso_1<0.1&&iso_2<1.5&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["dbriso1p75"]        = "(iso_1<0.1&&iso_2<1.75&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["dbriso1p0"]        = "(iso_1<0.1&&iso_2<1&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["dbriso2p0"]        = "(iso_1<0.1&&iso_2<2&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["dbriso2p5"]        = "(iso_1<0.1&&iso_2<2.5&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["mvaiso0"]        = "(iso_1<0.1&&iso_mva_new_2>0&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["mvaisom0p2"]        = "(iso_1<0.1&&iso_mva_new_2>-0.2&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["mvaisop0p2"]        = "(iso_1<0.1&&iso_mva_new_2>0.2&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["mvaisom0p1"]        = "(iso_1<0.1&&iso_mva_new_2>-0.1&&antiele_2&&antimu_2&&!leptonveto)";
      alias_map_["mvaisop0p1"]        = "(iso_1<0.1&&iso_mva_new_2>0.1&&antiele_2&&antimu_2&&!leptonveto)";

      alias_map_["inclusivenolv"]         = "(iso_1<0.1&&iso_2<1.5 && antiele_2 && antimu_2)";
      alias_map_["btag"] = "(n_jets<=1 && n_bjets>=1)";
      alias_map_["nobtag"] = "n_bjets==0";
      // Categories for iso studies
      alias_map_["incnotauiso"]      = "(iso_1<0.1 && antiele_2 && antimu_2 && !leptonveto&&"+alias_map_["notwoprong"]+")";
      alias_map_["dbloose"]          = "(db_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["dbmedium"]         = "(db_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["dbtight"]         = "(db_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["btagdbloose"]          = "(db_loose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagdbmedium"]         = "(db_medium_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagdbtight"]         = "(db_tight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["nobtagdbloose"]          = "(db_loose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagdbmedium"]         = "(db_medium_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagdbtight"]         = "(db_tight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["puwloose"]          = "(puw_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["puwmedium"]         = "(puw_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["puwtight"]         = "(puw_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldvloose"]    = "(mvadbold_vloose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldloose"]     = "(mvadbold_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldmedium"]    = "(mvadbold_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldtight"]     = "(mvadbold_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldvtight"]     = "(mvadbold_vtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldvvtight"]    = "(mvadbold_vvtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["btagmvadboldvloose"]    = "(mva_olddm_vloose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldloose"]     = "(mva_olddm_loose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldmedium"]    = "(mva_olddm_medium_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldtight"]     = "(mva_olddm_tight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldvtight"]     = "(mva_olddm_vtight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldvvtight"]    = "(mva_olddm_vvtight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["nobtagmvadboldvloose"]    = "(mva_olddm_vloose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldloose"]     = "(mva_olddm_loose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldmedium"]    = "(mva_olddm_medium_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldtight"]     = "(mva_olddm_tight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldvtight"]     = "(mva_olddm_vtight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldvvtight"]    = "(mva_olddm_vvtight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["mvadbnewvloose"]    = "(mvadbnew_vloose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewloose"]     = "(mvadbnew_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewmedium"]    = "(mvadbnew_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewtight"]     = "(mvadbnew_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewvtight"]     = "(mvadbnew_vtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewvvtight"]    = "(mvadbnew_vvtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldvloose"]    = "(mvapwold_vloose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldloose"]     = "(mvapwold_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldmedium"]    = "(mvapwold_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldtight"]     = "(mvapwold_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldvtight"]     = "(mvapwold_vtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldvvtight"]    = "(mvapwold_vvtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewvloose"]    = "(mvapwnew_vloose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewloose"]     = "(mvapwnew_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewmedium"]    = "(mvapwnew_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewtight"]     = "(mvapwnew_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewvtight"]     = "(mvapwnew_vtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewvvtight"]    = "(mvapwnew_vvtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["incnoiso"]         = "(iso_2<1.5 && antiele_2 && antimu_2 && !leptonveto)";
      alias_map_["incnoisonolv"]     = "(iso_2<1.5 && antiele_2 && antimu_2)";
      alias_map_["db03iso0p1"]            = "(iso_1_db03<0.1&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p11"]            = "(iso_1_db03<0.11&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p12"]            = "(iso_1_db03<0.12&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p13"]            = "(iso_1_db03<0.13&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p14"]            = "(iso_1_db03<0.14&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p15"]            = "(iso_1_db03<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p16"]            = "(iso_1_db03<0.16&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p17"]            = "(iso_1_db03<0.17&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p18"]            = "(iso_1_db03<0.18&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p09"]           = "(iso_1_db03<0.09&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p08"]           = "(iso_1_db03<0.08&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03iso0p07"]           = "(iso_1_db03<0.07&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p1"]            = "(iso_1_ea03<0.1&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p11"]            = "(iso_1_ea03<0.11&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p12"]            = "(iso_1_ea03<0.12&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p13"]            = "(iso_1_ea03<0.13&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p14"]            = "(iso_1_ea03<0.14&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p15"]            = "(iso_1_ea03<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p16"]            = "(iso_1_ea03<0.16&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p17"]            = "(iso_1_ea03<0.17&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p18"]            = "(iso_1_ea03<0.18&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p09"]           = "(iso_1_ea03<0.09&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p08"]           = "(iso_1_ea03<0.08&&"+alias_map_["incnoiso"]+")";
      alias_map_["ea03iso0p07"]           = "(iso_1_ea03<0.07&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p1"]            = "(iso_1_db03allch<0.1&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p11"]            = "(iso_1_db03allch<0.11&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p12"]            = "(iso_1_db03allch<0.12&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p13"]            = "(iso_1_db03allch<0.13&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p14"]            = "(iso_1_db03allch<0.14&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p15"]            = "(iso_1_db03allch<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p16"]            = "(iso_1_db03allch<0.16&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p17"]            = "(iso_1_db03allch<0.17&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p18"]            = "(iso_1_db03allch<0.18&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p09"]           = "(iso_1_db03allch<0.09&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p08"]           = "(iso_1_db03allch<0.08&&"+alias_map_["incnoiso"]+")";
      alias_map_["db03allchiso0p07"]           = "(iso_1_db03allch<0.07&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p1"]            = "(iso_1_db04allch<0.1&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p11"]            = "(iso_1_db04allch<0.11&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p12"]            = "(iso_1_db04allch<0.12&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p13"]            = "(iso_1_db04allch<0.13&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p14"]            = "(iso_1_db04allch<0.14&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p15"]            = "(iso_1_db04allch<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p16"]            = "(iso_1_db04allch<0.16&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p17"]            = "(iso_1_db04allch<0.17&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p18"]            = "(iso_1_db04allch<0.18&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p09"]           = "(iso_1_db04allch<0.09&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p08"]           = "(iso_1_db04allch<0.08&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04allchiso0p07"]           = "(iso_1_db04allch<0.07&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p1"]            = "(iso_1_db04<0.1&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p11"]            = "(iso_1_db04<0.11&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p12"]            = "(iso_1_db04<0.12&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p13"]            = "(iso_1_db04<0.13&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p14"]            = "(iso_1_db04<0.14&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p15"]            = "(iso_1_db04<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p16"]            = "(iso_1_db04<0.16&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p17"]            = "(iso_1_db04<0.17&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p18"]            = "(iso_1_db04<0.18&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p09"]           = "(iso_1_db04<0.09&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p08"]           = "(iso_1_db04<0.08&&"+alias_map_["incnoiso"]+")";
      alias_map_["db04iso0p07"]           = "(iso_1_db04<0.07&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p1"]            = "(iso_1_trk03<0.1&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p11"]            = "(iso_1_trk03<0.11&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p12"]            = "(iso_1_trk03<0.12&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p13"]            = "(iso_1_trk03<0.13&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p14"]            = "(iso_1_trk03<0.14&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p15"]            = "(iso_1_trk03<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p16"]            = "(iso_1_trk03<0.16&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p17"]            = "(iso_1_trk03<0.17&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p18"]            = "(iso_1_trk03<0.18&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p09"]           = "(iso_1_trk03<0.09&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p08"]           = "(iso_1_trk03<0.08&&"+alias_map_["incnoiso"]+")";
      alias_map_["trk03iso0p07"]           = "(iso_1_trk03<0.07&&"+alias_map_["incnoiso"]+")";
      alias_map_["qcd_loose_shape"]         = "(iso_1>0.2 && iso_1<0.5 && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto)";
      alias_map_["qcd_vloose_shape"]         = "(iso_1>0.2 && iso_1<0.5 && mva_olddm_tight_2<0.5 && antiele_2 && antimu_2 && !leptonveto)";
      alias_map_["vbf"] = "(n_jets>=2 && n_jetsingap==0 && mjj>500 && jdeta>3.5)";
      alias_map_["1jet"] = "(!("+alias_map_["vbf"]+")"+"&& n_jets>=1 && n_bjets==0)";

      alias_map_["btagnotwoprong"] = "(n_jets<=1 && n_bjets>=1&&"+alias_map_["notwoprong"]+")";
      alias_map_["btagpt20"] = "(n_lowpt_jets<=1 && n_bjets>=1)";

      alias_map_["nobtagnotwoprong"] = "(n_bjets==0&&"+alias_map_["notwoprong"]+")";
      //for making CSV control plot
      alias_map_["prebtag"] = "(n_jets<=1 && n_prebjets>=1)";
      //MSSM update analysis style categories:
      alias_map_["btaglow"] = "(n_jets<=1 && n_bjets>=1 && pt_2>30 && pt_2<45)";
      alias_map_["btaghigh"] = "(n_jets<=1 && n_bjets>=1 && pt_2>45)";
      alias_map_["btaghighnotwoprong"] = "(n_jets<=1 && n_bjets>=1 && pt_2>45&&"+alias_map_["notwoprong"]+")";
      alias_map_["nobtaglow"] = "n_bjets==0 && pt_2>30 && pt_2<45";
      alias_map_["nobtagmed"] = "n_bjets==0 && pt_2>45 && pt_2<60";
      alias_map_["nobtaghigh"] = "n_bjets==0 && pt_2>60";
      alias_map_["nobtaghighnotwoprong"] = "(n_bjets==0 && pt_2>60 &&"+alias_map_["notwoprong"]+")";
      alias_map_["2jet0tag"] = "(n_lowpt_jets>=2 && n_bjets==0)";
      alias_map_["2jet1tag"] = "(n_lowpt_jets>=2 && n_bjets==1)";
      alias_map_["2jet2tag"] = "(n_lowpt_jets>=2 && n_bjets>=2)";
      alias_map_["2jet0taghigh"] = "(n_jets>=2 && n_bjets==0)";
      alias_map_["2jet1taghigh"] = "(n_jets>=2 && n_bjets==1)";
      alias_map_["2jet2taghigh"] = "(n_jets>=2 && n_bjets>=2)";
    } else if (ch_ == channel::tt) {
      alias_map_["btag"] = "(n_jets<=1 && n_loose_bjets>=1)";
      alias_map_["nobtag"] = "n_loose_bjets==0";
      if(year_.find("6")!=year_.npos) alias_map_["btag"] = "(n_jets<=1 && n_bjets>=1)";
      if(year_.find("6")!=year_.npos) alias_map_["nobtag"] = "n_bjets==0";
      //alias_map_["btag"] = "(n_jets<=1 && n_bjets>=1)";
      //alias_map_["nobtag"] = "n_bjets==0";
      alias_map_["notwoprong"]      ="(tau_decay_mode_1!=5&&tau_decay_mode_2!=5&&tau_decay_mode_1!=6&&tau_decay_mode_2!=6)";
      alias_map_["incnotauiso"]          = "antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto";
      alias_map_["incvlelm"]         = "(iso_1<1&&iso_2<1 && antie_vloose_1>0 && antimu_loose_1>0 && antie_vloose_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["incvletm"]         = "(iso_1<1&&iso_2<1 && antie_vloose_1>0 && antimu_tight_1>0 && antie_vloose_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["inclelm"]         = "(iso_1<1&&iso_2<1 && antie_loose_1>0 && antimu_loose_1>0 && antie_loose_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["incletm"]         = "(iso_1<1&&iso_2<1 && antie_loose_1>0 && antimu_tight_1>0 && antie_loose_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["incmelm"]         = "(iso_1<1&&iso_2<1 && antie_medium_1>0 && antimu_loose_1>0 && antie_medium_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["incmetm"]         = "(iso_1<1&&iso_2<1 && antie_medium_1>0 && antimu_tight_1>0 && antie_medium_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["inctelm"]         = "(iso_1<1&&iso_2<1 && antie_tight_1>0 && antimu_loose_1>0 && antie_tight_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["inctetm"]         = "(iso_1<1&&iso_2<1 && antie_tight_1>0 && antimu_tight_1>0 && antie_tight_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["incvtelm"]         = "(iso_1<1&&iso_2<1 && antie_vtight_1>0 && antimu_loose_1>0 && antie_vtight_2>0 && antimu_loose_2>0 && !leptonveto)";
      alias_map_["incvtetm"]         = "(iso_1<1&&iso_2<1 && antie_vtight_1>0 && antimu_tight_1>0 && antie_vtight_2>0 && antimu_tight_2>0 && !leptonveto)";
      alias_map_["dbloose"]          = "(db_loose_1&&db_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["dbmedium"]         = "(db_medium_1&&db_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["dbtight"]         = "(db_tight_1&&db_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["puwloose"]          = "(puw_loose_1&&puw_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["btagdbloose"]          = "(db_loose_1&&db_loose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagdbmedium"]         = "(db_medium_1&&db_medium_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagdbtight"]         = "(db_tight_1&&db_tight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["nobtagdbloose"]          = "(db_loose_1&&db_loose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagdbmedium"]         = "(db_medium_1&&db_medium_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagdbtight"]         = "(db_tight_1&&db_tight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["puwmedium"]         = "(puw_medium_1&&puw_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["puwtight"]         = "(puw_tight_1&&puw_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["btagmvadboldvloose"]    = "(mva_olddm_vloose_1&&mva_olddm_vloose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldloose"]     = "(mva_olddm_loose_1&&mva_olddm_loose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldmedium"]    = "(mva_olddm_medium_1&&mva_olddm_medium_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldtight"]     = "(mva_olddm_tight_1&&mva_olddm_tight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldvtight"]     = "(mva_olddm_vtight_1&&mva_olddm_vtight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["btagmvadboldvvtight"]    = "(mva_olddm_vvtight_1&&mva_olddm_vvtight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["btag"]+")";
      alias_map_["nobtagmvadboldvloose"]    = "(mva_olddm_vloose_1&&mva_olddm_vloose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldloose"]     = "(mva_olddm_loose_1&&mva_olddm_loose_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldmedium"]    = "(mva_olddm_medium_1&&mva_olddm_medium_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldtight"]     = "(mva_olddm_tight_1&&mva_olddm_tight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldvtight"]     = "(mva_olddm_vtight_1&&mva_olddm_vtight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmvadboldvvtight"]    = "(mva_olddm_vvtight_1&&mva_olddm_vvtight_2&&"+alias_map_["incnotauiso"]+"&&"+alias_map_["nobtag"]+")";
      alias_map_["mvadbnewvloose"]    = "(mvadbnew_vloose_1&&mvadbnew_vloose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewloose"]     = "(mvadbnew_loose_1&&mvadbnew_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewmedium"]    = "(mvadbnew_medium_1&&mvadbnew_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewtight"]     = "(mvadbnew_tight_1&&mvadbnew_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewvtight"]     = "(mvadbnew_vtight_1&&mvadbnew_vtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadbnewvvtight"]    = "(mvadbnew_vvtight_1&&mvadbnew_vvtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldvloose"]    = "(mvadbold_vloose_1&&mvadbold_vloose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldloose"]     = "(mvadbold_loose_1&&mvadbold_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldmedium"]    = "(mvadbold_medium_1&&mvadbold_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldtight"]     = "(mvadbold_tight_1&&mvadbold_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldvtight"]     = "(mvadbold_vtight_1&&mvadbold_vtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvadboldvvtight"]    = "(mvadbold_vvtight_1&&mvadbold_vvtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldvloose"]    = "(mvapwold_vloose_1&&mvapwold_vloose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldloose"]     = "(mvapwold_loose_1&&mvapwold_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldmedium"]    = "(mvapwold_medium_1&&mvapwold_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldtight"]     = "(mvapwold_tight_1&&mvapwold_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldvtight"]     = "(mvapwold_vtight_1&&mvapwold_vtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwoldvvtight"]    = "(mvapwold_vvtight_1&&mvapwold_vvtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewvloose"]    = "(mvapwnew_vloose_1&&mvapwnew_vloose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewloose"]     = "(mvapwnew_loose_1&&mvapwnew_loose_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewmedium"]    = "(mvapwnew_medium_1&&mvapwnew_medium_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewtight"]     = "(mvapwnew_tight_1&&mvapwnew_tight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewvtight"]     = "(mvapwnew_vtight_1&&mvapwnew_vtight_2&&"+alias_map_["incnotauiso"]+")";
      alias_map_["mvapwnewvvtight"]    = "(mvapwnew_vvtight_2&&mvapwnew_vvtight_2&&"+alias_map_["incnotauiso"]+")";


      // SM Categories
      alias_map_["inclusive"]         = "1";
     // alias_map_["baseline"]          = "1";
      
      alias_map_["baseline"]          = "mva_olddm_vtight_1>0.5 && mva_olddm_vtight_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto";
      if(year_.find("6")!=year_.npos) alias_map_["baseline"]          = "mva_olddm_tight_1>0.5 && mva_olddm_tight_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto";
      alias_map_["tt_qcd_norm"]       = "mva_olddm_medium_1>0.5 && mva_olddm_loose_2>0.5 &&mva_olddm_vtight_2<0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto";
      if(year_.find("6")!=year_.npos) alias_map_["tt_qcd_norm"] = "mva_olddm_tight_1>0.5 && mva_olddm_medium_2>0.5 &&mva_olddm_tight_2<0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto";

      alias_map_["inclusivenolv"]         = "iso_1<1.0 && iso_2<1.0 && antiele_1 && antimu_1 && antiele_2 && antimu_2";
      //alias_map_["qcd_loose_shape"]   = "iso_1>1.0 && iso_2>1.0 && antiele_1 && antimu_1 && antiele_2 && antimu_2";
      alias_map_["qcd_loose_shape"]   = "iso_1>1.0 && iso_2>1.0 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto";

      alias_map_["btagnotwoprong"] = "(n_jets<=1 && n_bjets>=1&&"+alias_map_["notwoprong"]+")";

      alias_map_["nobtagnotwoprong"] = "(n_bjets==0 &&"+alias_map_["notwoprong"]+")";
      //for making CSV control plot
      alias_map_["prebtag"] = "(n_jets<=1 && n_prebjets>=1)";
      //MSSM update analysis style categories:
      alias_map_["btaglow"] = "(pt_2<60 &&"+alias_map_["btag"]+")";
      alias_map_["btaghigh"] = "(pt_2>60 &&"+alias_map_["btag"]+")";
      alias_map_["nobtaglow"] = "(pt_2<60 &&"+alias_map_["nobtag"]+")";
      alias_map_["nobtagmed"] = "(pt_2>60 && pt_2<80 &&"+alias_map_["nobtag"]+")";
      alias_map_["nobtaghigh"] = "(pt_2>80 &&"+alias_map_["nobtag"]+")";
    } else if (ch_ == channel::em) {
      // SM Categories
      alias_map_["inclusive"]         = "1";
      alias_map_["baseline"]          = "iso_1<0.15 && iso_2<0.15 && !leptonveto";
      if(year_.find("6")!=year_.npos) alias_map_["baseline"]          = "iso_1<0.15 && iso_2<0.2 && !leptonveto";
      alias_map_["inclusivenolv"]         = "iso_1<0.15 && iso_2<0.15";
      //Categories can be added using inclusive alias as follows:
      alias_map_["vbf"] = "(n_jets>=2 && n_jetsingap==0 && mjj>500 && jdeta>3.5)";
      alias_map_["1jet"] = "(!("+alias_map_["vbf"]+")"+"&& n_jets>=1 && n_bjets==0)";
      alias_map_["incnoiso"]         = "!leptonveto";
      alias_map_["incnoisowmu"]         = "iso_2<0.15 && !leptonveto";
      alias_map_["incnoisowe"]         = "iso_1<0.15 && !leptonveto";
      alias_map_["incnoisonolvwmu"]         = "iso_2<0.15";
      alias_map_["incnoisonolvwe"]         = "iso_1<0.15";
      alias_map_["edb03iso0p15"]            = "(iso_1_db03<0.15&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p16"]            = "(iso_1_db03<0.16&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p17"]            = "(iso_1_db03<0.17&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p18"]            = "(iso_1_db03<0.18&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p19"]            = "(iso_1_db03<0.18&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p20"]            = "(iso_1_db03<0.19&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p21"]            = "(iso_1_db03<0.21&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p1"]           = "(iso_1_db03<0.1&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p11"]           = "(iso_1_db03<0.11&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p12"]           = "(iso_1_db03<0.12&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p13"]           = "(iso_1_db03<0.13&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03iso0p14"]           = "(iso_1_db03<0.14&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p15"]            = "(iso_2_db03<0.15&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p16"]            = "(iso_2_db03<0.16&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p17"]            = "(iso_2_db03<0.17&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p18"]            = "(iso_2_db03<0.18&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p19"]            = "(iso_2_db03<0.19&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p20"]            = "(iso_2_db03<0.20&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p21"]            = "(iso_2_db03<0.21&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p1"]           = "(iso_2_db03<0.1&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p11"]           = "(iso_2_db03<0.11&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p12"]           = "(iso_2_db03<0.12&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p13"]           = "(iso_2_db03<0.13&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03iso0p14"]           = "(iso_2_db03<0.14&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p15"]            = "(iso_2_db04<0.15&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p16"]            = "(iso_2_db04<0.16&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p17"]            = "(iso_2_db04<0.17&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p18"]            = "(iso_2_db04<0.18&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p19"]            = "(iso_2_db04<0.19&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p20"]            = "(iso_2_db04<0.20&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p21"]            = "(iso_2_db04<0.21&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p1"]           = "(iso_2_db04<0.1&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p11"]           = "(iso_2_db04<0.11&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p12"]           = "(iso_2_db04<0.12&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p13"]           = "(iso_2_db04<0.13&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04iso0p14"]           = "(iso_2_db04<0.14&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p15"]            = "(iso_2_db04<0.15&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p16"]            = "(iso_2_db04<0.16&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p17"]            = "(iso_2_db04<0.17&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p18"]            = "(iso_2_db04<0.18&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p19"]            = "(iso_2_db04<0.19&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p20"]            = "(iso_2_db04<0.20&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p21"]            = "(iso_2_db04<0.21&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p1"]           = "(iso_2_db04<0.1&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p11"]           = "(iso_2_db04<0.11&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p12"]           = "(iso_2_db04<0.12&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p13"]           = "(iso_2_db04<0.13&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04iso0p14"]           = "(iso_2_db04<0.14&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p14"]           = "(iso_2_trk03<0.14&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p15"]            = "(iso_2_trk03<0.15&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p16"]            = "(iso_2_trk03<0.16&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p17"]            = "(iso_2_trk03<0.17&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p18"]            = "(iso_2_trk03<0.18&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p19"]            = "(iso_2_trk03<0.19&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p20"]            = "(iso_2_trk03<0.20&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p21"]            = "(iso_2_trk03<0.21&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p1"]           = "(iso_2_trk03<0.1&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p11"]           = "(iso_2_trk03<0.11&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p12"]           = "(iso_2_trk03<0.12&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mtrk03iso0p13"]           = "(iso_2_trk03<0.13&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p15"]            = "(iso_1_ea03<0.15&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p16"]            = "(iso_1_ea03<0.16&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p17"]            = "(iso_1_ea03<0.17&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p18"]            = "(iso_1_ea03<0.18&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p19"]            = "(iso_1_ea03<0.19&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p20"]            = "(iso_1_ea03<0.20&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p21"]            = "(iso_1_ea03<0.21&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p1"]           = "(iso_1_ea03<0.1&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p11"]           = "(iso_1_ea03<0.11&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p12"]           = "(iso_1_ea03<0.12&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p13"]           = "(iso_1_ea03<0.13&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["eea03iso0p14"]           = "(iso_1_ea03<0.14&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p15"]            = "(iso_2_ea03<0.15&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p16"]            = "(iso_2_ea03<0.16&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p17"]            = "(iso_2_ea03<0.17&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p18"]            = "(iso_2_ea03<0.18&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p19"]            = "(iso_2_ea03<0.19&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p20"]            = "(iso_2_ea03<0.20&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p21"]            = "(iso_2_ea03<0.21&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p1"]           = "(iso_2_ea03<0.1&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p11"]           = "(iso_2_ea03<0.11&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p12"]           = "(iso_2_ea03<0.12&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p13"]           = "(iso_2_ea03<0.13&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mea03iso0p14"]           = "(iso_2_ea03<0.14&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p15"]            = "(iso_1_db03allch<0.15&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p16"]            = "(iso_1_db03allch<0.16&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p17"]            = "(iso_1_db03allch<0.17&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p18"]            = "(iso_1_db03allch<0.18&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p19"]            = "(iso_1_db03allch<0.19&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p20"]            = "(iso_1_db03allch<0.20&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p21"]            = "(iso_1_db03allch<0.21&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p1"]           = "(iso_1_db03allch<0.1&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p11"]           = "(iso_1_db03allch<0.11&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p12"]           = "(iso_1_db03allch<0.12&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p13"]           = "(iso_1_db03allch<0.13&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb03allchiso0p14"]           = "(iso_1_db03allch<0.14&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p15"]            = "(iso_2_db03allch<0.15&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p16"]            = "(iso_2_db03allch<0.16&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p17"]            = "(iso_2_db03allch<0.17&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p18"]            = "(iso_2_db03allch<0.18&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p19"]            = "(iso_2_db03allch<0.19&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p20"]            = "(iso_2_db03allch<0.20&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p21"]            = "(iso_2_db03allch<0.21&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p1"]           = "(iso_2_db03allch<0.1&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p11"]           = "(iso_2_db03allch<0.11&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p12"]           = "(iso_2_db03allch<0.12&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p13"]           = "(iso_2_db03allch<0.13&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb03allchiso0p14"]           = "(iso_2_db03allch<0.14&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p15"]            = "(iso_1_db04allch<0.15&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p16"]            = "(iso_1_db04allch<0.16&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p17"]            = "(iso_1_db04allch<0.17&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p18"]            = "(iso_1_db04allch<0.18&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p19"]            = "(iso_1_db04allch<0.19&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p20"]            = "(iso_1_db04allch<0.20&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p21"]            = "(iso_1_db04allch<0.21&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p1"]           = "(iso_1_db04allch<0.1&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p11"]           = "(iso_1_db04allch<0.11&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p12"]           = "(iso_1_db04allch<0.12&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p13"]           = "(iso_1_db04allch<0.13&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["edb04allchiso0p14"]           = "(iso_1_db04allch<0.14&&iso_2<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p15"]            = "(iso_2_db04allch<0.15&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p16"]            = "(iso_2_db04allch<0.16&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p17"]            = "(iso_2_db04allch<0.17&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p18"]            = "(iso_2_db04allch<0.18&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p19"]            = "(iso_2_db04allch<0.19&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p20"]            = "(iso_2_db04allch<0.20&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p21"]            = "(iso_2_db04allch<0.21&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p1"]           = "(iso_2_db04allch<0.1&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p11"]           = "(iso_2_db04allch<0.11&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p12"]           = "(iso_2_db04allch<0.12&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p13"]           = "(iso_2_db04allch<0.13&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";
      alias_map_["mdb04allchiso0p14"]           = "(iso_2_db04allch<0.14&&iso_1<0.15&&"+alias_map_["incnoiso"]+")";

      //alias_map_["qcd_loose_shape"]         = "(iso_1>0.2&&iso_1<0.5  && iso_2>0.2&&iso_2<0.5)";
      alias_map_["qcd_loose_shape"]         = "(iso_1>0.2&&iso_1<0.5  && iso_2>0.2&&iso_2<0.5 && !leptonveto)";
      alias_map_["btag"] = "(n_jets<=1 && n_bjets>=1)";
      alias_map_["nobtag"] = "n_bjets==0";
      alias_map_["ttcontrolalt"] = "(n_jets>=1 && n_bjets>=1 && pzeta<-50)";
      alias_map_["ttcontrol"] = "(pzeta<-20 && met>80)";
      //for making CSV control plot
      alias_map_["prebtag"] = "(n_jets<=1 && n_prebjets>=1)";
    } else if (ch_ == channel::zmm || ch_ == channel::zee) {
      alias_map_["inclusive"]         = "1";
      alias_map_["btag"] = "(n_jets<=1 && n_bjets>=1)";
      alias_map_["nobtag"] = "n_bjets==0";
      alias_map_["baseline"]         = "(iso_1<0.1 && iso_2<0.1)";
      if(ch_ == channel::zmm) alias_map_["baseline"] = "(iso_1<0.15&&iso_2<0.15)";
    } else if (ch_ == channel::tpzmm || ch_ == channel::tpzee) {
      alias_map_["inclusive"]         = "1";
      alias_map_["baseline"]         = "(iso_1<0.1 && iso_2<0.1)";
    } else if (ch_ == channel::wmnu) {
      alias_map_["inclusive"]         = "1";
      alias_map_["baseline"]         = "(iso_1<0.1)";
    }
    
    // Selection control regions
    alias_map_["sel"]                    = "1";
    alias_map_["w_sdb"]                  = "mt_1>70.";
    alias_map_["w_sdb_os"]               = "os";
    alias_map_["w_vbf_sdb"]              = "mt_1>60. && mt_1<120.";
    alias_map_["w_os"]                   = "os";
    alias_map_["w_vbf_os"]               = "os";
    alias_map_["w_ss"]                   = "!os";
    alias_map_["w_shape_os"]             = "os"; 

    if (ch_ == channel::em) {
    // Inclusive region for e-mu fake rate normalisation
      // Sideband region for e-mu SS fake shape estimation
      alias_map_["ss"]                        = "!os";
      alias_map_["sel"]                       = "1";
    }

//    alias_map_["ZTT_Shape_Sample"]  = "DYJetsToLL_M-50-LO";

    // Samples to combine for diboson contribution
    samples_alias_map_["vv_samples"] = {
     "T-tW", "Tbar-tW", "T-t","Tbar-t",
     "WWTo1L1Nu2Q",
     "VVTo2L2Nu","ZZTo2L2Q","ZZTo4L",
     "WZTo2L2Q","WZJetsToLLLNu","WZTo1L3Nu","WZTo1L1Nu2Q"
    };

    if(year_.find("6")!=year_.npos){
      samples_alias_map_["vv_samples"] = {
       "T-tW", "Tbar-tW","Tbar-t","T-t",
       "WWTo1L1Nu2Q","WZJToLLLNu",
       "VVTo2L2Nu","ZZTo2L2Q","ZZTo4L",
       "WZTo2L2Q","WZTo1L3Nu","WZTo1L1Nu2Q"
      };
    } 
/*  if(ch_==channel::em){
    samples_alias_map_["vv_samples"] = {
//     "WZJetsTo3LNu",
     "T-tW", "Tbar-tW", "WWinclusive","WZinclusive", "ZZinclusive"//,"WWTo2L2Nu","WWTo4Q","WZTo1L1Nu2Q","ZZTo4L"
    //,"WWTo2L2Nu","WWTo4Q","WZTo1L1Nu2Q","ZZTo4L"
    };
  }*/

    samples_alias_map_["wgam_samples"] = {
        "WGToLNuG","WGstarToLNuEE","WGstarToLNuMuMu"
     };


    samples_alias_map_["top_samples"] = {
     "TT-ext"
    };
    if(year_.find("6")!=year_.npos){
      samples_alias_map_["top_samples"] = {
       "TT"
      };
    }

 
   samples_alias_map_["ztt_shape_samples"]={
    "DYJetsToLL_M-50-LO-ext",
    "DY1JetsToLL_M-50-LO","DY2JetsToLL_M-50-LO",
    "DY3JetsToLL_M-50-LO","DY4JetsToLL_M-50-LO"/*,
    /"DYJetsToLL_M-150-LO"*/,"DYJetsToLL_M-10-ext"
   };
/*   if(year_.find("6")!=year_.npos){
   samples_alias_map_["ztt_shape_samples"]={
      "DYJetsToLL_M-LO"
     };*/
   if(year_.find("6")!=year_.npos){
   samples_alias_map_["ztt_shape_samples"]={
      "DYJetsToLL-LO-ext","DY1JetsToLL-LO","DY2JetsToLL-LO","DY3JetsToLL-LO",
      "DY4JetsToLL-LO","DYJetsToLL_M-10-50-LO"
     };
  }

   samples_alias_map_["data_samples"] = {
   "SingleElectron-2015D-prompt"
   };
   if(ch_==channel::et || ch_==channel::zee|| ch_==channel::tpzee){
     samples_alias_map_["data_samples"] = {
      "SingleElectron-2015D"
     };
    if(year_.find("6")!=year_.npos){
       samples_alias_map_["data_samples"] = {
        "SingleElectronB","SingleElectronC","SingleElectronD"
       };
    }
   }
   if(ch_==channel::mt || ch_==channel::zmm || ch_==channel::tpzmm || ch_==channel::wmnu){
     samples_alias_map_["data_samples"] = {
     "SingleMuon-2015D"
     };
    if(year_.find("6")!=year_.npos){
     samples_alias_map_["data_samples"] = {
       "SingleMuonB","SingleMuonC","SingleMuonD"
       };
    }
   }
   if(ch_==channel::tt){
     samples_alias_map_["data_samples"] = {
      "Tau-2015D"
     };
    if(year_.find("6")!=year_.npos){ 
      samples_alias_map_["data_samples"] = {
        "TauB","TauC","TauD"
      };
    } 
   }
   if(ch_==channel::em){
    samples_alias_map_["data_samples"] = {
      "MuonEG-2015D"
    };
    if(year_.find("6")!=year_.npos){
      samples_alias_map_["data_samples"] = {
        "MuonEGB","MuonEGC","MuonEGD"
      };
    }
   }




   samples_alias_map_["ztt_samples"]={
     "DYJetsToLL_M-50-LO-ext",
     "DY1JetsToLL_M-50-LO","DY2JetsToLL_M-50-LO",
     "DY3JetsToLL_M-50-LO","DY4JetsToLL_M-50-LO"/*,
     "DYJetsToLL_M-150-LO"*/,"DYJetsToLL_M-10-ext"
   };
 /*  if(year_.find("6")!=year_.npos){
     samples_alias_map_["ztt_samples"]={
       "DYJetsToLL-LO"
     };*/
   if(year_.find("6")!=year_.npos){
   samples_alias_map_["ztt_samples"]={
      "DYJetsToLL-LO-ext","DY1JetsToLL-LO","DY2JetsToLL-LO","DY3JetsToLL-LO",
      "DY4JetsToLL-LO","DYJetsToLL_M-10-50-LO"
     };
   }

  
 if(ch_!=channel::em){
  samples_alias_map_["qcd_sub_samples"] = {
   "DYJetsToLL_M-50-LO-ext",/*"DYJetsToLL_M-150-LO",*/"DYJetsToLL_M-10-ext",
   "T-tW", "Tbar-tW", "T-t","Tbar-t",
   "WWTo1L1Nu2Q",
   "VVTo2L2Nu","ZZTo2L2Q","ZZTo4L",
   "WZTo2L2Q","WZJetsToLLLNu","WZTo1L3Nu","WZTo1L1Nu2Q",
   "WJetsToLNu-LO","TT-ext",
   "DY1JetsToLL_M-50-LO","DY2JetsToLL_M-50-LO",
   "DY3JetsToLL_M-50-LO","DY4JetsToLL_M-50-LO",
   "W1JetsToLNu-LO","W2JetsToLNu-LO",
   "W3JetsToLNu-LO","W4JetsToLNu-LO"
   };

  /*if(year_.find("6")!=year_.npos){
     samples_alias_map_["qcd_sub_samples"] = {
        "DYJetsToLL-LO",
        "T-tW","T-t", "Tbar-tW","Tbar-t",
        "WWTo1L1Nu2Q","VVTo2L2Nu",
        "ZZTo2L2Q","ZZTo4L",
        "WZTo2L2Q","WZTo1L3Nu","WZTo1L1Nu2Q",
        "TT","WJetsToLNu-LO",
        };*/
  if(year_.find("6")!=year_.npos){
     samples_alias_map_["qcd_sub_samples"] = {
      "DYJetsToLL-LO-ext","DY1JetsToLL-LO","DY2JetsToLL-LO","DY3JetsToLL-LO",
      "DY4JetsToLL-LO","DYJetsToLL_M-10-50-LO","ZZTo4L",
        "T-tW","T-t","Tbar-tW","Tbar-t",
        "WWTo1L1Nu2Q","VVTo2L2Nu",
        "ZZTo2L2Q","WZJToLLLNu",
        "WZTo2L2Q","WZTo1L3Nu","WZTo1L1Nu2Q",
        "TT",
        "WJetsToLNu-LO","W1JetsToLNu-LO","W2JetsToLNu-LO",
        "W3JetsToLNu-LO","W4JetsToLNu-LO"
        };
  }
  }
  

 if(ch_==channel::em){
  samples_alias_map_["qcd_sub_samples"] = {
   "DYJetsToLL_M-50-LO-ext",/*"DYJetsToLL_M-150-LO",*/"DYJetsToLL_M-10-ext",
   "T-tW", "Tbar-tW", "T-t","Tbar-t",
   "WWTo1L1Nu2Q","VVTo2L2Nu",
   "ZZTo2L2Q","ZZTo4L",
   "WZTo2L2Q","WZJetsToLLLNu","WZTo1L3Nu","WZTo1L1Nu2Q",
   "TT-ext","WJetsToLNu-LO",
   "DY1JetsToLL_M-50-LO","DY2JetsToLL_M-50-LO",
   "DY3JetsToLL_M-50-LO","DY4JetsToLL_M-50-LO",
   "W1JetsToLNu-LO","W2JetsToLNu-LO",
   "W3JetsToLNu-LO","W4JetsToLNu-LO",
   "WGToLNuG","WGstarToLNuEE","WGstarToLNuMuMu"
   };
  if(year_.find("6")!=year_.npos){
  samples_alias_map_["qcd_sub_samples"] = {
   "DYJetsToLL-LO-ext",
   "DY1JetsToLL-LO","DY2JetsToLL-LO",
   "DY3JetsToLL-LO","DY4JetsToLL-LO",
   "DYJetsToLL_M-10-50-LO",
   "T-tW", "Tbar-tW", "Tbar-t", "T-t",
   "WWTo1L1Nu2Q","VVTo2L2Nu", "ZZTo4L",
   "ZZTo2L2Q","WZJToLLLNu",
   "WZTo2L2Q","WZTo1L3Nu","WZTo1L1Nu2Q",
   "TT","WJetsToLNu-LO",
   "W1JetsToLNu-LO","W2JetsToLNu-LO",
   "W3JetsToLNu-LO","W4JetsToLNu-LO"
   };
  }
  }

  samples_alias_map_["w_sub_samples"] = {
   "DYJetsToLL_M-50-LO-ext",/*"DYJetsToLL_M-150-LO",*/"DYJetsToLL_M-10-ext",
   "DY1JetsToLL_M-50-LO","DY2JetsToLL_M-50-LO",
   "DY3JetsToLL_M-50-LO","DY4JetsToLL_M-50-LO",
   "T-tW", "Tbar-tW", "T-t","Tbar-t",
   "WWTo1L1Nu2Q","VVTo2L2Nu",
   "ZZTo2L2Q","ZZTo4L",
   "WZTo2L2Q","WZJetsToLLLNu","WZTo1L3Nu","WZTo1L1Nu2Q",
   "TT-ext"
   };
   /*if(year_.find("6")!=year_.npos){
    samples_alias_map_["w_sub_samples"] = {
     "DYJetsToLL-LO",
     "T-tW", "Tbar-tW","T-t", "Tbar-t",
     "WWTo1L1Nu2Q","VVTo2L2Nu",
     "ZZTo2L2Q","ZZTo4L",
     "WZTo2L2Q","WZTo1L3Nu","WZTo1L1Nu2Q",
     "TT"
     };*/
   if(year_.find("6")!=year_.npos){
    samples_alias_map_["w_sub_samples"] = {
      "DYJetsToLL-LO-ext","DY1JetsToLL-LO","DY2JetsToLL-LO","DY3JetsToLL-LO",
      "DY4JetsToLL-LO","DYJetsToLL_M-10-50-LO",
     "T-tW", "Tbar-tW", "Tbar-t","T-t",
     "WWTo1L1Nu2Q","VVTo2L2Nu","ZZTo4L",
     "ZZTo2L2Q","WZJToLLLNu",
     "WZTo2L2Q","WZTo1L3Nu","WZTo1L1Nu2Q",
     "TT"
     };
   } 

samples_alias_map_["wjets_samples"] = {
  "WJetsToLNu-LO",
  "W1JetsToLNu-LO","W2JetsToLNu-LO",
  "W3JetsToLNu-LO","W4JetsToLNu-LO"
 };
if(ch_==channel::em) {
    samples_alias_map_["wjets_samples"] = {
      "WJetsToLNu-LO",
      "W1JetsToLNu-LO","W2JetsToLNu-LO",
      "W3JetsToLNu-LO","W4JetsToLNu-LO"
     };
}
/*if(year_.find("6")!=year_.npos){
  samples_alias_map_["wjets_samples"] = {
    "WJetsToLNu-LO"
  };*/
if(year_.find("6")!=year_.npos){
  samples_alias_map_["wjets_samples"] = {
    "WJetsToLNu-LO","W1JetsToLNu-LO","W2JetsToLNu-LO",
    "W3JetsToLNu-LO","W4JetsToLNu-LO"
  };
}

sample_names_={};
push_back(sample_names_,this->ResolveSamplesAlias("ztt_samples"));
push_back(sample_names_,this->ResolveSamplesAlias("vv_samples"));
if(ch_==channel::em){
  push_back(sample_names_,this->ResolveSamplesAlias("wgam_samples"));
}
push_back(sample_names_,this->ResolveSamplesAlias("wjets_samples"));
push_back(sample_names_,this->ResolveSamplesAlias("top_samples"));
push_back(sample_names_,this->ResolveSamplesAlias("data_samples"));

 }

  void HTTRun2Analysis::SetQCDRatio(double const& ratio){
    qcd_os_ss_factor_ = ratio;
  }

  void HTTRun2Analysis::AddSample(std::string const& sample) {
    sample_names_.push_back(sample);
  }

  void HTTRun2Analysis::ParseParamFile(std::string const& file) {
    SimpleParamParser parser;
    parser.ParseFile(file);
    std::cout << "[HTTRun2Analysis::ParseParamFile] Extracting sample info from file " << file << std::endl;
    lumi_ = parser.GetParam<double>("LUMI_"+year_+"_"+Channel2String(ch_));
    std::cout << "[HTTRun2Analysis::ParseParamFile] Integrated luminosity set to " << lumi_ << " /pb" << std::endl;
//    if (verbosity_ > 1) std::cout << boost::format("%-25s %152 %15.3f %15.3f %15.3f\n") % "Sample" % "Events" % "Cross Section" % "Sample Lumi" % "Rel. Lumi";
    if (verbosity_ > 2) std::cout << "-----------------------------------------------------------------------------------------\n";
    for (auto sample : sample_names_) {
      std::string lookup = sample;
      if (sample.find("Special") != sample.npos) {
        size_t found = sample.find('_');
        found = sample.find('_',found+1);
        lookup = sample.substr(found+1);
      }
      if (parser.ParamExists("XS_"+sample) && parser.ParamExists("EVT_"+sample)) {
        double evt = parser.GetParam<double>("EVT_"+sample);
        double xs = parser.GetParam<double>("XS_"+sample);
        if (xs <= 0) continue;
        sample_info_[sample] = std::make_pair(evt, xs);
      //  if (verbosity_ > 1) std::cout << boost::format("%-25s %15i %15.3f %15.3f %15.3f\n") % sample % unsigned(evt+0.5) % xs % (evt/xs) % (evt/(xs*lumi_));
      }
    }
  }

  void HTTRun2Analysis::AddSMSignalSamples(std::vector<std::string> masses) {
    for (auto m : masses) {
        sample_names_.push_back("GluGluHToTauTau_M-"+m);
        sample_names_.push_back("VBFHToTauTau_M-"+m);
        sample_names_.push_back("WplusHToTauTau_M-"+m);
        sample_names_.push_back("WminusHToTauTau_M-"+m);
        sample_names_.push_back("ZHToTauTau_M-"+m);
        sample_names_.push_back("TTHToTauTau_M-"+m);
    }
  }

  void HTTRun2Analysis::AddMSSMSignalSamples(std::vector<std::string> masses) {
    for (auto m : masses) {
        sample_names_.push_back("SUSYGluGluToHToTauTau_M-"+m);
        sample_names_.push_back("SUSYGluGluToBBHToTauTau_M-"+m);
    }
  }
 
  void HTTRun2Analysis::AddMSSMSignalSamplesBBH(std::vector<std::string> masses){
  for (auto m : masses) {
        sample_names_.push_back("SUSYGluGluToBBHToTauTau_M-"+m);
    }
  }

  void HTTRun2Analysis::AddMSSMSignalSamplesGGH(std::vector<std::string> masses){
  for (auto m : masses) {
        sample_names_.push_back("SUSYGluGluToHToTauTau_M-"+m);
    }
  }


  void HTTRun2Analysis::AddHhhSignalSamples(std::vector<std::string> masses) {
    for (auto m : masses) {
        sample_names_.push_back("GluGluToRadionToHHTo2B2Tau_M-"+m);
    }
  }


  void HTTRun2Analysis::ReadTrees(std::string const& folder, std::string const& fallback_folder) {
    std::cout << "[HTTRun2Analysis::ReadTrees] Reading input files..." << std::endl;
    std::vector<std::string> result_summary;
    for (auto name : sample_names_) {
      // The input file is folder + sample name + channel + year
      std::string input_filename = folder+"/"+name+"_"+Channel2String(ch_)+"_"+year_+".root";
      std::string label = name;
      TFile *tmp_file = nullptr;
      if (boost::filesystem::exists(input_filename)) tmp_file = TFile::Open(input_filename.c_str());
      if (!tmp_file && fallback_folder != "") {
        if (verbosity_ > 2) std::cout << "[HTTRun2Analysis::ReadTrees] " << input_filename << " not found, trying fallback folder" << std::endl;
        input_filename = fallback_folder+"/"+name+"_"+Channel2String(ch_)+"_"+year_+".root";
        if (boost::filesystem::exists(input_filename)) tmp_file = TFile::Open(input_filename.c_str());
      }
      if (!tmp_file) {
        std::cout << "[HTTRun2Analysis::ReadTrees] Warning: " << input_filename << " cannot be opened" << std::endl;
        continue;
      }
      if (verbosity_ > 2) result_summary.push_back((boost::format("%-70s %s %-30s\n") % input_filename % "-->" % label).str());
      gDirectory->cd("/");
      TTree *tmp_tree = dynamic_cast<TTree*>(gDirectory->Get("ntuple"));
      if (!tmp_tree) {
        std::cerr << "[HTTRun2Analysis::ReadTrees] Warning: Unable to extract TTree from file " << input_filename << std::endl;
        continue;        
      }
      tmp_tree->SetEstimate(100000);
      tfiles_[label] = tmp_file;
      ttrees_[label] = tmp_tree;
    }
    for (auto str : result_summary) std::cout << str;
  }

  double HTTRun2Analysis::GetLumiScale(std::string const& sample) {
    auto it = sample_info_.find(sample);
    if (it != sample_info_.end()) {
      double evt = it->second.first;
      double xs = it->second.second;
      return ((xs*lumi_)/evt);
    } else if(sample.find("SingleMuon")==sample.npos && sample.find("SingleEle")==sample.npos && sample.find("MuonEG")==sample.npos && sample.find("Tau")==sample.npos){
      std::cout << "[HTTRun2Analysis::GetLumiScale] Warning: lumi scale not found for sample " << sample << std::endl;
      return 1.0;
    } else {
      return 1.0;
    }
  }

  double HTTRun2Analysis::GetLumiScaleFixedXS(std::string const& sample, double xs) {
      auto it = sample_info_.find(sample);
    if (it != sample_info_.end()) {
      double evt = it->second.first;;
      return ((xs*lumi_)/evt);
    } else {
      return 1.0;
    }
  }


  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateData(unsigned /*method*/, std::string var, std::string sel, std::string cat, std::string wt) {
    cat += "&&" + alias_map_["baseline"];
    auto data_norm = this->GetRate(this->ResolveSamplesAlias("data_samples"), sel, cat, wt);
    TH1F data_hist = this->GetShape(var, this->ResolveSamplesAlias("data_samples"), sel, cat, wt);
    SetNorm(&data_hist, data_norm.first);
    return std::make_pair(data_hist, data_norm);
  }



  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateZTT(unsigned /*method*/, std::string var, std::string sel, std::string cat, std::string wt) {
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateZTT] --------------------------------------------------------\n";
    cat += "&&" + alias_map_["baseline"];
    Value ztt_norm;
    std::string inclusive = this->ResolveAlias("inclusive") + "&&" + alias_map_["baseline"]; 
    if(verbosity_) this->GetRateViaRefEfficiency(this->ResolveSamplesAlias("ztt_samples"), this->ResolveSamplesAlias("ztt_samples"), sel, inclusive, sel, cat, wt);
    std::vector<std::string> ztt_samples = this->ResolveSamplesAlias("ztt_samples");
    if (verbosity_) {
      std::cout << "ztt_samples: ";
      for (unsigned i = 0; i < ztt_samples.size(); ++i) {
        std::cout << ztt_samples[i];
        if (i != ztt_samples.size()-1) std::cout << ", ";
      }
      std::cout << std::endl;
    }
//    std::vector<std::string> ztt_shape_samples = this->ResolveSamplesAlias("ztt_samples");
    ztt_norm = this->GetLumiScaledRate(ztt_samples, sel, cat, wt) ;
    TH1F ztt_hist = this->GetLumiScaledShape(var, ztt_samples, sel, cat, wt);
//    if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
 //     % this->ResolveAlias("ZTT_Shape_Sample") % sel % cat % wt;
    SetNorm(&ztt_hist, ztt_norm.first);
    return std::make_pair(ztt_hist, ztt_norm);
  }
  
  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateZL(unsigned /*method*/, std::string var, std::string sel, std::string cat, std::string wt) {
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateZL --------------------------------------------------------\n";
    cat += "&&" + alias_map_["baseline"];
    Value zl_norm;
    //ztt_norm = this->GetRateViaRefEfficiency(this->ResolveAlias("ZTT_Eff_Sample"), "DYJetsToLL", "os", this->ResolveAlias("inclusive"), sel, cat, wt);
    std::vector<std::string> zl_samples = this->ResolveSamplesAlias("ztt_samples");
    zl_norm = this->GetLumiScaledRate(zl_samples, sel, cat, wt) ;
    if (verbosity_) {
      std::cout << "zl_samples: ";
      for (unsigned i = 0; i < zl_samples.size(); ++i) {
        std::cout << zl_samples[i];
        if (i != zl_samples.size()-1) std::cout << ", ";
      }
      std::cout << std::endl;
    }
    TH1F zl_hist = this->GetLumiScaledShape(var, zl_samples, sel, cat, wt);
//    if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
 //     % "DYJetsToLL" % sel % cat % wt;
    SetNorm(&zl_hist, zl_norm.first);
    return std::make_pair(zl_hist, zl_norm);
  }
  
  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateZJ(unsigned /*method*/, std::string var, std::string sel, std::string cat, std::string wt) {
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateZJ --------------------------------------------------------\n";
    cat += "&&" + alias_map_["baseline"];
    Value zj_norm;
    //ztt_norm = this->GetRateViaRefEfficiency(this->ResolveAlias("ZTT_Eff_Sample"), "DYJetsToLL", "os", this->ResolveAlias("inclusive"), sel, cat, wt);
    std::vector<std::string> zj_samples = this->ResolveSamplesAlias("ztt_samples");
    if (verbosity_) {
      std::cout << "zj_samples: ";
      for (unsigned i = 0; i < zj_samples.size(); ++i) {
        std::cout << zj_samples[i];
        if (i != zj_samples.size()-1) std::cout << ", ";
      }
      std::cout << std::endl;
    }
    zj_norm = this->GetLumiScaledRate(zj_samples, sel, cat, wt) ;
    TH1F zj_hist = this->GetLumiScaledShape(var, zj_samples, sel, cat, wt);
//    if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
 //     % "zj_samples" % sel % cat % wt;
    SetNorm(&zj_hist, zj_norm.first);
    return std::make_pair(zj_hist, zj_norm);
  }
  
  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateZLL(unsigned /*method*/, std::string var, std::string sel, std::string cat, std::string wt) {
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateZLL --------------------------------------------------------\n";
    cat += "&&" + alias_map_["baseline"];
    std::string inclusive = this->ResolveAlias("inclusive") + "&&" + alias_map_["baseline"]; 
    Value zl_norm;
    if(verbosity_) this->GetRateViaRefEfficiency(this->ResolveSamplesAlias("ztt_samples"), this->ResolveSamplesAlias("ztt_samples"), sel, inclusive, sel, cat, wt);
    std::vector<std::string> zll_samples = this->ResolveSamplesAlias("ztt_samples");
    if (verbosity_) {
      std::cout << "zll_samples: ";
      for (unsigned i = 0; i < zll_samples.size(); ++i) {
        std::cout << zll_samples[i];
        if (i != zll_samples.size()-1) std::cout << ", ";
      }
      std::cout << std::endl;
    }
    zl_norm = this->GetLumiScaledRate(zll_samples, sel, cat, wt) ;
    TH1F zl_hist = this->GetLumiScaledShape(var, zll_samples, sel, cat, wt);
//    if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
 //     % "DYJetsToLL" % sel % cat % wt;
    SetNorm(&zl_hist, zl_norm.first);
    return std::make_pair(zl_hist, zl_norm);
  }


  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateTOP(unsigned /*method*/, std::string var, std::string sel, std::string cat, std::string wt) {
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateTOP] --------------------------------------------------------\n";
    cat += "&&" + alias_map_["baseline"];
    std::vector<std::string> top_samples = this->ResolveSamplesAlias("top_samples");
    if (verbosity_) {
      std::cout << "top_samples: ";
      for (unsigned i = 0; i < top_samples.size(); ++i) {
        std::cout << top_samples[i];
        if (i != top_samples.size()-1) std::cout << ", ";
      }
      std::cout << std::endl;
    }
    auto top_norm = this->GetLumiScaledRate(top_samples, sel, cat, wt);
    std::string top_shape_cat = cat;
    TH1F top_hist = this->GetLumiScaledShape(var, top_samples, sel, top_shape_cat, wt);
    // TH1F top_hist = this->GetLumiScaledShape(var, top_shape_sample, sel, top_shape_cat, wt);
    if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
      % "top_samples" % sel % top_shape_cat % wt;
    SetNorm(&top_hist, top_norm.first);
    return std::make_pair(top_hist, top_norm);
  }

  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateVV(unsigned /*method*/, std::string var, std::string sel, std::string cat, std::string wt) {
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateVV] ---------------------------------------------------------\n";
    cat += "&&" + alias_map_["baseline"];
    std::vector<std::string> vv_samples = this->ResolveSamplesAlias("vv_samples");
    auto vv_norm = this->GetLumiScaledRate(vv_samples, sel, cat, wt);
    std::string vv_shape_cat = cat;
    TH1F vv_hist = this->GetLumiScaledShape(var, vv_samples, sel, vv_shape_cat, wt);
    if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
      % "vv_samples" % sel % vv_shape_cat % wt;
    SetNorm(&vv_hist, vv_norm.first);
    return std::make_pair(vv_hist, vv_norm);
  }
  
  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateWGamma(unsigned /*method*/, std::string var, std::string sel, std::string cat, std::string wt) {
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateWGamma] ---------------------------------------------------------\n";
    cat += "&&" + alias_map_["baseline"];
    std::vector<std::string> wgam_samples = this->ResolveSamplesAlias("wgam_samples");
    auto wgam_norm = this->GetLumiScaledRate(wgam_samples, sel, cat, wt);
    std::string wgam_shape_cat = cat;
    TH1F wgam_hist = this->GetLumiScaledShape(var, wgam_samples, sel, wgam_shape_cat, wt);
    if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
      % "wgam_samples" % sel % wgam_shape_cat % wt;
    SetNorm(&wgam_hist, wgam_norm.first);
    return std::make_pair(wgam_hist, wgam_norm);
  }

  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateW(unsigned method, std::string var, std::string sel, std::string cat, std::string wt) {
    cat += "&&" + alias_map_["baseline"];
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateW] ----------------------------------------------------------\n";
    std::vector<std::string> w_sub_samples = this->ResolveSamplesAlias("w_sub_samples");
    //std::string cat_nobtag = "n_jets <=1 && "+alias_map_["baseline"];
    std::string cat_nobtag = "n_jets <=1 && n_lowpt_jets>=1 && "+alias_map_["baseline"];
    std::string w_extrap_cat = cat;
    std::string w_extrap_cat_nobtag = cat_nobtag;
    std::string w_extrp_sdb_sel = this->ResolveAlias("w_os")+" && "+this->ResolveAlias("w_sdb");
    std::string w_extrp_sig_sel = this->ResolveAlias("w_os")+" && "+this->ResolveAlias("sel");
    std::string w_sdb_sel = this->ResolveAlias("w_sdb_os")+" && "+this->ResolveAlias("w_sdb");
    std::string w_sdb_sel_osss = this->ResolveAlias("w_sdb");
    
    Value w_norm;
    std::vector<std::string> wjets_samples = this->ResolveSamplesAlias("wjets_samples");
    if(method == 8 || method == 9 || method == 15) {
      w_norm = this->GetLumiScaledRate(wjets_samples, sel, cat, wt);
    } else if(method == 10 || method == 11){
     w_norm = this->GetRateViaWMethod(wjets_samples, w_extrap_cat, w_extrp_sdb_sel, w_extrp_sig_sel, 
        this->ResolveSamplesAlias("data_samples"), cat, w_sdb_sel, w_sub_samples, wt, ValueFnMap());
    } else if(method == 12 || method == 13 || method == 14){
     w_norm = this->GetRateViaWOSSSMethod(wjets_samples, w_extrap_cat, w_extrp_sdb_sel, w_extrp_sig_sel, 
        this->ResolveSamplesAlias("data_samples"), cat, w_sdb_sel_osss, w_sub_samples, !do_ss_, wt, ValueFnMap());
    } else if(method == 16){
     Value w_norm_nobtag = this->GetRateViaWOSSSMethod(wjets_samples, w_extrap_cat_nobtag, w_extrp_sdb_sel, w_extrp_sig_sel,
        this->ResolveSamplesAlias("data_samples"),cat_nobtag, w_sdb_sel_osss,w_sub_samples,!do_ss_, wt, ValueFnMap());
        std::string btag_extrap_sel;
        btag_extrap_sel = do_ss_ ?  "!os && "+this->ResolveAlias("sel") : "os &&"+this->ResolveAlias("sel");
        Value btag_norm = GetRate(wjets_samples, btag_extrap_sel, cat, wt);
        Value nobtag_norm = GetRate(wjets_samples,btag_extrap_sel,cat_nobtag,wt);
        Value btag_extrap_factor = ValueDivide(btag_norm,nobtag_norm);
        if (verbosity_) PrintValue("ExtrapFactor into full b-tag definition",btag_extrap_factor);
        w_norm = ValueProduct(w_norm_nobtag,btag_extrap_factor);
   }

    std::string w_shape_cat = cat;
    if (method == 14) w_shape_cat = "n_jets<=1&&n_loose_bjets>=1&&"+alias_map_["baseline"];
    std::string w_shape_sel = this->ResolveAlias("w_shape_os") + " && " + this->ResolveAlias("sel");
    TH1F w_hist = this->GetShape(var, wjets_samples, w_shape_sel, w_shape_cat, wt);
    //if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
     // % wjets_samples % w_shape_sel % w_shape_cat % wt;
    if(verbosity_ > 1){
      std::string bin_delim="";
      if(var.find("[")!=std::string::npos) bin_delim = "[";
      if(var.find("(")!=std::string::npos) bin_delim = "(";
      this->KolmogorovTest(var.substr(0, var.find(bin_delim, 0)), wjets_samples.at(0), w_shape_sel, cat,wjets_samples.at(0), w_shape_sel, w_shape_cat,wt);
      TH1F default_w_hist = this->GetShape(var,wjets_samples,w_shape_sel, cat, wt);
      std::cout << "ROOT KS test: "<< default_w_hist.KolmogorovTest(&w_hist) <<std::endl;
    }

    SetNorm(&w_hist, w_norm.first);
    return std::make_pair(w_hist, w_norm);
  }

  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateQCD(unsigned method, std::string var, std::string /*sel*/, std::string cat, std::string wt) {
    if (verbosity_) std::cout << "[HTTRun2Analysis::GenerateQCD] --------------------------------------------------------\n";
    std::string maincat = cat + "&&" + alias_map_["baseline"];
    std::string ttqcdcat = cat + "&&" + alias_map_["incnotauiso"];
    Value qcd_norm;
    TH1F qcd_hist;
//    if (ch_ != channel::em) {
      std::vector<std::string> qcd_sub_samples;
      qcd_sub_samples = this->ResolveSamplesAlias("qcd_sub_samples");
      
      std::vector<std::string> w_sub_samples;
       if(ch_ != channel::em){
          w_sub_samples = this->ResolveSamplesAlias("w_sub_samples");
       }
      std::string qcd_sdb_sel = "!os && " + this->ResolveAlias("sel");
      std::string w_extrp_sdb_sel = this->ResolveAlias("w_ss")+" && "+this->ResolveAlias("w_sdb");
      std::string w_extrp_sig_sel = this->ResolveAlias("w_ss")+" && "+this->ResolveAlias("sel");
      std::string w_sdb_sel = "!os && "+this->ResolveAlias("w_sdb");
      std::string w_sdb_sel_osss = this->ResolveAlias("w_sdb");
      std::string qcd_cat = maincat;
      std::string qcd_sdb_cat = maincat;
      std::string qcd_extrap_sel = "!os && " + this->ResolveAlias("sel");
      //Default QCD method for fully hadronic is different...
      if(ch_ == channel::tt && method == 8){
        qcd_sdb_sel =  do_ss_ ? "!os &&" + this->ResolveAlias("sel") : "os && " + this->ResolveAlias("sel");
        qcd_sdb_cat = ttqcdcat + "&&" + alias_map_["tt_qcd_norm"];
      }
      
      //Work out the W norm in the SS region separately if it is data driven
      Value w_ss_norm;
      std::vector<std::string> wjets_samples = this->ResolveSamplesAlias("wjets_samples");
      if(ch_ != channel::em) {
         if(method == 10 || method == 11) { 
          w_ss_norm = this->GetRateViaWMethod(wjets_samples, qcd_cat, w_extrp_sdb_sel, w_extrp_sig_sel, 
              this->ResolveSamplesAlias("data_samples"), qcd_cat, w_sdb_sel, w_sub_samples, wt, ValueFnMap());
         }else if(method == 12 || method == 13 || method == 14){
          w_ss_norm = this->GetRateViaWOSSSMethod(wjets_samples, qcd_cat, w_extrp_sdb_sel, w_extrp_sig_sel, 
              this->ResolveSamplesAlias("data_samples"), qcd_cat, w_sdb_sel_osss, w_sub_samples, false, wt, ValueFnMap());
         }else if(method==16){
          //std::string qcd_cat_nobtag = "n_jets<=1 && "+alias_map_["baseline"];
          std::string qcd_cat_nobtag = "n_jets<=1 && n_lowpt_jets>= 1 && "+alias_map_["baseline"];
          Value w_ss_norm_nobtag = this->GetRateViaWOSSSMethod(wjets_samples, qcd_cat_nobtag, w_extrp_sdb_sel, w_extrp_sig_sel,
             this->ResolveSamplesAlias("data_samples"),qcd_cat_nobtag, w_sdb_sel_osss,w_sub_samples,false, wt, ValueFnMap());
         std::string os_sel_btag = "!os && " + this->ResolveAlias("sel");
         Value btag_norm = GetRate(wjets_samples, os_sel_btag, qcd_cat, wt);
         Value nobtag_norm = GetRate(wjets_samples,os_sel_btag,qcd_cat_nobtag,wt);
         Value btag_extrap_factor = ValueDivide(btag_norm,nobtag_norm);
         if (verbosity_) PrintValue("ExtrapFactor into full b-tag definition",btag_extrap_factor);
         w_ss_norm = ValueProduct(btag_extrap_factor,w_ss_norm_nobtag);
        } else {
          w_ss_norm = std::make_pair(0.0,0.0);
        }
      } else {
        w_ss_norm = std::make_pair(0.0,0.0);
      }
     
     std::map<std::string,std::function<Value()>> wjets_ss_vals;
     for(unsigned i=0; i< wjets_samples.size() ; ++i){ //We have wjets_samples.size() w+jets samples to subtract, but only one data driven norm
      wjets_ss_vals[wjets_samples.at(i)] = [&]()->HTTRun2Analysis::Value{return std::make_pair(w_ss_norm.first/(wjets_samples.size()),w_ss_norm.second/std::sqrt(wjets_samples.size()));};
    }
       
      if(ch_ != channel::tt){
        if(method == 8 || method == 9) {
          qcd_norm = this->GetRateViaQCDMethod(std::make_pair(qcd_os_ss_factor_,0.), this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_cat, qcd_sub_samples, wt, ValueFnMap());
        } else if (method == 10 || method == 11 || method == 12 || method == 13 || method == 14 || method == 16) {
          qcd_norm = this->GetRateViaQCDMethod(std::make_pair(qcd_os_ss_factor_,0.), this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_cat, qcd_sub_samples, wt,//{
          wjets_ss_vals);
        } else if (method == 15 && !do_ss_){
         qcd_norm = this->GetRateViaQCDMethod(std::make_pair(1.0,0.), this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_cat, qcd_sub_samples, wt+"*wt_em_qcd",ValueFnMap());
        } else if (method ==15 && do_ss_){
         qcd_norm = this->GetRateViaQCDMethod(std::make_pair(1.0,0.), this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_cat, qcd_sub_samples, wt,ValueFnMap());
        }
      } else {
        if(method == 8){
          qcd_norm = this->GetRateViaTauTauQCDMethod(this->ResolveSamplesAlias("data_samples"),qcd_sdb_cat,qcd_cat,qcd_extrap_sel,qcd_sdb_sel,qcd_sdb_cat, qcd_sub_samples, wt, ValueFnMap());
        } else if (method == 10){
          qcd_norm = this->GetRateViaTauTauQCDMethod(this->ResolveSamplesAlias("data_samples"),qcd_sdb_cat,qcd_cat,qcd_extrap_sel,qcd_sdb_sel,qcd_sdb_cat, qcd_sub_samples, wt,//{
          wjets_ss_vals);
       }
     }
      if (qcd_norm.first <= 0.0) {
        double default_rate = 0.0000001;
        std::cout << "[HTTRun2Analysis::GenerateQCD] Warning, QCD rate is negative (" 
          << qcd_norm.first << "), setting to " << default_rate << " and maintaining error" << std::endl;
        qcd_norm.first = default_rate;
      }
      std::string qcd_shape_cat = cat;
      if (ch_ != channel::tt){
        if (method == 8) {
          qcd_shape_cat += "&&" + alias_map_["baseline"];
          qcd_hist = this->GetShapeViaQCDMethod(var, this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_shape_cat, qcd_sub_samples, wt, ValueFnMap());
        } else if(method == 10 || method==12 || method == 14) {
          qcd_shape_cat += "&&" + alias_map_["baseline"];        
            if(method == 14) qcd_shape_cat = "n_jets<=1 && n_loose_bjets>=1 &&" +alias_map_["baseline"];
          qcd_hist = this->GetShapeViaQCDMethod(var, this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_shape_cat, qcd_sub_samples, wt, //{
          wjets_ss_vals);
        } else if (method == 15 &&  !do_ss_){
          qcd_shape_cat += "&&" + alias_map_["baseline"];
          qcd_hist = this->GetShapeViaQCDMethod(var, this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_shape_cat, qcd_sub_samples, wt+"*wt_em_qcd",ValueFnMap());
        } else if (method == 15 && do_ss_){
          qcd_shape_cat += "&&" + alias_map_["baseline"];
          qcd_hist = this->GetShapeViaQCDMethod(var, this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_shape_cat, qcd_sub_samples, wt,ValueFnMap());
        } else {
          if (method == 9 || method == 11 ||method == 13) {
              qcd_shape_cat += "&&" + this->ResolveAlias("qcd_loose_shape");
          }
          qcd_hist = this->GetShape(var, this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_shape_cat, wt);
//        if (verbosity_) std::cout << "Shape: " << boost::format("%s,'%s','%s','%s'\n")
 //         % this->ResolveSamplesAlias("data_samples") % qcd_sdb_sel % qcd_shape_cat % wt;
        }
     } else {
     if (method == 8){ 
       qcd_shape_cat += "&&" + alias_map_["tt_qcd_norm"];
       qcd_hist = this->GetShapeViaQCDMethod(var, this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_shape_cat, qcd_sub_samples, wt, ValueFnMap());
     } else if (method ==10){
       qcd_shape_cat += "&&" + alias_map_["tt_qcd_norm"];
          qcd_hist = this->GetShapeViaQCDMethod(var, this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, qcd_shape_cat, qcd_sub_samples, wt, //{
          wjets_ss_vals);
       }
    }

    if(verbosity_ > 1) {
      std::string bin_delim="";
      if(var.find("[")!=std::string::npos) bin_delim = "[";
      if(var.find("(")!=std::string::npos) bin_delim = "(";
      std::string default_cat = cat;
      default_cat  += "&&" +alias_map_["baseline"];
      this->KolmogorovTest(var.substr(0, var.find(bin_delim, 0)), (this->ResolveSamplesAlias("data_samples")).at(0), qcd_sdb_sel, default_cat,(this->ResolveSamplesAlias("data_samples")).at(0), qcd_sdb_sel, qcd_shape_cat,wt); //Not exactly equivalent to comparison with default norm as below, but good to give an idea
    //TH1F default_qcd_hist = this->GetShapeViaQCDMethod(var, this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, cat, qcd_sub_samples, wt, ValueFnMap());
      TH1F default_qcd_hist = this->GetShapeViaQCDMethod(var,this->ResolveSamplesAlias("data_samples"), qcd_sdb_sel, default_cat, qcd_sub_samples, wt, {
        {wjets_samples.at(0), [&]()->HTTRun2Analysis::Value {
          return w_ss_norm;}
         }
        });
      std::cout << "ROOT KS test: "<< default_qcd_hist.KolmogorovTest(&qcd_hist) <<std::endl;
    }

    SetNorm(&qcd_hist, qcd_norm.first);
    return std::make_pair(qcd_hist, qcd_norm);
   }

  HTTRun2Analysis::HistValuePair HTTRun2Analysis::GenerateSignal(std::string sample, std::string var, std::string sel, std::string cat, std::string wt, double xs) {
    cat += "&&" + alias_map_["baseline"];
    Value signal_norm;
    if (xs > 0) {
      if (verbosity_ > 2) std::cout << "[HTTRun2Analysis::GenerateSignal] " << sample << " scaled to lumi using cross section " << xs << " pb" << std::endl;
      signal_norm = GetRate(sample, sel, cat, wt);
      signal_norm = ValueProduct(signal_norm, std::make_pair(this->GetLumiScaleFixedXS(sample, xs), 0.0));
    } else {
      signal_norm = GetLumiScaledRate(sample, sel, cat, wt);
    }
    TH1F signal_shape = this->GetShape(var, sample, sel, cat, wt);
    SetNorm(&signal_shape, signal_norm.first);
    return std::make_pair(signal_shape, signal_norm);
  }

  void HTTRun2Analysis::FillSMSignal(HistValueMap & hmap, 
                    std::vector<std::string> const& masses,
                    std::string const& var,
                    std::string const& sel,
                    std::string const& cat,
                    std::string const& wt,
                    std::string const& infix,
                    std::string const& postfix,
                    double fixed_xs) {
    for (auto const& m : masses) {
        hmap["ggH"+infix+m+postfix] = this->GenerateSignal("GluGluHToTauTau_M-"+m,    var, sel, cat, wt, fixed_xs);
        hmap["qqH"+infix+m+postfix] = this->GenerateSignal("VBFHToTauTau_M-"+m,        var, sel, cat, wt, fixed_xs);
        hmap["WplusH"+infix+m+postfix] = this->GenerateSignal("WplusHToTauTau_M-"+m,        var, sel, cat, wt, fixed_xs);
        hmap["WminusH"+infix+m+postfix] = this->GenerateSignal("WminusHToTauTau_M-"+m,        var, sel, cat, wt, fixed_xs);
        hmap["ZH"+infix+m+postfix] = this->GenerateSignal("ZHToTauTau_M-"+m,        var, sel, cat, wt, fixed_xs);
        hmap["TTH"+infix+m+postfix] = this->GenerateSignal("TTHToTauTau_M-"+m,        var, sel, cat, wt, fixed_xs);
    }
  }
  

  void HTTRun2Analysis::FillMSSMSignal(HistValueMap & hmap, 
                    std::vector<std::string> const& masses,
                    std::string const& var,
                    std::string const& sel,
                    std::string const& cat,
                    std::string const& wt,
                    std::string const& infix,
                    std::string const& postfix,
                    double fixed_xs) {
    for (auto const& m : masses) {
        hmap["ggH"+infix+m+postfix] = this->GenerateSignal("SUSYGluGluToHToTauTau_M-"+m, var, sel, cat, wt, fixed_xs);
        hmap["bbH"+infix+m+postfix] = this->GenerateSignal("SUSYGluGluToBBHToTauTau_M-"+m, var, sel, cat, wt, fixed_xs);
    }
  }

  void HTTRun2Analysis::FillMSSMSignalBBH(HistValueMap & hmap, 
                    std::vector<std::string> const& masses,
                    std::string const& var,
                    std::string const& sel,
                    std::string const& cat,
                    std::string const& wt,
                    std::string const& infix,
                    std::string const& postfix,
                    double fixed_xs) {
    for (auto const& m : masses) {
        hmap["bbH"+infix+m+postfix] = this->GenerateSignal("SUSYGluGluToBBHToTauTau_M-"+m, var, sel, cat, wt, fixed_xs);
    }
  }

  void HTTRun2Analysis::FillMSSMSignalGGH(HistValueMap & hmap, 
                    std::vector<std::string> const& masses,
                    std::string const& var,
                    std::string const& sel,
                    std::string const& cat,
                    std::string const& wt,
                    std::string const& infix,
                    std::string const& postfix,
                    double fixed_xs) {
    for (auto const& m : masses) {
        hmap["ggH"+infix+m+postfix] = this->GenerateSignal("SUSYGluGluToHToTauTau_M-"+m, var, sel, cat, wt, fixed_xs);
    }
  }




  void HTTRun2Analysis::FillHhhSignal(HistValueMap & hmap, 
                    std::vector<std::string> const& masses,
                    std::string const& var,
                    std::string const& sel,
                    std::string const& cat,
                    std::string const& wt,
                    std::string const& infix,
                    std::string const& postfix,
                    double fixed_xs) {
    for (auto const& m : masses) {
        hmap["ggH"+infix+m+postfix] = this->GenerateSignal("GluGluToRadionToHHTo2B2Tau_M-"+m, var, sel, cat, wt, fixed_xs);
    }
  }


  void HTTRun2Analysis::FillHistoMap(HistValueMap & hmap, unsigned method,
                        std::string var,
                        std::string sel,
                        std::string cat,
                        std::string wt,
                        std::string postfix) {
    Value total_bkr;
    // Data
    auto data_pair = this->GenerateData(method, var, sel, cat, wt);
    PrintValue("data_obs"+postfix, data_pair.second);
    hmap["data_obs"+postfix] = data_pair;
    //Splitting TT into TTT and TTJ - let's disable this for now
    std::string ttt_sel, ttj_sel;
    ttt_sel = sel+"&&"+this->ResolveAlias("ztt_sel");
    ttj_sel = sel+"&&!"+this->ResolveAlias("ztt_sel");
    if(ch_ == channel::zmm || ch_==channel::zee){
      ttt_sel = "0";
      ttj_sel = sel;
    }
    auto topt_pair = this->GenerateTOP(method, var, ttt_sel, cat, wt);
    auto topj_pair = this->GenerateTOP(method, var, ttj_sel, cat, wt);
    PrintValue("TTT"+postfix, topt_pair.second);
    PrintValue("TTJ"+postfix, topj_pair.second);
    hmap["TTT"+postfix] = topt_pair;
    hmap["TTJ"+postfix] = topj_pair;
    //auto top_pair = this->GenerateTOP(method, var, sel, cat, wt);
    std::string top_map_label = "TT";
    Value tt_norm = ValueAdd(topt_pair.second, topj_pair.second);
    TH1F tt_hist = topt_pair.first;
    PrintValue(top_map_label+postfix, tt_norm);
    tt_hist.Add(&topj_pair.first);
    /*PrintValue(top_map_label+postfix, top_pair.second);
    total_bkr = ValueAdd(total_bkr, top_pair.second);*/
    total_bkr = ValueAdd(total_bkr, tt_norm);
    hmap[top_map_label+postfix] = std::make_pair(tt_hist,tt_norm);
    TH1F total_hist = hmap[top_map_label+postfix].first; 
    // Diboson
    auto vvt_pair = this->GenerateVV(method, var, ttt_sel, cat, wt);
    auto vvj_pair = this->GenerateVV(method, var, ttj_sel, cat, wt);
    PrintValue("VVT"+postfix, vvt_pair.second);
    PrintValue("VVJ"+postfix, vvj_pair.second);
    hmap["VVT"+postfix] = vvt_pair;
    hmap["VVJ"+postfix] = vvj_pair;
    //auto vv_pair = this->GenerateVV(method, var, sel, cat, wt);
    std::string vv_map_label =  "VV";
    Value vv_norm = ValueAdd(vvt_pair.second, vvj_pair.second);
    TH1F vv_hist = vvt_pair.first;
    PrintValue(vv_map_label+postfix, vv_norm);
    vv_hist.Add(&vvj_pair.first);
    /*PrintValue(vv_map_label+postfix, vv_pair.second);
    total_bkr = ValueAdd(total_bkr, vv_pair.second);*/
    total_bkr = ValueAdd(total_bkr, vv_norm);
    hmap[vv_map_label+postfix] = std::make_pair(vv_hist, vv_norm);
    total_hist.Add(&hmap[vv_map_label+postfix].first,1.0);
    // Z->ll
    if (ch_ != channel::em && ch_!= channel::zee && ch_!= channel::zmm && ch_!= channel::tpzee && ch_!= channel::tpzmm && ch_!=channel::wmnu) {
      std::string zl_sel, zj_sel;
      zl_sel= sel+"&&"+this->ResolveAlias("zl_sel");
      zj_sel= sel+"&&"+this->ResolveAlias("zj_sel");
      auto zl_pair = this->GenerateZL(method, var, zl_sel, cat, wt);
      auto zj_pair = this->GenerateZJ(method, var, zj_sel, cat, wt);
      Value zll_norm = ValueAdd(zl_pair.second, zj_pair.second);
      TH1F zll_hist = zl_pair.first;
      zll_hist.Add(&zj_pair.first);
      PrintValue("ZLL"+postfix, zll_norm);
      PrintValue("ZL"+postfix, zl_pair.second);
      PrintValue("ZJ"+postfix, zj_pair.second);
      total_bkr = ValueAdd(total_bkr, zll_norm);
      hmap["ZLL"+postfix] = std::make_pair(zll_hist, zll_norm);
      hmap["ZL"+postfix]  = zl_pair;
      hmap["ZJ"+postfix]  = zj_pair;
      total_hist.Add(&hmap["ZLL"+postfix].first,1.0);
    } else {
      std::string zll_sel;
      zll_sel = sel+"&&"+this->ResolveAlias("zll_sel");
      auto zll_pair = this->GenerateZLL(method, var, zll_sel, cat, wt);
      std::string zll_map_label = "ZLL";
      PrintValue(zll_map_label+postfix, zll_pair.second);
      total_bkr = ValueAdd(total_bkr, zll_pair.second);
      hmap[zll_map_label+postfix] = zll_pair;
      total_hist.Add(&hmap[zll_map_label+postfix].first,1.0);
    }
    // Z->tautau
    if(ch_!= channel::zee && ch_!= channel::zmm && ch_!= channel::tpzee && ch_!= channel::tpzmm && ch_!=channel::wmnu) {
      std::string ztt_sel;
      ztt_sel = sel+"&&"+this->ResolveAlias("ztt_sel");
      auto ztt_pair = this->GenerateZTT(method, var, ztt_sel, cat, wt);
      std::string ztt_map_label = "ZTT";
      //std::string ztt_map_label = (ch_ == channel::em) ? "Ztt" : "ZTT";
      PrintValue(ztt_map_label+postfix, ztt_pair.second);
      total_bkr = ValueAdd(total_bkr, ztt_pair.second);
      hmap[ztt_map_label+postfix] = ztt_pair;
      total_hist.Add(&hmap[ztt_map_label+postfix].first,1.0);
    }    // W+jets
    if (ch_ != channel::em) {
      auto w_pair = this->GenerateW(method, var, sel, cat, wt);
      PrintValue("W"+postfix, w_pair.second);
      total_bkr = ValueAdd(total_bkr, w_pair.second);
      hmap["W"+postfix] = w_pair;
      total_hist.Add(&hmap["W"+postfix].first,1.0);
    }
    else if (ch_ == channel::em) {
      auto w_pair = this->GenerateW(method, var, sel, cat, wt);
      PrintValue("WJets"+postfix, w_pair.second);
      total_bkr = ValueAdd(total_bkr, w_pair.second);
      hmap["WJets"+postfix] = w_pair;
      total_hist.Add(&hmap["WJets"+postfix].first,1.0);
      
      auto wgam_pair = this->GenerateWGamma(method, var, sel, cat, wt);
      PrintValue("WGam"+postfix, wgam_pair.second);
      total_bkr = ValueAdd(total_bkr, wgam_pair.second);
      hmap["WGam"+postfix] = wgam_pair;
      total_hist.Add(&hmap["WGam"+postfix].first,1.0);
      
      TH1F total_W = hmap["WJets"+postfix].first; 
      Value total_W_norm;
      total_W_norm = ValueAdd(total_W_norm, w_pair.second);
      total_W_norm = ValueAdd(total_W_norm, wgam_pair.second);
      total_W.Add(&hmap["WGam"+postfix].first,1.0);
      hmap["W"+postfix] = std::make_pair(total_W,total_W_norm);
      PrintValue("W"+postfix, total_W_norm);
    }
    // QCD/Fakes
    if(ch_!= channel::tpzee && ch_!= channel::tpzmm && ch_!=channel::wmnu) {
      auto qcd_pair = this->GenerateQCD(method, var, sel, cat, wt);
      std::string qcd_map_label = "QCD";
      //std::string qcd_map_label = (ch_ == channel::em) ? "Fakes" : "QCD";
      PrintValue(qcd_map_label+postfix, qcd_pair.second);
      total_bkr = ValueAdd(total_bkr, qcd_pair.second);
      hmap[qcd_map_label+postfix] = qcd_pair;
      total_hist.Add(&hmap["QCD"+postfix].first,1.0);
    }
    // Print the total background yield
    PrintValue("Total"+postfix, total_bkr);
    //Until there is data, fill the data with sum of the backgrounds
    hmap["total_bkg"+postfix] = std::make_pair(total_hist,total_bkr);
//    hmap["data_obs"+postfix] = std::make_pair(total_hist,total_bkr);
    return;
  }

  std::string HTTRun2Analysis::BuildCutString(std::string const& selection,
      std::string const& category,
      std::string const& weight) {
    std::string full_selection;
    if (weight != "" && (selection != "" || category != "")) full_selection += "( ";

    if (selection != "")                    full_selection += ("(" + selection + ")");
    if (selection != "" && category != "")  full_selection += " && ";
    if (category != "")                     full_selection += ("(" + category + ")");
    if (weight != "" && (selection != "" || category != "")) full_selection += " ) * ";
    if (weight != "") full_selection += ("("+weight+")");
    return full_selection;                                      
  }

  std::string HTTRun2Analysis::BuildVarString(std::string const& variable) {
    std::string full_variable = variable;
    if (full_variable.find_last_of("(") != full_variable.npos 
        && full_variable.find("[") == full_variable.npos
        && full_variable.find("]") == full_variable.npos) {
      full_variable.insert(full_variable.find_last_of("("),">>htemp");
    }
    return full_variable;
  }


  TH1F HTTRun2Analysis::GetShape(std::string const& variable,
                                       std::string const& sample, 
                                       std::string const& selection, 
                                       std::string const& category, 
                                       std::string const& weight) {
    TH1::SetDefaultSumw2(true);
    std::string full_variable = BuildVarString(variable);
    std::size_t begin_var = full_variable.find("[");
    std::size_t end_var   = full_variable.find("]");
    TH1F *htemp = nullptr;
    if (begin_var != full_variable.npos && end_var != full_variable.npos) {
      std::string binning = full_variable.substr(begin_var+1, end_var-begin_var-1);
      std::vector<std::string> string_vec;
      boost::split(string_vec, binning, boost::is_any_of(","));
      std::vector<double> bin_vec;
      for (auto str : string_vec) bin_vec.push_back(boost::lexical_cast<double>(str));
      TH1::AddDirectory(true);
      htemp = new TH1F("htemp","htemp", bin_vec.size()-1, &(bin_vec[0]));
      TH1::AddDirectory(false);
      full_variable.erase(begin_var, full_variable.npos);
      full_variable += ">>htemp";
    }
    std::string full_selection = BuildCutString(selection, category, weight);
    // std::cout << full_selection << std::endl;
    // std::cout << full_variable << std::endl;
    TH1::AddDirectory(true);
    //In the case of an empty sample, can return htemp. This is only created above if the []
    //binning is used, so we have to create an empty version in the case of () binning.
    if(ttrees_[sample]->GetEntries() == 0) {
      std::size_t begin_var = full_variable.find_last_of("(");
      std::size_t end_var   = full_variable.find_last_of(")");
      if (begin_var != full_variable.npos && end_var != full_variable.npos) {
        std::string binning = full_variable.substr(begin_var+1, end_var-begin_var-1);
        std::vector<std::string> string_vec;
        boost::split(string_vec, binning, boost::is_any_of(","));
        std::vector<double> bin_vec;
        for (auto str : string_vec) bin_vec.push_back(boost::lexical_cast<double>(str));
        htemp = new TH1F("htemp","htemp", bin_vec[0],bin_vec[1],bin_vec[bin_vec.size()-1]);
      }
      return (*htemp);
    }
    ttrees_[sample]->Draw(full_variable.c_str(), full_selection.c_str(), "goff");
    TH1::AddDirectory(false);
    htemp = (TH1F*)gDirectory->Get("htemp");
    TH1F result = (*htemp);
    gDirectory->Delete("htemp;*");
    auto rate = GetRate(sample, selection, category, weight);
    SetNorm(&result, rate.first);
    if(result.Integral(1,result.GetNbinsX()) == 0) std::cout<<"Warning - no shape for sample "<<sample<<std::endl;
    return result;
  }

  TH1F HTTRun2Analysis::GetShape(std::string const& variable,
                                       std::vector<std::string> const& samples, 
                                       std::string const& selection, 
                                       std::string const& category, 
                                       std::string const& weight) {

    TH1F result = GetShape(variable, samples.at(0), selection, category, weight);
    if (samples.size() > 1) {
      for (unsigned i = 1; i < samples.size(); ++i) {
        TH1F tmp = GetShape(variable, samples.at(i), selection, category, weight);
        result.Add(&tmp);
      }
    }
    return result;
  }


  TH1F HTTRun2Analysis::GetLumiScaledShape(std::string const& variable,
                                       std::string const& sample, 
                                       std::string const& selection, 
                                       std::string const& category, 
                                       std::string const& weight) {
    TH1F result = GetShape(variable, sample, selection, category, weight);
    result.Scale(GetLumiScale(sample));
    return result;
  }

  TH1F HTTRun2Analysis::GetLumiScaledShape(std::string const& variable,
                                       std::vector<std::string> const& samples, 
                                       std::string const& selection, 
                                       std::string const& category, 
                                       std::string const& weight) {

    TH1F result = GetLumiScaledShape(variable, samples.at(0), selection, category, weight);
    if (samples.size() > 1) {
      for (unsigned i = 1; i < samples.size(); ++i) {
        TH1F tmp = GetLumiScaledShape(variable, samples.at(i), selection, category, weight);
        result.Add(&tmp);
      }
    }
    return result;
  }


  std::pair<double, double> HTTRun2Analysis::GetRate(std::string const& sample, 
                                      std::string const& selection, 
                                      std::string const& category, 
                                      std::string const& weight) {
    if(verbosity_>2){ std::cout << "--GetRate-- Sample:\"" << sample << "\" Selection:\"" << selection << "\" Category:\"" 
      << category << "\" Weight:\"" << weight << "\"" << std::endl;}
    std::string full_selection = BuildCutString(selection, category, weight);
    TH1::AddDirectory(true);
    //If the tree is empty, return 0
    if(ttrees_[sample]->GetEntries() == 0) return std::make_pair(0,0);
    ttrees_[sample]->Draw("0.5>>htemp(1,0,1)", full_selection.c_str(), "goff");
    TH1::AddDirectory(false);
    TH1F *htemp = (TH1F*)gDirectory->Get("htemp");
    auto result = std::make_pair(Integral(htemp), Error(htemp));
    gDirectory->Delete("htemp;*");
    return result;
  }

  std::pair<double, double> HTTRun2Analysis::GetRate(std::vector<std::string> const& samples, 
                                      std::string const& selection, 
                                      std::string const& category, 
                                      std::string const& weight) {
    auto result = GetRate(samples.at(0),selection,category,weight);
    double err_sqr = result.second*result.second;
    if(samples.size() > 1){
      for(unsigned i = 1; i < samples.size(); ++i){
        auto tmp = GetRate(samples.at(i),selection,category,weight);
        result.first += tmp.first;
        err_sqr += (tmp.second*tmp.second);
      }
    }
    result.second = sqrt(err_sqr);
    return result;
  }


  std::pair<double, double> HTTRun2Analysis::GetLumiScaledRate(std::string const& sample, 
                                      std::string const& selection, 
                                      std::string const& category, 
                                      std::string const& weight) {
    auto result = GetRate(sample, selection, category, weight);
    double sf = GetLumiScale(sample);
    result.first *= sf;
    result.second *= sf;
    if(verbosity_ >0 ){
      PrintValue(sample,result);
    }
    return result;
  }
  std::pair<double, double> HTTRun2Analysis::GetLumiScaledRate(std::vector<std::string> const& samples, 
                                      std::string const& selection, 
                                      std::string const& category, 
                                      std::string const& weight) {
    auto result = GetLumiScaledRate(samples.at(0), selection, category, weight);
    double err_sqr = result.second * result.second;
    if (samples.size() > 1) {
      for (unsigned i = 1; i < samples.size(); ++i) {
        auto tmp = GetLumiScaledRate(samples.at(i), selection, category, weight);
        result.first += tmp.first;
        err_sqr += (tmp.second * tmp.second);
      }
    }
    result.second = sqrt(err_sqr);
    return result;
  }

  std::pair<double, double> HTTRun2Analysis::SampleEfficiency(std::vector<std::string> const& samples, 
                          std::string const& ref_selection, 
                          std::string const& ref_category,
                          std::string const& target_selection, 
                          std::string const& target_category,  
                          std::string const& weight) {
    auto num = GetRate(samples, target_selection, target_category, weight);
    auto den = GetRate(samples, ref_selection, ref_category, weight);
    double num_eff = std::pow(num.first / num.second, 2.0) ;
    unsigned num_eff_rounded = unsigned(num_eff+0.5);
    double den_eff = std::pow(den.first / den.second, 2.0) ;
    unsigned den_eff_rounded = unsigned(den_eff+0.5);
    double eff = num.first / den.first;
    TEfficiency teff;
    double eff_err_up   = teff.ClopperPearson(den_eff_rounded,num_eff_rounded,0.683,1)-(num_eff/den_eff);
    double eff_err_down = (num_eff/den_eff)-teff.ClopperPearson(den_eff_rounded,num_eff_rounded,0.683,0);
    double eff_err = (eff_err_up/(num_eff/den_eff)) * eff;
    if (num.first == 0.0) {
      std::cout << "[HTTRun2Analysis::SampleEfficiency] Numerator is zero, setting error to zero" << std::endl;
      eff_err = 0.0;
    }
    auto result = std::make_pair(eff, eff_err);
    if (verbosity_ > 0) {
      std::cout << "[HTTRun2Analysis::SampleEfficiency]" << std::endl;
      std::cout << "Numerator:   "; 
      if (samples.size() > 1) {
        for (unsigned i = 1; i < samples.size(); ++i) {
            std::cout << samples[i] << " ";
        }
      }
      std::cout <<  target_selection << " " << target_category << " " << weight << std::endl;
      std::cout << "Denominator: "; 
      if (samples.size() > 1) {
        for (unsigned i = 1; i < samples.size(); ++i) {
            std::cout << samples[i] << " ";
        }
      }
      std::cout <<  ref_selection << " " << ref_category << " " << weight << std::endl;
      PrintValue("Numerator",num);
      PrintValue("Denominator",den);
      std::cout << "Effective Numerator:   " << num_eff_rounded << std::endl;
      std::cout << "Effective Denominator: " << den_eff_rounded << std::endl;
      std::cout << "Error down (relative): " << eff_err_down/(num_eff/den_eff) << std::endl;
      std::cout << "Error up   (relative): " << eff_err_up/(num_eff/den_eff) << std::endl;
    }
    return result;
  }
  
  std::pair<double, double> HTTRun2Analysis::SampleEfficiency(std::string const& sample, 
                          std::string const& ref_selection, 
                          std::string const& ref_category,
                          std::string const& target_selection, 
                          std::string const& target_category,  
                          std::string const& weight) {
    auto num = GetRate(sample, target_selection, target_category, weight);
    auto den = GetRate(sample, ref_selection, ref_category, weight);
    double num_eff = std::pow(num.first / num.second, 2.0) ;
    unsigned num_eff_rounded = unsigned(num_eff+0.5);
    double den_eff = std::pow(den.first / den.second, 2.0) ;
    unsigned den_eff_rounded = unsigned(den_eff+0.5);
    double eff = num.first / den.first;
    TEfficiency teff;
    double eff_err_up   = teff.ClopperPearson(den_eff_rounded,num_eff_rounded,0.683,1)-(num_eff/den_eff);
    double eff_err_down = (num_eff/den_eff)-teff.ClopperPearson(den_eff_rounded,num_eff_rounded,0.683,0);
    double eff_err = (eff_err_up/(num_eff/den_eff)) * eff;
    if (num.first == 0.0) {
      std::cout << "[HTTRun2Analysis::SampleEfficiency] Numerator is zero, setting error to zero" << std::endl;
      eff_err = 0.0;
    }
    auto result = std::make_pair(eff, eff_err);
    if (verbosity_ > 0) {
      std::cout << "[HTTRun2Analysis::SampleEfficiency]" << std::endl;
      std::cout << "Numerator:   " << boost::format("%s,'%s','%s','%s'\n") % sample % target_selection
                % target_category % weight;
      std::cout << "Denominator: " << boost::format("%s,'%s','%s','%s'\n") % sample % ref_selection
                % ref_category % weight;
      PrintValue("Numerator",num);
      PrintValue("Denominator",den);
      std::cout << "Effective Numerator:   " << num_eff_rounded << std::endl;
      std::cout << "Effective Denominator: " << den_eff_rounded << std::endl;
      std::cout << "Error down (relative): " << eff_err_down/(num_eff/den_eff) << std::endl;
      std::cout << "Error up   (relative): " << eff_err_up/(num_eff/den_eff) << std::endl;
    }
    return result;
  }


  HTTRun2Analysis::Value HTTRun2Analysis::SampleRatio(std::vector<std::string> const& sample, 
                          std::string const& ref_selection, 
                          std::string const& ref_category,
                          std::string const& target_selection, 
                          std::string const& target_category,  
                          std::string const& weight) {
    Value num = GetRate(sample, target_selection, target_category, weight);
    Value den = GetRate(sample, ref_selection, ref_category, weight);
    return ValueDivide(num, den);
  }

  HTTRun2Analysis::Value HTTRun2Analysis::DataSampleRatio(std::vector<std::string> const& sample, 
                          std::string const& ref_selection, 
                          std::string const& ref_category,
                          std::string const& target_selection, 
                          std::string const& target_category,  
                          std::string const& weight,
                          std::vector<std::string> const& subtr_samples,
                          std::map<std::string, std::function<Value()>> func_dict) {
    Value num = GetRate(sample, target_selection, target_category, weight);
    Value total_bkg_num;
    for (unsigned i = 0; i < subtr_samples.size(); ++i) {
      Value bkr;
      if (func_dict.count(subtr_samples[i])) {
        bkr = ((*func_dict.find(subtr_samples[i])).second)(); // find and evaluate function
      } else {
        bkr = GetLumiScaledRate(subtr_samples[i], target_selection, target_category, weight);
      }
      if (verbosity_) PrintValue("-"+subtr_samples[i], bkr);
      double new_err = std::sqrt((total_bkg_num.second * total_bkg_num.second) + (bkr.second * bkr.second));
      total_bkg_num.first += bkr.first;
      total_bkg_num.second = new_err;
    }
    double ratio_num_err = std::sqrt((total_bkg_num.second * total_bkg_num.second) + (num.second * num.second));
    Value bkg_sub_num(num.first - total_bkg_num.first, ratio_num_err);
    Value den = GetRate(sample, ref_selection, ref_category, weight);
    Value total_bkg_den;
    for (unsigned i = 0; i < subtr_samples.size(); ++i) {
      Value bkr;
      if (func_dict.count(subtr_samples[i])) {
        bkr = ((*func_dict.find(subtr_samples[i])).second)(); // find and evaluate function
      } else {
        bkr = GetLumiScaledRate(subtr_samples[i], ref_selection, ref_category, weight);
      }
      if (verbosity_) PrintValue("-"+subtr_samples[i], bkr);
      double new_err = std::sqrt((total_bkg_den.second * total_bkg_den.second) + (bkr.second * bkr.second));
      total_bkg_den.first += bkr.first;
      total_bkg_den.second = new_err;
    }
    double ratio_den_err = std::sqrt((total_bkg_den.second * total_bkg_den.second) + (den.second * den.second));
    Value bkg_sub_den(den.first - total_bkg_den.first, ratio_den_err);
    return ValueDivide(bkg_sub_num, bkg_sub_den);
  }

  HTTRun2Analysis::Value HTTRun2Analysis::GetRateViaRefEfficiency(std::vector<std::string> const& target_samples, 
                          std::vector<std::string> const& ref_samples,
                          std::string const& ref_selection, 
                          std::string const& ref_category,
                          std::string const& target_selection, 
                          std::string const& target_category,  
                          std::string const& weight) {
    if (verbosity_) {
      std::cout << "[HTTRun2Analysis::GetRateViaRefEfficiency]\n";
      std::cout << "ReferenceRate:   " << ref_selection << " "; 
      if (ref_samples.size() > 1) {
        for (unsigned i = 1; i < ref_samples.size(); ++i) {
            std::cout << ref_samples[i] << " ";
        }
      }
      std::cout <<  ref_category << " " << weight << std::endl;;
      std::cout << "Efficiency:   " << target_selection << " " ;
      if (target_samples.size() > 1) {
        for (unsigned i = 1; i < target_samples.size(); ++i) {
            std::cout << target_samples[i] << " " ;
        }
      }
      std::cout <<  target_category << " " << weight << std::endl;;
    }
    auto ref_rate = GetLumiScaledRate(ref_samples, ref_selection, ref_category, weight);
    if (verbosity_) PrintValue("ReferenceRate", ref_rate);
    auto target_eff = SampleEfficiency(target_samples, ref_selection, ref_category, target_selection, target_category, weight);
    if (verbosity_) PrintValue("Efficiency", target_eff);
    return ValueProduct(ref_rate, target_eff);
  }

  HTTRun2Analysis::Value HTTRun2Analysis::GetRateViaRefEfficiency(std::string const& target_sample, 
                          std::string const& ref_sample,
                          std::string const& ref_selection, 
                          std::string const& ref_category,
                          std::string const& target_selection, 
                          std::string const& target_category,  
                          std::string const& weight) {
    if (verbosity_) {
      std::cout << "[HTTRun2Analysis::GetRateViaRefEfficiency]\n";
      std::cout << "ReferenceRate:   " << boost::format("%s,'%s','%s','%s'\n") % ref_sample % ref_selection 
                % ref_category % weight;
      std::cout << "Efficiency:      " << boost::format("%s,'%s','%s','%s'\n") %target_sample % target_selection 
                % target_category % weight;
    }
    auto ref_rate = GetLumiScaledRate(ref_sample, ref_selection, ref_category, weight);
    if (verbosity_) PrintValue("ReferenceRate", ref_rate);
    auto target_eff = SampleEfficiency(target_sample, ref_selection, ref_category, target_selection, target_category, weight);
    if (verbosity_) PrintValue("Efficiency", target_eff);
    return ValueProduct(ref_rate, target_eff);
  }

  HTTRun2Analysis::Value HTTRun2Analysis::GetRateViaFakesMethod(std::string const& sel,
                              std::string const& cat, 
                              std::string const& wt) {
    auto e_fakes = this->GetRate("Special_20_Data", sel, cat, wt);
    auto m_fakes = this->GetRate("Special_21_Data", sel, cat, wt);
    auto em_fakes = this->GetRate("Special_22_Data", sel, cat, wt);
    if (verbosity_) PrintValue("Electron Fakes", e_fakes);
    if (verbosity_) PrintValue("Muon Fakes", m_fakes);
    if (verbosity_) PrintValue("Double Fakes", em_fakes);
    auto norm = ValueAdd(e_fakes, m_fakes);
    norm = ValueSubtract(norm, em_fakes);
    return norm;
  }


  HTTRun2Analysis::Value HTTRun2Analysis::GetRateViaWMethod(std::vector<std::string> const& w_sample,
                          std::string const& ratio_cat,
                          std::string const& ratio_control_sel,
                          std::string const& ratio_signal_sel,
                          std::vector<std::string> const& data_sample,
                          std::string const& cat,
                          std::string const& control_sel,
                          std::vector<std::string> const& sub_samples,
                          std::string const& wt,
                          std::map<std::string, std::function<Value()>> dict
                          ) {
    if (verbosity_) {
      std::cout << "[HTTRun2Analysis::GetRateViaWMethod]\n";
      //std::cout << "ExtrapFactor:   " << boost::format("%s,'%s'/'%s','%s','%s'\n") % w_sample % ratio_signal_sel 
       //         % ratio_control_sel % ratio_cat % wt;
//      std::cout << "Sideband:       " << boost::format("%s,'%s','%s','%s'\n") % data_sample % control_sel % cat % wt;
    }
    Value ratio = SampleRatio(w_sample, ratio_control_sel, ratio_cat, ratio_signal_sel, ratio_cat, wt); 
    Value data_control = GetRate(data_sample, control_sel, cat, wt);
 //   if (verbosity_) PrintValue(data_sample, data_control);
    Value total_bkg;
    for (unsigned i = 0; i < sub_samples.size(); ++i) {
      Value bkr;
      if (dict.count(sub_samples[i])) {
        bkr = ((*dict.find(sub_samples[i])).second)(); // find and evaluate function
      } else {
        bkr = GetLumiScaledRate(sub_samples[i], control_sel, cat, wt);
      }
      if (verbosity_) PrintValue("-"+sub_samples[i], bkr);
      double new_err = std::sqrt((total_bkg.second * total_bkg.second) + (bkr.second * bkr.second));
      total_bkg.first += bkr.first;
      total_bkg.second = new_err;
    }
    if (verbosity_) PrintValue("TotalBkg", total_bkg);
    double w_control_err = std::sqrt((total_bkg.second * total_bkg.second) + (data_control.second * data_control.second));
    Value w_control(data_control.first - total_bkg.first, w_control_err);
    if (verbosity_) PrintValue("WSideband", w_control);
    if (verbosity_) PrintValue("ExtrapFactor", ratio);
    Value w_signal = ValueProduct(w_control, ratio);
    return w_signal;
  }


  HTTRun2Analysis::Value HTTRun2Analysis::GetRateViaWOSSSMethod(std::vector<std::string> const& w_sample,
                          std::string const& ratio_cat,
                          std::string const& ratio_control_sel,
                          std::string const& ratio_signal_sel,
                          std::vector<std::string> const& data_sample,
                          std::string const& cat,
                          std::string const& control_sel,
                          std::vector<std::string> const& sub_samples,
                          bool get_os,
                          std::string const& wt,
                          std::map<std::string, std::function<Value()>> dict
                          ) {
    if (verbosity_) {
      std::cout << "[HTTRun2Analysis::GetRateViaWOSSSMethod]\n";
//      std::cout << "ExtrapFactor:   " << boost::format("%s,'%s'/'%s','%s','%s'\n") % w_sample % ratio_signal_sel 
 //               % ratio_control_sel % ratio_cat % wt;
//      std::cout << "Sideband:       " << boost::format("%s,'%s','%s','%s'\n") % data_sample % control_sel % cat % wt;
    }
    std::string os_ctr_sel = "os && "+control_sel;
    std::string ss_ctr_sel = "!os && "+control_sel;
    Value w_os_ss_ratio = SampleRatio(w_sample, "!os", ratio_cat, "os",ratio_cat,wt);
    Value ratio = SampleRatio(w_sample, ratio_control_sel, ratio_cat, ratio_signal_sel, ratio_cat, wt);
    Value data_control_os = GetRate(data_sample, os_ctr_sel, cat, wt);
    Value data_control_ss = GetRate(data_sample, ss_ctr_sel, cat, wt);
    if(verbosity_) PrintValue("DataControlOS", data_control_os);
    if(verbosity_) PrintValue("DataControlSS", data_control_ss);
 //   if (verbosity_) PrintValue(data_sample, data_control);
    Value total_bkg_os;
    Value total_bkg_ss;
    for (unsigned i = 0; i < sub_samples.size(); ++i) {
      Value bkr_os;
      Value bkr_ss;
      if (dict.count(sub_samples[i])) {
        bkr_os = ((*dict.find(sub_samples[i])).second)(); // find and evaluate function
        bkr_ss = ((*dict.find(sub_samples[i])).second)(); // find and evaluate function
      } else {
        bkr_os = GetLumiScaledRate(sub_samples[i], os_ctr_sel, cat, wt);
        bkr_ss = GetLumiScaledRate(sub_samples[i], ss_ctr_sel, cat, wt);
      }
      if (verbosity_) PrintValue("-os"+sub_samples[i], bkr_os);
      if (verbosity_) PrintValue("-ss"+sub_samples[i], bkr_ss);
      double new_err_os = std::sqrt((total_bkg_os.second * total_bkg_os.second) + (bkr_os.second * bkr_os.second));
      double new_err_ss = std::sqrt((total_bkg_ss.second * total_bkg_ss.second) + (bkr_ss.second * bkr_ss.second));
      total_bkg_os.first += bkr_os.first;
      total_bkg_ss.first += bkr_ss.first;
      total_bkg_os.second = new_err_os;
      total_bkg_ss.second = new_err_ss;
    }
    if (verbosity_) PrintValue("TotalBkgOS", total_bkg_os);
    if (verbosity_) PrintValue("TotalBkgSS", total_bkg_ss);
//    double w_control_err = std::sqrt((total_bkg_os.second * total_bkg_os.second) + (data_control_os.second * data_control_os.second) + (data_control_ss.second *data_control_ss.second) + (total_bkg_ss.second * total_bkg_ss.second));
    double YOSterm = std::sqrt((total_bkg_os.second * total_bkg_os.second)+(data_control_os.second*data_control_os.second));
    double YSSterm = std::sqrt((total_bkg_os.second * total_bkg_os.second)+(data_control_os.second*data_control_os.second))*qcd_os_ss_factor_;
    double RQterm = 0.1*qcd_os_ss_factor_*(std::sqrt((total_bkg_ss.second* total_bkg_ss.second)+(data_control_ss.second*data_control_ss.second))*w_os_ss_ratio.first-std::sqrt((total_bkg_os.second * total_bkg_os.second)+(data_control_os.second*data_control_os.second)))/(w_os_ss_ratio.first-qcd_os_ss_factor_);
    double RWterm = 0.08*w_os_ss_ratio.first*(data_control_ss.first - total_bkg_ss.first);
    if(get_os){
      YOSterm *= w_os_ss_ratio.first;
      YSSterm *= w_os_ss_ratio.first;
      RQterm *= w_os_ss_ratio.first;
      RWterm = 0.08*w_os_ss_ratio.first*(qcd_os_ss_factor_*qcd_os_ss_factor_*std::sqrt((total_bkg_ss.second* total_bkg_ss.second)+(data_control_ss.second*data_control_ss.second))-qcd_os_ss_factor_*std::sqrt((total_bkg_os.second * total_bkg_os.second)+(data_control_os.second*data_control_os.second)))/(w_os_ss_ratio.first-qcd_os_ss_factor_);
    }
    double w_control_err = (std::sqrt(YOSterm*YOSterm + YSSterm*YSSterm + RQterm*RQterm + RWterm*RWterm))/(w_os_ss_ratio.first - qcd_os_ss_factor_);
    double w_control_first = ((data_control_os.first - total_bkg_os.first) - (data_control_ss.first - total_bkg_ss.first)*qcd_os_ss_factor_)/(w_os_ss_ratio.first-qcd_os_ss_factor_);
    

    if(get_os) w_control_first*=w_os_ss_ratio.first;
    Value w_control(w_control_first, w_control_err);
    if (verbosity_) PrintValue("WSideband", w_control);
    if (verbosity_) PrintValue("ExtrapFactor", ratio);
    if (verbosity_) PrintValue("W OS/SS Ratio", w_os_ss_ratio);
    Value w_signal = ValueProduct(w_control, ratio);
    return w_signal;
  }


  HTTRun2Analysis::Value HTTRun2Analysis::GetRateViaTauTauQCDMethod(std::vector<std::string> const& data_sample,
                          std::string const& ratio_control_cat,
                          std::string const& ratio_signal_cat,
                          std::string const& ratio_sel,
                          std::string const& control_sel,
                          std::string const& cat,
                          std::vector<std::string> const& sub_samples,
                          std::string const& wt,
                          std::map<std::string, std::function<Value()>> dict
                          ) {
    if (verbosity_) {
      std::cout << "[HTTRun2Analysis::GetRateViaTauTauQCDMethod]\n";
      std::cout << "ExtrapFactor:   " << boost::format("'%s'/'%s','%s','%s'\n") % ratio_signal_cat 
                % ratio_control_cat % ratio_sel % wt;
      std::cout << "ExtrapFactorOS:   " << boost::format("'%s'/'%s','%s','%s'\n") % ratio_signal_cat
                % cat % control_sel % wt;
      std::cout << "Sideband:       " << boost::format("'%s','%s','%s'\n") %  control_sel % cat % wt;
    }

    Value ratio = DataSampleRatio(data_sample, ratio_sel, ratio_control_cat, ratio_sel, ratio_signal_cat, wt, sub_samples, dict);
    Value ratio_os = DataSampleRatio(data_sample, control_sel, cat, control_sel, ratio_signal_cat, wt, sub_samples, dict);
    Value data_control = GetRate(data_sample, control_sel, cat, wt);
 //   if (verbosity_) PrintValue(data_sample, data_control);
    Value total_bkg;
    for (unsigned i = 0; i < sub_samples.size(); ++i) {
      Value bkr;
      if (dict.count(sub_samples[i])) {
        bkr = ((*dict.find(sub_samples[i])).second)(); // find and evaluate function
      } else {
        bkr = GetLumiScaledRate(sub_samples[i], control_sel, cat, wt);
      }
      if (verbosity_) PrintValue("-"+sub_samples[i], bkr);
      double new_err = std::sqrt((total_bkg.second * total_bkg.second) + (bkr.second * bkr.second));
      total_bkg.first += bkr.first;
      total_bkg.second = new_err;
    }
    if (verbosity_) PrintValue("TotalBkg", total_bkg);
    double qcd_control_err = std::sqrt((total_bkg.second * total_bkg.second) + (data_control.second * data_control.second));
    Value qcd_control(data_control.first - total_bkg.first, qcd_control_err);
    if (verbosity_) PrintValue("QCDSideband", qcd_control);
    if (verbosity_) PrintValue("ExtrapFactor", ratio);
    if (verbosity_) PrintValue("ExtrapFactorOS", ratio_os);
    Value qcd_signal = ValueProduct(qcd_control, ratio);
    return qcd_signal;
  }

  HTTRun2Analysis::Value HTTRun2Analysis::GetRateViaQCDMethod(HTTRun2Analysis::Value const& ratio,
                          std::vector<std::string> const& data_sample,
                          std::string const& control_selection,
                          std::string const& category,
                          std::vector<std::string> const& sub_samples,
                          std::string const& weight,
                          std::map<std::string, std::function<Value()>> dict
                          ) {
    if (verbosity_) {
      std::cout << "[HTTRun2Analysis::GetRateViaQCDMethod]\n";
//      std::cout << "Sideband:       " << boost::format("%s,'%s','%s','%s'\n") % data_sample % control_selection % category % weight;
    }
    Value data_control = GetRate(data_sample, control_selection, category, weight);
//    if (verbosity_) PrintValue(data_sample, data_control);
    Value total_bkg;
    for (unsigned i = 0; i < sub_samples.size(); ++i) {
      Value bkr;
      if (dict.count(sub_samples[i])) {
        bkr = ((*dict.find(sub_samples[i])).second)(); // find and evaluate function
      } else {
        bkr = GetLumiScaledRate(sub_samples[i], control_selection, category, weight);
      }
      if (verbosity_) PrintValue("-"+sub_samples[i], bkr);
      double new_err = std::sqrt((total_bkg.second * total_bkg.second) + (bkr.second * bkr.second));
      total_bkg.first += bkr.first;
      total_bkg.second = new_err;
    }
    if (verbosity_) PrintValue("TotalBkg", total_bkg);
    double qcd_control_err = std::sqrt((total_bkg.second * total_bkg.second) + (data_control.second * data_control.second));
    Value qcd_control(data_control.first - total_bkg.first, qcd_control_err);
    if (verbosity_) PrintValue("QCDSideband", qcd_control);
    Value qcd_signal = ValueProduct(qcd_control, ratio);
    if (verbosity_) PrintValue("OS/SS Factor", ratio);
    return qcd_signal;
  }

  TH1F HTTRun2Analysis::GetShapeViaQCDMethod(std::string const& variable,
                          std::vector<std::string> const& data_sample,
                          std::string const& selection,
                          std::string const& category,
                          std::vector<std::string> const& sub_samples,
                          std::string const& weight,
                          std::map<std::string, std::function<Value()>> dict
                          ) {
    if (verbosity_) {
      std::cout << "[HTTRun2Analysis::GetShapeViaQCDMethod]\n";
//      std::cout << "Sideband:       " << boost::format("%s,'%s','%s','%s'\n") % data_sample % selection % category % weight;

    }
    TH1F result = GetLumiScaledShape(variable, data_sample, selection, category, weight);
    for (unsigned i = 0; i < sub_samples.size(); ++i) {
      if (dict.count(sub_samples[i])) {
        Value bkr_rate = ((*dict.find(sub_samples[i])).second)(); // find and evaluate function
        TH1F tmp = GetShape(variable, sub_samples.at(i), selection, category, weight);
        SetNorm(&tmp, bkr_rate.first);
        result.Add(&tmp, -1.);
      } else {
      TH1F tmp = GetLumiScaledShape(variable, sub_samples[i], selection, category, weight);
      result.Add(&tmp, -1.);
     }
    }
    return result;
  }


  TH1F HTTRun2Analysis::GetShapeViaFakesMethod(std::string const& var,
                              std::string const& sel,
                              std::string const& cat,
                              std::string const& wt) {
    TH1F fr_e_hist =  this->GetShape(var, "Special_20_Data", sel, cat, wt);
    TH1F fr_m_hist =  this->GetShape(var, "Special_21_Data", sel, cat, wt);
    TH1F fr_em_hist = this->GetShape(var, "Special_22_Data", sel, cat, wt);
    fr_e_hist.Add(&fr_m_hist,1.0);
    fr_e_hist.Add(&fr_em_hist,-1.0);
    return fr_e_hist;
  }

  HTTRun2Analysis::Value HTTRun2Analysis::ValueProduct(Value const& p1, Value const& p2) {
    if (p1.first == 0.0 || p2.first == 0.0) {
      //if (verbosity_ > 0) std::cout << "[HTTRun2Analysis::ValueProduct] At least one value is zero, returning 0.0 +/- 0.0" << std::endl;
      return std::make_pair(0.0, 0.0);
    }
    double f = p1.first * p2.first;
    double a_sqrd = std::pow(p1.second / p1.first, 2.0);
    double b_sqrd = std::pow(p2.second / p2.first, 2.0);
    double f_err = f * sqrt( a_sqrd + b_sqrd );
    return std::make_pair(f, f_err);
  }
  HTTRun2Analysis::Value HTTRun2Analysis::ValueDivide(Value const& p1, Value const& p2) {
    if (p1.first == 0.0 && p2.first == 0.0) {
      std::cout << "[HTTRun2Analysis::ValueDivide] Numerator and denominator both zero, returning 0.0 +/- 0.0" << std::endl;
      return std::make_pair(0.0, 0.0);
    }
    if (p1.first == 0.0) {
      //if (verbosity_ > 0) std::cout << "[HTTRun2Analysis::ValueDivide] Numerator is zero, returning 0.0 +/- 0.0" << std::endl;
      return std::make_pair(0.0, 0.0);
    }
    if (p2.first == 0.0) {
      std::cout << "[HTTRun2Analysis::ValueDivide] Denominator is zero, returning 0.0 +/- 0.0" << std::endl;
      return std::make_pair(0.0, 0.0);
    }
    double f = p1.first / p2.first;
    double a_sqrd = std::pow(p1.second / p1.first, 2.0);
    double b_sqrd = std::pow(p2.second / p2.first, 2.0);
    double f_err = f * sqrt( a_sqrd + b_sqrd );
    return std::make_pair(f, f_err);
  }
  HTTRun2Analysis::Value HTTRun2Analysis::ValueAdd(Value const& p1, Value const& p2) {
    double f = p1.first + p2.first;
    double f_err = sqrt( p1.second*p1.second + p2.second*p2.second );
    return std::make_pair(f, f_err);
  }
  HTTRun2Analysis::Value HTTRun2Analysis::ValueSubtract(Value const& p1, Value const& p2) {
    double f = p1.first - p2.first;
    double f_err = sqrt( p1.second*p1.second + p2.second*p2.second );
    return std::make_pair(f, f_err);
  }

  void HTTRun2Analysis::PrintValue(std::string const& label, Value const& val) {
    //std::cout << (boost::format("%-45s %10.3f +/-   %10.3f  (%.4f)") % (label+":") % val.first % val.second % (val.second/val.first)) << std::endl;
    std::cout << (boost::format("%-45s %10.5f +/-   %10.5f  (%.4f)") % (label+":") % val.first % val.second % (val.second/val.first)) << std::endl;
  }

  std::string HTTRun2Analysis::ResolveAlias(std::string const& al) {
    if (alias_map_.count(al)) {
      return alias_map_.find(al)->second;
    } else {
      std::cerr << "Warning in <HTTRun2Analysis::ResolveAlias>: No alias for " << al << " found" << std::endl;
      return al;
    }
  }

  std::vector<std::string> HTTRun2Analysis::ResolveSamplesAlias(std::string const& al) {
    if (samples_alias_map_.count(al)) {
      return samples_alias_map_.find(al)->second;
    } else {
      std::cerr << "Warning in <HTTRun2Analysis::ResolveSamplesAlias>: No alias for " << al << " found" << std::endl;
      return std::vector<std::string>();
    }
  }

  void HTTRun2Analysis::SetAlias(std::string const& al, std::string const& val) {
    alias_map_[al] = val;
  }

  double HTTRun2Analysis::KolmogorovTest(std::string const& variable,
      std::string const& sample1, std::string const& selection1, std::string const& category1,
      std::string const& sample2, std::string const& selection2, std::string const& category2,
      std::string const& weight) {
    std::cout << "[HTTRun2Analysis::KolmogorovTest] Calculating statistic for shapes:" << std::endl;
    std::cout << "[1] " << boost::format("%s,'%s','%s','%s'\n") % sample1 % selection1 % category1 % weight;
    std::cout << "[2] " << boost::format("%s,'%s','%s','%s'\n") % sample2 % selection2 % category2 % weight;
    std::string full1 = BuildCutString(selection1, category1, weight);
    std::string full2 = BuildCutString(selection2, category2, weight);
    TH1::AddDirectory(true);
    ttrees_[sample1]->Draw(">>elist", full1.c_str(), "entrylist");
    TEntryList *elist1 = (TEntryList*)gDirectory->Get("elist");
    TH1::AddDirectory(false);
    unsigned entries1 = elist1->GetN();
    double x1;
    double wt1;
    std::vector<std::pair<double,double>> a(entries1,std::make_pair(0.,0.));
    ttrees_[sample1]->SetBranchAddress(variable.c_str(), &x1);
    ttrees_[sample1]->SetBranchAddress(weight.c_str(), &wt1);
    for (unsigned i = 0; i < entries1; ++i) {
      unsigned entry = elist1->GetEntry(i);
      ttrees_[sample1]->GetEntry(entry);
      a[i].first = x1;
      a[i].second = wt1;
    }
    gDirectory->Delete("elist;*");
    TH1::AddDirectory(true);
    ttrees_[sample2]->Draw(">>elist", full2.c_str(), "entrylist");
    TEntryList *elist2 = (TEntryList*)gDirectory->Get("elist");
    TH1::AddDirectory(false);
    unsigned entries2 = elist2->GetN();
    double x2;
    double wt2;
    std::vector<std::pair<double,double>> b(entries2,std::make_pair(0.,0.));
    ttrees_[sample2]->SetBranchAddress(variable.c_str(), &x2);
    ttrees_[sample2]->SetBranchAddress(weight.c_str(), &wt2);
    for (unsigned i = 0; i < entries2; ++i) {
      unsigned entry = elist2->GetEntry(i);
      ttrees_[sample2]->GetEntry(entry);
      b[i].first = x2;
      b[i].second = wt2;
    }
    gDirectory->Delete("elist;*");
    std::sort(a.begin(), a.end(), [](const std::pair<double, double>& first, const std::pair<double, double>& second)
      {
        return first.first < second.first;
      });
    std::sort(b.begin(), b.end(), [](const std::pair<double, double>& first, const std::pair<double, double>& second)
      {
        return first.first < second.first;
      });
    double prob = -1;
    int na = a.size();
    int nb = b.size();
    double rna = 0.;
    double rnb = 0.;
    for (unsigned i = 0; i < a.size(); ++i) rna += a[i].second;
    for (unsigned i = 0; i < b.size(); ++i) rnb += b[i].second;

    /* 
      The implementation below shamelessly stolen from:
      http://root.cern.ch/root/html/TMath.html#TMath:KolmogorovTest
      and modified to take into account weighted events.  Quite possible
      that this is not statistically valid, but seems to behave well enough
      when weights are close to unity.
    */
    // double rna = a.size();
    // double rnb = b.size();
    std::cout << "Entries: (1) " << rna << " (2) " << rnb << std::endl;
    // double sa  = 1./rna;
    // double sb  = 1./rnb;
    double rdiff = 0;
    double rdmax = 0;
    int ia = 0;
    int ib = 0;
    bool ok = false;
    for (int i=0;i<na+nb;i++) {
      if (a[ia].first < b[ib].first) {
         rdiff -= a[ia].second/rna;
         ia++;
         if (ia >= na) {ok = true; break;}
      } else if (a[ia].first > b[ib].first) {
         rdiff += b[ib].second/rnb;
         ib++;
         if (ib >= nb) {ok = true; break;}
      } else {
         // special cases for the ties
         double x = a[ia].first;
         while(a[ia].first == x && ia < na) {
            rdiff -= a[ia].second/rna;
            ia++;
         }
         while(b[ib].first == x && ib < nb) {
            rdiff += b[ib].second/rnb;
            ib++;
         }
         if (ia >= na) {ok = true; break;}
         if (ib >= nb) {ok = true; break;}
      }
      rdmax = TMath::Max(rdmax,TMath::Abs(rdiff));
    }
    if (ok) {
      rdmax = TMath::Max(rdmax,TMath::Abs(rdiff));
      Double_t z = rdmax * TMath::Sqrt(rna*rnb/(rna+rnb));
      prob = TMath::KolmogorovProb(z);
    }
    std::cout << " Kolmogorov Probability = " << prob << ", rmax=" << rdmax << std::endl;
    return prob;
  }



}

