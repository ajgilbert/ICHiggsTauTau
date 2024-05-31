#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/interface/HTTCategories.h"
#include "UserCode/ICHiggsTauTau/interface/PFCandidate.hh"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "Utilities/interface/FnRootTools.h"
#include "UserCode/ICHiggsTauTau/interface/L1TObject.hh"

#include "TMVA/Reader.h"
#include "TVector3.h"
#include "boost/format.hpp"
#include "TMath.h"
#include "TLorentzVector.h"


double CPWeight(double x, double sm, double ps, double mm ){
        x*=M_PI/180; //convert to radians
        return cos(x)*cos(x)*sm + sin(x)*sin(x)*ps + 2*cos(x)*sin(x)*(mm-sm/2-ps/2);

}

namespace ic {

  void HTTCategories::CorrectMETForShift(ic::Met * met, ROOT::Math::PxPyPzEVector const& shift) {
    double metx = met->vector().px() - shift.px();
    double mety = met->vector().py() - shift.py();
    double metet = sqrt(metx*metx + mety*mety);
    ROOT::Math::PxPyPzEVector new_met(metx, mety, 0, metet);
    met->set_vector(ROOT::Math::PtEtaPhiEVector(new_met));
  }

  HTTCategories::HTTCategories(std::string const& name) : ModuleBase(name), 
      channel_(channel::et), 
      era_(era::data_2012_rereco),
      strategy_(strategy::paper2013) {
      ditau_label_ = "emtauCandidates";
      jets_label_ = "pfJetsPFlow";
      met_label_ = "pfMVAMetNoLeptons";
      fs_ = NULL;
      write_tree_ = true;
      make_sync_ntuple_ = false;
      make_mva_ntuple_ = false;
      sync_output_name_ = "SYNC.root";
      is_embedded_=false;
      is_data_=false;
      systematic_shift_ = false;
      do_qcd_scale_wts_ = false;
      do_mssm_higgspt_ = false;
      do_faketaus_ = false;
      trg_applied_in_mc_ = true;
}

  HTTCategories::~HTTCategories() {
	;
  }

  int HTTCategories::PreAnalysis() {
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "HTTCategories" << std::endl;
    std::cout << "-------------------------------------" << std::endl;    
      std::cout << boost::format(param_fmt()) % "channel"         % Channel2String(channel_);
      std::cout << boost::format(param_fmt()) % "strategy"        % Strategy2String(strategy_);
      std::cout << boost::format(param_fmt()) % "era"             % Era2String(era_);
      std::cout << boost::format(param_fmt()) % "dilepton_label"  % ditau_label_;
      std::cout << boost::format(param_fmt()) % "met_label"       % met_label_;
      std::cout << boost::format(param_fmt()) % "jets_label"      % jets_label_;
      std::cout << boost::format(param_fmt()) % "write_tree"      % write_tree_;
      std::cout << boost::format(param_fmt()) % "make_sync_ntuple" % make_sync_ntuple_;
      std::cout << boost::format(param_fmt()) % "make_mva_ntuple" % make_mva_ntuple_;

    std::string MVADM2017;
    if(era_ == era::data_2016 || era_ == era::data_2017 || era_ == era::data_2018){ 
         MVADM2017 = "MVADM2017v1";

    } else if(era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP || era_ == era::data_2017UL || era_ == era::data_2018UL){ 
         MVADM2017 = "MVADM2017v2";
    } else {
         MVADM2017 = " ";
    }

 
    if (!is_data_) {
      // initialize IP corrector
      std::string cmsswBase = (getenv ("CMSSW_BASE"));
      TString ip_corr_filename = "";
      TString ip_corr_filename_elecs = "";
      if(is_embedded_) {
        if(era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) ip_corr_filename = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ip_embed_2016.root";
        if(era_ == era::data_2017 || era_ == era::data_2017UL) ip_corr_filename = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ip_embed_2017.root";
        else                       ip_corr_filename = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ip_embed_2018.root";

        if(era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) ip_corr_filename_elecs = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ipEle_embed_2016.root";
        if(era_ == era::data_2017 || era_ == era::data_2017UL) ip_corr_filename_elecs = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ipEle_embed_2017.root";
        else                       ip_corr_filename_elecs = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ipEle_embed_2018.root";
      } else {
        if(era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) ip_corr_filename = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ip_2016.root";
        if(era_ == era::data_2017 || era_ == era::data_2017UL) ip_corr_filename = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ip_2017.root";
        else                       ip_corr_filename = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ip_2018.root";

        if(era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) ip_corr_filename_elecs = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ipEle_2016.root";
        if(era_ == era::data_2017 || era_ == era::data_2017UL) ip_corr_filename_elecs = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ipEle_2017.root";
        else                       ip_corr_filename_elecs = TString(cmsswBase)+"/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/input/ip_corrs/ipEle_2018.root";
      }
      ipCorrector.Init(ip_corr_filename);
      ipCorrectorEle.Init(ip_corr_filename_elecs);
    }


    rand = new TRandom3(0);

    // for smearing tests
    //TFile f_smear("smear_file_2018_zmm_formula.root");
    //fns_["smear_func_ptlow30_etalow0p0"] = (TF1*)f_smear.Get(
    //  "smear_func_ptlow30_etalow0p0"
    //);
    //fns_["smear_func_ptlow30_etalow1p1"] = (TF1*)f_smear.Get(
    //  "smear_func_ptlow30_etalow1p1"
    //);
    //fns_["smear_func_ptlow30_etalow2p1"] = (TF1*)f_smear.Get(
    //  "smear_func_ptlow30_etalow2p1"
    //);
    //fns_["smear_func_ptlow30_etalow3p0"] = (TF1*)f_smear.Get(
    //  "smear_func_ptlow30_etalow3p0"
    //);
    //fns_["smear_func_ptlow60_etalow0p0"] = (TF1*)f_smear.Get(
    //  "smear_func_ptlow60_etalow0p0"
    //);
    //fns_["smear_func_ptlow60_etalow1p1"] = (TF1*)f_smear.Get(
    //  "smear_func_ptlow60_etalow1p1"
    //);
    //fns_["smear_func_ptlow60_etalow2p1"] = (TF1*)f_smear.Get(
    //  "smear_func_ptlow60_etalow2p1"
    //);
    //fns_["smear_func_ptlow60_etalow3p0"] = (TF1*)f_smear.Get(
    //  "smear_func_ptlow60_etalow3p0"
    //);
    //f_smear.Close();
    if (fs_ && write_tree_) {
      outtree_ = fs_->make<TTree>("ntuple","ntuple");
      
      outtree_->Branch("event",             &event_);
      outtree_->Branch("run",               &run_);
      outtree_->Branch("lumi",               &lumi_);
      outtree_->Branch("npu",               &n_pu_, "n_pu/F");
      outtree_->Branch("rho",               &rho_, "rho/F");
      outtree_->Branch("puweight",          &pu_weight_, "pu_weight/F");
      outtree_->Branch("wt",                &wt_.var_double);
      outtree_->Branch("wt_dysoup",         &wt_dysoup_);
      outtree_->Branch("wt_dysoup_NLO",         &wt_dysoup_NLO_);
      outtree_->Branch("wt_wsoup",         &wt_wsoup_);
      // btag weights (for event reweighting method)
      outtree_->Branch("wt_btag",           &wt_btag_);
      if(!systematic_shift_) {
        outtree_->Branch("wt_btag_up",        &wt_btag_up_);
        outtree_->Branch("wt_btag_down",      &wt_btag_down_);
        outtree_->Branch("wt_btag_real_down", &wt_btag_real_down_);
        outtree_->Branch("wt_btag_real_up",   &wt_btag_real_up_);
        outtree_->Branch("wt_btag_fake_down", &wt_btag_fake_down_);
        outtree_->Branch("wt_btag_fake_up",   &wt_btag_fake_up_);
      }
      //
      outtree_->Branch("trigweight_1", &trigweight_1_, "trigweight_1/F");
      outtree_->Branch("trigweight_2", &trigweight_2_, "trigweight_2/F");
      outtree_->Branch("xtrg_sf", &xtrg_sf_);
      outtree_->Branch("single_l_sf", &single_l_sf_);
      outtree_->Branch("OR_notrig", OR_notrig_);
      outtree_->Branch("xtrg_notrig", xtrg_notrig_);
      outtree_->Branch("idisoweight_1", &idisoweight_1_, "idisoweight_1/F");
      outtree_->Branch("idisoweight_2", &idisoweight_2_, "idisoweight_2/F");
      outtree_->Branch("trackingweight_1", &trackingweight_1_, "trackingweight_1/F");
      outtree_->Branch("trackingweight_2", &trackingweight_2_, "trackingweight_1/F");
      outtree_->Branch("wt_emb_sel_kit", &wt_emb_sel_kit_, "wt_emb_sel_kit/F");
      outtree_->Branch("wt_emb_sel_trg_kit", &wt_emb_sel_trg_kit_, "wt_emb_sel_trg_kit/F");
      outtree_->Branch("effweight", &effweight_, "effweight/F");
      outtree_->Branch("wt_quarkmass", &wt_quarkmass_);
      outtree_->Branch("wt_fullquarkmass", & wt_fullquarkmass_);
      outtree_->Branch("wt_tau_id_extra", &wt_tau_id_extra_);
      outtree_->Branch("tauFlag_1", &tauFlag_1_);
      outtree_->Branch("tauFlag_2", &tauFlag_2_);
      outtree_->Branch("tau_leg_SF", &tau_leg_SF_);

      outtree_->Branch("ip_mag_1", &ip_mag_1_);
      outtree_->Branch("ip_mag_2", &ip_mag_2_);
      outtree_->Branch("ip_sig_1", &ip_sig_1_);
      outtree_->Branch("ip_sig_2", &ip_sig_2_);

      if(!systematic_shift_) {
        outtree_->Branch("ip_sig_1_raw",&ip_sig_1_raw_);
        outtree_->Branch("ip_sig_2_raw",&ip_sig_2_raw_);
        outtree_->Branch("ip_sig_1_up",&ip_sig_1_up_);
        outtree_->Branch("ip_sig_2_up",&ip_sig_2_up_);
        outtree_->Branch("ip_sig_1_down",&ip_sig_1_down_);
        outtree_->Branch("ip_sig_2_down",&ip_sig_2_down_);
      }

      outtree_->Branch("q_tot_1", &q_tot_1_);
      outtree_->Branch("q_tot_2", &q_tot_2_);
      outtree_->Branch("wt_tau_id_mvadm", & wt_tau_id_mvadm_);
      outtree_->Branch("wt_tau_trg_ic",    &wt_tau_trg_ic_);
      outtree_->Branch("wt_tau_trg_mvadm",    &wt_tau_trg_mvadm_);
      outtree_->Branch("wt_tau_trg_mssm_doubleonly",    &wt_tau_trg_mssm_doubleonly_);
      outtree_->Branch("wt_tau_trg_mssm",    &wt_tau_trg_mssm_);
      outtree_->Branch("wt_tau_trg_mssm_dm0_up",    &wt_tau_trg_mssm_dm0_up_);
      outtree_->Branch("wt_tau_trg_mssm_dm0_down",    &wt_tau_trg_mssm_dm0_down_);
      outtree_->Branch("wt_tau_trg_mssm_dm1_up",    &wt_tau_trg_mssm_dm1_up_);
      outtree_->Branch("wt_tau_trg_mssm_dm1_down",    &wt_tau_trg_mssm_dm1_down_);
      outtree_->Branch("wt_tau_trg_mssm_dm10_up",    &wt_tau_trg_mssm_dm10_up_);
      outtree_->Branch("wt_tau_trg_mssm_dm10_down",    &wt_tau_trg_mssm_dm10_down_);
      outtree_->Branch("wt_tau_trg_mssm_dm11_up",    &wt_tau_trg_mssm_dm11_up_);
      outtree_->Branch("wt_tau_trg_mssm_dm11_down",    &wt_tau_trg_mssm_dm11_down_);
      outtree_->Branch("wt_tau_trg_mssm_highpt_dm0_up",    &wt_tau_trg_mssm_highpt_dm0_up_);
      outtree_->Branch("wt_tau_trg_mssm_highpt_dm0_down",    &wt_tau_trg_mssm_highpt_dm0_down_);
      outtree_->Branch("wt_tau_trg_mssm_highpt_dm1_up",    &wt_tau_trg_mssm_highpt_dm1_up_);
      outtree_->Branch("wt_tau_trg_mssm_highpt_dm1_down",    &wt_tau_trg_mssm_highpt_dm1_down_);
      outtree_->Branch("wt_tau_trg_mssm_highpt_dm10_up",    &wt_tau_trg_mssm_highpt_dm10_up_);
      outtree_->Branch("wt_tau_trg_mssm_highpt_dm10_down",    &wt_tau_trg_mssm_highpt_dm10_down_);
      outtree_->Branch("wt_tau_trg_mssm_highpt_dm11_up",    &wt_tau_trg_mssm_highpt_dm11_up_);
      outtree_->Branch("wt_tau_trg_mssm_highpt_dm11_down",    &wt_tau_trg_mssm_highpt_dm11_down_);
      outtree_->Branch("wt_tau_trg_mssm_singletau_up",    &wt_tau_trg_mssm_singletau_up_);
      outtree_->Branch("wt_tau_trg_mssm_singletau_down",    &wt_tau_trg_mssm_singletau_down_);
      outtree_->Branch("wt_tau_trg_mssm_singlelep_up",    &wt_tau_trg_mssm_singlelep_up_);
      outtree_->Branch("wt_tau_trg_mssm_singlelep_down",    &wt_tau_trg_mssm_singlelep_down_);
      outtree_->Branch("wt_tau_trg_mssm_crosslep_up",    &wt_tau_trg_mssm_crosslep_up_);
      outtree_->Branch("wt_tau_trg_mssm_crosslep_down",    &wt_tau_trg_mssm_crosslep_down_);

      outtree_->Branch("wt_tau_id_mssm",    &wt_tau_id_mssm_);
      outtree_->Branch("wt_tau_id_mssm_bin5_up",    &wt_tau_id_mssm_bin5_up_);
      outtree_->Branch("wt_tau_id_mssm_bin5_down",    &wt_tau_id_mssm_bin5_down_);
      outtree_->Branch("wt_tau_id_mssm_bin6_up",    &wt_tau_id_mssm_bin6_up_);
      outtree_->Branch("wt_tau_id_mssm_bin6_down",    &wt_tau_id_mssm_bin6_down_);

      outtree_->Branch("wt_tau_id_dm_vs_pt",    &wt_tau_sf_id_dm_vs_pt_);
      outtree_->Branch("wt_tau_id_dm",    &wt_tau_sf_id_dm_);
      outtree_->Branch("wt_tau_id_pt",    &wt_tau_sf_id_pt_);

      outtree_->Branch("tau_idisoweight_1", &tau_idisoweight_1_);
      outtree_->Branch("tau_idisoweight_2", &tau_idisoweight_2_);

      outtree_->Branch("tau_idisoweight_ratio_1_up", &tau_idisoweight_ratio_1_up_);
      outtree_->Branch("tau_idisoweight_ratio_2_up", &tau_idisoweight_ratio_2_up_);

      outtree_->Branch("tau_idisoweight_ratio_1_down", &tau_idisoweight_ratio_1_down_);
      outtree_->Branch("tau_idisoweight_ratio_2_down", &tau_idisoweight_ratio_2_down_);

      outtree_->Branch("tau_idisoweight_vvvloose_1", &tau_idisoweight_vvvloose_1_);
      outtree_->Branch("tau_idisoweight_vvvloose_2", &tau_idisoweight_vvvloose_2_);

      outtree_->Branch("tau_idisoweight_vvloose_1", &tau_idisoweight_vvloose_1_);
      outtree_->Branch("tau_idisoweight_vvloose_2", &tau_idisoweight_vvloose_2_);

      outtree_->Branch("tau_idisoweight_vloose_1", &tau_idisoweight_vloose_1_);
      outtree_->Branch("tau_idisoweight_vloose_2", &tau_idisoweight_vloose_2_);

      outtree_->Branch("tau_idisoweight_loose_1", &tau_idisoweight_loose_1_);
      outtree_->Branch("tau_idisoweight_loose_2", &tau_idisoweight_loose_2_);

      outtree_->Branch("tau_idisoweight_medium_1", &tau_idisoweight_medium_1_);
      outtree_->Branch("tau_idisoweight_medium_2", &tau_idisoweight_medium_2_);

      outtree_->Branch("tau_idisoweight_tight_1", &tau_idisoweight_tight_1_);
      outtree_->Branch("tau_idisoweight_tight_2", &tau_idisoweight_tight_2_);

      outtree_->Branch("tau_idisoweight_ratio_1_uncert0_up", &tau_idisoweight_ratio_1_uncert0_up_);
      outtree_->Branch("tau_idisoweight_ratio_2_uncert0_up", &tau_idisoweight_ratio_2_uncert0_up_);

      outtree_->Branch("tau_idisoweight_ratio_1_uncert0_down", &tau_idisoweight_ratio_1_uncert0_down_);
      outtree_->Branch("tau_idisoweight_ratio_2_uncert0_down", &tau_idisoweight_ratio_2_uncert0_down_);

      outtree_->Branch("tau_idisoweight_ratio_1_uncert1_up", &tau_idisoweight_ratio_1_uncert1_up_);
      outtree_->Branch("tau_idisoweight_ratio_2_uncert1_up", &tau_idisoweight_ratio_2_uncert1_up_);

      outtree_->Branch("tau_idisoweight_ratio_1_uncert1_down", &tau_idisoweight_ratio_1_uncert1_down_);
      outtree_->Branch("tau_idisoweight_ratio_2_uncert1_down", &tau_idisoweight_ratio_2_uncert1_down_);

      outtree_->Branch("tau_idisoweight_ratio_1_syst_all_eras_up", &tau_idisoweight_ratio_1_syst_all_eras_up_);
      outtree_->Branch("tau_idisoweight_ratio_2_syst_all_eras_up", &tau_idisoweight_ratio_2_syst_all_eras_up_);

      outtree_->Branch("tau_idisoweight_ratio_1_syst_all_eras_down", &tau_idisoweight_ratio_1_syst_all_eras_down_);
      outtree_->Branch("tau_idisoweight_ratio_2_syst_all_eras_down", &tau_idisoweight_ratio_2_syst_all_eras_down_);

      outtree_->Branch("tau_idisoweight_ratio_1_syst_year_up", &tau_idisoweight_ratio_1_syst_year_up_);
      outtree_->Branch("tau_idisoweight_ratio_2_syst_year_up", &tau_idisoweight_ratio_2_syst_year_up_);

      outtree_->Branch("tau_idisoweight_ratio_1_syst_year_down", &tau_idisoweight_ratio_1_syst_year_down_);
      outtree_->Branch("tau_idisoweight_ratio_2_syst_year_down", &tau_idisoweight_ratio_2_syst_year_down_);

      outtree_->Branch("tau_idisoweight_ratio_1_syst_dm_year_up", &tau_idisoweight_ratio_1_syst_dm_year_up_);
      outtree_->Branch("tau_idisoweight_ratio_2_syst_dm_year_up", &tau_idisoweight_ratio_2_syst_dm_year_up_);

      outtree_->Branch("tau_idisoweight_ratio_1_syst_dm_year_down", &tau_idisoweight_ratio_1_syst_dm_year_down_);
      outtree_->Branch("tau_idisoweight_ratio_2_syst_dm_year_down", &tau_idisoweight_ratio_2_syst_dm_year_down_);

      outtree_->Branch("tau_idisoweight_ratio_1_highpt_up", &tau_idisoweight_ratio_1_highpt_up_);
      outtree_->Branch("tau_idisoweight_ratio_2_highpt_up", &tau_idisoweight_ratio_2_highpt_up_);

      outtree_->Branch("tau_idisoweight_ratio_1_highpt_down", &tau_idisoweight_ratio_1_highpt_down_);
      outtree_->Branch("tau_idisoweight_ratio_2_highpt_down", &tau_idisoweight_ratio_2_highpt_down_);

      outtree_->Branch("wt_mg_nnlops", & wt_mg_nnlops_);
      outtree_->Branch("wt_ph_nnlops", & wt_ph_nnlops_);
      if(!systematic_shift_){
        outtree_->Branch("wt_tau_id_lowpt_mvadm0_up", &wt_tau_id_lowpt_mvadm0_up_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm1_up", &wt_tau_id_lowpt_mvadm1_up_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm2_up", &wt_tau_id_lowpt_mvadm2_up_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm10_up", &wt_tau_id_lowpt_mvadm10_up_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm11_up", &wt_tau_id_lowpt_mvadm11_up_);
        outtree_->Branch("wt_tau_id_highpt_mvadm0_up", &wt_tau_id_highpt_mvadm0_up_);
        outtree_->Branch("wt_tau_id_highpt_mvadm1_up", &wt_tau_id_highpt_mvadm1_up_);
        outtree_->Branch("wt_tau_id_highpt_mvadm2_up", &wt_tau_id_highpt_mvadm2_up_);
        outtree_->Branch("wt_tau_id_highpt_mvadm10_up", &wt_tau_id_highpt_mvadm10_up_);
        outtree_->Branch("wt_tau_id_highpt_mvadm11_up", &wt_tau_id_highpt_mvadm11_up_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm0_down", &wt_tau_id_lowpt_mvadm0_down_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm1_down", &wt_tau_id_lowpt_mvadm1_down_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm2_down", &wt_tau_id_lowpt_mvadm2_down_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm10_down", &wt_tau_id_lowpt_mvadm10_down_);
        outtree_->Branch("wt_tau_id_lowpt_mvadm11_down", &wt_tau_id_lowpt_mvadm11_down_);
        outtree_->Branch("wt_tau_id_highpt_mvadm0_down", &wt_tau_id_highpt_mvadm0_down_);
        outtree_->Branch("wt_tau_id_highpt_mvadm1_down", &wt_tau_id_highpt_mvadm1_down_);
        outtree_->Branch("wt_tau_id_highpt_mvadm2_down", &wt_tau_id_highpt_mvadm2_down_);
        outtree_->Branch("wt_tau_id_highpt_mvadm10_down", &wt_tau_id_highpt_mvadm10_down_);
        outtree_->Branch("wt_tau_id_highpt_mvadm11_down", &wt_tau_id_highpt_mvadm11_down_);

        outtree_->Branch("wt_tau_id_dm0_up", &wt_tau_id_dm0_up_);
        outtree_->Branch("wt_tau_id_dm1_up", &wt_tau_id_dm1_up_);
        outtree_->Branch("wt_tau_id_dm10_up", &wt_tau_id_dm10_up_);
        outtree_->Branch("wt_tau_id_dm11_up"  , &wt_tau_id_dm11_up_);  
        outtree_->Branch("wt_tau_id_dm0_down", &wt_tau_id_dm0_down_);
        outtree_->Branch("wt_tau_id_dm1_down", &wt_tau_id_dm1_down_);
        outtree_->Branch("wt_tau_id_dm10_down", &wt_tau_id_dm10_down_);
        outtree_->Branch("wt_tau_id_dm11_down", &wt_tau_id_dm11_down_);
        outtree_->Branch("wt_tau_id_pt_bin1_up", &wt_tau_id_pt_bin1_up_);
        outtree_->Branch("wt_tau_id_pt_bin2_up", &wt_tau_id_pt_bin2_up_);
        outtree_->Branch("wt_tau_id_pt_bin3_up", &wt_tau_id_pt_bin3_up_);
        outtree_->Branch("wt_tau_id_pt_bin4_up", &wt_tau_id_pt_bin4_up_);
        outtree_->Branch("wt_tau_id_pt_bin5_up", &wt_tau_id_pt_bin5_up_);
        outtree_->Branch("wt_tau_id_pt_bin1_down", &wt_tau_id_pt_bin1_down_);
        outtree_->Branch("wt_tau_id_pt_bin2_down", &wt_tau_id_pt_bin2_down_);
        outtree_->Branch("wt_tau_id_pt_bin3_down", &wt_tau_id_pt_bin3_down_);
        outtree_->Branch("wt_tau_id_pt_bin4_down", &wt_tau_id_pt_bin4_down_);
        outtree_->Branch("wt_tau_id_pt_bin5_down", &wt_tau_id_pt_bin5_down_);
        outtree_->Branch("wt_tau_trg_dm0_up", &wt_tau_trg_dm0_up_);
        outtree_->Branch("wt_tau_trg_dm1_up", &wt_tau_trg_dm1_up_);
        outtree_->Branch("wt_tau_trg_dm10_up", &wt_tau_trg_dm10_up_);
        outtree_->Branch("wt_tau_trg_dm11_up", &wt_tau_trg_dm11_up_);
        outtree_->Branch("wt_tau_trg_dm0_down", &wt_tau_trg_dm0_down_);
        outtree_->Branch("wt_tau_trg_dm1_down", &wt_tau_trg_dm1_down_);
        outtree_->Branch("wt_tau_trg_dm10_down", &wt_tau_trg_dm10_down_);
        outtree_->Branch("wt_tau_trg_dm11_down", &wt_tau_trg_dm11_down_);

        outtree_->Branch("wt_tau_trg_mvadm0_up",    &wt_tau_trg_mvadm0_up_);
        outtree_->Branch("wt_tau_trg_mvadm1_up",    &wt_tau_trg_mvadm1_up_);
        outtree_->Branch("wt_tau_trg_mvadm2_up",    &wt_tau_trg_mvadm2_up_);
        outtree_->Branch("wt_tau_trg_mvadm10_up",   &wt_tau_trg_mvadm10_up_);
        outtree_->Branch("wt_tau_trg_mvadm11_up",   &wt_tau_trg_mvadm11_up_);
        outtree_->Branch("wt_tau_trg_mvadm0_down",  &wt_tau_trg_mvadm0_down_);
	outtree_->Branch("wt_tau_trg_mvadm1_down",  &wt_tau_trg_mvadm1_down_);
        outtree_->Branch("wt_tau_trg_mvadm2_down",  &wt_tau_trg_mvadm2_down_);
        outtree_->Branch("wt_tau_trg_mvadm10_down", &wt_tau_trg_mvadm10_down_);
        outtree_->Branch("wt_tau_trg_mvadm11_down", &wt_tau_trg_mvadm11_down_);

        outtree_->Branch("wt_ue_up", & wt_ue_up_);
        outtree_->Branch("wt_ue_down", & wt_ue_down_);

      }
     
      outtree_->Branch("wt_prefire", &wt_prefire_);
      outtree_->Branch("wt_zpt",&wt_zpt_);
      outtree_->Branch("wt_zpt_NLO",&wt_zpt_NLO_);

      if(!systematic_shift_) {
        outtree_->Branch("wt_prefire_up", &wt_prefire_up_);
        outtree_->Branch("wt_prefire_down", &wt_prefire_down_);
      }

      if(do_mssm_higgspt_){
        outtree_->Branch("wt_ggh_t", &wt_ggh_t_);
        outtree_->Branch("wt_ggh_b", &wt_ggh_b_);
        outtree_->Branch("wt_ggh_i", &wt_ggh_i_);
        outtree_->Branch("wt_ggH_t", &wt_ggH_t_);
        outtree_->Branch("wt_ggH_b", &wt_ggH_b_);
        outtree_->Branch("wt_ggH_i", &wt_ggH_i_);
        outtree_->Branch("wt_ggA_t", &wt_ggA_t_);
        outtree_->Branch("wt_ggA_b", &wt_ggA_b_);
        outtree_->Branch("wt_ggA_i", &wt_ggA_i_);   

        outtree_->Branch("wt_ggh_t_scale_up", &wt_ggh_t_scale_up_);
        outtree_->Branch("wt_ggh_b_scale_up", &wt_ggh_b_scale_up_);
        outtree_->Branch("wt_ggh_i_scale_up", &wt_ggh_i_scale_up_);
        outtree_->Branch("wt_ggA_t_scale_up", &wt_ggA_t_scale_up_);
        outtree_->Branch("wt_ggA_b_scale_up", &wt_ggA_b_scale_up_);
        outtree_->Branch("wt_ggA_i_scale_up", &wt_ggA_i_scale_up_); 

        outtree_->Branch("wt_ggh_t_scale_down", &wt_ggh_t_scale_down_);
        outtree_->Branch("wt_ggh_b_scale_down", &wt_ggh_b_scale_down_);
        outtree_->Branch("wt_ggh_i_scale_down", &wt_ggh_i_scale_down_);
        outtree_->Branch("wt_ggA_t_scale_down", &wt_ggA_t_scale_down_);
        outtree_->Branch("wt_ggA_b_scale_down", &wt_ggA_b_scale_down_);
        outtree_->Branch("wt_ggA_i_scale_down", &wt_ggA_i_scale_down_);

        outtree_->Branch("wt_ggh_t_hdamp_up", &wt_ggh_t_hdamp_up_);
        outtree_->Branch("wt_ggh_b_hdamp_up", &wt_ggh_b_hdamp_up_);
        outtree_->Branch("wt_ggh_i_hdamp_up", &wt_ggh_i_hdamp_up_);
        outtree_->Branch("wt_ggA_t_hdamp_up", &wt_ggA_t_hdamp_up_);
        outtree_->Branch("wt_ggA_b_hdamp_up", &wt_ggA_b_hdamp_up_);
        outtree_->Branch("wt_ggA_i_hdamp_up", &wt_ggA_i_hdamp_up_);

        outtree_->Branch("wt_ggh_t_hdamp_down", &wt_ggh_t_hdamp_down_);
        outtree_->Branch("wt_ggh_b_hdamp_down", &wt_ggh_b_hdamp_down_);
        outtree_->Branch("wt_ggh_i_hdamp_down", &wt_ggh_i_hdamp_down_);
        outtree_->Branch("wt_ggA_t_hdamp_down", &wt_ggA_t_hdamp_down_);
        outtree_->Branch("wt_ggA_b_hdamp_down", &wt_ggA_b_hdamp_down_);
        outtree_->Branch("wt_ggA_i_hdamp_down", &wt_ggA_i_hdamp_down_);
      }
      
      // fake-factor weights
      if (do_ff_weights_ && (channel_ == channel::et || channel_ == channel::mt || channel_ == channel::tt)){
           
        if (strategy_ == strategy::smsummer16 || strategy_ == strategy::cpsummer16 || strategy_ == strategy::legacy16  || strategy_ == strategy::cpdecays16 || strategy_ == strategy::cpsummer17 || strategy_ == strategy::cpdecays17 || strategy_ == strategy::cpdecays18) {
          outtree_->Branch("wt_ff_1"  , &wt_ff_1_); 
          outtree_->Branch("wt_ff_dmbins_1"  , &wt_ff_dmbins_1_);

          outtree_->Branch("wt_ff_mssm_1"  , &wt_ff_mssm_1_);
 
          outtree_->Branch("w_frac",             &w_frac_);
          outtree_->Branch("qcd_frac",             &qcd_frac_);
          outtree_->Branch("ttbar_frac",             &ttbar_frac_);
          if(do_ff_systematics_ && !systematic_shift_){
            outtree_->Branch("wt_ff_qcd_1"  , &wt_ff_qcd_1_);
            outtree_->Branch("wt_ff_wjets_syst_up_1"  , &wt_ff_wjets_syst_up_1_);
            outtree_->Branch("wt_ff_wjets_syst_down_1"  , &wt_ff_wjets_syst_down_1_);
            outtree_->Branch("wt_ff_ttbar_syst_up_1"  , &wt_ff_ttbar_syst_up_1_);
            outtree_->Branch("wt_ff_ttbar_syst_down_1"  , &wt_ff_ttbar_syst_down_1_);
            outtree_->Branch("wt_ff_qcd_syst_up_1"  , &wt_ff_qcd_syst_up_1_);
            outtree_->Branch("wt_ff_qcd_syst_down_1"  , &wt_ff_qcd_syst_down_1_);
            outtree_->Branch("wt_ff_qcd_met_up_1"  , &wt_ff_qcd_met_up_1_);
            outtree_->Branch("wt_ff_qcd_met_down_1"  , &wt_ff_qcd_met_down_1_);
            outtree_->Branch("wt_ff_ttbar_met_up_1"  , &wt_ff_ttbar_met_up_1_);
            outtree_->Branch("wt_ff_ttbar_met_down_1"  , &wt_ff_ttbar_met_down_1_);

            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_up_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_down_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_up_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_down_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_up_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_down_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_up_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_down_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_up_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_down_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_up_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_down_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm1_up_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm1_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm1_down_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm1_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm1_up_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm1_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm1_down_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm1_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm1_up_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm1_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm1_down_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm1_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm2_up_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm2_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm2_down_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm2_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm2_up_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm2_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm2_down_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm2_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm2_up_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm2_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm2_down_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm2_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm10_up_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm10_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm10_down_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm10_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm10_up_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm10_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm10_down_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm10_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm10_up_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm10_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm10_down_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm10_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm11_up_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm11_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet0_mvadm11_down_1"  , &wt_ff_qcd_stat_unc1_njet0_mvadm11_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm11_up_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm11_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet1_mvadm11_down_1"  , &wt_ff_qcd_stat_unc1_njet1_mvadm11_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm11_up_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm11_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc1_njet2_mvadm11_down_1"  , &wt_ff_qcd_stat_unc1_njet2_mvadm11_down_1_);            
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_lt3_up_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_lt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_lt3_down_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_lt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_up_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_down_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_up_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_down_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_up_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_down_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_up_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_down_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_up_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_down_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm1_up_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm1_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm1_down_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm1_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm1_up_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm1_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm1_down_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm1_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm1_up_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm1_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm1_down_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm1_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm2_up_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm2_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm2_down_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm2_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm2_up_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm2_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm2_down_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm2_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm2_up_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm2_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm2_down_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm2_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm10_up_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm10_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm10_down_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm10_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm10_up_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm10_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm10_down_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm10_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm10_up_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm10_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm10_down_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm10_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm11_up_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm11_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet0_mvadm11_down_1"  , &wt_ff_qcd_stat_unc2_njet0_mvadm11_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm11_up_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm11_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet1_mvadm11_down_1"  , &wt_ff_qcd_stat_unc2_njet1_mvadm11_down_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm11_up_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm11_up_1_);
            outtree_->Branch("wt_ff_qcd_stat_unc2_njet2_mvadm11_down_1"  , &wt_ff_qcd_stat_unc2_njet2_mvadm11_down_1_);

            outtree_->Branch("wt_ff_dmbins_qcd_1"  , &wt_ff_dmbins_qcd_1_);
            outtree_->Branch("wt_ff_dmbins_wjets_1"  , &wt_ff_dmbins_wjets_1_);
            outtree_->Branch("wt_ff_dmbins_ttbar_1"  , &wt_ff_dmbins_ttbar_1_);
            outtree_->Branch("wt_ff_dmbins_wjets_syst_up_1"  , &wt_ff_dmbins_wjets_syst_up_1_);
            outtree_->Branch("wt_ff_dmbins_wjets_syst_down_1"  , &wt_ff_dmbins_wjets_syst_down_1_);
            outtree_->Branch("wt_ff_dmbins_ttbar_syst_up_1"  , &wt_ff_dmbins_ttbar_syst_up_1_);
            outtree_->Branch("wt_ff_dmbins_ttbar_syst_down_1"  , &wt_ff_dmbins_ttbar_syst_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_syst_up_1"  , &wt_ff_dmbins_qcd_syst_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_syst_down_1"  , &wt_ff_dmbins_qcd_syst_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_met_up_1"  , &wt_ff_dmbins_qcd_met_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_met_down_1"  , &wt_ff_dmbins_qcd_met_down_1_);
            outtree_->Branch("wt_ff_dmbins_ttbar_met_up_1"  , &wt_ff_dmbins_ttbar_met_up_1_);
            outtree_->Branch("wt_ff_dmbins_ttbar_met_down_1"  , &wt_ff_dmbins_ttbar_met_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_up_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_down_1"  , &wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_down_1_);            
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_down_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_up_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_up_1_);
            outtree_->Branch("wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_down_1"  , &wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_down_1_);

            outtree_->Branch("wt_ff_mssm_qcd_syst_up_1"  , &wt_ff_mssm_qcd_syst_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_syst_down_1"  , &wt_ff_mssm_qcd_syst_down_1_);
            outtree_->Branch("wt_ff_mssm_wjets_syst_up_1"  , &wt_ff_mssm_wjets_syst_up_1_);
            outtree_->Branch("wt_ff_mssm_wjets_syst_down_1"  , &wt_ff_mssm_wjets_syst_down_1_);
            outtree_->Branch("wt_ff_mssm_wjets_syst_extrap_up_1"  , &wt_ff_mssm_wjets_syst_extrap_up_1_);
            outtree_->Branch("wt_ff_mssm_wjets_syst_extrap_down_1"  , &wt_ff_mssm_wjets_syst_extrap_down_1_);
            outtree_->Branch("wt_ff_mssm_ttbar_syst_up_1"  , &wt_ff_mssm_ttbar_syst_up_1_);
            outtree_->Branch("wt_ff_mssm_ttbar_syst_down_1"  , &wt_ff_mssm_ttbar_syst_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_up_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_down_1"  , &wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_up_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_down_1"  , &wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_down_1_);

            outtree_->Branch("wt_ff_mssm_syst_alt_func_up_1"  , &wt_ff_mssm_syst_alt_func_up_1_);
            outtree_->Branch("wt_ff_mssm_syst_alt_func_down_1"  , &wt_ff_mssm_syst_alt_func_down_1_);
            outtree_->Branch("wt_ff_mssm_qcd_syst_met_closure_up_1"  , &wt_ff_mssm_qcd_syst_met_closure_up_1_);
            outtree_->Branch("wt_ff_mssm_qcd_syst_met_closure_down_1"  , &wt_ff_mssm_qcd_syst_met_closure_down_1_);

            if(channel_ == channel::tt) {
              outtree_->Branch("wt_ff_mssm_qcd_stat_dR_unc1_up_1"  , &wt_ff_mssm_qcd_stat_dR_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_dR_unc1_down_1"  , &wt_ff_mssm_qcd_stat_dR_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_dR_unc2_up_1"  , &wt_ff_mssm_qcd_stat_dR_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_dR_unc2_down_1"  , &wt_ff_mssm_qcd_stat_dR_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_pt_unc1_up_1"  , &wt_ff_mssm_qcd_stat_pt_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_pt_unc1_down_1"  , &wt_ff_mssm_qcd_stat_pt_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_pt_unc2_up_1"  , &wt_ff_mssm_qcd_stat_pt_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_pt_unc2_down_1"  , &wt_ff_mssm_qcd_stat_pt_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_syst_pt_2_closure_up_1"  , &wt_ff_mssm_qcd_syst_pt_2_closure_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_syst_pt_2_closure_down_1"  , &wt_ff_mssm_qcd_syst_pt_2_closure_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_syst_dR_closure_up_1"  , &wt_ff_mssm_qcd_syst_dr_closure_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_syst_dR_closure_down_1"  , &wt_ff_mssm_qcd_syst_dr_closure_down_1_);
            }

            if(channel_ == channel::et || channel_ == channel::mt) {

              outtree_->Branch("wt_ff_wjets_1"  , &wt_ff_wjets_1_);
              outtree_->Branch("wt_ff_ttbar_1"  , &wt_ff_ttbar_1_);

              outtree_->Branch("wt_ff_mssm_qcd_1"  , &wt_ff_mssm_qcd_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_1"  , &wt_ff_mssm_ttbar_1_);
              outtree_->Branch("wt_ff_mssm_wjets_1"  , &wt_ff_mssm_wjets_1_);

              outtree_->Branch("ff_frac_mssm_wjets"  , &ff_frac_mssm_wjets_);
              outtree_->Branch("ff_frac_mssm_ttbar"  , &ff_frac_mssm_ttbar_);
              outtree_->Branch("ff_frac_mssm_qcd"  , &ff_frac_mssm_qcd_);

              outtree_->Branch("wt_ff_mssm_qcd_stat_ss_njets0_unc1_up_1"  , &wt_ff_mssm_qcd_stat_ss_njets0_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_ss_njets0_unc1_down_1"  , &wt_ff_mssm_qcd_stat_ss_njets0_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_ss_njets1_unc1_up_1"  , &wt_ff_mssm_qcd_stat_ss_njets1_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_ss_njets1_unc1_down_1"  , &wt_ff_mssm_qcd_stat_ss_njets1_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_os_unc1_up_1"  , &wt_ff_mssm_qcd_stat_os_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_os_unc1_down_1"  , &wt_ff_mssm_qcd_stat_os_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_ss_njets0_unc2_up_1"  , &wt_ff_mssm_qcd_stat_ss_njets0_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_ss_njets0_unc2_down_1"  , &wt_ff_mssm_qcd_stat_ss_njets0_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_ss_njets1_unc2_up_1"  , &wt_ff_mssm_qcd_stat_ss_njets1_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_ss_njets1_unc2_down_1"  , &wt_ff_mssm_qcd_stat_ss_njets1_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_os_unc2_up_1"  , &wt_ff_mssm_qcd_stat_os_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_os_unc2_down_1"  , &wt_ff_mssm_qcd_stat_os_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_up_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_down_1"  , &wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_up_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_down_1"  , &wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_met_njets0_unc1_up_1"  , &wt_ff_mssm_wjets_stat_met_njets0_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_met_njets0_unc1_down_1"  , &wt_ff_mssm_wjets_stat_met_njets0_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_met_njets1_unc1_up_1"  , &wt_ff_mssm_wjets_stat_met_njets1_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_met_njets1_unc1_down_1"  , &wt_ff_mssm_wjets_stat_met_njets1_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_up_1"  , &wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_down_1"  , &wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_up_1"  , &wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_down_1"  , &wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_extrap_unc1_up_1"  , &wt_ff_mssm_wjets_stat_extrap_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_extrap_unc1_down_1"  , &wt_ff_mssm_wjets_stat_extrap_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_met_njets0_unc2_up_1"  , &wt_ff_mssm_wjets_stat_met_njets0_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_met_njets0_unc2_down_1"  , &wt_ff_mssm_wjets_stat_met_njets0_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_met_njets1_unc2_up_1"  , &wt_ff_mssm_wjets_stat_met_njets1_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_met_njets1_unc2_down_1"  , &wt_ff_mssm_wjets_stat_met_njets1_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_up_1"  , &wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_down_1"  , &wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_up_1"  , &wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_down_1"  , &wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_extrap_unc2_up_1"  , &wt_ff_mssm_wjets_stat_extrap_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_stat_extrap_unc2_down_1"  , &wt_ff_mssm_wjets_stat_extrap_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_up_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_down_1"  , &wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_l_pt_unc1_up_1"  , &wt_ff_mssm_ttbar_stat_l_pt_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_l_pt_unc1_down_1"  , &wt_ff_mssm_ttbar_stat_l_pt_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_l_pt_unc2_up_1"  , &wt_ff_mssm_ttbar_stat_l_pt_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_l_pt_unc2_down_1"  , &wt_ff_mssm_ttbar_stat_l_pt_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_met_unc1_up_1"  , &wt_ff_mssm_ttbar_stat_met_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_met_unc1_down_1"  , &wt_ff_mssm_ttbar_stat_met_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_met_unc2_up_1"  , &wt_ff_mssm_ttbar_stat_met_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_stat_met_unc2_down_1"  , &wt_ff_mssm_ttbar_stat_met_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_frac_up_1"  , &wt_ff_mssm_wjets_frac_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_frac_down_1"  , &wt_ff_mssm_wjets_frac_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_frac_up_1"  , &wt_ff_mssm_ttbar_frac_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_frac_down_1"  , &wt_ff_mssm_ttbar_frac_down_1_);
              outtree_->Branch("wt_ff_mssm_low_pt_up_1"  , &wt_ff_mssm_low_pt_up_1_);
              outtree_->Branch("wt_ff_mssm_low_pt_down_1"  , &wt_ff_mssm_low_pt_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_l_pt_unc1_up_1"  , &wt_ff_mssm_qcd_stat_l_pt_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_l_pt_unc1_down_1"  , &wt_ff_mssm_qcd_stat_l_pt_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_l_pt_unc2_up_1"  , &wt_ff_mssm_qcd_stat_l_pt_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_l_pt_unc2_down_1"  , &wt_ff_mssm_qcd_stat_l_pt_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_iso_unc1_up_1"  , &wt_ff_mssm_qcd_stat_iso_unc1_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_iso_unc1_down_1"  , &wt_ff_mssm_qcd_stat_iso_unc1_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_iso_unc2_up_1"  , &wt_ff_mssm_qcd_stat_iso_unc2_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_stat_iso_unc2_down_1"  , &wt_ff_mssm_qcd_stat_iso_unc2_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_syst_iso_up_1"  , &wt_ff_mssm_qcd_syst_iso_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_syst_iso_down_1"  , &wt_ff_mssm_qcd_syst_iso_down_1_);

              outtree_->Branch("wt_ff_mssm_ttbar_syst_met_closure_up_1"  , &wt_ff_mssm_ttbar_syst_met_closure_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_syst_met_closure_down_1"  , &wt_ff_mssm_ttbar_syst_met_closure_down_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_syst_l_pt_closure_up_1"  , &wt_ff_mssm_ttbar_syst_l_pt_closure_up_1_);
              outtree_->Branch("wt_ff_mssm_ttbar_syst_l_pt_closure_down_1"  , &wt_ff_mssm_ttbar_syst_l_pt_closure_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_syst_met_closure_up_1"  , &wt_ff_mssm_wjets_syst_met_closure_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_syst_met_closure_down_1"  , &wt_ff_mssm_wjets_syst_met_closure_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_syst_l_pt_closure_up_1"  , &wt_ff_mssm_wjets_syst_l_pt_closure_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_syst_l_pt_closure_down_1"  , &wt_ff_mssm_wjets_syst_l_pt_closure_down_1_);
              outtree_->Branch("wt_ff_mssm_wjets_syst_bkg_up_1"  , &wt_ff_mssm_wjets_syst_bkg_up_1_);
              outtree_->Branch("wt_ff_mssm_wjets_syst_bkg_down_1"  , &wt_ff_mssm_wjets_syst_bkg_down_1_);
              outtree_->Branch("wt_ff_mssm_qcd_syst_bkg_up_1"  , &wt_ff_mssm_qcd_syst_bkg_up_1_);
              outtree_->Branch("wt_ff_mssm_qcd_syst_bkg_down_1"  , &wt_ff_mssm_qcd_syst_bkg_down_1_);

              outtree_->Branch("wt_ff_qcd_l_pt_up_1"  , &wt_ff_qcd_l_pt_up_1_);
              outtree_->Branch("wt_ff_qcd_l_pt_down_1"  , &wt_ff_qcd_l_pt_down_1_);
              outtree_->Branch("wt_ff_wjets_met_up_1"  , &wt_ff_wjets_met_up_1_);
              outtree_->Branch("wt_ff_wjets_met_down_1"  , &wt_ff_wjets_met_down_1_);
              outtree_->Branch("wt_ff_wjets_l_pt_up_1"  , &wt_ff_wjets_l_pt_up_1_);
              outtree_->Branch("wt_ff_wjets_l_pt_down_1"  , &wt_ff_wjets_l_pt_down_1_);

              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_up_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_down_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_up_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_down_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_up_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_down_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_up_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_down_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_up_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_down_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_up_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_down_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm1_up_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm1_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm1_down_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm1_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm1_up_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm1_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm1_down_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm1_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm1_up_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm1_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm1_down_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm1_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm2_up_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm2_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm2_down_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm2_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm2_up_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm2_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm2_down_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm2_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm2_up_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm2_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm2_down_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm2_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm10_up_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm10_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm10_down_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm10_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm10_up_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm10_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm10_down_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm10_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm10_up_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm10_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm10_down_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm10_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm11_up_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm11_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet0_mvadm11_down_1"  , &wt_ff_wjets_stat_unc1_njet0_mvadm11_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm11_up_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm11_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet1_mvadm11_down_1"  , &wt_ff_wjets_stat_unc1_njet1_mvadm11_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm11_up_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm11_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc1_njet2_mvadm11_down_1"  , &wt_ff_wjets_stat_unc1_njet2_mvadm11_down_1_);            
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_up_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_down_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_up_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_down_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_up_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_down_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_up_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_down_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_up_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_down_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_up_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_down_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm1_up_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm1_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm1_down_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm1_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm1_up_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm1_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm1_down_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm1_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm1_up_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm1_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm1_down_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm1_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm2_up_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm2_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm2_down_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm2_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm2_up_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm2_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm2_down_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm2_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm2_up_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm2_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm2_down_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm2_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm10_up_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm10_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm10_down_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm10_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm10_up_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm10_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm10_down_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm10_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm10_up_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm10_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm10_down_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm10_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm11_up_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm11_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet0_mvadm11_down_1"  , &wt_ff_wjets_stat_unc2_njet0_mvadm11_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm11_up_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm11_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet1_mvadm11_down_1"  , &wt_ff_wjets_stat_unc2_njet1_mvadm11_down_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm11_up_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm11_up_1_);
              outtree_->Branch("wt_ff_wjets_stat_unc2_njet2_mvadm11_down_1"  , &wt_ff_wjets_stat_unc2_njet2_mvadm11_down_1_);

              outtree_->Branch("wt_ff_dmbins_qcd_l_pt_up_1"  , &wt_ff_dmbins_qcd_l_pt_up_1_);
              outtree_->Branch("wt_ff_dmbins_qcd_l_pt_down_1"  , &wt_ff_dmbins_qcd_l_pt_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_met_up_1"  , &wt_ff_dmbins_wjets_met_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_met_down_1"  , &wt_ff_dmbins_wjets_met_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_l_pt_up_1"  , &wt_ff_dmbins_wjets_l_pt_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_l_pt_down_1"  , &wt_ff_dmbins_wjets_l_pt_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_up_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_down_1"  , &wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_down_1_);            
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_down_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_up_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_up_1_);
              outtree_->Branch("wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_down_1"  , &wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_down_1_);
            }

          }
         }
      }
      
     if(do_qcd_scale_wts_){
       outtree_->Branch("wt_mur1_muf1",    &scale1_);
       outtree_->Branch("wt_mur1_muf2",    &scale2_);
       outtree_->Branch("wt_mur1_muf0p5",  &scale3_);
       outtree_->Branch("wt_mur2_muf1",    &scale4_);
       outtree_->Branch("wt_mur2_muf2",    &scale5_);
       outtree_->Branch("wt_mur2_muf0p5",  &scale6_);
       outtree_->Branch("wt_mur0p5_muf1",  &scale7_);
       outtree_->Branch("wt_mur0p5_muf2",  &scale8_);
       outtree_->Branch("wt_mur0p5_muf0p5",&scale9_);
     }
      
      outtree_->Branch("os",                &os_);
      outtree_->Branch("m_sv",              &m_sv_.var_double);
      outtree_->Branch("mt_sv",             &mt_sv_.var_double);
      outtree_->Branch("m_vis",             &m_vis_.var_double);
      outtree_->Branch("pt_vis",            &pt_vis_.var_double);
      outtree_->Branch("pt_h",              &pt_h_.var_double);
      outtree_->Branch("eta_h",             &eta_h_);
      outtree_->Branch("pt_tt",             &pt_tt_.var_double);
      outtree_->Branch("pt_tt_inc_met",             pt_tt_inc_met_);
      outtree_->Branch("mt_tot",            &mt_tot_.var_double);
      outtree_->Branch("mt_lep",            &mt_lep_.var_double);
      outtree_->Branch("mt_2",              &mt_2_.var_double);
      outtree_->Branch("mt_1",              &mt_1_.var_double);
      outtree_->Branch("pf_mt_1"    , &pf_mt_1_     );
      outtree_->Branch("pzeta",             &pzeta_.var_double);
      outtree_->Branch("iso_1",             &iso_1_.var_double);
      outtree_->Branch("iso_2",             &iso_2_.var_double);
      outtree_->Branch("leptonveto",        &lepton_veto_);
      outtree_->Branch("dilepton_veto",     &dilepton_veto_);
      outtree_->Branch("extraelec_veto",    &extraelec_veto_);
      outtree_->Branch("extramuon_veto",    &extramuon_veto_);
      outtree_->Branch("met",               &met_.var_double);
      outtree_->Branch("pf_met",               &pf_met_);
      outtree_->Branch("n_jets",            &n_jets_);
      outtree_->Branch("n_bjets",           &n_bjets_);
      outtree_->Branch("n_deepbjets",       &n_deepbjets_);
      outtree_->Branch("deepbpt_1",       &deepbpt_1_);
      outtree_->Branch("deepbeta_1",       &deepbeta_1_);
      outtree_->Branch("deep_probb_1",       &deep_probb_1_);
      outtree_->Branch("deep_probbb_1",       &deep_probbb_1_);
      outtree_->Branch("deep_problepb_1",       &deep_problepb_1_);
      outtree_->Branch("deepcsv_probb_1",       &deepcsv_probb_1_);
      outtree_->Branch("deepcsv_probbb_1",       &deepcsv_probbb_1_);

      outtree_->Branch("n_prebjets",        &n_prebjets_);
      outtree_->Branch("n_loose_bjets",     &n_loose_bjets_);
      outtree_->Branch("n_btag",            &n_btag_);
      outtree_->Branch("n_loose_btag",      &n_loose_btag_);
      outtree_->Branch("mjj",               &mjj_.var_double);
      outtree_->Branch("jdeta",             &jdeta_.var_double);
      outtree_->Branch("jdphi",             &jdphi_);
      outtree_->Branch("dphi_jtt",          &dphi_jtt_);
      //outtree_->Branch("shifted_dphi_jtt",  &shifted_dphi_jtt_);
      //outtree_->Branch("dphi_jtt_smear",    &dphi_jtt_smear_);
      //outtree_->Branch("shifted_dphi_jtt_smear",    &shifted_dphi_jtt_smear_);
      //outtree_->Branch("smear_jet_delta",   &smear_jet_delta_);
      //outtree_->Branch("residual_pt",       &residual_pt_);
      //outtree_->Branch("residual_phi",      &residual_phi_);
      //outtree_->Branch("dphi_residual_tt",  &dphi_residual_tt_);
      outtree_->Branch("dijetpt",           &dijetpt_);
      outtree_->Branch("mt_1_nomu"    , &mt_1_nomu_     );
      outtree_->Branch("mu_pt", &mu_pt_);
      outtree_->Branch("sjdphi",             &sjdphi_);
      outtree_->Branch("n_lowpt_jets",      &n_lowpt_jets_);
      outtree_->Branch("pt_2",              &pt_2_.var_double);
      outtree_->Branch("pt_1",              &pt_1_.var_double);
      outtree_->Branch("eta_1",             &eta_1_.var_double);
      outtree_->Branch("eta_2",             &eta_2_.var_double);
      outtree_->Branch("gen_match_1", &gen_match_1_);
      outtree_->Branch("gen_match_2", &gen_match_2_);
      outtree_->Branch("gen_match_1_pt", &gen_match_1_pt_);
      outtree_->Branch("gen_match_2_pt", &gen_match_2_pt_);
      outtree_->Branch("gen_match_undecayed_1_pt", &gen_match_undecayed_1_pt_);
      outtree_->Branch("gen_match_undecayed_2_pt", &gen_match_undecayed_2_pt_);
      outtree_->Branch("gen_sjdphi", &gen_sjdphi_);
      outtree_->Branch("ngenjets" , &ngenjets_);
      outtree_->Branch("ngenjets20" , &ngenjets20_);
      outtree_->Branch("genjetpt_1" , &genjetpt_1_);
      outtree_->Branch("genjetpt_2" , &genjetpt_2_);
      outtree_->Branch("genjeteta_1" , &genjeteta_1_);
      outtree_->Branch("genjeteta_2" , &genjeteta_2_);
      outtree_->Branch("genM", &gen_m_);
      outtree_->Branch("genpT", &gen_pt_);
      outtree_->Branch("m_1", &m_1_, "m_1/F");
      outtree_->Branch("m_2", &m_2_, "m_2/F");

      outtree_->Branch("wt_cp_sm", &wt_cp_sm_);
      outtree_->Branch("wt_cp_ps", &wt_cp_ps_);
      outtree_->Branch("wt_cp_mm", &wt_cp_mm_);

      outtree_->Branch("wt_cp_sm_alt"       , &wt_cp_sm_alt_);
      outtree_->Branch("wt_cp_ps_alt"       , &wt_cp_ps_alt_);
      outtree_->Branch("wt_cp_mm_alt"       , &wt_cp_mm_alt_);

      outtree_->Branch("gen_phitt", &gen_phitt_);

      outtree_->Branch("wt_cp_prod_sm",&wt_cp_prod_sm_);
      outtree_->Branch("wt_cp_prod_ps",&wt_cp_prod_ps_);
      outtree_->Branch("wt_cp_prod_mm",&wt_cp_prod_mm_);

      outtree_->Branch("mva_dm_1", &tau_mva_decay_mode_1_);
      outtree_->Branch("mva_dm_2", &tau_mva_decay_mode_2_);

      outtree_->Branch("a1_flag", &a1_flag_);
      outtree_->Branch("a1_flag_2", &a1_flag_2_);
      outtree_->Branch("aco_angle_1", &aco_angle_1_);
      outtree_->Branch("aco_angle_2", &aco_angle_2_);
      outtree_->Branch("aco_angle_3", &aco_angle_3_);
      outtree_->Branch("aco_angle_4", &aco_angle_4_);
      outtree_->Branch("aco_angle_5", &aco_angle_5_);
      outtree_->Branch("aco_angle_6", &aco_angle_6_);
      outtree_->Branch("pv_angle", &pv_angle_);

      outtree_->Branch("sv_x_2", &sv_x_2_);
      outtree_->Branch("sv_y_2", &sv_y_2_);
      outtree_->Branch("sv_z_2", &sv_z_2_);
      outtree_->Branch("sv_mag_2", &sv_mag_2_);
      outtree_->Branch("sv_dphi_2", &sv_dphi_2_);
      outtree_->Branch("sv_deta_2", &sv_deta_2_);
      outtree_->Branch("hasSV_2", &hasSV_2_);

      outtree_->Branch("aco_angle_rand", &aco_angle_rand_);
      outtree_->Branch("aco_sign", &aco_sign_);
      outtree_->Branch("aco_sign_rand", &aco_sign_rand_);
      outtree_->Branch("lead_pt_1", &lead_pt_1_);
      outtree_->Branch("lead_pt_2", &lead_pt_2_);
      outtree_->Branch("alpha1_1", &alpha1_1_);
      outtree_->Branch("alpha1_2", &alpha1_2_);
      outtree_->Branch("alpha2_1", &alpha2_1_);
      outtree_->Branch("alpha2_2", &alpha2_2_);

      outtree_->Branch("mass0",         &mass0_);
      outtree_->Branch("mass1",         &mass1_);
      outtree_->Branch("mass2",         &mass2_);
      outtree_->Branch("y_1_1", &y_1_1_);
      outtree_->Branch("y_1_2", &y_1_2_);
      outtree_->Branch("y_2_2", &y_2_2_);
      outtree_->Branch("y_3_2", &y_3_2_);
      outtree_->Branch("y_4_2", &y_4_2_);
      outtree_->Branch("rand", &rand_);

      outtree_->Branch("IC_BDT_max_score", &IC_BDT_max_score_);
      outtree_->Branch("IC_BDT_max_index", &IC_BDT_max_index_);

      // moved these here from !(systematics_shifts) because want to use
      // them in BDT
      outtree_->Branch("phi_1",             &phi_1_.var_double);
      outtree_->Branch("phi_2",             &phi_2_.var_double);
      outtree_->Branch("dphi",              &dphi_);
      outtree_->Branch("dR",              &dR_);
      outtree_->Branch("jphi_1",            &jphi_1_, "jphi_1/F");
      outtree_->Branch("jphi_2",            &jphi_2_, "jphi_2/F");
      outtree_->Branch("bpt_1",             &bpt_1_.var_float);
      outtree_->Branch("bpt_2",             &bpt_2_.var_float);
      outtree_->Branch("beta_1",            &beta_1_.var_float);
      outtree_->Branch("beta_2",            &beta_2_.var_float);
      outtree_->Branch("bcsv_1",            &bcsv_1_.var_float);
      outtree_->Branch("bcsv_2",            &bcsv_2_);
      outtree_->Branch("uncorrmet", &uncorrmet_.var_float, "met/F");
      outtree_->Branch("uncorrmetphi", &uncorrmet_phi_.var_float, "met_phi/F");
      outtree_->Branch("met_dphi_1",             &met_dphi_1_);
      outtree_->Branch("met_dphi_2",             &met_dphi_2_);
      outtree_->Branch("newmet_dphi_1",             &newmet_dphi_1_);
      outtree_->Branch("newmet_dphi_2",             &newmet_dphi_2_);
      outtree_->Branch("fake_met_dphi_1",             &fake_met_dphi_1_);
      outtree_->Branch("fake_met_dphi_2",             &fake_met_dphi_2_);
      outtree_->Branch("newmet",             &newmet_);
      outtree_->Branch("fake_met",             &fake_met_);
      outtree_->Branch("gen_met",             &gen_met_);

      outtree_->Branch("fake_tau_met_dphi_1",             &fake_tau_met_dphi_1_);
      outtree_->Branch("fake_tau_met_dphi_2",             &fake_tau_met_dphi_2_);
      outtree_->Branch("fake_tau_met",             &fake_tau_met_);
      outtree_->Branch("gen_tau_met",             &gen_tau_met_);

      outtree_->Branch("metx", &metx_);
      outtree_->Branch("mety", &mety_);
      outtree_->Branch("gen_tau_metx", &gen_tau_metx_);
      outtree_->Branch("gen_tau_mety", &gen_tau_mety_);
      outtree_->Branch("gen_metx", &gen_metx_);
      outtree_->Branch("gen_mety", &gen_mety_);
      outtree_->Branch("fake_tau_metx", &fake_tau_metx_);
      outtree_->Branch("fake_tau_mety", &fake_tau_mety_);

      outtree_->Branch("qcd_frac_score",             &qcd_frac_score_);
      outtree_->Branch("w_frac_score",             &w_frac_score_);

      // NN Tau IDs
      outtree_->Branch("deepTauVsEle_iso_1",       &deepTauVsEle_iso_1_);
      outtree_->Branch("deepTauVsEle_iso_2",       &deepTauVsEle_iso_2_);
      outtree_->Branch("deepTauVsMu_iso_1",        &deepTauVsMu_iso_1_);
      outtree_->Branch("deepTauVsMu_iso_2",        &deepTauVsMu_iso_2_);

      outtree_->Branch("deepTauVsJets_vvvloose_1", &deepTauVsJets_vvvloose_1_);
      outtree_->Branch("deepTauVsJets_vvvloose_2", &deepTauVsJets_vvvloose_2_);
      outtree_->Branch("deepTauVsJets_vvloose_1",  &deepTauVsJets_vvloose_1_);
      outtree_->Branch("deepTauVsJets_vvloose_2",  &deepTauVsJets_vvloose_2_);
      outtree_->Branch("deepTauVsJets_vloose_1",   &deepTauVsJets_vloose_1_);
      outtree_->Branch("deepTauVsJets_vloose_2",   &deepTauVsJets_vloose_2_);
      outtree_->Branch("deepTauVsJets_loose_1",    &deepTauVsJets_loose_1_);
      outtree_->Branch("deepTauVsJets_loose_2",    &deepTauVsJets_loose_2_);
      outtree_->Branch("deepTauVsJets_medium_1",   &deepTauVsJets_medium_1_);
      outtree_->Branch("deepTauVsJets_medium_2",   &deepTauVsJets_medium_2_);
      outtree_->Branch("deepTauVsJets_tight_1",    &deepTauVsJets_tight_1_);
      outtree_->Branch("deepTauVsJets_tight_2",    &deepTauVsJets_tight_2_);
      outtree_->Branch("deepTauVsJets_vtight_1",   &deepTauVsJets_vtight_1_);
      outtree_->Branch("deepTauVsJets_vtight_2",   &deepTauVsJets_vtight_2_);
      outtree_->Branch("deepTauVsJets_vvtight_1",  &deepTauVsJets_vvtight_1_);
      outtree_->Branch("deepTauVsJets_vvtight_2",  &deepTauVsJets_vvtight_2_);

      outtree_->Branch("deepTauVsEle_vvvloose_1",  &deepTauVsEle_vvvloose_1_);
      outtree_->Branch("deepTauVsEle_vvvloose_2",  &deepTauVsEle_vvvloose_2_);
      outtree_->Branch("deepTauVsEle_vvloose_1",   &deepTauVsEle_vvloose_1_);
      outtree_->Branch("deepTauVsEle_vvloose_2",   &deepTauVsEle_vvloose_2_);
      outtree_->Branch("deepTauVsEle_vloose_1",    &deepTauVsEle_vloose_1_);
      outtree_->Branch("deepTauVsEle_vloose_2",    &deepTauVsEle_vloose_2_);
      outtree_->Branch("deepTauVsEle_loose_1",     &deepTauVsEle_loose_1_);
      outtree_->Branch("deepTauVsEle_loose_2",     &deepTauVsEle_loose_2_);
      outtree_->Branch("deepTauVsEle_medium_1",    &deepTauVsEle_medium_1_);
      outtree_->Branch("deepTauVsEle_medium_2",    &deepTauVsEle_medium_2_);
      outtree_->Branch("deepTauVsEle_tight_1",     &deepTauVsEle_tight_1_);
      outtree_->Branch("deepTauVsEle_tight_2",     &deepTauVsEle_tight_2_);
      outtree_->Branch("deepTauVsEle_vtight_1",    &deepTauVsEle_vtight_1_);
      outtree_->Branch("deepTauVsEle_vtight_2",    &deepTauVsEle_vtight_2_);
      outtree_->Branch("deepTauVsEle_vvtight_1",   &deepTauVsEle_vvtight_1_);
      outtree_->Branch("deepTauVsEle_vvtight_2",   &deepTauVsEle_vvtight_2_);
      outtree_->Branch("deepTauVsMu_vvvloose_1",  &deepTauVsMu_vvvloose_1_);
      outtree_->Branch("deepTauVsMu_vvvloose_2",  &deepTauVsMu_vvvloose_2_);
      outtree_->Branch("deepTauVsMu_vvloose_1",   &deepTauVsMu_vvloose_1_);
      outtree_->Branch("deepTauVsMu_vvloose_2",   &deepTauVsMu_vvloose_2_);
      outtree_->Branch("deepTauVsMu_vloose_1",    &deepTauVsMu_vloose_1_);
      outtree_->Branch("deepTauVsMu_vloose_2",    &deepTauVsMu_vloose_2_);
      outtree_->Branch("deepTauVsMu_loose_1",     &deepTauVsMu_loose_1_);
      outtree_->Branch("deepTauVsMu_loose_2",     &deepTauVsMu_loose_2_);
      outtree_->Branch("deepTauVsMu_medium_1",    &deepTauVsMu_medium_1_);
      outtree_->Branch("deepTauVsMu_medium_2",    &deepTauVsMu_medium_2_);
      outtree_->Branch("deepTauVsMu_tight_1",     &deepTauVsMu_tight_1_);
      outtree_->Branch("deepTauVsMu_tight_2",     &deepTauVsMu_tight_2_);
      outtree_->Branch("deepTauVsMu_vtight_1",    &deepTauVsMu_vtight_1_);
      outtree_->Branch("deepTauVsMu_vtight_2",    &deepTauVsMu_vtight_2_);
      outtree_->Branch("deepTauVsMu_vvtight_1",   &deepTauVsMu_vvtight_1_);
      outtree_->Branch("deepTauVsMu_vvtight_2",   &deepTauVsMu_vvtight_2_);

      // DeepTau v2p5
      outtree_->Branch("deepTauVsJetsV2p5_vvvloose_1", &deepTauVsJetsV2p5_vvvloose_1_);
      outtree_->Branch("deepTauVsJetsV2p5_vvvloose_2", &deepTauVsJetsV2p5_vvvloose_2_);
      outtree_->Branch("deepTauVsJetsV2p5_vvloose_1",  &deepTauVsJetsV2p5_vvloose_1_);
      outtree_->Branch("deepTauVsJetsV2p5_vvloose_2",  &deepTauVsJetsV2p5_vvloose_2_);
      outtree_->Branch("deepTauVsJetsV2p5_vloose_1",   &deepTauVsJetsV2p5_vloose_1_);
      outtree_->Branch("deepTauVsJetsV2p5_vloose_2",   &deepTauVsJetsV2p5_vloose_2_);
      outtree_->Branch("deepTauVsJetsV2p5_loose_1",    &deepTauVsJetsV2p5_loose_1_);
      outtree_->Branch("deepTauVsJetsV2p5_loose_2",    &deepTauVsJetsV2p5_loose_2_);
      outtree_->Branch("deepTauVsJetsV2p5_medium_1",   &deepTauVsJetsV2p5_medium_1_);
      outtree_->Branch("deepTauVsJetsV2p5_medium_2",   &deepTauVsJetsV2p5_medium_2_);
      outtree_->Branch("deepTauVsJetsV2p5_tight_1",    &deepTauVsJetsV2p5_tight_1_);
      outtree_->Branch("deepTauVsJetsV2p5_tight_2",    &deepTauVsJetsV2p5_tight_2_);
      outtree_->Branch("deepTauVsJetsV2p5_vtight_1",   &deepTauVsJetsV2p5_vtight_1_);
      outtree_->Branch("deepTauVsJetsV2p5_vtight_2",   &deepTauVsJetsV2p5_vtight_2_);
      outtree_->Branch("deepTauVsJetsV2p5_vvtight_1",  &deepTauVsJetsV2p5_vvtight_1_);
      outtree_->Branch("deepTauVsJetsV2p5_vvtight_2",  &deepTauVsJetsV2p5_vvtight_2_);
      outtree_->Branch("deepTauVsEleV2p5_vvloose_1",   &deepTauVsEleV2p5_vvloose_1_);
      outtree_->Branch("deepTauVsEleV2p5_vvloose_2",   &deepTauVsEleV2p5_vvloose_2_);
      outtree_->Branch("deepTauVsEleV2p5_tight_1",     &deepTauVsEleV2p5_tight_1_);
      outtree_->Branch("deepTauVsEleV2p5_tight_2",     &deepTauVsEleV2p5_tight_2_);
      outtree_->Branch("deepTauVsMuV2p5_tight_1",      &deepTauVsMuV2p5_tight_1_);
      outtree_->Branch("deepTauVsMuV2p5_tight_2",      &deepTauVsMuV2p5_tight_2_);
      outtree_->Branch("deepTauVsMuV2p5_vloose_1",     &deepTauVsMuV2p5_vloose_1_);
      outtree_->Branch("deepTauVsMuV2p5_vloose_2",     &deepTauVsMuV2p5_vloose_2_);
      outtree_->Branch("iso_2_V2p5",    &iso_2_V2p5_);

      outtree_->Branch("tau_decay_mode_2",    &tau_decay_mode_2_);
      outtree_->Branch("tau_decay_mode_1",    &tau_decay_mode_1_);
      outtree_->Branch("trg_singleelectron",    &trg_singleelectron_);
      outtree_->Branch("trg_singlemuon",    &trg_singlemuon_);
      outtree_->Branch("trg_singlemuon24",    &trg_singlemuon24_);
      outtree_->Branch("trg_singlemuon27",    &trg_singlemuon27_);
      outtree_->Branch("trg_doubletau",    &trg_doubletau_);
      outtree_->Branch("trg_doubletau_mssm",    &trg_doubletau_mssm_);

      outtree_->Branch("trg_tt_monitoring_1",    &trg_tt_monitoring_1_);
      outtree_->Branch("trg_tt_monitoring_2",    &trg_tt_monitoring_2_);
      outtree_->Branch("trg_tt_monitoring_3",    &trg_tt_monitoring_3_);
      outtree_->Branch("trg_tt_monitoring_4",    &trg_tt_monitoring_4_);
      outtree_->Branch("trg_tt_monitoring_5",    &trg_tt_monitoring_5_);
      outtree_->Branch("trg_tt_monitoring_6",    &trg_tt_monitoring_6_);
      outtree_->Branch("trg_tt_monitoring_l1",    &trg_tt_monitoring_l1_);

      outtree_->Branch("trg_vbfdoubletau",    &trg_vbfdoubletau_);
      outtree_->Branch("trg_muonelectron",    &trg_muonelectron_);
      outtree_->Branch("trg_singletau_1",    &trg_singletau_1_);
      outtree_->Branch("trg_singletau_2",    &trg_singletau_2_);
      outtree_->Branch("trg_mutaucross",    &trg_mutaucross_);
      outtree_->Branch("trg_etaucross",    &trg_etaucross_);
      outtree_->Branch("trg_etau_matched", &trg_etau_matched_);
      outtree_->Branch("jpt_1",             &jpt_1_.var_double);
      outtree_->Branch("jpt_2",             &jpt_2_.var_double);
      outtree_->Branch("jnemfrac_1",             &jnemfrac_1_);
      outtree_->Branch("jeta_1",            &jeta_1_.var_double);
      outtree_->Branch("jeta_2",            &jeta_2_.var_double);
      outtree_->Branch("jmva_1",             &jmva_1_);
      outtree_->Branch("jmva_2",             &jmva_2_);

      if (channel_ == channel::em) {
        outtree_->Branch("wt_em_qcd",         &wt_em_qcd_);
        if(!systematic_shift_) {
          outtree_->Branch("wt_em_qcd_extrapup",      &wt_em_qcd_extrapup_);
          outtree_->Branch("wt_em_qcd_extrapdown",    &wt_em_qcd_extrapdown_);
          outtree_->Branch("wt_em_qcd_njets0_unc1_up",    &wt_em_qcd_njets0_unc1_up_);
          outtree_->Branch("wt_em_qcd_njets0_unc1_down",  &wt_em_qcd_njets0_unc1_down_);
          outtree_->Branch("wt_em_qcd_njets0_unc2_up",    &wt_em_qcd_njets0_unc2_up_);
          outtree_->Branch("wt_em_qcd_njets0_unc2_down",  &wt_em_qcd_njets0_unc2_down_);
          outtree_->Branch("wt_em_qcd_njets1_unc1_up",    &wt_em_qcd_njets1_unc1_up_);
          outtree_->Branch("wt_em_qcd_njets1_unc1_down",  &wt_em_qcd_njets1_unc1_down_);
          outtree_->Branch("wt_em_qcd_njets1_unc2_up",    &wt_em_qcd_njets1_unc2_up_);
          outtree_->Branch("wt_em_qcd_njets1_unc2_down",  &wt_em_qcd_njets1_unc2_down_);
          outtree_->Branch("wt_em_qcd_njets2_unc1_up",    &wt_em_qcd_njets2_unc1_up_);
          outtree_->Branch("wt_em_qcd_njets2_unc1_down",  &wt_em_qcd_njets2_unc1_down_);
          outtree_->Branch("wt_em_qcd_njets2_unc2_up",    &wt_em_qcd_njets2_unc2_up_);
          outtree_->Branch("wt_em_qcd_njets2_unc2_down",  &wt_em_qcd_njets2_unc2_down_);
        }
      }
 
      outtree_->Branch("jet_flav_1", &jet_flav_1_);
      outtree_->Branch("jet_flav_2", &jet_flav_2_);
      outtree_->Branch("jet_flav_3", &jet_flav_3_);
      outtree_->Branch("jet_pt_1", &jet_pt_1_);
      outtree_->Branch("jet_pt_2", &jet_pt_2_);
      if(do_faketaus_){
        outtree_->Branch("tau_dm_1", &tau_dm_1_);
        outtree_->Branch("tau_pt_1", &tau_pt_1_);
      }

      outtree_->Branch("q_1", &q_1_);
      outtree_->Branch("q_2", &q_2_);
      outtree_->Branch("wt_zpt_down",       &wt_zpt_down_);
      outtree_->Branch("wt_zpt_embed_ic",       &wt_zpt_embed_ic_);
        
      outtree_->Branch("wt_tquark_up",      &wt_tquark_up_);
      outtree_->Branch("wt_tquark_down",    &wt_tquark_down_);
      outtree_->Branch("wt_tquark_alt",      &wt_tquark_alt_);
                                                                
      //Variables needed for control plots need only be generated for central systematics
      if(!systematic_shift_) {
        outtree_->Branch("wt_tau_fake_up",    &wt_tau_fake_up_);
        outtree_->Branch("wt_tau_fake_down",  &wt_tau_fake_down_);
        outtree_->Branch("wt_zpt_up",         &wt_zpt_up_);
        outtree_->Branch("wt_zpt_up_NLO",         &wt_zpt_up_NLO_);
        outtree_->Branch("wt_efake_rate_up", &wt_efake_rate_up_);
        outtree_->Branch("wt_efake_rate_down", &wt_efake_rate_down_);
        outtree_->Branch("wt_mfake_rate_up", &wt_mfake_rate_up_);
        outtree_->Branch("wt_mfake_rate_down", &wt_mfake_rate_down_);
        outtree_->Branch("n_vtx",             &n_vtx_);
        outtree_->Branch("met_phi",           &met_phi_.var_double);
        outtree_->Branch("wt_qcdscale_up"   , &wt_qcdscale_up_    );
        outtree_->Branch("wt_qcdscale_down" , &wt_qcdscale_down_ );
        outtree_->Branch("metcov00", &metCov00_, "metCov00/F");
        outtree_->Branch("metcov01", &metCov01_, "metCov01/F");
        outtree_->Branch("metcov10", &metCov10_, "metCov10/F");
        outtree_->Branch("metcov11", &metCov11_, "metCov11/F");
        outtree_->Branch("wt_ps_isr_up" , &wt_ps_isr_up_ );
        outtree_->Branch("wt_ps_isr_down" , &wt_ps_isr_down_ );
        outtree_->Branch("wt_ps_fsr_up" , &wt_ps_fsr_up_ );
        outtree_->Branch("wt_ps_fsr_down" , &wt_ps_fsr_down_ );

        outtree_->Branch("d0_1", &d0_1_.var_float, "d0_1/F");
        outtree_->Branch("dZ_1", &dz_1_.var_float, "dz_1/F");
        outtree_->Branch("d0_2", &d0_2_.var_float, "d0_2/F");
        outtree_->Branch("dZ_2", &dz_2_.var_float, "dz_2/F");

        if (channel_ == channel::em) {
          outtree_->Branch("pzetavis",          &pzetavis_.var_double);
          outtree_->Branch("pzetamiss",         &pzetamiss_.var_double);
          outtree_->Branch("mt_ll",             &mt_ll_);
          outtree_->Branch("dz_1",              &dz_1_.var_double);
          outtree_->Branch("dz_2",              &dz_2_.var_double);
        } 
      }
    }
    if(make_sync_ntuple_) {
      lOFile = new TFile(sync_output_name_.c_str(), "RECREATE");
      lOFile->cd();
      synctree_ = new TTree("TauCheck", "TauCheck");

      synctree_->Branch("run", &run_, "run/I");
      synctree_->Branch("lumi", &lumi_, "lumi/I");
      synctree_->Branch("event", &event_, "event/l");

      synctree_->Branch("npv", &n_vtx_, "n_vtx/I");
      synctree_->Branch("npu", &n_pu_, "n_pu/F");
      synctree_->Branch("rho", &rho_, "rho/F");

      synctree_->Branch("puweight", &pu_weight_, "pu_weight/F");

     if(channel_==channel::et) synctree_->Branch("trigweight_1", &et_trg_);
     else synctree_->Branch("trigweight_1", &trigweight_1_, "trigweight_1/F");

      synctree_->Branch("trigweight_2", &et_trg_cross_);
      synctree_->Branch("trigweight", &et_trg_or_);
      synctree_->Branch("idisoweight_1", &idisoweight_1_);
      if(channel_==channel::et) synctree_->Branch("idisoweight_2", &mvadm_idiso_et_);
      else synctree_->Branch("idisoweight_2", &idisoweight_2_, "idisoweight_2/F");
      synctree_->Branch("embweight", &wt_embedding_yield_);
      synctree_->Branch("trkeffweight", &trackingweight_1_, "trackingweight_1/F");
      synctree_->Branch("trackingweight_1", &trackingweight_1_, "trackingweight_1/F");
      synctree_->Branch("trackingweight_2", &trackingweight_2_, "trackingweight_2/F");
      synctree_->Branch("effweight", &effweight_, "effweight/F");
      synctree_->Branch("weight", &wt_.var_float, "wt/F");

      synctree_->Branch("m_vis", &m_vis_.var_float, "m_vis/F");
      synctree_->Branch("mt_sv", &mt_sv_.var_float, "mt_sv/F");
      synctree_->Branch("mt_fast", &m_sv_.var_float, "mt_fast/F");
      synctree_->Branch("m_sv",   &m_sv_.var_float, "m_sv/F");
      synctree_->Branch("m_fast", &m_sv_.var_float, "m_fast/F");
      synctree_->Branch("pt_sv", &pt_h_.var_float, "pt_h/F");
      synctree_->Branch("eta_sv", &eta_h_, "eta_h/F");
      synctree_->Branch("phi_sv", &phi_h_, "phi_h/F");

      // Lepton 1 properties
      if ((strategy_ == strategy::cpsummer17 || strategy_ == strategy::cpdecays17 || strategy_ == strategy::cpdecays18) && (channel_ == channel::et || channel_ == channel::em)) { 
        synctree_->Branch("eCorr_1",           &E_1_);
      }
      synctree_->Branch("pt_1", &pt_1_.var_float, "pt_1/F");
      synctree_->Branch("raw_pt_1", &raw_pt_1_);
      synctree_->Branch("phi_1", &phi_1_.var_float, "phi_1/F");
      synctree_->Branch("eta_1", &eta_1_.var_float, "eta_1/F");
      synctree_->Branch("m_1", &m_1_, "m_1/F");
      synctree_->Branch("q_1", &q_1_, "q_1/I");
      synctree_->Branch("iso_1", &iso_1_.var_float, "iso_1/F");
      synctree_->Branch("iso_2", &iso_2_.var_float, "iso_2/F");
      synctree_->Branch("d0_1", &d0_1_.var_float, "d0_1/F");
      synctree_->Branch("dZ_1", &dz_1_.var_float, "dz_1/F");
      synctree_->Branch("mt_1", &mt_1_.var_float, "mt_1/F");

      // repeated twice in case want different branch names for sync
      synctree_->Branch("tau_decay_mode_1", & tau_decay_mode_1_, "tau_decay_mode_1/I");
      synctree_->Branch("tau_decay_mode_2", & tau_decay_mode_2_, "tau_decay_mode_2/I");
      synctree_->Branch("dm_1", & tau_decay_mode_1_, "dm_1/I");
      synctree_->Branch("dm_2", & tau_decay_mode_2_, "dm_2/I");
      // MVA tau decay modes (for decays analysis specifically)
      synctree_->Branch("mva_dm_1",         & tau_mva_decay_mode_1_, "mva_dm_1/I");
      synctree_->Branch("mva_dm_2",         & tau_mva_decay_mode_2_, "mva_dm_2/I");
      synctree_->Branch("dmMVA_1",         & tau_mva_decay_mode_1_, "dmMVA_1/I");
      synctree_->Branch("dmMVA_2",         & tau_mva_decay_mode_2_, "dmMVA_2/I");

      // deeptau variables
      synctree_->Branch("deepTauVsJetsRaw_1", & iso_1_, "deepTauVsJetsRaw_1/D");
      synctree_->Branch("deepTauVsJetsRaw_2", & iso_2_, "deepTauVsJetsRaw_2/D");
      synctree_->Branch("deepTauVsEleRaw_1",  & deepTauVsEle_iso_1_, "deepTauVsEleRaw_1/D");
      synctree_->Branch("deepTauVsEleRaw_2",  & deepTauVsEle_iso_2_, "deepTauVsEleRaw_2/D");
      synctree_->Branch("deepTauVsMuRaw_1",   & deepTauVsMu_iso_1_, "deepTauVsMuRaw_1/D");
      synctree_->Branch("deepTauVsMuRaw_2",   & deepTauVsMu_iso_2_, "deepTauVsMuRaw_2/D");

      // tauspinner weights
      synctree_->Branch("tauspinnerH",      & wt_cp_sm_, "tauspinnerH/D");
      synctree_->Branch("tauspinnerA",      & wt_cp_ps_, "tauspinnerA/D");
      synctree_->Branch("tauspinnerMaxMix", & wt_cp_mm_, "tauspinnerMaxMix/D");

      synctree_->Branch("pv_angle", &pv_angle_);
      // acoplanarity angles
      if (channel_ == channel::mt || channel_ == channel::et ) {
        // 00: mu + pi (IP)
        // 01: mu + rho/a1 (DP)
        synctree_->Branch("acotautau_00",      & aco_angle_6_, "acotautau_00/D");
        synctree_->Branch("acotautau_01",      & aco_angle_5_, "acotautau_01/D");
        synctree_->Branch("gen_match_2_pt",      &gen_match_2_pt_);
      }
      else if (channel_ == channel::tt) {
        // 00: rho + pi (IP)
        // 01: rho + rho (DP)
        synctree_->Branch("acotautau_00",      & aco_angle_6_, "acotautau_00/D");
        synctree_->Branch("acotautau_01",      & aco_angle_1_, "acotautau_01/D");
        synctree_->Branch("gen_match_1_pt",      &gen_match_1_pt_);
        synctree_->Branch("gen_match_2_pt",      &gen_match_2_pt_);

      }

      // IP variables
      synctree_->Branch("ip_sig_1", &ip_sig_1_);
      synctree_->Branch("ip_sig_2", &ip_sig_2_);
      synctree_->Branch("ip_mag_1", &ip_mag_1_);
      synctree_->Branch("ip_mag_2", &ip_mag_2_);
      synctree_->Branch("ipx_1", &ipx_1_);
      synctree_->Branch("ipy_1", &ipy_1_);
      synctree_->Branch("ipz_1", &ipz_1_);
      synctree_->Branch("ipx_2", &ipx_2_);
      synctree_->Branch("ipy_2", &ipy_2_);
      synctree_->Branch("ipz_2", &ipz_2_);


      synctree_->Branch("ipx_uncorr_1", &ipx_uncorr_1_);
      synctree_->Branch("ipy_uncorr_1", &ipy_uncorr_1_);
      synctree_->Branch("ipz_uncorr_1", &ipz_uncorr_1_);
      synctree_->Branch("ipx_uncorr_2", &ipx_uncorr_2_);
      synctree_->Branch("ipy_uncorr_2", &ipy_uncorr_2_);
      synctree_->Branch("ipz_uncorr_2", &ipz_uncorr_2_);
      synctree_->Branch("ip_sig_uncorr_1", &ip_sig_uncorr_1_);
      synctree_->Branch("ip_sig_uncorr_2", &ip_sig_uncorr_2_);

      synctree_->Branch("gen_ipx_1", &gen_ipx_1_);
      synctree_->Branch("gen_ipy_1", &gen_ipy_1_);
      synctree_->Branch("gen_ipz_1", &gen_ipz_1_);
      synctree_->Branch("gen_ipx_2", &gen_ipx_2_);
      synctree_->Branch("gen_ipy_2", &gen_ipy_2_);
      synctree_->Branch("gen_ipz_2", &gen_ipz_2_);

      // refitted vertex parameters
      synctree_->Branch("pvx", &pvx_);
      synctree_->Branch("pvy", &pvy_);
      synctree_->Branch("pvz", &pvz_);

      // bool for refitted vertex
      synctree_->Branch("use_refitted_vertex", & use_refitted_vertex_);

      synctree_->Branch("tDeepTau2017v2p1VSjetraw_1",      &iso_1_);
      synctree_->Branch("tDeepTau2017v2p1VSjetraw_2",      &iso_2_);
      synctree_->Branch("tDeepTau2017v2p1VSeraw_1",       &deepTauVsEle_iso_1_);
      synctree_->Branch("tDeepTau2017v2p1VSeraw_2",       &deepTauVsEle_iso_2_);
      synctree_->Branch("tDeepTau2017v2p1VSmuraw_1",        &deepTauVsMu_iso_1_);
      synctree_->Branch("tDeepTau2017v2p1VSmuraw_2",        &deepTauVsMu_iso_2_);
      synctree_->Branch("tVVVLooseDeepTau2017v2p1VSjet_1", &deepTauVsJets_vvvloose_1_);
      synctree_->Branch("tVVVLooseDeepTau2017v2p1VSjet_2", &deepTauVsJets_vvvloose_2_);
      synctree_->Branch("tVVLooseDeepTau2017v2p1VSjet_1",  &deepTauVsJets_vvloose_1_);
      synctree_->Branch("tVVLooseDeepTau2017v2p1VSjet_2",  &deepTauVsJets_vvloose_2_);
      synctree_->Branch("tVLooseDeepTau2017v2p1VSjet_1",   &deepTauVsJets_vloose_1_);
      synctree_->Branch("tVLooseDeepTau2017v2p1VSjet_2",   &deepTauVsJets_vloose_2_);
      synctree_->Branch("tLooseDeepTau2017v2p1VSjet_1",    &deepTauVsJets_loose_1_);
      synctree_->Branch("tLooseDeepTau2017v2p1VSjet_2",    &deepTauVsJets_loose_2_);
      synctree_->Branch("tMediumDeepTau2017v2p1VSjet_1",   &deepTauVsJets_medium_1_);
      synctree_->Branch("tMediumDeepTau2017v2p1VSjet_2",   &deepTauVsJets_medium_2_);
      synctree_->Branch("tTightDeepTau2017v2p1VSjet_1",    &deepTauVsJets_tight_1_);
      synctree_->Branch("tTightDeepTau2017v2p1VSjet_2",    &deepTauVsJets_tight_2_);
      synctree_->Branch("tVTightDeepTau2017v2p1VSjet_1",   &deepTauVsJets_vtight_1_);
      synctree_->Branch("tVTightDeepTau2017v2p1VSjet_2",   &deepTauVsJets_vtight_2_);
      synctree_->Branch("tVVTightDeepTau2017v2p1VSjet_1",  &deepTauVsJets_vvtight_1_);
      synctree_->Branch("tVVTightDeepTau2017v2p1VSjet_2",  &deepTauVsJets_vvtight_2_);
      synctree_->Branch("tVVVLooseDeepTau2017v2p1VSe_1",  &deepTauVsEle_vvvloose_1_);
      synctree_->Branch("tVVVLooseDeepTau2017v2p1VSe_2",  &deepTauVsEle_vvvloose_2_);
      synctree_->Branch("tVVLooseDeepTau2017v2p1VSe_1",   &deepTauVsEle_vvloose_1_);
      synctree_->Branch("tVVLooseDeepTau2017v2p1VSe_2",   &deepTauVsEle_vvloose_2_);
      synctree_->Branch("tVLooseDeepTau2017v2p1VSe_1",    &deepTauVsEle_vloose_1_);
      synctree_->Branch("tVLooseDeepTau2017v2p1VSe_2",    &deepTauVsEle_vloose_2_);
      synctree_->Branch("tLooseDeepTau2017v2p1VSe_1",     &deepTauVsEle_loose_1_);
      synctree_->Branch("tLooseDeepTau2017v2p1VSe_2",     &deepTauVsEle_loose_2_);
      synctree_->Branch("tMediumDeepTau2017v2p1VSe_1",    &deepTauVsEle_medium_1_);
      synctree_->Branch("tMediumDeepTau2017v2p1VSe_2",    &deepTauVsEle_medium_2_);
      synctree_->Branch("tTightDeepTau2017v2p1VSe_1",     &deepTauVsEle_tight_1_);
      synctree_->Branch("tTightDeepTau2017v2p1VSe_2",     &deepTauVsEle_tight_2_);
      synctree_->Branch("tVTightDeepTau2017v2p1VSe_1",    &deepTauVsEle_vtight_1_);
      synctree_->Branch("tVTightDeepTau2017v2p1VSe_2",    &deepTauVsEle_vtight_2_);
      synctree_->Branch("tVVTightDeepTau2017v2p1VSe_1",   &deepTauVsEle_vvtight_1_);
      synctree_->Branch("tVVTightDeepTau2017v2p1VSe_2",   &deepTauVsEle_vvtight_2_);
      synctree_->Branch("tVVVLooseDeepTau2017v2p1VSmu_1",  &deepTauVsMu_vvvloose_1_);
      synctree_->Branch("tVVVLooseDeepTau2017v2p1VSmu_2",  &deepTauVsMu_vvvloose_2_);
      synctree_->Branch("tVVLooseDeepTau2017v2p1VSmu_1",   &deepTauVsMu_vvloose_1_);
      synctree_->Branch("tVVLooseDeepTau2017v2p1VSmu_2",   &deepTauVsMu_vvloose_2_);
      synctree_->Branch("tVLooseDeepTau2017v2p1VSmu_1",    &deepTauVsMu_vloose_1_);
      synctree_->Branch("tVLooseDeepTau2017v2p1VSmu_2",    &deepTauVsMu_vloose_2_);
      synctree_->Branch("tLooseDeepTau2017v2p1VSmu_1",     &deepTauVsMu_loose_1_);
      synctree_->Branch("tLooseDeepTau2017v2p1VSmu_2",     &deepTauVsMu_loose_2_);
      synctree_->Branch("tMediumDeepTau2017v2p1VSmu_1",    &deepTauVsMu_medium_1_);
      synctree_->Branch("tMediumDeepTau2017v2p1VSmu_2",    &deepTauVsMu_medium_2_);
      synctree_->Branch("tTightDeepTau2017v2p1VSmu_1",    &deepTauVsMu_tight_1_);
      synctree_->Branch("tTightDeepTau2017v2p1VSmu_2",    &deepTauVsMu_tight_2_);
      synctree_->Branch("tVTightDeepTau2017v2p1VSmu_1",    &deepTauVsMu_vtight_1_);
      synctree_->Branch("tVTightDeepTau2017v2p1VSmu_2",    &deepTauVsMu_vtight_2_);
      synctree_->Branch("tVVTightDeepTau2017v2p1VSmu_1",    &deepTauVsMu_vvtight_1_);
      synctree_->Branch("tVVTightDeepTau2017v2p1VSmu_2",    &deepTauVsMu_vvtight_2_);

      // DeepTau v2p5
      synctree_->Branch("tVVVLooseDeepTau2018v2p5VSjet_1", &deepTauVsJetsV2p5_vvvloose_1_);
      synctree_->Branch("tVVVLooseDeepTau2018v2p5VSjet_2", &deepTauVsJetsV2p5_vvvloose_2_);
      synctree_->Branch("tVVLooseDeepTau2018v2p5VSjet_1",  &deepTauVsJetsV2p5_vvloose_1_);
      synctree_->Branch("tVVLooseDeepTau2018v2p5VSjet_2",  &deepTauVsJetsV2p5_vvloose_2_);
      synctree_->Branch("tVLooseDeepTau2018v2p5VSjet_1",   &deepTauVsJetsV2p5_vloose_1_);
      synctree_->Branch("tVLooseDeepTau2018v2p5VSjet_2",   &deepTauVsJetsV2p5_vloose_2_);
      synctree_->Branch("tLooseDeepTau2018v2p5VSjet_1",    &deepTauVsJetsV2p5_loose_1_);
      synctree_->Branch("tLooseDeepTau2018v2p5VSjet_2",    &deepTauVsJetsV2p5_loose_2_);
      synctree_->Branch("tMediumDeepTau2018v2p5VSjet_1",   &deepTauVsJetsV2p5_medium_1_);
      synctree_->Branch("tMediumDeepTau2018v2p5VSjet_2",   &deepTauVsJetsV2p5_medium_2_);
      synctree_->Branch("tTightDeepTau2018v2p5VSjet_1",    &deepTauVsJetsV2p5_tight_1_);
      synctree_->Branch("tTightDeepTau2018v2p5VSjet_2",    &deepTauVsJetsV2p5_tight_2_);
      synctree_->Branch("tVTightDeepTau2018v2p5VSjet_1",   &deepTauVsJetsV2p5_vtight_1_);
      synctree_->Branch("tVTightDeepTau2018v2p5VSjet_2",   &deepTauVsJetsV2p5_vtight_2_);
      synctree_->Branch("tVVTightDeepTau2018v2p5VSjet_1",  &deepTauVsJetsV2p5_vvtight_1_);
      synctree_->Branch("tVVTightDeepTau2018v2p5VSjet_2",  &deepTauVsJetsV2p5_vvtight_2_);
      synctree_->Branch("tVVLooseDeepTau2018v2p5VSe_1",   &deepTauVsEleV2p5_vvloose_1_);
      synctree_->Branch("tVVLooseDeepTau2018v2p5VSe_2",   &deepTauVsEleV2p5_vvloose_2_);
      synctree_->Branch("tTightDeepTau2018v2p5VSe_1",     &deepTauVsEleV2p5_tight_1_);
      synctree_->Branch("tTightDeepTau2018v2p5VSe_2",     &deepTauVsEleV2p5_tight_2_);
      synctree_->Branch("tTightDeepTau2018v2p5VSmu_1",      &deepTauVsMuV2p5_tight_1_);
      synctree_->Branch("tTightDeepTau2018v2p5VSmu_2",      &deepTauVsMuV2p5_tight_2_);
      synctree_->Branch("tVLooseDeepTau2018v2p5VSmu_1",     &deepTauVsMuV2p5_vloose_1_);
      synctree_->Branch("tVLooseDeepTau2018v2p5VSmu_2",     &deepTauVsMuV2p5_vloose_2_);
      synctree_->Branch("tDeepTau2018v2p5VSjetraw_1",    &iso_1_V2p5_);
      synctree_->Branch("tDeepTau2018v2p5VSjetraw_2",    &iso_2_V2p5_);


      synctree_->Branch("pt_2", &pt_2_.var_float, "pt_2/F");
      synctree_->Branch("phi_2", &phi_2_.var_float, "phi_2/F");
      synctree_->Branch("eta_2", &eta_2_.var_float, "eta_2/F");
      synctree_->Branch("m_2", &m_2_, "lM2/F");
      synctree_->Branch("q_2", &q_2_, "lq2/I");
      synctree_->Branch("iso_2", &iso_2_);
      synctree_->Branch("d0_2", &d0_2_.var_float, "d0_2/F");
      synctree_->Branch("dZ_2", &dz_2_.var_float, "dz_2/F");
      synctree_->Branch("mt_2", &mt_2_.var_float, "mt_2/F");

      synctree_->Branch("os", &os_, "os/O");
      synctree_->Branch("dilepton_veto", &dilepton_veto_, "dilepton_veto/O");
      synctree_->Branch("extraelec_veto", &extraelec_veto_, "extraelec_veto/O");
      synctree_->Branch("extramuon_veto", &extramuon_veto_, "extramuon_veto/O");
      synctree_->Branch("gen_match_1", &gen_match_1_, "gen_match_1/i");
      synctree_->Branch("gen_match_2", &gen_match_2_,"gen_match_2/i");

      synctree_->Branch("met", &met_.var_float, "met/F");
      synctree_->Branch("metphi", &met_phi_.var_float, "met_phi/F");
      synctree_->Branch("uncorrmet", &uncorrmet_.var_float, "met/F");
      synctree_->Branch("uncorrmetphi", &uncorrmet_phi_.var_float, "met_phi/F");
      synctree_->Branch("metcov00", &metCov00_, "metCov00/F");
      synctree_->Branch("metcov01", &metCov01_, "metCov01/F");
      synctree_->Branch("metcov10", &metCov10_, "metCov10/F");
      synctree_->Branch("metcov11", &metCov11_, "metCov11/F");
      synctree_->Branch("genpX", &gen_px_, "genpX/F");
      synctree_->Branch("genpY", &gen_py_, "genpY/F");
      synctree_->Branch("vispX", &vis_px_, "vispX/F");
      synctree_->Branch("vispY", &vis_py_, "vispY/F");

      synctree_->Branch("pt_tt", &pt_tt_.var_float, "pt_tt/F");
      synctree_->Branch("mt_tot", &mt_tot_.var_float, "mt_tot/F");

      synctree_->Branch("pzetavis", &pzetavis_.var_float, "pzetavis/F");
      synctree_->Branch("pzetamiss", &pzetamiss_.var_float, "pzetamiss/F");

      synctree_->Branch("njets", &n_jets_, "n_jets/I");
      synctree_->Branch("njetspt20", &n_lowpt_jets_, "n_lowpt_jets/I");
      synctree_->Branch("jpt_1", &jpt_1_.var_float, "jpt_1/F");
      synctree_->Branch("jeta_1", &jeta_1_.var_float, "jeta_1/F");
      synctree_->Branch("jphi_1", &jphi_1_, "jphi_1/F");
      synctree_->Branch("jmva_1", &jmva_1_, "jmva_1/F");
      synctree_->Branch("jpt_2", &jpt_2_.var_float, "jpt_2/F");
      synctree_->Branch("jeta_2", &jeta_2_.var_float, "jeta_2/F");
      synctree_->Branch("jphi_2", &jphi_2_, "jphi_2/F");
      synctree_->Branch("jmva_2", &jmva_2_, "jmva_2/F");
      synctree_->Branch("mjj", &mjj_.var_float, "mjj/F");
      synctree_->Branch("jdeta", &jdeta_.var_float, "jdeta/F");
      synctree_->Branch("jdphi", &jdphi_, "jdphi/F");

      // Number of b-tagging jets passing above selections
      //synctree_->Branch("nbtag", &n_bjets_, "n_bjets/I");
      synctree_->Branch("nbtag", &n_deepbjets_);
      synctree_->Branch("nloosebtag",     &n_loose_bjets_);
      synctree_->Branch("bpt_1", &bpt_1_.var_float, "bpt_1/F");
      synctree_->Branch("beta_1", &beta_1_.var_float, "beta_1/F");
      synctree_->Branch("bphi_1", &bphi_1_, "bphi_1/F");
      synctree_->Branch("bmva_1", &bmva_1_, "bmva_1/F");
      synctree_->Branch("bcsv_1", &bcsv_1_.var_float, "bcsv_1/F");
      synctree_->Branch("bpt_2", &bpt_2_.var_float, "bpt_2/F");
      synctree_->Branch("beta_2", &beta_2_.var_float, "beta_2/F");
      synctree_->Branch("bphi_2", &bphi_2_, "bphi_2/F");
      synctree_->Branch("bmva_2", &bmva_2_, "bmva_2/F");
      synctree_->Branch("bcsv_2", &bcsv_2_, "bcsv_2/F");
      
      synctree_->Branch("trg_singleelectron",    &trg_singleelectron_);
      synctree_->Branch("trg_singlemuon",    &trg_singlemuon_);
      synctree_->Branch("trg_doubletau",    &trg_doubletau_);

      synctree_->Branch("trg_tt_monitoring_1",    &trg_tt_monitoring_1_);
      synctree_->Branch("trg_tt_monitoring_2",    &trg_tt_monitoring_2_);
      synctree_->Branch("trg_tt_monitoring_3",    &trg_tt_monitoring_3_);
      synctree_->Branch("trg_tt_monitoring_4",    &trg_tt_monitoring_4_);
      synctree_->Branch("trg_tt_monitoring_5",    &trg_tt_monitoring_5_);
      synctree_->Branch("trg_tt_monitoring_6",    &trg_tt_monitoring_6_);

      synctree_->Branch("trg_muonelectron",    &trg_muonelectron_);
      synctree_->Branch("trg_singletau_1",    &trg_singletau_1_);
      synctree_->Branch("trg_singletau_2",    &trg_singletau_2_);
      synctree_->Branch("trg_mutaucross", &trg_mutaucross_);
      synctree_->Branch("trg_etaucross", &trg_etaucross_);
      
      synctree_->Branch("generatorWeight", &generatorWeight_);
      synctree_->Branch("muonEffTrgWeight", &muonEffTrgWeight_);
      synctree_->Branch("muonEffIDWeight_1", &muonEffIDWeight_1_);
      synctree_->Branch("muonEffIDWeight_2", &muonEffIDWeight_2_);
      synctree_->Branch("et_triggerweight_ic", &et_triggerweight_ic_);
      synctree_->Branch("mt_triggerweight_ic", &mt_triggerweight_ic_);
      synctree_->Branch("idisoWeight_1", &idisoweight_1_);
      synctree_->Branch("tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2", &tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2_);
      synctree_->Branch("tauIDScaleFactorWeight_highpt_deeptauid_2", &tauIDScaleFactorWeight_highpt_deeptauid_2_);

      synctree_->Branch("wt_ggh_t", &wt_ggh_t_);
      synctree_->Branch("wt_ggh_b", &wt_ggh_b_);
      synctree_->Branch("wt_ggh_i", &wt_ggh_i_);
      synctree_->Branch("wt_ggH_t", &wt_ggH_t_);
      synctree_->Branch("wt_ggH_b", &wt_ggH_b_);
      synctree_->Branch("wt_ggH_i", &wt_ggH_i_);
      synctree_->Branch("wt_ggA_t", &wt_ggA_t_);
      synctree_->Branch("wt_ggA_b", &wt_ggA_b_);
      synctree_->Branch("wt_ggA_i", &wt_ggA_i_);

    }

    if(make_mva_ntuple_) {
      lOFile = new TFile(mva_output_name_.c_str(), "RECREATE");
      lOFile->cd();
      mvatree_ = new TTree("ntuple", "ntuple");
      mvatree_->Branch("genM", &gen_m_);
      mvatree_->Branch("genpT", &gen_pt_);
      mvatree_->Branch("gen_nu_p_1",   &gen_nu_p_1_);
      mvatree_->Branch("gen_nu_p_2",   &gen_nu_p_2_);
      mvatree_->Branch("gen_nu_p_1",   &gen_nu_p_1_);
      mvatree_->Branch("gen_nu_p_2",   &gen_nu_p_2_);
      mvatree_->Branch("gen_nu_px_1",   &gen_nu_px_1_);
      mvatree_->Branch("gen_nu_px_2",   &gen_nu_px_2_);
      mvatree_->Branch("gen_nu_py_1",   &gen_nu_py_1_);
      mvatree_->Branch("gen_nu_py_2",   &gen_nu_py_2_);
      mvatree_->Branch("gen_nu_pz_1",   &gen_nu_pz_1_);
      mvatree_->Branch("gen_nu_pz_2",   &gen_nu_pz_2_);
      mvatree_->Branch("gen_nu_phi_1", &gen_nu_phi_1_);
      mvatree_->Branch("gen_nu_phi_2", &gen_nu_phi_2_);
      mvatree_->Branch("gen_nu_eta_1", &gen_nu_eta_1_);
      mvatree_->Branch("gen_nu_eta_2", &gen_nu_eta_2_);

      mvatree_->Branch("gen_vis_p_1",   &gen_vis_p_1_);
      mvatree_->Branch("gen_vis_p_2",   &gen_vis_p_2_);
      mvatree_->Branch("gen_vis_px_1",  &gen_vis_px_1_);
      mvatree_->Branch("gen_vis_px_2",  &gen_vis_px_2_);
      mvatree_->Branch("gen_vis_py_1",  &gen_vis_py_1_);
      mvatree_->Branch("gen_vis_py_2",  &gen_vis_py_2_);
      mvatree_->Branch("gen_vis_pz_1",  &gen_vis_pz_1_);
      mvatree_->Branch("gen_vis_pz_2",  &gen_vis_pz_2_);
      mvatree_->Branch("gen_vis_E_1",   &gen_vis_E_1_);
      mvatree_->Branch("gen_vis_E_2",   &gen_vis_E_2_);
      mvatree_->Branch("gen_vis_phi_1", &gen_vis_phi_1_);
      mvatree_->Branch("gen_vis_phi_2", &gen_vis_phi_2_);
      mvatree_->Branch("gen_vis_eta_1", &gen_vis_eta_1_);
      mvatree_->Branch("gen_vis_eta_2", &gen_vis_eta_2_);

      mvatree_->Branch("gen_neutral_p_1",   &gen_neutral_p_1_);
      mvatree_->Branch("gen_neutral_p_2",   &gen_neutral_p_2_);
      mvatree_->Branch("gen_neutral_px_1",   &gen_neutral_px_1_);
      mvatree_->Branch("gen_neutral_px_2",   &gen_neutral_px_2_);
      mvatree_->Branch("gen_neutral_py_1",   &gen_neutral_py_1_);
      mvatree_->Branch("gen_neutral_py_2",   &gen_neutral_py_2_);
      mvatree_->Branch("gen_neutral_pz_1",   &gen_neutral_pz_1_);
      mvatree_->Branch("gen_neutral_pz_2",   &gen_neutral_pz_2_);
      mvatree_->Branch("gen_neutral_E_1",   &gen_neutral_E_1_);
      mvatree_->Branch("gen_neutral_E_2",   &gen_neutral_E_2_);
      mvatree_->Branch("gen_neutral_phi_1", &gen_neutral_phi_1_);
      mvatree_->Branch("gen_neutral_phi_2", &gen_neutral_phi_2_);
      mvatree_->Branch("gen_neutral_eta_1", &gen_neutral_eta_1_);
      mvatree_->Branch("gen_neutral_eta_2", &gen_neutral_eta_2_);

      mvatree_->Branch("os",                &os_);
      mvatree_->Branch("q1",                &q1_);
      mvatree_->Branch("q2",                &q2_);
      mvatree_->Branch("pt_vis",            &pt_vis_.var_double);
      mvatree_->Branch("pt_tt",             &pt_tt_.var_double);
      mvatree_->Branch("m_vis",             &m_vis_.var_double);

      mvatree_->Branch("pt_2",              &pt_2_.var_double);
      mvatree_->Branch("pt_1",              &pt_1_.var_double);
      mvatree_->Branch("iso_1", &iso_1_.var_float, "iso_1/F");
      mvatree_->Branch("met",               &met_.var_double);
      mvatree_->Branch("wt_cp_sm", &wt_cp_sm_);
      mvatree_->Branch("wt_cp_ps", &wt_cp_ps_);
      mvatree_->Branch("wt_cp_mm", &wt_cp_mm_);
      mvatree_->Branch("mva_dm_1", &tau_mva_decay_mode_1_);
      mvatree_->Branch("mva_dm_2", &tau_mva_decay_mode_2_);
      mvatree_->Branch("tau_decay_mode_2",    &tau_decay_mode_2_);
      mvatree_->Branch("tau_decay_mode_1",    &tau_decay_mode_1_);
      mvatree_->Branch("trg_singleelectron",    &trg_singleelectron_);
      mvatree_->Branch("trg_singlemuon",    &trg_singlemuon_);
      mvatree_->Branch("trg_doubletau",    &trg_doubletau_);
      mvatree_->Branch("deepTauVsJets_medium_1",   &deepTauVsJets_medium_1_);
      mvatree_->Branch("deepTauVsJets_medium_2",   &deepTauVsJets_medium_2_);
      mvatree_->Branch("deepTauVsEle_vvloose_1",   &deepTauVsEle_vvloose_1_);
      mvatree_->Branch("deepTauVsEle_vvloose_2",   &deepTauVsEle_vvloose_2_);
      mvatree_->Branch("deepTauVsEle_tight_1",   &deepTauVsEle_tight_1_);
      mvatree_->Branch("deepTauVsEle_tight_2",   &deepTauVsEle_tight_2_);
      mvatree_->Branch("deepTauVsMu_vloose_1",   &deepTauVsMu_vloose_1_);
      mvatree_->Branch("deepTauVsMu_vloose_2",   &deepTauVsMu_vloose_2_);
      mvatree_->Branch("deepTauVsMu_tight_1",   &deepTauVsMu_tight_1_);
      mvatree_->Branch("deepTauVsMu_tight_2",   &deepTauVsMu_tight_2_);
      mvatree_->Branch("aco_angle_1", &aco_angle_1_);
      mvatree_->Branch("aco_angle_2", &aco_angle_2_);
      mvatree_->Branch("aco_angle_3", &aco_angle_3_);
      mvatree_->Branch("aco_angle_4", &aco_angle_4_);
      mvatree_->Branch("aco_angle_5", &aco_angle_5_);
      mvatree_->Branch("aco_angle_6", &aco_angle_6_);
      mvatree_->Branch("aco_angle_7", &aco_angle_7_);
      mvatree_->Branch("pv_angle", &pv_angle_);
      mvatree_->Branch("mass0",         &mass0_);
      mvatree_->Branch("mass1",         &mass1_);
      mvatree_->Branch("mass2",         &mass2_);
      mvatree_->Branch("ip_sig_1", &ip_sig_1_);
      mvatree_->Branch("ip_sig_2", &ip_sig_2_);
      mvatree_->Branch("rand", &rand_);
      mvatree_->Branch("tauFlag_1", &tauFlag_1_);
      mvatree_->Branch("tauFlag_2", &tauFlag_2_);
      mvatree_->Branch("gen_phitt", &gen_phitt_);
      mvatree_->Branch("y_1_1", &y_1_1_);
      mvatree_->Branch("y_1_2", &y_1_2_);
      mvatree_->Branch("y_2_2", &y_2_2_);
      mvatree_->Branch("y_3_2", &y_3_2_);
      mvatree_->Branch("y_4_2", &y_4_2_);
      mvatree_->Branch("metcov00", &metCov00_, "metCov00/F");
      mvatree_->Branch("metcov01", &metCov01_, "metCov01/F");
      mvatree_->Branch("metcov10", &metCov10_, "metCov10/F");
      mvatree_->Branch("metcov11", &metCov11_, "metCov11/F");
      mvatree_->Branch("pi_px_1", &pi_px_1_);
      mvatree_->Branch("pi_py_1", &pi_py_1_);
      mvatree_->Branch("pi_pz_1", &pi_pz_1_);
      mvatree_->Branch("pi_E_1", &pi_E_1_);
      mvatree_->Branch("pi2_px_1", &pi2_px_1_);
      mvatree_->Branch("pi2_py_1", &pi2_py_1_);
      mvatree_->Branch("pi2_pz_1", &pi2_pz_1_);
      mvatree_->Branch("pi2_E_1", &pi2_E_1_);
      mvatree_->Branch("pi3_px_1", &pi3_px_1_);
      mvatree_->Branch("pi3_py_1", &pi3_py_1_);
      mvatree_->Branch("pi3_pz_1", &pi3_pz_1_);
      mvatree_->Branch("pi3_E_1", &pi3_E_1_);
      mvatree_->Branch("pi_px_2", &pi_px_2_);
      mvatree_->Branch("pi_py_2", &pi_py_2_);
      mvatree_->Branch("pi_pz_2", &pi_pz_2_);
      mvatree_->Branch("pi_E_2", &pi_E_2_);
      mvatree_->Branch("pi2_px_2", &pi2_px_2_);
      mvatree_->Branch("pi2_py_2", &pi2_py_2_);
      mvatree_->Branch("pi2_pz_2", &pi2_pz_2_);
      mvatree_->Branch("pi2_E_2", &pi2_E_2_);
      mvatree_->Branch("pi3_px_2", &pi3_px_2_);
      mvatree_->Branch("pi3_py_2", &pi3_py_2_);
      mvatree_->Branch("pi3_pz_2", &pi3_pz_2_);
      mvatree_->Branch("pi3_E_2", &pi3_E_2_);
      mvatree_->Branch("pi0_px_1", &pi0_px_1_);
      mvatree_->Branch("pi0_py_1", &pi0_py_1_);
      mvatree_->Branch("pi0_pz_1", &pi0_pz_1_);
      mvatree_->Branch("pi0_E_1", &pi0_E_1_);
      mvatree_->Branch("pi0_px_2", &pi0_px_2_);
      mvatree_->Branch("pi0_py_2", &pi0_py_2_);
      mvatree_->Branch("pi0_pz_2", &pi0_pz_2_);
      mvatree_->Branch("pi0_E_2", &pi0_E_2_);
      mvatree_->Branch("gam1_px_1", &gam1_px_1_);
      mvatree_->Branch("gam1_py_1", &gam1_py_1_);
      mvatree_->Branch("gam1_pz_1", &gam1_pz_1_);
      mvatree_->Branch("gam1_E_1", &gam1_E_1_);
      mvatree_->Branch("gam1_px_2", &gam1_px_2_);
      mvatree_->Branch("gam1_py_2", &gam1_py_2_);
      mvatree_->Branch("gam1_pz_2", &gam1_pz_2_);
      mvatree_->Branch("gam1_E_2", &gam1_E_2_);
      mvatree_->Branch("gam2_px_1", &gam2_px_1_);
      mvatree_->Branch("gam2_py_1", &gam2_py_1_);
      mvatree_->Branch("gam2_pz_1", &gam2_pz_1_);
      mvatree_->Branch("gam2_E_1", &gam2_E_1_);
      mvatree_->Branch("gam2_px_2", &gam2_px_2_);
      mvatree_->Branch("gam2_py_2", &gam2_py_2_);
      mvatree_->Branch("gam2_pz_2", &gam2_pz_2_);
      mvatree_->Branch("gam2_E_2", &gam2_E_2_);
      mvatree_->Branch("gam3_px_1", &gam3_px_1_);
      mvatree_->Branch("gam3_py_1", &gam3_py_1_);
      mvatree_->Branch("gam3_pz_1", &gam3_pz_1_);
      mvatree_->Branch("gam3_E_1", &gam3_E_1_);
      mvatree_->Branch("gam3_px_2", &gam3_px_2_);
      mvatree_->Branch("gam3_py_2", &gam3_py_2_);
      mvatree_->Branch("gam3_pz_2", &gam3_pz_2_);
      mvatree_->Branch("gam3_E_2", &gam3_E_2_);
      mvatree_->Branch("gam4_px_1", &gam4_px_1_);
      mvatree_->Branch("gam4_py_1", &gam4_py_1_);
      mvatree_->Branch("gam4_pz_1", &gam4_pz_1_);
      mvatree_->Branch("gam4_E_1", &gam4_E_1_);
      mvatree_->Branch("gam4_px_2", &gam4_px_2_);
      mvatree_->Branch("gam4_py_2", &gam4_py_2_);
      mvatree_->Branch("gam4_pz_2", &gam4_pz_2_);
      mvatree_->Branch("gam4_E_2", &gam4_E_2_);
      mvatree_->Branch("metx", &metx_);
      mvatree_->Branch("mety", &mety_);
      mvatree_->Branch("gen_tau_metx", &gen_tau_metx_);
      mvatree_->Branch("gen_tau_mety", &gen_tau_mety_);
      mvatree_->Branch("gen_metx", &gen_metx_);
      mvatree_->Branch("gen_mety", &gen_mety_);
      mvatree_->Branch("ip_x_1", &ip_x_1_);
      mvatree_->Branch("ip_y_1", &ip_y_1_);
      mvatree_->Branch("ip_z_1", &ip_z_1_);
      mvatree_->Branch("ip_x_2", &ip_x_2_);
      mvatree_->Branch("ip_y_2", &ip_y_2_);
      mvatree_->Branch("ip_z_2", &ip_z_2_);
      mvatree_->Branch("sv_x_1", &sv_x_1_);
      mvatree_->Branch("sv_y_1", &sv_y_1_);
      mvatree_->Branch("sv_z_1", &sv_z_1_);
      mvatree_->Branch("sv_x_2", &sv_x_2_);
      mvatree_->Branch("sv_y_2", &sv_y_2_);
      mvatree_->Branch("sv_z_2", &sv_z_2_);
      mvatree_->Branch("ipcov00_1", &ipcov00_1_);
      mvatree_->Branch("ipcov01_1", &ipcov01_1_);
      mvatree_->Branch("ipcov02_1", &ipcov02_1_);
      mvatree_->Branch("ipcov10_1", &ipcov10_1_);
      mvatree_->Branch("ipcov11_1", &ipcov11_1_);
      mvatree_->Branch("ipcov12_1", &ipcov12_1_);
      mvatree_->Branch("ipcov20_1", &ipcov20_1_);
      mvatree_->Branch("ipcov21_1", &ipcov21_1_);
      mvatree_->Branch("ipcov22_1", &ipcov22_1_);
      mvatree_->Branch("ipcov00_2", &ipcov00_2_);
      mvatree_->Branch("ipcov01_2", &ipcov01_2_);
      mvatree_->Branch("ipcov02_2", &ipcov02_2_);
      mvatree_->Branch("ipcov10_2", &ipcov10_2_);
      mvatree_->Branch("ipcov11_2", &ipcov11_2_);
      mvatree_->Branch("ipcov12_2", &ipcov12_2_);
      mvatree_->Branch("ipcov20_2", &ipcov20_2_);
      mvatree_->Branch("ipcov21_2", &ipcov21_2_);
      mvatree_->Branch("ipcov22_2", &ipcov22_2_);
      mvatree_->Branch("svcov00_1", &svcov00_1_);
      mvatree_->Branch("svcov01_1", &svcov01_1_);
      mvatree_->Branch("svcov02_1", &svcov02_1_);
      mvatree_->Branch("svcov10_1", &svcov10_1_);
      mvatree_->Branch("svcov11_1", &svcov11_1_);
      mvatree_->Branch("svcov12_1", &svcov12_1_);
      mvatree_->Branch("svcov20_1", &svcov20_1_);
      mvatree_->Branch("svcov21_1", &svcov21_1_);
      mvatree_->Branch("svcov22_1", &svcov22_1_);
      mvatree_->Branch("svcov00_2", &svcov00_2_);
      mvatree_->Branch("svcov01_2", &svcov01_2_);
      mvatree_->Branch("svcov02_2", &svcov02_2_);
      mvatree_->Branch("svcov10_2", &svcov10_2_);
      mvatree_->Branch("svcov11_2", &svcov11_2_);
      mvatree_->Branch("svcov12_2", &svcov12_2_);
      mvatree_->Branch("svcov20_2", &svcov20_2_);
      mvatree_->Branch("svcov21_2", &svcov21_2_);
      mvatree_->Branch("svcov22_2", &svcov22_2_);
      mvatree_->Branch("pvx", &pvx_);
      mvatree_->Branch("pvy", &pvy_);
      mvatree_->Branch("pvz", &pvz_);
      mvatree_->Branch("q_pi_1", &q_pi_1_);
      mvatree_->Branch("q_pi2_1", &q_pi2_1_);
      mvatree_->Branch("q_pi3_1", &q_pi3_1_);
      mvatree_->Branch("q_pi_2", &q_pi_2_);
      mvatree_->Branch("q_pi2_2", &q_pi2_2_);
      mvatree_->Branch("q_pi3_2", &q_pi3_2_);	  
    }
    return 0;
  }

  int HTTCategories::Execute(TreeEvent *event) {
    
    if (event->Exists("trg_singleelectron")) trg_singleelectron_ = event->Get<bool>("trg_singleelectron");
    if (event->Exists("trg_singlemuon"))     trg_singlemuon_     = event->Get<bool>("trg_singlemuon");
    if (event->Exists("trg_singlemuon24"))     trg_singlemuon24_     = event->Get<bool>("trg_singlemuon24");
    else trg_singlemuon24_ = false;
    if (event->Exists("trg_singlemuon27"))     trg_singlemuon27_     = event->Get<bool>("trg_singlemuon27");
    else trg_singlemuon27_ = false;
    if (event->Exists("trg_doubletau"))      trg_doubletau_      = event->Get<bool>("trg_doubletau");

    if (event->Exists("trg_tt_monitoring_1"))trg_tt_monitoring_1_ = event->Get<bool>("trg_tt_monitoring_1");
    if (event->Exists("trg_tt_monitoring_2"))trg_tt_monitoring_2_ = event->Get<bool>("trg_tt_monitoring_2");
    if (event->Exists("trg_tt_monitoring_3"))trg_tt_monitoring_3_ = event->Get<bool>("trg_tt_monitoring_3");
    if (event->Exists("trg_tt_monitoring_4"))trg_tt_monitoring_4_ = event->Get<bool>("trg_tt_monitoring_4");
    if (event->Exists("trg_tt_monitoring_5"))trg_tt_monitoring_5_ = event->Get<bool>("trg_tt_monitoring_5");
    if (event->Exists("trg_tt_monitoring_6"))trg_tt_monitoring_6_ = event->Get<bool>("trg_tt_monitoring_6");
    if (event->Exists("trg_tt_monitoring_l1"))trg_tt_monitoring_l1_ = event->Get<float>("trg_tt_monitoring_l1");

    if (event->Exists("trg_doubletau_mssm"))      trg_doubletau_mssm_      = event->Get<bool>("trg_doubletau_mssm");
    if (event->Exists("trg_vbfdoubletau"))   trg_vbfdoubletau_   = event->Get<bool>("trg_vbfdoubletau");
    if (event->Exists("trg_muonelectron"))   trg_muonelectron_   = event->Get<bool>("trg_muonelectron");
    if (event->Exists("trg_singletau_1"))    trg_singletau_1_      = event->Get<bool>("trg_singletau_1");
    if (event->Exists("trg_singletau_2"))    trg_singletau_2_      = event->Get<bool>("trg_singletau_2");
    if (event->Exists("trg_mutaucross"))     trg_mutaucross_ = event->Get<bool>("trg_mutaucross");
    if (event->Exists("trg_etaucross"))     trg_etaucross_ = event->Get<bool>("trg_etaucross");
 
    if(!trg_applied_in_mc_) {
      trg_singleelectron_ = true;
      trg_singlemuon_     = true;
      trg_doubletau_      = true;
      trg_tt_monitoring_1_ = true;
      trg_tt_monitoring_2_ = true;
      trg_tt_monitoring_3_ = true;
      trg_tt_monitoring_4_ = true;
      trg_tt_monitoring_5_ = true;
      trg_tt_monitoring_6_ = true;
      trg_vbfdoubletau_   = true;
      trg_muonelectron_   = true;
      trg_singletau_1_      = true;
      trg_singletau_2_      = true; 
      trg_mutaucross_ = true;
      trg_etaucross_ = true;
    } 


    std::string MVADM2017;
    if(era_ == era::data_2016 || era_ == era::data_2017 || era_ == era::data_2018){ 
         MVADM2017 = "MVADM2017v1";

    } else if(era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP || era_ == era::data_2017UL || era_ == era::data_2018UL){ 
         MVADM2017 = "MVADM2017v2";
    } else {
         MVADM2017 = " ";
    }

    w_frac_score_ = event->Exists("w_frac_score") ? event->Get<double>("w_frac_score") : -1.; 
    qcd_frac_score_ = event->Exists("qcd_frac_score") ? event->Get<double>("qcd_frac_score") : -1.;

    w_frac_ = event->Exists("w_frac") ? event->Get<double>("w_frac") : -1.; 
    qcd_frac_ = event->Exists("qcd_frac") ? event->Get<double>("qcd_frac") : -1.; 
    ttbar_frac_ = event->Exists("ttbar_frac") ? event->Get<double>("ttbar_frac") : -1.; 

    // Get the objects we need from the event
    EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
    
    wt_ = {eventInfo->total_weight(), static_cast<float>(eventInfo->total_weight())};
    wt_wsoup_ = eventInfo->weight_defined("wsoup") ? eventInfo->weight("wsoup") : 1.0;
    wt_dysoup_ = eventInfo->weight_defined("dysoup") ? eventInfo->weight("dysoup") : 1.0;
    wt_dysoup_NLO_ = eventInfo->weight_defined("dysoup_NLO") ? eventInfo->weight("dysoup_NLO") : 1.0;

    //std::cout << (unsigned long long) eventInfo->event() << std::endl; 
    //eventInfo->print_weights();
    //eventInfo->print_all_weights();
   
    tau_idisoweight_1_ = event->Exists("tau_idisoweight_1") ? event->Get<double>("tau_idisoweight_1") : 1.0;
    tau_idisoweight_2_ = event->Exists("tau_idisoweight_2") ? event->Get<double>("tau_idisoweight_2") : 1.0;

    tau_idisoweight_ratio_1_up_ = event->Exists("tau_idisoweight_ratio_1_up") ? event->Get<double>("tau_idisoweight_ratio_1_up") : 1.0;
    tau_idisoweight_ratio_2_up_ = event->Exists("tau_idisoweight_ratio_2_up") ? event->Get<double>("tau_idisoweight_ratio_2_up") : 1.0;

    tau_idisoweight_ratio_1_down_ = event->Exists("tau_idisoweight_ratio_1_down") ? event->Get<double>("tau_idisoweight_ratio_1_down") : 1.0;
    tau_idisoweight_ratio_2_down_ = event->Exists("tau_idisoweight_ratio_2_down") ? event->Get<double>("tau_idisoweight_ratio_2_down") : 1.0;

    tau_idisoweight_ratio_1_uncert0_up_ = event->Exists("tau_idisoweight_ratio_1_uncert0_up") ? event->Get<double>("tau_idisoweight_ratio_1_uncert0_up") : 1.0;
    tau_idisoweight_ratio_2_uncert0_up_ = event->Exists("tau_idisoweight_ratio_2_uncert0_up") ? event->Get<double>("tau_idisoweight_ratio_2_uncert0_up") : 1.0;

    tau_idisoweight_ratio_1_uncert0_down_ = event->Exists("tau_idisoweight_ratio_1_uncert0_down") ? event->Get<double>("tau_idisoweight_ratio_1_uncert0_down") : 1.0;
    tau_idisoweight_ratio_2_uncert0_down_ = event->Exists("tau_idisoweight_ratio_2_uncert0_down") ? event->Get<double>("tau_idisoweight_ratio_2_uncert0_down") : 1.0;

    tau_idisoweight_ratio_1_uncert1_up_ = event->Exists("tau_idisoweight_ratio_1_uncert1_up") ? event->Get<double>("tau_idisoweight_ratio_1_uncert1_up") : 1.0;
    tau_idisoweight_ratio_2_uncert1_up_ = event->Exists("tau_idisoweight_ratio_2_uncert1_up") ? event->Get<double>("tau_idisoweight_ratio_2_uncert1_up") : 1.0;

    tau_idisoweight_ratio_1_uncert1_down_ = event->Exists("tau_idisoweight_ratio_1_uncert1_down") ? event->Get<double>("tau_idisoweight_ratio_1_uncert1_down") : 1.0;
    tau_idisoweight_ratio_2_uncert1_down_ = event->Exists("tau_idisoweight_ratio_2_uncert1_down") ? event->Get<double>("tau_idisoweight_ratio_2_uncert1_down") : 1.0;

    tau_idisoweight_ratio_1_syst_all_eras_up_ = event->Exists("tau_idisoweight_ratio_1_syst_all_eras_up") ? event->Get<double>("tau_idisoweight_ratio_1_syst_all_eras_up") : 1.0;
    tau_idisoweight_ratio_2_syst_all_eras_up_ = event->Exists("tau_idisoweight_ratio_2_syst_all_eras_up") ? event->Get<double>("tau_idisoweight_ratio_2_syst_all_eras_up") : 1.0;

    tau_idisoweight_ratio_1_syst_all_eras_down_ = event->Exists("tau_idisoweight_ratio_1_syst_all_eras_down") ? event->Get<double>("tau_idisoweight_ratio_1_syst_all_eras_down") : 1.0;
    tau_idisoweight_ratio_2_syst_all_eras_down_ = event->Exists("tau_idisoweight_ratio_2_syst_all_eras_down") ? event->Get<double>("tau_idisoweight_ratio_2_syst_all_eras_down") : 1.0;

    tau_idisoweight_ratio_1_syst_year_up_ = event->Exists("tau_idisoweight_ratio_1_syst_year_up") ? event->Get<double>("tau_idisoweight_ratio_1_syst_year_up") : 1.0;
    tau_idisoweight_ratio_2_syst_year_up_ = event->Exists("tau_idisoweight_ratio_2_syst_year_up") ? event->Get<double>("tau_idisoweight_ratio_2_syst_year_up") : 1.0;

    tau_idisoweight_ratio_1_syst_year_down_ = event->Exists("tau_idisoweight_ratio_1_syst_year_down") ? event->Get<double>("tau_idisoweight_ratio_1_syst_year_down") : 1.0;
    tau_idisoweight_ratio_2_syst_year_down_ = event->Exists("tau_idisoweight_ratio_2_syst_year_down") ? event->Get<double>("tau_idisoweight_ratio_2_syst_year_down") : 1.0;

    tau_idisoweight_ratio_1_syst_dm_year_up_ = event->Exists("tau_idisoweight_ratio_1_syst_dm_year_up") ? event->Get<double>("tau_idisoweight_ratio_1_syst_dm_year_up") : 1.0;
    tau_idisoweight_ratio_2_syst_dm_year_up_ = event->Exists("tau_idisoweight_ratio_2_syst_dm_year_up") ? event->Get<double>("tau_idisoweight_ratio_2_syst_dm_year_up") : 1.0;

    tau_idisoweight_ratio_1_syst_dm_year_down_ = event->Exists("tau_idisoweight_ratio_1_syst_dm_year_down") ? event->Get<double>("tau_idisoweight_ratio_1_syst_dm_year_down") : 1.0;
    tau_idisoweight_ratio_2_syst_dm_year_down_ = event->Exists("tau_idisoweight_ratio_2_syst_dm_year_down") ? event->Get<double>("tau_idisoweight_ratio_2_syst_dm_year_down") : 1.0;

    tau_idisoweight_ratio_1_highpt_up_ = event->Exists("tau_idisoweight_ratio_1_highpt_up") ? event->Get<double>("tau_idisoweight_ratio_1_highpt_up") : 1.0;
    tau_idisoweight_ratio_2_highpt_up_ = event->Exists("tau_idisoweight_ratio_2_highpt_up") ? event->Get<double>("tau_idisoweight_ratio_2_highpt_up") : 1.0;
    
    tau_idisoweight_ratio_1_highpt_down_ = event->Exists("tau_idisoweight_ratio_1_highpt_down") ? event->Get<double>("tau_idisoweight_ratio_1_highpt_down") : 1.0;
    tau_idisoweight_ratio_2_highpt_down_ = event->Exists("tau_idisoweight_ratio_2_highpt_down") ? event->Get<double>("tau_idisoweight_ratio_2_highpt_down") : 1.0;

    tau_idisoweight_ratio_1_down_ = event->Exists("tau_idisoweight_ratio_1_down") ? event->Get<double>("tau_idisoweight_ratio_1_down") : 1.0;
    tau_idisoweight_ratio_2_down_ = event->Exists("tau_idisoweight_ratio_2_down") ? event->Get<double>("tau_idisoweight_ratio_2_down") : 1.0;

    tau_idisoweight_vvvloose_1_ = event->Exists("tau_idisoweight_vvvloose_1") ? event->Get<double>("tau_idisoweight_vvvloose_1") : 1.0;
    tau_idisoweight_vvvloose_2_ = event->Exists("tau_idisoweight_vvvloose_2") ? event->Get<double>("tau_idisoweight_vvvloose_2") : 1.0;

    tau_idisoweight_vvloose_1_ = event->Exists("tau_idisoweight_vvloose_1") ? event->Get<double>("tau_idisoweight_vvloose_1") : 1.0;
    tau_idisoweight_vvloose_2_ = event->Exists("tau_idisoweight_vvloose_2") ? event->Get<double>("tau_idisoweight_vvloose_2") : 1.0;

    tau_idisoweight_vloose_1_ = event->Exists("tau_idisoweight_vloose_1") ? event->Get<double>("tau_idisoweight_vloose_1") : 1.0;
    tau_idisoweight_vloose_2_ = event->Exists("tau_idisoweight_vloose_2") ? event->Get<double>("tau_idisoweight_vloose_2") : 1.0;

    tau_idisoweight_loose_1_ = event->Exists("tau_idisoweight_loose_1") ? event->Get<double>("tau_idisoweight_loose_1") : 1.0;
    tau_idisoweight_loose_2_ = event->Exists("tau_idisoweight_loose_2") ? event->Get<double>("tau_idisoweight_loose_2") : 1.0;

    tau_idisoweight_medium_1_ = event->Exists("tau_idisoweight_medium_1") ? event->Get<double>("tau_idisoweight_medium_1") : 1.0;
    tau_idisoweight_medium_2_ = event->Exists("tau_idisoweight_medium_2") ? event->Get<double>("tau_idisoweight_medium_2") : 1.0;

    tau_idisoweight_tight_1_ = event->Exists("tau_idisoweight_tight_1") ? event->Get<double>("tau_idisoweight_tight_1") : 1.0;
    tau_idisoweight_tight_2_ = event->Exists("tau_idisoweight_tight_2") ? event->Get<double>("tau_idisoweight_tight_2") : 1.0;
   
    wt_tau_id_dm0_up_ =  (event->Exists("wt_tau_id_dm0_up")) ? event->Get<double>("wt_tau_id_dm0_up") : 1.;
    wt_tau_id_dm1_up_ =  (event->Exists("wt_tau_id_dm1_up")) ? event->Get<double>("wt_tau_id_dm1_up") : 1.;
    wt_tau_id_dm10_up_ =  (event->Exists("wt_tau_id_dm10_up")) ? event->Get<double>("wt_tau_id_dm10_up") : 1.;
    wt_tau_id_dm11_up_ =  (event->Exists("wt_tau_id_dm11_up")) ? event->Get<double>("wt_tau_id_dm11_up") : 1.;
    wt_tau_id_dm0_down_ =  (event->Exists("wt_tau_id_dm0_down")) ? event->Get<double>("wt_tau_id_dm0_down") : 1.;
    wt_tau_id_dm1_down_ =  (event->Exists("wt_tau_id_dm1_down")) ? event->Get<double>("wt_tau_id_dm1_down") : 1.;
    wt_tau_id_dm10_down_ =  (event->Exists("wt_tau_id_dm10_down")) ? event->Get<double>("wt_tau_id_dm10_down") : 1.;
    wt_tau_id_dm11_down_ =  (event->Exists("wt_tau_id_dm11_down")) ? event->Get<double>("wt_tau_id_dm11_down") : 1.;

    wt_tau_id_pt_bin1_up_ =  (event->Exists("wt_tau_id_pt_bin1_up")) ? event->Get<double>("wt_tau_id_pt_bin1_up") : 1.;
    wt_tau_id_pt_bin2_up_ =  (event->Exists("wt_tau_id_pt_bin2_up")) ? event->Get<double>("wt_tau_id_pt_bin2_up") : 1.;
    wt_tau_id_pt_bin3_up_ =  (event->Exists("wt_tau_id_pt_bin3_up")) ? event->Get<double>("wt_tau_id_pt_bin3_up") : 1.;
    wt_tau_id_pt_bin4_up_ =  (event->Exists("wt_tau_id_pt_bin4_up")) ? event->Get<double>("wt_tau_id_pt_bin4_up") : 1.;
    wt_tau_id_pt_bin5_up_ =  (event->Exists("wt_tau_id_pt_bin5_up")) ? event->Get<double>("wt_tau_id_pt_bin5_up") : 1.;
    wt_tau_id_pt_bin1_down_ =  (event->Exists("wt_tau_id_pt_bin1_down")) ? event->Get<double>("wt_tau_id_pt_bin1_down") : 1.;
    wt_tau_id_pt_bin2_down_ =  (event->Exists("wt_tau_id_pt_bin2_down")) ? event->Get<double>("wt_tau_id_pt_bin2_down") : 1.;
    wt_tau_id_pt_bin3_down_ =  (event->Exists("wt_tau_id_pt_bin3_down")) ? event->Get<double>("wt_tau_id_pt_bin3_down") : 1.;
    wt_tau_id_pt_bin4_down_ =  (event->Exists("wt_tau_id_pt_bin4_down")) ? event->Get<double>("wt_tau_id_pt_bin4_down") : 1.;
    wt_tau_id_pt_bin5_down_ =  (event->Exists("wt_tau_id_pt_bin5_down")) ? event->Get<double>("wt_tau_id_pt_bin5_down") : 1.;

    wt_tau_trg_dm0_up_ =  (event->Exists("wt_tau_trg_dm0_up")) ? event->Get<double>("wt_tau_trg_dm0_up") : 1.;
    wt_tau_trg_dm1_up_ =  (event->Exists("wt_tau_trg_dm1_up")) ? event->Get<double>("wt_tau_trg_dm1_up") : 1.;
    wt_tau_trg_dm10_up_ =  (event->Exists("wt_tau_trg_dm10_up")) ? event->Get<double>("wt_tau_trg_dm10_up") : 1.;
    wt_tau_trg_dm11_up_ =  (event->Exists("wt_tau_trg_dm11_up")) ? event->Get<double>("wt_tau_trg_dm11_up") : 1.;
    wt_tau_trg_dm0_down_ =  (event->Exists("wt_tau_trg_dm0_down")) ? event->Get<double>("wt_tau_trg_dm0_down") : 1.;
    wt_tau_trg_dm1_down_ =  (event->Exists("wt_tau_trg_dm1_down")) ? event->Get<double>("wt_tau_trg_dm1_down") : 1.;
    wt_tau_trg_dm10_down_ =  (event->Exists("wt_tau_trg_dm10_down")) ? event->Get<double>("wt_tau_trg_dm10_down") : 1.;
    wt_tau_trg_dm11_down_ =  (event->Exists("wt_tau_trg_dm11_down")) ? event->Get<double>("wt_tau_trg_dm11_down") : 1.;

    wt_tau_trg_mvadm0_up_ =     (event->Exists("wt_tau_trg_mvadm0_up")) ?    event->Get<double>("wt_tau_trg_mvadm0_up") : 1.;
    wt_tau_trg_mvadm1_up_ =     (event->Exists("wt_tau_trg_mvadm1_up")) ?    event->Get<double>("wt_tau_trg_mvadm1_up") : 1.;
    wt_tau_trg_mvadm2_up_ =     (event->Exists("wt_tau_trg_mvadm2_up")) ?    event->Get<double>("wt_tau_trg_mvadm2_up") : 1.;
    wt_tau_trg_mvadm10_up_ =    (event->Exists("wt_tau_trg_mvadm10_up")) ?   event->Get<double>("wt_tau_trg_mvadm10_up") : 1.;
    wt_tau_trg_mvadm11_up_ =    (event->Exists("wt_tau_trg_mvadm11_up")) ?   event->Get<double>("wt_tau_trg_mvadm11_up") : 1.;
    wt_tau_trg_mvadm0_down_ =   (event->Exists("wt_tau_trg_mvadm0_down")) ?  event->Get<double>("wt_tau_trg_mvadm0_down") : 1.;
    wt_tau_trg_mvadm1_down_ =   (event->Exists("wt_tau_trg_mvadm1_down")) ?  event->Get<double>("wt_tau_trg_mvadm1_down") : 1.;
    wt_tau_trg_mvadm2_down_ =   (event->Exists("wt_tau_trg_mvadm2_down")) ?  event->Get<double>("wt_tau_trg_mvadm2_down") : 1.;
    wt_tau_trg_mvadm10_down_ =  (event->Exists("wt_tau_trg_mvadm10_down")) ? event->Get<double>("wt_tau_trg_mvadm10_down") : 1.;
    wt_tau_trg_mvadm11_down_ =  (event->Exists("wt_tau_trg_mvadm11_down")) ? event->Get<double>("wt_tau_trg_mvadm11_down") : 1.;

    wt_tau_trg_ic_ =     (event->Exists("wt_tau_trg_ic")) ?    event->Get<double>("wt_tau_trg_ic") : 1.;
    wt_tau_trg_mvadm_ =     (event->Exists("wt_tau_trg_mvadm")) ?    event->Get<double>("wt_tau_trg_mvadm") : 1.;

    wt_tau_id_mvadm_                  = (event->Exists("wt_tau_id_mvadm")) ? event->Get<double>("wt_tau_id_mvadm") : 1.0;
    wt_tau_id_lowpt_mvadm0_up_     = (event->Exists("wt_tau_id_lowpt_mvadm0_up")) ? event->Get<double>("wt_tau_id_lowpt_mvadm0_up") : 1.0;
    wt_tau_id_lowpt_mvadm1_up_     = (event->Exists("wt_tau_id_lowpt_mvadm1_up")) ? event->Get<double>("wt_tau_id_lowpt_mvadm1_up") : 1.0;
    wt_tau_id_lowpt_mvadm2_up_     = (event->Exists("wt_tau_id_lowpt_mvadm2_up")) ? event->Get<double>("wt_tau_id_lowpt_mvadm2_up") : 1.0;
    wt_tau_id_lowpt_mvadm10_up_    = (event->Exists("wt_tau_id_lowpt_mvadm10_up")) ? event->Get<double>("wt_tau_id_lowpt_mvadm10_up") : 1.0;
    wt_tau_id_lowpt_mvadm11_up_    = (event->Exists("wt_tau_id_lowpt_mvadm11_up")) ? event->Get<double>("wt_tau_id_lowpt_mvadm11_up") : 1.0;
    wt_tau_id_highpt_mvadm0_up_    = (event->Exists("wt_tau_id_highpt_mvadm0_up")) ? event->Get<double>("wt_tau_id_highpt_mvadm0_up") : 1.0;
    wt_tau_id_highpt_mvadm1_up_    = (event->Exists("wt_tau_id_highpt_mvadm1_up")) ? event->Get<double>("wt_tau_id_highpt_mvadm1_up") : 1.0;
    wt_tau_id_highpt_mvadm2_up_    = (event->Exists("wt_tau_id_highpt_mvadm2_up")) ? event->Get<double>("wt_tau_id_highpt_mvadm2_up") : 1.0;
    wt_tau_id_highpt_mvadm10_up_   = (event->Exists("wt_tau_id_highpt_mvadm10_up")) ? event->Get<double>("wt_tau_id_highpt_mvadm10_up") : 1.0;
    wt_tau_id_highpt_mvadm11_up_   = (event->Exists("wt_tau_id_highpt_mvadm11_up")) ? event->Get<double>("wt_tau_id_highpt_mvadm11_up") : 1.0;
    wt_tau_id_lowpt_mvadm0_down_   = (event->Exists("wt_tau_id_lowpt_mvadm0_down")) ? event->Get<double>("wt_tau_id_lowpt_mvadm0_down") : 1.0;
    wt_tau_id_lowpt_mvadm1_down_   = (event->Exists("wt_tau_id_lowpt_mvadm1_down")) ? event->Get<double>("wt_tau_id_lowpt_mvadm1_down") : 1.0;
    wt_tau_id_lowpt_mvadm2_down_   = (event->Exists("wt_tau_id_lowpt_mvadm2_down")) ? event->Get<double>("wt_tau_id_lowpt_mvadm2_down") : 1.0;
    wt_tau_id_lowpt_mvadm10_down_  = (event->Exists("wt_tau_id_lowpt_mvadm10_down")) ? event->Get<double>("wt_tau_id_lowpt_mvadm10_down") : 1.0;
    wt_tau_id_lowpt_mvadm11_down_  = (event->Exists("wt_tau_id_lowpt_mvadm11_down")) ? event->Get<double>("wt_tau_id_lowpt_mvadm11_down") : 1.0;
    wt_tau_id_highpt_mvadm0_down_  = (event->Exists("wt_tau_id_highpt_mvadm0_down")) ? event->Get<double>("wt_tau_id_highpt_mvadm0_down") : 1.0;
    wt_tau_id_highpt_mvadm1_down_  = (event->Exists("wt_tau_id_highpt_mvadm1_down")) ? event->Get<double>("wt_tau_id_highpt_mvadm1_down") : 1.0;
    wt_tau_id_highpt_mvadm2_down_  = (event->Exists("wt_tau_id_highpt_mvadm2_down")) ? event->Get<double>("wt_tau_id_highpt_mvadm2_down") : 1.0;
    wt_tau_id_highpt_mvadm10_down_ = (event->Exists("wt_tau_id_highpt_mvadm10_down")) ? event->Get<double>("wt_tau_id_highpt_mvadm10_down") : 1.0;
    wt_tau_id_highpt_mvadm11_down_ = (event->Exists("wt_tau_id_highpt_mvadm11_down")) ? event->Get<double>("wt_tau_id_highpt_mvadm11_down") : 1.0;

    wt_tau_trg_mssm_doubleonly_ =     (event->Exists("wt_tau_trg_mssm_doubleonly")) ?    event->Get<double>("wt_tau_trg_mssm_doubleonly") : 1.;
    wt_tau_trg_mssm_ =     (event->Exists("wt_tau_trg_mssm")) ?    event->Get<double>("wt_tau_trg_mssm") : 1.;
    wt_tau_trg_mssm_dm0_up_ =     (event->Exists("wt_tau_trg_mssm_lowpt_dm0_up")) ?    event->Get<double>("wt_tau_trg_mssm_lowpt_dm0_up") : 1.;
    wt_tau_trg_mssm_dm0_down_ =     (event->Exists("wt_tau_trg_mssm_lowpt_dm0_down")) ?    event->Get<double>("wt_tau_trg_mssm_lowpt_dm0_down") : 1.;
    wt_tau_trg_mssm_dm1_up_ =     (event->Exists("wt_tau_trg_mssm_lowpt_dm1_up")) ?    event->Get<double>("wt_tau_trg_mssm_lowpt_dm1_up") : 1.;
    wt_tau_trg_mssm_dm1_down_ =     (event->Exists("wt_tau_trg_mssm_lowpt_dm1_down")) ?    event->Get<double>("wt_tau_trg_mssm_lowpt_dm1_down") : 1.;
    wt_tau_trg_mssm_dm10_up_ =     (event->Exists("wt_tau_trg_mssm_lowpt_dm10_up")) ?    event->Get<double>("wt_tau_trg_mssm_lowpt_dm10_up") : 1.;
    wt_tau_trg_mssm_dm10_down_ =     (event->Exists("wt_tau_trg_mssm_lowpt_dm10_down")) ?    event->Get<double>("wt_tau_trg_mssm_lowpt_dm10_down") : 1.;
    wt_tau_trg_mssm_dm11_up_ =     (event->Exists("wt_tau_trg_mssm_lowpt_dm11_up")) ?    event->Get<double>("wt_tau_trg_mssm_lowpt_dm11_up") : 1.;
    wt_tau_trg_mssm_dm11_down_ =     (event->Exists("wt_tau_trg_mssm_lowpt_dm11_down")) ?    event->Get<double>("wt_tau_trg_mssm_lowpt_dm11_down") : 1.;
    wt_tau_trg_mssm_highpt_dm0_up_ =     (event->Exists("wt_tau_trg_mssm_highpt_dm0_up")) ?    event->Get<double>("wt_tau_trg_mssm_highpt_dm0_up") : 1.;
    wt_tau_trg_mssm_highpt_dm0_down_ =     (event->Exists("wt_tau_trg_mssm_highpt_dm0_down")) ?    event->Get<double>("wt_tau_trg_mssm_highpt_dm0_down") : 1.;
    wt_tau_trg_mssm_highpt_dm1_up_ =     (event->Exists("wt_tau_trg_mssm_highpt_dm1_up")) ?    event->Get<double>("wt_tau_trg_mssm_highpt_dm1_up") : 1.;
    wt_tau_trg_mssm_highpt_dm1_down_ =     (event->Exists("wt_tau_trg_mssm_highpt_dm1_down")) ?    event->Get<double>("wt_tau_trg_mssm_highpt_dm1_down") : 1.;
    wt_tau_trg_mssm_highpt_dm10_up_ =     (event->Exists("wt_tau_trg_mssm_highpt_dm10_up")) ?    event->Get<double>("wt_tau_trg_mssm_highpt_dm10_up") : 1.;
    wt_tau_trg_mssm_highpt_dm10_down_ =     (event->Exists("wt_tau_trg_mssm_highpt_dm10_down")) ?    event->Get<double>("wt_tau_trg_mssm_highpt_dm10_down") : 1.;
    wt_tau_trg_mssm_highpt_dm11_up_ =     (event->Exists("wt_tau_trg_mssm_highpt_dm11_up")) ?    event->Get<double>("wt_tau_trg_mssm_highpt_dm11_up") : 1.;
    wt_tau_trg_mssm_highpt_dm11_down_ =     (event->Exists("wt_tau_trg_mssm_highpt_dm11_down")) ?    event->Get<double>("wt_tau_trg_mssm_highpt_dm11_down") : 1.;
    wt_tau_trg_mssm_singletau_up_ =     (event->Exists("wt_tau_trg_mssm_singletau_up")) ?    event->Get<double>("wt_tau_trg_mssm_singletau_up") : 1.;
    wt_tau_trg_mssm_singletau_down_ =     (event->Exists("wt_tau_trg_mssm_singletau_down")) ?    event->Get<double>("wt_tau_trg_mssm_singletau_down") : 1.;
    wt_tau_trg_mssm_singlelep_up_ =     (event->Exists("wt_tau_trg_mssm_singlelep_up")) ?    event->Get<double>("wt_tau_trg_mssm_singlelep_up") : 1.;
    wt_tau_trg_mssm_singlelep_down_ =     (event->Exists("wt_tau_trg_mssm_singlelep_down")) ?    event->Get<double>("wt_tau_trg_mssm_singlelep_down") : 1.;
    wt_tau_trg_mssm_crosslep_up_ =     (event->Exists("wt_tau_trg_mssm_crosslep_up")) ?    event->Get<double>("wt_tau_trg_mssm_crosslep_up") : 1.;
    wt_tau_trg_mssm_crosslep_down_ =     (event->Exists("wt_tau_trg_mssm_crosslep_down")) ?    event->Get<double>("wt_tau_trg_mssm_crosslep_down") : 1.;

    wt_tau_id_mssm_ =     (event->Exists("wt_tau_id_mssm")) ?    event->Get<double>("wt_tau_id_mssm") : 1.;
    wt_tau_id_mssm_bin5_up_ =     (event->Exists("wt_tau_id_mssm_bin5_up")) ?    event->Get<double>("wt_tau_id_mssm_bin5_up") : 1.;
    wt_tau_id_mssm_bin5_down_ =     (event->Exists("wt_tau_id_mssm_bin5_down")) ?    event->Get<double>("wt_tau_id_mssm_bin5_down") : 1.;
    wt_tau_id_mssm_bin6_up_ =     (event->Exists("wt_tau_id_mssm_bin6_up")) ?    event->Get<double>("wt_tau_id_mssm_bin6_up") : 1.;
    wt_tau_id_mssm_bin6_down_ =     (event->Exists("wt_tau_id_mssm_bin6_down")) ?    event->Get<double>("wt_tau_id_mssm_bin6_down") : 1.;

    wt_tau_sf_id_dm_vs_pt_ = (event->Exists("tau_sf_id_dm_vs_pt")) ?    event->Get<double>("tau_sf_id_dm_vs_pt") : 1.;
    wt_tau_sf_id_dm_ = (event->Exists("tau_sf_id_dm")) ?    event->Get<double>("tau_sf_id_dm") : 1.;
    wt_tau_sf_id_pt_ = (event->Exists("tau_sf_id_pt")) ?    event->Get<double>("tau_sf_id_pt") : 1.;

    wt_zpt_ = (event->Exists("zpt_sf")) ?    event->Get<double>("zpt_sf") : 1.;
    wt_zpt_NLO_ = (event->Exists("zpt_sf_NLO")) ?    event->Get<double>("zpt_sf_NLO") : 1.;

    if(do_mssm_higgspt_){
      wt_ggh_t_ = event->Exists("wt_ggh_t") ? event->Get<double>("wt_ggh_t") : 1.0;
      wt_ggh_b_ = event->Exists("wt_ggh_b") ? event->Get<double>("wt_ggh_b") : 1.0;
      wt_ggh_i_ = event->Exists("wt_ggh_i") ? event->Get<double>("wt_ggh_i") : 1.0;
      wt_ggH_t_ = event->Exists("wt_ggH_t") ? event->Get<double>("wt_ggH_t") : 1.0;
      wt_ggH_b_ = event->Exists("wt_ggH_b") ? event->Get<double>("wt_ggH_b") : 1.0;
      wt_ggH_i_ = event->Exists("wt_ggH_i") ? event->Get<double>("wt_ggH_i") : 1.0;
      wt_ggA_t_ = event->Exists("wt_ggA_t") ? event->Get<double>("wt_ggA_t") : 1.0;
      wt_ggA_b_ = event->Exists("wt_ggA_b") ? event->Get<double>("wt_ggA_b") : 1.0;
      wt_ggA_i_ = event->Exists("wt_ggA_i") ? event->Get<double>("wt_ggA_i") : 1.0;  

      wt_ggh_t_scale_up_ = event->Exists("wt_ggh_t_scale_up") ? event->Get<double>("wt_ggh_t_scale_up") : 1.0;
      wt_ggh_b_scale_up_ = event->Exists("wt_ggh_b_scale_up") ? event->Get<double>("wt_ggh_b_scale_up") : 1.0;
      wt_ggh_i_scale_up_ = event->Exists("wt_ggh_i_scale_up") ? event->Get<double>("wt_ggh_i_scale_up") : 1.0;
      wt_ggA_t_scale_up_ = event->Exists("wt_ggA_t_scale_up") ? event->Get<double>("wt_ggA_t_scale_up") : 1.0;
      wt_ggA_b_scale_up_ = event->Exists("wt_ggA_b_scale_up") ? event->Get<double>("wt_ggA_b_scale_up") : 1.0;
      wt_ggA_i_scale_up_ = event->Exists("wt_ggA_i_scale_up") ? event->Get<double>("wt_ggA_i_scale_up") : 1.0;

      wt_ggh_t_scale_down_ = event->Exists("wt_ggh_t_scale_down") ? event->Get<double>("wt_ggh_t_scale_down") : 1.0;
      wt_ggh_b_scale_down_ = event->Exists("wt_ggh_b_scale_down") ? event->Get<double>("wt_ggh_b_scale_down") : 1.0;
      wt_ggh_i_scale_down_ = event->Exists("wt_ggh_i_scale_down") ? event->Get<double>("wt_ggh_i_scale_down") : 1.0;
      wt_ggA_t_scale_down_ = event->Exists("wt_ggA_t_scale_down") ? event->Get<double>("wt_ggA_t_scale_down") : 1.0;
      wt_ggA_b_scale_down_ = event->Exists("wt_ggA_b_scale_down") ? event->Get<double>("wt_ggA_b_scale_down") : 1.0;
      wt_ggA_i_scale_down_ = event->Exists("wt_ggA_i_scale_down") ? event->Get<double>("wt_ggA_i_scale_down") : 1.0;

      wt_ggh_t_hdamp_up_ = event->Exists("wt_ggh_t_hdamp_up") ? event->Get<double>("wt_ggh_t_hdamp_up") : 1.0;
      wt_ggh_b_hdamp_up_ = event->Exists("wt_ggh_b_hdamp_up") ? event->Get<double>("wt_ggh_b_hdamp_up") : 1.0;
      wt_ggh_i_hdamp_up_ = event->Exists("wt_ggh_i_hdamp_up") ? event->Get<double>("wt_ggh_i_hdamp_up") : 1.0;
      wt_ggA_t_hdamp_up_ = event->Exists("wt_ggA_t_hdamp_up") ? event->Get<double>("wt_ggA_t_hdamp_up") : 1.0;
      wt_ggA_b_hdamp_up_ = event->Exists("wt_ggA_b_hdamp_up") ? event->Get<double>("wt_ggA_b_hdamp_up") : 1.0;
      wt_ggA_i_hdamp_up_ = event->Exists("wt_ggA_i_hdamp_up") ? event->Get<double>("wt_ggA_i_hdamp_up") : 1.0;

      wt_ggh_t_hdamp_down_ = event->Exists("wt_ggh_t_hdamp_down") ? event->Get<double>("wt_ggh_t_hdamp_down") : 1.0;
      wt_ggh_b_hdamp_down_ = event->Exists("wt_ggh_b_hdamp_down") ? event->Get<double>("wt_ggh_b_hdamp_down") : 1.0;
      wt_ggh_i_hdamp_down_ = event->Exists("wt_ggh_i_hdamp_down") ? event->Get<double>("wt_ggh_i_hdamp_down") : 1.0;
      wt_ggA_t_hdamp_down_ = event->Exists("wt_ggA_t_hdamp_down") ? event->Get<double>("wt_ggA_t_hdamp_down") : 1.0;
      wt_ggA_b_hdamp_down_ = event->Exists("wt_ggA_b_hdamp_down") ? event->Get<double>("wt_ggA_b_hdamp_down") : 1.0;
      wt_ggA_i_hdamp_down_ = event->Exists("wt_ggA_i_hdamp_down") ? event->Get<double>("wt_ggA_i_hdamp_down") : 1.0; 
    }
    wt_quarkmass_ = 1.0;
    wt_fullquarkmass_ =1.0;
    if(event->Exists("wt_quarkmass")) wt_quarkmass_ = event->Get<double>("wt_quarkmass");
    if(event->Exists("wt_fullquarkmass")) wt_fullquarkmass_ = event->Get<double>("wt_fullquarkmass");

    wt_ue_up_    = event->Exists("wt_ue_up") ? event->Get<double>("wt_ue_up") : 1.0;
    wt_ue_down_  = event->Exists("wt_ue_down") ? event->Get<double>("wt_ue_down") : 1.0;

    wt_prefire_ = eventInfo->weight_defined("wt_prefire") ? eventInfo->weight("wt_prefire") : 1.0;
    wt_prefire_up_ = eventInfo->weight_defined("wt_prefire_up") ? eventInfo->weight("wt_prefire_up") : 1.0;
    wt_prefire_down_ = eventInfo->weight_defined("wt_prefire_down") ? eventInfo->weight("wt_prefire_down") : 1.0;

    wt_tau_id_extra_ = event->Exists("wt_tau_id_extra") ? event->Get<double>("wt_tau_id_extra") : 1.0;
    wt_mg_nnlops_ = event->Exists("wt_mg_nnlops") ? event->Get<double>("wt_mg_nnlops") : 1.0;
    wt_ph_nnlops_ = event->Exists("wt_ph_nnlops") ? event->Get<double>("wt_ph_nnlops") : 1.0;
    
    run_ = eventInfo->run();
    event_ = (unsigned long long) eventInfo->event();
    lumi_ = eventInfo->lumi_block();

    // temporary fix to removed duplicated events
    if (run_==prev_run_ && event_==prev_event_ && lumi_==prev_lumi_) { 
      return 1;
    }
    prev_run_ = run_;
    prev_event_ = event_;
    prev_lumi_ = lumi_;

    rand->SetSeed(event_);
    rand_ = rand->Uniform();
  
    // fake-factor weights        
    if (do_ff_weights_ && (channel_ == channel::et || channel_ == channel::mt || channel_ == channel::tt)){
      
      if(event->Exists("wt_ff_1")) wt_ff_1_ = event->Get<double>("wt_ff_1");
      if(event->Exists("wt_ff_dmbins_1")) wt_ff_dmbins_1_ = event->Get<double>("wt_ff_dmbins_1");

      if(channel_ == channel::et || channel_ == channel::mt) {
        if(event->Exists("wt_ff_qcd_1")) wt_ff_qcd_1_ = event->Get<double>("wt_ff_qcd_1");
        if(event->Exists("wt_ff_wjets_1")) wt_ff_wjets_1_ = event->Get<double>("wt_ff_wjets_1");
        if(event->Exists("wt_ff_ttbar_1")) wt_ff_ttbar_1_ = event->Get<double>("wt_ff_ttbar_1");

        if(event->Exists("wt_ff_dmbins_qcd_1")) wt_ff_dmbins_qcd_1_ = event->Get<double>("wt_ff_dmbins_qcd_1");
        if(event->Exists("wt_ff_dmbins_wjets_1")) wt_ff_dmbins_wjets_1_ = event->Get<double>("wt_ff_dmbins_wjets_1");
        if(event->Exists("wt_ff_dmbins_ttbar_1")) wt_ff_dmbins_ttbar_1_ = event->Get<double>("wt_ff_dmbins_ttbar_1");
      }

      wt_ff_mssm_1_ = event->Exists("wt_ff_mssm_1") ? event->Get<double>("wt_ff_mssm_1") : 0.0;

      if(channel_ == channel::et || channel_ == channel::mt) {
        wt_ff_mssm_qcd_1_ = event->Exists("wt_ff_mssm_qcd_1") ? event->Get<double>("wt_ff_mssm_qcd_1") : 0.0;
        wt_ff_mssm_ttbar_1_ = event->Exists("wt_ff_mssm_ttbar_1") ? event->Get<double>("wt_ff_mssm_ttbar_1") : 0.0;
        wt_ff_mssm_wjets_1_ = event->Exists("wt_ff_mssm_wjets_1") ? event->Get<double>("wt_ff_mssm_wjets_1") : 0.0;
        ff_frac_mssm_wjets_ = event->Exists("ff_frac_mssm_wjets") ? event->Get<double>("ff_frac_mssm_wjets") : 0.0;
        ff_frac_mssm_qcd_ = event->Exists("ff_frac_mssm_qcd") ? event->Get<double>("ff_frac_mssm_qcd") : 0.0;
        ff_frac_mssm_ttbar_ = event->Exists("ff_frac_mssm_ttbar") ? event->Get<double>("ff_frac_mssm_ttbar") : 0.0;
      }

      if(do_ff_systematics_){

        wt_ff_wjets_syst_up_1_ = event->Exists("wt_ff_wjets_syst_up_1") ? event->Get<double>("wt_ff_wjets_syst_up_1") : 0.0;
        wt_ff_wjets_syst_down_1_ = event->Exists("wt_ff_wjets_syst_down_1") ? event->Get<double>("wt_ff_wjets_syst_down_1") : 0.0;
        wt_ff_wjets_met_up_1_ = event->Exists("wt_ff_wjets_met_up_1") ? event->Get<double>("wt_ff_wjets_met_up_1") : 0.0;
        wt_ff_wjets_met_down_1_ = event->Exists("wt_ff_wjets_met_down_1") ? event->Get<double>("wt_ff_wjets_met_down_1") : 0.0;
        wt_ff_wjets_l_pt_up_1_ = event->Exists("wt_ff_wjets_l_pt_up_1") ? event->Get<double>("wt_ff_wjets_l_pt_up_1") : 0.0;
        wt_ff_wjets_l_pt_down_1_ = event->Exists("wt_ff_wjets_l_pt_down_1") ? event->Get<double>("wt_ff_wjets_l_pt_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_lt3_down_1") : 0.0;    
        wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_lt3_down_1") : 0.0;    
        wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_lt3_down_1") : 0.0;    
        wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm0_sig_gt3_down_1") : 0.0;    
        wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm0_sig_gt3_down_1") : 0.0;    
        wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm0_sig_gt3_down_1") : 0.0;    
        wt_ff_wjets_stat_unc1_njet0_mvadm1_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm1_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm1_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm1_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm1_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm1_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet1_mvadm1_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm1_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm1_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet1_mvadm1_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm1_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm1_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet2_mvadm1_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm1_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm1_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet2_mvadm1_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm1_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm1_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm2_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm2_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm2_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm2_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm2_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm2_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet1_mvadm2_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm2_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm2_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet1_mvadm2_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm2_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm2_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet2_mvadm2_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm2_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm2_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet2_mvadm2_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm2_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm2_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm10_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm10_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm10_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm10_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm10_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm10_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet1_mvadm10_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm10_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm10_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet1_mvadm10_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm10_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm10_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet2_mvadm10_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm10_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm10_up_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet2_mvadm10_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm10_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm10_down_1") : 0.0;
        wt_ff_wjets_stat_unc1_njet0_mvadm11_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm11_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm11_up_1") : 0.0; 
        wt_ff_wjets_stat_unc1_njet0_mvadm11_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet0_mvadm11_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet0_mvadm11_down_1") : 0.0; 
        wt_ff_wjets_stat_unc1_njet1_mvadm11_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm11_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm11_up_1") : 0.0; 
        wt_ff_wjets_stat_unc1_njet1_mvadm11_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet1_mvadm11_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet1_mvadm11_down_1") : 0.0; 
        wt_ff_wjets_stat_unc1_njet2_mvadm11_up_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm11_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm11_up_1") : 0.0; 
        wt_ff_wjets_stat_unc1_njet2_mvadm11_down_1_ = event->Exists("wt_ff_wjets_stat_unc1_njet2_mvadm11_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc1_njet2_mvadm11_down_1") : 0.0;

        wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_lt3_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_lt3_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_lt3_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm0_sig_gt3_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm0_sig_gt3_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm0_sig_gt3_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm1_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm1_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm1_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm1_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm1_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm1_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm1_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm1_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm1_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm1_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm1_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm1_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm1_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm1_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm1_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm1_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm1_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm1_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm2_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm2_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm2_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm2_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm2_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm2_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm2_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm2_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm2_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm2_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm2_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm2_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm2_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm2_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm2_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm2_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm2_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm2_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm10_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm10_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm10_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm10_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm10_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm10_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm10_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm10_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm10_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm10_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm10_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm10_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm10_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm10_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm10_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm10_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm10_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm10_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm11_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm11_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm11_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet0_mvadm11_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet0_mvadm11_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet0_mvadm11_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm11_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm11_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm11_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet1_mvadm11_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet1_mvadm11_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet1_mvadm11_down_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm11_up_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm11_up_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm11_up_1") : 0.0;
        wt_ff_wjets_stat_unc2_njet2_mvadm11_down_1_ = event->Exists("wt_ff_wjets_stat_unc2_njet2_mvadm11_down_1") ? event->Get<double>("wt_ff_wjets_stat_unc2_njet2_mvadm11_down_1") : 0.0;

        wt_ff_qcd_syst_up_1_ = event->Exists("wt_ff_qcd_syst_up_1") ? event->Get<double>("wt_ff_qcd_syst_up_1") : 0.0;
        wt_ff_qcd_syst_down_1_ = event->Exists("wt_ff_qcd_syst_down_1") ? event->Get<double>("wt_ff_qcd_syst_down_1") : 0.0;
        wt_ff_qcd_met_up_1_ = event->Exists("wt_ff_qcd_met_up_1") ? event->Get<double>("wt_ff_qcd_met_up_1") : 0.0;
        wt_ff_qcd_met_down_1_ = event->Exists("wt_ff_qcd_met_down_1") ? event->Get<double>("wt_ff_qcd_met_down_1") : 0.0;
        wt_ff_qcd_l_pt_up_1_ = event->Exists("wt_ff_qcd_l_pt_up_1") ? event->Get<double>("wt_ff_qcd_l_pt_up_1") : 0.0;
        wt_ff_qcd_l_pt_down_1_ = event->Exists("wt_ff_qcd_l_pt_down_1") ? event->Get<double>("wt_ff_qcd_l_pt_down_1") : 0.0;
        wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_lt3_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_lt3_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_lt3_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm0_sig_gt3_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm0_sig_gt3_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm0_sig_gt3_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm1_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm1_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm1_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm1_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm1_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm1_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm1_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm1_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm1_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm1_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm1_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm1_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm1_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm1_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm1_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm1_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm1_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm1_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm2_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm2_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm2_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm2_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm2_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm2_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm2_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm2_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm2_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm2_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm2_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm2_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm2_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm2_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm2_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm2_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm2_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm2_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm10_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm10_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm10_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm10_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm10_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm10_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm10_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm10_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm10_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm10_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm10_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm10_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm10_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm10_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm10_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm10_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm10_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm10_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm11_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm11_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm11_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet0_mvadm11_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet0_mvadm11_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet0_mvadm11_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm11_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm11_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm11_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet1_mvadm11_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet1_mvadm11_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet1_mvadm11_down_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm11_up_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm11_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm11_up_1") : 0.0; 
        wt_ff_qcd_stat_unc1_njet2_mvadm11_down_1_ = event->Exists("wt_ff_qcd_stat_unc1_njet2_mvadm11_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc1_njet2_mvadm11_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_lt3_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_lt3_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_lt3_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm0_sig_gt3_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm0_sig_gt3_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm0_sig_gt3_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm1_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm1_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm1_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm1_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm1_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm1_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm1_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm1_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm1_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm1_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm1_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm1_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm1_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm1_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm1_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm1_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm1_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm1_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm2_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm2_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm2_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm2_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm2_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm2_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm2_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm2_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm2_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm2_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm2_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm2_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm2_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm2_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm2_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm2_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm2_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm2_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm10_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm10_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm10_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm10_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm10_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm10_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm10_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm10_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm10_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm10_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm10_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm10_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm10_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm10_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm10_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm10_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm10_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm10_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm11_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm11_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm11_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet0_mvadm11_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet0_mvadm11_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet0_mvadm11_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm11_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm11_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm11_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet1_mvadm11_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet1_mvadm11_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet1_mvadm11_down_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm11_up_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm11_up_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm11_up_1") : 0.0;
        wt_ff_qcd_stat_unc2_njet2_mvadm11_down_1_ = event->Exists("wt_ff_qcd_stat_unc2_njet2_mvadm11_down_1") ? event->Get<double>("wt_ff_qcd_stat_unc2_njet2_mvadm11_down_1") : 0.0;


        wt_ff_dmbins_ttbar_syst_up_1_ = event->Exists("wt_ff_dmbins_ttbar_syst_up_1") ? event->Get<double>("wt_ff_dmbins_ttbar_syst_up_1") : 0.0;
        wt_ff_dmbins_ttbar_syst_down_1_ = event->Exists("wt_ff_dmbins_ttbar_syst_down_1") ? event->Get<double>("wt_ff_dmbins_ttbar_syst_down_1") : 0.0;
        wt_ff_dmbins_wjets_syst_up_1_ = event->Exists("wt_ff_dmbins_wjets_syst_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_syst_up_1") : 0.0;
        wt_ff_dmbins_wjets_syst_down_1_ = event->Exists("wt_ff_dmbins_wjets_syst_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_syst_down_1") : 0.0;
        wt_ff_dmbins_wjets_met_up_1_ = event->Exists("wt_ff_dmbins_wjets_met_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_met_up_1") : 0.0;
        wt_ff_dmbins_wjets_met_down_1_ = event->Exists("wt_ff_dmbins_wjets_met_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_met_down_1") : 0.0;
        wt_ff_dmbins_wjets_l_pt_up_1_ = event->Exists("wt_ff_dmbins_wjets_l_pt_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_l_pt_up_1") : 0.0;
        wt_ff_dmbins_wjets_l_pt_down_1_ = event->Exists("wt_ff_dmbins_wjets_l_pt_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_l_pt_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet0_dm0_down_1") : 0.0;    
        wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet1_dm0_down_1") : 0.0;    
        wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet2_dm0_down_1") : 0.0;    
        wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet0_dm1_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet1_dm1_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet2_dm1_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet0_dm10_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet1_dm10_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet2_dm10_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_up_1") : 0.0; 
        wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet0_dm11_down_1") : 0.0; 
        wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_up_1") : 0.0; 
        wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet1_dm11_down_1") : 0.0; 
        wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_up_1") : 0.0; 
        wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc1_njet2_dm11_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet0_dm0_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet1_dm0_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet2_dm0_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet0_dm1_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet1_dm1_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet2_dm1_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet0_dm10_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet1_dm10_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet2_dm10_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet0_dm11_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet1_dm11_down_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_up_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_up_1") : 0.0;
        wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_down_1_ = event->Exists("wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_wjets_stat_unc2_njet2_dm11_down_1") : 0.0;
        wt_ff_dmbins_qcd_syst_up_1_ = event->Exists("wt_ff_dmbins_qcd_syst_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_syst_up_1") : 0.0;
        wt_ff_dmbins_qcd_syst_down_1_ = event->Exists("wt_ff_dmbins_qcd_syst_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_syst_down_1") : 0.0;
        wt_ff_dmbins_qcd_met_up_1_ = event->Exists("wt_ff_dmbins_qcd_met_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_met_up_1") : 0.0;
        wt_ff_dmbins_qcd_met_down_1_ = event->Exists("wt_ff_dmbins_qcd_met_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_met_down_1") : 0.0;
        wt_ff_dmbins_qcd_l_pt_up_1_ = event->Exists("wt_ff_dmbins_qcd_l_pt_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_l_pt_up_1") : 0.0;
        wt_ff_dmbins_qcd_l_pt_down_1_ = event->Exists("wt_ff_dmbins_qcd_l_pt_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_l_pt_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet0_dm0_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet1_dm0_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet2_dm0_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet0_dm1_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet1_dm1_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet2_dm1_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet0_dm10_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet1_dm10_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet2_dm10_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet0_dm11_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet1_dm11_down_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_up_1") : 0.0; 
        wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc1_njet2_dm11_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet0_dm0_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet1_dm0_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet2_dm0_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet0_dm1_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet1_dm1_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet2_dm1_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet0_dm10_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet1_dm10_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet2_dm10_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet0_dm11_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet1_dm11_down_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_up_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_up_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_up_1") : 0.0;
        wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_down_1_ = event->Exists("wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_down_1") ? event->Get<double>("wt_ff_dmbins_qcd_stat_unc2_njet2_dm11_down_1") : 0.0;

        wt_ff_mssm_qcd_syst_up_1_ = event->Exists("wt_ff_mssm_qcd_syst_up_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_up_1") : 0.0;
        wt_ff_mssm_qcd_syst_down_1_ = event->Exists("wt_ff_mssm_qcd_syst_down_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_down_1") : 0.0;
        wt_ff_mssm_wjets_syst_up_1_ = event->Exists("wt_ff_mssm_wjets_syst_up_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_up_1") : 0.0;
        wt_ff_mssm_wjets_syst_down_1_ = event->Exists("wt_ff_mssm_wjets_syst_down_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_down_1") : 0.0;
        wt_ff_mssm_wjets_syst_extrap_up_1_ = event->Exists("wt_ff_mssm_wjets_syst_extrap_up_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_extrap_up_1") : 0.0;
        wt_ff_mssm_wjets_syst_extrap_down_1_ = event->Exists("wt_ff_mssm_wjets_syst_extrap_down_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_extrap_down_1") : 0.0;
        wt_ff_mssm_ttbar_syst_up_1_ = event->Exists("wt_ff_mssm_ttbar_syst_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_syst_up_1") : 0.0;
        wt_ff_mssm_ttbar_syst_down_1_ = event->Exists("wt_ff_mssm_ttbar_syst_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_syst_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc1_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc2_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_low_unc3_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc1_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc2_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_med_unc3_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc1_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc2_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet0_jet_pt_high_unc3_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc1_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc2_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_low_unc3_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc1_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc2_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_med_unc3_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc1_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc2_down_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_up_1") : 0.0;
        wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_njet1_jet_pt_high_unc3_down_1") : 0.0;

        if(channel_ == channel::tt) {
          wt_ff_mssm_qcd_stat_dR_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_dR_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_dR_unc1_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_dR_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_dR_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_dR_unc1_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_dR_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_dR_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_dR_unc2_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_dR_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_dR_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_dR_unc2_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_pt_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_pt_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_pt_unc1_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_pt_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_pt_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_pt_unc1_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_pt_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_pt_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_pt_unc2_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_pt_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_pt_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_pt_unc2_down_1") : 0.0;
        }

        if(channel_ == channel::et || channel_ == channel::mt) {
          wt_ff_mssm_qcd_stat_ss_njets0_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_ss_njets0_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_ss_njets0_unc1_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_ss_njets0_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_ss_njets0_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_ss_njets0_unc1_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_ss_njets1_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_ss_njets1_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_ss_njets1_unc1_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_ss_njets1_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_ss_njets1_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_ss_njets1_unc1_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_os_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_os_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_os_unc1_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_os_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_os_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_os_unc1_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_ss_njets0_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_ss_njets0_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_ss_njets0_unc2_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_ss_njets0_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_ss_njets0_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_ss_njets0_unc2_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_ss_njets1_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_ss_njets1_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_ss_njets1_unc2_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_ss_njets1_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_ss_njets1_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_ss_njets1_unc2_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_os_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_os_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_os_unc2_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_os_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_os_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_os_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc3_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_low_unc4_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc3_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_med_unc4_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc3_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet0_jet_pt_high_unc4_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc3_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_low_unc4_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc3_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_med_unc4_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc3_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_njet1_jet_pt_high_unc4_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_met_njets0_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_met_njets0_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_met_njets0_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_met_njets0_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_met_njets0_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_met_njets0_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_met_njets1_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_met_njets1_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_met_njets1_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_met_njets1_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_met_njets1_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_met_njets1_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_l_pt_njets0_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_l_pt_njets1_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_extrap_unc1_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_extrap_unc1_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_extrap_unc1_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_extrap_unc1_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_extrap_unc1_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_extrap_unc1_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_met_njets0_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_met_njets0_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_met_njets0_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_met_njets0_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_met_njets0_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_met_njets0_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_met_njets1_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_met_njets1_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_met_njets1_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_met_njets1_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_met_njets1_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_met_njets1_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_l_pt_njets0_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_l_pt_njets1_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_stat_extrap_unc2_up_1_ = event->Exists("wt_ff_mssm_wjets_stat_extrap_unc2_up_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_extrap_unc2_up_1") : 0.0;
          wt_ff_mssm_wjets_stat_extrap_unc2_down_1_ = event->Exists("wt_ff_mssm_wjets_stat_extrap_unc2_down_1") ? event->Get<double>("wt_ff_mssm_wjets_stat_extrap_unc2_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_low_unc1_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_low_unc2_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_low_unc3_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_med_unc1_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_med_unc2_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_med_unc3_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_high_unc1_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_high_unc2_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_jet_pt_high_unc3_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_l_pt_unc1_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_l_pt_unc1_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_l_pt_unc1_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_l_pt_unc1_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_l_pt_unc1_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_l_pt_unc1_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_l_pt_unc2_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_l_pt_unc2_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_l_pt_unc2_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_l_pt_unc2_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_l_pt_unc2_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_l_pt_unc2_down_1") : 0.0;

          wt_ff_mssm_ttbar_stat_met_unc1_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_met_unc1_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_met_unc1_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_met_unc1_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_met_unc1_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_met_unc1_down_1") : 0.0;
          wt_ff_mssm_ttbar_stat_met_unc2_up_1_ = event->Exists("wt_ff_mssm_ttbar_stat_met_unc2_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_met_unc2_up_1") : 0.0;
          wt_ff_mssm_ttbar_stat_met_unc2_down_1_ = event->Exists("wt_ff_mssm_ttbar_stat_met_unc2_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_stat_met_unc2_down_1") : 0.0;
          wt_ff_mssm_wjets_frac_up_1_ = event->Exists("wt_ff_mssm_wjets_frac_up_1") ? event->Get<double>("wt_ff_mssm_wjets_frac_up_1") : 0.0;
          wt_ff_mssm_wjets_frac_down_1_ = event->Exists("wt_ff_mssm_wjets_frac_down_1") ? event->Get<double>("wt_ff_mssm_wjets_frac_down_1") : 0.0;
          wt_ff_mssm_ttbar_frac_up_1_ = event->Exists("wt_ff_mssm_ttbar_frac_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_frac_up_1") : 0.0;
          wt_ff_mssm_ttbar_frac_down_1_ = event->Exists("wt_ff_mssm_ttbar_frac_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_frac_down_1") : 0.0;
          wt_ff_mssm_low_pt_up_1_ = event->Exists("wt_ff_mssm_low_pt_up_1") ? event->Get<double>("wt_ff_mssm_low_pt_up_1") : 0.0;
          wt_ff_mssm_low_pt_down_1_ = event->Exists("wt_ff_mssm_low_pt_down_1") ? event->Get<double>("wt_ff_mssm_low_pt_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_l_pt_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_l_pt_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_l_pt_unc1_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_l_pt_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_l_pt_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_l_pt_unc1_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_l_pt_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_l_pt_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_l_pt_unc2_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_l_pt_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_l_pt_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_l_pt_unc2_down_1") : 0.0;

          wt_ff_mssm_qcd_stat_iso_unc1_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_iso_unc1_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_iso_unc1_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_iso_unc1_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_iso_unc1_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_iso_unc1_down_1") : 0.0;
          wt_ff_mssm_qcd_stat_iso_unc2_up_1_ = event->Exists("wt_ff_mssm_qcd_stat_iso_unc2_up_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_iso_unc2_up_1") : 0.0;
          wt_ff_mssm_qcd_stat_iso_unc2_down_1_ = event->Exists("wt_ff_mssm_qcd_stat_iso_unc2_down_1") ? event->Get<double>("wt_ff_mssm_qcd_stat_iso_unc2_down_1") : 0.0;
          wt_ff_mssm_qcd_syst_iso_up_1_ = event->Exists("wt_ff_mssm_qcd_syst_iso_up_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_iso_up_1") : 0.0;
          wt_ff_mssm_qcd_syst_iso_down_1_ = event->Exists("wt_ff_mssm_qcd_syst_iso_down_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_iso_down_1") : 0.0;
        }

        wt_ff_mssm_syst_alt_func_up_1_ = event->Exists("wt_ff_mssm_syst_alt_func_up_1") ? event->Get<double>("wt_ff_mssm_syst_alt_func_up_1") : 0.0;
        wt_ff_mssm_syst_alt_func_down_1_ = event->Exists("wt_ff_mssm_syst_alt_func_down_1") ? event->Get<double>("wt_ff_mssm_syst_alt_func_down_1") : 0.0;
        wt_ff_mssm_qcd_syst_met_closure_up_1_ = event->Exists("wt_ff_mssm_qcd_syst_met_closure_up_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_met_closure_up_1") : 0.0;
        wt_ff_mssm_qcd_syst_met_closure_down_1_ = event->Exists("wt_ff_mssm_qcd_syst_met_closure_down_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_met_closure_down_1") : 0.0;
        wt_ff_mssm_ttbar_syst_met_closure_up_1_ = event->Exists("wt_ff_mssm_ttbar_syst_met_closure_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_syst_met_closure_up_1") : 0.0;
        wt_ff_mssm_ttbar_syst_met_closure_down_1_ = event->Exists("wt_ff_mssm_ttbar_syst_met_closure_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_syst_met_closure_down_1") : 0.0;
        wt_ff_mssm_ttbar_syst_l_pt_closure_up_1_ = event->Exists("wt_ff_mssm_ttbar_syst_l_pt_closure_up_1") ? event->Get<double>("wt_ff_mssm_ttbar_syst_l_pt_closure_up_1") : 0.0;
        wt_ff_mssm_ttbar_syst_l_pt_closure_down_1_ = event->Exists("wt_ff_mssm_ttbar_syst_l_pt_closure_down_1") ? event->Get<double>("wt_ff_mssm_ttbar_syst_l_pt_closure_down_1") : 0.0;
        wt_ff_mssm_wjets_syst_met_closure_up_1_ = event->Exists("wt_ff_mssm_wjets_syst_met_closure_up_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_met_closure_up_1") : 0.0;
        wt_ff_mssm_wjets_syst_met_closure_down_1_ = event->Exists("wt_ff_mssm_wjets_syst_met_closure_down_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_met_closure_down_1") : 0.0;
        wt_ff_mssm_wjets_syst_l_pt_closure_up_1_ = event->Exists("wt_ff_mssm_wjets_syst_l_pt_closure_up_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_l_pt_closure_up_1") : 0.0;
        wt_ff_mssm_wjets_syst_l_pt_closure_down_1_ = event->Exists("wt_ff_mssm_wjets_syst_l_pt_closure_down_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_l_pt_closure_down_1") : 0.0;
        wt_ff_mssm_wjets_syst_bkg_up_1_ = event->Exists("wt_ff_mssm_wjets_syst_bkg_up_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_bkg_up_1") : 0.0;
        wt_ff_mssm_wjets_syst_bkg_down_1_ = event->Exists("wt_ff_mssm_wjets_syst_bkg_down_1") ? event->Get<double>("wt_ff_mssm_wjets_syst_bkg_down_1") : 0.0;
        wt_ff_mssm_qcd_syst_bkg_up_1_ = event->Exists("wt_ff_mssm_qcd_syst_bkg_up_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_bkg_up_1") : 0.0;
        wt_ff_mssm_qcd_syst_bkg_down_1_ = event->Exists("wt_ff_mssm_qcd_syst_bkg_down_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_bkg_down_1") : 0.0;
        wt_ff_mssm_qcd_syst_pt_2_closure_up_1_ = event->Exists("wt_ff_mssm_qcd_syst_pt_2_closure_up_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_pt_2_closure_up_1") : 0.0;
        wt_ff_mssm_qcd_syst_pt_2_closure_down_1_ = event->Exists("wt_ff_mssm_qcd_syst_pt_2_closure_down_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_pt_2_closure_down_1") : 0.0;
        wt_ff_mssm_qcd_syst_dr_closure_up_1_ = event->Exists("wt_ff_mssm_qcd_syst_dr_closure_up_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_dr_closure_up_1") : 0.0;
        wt_ff_mssm_qcd_syst_dr_closure_down_1_ = event->Exists("wt_ff_mssm_qcd_syst_dr_closure_down_1") ? event->Get<double>("wt_ff_mssm_qcd_syst_dr_closure_down_1") : 0.0;
      }

    }
    
    
   if(do_qcd_scale_wts_){
     // note some of these labels may be generator dependent so need to make sure you check before using them
      
     if(eventInfo->weight_defined("1001")) scale1_ = eventInfo->weight("1001"); else scale1_=1.0;
     if(eventInfo->weight_defined("1002")) scale2_ = eventInfo->weight("1002"); else scale2_=1.0;
     if(eventInfo->weight_defined("1003")) scale3_ = eventInfo->weight("1003"); else scale3_=1.0;
     if(eventInfo->weight_defined("1004")) scale4_ = eventInfo->weight("1004"); else scale4_=1.0;
     if(eventInfo->weight_defined("1005")) scale5_ = eventInfo->weight("1005"); else scale5_=1.0;
     if(eventInfo->weight_defined("1006")) scale6_ = eventInfo->weight("1006"); else scale6_=1.0;
     if(eventInfo->weight_defined("1007")) scale7_ = eventInfo->weight("1007"); else scale7_=1.0;
     if(eventInfo->weight_defined("1008")) scale8_ = eventInfo->weight("1008"); else scale8_=1.0;
     if(eventInfo->weight_defined("1009")) scale9_ = eventInfo->weight("1009"); else scale9_=1.0;    

     // W-jets weights are numbered 1-9 - be careful this doesn't overwrite scale weights for other samples with some other weight!
     if(eventInfo->weight_defined("1")) scale1_ = eventInfo->weight("1"); else scale1_=1.0;
     if(eventInfo->weight_defined("2")) scale2_ = eventInfo->weight("2"); else scale2_=1.0;
     if(eventInfo->weight_defined("3")) scale3_ = eventInfo->weight("3"); else scale3_=1.0;
     if(eventInfo->weight_defined("4")) scale4_ = eventInfo->weight("4"); else scale4_=1.0;
     if(eventInfo->weight_defined("5")) scale5_ = eventInfo->weight("5"); else scale5_=1.0;
     if(eventInfo->weight_defined("6")) scale6_ = eventInfo->weight("6"); else scale6_=1.0;
     if(eventInfo->weight_defined("7")) scale7_ = eventInfo->weight("7"); else scale7_=1.0;
     if(eventInfo->weight_defined("8")) scale8_ = eventInfo->weight("8"); else scale8_=1.0;
     if(eventInfo->weight_defined("9")) scale9_ = eventInfo->weight("9"); else scale9_=1.0;
   }
    
    std::vector<PileupInfo *> puInfo;
    float true_int = -1;

    if (event->Exists("pileupInfo") && !is_data_ && !is_embedded_ ) {
     puInfo = event->GetPtrVec<PileupInfo>("pileupInfo");
      for (unsigned i = 0; i < puInfo.size(); ++i) {
        if (puInfo[i]->bunch_crossing() == 0)
          true_int = puInfo[i]->true_num_interactions();
      }
    }
    n_pu_ = true_int;
    rho_ = eventInfo->jet_rho();
    if(event->Exists("gen_match_1")) gen_match_1_ = MCOrigin2UInt(event->Get<ic::mcorigin>("gen_match_1"));
    if(event->Exists("gen_match_2")) gen_match_2_ = MCOrigin2UInt(event->Get<ic::mcorigin>("gen_match_2"));
    if(event->Exists("gen_match_1_pt")) gen_match_1_pt_ = event->Get<double>("gen_match_1_pt");
    if(event->Exists("gen_match_2_pt")) gen_match_2_pt_ = event->Get<double>("gen_match_2_pt");

    if(event->Exists("gen_match_undecayed_1_pt")) gen_match_undecayed_1_pt_ = event->Get<double>("gen_match_undecayed_1_pt");
    if(event->Exists("gen_match_undecayed_2_pt")) gen_match_undecayed_2_pt_ = event->Get<double>("gen_match_undecayed_2_pt");

    if(event->Exists("ngenjets")) ngenjets_ = event->Get<unsigned>("ngenjets");
    if(event->Exists("ngenjets20")) ngenjets20_ = event->Get<unsigned>("ngenjets20");
    if(event->Exists("genjetpt_1")) genjetpt_1_ = event->Get<double>("genjetpt_1");
    if(event->Exists("genjetpt_2")) genjetpt_2_ = event->Get<double>("genjetpt_2");
    if(event->Exists("genjeteta_1")) genjeteta_1_ = event->Get<double>("genjeteta_1");
    if(event->Exists("genjeteta_2")) genjeteta_2_ = event->Get<double>("genjeteta_2");
    if(event->Exists("gen_sjdphi")) gen_sjdphi_ = event->Get<double>("gen_sjdphi");
    if(event->Exists("tauFlag1")) tauFlag_1_ = event->Get<int>("tauFlag1");
    if(event->Exists("tauFlag2")) tauFlag_2_ = event->Get<int>("tauFlag2");
   
    wt_tau_fake_up_ = 1.0;
    wt_tau_fake_down_ = 1.0;
    wt_tquark_up_ = 1.0;
    wt_tquark_down_ = 1.0;
    wt_tquark_alt_ = 1.0;
    wt_zpt_up_ = 1.0;
    wt_zpt_up_NLO_ = 1.0;
    wt_zpt_down_ = 1.0;
    wt_zpt_down_NLO_ = 1.0;
    wt_zpt_embed_ic_ = 1.0;
    wt_em_qcd_ = 1.0;
    wt_efake_rate_up_ = 1.0;
    wt_efake_rate_down_ = 1.0;
    wt_mfake_rate_up_ = 1.0;
    wt_mfake_rate_down_ = 1.0;
    if (event->Exists("wt_tau_fake_up"))    wt_tau_fake_up_   = event->Get<double>("wt_tau_fake_up");
    if (event->Exists("wt_tau_fake_down"))  wt_tau_fake_down_ = event->Get<double>("wt_tau_fake_down");
    if (event->Exists("wt_tquark_up"))      wt_tquark_up_   = event->Get<double>("wt_tquark_up");
    if (event->Exists("wt_tquark_down"))    wt_tquark_down_ = event->Get<double>("wt_tquark_down");
    if (event->Exists("wt_tquark_alt"))      wt_tquark_alt_   = event->Get<double>("wt_tquark_alt");
    if (event->Exists("wt_zpt_up"))         wt_zpt_up_   = event->Get<double>("wt_zpt_up");
    if (event->Exists("wt_zpt_down"))       wt_zpt_down_ = event->Get<double>("wt_zpt_down");
    if (event->Exists("wt_zpt_up_NLO"))         wt_zpt_up_NLO_   = event->Get<double>("wt_zpt_up_NLO");
    if (event->Exists("wt_zpt_down_NLO"))       wt_zpt_down_NLO_ = event->Get<double>("wt_zpt_down_NLO");
    if (event->Exists("wt_zpt_embed_ic"))      wt_zpt_embed_ic_ = event->Get<double>("wt_zpt_embed_ic");
    if (event->Exists("wt_efake_rate_up"))  wt_efake_rate_up_   = event->Get<double>("wt_efake_rate_up");
    if (event->Exists("wt_efake_rate_down")) wt_efake_rate_down_ = event->Get<double>("wt_efake_rate_down");
    if (event->Exists("wt_mfake_rate_up"))  wt_mfake_rate_up_   = event->Get<double>("wt_mfake_rate_up");
    if (event->Exists("wt_mfake_rate_down")) wt_mfake_rate_down_ = event->Get<double>("wt_mfake_rate_down");
    if (event->Exists("wt_em_qcd"))         wt_em_qcd_ = event->Get<double>("wt_em_qcd");
    if (event->Exists("wt_em_qcd_extrapup"))      wt_em_qcd_extrapup_ = event->Get<double>("wt_em_qcd_extrapup");
    if (event->Exists("wt_em_qcd_extrapdown"))    wt_em_qcd_extrapdown_ = event->Get<double>("wt_em_qcd_extrapdown");

    wt_em_qcd_njets0_unc1_up_   = (event->Exists("wt_em_qcd_njets0_unc1_up"  )) ? event->Get<double>("wt_em_qcd_njets0_unc1_up"  ) : 1.0; 
    wt_em_qcd_njets0_unc1_down_ = (event->Exists("wt_em_qcd_njets0_unc1_down")) ? event->Get<double>("wt_em_qcd_njets0_unc1_down") : 1.0; 
    wt_em_qcd_njets0_unc2_up_   = (event->Exists("wt_em_qcd_njets0_unc2_up"  )) ? event->Get<double>("wt_em_qcd_njets0_unc2_up"  ) : 1.0; 
    wt_em_qcd_njets0_unc2_down_ = (event->Exists("wt_em_qcd_njets0_unc2_down")) ? event->Get<double>("wt_em_qcd_njets0_unc2_down") : 1.0; 
    wt_em_qcd_njets1_unc1_up_   = (event->Exists("wt_em_qcd_njets1_unc1_up"  )) ? event->Get<double>("wt_em_qcd_njets1_unc1_up"  ) : 1.0; 
    wt_em_qcd_njets1_unc1_down_ = (event->Exists("wt_em_qcd_njets1_unc1_down")) ? event->Get<double>("wt_em_qcd_njets1_unc1_down") : 1.0; 
    wt_em_qcd_njets1_unc2_up_   = (event->Exists("wt_em_qcd_njets1_unc2_up"  )) ? event->Get<double>("wt_em_qcd_njets1_unc2_up"  ) : 1.0; 
    wt_em_qcd_njets1_unc2_down_ = (event->Exists("wt_em_qcd_njets1_unc2_down")) ? event->Get<double>("wt_em_qcd_njets1_unc2_down") : 1.0; 
    wt_em_qcd_njets2_unc1_up_   = (event->Exists("wt_em_qcd_njets2_unc1_up"  )) ? event->Get<double>("wt_em_qcd_njets2_unc1_up"  ) : 1.0; 
    wt_em_qcd_njets2_unc1_down_ = (event->Exists("wt_em_qcd_njets2_unc1_down")) ? event->Get<double>("wt_em_qcd_njets2_unc1_down") : 1.0; 
    wt_em_qcd_njets2_unc2_up_   = (event->Exists("wt_em_qcd_njets2_unc2_up"  )) ? event->Get<double>("wt_em_qcd_njets2_unc2_up"  ) : 1.0; 
    wt_em_qcd_njets2_unc2_down_ = (event->Exists("wt_em_qcd_njets2_unc2_down")) ? event->Get<double>("wt_em_qcd_njets2_unc2_down") : 1.0;

    
    if (!is_embedded_ && eventInfo->weight_defined("pileup")) pu_weight_ = eventInfo->weight("pileup"); else pu_weight_ = 0.0;
    if (event->Exists("trigweight_1")) trigweight_1_ = event->Get<double>("trigweight_1"); else trigweight_1_ = 1.0;
    if (event->Exists("trigweight_2")) trigweight_2_ = event->Get<double>("trigweight_2"); else trigweight_2_ = 1.0;
    if (event->Exists("tau_leg_SF")) tau_leg_SF_ = event->Get<double>("tau_leg_SF"); else tau_leg_SF_ = 1.0;
    if (event->Exists("xtrg_sf")) xtrg_sf_ = event->Get<double>("xtrg_sf"); else xtrg_sf_ = 1.0;
    if (event->Exists("single_l_sf")) single_l_sf_ = event->Get<double>("single_l_sf"); else single_l_sf_ = 1.0;
    if (event->Exists("xtrg_notrig")) xtrg_notrig_ = event->Get<double>("xtrg_notrig"); else xtrg_notrig_ = 1.0;
    if (event->Exists("OR_notrig")) OR_notrig_ = event->Get<double>("OR_notrig"); else OR_notrig_ = 1.0; 
    if (event->Exists("idisoweight_1")) idisoweight_1_ = event->Get<double>("idisoweight_1"); else idisoweight_1_ = 0.0;
    if (event->Exists("idisoweight_2")) idisoweight_2_ = event->Get<double>("idisoweight_2"); else idisoweight_2_ = 1.0;
    if (event->Exists("trackingweight_1")) trackingweight_1_ = event->Get<double>("trackingweight_1"); else trackingweight_1_ = 1.0;
    if (event->Exists("trackingweight_2")) trackingweight_2_ = event->Get<double>("trackingweight_2"); else trackingweight_2_ = 1.0;
    if (event->Exists("wt_emb_sel_kit")) wt_emb_sel_kit_ = event->Get<double>("wt_emb_sel_kit"); else wt_emb_sel_kit_ = 1.0;
    if (event->Exists("wt_emb_sel_trg_kit")) wt_emb_sel_trg_kit_ = event->Get<double>("wt_emb_sel_trg_kit"); else wt_emb_sel_trg_kit_ = 1.0;
    if (event->Exists("wt_embedding_yield")) wt_embedding_yield_ = event->Get<double>("wt_embedding_yield"); else wt_embedding_yield_ = 0.0;
    if (eventInfo->weight_defined("lepton")) effweight_ = eventInfo->weight("lepton"); else effweight_ = 0.0;
    et_trg_= (event->Exists("et_trg_single")) ? event->Get<double>("et_trg_single") : 0;
    et_trg_cross_= (event->Exists("et_trg_cross")) ? event->Get<double>("et_trg_cross") : 0;
    et_trg_or_= (event->Exists("et_trg_or")) ? event->Get<double>("et_trg_or") : 0;


    if (eventInfo->weight_defined("wt_embedding")) generatorWeight_ = std::min(eventInfo->weight("wt_embedding"),1.0); else generatorWeight_ = wt_.var_double/fabs(wt_.var_double);
    if (eventInfo->weight_defined("muonEffTrgWeight")) muonEffTrgWeight_= eventInfo->weight("muonEffTrgWeight"); else muonEffTrgWeight_ = 0.0;
    if (eventInfo->weight_defined("muonEffIDWeight_1")) muonEffIDWeight_1_ = eventInfo->weight("muonEffIDWeight_1"); else muonEffIDWeight_1_ = 0.0;
    if (eventInfo->weight_defined("muonEffIDWeight_2")) muonEffIDWeight_2_ = eventInfo->weight("muonEffIDWeight_2"); else muonEffIDWeight_2_ = 0.0;
    if (eventInfo->weight_defined("et_triggerweight_ic")) et_triggerweight_ic_ = eventInfo->weight("et_triggerweight_ic"); else et_triggerweight_ic_ = 0.0;
    if (eventInfo->weight_defined("mt_triggerweight_ic")) mt_triggerweight_ic_ = eventInfo->weight("mt_triggerweight_ic"); else mt_triggerweight_ic_ = 0.0;
    if (eventInfo->weight_defined("tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2")) tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2_ = eventInfo->weight("tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2"); else tauIDScaleFactorWeight_medium_DeepTau2017v2p1VSjet_2_ = 0.0;
    if (eventInfo->weight_defined("tauIDScaleFactorWeight_highpt_deeptauid_2")) tauIDScaleFactorWeight_highpt_deeptauid_2_ = eventInfo->weight("tauIDScaleFactorWeight_highpt_deeptauid_2"); else tauIDScaleFactorWeight_highpt_deeptauid_2_ = 0.0;


    mvadm_idiso_et_ =  (event->Exists("wt_tau_id_mvadm_sync")) ? event->Get<double>("wt_tau_id_mvadm_sync") : 1.0; 
    std::vector<CompositeCandidate *> const& ditau_vec = event->GetPtrVec<CompositeCandidate>(ditau_label_);
    CompositeCandidate const* ditau = ditau_vec.at(0);
    Candidate const* lep1 = ditau->GetCandidate("lepton1");
    Candidate const* lep2 = ditau->GetCandidate("lepton2");
    Met const* mets = NULL;
    mets = event->GetPtr<Met>(met_label_);

    std::vector<PFJet*> jets = event->GetPtrVec<PFJet>(jets_label_);
    std::vector<PFJet*> uncleaned_jets = event->GetPtrVec<PFJet>("ak4PFJetsCHSUnFiltered");
    std::sort(jets.begin(), jets.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
    std::vector<PFJet*> lowpt_jets = jets;
    ic::erase_if(jets,!boost::bind(MinPtMaxEta, _1, 30.0, 4.7));
    ic::erase_if(lowpt_jets,!boost::bind(MinPtMaxEta, _1, 20.0, 4.7));
    std::vector<PFJet*> prebjets = lowpt_jets;

    double eta_cut = 2.4;
    if(era_ == era::data_2017 || era_ == era::data_2017UL || era_ == era::data_2018 || era_ == era::data_2018UL) eta_cut = 2.5;
    ic::erase_if(prebjets,!boost::bind(MinPtMaxEta, _1, 20.0, eta_cut));
    n_prebjets_ = prebjets.size();
    std::vector<PFJet*> bjets = prebjets;
    std::vector<PFJet*> loose_bjets = prebjets;
    // NEW
    std::vector<PFJet*> btag = prebjets;
    std::vector<PFJet*> loose_btag = prebjets;
    //
    std::vector<PFJet*> deepbjets = prebjets;
    std::string btag_label="s";
    std::string btag_label_extra ="";
    double btag_wp;
    double loose_btag_wp;
    double deepjet_wp;
    if (era_ == era::data_2017 || era_ == era::data_2017UL) {
      btag_wp = 0.4941;
      loose_btag_wp = 0.1522;
      btag_label = "pfDeepCSVJetTags:probb";
      btag_label_extra = "pfDeepCSVJetTags:probbb";

      deepjet_wp = 0.3033;
    }
    if (era_ == era::data_2018) {
      btag_wp = 0.4184;
      loose_btag_wp = 0.1241;
      btag_label = "pfDeepCSVJetTags:probb";
      btag_label_extra = "pfDeepCSVJetTags:probbb";

      deepjet_wp = 0.2770;
    }
    if (era_ == era::data_2018UL) {
      btag_wp = 0.4184;
      loose_btag_wp = 0.1241;
      btag_label = "pfDeepCSVJetTags:probb";
      btag_label_extra = "pfDeepCSVJetTags:probbb";

      deepjet_wp = 0.2783;
    }
    if (era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) {
      btag_wp = 0.6321;
      loose_btag_wp = 0.2217;
      btag_label = "pfDeepCSVJetTags:probb";
      btag_label_extra = "pfDeepCSVJetTags:probbb";

      deepjet_wp = 0.3093;
    }

    std::string deepjet_label_1 = "pfDeepFlavourJetTags:probb"; 
    std::string deepjet_label_2 = "pfDeepFlavourJetTags:probbb"; 
    std::string deepjet_label_3 = "pfDeepFlavourJetTags:problepb";

    auto filterDeepJetBTag = [deepjet_label_1, deepjet_label_2, deepjet_label_3, deepjet_wp] (PFJet* s1) -> bool {
      return s1->GetBDiscriminator(deepjet_label_1) + s1->GetBDiscriminator(deepjet_label_2) + s1->GetBDiscriminator(deepjet_label_3) > deepjet_wp;
    };


    auto sortRuleBTagSum = [btag_label, btag_label_extra] (PFJet* s1, PFJet* s2) -> bool {
      return s1->GetBDiscriminator(btag_label) + s1->GetBDiscriminator(btag_label_extra) > 
        s2->GetBDiscriminator(btag_label) + s2->GetBDiscriminator(btag_label_extra);
    };
    auto filterBTagSumTight = [btag_label, btag_label_extra, btag_wp] (PFJet* s1) -> bool {
      return s1->GetBDiscriminator(btag_label) + s1->GetBDiscriminator(btag_label_extra) > btag_wp;
    };
    auto filterBTagSumLoose = [btag_label, btag_label_extra, loose_btag_wp] (PFJet* s1) -> bool {
      return s1->GetBDiscriminator(btag_label) + s1->GetBDiscriminator(btag_label_extra) > loose_btag_wp;
    };
    ic::erase_if_not(loose_btag, filterBTagSumLoose);
    ic::erase_if_not(btag, filterBTagSumTight);

    ic::erase_if_not(loose_bjets, filterBTagSumLoose);
    // Instead of changing b-tag value in the promote/demote method we look for a map of bools
    // that say whether a jet should pass the WP or not
    if (event->Exists("retag_result") && false) { // set false as we use deepjet now
      auto const& retag_result = event->Get<std::map<std::size_t,bool>>("retag_result"); 
      ic::erase_if(bjets, !boost::bind(IsReBTagged, _1, retag_result));
    } else{ 
      ic::erase_if_not(bjets, filterBTagSumTight);
    }

    if (event->Exists("retag_result")) {
      auto const& retag_result = event->Get<std::map<std::size_t,bool>>("retag_result");
      ic::erase_if(deepbjets, !boost::bind(IsReBTagged, _1, retag_result));
    } else{
      ic::erase_if_not(deepbjets, filterDeepJetBTag);
    } 

    n_deepbjets_ = deepbjets.size();
    deepbpt_1_=-9999;
    deepbeta_1_=-9999;
    deep_probb_1_=-9999;
    deep_probbb_1_=-9999;
    deep_problepb_1_=-9999;
    deepcsv_probb_1_=-9999;
    deepcsv_probbb_1_=-9999;
    if(prebjets.size()>0) {
      if(n_deepbjets_>0){
        deepbpt_1_=prebjets[0]->pt();
        deepbeta_1_=prebjets[0]->eta();
      }
      deep_probb_1_ = prebjets[0]->GetBDiscriminator(deepjet_label_1);
      deep_probbb_1_ = prebjets[0]->GetBDiscriminator(deepjet_label_2);
      deep_problepb_1_ = prebjets[0]->GetBDiscriminator(deepjet_label_3);

      deepcsv_probb_1_ = prebjets[0]->GetBDiscriminator("pfDeepCSVJetTags:probb");
      deepcsv_probbb_1_ = prebjets[0]->GetBDiscriminator("pfDeepCSVJetTags:probbb");
    }


    // Btag weights
    wt_btag_           = event->Exists("btag_evt_weight") ? event->Get<double>("btag_evt_weight") : 1.;
    wt_btag_down_      = event->Exists("btag_evt_weight_down") ? event->Get<double>("btag_evt_weight_down") : 1.;
    wt_btag_up_        = event->Exists("btag_evt_weight_up") ? event->Get<double>("btag_evt_weight_up") : 1.;

    wt_btag_real_down_ = event->Exists("btag_evt_weight_realbtag_down") ? event->Get<double>("btag_evt_weight_realbtag_down") : 1.;
    wt_btag_real_up_   = event->Exists("btag_evt_weight_realbtag_up")   ? event->Get<double>("btag_evt_weight_realbtag_up") : 1.;
    wt_btag_fake_down_ = event->Exists("btag_evt_weight_fakebtag_down") ? event->Get<double>("btag_evt_weight_fakebtag_down") : 1.;
    wt_btag_fake_up_   = event->Exists("btag_evt_weight_fakebtag_up")   ? event->Get<double>("btag_evt_weight_fakebtag_up") : 1.;

    if (PairOppSign(ditau)) {
      os_ = true;
    } else {
      os_ = false;
    }

    q1_ = ditau->At(0)->charge();
    q2_ = ditau->At(1)->charge();

    //Fill extra lepton veto bools
    dilepton_veto_ = false;
    extraelec_veto_ = false;
    extramuon_veto_ = false;
    bool dilep_veto_=false;
    if(channel_ == channel::et) { 
        if(event->Exists("dielec_veto")) dilepton_veto_ = event->Get<bool>("dielec_veto");
        if(event->Exists("dimuon_veto")) dilep_veto_ = event->Get<bool>("dimuon_veto");
        if(event->Exists("extra_elec_veto")) extraelec_veto_ = event->Get<bool>("extra_elec_veto");
        if(event->Exists("extra_muon_veto")) extramuon_veto_ = event->Get<bool>("extra_muon_veto");
    }
    if(channel_ == channel::mt) { 
        if(event->Exists("dimuon_veto")) dilepton_veto_ = event->Get<bool>("dimuon_veto");
        if(event->Exists("dielec_veto")) dilep_veto_ = event->Get<bool>("dielec_veto");
        if(event->Exists("extra_elec_veto")) extraelec_veto_ = event->Get<bool>("extra_elec_veto");
        if(event->Exists("extra_muon_veto")) extramuon_veto_ = event->Get<bool>("extra_muon_veto");
    }
    if(channel_ == channel::em) { 
        if(event->Exists("extra_elec_veto")) extraelec_veto_ = event->Get<bool>("extra_elec_veto");
        if(event->Exists("extra_muon_veto")) extramuon_veto_ = event->Get<bool>("extra_muon_veto");
    }
    if(channel_ == channel::tt) {
        if(event->Exists("extra_elec_veto")) extraelec_veto_ = event->Get<bool>("extra_elec_veto");
        if(event->Exists("extra_muon_veto")) extramuon_veto_ = event->Get<bool>("extra_muon_veto");
    }
    if(channel_ == channel::zee || channel_ == channel::zmm) {
        if(event->Exists("extra_elec_veto")) extraelec_veto_ = event->Get<bool>("extra_elec_veto");
        if(event->Exists("extra_muon_veto")) extramuon_veto_ = event->Get<bool>("extra_muon_veto");
    }  

    lepton_veto_ = dilepton_veto_ || extraelec_veto_ || extramuon_veto_;
    
    if((strategy_==strategy::smsummer16 || strategy_ == strategy::cpsummer16 || strategy_ == strategy::legacy16 || strategy_ == strategy::cpdecays16 || strategy_ == strategy::cpsummer17 || strategy_ == strategy::cpdecays17 || strategy_ == strategy::cpdecays18) && !make_sync_ntuple_) dilepton_veto_ = dilep_veto_ || dilepton_veto_;

    n_vtx_ = eventInfo->good_vertices();
    if (event->Exists("svfitMass")) {
      m_sv_ = event->Get<double>("svfitMass");
    } else {
      m_sv_ = -9999;
    }

   if (event->Exists("svfitMT")){
     mt_sv_ = event->Get<double>("svfitMT");
   } else {
     mt_sv_ = -9999;
   }

    if (event->Exists("svfitHiggs")) {
      Candidate const& higgs = event->Get<Candidate>("svfitHiggs");
      pt_h_ = higgs.pt();
      eta_h_ = higgs.eta();
      phi_h_ = higgs.phi();
    } else {
      pt_h_ = -9999;
      eta_h_ = -9999;
      phi_h_ = -9999;
    }

    pt_tt_ = (ditau->vector() + mets->vector()).pt();

    Met* newmet = new Met();
    Candidate * newlep = new Candidate();

    ROOT::Math::PtEtaPhiEVector newvec;

    newvec.SetPxPyPzE(lep1->vector().Px(), lep1->vector().Py(),0.,lep1->pt());

    newlep->set_vector(newvec);
    newmet->set_vector(mets->vector()+newlep->vector());

    newmet_ = newmet->pt();
    newmet_dphi_1_=std::fabs(ROOT::Math::VectorUtil::DeltaPhi(newmet->vector(),lep1->vector()));
    newmet_dphi_2_=std::fabs(ROOT::Math::VectorUtil::DeltaPhi(newmet->vector(),lep2->vector()));   
    ROOT::Math::PtEtaPhiEVector fake_met_vec;  
    event->Exists("fake_met_vec") ? fake_met_vec = event->Get<ROOT::Math::PtEtaPhiEVector>("fake_met_vec") : fake_met_vec;
    fake_met_dphi_1_=std::fabs(ROOT::Math::VectorUtil::DeltaPhi(fake_met_vec,lep1->vector()));   
    fake_met_dphi_2_=std::fabs(ROOT::Math::VectorUtil::DeltaPhi(fake_met_vec,lep2->vector()));   


    ROOT::Math::PtEtaPhiEVector fake_tau_met_vec;
    event->Exists("fake_tau_met_vec") ? fake_tau_met_vec = event->Get<ROOT::Math::PtEtaPhiEVector>("fake_tau_met_vec") : fake_tau_met_vec;
    fake_tau_met_dphi_1_=std::fabs(ROOT::Math::VectorUtil::DeltaPhi(fake_tau_met_vec,lep1->vector()));
    fake_tau_met_dphi_2_=std::fabs(ROOT::Math::VectorUtil::DeltaPhi(fake_tau_met_vec,lep2->vector()));

    ROOT::Math::PtEtaPhiEVector gen_tau_met_vec;
    event->Exists("gen_tau_met_vec") ? gen_tau_met_vec = event->Get<ROOT::Math::PtEtaPhiEVector>("gen_tau_met_vec") : gen_tau_met_vec;
 
    ROOT::Math::PtEtaPhiEVector gen_met_vec;
    event->Exists("gen_met_vec") ? gen_met_vec = event->Get<ROOT::Math::PtEtaPhiEVector>("gen_met_vec") : gen_met_vec;

    fake_tau_metx_ = fake_tau_met_vec.Px();
    fake_tau_mety_ = fake_tau_met_vec.Py();
    gen_tau_metx_ = gen_tau_met_vec.Px();
    gen_tau_mety_ = gen_tau_met_vec.Py();
    gen_metx_ = gen_met_vec.Px();
    gen_mety_ = gen_met_vec.Py();

 
    event->Exists("fake_met") ? fake_met_ = event->Get<double>("fake_met") : 0.;
    event->Exists("gen_met") ? gen_met_ = event->Get<double>("gen_met") : 0.;

    event->Exists("fake_tau_met") ? fake_tau_met_ = event->Get<double>("fake_tau_met") : 0.;
    event->Exists("gen_tau_met") ? gen_tau_met_ = event->Get<double>("gen_tau_met") : 0.;

    pt_tt_inc_met_=-9999;
    if(channel_ == channel::zmm || channel_ == channel::zee) pt_tt_ = (ditau->vector()).pt(); 
    if(channel_ == channel::zmm || channel_ == channel::zee) pt_tt_inc_met_ = (ditau->vector()+mets->vector()).pt(); 
    m_vis_ = ditau->M();
    pt_vis_ = ditau->pt();

    mt_lep_ = MT(lep1,lep2);
    mt_ll_ = MT(ditau, mets);
    pzeta_ = PZeta(ditau, mets, 0.85);
    pzetamiss_ = PZeta(ditau, mets, 0.0);
    pzetavis_ = PZetaVis(ditau);

    met_dphi_1_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(mets->vector(),lep1->vector()));
    met_dphi_2_ = std::fabs(ROOT::Math::VectorUtil::DeltaPhi(mets->vector(),lep2->vector()));

    
    


    mt_1_ = MT(lep1, mets);
    mt_2_ = MT(lep2, mets);
    mt_tot_ = sqrt(pow(mt_lep_.var_double,2)+pow(mt_2_.var_double,2)+pow(mt_1_.var_double,2));

    // std::cout << lep1->pt() << "    " << lep2->pt() << "    " << lep1->eta() << "    " << lep2->eta() << std::endl;

    pt_1_ = lep1->pt();
    pt_2_ = lep2->pt();
    eta_1_ = lep1->eta();
    eta_2_ = lep2->eta();
    phi_1_ = lep1->phi();
    phi_2_ = lep2->phi();
    dphi_ = ROOT::Math::VectorUtil::DeltaPhi(lep1->vector(),lep2->vector());
    dR_ = std::fabs(ROOT::Math::VectorUtil::DeltaR(lep1->vector(),lep2->vector()));
    E_1_ = lep1->energy();
    E_2_ = lep2->energy();
    m_1_ = lep1->M();
    m_2_ = lep2->M();
    q_1_ = lep1->charge();
    q_2_ = lep2->charge();



    if(make_sync_ntuple_){
      event->Exists("genpX") ? gen_px_ = event->Get<double>("genpX") : 0.;
      event->Exists("genpY") ? gen_py_ = event->Get<double>("genpY") : 0.;
      event->Exists("vispX") ? vis_px_ = event->Get<double>("vispX") : 0.;
      event->Exists("vispY") ? vis_py_ = event->Get<double>("vispY") : 0.;
    }
    met_ = mets->vector().pt();
    met_phi_ = mets->vector().phi();
    event->Exists("genM") ? gen_m_ = event->Get<double>("genM") : 0.;
    event->Exists("genpT") ? gen_pt_ = event->Get<double>("genpT") : 0.;

    gen_nu_p_1_=-9999.;
    gen_nu_p_2_=-9999.;
    gen_nu_px_1_=-9999.;
    gen_nu_px_2_=-9999.;
    gen_nu_py_1_=-9999.;
    gen_nu_py_2_=-9999.;
    gen_nu_pz_1_=-9999.;
    gen_nu_pz_2_=-9999.;
    gen_nu_phi_1_=-9999.;
    gen_nu_phi_2_=-9999.;
    gen_nu_eta_1_=-9999.;
    gen_nu_eta_2_=-9999.;

    gen_vis_p_1_=-9999.;
    gen_vis_p_2_=-9999.;
    gen_vis_px_1_=-9999.;
    gen_vis_px_2_=-9999.;
    gen_vis_py_1_=-9999.;
    gen_vis_py_2_=-9999.;
    gen_vis_pz_1_=-9999.;
    gen_vis_pz_2_=-9999.;
    gen_vis_E_1_=-9999.;
    gen_vis_E_2_=-9999.;
    gen_vis_phi_1_=-9999.;
    gen_vis_phi_2_=-9999.;
    gen_vis_eta_1_=-9999.;
    gen_vis_eta_2_=-9999.;

    gen_neutral_p_1_=-9999.;
    gen_neutral_p_2_=-9999.;
    gen_neutral_px_1_=-9999.;
    gen_neutral_px_2_=-9999.;
    gen_neutral_py_1_=-9999.;
    gen_neutral_py_2_=-9999.;
    gen_neutral_pz_1_=-9999.;
    gen_neutral_pz_2_=-9999.;
    gen_neutral_E_1_=-9999.;
    gen_neutral_E_2_=-9999.;
    gen_neutral_phi_1_=-9999.;
    gen_neutral_phi_2_=-9999.;
    gen_neutral_eta_1_=-9999.;
    gen_neutral_eta_2_=-9999.;

    event->Exists("gen_nu_p_1")   ? gen_nu_p_1_   = event->Get<double>("gen_nu_p_1")   : -9999.;
    event->Exists("gen_nu_p_2")   ? gen_nu_p_2_   = event->Get<double>("gen_nu_p_2")   : -9999.;
    event->Exists("gen_nu_px_1")   ? gen_nu_px_1_   = event->Get<double>("gen_nu_px_1")   : -9999.;
    event->Exists("gen_nu_px_2")   ? gen_nu_px_2_   = event->Get<double>("gen_nu_px_2")   : -9999.;
    event->Exists("gen_nu_py_1")   ? gen_nu_py_1_   = event->Get<double>("gen_nu_py_1")   : -9999.;
    event->Exists("gen_nu_py_2")   ? gen_nu_py_2_   = event->Get<double>("gen_nu_py_2")   : -9999.;
    event->Exists("gen_nu_pz_1")   ? gen_nu_pz_1_   = event->Get<double>("gen_nu_pz_1")   : -9999.;
    event->Exists("gen_nu_pz_2")   ? gen_nu_pz_2_   = event->Get<double>("gen_nu_pz_2")   : -9999.;
    event->Exists("gen_nu_phi_1") ? gen_nu_phi_1_ = event->Get<double>("gen_nu_phi_1") : -9999.;
    event->Exists("gen_nu_phi_2") ? gen_nu_phi_2_ = event->Get<double>("gen_nu_phi_2") : -9999.;
    event->Exists("gen_nu_eta_1") ? gen_nu_eta_1_ = event->Get<double>("gen_nu_eta_1") : -9999.;
    event->Exists("gen_nu_eta_2") ? gen_nu_eta_2_ = event->Get<double>("gen_nu_eta_2") : -9999.;

    event->Exists("gen_vis_p_1")   ? gen_vis_p_1_   = event->Get<double>("gen_vis_p_1")   : -9999.;
    event->Exists("gen_vis_p_2")   ? gen_vis_p_2_   = event->Get<double>("gen_vis_p_2")   : -9999.;
    event->Exists("gen_vis_px_1")   ? gen_vis_px_1_   = event->Get<double>("gen_vis_px_1")   : -9999.;
    event->Exists("gen_vis_px_2")   ? gen_vis_px_2_   = event->Get<double>("gen_vis_px_2")   : -9999.;
    event->Exists("gen_vis_py_1")   ? gen_vis_py_1_   = event->Get<double>("gen_vis_py_1")   : -9999.;
    event->Exists("gen_vis_py_2")   ? gen_vis_py_2_   = event->Get<double>("gen_vis_py_2")   : -9999.;
    event->Exists("gen_vis_pz_1")   ? gen_vis_pz_1_   = event->Get<double>("gen_vis_pz_1")   : -9999.;
    event->Exists("gen_vis_pz_2")   ? gen_vis_pz_2_   = event->Get<double>("gen_vis_pz_2")   : -9999.;
    event->Exists("gen_vis_E_1")   ? gen_vis_E_1_   = event->Get<double>("gen_vis_E_1")   : -9999.;
    event->Exists("gen_vis_E_2")   ? gen_vis_E_2_   = event->Get<double>("gen_vis_E_2")   : -9999.;
    event->Exists("gen_vis_phi_1") ? gen_vis_phi_1_ = event->Get<double>("gen_vis_phi_1") : -9999.;
    event->Exists("gen_vis_phi_2") ? gen_vis_phi_2_ = event->Get<double>("gen_vis_phi_2") : -9999.;
    event->Exists("gen_vis_eta_1") ? gen_vis_eta_1_ = event->Get<double>("gen_vis_eta_1") : -9999.;
    event->Exists("gen_vis_eta_2") ? gen_vis_eta_2_ = event->Get<double>("gen_vis_eta_2") : -9999.;


    event->Exists("gen_neutral_p_1")   ? gen_neutral_p_1_   = event->Get<double>("gen_neutral_p_1")   : -9999.;
    event->Exists("gen_neutral_p_2")   ? gen_neutral_p_2_   = event->Get<double>("gen_neutral_p_2")   : -9999.;
    event->Exists("gen_neutral_px_1")   ? gen_neutral_px_1_   = event->Get<double>("gen_neutral_px_1")   : -9999.;
    event->Exists("gen_neutral_px_2")   ? gen_neutral_px_2_   = event->Get<double>("gen_neutral_px_2")   : -9999.;
    event->Exists("gen_neutral_py_1")   ? gen_neutral_py_1_   = event->Get<double>("gen_neutral_py_1")   : -9999.;
    event->Exists("gen_neutral_py_2")   ? gen_neutral_py_2_   = event->Get<double>("gen_neutral_py_2")   : -9999.;
    event->Exists("gen_neutral_pz_1")   ? gen_neutral_pz_1_   = event->Get<double>("gen_neutral_pz_1")   : -9999.;
    event->Exists("gen_neutral_pz_2")   ? gen_neutral_pz_2_   = event->Get<double>("gen_neutral_pz_2")   : -9999.;
    event->Exists("gen_neutral_E_1")   ? gen_neutral_E_1_   = event->Get<double>("gen_neutral_E_1")   : -9999.;
    event->Exists("gen_neutral_E_2")   ? gen_neutral_E_2_   = event->Get<double>("gen_neutral_E_2")   : -9999.;
    event->Exists("gen_neutral_phi_1") ? gen_neutral_phi_1_ = event->Get<double>("gen_neutral_phi_1") : -9999.;
    event->Exists("gen_neutral_phi_2") ? gen_neutral_phi_2_ = event->Get<double>("gen_neutral_phi_2") : -9999.;
    event->Exists("gen_neutral_eta_1") ? gen_neutral_eta_1_ = event->Get<double>("gen_neutral_eta_1") : -9999.;
    event->Exists("gen_neutral_eta_2") ? gen_neutral_eta_2_ = event->Get<double>("gen_neutral_eta_2") : -9999.;


    uncorrmet_ = met_;
    if (event->Exists("met_norecoil")) uncorrmet_ = event->Get<double>("met_norecoil");
    uncorrmet_phi_ = met_phi_;
    if (event->Exists("met_phi_norecoil")) uncorrmet_phi_ = event->Get<double>("met_phi_norecoil");

    std::vector<Met*> pfMet_vec = event->GetPtrVec<Met>("pfMetFromSlimmed");
    pf_met_ = pfMet_vec[0]->pt();

    pf_mt_1_ = MT(lep1, pfMet_vec[0]);

    metCov00_ = mets->xx_sig();
    metCov10_ = mets->yx_sig();
    metCov01_ = mets->xy_sig();
    metCov11_ = mets->yy_sig();
    
    if (channel_ == channel::et) {
      Electron const* elec = dynamic_cast<Electron const*>(lep1);
      Tau const* tau = dynamic_cast<Tau const*>(lep2);
      d0_1_ = elec->dxy_vertex();
      dz_1_ = elec->dz_vertex();
      d0_2_ = tau->lead_dxy_vertex();
      dz_2_ = tau->lead_dz_vertex();
     
      iso_1_ = PF03EAIsolationVal(elec, eventInfo->jet_rho());
      // add deepTau ID v2
      // Raw DNN scores
      iso_2_      = tau->HasTauID("byDeepTau2017v2p1VSjetraw")      ? tau->GetTauID("byDeepTau2017v2p1VSjetraw"):      0.;
      deepTauVsEle_iso_2_       = tau->HasTauID("byDeepTau2017v2p1VSeraw")        ? tau->GetTauID("byDeepTau2017v2p1VSeraw"):        0.;
      deepTauVsMu_iso_2_        = tau->HasTauID("byDeepTau2017v2p1VSmuraw")       ? tau->GetTauID("byDeepTau2017v2p1VSmuraw"):        0.;

      // Existing workpoints
      deepTauVsJets_vvvloose_2_ = tau->HasTauID("byVVVLooseDeepTau2017v2p1VSjet") ? tau->GetTauID("byVVVLooseDeepTau2017v2p1VSjet"): 0.;
      deepTauVsJets_vvloose_2_  = tau->HasTauID("byVVLooseDeepTau2017v2p1VSjet")  ? tau->GetTauID("byVVLooseDeepTau2017v2p1VSjet"):  0.;
      deepTauVsJets_vloose_2_   = tau->HasTauID("byVLooseDeepTau2017v2p1VSjet")   ? tau->GetTauID("byVLooseDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_loose_2_    = tau->HasTauID("byLooseDeepTau2017v2p1VSjet")    ? tau->GetTauID("byLooseDeepTau2017v2p1VSjet"):    0.;
      deepTauVsJets_medium_2_   = tau->HasTauID("byMediumDeepTau2017v2p1VSjet")   ? tau->GetTauID("byMediumDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_tight_2_    = tau->HasTauID("byTightDeepTau2017v2p1VSjet")    ? tau->GetTauID("byTightDeepTau2017v2p1VSjet"):    0.;
      deepTauVsJets_vtight_2_   = tau->HasTauID("byVTightDeepTau2017v2p1VSjet")   ? tau->GetTauID("byVTightDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_vvtight_2_  = tau->HasTauID("byVVTightDeepTau2017v2p1VSjet")  ? tau->GetTauID("byVVTightDeepTau2017v2p1VSjet"):  0.;

      deepTauVsEle_vvvloose_2_  = tau->HasTauID("byVVVLooseDeepTau2017v2p1VSe")   ? tau->GetTauID("byVVVLooseDeepTau2017v2p1VSe"):   0.;
      deepTauVsEle_vvloose_2_   = tau->HasTauID("byVVLooseDeepTau2017v2p1VSe")    ? tau->GetTauID("byVVLooseDeepTau2017v2p1VSe"):    0.;
      deepTauVsEle_vloose_2_    = tau->HasTauID("byVLooseDeepTau2017v2p1VSe")     ? tau->GetTauID("byVLooseDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_loose_2_     = tau->HasTauID("byLooseDeepTau2017v2p1VSe")      ? tau->GetTauID("byLooseDeepTau2017v2p1VSe"):      0.;
      deepTauVsEle_medium_2_    = tau->HasTauID("byMediumDeepTau2017v2p1VSe")     ? tau->GetTauID("byMediumDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_tight_2_     = tau->HasTauID("byTightDeepTau2017v2p1VSe")      ? tau->GetTauID("byTightDeepTau2017v2p1VSe"):      0.;
      deepTauVsEle_vtight_2_    = tau->HasTauID("byVTightDeepTau2017v2p1VSe")     ? tau->GetTauID("byVTightDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_vvtight_2_   = tau->HasTauID("byVVTightDeepTau2017v2p1VSe")    ? tau->GetTauID("byVVTightDeepTau2017v2p1VSe"):    0.;

      deepTauVsMu_vvvloose_2_   = tau->HasTauID("byVVVLooseDeepTau2017v2p1VSmu")  ? tau->GetTauID("byVVVLooseDeepTau2017v2p1VSmu"):   0.;
      deepTauVsMu_vvloose_2_    = tau->HasTauID("byVVLooseDeepTau2017v2p1VSmu")   ? tau->GetTauID("byVVLooseDeepTau2017v2p1VSmu"):    0.;
      deepTauVsMu_vloose_2_     = tau->HasTauID("byVLooseDeepTau2017v2p1VSmu")    ? tau->GetTauID("byVLooseDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_loose_2_      = tau->HasTauID("byLooseDeepTau2017v2p1VSmu")     ? tau->GetTauID("byLooseDeepTau2017v2p1VSmu"):      0.;
      deepTauVsMu_medium_2_     = tau->HasTauID("byMediumDeepTau2017v2p1VSmu")    ? tau->GetTauID("byMediumDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_tight_2_      = tau->HasTauID("byTightDeepTau2017v2p1VSmu")     ? tau->GetTauID("byTightDeepTau2017v2p1VSmu"):      0.;
      deepTauVsMu_vtight_2_     = tau->HasTauID("byVTightDeepTau2017v2p1VSmu")    ? tau->GetTauID("byVTightDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_vvtight_2_    = tau->HasTauID("byVVTightDeepTau2017v2p1VSmu")   ? tau->GetTauID("byVVTightDeepTau2017v2p1VSmu"):    0.;

      // DeepTau v2p5
      deepTauVsJetsV2p5_vvvloose_2_ = tau->HasTauID("byVVVLooseDeepTau2018v2p5VSjet") ? tau->GetTauID("byVVVLooseDeepTau2018v2p5VSjet"): 0.;
      deepTauVsJetsV2p5_vvloose_2_  = tau->HasTauID("byVVLooseDeepTau2018v2p5VSjet")  ? tau->GetTauID("byVVLooseDeepTau2018v2p5VSjet"):  0.;
      deepTauVsJetsV2p5_vloose_2_   = tau->HasTauID("byVLooseDeepTau2018v2p5VSjet")   ? tau->GetTauID("byVLooseDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_loose_2_    = tau->HasTauID("byLooseDeepTau2018v2p5VSjet")    ? tau->GetTauID("byLooseDeepTau2018v2p5VSjet"):    0.;
      deepTauVsJetsV2p5_medium_2_   = tau->HasTauID("byMediumDeepTau2018v2p5VSjet")   ? tau->GetTauID("byMediumDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_tight_2_    = tau->HasTauID("byTightDeepTau2018v2p5VSjet")    ? tau->GetTauID("byTightDeepTau2018v2p5VSjet"):    0.;
      deepTauVsJetsV2p5_vtight_2_   = tau->HasTauID("byVTightDeepTau2018v2p5VSjet")   ? tau->GetTauID("byVTightDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_vvtight_2_  = tau->HasTauID("byVVTightDeepTau2018v2p5VSjet")  ? tau->GetTauID("byVVTightDeepTau2018v2p5VSjet"):  0.;
      deepTauVsMuV2p5_vloose_2_     = tau->HasTauID("byVLooseDeepTau2018v2p5VSmu")    ? tau->GetTauID("byVLooseDeepTau2018v2p5VSmu"):     0.;
      deepTauVsMuV2p5_tight_2_      = tau->HasTauID("byTightDeepTau2018v2p5VSmu")     ? tau->GetTauID("byTightDeepTau2018v2p5VSmu"):      0.;
      deepTauVsEleV2p5_vvloose_2_   = tau->HasTauID("byVVLooseDeepTau2018v2p5VSe")    ? tau->GetTauID("byVVLooseDeepTau2018v2p5VSe"):    0.;
      deepTauVsEleV2p5_tight_2_     = tau->HasTauID("byTightDeepTau2018v2p5VSe")      ? tau->GetTauID("byTightDeepTau2018v2p5VSe"):      0.;
      iso_2_V2p5_      = tau->HasTauID("byDeepTau2018v2p5VSjetraw")      ? tau->GetTauID("byDeepTau2018v2p5VSjetraw"):      0.;

    }
    if (channel_ == channel::mt || channel_ == channel::mtmet) {
      Muon const* muon = dynamic_cast<Muon const*>(lep1);
      Tau const* tau = dynamic_cast<Tau const*>(lep2);
      d0_1_ = muon->dxy_vertex();
      dz_1_ = muon->dz_vertex();
      d0_2_ = tau->lead_dxy_vertex();
      dz_2_ = tau->lead_dz_vertex();
      iso_1_ = PF04IsolationVal(muon, 0.5, 0);
      // add deepTau ID v2
      // Raw DNN scores
      iso_2_      = tau->HasTauID("byDeepTau2017v2p1VSjetraw")      ? tau->GetTauID("byDeepTau2017v2p1VSjetraw"):      0.;
      deepTauVsEle_iso_2_       = tau->HasTauID("byDeepTau2017v2p1VSeraw")        ? tau->GetTauID("byDeepTau2017v2p1VSeraw"):        0.;
      deepTauVsMu_iso_2_        = tau->HasTauID("byDeepTau2017v2p1VSmuraw")       ? tau->GetTauID("byDeepTau2017v2p1VSmuraw"):        0.;

      // Existing workpoints
      deepTauVsJets_vvvloose_2_ = tau->HasTauID("byVVVLooseDeepTau2017v2p1VSjet") ? tau->GetTauID("byVVVLooseDeepTau2017v2p1VSjet"): 0.;
      deepTauVsJets_vvloose_2_  = tau->HasTauID("byVVLooseDeepTau2017v2p1VSjet")  ? tau->GetTauID("byVVLooseDeepTau2017v2p1VSjet"):  0.;
      deepTauVsJets_vloose_2_   = tau->HasTauID("byVLooseDeepTau2017v2p1VSjet")   ? tau->GetTauID("byVLooseDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_loose_2_    = tau->HasTauID("byLooseDeepTau2017v2p1VSjet")    ? tau->GetTauID("byLooseDeepTau2017v2p1VSjet"):    0.;
      deepTauVsJets_medium_2_   = tau->HasTauID("byMediumDeepTau2017v2p1VSjet")   ? tau->GetTauID("byMediumDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_tight_2_    = tau->HasTauID("byTightDeepTau2017v2p1VSjet")    ? tau->GetTauID("byTightDeepTau2017v2p1VSjet"):    0.;
      deepTauVsJets_vtight_2_   = tau->HasTauID("byVTightDeepTau2017v2p1VSjet")   ? tau->GetTauID("byVTightDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_vvtight_2_  = tau->HasTauID("byVVTightDeepTau2017v2p1VSjet")  ? tau->GetTauID("byVVTightDeepTau2017v2p1VSjet"):  0.;

      deepTauVsEle_vvvloose_2_  = tau->HasTauID("byVVVLooseDeepTau2017v2p1VSe")   ? tau->GetTauID("byVVVLooseDeepTau2017v2p1VSe"):   0.;
      deepTauVsEle_vvloose_2_   = tau->HasTauID("byVVLooseDeepTau2017v2p1VSe")    ? tau->GetTauID("byVVLooseDeepTau2017v2p1VSe"):    0.;
      deepTauVsEle_vloose_2_    = tau->HasTauID("byVLooseDeepTau2017v2p1VSe")     ? tau->GetTauID("byVLooseDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_loose_2_     = tau->HasTauID("byLooseDeepTau2017v2p1VSe")      ? tau->GetTauID("byLooseDeepTau2017v2p1VSe"):      0.;
      deepTauVsEle_medium_2_    = tau->HasTauID("byMediumDeepTau2017v2p1VSe")     ? tau->GetTauID("byMediumDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_tight_2_     = tau->HasTauID("byTightDeepTau2017v2p1VSe")      ? tau->GetTauID("byTightDeepTau2017v2p1VSe"):      0.;
      deepTauVsEle_vtight_2_    = tau->HasTauID("byVTightDeepTau2017v2p1VSe")     ? tau->GetTauID("byVTightDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_vvtight_2_   = tau->HasTauID("byVVTightDeepTau2017v2p1VSe")    ? tau->GetTauID("byVVTightDeepTau2017v2p1VSe"):    0.;

      deepTauVsMu_vvvloose_2_   = tau->HasTauID("byVVVLooseDeepTau2017v2p1VSmu")  ? tau->GetTauID("byVVVLooseDeepTau2017v2p1VSmu"):   0.;
      deepTauVsMu_vvloose_2_    = tau->HasTauID("byVVLooseDeepTau2017v2p1VSmu")   ? tau->GetTauID("byVVLooseDeepTau2017v2p1VSmu"):    0.;
      deepTauVsMu_vloose_2_     = tau->HasTauID("byVLooseDeepTau2017v2p1VSmu")    ? tau->GetTauID("byVLooseDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_loose_2_      = tau->HasTauID("byLooseDeepTau2017v2p1VSmu")     ? tau->GetTauID("byLooseDeepTau2017v2p1VSmu"):      0.;
      deepTauVsMu_medium_2_     = tau->HasTauID("byMediumDeepTau2017v2p1VSmu")    ? tau->GetTauID("byMediumDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_tight_2_      = tau->HasTauID("byTightDeepTau2017v2p1VSmu")     ? tau->GetTauID("byTightDeepTau2017v2p1VSmu"):      0.;
      deepTauVsMu_vtight_2_     = tau->HasTauID("byVTightDeepTau2017v2p1VSmu")    ? tau->GetTauID("byVTightDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_vvtight_2_    = tau->HasTauID("byVVTightDeepTau2017v2p1VSmu")   ? tau->GetTauID("byVVTightDeepTau2017v2p1VSmu"):    0.;

      // DeepTau v2p5
      deepTauVsJetsV2p5_vvvloose_2_ = tau->HasTauID("byVVVLooseDeepTau2018v2p5VSjet") ? tau->GetTauID("byVVVLooseDeepTau2018v2p5VSjet"): 0.;
      deepTauVsJetsV2p5_vvloose_2_  = tau->HasTauID("byVVLooseDeepTau2018v2p5VSjet")  ? tau->GetTauID("byVVLooseDeepTau2018v2p5VSjet"):  0.;
      deepTauVsJetsV2p5_vloose_2_   = tau->HasTauID("byVLooseDeepTau2018v2p5VSjet")   ? tau->GetTauID("byVLooseDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_loose_2_    = tau->HasTauID("byLooseDeepTau2018v2p5VSjet")    ? tau->GetTauID("byLooseDeepTau2018v2p5VSjet"):    0.;
      deepTauVsJetsV2p5_medium_2_   = tau->HasTauID("byMediumDeepTau2018v2p5VSjet")   ? tau->GetTauID("byMediumDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_tight_2_    = tau->HasTauID("byTightDeepTau2018v2p5VSjet")    ? tau->GetTauID("byTightDeepTau2018v2p5VSjet"):    0.;
      deepTauVsJetsV2p5_vtight_2_   = tau->HasTauID("byVTightDeepTau2018v2p5VSjet")   ? tau->GetTauID("byVTightDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_vvtight_2_  = tau->HasTauID("byVVTightDeepTau2018v2p5VSjet")  ? tau->GetTauID("byVVTightDeepTau2018v2p5VSjet"):  0.;
      deepTauVsMuV2p5_vloose_2_     = tau->HasTauID("byVLooseDeepTau2018v2p5VSmu")    ? tau->GetTauID("byVLooseDeepTau2018v2p5VSmu"):     0.;
      deepTauVsMuV2p5_tight_2_      = tau->HasTauID("byTightDeepTau2018v2p5VSmu")     ? tau->GetTauID("byTightDeepTau2018v2p5VSmu"):      0.;
      deepTauVsEleV2p5_vvloose_2_   = tau->HasTauID("byVVLooseDeepTau2018v2p5VSe")    ? tau->GetTauID("byVVLooseDeepTau2018v2p5VSe"):    0.;
      deepTauVsEleV2p5_tight_2_     = tau->HasTauID("byTightDeepTau2018v2p5VSe")      ? tau->GetTauID("byTightDeepTau2018v2p5VSe"):      0.;
      iso_2_V2p5_      = tau->HasTauID("byDeepTau2018v2p5VSjetraw")      ? tau->GetTauID("byDeepTau2018v2p5VSjetraw"):      0.;

    }
    if (channel_ == channel::em) {
      Electron const* elec = dynamic_cast<Electron const*>(lep1);
      Muon const* muon = dynamic_cast<Muon const*>(lep2);
      iso_1_ = PF03EAIsolationVal(elec, eventInfo->jet_rho()); //lepton_rho
      iso_2_ = PF04IsolationVal(muon, 0.5, 0);
      d0_1_ = elec->dxy_vertex();
      dz_1_ = elec->dz_vertex();
      d0_2_ = muon->dxy_vertex();
      dz_2_ = muon->dz_vertex();
      raw_pt_1_ = elec->ecalTrkEnergyPreCorr()/elec->energy()*pt_1_.var_float;
    }
    if (channel_ == channel::tt) {
      Tau const* tau1 = dynamic_cast<Tau const*>(lep1);
      Tau const* tau2 = dynamic_cast<Tau const*>(lep2);
      d0_1_ = tau1->lead_dxy_vertex();
      dz_1_ = tau1->lead_dz_vertex();
      d0_2_ = tau2->lead_dxy_vertex();
      dz_2_ = tau2->lead_dz_vertex();

      // add deepTau ID v2
      // Raw DNN scores
      iso_1_      = tau1->HasTauID("byDeepTau2017v2p1VSjetraw")      ? tau1->GetTauID("byDeepTau2017v2p1VSjetraw"):      0.;
      iso_2_      = tau2->HasTauID("byDeepTau2017v2p1VSjetraw")      ? tau2->GetTauID("byDeepTau2017v2p1VSjetraw"):      0.;
      deepTauVsEle_iso_1_       = tau1->HasTauID("byDeepTau2017v2p1VSeraw")        ? tau1->GetTauID("byDeepTau2017v2p1VSeraw"):        0.;
      deepTauVsEle_iso_2_       = tau2->HasTauID("byDeepTau2017v2p1VSeraw")        ? tau2->GetTauID("byDeepTau2017v2p1VSeraw"):        0.;
      deepTauVsMu_iso_1_       = tau1->HasTauID("byDeepTau2017v2p1VSmuraw")        ? tau1->GetTauID("byDeepTau2017v2p1VSmuraw"):        0.;
      deepTauVsMu_iso_2_       = tau2->HasTauID("byDeepTau2017v2p1VSmuraw")        ? tau2->GetTauID("byDeepTau2017v2p1VSmuraw"):        0.;

      // Existing workpoints
      deepTauVsJets_vvvloose_1_ = tau1->HasTauID("byVVVLooseDeepTau2017v2p1VSjet") ? tau1->GetTauID("byVVVLooseDeepTau2017v2p1VSjet"): 0.;
      deepTauVsJets_vvloose_1_  = tau1->HasTauID("byVVLooseDeepTau2017v2p1VSjet")  ? tau1->GetTauID("byVVLooseDeepTau2017v2p1VSjet"):  0.;
      deepTauVsJets_vloose_1_   = tau1->HasTauID("byVLooseDeepTau2017v2p1VSjet")   ? tau1->GetTauID("byVLooseDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_loose_1_    = tau1->HasTauID("byLooseDeepTau2017v2p1VSjet")    ? tau1->GetTauID("byLooseDeepTau2017v2p1VSjet"):    0.;
      deepTauVsJets_medium_1_   = tau1->HasTauID("byMediumDeepTau2017v2p1VSjet")   ? tau1->GetTauID("byMediumDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_tight_1_    = tau1->HasTauID("byTightDeepTau2017v2p1VSjet")    ? tau1->GetTauID("byTightDeepTau2017v2p1VSjet"):    0.;
      deepTauVsJets_vtight_1_   = tau1->HasTauID("byVTightDeepTau2017v2p1VSjet")   ? tau1->GetTauID("byVTightDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_vvtight_1_  = tau1->HasTauID("byVVTightDeepTau2017v2p1VSjet")  ? tau1->GetTauID("byVVTightDeepTau2017v2p1VSjet"):  0.;

      deepTauVsJets_vvvloose_2_ = tau2->HasTauID("byVVVLooseDeepTau2017v2p1VSjet") ? tau2->GetTauID("byVVVLooseDeepTau2017v2p1VSjet"): 0.;
      deepTauVsJets_vvloose_2_  = tau2->HasTauID("byVVLooseDeepTau2017v2p1VSjet")  ? tau2->GetTauID("byVVLooseDeepTau2017v2p1VSjet"):  0.;
      deepTauVsJets_vloose_2_   = tau2->HasTauID("byVLooseDeepTau2017v2p1VSjet")   ? tau2->GetTauID("byVLooseDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_loose_2_    = tau2->HasTauID("byLooseDeepTau2017v2p1VSjet")    ? tau2->GetTauID("byLooseDeepTau2017v2p1VSjet"):    0.;
      deepTauVsJets_medium_2_   = tau2->HasTauID("byMediumDeepTau2017v2p1VSjet")   ? tau2->GetTauID("byMediumDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_tight_2_    = tau2->HasTauID("byTightDeepTau2017v2p1VSjet")    ? tau2->GetTauID("byTightDeepTau2017v2p1VSjet"):    0.;
      deepTauVsJets_vtight_2_   = tau2->HasTauID("byVTightDeepTau2017v2p1VSjet")   ? tau2->GetTauID("byVTightDeepTau2017v2p1VSjet"):   0.;
      deepTauVsJets_vvtight_2_  = tau2->HasTauID("byVVTightDeepTau2017v2p1VSjet")  ? tau2->GetTauID("byVVTightDeepTau2017v2p1VSjet"):  0.;

      deepTauVsEle_vvvloose_1_  = tau1->HasTauID("byVVVLooseDeepTau2017v2p1VSe")   ? tau1->GetTauID("byVVVLooseDeepTau2017v2p1VSe"):   0.;
      deepTauVsEle_vvloose_1_   = tau1->HasTauID("byVVLooseDeepTau2017v2p1VSe")    ? tau1->GetTauID("byVVLooseDeepTau2017v2p1VSe"):    0.;
      deepTauVsEle_vloose_1_    = tau1->HasTauID("byVLooseDeepTau2017v2p1VSe")     ? tau1->GetTauID("byVLooseDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_loose_1_     = tau1->HasTauID("byLooseDeepTau2017v2p1VSe")      ? tau1->GetTauID("byLooseDeepTau2017v2p1VSe"):      0.;
      deepTauVsEle_medium_1_    = tau1->HasTauID("byMediumDeepTau2017v2p1VSe")     ? tau1->GetTauID("byMediumDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_tight_1_     = tau1->HasTauID("byTightDeepTau2017v2p1VSe")      ? tau1->GetTauID("byTightDeepTau2017v2p1VSe"):      0.;
      deepTauVsEle_vtight_1_    = tau1->HasTauID("byVTightDeepTau2017v2p1VSe")     ? tau1->GetTauID("byVTightDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_vvtight_1_   = tau1->HasTauID("byVVTightDeepTau2017v2p1VSe")    ? tau1->GetTauID("byVVTightDeepTau2017v2p1VSe"):    0.;

      deepTauVsEle_vvvloose_2_  = tau2->HasTauID("byVVVLooseDeepTau2017v2p1VSe")   ? tau2->GetTauID("byVVVLooseDeepTau2017v2p1VSe"):   0.;
      deepTauVsEle_vvloose_2_   = tau2->HasTauID("byVVLooseDeepTau2017v2p1VSe")    ? tau2->GetTauID("byVVLooseDeepTau2017v2p1VSe"):    0.;
      deepTauVsEle_vloose_2_    = tau2->HasTauID("byVLooseDeepTau2017v2p1VSe")     ? tau2->GetTauID("byVLooseDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_loose_2_     = tau2->HasTauID("byLooseDeepTau2017v2p1VSe")      ? tau2->GetTauID("byLooseDeepTau2017v2p1VSe"):      0.;
      deepTauVsEle_medium_2_    = tau2->HasTauID("byMediumDeepTau2017v2p1VSe")     ? tau2->GetTauID("byMediumDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_tight_2_     = tau2->HasTauID("byTightDeepTau2017v2p1VSe")      ? tau2->GetTauID("byTightDeepTau2017v2p1VSe"):      0.;
      deepTauVsEle_vtight_2_    = tau2->HasTauID("byVTightDeepTau2017v2p1VSe")     ? tau2->GetTauID("byVTightDeepTau2017v2p1VSe"):     0.;
      deepTauVsEle_vvtight_2_   = tau2->HasTauID("byVVTightDeepTau2017v2p1VSe")    ? tau2->GetTauID("byVVTightDeepTau2017v2p1VSe"):    0.;

      deepTauVsMu_vvvloose_1_  = tau1->HasTauID("byVVVLooseDeepTau2017v2p1VSmu")   ? tau1->GetTauID("byVVVLooseDeepTau2017v2p1VSmu"):   0.;
      deepTauVsMu_vvloose_1_   = tau1->HasTauID("byVVLooseDeepTau2017v2p1VSmu")    ? tau1->GetTauID("byVVLooseDeepTau2017v2p1VSmu"):    0.;
      deepTauVsMu_vloose_1_    = tau1->HasTauID("byVLooseDeepTau2017v2p1VSmu")     ? tau1->GetTauID("byVLooseDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_loose_1_     = tau1->HasTauID("byLooseDeepTau2017v2p1VSmu")      ? tau1->GetTauID("byLooseDeepTau2017v2p1VSmu"):      0.;
      deepTauVsMu_medium_1_    = tau1->HasTauID("byMediumDeepTau2017v2p1VSmu")     ? tau1->GetTauID("byMediumDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_tight_1_     = tau1->HasTauID("byTightDeepTau2017v2p1VSmu")      ? tau1->GetTauID("byTightDeepTau2017v2p1VSmu"):      0.;
      deepTauVsMu_vtight_1_    = tau1->HasTauID("byVTightDeepTau2017v2p1VSmu")     ? tau1->GetTauID("byVTightDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_vvtight_1_   = tau1->HasTauID("byVVTightDeepTau2017v2p1VSmu")    ? tau1->GetTauID("byVVTightDeepTau2017v2p1VSmu"):    0.;

      deepTauVsMu_vvvloose_2_  = tau2->HasTauID("byVVVLooseDeepTau2017v2p1VSmu")   ? tau2->GetTauID("byVVVLooseDeepTau2017v2p1VSmu"):   0.;
      deepTauVsMu_vvloose_2_   = tau2->HasTauID("byVVLooseDeepTau2017v2p1VSmu")    ? tau2->GetTauID("byVVLooseDeepTau2017v2p1VSmu"):    0.;
      deepTauVsMu_vloose_2_    = tau2->HasTauID("byVLooseDeepTau2017v2p1VSmu")     ? tau2->GetTauID("byVLooseDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_loose_2_     = tau2->HasTauID("byLooseDeepTau2017v2p1VSmu")      ? tau2->GetTauID("byLooseDeepTau2017v2p1VSmu"):      0.;
      deepTauVsMu_medium_2_    = tau2->HasTauID("byMediumDeepTau2017v2p1VSmu")     ? tau2->GetTauID("byMediumDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_tight_2_     = tau2->HasTauID("byTightDeepTau2017v2p1VSmu")      ? tau2->GetTauID("byTightDeepTau2017v2p1VSmu"):      0.;
      deepTauVsMu_vtight_2_    = tau2->HasTauID("byVTightDeepTau2017v2p1VSmu")     ? tau2->GetTauID("byVTightDeepTau2017v2p1VSmu"):     0.;
      deepTauVsMu_vvtight_2_   = tau2->HasTauID("byVVTightDeepTau2017v2p1VSmu")    ? tau2->GetTauID("byVVTightDeepTau2017v2p1VSmu"):    0.;

      // DeepTau v2p5
      deepTauVsJetsV2p5_vvvloose_1_ = tau1->HasTauID("byVVVLooseDeepTau2018v2p5VSjet") ? tau1->GetTauID("byVVVLooseDeepTau2018v2p5VSjet"): 0.;
      deepTauVsJetsV2p5_vvloose_1_  = tau1->HasTauID("byVVLooseDeepTau2018v2p5VSjet")  ? tau1->GetTauID("byVVLooseDeepTau2018v2p5VSjet"):  0.;
      deepTauVsJetsV2p5_vloose_1_   = tau1->HasTauID("byVLooseDeepTau2018v2p5VSjet")   ? tau1->GetTauID("byVLooseDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_loose_1_    = tau1->HasTauID("byLooseDeepTau2018v2p5VSjet")    ? tau1->GetTauID("byLooseDeepTau2018v2p5VSjet"):    0.;
      deepTauVsJetsV2p5_medium_1_   = tau1->HasTauID("byMediumDeepTau2018v2p5VSjet")   ? tau1->GetTauID("byMediumDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_tight_1_    = tau1->HasTauID("byTightDeepTau2018v2p5VSjet")    ? tau1->GetTauID("byTightDeepTau2018v2p5VSjet"):    0.;
      deepTauVsJetsV2p5_vtight_1_   = tau1->HasTauID("byVTightDeepTau2018v2p5VSjet")   ? tau1->GetTauID("byVTightDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_vvtight_1_  = tau1->HasTauID("byVVTightDeepTau2018v2p5VSjet")  ? tau1->GetTauID("byVVTightDeepTau2018v2p5VSjet"):  0.;
      deepTauVsMuV2p5_vloose_1_     = tau1->HasTauID("byVLooseDeepTau2018v2p5VSmu")    ? tau1->GetTauID("byVLooseDeepTau2018v2p5VSmu"):     0.;
      deepTauVsMuV2p5_tight_1_      = tau1->HasTauID("byTightDeepTau2018v2p5VSmu")     ? tau1->GetTauID("byTightDeepTau2018v2p5VSmu"):      0.;
      deepTauVsEleV2p5_vvloose_1_   = tau1->HasTauID("byVVLooseDeepTau2018v2p5VSe")    ? tau1->GetTauID("byVVLooseDeepTau2018v2p5VSe"):    0.;
      deepTauVsEleV2p5_tight_1_     = tau1->HasTauID("byTightDeepTau2018v2p5VSe")      ? tau1->GetTauID("byTightDeepTau2018v2p5VSe"):      0.;
      iso_1_V2p5_      = tau1->HasTauID("byDeepTau2018v2p5VSjetraw")      ? tau1->GetTauID("byDeepTau2018v2p5VSjetraw"):      0.;

      deepTauVsJetsV2p5_vvvloose_2_ = tau2->HasTauID("byVVVLooseDeepTau2018v2p5VSjet") ? tau2->GetTauID("byVVVLooseDeepTau2018v2p5VSjet"): 0.;
      deepTauVsJetsV2p5_vvloose_2_  = tau2->HasTauID("byVVLooseDeepTau2018v2p5VSjet")  ? tau2->GetTauID("byVVLooseDeepTau2018v2p5VSjet"):  0.;
      deepTauVsJetsV2p5_vloose_2_   = tau2->HasTauID("byVLooseDeepTau2018v2p5VSjet")   ? tau2->GetTauID("byVLooseDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_loose_2_    = tau2->HasTauID("byLooseDeepTau2018v2p5VSjet")    ? tau2->GetTauID("byLooseDeepTau2018v2p5VSjet"):    0.;
      deepTauVsJetsV2p5_medium_2_   = tau2->HasTauID("byMediumDeepTau2018v2p5VSjet")   ? tau2->GetTauID("byMediumDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_tight_2_    = tau2->HasTauID("byTightDeepTau2018v2p5VSjet")    ? tau2->GetTauID("byTightDeepTau2018v2p5VSjet"):    0.;
      deepTauVsJetsV2p5_vtight_2_   = tau2->HasTauID("byVTightDeepTau2018v2p5VSjet")   ? tau2->GetTauID("byVTightDeepTau2018v2p5VSjet"):   0.;
      deepTauVsJetsV2p5_vvtight_2_  = tau2->HasTauID("byVVTightDeepTau2018v2p5VSjet")  ? tau2->GetTauID("byVVTightDeepTau2018v2p5VSjet"):  0.;
      deepTauVsMuV2p5_vloose_2_     = tau2->HasTauID("byVLooseDeepTau2018v2p5VSmu")    ? tau2->GetTauID("byVLooseDeepTau2018v2p5VSmu"):     0.;
      deepTauVsMuV2p5_tight_2_      = tau2->HasTauID("byTightDeepTau2018v2p5VSmu")     ? tau2->GetTauID("byTightDeepTau2018v2p5VSmu"):      0.;
      deepTauVsEleV2p5_vvloose_2_   = tau2->HasTauID("byVVLooseDeepTau2018v2p5VSe")    ? tau2->GetTauID("byVVLooseDeepTau2018v2p5VSe"):    0.;
      deepTauVsEleV2p5_tight_2_     = tau2->HasTauID("byTightDeepTau2018v2p5VSe")      ? tau2->GetTauID("byTightDeepTau2018v2p5VSe"):      0.;
      iso_2_V2p5_      = tau2->HasTauID("byDeepTau2018v2p5VSjetraw")      ? tau2->GetTauID("byDeepTau2018v2p5VSjetraw"):      0.;
    }
    if (channel_ == channel::zee || channel_ == channel::tpzee) {
      Electron const* elec1 = dynamic_cast<Electron const*>(lep1);
      Electron const* elec2 = dynamic_cast<Electron const*>(lep2);
      iso_1_ = PF03EAIsolationVal(elec1, eventInfo->jet_rho()); //should be lepton_rho
      iso_2_ = PF03EAIsolationVal(elec2, eventInfo->jet_rho()); 
      d0_1_ = elec1->dxy_vertex();
      dz_1_ = elec1->dz_vertex();
      d0_2_ = elec2->dxy_vertex();
      dz_2_ = elec2->dz_vertex();
    }
    if (channel_ == channel::zmm || channel_ == channel::tpzmm) {
      Muon const* muon1 = dynamic_cast<Muon const*>(lep1);
      Muon const* muon2 = dynamic_cast<Muon const*>(lep2);

      iso_1_ = PF04IsolationVal(muon1, 0.5, 0);
      iso_2_ = PF04IsolationVal(muon2, 0.5, 0);
      d0_1_ = muon1->dxy_vertex();
      dz_1_ = muon1->dz_vertex();
      d0_2_ = muon2->dxy_vertex();
      dz_2_ = muon2->dz_vertex();
    }

    Tau const* tau1 = dynamic_cast<Tau const*>(lep1);
    Tau const* tau2 = dynamic_cast<Tau const*>(lep2);
    if (tau1) {
      tau_decay_mode_1_ = tau1->decay_mode();
    } else {
      tau_decay_mode_1_ = 0;
    }
    if (tau2) {
      tau_decay_mode_2_ = tau2->decay_mode();
    } else {
      tau_decay_mode_2_ = 0;
    }

    n_jets_ = jets.size();
    n_lowpt_jets_ = lowpt_jets.size();
    n_bjets_ = bjets.size();
    n_loose_bjets_ = loose_bjets.size();
    n_btag_ = btag.size();
    n_loose_btag_ = loose_btag.size();

    std::vector<Candidate *> leading_lepton;
    std::vector<Candidate *> subleading_lepton;
    leading_lepton.push_back(ditau->GetCandidate("lepton1"));
    subleading_lepton.push_back(ditau->GetCandidate("lepton2")); 
    std::vector<std::pair<ic::PFJet *, ic::Candidate *>> mu_matches = MatchByDR(uncleaned_jets, leading_lepton, 0.5, true, true);
    std::vector<std::pair<ic::PFJet *, ic::Candidate *>> tau_matches = MatchByDR(uncleaned_jets, subleading_lepton, 0.5, true, true);

    ROOT::Math::PxPyPzEVector shift_before; 
    ROOT::Math::PxPyPzEVector shift_after; 

    if(mu_matches.size() > 0) {
        jet_flav_1_ = (mu_matches.at(0)).first->parton_flavour();
        jet_pt_1_ = (mu_matches.at(0)).first->pt();

        shift_before += (mu_matches.at(0)).first->vector();
        shift_after += (mu_matches.at(0)).first->vector()*(mu_matches.at(0)).first->uncorrected_energy()/(mu_matches.at(0)).first->energy();
    } else {
      jet_flav_1_ = -9999;
      jet_pt_1_ = -9999;
    }
    if(tau_matches.size() > 0) {
        jet_flav_2_ = (tau_matches.at(0)).first->parton_flavour();
        jet_pt_2_ = (tau_matches.at(0)).first->pt();

        shift_before += (tau_matches.at(0)).first->vector();
        shift_after += (tau_matches.at(0)).first->vector()*(tau_matches.at(0)).first->uncorrected_energy()/(tau_matches.at(0)).first->energy();

    } else {
      jet_flav_2_ = -9999;
      jet_pt_2_ = -9999;
    }
    if(jets.size() > 0) {
        jet_flav_3_ = jets[0]->parton_flavour();
    } else jet_flav_3_ = -9999;
 
 
    if (n_lowpt_jets_ >= 1) {
      jpt_1_ = lowpt_jets[0]->pt();
      jnemfrac_1_=lowpt_jets[0]->neutral_em_energy_frac();
      jeta_1_ = lowpt_jets[0]->eta();
      jphi_1_ = lowpt_jets[0]->phi();
      jmva_1_ = lowpt_jets[0]->pu_id_mva_value();
      dphi_jtt_ =  ROOT::Math::VectorUtil::DeltaPhi(lowpt_jets[0]->vector(), ditau->vector());
      //shifted_dphi_jtt_ = (dphi_jtt_ < 0.) * (dphi_jtt_ + M_PI)
      //  + (dphi_jtt_ > 0.) * (dphi_jtt_ - M_PI);
      //// Make the following two respect to Z/mumu
      //residual_pt_ =  (mets->vector() + lowpt_jets[0]->vector() + ditau->vector()).pt();
      //residual_phi_ =  (mets->vector() + lowpt_jets[0]->vector() + ditau->vector()).phi();
      //dphi_residual_tt_ =  ROOT::Math::VectorUtil::DeltaPhi(
      //  (mets->vector() + lowpt_jets[0]->vector() + ditau->vector()), 
      //  ditau->vector()
      //);
    } else {
      jpt_1_ = -9999;
      jeta_1_ = -9999;
      jphi_1_ = -9999;
      jmva_1_ = -9999;
      dphi_jtt_ = -9999.;
      //shifted_dphi_jtt_ = -9999.;
      //residual_pt_ =  -9999.;
      //residual_phi_ = -9999.;
      //dphi_residual_tt_ = -9999.;
    }

    if (n_lowpt_jets_ >= 2) {
      jpt_2_ = lowpt_jets[1]->pt();
      jeta_2_ = lowpt_jets[1]->eta();
      jphi_2_ = lowpt_jets[1]->phi();
      jmva_2_ = lowpt_jets[1]->pu_id_mva_value();
      mjj_ = (lowpt_jets[0]->vector() + lowpt_jets[1]->vector()).M();
      jdeta_ = fabs(lowpt_jets[0]->eta() - lowpt_jets[1]->eta());
      jdphi_ =  ROOT::Math::VectorUtil::DeltaPhi(lowpt_jets[0]->vector(), lowpt_jets[1]->vector());
      dijetpt_ =  (lowpt_jets[0]->vector() + lowpt_jets[1]->vector()).pt();

      if(lowpt_jets[0]->eta() > lowpt_jets[1]->eta()){
        sjdphi_ =  ROOT::Math::VectorUtil::DeltaPhi(lowpt_jets[0]->vector(), lowpt_jets[1]->vector());
      }
      else{
        sjdphi_ =  ROOT::Math::VectorUtil::DeltaPhi(lowpt_jets[1]->vector(), lowpt_jets[0]->vector());
      }
      
    } else {
      jpt_2_ = -9999;
      jeta_2_ = -9999;
      jphi_2_ = -9999;
      mjj_ = -9999;
      jdeta_ = -9999;
      jdphi_ = -9999;
      dijetpt_ = -9999;
      jmva_2_ = -9999;
    }
    
    wt_qcdscale_up_   = eventInfo->weight_defined("1005") ? eventInfo->weight("1005") : 1.0;
    wt_qcdscale_down_ = eventInfo->weight_defined("1009") ? eventInfo->weight("1009") : 1.0;
    wt_ps_isr_up_   = eventInfo->weight_defined("genweight6") ? eventInfo->weight("genweight6") : 1.0;
    wt_ps_isr_down_ = eventInfo->weight_defined("genweight8") ? eventInfo->weight("genweight8") : 1.0;
    wt_ps_fsr_up_   = eventInfo->weight_defined("genweight7") ? eventInfo->weight("genweight7") : 1.0;
    wt_ps_fsr_down_ = eventInfo->weight_defined("genweight9") ? eventInfo->weight("genweight9") : 1.0;


    // full Run2 legacy
    // add medium WP (called it tight) btags
    if (n_btag_ >= 1) {
      bpt_1_ = btag[0]->pt();
      beta_1_ = btag[0]->eta();
      bphi_1_ = btag[0]->phi();
      bmva_1_ = btag[0]->pu_id_mva_value();
    
    } else {
      bpt_1_ = -9999;
      beta_1_ = -9999;
      bphi_1_ = -9999;
      bmva_1_ = -9999;
    }

    if (n_btag_ >= 2) {
      bpt_2_ = btag[1]->pt();
      beta_2_ = btag[1]->eta();
      bphi_2_ = btag[1]->phi();
      bmva_2_ = btag[1]->pu_id_mva_value();
    
    } else {
      bpt_2_ = -9999;
      beta_2_ = -9999;
      bphi_2_ = -9999;
      bmva_2_ = -9999;
    }
    // add loose tags as well
    if (n_loose_btag_ >= 1) {
      bpt_1_ = loose_btag[0]->pt();
      beta_1_ = loose_btag[0]->eta();
      bphi_1_ = loose_btag[0]->phi();
      bmva_1_ = loose_btag[0]->pu_id_mva_value();
    
    } else {
      bpt_1_ = -9999;
      beta_1_ = -9999;
      bphi_1_ = -9999;
      bmva_1_ = -9999;
    }

    if (n_loose_btag_ >= 2) {
      bpt_2_ = loose_btag[1]->pt();
      beta_2_ = loose_btag[1]->eta();
      bphi_2_ = loose_btag[1]->phi();
      bmva_2_ = loose_btag[1]->pu_id_mva_value();
    
    } else {
      bpt_2_ = -9999;
      beta_2_ = -9999;
      bphi_2_ = -9999;
      bmva_2_ = -9999;
    }

    if (prebjets.size() >= 1) {
      bcsv_1_ = prebjets[0]->GetBDiscriminator(btag_label);
    } else {
      bcsv_1_ = -9999;
    }
    if (prebjets.size() >= 2) {
      bcsv_2_ = prebjets[1]->GetBDiscriminator(btag_label);
    } else {
      bcsv_2_ = -9999;
    }

    
    if(do_faketaus_&&(channel_==channel::zmm||channel_==channel::em)){
      std::vector<Tau *> taus = event->GetPtrVec<Tau>("taus");
      std::sort(taus.begin(), taus.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));
      tau_pt_1_=-9999;
      tau_dm_1_=-9999;

      if(channel_==channel::zmm){
        std::vector<const ic::Candidate*> muons;
        muons.push_back(lep1);
        muons.push_back(lep2);

        Met *newmet = event->GetPtr<Met>("pfMET");
        double WZ_ratio = 80.4/91.2;
        TLorentzVector muon1_vec;
        TLorentzVector muon2_vec;

        muon1_vec.SetXYZM(muons[0]->vector().Px(),muons[0]->vector().Py(),muons[0]->vector().Pz(),0);
        muon2_vec.SetXYZM(muons[1]->vector().Px(),muons[1]->vector().Py(),muons[1]->vector().Pz(),0);
        TVector3 Zboost = (muon1_vec+muon2_vec).BoostVector();
        TLorentzVector Z_vec = muon1_vec+muon2_vec;
        double scale = sqrt( pow(Z_vec.Px(),2) + pow(Z_vec.Py(),2) + pow(Z_vec.Pz(),2) + pow(Z_vec.M()*WZ_ratio,2));
        TVector3 BoostBack(Z_vec.Px()/scale, Z_vec.Py()/scale, Z_vec.Pz()/scale); 
        muon1_vec.Boost(-Zboost);
        muon2_vec.Boost(-Zboost);
        muon1_vec.SetXYZM(muon1_vec.Px()*WZ_ratio,muon1_vec.Py()*WZ_ratio,muon1_vec.Pz()*WZ_ratio ,0);
        muon2_vec.SetXYZM(muon2_vec.Px()*WZ_ratio,muon2_vec.Py()*WZ_ratio,muon2_vec.Pz()*WZ_ratio ,0);
        muon1_vec.Boost(BoostBack);
        muon2_vec.Boost(BoostBack);
        ic::Candidate* muon1 = new ic::Candidate();
        ic::Candidate* muon2 = new ic::Candidate();
        muon1->set_pt(muon1_vec.Pt());
        muon1->set_eta(muon1_vec.Eta());
        muon1->set_phi(muon1_vec.Phi());
        muon1->set_energy(muon1_vec.E());
        muon2->set_pt(muon2_vec.Pt());
        muon2->set_eta(muon2_vec.Eta());
        muon2->set_phi(muon2_vec.Phi());
        muon2->set_energy(muon2_vec.E());
        ic::Candidate *not_removed_mu;
        if (eventInfo->event() % 2 == 0){
          newmet->set_vector(muon1->vector()+newmet->vector());
          not_removed_mu = muon2;
        } else {
          newmet->set_vector(muon2->vector()+newmet->vector());
          not_removed_mu = muon1;
        }
        mt_1_nomu_ = MT(not_removed_mu, newmet); 
        mu_pt_ = not_removed_mu->pt();

        ROOT::Math::PtEtaPhiEVector newvec;
        Met* modmet = new Met();
        newvec.SetPxPyPzE(not_removed_mu->vector().Px(), not_removed_mu->vector().Py(),0.,lep1->pt());
    
        newlep->set_vector(newvec);
        modmet->set_vector(newmet->vector()+newlep->vector());
    
        newmet_ = modmet->pt();
        newmet_dphi_1_=std::fabs(ROOT::Math::VectorUtil::DeltaPhi(modmet->vector(),lep1->vector()));
        newmet_dphi_2_=std::fabs(ROOT::Math::VectorUtil::DeltaPhi(modmet->vector(),lep2->vector()));

        std::vector<Candidate *> subleading_lepton;
        if(taus.size()>0) subleading_lepton.push_back(taus[0]);
        std::vector<std::pair<ic::PFJet *, ic::Candidate *>> tau_matches = MatchByDR(uncleaned_jets, subleading_lepton, 0.5, true, true);

        if(tau_matches.size() > 0) jet_pt_1_ = (tau_matches.at(0)).first->pt();
        else jet_pt_1_ = -9999;

      }

      if(taus.size()>0) { 
        deepTauVsEle_vvvloose_1_   = taus[0]->HasTauID("byVVVLooseDeepTau2017v2p1VSe")   ? taus[0]->GetTauID("byVVVLooseDeepTau2017v2p1VSe"):   0.;
        deepTauVsEle_vvloose_1_   = taus[0]->HasTauID("byVVLooseDeepTau2017v2p1VSe")   ? taus[0]->GetTauID("byVVLooseDeepTau2017v2p1VSe"):   0.;
        deepTauVsEle_tight_1_   = taus[0]->HasTauID("byTightDeepTau2017v2p1VSe")   ? taus[0]->GetTauID("byTightDeepTau2017v2p1VSe"):   0.;
        deepTauVsMu_vloose_1_      =    taus[0]->HasTauID("byVLooseDeepTau2017v2p1VSmu")   ? taus[0]->GetTauID("byVLooseDeepTau2017v2p1VSmu"):   0.;
        deepTauVsMu_tight_1_      =    taus[0]->HasTauID("byTightDeepTau2017v2p1VSmu")   ? taus[0]->GetTauID("byTightDeepTau2017v2p1VSmu"):   0.;
        deepTauVsJets_vvvloose_1_  = taus[0]->HasTauID("byVVVLooseDeepTau2017v2p1VSjet")   ? taus[0]->GetTauID("byVVVLooseDeepTau2017v2p1VSjet"):   0.;
        deepTauVsJets_medium_1_  = taus[0]->HasTauID("byMediumDeepTau2017v2p1VSjet")   ? taus[0]->GetTauID("byMediumDeepTau2017v2p1VSjet"):   0.;
        tau_mva_decay_mode_1_ = taus[0]->HasTauID(MVADM2017) ? taus[0]->GetTauID(MVADM2017) : 0.0;
      }
      else {
        deepTauVsEle_vvvloose_1_ = 0.;
        deepTauVsEle_vvloose_1_ = 0.;
        deepTauVsEle_tight_1_ = 0.;
        deepTauVsMu_vloose_1_ = 0.;
        deepTauVsMu_tight_1_ = 0.;
        deepTauVsJets_vvvloose_1_ = 0.;
        deepTauVsJets_medium_1_ = 0.;
        tau_mva_decay_mode_1_ = -1;
      } 

      if(taus.size()>0){
        tau_pt_1_ = taus[0]->pt();
        tau_dm_1_ = taus[0]->decay_mode();
      }

    }

    // cp stuff
    wt_cp_sm_=0.; wt_cp_ps_=0.; wt_cp_mm_=0.; 
    if(event->ExistsInTree("tauspinner")){
      EventInfo const* tauspinner = event->GetPtr<EventInfo>("tauspinner");
      wt_cp_sm_ = tauspinner->weight("wt_cp_0");
      wt_cp_ps_ = tauspinner->weight("wt_cp_0p5");
      wt_cp_mm_ = tauspinner->weight("wt_cp_0p25");

      if(tauspinner->weight_defined("wt_cp_0_alt")){
        wt_cp_sm_alt_ = tauspinner->weight("wt_cp_0_alt");
        wt_cp_ps_alt_ = tauspinner->weight("wt_cp_0p5_alt");
        wt_cp_mm_alt_ = tauspinner->weight("wt_cp_0p25_alt");
      }

      if(make_mva_ntuple_){
  
        gen_phitt_ = 0.;
        double max=0.;
        double step=0.01;
        for (double i=-90+step; i<=90.; i+=step) {
          double m = CPWeight(i, wt_cp_sm_, wt_cp_ps_, wt_cp_mm_ );
          if (m > max) { max = m; gen_phitt_ = i; }
        }
        if(wt_cp_sm_==1&&wt_cp_ps_==1&&wt_cp_mm_==1) gen_phitt_ = rand_*180. - 90. ;

      }

    }
    wt_cp_prod_sm_=0.; wt_cp_prod_ps_=0.; wt_cp_prod_mm_=0.; 
    if(eventInfo->weight_defined("sm_weight_nlo")) {        
      if(eventInfo->weight_defined("sm_weight_nlo")) wt_cp_prod_sm_ = eventInfo->weight("sm_weight_nlo");
      if(eventInfo->weight_defined("ps_weight_nlo")) wt_cp_prod_ps_ = eventInfo->weight("ps_weight_nlo");
      if(eventInfo->weight_defined("mm_weight_nlo")) wt_cp_prod_mm_ = eventInfo->weight("mm_weight_nlo");
    }

    // initalise here to add to sync
    pvx_   = -9999.;
    pvy_   = -9999.;
    pvz_   = -9999.;
    ipx_1_ = -9999.;
    ipy_1_ = -9999.;
    ipz_1_ = -9999.;
    ipx_2_ = -9999.;
    ipy_2_ = -9999.;
    ipz_2_ = -9999.;

    gen_ipx_1_ = -9999.;
    gen_ipy_1_ = -9999.;
    gen_ipz_1_ = -9999.;
    gen_ipx_2_ = -9999.;
    gen_ipy_2_ = -9999.;
    gen_ipz_2_ = -9999.;

    aco_angle_1_=-9999.;
    aco_angle_2_=-9999.;
    aco_angle_3_=-9999.;
    aco_angle_4_=-9999.;
    aco_angle_5_=-9999.;
    aco_angle_6_=-9999.;
    aco_angle_7_=-9999.;
    pv_angle_=-9999.;
    sv_x_2_=0.;
    sv_y_2_=-0.;
    sv_z_2_=0.;
    sv_mag_2_=0.;
    sv_dphi_2_=0.;
    sv_deta_2_=0.;
    hasSV_2_=false;
 
    y_1_1_=-9999.;
    y_1_2_=-9999.;
    y_2_2_=-9999.;
    y_3_2_=-9999.;
    y_4_2_=-9999.;

    pi_px_1_=0.;
    pi_py_1_=0.;
    pi_pz_1_=0.;
    pi_E_1_=0.;
    pi_px_2_=0.;
    pi_py_2_=0.;
    pi_pz_2_=0.;
    pi_E_2_=0.;
    pi2_px_1_=0.;
    pi2_py_1_=0.;
    pi2_pz_1_=0.;
    pi2_E_1_=0.;
    pi2_px_2_=0.;
    pi2_py_2_=0.;
    pi2_pz_2_=0.;
    pi2_E_2_=0.;
    pi3_px_1_=0.;
    pi3_py_1_=0.;
    pi3_pz_1_=0.;
    pi3_E_1_=0.;
    pi3_px_2_=0.;
    pi3_py_2_=0.;
    pi3_pz_2_=0.;
    pi3_E_2_=0.;

    pi0_px_1_=0.;
    pi0_py_1_=0.;
    pi0_pz_1_=0.;
    pi0_E_1_=0.;
    pi0_px_2_=0.;
    pi0_py_2_=0.;
    pi0_pz_2_=0.;
    pi0_E_2_=0.;
    gam1_px_1_=0.;
    gam1_py_1_=0.;
    gam1_pz_1_=0.;
    gam1_E_1_=0.;
    gam1_px_2_=0.;
    gam1_py_2_=0.;
    gam1_pz_2_=0.;
    gam1_E_2_=0.;
    gam2_px_1_=0.;
    gam2_py_1_=0.;
    gam2_pz_1_=0.;
    gam2_E_1_=0.;
    gam2_px_2_=0.;
    gam2_py_2_=0.;
    gam2_pz_2_=0.;
    gam2_E_2_=0.;
    gam3_px_1_=0.;
    gam3_py_1_=0.;
    gam3_pz_1_=0.;
    gam3_E_1_=0.;
    gam3_px_2_=0.;
    gam3_py_2_=0.;
    gam3_pz_2_=0.;
    gam3_E_2_=0.;
    gam4_px_1_=0.;
    gam4_py_1_=0.;
    gam4_pz_1_=0.;
    gam4_E_1_=0.;
    gam4_px_2_=0.;
    gam4_py_2_=0.;
    gam4_pz_2_=0.;
    gam4_E_2_=0.;
    metx_= mets->vector().Px();
    mety_= mets->vector().Py();
    ip_x_1_=0.;
    ip_y_1_=0.;
    ip_z_1_=0.;
    ip_x_2_=0.;
    ip_y_2_=0.;
    ip_z_2_=0.;

    sv_x_1_=0.;
    sv_y_1_=0.;
    sv_z_1_=0.;
    sv_x_2_=0.;
    sv_y_2_=0.;
    sv_z_2_=0.;

    ipcov00_1_=0;
    ipcov01_1_=0;
    ipcov02_1_=0;
    ipcov10_1_=0;
    ipcov11_1_=0;
    ipcov12_1_=0;
    ipcov20_1_=0;
    ipcov21_1_=0;
    ipcov22_1_=0;
    ipcov00_2_=0;
    ipcov01_2_=0;
    ipcov02_2_=0;
    ipcov10_2_=0;
    ipcov11_2_=0;
    ipcov12_2_=0;
    ipcov20_2_=0;
    ipcov21_2_=0;
    ipcov22_2_=0;

    svcov00_1_=0;
    svcov01_1_=0;
    svcov02_1_=0;
    svcov10_1_=0;
    svcov11_1_=0;
    svcov12_1_=0;
    svcov20_1_=0;
    svcov21_1_=0;
    svcov22_1_=0;
    svcov00_2_=0;
    svcov01_2_=0;
    svcov02_2_=0;
    svcov10_2_=0;
    svcov11_2_=0;
    svcov12_2_=0;
    svcov20_2_=0;
    svcov21_2_=0;
    svcov22_2_=0;
	
    q_pi_1_ = 0.;
    q_pi2_1_ = 0.;	
    q_pi3_1_ = 0.;	
    q_pi_2_ = 0.;
    q_pi2_2_ = 0.;	
    q_pi3_2_ = 0.;	

    use_refitted_vertex_ = false;

    TVector3 ipgen1(0.,0.,0.);
    TVector3 ipgen2(0.,0.,0.);

    if(event->Exists("gen_ip_1")) ipgen1 = event->Get<TVector3>("gen_ip_1");
    if(event->Exists("gen_ip_2")) ipgen2 = event->Get<TVector3>("gen_ip_2");

    a1_flag_ = false;
    if (channel_ == channel::tt && event->ExistsInTree("pfCandidates")) {

      Tau const* tau1 = dynamic_cast<Tau const*>(lep1);
      Tau const* tau2 = dynamic_cast<Tau const*>(lep2);

      std::vector<ic::Vertex*> & vertex_vec = event->GetPtrVec<ic::Vertex>("vertices");
      //std::vector<ic::Vertex*> & refit_vertex_vec = event->GetPtrVec<ic::Vertex>("refittedVertices");
      std::vector<ic::Vertex*> & refit_vertex_vec = event->GetPtrVec<ic::Vertex>("refittedVerticesBS");
      ic::Vertex* refit_vertex = vertex_vec[0];
      for(auto v : refit_vertex_vec) {
        if(v->id() == tau1->id()+tau2->id()) {
          refit_vertex = v;
          use_refitted_vertex_ = true;
        }
      }

      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pv_cov = refit_vertex->covariance();

      lead_pt_1_ =  tau1->lead_pt();
      lead_pt_2_ =  tau2->lead_pt();
     
      tau_mva_decay_mode_1_ = tau1->HasTauID(MVADM2017) ? tau1->GetTauID(MVADM2017) : 0.0;
      tau_mva_decay_mode_2_ = tau2->HasTauID(MVADM2017) ? tau2->GetTauID(MVADM2017) : 0.0;

      std::vector<ic::PFCandidate*> pfcands =  event->GetPtrVec<ic::PFCandidate>("pfCandidates");
      std::pair<ic::Candidate*, ic::Candidate*> rho1 = GetRho(tau1, pfcands);
      ic::Candidate *pi_tau1 = rho1.first;
      ic::Candidate *pi0_tau1 = rho1.second;


      std::pair<ic::Candidate*, ic::Candidate*> rho2 = GetRho(tau2, pfcands);
      ic::Candidate *pi_tau2 = rho2.first;
      ic::Candidate *pi0_tau2 = rho2.second;

      std::vector<ic::PFCandidate*> hads1 = GetTauHads(tau1, pfcands,0);
      std::vector<ic::PFCandidate*> gammas1 = GetTauGammas(tau1, pfcands, 1., 0);
      std::vector<ic::PFCandidate*> hads2 = GetTauHads(tau2, pfcands,0);
      std::vector<ic::PFCandidate*> gammas2 = GetTauGammas(tau2, pfcands, 1., 0);
      if(hads1.size()>2) hads1 = GetA1 (tau1, pfcands).first;
      if(hads2.size()>2) hads2 = GetA1 (tau2, pfcands).first;

      if(tau1->hasSV()&&hads1.size()==3) {
          q_pi_1_ = (double)hads1[0]->charge();
          q_pi2_1_ = (double)hads1[1]->charge();
          q_pi3_1_ = (double)hads1[2]->charge();		  
       }
	  
      if(tau2->hasSV()&&hads2.size()==3) {
          q_pi_2_ = (double)hads2[0]->charge();
          q_pi2_2_ = (double)hads2[1]->charge();
          q_pi3_2_ = (double)hads2[2]->charge();		  
      }	  

      if(tau_decay_mode_1_<3.) {
        pi_px_1_=pi_tau1->vector().Px();
        pi_py_1_=pi_tau1->vector().Py();
        pi_pz_1_=pi_tau1->vector().Pz();
        pi_E_1_=pi_tau1->vector().E();
      } else if(hads1.size()>2) {
        pi_px_1_=hads1[0]->vector().Px();
        pi_py_1_=hads1[0]->vector().Py();
        pi_pz_1_=hads1[0]->vector().Pz();
        pi_E_1_=hads1[0]->vector().E();

        pi2_px_1_=hads1[1]->vector().Px();
        pi2_py_1_=hads1[1]->vector().Py();
        pi2_pz_1_=hads1[1]->vector().Pz();
        pi2_E_1_=hads1[1]->vector().E();

        pi3_px_1_=hads1[2]->vector().Px();
        pi3_py_1_=hads1[2]->vector().Py();
        pi3_pz_1_=hads1[2]->vector().Pz();
        pi3_E_1_=hads1[2]->vector().E();

        if(tau1->hasSV()){
          sv_x_1_ = tau1->secondary_vertex().X() - refit_vertex->vx();
          sv_y_1_ = tau1->secondary_vertex().Y() - refit_vertex->vy();
          sv_z_1_ = tau1->secondary_vertex().Z() - refit_vertex->vz();

          ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > sv_cov = tau1->s_vtx_covariance()+pv_cov;

          svcov00_1_ = sv_cov(0,0);
          svcov01_1_ = sv_cov(0,1);
          svcov02_1_ = sv_cov(0,2);
          svcov10_1_ = sv_cov(1,0);
          svcov11_1_ = sv_cov(1,1);
          svcov12_1_ = sv_cov(1,2);
          svcov20_1_ = sv_cov(2,0);
          svcov21_1_ = sv_cov(2,1);
          svcov22_1_ = sv_cov(2,2);
        }
       
      }
      if(tau_decay_mode_2_<3.) {
        pi_px_2_=pi_tau2->vector().Px();
        pi_py_2_=pi_tau2->vector().Py();
        pi_pz_2_=pi_tau2->vector().Pz();
        pi_E_2_=pi_tau2->vector().E();
      } else if(hads2.size()>2) {
        pi_px_2_=hads2[0]->vector().Px();
        pi_py_2_=hads2[0]->vector().Py();
        pi_pz_2_=hads2[0]->vector().Pz();
        pi_E_2_=hads2[0]->vector().E();

        pi2_px_2_=hads2[1]->vector().Px();
        pi2_py_2_=hads2[1]->vector().Py();
        pi2_pz_2_=hads2[1]->vector().Pz();
        pi2_E_2_=hads2[1]->vector().E();

        pi3_px_2_=hads2[2]->vector().Px();
        pi3_py_2_=hads2[2]->vector().Py();
        pi3_pz_2_=hads2[2]->vector().Pz();
        pi3_E_2_=hads2[2]->vector().E();

        if(tau2->hasSV()) {
          sv_x_2_ = tau2->secondary_vertex().X() - refit_vertex->vx();
          sv_y_2_ = tau2->secondary_vertex().Y() - refit_vertex->vy();
          sv_z_2_ = tau2->secondary_vertex().Z() - refit_vertex->vz();

          ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > sv_cov = tau2->s_vtx_covariance()+pv_cov;

          svcov00_2_ = sv_cov(0,0);
          svcov01_2_ = sv_cov(0,1);
          svcov02_2_ = sv_cov(0,2);
          svcov10_2_ = sv_cov(1,0);
          svcov11_2_ = sv_cov(1,1);
          svcov12_2_ = sv_cov(1,2);
          svcov20_2_ = sv_cov(2,0);
          svcov21_2_ = sv_cov(2,1);
          svcov22_2_ = sv_cov(2,2);
        }
      }
      if(gammas1.size()>0 && tau_decay_mode_1_<2){
        pi0_px_1_=pi0_tau1->vector().Px();
        pi0_py_1_=pi0_tau1->vector().Py();
        pi0_pz_1_=pi0_tau1->vector().Pz();
        pi0_E_1_=pi0_tau1->vector().E();
        gam1_px_1_=gammas1[0]->vector().Px();
        gam1_py_1_=gammas1[0]->vector().Py();
        gam1_pz_1_=gammas1[0]->vector().Pz();
        gam1_E_1_=gammas1[0]->vector().E();
        if(gammas1.size()>1) gam2_px_1_= gammas1[1]->vector().Px();
        if(gammas1.size()>1) gam2_py_1_= gammas1[1]->vector().Py();
        if(gammas1.size()>1) gam2_pz_1_= gammas1[1]->vector().Pz();
        if(gammas1.size()>1) gam2_E_1_= gammas1[1]->vector().E();
        if(gammas1.size()>2) gam3_px_1_= gammas1[2]->vector().Px();
        if(gammas1.size()>2) gam3_py_1_= gammas1[2]->vector().Py();
        if(gammas1.size()>2) gam3_pz_1_= gammas1[2]->vector().Pz();
        if(gammas1.size()>2) gam3_E_1_= gammas1[2]->vector().E();
        if(gammas1.size()>3) gam4_px_1_= gammas1[3]->vector().Px();
        if(gammas1.size()>3) gam4_py_1_= gammas1[3]->vector().Py();
        if(gammas1.size()>3) gam4_pz_1_= gammas1[3]->vector().Pz();
        if(gammas1.size()>3) gam4_E_1_= gammas1[3]->vector().E();
      }
      if(gammas2.size()>0 && tau_decay_mode_2_<2) {
        pi0_px_2_=pi0_tau2->vector().Px();
        pi0_py_2_=pi0_tau2->vector().Py();
        pi0_pz_2_=pi0_tau2->vector().Pz();
        pi0_E_2_=pi0_tau2->vector().E();

        gam1_px_2_=gammas2[0]->vector().Px();
        gam1_py_2_=gammas2[0]->vector().Py();
        gam1_pz_2_=gammas2[0]->vector().Pz();
        gam1_E_2_=gammas2[0]->vector().E();
        if(gammas2.size()>1) gam2_px_2_= gammas2[1]->vector().Px();
        if(gammas2.size()>1) gam2_py_2_= gammas2[1]->vector().Py();
        if(gammas2.size()>1) gam2_pz_2_= gammas2[1]->vector().Pz();
        if(gammas2.size()>1) gam2_E_2_= gammas2[1]->vector().E();
        if(gammas2.size()>2) gam3_px_2_= gammas2[2]->vector().Px();
        if(gammas2.size()>2) gam3_py_2_= gammas2[2]->vector().Py();
        if(gammas2.size()>2) gam3_pz_2_= gammas2[2]->vector().Pz();
        if(gammas2.size()>2) gam3_E_2_= gammas2[2]->vector().E();
        if(gammas2.size()>3) gam4_px_2_= gammas2[3]->vector().Px();
        if(gammas2.size()>3) gam4_py_2_= gammas2[3]->vector().Py();
        if(gammas2.size()>3) gam4_pz_2_= gammas2[3]->vector().Pz();
        if(gammas2.size()>3) gam4_E_2_= gammas2[3]->vector().E();
      }

      TLorentzVector lvec1;
      TLorentzVector lvec2;
      TLorentzVector lvec3;
      TLorentzVector lvec4;
      TLorentzVector pvtosv;

      std::pair<TVector3,double> ipandsig_1 = IPAndSignificance(tau1, refit_vertex, pfcands);
      std::pair<TVector3,double> ipandsig_2 = IPAndSignificance(tau2, refit_vertex, pfcands);
      ip_mag_1_ = ipandsig_1.first.Mag();
      ip_sig_1_ = ipandsig_1.second;
      ip_mag_2_ = ipandsig_2.first.Mag();
      ip_sig_2_ = ipandsig_2.second;

      TVector3 ip_corrected_1 = ipandsig_1.first;
      TVector3 ip_corrected_2 = ipandsig_2.first;

      std::pair<TVector3,ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> ipandcov_1 = IPAndCov(tau1, refit_vertex, pfcands);
      std::pair<TVector3,ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> ipandcov_2 = IPAndCov(tau2, refit_vertex, pfcands);
      CovMatrix cov_corrected_1 = ipandcov_1.second;
      CovMatrix cov_corrected_2 = ipandcov_2.second;


      // add this part for the mt channel as well!!!!!
      if(!is_data_) {

        ip_sig_1_raw_= ip_sig_1_;
        ip_sig_2_raw_= ip_sig_2_;

        ip_corrected_1 = ipCorrector.correctIp(
                         ipandsig_1.first,
                         ipgen1,
                         fabs(lep1->eta())
        );

        ip_corrected_2 = ipCorrector.correctIp(
                         ipandsig_2.first,
                         ipgen2,
                         fabs(lep2->eta())
        );

      
        // Correct IP covariance matrix
        cov_corrected_1 = ipCorrector.correctIpCov(
                    ipandcov_1.second, 
                    fabs(lep1->eta())
        );

        cov_corrected_2 = ipCorrector.correctIpCov(
                    ipandcov_2.second,
                    fabs(lep2->eta())
        );


        ROOT::Math::SVector<double, 3> ip_svec_1;
        ip_svec_1(0) = ip_corrected_1.X();
        ip_svec_1(1) = ip_corrected_1.Y();
        ip_svec_1(2) = ip_corrected_1.Z();
        ip_svec_1 = ip_svec_1.Unit();
        ip_sig_1_ = ip_corrected_1.Mag()/sqrt(ROOT::Math::Dot( ip_svec_1, cov_corrected_1 * ip_svec_1));


        ROOT::Math::SVector<double, 3> ip_svec_2;
        ip_svec_2(0) = ip_corrected_2.X();
        ip_svec_2(1) = ip_corrected_2.Y();
        ip_svec_2(2) = ip_corrected_2.Z();
        ip_svec_2 = ip_svec_2.Unit();
        ip_sig_2_ = ip_corrected_2.Mag()/sqrt(ROOT::Math::Dot( ip_svec_2, cov_corrected_2 * ip_svec_2));

        ip_sig_1_down_= ip_sig_1_ -0.25*(ip_sig_1_-ip_sig_1_raw_);
        ip_sig_1_up_  = ip_sig_1_ +0.25*(ip_sig_1_-ip_sig_1_raw_);
        ip_sig_2_down_= ip_sig_2_ -0.25*(ip_sig_2_-ip_sig_2_raw_);
        ip_sig_2_up_  = ip_sig_2_ +0.25*(ip_sig_2_-ip_sig_2_raw_);

      } else {
        ip_sig_1_raw_= ip_sig_1_;
        ip_sig_2_raw_= ip_sig_2_;
        ip_sig_1_up_= ip_sig_1_;
        ip_sig_2_up_= ip_sig_2_;
        ip_sig_1_down_= ip_sig_1_;
        ip_sig_2_down_= ip_sig_2_;
      }

      ip_x_1_ = ip_corrected_1.X();
      ip_y_1_ = ip_corrected_1.Y();
      ip_z_1_ = ip_corrected_1.Z();
      ip_x_2_ = ip_corrected_2.X();
      ip_y_2_ = ip_corrected_2.Y();
      ip_z_2_ = ip_corrected_2.Z();

      ipcov00_1_ = cov_corrected_1(0,0);
      ipcov01_1_ = cov_corrected_1(0,1);
      ipcov02_1_ = cov_corrected_1(0,2);
      ipcov10_1_ = cov_corrected_1(1,0);
      ipcov11_1_ = cov_corrected_1(1,1);
      ipcov12_1_ = cov_corrected_1(1,2);
      ipcov20_1_ = cov_corrected_1(2,0);
      ipcov21_1_ = cov_corrected_1(2,1);
      ipcov22_1_ = cov_corrected_1(2,2);

      ipcov00_2_ = cov_corrected_2(0,0);
      ipcov01_2_ = cov_corrected_2(0,1);
      ipcov02_2_ = cov_corrected_2(0,2);
      ipcov10_2_ = cov_corrected_2(1,0);
      ipcov11_2_ = cov_corrected_2(1,1);
      ipcov12_2_ = cov_corrected_2(1,2);
      ipcov20_2_ = cov_corrected_2(2,0);
      ipcov21_2_ = cov_corrected_2(2,1);
      ipcov22_2_ = cov_corrected_2(2,2);

      // uncorrected IP
      //TVector3 ip1 = (ipandsig_1.first).Unit();
      //TVector3 ip2 = (ipandsig_2.first).Unit();

      // corrected IP
      TVector3 ip1 = ip_corrected_1.Unit();
      TVector3 ip2 = ip_corrected_2.Unit();

      std::vector<ic::PFCandidate*> charged_cands_1 = GetTauHads(tau1, pfcands,0);
      std::vector<ic::PFCandidate*> charged_cands_2 = GetTauHads(tau2, pfcands,0);
      q_tot_1_=0., q_tot_2_=0.;

      for(auto c : charged_cands_1) q_tot_1_+=c->charge();
      for(auto c : charged_cands_2) q_tot_2_+=c->charge();

      auto primary_vtx = refit_vertex;

      if(tau_decay_mode_1_==0&&tau_decay_mode_2_==0) {

        TLorentzVector pvtosv1(
                tau1->svx() - primary_vtx->vx(),
                tau1->svy() - primary_vtx->vy(),
                tau1->svz() - primary_vtx->vz(),
                0.);
        TLorentzVector pvtosv2(
                tau2->svx() - primary_vtx->vx(),
                tau2->svy() - primary_vtx->vy(),
                tau2->svz() - primary_vtx->vz(),
                0.);

        lvec3 = ConvertToLorentz(tau1->vector()); 
        lvec4 = ConvertToLorentz(tau2->vector()); 

        lvec1 = TLorentzVector(ip1, 0.);
        lvec2 = TLorentzVector(ip2, 0.);

        aco_angle_6_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
        //double test = IPAcoAngle_Test(lvec1, lvec2, lvec3, lvec4,false); 
        //std::cout << " " << std::endl;
        //std::cout << ip_x_1_ << " " << ip_y_1_ << " " << ip_z_1_ << " " << lvec3.X() << " " << lvec3.Y() << " " << lvec3.Z() << " " << lvec3.E() << std::endl;
        //std::cout << ip_x_1_ << " " << ip_y_1_ << " " << ip_z_1_ << " " << pi_px_1_ << " " << pi_py_1_  << " " << pi_pz_1_  << " " << pi_E_1_  << std::endl;
        //std::cout << ip_x_2_ << " " << ip_y_2_ << " " << ip_z_2_ << " " << lvec4.X() << " " << lvec4.Y() << " " << lvec4.Z() << " " << lvec4.E() << std::endl;
        //std::cout << ip_x_2_ << " " << ip_y_2_ << " " << ip_z_2_ << " " << pi_px_2_ << " " << pi_py_2_  << " " << pi_pz_2_  << " " << pi_E_2_  << std::endl;

        //std::cout << test << std::endl;
        aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

        if(event_ % 2) {
          aco_angle_rand_ = aco_angle_6_;
          aco_sign_rand_ = aco_sign_;
        } else {
          aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
          aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
        }

        alpha1_1_ = AlphaAngle(lvec3.Vect(), ip1);
        alpha1_2_ = AlphaAngle(lvec4.Vect(), ip2);

      }

      if((tau_decay_mode_1_==1&&tau_decay_mode_2_==0) || (tau_decay_mode_1_==0&&tau_decay_mode_2_==1)) {
               
        ic::Candidate *pi;
        ic::Candidate *pi0;
        ic::Tau const *tau;
        ic::Tau const *tau_1;

        TVector3 ip;
        TVector3 ip_alt;
        if(tau_decay_mode_1_==1&&tau_decay_mode_2_==0) {
          pi = pi_tau1;
          pi0 = pi0_tau1;
          tau = tau2;
          tau_1 = tau1;
          ip = ip2;
          ip_alt = ip1;
        } else {
          pi = pi_tau2;
          pi0 = pi0_tau2;
          tau = tau1;
          tau_1 = tau2;
          ip = ip1;
          ip_alt = ip2;
        }
        double cp_sign_ = YRho(std::vector<Candidate*>({pi, pi0}),TVector3());
        y_1_1_ = YRho(std::vector<Candidate*>({pi, pi0}),TVector3());
        TLorentzVector pvtosv1(
                tau->svx() - primary_vtx->vx(),
                tau->svy() - primary_vtx->vy(),
                tau->svz() - primary_vtx->vz(),
                0.);
        lvec1 = ConvertToLorentz(pi0->vector()); //pi zero from rho
        lvec3 = ConvertToLorentz(pi->vector()); //pi charge from rho
        lvec4 = ConvertToLorentz(tau->vector()); //pi charge from tau

        //TVector3 ip = (pvtosv1.Vect() - pvtosv1.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
        lvec2 = TLorentzVector(ip, 0.);

        TLorentzVector lvec3_1 = ConvertToLorentz(pi->vector());

        TLorentzVector pvtosv2(
                tau_1->svx() - primary_vtx->vx(),
                tau_1->svy() - primary_vtx->vy(),
                tau_1->svz() - primary_vtx->vz(),
                0.);

        //TVector3 ip1 = (pvtosv2.Vect() - pvtosv2.Vect().Dot(lvec3_1.Vect().Unit())*lvec3_1.Vect().Unit()).Unit();
        TLorentzVector lvec1_1 = TLorentzVector(ip_alt, 0.);


        aco_angle_5_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
        aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);
        if (cp_sign_<0) {
          if (aco_angle_5_<M_PI)  aco_angle_5_ = aco_angle_5_+M_PI;
          else                    aco_angle_5_ = aco_angle_5_-M_PI;
        }

        if(event_ % 2) {
          aco_angle_rand_ = aco_angle_5_;
          aco_sign_rand_ = aco_sign_;
        } else {
          aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
          aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
          if (cp_sign_<0) {
            if (aco_angle_rand_<M_PI) aco_angle_rand_ = aco_angle_rand_+M_PI;
            else                      aco_angle_rand_ = aco_angle_rand_-M_PI;
          }
        }

        aco_angle_6_ = IPAcoAngle(lvec1_1, lvec2, lvec3_1, lvec4,false);
        if (cp_sign_<0) {
          if (aco_angle_6_<M_PI)  aco_angle_6_ = aco_angle_6_+M_PI;
          else                    aco_angle_6_ = aco_angle_6_-M_PI;
        }

        alpha1_1_ = AlphaAngle(lvec4.Vect(), ip);
        alpha2_2_ = AlphaAngleRho(lvec3.Vect(), lvec1.Vect());

      }

      if(tau_decay_mode_1_==1&&tau_decay_mode_2_==1){
        lvec1 = ConvertToLorentz(pi0_tau1->vector());
        lvec2 = ConvertToLorentz(pi0_tau2->vector());
        lvec3 = ConvertToLorentz(pi_tau1->vector());
        lvec4 = ConvertToLorentz(pi_tau2->vector());
        double cp_sign_ = YRho(std::vector<Candidate*>({pi0_tau1, pi_tau1}),TVector3())
            *YRho(std::vector<Candidate*>({pi0_tau2, pi_tau2}),TVector3());
        y_1_1_ = YRho(std::vector<Candidate*>({pi0_tau1, pi_tau1}),TVector3());
        y_1_2_ = YRho(std::vector<Candidate*>({pi0_tau2, pi_tau2}),TVector3());

        //double cp_sign_1_ = YRho(std::vector<Candidate*>({pi0_tau1, pi_tau1}),TVector3()); 
        //double cp_sign_2_ = YRho(std::vector<Candidate*>({pi0_tau2, pi_tau2}),TVector3()); 

        aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
        aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

        TLorentzVector lvec1_1 = TLorentzVector(ip1, 0.);
        TLorentzVector lvec2_1 = TLorentzVector(ip2, 0.);

        aco_angle_6_ = IPAcoAngle(lvec1_1, lvec2_1, lvec3, lvec4,false);
        aco_angle_5_ = IPAcoAngle(lvec1_1, lvec2, lvec3, lvec4,false);
        aco_angle_7_ = IPAcoAngle(lvec1, lvec2_1, lvec3, lvec4,false);

        if (cp_sign_<0) {
          if (aco_angle_6_<M_PI)  aco_angle_6_ = aco_angle_6_+M_PI;
          else                    aco_angle_6_ = aco_angle_6_-M_PI;
        }

        if (cp_sign_<0) {
          if (aco_angle_1_<M_PI) aco_angle_1_ += M_PI;
          else                   aco_angle_1_ -= M_PI;
        }

        if (cp_sign_<0) {
          if (aco_angle_5_<M_PI) aco_angle_5_ += M_PI;
          else                   aco_angle_5_ -= M_PI;
        }
        if (cp_sign_<0) {
          if (aco_angle_7_<M_PI) aco_angle_7_ += M_PI;
          else                   aco_angle_7_ -= M_PI;
        }

        if(event_ % 2) {
          aco_angle_rand_ = aco_angle_1_;
          aco_sign_rand_ = aco_sign_;
        } else {
          aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
          aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
          if (cp_sign_<0) {
            if (aco_angle_rand_<M_PI) aco_angle_rand_ = aco_angle_rand_+M_PI;
            else                      aco_angle_rand_ = aco_angle_rand_-M_PI;
          }
        } 

        alpha2_1_ = AlphaAngleRho(lvec3.Vect(), lvec1.Vect());
        alpha2_2_ = AlphaAngleRho(lvec4.Vect(), lvec2.Vect()); 
      }

      else if((tau_decay_mode_1_==0&&tau_decay_mode_2_>=10) || (tau_decay_mode_1_>=10&&tau_decay_mode_2_==0)){


        TVector3 ip;

        ic::Tau const *tau_1;
        std::vector<ic::PFCandidate*> a1_daughters;
        if(tau_decay_mode_1_>=10) {
          a1_daughters  = GetA1(tau1, pfcands).first;
          tau_1 = tau2;
          ip = ip2;
        } else {
          a1_daughters  = GetA1(tau2, pfcands).first;
          tau_1 = tau1;
          ip = ip1;
        }

        TLorentzVector pvtosv1(
                tau_1->svx() - primary_vtx->vx(),
                tau_1->svy() - primary_vtx->vy(),
                tau_1->svz() - primary_vtx->vz(),
                0.);
        lvec3 = ConvertToLorentz(tau_1->vector()); //pi charge from tau

//        TVector3 ip = (pvtosv1.Vect() - pvtosv1.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();

        lvec1 = TLorentzVector(ip, 0.);

        if (a1_daughters.size()>2){

            a1_flag_ = true;
            if(a1_daughters.size()==3) a1_flag_2_ = true;

          

            lvec2 = ConvertToLorentz(a1_daughters[0]->vector()); //pi zero from rho
            lvec4 = ConvertToLorentz(a1_daughters[1]->vector()); //pi charge from rho
            
            aco_angle_5_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);

            aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

            double cp_sign_ = YRho(std::vector<Candidate*>({a1_daughters[0], a1_daughters[1]}),TVector3());        
            if (cp_sign_<0) {
              if (aco_angle_5_<M_PI)  aco_angle_5_ = aco_angle_5_+M_PI;
              else                    aco_angle_5_ = aco_angle_5_-M_PI;
            }
            if(event_ % 2) {
              aco_angle_rand_ = aco_angle_5_;
              aco_sign_rand_ = aco_sign_;
            } else {
              aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
              aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
              if (cp_sign_<0) {
                if (aco_angle_rand_<M_PI) aco_angle_rand_ = aco_angle_rand_+M_PI;
                else                      aco_angle_rand_ = aco_angle_rand_-M_PI;
              }
            }

            TLorentzVector lvec5 = ConvertToLorentz(a1_daughters[2]->vector());

          }

          alpha1_1_ = AlphaAngle(lvec3.Vect(), ip);
          alpha2_2_ = AlphaAngleRho(lvec4.Vect(), lvec2.Vect());

      }


      else if((tau_decay_mode_1_==1&&tau_decay_mode_2_>=10) || (tau_decay_mode_1_>=10&&tau_decay_mode_2_==1)){
      
        ic::Candidate *pi;
        ic::Candidate *pi0;

        std::vector<ic::PFCandidate*> a1_daughters;
        std::pair<ic::Candidate*, ic::Candidate*> rho_daughters;
        if(tau_decay_mode_1_>9){
          a1_daughters  = GetA1(tau1, pfcands).first;
          rho_daughters = GetRho(tau2, pfcands);
          pi0 = pi0_tau2;
          pi = pi_tau2;
        }
        if(tau_decay_mode_2_>9){
          a1_daughters  = GetA1(tau2, pfcands).first;
          rho_daughters = GetRho(tau1, pfcands);
          pi0 = pi0_tau1;
          pi = pi_tau1;

        }
        double  cp_sign_1_=0;
        double  cp_sign_2_=0;
        double  cp_sign_3_=0;
        double  cp_sign_4_=0;

        if (a1_daughters.size()>2){
            a1_flag_ = true;
            if(a1_daughters.size()==3) a1_flag_2_ = true;
            mass0_ = (a1_daughters[0]->vector() + a1_daughters[1]->vector() + a1_daughters[2]->vector()).M();
            mass1_ = (a1_daughters[0]->vector() + a1_daughters[1]->vector()).M();
            mass2_ = (a1_daughters[0]->vector() + a1_daughters[2]->vector()).M();

            //aco_angle_1_ = AcoplanarityAngle(std::vector<Candidate*> ({rho_daughters.first,rho_daughters.second}), std::vector<Candidate*> ({a1_daughters[0],a1_daughters[1]}));

            lvec1 = ConvertToLorentz(pi0->vector());
            lvec3 = ConvertToLorentz(pi->vector());
            lvec2 = ConvertToLorentz(a1_daughters[0]->vector()); //pi zero from rho
            lvec4 = ConvertToLorentz(a1_daughters[1]->vector()); //pi charge from rho

            aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
            aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

            //aco_angle_2_ treats rho like neutral pion 
            lvec2 = ConvertToLorentz(a1_daughters[0]->vector()+a1_daughters[1]->vector());
            lvec4 = ConvertToLorentz(a1_daughters[2]->vector());
            aco_angle_2_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);

            //aco_angle_3 is same as aco_angle_1 but makes other choise for a1->rho
            lvec2 = ConvertToLorentz(a1_daughters[0]->vector()); //pi charge from rho 
            lvec4 = ConvertToLorentz(a1_daughters[2]->vector()); //pi zero from rho
            aco_angle_3_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);

            //aco_angle_4 is same as aco_angle_2 but makes other choise for a1->rho
            lvec2 = ConvertToLorentz(a1_daughters[0]->vector()+a1_daughters[2]->vector());
            lvec4 = ConvertToLorentz(a1_daughters[1]->vector());
            aco_angle_4_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);

            Candidate* rho_1  = new Candidate();
            Candidate* rho_2  = new Candidate();
            rho_1->set_vector(a1_daughters[0]->vector()+a1_daughters[1]->vector());
            rho_2->set_vector(a1_daughters[0]->vector()+a1_daughters[2]->vector());
            //aco_angle_3_ = AcoplanarityAngle(std::vector<Candidate*> ({rho_daughters.first,rho_daughters.second}), std::vector<Candidate*> ({rho_1,a1_daughters[2]}));
            //aco_angle_4_ = AcoplanarityAngle(std::vector<Candidate*> ({rho_daughters.first,rho_daughters.second}), std::vector<Candidate*> ({rho_2,a1_daughters[1]}));

            y_1_1_ = YRho(std::vector<Candidate*>({rho_daughters.first, rho_daughters.second}),TVector3());
            y_1_2_ = YRho(std::vector<Candidate*>({a1_daughters[0], a1_daughters[1]}),TVector3());
            y_3_2_ = YRho(std::vector<Candidate*>({a1_daughters[0], a1_daughters[2]}),TVector3());
            y_2_2_ = YA1(std::vector<Candidate*>({rho_1, a1_daughters[2]}),TVector3());
            y_4_2_ = YA1(std::vector<Candidate*>({rho_2, a1_daughters[1]}),TVector3());

            cp_sign_1_ = y_1_1_*y_1_2_;
            cp_sign_2_ = y_1_1_*y_2_2_;
            cp_sign_3_ = y_1_1_*y_3_2_;
            cp_sign_4_ = y_1_1_*y_4_2_;
        }   

        if (cp_sign_1_<0) {
          if (aco_angle_1_<M_PI) aco_angle_1_ += M_PI;
          else                   aco_angle_1_ -= M_PI;
        }
        if (cp_sign_2_<0) {
          if (aco_angle_2_<M_PI) aco_angle_2_ += M_PI;
          else                   aco_angle_2_ -= M_PI;
        }
        if (cp_sign_3_<0) { 
          if (aco_angle_3_<M_PI) aco_angle_3_ += M_PI;
          else                   aco_angle_3_ -= M_PI;
        }
        if (cp_sign_4_<0) {
          if (aco_angle_4_<M_PI) aco_angle_4_ += M_PI;
          else                   aco_angle_4_ -= M_PI;
        }

        if(event_ % 2) {
          aco_angle_rand_ = aco_angle_1_;
          aco_sign_rand_ = aco_sign_;
        } else {
          aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
          aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
          if (cp_sign_1_<0) {
            if (aco_angle_rand_<M_PI) aco_angle_rand_ = aco_angle_rand_+M_PI;
            else                      aco_angle_rand_ = aco_angle_rand_-M_PI;
          }
        }

        alpha2_1_ = AlphaAngleRho(lvec3.Vect(), lvec1.Vect());
        alpha2_2_ = AlphaAngleRho(lvec4.Vect(), lvec2.Vect());

      } else if(tau_decay_mode_1_>=10&&tau_decay_mode_2_>=10){

        std::vector<ic::PFCandidate*> a1_daughters_1 = GetA1(tau1, pfcands).first;
        std::vector<ic::PFCandidate*> a1_daughters_2 = GetA1(tau2, pfcands).first;

        if(tau1->hasSV()&&tau2->hasSV()&&a1_daughters_1.size()==3&&a1_daughters_2.size()==3) {

          Tau const *tau_neg;
          Tau const *tau_pos;

          if(tau1->charge()<tau2->charge()) {
            tau_neg = tau1;
            tau_pos = tau2;
          } else {
            tau_neg = tau2;
            tau_pos = tau1;
          }

          std::vector<ic::PFCandidate*> a1_daughters_neg = GetA1(tau_neg, pfcands).first;
          std::vector<ic::PFCandidate*> a1_daughters_pos = GetA1(tau_pos, pfcands).first;
          
          TVector3 svminuspv_1(
                tau_neg->secondary_vertex().X() - primary_vtx->vx(),
                tau_neg->secondary_vertex().Y() - primary_vtx->vy(),
                tau_neg->secondary_vertex().Z() - primary_vtx->vz()
          );
          TVector3 svminuspv_2(
                tau_pos->secondary_vertex().X() - primary_vtx->vx(),
                tau_pos->secondary_vertex().Y() - primary_vtx->vy(),
                tau_pos->secondary_vertex().Z() - primary_vtx->vz()
          ); 

          std::vector<TLorentzVector> pis_1 = { 
                ConvertToLorentz(a1_daughters_neg[0]->vector()),
                ConvertToLorentz(a1_daughters_neg[1]->vector()),
                ConvertToLorentz(a1_daughters_neg[2]->vector())
          };

          std::vector<TLorentzVector> pis_2 = {
                ConvertToLorentz(a1_daughters_pos[0]->vector()),
                ConvertToLorentz(a1_daughters_pos[1]->vector()),
                ConvertToLorentz(a1_daughters_pos[2]->vector())
          };

          std::vector<double> charges_1 = {
                (double)a1_daughters_neg[0]->charge(),
                (double)a1_daughters_neg[1]->charge(),
                (double)a1_daughters_neg[2]->charge()
          };

          std::vector<double> charges_2 = {
                (double)a1_daughters_pos[0]->charge(),
                (double)a1_daughters_pos[1]->charge(),
                (double)a1_daughters_pos[2]->charge()
          };

          pv_angle_=PolarimetricA1A1(svminuspv_1, svminuspv_2, ConvertToLorentz(a1_daughters_neg[0]->vector()+a1_daughters_neg[1]->vector()+a1_daughters_neg[2]->vector()), ConvertToLorentz(a1_daughters_pos[0]->vector()+a1_daughters_pos[1]->vector()+a1_daughters_pos[2]->vector()), pis_1, pis_2, charges_1, charges_2);


        } else pv_angle_=-9999;

        if(a1_daughters_1.size()>2 && a1_daughters_2.size()>2) {
          a1_flag_ = true;
          if(a1_daughters_1.size()==3 && a1_daughters_2.size()==3) a1_flag_2_ = true;
          Candidate* rho_1  = new Candidate();
          Candidate* rho_2  = new Candidate();
          rho_1->set_vector(a1_daughters_1[0]->vector()+a1_daughters_1[1]->vector());
          rho_2->set_vector(a1_daughters_2[0]->vector()+a1_daughters_2[1]->vector());

          //aco_angle_1_ = AcoplanarityAngle(std::vector<Candidate*> ({a1_daughters_1[0],a1_daughters_1[1]}), std::vector<Candidate*> ({a1_daughters_2[0],a1_daughters_2[1]}));

          aco_angle_2_ = AcoplanarityAngle(std::vector<Candidate*> ({a1_daughters_1[0],a1_daughters_1[1]}), std::vector<Candidate*> ({rho_2,a1_daughters_2[2]}));

          aco_angle_3_ = AcoplanarityAngle(std::vector<Candidate*> ({a1_daughters_2[0],a1_daughters_2[1]}), std::vector<Candidate*> ({rho_1,a1_daughters_1[2]}));

          aco_angle_4_ = AcoplanarityAngle(std::vector<Candidate*> ({rho_1,a1_daughters_1[2]}), std::vector<Candidate*> ({rho_2,a1_daughters_2[2]})); 

          lvec1 = ConvertToLorentz(a1_daughters_1[0]->vector());
          lvec3 = ConvertToLorentz(a1_daughters_1[1]->vector());
          lvec2 = ConvertToLorentz(a1_daughters_2[0]->vector()); //pi zero from rho
          lvec4 = ConvertToLorentz(a1_daughters_2[1]->vector()); //pi charge from rho

          aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
          aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

          double yrho_1_ = YRho(std::vector<Candidate*>({a1_daughters_1[0],a1_daughters_1[1]}),TVector3());
          double yrho_2_ = YRho(std::vector<Candidate*>({a1_daughters_2[0],a1_daughters_2[1]}),TVector3());

          double ya1_1_ = YA1(std::vector<Candidate*>({rho_1, a1_daughters_1[2]}),TVector3());
          double ya1_2_ = YA1(std::vector<Candidate*>({rho_2, a1_daughters_2[2]}),TVector3());

          double cp_sign_1_=yrho_1_*yrho_2_;
          double cp_sign_2_=yrho_1_*ya1_2_;
          double cp_sign_3_=yrho_2_*ya1_1_;
          double cp_sign_4_=ya1_1_*ya1_2_;

          y_1_1_ = yrho_1_; 
          y_1_2_ = yrho_2_;

          //y_3_1_ =   
          //y_2_1_ =   
          //y_4_1_ = 

          //y_3_2_ = 
          //y_2_2_ = 
          //y_4_2_ = 

          //  y_1_1_ = YRho(std::vector<Candidate*>({rho_daughters.first, rho_daughters.second}),TVector3());
          //  y_1_2_ = YRho(std::vector<Candidate*>({a1_daughters[0], a1_daughters[1]}),TVector3());
          //  y_3_2_ = YRho(std::vector<Candidate*>({a1_daughters[0], a1_daughters[2]}),TVector3());
          //  y_2_2_ = YA1(std::vector<Candidate*>({rho_1, a1_daughters[2]}),TVector3());
          //  y_4_2_ = YA1(std::vector<Candidate*>({rho_2, a1_daughters[1]}),TVector3());

          if (cp_sign_1_<0) {
            if (aco_angle_1_<M_PI) aco_angle_1_ += M_PI;
            else                   aco_angle_1_ -= M_PI;
          }
          if (cp_sign_2_<0) {
            if (aco_angle_2_<M_PI) aco_angle_2_ += M_PI;
            else                   aco_angle_2_ -= M_PI;
          }
          if (cp_sign_3_<0) {
            if (aco_angle_3_<M_PI) aco_angle_3_ += M_PI;
            else                   aco_angle_3_ -= M_PI;
          }
          if (cp_sign_4_<0) {
            if (aco_angle_4_<M_PI) aco_angle_4_ += M_PI;
            else                   aco_angle_4_ -= M_PI;
          }

          if(event_ % 2) {
            aco_angle_rand_ = aco_angle_1_;
            aco_sign_rand_ = aco_sign_;
          } else {
            aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
            aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
            if (cp_sign_1_<0) {
              if (aco_angle_rand_<M_PI) aco_angle_rand_ = aco_angle_rand_+M_PI;
              else                      aco_angle_rand_ = aco_angle_rand_-M_PI;
            }
          }

          alpha2_1_ = AlphaAngleRho(lvec3.Vect(), lvec1.Vect());
          alpha2_2_ = AlphaAngleRho(lvec4.Vect(), lvec2.Vect());
        }
      }

      // save IP for sync 
      ipx_1_ = ip1.X();
      ipy_1_ = ip1.Y();
      ipz_1_ = ip1.Z();
      ipx_2_ = ip2.X();
      ipy_2_ = ip2.Y();
      ipz_2_ = ip2.Z();

      // gen IP
      gen_ipx_1_ = ipgen1.X();
      gen_ipy_1_ = ipgen1.Y();
      gen_ipz_1_ = ipgen1.Z();
      gen_ipx_2_ = ipgen2.X();
      gen_ipy_2_ = ipgen2.Y();
      gen_ipz_2_ = ipgen2.Z();

      // primary vertex
      pvx_ = primary_vtx->vx();
      pvy_ = primary_vtx->vy();
      pvz_ = primary_vtx->vz();

    }
    else if (channel_ == channel::mt && event->ExistsInTree("pfCandidates")) {
      Muon const* muon1 = dynamic_cast<Muon const*>(lep1);
      Tau const* tau2 = dynamic_cast<Tau const*>(lep2);

      std::vector<ic::Vertex*> & vertex_vec = event->GetPtrVec<ic::Vertex>("vertices");
      //std::vector<ic::Vertex*> & refit_vertex_vec = event->GetPtrVec<ic::Vertex>("refittedVertices");
      std::vector<ic::Vertex*> & refit_vertex_vec = event->GetPtrVec<ic::Vertex>("refittedVerticesBS");
      ic::Vertex* refit_vertex = new ic::Vertex();
      if(vertex_vec.size()>0) refit_vertex = vertex_vec[0];
      for(auto v : refit_vertex_vec) {
        if(v->id() == muon1->id()+tau2->id()) {
          refit_vertex = v;
          use_refitted_vertex_ = true;
        }
      }

      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pv_cov = refit_vertex->covariance();

      lead_pt_2_ =  tau2->lead_pt();

      tau_mva_decay_mode_2_ = tau2->HasTauID(MVADM2017) ? tau2->GetTauID(MVADM2017) : 0.0;


      std::vector<ic::PFCandidate*> pfcands =  event->GetPtrVec<ic::PFCandidate>("pfCandidates");
      std::pair<ic::Candidate*, ic::Candidate*> rho = GetRho(tau2, pfcands);
      ic::Candidate *pi_tau2 = rho.first;
      ic::Candidate *pi0_tau2 = rho.second;
     
      std::vector<ic::PFCandidate*> hads2 = GetTauHads(tau2, pfcands,0);
      std::vector<ic::PFCandidate*> gammas2 = GetTauGammas(tau2, pfcands, 1., 0);
      if(hads2.size()>2) hads2 = GetA1 (tau2, pfcands).first;

      pi_px_1_=lep1->vector().Px();
      pi_py_1_=lep1->vector().Py();
      pi_pz_1_=lep1->vector().Pz();
      pi_E_1_=lep1->vector().E();
      if(tau_decay_mode_2_<3.) {
        pi_px_2_=pi_tau2->vector().Px();
        pi_py_2_=pi_tau2->vector().Py();
        pi_pz_2_=pi_tau2->vector().Pz();
        pi_E_2_=pi_tau2->vector().E();
      } else if(hads2.size()>2) {
        pi_px_2_=hads2[0]->vector().Px();
        pi_py_2_=hads2[0]->vector().Py();
        pi_pz_2_=hads2[0]->vector().Pz();
        pi_E_2_=hads2[0]->vector().E();

        pi2_px_2_=hads2[1]->vector().Px();
        pi2_py_2_=hads2[1]->vector().Py();
        pi2_pz_2_=hads2[1]->vector().Pz();
        pi2_E_2_=hads2[1]->vector().E();

        pi3_px_2_=hads2[2]->vector().Px();
        pi3_py_2_=hads2[2]->vector().Py();
        pi3_pz_2_=hads2[2]->vector().Pz();
        pi3_E_2_=hads2[2]->vector().E();

        if(tau2->hasSV()) {
          sv_x_2_ = tau2->secondary_vertex().X() - refit_vertex->vx();
          sv_y_2_ = tau2->secondary_vertex().Y() - refit_vertex->vy();
          sv_z_2_ = tau2->secondary_vertex().Z() - refit_vertex->vz();

          ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > sv_cov = tau2->s_vtx_covariance()+pv_cov;

          svcov00_2_ = sv_cov(0,0);
          svcov01_2_ = sv_cov(0,1);
          svcov02_2_ = sv_cov(0,2);
          svcov10_2_ = sv_cov(1,0);
          svcov11_2_ = sv_cov(1,1);
          svcov12_2_ = sv_cov(1,2);
          svcov20_2_ = sv_cov(2,0);
          svcov21_2_ = sv_cov(2,1);
          svcov22_2_ = sv_cov(2,2);
        }
      }

      if(gammas2.size()>0 && tau_decay_mode_2_<2) {
        pi0_px_2_=pi0_tau2->vector().Px();
        pi0_py_2_=pi0_tau2->vector().Py();
        pi0_pz_2_=pi0_tau2->vector().Pz();
        pi0_E_2_=pi0_tau2->vector().E();

        gam1_px_2_=gammas2[0]->vector().Px();
        gam1_py_2_=gammas2[0]->vector().Py();
        gam1_pz_2_=gammas2[0]->vector().Pz();
        gam1_E_2_=gammas2[0]->vector().E();
        if(gammas2.size()>1) gam2_px_2_= gammas2[1]->vector().Px();
        if(gammas2.size()>1) gam2_py_2_= gammas2[1]->vector().Py();
        if(gammas2.size()>1) gam2_pz_2_= gammas2[1]->vector().Pz();
        if(gammas2.size()>1) gam2_E_2_= gammas2[1]->vector().E();
        if(gammas2.size()>2) gam3_px_2_= gammas2[2]->vector().Px();
        if(gammas2.size()>2) gam3_py_2_= gammas2[2]->vector().Py();
        if(gammas2.size()>2) gam3_pz_2_= gammas2[2]->vector().Pz();
        if(gammas2.size()>2) gam3_E_2_= gammas2[2]->vector().E();
        if(gammas2.size()>3) gam4_px_2_= gammas2[3]->vector().Px();
        if(gammas2.size()>3) gam4_py_2_= gammas2[3]->vector().Py();
        if(gammas2.size()>3) gam4_pz_2_= gammas2[3]->vector().Pz();
        if(gammas2.size()>3) gam4_E_2_= gammas2[3]->vector().E();
      }
 
      TLorentzVector lvec1;
      TLorentzVector lvec2;
      TLorentzVector lvec3;
      TLorentzVector lvec4;
      TLorentzVector pvtosv1;

      auto primary_vtx = refit_vertex;

      pvtosv1.SetXYZT(
              muon1->vx() - refit_vertex->vx(),
              muon1->vy() - refit_vertex->vy(),
              muon1->vz() - refit_vertex->vz(),
              0.);

      std::pair<TVector3,double> ipandsig_1 = IPAndSignificanceMuon(muon1, refit_vertex);
      std::pair<TVector3,double> ipandsig_2 = IPAndSignificance(tau2, refit_vertex, pfcands);
      ip_mag_1_ = ipandsig_1.first.Mag();
      ip_sig_1_ = ipandsig_1.second;
      ip_mag_2_ = ipandsig_2.first.Mag();
      ip_sig_2_ = ipandsig_2.second;

      TVector3 ip_corrected_1 = ipandsig_1.first;
      TVector3 ip_corrected_2 = ipandsig_2.first;

      std::pair<TVector3,ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> ipandcov_1 = IPAndCovMuon(muon1, refit_vertex);
      std::pair<TVector3,ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> ipandcov_2 = IPAndCov(tau2, refit_vertex, pfcands);
      CovMatrix cov_corrected_1 = ipandcov_1.second;
      CovMatrix cov_corrected_2 = ipandcov_2.second;

      if(!is_data_) {

        ip_sig_1_raw_= ip_sig_1_;
        ip_sig_2_raw_= ip_sig_2_;


        ip_corrected_1 = ipCorrector.correctIp(
                         ipandsig_1.first,
                         ipgen1,
                         fabs(lep1->eta())
        );

        ip_corrected_2 = ipCorrector.correctIp(
                         ipandsig_2.first,
                         ipgen2,
                         fabs(lep2->eta())
        );


        // Correct IP covariance matrix
        cov_corrected_1 = ipCorrector.correctIpCov(
                    ipandcov_1.second, 
                    fabs(lep1->eta())
        );

        cov_corrected_2 = ipCorrector.correctIpCov(
                    ipandcov_2.second,
                    fabs(lep2->eta())
        );

        ROOT::Math::SVector<double, 3> ip_svec_1;
        ip_svec_1(0) = ip_corrected_1.X();
        ip_svec_1(1) = ip_corrected_1.Y();
        ip_svec_1(2) = ip_corrected_1.Z();
        ip_svec_1 = ip_svec_1.Unit();
        ip_sig_1_ = ip_corrected_1.Mag()/sqrt(ROOT::Math::Dot( ip_svec_1, cov_corrected_1 * ip_svec_1));

        ROOT::Math::SVector<double, 3> ip_svec_2;
        ip_svec_2(0) = ip_corrected_2.X();
        ip_svec_2(1) = ip_corrected_2.Y();
        ip_svec_2(2) = ip_corrected_2.Z();
        ip_svec_2 = ip_svec_2.Unit();
        ip_sig_2_ = ip_corrected_2.Mag()/sqrt(ROOT::Math::Dot( ip_svec_2, cov_corrected_2 * ip_svec_2));

        ip_sig_1_down_= ip_sig_1_ -0.2*(ip_sig_1_-ip_sig_1_raw_);
        ip_sig_1_up_  = ip_sig_1_ +0.2*(ip_sig_1_-ip_sig_1_raw_);
        ip_sig_2_down_= ip_sig_2_ -0.2*(ip_sig_2_-ip_sig_2_raw_);
        ip_sig_2_up_  = ip_sig_2_ +0.2*(ip_sig_2_-ip_sig_2_raw_);

      } else {
        ip_sig_1_raw_= ip_sig_1_;
        ip_sig_2_raw_= ip_sig_2_;
        ip_sig_1_up_= ip_sig_1_;
        ip_sig_2_up_= ip_sig_2_;
        ip_sig_1_down_= ip_sig_1_;
        ip_sig_2_down_= ip_sig_2_;
      }

      ip_x_1_ = ip_corrected_1.X();
      ip_y_1_ = ip_corrected_1.Y();
      ip_z_1_ = ip_corrected_1.Z();
      ip_x_2_ = ip_corrected_2.X();
      ip_y_2_ = ip_corrected_2.Y();
      ip_z_2_ = ip_corrected_2.Z();

      ipcov00_1_ = cov_corrected_1(0,0);
      ipcov01_1_ = cov_corrected_1(0,1);
      ipcov02_1_ = cov_corrected_1(0,2);
      ipcov10_1_ = cov_corrected_1(1,0);
      ipcov11_1_ = cov_corrected_1(1,1);
      ipcov12_1_ = cov_corrected_1(1,2);
      ipcov20_1_ = cov_corrected_1(2,0);
      ipcov21_1_ = cov_corrected_1(2,1);
      ipcov22_1_ = cov_corrected_1(2,2);

      ipcov00_2_ = cov_corrected_2(0,0);
      ipcov01_2_ = cov_corrected_2(0,1);
      ipcov02_2_ = cov_corrected_2(0,2);
      ipcov10_2_ = cov_corrected_2(1,0);
      ipcov11_2_ = cov_corrected_2(1,1);
      ipcov12_2_ = cov_corrected_2(1,2);
      ipcov20_2_ = cov_corrected_2(2,0);
      ipcov21_2_ = cov_corrected_2(2,1);
      ipcov22_2_ = cov_corrected_2(2,2);

      // uncorrected IP
      //TVector3 ip1 = (ipandsig_1.first).Unit();
      //TVector3 ip2 = (ipandsig_2.first).Unit();
      
      // corrected IP
      TVector3 ip1 = ip_corrected_1.Unit();
      TVector3 ip2 = ip_corrected_2.Unit();

      lvec1 = TLorentzVector(ip1, 0.);
      lvec3 = ConvertToLorentz(muon1->vector());


      if(tau_decay_mode_2_==0) {

        TLorentzVector pvtosv2(
                tau2->svx() - primary_vtx->vx(),
                tau2->svy() - primary_vtx->vy(),
                tau2->svz() - primary_vtx->vz(),
                0.);
        lvec4 = ConvertToLorentz(tau2->vector()); 

        //TVector3 ip2 = (pvtosv2.Vect() - pvtosv2.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
        lvec2 = TLorentzVector(ip2, 0.);

        aco_angle_6_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
        aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

        if(event_ % 2) {
          aco_angle_rand_ = aco_angle_6_;
          aco_sign_rand_ = aco_sign_;
        } else {
          aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
          aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
        }

        alpha1_1_ = AlphaAngle(lvec3.Vect(), ip1);
        alpha1_2_ = AlphaAngle(lvec4.Vect(), ip2);
      }

      if(tau_decay_mode_2_==1) {
               
        ic::Candidate *pi;
        ic::Candidate *pi0;

        pi = pi_tau2;
        pi0 = pi0_tau2;

        lvec2 = ConvertToLorentz(pi0->vector()); //pi zero from rho
        lvec4 = ConvertToLorentz(pi->vector()); //pi charge from rho

        double cp_sign_ = YRho(std::vector<Candidate*>({pi, pi0}),TVector3());
        y_1_1_ = YRho(std::vector<Candidate*>({pi, pi0}),TVector3());
        TLorentzVector pvtosv2(
                tau2->svx() - primary_vtx->vx(),
                tau2->svy() - primary_vtx->vy(),
                tau2->svz() - primary_vtx->vz(),
                0.);
        
        TLorentzVector lvec4_2 = ConvertToLorentz(tau2->vector());

        //TVector3 ip2 = (pvtosv2.Vect() - pvtosv2.Vect().Dot(lvec4_2.Vect().Unit())*lvec4_2.Vect().Unit()).Unit();
        TLorentzVector lvec2_2 = TLorentzVector(ip2, 0.);

        aco_angle_5_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
        aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);
        if (cp_sign_<0) {
          if (aco_angle_5_<M_PI)  aco_angle_5_ = aco_angle_5_+M_PI;
          else                    aco_angle_5_ = aco_angle_5_-M_PI;
        }
        aco_angle_6_ = IPAcoAngle(lvec1, lvec2_2, lvec3, lvec4_2,false);
        if (cp_sign_<0) {
          if (aco_angle_6_<M_PI)  aco_angle_6_ = aco_angle_6_+M_PI;
          else                    aco_angle_6_ = aco_angle_6_-M_PI;
        }

        if(event_ % 2) {
          aco_angle_rand_ = aco_angle_5_;
          aco_sign_rand_ = aco_sign_;
        } else {
          aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
          aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
          if (cp_sign_<0) {
            if (aco_angle_rand_<M_PI) aco_angle_rand_ = aco_angle_rand_+M_PI;
            else                      aco_angle_rand_ = aco_angle_rand_-M_PI;
          }
        }

        alpha1_1_ = AlphaAngle(lvec3.Vect(), ip1);
        alpha1_2_ = AlphaAngle(lvec4_2.Vect(), ip2);
        alpha2_2_ = AlphaAngleRho(lvec4.Vect(), lvec2.Vect());
      }

      else if(tau_decay_mode_2_>=10){

        std::vector<ic::PFCandidate*> a1_daughters;
        std::pair<ic::Candidate*, ic::Candidate*> rho_daughters;
        a1_daughters  = GetA1(tau2, pfcands).first;

        if (a1_daughters.size()>2){

            if(tau2->hasSV()) {

              TVector3 svminuspv_2(
                    tau2->secondary_vertex().X() - primary_vtx->vx(),
                    tau2->secondary_vertex().Y() - primary_vtx->vy(),
                    tau2->secondary_vertex().Z() - primary_vtx->vz()
              );

              sv_x_2_=svminuspv_2.X();
              sv_y_2_=svminuspv_2.Y();
              sv_z_2_=svminuspv_2.Z();
              sv_mag_2_=svminuspv_2.Mag();
              sv_dphi_2_=svminuspv_2.DeltaPhi (ConvertToTVector3(tau2->vector()));
              sv_deta_2_=svminuspv_2.Eta()-tau2->eta();
              hasSV_2_=true;
            }

            lvec2 = ConvertToLorentz(a1_daughters[0]->vector()); //pi zero from rho
            lvec4 = ConvertToLorentz(a1_daughters[1]->vector()); //pi charge from rho
            aco_angle_5_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
            aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

            double cp_sign_ = YRho(std::vector<Candidate*>({a1_daughters[0], a1_daughters[1]}),TVector3());

            if (cp_sign_<0) {
              if (aco_angle_5_<M_PI)  aco_angle_5_ = aco_angle_5_+M_PI;
              else                    aco_angle_5_ = aco_angle_5_-M_PI;
            }

            if(event_ % 2) {
              aco_angle_rand_ = aco_angle_5_;
              aco_sign_rand_ = aco_sign_;
            } else {
              aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
              aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
              if (cp_sign_<0) {
                if (aco_angle_rand_<M_PI) aco_angle_rand_ = aco_angle_rand_+M_PI;
                else                      aco_angle_rand_ = aco_angle_rand_-M_PI;
              }
            }

          }
          alpha1_1_ = AlphaAngle(lvec3.Vect(), ip1);
          alpha2_2_ = AlphaAngleRho(lvec4.Vect(), lvec2.Vect());
      }
     
      // save IP for sync 
      ipx_1_ = ip1.X();
      ipy_1_ = ip1.Y();
      ipz_1_ = ip1.Z();
      ipx_2_ = ip2.X();
      ipy_2_ = ip2.Y();
      ipz_2_ = ip2.Z();

      // gen IP
      gen_ipx_1_ = ipgen1.X();
      gen_ipy_1_ = ipgen1.Y();
      gen_ipz_1_ = ipgen1.Z();
      gen_ipx_2_ = ipgen2.X();
      gen_ipy_2_ = ipgen2.Y();
      gen_ipz_2_ = ipgen2.Z();

      // primary vertex
      pvx_ = primary_vtx->vx();
      pvy_ = primary_vtx->vy();
      pvz_ = primary_vtx->vz();
      
    }
    else if (channel_ == channel::et && event->ExistsInTree("pfCandidates")) {
      Electron const* ele1 = dynamic_cast<Electron const*>(lep1);
      Tau const* tau2 = dynamic_cast<Tau const*>(lep2);

      std::vector<ic::Vertex*> & vertex_vec = event->GetPtrVec<ic::Vertex>("vertices");
      //std::vector<ic::Vertex*> & refit_vertex_vec = event->GetPtrVec<ic::Vertex>("refittedVertices");
      std::vector<ic::Vertex*> & refit_vertex_vec = event->GetPtrVec<ic::Vertex>("refittedVerticesBS");
      ic::Vertex* refit_vertex = vertex_vec[0];
      for(auto v : refit_vertex_vec) {
        if(v->id() == ele1->id()+tau2->id()) {
          refit_vertex = v;
          use_refitted_vertex_ = true;
        }
      }

      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pv_cov = refit_vertex->covariance();

      lead_pt_2_ =  tau2->lead_pt();
      raw_pt_1_ = ele1->ecalTrkEnergyPreCorr()/ele1->energy()*pt_1_.var_float;

      tau_mva_decay_mode_2_ = tau2->HasTauID(MVADM2017) ? tau2->GetTauID(MVADM2017) : 0.0;

      std::vector<ic::PFCandidate*> pfcands =  event->GetPtrVec<ic::PFCandidate>("pfCandidates");
      std::pair<ic::Candidate*, ic::Candidate*> rho = GetRho(tau2, pfcands);
      ic::Candidate *pi_tau2 = rho.first;
      ic::Candidate *pi0_tau2 = rho.second;

      std::vector<ic::PFCandidate*> hads2 = GetTauHads(tau2, pfcands,0);
      std::vector<ic::PFCandidate*> gammas2 = GetTauGammas(tau2, pfcands, 1., 0);
      if(hads2.size()>2) hads2 = GetA1 (tau2, pfcands).first;

      pi_px_1_=lep1->vector().Px();
      pi_py_1_=lep1->vector().Py();
      pi_pz_1_=lep1->vector().Pz();
      pi_E_1_=lep1->vector().E();

      if(tau_decay_mode_2_<3.) {
        pi_px_2_=pi_tau2->vector().Px();
        pi_py_2_=pi_tau2->vector().Py();
        pi_pz_2_=pi_tau2->vector().Pz();
        pi_E_2_=pi_tau2->vector().E();
      } else if(hads2.size()>2) {
        pi_px_2_=hads2[0]->vector().Px();
        pi_py_2_=hads2[0]->vector().Py();
        pi_pz_2_=hads2[0]->vector().Pz();
        pi_E_2_=hads2[0]->vector().E();

        pi2_px_2_=hads2[1]->vector().Px();
        pi2_py_2_=hads2[1]->vector().Py();
        pi2_pz_2_=hads2[1]->vector().Pz();
        pi2_E_2_=hads2[1]->vector().E();

        pi3_px_2_=hads2[2]->vector().Px();
        pi3_py_2_=hads2[2]->vector().Py();
        pi3_pz_2_=hads2[2]->vector().Pz();
        pi3_E_2_=hads2[2]->vector().E();

        if(tau2->hasSV()) {
          sv_x_2_ = tau2->secondary_vertex().X() - refit_vertex->vx();
          sv_y_2_ = tau2->secondary_vertex().Y() - refit_vertex->vy();
          sv_z_2_ = tau2->secondary_vertex().Z() - refit_vertex->vz();

          ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > sv_cov = tau2->s_vtx_covariance()+pv_cov;

          svcov00_2_ = sv_cov(0,0);
          svcov01_2_ = sv_cov(0,1);
          svcov02_2_ = sv_cov(0,2);
          svcov10_2_ = sv_cov(1,0);
          svcov11_2_ = sv_cov(1,1);
          svcov12_2_ = sv_cov(1,2);
          svcov20_2_ = sv_cov(2,0);
          svcov21_2_ = sv_cov(2,1);
          svcov22_2_ = sv_cov(2,2);
        }
      }

      if(gammas2.size()>0 && tau_decay_mode_2_<2) {
        pi0_px_2_=pi0_tau2->vector().Px();
        pi0_py_2_=pi0_tau2->vector().Py();
        pi0_pz_2_=pi0_tau2->vector().Pz();
        pi0_E_2_=pi0_tau2->vector().E();
      
        gam1_px_2_=gammas2[0]->vector().Px();
        gam1_py_2_=gammas2[0]->vector().Py();
        gam1_pz_2_=gammas2[0]->vector().Pz();
        gam1_E_2_=gammas2[0]->vector().E(); 
        if(gammas2.size()>1) gam2_px_2_= gammas2[1]->vector().Px();
        if(gammas2.size()>1) gam2_py_2_= gammas2[1]->vector().Py();
        if(gammas2.size()>1) gam2_pz_2_= gammas2[1]->vector().Pz();
        if(gammas2.size()>1) gam2_E_2_= gammas2[1]->vector().E();
        if(gammas2.size()>2) gam3_px_2_= gammas2[2]->vector().Px();
        if(gammas2.size()>2) gam3_py_2_= gammas2[2]->vector().Py();
        if(gammas2.size()>2) gam3_pz_2_= gammas2[2]->vector().Pz();
        if(gammas2.size()>2) gam3_E_2_= gammas2[2]->vector().E();
        if(gammas2.size()>3) gam4_px_2_= gammas2[3]->vector().Px();
        if(gammas2.size()>3) gam4_py_2_= gammas2[3]->vector().Py();
        if(gammas2.size()>3) gam4_pz_2_= gammas2[3]->vector().Pz();
        if(gammas2.size()>3) gam4_E_2_= gammas2[3]->vector().E();
      }
 
      TLorentzVector lvec1;
      TLorentzVector lvec2;
      TLorentzVector lvec3;
      TLorentzVector lvec4;
      TLorentzVector pvtosv;

      auto primary_vtx = refit_vertex;

      std::pair<TVector3,double> ipandsig_1 = IPAndSignificanceElectron(ele1, refit_vertex);
      std::pair<TVector3,double> ipandsig_2 = IPAndSignificance(tau2, refit_vertex, pfcands);
      ip_mag_1_ = ipandsig_1.first.Mag();
      ip_sig_1_ = ipandsig_1.second;
      ip_mag_2_ = ipandsig_2.first.Mag();
      ip_sig_2_ = ipandsig_2.second;

      TVector3 ip_corrected_1 = ipandsig_1.first;
      TVector3 ip_corrected_2 = ipandsig_2.first;


      ipx_uncorr_1_ = ipandsig_1.first.X();
      ipy_uncorr_1_ = ipandsig_1.first.Y();
      ipz_uncorr_1_ = ipandsig_1.first.Z();
      ipx_uncorr_2_ = ipandsig_2.first.X();
      ipy_uncorr_2_ = ipandsig_2.first.Y();
      ipz_uncorr_2_ = ipandsig_2.first.Z();
      ip_sig_uncorr_1_ = ipandsig_1.second;
      ip_sig_uncorr_2_ = ipandsig_1.second;

      std::pair<TVector3,ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> ipandcov_1 = IPAndCovElec(ele1, refit_vertex);
      std::pair<TVector3,ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> ipandcov_2 = IPAndCov(tau2, refit_vertex, pfcands);
      CovMatrix cov_corrected_1 = ipandcov_1.second;
      CovMatrix cov_corrected_2 = ipandcov_2.second;

      if(!is_data_) {

        ip_sig_1_raw_= ip_sig_1_;
        ip_sig_2_raw_= ip_sig_2_;

        ip_corrected_1 = ipCorrectorEle.correctIp(
                         ipandsig_1.first,
                         ipgen1,
                         fabs(lep1->eta())
        );

        ip_corrected_2 = ipCorrector.correctIp(
                         ipandsig_2.first,
                         ipgen2,
                         fabs(lep2->eta())
        );

        // Correct IP covariance matrix
        cov_corrected_1 = ipCorrectorEle.correctIpCov(
                    ipandcov_1.second,
                    fabs(lep1->eta())
        );

        cov_corrected_2 = ipCorrector.correctIpCov(
                    ipandcov_2.second,
                    fabs(lep2->eta())
        );

        ROOT::Math::SVector<double, 3> ip_svec_1;
        ip_svec_1(0) = ip_corrected_1.X();
        ip_svec_1(1) = ip_corrected_1.Y();
        ip_svec_1(2) = ip_corrected_1.Z();
        ip_svec_1 = ip_svec_1.Unit();
        ip_sig_1_ = ip_corrected_1.Mag()/sqrt(ROOT::Math::Dot( ip_svec_1, cov_corrected_1 * ip_svec_1));

        ROOT::Math::SVector<double, 3> ip_svec_2;
        ip_svec_2(0) = ip_corrected_2.X();
        ip_svec_2(1) = ip_corrected_2.Y();
        ip_svec_2(2) = ip_corrected_2.Z();
        ip_svec_2 = ip_svec_2.Unit();
        ip_sig_2_ = ip_corrected_2.Mag()/sqrt(ROOT::Math::Dot( ip_svec_2, cov_corrected_2 * ip_svec_2));

        ip_sig_1_down_= ip_sig_1_ -0.2*(ip_sig_1_-ip_sig_1_raw_);
        ip_sig_1_up_  = ip_sig_1_ +0.2*(ip_sig_1_-ip_sig_1_raw_);
        ip_sig_2_down_= ip_sig_2_ -0.2*(ip_sig_2_-ip_sig_2_raw_);
        ip_sig_2_up_  = ip_sig_2_ +0.2*(ip_sig_2_-ip_sig_2_raw_);

      } else {
        ip_sig_1_raw_= ip_sig_1_;
        ip_sig_2_raw_= ip_sig_2_;
        ip_sig_1_up_= ip_sig_1_;
        ip_sig_2_up_= ip_sig_2_;
        ip_sig_1_down_= ip_sig_1_;
        ip_sig_2_down_= ip_sig_2_;
      }

      ip_x_1_ = ip_corrected_1.X();
      ip_y_1_ = ip_corrected_1.Y();
      ip_z_1_ = ip_corrected_1.Z();
      ip_x_2_ = ip_corrected_2.X();
      ip_y_2_ = ip_corrected_2.Y();
      ip_z_2_ = ip_corrected_2.Z();

      ipcov00_1_ = cov_corrected_1(0,0); 
      ipcov01_1_ = cov_corrected_1(0,1); 
      ipcov02_1_ = cov_corrected_1(0,2); 
      ipcov10_1_ = cov_corrected_1(1,0); 
      ipcov11_1_ = cov_corrected_1(1,1); 
      ipcov12_1_ = cov_corrected_1(1,2); 
      ipcov20_1_ = cov_corrected_1(2,0); 
      ipcov21_1_ = cov_corrected_1(2,1); 
      ipcov22_1_ = cov_corrected_1(2,2); 

      ipcov00_2_ = cov_corrected_2(0,0);
      ipcov01_2_ = cov_corrected_2(0,1);
      ipcov02_2_ = cov_corrected_2(0,2);
      ipcov10_2_ = cov_corrected_2(1,0);
      ipcov11_2_ = cov_corrected_2(1,1);
      ipcov12_2_ = cov_corrected_2(1,2);
      ipcov20_2_ = cov_corrected_2(2,0);
      ipcov21_2_ = cov_corrected_2(2,1);
      ipcov22_2_ = cov_corrected_2(2,2);

      // corrected IP
      TVector3 ip1 = ip_corrected_1.Unit();
      TVector3 ip2 = ip_corrected_2.Unit();

      lvec1 = TLorentzVector(ip1, 0.);
      lvec3 = ConvertToLorentz(ele1->vector());


      if(tau_decay_mode_2_==0) {

        TLorentzVector pvtosv2(
                tau2->svx() - primary_vtx->vx(),
                tau2->svy() - primary_vtx->vy(),
                tau2->svz() - primary_vtx->vz(),
                0.);
        lvec4 = ConvertToLorentz(tau2->vector());

        //TVector3 ip2 = (pvtosv2.Vect() - pvtosv2.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
        lvec2 = TLorentzVector(ip2, 0.);

        aco_angle_6_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
        aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

        if(event_ % 2) {
          aco_angle_rand_ = aco_angle_6_;
          aco_sign_rand_ = aco_sign_;
        } else {
          aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
          aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
        }

        alpha1_1_ = AlphaAngle(lvec3.Vect(), ip1);
        alpha1_2_ = AlphaAngle(lvec4.Vect(), ip2);
      }

      if(tau_decay_mode_2_==1) {
               
        ic::Candidate *pi;
        ic::Candidate *pi0;
        
        pi = pi_tau2;
        pi0 = pi0_tau2;
        
        lvec2 = ConvertToLorentz(pi0->vector()); //pi zero from rho
        lvec4 = ConvertToLorentz(pi->vector()); //pi charge from rho
        
        double cp_sign_ = YRho(std::vector<Candidate*>({pi, pi0}),TVector3());
        y_1_1_ = YRho(std::vector<Candidate*>({pi, pi0}),TVector3());
        TLorentzVector pvtosv2(
                tau2->svx() - primary_vtx->vx(),
                tau2->svy() - primary_vtx->vy(),
                tau2->svz() - primary_vtx->vz(),
                0.);
        
        TLorentzVector lvec4_2 = ConvertToLorentz(tau2->vector());
        
        //TVector3 ip2 = (pvtosv2.Vect() - pvtosv2.Vect().Dot(lvec4_2.Vect().Unit())*lvec4_2.Vect().Unit()).Unit();
        TLorentzVector lvec2_2 = TLorentzVector(ip2, 0.);
        
        aco_angle_5_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
        aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);
        if (cp_sign_<0) {
          if (aco_angle_5_<M_PI)  aco_angle_5_ = aco_angle_5_+M_PI;
          else                    aco_angle_5_ = aco_angle_5_-M_PI;
        }
        aco_angle_6_ = IPAcoAngle(lvec1, lvec2_2, lvec3, lvec4_2,false);
        if (cp_sign_<0) {
          if (aco_angle_6_<M_PI)  aco_angle_6_ = aco_angle_6_+M_PI;
          else                    aco_angle_6_ = aco_angle_6_-M_PI;
        }
        
        alpha1_1_ = AlphaAngle(lvec3.Vect(), ip1);
        alpha1_2_ = AlphaAngle(lvec4_2.Vect(), ip2);
        alpha2_2_ = AlphaAngleRho(lvec4.Vect(), lvec2.Vect());
      }

      else if(tau_decay_mode_2_>=10){

        std::vector<ic::PFCandidate*> a1_daughters;
        std::pair<ic::Candidate*, ic::Candidate*> rho_daughters;
        a1_daughters  = GetA1(tau2, pfcands).first;

        if (a1_daughters.size()>2){

            lvec2 = ConvertToLorentz(a1_daughters[0]->vector()); //pi zero from rho
            lvec4 = ConvertToLorentz(a1_daughters[1]->vector()); //pi charge from rho
            aco_angle_5_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
            aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

            double cp_sign_ = YRho(std::vector<Candidate*>({a1_daughters[0], a1_daughters[1]}),TVector3());

            if (cp_sign_<0) {
              if (aco_angle_5_<M_PI)  aco_angle_5_ = aco_angle_5_+M_PI;
              else                    aco_angle_5_ = aco_angle_5_-M_PI;
            }

            if(event_ % 2) {
              aco_angle_rand_ = aco_angle_5_;
              aco_sign_rand_ = aco_sign_;
            } else {
              aco_angle_rand_ = IPAcoAngle(lvec2, lvec1, lvec4, lvec3,false);
              aco_sign_rand_ = IPAcoSign(lvec2, lvec1, lvec4, lvec3,false);
              if (cp_sign_<0) {
                if (aco_angle_rand_<M_PI) aco_angle_rand_ = aco_angle_rand_+M_PI;
                else                      aco_angle_rand_ = aco_angle_rand_-M_PI;
              }
            }

          }
          alpha1_1_ = AlphaAngle(lvec3.Vect(), ip1);
          alpha2_2_ = AlphaAngleRho(lvec4.Vect(), lvec2.Vect());
      }
    
      // save IP for sync 
      ipx_1_ = ip_corrected_1.X();
      ipy_1_ = ip_corrected_1.Y();
      ipz_1_ = ip_corrected_1.Z();
      ipx_2_ = ip_corrected_2.X();
      ipy_2_ = ip_corrected_2.Y();
      ipz_2_ = ip_corrected_2.Z();

      // gen IP
      gen_ipx_1_ = ipgen1.X();
      gen_ipy_1_ = ipgen1.Y();
      gen_ipz_1_ = ipgen1.Z();
      gen_ipx_2_ = ipgen2.X();
      gen_ipy_2_ = ipgen2.Y();
      gen_ipz_2_ = ipgen2.Z();

      // primary vertex
      pvx_ = primary_vtx->vx();
      pvy_ = primary_vtx->vy();
      pvz_ = primary_vtx->vz();

    }
      
    // signal background event classification
    IC_BDT_max_score_ = event->Exists("IC_BDT_max_score") ? event->Get<float>("IC_BDT_max_score") : -999.0;
    IC_BDT_max_index_ = event->Exists("IC_BDT_max_index") ? event->Get<int>("IC_BDT_max_index") : -999;

    // smear JERphi things
    //gRandom->SetSeed((int)((lep1->eta()+2.5)*100000 + (lep1->phi()+4)*1000));
    //double smear_jphi_1 = 0.;
    //double smear_jphi_2 = 0.;
    //shifted_dphi_jtt_smear_ = -9999.;
    //dphi_jtt_smear_ = -9999.;
    //sjdphi_smear_ = -9999.;

    //// set dphi_jtt_smear_ equal to dphi_jtt_ for data, shift for MC.
    //if (n_lowpt_jets_ >= 1) {
    //  dphi_jtt_smear_ =  ROOT::Math::VectorUtil::DeltaPhi(lowpt_jets[0]->vector(), ditau->vector());
    //  if(!is_data_) {
    //    // do jer shifts for dphi here
    //    if (lowpt_jets[0]->vector().pt() > 60.) {
    //      if (fabs(lowpt_jets[0]->vector().eta()) > 0. && fabs(lowpt_jets[0]->vector().eta()) < 1.1) {
    //        smear_jphi_1 = fns_["smear_func_ptlow60_etalow0p0"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[0]->vector().eta()) > 1.1 && fabs(lowpt_jets[0]->vector().eta()) < 2.1){
    //        smear_jphi_1 = fns_["smear_func_ptlow60_etalow1p1"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[0]->vector().eta()) > 2.1 && fabs(lowpt_jets[0]->vector().eta()) < 3.0){
    //        smear_jphi_1 = fns_["smear_func_ptlow60_etalow2p1"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[0]->vector().eta()) > 3. && fabs(lowpt_jets[0]->vector().eta()) < 5.){
    //        smear_jphi_1 = fns_["smear_func_ptlow60_etalow3p0"]->GetRandom();
    //      }
    //    }
    //    else if (lowpt_jets[0]->vector().pt() > 30.) {
    //      if (fabs(lowpt_jets[0]->vector().eta()) > 0. && fabs(lowpt_jets[0]->vector().eta()) < 1.1) {
    //        smear_jphi_1 = fns_["smear_func_ptlow30_etalow0p0"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[0]->vector().eta()) > 1.1 && fabs(lowpt_jets[0]->vector().eta()) < 2.1){
    //        smear_jphi_1 = fns_["smear_func_ptlow30_etalow1p1"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[0]->vector().eta()) > 2.1 && fabs(lowpt_jets[0]->vector().eta()) < 3.0){
    //        smear_jphi_1 = fns_["smear_func_ptlow30_etalow2p1"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[0]->vector().eta()) > 3. && fabs(lowpt_jets[0]->vector().eta()) < 5.){
    //        smear_jphi_1 = fns_["smear_func_ptlow30_etalow3p0"]->GetRandom();
    //      }
    //    }
    //  }
    //  // shifted dphi_jtt_ is already wrapped correctly so can just subtract
    //  // smear_phi_1 i think
    //  shifted_dphi_jtt_smear_ = boundPhi(shifted_dphi_jtt_ - smear_jphi_1);
    //  dphi_jtt_smear_ = (shifted_dphi_jtt_smear_ > 0.) * (shifted_dphi_jtt_smear_ - M_PI)
    //    + (shifted_dphi_jtt_smear_ < 0.) * (shifted_dphi_jtt_smear_ + M_PI);
    //  smear_jet_delta_ = smear_jphi_1;
    //}

    //// if two jets smear lowpt_jets[1] for dphi_jj variables
    //if (n_lowpt_jets_ >= 2) {
    //  if(lowpt_jets[0]->eta() > lowpt_jets[1]->eta()){
    //    sjdphi_smear_ =  ROOT::Math::VectorUtil::DeltaPhi(lowpt_jets[0]->vector(), lowpt_jets[1]->vector());
    //  }
    //  else{
    //    sjdphi_smear_ =  ROOT::Math::VectorUtil::DeltaPhi(lowpt_jets[1]->vector(), lowpt_jets[0]->vector());
    //  }
    //  if(!is_data_) {
    //    // do jer shifts for lowpt_jets[1] here, jet[0] smear value is already
    //    // obtained above
    //    if (lowpt_jets[1]->vector().pt() > 60.) {
    //      if (fabs(lowpt_jets[1]->vector().eta()) > 0. && fabs(lowpt_jets[1]->vector().eta()) < 1.1) {
    //        smear_jphi_2 = fns_["smear_func_ptlow60_etalow0p0"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[1]->vector().eta()) > 1.1 && fabs(lowpt_jets[1]->vector().eta()) < 2.1){
    //        smear_jphi_2 = fns_["smear_func_ptlow60_etalow1p1"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[1]->vector().eta()) > 2.1 && fabs(lowpt_jets[1]->vector().eta()) < 3.0){
    //        smear_jphi_2 = fns_["smear_func_ptlow60_etalow2p1"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[1]->vector().eta()) > 3. && fabs(lowpt_jets[1]->vector().eta()) < 5.){
    //        smear_jphi_2 = fns_["smear_func_ptlow60_etalow3p0"]->GetRandom();
    //      }
    //    }
    //    else if (lowpt_jets[1]->vector().pt() > 30.) {
    //      if (fabs(lowpt_jets[1]->vector().eta()) > 0. && fabs(lowpt_jets[1]->vector().eta()) < 1.1) {
    //        smear_jphi_2 = fns_["smear_func_ptlow30_etalow0p0"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[1]->vector().eta()) > 1.1 && fabs(lowpt_jets[1]->vector().eta()) < 2.1){
    //        smear_jphi_2 = fns_["smear_func_ptlow30_etalow1p1"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[1]->vector().eta()) > 2.1 && fabs(lowpt_jets[1]->vector().eta()) < 3.0){
    //        smear_jphi_2 = fns_["smear_func_ptlow30_etalow2p1"]->GetRandom();
    //      }
    //      else if (fabs(lowpt_jets[1]->vector().eta()) > 3. && fabs(lowpt_jets[1]->vector().eta()) < 5.){
    //        smear_jphi_2 = fns_["smear_func_ptlow30_etalow3p0"]->GetRandom();
    //      }
    //    }
    //  }
    //  // calculate smeared sjdphi
    //  if(lowpt_jets[0]->eta() > lowpt_jets[1]->eta()){
    //    sjdphi_smear_ =  boundPhi(sjdphi_ + smear_jphi_2 - smear_jphi_1);
    //  }
    //  else{
    //    sjdphi_smear_ =  boundPhi(sjdphi_ + smear_jphi_1 - smear_jphi_2);
    //  }
    //}

    // add 
    // add bools for etau cross trigger SF measurment
    if(channel_ == channel::mt) {
      trg_etau_matched_ = false;
      std::vector<ic::L1TObject*> l1taus = event->GetPtrVec<ic::L1TObject>("L1Taus");
      std::vector<Candidate const *> match_tau;
      match_tau.push_back(lep2);
      std::pair<bool, unsigned> match_etau;
      if(strategy_ == strategy::cpdecays17 || (is_data_ && strategy_ == strategy::cpdecays18 && run_ < 317509)){
        std::string trg_etau_label = "triggerObjectsIsoMu20Tau27";
        std::string filter_etau_Data_MC = "hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched";
        std::string filter_etau_Embed = "hltL1sMu18erTau24erIorMu20erTau24er";

        std::vector<TriggerObject *> const& objs_etau = event->GetPtrVec<TriggerObject>(trg_etau_label);

        if (is_embedded_)
            match_etau = IsFilterMatchedWithIndex(lep2, objs_etau, filter_etau_Embed, 0.5);
        else
            match_etau = IsFilterMatchedWithIndex(lep2, objs_etau, filter_etau_Data_MC, 0.5);

        //HLT pt cut - only for MC and for embedding we only have L1 
        if (match_etau.first && !is_embedded_)  trg_etau_matched_ = objs_etau[match_etau.second]->vector().Pt()>30;

        //L1 cut
        std::vector<ic::L1TObject*> passed_l1_taus_etau;
        for(unsigned ta=0; ta<l1taus.size(); ++ta){
            if(l1taus[ta]->isolation()!=0 && l1taus[ta]->vector().Pt() >= 26.) passed_l1_taus_etau.push_back(l1taus[ta]);
        }
        bool match_l1_part_etau = (MatchByDR(match_tau,passed_l1_taus_etau,0.5,true,true)).size() == 1;
        trg_etau_matched_ = trg_etau_matched_ && match_l1_part_etau;

      }
      else if(strategy_ == strategy::cpdecays18){
        std::string trg_etau_label = "triggerObjectsIsoMu20TauHPS27";
        std::string filter_etau_Data_MC = "hltHpsSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched";
        std::string filter_etau_Embed = "hltL1sBigORMu18erTauXXer2p1";
        //if(run_< 317509 && is_data_) {
        //  std::string trg_etau_label = "triggerObjectsIsoMu20Tau27";
        //  std::string filter_etau_Data_MC = "hltSelectedPFTau27LooseChargedIsolationAgainstMuonL1HLTMatched";
        //} 


        std::vector<TriggerObject *> const& objs_etau = event->GetPtrVec<TriggerObject>(trg_etau_label);

        if (is_embedded_)
            match_etau = IsFilterMatchedWithIndex(lep2, objs_etau, filter_etau_Embed, 0.5);
        else
            match_etau = IsFilterMatchedWithIndex(lep2, objs_etau, filter_etau_Data_MC, 0.5);

        //HLT pt cut
        if (match_etau.first && !is_embedded_)  trg_etau_matched_ = objs_etau[match_etau.second]->vector().Pt()>30;
        else if(is_embedded_)  trg_etau_matched_ = match_etau.first;       

        //L1 cut
        std::vector<ic::L1TObject*> passed_l1_taus_etau;
        for(unsigned ta=0; ta<l1taus.size(); ++ta){
            if(l1taus[ta]->isolation()!=0 && l1taus[ta]->vector().Pt() >= 26.) passed_l1_taus_etau.push_back(l1taus[ta]);
        }
        bool match_l1_part_etau = (MatchByDR(match_tau,passed_l1_taus_etau,0.5,true,true)).size() == 1;
        trg_etau_matched_ = trg_etau_matched_ && match_l1_part_etau;
      }
    }   

    if (ip_sig_1_ != ip_sig_1_) ip_sig_1_=0.;
    if (ip_sig_2_ != ip_sig_2_) ip_sig_2_=0.;

    aco_angle_1_ = aco_angle_1_!=aco_angle_1_ ? -9999. : aco_angle_1_;
    aco_angle_2_ = aco_angle_2_!=aco_angle_2_ ? -9999. : aco_angle_2_;
    aco_angle_3_ = aco_angle_3_!=aco_angle_3_ ? -9999. : aco_angle_3_;
    aco_angle_4_ = aco_angle_4_!=aco_angle_4_ ? -9999. : aco_angle_4_;
    aco_angle_5_ = aco_angle_5_!=aco_angle_5_ ? -9999. : aco_angle_5_;
    aco_angle_6_ = aco_angle_6_!=aco_angle_6_ ? -9999. : aco_angle_6_;
    aco_angle_7_ = aco_angle_7_!=aco_angle_7_ ? -9999. : aco_angle_7_;
 
    if (write_tree_ && fs_) outtree_->Fill();
    if (make_sync_ntuple_) synctree_->Fill();
    if (make_mva_ntuple_) mvatree_->Fill();

    return 0;
  }

  int HTTCategories::PostAnalysis() {
    if(make_sync_ntuple_) {   
      lOFile->cd();
      synctree_->Write();
      lOFile->Close();
    }
    if(make_mva_ntuple_) {
      lOFile->cd();
      mvatree_->Write();
      lOFile->Close();
    }
    return 0;
  }

  void HTTCategories::PrintInfo() {
    ;
  }
}
