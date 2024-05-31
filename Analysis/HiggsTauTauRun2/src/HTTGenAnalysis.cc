#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/interface/HTTGenAnalysis.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/PolarimetricA1.h"

TVector3 GenIP (ic::GenParticle *h, ic::GenParticle *t) {

  TVector3 pvtosv(
           t->vtx().vx() - h->vtx().vx(),
           t->vtx().vy() - h->vtx().vy(),
           t->vtx().vz() - h->vtx().vz()
           );

  TVector3 momenta(
                  t->vector().Px(),
                  t->vector().Py(),
                  t->vector().Pz()
                  );

  double projection=pvtosv*momenta/momenta.Mag2();
  TVector3 gen_ip=pvtosv-momenta*projection;

  return gen_ip;
}

std::vector<ic::GenParticle> FamilyTree (std::vector<ic::GenParticle> &v, ic::GenParticle p, std::vector<ic::GenParticle*> gen_particles, unsigned &outputID){ 
  if(p.daughters().size() == 0){
    unsigned ID = std::fabs(p.pdgid());
    if(ID == 11) outputID = 11;
    else if(ID == 13) outputID = 13;
    if(!(ID==12 || ID==14 || ID==16)){
      v.push_back(p);
    }
  }
  else{
    for(size_t i = 0; i < p.daughters().size(); ++i ){
      ic::GenParticle d = *gen_particles[p.daughters().at(i)];
      FamilyTree(v,d, gen_particles, outputID);
    }
  }
  return v;
}

struct swap_labels{
  bool operator() (std::string a,std::string b) {
    if(a=="t" && b!="t") return false;
    else if(a=="m" &&  (b!="m" && b!="t")) return false;
    else return true;
  }
};

struct PtComparator{
  bool operator() (ic::Candidate a, ic::Candidate b) {
    return (a.vector().Pt() > b.vector().Pt());
  }
};

struct PtComparatorGenPart{
  bool operator() (ic::GenParticle *a, ic::GenParticle *b) {
    return (a->vector().Pt() > b->vector().Pt());
  }
};

  
namespace ic {

  HTTGenAnalysis::HTTGenAnalysis(std::string const& name) : ModuleBase(name), era_(era::data_2018) {
    fs_ = NULL;
    bbtag_eff_ = nullptr;
  }

  HTTGenAnalysis::~HTTGenAnalysis() {
    ;
  }

  int HTTGenAnalysis::PreAnalysis() {
    TFile f("input/mssm_higgspt/higgs_pt_v2_mssm_mode.root");
    mssm_w_ = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
    f.Close();
    std::string mass_str = mssm_mass_;
    if(mssm_w_->function(("h_"+mass_str+"_t_ratio").c_str())){
      fns_["h_t_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("h_"+mass_str+"_t_ratio").c_str())->functor(mssm_w_->argSet("h_pt")));        
      fns_["h_b_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("h_"+mass_str+"_b_ratio").c_str())->functor(mssm_w_->argSet("h_pt")));
      fns_["h_i_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("h_"+mass_str+"_i_ratio").c_str())->functor(mssm_w_->argSet("h_pt")));
      fns_["H_t_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("H_"+mass_str+"_t_ratio").c_str())->functor(mssm_w_->argSet("h_pt")));
      fns_["H_b_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("H_"+mass_str+"_b_ratio").c_str())->functor(mssm_w_->argSet("h_pt")));
      fns_["H_i_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("H_"+mass_str+"_i_ratio").c_str())->functor(mssm_w_->argSet("h_pt")));
      fns_["A_t_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("A_"+mass_str+"_t_ratio").c_str())->functor(mssm_w_->argSet("h_pt")));
      fns_["A_b_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("A_"+mass_str+"_b_ratio").c_str())->functor(mssm_w_->argSet("h_pt")));
      fns_["A_i_ratio"] = std::shared_ptr<RooFunctor>(mssm_w_->function(("A_"+mass_str+"_i_ratio").c_str())->functor(mssm_w_->argSet("h_pt"))); 
    }

    if(make_mva_ntuple_) {

      lOFile = new TFile(mva_output_name_.c_str(), "RECREATE");
      lOFile->cd();
      mvatree_ = new TTree("ntuple", "ntuple");

      mvatree_->Branch("rand", &rand_);

      mvatree_->Branch("wt_cp_sm"       , &wt_cp_sm_);
      mvatree_->Branch("wt_cp_ps"       , &wt_cp_ps_);
      mvatree_->Branch("wt_cp_mm"       , &wt_cp_mm_);

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

      mvatree_->Branch("nu_px_1", &nu_px_1_);
      mvatree_->Branch("nu_py_1", &nu_py_1_);
      mvatree_->Branch("nu_pz_1", &nu_pz_1_);
      mvatree_->Branch("nu_E_1",  &nu_E_1_);
      mvatree_->Branch("nu_px_2", &nu_px_2_);
      mvatree_->Branch("nu_py_2", &nu_py_2_);
      mvatree_->Branch("nu_pz_2", &nu_pz_2_);
      mvatree_->Branch("nu_E_2",  &nu_E_2_);

      mvatree_->Branch("sv_x_1", &sv_x_1_);
      mvatree_->Branch("sv_y_1", &sv_y_1_);
      mvatree_->Branch("sv_z_1", &sv_z_1_);
      mvatree_->Branch("sv_x_2", &sv_x_2_);
      mvatree_->Branch("sv_y_2", &sv_y_2_);
      mvatree_->Branch("sv_z_2", &sv_z_2_);

      mvatree_->Branch("ip_x_1", &ip_x_1_);
      mvatree_->Branch("ip_y_1", &ip_y_1_);
      mvatree_->Branch("ip_z_1", &ip_z_1_);
      mvatree_->Branch("ip_x_2", &ip_x_2_);
      mvatree_->Branch("ip_y_2", &ip_y_2_);
      mvatree_->Branch("ip_z_2", &ip_z_2_);
      
      mvatree_->Branch("tau1_charge", &tau1_charge);
      mvatree_->Branch("tau2_charge", &tau2_charge);

      mvatree_->Branch("dm_1", &tauFlag_1_);
      mvatree_->Branch("dm_2", &tauFlag_2_);

      mvatree_->Branch("metx"         , &metx_         );
      mvatree_->Branch("mety"         , &mety_         );

      mvatree_->Branch("taupos_px", &taupos_px_);
      mvatree_->Branch("taupos_py", &taupos_py_);
      mvatree_->Branch("taupos_pz", &taupos_pz_);
      mvatree_->Branch("taupos_E", &taupos_E_);
      mvatree_->Branch("tauneg_px", &tauneg_px_);
      mvatree_->Branch("tauneg_py", &tauneg_py_);
      mvatree_->Branch("tauneg_pz", &tauneg_pz_);
      mvatree_->Branch("tauneg_E", &tauneg_E_);
      mvatree_->Branch("taupos_polvec_x", &taupos_polvec_x_);
      mvatree_->Branch("taupos_polvec_y", &taupos_polvec_y_);
      mvatree_->Branch("taupos_polvec_z", &taupos_polvec_z_);
      mvatree_->Branch("taupos_polvec_E", &taupos_polvec_E_);
      mvatree_->Branch("tauneg_polvec_x", &tauneg_polvec_x_);
      mvatree_->Branch("tauneg_polvec_y", &tauneg_polvec_y_);
      mvatree_->Branch("tauneg_polvec_z", &tauneg_polvec_z_);
      mvatree_->Branch("tauneg_polvec_E", &tauneg_polvec_E_);

      mvatree_->Branch("reco_sv_x_1", &reco_sv_x_1_);
      mvatree_->Branch("reco_sv_y_1", &reco_sv_y_1_);
      mvatree_->Branch("reco_sv_z_1", &reco_sv_z_1_);
      mvatree_->Branch("reco_ip_x_1", &reco_ip_x_1_);
      mvatree_->Branch("reco_ip_y_1", &reco_ip_y_1_);
      mvatree_->Branch("reco_ip_z_1", &reco_ip_z_1_);
      mvatree_->Branch("reco_dm_1",  &reco_dm_1_);
      mvatree_->Branch("reco_metx"         , &reco_metx_         );
      mvatree_->Branch("reco_mety"         , &reco_mety_         );
      mvatree_->Branch("reco_pi_px_1",  &reco_pi_px_1_);
      mvatree_->Branch("reco_pi_py_1",  &reco_pi_py_1_);
      mvatree_->Branch("reco_pi_pz_1",  &reco_pi_pz_1_);
      mvatree_->Branch("reco_pi_E_1",   &reco_pi_E_1_);
      mvatree_->Branch("reco_pi2_px_1",  &reco_pi2_px_1_);
      mvatree_->Branch("reco_pi2_py_1",  &reco_pi2_py_1_);
      mvatree_->Branch("reco_pi2_pz_1",  &reco_pi2_pz_1_);
      mvatree_->Branch("reco_pi2_E_1",   &reco_pi2_E_1_);
      mvatree_->Branch("reco_pi3_px_1",  &reco_pi3_px_1_);
      mvatree_->Branch("reco_pi3_py_1",  &reco_pi3_py_1_);
      mvatree_->Branch("reco_pi3_pz_1",  &reco_pi3_pz_1_);
      mvatree_->Branch("reco_pi3_E_1",   &reco_pi3_E_1_);
      mvatree_->Branch("reco_pi0_px_1", &reco_pi0_px_1_);
      mvatree_->Branch("reco_pi0_py_1", &reco_pi0_py_1_);
      mvatree_->Branch("reco_pi0_pz_1", &reco_pi0_pz_1_);
      mvatree_->Branch("reco_pi0_E_1",  &reco_pi0_E_1_);

      mvatree_->Branch("reco_sv_x_2", &reco_sv_x_2_);
      mvatree_->Branch("reco_sv_y_2", &reco_sv_y_2_);
      mvatree_->Branch("reco_sv_z_2", &reco_sv_z_2_);
      mvatree_->Branch("reco_ip_x_2", &reco_ip_x_2_);
      mvatree_->Branch("reco_ip_y_2", &reco_ip_y_2_);
      mvatree_->Branch("reco_ip_z_2", &reco_ip_z_2_);
      mvatree_->Branch("reco_dm_2",  &reco_dm_2_);
      mvatree_->Branch("reco_pi_px_2",  &reco_pi_px_2_);
      mvatree_->Branch("reco_pi_py_2",  &reco_pi_py_2_);
      mvatree_->Branch("reco_pi_pz_2",  &reco_pi_pz_2_);
      mvatree_->Branch("reco_pi_E_2",   &reco_pi_E_2_);
      mvatree_->Branch("reco_pi2_px_2",  &reco_pi2_px_2_);
      mvatree_->Branch("reco_pi2_py_2",  &reco_pi2_py_2_);
      mvatree_->Branch("reco_pi2_pz_2",  &reco_pi2_pz_2_);
      mvatree_->Branch("reco_pi2_E_2",   &reco_pi2_E_2_);
      mvatree_->Branch("reco_pi3_px_2",  &reco_pi3_px_2_);
      mvatree_->Branch("reco_pi3_py_2",  &reco_pi3_py_2_);
      mvatree_->Branch("reco_pi3_pz_2",  &reco_pi3_pz_2_);
      mvatree_->Branch("reco_pi3_E_2",   &reco_pi3_E_2_);
      mvatree_->Branch("reco_pi0_px_2", &reco_pi0_px_2_);
      mvatree_->Branch("reco_pi0_py_2", &reco_pi0_py_2_);
      mvatree_->Branch("reco_pi0_pz_2", &reco_pi0_pz_2_);
      mvatree_->Branch("reco_pi0_E_2",  &reco_pi0_E_2_);	  

    }

      
    rand = new TRandom3(0);
    if(fs_){  
      outtree_ = fs_->make<TTree>("gen_ntuple","gen_ntuple");
      outtree_->Branch("event"       , &event_       );
      outtree_->Branch("wt"       , &wt_       );
      outtree_->Branch("hasTaus"       , &hasTaus_       );
      outtree_->Branch("gen_wt"       , &gen_wt_       );
      outtree_->Branch("wt_nnlops"       , &wt_nnlops_       );
      outtree_->Branch("wt_z_pol", &wt_z_pol_);
      outtree_->Branch("wt_stitch"       , &wt_stitch_       );
      outtree_->Branch("wt_topmass"       , &wt_topmass_       );
      outtree_->Branch("wt_topmass_2"     , &wt_topmass_2_  );
      outtree_->Branch("wt_ps_down", &wt_ps_down_);
      outtree_->Branch("wt_ps_up", &wt_ps_up_);
      outtree_->Branch("wt_ue_down", &wt_ue_down_);
      outtree_->Branch("wt_ue_up", &wt_ue_up_);
      outtree_->Branch("npNLO", &npNLO_);
      outtree_->Branch("cand_1"       , &cand_1_       );
      outtree_->Branch("cand_2"       , &cand_2_       );
      outtree_->Branch("match_1"       , &match_1_       );
      outtree_->Branch("match_2"       , &match_2_       );

      outtree_->Branch("wt_cp_sm"       , &wt_cp_sm_);
      outtree_->Branch("wt_cp_ps"       , &wt_cp_ps_);
      outtree_->Branch("wt_cp_mm"       , &wt_cp_mm_);
      outtree_->Branch("WTp"       , &WTp_);
      outtree_->Branch("WTm"       , &WTm_);

      outtree_->Branch("wt_cp_sm_alt"       , &wt_cp_sm_alt_);
      outtree_->Branch("wt_cp_ps_alt"       , &wt_cp_ps_alt_);
      outtree_->Branch("wt_cp_mm_alt"       , &wt_cp_mm_alt_);

      outtree_->Branch("wt_cp_prod_sm",&wt_cp_prod_sm_);
      outtree_->Branch("wt_cp_prod_ps",&wt_cp_prod_ps_);
      outtree_->Branch("wt_cp_prod_mm",&wt_cp_prod_mm_);

      outtree_->Branch("wt_vlq_sm",&wt_vlq_sm_);
      outtree_->Branch("wt_vlq_full",&wt_vlq_full_);

      outtree_->Branch("Ediff_1"       , &Ediff_1_);
      outtree_->Branch("Ediff_2"       , &Ediff_2_);

      outtree_->Branch("fakea1_dR", &fakea1_dR_);
      outtree_->Branch("rho_dR", &rho_dR_);
      outtree_->Branch("tauFlag_1", &tauFlag_1_);
      outtree_->Branch("tauFlag_2", &tauFlag_2_);

      outtree_->Branch("wt_vlq_mur1_muf0p5",    &scale_vlq_1_); //1034
      outtree_->Branch("wt_vlq_mur1_muf2",    &scale_vlq_2_); //1019
      outtree_->Branch("wt_vlq_mur2_muf1",    &scale_vlq_3_);
      outtree_->Branch("wt_vlq_mur0p5_muf1",    &scale_vlq_4_);
      outtree_->Branch("wt_vlq_mur2_muf2",    &scale_vlq_5_);
      outtree_->Branch("wt_vlq_mur0p5_muf2",    &scale_vlq_6_);
      outtree_->Branch("wt_vlq_mur2_muf0p5",    &scale_vlq_7_);
      outtree_->Branch("wt_vlq_mur0p5_muf0p5",    &scale_vlq_8_);
      outtree_->Branch("dm_1", &tauFlag_1_);
      outtree_->Branch("dm_2", &tauFlag_2_);

      if(do_theory_uncert_){
        outtree_->Branch("wt_mur1_muf1",    &scale1_);
        outtree_->Branch("wt_mur1_muf2",    &scale2_);
        outtree_->Branch("wt_mur1_muf0p5",  &scale3_);
        outtree_->Branch("wt_mur2_muf1",    &scale4_);
        outtree_->Branch("wt_mur2_muf2",    &scale5_);
        outtree_->Branch("wt_mur2_muf0p5",  &scale6_);
        outtree_->Branch("wt_mur0p5_muf1",  &scale7_);
        outtree_->Branch("wt_mur0p5_muf2",  &scale8_);
        outtree_->Branch("wt_mur0p5_muf0p5",&scale9_);

        outtree_->Branch("wt_pdf_0", &pdfweight_0_);
        outtree_->Branch("wt_pdf_1", &pdfweight_1_);
        outtree_->Branch("wt_pdf_2", &pdfweight_2_);
        outtree_->Branch("wt_pdf_3", &pdfweight_3_);
        outtree_->Branch("wt_pdf_4", &pdfweight_4_);
        outtree_->Branch("wt_pdf_5", &pdfweight_5_);
        outtree_->Branch("wt_pdf_6", &pdfweight_6_);
        outtree_->Branch("wt_pdf_7", &pdfweight_7_);
        outtree_->Branch("wt_pdf_8", &pdfweight_8_);
        outtree_->Branch("wt_pdf_9", &pdfweight_9_);
        outtree_->Branch("wt_pdf_10", &pdfweight_10_);
        outtree_->Branch("wt_pdf_11", &pdfweight_11_);
        outtree_->Branch("wt_pdf_12", &pdfweight_12_);
        outtree_->Branch("wt_pdf_13", &pdfweight_13_);
        outtree_->Branch("wt_pdf_14", &pdfweight_14_);
        outtree_->Branch("wt_pdf_15", &pdfweight_15_);
        outtree_->Branch("wt_pdf_16", &pdfweight_16_);
        outtree_->Branch("wt_pdf_17", &pdfweight_17_);
        outtree_->Branch("wt_pdf_18", &pdfweight_18_);
        outtree_->Branch("wt_pdf_19", &pdfweight_19_);
        outtree_->Branch("wt_pdf_20", &pdfweight_20_);
        outtree_->Branch("wt_pdf_21", &pdfweight_21_);
        outtree_->Branch("wt_pdf_22", &pdfweight_22_);
        outtree_->Branch("wt_pdf_23", &pdfweight_23_);
        outtree_->Branch("wt_pdf_24", &pdfweight_24_);
        outtree_->Branch("wt_pdf_25", &pdfweight_25_);
        outtree_->Branch("wt_pdf_26", &pdfweight_26_);
        outtree_->Branch("wt_pdf_27", &pdfweight_27_);
        outtree_->Branch("wt_pdf_28", &pdfweight_28_);
        outtree_->Branch("wt_pdf_29", &pdfweight_29_);
        outtree_->Branch("wt_pdf_30", &pdfweight_30_);

        outtree_->Branch("wt_h_tb",    &wt_h_tb_);
        outtree_->Branch("wt_h_t",     &wt_h_t_);
        outtree_->Branch("wt_h_b",     &wt_h_b_);
        outtree_->Branch("wt_A_tb",    &wt_A_tb_);
        outtree_->Branch("wt_A_t",     &wt_A_t_);
        outtree_->Branch("wt_A_b",     &wt_A_b_);
        outtree_->Branch("wt_H_tb",    &wt_H_tb_);
        outtree_->Branch("wt_H_t",     &wt_H_t_);
        outtree_->Branch("wt_H_b",     &wt_H_b_);

        outtree_->Branch("wt_hfact_nom",     &wt_hfact_nom_);
        outtree_->Branch("wt_hfact_up",     &wt_hfact_up_);
        outtree_->Branch("wt_hfact_down",     &wt_hfact_down_);

        outtree_->Branch("wt_lhe_nominal",     &wt_lhe_nominal_);
      }
      outtree_->Branch("tau_pt_1_tt"        , &tau_pt_1_tt_        );
      outtree_->Branch("tau_pt_1_mt"        , &tau_pt_1_mt_        );
      outtree_->Branch("tau_pt_1_et"        , &tau_pt_1_et_        );
      outtree_->Branch("tau_pt_1_sf"        , &tau_pt_1_sf_        );
      outtree_->Branch("gen_tau_pt_1"        , &gen_tau_pt_1_        );
      outtree_->Branch("gen_tau_eta_1"        , &gen_tau_eta_1_        );
      outtree_->Branch("gen_tau_dm_1"        , &gen_tau_dm_1_        );
      outtree_->Branch("passed"       ,&passed_       );
      outtree_->Branch("pt_1"        , &pt_1_        );
      outtree_->Branch("pt_2"        , &pt_2_        );
      outtree_->Branch("eta_1"       , &eta_1_       );
      outtree_->Branch("eta_2"       , &eta_2_       );
      outtree_->Branch("phi_1"       , &phi_1_       );
      outtree_->Branch("phi_2"       , &phi_2_       );
      outtree_->Branch("met"         , &met_         );
      outtree_->Branch("m_vis"       , &m_vis_       );
      outtree_->Branch("pt_tt"       , &pt_tt_       );
      outtree_->Branch("mt_tot"       , &mt_tot_       );
      outtree_->Branch("mass"       , &mass_       );
      outtree_->Branch("wtzpt"       , &wtzpt_       );
      outtree_->Branch("mt_1"        , &mt_1_        );
      outtree_->Branch("mt_2"        , &mt_2_        );
      outtree_->Branch("pzeta"       , &pzeta_       );
      outtree_->Branch("n_bjets"     , &n_bjets_     );
      outtree_->Branch("n_bjets_pt25"     , &n_bjets_pt25_     );
      outtree_->Branch("n_bjets_eta2p5"     , &n_bjets_eta2p5_     );
      outtree_->Branch("n_bjets_noscale"     , &n_bjets_noscale_);
      outtree_->Branch("n_bjets_eta2p5_noscale"     , &n_bjets_eta2p5_noscale_);
      outtree_->Branch("n_bpartons"     , &n_bpartons_     );
      outtree_->Branch("n_jets"      , &n_jets_      );
      outtree_->Branch("n_jets_nofilter"      , &n_jets_nofilter_);
      outtree_->Branch("n_jetsingap" , &n_jetsingap_ );
      outtree_->Branch("jpt_1"       , &jpt_1_       );
      outtree_->Branch("jpt_2"       , &jpt_2_       );
      outtree_->Branch("jpt_3"       , &jpt_3_       );
      outtree_->Branch("bpt_1"       , &bpt_1_       );
      outtree_->Branch("jeta_1"      , &jeta_1_      );
      outtree_->Branch("jeta_2"      , &jeta_2_      );
      outtree_->Branch("jphi_1"      , &jphi_1_      );
      outtree_->Branch("jphi_2"      , &jphi_2_      );
      outtree_->Branch("mjj"         , &mjj_         );
      outtree_->Branch("largest_gen_mjj",   &largest_gen_mjj_);
      outtree_->Branch("jdeta"       , &jdeta_       );
      outtree_->Branch("higgsDecay"  , &decayType    );
      outtree_->Branch("genpt_1"     , &genpt_1_        );
      outtree_->Branch("genpt_2"     , &genpt_2_        );
      outtree_->Branch("geneta_2"    , &geneta_2_       );
      outtree_->Branch("geneta_1"    , &geneta_1_       );
      outtree_->Branch("HiggsPt"     , &HiggsPt_     );
      outtree_->Branch("HpT"     , &pT_A_     );
      outtree_->Branch("partons"     , &partons_);
      outtree_->Branch("partons_lhe"     , &partons_lhe_);
      outtree_->Branch("quarks_initial"     , &quarks_initial_);
      outtree_->Branch("quarks_final"     , &quarks_final_);
      outtree_->Branch("gluons_initial"     , &gluons_initial_);
      outtree_->Branch("gluons_final"     , &gluons_final_);
      outtree_->Branch("parton_pt"     , &parton_pt_);
      outtree_->Branch("parton_pt_2"     , &parton_pt_2_);
      outtree_->Branch("parton_pt_3"     , &parton_pt_3_);
      outtree_->Branch("parton_mjj",    &parton_mjj_);
      outtree_->Branch("parton_HpT", &parton_HpT_);
      outtree_->Branch("parton_Zmass", &parton_Zmass_);
      outtree_->Branch("D0"     , &D0_);
      outtree_->Branch("D0star"     , &D0star_);
      outtree_->Branch("DCP"     , &DCP_);
      outtree_->Branch("sjdphi"     , &sjdphi_);
      outtree_->Branch("spjdphi"     , &spjdphi_);
      outtree_->Branch("ysep"     , &ysep_);
      outtree_->Branch("n_pjets"     , &n_pjets_);
      outtree_->Branch("n_pu",      &n_pu_);
      outtree_->Branch("lead_b_pt", &lead_b_pt_);
      outtree_->Branch("lead_b_eta", &lead_b_eta_); 
      outtree_->Branch("aco_sign", &aco_sign_);
      outtree_->Branch("aco_sign_smear", &aco_sign_smear_);
      outtree_->Branch("aco_angle_smear", &aco_angle_smear_);
      outtree_->Branch("aco_angle_1", &aco_angle_1_);
      outtree_->Branch("pv_aco_angle", &pv_aco_angle_);
      outtree_->Branch("pv_lin_angle", &pv_lin_angle_);


      outtree_->Branch("costheta_1", &costheta_1_);
      outtree_->Branch("costheta_tf_1", &costheta_tf_1_);
      outtree_->Branch("costheta_2", &costheta_2_);
      outtree_->Branch("costheta_tf_2", &costheta_tf_2_);
      
      outtree_->Branch("q_2", &q_2_);
      outtree_->Branch("q_1", &q_1_);

      outtree_->Branch("aco_angle_2", &aco_angle_2_);
      outtree_->Branch("aco_angle_3", &aco_angle_3_);
      outtree_->Branch("aco_angle_4", &aco_angle_4_);
      outtree_->Branch("aco_lin_angle_1", &aco_lin_angle_1_);
      outtree_->Branch("angle_1", &angle_1_);
      outtree_->Branch("spin_var", &spin_var_);
      outtree_->Branch("spin_var_vis", &spin_var_vis_);
      outtree_->Branch("mass_undecayed_taus", &mass_undecayed_taus_);

      outtree_->Branch("t_E_frac_1", &t_E_frac_1_);
      outtree_->Branch("t_E_frac_2", &t_E_frac_2_);
      outtree_->Branch("t_c_1", &t_c_1_);
      outtree_->Branch("t_c_2", &t_c_2_);

      outtree_->Branch("p_aco_angle_1", &p_aco_angle_1_);
      outtree_->Branch("p_aco_angle_2", &p_aco_angle_2_);
      outtree_->Branch("p_aco_angle_3", &p_aco_angle_3_);
      outtree_->Branch("p_aco_angle_4", &p_aco_angle_4_);
      outtree_->Branch("dR1", &dR1_);
      outtree_->Branch("dR2", &dR2_); 
      outtree_->Branch("aco_angle_5", &aco_angle_5_);
      outtree_->Branch("aco_angle_6", &aco_angle_6_);
      outtree_->Branch("cp_sign_1",     &cp_sign_1_);
      outtree_->Branch("cp_sign_2",     &cp_sign_2_);
      outtree_->Branch("cp_sign_3",     &cp_sign_3_);
      outtree_->Branch("cp_sign_4",     &cp_sign_4_);
      outtree_->Branch("cp_sign_5",     &cp_sign_5_);
      outtree_->Branch("cp_sign_6",     &cp_sign_6_);
      outtree_->Branch("rand", &rand_);
      outtree_->Branch("y_1_1",     &y_1_1_);
      outtree_->Branch("y_2_1",     &y_3_1_);
      outtree_->Branch("y_3_1",     &y_4_1_);
      outtree_->Branch("y_4_1",     &y_4_1_);
      outtree_->Branch("y_1_2",     &y_1_2_);
      outtree_->Branch("y_2_2",     &y_3_2_);
      outtree_->Branch("y_3_2",     &y_4_2_);
      outtree_->Branch("y_4_2",     &y_4_2_);

      outtree_->Branch("a1_mass", &a1_mass_);
      outtree_->Branch("E_1_1", &E_1_1_);
      outtree_->Branch("E_2_1", &E_2_1_);
      outtree_->Branch("E_2_2", &E_1_2_);
      outtree_->Branch("E_2_2", &E_2_2_);
      outtree_->Branch("E_3_2", &E_3_2_);


      outtree_->Branch("cp_channel",    &cp_channel_);
      outtree_->Branch("pi0_E_1", &pi0_E_1_);
      outtree_->Branch("pi0_E_2", &pi0_E_2_);
      outtree_->Branch("pi_E_1", &pi_E_1_);
      outtree_->Branch("nu_E_1",  &nu_E_1_);
      outtree_->Branch("pi_E_2", &pi_E_2_);
      outtree_->Branch("mass1", &mass1_);
      outtree_->Branch("mass2", &mass2_);    
      outtree_->Branch("pt_rho_pi", &pt_rho_pi_);
      outtree_->Branch("pt_rho_pi0", &pt_rho_pi0_);
      outtree_->Branch("pt_other_pi0", &pt_other_pi0_); 
 
      outtree_->Branch("wt_ggh_t", &wt_ggh_t_);
      outtree_->Branch("wt_ggh_b", &wt_ggh_b_);
      outtree_->Branch("wt_ggh_i", &wt_ggh_i_);
      outtree_->Branch("wt_ggH_t", &wt_ggH_t_);
      outtree_->Branch("wt_ggH_b", &wt_ggH_b_);
      outtree_->Branch("wt_ggH_i", &wt_ggH_i_);
      outtree_->Branch("wt_ggA_t", &wt_ggA_t_);
      outtree_->Branch("wt_ggA_b", &wt_ggA_b_);
      outtree_->Branch("wt_ggA_i", &wt_ggA_i_);
      outtree_->Branch("wt_ggA_i", &wt_ggA_i_);
      outtree_->Branch("pT_A", &pT_A_);

      outtree_->Branch("reco_pvx", &reco_pvx_);
      outtree_->Branch("reco_pvy", &reco_pvy_);
      outtree_->Branch("reco_pvz", &reco_pvz_);

      outtree_->Branch("reco_bs_pvx", &reco_bs_pvx_);
      outtree_->Branch("reco_bs_pvy", &reco_bs_pvy_);
      outtree_->Branch("reco_bs_pvz", &reco_bs_pvz_);

      outtree_->Branch("gen_pvx",  &gen_pvx_);
      outtree_->Branch("gen_pvy",  &gen_pvy_);
      outtree_->Branch("gen_pvz",  &gen_pvz_);

      outtree_->Branch("wt_ps_isr_up" , &wt_ps_isr_up_ );
      outtree_->Branch("wt_ps_isr_down" , &wt_ps_isr_down_ );
      outtree_->Branch("wt_ps_fsr_up" , &wt_ps_fsr_up_ );
      outtree_->Branch("wt_ps_fsr_down" , &wt_ps_fsr_down_ );

      outtree_->Branch("count_taus" , &count_taus_ );

      outtree_->Branch("gammas_deta_1", &gammas_deta_1_);
      outtree_->Branch("gammas_dphi_1", &gammas_dphi_1_);
      outtree_->Branch("gammas_deta_2", &gammas_deta_2_);
      outtree_->Branch("gammas_dphi_2", &gammas_dphi_2_);
     
      outtree_->Branch("pi0s_deta_1", &pi0s_deta_1_);
      outtree_->Branch("pi0s_dphi_1", &pi0s_dphi_1_);
      outtree_->Branch("pi0s_deta_2", &pi0s_deta_2_);
      outtree_->Branch("pi0s_dphi_2", &pi0s_dphi_2_);

      outtree_->Branch("rho_deta_1", &rho_deta_1_);
      outtree_->Branch("rho_dphi_1", &rho_dphi_1_);
      outtree_->Branch("rho_deta_2", &rho_deta_2_);
      outtree_->Branch("rho_dphi_2", &rho_dphi_2_);
 
    }
    count_ee_ = 0;
    count_em_ = 0;
    count_et_ = 0;
    count_mm_ = 0;
    count_mt_ = 0;
    count_tt_ = 0;


    /* topmass_wts_toponly_ = GetFromTFile<TH1F>("input/ggh_weights/quarkmass_uncerts_hnnlo.root","/","nom"); */
   
    /* ps_0jet_up_ = GetFromTFile<TH1D>("input/ggh_weights/MG_ps_uncerts.root","/","ps_0jet_up"); */
    /* ps_0jet_down_ = GetFromTFile<TH1D>("input/ggh_weights/MG_ps_uncerts.root","/","ps_0jet_down"); */
    /* ps_1jet_up_ = GetFromTFile<TH1D>("input/ggh_weights/MG_ps_uncerts.root","/","ps_1jet_up"); */
    /* ps_1jet_down_ = GetFromTFile<TH1D>("input/ggh_weights/MG_ps_uncerts.root","/","ps_1jet_down"); */
    /* ps_2jet_up_ = GetFromTFile<TH1D>("input/ggh_weights/MG_ps_uncerts.root","/","ps_2jet_up"); */
    /* ps_2jet_down_ = GetFromTFile<TH1D>("input/ggh_weights/MG_ps_uncerts.root","/","ps_2jet_down"); */
    /* ps_3jet_up_ = GetFromTFile<TH1D>("input/ggh_weights/MG_ps_uncerts.root","/","ps_3jet_up"); */
    /* ps_3jet_down_ = GetFromTFile<TH1D>("input/ggh_weights/MG_ps_uncerts.root","/","ps_3jet_down"); */  

    std::string csv_file_path = "";

    if (era_==era::data_2016)       csv_file_path = "./input/btag_sf/DeepJet_2016LegacySF_V1.csv";
    else if (era_==era::data_2017)  csv_file_path = "./input/btag_sf/DeepFlavour_94XSF_V4_B_F.csv";
    else if (era_==era::data_2018)  csv_file_path = "./input/btag_sf/DeepJet_102XSF_V2.csv";
    
    //calib  = new const BTagCalibration("deepjet",csv_file_path);
    
    
    //reader_comb = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"});
    //reader_comb->load(*calib,BTagEntry::FLAV_B,"comb"); 
    
    
    
    std::string file = "input/ggh_weights/NNLOPS_reweight.root";
    ggh_weights_ = new TFile(file.c_str());
    ggh_weights_->cd();
    
    ggh_ph_0jet_ = (TGraph*)gDirectory->Get("gr_NNLOPSratio_pt_powheg_0jet");
    ggh_ph_1jet_ = (TGraph*)gDirectory->Get("gr_NNLOPSratio_pt_powheg_1jet");
    ggh_ph_2jet_ = (TGraph*)gDirectory->Get("gr_NNLOPSratio_pt_powheg_2jet");
    ggh_ph_3jet_ = (TGraph*)gDirectory->Get("gr_NNLOPSratio_pt_powheg_3jet");
    
    ggh_weights_->Close();
    return 0;
  }

  int HTTGenAnalysis::Execute(TreeEvent *event) {
    
    EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
    event_ = (unsigned long long) eventInfo->event();
    wt_ = 1;
    rand->SetSeed(event_);
    rand_ = rand->Uniform();   
    hasTaus_=false; 
    wt_ = eventInfo->total_weight();
  
    if (eventInfo->weight_defined("gen_weight")) gen_wt_ = eventInfo->weight("gen_weight");
  
    //std::cout << "================event " << event_ << "===================" << std::endl;
  
    wt_cp_sm_=0.0; 
    wt_cp_ps_=0.0; 
    wt_cp_mm_=0.0; 
  
    // test tauspinner weights
    if(event->ExistsInTree("tauspinner")){
      EventInfo const* tauspinner = event->GetPtr<EventInfo>("tauspinner");
  
      if (tauspinner->weight_defined("wt_z_pol")) wt_z_pol_ = tauspinner->weight("wt_z_pol"); else wt_z_pol_=1.0;
  
      if (tauspinner->weight_defined("wt_cp_0")){
        wt_cp_sm_ = tauspinner->weight("wt_cp_0");
        wt_cp_ps_ = tauspinner->weight("wt_cp_0p5"); 
        wt_cp_mm_ = tauspinner->weight("wt_cp_0p25");
  
        //wt_cp_sm_alt_ = tauspinner->weight("wt_cp_0_alt");
        //wt_cp_ps_alt_ = tauspinner->weight("wt_cp_0p5_alt");
        //wt_cp_mm_alt_ = tauspinner->weight("wt_cp_0p25_alt");
        //WTm_ = tauspinner->weight("WTm");
        //WTp_ = tauspinner->weight("WTp");
      }
  
    }
  
  
    wt_cp_prod_sm_=0.0;
    wt_cp_prod_ps_=0.0;
    wt_cp_prod_mm_=0.0;
  
    if (eventInfo->weight_defined("sm_weight_nlo")){
      if(eventInfo->weight_defined("sm_weight_nlo")) wt_cp_prod_sm_ = eventInfo->weight("sm_weight_nlo");
      if(eventInfo->weight_defined("ps_weight_nlo")) wt_cp_prod_ps_ = eventInfo->weight("ps_weight_nlo");
      if(eventInfo->weight_defined("mm_weight_nlo")) wt_cp_prod_mm_ = eventInfo->weight("mm_weight_nlo");
    } else {
      if(eventInfo->weight_defined("sm_weight")) wt_cp_prod_sm_ = eventInfo->weight("sm_weight");
      if(eventInfo->weight_defined("ps_weight")) wt_cp_prod_ps_ = eventInfo->weight("ps_weight");
      if(eventInfo->weight_defined("mm_weight")) wt_cp_prod_mm_ = eventInfo->weight("mm_weight");
   }
  
  
    wt_ps_isr_up_   = eventInfo->weight_defined("genweight6") ? eventInfo->weight("genweight6") : 1.0;
    wt_ps_isr_down_ = eventInfo->weight_defined("genweight8") ? eventInfo->weight("genweight8") : 1.0;
    wt_ps_fsr_up_   = eventInfo->weight_defined("genweight7") ? eventInfo->weight("genweight7") : 1.0;
    wt_ps_fsr_down_ = eventInfo->weight_defined("genweight9") ? eventInfo->weight("genweight9") : 1.0;
  
   scale_vlq_1_ = eventInfo->weight_defined("1034") ? eventInfo->weight("1034") : 1.0; // 1 0.5
   scale_vlq_2_ = eventInfo->weight_defined("1019") ? eventInfo->weight("1019") : 1.0; // 1 2
   scale_vlq_3_ = eventInfo->weight_defined("1009") ? eventInfo->weight("1009") : 1.0; // 2 1
   scale_vlq_4_ = eventInfo->weight_defined("1014") ? eventInfo->weight("1014") : 1.0; // 0.5 1
   scale_vlq_5_ = eventInfo->weight_defined("1024") ? eventInfo->weight("1024") : 1.0; // 2 2
   scale_vlq_6_ = eventInfo->weight_defined("1029") ? eventInfo->weight("1029") : 1.0; // 0.5 2
   scale_vlq_7_ = eventInfo->weight_defined("1039") ? eventInfo->weight("1039") : 1.0; // 2 0.5
   scale_vlq_8_ = eventInfo->weight_defined("1044") ? eventInfo->weight("1044") : 1.0; // 0.5 0.5
  
   wt_vlq_sm_ = eventInfo->weight_defined("sm") ? eventInfo->weight("sm") : 1.0;
   wt_vlq_full_ = eventInfo->weight_defined("full") ? eventInfo->weight("full") : 1.0;
  
   if(do_theory_uncert_){
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


      if(eventInfo->weight_defined("5000")) pdfweight_0_ = eventInfo->weight("5000"); else if(eventInfo->weight_defined("10000")) pdfweight_0_ = eventInfo->weight("10000"); else pdfweight_0_=1.0;
      if(eventInfo->weight_defined("5001")) pdfweight_1_ = eventInfo->weight("5001"); else if(eventInfo->weight_defined("10001")) pdfweight_1_ = eventInfo->weight("10001"); else pdfweight_1_=1.0;
      if(eventInfo->weight_defined("5002")) pdfweight_2_ = eventInfo->weight("5002"); else if(eventInfo->weight_defined("10002")) pdfweight_2_ = eventInfo->weight("10002"); else pdfweight_2_=1.0;
      if(eventInfo->weight_defined("5003")) pdfweight_3_ = eventInfo->weight("5003"); else if(eventInfo->weight_defined("10003")) pdfweight_3_ = eventInfo->weight("10003"); else pdfweight_3_=1.0;
      if(eventInfo->weight_defined("5004")) pdfweight_4_ = eventInfo->weight("5004"); else if(eventInfo->weight_defined("10004")) pdfweight_4_ = eventInfo->weight("10004"); else pdfweight_4_=1.0;
      if(eventInfo->weight_defined("5005")) pdfweight_5_ = eventInfo->weight("5005"); else if(eventInfo->weight_defined("10005")) pdfweight_5_ = eventInfo->weight("10005"); else pdfweight_5_=1.0;
      if(eventInfo->weight_defined("5006")) pdfweight_6_ = eventInfo->weight("5006"); else if(eventInfo->weight_defined("10006")) pdfweight_6_ = eventInfo->weight("10006"); else pdfweight_6_=1.0;
      if(eventInfo->weight_defined("5007")) pdfweight_7_ = eventInfo->weight("5007"); else if(eventInfo->weight_defined("10007")) pdfweight_7_ = eventInfo->weight("10007"); else pdfweight_7_=1.0;
      if(eventInfo->weight_defined("5008")) pdfweight_8_ = eventInfo->weight("5008"); else if(eventInfo->weight_defined("10008")) pdfweight_8_ = eventInfo->weight("10008"); else pdfweight_8_=1.0;
      if(eventInfo->weight_defined("5009")) pdfweight_9_ = eventInfo->weight("5009"); else if(eventInfo->weight_defined("10009")) pdfweight_9_ = eventInfo->weight("10009"); else pdfweight_9_=1.0;
      if(eventInfo->weight_defined("5010")) pdfweight_10_ = eventInfo->weight("5010"); else if(eventInfo->weight_defined("10010")) pdfweight_10_ = eventInfo->weight("10010"); else pdfweight_10_=1.0;
      if(eventInfo->weight_defined("5011")) pdfweight_11_ = eventInfo->weight("5011"); else if(eventInfo->weight_defined("10011")) pdfweight_11_ = eventInfo->weight("10011"); else pdfweight_11_=1.0;
      if(eventInfo->weight_defined("5012")) pdfweight_12_ = eventInfo->weight("5012"); else if(eventInfo->weight_defined("10012")) pdfweight_12_ = eventInfo->weight("10012"); else pdfweight_12_=1.0;
      if(eventInfo->weight_defined("5013")) pdfweight_13_ = eventInfo->weight("5013"); else if(eventInfo->weight_defined("10013")) pdfweight_13_ = eventInfo->weight("10013"); else pdfweight_13_=1.0;
      if(eventInfo->weight_defined("5014")) pdfweight_14_ = eventInfo->weight("5014"); else if(eventInfo->weight_defined("10014")) pdfweight_14_ = eventInfo->weight("10014"); else pdfweight_14_=1.0;
      if(eventInfo->weight_defined("5015")) pdfweight_15_ = eventInfo->weight("5015"); else if(eventInfo->weight_defined("10015")) pdfweight_15_ = eventInfo->weight("10015"); else pdfweight_15_=1.0;
      if(eventInfo->weight_defined("5016")) pdfweight_16_ = eventInfo->weight("5016"); else if(eventInfo->weight_defined("10016")) pdfweight_16_ = eventInfo->weight("10016"); else pdfweight_16_=1.0;
      if(eventInfo->weight_defined("5017")) pdfweight_17_ = eventInfo->weight("5017"); else if(eventInfo->weight_defined("10017")) pdfweight_17_ = eventInfo->weight("10017"); else pdfweight_17_=1.0;
      if(eventInfo->weight_defined("5018")) pdfweight_18_ = eventInfo->weight("5018"); else if(eventInfo->weight_defined("10018")) pdfweight_18_ = eventInfo->weight("10018"); else pdfweight_18_=1.0;
      if(eventInfo->weight_defined("5019")) pdfweight_19_ = eventInfo->weight("5019"); else if(eventInfo->weight_defined("10019")) pdfweight_19_ = eventInfo->weight("10019"); else pdfweight_19_=1.0;
      if(eventInfo->weight_defined("5020")) pdfweight_20_ = eventInfo->weight("5020"); else if(eventInfo->weight_defined("10020")) pdfweight_20_ = eventInfo->weight("10020"); else pdfweight_20_=1.0;
      if(eventInfo->weight_defined("5021")) pdfweight_21_ = eventInfo->weight("5021"); else if(eventInfo->weight_defined("10021")) pdfweight_21_ = eventInfo->weight("10021"); else pdfweight_21_=1.0;
      if(eventInfo->weight_defined("5022")) pdfweight_22_ = eventInfo->weight("5022"); else if(eventInfo->weight_defined("10022")) pdfweight_22_ = eventInfo->weight("10022"); else pdfweight_22_=1.0;
      if(eventInfo->weight_defined("5023")) pdfweight_23_ = eventInfo->weight("5023"); else if(eventInfo->weight_defined("10023")) pdfweight_23_ = eventInfo->weight("10023"); else pdfweight_23_=1.0;
      if(eventInfo->weight_defined("5024")) pdfweight_24_ = eventInfo->weight("5024"); else if(eventInfo->weight_defined("10024")) pdfweight_24_ = eventInfo->weight("10024"); else pdfweight_24_=1.0;
      if(eventInfo->weight_defined("5025")) pdfweight_25_ = eventInfo->weight("5025"); else if(eventInfo->weight_defined("10025")) pdfweight_25_ = eventInfo->weight("10025"); else pdfweight_25_=1.0;
      if(eventInfo->weight_defined("5026")) pdfweight_26_ = eventInfo->weight("5026"); else if(eventInfo->weight_defined("10026")) pdfweight_26_ = eventInfo->weight("10026"); else pdfweight_26_=1.0;
      if(eventInfo->weight_defined("5027")) pdfweight_27_ = eventInfo->weight("5027"); else if(eventInfo->weight_defined("10027")) pdfweight_27_ = eventInfo->weight("10027"); else pdfweight_27_=1.0;
      if(eventInfo->weight_defined("5028")) pdfweight_28_ = eventInfo->weight("5028"); else if(eventInfo->weight_defined("10028")) pdfweight_28_ = eventInfo->weight("10028"); else pdfweight_28_=1.0;
      if(eventInfo->weight_defined("5029")) pdfweight_29_ = eventInfo->weight("5029"); else if(eventInfo->weight_defined("10029")) pdfweight_29_ = eventInfo->weight("10029"); else pdfweight_29_=1.0;
      if(eventInfo->weight_defined("5030")) pdfweight_30_ = eventInfo->weight("5030"); else if(eventInfo->weight_defined("10030")) pdfweight_30_ = eventInfo->weight("10030"); else pdfweight_30_=1.0;

      if(eventInfo->weight_defined("h_tb")) wt_h_tb_ = eventInfo->weight("h_tb"); else wt_h_tb_=1.0;
      if(eventInfo->weight_defined("h_t"))  wt_h_t_ = eventInfo->weight("h_t"); else wt_h_t_=1.0;
      if(eventInfo->weight_defined("h_b"))  wt_h_b_ = eventInfo->weight("h_b"); else wt_h_b_=1.0;
      if(eventInfo->weight_defined("A_tb")) wt_A_tb_ = eventInfo->weight("A_tb"); else wt_A_tb_=1.0;
      if(eventInfo->weight_defined("A_t"))  wt_A_t_ = eventInfo->weight("A_t"); else wt_A_t_=1.0;
      if(eventInfo->weight_defined("A_b"))  wt_A_b_ = eventInfo->weight("A_b"); else wt_A_b_=1.0;
      if(eventInfo->weight_defined("H_tb")) wt_H_tb_ = eventInfo->weight("H_tb"); else wt_H_tb_=1.0;
      if(eventInfo->weight_defined("H_t"))  wt_H_t_ = eventInfo->weight("H_t"); else wt_H_t_=1.0;
      if(eventInfo->weight_defined("H_b"))  wt_H_b_ = eventInfo->weight("H_b"); else wt_H_b_=1.0;

      if(eventInfo->weight_defined("lhe_nominal")) wt_lhe_nominal_ = eventInfo->weight("lhe_nominal"); else wt_lhe_nominal_=1.0;

      if(eventInfo->weight_defined("hfact_nominal")) wt_hfact_nom_ = eventInfo->weight("hfact_nominal"); else wt_hfact_nom_=1.0;
      if(eventInfo->weight_defined("hfact_up")) wt_hfact_up_ = eventInfo->weight("hfact_up"); else wt_hfact_up_=1.0;
      if(eventInfo->weight_defined("hfact_down")) wt_hfact_down_ = eventInfo->weight("hfact_down"); else wt_hfact_down_=1.0;


      // 1005 = up
      // 1009 = down 

      // For MG5 samples:
      //  <weight MUF="1.0" MUR="2.0" PDF="292200" id="1002"> MUR=2.0  </weight>
      //  <weight MUF="1.0" MUR="0.5" PDF="292200" id="1003"> MUR=0.5  </weight>
      //  <weight MUF="2.0" MUR="1.0" PDF="292200" id="1004"> MUF=2.0  </weight>
      //  <weight MUF="2.0" MUR="2.0" PDF="292200" id="1005"> MUR=2.0 MUF=2.0  </weight>
      //  <weight MUF="2.0" MUR="0.5" PDF="292200" id="1006"> MUR=0.5 MUF=2.0  </weight>
      //  <weight MUF="0.5" MUR="1.0" PDF="292200" id="1007"> MUF=0.5  </weight>
      //  <weight MUF="0.5" MUR="2.0" PDF="292200" id="1008"> MUR=2.0 MUF=0.5  </weight>
      //  <weight MUF="0.5" MUR="0.5" PDF="292200" id="1009"> MUR=0.5 MUF=0.5  </weight>
      //

      // for NNLOPS:
      //  <weight id="1001"> muR=1 muF=1 </weight>
      //  <weight id="1002"> muR=1 muF=2 </weight>
      //  <weight id="1003"> muR=1 muF=0.5 </weight>
      //  <weight id="1004"> muR=2 muF=1 </weight>
      //  <weight id="1005"> muR=2 muF=2 </weight>
      //  <weight id="1006"> muR=2 muF=0.5 </weight>
      //  <weight id="1007"> muR=0.5 muF=1 </weight>
      //  <weight id="1008"> muR=0.5 muF=2 </weight>
      //  <weight id="1009"> muR=0.5 muF=0.5 </weight>

    }
    
    if (event->Exists("D0")) D0_ = event->Get<float>("D0");
    else D0_ = -9999;
    if (event->Exists("DCP")) DCP_ = event->Get<float>("DCP");
    else DCP_ = -9999; 
    if(D0_!=-9999&&DCP_!=-9999) D0star_ = D0_*DCP_/fabs(DCP_);
    else D0star_ = -9999;
    
    std::vector<ic::GenParticle*> gen_particles = event->GetPtrVec<ic::GenParticle>("genParticles");
    std::vector<ic::GenJet*> gen_jets;
    if(event->ExistsInTree("genJets")) gen_jets = event->GetPtrVec<ic::GenJet>("genJets");

    std::vector<ic::GenParticle> higgs_products;
    std::vector<ic::GenParticle> gen_taus;
    ic::Candidate met; 
    std::vector<ic::GenParticle> prompt_leptons;
    std::vector<std::string> decay_types;

    double pT=0;
    HiggsPt_=-9999;
    partons_lhe_=0;
    partons_=0;
    parton_mjj_=-9999;
    parton_HpT_=-9999;
    parton_Zmass_=-9999;
    quarks_initial_=0;
    quarks_final_=0;
    gluons_initial_=0;
    gluons_final_=0;

    //double higgs_eta = 0;
    std::vector<double> parton_pt_vec = {};
    bool lhe_exists = event->ExistsInTree("lheParticles");
    std::vector<GenParticle*> outparts;
    ROOT::Math::PtEtaPhiEVector gen_boson_lhe;
    if(lhe_exists){
      std::vector<GenParticle*> const& lhe_parts = event->GetPtrVec<GenParticle>("lheParticles");
      parton_pt_=-9999;
      double largest_b_pt=-1;
      double lead_b_eta=-9999;
      for(unsigned i = 0; i< lhe_parts.size(); ++i){
           unsigned id = abs(lhe_parts[i]->pdgid());

           if(lhe_parts[i]->status() != 1) {
             if (id >= 1 && id <=6) quarks_initial_++;
             if (id==21) gluons_initial_++;   
           }
           if(lhe_parts[i]->status() != 1) continue;

           if(id==5 && lhe_parts[i]->pt()>largest_b_pt) {
             largest_b_pt = lhe_parts[i]->pt();
             lead_b_eta=lhe_parts[i]->eta();
           }
           if(id==25) parton_HpT_ = lhe_parts[i]->pt();
           if(id==11 ||id==13 || id==15) gen_boson_lhe+=lhe_parts[i]->vector();
           if (id >= 1 && id <=6) quarks_final_++; 
           if (id==21) gluons_final_++; 
           if ((id >= 1 && id <=6) || id == 21){ 
             outparts.push_back(lhe_parts[i]);
             partons_++;
             parton_pt_vec.push_back(lhe_parts[i]->pt());
             if(lhe_parts[i]->pt()>=10) partons_lhe_++; 
        }
      }
    lead_b_pt_ = largest_b_pt;
    lead_b_eta_ = lead_b_eta; 
    } else {
        partons_=0;
        for (unsigned i = 0; i < gen_particles.size(); ++i) {
          if (!(gen_particles[i]->statusFlags()[IsHardProcess])) continue;
          if (gen_particles[i]->status()!=23) continue;
          if(gen_particles[gen_particles[i]->mothers()[0]]->pdgid() == 2212) continue;
          unsigned id = abs(gen_particles[i]->pdgid());
          if (id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 6 || id == 21 ) {
            partons_++;
            parton_pt_vec.push_back(gen_particles[i]->pt());
          }
        }
    }
    std::sort(outparts.begin(),outparts.end(),PtComparatorGenPart());
    if(outparts.size()>1) parton_mjj_ = (outparts[0]->vector()+outparts[1]->vector()).M();
    else parton_mjj_ = -9999;
    if(outparts.size()>2){
      double parton_mjj_2 = (outparts[0]->vector()+outparts[2]->vector()).M(); 
      double parton_mjj_3 = (outparts[1]->vector()+outparts[2]->vector()).M();  
      if(parton_mjj_ < std::max(parton_mjj_2, parton_mjj_3)) parton_mjj_ = std::max(parton_mjj_2, parton_mjj_3);
   
    }
    std::sort(parton_pt_vec.begin(),parton_pt_vec.end());
    std::reverse(parton_pt_vec.begin(),parton_pt_vec.end());
    if (parton_pt_vec.size()>0) parton_pt_ = parton_pt_vec[0];
    else parton_pt_ = -9999;
    if (parton_pt_vec.size()>1) parton_pt_2_ = parton_pt_vec[1];
    else parton_pt_2_ = -9999;
    if (parton_pt_vec.size()>2) parton_pt_3_ = parton_pt_vec[2];
    else parton_pt_3_ = -9999;

    npNLO_ = eventInfo->npNLO();
    //if(npNLO_<0) npNLO_ = 2; 
    double n_inc_ = 3089015.;
    double n2_    = 14254055;
    double f2_   = 0.279662;
    if(npNLO_>=2) wt_stitch_ = (n_inc_*f2_) / ( (n_inc_*f2_) + n2_ );
    else wt_stitch_=1.;
    
    ROOT::Math::PtEtaPhiEVector gen_boson;

    for(unsigned i=0; i<gen_particles.size(); ++i){
      if((gen_particles[i]->statusFlags()[FromHardProcessBeforeFSR] || gen_particles[i]->statusFlags()[IsLastCopy]) && gen_particles[i]->pdgid() == 25) {
          HiggsPt_ = gen_particles[i]->pt();
           //higgs_eta = gen_particles[i]->eta();
           //wt_topmass_ = topmass_wts_toponly_.GetBinContent(topmass_wts_toponly_.FindBin(HiggsPt_))*1.006;
           //wt_topmass_2_ = topmass_wts_.GetBinContent(topmass_wts_.FindBin(HiggsPt_))*0.985; //*sm = 0.985022, mix= 0.985167 ps=0.985076 -> all = 0.985 to 3dp so use thsi number
      }

      
      ic::GenParticle part = *gen_particles[i];
      ic::GenParticle higgs_product;
      unsigned genID = std::fabs(part.pdgid());
      bool status_flag_t = part.statusFlags().at(0);
      bool status_flag_tlc = part.statusFlags().at(13);
      bool status_hard_process = part.statusFlags().at(7);
      
      if ( (genID >= 11 && genID <= 16 && part.statusFlags()[FromHardProcess] && abs(part.status())==1) || part.statusFlags()[IsDirectHardProcessTauDecayProduct]) gen_boson+=part.vector();
      if (!lhe_exists && status_hard_process &&(genID == 1 || genID == 2 || genID == 3 || genID == 4 || genID == 5 || genID == 6 || genID == 21) && gen_particles[part.mothers().at(0)]->pdgid() != 2212 ) partons_++;

      //if(genID==25 || genID==36 || genID==35 ) std::cout << genID << std::endl;
 
      if(genID==36 && gen_particles[i]->statusFlags()[IsLastCopy]){
        pT = gen_particles[i]->vector().Pt();
        pT_A_ = pT;
      }
      if((genID==23||genID==25||genID==35||genID==36) && gen_particles[i]->statusFlags()[IsLastCopy]){
        pT = gen_particles[i]->vector().Pt();
        pT_A_ = pT;
      }

      std::vector<ic::Electron*> reco_electrons = {};//event->GetPtrVec<ic::Electron>("electrons");
      std::vector<ic::Muon*> reco_muons = {};//event->GetPtrVec<ic::Muon>("muons");     

      // add tau_neutrinos 4-vectors to get gen met
      if(genID == 12 || genID == 14 || genID == 16){
        met.set_vector(met.vector() + part.vector());
        continue;
      }
      if(channel_str_=="zmm") {
        if(!(genID == 13 && gen_particles[i]->statusFlags()[IsPrompt] && gen_particles[i]->statusFlags()[IsLastCopy])) continue;
        higgs_products.push_back(*(gen_particles[i]));
        decay_types.push_back("m");
        std::vector<ic::GenParticle *> match_muons = {gen_particles[i]};
        if(fabs(gen_particles[i]->eta())<2.4 && gen_particles[i]->pt()>20.) {
          if(decay_types.size()==1){
            cand_1_ = true;
            match_1_ = (MatchByDR(match_muons,reco_muons,0.5,true,true).size()>0);
          } else if (decay_types.size()==2) {
            cand_2_ = true;
            match_2_ = (MatchByDR(match_muons,reco_muons,0.5,true,true).size()>0);
          }
        }
      }
      if(channel_str_=="zee") {
        if(!(genID == 11 && gen_particles[i]->statusFlags()[IsPrompt] && gen_particles[i]->statusFlags()[IsLastCopy])) continue;
        higgs_products.push_back(*(gen_particles[i]));
        decay_types.push_back("e");
        std::vector<ic::GenParticle *> match_elecs = {gen_particles[i]};
        if(fabs(gen_particles[i]->eta())<2.5 && gen_particles[i]->pt()>20.) {
          if(decay_types.size()==1){
            cand_1_ = true;
            match_1_ = (MatchByDR(match_elecs,reco_electrons,0.5,true,true).size()>0);
          } else if (decay_types.size()==2) {
            cand_2_ = true;
            match_2_ = (MatchByDR(match_elecs,reco_electrons,0.5,true,true).size()>0); 
          }
        }
      }
      if(!(genID == 15 && status_flag_t && status_flag_tlc)) continue;
      hasTaus_=true;
      gen_taus.push_back(part);
      std::vector<ic::GenParticle> family;
      unsigned outputID = 15;
      FamilyTree(family, part, gen_particles, outputID);
      if(family.size()==1 && (outputID ==11 || outputID ==13)){
        higgs_products.push_back(family[0]);
        if (outputID == 11) {decay_types.push_back("e");}  
        else if (outputID == 13) {decay_types.push_back("m");}
      } else {
        decay_types.push_back("t");  
        ic::GenParticle had_tau;
        int charge = 0;
        int pdgid = 15;
        for(unsigned j=0; j<family.size(); ++j){
          had_tau.set_vector(had_tau.vector() + family[j].vector());
          charge += family[j].charge();
        }
        pdgid = 15*charge;
        had_tau.set_charge(charge);
        had_tau.set_pdgid(pdgid);
        higgs_products.push_back(had_tau);
      }
    }

    mass_=gen_boson.M();
    parton_Zmass_=gen_boson_lhe.M();

    std::sort(higgs_products.begin(),higgs_products.end(),PtComparator());
    std::sort(gen_taus.begin(),gen_taus.end(),PtComparator());
    
    if(gen_taus.size()>=2){
      genpt_1_ = gen_taus[0].vector().Pt();
      genpt_2_ = gen_taus[1].vector().Pt();
      geneta_1_ = gen_taus[0].vector().Rapidity();
      geneta_2_ = gen_taus[1].vector().Rapidity();
    } else {
      genpt_1_ =  -9999;
      genpt_2_ =  -9999;
      geneta_1_ = -9999;
      geneta_2_ = -9999; 
    }
    
    std::vector<ic::GenParticle> electrons;
    std::vector<ic::GenParticle> muons;
    std::vector<ic::GenParticle> taus;
    
    double min_tau_pt[2];
    min_tau_pt[0]     = min_tau1_pt_;
    min_tau_pt[1]     = min_tau2_pt_;
    
    for(unsigned i=0; i<higgs_products.size(); ++i){
      unsigned ID = std::fabs(higgs_products[i].pdgid());
      double eta = std::fabs(higgs_products[i].vector().Rapidity());
      double pt = higgs_products[i].vector().Pt();
      if(ID == 11){
        if(pt > min_e_pt_ && eta < max_e_eta_) electrons.push_back(higgs_products[i]);  
      } else if(ID == 13){
        if(pt > min_mu_pt_ && eta < max_mu_eta_) muons.push_back(higgs_products[i]);  
      } else if(ID == 15){
        if(pt > min_tau_pt[i] && eta < max_tau_eta_) taus.push_back(higgs_products[i]);
      }
    }

    //size of decay_types vector should always be 2 but added this if statement just to be sure
    decayType = "";
    std::sort(decay_types.begin(),decay_types.end(),swap_labels());
    for(unsigned i=0; i< decay_types.size(); ++i){
      decayType += decay_types[i];
    }
    
    if(decayType == "ee") count_ee_++;
    if(decayType == "em") count_em_++;
    if(decayType == "et") count_et_++;
    if(decayType == "mm") count_mm_++;
    if(decayType == "mt") count_mt_++;
    if(decayType == "tt") count_tt_++;

    pt_1_ = -9999.;
    pt_2_ = -9999.;
    ic::Candidate lep1;
    ic::Candidate lep2;
    passed_ = false;
    
    if(channel_str_ == "ee"){
      if(electrons.size() == 2){
        lep1 = electrons[0];
        lep2 = electrons[1];
        passed_ = true;
      }
    } else if(channel_str_ == "mm"){
      if(muons.size() == 2){
        lep1 = muons[0];
        lep2 = muons[1];
        passed_ = true;
      }
    } else if(channel_str_ == "zmm"){
      if(muons.size() == 2){
        lep1 = muons[0];
        lep2 = muons[1];
        passed_ = true;
      }
    } else if(channel_str_ == "em"){
      if(electrons.size() == 1 && muons.size() == 1){
        lep1 = electrons[0];
        lep2 = muons[0];
        passed_ = true;
      }
    } else if(channel_str_ == "et"){
      if(electrons.size() == 1 && taus.size() == 1){
        lep1 = electrons[0];
        lep2 = taus[0];
        passed_ = true;
      }
    } else if(channel_str_ == "mt"){
      if(muons.size() == 1 && taus.size() == 1){
        lep1 = muons[0];
        lep2 = taus[0];
        passed_ = true;
      }
    } else if(channel_str_ == "tt"){
      if(taus.size() == 2){
        lep1 = taus[0];
        lep2 = taus[1];
        passed_ = true;
      }
    }

    std::vector<ic::Vertex*> gen_vertices;
    if(event->ExistsInTree("genVertices")) gen_vertices = event->GetPtrVec<ic::Vertex>("genVertices");

    std::vector<GenJet> gen_tau_jets = BuildTauJets(gen_particles, false,true);
    std::vector<GenJet *> gen_tau_jets_ptr;
    for (auto & x : gen_tau_jets) gen_tau_jets_ptr.push_back(&x);
    ic::erase_if(gen_tau_jets_ptr, !boost::bind(MinPtMaxEta, _1, 15.0, 999.));
    std::sort(gen_tau_jets_ptr.begin(), gen_tau_jets_ptr.end(), bind(&Candidate::pt, _1) > bind(&Candidate::pt, _2));

    tau_pt_1_tt_=-9999; 
    tau_pt_1_mt_=-9999; 
    tau_pt_1_et_=-9999; 
    tau_pt_1_sf_=-9999; 
    gen_tau_pt_1_=-9999;
    gen_tau_dm_1_=-1;
    gen_tau_eta_1_=-9999;

    tau_pt_2_tt_=-9999;
    gen_tau_pt_2_=-9999;
    gen_tau_dm_2_=-1;
    gen_tau_eta_2_=-9999;

    // match gen tau jets to taus passing selections for et, mt, and tt channels
//    if(event->ExistsInTree("taus") && gen_tau_jets_ptr.size()>=1) {
    if(event->ExistsInTree("taus") && gen_tau_jets_ptr.size()==2) {
//      std::vector<GenJet *> lead_gen_tau = {gen_tau_jets_ptr[0]};
      std::vector<GenJet *> lead_gen_tau = {gen_tau_jets_ptr[0]};
      std::vector<GenJet *> sublead_gen_tau = {gen_tau_jets_ptr[1]};

      std::vector<Tau *> taus = event->GetPtrVec<Tau>("taus");
      ic::erase_if(taus,!boost::bind(MinPtMaxEta, _1, 20.0, 2.3));         
      std::vector<Tau *> taus_tt;
      std::vector<Tau *> taus_mt;
      std::vector<Tau *> taus_et;
      std::vector<Tau *> taus_sf;

      gen_tau_pt_1_ = lead_gen_tau[0]->pt();
      gen_tau_eta_1_ = lead_gen_tau[0]->eta();
      gen_tau_dm_1_ = lead_gen_tau[0]->flavour();

      gen_tau_pt_2_ = sublead_gen_tau[0]->pt();
      gen_tau_eta_2_ = sublead_gen_tau[0]->eta();
      gen_tau_dm_2_ = sublead_gen_tau[0]->flavour();

      for(auto t : taus) {
        if(t->GetTauID("byMediumDeepTau2017v2p1VSjet") > 0.5 && t->GetTauID("byVVLooseDeepTau2017v2p1VSe") > 0.5 && t->GetTauID("byVLooseDeepTau2017v2p1VSmu") > 0.5 && t->GetTauID("decayModeFindingNewDMs") > 0.5 && (t->decay_mode()<2 || t->decay_mode()>9) && fabs(t->charge()) == 1 && fabs(t->lead_dz_vertex()) < 0.2) {
        taus_tt.push_back(t);
        if(t->GetTauID("byTightDeepTau2017v2p1VSe") > 0.5 )  taus_et.push_back(t);
        if(t->GetTauID("byTightDeepTau2017v2p1VSmu") > 0.5 ) taus_mt.push_back(t);
        if(t->GetTauID("byTightDeepTau2017v2p1VSmu") > 0.5 && t->GetTauID("byVLooseDeepTau2017v2p1VSe") > 0.5 ) taus_sf.push_back(t);
        }
      }
      std::vector<std::pair<ic::GenJet *, ic::Tau *>> tt_matches_lead =  MatchByDR(lead_gen_tau,taus_tt,0.5,true,true); 
      std::vector<std::pair<ic::GenJet *, ic::Tau *>> tt_matches_sublead =  MatchByDR(sublead_gen_tau,taus_tt,0.5,true,true);

      std::vector<std::pair<ic::GenJet *, ic::Tau *>> mt_matches =  MatchByDR(lead_gen_tau,taus_mt,0.5,true,true); 
      std::vector<std::pair<ic::GenJet *, ic::Tau *>> et_matches =  MatchByDR(lead_gen_tau,taus_et,0.5,true,true); 
      std::vector<std::pair<ic::GenJet *, ic::Tau *>> sf_matches =  MatchByDR(lead_gen_tau,taus_sf,0.5,true,true); 
 
      if(tt_matches_lead.size()>0) tau_pt_1_tt_=tt_matches_lead[0].second->pt();
      if(tt_matches_sublead.size()>0) tau_pt_2_tt_=tt_matches_sublead[0].second->pt();
      if(mt_matches.size()>0) tau_pt_1_mt_=mt_matches[0].second->pt();
      if(et_matches.size()>0) tau_pt_1_et_=et_matches[0].second->pt();
      if(sf_matches.size()>0) tau_pt_1_sf_=sf_matches[0].second->pt();
 
    }

    cp_sign_1_ = 999;
    cp_sign_2_ = 999;
    cp_sign_3_ = 999;
    cp_sign_4_ = 999;
    aco_angle_1_ = 9999.;
    aco_angle_2_ = 9999.;
    aco_angle_3_ = 9999.;
    aco_angle_4_ = 9999.;
    cp_channel_=-1;
    Ediff_1_ = -9999.;
    Ediff_2_ = -9999.;  
    Pfrac_1_=-9999.; 
    Pfrac_2_=-9999.;

    pt_tt_ = (met.vector()+lep1.vector()+lep2.vector()).Pt();
    mass_ = (met.vector()+lep1.vector()+lep2.vector()).M();


    if(passed_){
      pt_1_  = lep1.vector().Pt();
      pt_2_  = lep2.vector().Pt();
      eta_1_ = lep1.vector().Rapidity();
      eta_2_ = lep2.vector().Rapidity();
      phi_1_ = lep1.vector().Phi();
      phi_2_ = lep2.vector().Phi();
      met_   = met.vector().Pt();
      pt_tt_ = (met.vector()+lep1.vector()+lep2.vector()).Pt();
      //mass_ = (met.vector()+lep1.vector()+lep2.vector()).M();
      //wtzpt_ = z_pt_weights_sm_.GetBinContent(z_pt_weights_sm_.FindBin(mass_,pt_tt_));
      m_vis_ = (lep1.vector()+lep2.vector()).M();
      if(!(channel_str_ == "zmm") && !(channel_str_ == "zee")){
        mt_1_ = MT(&lep1, &met);
        mt_2_ = MT(&lep2, &met);
        //double mt_lep = MT(&lep1, &lep1);
        //mt_tot_ = sqrt(mt_lep*mt_lep + mt_1_*mt_1_ + mt_2_*mt_2_);
      }

      ic::CompositeCandidate *ditau = new ic::CompositeCandidate();
      ditau->AddCandidate("lep1",&lep1);
      ditau->AddCandidate("lep2",&lep2);
      pzeta_ = PZeta(ditau, &met , 0.85);
    } else {
      pt_1_  = -9999;
      pt_2_  = -9999;
      eta_1_ = -9999;
      eta_2_ = -9999;
      phi_1_ = -9999;
      phi_2_ = -9999;
      met_   = -9999;
      //mass_= -9999;
      m_vis_ = -9999;
      mt_1_ = -9999;
      mt_2_ = -9999;
      pzeta_ = -9999;
    }
    
    std::vector<ic::GenJet*> filtered_jets;
    std::vector<ic::GenJet*> bjets;
   
    n_bpartons_=0;
 
    std::vector<GenParticle *> sel_bquarks;
    for (unsigned i=0; i < gen_particles.size(); ++i){
      std::vector<bool> status_flags = gen_particles[i]->statusFlags();
      unsigned id = abs(gen_particles[i]->pdgid());  
      //if(id == 5 && status_flags[FromHardProcess] && status_flags[IsLastCopy] && gen_particles[i]->vector().Pt()>0){
      if(id == 5 && status_flags[FromHardProcess] && status_flags[IsLastCopy] && gen_particles[i]->vector().Pt()>0) n_bpartons_++;
      if(id == 5 && status_flags[IsLastCopy] && gen_particles[i]->vector().Pt()>5){
        sel_bquarks.push_back(gen_particles[i]);
      }
    }
     
    for(unsigned i=0; i<gen_jets.size(); ++i){
      ic::GenJet *jet = gen_jets[i];
      double jetPt = jet->vector().Pt();
      double jetEta = std::fabs(jet->vector().Rapidity());
      if(jetPt > min_jet_pt_ && jetEta < max_jet_eta_) filtered_jets.push_back(jet); 
      if(jetPt > 20 && jetEta < 2.5){
          bool MatchedToB = false;
          for(unsigned j=0; j<sel_bquarks.size(); ++j) if(DRLessThan(std::make_pair(jet, sel_bquarks[j]),0.5)) MatchedToB = true;
          if(MatchedToB) bjets.push_back(jet); 
      }
    }
    
    for(unsigned i=0; i<bjets.size(); ++i){
      ic::GenJet *jet = bjets[i];
      bool MatchedToPrompt = false;
      for(unsigned j=0; j<higgs_products.size(); ++j){
        if(DRLessThan(std::make_pair(jet, &higgs_products[j]),0.5)) MatchedToPrompt = true;
      }
      //remove jets that are matched to Higgs decay products
      if(MatchedToPrompt) bjets.erase (bjets.begin()+i);
    }
   
    std::vector<ic::GenJet*> bjets_eta2p4 = bjets;
    std::vector<ic::GenJet*> bjets_25 = bjets;

    ic::erase_if(bjets_eta2p4,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));
    ic::erase_if(bjets_25,!boost::bind(MinPtMaxEta, _1, 25.0, 2.4));
 
    n_bjets_pt25_ = bjets_25.size();
    n_bjets_eta2p5_noscale_ = bjets.size();
    n_bjets_noscale_ = bjets_eta2p4.size();
    bpt_1_=-9999.;
    if(n_bjets_noscale_>0) bpt_1_ = bjets_eta2p4[0]->pt();

    if(bbtag_eff_ != nullptr) {
      for(unsigned i=0; i<bjets.size(); ++i){
        ic::GenJet *jet = bjets[i];
        bool fail = false;
        double eff = bbtag_eff_->GetBinContent(bbtag_eff_->GetXaxis()->FindBin(jet->pt()),bbtag_eff_->GetYaxis()->FindBin(fabs(jet->eta()))); 
        //double sf = reader_comb->eval_auto_bounds("central",BTagEntry::FLAV_B, jet->eta(), jet->pt());
        double sf=1.;
        //std::cout << eff << "    " << sf << std::endl;

        double r = rand->Uniform();   
        if(r>eff*sf) fail = true;
        //remove jets that do not pass btagging
        if(fail) bjets.erase(bjets.begin()+i);
      }
    }

    bjets_eta2p4 = bjets;
    ic::erase_if(bjets_eta2p4,!boost::bind(MinPtMaxEta, _1, 20.0, 2.4));

    n_bjets_ = bjets_eta2p4.size();
    n_bjets_eta2p5_ = bjets.size();
   
 
    for(unsigned i=0; i<filtered_jets.size(); ++i){
      ic::GenJet *jet = filtered_jets[i];
      bool MatchedToPrompt = false;
      for(unsigned j=0; j<higgs_products.size(); ++j){
        if(DRLessThan(std::make_pair(jet, &higgs_products[j]),0.5)) MatchedToPrompt = true;
      }
      //remove jets that are matched to Higgs decay products
      if(MatchedToPrompt) filtered_jets.erase (filtered_jets.begin()+i);
    }

    for(unsigned i=0; i<gen_jets.size(); ++i){
      ic::GenJet *jet = gen_jets[i];
      bool MatchedToPrompt = false;
      for(unsigned j=0; j<higgs_products.size(); ++j){
        if(DRLessThan(std::make_pair(jet, &higgs_products[j]),0.5)) MatchedToPrompt = true;
      }
      //remove jets that are matched to Higgs decay products
            if(MatchedToPrompt) gen_jets.erase (gen_jets.begin()+i);
    }

    largest_gen_mjj_=-9999;
    if(gen_jets.size()>1){
      for(unsigned i=0; i<gen_jets.size()-1;++i){
        for(unsigned j=i+1; j<gen_jets.size();++j){
          double mass = (gen_jets[i]->vector()+gen_jets[j]->vector()).M();
         if(mass>largest_gen_mjj_) largest_gen_mjj_ = mass;
         }
      }
    }

    n_jets_ = filtered_jets.size();
    jpt_1_       = -9999;
    jeta_1_      = -9999;
    jphi_1_      = -9999;
    jpt_2_       = -9999;
    jpt_3_       = -9999;
    jeta_2_      = -9999;
    jphi_2_      = -9999;
    mjj_         = -9999;
    jdeta_       = -9999;
    n_jetsingap_ = 9999;
    if(n_jets_ > 0){
      jpt_1_  = filtered_jets[0]->vector().Pt();
      jeta_1_ = filtered_jets[0]->vector().Rapidity();
      jphi_1_ = filtered_jets[0]->vector().Phi();
    }
    if(n_jets_ > 1){
      jpt_2_  = filtered_jets[1]->vector().Pt();
      jeta_2_ = filtered_jets[1]->vector().Rapidity();
      jphi_2_ = filtered_jets[1]->vector().Phi();
      mjj_   = (filtered_jets[0]->vector()+filtered_jets[1]->vector()).M();
      jdeta_  = std::fabs(filtered_jets[0]->vector().Rapidity() - filtered_jets[1]->vector().Rapidity());
      double max_jeta = std::max(jeta_1_,jeta_2_);
      double min_jeta = std::min(jeta_1_,jeta_2_);
      n_jetsingap_ = 0;
      for(unsigned i=2; i<n_jets_; ++i){
         double jeta_3 = filtered_jets[i]->vector().Rapidity();
         if(jeta_3 > min_jeta && jeta_3 < max_jeta) n_jetsingap_++;    
      }
    }
    if(n_jets_ > 2){
      jpt_3_  = filtered_jets[2]->vector().Pt();      
    }

    n_pjets_=0;    
    if(filtered_jets.size()>=2){
      if(filtered_jets[0]->eta() > filtered_jets[1]->eta()){
        sjdphi_ =  ROOT::Math::VectorUtil::DeltaPhi(filtered_jets[0]->vector(), filtered_jets[1]->vector());
      }
      else{
        sjdphi_ =  ROOT::Math::VectorUtil::DeltaPhi(filtered_jets[1]->vector(), filtered_jets[0]->vector());
      }
     // // sort jets higher and lower than higgs eta_1
     // std::vector<GenJet> jets_high;
     // std::vector<GenJet> jets_low;
     // for (unsigned i=0; i<filtered_jets.size(); ++i){
     //   if (filtered_jets[i]->eta() > higgs_eta) jets_high.push_back(filtered_jets[i]);
     //   else jets_low.push_back(filtered_jets[i]);
     // }
     // if(jets_low.size()>0) n_pjets_++;
     // if(jets_high.size()>0) n_pjets_++;
     // Candidate pseudo_jet_a;
     // Candidate pseudo_jet_b;
     // for (auto j : jets_low) pseudo_jet_a.set_vector(pseudo_jet_a.vector()+j.vector());
     // for (auto j : jets_high) pseudo_jet_b.set_vector(pseudo_jet_b.vector()+j.vector());
     // spjdphi_ =  ROOT::Math::VectorUtil::DeltaPhi(pseudo_jet_a.vector(),pseudo_jet_b.vector());
     // for (unsigned i=0; i<filtered_jets.size(); ++i){
     //   double dEta = std::fabs(higgs_eta - filtered_jets[i]->eta());
     //   if(i==0 || dEta<ysep_) ysep_ = dEta;
     // }
    } else {
      sjdphi_ = -9999;
      spjdphi_ = -9999;
    }

    /* wt_ps_down_ = 1.0; */
    /* wt_ps_up_ = 1.0; */
    /* if(n_jets_==0){ */
    /*   wt_ps_up_ =  ps_0jet_up_  .GetBinContent(ps_0jet_up_  .FindBin(HiggsPt_)); */  
    /*   wt_ps_down_ =  ps_0jet_down_.GetBinContent(ps_0jet_down_.FindBin(HiggsPt_)); */   
    /* } */
    /* if(n_jets_==1){ */ 
    /*   wt_ps_up_ =  ps_1jet_up_ .GetBinContent(ps_1jet_up_  .FindBin(HiggsPt_)); */  
    /*   wt_ps_down_ =  ps_1jet_down_.GetBinContent(ps_1jet_down_.FindBin(HiggsPt_)); */   
    /* } */
    /* if(n_jets_==2){ */
    /*   wt_ps_up_ =  ps_2jet_up_  .GetBinContent(ps_2jet_up_  .FindBin(HiggsPt_)); */  
    /*   wt_ps_down_ =  ps_2jet_down_.GetBinContent(ps_2jet_down_.FindBin(HiggsPt_)); */     
    /* } */
    /* if(n_jets_>2){ */
    /*   wt_ps_up_ =  ps_3jet_up_  .GetBinContent(ps_3jet_up_  .FindBin(HiggsPt_)); */
    /*   wt_ps_down_ =  ps_3jet_down_.GetBinContent(ps_3jet_down_.FindBin(HiggsPt_)); */
    /* } */
    //wt_ue_up_ = ue_up_  .GetBinContent(ue_up_  .FindBin(n_jets_)); 
    //wt_ue_down_ = ue_down_  .GetBinContent(ue_down_  .FindBin(n_jets_)); 

    if (n_jets_      == 0) wt_nnlops_ = ggh_ph_0jet_->Eval(std::min(pT_A_, 125.));
    else if (n_jets_ == 1) wt_nnlops_ = ggh_ph_1jet_->Eval(std::min(pT_A_, 625.));
    else if (n_jets_ == 2) wt_nnlops_ = ggh_ph_2jet_->Eval(std::min(pT_A_, 800.));
    else if (n_jets_ >= 3) wt_nnlops_ = ggh_ph_3jet_->Eval(std::min(pT_A_, 925.));
    else wt_nnlops_ = 1.;


    std::vector<PileupInfo *> puInfo;
    float true_int = -1;
    
    if(event->ExistsInTree("pileupInfo")){
      puInfo = event->GetPtrVec<PileupInfo>("pileupInfo");
        for (unsigned i = 0; i < puInfo.size(); ++i) {
          if (puInfo[i]->bunch_crossing() == 0)
            true_int = puInfo[i]->true_num_interactions();
        }
    
      n_pu_ = true_int;
    }
    
    auto args = std::vector<double>{pT};
    if(fns_["h_t_ratio"]){
      wt_ggh_t_ = fns_["h_t_ratio"]->eval(args.data());        
      wt_ggh_b_ = fns_["h_b_ratio"]->eval(args.data());
      wt_ggh_i_ = fns_["h_i_ratio"]->eval(args.data());
      wt_ggH_t_ = fns_["H_t_ratio"]->eval(args.data());
      wt_ggH_b_ = fns_["H_b_ratio"]->eval(args.data());
      wt_ggH_i_ = fns_["H_i_ratio"]->eval(args.data());
      wt_ggA_t_ = fns_["A_t_ratio"]->eval(args.data());
      wt_ggA_b_ = fns_["A_b_ratio"]->eval(args.data());
      wt_ggA_i_ = fns_["A_i_ratio"]->eval(args.data());
    }

    /*std::vector<std::pair<GenParticle*, GenParticle*>> rho_daughters;
    std::vector<std::pair<GenParticle*, GenParticle*>> prho_daughters;
    std::vector<std::pair<GenParticle*, GenParticle*>> l_daughters;
    std::vector<std::vector<GenParticle*>> a1_daughters;
    std::vector<std::vector<GenParticle*>> fakea1_daughters;
    std::vector<std::vector<GenParticle*>> charged_daughters;
    std::vector<std::vector<GenParticle*>> a1_2_daughters;
    std::vector<Candidate*> tau_rhos;
    std::vector<GenParticle*> pi_daughters;
    std::vector<GenParticle*> kaon_daughters;
    count_taus_=0;   
    for (unsigned i = 0; i < gen_particles.size(); ++i) {
      if(std::fabs(gen_particles[i]->pdgid()) == 15 && gen_particles[i]->statusFlags()[IsLastCopy] 
        && gen_particles[i]->statusFlags()[FromHardProcess]){

        GenParticle* tau = gen_particles[i];
        ++count_taus_; //expect two taus for tt channel
        bool foundRho = false;
        bool foundPi = false;
        bool foundLep = false;
        bool foundA1 = false;
        bool foundA1_2 = false;
        bool found3Pi1P0 = false;
        bool foundK = false;
        unsigned count_pi0 = 0;
        unsigned count_pi = 0;
        unsigned count_K = 0;
        unsigned count_hadr = 0;
        unsigned count_gamma =0;
        unsigned count_lep=0;
        unsigned count_charged=0;

        std::pair<GenParticle*,GenParticle*> rho = std::make_pair(new GenParticle(), new GenParticle());
        std::pair<GenParticle*,GenParticle*> prho = std::make_pair(new GenParticle(), new GenParticle());
        //std::pair<GenParticle*,GenParticle*> lep = std::make_pair(new GenParticle(), new GenParticle());
        std::vector<GenParticle*> charged_vec = {};
        std::vector<GenParticle*> a1 = {};
        std::vector<GenParticle*> fakea1 = {};
        std::vector<GenParticle*> a1_2 = {};
        GenParticle* pi = new GenParticle();
        GenParticle* kaon = new GenParticle();
        std::vector<int> daughters = gen_particles[i]->daughters();
        charged_vec.push_back(gen_particles[i]);
        for (unsigned d : daughters){
          unsigned daughter_id = fabs(gen_particles[d]->pdgid());
          int daughter_charge = gen_particles[d]->charge();
          if(daughter_id == 12 || daughter_id == 14 || daughter_id == 16) continue;
          if(daughter_id == 11 || daughter_id == 13) {
            ++count_lep;
            //lep.first = gen_particles[d];
            pi = gen_particles[d];
            //lep.second = tau;
            prho.first = pi;
            prho.second = tau;

            continue;
          }
          ++count_hadr;
          if(daughter_charge!=0) {
            ++count_charged;
            charged_vec.push_back(gen_particles[d]);
          }
          if (daughter_id == 22) {
            ++count_gamma;
            continue;
          }
          if (daughter_id == 211) {
            //for (auto d: daughters) 
            //    std::cout << "daughter Ids: " << gen_particles[d]->pdgid() << std::endl;
            ++count_pi;
            //std::cout << "pi charge  " << gen_particles[d]->vector() << std::endl;
            rho.first = gen_particles[d];
            pi = gen_particles[d];
            prho.first = pi;
            prho.second = tau;
            a1.push_back(gen_particles[d]);
            fakea1.push_back(gen_particles[d]);
            a1_2.push_back(gen_particles[d]);
            continue;
          }
          //if (daughter_id == 130 || daughter_id == 310 
          //        || daughter_id == 311 || daughter_id == 321) {
          if (daughter_id == 321) { // K^\pm
            //std::cout << "found 321 kaon" << std::endl;
            //for (auto d: daughters) 
            //    std::cout << "daughter Ids: " << gen_particles[d]->pdgid() << std::endl;
            ++count_K;
            kaon = gen_particles[d];
            prho.first = kaon;
            prho.second = tau;
            continue;
          }
          if (daughter_id == 111) {
            ++count_pi0;
            rho.second = gen_particles[d];
            a1_2.push_back(gen_particles[d]); 
            fakea1.push_back(gen_particles[d]);
            //std::cout << "pi zero " << gen_particles[d]->vector() << std::endl;
            continue;
          }
          if (daughter_id == 213) {
            //std::cout << "found 213 rho " << daughter_id << std::endl;
            //for (auto d: daughters) 
            //    std::cout << "daughter Ids: " << gen_particles[d]->pdgid() << std::endl;
            count_pi0 = 0;
            count_pi = 0;
            count_K = 0;
            count_hadr=0;
            count_gamma=0;
            std::vector<int> gdaughters = gen_particles[d]->daughters(); 
            for (unsigned g : gdaughters){
              unsigned gdaughter_id = fabs(gen_particles[g]->pdgid());
              ++count_hadr;
              if (gdaughter_id == 22) {
                ++count_gamma;
                continue;
              }
              if (gdaughter_id == 211) {
                ++count_pi;
                rho.first = gen_particles[g];
                continue;
              }
              if (gdaughter_id == 111) {
                ++count_pi0;
                rho.second = gen_particles[g];
                continue;
              }
            }
            foundRho = true;
            break;
          }
        }
        foundRho = foundRho || (count_hadr-count_gamma==2 && count_pi==1 && count_pi0==1 && count_K==0);
        foundPi = (count_hadr-count_gamma==1 && count_pi==1 && count_pi0==0 && count_K==0);
        foundLep = count_lep==1;
        foundA1 = foundA1 || (count_hadr-count_gamma==3 && count_pi==3 && count_pi0==0);
        found3Pi1P0 = found3Pi1P0 || (count_hadr-count_gamma==4 && count_pi==3 && count_pi0==1);
        foundA1_2 = foundA1_2 || (count_hadr-count_gamma==3 && count_pi==1 && count_pi0==2);
        foundK = count_K==1 && count_pi0==0 && count_pi==0;
        if (foundK) {
          kaon_daughters.push_back(kaon);
          prho_daughters.push_back(prho);
        }
        if(foundRho) {
          rho_daughters.push_back(rho);
          Candidate * tauvis = new Candidate();
          tauvis->set_vector(rho.first->vector()+rho.second->vector());
          tau_rhos.push_back(tauvis);
        }
        if((foundPi || foundLep)) {
          prho_daughters.push_back(prho);
          pi_daughters.push_back(pi);
        }
        //if(foundLep) {
        //  //l_daughters.push_back(lep);
        //  prho_daughters.push_back(prho);
        //  pi_daughters.push_back(pi);
        //}
        if(foundA1){
          a1_daughters.push_back(a1);
        }
        if(found3Pi1P0){
          fakea1_daughters.push_back(fakea1);
        }
        if(count_charged==3) {
          charged_daughters.push_back(charged_vec);
        }

        if(foundA1_2){
          a1_2_daughters.push_back(a1_2); 
        }
      }
    }

    TLorentzVector lvec1;
    TLorentzVector lvec2;
    TLorentzVector lvec3;
    TLorentzVector lvec4;
    TVector3 ip;
    if (rho_daughters.size()==1 && pi_daughters.size()==1){
        cp_channel_=2;
        cp_sign_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3());

        ip = GenIP(pi_daughters[0], prho_daughters[0].second);

        lvec1 = ConvertToLorentz(rho_daughters[0].second->vector()); //pi zero from rho
        lvec3 = ConvertToLorentz(rho_daughters[0].first->vector()); //pi charge from rho
        lvec4 = ConvertToLorentz(pi_daughters[0]->vector()); //pi charge from tau

        lvec2 = TLorentzVector(ip, 0.);
    }
    else if (rho_daughters.size()==1 && kaon_daughters.size()==1){
        cp_channel_=12;
        cp_sign_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3());

        ip = GenIP(kaon_daughters[0], prho_daughters[0].second);

        lvec1 = ConvertToLorentz(rho_daughters[0].second->vector()); //pi zero from rho
        lvec3 = ConvertToLorentz(rho_daughters[0].first->vector()); //pi charge from rho
        lvec4 = ConvertToLorentz(kaon_daughters[0]->vector()); //pi charge from tau

        lvec2 = TLorentzVector(ip, 0.);
    }
    
    else if (rho_daughters.size()==2){
        cp_channel_=3;  
        if(rho_daughters[1].first->charge()<0) {
          lvec1 = ConvertToLorentz(rho_daughters[0].second->vector());   
          lvec2 = ConvertToLorentz(rho_daughters[1].second->vector());
          lvec3 = ConvertToLorentz(rho_daughters[0].first->vector());   
          lvec4 = ConvertToLorentz(rho_daughters[1].first->vector());
        } else {
          lvec1 = ConvertToLorentz(rho_daughters[1].second->vector());
          lvec2 = ConvertToLorentz(rho_daughters[0].second->vector());
          lvec3 = ConvertToLorentz(rho_daughters[1].first->vector());
          lvec4 = ConvertToLorentz(rho_daughters[0].first->vector());
        }
        cp_sign_1_ = YRho(std::vector<GenParticle*>(
                    {rho_daughters[0].first, rho_daughters[0].second}),TVector3())*
                    YRho(std::vector<GenParticle*>(
                    {rho_daughters[1].first, rho_daughters[1].second}),TVector3());
    }
    if(cp_channel_==2 || cp_channel_==5 || cp_channel_==12){
      aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);    
      if (cp_sign_1_>0) {
        if (aco_angle_1_<M_PI) aco_angle_1_ += M_PI;
        else                   aco_angle_1_ -= M_PI;
      }  
    }

    if(cp_channel_==1 || cp_channel_==11) {
      aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
    }

    if(cp_channel_!=-1 && cp_channel_==3){
      aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);    

      if (cp_sign_1_>0) {
        if (aco_angle_1_<M_PI) aco_angle_1_ += M_PI;
        else                   aco_angle_1_ -= M_PI;
      }  
    }*/
    

    // =========================
    
    std::vector<std::pair<GenParticle*, GenParticle*>> rho_daughters;
    std::vector<std::pair<GenParticle*, GenParticle*>> prho_daughters;
    std::vector<std::pair<GenParticle*, GenParticle*>> l_daughters;
   
    std::vector<std::vector<GenParticle *>> pis = {{}, {}};
    std::vector<std::vector<GenParticle *>> pi0s = {{}, {}};

 

    std::vector<GenParticle*> tau_neutrinos;
    std::vector<std::vector<GenParticle*>> a1_daughters;
    std::vector<std::vector<GenParticle*>> fakea1_daughters;
    std::vector<std::vector<GenParticle*>> charged_daughters;
    std::vector<std::vector<GenParticle*>> a1_2_daughters;
    std::vector<Candidate*> tau_rhos;
    std::vector<GenParticle*> pi_daughters;
    std::vector<GenParticle*> kaon_daughters;
    std::vector<GenParticle*> undecayed_taus;
    unsigned count_taus=0;   
    for (unsigned i = 0; i < gen_particles.size(); ++i) {
      if(std::fabs(gen_particles[i]->pdgid()) == 15 && gen_particles[i]->statusFlags()[IsLastCopy]){
        GenParticle* tau = gen_particles[i];
        undecayed_taus.push_back(tau);
        bool foundRho = false;
        bool foundPi = false;
        bool foundLep = false;
        bool foundA1 = false;
        bool foundA1_2 = false;
        bool found3Pi1P0 = false;
        unsigned count_pi0 = 0;
        unsigned count_pi = 0;
        unsigned count_K = 0;
        unsigned count_hadr = 0;
        unsigned count_gamma =0;
        unsigned count_lep=0;
        unsigned count_charged=0;

        std::pair<GenParticle*,GenParticle*> rho = std::make_pair(new GenParticle(), new GenParticle());
        std::pair<GenParticle*,GenParticle*> prho = std::make_pair(new GenParticle(), new GenParticle());
        std::pair<GenParticle*,GenParticle*> lep = std::make_pair(new GenParticle(), new GenParticle());
        std::vector<GenParticle*> charged_vec = {};
        std::vector<GenParticle*> a1 = {};
        std::vector<GenParticle*> fakea1 = {};
        std::vector<GenParticle*> a1_2 = {};
        GenParticle* pi = new GenParticle();
        GenParticle* kaon = new GenParticle();
        std::vector<int> daughters = gen_particles[i]->daughters();
        charged_vec.push_back(gen_particles[i]);
        //std::cout << "tau: " << gen_particles[i] << std::endl;
        for (unsigned d : daughters){
          std::vector<bool> status_flags = gen_particles[d]->statusFlags(); 
          unsigned daughter_id = fabs(gen_particles[d]->pdgid());
          //std::cout << "ID in other function: " << gen_particles[d]->pdgid() << std::endl;
          int daughter_charge = gen_particles[d]->charge();
          if(daughter_id == 12 || daughter_id == 14 || daughter_id == 16) {
            if(status_flags[IsDirectPromptTauDecayProduct] && daughter_id==16) tau_neutrinos.push_back(gen_particles[d]);
            //tau_neutrinos.push_back(gen_particles[d]);
            //std::cout << "neutrino " << gen_particles[d]->pdgid() << std::endl;
            continue;
          }
          if(daughter_id == 11 || daughter_id == 13) {
            ++count_lep;
            lep.first = gen_particles[d];
            //pi = gen_particles[d];
            lep.second = tau;
            prho.first = pi;
            prho.second = tau;

            continue;
          }
          ++count_hadr;
          if(daughter_charge!=0) {
            ++count_charged;
            charged_vec.push_back(gen_particles[d]);
          }
          if (daughter_id == 22) {
            ++count_gamma;
            continue;
          }
          if (daughter_id == 211) {
            ++count_pi;
            //std::cout << "pi charge  " << gen_particles[d]->vector() << std::endl;
            rho.first = gen_particles[d];
            pi = gen_particles[d];
            pis[count_taus].push_back(gen_particles[d]);
            prho.first = pi;
            prho.second = tau;
            a1.push_back(gen_particles[d]);
            fakea1.push_back(gen_particles[d]);
            a1_2.push_back(gen_particles[d]);
            continue;
          }
          if (daughter_id == 321) { // K^\pm
            ++count_K;
            kaon = gen_particles[d];
            prho.first = kaon;
            prho.second = tau;
            continue;
          }
          if (daughter_id == 111) {
            ++count_pi0;
            pi0s[count_taus].push_back(gen_particles[d]);
            rho.second = gen_particles[d];
            a1_2.push_back(gen_particles[d]); 
            fakea1.push_back(gen_particles[d]);
            //std::cout << "pi zero " << gen_particles[d]->vector() << std::endl;
            continue;
          }
        }
        foundRho = foundRho || (count_hadr-count_gamma==2 && count_pi==1 && count_pi0==1);
        foundPi = (count_hadr-count_gamma==1 && count_pi==1 && count_pi0==0);
        foundLep = count_lep>0;
        foundA1 = (count_hadr-count_gamma==3 && count_pi==3 && count_pi0==0);
        found3Pi1P0 = (count_hadr-count_gamma==4 && count_pi==3 && count_pi0==1);
        foundA1_2 = (count_hadr-count_gamma==3 && count_pi==1 && count_pi0==2);
        if (count_hadr==1 && count_K==1 && count_pi0==0 && count_pi==0) {
          kaon_daughters.push_back(kaon);
          prho_daughters.push_back(prho);
        }
        if(foundRho&&!foundLep) {
          rho_daughters.push_back(rho);
          Candidate * tauvis = new Candidate();
          tauvis->set_vector(rho.first->vector()+rho.second->vector());
          tau_rhos.push_back(tauvis);
        }
        if(foundPi&&!foundLep) {
          prho_daughters.push_back(prho);
          pi_daughters.push_back(pi);
        }
        if(foundLep) {
          l_daughters.push_back(lep);
        //  prho_daughters.push_back(prho);
        }
        if(foundA1&&!foundLep){
          a1_daughters.push_back(a1);
        }
        if(found3Pi1P0&&!foundLep){
          fakea1_daughters.push_back(fakea1);
        }
        if(count_charged==3&&!foundLep) {
          charged_daughters.push_back(charged_vec);
        }

        if(foundA1_2&&!foundLep){
          a1_2_daughters.push_back(a1_2); 
        }
      count_taus++;
      if(count_taus==2) break; 
//      std::cout << count_taus << std::endl;
      }
    }
  
    spin_var_ = -9999.;
    mass_undecayed_taus_=-9999;
    if(taus.size()==2) {
      spin_var_ = SpinVar(undecayed_taus[0], undecayed_taus[1]);
      mass_undecayed_taus_ = (undecayed_taus[0]->vector()+undecayed_taus[1]->vector()).M();
    }

    fakea1_dR_ = -1;
    rho_dR_ = -1;

    if(rho_daughters.size()==1) {
      rho_dR_ = ROOT::Math::VectorUtil::DeltaR(rho_daughters[0].first->vector(),rho_daughters[0].second->vector());
    }


    if(fakea1_daughters.size()==1) {
      ROOT::Math::PtEtaPhiEVector charged_vec;
      ROOT::Math::PtEtaPhiEVector neutral_vec;
      for (auto p : fakea1_daughters[0]) {
        if (std::fabs(p->pdgid()) == 111) neutral_vec = p->vector();
        else charged_vec += p->vector();
      }
      fakea1_dR_ = ROOT::Math::VectorUtil::DeltaR(neutral_vec,charged_vec);
    }


    // fill vars for MVA training ntuple if option requested
    
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

    nu_px_1_=0.;
    nu_py_1_=0.;
    nu_pz_1_=0.;
    nu_E_1_=0.;
    nu_px_2_=0.;
    nu_py_2_=0.;
    nu_pz_2_=0.;
    nu_E_2_=0.;

    sv_x_1_=0.;
    sv_y_1_=0.;
    sv_z_1_=0.;
    sv_x_2_=0.;
    sv_y_2_=0.;
    sv_z_2_=0.;

    metx_   = met.vector().Px();
    mety_   = met.vector().Py();

    tau1_charge=0, tau2_charge=0;
    tauFlag_1_=-1, tauFlag_2_=-1;


    if(gen_vertices.size()>0 && tau_neutrinos.size()==2) {

      if (rho_daughters.size()>=1) {
        pi_px_1_ = rho_daughters[0].first->vector().Px();
        pi_py_1_ = rho_daughters[0].first->vector().Py();
        pi_pz_1_ = rho_daughters[0].first->vector().Pz();
        pi_E_1_ =  rho_daughters[0].first->vector().E();

        pi0_px_1_ = rho_daughters[0].second->vector().Px();
        pi0_py_1_ = rho_daughters[0].second->vector().Py();
        pi0_pz_1_ = rho_daughters[0].second->vector().Pz();
        pi0_E_1_ =  rho_daughters[0].second->vector().E();

        sv_x_1_ = rho_daughters[0].first->vtx().vx() - gen_vertices[0]->vx();
        sv_y_1_ = rho_daughters[0].first->vtx().vy() - gen_vertices[0]->vy();
        sv_z_1_ = rho_daughters[0].first->vtx().vz() - gen_vertices[0]->vz();
        tauFlag_1_=1;
        TVector3 gen_ip_1 = GetGenImpactParam (*(gen_vertices[0]), rho_daughters[0].first->vtx(), rho_daughters[0].first->vector());
        ip_x_1_ = gen_ip_1.X();
        ip_y_1_ = gen_ip_1.Y();
        ip_z_1_ = gen_ip_1.Z();

 
        tau1_charge = rho_daughters[0].first->charge();
        for (auto n: tau_neutrinos) {
          if(tau1_charge*n->pdgid()>=0) continue; 
          nu_px_1_ = n->vector().Px();
          nu_py_1_ = n->vector().Py();
          nu_pz_1_ = n->vector().Pz();
          nu_E_1_ = n->vector().E();
          break;
        }

        if (pi_daughters.size()>=1) {
          pi_px_2_ = pi_daughters[0]->vector().Px();
          pi_py_2_ = pi_daughters[0]->vector().Py();
          pi_pz_2_ = pi_daughters[0]->vector().Pz();
          pi_E_2_ =  pi_daughters[0]->vector().E();
          sv_x_2_ = pi_daughters[0]->vtx().vx() - gen_vertices[0]->vx();
          sv_y_2_ = pi_daughters[0]->vtx().vy() - gen_vertices[0]->vy();
          sv_z_2_ = pi_daughters[0]->vtx().vz() - gen_vertices[0]->vz();
          tauFlag_2_=0;

          TVector3 gen_ip_2 = GetGenImpactParam (*(gen_vertices[0]), pi_daughters[0]->vtx(), pi_daughters[0]->vector());
          ip_x_2_ = gen_ip_2.X();
          ip_y_2_ = gen_ip_2.Y();
          ip_z_2_ = gen_ip_2.Z();

          tau2_charge = pi_daughters[0]->charge();
          for (auto n: tau_neutrinos) {
            if(tau2_charge*n->pdgid()>=0) continue;               
            nu_px_2_ = n->vector().Px();
            nu_py_2_ = n->vector().Py();
            nu_pz_2_ = n->vector().Pz();
            nu_E_2_ = n->vector().E();
            break;
          }
        }

        if (a1_daughters.size()>=1) {
          std::vector<GenParticle*> a1_daughters_sorted_1 = SortA1Products(a1_daughters[0]);
          pi_px_2_ = a1_daughters_sorted_1[0]->vector().Px();
          pi_py_2_ = a1_daughters_sorted_1[0]->vector().Py();
          pi_pz_2_ = a1_daughters_sorted_1[0]->vector().Pz();
          pi_E_2_ =  a1_daughters_sorted_1[0]->vector().E();
          pi2_px_2_ = a1_daughters_sorted_1[1]->vector().Px();
          pi2_py_2_ = a1_daughters_sorted_1[1]->vector().Py();
          pi2_pz_2_ = a1_daughters_sorted_1[1]->vector().Pz();
          pi2_E_2_ =  a1_daughters_sorted_1[1]->vector().E();
          pi3_px_2_ = a1_daughters_sorted_1[2]->vector().Px();
          pi3_py_2_ = a1_daughters_sorted_1[2]->vector().Py();
          pi3_pz_2_ = a1_daughters_sorted_1[2]->vector().Pz();
          pi3_E_2_ =  a1_daughters_sorted_1[2]->vector().E();

          sv_x_2_ = a1_daughters_sorted_1[0]->vtx().vx() - gen_vertices[0]->vx();
          sv_y_2_ = a1_daughters_sorted_1[0]->vtx().vy() - gen_vertices[0]->vy();
          sv_z_2_ = a1_daughters_sorted_1[0]->vtx().vz() - gen_vertices[0]->vz();
          tauFlag_2_=10;

          tau2_charge = a1_daughters_sorted_1[0]->charge()+a1_daughters_sorted_1[1]->charge()+a1_daughters_sorted_1[2]->charge();
          for (auto n: tau_neutrinos) {
            if(tau2_charge*n->pdgid()>=0) continue;
            nu_px_2_ = n->vector().Px();
            nu_py_2_ = n->vector().Py();
            nu_pz_2_ = n->vector().Pz();
            nu_E_2_ = n->vector().E();
            //std::cout << (n->vector()+a1_daughters_sorted_1[0]->vector()+a1_daughters_sorted_1[1]->vector()+a1_daughters_sorted_1[2]->vector()).M() << std::endl;
            break;
          }

        }

        if (rho_daughters.size()>1) {
          pi_px_2_ = rho_daughters[1].first->vector().Px();
          pi_py_2_ = rho_daughters[1].first->vector().Py();
          pi_pz_2_ = rho_daughters[1].first->vector().Pz();
          pi_E_2_ =  rho_daughters[1].first->vector().E();

          pi0_px_2_ = rho_daughters[1].second->vector().Px();
          pi0_py_2_ = rho_daughters[1].second->vector().Py();
          pi0_pz_2_ = rho_daughters[1].second->vector().Pz();
          pi0_E_2_ =  rho_daughters[1].second->vector().E();

          sv_x_2_ = rho_daughters[1].first->vtx().vx() - gen_vertices[0]->vx();
          sv_y_2_ = rho_daughters[1].first->vtx().vy() - gen_vertices[0]->vy();
          sv_z_2_ = rho_daughters[1].first->vtx().vz() - gen_vertices[0]->vz();
          tauFlag_2_=1;

          TVector3 gen_ip_2 = GetGenImpactParam (*(gen_vertices[0]), rho_daughters[1].first->vtx(), rho_daughters[1].first->vector());
          ip_x_2_ = gen_ip_2.X();
          ip_y_2_ = gen_ip_2.Y();
          ip_z_2_ = gen_ip_2.Z();

          tau2_charge = rho_daughters[1].first->charge();
          for (auto n: tau_neutrinos) {
            if(tau2_charge*n->pdgid()>=0) continue;
            nu_px_2_ = n->vector().Px();
            nu_py_2_ = n->vector().Py();
            nu_pz_2_ = n->vector().Pz();
            nu_E_2_ = n->vector().E();
            //std::cout << (n->vector()+rho_daughters[1].first->vector()+rho_daughters[1].second->vector()).M() << std::endl;
            break;
          }
        }
      } else if (a1_daughters.size()>=1) {

          std::vector<GenParticle*> a1_daughters_sorted_1 = SortA1Products(a1_daughters[0]);
          pi_px_1_ = a1_daughters_sorted_1[0]->vector().Px();
          pi_py_1_ = a1_daughters_sorted_1[0]->vector().Py();
          pi_pz_1_ = a1_daughters_sorted_1[0]->vector().Pz();
          pi_E_1_ =  a1_daughters_sorted_1[0]->vector().E();
          pi2_px_1_ = a1_daughters_sorted_1[1]->vector().Px();
          pi2_py_1_ = a1_daughters_sorted_1[1]->vector().Py();
          pi2_pz_1_ = a1_daughters_sorted_1[1]->vector().Pz();
          pi2_E_1_ =  a1_daughters_sorted_1[1]->vector().E();
          pi3_px_1_ = a1_daughters_sorted_1[2]->vector().Px();
          pi3_py_1_ = a1_daughters_sorted_1[2]->vector().Py();
          pi3_pz_1_ = a1_daughters_sorted_1[2]->vector().Pz();
          pi3_E_1_ =  a1_daughters_sorted_1[2]->vector().E();

          sv_x_1_ = a1_daughters_sorted_1[0]->vtx().vx() - gen_vertices[0]->vx();
          sv_y_1_ = a1_daughters_sorted_1[0]->vtx().vy() - gen_vertices[0]->vy();
          sv_z_1_ = a1_daughters_sorted_1[0]->vtx().vz() - gen_vertices[0]->vz();
          tauFlag_1_=10;

          tau1_charge = a1_daughters_sorted_1[0]->charge()+a1_daughters_sorted_1[1]->charge()+a1_daughters_sorted_1[2]->charge();
          for (auto n: tau_neutrinos) {
            if(tau1_charge*n->pdgid()>=0) continue;
            nu_px_1_ = n->vector().Px();
            nu_py_1_ = n->vector().Py();
            nu_pz_1_ = n->vector().Pz();
            nu_E_1_ = n->vector().E();
            //std::cout << (n->vector()+a1_daughters_sorted_1[0]->vector()+a1_daughters_sorted_1[1]->vector()+a1_daughters_sorted_1[2]->vector()).M() << std::endl;
            break;
          }

          if (pi_daughters.size()>=1) {
            pi_px_2_ = pi_daughters[0]->vector().Px();
            pi_py_2_ = pi_daughters[0]->vector().Py();
            pi_pz_2_ = pi_daughters[0]->vector().Pz();
            pi_E_2_ =  pi_daughters[0]->vector().E();
            sv_x_2_ = pi_daughters[0]->vtx().vx() - gen_vertices[0]->vx();
            sv_y_2_ = pi_daughters[0]->vtx().vy() - gen_vertices[0]->vy();
            sv_z_2_ = pi_daughters[0]->vtx().vz() - gen_vertices[0]->vz();
            tauFlag_2_=0;

            TVector3 gen_ip_2 = GetGenImpactParam (*(gen_vertices[0]), pi_daughters[0]->vtx(), pi_daughters[0]->vector());
            ip_x_2_ = gen_ip_2.X();
            ip_y_2_ = gen_ip_2.Y();
            ip_z_2_ = gen_ip_2.Z();

            tau2_charge = pi_daughters[0]->charge();
            for (auto n: tau_neutrinos) {
              if(tau2_charge*n->pdgid()>=0) continue;
              nu_px_2_ = n->vector().Px();
              nu_py_2_ = n->vector().Py();
              nu_pz_2_ = n->vector().Pz();
              nu_E_2_ = n->vector().E();
              //std::cout << (n->vector()+pi_daughters[0]->vector()).M() << std::endl;
              break;
            }
          }

          if (a1_daughters.size()>1) {
            std::vector<GenParticle*> a1_daughters_sorted_2 = SortA1Products(a1_daughters[1]);
            pi_px_2_ = a1_daughters_sorted_2[0]->vector().Px();
            pi_py_2_ = a1_daughters_sorted_2[0]->vector().Py();
            pi_pz_2_ = a1_daughters_sorted_2[0]->vector().Pz();
            pi_E_2_ =  a1_daughters_sorted_2[0]->vector().E();
            pi2_px_2_ = a1_daughters_sorted_2[1]->vector().Px();
            pi2_py_2_ = a1_daughters_sorted_2[1]->vector().Py();
            pi2_pz_2_ = a1_daughters_sorted_2[1]->vector().Pz();
            pi2_E_2_ =  a1_daughters_sorted_2[1]->vector().E();
            pi3_px_2_ = a1_daughters_sorted_2[2]->vector().Px();
            pi3_py_2_ = a1_daughters_sorted_2[2]->vector().Py();
            pi3_pz_2_ = a1_daughters_sorted_2[2]->vector().Pz();
            pi3_E_2_ =  a1_daughters_sorted_2[2]->vector().E();

            sv_x_2_ = a1_daughters_sorted_2[0]->vtx().vx() - gen_vertices[0]->vx();
            sv_y_2_ = a1_daughters_sorted_2[0]->vtx().vy() - gen_vertices[0]->vy();
            sv_z_2_ = a1_daughters_sorted_2[0]->vtx().vz() - gen_vertices[0]->vz();
            tauFlag_2_=10;

            tau2_charge = a1_daughters_sorted_2[0]->charge()+a1_daughters_sorted_2[1]->charge()+a1_daughters_sorted_2[2]->charge();
            for (auto n: tau_neutrinos) {
              if(tau2_charge*n->pdgid()>=0) continue;
              nu_px_2_ = n->vector().Px();
              nu_py_2_ = n->vector().Py();
              nu_pz_2_ = n->vector().Pz();
              nu_E_2_ = n->vector().E();
              //std::cout << (n->vector()+a1_daughters_sorted_2[0]->vector()+a1_daughters_sorted_2[1]->vector()+a1_daughters_sorted_2[2]->vector()).M() << std::endl;
              break;
            }

          }
      } else if (pi_daughters.size()>=1) {

          pi_px_1_ = pi_daughters[0]->vector().Px();
          pi_py_1_ = pi_daughters[0]->vector().Py();
          pi_pz_1_ = pi_daughters[0]->vector().Pz();
          pi_E_1_ =  pi_daughters[0]->vector().E();
          tauFlag_1_=0;

          tau1_charge = pi_daughters[0]->charge();
          for (auto n: tau_neutrinos) {
            if(tau1_charge*n->pdgid()>=0) continue;
            nu_px_1_ = n->vector().Px();
            nu_py_1_ = n->vector().Py();
            nu_pz_1_ = n->vector().Pz();
            nu_E_1_ = n->vector().E();
            //std::cout << (n->vector()+pi_daughters[0]->vector()).M() << std::endl;
            break;
          }

          sv_x_1_ = pi_daughters[0]->vtx().vx() - gen_vertices[0]->vx();
          sv_y_1_ = pi_daughters[0]->vtx().vy() - gen_vertices[0]->vy();
          sv_z_1_ = pi_daughters[0]->vtx().vz() - gen_vertices[0]->vz();

          TVector3 gen_ip_1 = GetGenImpactParam (*(gen_vertices[0]), pi_daughters[0]->vtx(), pi_daughters[0]->vector());
          ip_x_1_ = gen_ip_1.X();
          ip_y_1_ = gen_ip_1.Y();
          ip_z_1_ = gen_ip_1.Z();

          if (pi_daughters.size()>1) {
            pi_px_2_ = pi_daughters[1]->vector().Px();
            pi_py_2_ = pi_daughters[1]->vector().Py();
            pi_pz_2_ = pi_daughters[1]->vector().Pz();
            pi_E_2_ =  pi_daughters[1]->vector().E();

            tau2_charge = pi_daughters[1]->charge();
            for (auto n: tau_neutrinos) {
              if(tau2_charge*n->pdgid()>=0) continue;
              nu_px_2_ = n->vector().Px();
              nu_py_2_ = n->vector().Py();
              nu_pz_2_ = n->vector().Pz();
              nu_E_2_ = n->vector().E();
              //std::cout << (n->vector()+pi_daughters[1]->vector()).M() << std::endl;
              break;
            }

            sv_x_2_ = pi_daughters[1]->vtx().vx() - gen_vertices[0]->vx();
            sv_y_2_ = pi_daughters[1]->vtx().vy() - gen_vertices[0]->vy();
            sv_z_2_ = pi_daughters[1]->vtx().vz() - gen_vertices[0]->vz();
            tauFlag_2_=0;

            TVector3 gen_ip_2 = GetGenImpactParam (*(gen_vertices[0]), pi_daughters[1]->vtx(), pi_daughters[1]->vector());
            ip_x_2_ = gen_ip_2.X();
            ip_y_2_ = gen_ip_2.Y();
            ip_z_2_ = gen_ip_2.Z();

          }

      }
    }
    

    TLorentzVector lvec1;
    TLorentzVector lvec2;
    TLorentzVector lvec3;
    TLorentzVector lvec4;

    if(a1_2_daughters.size()==1 && rho_daughters.size()==1){
      cp_channel_=6;
      if (a1_2_daughters[0][1]->charge()!=0) {
        auto temp = a1_2_daughters[0][1];
        a1_2_daughters[0][1] = a1_2_daughters[0][0];
        a1_2_daughters[0][0] = temp;
      }
      if (a1_2_daughters[0][2]->charge()!=0) {
        auto temp = a1_2_daughters[0][2];
        a1_2_daughters[0][2] = a1_2_daughters[0][0];
        a1_2_daughters[0][0] = temp;
      }
      mass1_ = (a1_2_daughters[0][0]->vector() + a1_2_daughters[0][1]->vector()).M();
      mass2_ = (a1_2_daughters[0][0]->vector() + a1_2_daughters[0][2]->vector()).M();
      if(std::fabs(mass2_-0.775) < std::fabs(mass1_-0.775)) {
        double temp_mass = mass1_;
        mass1_ = mass2_;
        mass2_ = temp_mass;
        auto temp = a1_2_daughters[0][1];
        a1_2_daughters[0][1] = a1_2_daughters[0][2];
        a1_2_daughters[0][2] = temp;
      }

      std::pair<GenParticle*,GenParticle*> p_rho = std::make_pair(new GenParticle(), new GenParticle());
      std::pair<GenParticle*,GenParticle*> p_rho_1 = std::make_pair(new GenParticle(), new GenParticle());
      p_rho_1.first = a1_2_daughters[0][0];
      p_rho_1.second = a1_2_daughters[0][1];
      p_rho.first = a1_2_daughters[0][0];
      p_rho.second = a1_2_daughters[0][1];
      p_rho.second->set_vector(a1_2_daughters[0][1]->vector()+a1_2_daughters[0][2]->vector());

      pt_rho_pi_ = a1_2_daughters[0][0]->pt();
      pt_rho_pi0_ = a1_2_daughters[0][1]->pt();
      pt_other_pi0_ = a1_2_daughters[0][2]->pt();

      GenParticle* rho_1  = new GenParticle();
      GenParticle* rho_2  = new GenParticle();
     
      //rho_1->set_vector(a1_2_daughters[0][0]->vector()+a1_2_daughters[0][1]->vector());

      lvec1 = ConvertToLorentz(rho_daughters[0].second->vector());
      lvec2 = ConvertToLorentz(p_rho_1.second->vector());
      lvec3 = ConvertToLorentz(rho_daughters[0].first->vector());
      lvec4 = ConvertToLorentz(p_rho_1.first->vector());
      cp_sign_1_ = YRho(std::vector<GenParticle*>(
                    {rho_daughters[0].first, rho_daughters[0].second}),TVector3())*
                    YRho(std::vector<GenParticle*>(
                    {p_rho_1.first, p_rho_1.second}),TVector3());

      aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);

      if (cp_sign_1_<0) {
        if (aco_angle_1_<M_PI) aco_angle_1_ += M_PI;
        else                   aco_angle_1_ -= M_PI;
      }

      lvec1 = ConvertToLorentz(rho_daughters[0].second->vector());
      lvec2 = ConvertToLorentz(p_rho.second->vector());
      lvec3 = ConvertToLorentz(rho_daughters[0].first->vector());
      lvec4 = ConvertToLorentz(p_rho.first->vector());
      cp_sign_2_ = YRho(std::vector<GenParticle*>(
                    {rho_daughters[0].first, rho_daughters[0].second}),TVector3())*
                    YRho(std::vector<GenParticle*>(
                    {p_rho.first, p_rho.second}),TVector3());

      aco_angle_2_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);

      if (cp_sign_2_<0) {
        if (aco_angle_2_<M_PI) aco_angle_2_ += M_PI;
        else                   aco_angle_2_ -= M_PI;
      }



      //aco_angle_2_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({rho_1,a1_2_daughters[0][2]}));
      //cp_sign_2_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YA1(std::vector<GenParticle*>({rho_1, a1_2_daughters[0][2]}),TVector3());
      ////std::cout << aco_angle_2_ << std::endl;

      //aco_angle_1_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_2_daughters[0][0],a1_2_daughters[0][1]}));
      //aco_angle_3_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_2_daughters[0][0],a1_2_daughters[0][2]}));
      ////std::cout << a1_2_daughters[0][0]->charge() << "    " << a1_2_daughters[0][1]->charge() << "    " << a1_2_daughters[0][2]->charge() << "    " << aco_angle_1_ << std::endl;
      //aco_angle_5_ = aco_angle_1_;
      //aco_angle_6_ = aco_angle_3_;

      //GenParticle* rho_1  = new GenParticle();
      //GenParticle* rho_2  = new GenParticle();
      //rho_1->set_vector(a1_2_daughters[0][0]->vector()+a1_2_daughters[0][1]->vector());
      //rho_2->set_vector(a1_2_daughters[0][0]->vector()+a1_2_daughters[0][2]->vector());
      //aco_angle_2_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({rho_1,a1_2_daughters[0][2]}));
      //aco_angle_4_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({rho_2,a1_2_daughters[0][1]}));


      //cp_sign_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YRho(std::vector<GenParticle*>({a1_2_daughters[0][0], a1_2_daughters[0][1]}),TVector3());
      //cp_sign_3_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YRho(std::vector<GenParticle*>({a1_2_daughters[0][0], a1_2_daughters[0][2]}),TVector3());

      //cp_sign_2_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YA1(std::vector<GenParticle*>({rho_1, a1_2_daughters[0][2]}),TVector3());
      //cp_sign_4_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YA1(std::vector<GenParticle*>({rho_2, a1_2_daughters[0][1]}),TVector3()); 

      //mass1_ = (a1_2_daughters[0][0]->vector() + a1_2_daughters[0][1]->vector()).M();
      //mass2_ = (a1_2_daughters[0][0]->vector() + a1_2_daughters[0][2]->vector()).M();


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

    
      a1_daughters = a1_2_daughters;

      dR1_ = ROOT::Math::VectorUtil::DeltaR(a1_daughters[0][0]->vector(), a1_daughters[0][1]->vector());
      if(rand_<0.5) dR2_ = ROOT::Math::VectorUtil::DeltaR(a1_daughters[0][0]->vector(), a1_daughters[0][2]->vector());
      else dR2_ = ROOT::Math::VectorUtil::DeltaR(a1_daughters[0][1]->vector(), a1_daughters[0][2]->vector()); // pick this randomly so that it more closly follows what would happen for 3 prong channel

      // designed to mimic what is one for the 3 prong case where each of neutral particle are assigned a chare at random
      if(rand_<0.5){
        a1_daughters[0][0]->set_charge(-1);
        a1_daughters[0][1]->set_charge(1);
      } else {
        a1_daughters[0][0]->set_charge(1);
        a1_daughters[0][1]->set_charge(-1);
      }
      a1_daughters[0] = SortA1Products(a1_daughters[0]);

      //mass1_ = (a1_daughters[0][0]->vector() + a1_daughters[0][1]->vector()).M();
      //mass2_ = (a1_daughters[0][0]->vector() + a1_daughters[0][2]->vector()).M();

      //std::cout << a1_daughters[0][0]->charge() << "    " << a1_daughters[0][1]->charge() << "    " << a1_daughters[0][2]->charge() << "    " << mass1_ << "    " << mass2_ << std::endl;

      p_aco_angle_1_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][0],a1_daughters[0][1]}));
      p_aco_angle_3_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][0],a1_daughters[0][2]}));

      rho_1  = new GenParticle();
      rho_2  = new GenParticle();
      rho_1->set_vector(a1_daughters[0][0]->vector()+a1_daughters[0][1]->vector());
      rho_2->set_vector(a1_daughters[0][0]->vector()+a1_daughters[0][2]->vector());
      p_aco_angle_2_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({rho_1,a1_daughters[0][2]}));
      p_aco_angle_4_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({rho_2,a1_daughters[0][1]}));

      double cp_sign_1 = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YRho(std::vector<GenParticle*>({a1_daughters[0][0], a1_daughters[0][1]}),TVector3());
      double cp_sign_3 = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YRho(std::vector<GenParticle*>({a1_daughters[0][0], a1_daughters[0][2]}),TVector3());

      double cp_sign_2 = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YA1(std::vector<GenParticle*>({rho_1, a1_daughters[0][2]}),TVector3());
      double cp_sign_4 = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YA1(std::vector<GenParticle*>({rho_2, a1_daughters[0][1]}),TVector3());

      if (cp_sign_1<0) {
        if (p_aco_angle_1_<M_PI) p_aco_angle_1_ += M_PI;
        else                   p_aco_angle_1_ -= M_PI;
      }
      if (cp_sign_2<0) {
        if (p_aco_angle_2_<M_PI) p_aco_angle_2_ += M_PI;
        else                   p_aco_angle_2_ -= M_PI;
      }
      if (cp_sign_3<0) {
        if (p_aco_angle_3_<M_PI) p_aco_angle_3_ += M_PI;
        else                   p_aco_angle_3_ -= M_PI;
      }
      if (cp_sign_4<0) {
        if (p_aco_angle_4_<M_PI) p_aco_angle_4_ += M_PI;
        else                   p_aco_angle_4_ -= M_PI;
      }

      a1_daughters = {};

    }

    if(a1_daughters.size()==1){
      mass1_ = (a1_daughters[0][0]->vector() + a1_daughters[0][1]->vector()).M();
      mass2_ = (a1_daughters[0][0]->vector() + a1_daughters[0][2]->vector()).M();
      a1_mass_ = (a1_daughters[0][0]->vector() + a1_daughters[0][1]->vector()+a1_daughters[0][2]->vector()).M();

    }
    if(a1_daughters.size()==1 && rho_daughters.size()==1){
      cp_channel_=4;
      cp_sign_1_ = 999;
      cp_sign_2_ = 999;
      cp_sign_3_ = 999;
      cp_sign_4_ = 999;
      aco_angle_1_ = 9999.;
      aco_angle_2_ = 9999.;
      aco_angle_3_ = 9999.;
      aco_angle_4_ = 9999.;
      a1_daughters[0] = SortA1Products(a1_daughters[0]);
      mass1_ = (a1_daughters[0][0]->vector() + a1_daughters[0][1]->vector()).M();
      mass2_ = (a1_daughters[0][0]->vector() + a1_daughters[0][2]->vector()).M();

      aco_angle_1_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][0],a1_daughters[0][1]}));
      aco_angle_3_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][0],a1_daughters[0][2]}));
      //aco_angle_1_ = AcoplanarityAngle2(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][0],a1_daughters[0][1]}));
      //aco_angle_3_ = AcoplanarityAngle2(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][0],a1_daughters[0][2]}));


      GenParticle* rho_1  = new GenParticle();
      GenParticle* rho_2  = new GenParticle();
      rho_1->set_vector(a1_daughters[0][0]->vector()+a1_daughters[0][1]->vector());
      rho_2->set_vector(a1_daughters[0][0]->vector()+a1_daughters[0][2]->vector());
      aco_angle_2_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({rho_1,a1_daughters[0][2]}));
      aco_angle_4_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({rho_2,a1_daughters[0][1]}));

      //aco_angle_2_ = AcoplanarityAngle2(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][2],rho_1}));
      //aco_angle_4_ = AcoplanarityAngle2(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][1],rho_1}));

      cp_sign_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YRho(std::vector<GenParticle*>({a1_daughters[0][0], a1_daughters[0][1]}),TVector3());
      cp_sign_3_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YRho(std::vector<GenParticle*>({a1_daughters[0][0], a1_daughters[0][2]}),TVector3()); 

      cp_sign_2_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YA1(std::vector<GenParticle*>({rho_1, a1_daughters[0][2]}),TVector3());
      cp_sign_4_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YA1(std::vector<GenParticle*>({rho_2, a1_daughters[0][1]}),TVector3());

      //aco_angle_5_ = AcoplanarityAngle2(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][1],a1_daughters[0][2]}));
      aco_angle_5_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][1],a1_daughters[0][2]}));
      cp_sign_5_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YRho(std::vector<GenParticle*>({a1_daughters[0][1], a1_daughters[0][2]}),TVector3());
      GenParticle* rho_3  = new GenParticle();
      rho_3->set_vector(a1_daughters[0][1]->vector()+a1_daughters[0][2]->vector());
      //aco_angle_6_ = AcoplanarityAngle2(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({a1_daughters[0][0],rho_3}));
      aco_angle_6_ = AcoplanarityAngle(std::vector<GenParticle*> ({rho_daughters[0].first,rho_daughters[0].second}), std::vector<GenParticle*> ({rho_3,a1_daughters[0][0]}));
      cp_sign_6_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3())*YA1(std::vector<GenParticle*>({rho_3, a1_daughters[0][0]}),TVector3());
   
      y_1_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3());
      y_2_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3());
      y_3_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3());
      y_4_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3());

      y_1_2_ = YRho(std::vector<GenParticle*>({a1_daughters[0][0], a1_daughters[0][1]}),TVector3());
      y_2_2_ = YRho(std::vector<GenParticle*>({a1_daughters[0][0], a1_daughters[0][2]}),TVector3());
      y_3_2_ = YA1(std::vector<GenParticle*>({rho_1, a1_daughters[0][2]}),TVector3());
      y_4_2_ = YA1(std::vector<GenParticle*>({rho_2, a1_daughters[0][1]}),TVector3());

      a1_mass_ = (a1_daughters[0][0]->vector() + a1_daughters[0][1]->vector()+a1_daughters[0][2]->vector()).M();
      E_1_1_ = rho_daughters[0].first->energy();
      E_2_1_ = rho_daughters[0].second->energy();
      E_1_2_ = a1_daughters[0][0]->energy();
      E_2_2_ = a1_daughters[0][1]->energy();
      E_3_2_ = a1_daughters[0][2]->energy();


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
      if (cp_sign_5_<0) {
        if (aco_angle_5_<M_PI) aco_angle_5_ += M_PI;
        else                   aco_angle_5_ -= M_PI;
      }
      if (cp_sign_6_<0) {
        if (aco_angle_6_<M_PI) aco_angle_6_ += M_PI;
        else                   aco_angle_6_ -= M_PI;
      }
  
    }
 
    if(prho_daughters.size()==1){
      Pfrac_1_ = prho_daughters[0].first->vector().P()/prho_daughters[0].second->vector().P();
    }
    if(prho_daughters.size()==2){
      Pfrac_2_ = prho_daughters[1].first->vector().P()/prho_daughters[1].second->vector().P();
    }

    if(prho_daughters.size()==2){
      //cp_channel_=1;
    }

    if (pi_daughters.size()==2 && prho_daughters.size()==2){
        cp_channel_=1;

        TLorentzVector pvtosv1(
                pi_daughters[0]->vtx().vx() - prho_daughters[0].second->vtx().vx(),
                pi_daughters[0]->vtx().vy() - prho_daughters[0].second->vtx().vy(),
                pi_daughters[0]->vtx().vz() - prho_daughters[0].second->vtx().vz(),
                0.);

        TLorentzVector pvtosv2(
                pi_daughters[1]->vtx().vx() - prho_daughters[1].second->vtx().vx(),
                pi_daughters[1]->vtx().vy() - prho_daughters[1].second->vtx().vy(),
                pi_daughters[1]->vtx().vz() - prho_daughters[1].second->vtx().vz(),
                0.);

        if(pi_daughters[0]->charge()<0) {
          lvec4 = ConvertToLorentz(pi_daughters[1]->vector()); 
          lvec3 = ConvertToLorentz(pi_daughters[0]->vector()); 

          TVector3 ip1 = (pvtosv1.Vect() - pvtosv1.Vect().Dot(lvec3.Vect().Unit())*lvec3.Vect().Unit()).Unit();
          lvec1 = TLorentzVector(ip1, 0.);

          TVector3 ip2 = (pvtosv2.Vect() - pvtosv2.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
          lvec2 = TLorentzVector(ip2, 0.);
        }
    }

    if (rho_daughters.size()==1 && pi_daughters.size()==1){
        cp_channel_=2;
        cp_sign_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3());

        pi_px_1_ = pi_daughters[0]->vector().Px();
        pi_py_1_ = pi_daughters[0]->vector().Py();
        pi_pz_1_ = pi_daughters[0]->vector().Pz();
        pi_E_1_ =  pi_daughters[0]->vector().E();

        pi_px_2_ = rho_daughters[0].first->vector().Px();
        pi_py_2_ = rho_daughters[0].first->vector().Py();
        pi_pz_2_ = rho_daughters[0].first->vector().Pz();
        pi_E_2_ =  rho_daughters[0].first->vector().E();
        pi0_px_2_ = rho_daughters[0].second->vector().Px();
        pi0_py_2_ = rho_daughters[0].second->vector().Py();
        pi0_pz_2_ = rho_daughters[0].second->vector().Pz();
        pi0_E_2_ =  rho_daughters[0].second->vector().E();

        for (unsigned i = 0; i < gen_vertices.size(); i++) {
          gen_pvx_ = gen_vertices[0]->vx();
          gen_pvy_ = gen_vertices[0]->vy();
          gen_pvz_ = gen_vertices[0]->vz();
        }
        TLorentzVector pvtosv(
                pi_daughters[0]->vtx().vx() - prho_daughters[0].second->vtx().vx(), 
                pi_daughters[0]->vtx().vy() - prho_daughters[0].second->vtx().vy(), 
                pi_daughters[0]->vtx().vz() - prho_daughters[0].second->vtx().vz(), 
                0.);

        lvec1 = ConvertToLorentz(rho_daughters[0].second->vector()); //pi zero from rho
        lvec3 = ConvertToLorentz(rho_daughters[0].first->vector()); //pi charge from rho
        lvec4 = ConvertToLorentz(pi_daughters[0]->vector()); //pi charge from tau

        TVector3 ip = (pvtosv.Vect() - pvtosv.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
        lvec2 = TLorentzVector(ip, 0.);
    }
    else if (rho_daughters.size()==1 && kaon_daughters.size()==1){
        cp_channel_=12;
        cp_sign_1_ = YRho(std::vector<GenParticle*>({rho_daughters[0].first, rho_daughters[0].second}),TVector3());

        for (unsigned i = 0; i < gen_vertices.size(); i++) {
          gen_pvx_ = gen_vertices[0]->vx();
          gen_pvy_ = gen_vertices[0]->vy();
          gen_pvz_ = gen_vertices[0]->vz();
        }
        TLorentzVector pvtosv(
                kaon_daughters[0]->vtx().vx() - prho_daughters[0].second->vtx().vx(), 
                kaon_daughters[0]->vtx().vy() - prho_daughters[0].second->vtx().vy(), 
                kaon_daughters[0]->vtx().vz() - prho_daughters[0].second->vtx().vz(), 
                0.);

        lvec1 = ConvertToLorentz(rho_daughters[0].second->vector()); //pi zero from rho
        lvec3 = ConvertToLorentz(rho_daughters[0].first->vector()); //pi charge from rho
        lvec4 = ConvertToLorentz(kaon_daughters[0]->vector()); //pi charge from tau

        TVector3 ip = (pvtosv.Vect() - pvtosv.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
        lvec2 = TLorentzVector(ip, 0.);
    }

    else if (rho_daughters.size()==1 && l_daughters.size()==1){
        cp_channel_=5;
        std::vector<ic::Vertex*> vertex_vec = {};
        if(event->ExistsInTree("vertices")) vertex_vec = event->GetPtrVec<ic::Vertex>("vertices"); // reco
        for (unsigned i = 0; i < vertex_vec.size(); i++) {
          reco_pvx_ = vertex_vec[0]->vx();
          reco_pvy_ = vertex_vec[0]->vy();
          reco_pvz_ = vertex_vec[0]->vz();
        }
        for (unsigned i = 0; i < gen_vertices.size(); i++) {
          gen_pvx_ = gen_vertices[0]->vx();
          gen_pvy_ = gen_vertices[0]->vy();
          gen_pvz_ = gen_vertices[0]->vz();
        }
        TLorentzVector pvtosv(
                l_daughters[0].first->vtx().vx() - l_daughters[0].second->vtx().vx(), 
                l_daughters[0].first->vtx().vy() - l_daughters[0].second->vtx().vy(), 
                l_daughters[0].first->vtx().vz() - l_daughters[0].second->vtx().vz(), 
                0.);

        if(l_daughters[0].first->charge()<0) {

          lvec1 = ConvertToLorentz(rho_daughters[0].second->vector());
          lvec3 = ConvertToLorentz(rho_daughters[0].first->vector());
          lvec4 = ConvertToLorentz(l_daughters[0].first->vector());
          TVector3 ip = (pvtosv.Vect() - pvtosv.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
          lvec2 = TLorentzVector(ip, 0.);
        } else {
          lvec2 = ConvertToLorentz(rho_daughters[0].second->vector());
          lvec4 = ConvertToLorentz(rho_daughters[0].first->vector());
          lvec3 = ConvertToLorentz(l_daughters[0].first->vector());
          TVector3 ip = (pvtosv.Vect() - pvtosv.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
          lvec1 = TLorentzVector(ip, 0.);
        }
        cp_sign_1_ = YRho(std::vector<GenParticle*>(
                    {rho_daughters[0].first, rho_daughters[0].second}),TVector3());
    }

    if (a1_daughters.size()==1 && l_daughters.size()==1){

      a1_daughters[0] = SortA1Products(a1_daughters[0]);

        cp_channel_=51;
        std::vector<ic::Vertex*> vertex_vec = {};
        if(event->ExistsInTree("vertices")) vertex_vec = event->GetPtrVec<ic::Vertex>("vertices"); // reco
        for (unsigned i = 0; i < vertex_vec.size(); i++) {
          reco_pvx_ = vertex_vec[0]->vx();
          reco_pvy_ = vertex_vec[0]->vy();
          reco_pvz_ = vertex_vec[0]->vz();
        }
        for (unsigned i = 0; i < gen_vertices.size(); i++) {
          gen_pvx_ = gen_vertices[0]->vx();
          gen_pvy_ = gen_vertices[0]->vy();
          gen_pvz_ = gen_vertices[0]->vz();
        }
        TLorentzVector pvtosv(
                l_daughters[0].first->vtx().vx() - l_daughters[0].second->vtx().vx(),
                l_daughters[0].first->vtx().vy() - l_daughters[0].second->vtx().vy(),
                l_daughters[0].first->vtx().vz() - l_daughters[0].second->vtx().vz(),
                0.);

        if(l_daughters[0].first->charge()<0) {

          lvec1 = ConvertToLorentz(a1_daughters[0][1]->vector());
          lvec3 = ConvertToLorentz(a1_daughters[0][0]->vector());
          lvec4 = ConvertToLorentz(l_daughters[0].first->vector());
          TVector3 ip = (pvtosv.Vect() - pvtosv.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
          lvec2 = TLorentzVector(ip, 0.);
        } else {
          lvec2 = ConvertToLorentz(a1_daughters[0][1]->vector());
          lvec4 = ConvertToLorentz(a1_daughters[0][0]->vector());
          lvec3 = ConvertToLorentz(l_daughters[0].first->vector());
          TVector3 ip = (pvtosv.Vect() - pvtosv.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
          lvec1 = TLorentzVector(ip, 0.);
        }
        cp_sign_1_ = YRho(std::vector<GenParticle*>(
                    {a1_daughters[0][0], a1_daughters[0][1]}),TVector3());
    }

    if (rho_daughters.size()==2){
        cp_channel_=3;  
        if(rho_daughters[1].first->charge()<0) {
          lvec1 = ConvertToLorentz(rho_daughters[0].second->vector());   
          lvec2 = ConvertToLorentz(rho_daughters[1].second->vector());
          lvec3 = ConvertToLorentz(rho_daughters[0].first->vector());   
          lvec4 = ConvertToLorentz(rho_daughters[1].first->vector());
        } else {
          lvec1 = ConvertToLorentz(rho_daughters[1].second->vector());
          lvec2 = ConvertToLorentz(rho_daughters[0].second->vector());
          lvec3 = ConvertToLorentz(rho_daughters[1].first->vector());
          lvec4 = ConvertToLorentz(rho_daughters[0].first->vector());
        }
        cp_sign_1_ = YRho(std::vector<GenParticle*>(
                    {rho_daughters[0].first, rho_daughters[0].second}),TVector3())*
                    YRho(std::vector<GenParticle*>(
                    {rho_daughters[1].first, rho_daughters[1].second}),TVector3());
        Ediff_1_ = (rho_daughters[0].first->vector().E() - rho_daughters[0].second->vector().E())/
            (rho_daughters[0].first->vector().E() + rho_daughters[0].second->vector().E());
        Ediff_2_ = (rho_daughters[1].first->vector().E() - rho_daughters[1].second->vector().E())/
            (rho_daughters[1].first->vector().E() + rho_daughters[1].second->vector().E());
    }
    if(cp_channel_==2 || cp_channel_==5 || cp_channel_==12){
      aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);    
      if (cp_sign_1_>0) {
        if (aco_angle_1_<M_PI) aco_angle_1_ += M_PI;
        else                   aco_angle_1_ -= M_PI;
      }  
    }

    if(cp_channel_==1 || cp_channel_==11) {
      aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
    }

    if(cp_channel_!=-1 && cp_channel_==3){
      aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);    

      if (cp_sign_1_>0) {
        if (aco_angle_1_<M_PI) aco_angle_1_ += M_PI;
        else                   aco_angle_1_ -= M_PI;
      }  

      pi_E_1_ = rho_daughters[0].first->vector().E();
      pi_E_2_ = rho_daughters[1].first->vector().E();
      pi0_E_1_ = rho_daughters[0].second->vector().E();
      pi0_E_2_ = rho_daughters[1].second->vector().E();

    }
    
    if(charged_daughters.size()==2){
      cp_channel_=7;
      int taucharge1 = -charged_daughters[0][0]->charge(), taucharge2=-charged_daughters[1][0]->charge();
      for(auto x : charged_daughters[0]) taucharge1+=x->charge();
      for(auto x : charged_daughters[1]) taucharge2+=x->charge();
      // first tau
      auto os_pi_1 = charged_daughters[0][1];
      auto ss_pi_1_1 = charged_daughters[0][2];
      auto ss_pi_2_1 = charged_daughters[0][3]; 
      if(charged_daughters[0][2]->charge()!=charged_daughters[0][1]->charge() && charged_daughters[0][2]->charge()!=charged_daughters[0][3]->charge()) {
        os_pi_1 = charged_daughters[0][2];
        ss_pi_1_1 = charged_daughters[0][1];
        ss_pi_2_1 = charged_daughters[0][3];
      }
      if(charged_daughters[0][3]->charge()!=charged_daughters[0][1]->charge() && charged_daughters[0][3]->charge()!=charged_daughters[0][2]->charge()) {
        os_pi_1 = charged_daughters[0][3];
        ss_pi_1_1 = charged_daughters[0][1];
        ss_pi_2_1 = charged_daughters[0][2];
      }
      auto os_pi_2 = charged_daughters[1][1];
      auto ss_pi_1_2 = charged_daughters[1][2];
      auto ss_pi_2_2 = charged_daughters[1][3];
      if(charged_daughters[1][2]->charge()!=charged_daughters[1][1]->charge() && charged_daughters[1][2]->charge()!=charged_daughters[1][3]->charge()) {
        os_pi_2 = charged_daughters[1][2];
        ss_pi_1_2 = charged_daughters[1][1];
        ss_pi_2_2 = charged_daughters[1][3];
      }
      if(charged_daughters[1][3]->charge()!=charged_daughters[1][1]->charge() && charged_daughters[1][3]->charge()!=charged_daughters[1][2]->charge()) {
        os_pi_2 = charged_daughters[1][3];
        ss_pi_1_2 = charged_daughters[1][1];
        ss_pi_2_2 = charged_daughters[1][2];
      }
      TVector3 boost1 = ConvertToLorentz(charged_daughters[0][0]->vector()).BoostVector();
      TVector3 boost2 = ConvertToLorentz(charged_daughters[1][0]->vector()).BoostVector();
      TVector3 boost = (ConvertToLorentz(charged_daughters[0][0]->vector()+charged_daughters[1][0]->vector())).BoostVector();
      std::vector<TLorentzVector> TauA1andProd1 = {
        ConvertToLorentz(charged_daughters[0][0]->vector()),
        ConvertToLorentz(os_pi_1->vector()), 
        ConvertToLorentz(ss_pi_1_1->vector()), 
        ConvertToLorentz(ss_pi_2_1->vector()) 
      };
      std::vector<TLorentzVector> TauA1andProd2 = {
        ConvertToLorentz(charged_daughters[1][0]->vector()),
        ConvertToLorentz(os_pi_2->vector()),
        ConvertToLorentz(ss_pi_1_2->vector()),
        ConvertToLorentz(ss_pi_2_2->vector())
      };

      PolarimetricA1  a1pol_1;
      a1pol_1.Configure(TauA1andProd1,taucharge1);
      auto l_h1 = a1pol_1.PVC();
      PolarimetricA1  a1pol_2;
      a1pol_2.Configure(TauA1andProd2,taucharge2);
      auto l_h2 = a1pol_2.PVC();

      l_h1.Boost(-boost);
      l_h2.Boost(-boost);
 
      auto h1 = l_h1.Vect().Unit();
      auto h2 = l_h2.Vect().Unit();

      auto l_n1 = ConvertToLorentz(charged_daughters[0][0]->vector());
      auto l_n2 = ConvertToLorentz(charged_daughters[1][0]->vector());

      l_n1.Boost(-boost);
      l_n2.Boost(-boost);

      TVector3 n1 = l_n1.Vect().Unit();
      TVector3 n2 = l_n2.Vect().Unit();

      auto k1 = h1.Cross(n1).Unit();   
      auto k2 = h2.Cross(n2).Unit();  
                                      
      aco_angle_1_ = acos(k1.Dot(k2));
      double sign = h1.Cross(h2).Dot(n1);
      if(sign>0) aco_angle_1_ = 2*M_PI-aco_angle_1_;
                                      
    }                                 


    aco_angle_smear_=-9999;
    aco_sign_smear_=-9999;
    aco_sign_=-9999;

    if (pi_daughters.size()==2 && prho_daughters.size()==2){
        cp_channel_=101;

        TLorentzVector pvtosv1(
                pi_daughters[0]->vtx().vx() - prho_daughters[0].second->vtx().vx(),
                pi_daughters[0]->vtx().vy() - prho_daughters[0].second->vtx().vy(),
                pi_daughters[0]->vtx().vz() - prho_daughters[0].second->vtx().vz(),
                0.);

        TLorentzVector pvtosv2(
                pi_daughters[1]->vtx().vx() - prho_daughters[1].second->vtx().vx(),
                pi_daughters[1]->vtx().vy() - prho_daughters[1].second->vtx().vy(),
                pi_daughters[1]->vtx().vz() - prho_daughters[1].second->vtx().vz(),
                0.);

        TLorentzVector lvec3 = ConvertToLorentz(pi_daughters[0]->vector()); //pi charge from lep
        TLorentzVector lvec4 = ConvertToLorentz(pi_daughters[1]->vector()); //pi charge from tau

        TVector3 ip1 = (pvtosv1.Vect() - pvtosv1.Vect().Dot(lvec3.Vect().Unit())*lvec3.Vect().Unit()).Unit();
        TLorentzVector lvec1 = TLorentzVector(ip1, 0.);

        TVector3 ip2 = (pvtosv2.Vect() - pvtosv2.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
        TLorentzVector lvec2 = TLorentzVector(ip2, 0.);

        aco_angle_1_ = IPAcoAngle(lvec1, lvec2, lvec3, lvec4,false);
        aco_lin_angle_1_ = IPAcoLinAngle(lvec1, lvec2, lvec3, lvec4,false);
        angle_1_ = IPAngle(lvec1, lvec2, lvec3, lvec4,false);
        aco_sign_ = IPAcoSign(lvec1, lvec2, lvec3, lvec4,false);

        t_E_frac_1_=pi_daughters[0]->energy()/prho_daughters[0].second->energy();
        t_E_frac_2_=pi_daughters[1]->energy()/prho_daughters[1].second->energy();
        std::pair<double,double> c1c2_pair = CTauAngle(lvec1, lvec2, lvec3, lvec4,false);
        std::pair<double,double> a1a2_pair = IPAcoAngleSep(lvec1, lvec2, lvec3, lvec4,false);
        (void)a1a2_pair;
        t_c_1_ = c1c2_pair.first;
        t_c_2_ = c1c2_pair.second;

        rand->SetSeed((int)((pi_daughters[0]->eta()+5)*100000 + (pi_daughters[0]->phi()+4)*1000));

        double rx = rand->Gaus(0,0.001);
        double ry = rand->Gaus(0,0.001);
        double rz = rand->Gaus(0,0.0015);

        TLorentzVector pvtosv1_smear(
                pi_daughters[0]->vtx().vx() - prho_daughters[0].second->vtx().vx() -rx,
                pi_daughters[0]->vtx().vy() - prho_daughters[0].second->vtx().vy() -ry,
                pi_daughters[0]->vtx().vz() - prho_daughters[0].second->vtx().vz() -rz,
                0.);

        TLorentzVector pvtosv2_smear(
                pi_daughters[1]->vtx().vx() - prho_daughters[1].second->vtx().vx() -rx,
                pi_daughters[1]->vtx().vy() - prho_daughters[1].second->vtx().vy() -ry,
                pi_daughters[1]->vtx().vz() - prho_daughters[1].second->vtx().vz() -rz,
                0.);

        TVector3 ip1_smear = (pvtosv1_smear.Vect() - pvtosv1_smear.Vect().Dot(lvec3.Vect().Unit())*lvec3.Vect().Unit()).Unit();
        TLorentzVector lvec1_smear = TLorentzVector(ip1_smear, 0.);

        TVector3 ip2_smear = (pvtosv2_smear.Vect() - pvtosv2_smear.Vect().Dot(lvec4.Vect().Unit())*lvec4.Vect().Unit()).Unit();
        TLorentzVector lvec2_smear = TLorentzVector(ip2_smear, 0.);

        aco_angle_smear_ = IPAcoAngle(lvec1_smear, lvec2_smear, lvec3, lvec4,false);
        aco_sign_smear_ = IPAcoSign(lvec1_smear, lvec2_smear, lvec3, lvec4,false);


    }

    std::vector<ic::Vertex*> vertex_vec = {}; 
    if(event->ExistsInTree("vertices")) vertex_vec = event->GetPtrVec<ic::Vertex>("vertices"); 
    // std::vector<ic::Vertex*> & vertex_bs_vec = event->GetPtrVec<ic::Vertex>("verticesBS"); 
    std::vector<ic::Vertex*> vertex_bs_vec = {};
    if(event->ExistsInTree("refittedVerticesBS")) vertex_bs_vec = event->GetPtrVec<ic::Vertex>("refittedVerticesBS");

 
    if (vertex_vec.size() > 0) {
      reco_pvx_ = vertex_vec[0]->vx();
      reco_pvy_ = vertex_vec[0]->vy();
      reco_pvz_ = vertex_vec[0]->vz();
    }
    if (vertex_bs_vec.size() > 0) {
      reco_bs_pvx_ = vertex_bs_vec[0]->vx();
      reco_bs_pvy_ = vertex_bs_vec[0]->vy();
      reco_bs_pvz_ = vertex_bs_vec[0]->vz();
    } else if(vertex_vec.size() > 0) {
      reco_bs_pvx_ = vertex_vec[0]->vx();
      reco_bs_pvy_ = vertex_vec[0]->vy();
      reco_bs_pvz_ = vertex_vec[0]->vz();
    }
    if (gen_vertices.size() > 0) {
      gen_pvx_ = gen_vertices[0]->vx();
      gen_pvy_ = gen_vertices[0]->vy();
      gen_pvz_ = gen_vertices[0]->vz();
    }

    reco_pi_px_1_=-9999, reco_pi_py_1_=-9999, reco_pi_pz_1_=-9999, reco_pi_E_1_=-9999;
    reco_pi2_px_1_=-9999, reco_pi2_py_1_=-9999, reco_pi2_pz_1_=-9999, reco_pi2_E_1_=-9999;
    reco_pi3_px_1_=-9999, reco_pi3_py_1_=-9999, reco_pi3_pz_1_=-9999, reco_pi3_E_1_=-9999;
    reco_pi0_px_1_=-9999, reco_pi0_py_1_=-9999, reco_pi0_pz_1_=-9999, reco_pi0_E_1_=-9999;
    reco_sv_x_1_=-9999, reco_sv_y_1_=-9999, reco_sv_z_1_=-9999; 
    reco_ip_x_1_=-9999, reco_ip_y_1_=-9999, reco_ip_z_1_=-9999;
    reco_metx_=-9999, reco_mety_=-9999; 
    reco_dm_1_=-9999; 
    // Adding second tau info	
    reco_pi_px_2_=-9999, reco_pi_py_2_=-9999, reco_pi_pz_2_=-9999, reco_pi_E_2_=-9999;
    reco_pi2_px_2_=-9999, reco_pi2_py_2_=-9999, reco_pi2_pz_2_=-9999, reco_pi2_E_2_=-9999;
    reco_pi3_px_2_=-9999, reco_pi3_py_2_=-9999, reco_pi3_pz_2_=-9999, reco_pi3_E_2_=-9999;
    reco_pi0_px_2_=-9999, reco_pi0_py_2_=-9999, reco_pi0_pz_2_=-9999, reco_pi0_E_2_=-9999;
    reco_sv_x_2_=-9999, reco_sv_y_2_=-9999, reco_sv_z_2_=-9999; 
    reco_ip_x_2_=-9999, reco_ip_y_2_=-9999, reco_ip_z_2_=-9999;
    reco_dm_2_=-9999; 

    //if(event->ExistsInTree("taus") && gen_tau_jets_ptr.size()>=1) {
    if(event->ExistsInTree("taus") && gen_tau_jets_ptr.size()>=2) {
      std::vector<GenJet *> lead_gen_tau;
      std::vector<GenJet *> sublead_gen_tau;

      if(tau1_charge == gen_tau_jets_ptr[0]->charge() && tau2_charge == gen_tau_jets_ptr[1]->charge()){
         lead_gen_tau = {gen_tau_jets_ptr[0]};
         sublead_gen_tau = {gen_tau_jets_ptr[1]};
      } else {
         lead_gen_tau = {gen_tau_jets_ptr[1]};
         sublead_gen_tau = {gen_tau_jets_ptr[0]};
      }

      std::vector<Tau *> taus = event->GetPtrVec<Tau>("taus");
      ic::erase_if(taus,!boost::bind(MinPtMaxEta, _1, 20.0, 2.3));         
      std::vector<Tau *> taus_tt;

      for(auto t : taus) {
        if(t->GetTauID("byVVVLooseDeepTau2017v2p1VSjet") > 0.5 && t->GetTauID("byVVLooseDeepTau2017v2p1VSe") > 0.5 && t->GetTauID("byVLooseDeepTau2017v2p1VSmu") > 0.5 && t->GetTauID("decayModeFindingNewDMs") > 0.5 && (t->decay_mode()<2 || t->decay_mode()>9) && fabs(t->charge()) == 1 && fabs(t->lead_dz_vertex()) < 0.2) {
        taus_tt.push_back(t);
        }
      }
      std::vector<std::pair<ic::GenJet *, ic::Tau *>> tt_matches_lead =  MatchByDR(lead_gen_tau,taus_tt,0.5,true,true); 
      std::vector<std::pair<ic::GenJet *, ic::Tau *>> tt_matches_sublead =  MatchByDR(sublead_gen_tau,taus_tt,0.5,true,true);

      if (tt_matches_lead.size()>0 && tt_matches_sublead.size()>0) {

        std::vector<ic::Vertex*> & vertex_vec = event->GetPtrVec<ic::Vertex>("vertices");
        ic::Vertex* refit_vertex = vertex_vec[0];

        ic::Tau *reco_tau_lead = tt_matches_lead[0].second;
        ic::Tau *reco_tau_sublead = tt_matches_sublead[0].second;

        std::vector<ic::PFCandidate*> pfcands =  event->GetPtrVec<ic::PFCandidate>("pfCandidates");

        std::vector<Met*> reco_met_vec = event->GetPtrVec<Met>("puppiMet");
        Met *reco_met = reco_met_vec.at(0);
        reco_metx_ = reco_met->vector().Px();
        reco_mety_ = reco_met->vector().Py();
      
        int tau_decay_mode_1_ = reco_tau_lead->decay_mode();
        int tau_decay_mode_2_ = reco_tau_sublead->decay_mode();
 
        reco_dm_1_ = reco_tau_lead->HasTauID("MVADM2017v1") ? reco_tau_lead->GetTauID("MVADM2017v1") : 0.0;
        reco_dm_2_ = reco_tau_sublead->HasTauID("MVADM2017v1") ? reco_tau_sublead->GetTauID("MVADM2017v1") : 0.0;

        if(tau_decay_mode_1_==0&&reco_dm_1_>0) reco_dm_1_=-1;
        if(tau_decay_mode_2_==0&&reco_dm_2_>0) reco_dm_2_=-1;

        if(tau_decay_mode_1_<3){
          std::pair<TVector3,double> ipandsig_1 = IPAndSignificance(reco_tau_lead, refit_vertex, pfcands);
          TVector3 ip_1 = ipandsig_1.first;
          reco_ip_x_1_ = ip_1.X();
          reco_ip_y_1_ = ip_1.Y();
          reco_ip_z_1_ = ip_1.Z();
        }
        if(tau_decay_mode_2_<3){
          std::pair<TVector3,double> ipandsig_2 = IPAndSignificance(reco_tau_sublead, refit_vertex, pfcands);
          TVector3 ip_2 = ipandsig_2.first;
          reco_ip_x_2_ = ip_2.X();
          reco_ip_y_2_ = ip_2.Y();
          reco_ip_z_2_ = ip_2.Z();
        }

        if(tau_decay_mode_1_<3) {
          reco_pi_px_1_=reco_tau_lead->vector().Px();
          reco_pi_py_1_=reco_tau_lead->vector().Py();
          reco_pi_pz_1_=reco_tau_lead->vector().Pz();
          reco_pi_E_1_ =reco_tau_lead->vector().E(); 
        }
        if(tau_decay_mode_2_<3) {
          reco_pi_px_2_=reco_tau_sublead->vector().Px();
          reco_pi_py_2_=reco_tau_sublead->vector().Py();
          reco_pi_pz_2_=reco_tau_sublead->vector().Pz();
          reco_pi_E_2_ =reco_tau_sublead->vector().E();
        }


        if(tau_decay_mode_1_==1) {
          std::pair<ic::Candidate*, ic::Candidate*> rho1 = GetRho(reco_tau_lead, pfcands);
          ic::Candidate *pi_tau1 = rho1.first;
          ic::Candidate *pi0_tau1 = rho1.second;
          reco_pi_px_1_=pi_tau1->vector().Px();
          reco_pi_py_1_=pi_tau1->vector().Py();
          reco_pi_pz_1_=pi_tau1->vector().Pz();
          reco_pi_E_1_=pi_tau1->vector().E(); 
          reco_pi0_px_1_=pi0_tau1->vector().Px();
          reco_pi0_py_1_=pi0_tau1->vector().Py();
          reco_pi0_pz_1_=pi0_tau1->vector().Pz();
          reco_pi0_E_1_=pi0_tau1->vector().E();
        }
        if(tau_decay_mode_2_==1) {
          std::pair<ic::Candidate*, ic::Candidate*> rho2 = GetRho(reco_tau_sublead, pfcands);
          ic::Candidate *pi_tau2 = rho2.first;
          ic::Candidate *pi0_tau2 = rho2.second;
          reco_pi_px_2_=pi_tau2->vector().Px();
          reco_pi_py_2_=pi_tau2->vector().Py();
          reco_pi_pz_2_=pi_tau2->vector().Pz();
          reco_pi_E_2_=pi_tau2->vector().E();
          reco_pi0_px_2_=pi0_tau2->vector().Px();
          reco_pi0_py_2_=pi0_tau2->vector().Py();
          reco_pi0_pz_2_=pi0_tau2->vector().Pz();
          reco_pi0_E_2_=pi0_tau2->vector().E();
        }

        if(tau_decay_mode_1_==10 || tau_decay_mode_1_==11){
          if(reco_tau_lead->hasSV()){
            reco_sv_x_1_ = reco_tau_lead->secondary_vertex().X() - refit_vertex->vx();
            reco_sv_y_1_ = reco_tau_lead->secondary_vertex().Y() - refit_vertex->vy();
            reco_sv_z_1_ = reco_tau_lead->secondary_vertex().Z() - refit_vertex->vz();
          }
          std::vector<ic::PFCandidate*> hads1 = GetTauHads(reco_tau_lead, pfcands,0);
          if(hads1.size()>2){
             hads1 = GetA1 (reco_tau_lead, pfcands).first;
             reco_pi_px_1_=hads1[0]->vector().Px();
             reco_pi_py_1_=hads1[0]->vector().Py();
             reco_pi_pz_1_=hads1[0]->vector().Pz();
             reco_pi_E_1_=hads1[0]->vector().E();
             reco_pi2_px_1_=hads1[1]->vector().Px();
             reco_pi2_py_1_=hads1[1]->vector().Py();
             reco_pi2_pz_1_=hads1[1]->vector().Pz();
             reco_pi2_E_1_=hads1[1]->vector().E();
             reco_pi3_px_1_=hads1[2]->vector().Px();
             reco_pi3_py_1_=hads1[2]->vector().Py();
             reco_pi3_pz_1_=hads1[2]->vector().Pz();
             reco_pi3_E_1_=hads1[2]->vector().E();
          }
          if(tau_decay_mode_1_==11){
             ic::Candidate* gamma1 = GetA1(reco_tau_lead, pfcands).second;
             reco_pi0_px_1_=gamma1->vector().Px();
             reco_pi0_py_1_=gamma1->vector().Py();
             reco_pi0_pz_1_=gamma1->vector().Pz();
             reco_pi0_E_1_= gamma1->vector().E();
          } 
        }
        if(tau_decay_mode_2_==10 || tau_decay_mode_2_==11){
          if(reco_tau_sublead->hasSV()){
            reco_sv_x_2_ = reco_tau_sublead->secondary_vertex().X() - refit_vertex->vx();
            reco_sv_y_2_ = reco_tau_sublead->secondary_vertex().Y() - refit_vertex->vy();
            reco_sv_z_2_ = reco_tau_sublead->secondary_vertex().Z() - refit_vertex->vz();
          }
          std::vector<ic::PFCandidate*> hads2 = GetTauHads(reco_tau_sublead, pfcands,0);
          if(hads2.size()>2){
             hads2 = GetA1 (reco_tau_sublead, pfcands).first;
             reco_pi_px_2_=hads2[0]->vector().Px();
             reco_pi_py_2_=hads2[0]->vector().Py();
             reco_pi_pz_2_=hads2[0]->vector().Pz();
             reco_pi_E_2_=hads2[0]->vector().E();
             reco_pi2_px_2_=hads2[1]->vector().Px();
             reco_pi2_py_2_=hads2[1]->vector().Py();
             reco_pi2_pz_2_=hads2[1]->vector().Pz();
             reco_pi2_E_2_=hads2[1]->vector().E();
             reco_pi3_px_2_=hads2[2]->vector().Px();
             reco_pi3_py_2_=hads2[2]->vector().Py();
             reco_pi3_pz_2_=hads2[2]->vector().Pz();
             reco_pi3_E_2_=hads2[2]->vector().E();
          }
          if(tau_decay_mode_2_==11){
             ic::Candidate* gamma2 = GetA1(reco_tau_sublead, pfcands).second;
             reco_pi0_px_2_=gamma2->vector().Px();
             reco_pi0_py_2_=gamma2->vector().Py();
             reco_pi0_pz_2_=gamma2->vector().Pz();
             reco_pi0_E_2_= gamma2->vector().E();
          }
        }
      }
    }

    pv_aco_angle_=-9999;
    pv_lin_angle_=-9999;
    
    costheta_1_=-9999;
    costheta_tf_1_=-9999;
    costheta_2_=-9999;
    costheta_tf_2_=-9999;
    q_1_=-9999;
    q_2_=-9999;
    
    if(undecayed_taus.size()>1) {
    
      q_1_ = undecayed_taus[0]->charge();
      q_2_ = undecayed_taus[1]->charge();
    
      TLorentzVector l_h1= PolarimetricVector(pis[0], pi0s[0], undecayed_taus[0],false);
      TLorentzVector l_h2= PolarimetricVector(pis[1], pi0s[1], undecayed_taus[1],false);
      TLorentzVector l_tau_1 = ConvertToLorentz(undecayed_taus[0]->vector());
      TLorentzVector l_tau_2 = ConvertToLorentz(undecayed_taus[1]->vector());
    
      TVector3 boost = ConvertToLorentz(undecayed_taus[0]->vector()+undecayed_taus[1]->vector()).BoostVector();
      TVector3 boost1 = ConvertToLorentz(undecayed_taus[0]->vector()).BoostVector();
      TVector3 boost2 = ConvertToLorentz(undecayed_taus[1]->vector()).BoostVector();
   
      TVector3 boson = (ConvertToLorentz(undecayed_taus[0]->vector()+undecayed_taus[1]->vector())).Vect();
 
      ROOT::Math::PtEtaPhiEVector vis_tau_1, vis_tau_2, vis_system;
      for (auto p : pis[0]) vis_tau_1+=p->vector();
      for (auto p : pis[1]) vis_tau_2+=p->vector();
      for (auto p : pi0s[0]) vis_tau_1+=p->vector();
      for (auto p : pi0s[1]) vis_tau_2+=p->vector(); 
   
      // checks gamma seperation in single pi0 events
      gammas_deta_1_=-9999;
      gammas_dphi_1_=-9999;
      gammas_deta_2_=-9999;
      gammas_dphi_2_=-9999;

      pi0s_deta_1_=-9999;
      pi0s_dphi_1_=-9999;
      pi0s_deta_2_=-9999;
      pi0s_dphi_2_=-9999;

      rho_deta_1_=-9999;
      rho_dphi_1_=-9999;
      rho_deta_2_=-9999;
      rho_dphi_2_=-9999;


      if (pi0s.size()>0 && pi0s[0].size()==1 && pis[0].size()==1) {
        rho_deta_1_ = pi0s[0][0]->eta() - pis[0][0]->eta();
        rho_dphi_1_ = ROOT::Math::VectorUtil::DeltaPhi(pi0s[0][0]->vector(), pis[0][0]->vector());
        auto gammas = pi0s[0][0]->daughters();
        if (gammas.size()==2) {
          if(gen_particles[gammas[0]]->pdgid()==22 && gen_particles[gammas[0]]->pdgid()==22) {
            gammas_deta_1_ = gen_particles[gammas[0]]->eta() - gen_particles[gammas[1]]->eta();
            gammas_dphi_1_ = ROOT::Math::VectorUtil::DeltaPhi(gen_particles[gammas[0]]->vector(), gen_particles[gammas[1]]->vector()); 
          }
        }
      }
      if (pi0s.size()>1 && pi0s[1].size()==1 && pis[1].size()==1) {
        rho_deta_2_ = pi0s[1][0]->eta() - pis[1][0]->eta();
        rho_dphi_2_ = ROOT::Math::VectorUtil::DeltaPhi(pi0s[1][0]->vector(), pis[1][0]->vector());
        auto gammas = pi0s[1][0]->daughters();
        if (gammas.size()==2) {
          if(gen_particles[gammas[0]]->pdgid()==22 && gen_particles[gammas[0]]->pdgid()==22) {
            gammas_deta_2_ = gen_particles[gammas[0]]->eta() - gen_particles[gammas[1]]->eta();
            gammas_dphi_2_ = ROOT::Math::VectorUtil::DeltaPhi(gen_particles[gammas[0]]->vector(), gen_particles[gammas[1]]->vector());
          }
        }
      }

      if (pi0s.size()>0 && pi0s[0].size()==2 && pis[0].size()==1) {
        pi0s_deta_1_ = pi0s[0][0]->eta() - pi0s[0][1]->eta();
        pi0s_dphi_1_ = ROOT::Math::VectorUtil::DeltaPhi(pi0s[0][0]->vector(), pi0s[0][1]->vector());
      }

      if (pi0s.size()>1 && pi0s[1].size()==2 && pis[1].size()==1) {
        pi0s_deta_2_ = pi0s[1][0]->eta() - pi0s[1][1]->eta();
        pi0s_dphi_2_ = ROOT::Math::VectorUtil::DeltaPhi(pi0s[1][0]->vector(), pi0s[1][1]->vector());
      }
 
      vis_system=vis_tau_1+vis_tau_2;
    
      TVector3 vis_boost = ConvertToLorentz(vis_system).BoostVector();
    
      TVector3 h1, h2, n1, n2, k1, k2;
      double pv_sign=0.;
    
      l_h1= PolarimetricVector(pis[0], pi0s[0], undecayed_taus[0], false);
      l_h2= PolarimetricVector(pis[1], pi0s[1], undecayed_taus[1], false);
      l_tau_1 = ConvertToLorentz(undecayed_taus[0]->vector());
      l_tau_2 = ConvertToLorentz(undecayed_taus[1]->vector());

     
      taupos_px_=-9999;
      taupos_py_=-9999;
      taupos_pz_=-9999;
      taupos_E_=-9999;
      tauneg_px_=-9999;
      tauneg_py_=-9999;
      tauneg_pz_=-9999;
      tauneg_E_=-9999;
      taupos_polvec_x_=-9999;
      taupos_polvec_y_=-9999;
      taupos_polvec_z_=-9999;
      taupos_polvec_E_=-9999;
      tauneg_polvec_x_=-9999;
      tauneg_polvec_y_=-9999;
      tauneg_polvec_z_=-9999;
      tauneg_polvec_E_=-9999;

      if(undecayed_taus[0]->charge()>0&&undecayed_taus[1]->charge()<0) {
        taupos_px_=l_tau_1.X();
        taupos_py_=l_tau_1.Y();
        taupos_pz_=l_tau_1.Z();
        taupos_E_=l_tau_1.T();
        tauneg_px_=l_tau_2.X();
        tauneg_py_=l_tau_2.Y();
        tauneg_pz_=l_tau_2.Z();
        tauneg_E_=l_tau_2.T();

        taupos_polvec_x_=l_h1.X();
        taupos_polvec_y_=l_h1.Y();
        taupos_polvec_z_=l_h1.Z();
        taupos_polvec_E_=l_h1.T();
        tauneg_polvec_x_=l_h2.X();
        tauneg_polvec_y_=l_h2.Y();
        tauneg_polvec_z_=l_h2.Z();
        tauneg_polvec_E_=l_h2.T();
      } else if(undecayed_taus[0]->charge()<0&&undecayed_taus[1]->charge()>0) {
        taupos_px_=l_tau_2.X();
        taupos_py_=l_tau_2.Y();
        taupos_pz_=l_tau_2.Z();
        taupos_E_=l_tau_2.T();
        tauneg_px_=l_tau_1.X();
        tauneg_py_=l_tau_1.Y();
        tauneg_pz_=l_tau_1.Z();
        tauneg_E_=l_tau_1.T();

        taupos_polvec_x_=l_h2.X();
        taupos_polvec_y_=l_h2.Y();
        taupos_polvec_z_=l_h2.Z();
        taupos_polvec_E_=l_h2.T();
        tauneg_polvec_x_=l_h1.X();
        tauneg_polvec_y_=l_h1.Y();
        tauneg_polvec_z_=l_h1.Z();
        tauneg_polvec_E_=l_h1.T();
      } 

 
      l_h1.Boost(-boost);
      l_h2.Boost(-boost);
    
      l_tau_1.Boost(-boost);
      l_tau_2.Boost(-boost);

      h1 = l_h1.Vect().Unit();
      h2 = l_h2.Vect().Unit();

      n1 = l_tau_1.Vect().Unit();
      n2 = l_tau_2.Vect().Unit();

      k1 = h1.Cross(n1).Unit();
      k2 = h2.Cross(n2).Unit();

      pv_lin_angle_ = acos(h1.Dot(h2));
      pv_aco_angle_ = acos(k1.Dot(k2));
      pv_sign = h1.Cross(h2).Dot(n1);
      if(pv_sign>0) pv_aco_angle_ = 2*M_PI-pv_aco_angle_;
    
      costheta_1_=h1.Dot(n1);
      costheta_2_=h2.Dot(n2);
    
      //bost back to lab frame, then to tau frame
      l_h1.Boost(boost);
      l_h2.Boost(boost);
      l_h1.Boost(-boost1);
      l_h2.Boost(-boost2);
    
      TVector3 h1_tf = l_h1.Vect().Unit();
      TVector3 h2_tf = l_h2.Vect().Unit();
    
      costheta_tf_1_=h1_tf.Dot(n1);
      costheta_tf_2_=h2_tf.Dot(n2);

 
    }



    if(fs_) outtree_->Fill();
    if(make_mva_ntuple_ && tauFlag_1_>=0 && tauFlag_2_>=0) mvatree_->Fill();
    return 0;
  }
  int HTTGenAnalysis::PostAnalysis() {
    std::cout << "ee count = " << count_ee_ << std::endl;
    std::cout << "em count = " << count_em_ << std::endl;
    std::cout << "et count = " << count_et_ << std::endl;
    std::cout << "mm count = " << count_mm_ << std::endl;
    std::cout << "mt count = " << count_mt_ << std::endl;
    std::cout << "tt count = " << count_tt_ << std::endl;
    std::cout << "count taus = " << count_taus_ << std::endl;

    if(make_mva_ntuple_) {
      lOFile->cd();
      mvatree_->Write();
      lOFile->Close();
    }

    return 0;
  }

  void HTTGenAnalysis::PrintInfo() {
    ;
  }

}
