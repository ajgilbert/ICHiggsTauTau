#ifndef ICHiggsTauTau_HiggsTauTauRun2_HTTGenAnalysis_h
#define ICHiggsTauTau_HiggsTauTauRun2_HTTGenAnalysis_h

#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/TreeEvent.h"
#include "UserCode/ICHiggsTauTau/Analysis/Core/interface/ModuleBase.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/HistoSet.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/interface/HTTConfig.h"
#include "UserCode/ICHiggsTauTau/interface/GenParticle.hh"
#include "UserCode/ICHiggsTauTau/interface/GenJet.hh"
#include "UserCode/ICHiggsTauTau/interface/Candidate.hh"
#include "UserCode/ICHiggsTauTau/interface/Electron.hh"
#include "UserCode/ICHiggsTauTau/interface/Muon.hh"
#include "UserCode/ICHiggsTauTau/interface/CompositeCandidate.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "TRandom3.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"
#include <string>
#include "Utilities/interface/FnRootTools.h"

#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/BTagCalibrationStandalone.h"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/interface/BTagWeightRun2.h"

namespace ic {
  
class HTTGenAnalysis : public ModuleBase {
 private:
  CLASS_MEMBER(HTTGenAnalysis, fwlite::TFileService*, fs)
  CLASS_MEMBER(HTTGenAnalysis, std::string, channel_str)
  CLASS_MEMBER(HTTGenAnalysis, double, min_jet_pt    )
  CLASS_MEMBER(HTTGenAnalysis, double, max_jet_eta  )
  CLASS_MEMBER(HTTGenAnalysis, double, min_e_pt     )
  CLASS_MEMBER(HTTGenAnalysis, double, min_mu_pt    )
  CLASS_MEMBER(HTTGenAnalysis, double, min_tau1_pt  )
  CLASS_MEMBER(HTTGenAnalysis, double, min_tau2_pt  )
  CLASS_MEMBER(HTTGenAnalysis, double, max_e_eta    )
  CLASS_MEMBER(HTTGenAnalysis, double, max_mu_eta   )
  CLASS_MEMBER(HTTGenAnalysis, double, max_tau_eta  )
  CLASS_MEMBER(HTTGenAnalysis, bool, do_theory_uncert)
  CLASS_MEMBER(HTTGenAnalysis, TH2F*,bbtag_eff)
  CLASS_MEMBER(HTTGenAnalysis, std::string, mssm_mass  )
  CLASS_MEMBER(HTTGenAnalysis, bool, make_mva_ntuple)
  CLASS_MEMBER(HTTGenAnalysis, std::string, mva_output_name)
  CLASS_MEMBER(HTTGenAnalysis, ic::era, era)
  
  TTree *outtree_;
  TRandom3  *rand;
  TTree *mvatree_;
  TFile *lOFile; 
 
  unsigned long long event_;

  unsigned n_rho_;
  unsigned n_tot_;
  
  unsigned count_ee_;
  unsigned count_em_;
  unsigned count_et_;
  unsigned count_mm_;
  unsigned count_mt_;
  unsigned count_tt_;
  unsigned count_taus_;
  
  bool passed_;
  double pt_1_;
  double pt_2_;
  double eta_1_;
  double eta_2_;
  double phi_1_;
  double phi_2_;
  double met_;
  double mt_1_;
  double mt_2_;
  double pzeta_;
  double n_bjets_;
  double n_bjets_eta2p5_;
  double n_bjets_noscale_;
  double n_bjets_eta2p5_noscale_;
  double n_bjets_pt25_;
  int n_bpartons_;

  unsigned n_jets_nofilter_;
  unsigned n_jets_;
  unsigned n_jetsingap_;
  double jpt_1_;  
  double jpt_2_;  
  double jpt_3_;
  double bpt_1_;  
  double jeta_1_; 
  double jeta_2_; 
  double jphi_1_;
  double jphi_2_;
  double mjj_;
  double pmjj_;
  double mjj2_;
  double minmjj_;
  double maxmjj_;
  double jdR_;
  double pjdR_;
  double jdR2_;
  double jdeta_;
  double m_vis_;
  double pt_tt_;
  double mt_tot_;
  double wtzpt_;
  double mass_;
  double wt_;
  double gen_wt_;
  bool hasTaus_;
  double wt_z_pol_;
  double wt_cp_sm_;
  double wt_cp_ps_;
  double wt_cp_mm_;
  double wt_cp_sm_alt_;
  double wt_cp_ps_alt_;
  double wt_cp_mm_alt_;
  double wt_vlq_sm_;
  double wt_vlq_full_;
  double WTm_, WTp_;
  double wt_cp_prod_sm_;
  double wt_cp_prod_ps_;
  double wt_cp_prod_mm_;

  double HiggsPt_;
  std::string decayType;
  double wt_ggh_pt_           ;
  double wt_ggh_pt_up_        ;
  double wt_ggh_pt_down_      ;
  double wt_ggh_pt_herwig_    ;
  double wt_ggh_pt_amc_       ;
  double wt_ggh_pt_pythiaup_  ;
  double wt_ggh_pt_pythiadown_;
  double wt_ggh_pt_scalehigh_ ;
  double wt_ggh_pt_scalelow_  ;
  double genpt_1_;
  double genpt_2_;
  double geneta_1_;
  double geneta_2_;
  double n_jets_offline_;
  double n_bjets_offline_;
  
  unsigned partons_lhe_;
  unsigned partons_;
  double parton_pt_;
  double parton_pt_2_;
  double parton_pt_3_;
  double parton_mjj_;
  double parton_HpT_;
  double parton_Zmass_;
  float D0_;
  float DCP_;
  float D0star_;
  double sjdphi_;
  double sjdphi2_;
  double sjdphi3_;

  double spjdphi_;
  double ysep_;
  unsigned n_pjets_;


  float n_pu_;
  
  double pT_A_;
 
  TH1F topmass_wts_;
  TH1F topmass_wts_toponly_;

  TH1D ps_0jet_up_;
  TH1D ps_0jet_down_;
  TH1D ps_1jet_up_;
  TH1D ps_1jet_down_;
  TH1D ps_2jet_up_;
  TH1D ps_2jet_down_;
  TH1D ps_3jet_up_;
  TH1D ps_3jet_down_;
  TH1D ue_down_;
  TH1D ue_up_;

  double wt_topmass_;
  double wt_topmass_2_;
  double wt_ps_down_;
  double wt_ps_up_;
  double wt_ue_down_;
  double wt_ue_up_;
  double scale1_;
  double scale2_;
  double scale3_;
  double scale4_;
  double scale5_;
  double scale6_;
  double scale7_;
  double scale8_;
  double scale9_;

  double scale_vlq_1_;
  double scale_vlq_2_;
  double scale_vlq_3_;
  double scale_vlq_4_;
  double scale_vlq_5_;
  double scale_vlq_6_;
  double scale_vlq_7_;
  double scale_vlq_8_;
 
  double pdfweight_0_;
  double pdfweight_1_;
  double pdfweight_2_;
  double pdfweight_3_;
  double pdfweight_4_;
  double pdfweight_5_;
  double pdfweight_6_;
  double pdfweight_7_;
  double pdfweight_8_;
  double pdfweight_9_;
  double pdfweight_10_;
  double pdfweight_11_;
  double pdfweight_12_;
  double pdfweight_13_;
  double pdfweight_14_;
  double pdfweight_15_;
  double pdfweight_16_;
  double pdfweight_17_;
  double pdfweight_18_;
  double pdfweight_19_;
  double pdfweight_20_;
  double pdfweight_21_;
  double pdfweight_22_;
  double pdfweight_23_;
  double pdfweight_24_;
  double pdfweight_25_;
  double pdfweight_26_;
  double pdfweight_27_;
  double pdfweight_28_;
  double pdfweight_29_;
  double pdfweight_30_;

  double wt_lhe_nominal_;
  double wt_h_tb_;
  double wt_h_t_;
  double wt_h_b_;
  double wt_A_tb_;
  double wt_A_t_;
  double wt_A_b_;
  double wt_H_tb_;
  double wt_H_t_;
  double wt_H_b_; 

  double wt_hfact_nom_;
  double wt_hfact_up_;
  double wt_hfact_down_;

  std::shared_ptr<RooWorkspace> mssm_w_;
  std::map<std::string, std::shared_ptr<RooFunctor>> fns_;
  
  double wt_ggh_t_;
  double wt_ggh_b_;
  double wt_ggh_i_;
  double wt_ggH_t_;
  double wt_ggH_b_;
  double wt_ggH_i_;
  double wt_ggA_t_;
  double wt_ggA_b_;
  double wt_ggA_i_;
  double wt_stitch_;
  int npNLO_;

  bool cand_1_;
  bool cand_2_;
  bool match_1_;
  bool match_2_;
  
  TH2D z_pt_weights_sm_;

  double aco_lin_angle_1_;
  double angle_1_;
  double spin_var_;
  double spin_var_vis_;
  double mass_undecayed_taus_;

  double t_E_frac_1_;
  double t_E_frac_2_;
  double t_c_1_;
  double t_c_2_;

  double p_aco_angle_1_;
  double p_aco_angle_2_;
  double p_aco_angle_3_;
  double p_aco_angle_4_;
  double dR1_, dR2_;

  double aco_sign_;
  double aco_sign_smear_;
  double aco_angle_smear_;
  double pv_lin_angle_;
  double pv_aco_angle_;
  double costheta_1_;
  double costheta_tf_1_;
  double costheta_2_;
  double costheta_tf_2_;

  double aco_angle_1_;
  double aco_angle_2_;
  double aco_angle_3_;
  double aco_angle_4_;
  double aco_angle_5_;
  double aco_angle_6_;
  double cp_sign_1_;
  double cp_sign_2_;
  double cp_sign_3_;
  double cp_sign_4_;
  double cp_sign_5_;
  double cp_sign_6_;
  double mass1_, mass2_; 
  double a1_mass_;
  int q_1_;
  int q_2_;

  double pt_rho_pi_ = -1;
  double pt_rho_pi0_ = -1;
  double pt_other_pi0_ =-1;

  double E_1_1_;
  double E_2_1_;
  double E_1_2_;
  double E_2_2_;
  double E_3_2_;

  double y_1_1_;
  double y_2_1_;
  double y_3_1_;
  double y_4_1_;

  double y_1_2_;
  double y_2_2_;
  double y_3_2_;
  double y_4_2_;
  double rand_;
 
  int cp_channel_;
  
  double ip_dxy_res_1_; 
  double ip_dxy_res_2_; 
  double ip_dx_res_1_; 
  double ip_dx_res_2_; 
  double ip_dy_res_1_; 
  double ip_dy_res_2_;
  double ip_dz_res_1_;
  double ip_dz_res_2_;

  double Ediff_1_;
  double Ediff_2_;
  double Pfrac_1_; 
  double Pfrac_2_;

  double reco_pvx_;
  double reco_pvy_;
  double reco_pvz_;

  double reco_bs_pvx_;
  double reco_bs_pvy_;
  double reco_bs_pvz_;
  double gen_pvx_;
  double gen_pvy_;
  double gen_pvz_;

  double fakea1_dR_ = -1;
  double rho_dR_ = -1;
  double tauFlag_1_ = -1;
  double tauFlag_2_ = -1;

  double largest_gen_mjj_=-9999;

  double lead_b_eta_=-9999, lead_b_pt_=-9999;

  double wt_ps_isr_up_  ;
  double wt_ps_isr_down_;
  double wt_ps_fsr_up_;
  double wt_ps_fsr_down_;

  float tau_pt_1_tt_, tau_pt_1_mt_, tau_pt_1_et_, tau_pt_1_sf_, gen_tau_pt_1_, gen_tau_eta_1_,tau_pt_2_tt_, tau_pt_2_mt_, tau_pt_2_et_, tau_pt_2_sf_, gen_tau_pt_2_, gen_tau_eta_2_;
  int gen_tau_dm_1_, gen_tau_dm_2_;
  float pi_px_1_, pi_py_1_, pi_pz_1_, pi_E_1_, pi_px_2_, pi_py_2_, pi_pz_2_, pi_E_2_;
  float pi2_px_1_, pi2_py_1_, pi2_pz_1_, pi2_E_1_, pi2_px_2_, pi2_py_2_, pi2_pz_2_, pi2_E_2_;
  float pi3_px_1_, pi3_py_1_, pi3_pz_1_, pi3_E_1_, pi3_px_2_, pi3_py_2_, pi3_pz_2_, pi3_E_2_;
  float pi0_px_1_, pi0_py_1_, pi0_pz_1_, pi0_E_1_, pi0_px_2_, pi0_py_2_, pi0_pz_2_, pi0_E_2_;
  float nu_px_1_, nu_py_1_, nu_pz_1_, nu_E_1_, nu_px_2_, nu_py_2_, nu_pz_2_, nu_E_2_;
  float sv_x_1_, sv_y_1_, sv_z_1_, sv_x_2_, sv_y_2_, sv_z_2_;
  float ip_x_1_, ip_y_1_, ip_z_1_, ip_x_2_, ip_y_2_, ip_z_2_;
  float metx_, mety_;
  int tau1_charge, tau2_charge;

  double reco_pi_px_1_, reco_pi_py_1_, reco_pi_pz_1_, reco_pi_E_1_;
  double reco_pi2_px_1_, reco_pi2_py_1_, reco_pi2_pz_1_, reco_pi2_E_1_;
  double reco_pi3_px_1_, reco_pi3_py_1_, reco_pi3_pz_1_, reco_pi3_E_1_;
  double reco_pi0_px_1_, reco_pi0_py_1_, reco_pi0_pz_1_, reco_pi0_E_1_;
  double reco_sv_x_1_, reco_sv_y_1_, reco_sv_z_1_;
  double reco_ip_x_1_, reco_ip_y_1_, reco_ip_z_1_;
  double reco_metx_, reco_mety_;
  int reco_dm_1_;
 
  double reco_pi_px_2_, reco_pi_py_2_, reco_pi_pz_2_, reco_pi_E_2_;
  double reco_pi2_px_2_, reco_pi2_py_2_, reco_pi2_pz_2_, reco_pi2_E_2_;
  double reco_pi3_px_2_, reco_pi3_py_2_, reco_pi3_pz_2_, reco_pi3_E_2_;
  double reco_pi0_px_2_, reco_pi0_py_2_, reco_pi0_pz_2_, reco_pi0_E_2_;
  double reco_sv_x_2_, reco_sv_y_2_, reco_sv_z_2_;
  double reco_ip_x_2_, reco_ip_y_2_, reco_ip_z_2_;
  int reco_dm_2_;

  double gammas_deta_1_;
  double gammas_dphi_1_;
  double gammas_deta_2_;
  double gammas_dphi_2_;

  double pi0s_deta_1_;
  double pi0s_dphi_1_;
  double pi0s_deta_2_;
  double pi0s_dphi_2_;

  double rho_deta_1_;
  double rho_dphi_1_;
  double rho_deta_2_;
  double rho_dphi_2_;
 
  const BTagCalibration *calib;
  BTagCalibrationReader* reader_incl;
  BTagCalibrationReader* reader_comb;

  TFile *ggh_weights_;
  TGraph *ggh_ph_0jet_;
  TGraph *ggh_ph_1jet_;
  TGraph *ggh_ph_2jet_;
  TGraph *ggh_ph_3jet_;
  double wt_nnlops_;

  unsigned quarks_initial_;
  unsigned quarks_final_;
  unsigned gluons_initial_;
  unsigned gluons_final_;

  float taupos_px_;
  float taupos_py_;
  float taupos_pz_;
  float taupos_E_;
  float tauneg_px_;
  float tauneg_py_;
  float tauneg_pz_;
  float tauneg_E_;
  float taupos_polvec_x_;
  float taupos_polvec_y_;
  float taupos_polvec_z_;
  float taupos_polvec_E_;
  float tauneg_polvec_x_;
  float tauneg_polvec_y_;
  float tauneg_polvec_z_;
  float tauneg_polvec_E_;
 
 public:
  HTTGenAnalysis(std::string const& name);
  virtual ~HTTGenAnalysis();

  virtual int PreAnalysis();
  virtual int Execute(TreeEvent *event);
  virtual int PostAnalysis();
  virtual void PrintInfo();
  

};

}

#endif
