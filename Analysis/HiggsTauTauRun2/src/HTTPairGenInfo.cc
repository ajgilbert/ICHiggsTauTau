#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/interface/HTTPairGenInfo.h"
#include "UserCode/ICHiggsTauTau/interface/PFJet.hh"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include <boost/functional/hash.hpp>
#include "boost/algorithm/string.hpp"
#include "UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/interface/HTTConfig.h"
#include "boost/lexical_cast.hpp"
#include "boost/format.hpp"


namespace ic {

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

  HTTPairGenInfo::HTTPairGenInfo(std::string const& name) : ModuleBase(name), channel_(channel::mt), era_(era::data_2018) {
    ditau_label_ = "ditau";
    fs_ = NULL;
    write_plots_ = false;
    hists_.resize(1);
    ngenjets_=false;
    is_embedded_=false;
  }

  HTTPairGenInfo::~HTTPairGenInfo() {
    ;
  }

  int HTTPairGenInfo::PreAnalysis() {
  if (fs_ && write_plots_) {
      hists_[0] = new Dynamic2DHistoSet(fs_->mkdir("httpairgeninfo"));
      for (unsigned i = 0; i < hists_.size(); ++i) {
        hists_[i]->Create("relpt_vs_drlep_lead", 50, 0, 10, 20, 0, 0.5);
        hists_[i]->Create("relpt_vs_drtau_lead", 50, 0, 10, 20, 0, 0.5);
        hists_[i]->Create("relpt_vs_drlep_sublead", 50, 0, 10, 20, 0, 0.5);
        hists_[i]->Create("relpt_vs_drtau_sublead", 50, 0, 10, 20, 0, 0.5);
      }
    }

    if(is_embedded_&& (channel_==channel::em || channel_==channel::et || channel_==channel::mt || channel_==channel::tt) && (era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP || era_ == era::data_2017 || era_ == era::data_2017UL || era_ == era::data_2018 || era_ == era::data_2018UL)) {

      // Retrieve file with TGraphs of corrections
      std::string file = "input/recoilfits/embed_fake_met_corrections.root";
      embed_corr_ = new TFile(file.c_str());
      embed_corr_->cd();
      std::string graph_name = "rms_corr";
      if(channel_==channel::em) graph_name+="_em";
      if(channel_==channel::et) graph_name+="_et";
      if(channel_==channel::mt) graph_name+="_mt";
      if(channel_==channel::tt) graph_name+="_tt";

      
      if(era_ == era::data_2016 ) graph_name+="_2016";
      if(era_ == era::data_2017) graph_name+="_2017";
      if(era_ == era::data_2018) graph_name+="_2018"; 
      if(era_ == era::data_2016UL_preVFP ) graph_name+="_2016UL_preVFP";
      if(era_ == era::data_2016UL_postVFP) graph_name+="_2016UL_postVFP";
      if(era_ == era::data_2017UL) graph_name+="_2017UL";
      if(era_ == era::data_2018UL) graph_name+="_2018UL";





      gr_met_corr_ = (TGraph*)gDirectory->Get(graph_name.c_str());
      std::cout << "Input TGraph for embedded MET correction: " << graph_name << std::endl;
      //gr_met_corr_up_ = (TGraph*)gDirectory->Get((graph_name+"_up").c_str());
      //gr_met_corr_down_ = (TGraph*)gDirectory->Get((graph_name+"_down").c_str());

      embed_corr_->Close();
    }

    std::string scale_str ="1.";
    std::string scale_str_2 ="1.";

    if(is_embedded_) {
      //if(channel_==channel::tt) {
      //  if(era_ == era::data_2016) scale_str = "(x-y)*0.934 + y*(1.+0.028)";
      //  if(era_ == era::data_2017) scale_str = "(x-y)*0.866 + y*(1.+0.022)";
      //  //if(era_ == era::data_2018) scale_str = "(x-y)*0.891 + y*(1.+0.038)";
      //  if(era_ == era::data_2018) scale_str = "(x-y)*0.883 + y*(1.+0.038)";
      //} else if (channel_==channel::et){
      //  if(era_ == era::data_2016) scale_str = "(x-y)*0.972 + y*(1.+-0.003)";
      //  if(era_ == era::data_2017) scale_str = "(x-y)*0.911 + y*(1.+0.014)";
      //  //if(era_ == era::data_2018) scale_str = "(x-y)*0.960 + y*(1.+0.018)";
      //  if(era_ == era::data_2018) scale_str = "(x-y)*0.952 + y*(1.+0.018)";
      //} else if (channel_==channel::mt){
      //  if(era_ == era::data_2016) scale_str = "(x-y)*0.973 + y*(1.+-0.010)";
      //  if(era_ == era::data_2017) scale_str = "(x-y)*0.906 + y*(1.+0.008)";
      //  //if(era_ == era::data_2018) scale_str = "(x-y)*0.956 + y*(1.+0.016)";
      //  if(era_ == era::data_2018) scale_str = "(x-y)*0.950 + y*(1.+0.016)";
      //} else if (channel_==channel::em){
      //  if(era_ == era::data_2016) scale_str = "(x-y)*1.006 + y*(1.+0.003)";
      //  if(era_ == era::data_2017) scale_str = "(x-y)*0.926 + y*(1.+0.010)";
      //  //if(era_ == era::data_2018) scale_str = "(x-y)*1.008 + y*(1.+0.019)";
      //  if(era_ == era::data_2018) scale_str = "(x-y)*0.998 + y*(1.+0.019)";
      //}

      // so called v3 of corrections where em has no applied correction and is used to "calibrate" other channels
      if(channel_==channel::tt) {
        if(era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) scale_str = "(x-y)*0.929 + y*(1.+0.005)";
        if(era_ == era::data_2017 || era_ == era::data_2017UL) scale_str = "(x-y)*0.935 + y*(1.+-0.002)";
        if(era_ == era::data_2018 || era_ == era::data_2018UL) scale_str = "(x-y)*0.885 + y*(1.+-0.004)";
      } else if (channel_==channel::et){
        if(era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) scale_str = "(x-y)*0.966 + y*(1.+-0.003)";
        if(era_ == era::data_2017 || era_ == era::data_2017UL) scale_str = "(x-y)*0.984 + y*(1.+-0.002)";
        if(era_ == era::data_2018 || era_ == era::data_2018UL) scale_str = "(x-y)*0.954 + y*(1.+-0.007)";
      } else if (channel_==channel::mt){
        if(era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) scale_str = "(x-y)*0.967 + y*(1.+-0.007)";
        if(era_ == era::data_2017 || era_ == era::data_2017UL) scale_str = "(x-y)*0.979 + y*(1.+-0.000)";
        if(era_ == era::data_2018 || era_ == era::data_2018UL) scale_str = "(x-y)*0.952 + y*(1.+0.001)";
      } else if (channel_==channel::em){
        if(era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP) scale_str = "(x-y)*1.000 + y*(1.+0.000)";
        if(era_ == era::data_2017 || era_ == era::data_2017UL) scale_str = "(x-y)*1.000 + y*(1.+0.000)";
        if(era_ == era::data_2018 || era_ == era::data_2018UL) scale_str = "(x-y)*1.000 + y*(1.+0.000)";
      }
    } //else {
      //no njets bins version
      //scale_str  ="(x-y)*min(0.91386 + 0.00323*z, 1.06873 + 0.00014*z) + y";
      //scale_str_2="(x-y)*min(0.91252 + 0.00333*z, 1.06925 + 0.00022*z) + y";
    //}

    std::cout << "smear+scale function for met shift:" << std::endl;
    std::cout << scale_str << std::endl;
    //std::cout << scale_str_2 << std::endl;
    func = new TF2("func",scale_str.c_str());
    //func2 = new TF2("func2",scale_str_2.c_str());

    return 0;
  }

  int HTTPairGenInfo::Execute(TreeEvent *event) {


    std::vector<CompositeCandidate *> const& ditau_vec = event->GetPtrVec<CompositeCandidate>(ditau_label_);
    CompositeCandidate const* ditau = ditau_vec.at(0);

    std::vector<Candidate *> leading_lepton;
    std::vector<Candidate *> subleading_lepton;
    leading_lepton.push_back(ditau->GetCandidate("lepton1"));
    subleading_lepton.push_back(ditau->GetCandidate("lepton2")); 

    std::vector<GenParticle *> const& particles = event->GetPtrVec<GenParticle>("genParticles");
    std::vector<GenParticle *> sel_particles;
    std::vector<GenParticle *> undecayed_taus;
    
    double gen_match_undecayed_1_pt = -1;
    double gen_match_undecayed_2_pt = -1;
    double gen_match_undecayed_1_eta = -1;
    double gen_match_undecayed_2_eta = -1;
    double gen_met=0.;
    GenParticle *gen_match_undecayed_1_;
    GenParticle *gen_match_undecayed_2_;


    TVector3 gen_ip_1(0.,0.,0.);
    TVector3 gen_ip_2(0.,0.,0.);
    bool foundboson=false;
    GenParticle *h = new GenParticle();
    double bosonpT = -1.;

    ROOT::Math::PtEtaPhiEVector neutrinos; 
    ROOT::Math::PtEtaPhiEVector tau_neutrinos; 
    ROOT::Math::PtEtaPhiEVector bosons; 
    ROOT::Math::PtEtaPhiEVector gen_neutral_1;
    ROOT::Math::PtEtaPhiEVector gen_neutral_2;

 
    for (unsigned i=0; i < particles.size(); ++i){
      std::vector<bool> status_flags_start = particles[i]->statusFlags();
      if ( ((abs(particles[i]->pdgid()) == 11 )||(abs(particles[i]->pdgid()) == 13 /*&& particles[i]->status()==1*/)) && particles[i]->pt() > 8. && (status_flags_start[IsPrompt] || status_flags_start[IsDirectPromptTauDecayProduct] /*|| status_flags_start[IsDirectHadronDecayProduct]*/)){
        sel_particles.push_back(particles[i]);
      }
      if ( ((abs(particles[i]->pdgid()) == 12 )||(abs(particles[i]->pdgid()) == 14)||(abs(particles[i]->pdgid()) == 16)) && (status_flags_start[IsPrompt] || status_flags_start[IsDirectPromptTauDecayProduct] /*|| status_flags_start[IsDirectHadronDecayProduct]*/)) neutrinos+=particles[i]->vector();
      if ( ((abs(particles[i]->pdgid()) == 12 )||(abs(particles[i]->pdgid()) == 14)||(abs(particles[i]->pdgid()) == 16)) && (status_flags_start[IsDirectPromptTauDecayProduct])) tau_neutrinos+=particles[i]->vector();
      if(channel_!=channel::zmm&&status_flags_start[IsPrompt] && status_flags_start[IsLastCopy] && abs(particles[i]->pdgid()) == 15) undecayed_taus.push_back(particles[i]);
      if(channel_==channel::zmm&&status_flags_start[IsPrompt] && status_flags_start[IsLastCopy] && abs(particles[i]->pdgid()) == 13) undecayed_taus.push_back(particles[i]);
      if(status_flags_start[IsLastCopy] && (abs(particles[i]->pdgid()) == 23 || abs(particles[i]->pdgid()) == 24 || abs(particles[i]->pdgid()) == 25 || abs(particles[i]->pdgid()) == 35 || abs(particles[i]->pdgid()) == 6) ) {h = particles[i]; foundboson=true; }
      if(status_flags_start[IsLastCopy] && (abs(particles[i]->pdgid()) == 23 || abs(particles[i]->pdgid()) == 24 || abs(particles[i]->pdgid()) == 25 || abs(particles[i]->pdgid()) == 35) ) {bosons+=particles[i]->vector(); }

    }
    bosonpT = bosons.Pt();


    gen_met=neutrinos.Pt();
    double gen_tau_met=tau_neutrinos.Pt();

    event->Add("gen_met",gen_met);
    event->Add("gen_tau_met",gen_tau_met);
    Met *mets = NULL;
    mets = event->GetPtr<Met>("pfMET");

    double fake_met = (mets->vector() - neutrinos).Pt();
    double fake_tau_met = (mets->vector() - tau_neutrinos).Pt();
    event->Add("fake_met",fake_met);
    event->Add("fake_tau_met",fake_tau_met);
    ROOT::Math::PtEtaPhiEVector  fake_met_vec = (mets->vector() - neutrinos);
    ROOT::Math::PtEtaPhiEVector  fake_tau_met_vec = (mets->vector() - tau_neutrinos);
    
    event->Add("fake_met_vec",fake_met_vec);
    event->Add("fake_tau_met_vec",fake_tau_met_vec);
    event->Add("gen_tau_met_vec",tau_neutrinos);
    event->Add("gen_met_vec",neutrinos);


    if(is_embedded_&& (channel_==channel::em || channel_==channel::et || channel_==channel::mt || channel_==channel::tt) && (era_ == era::data_2016 || era_ == era::data_2016UL_preVFP || era_ == era::data_2016UL_postVFP || era_ == era::data_2017 || era_ == era::data_2017UL || era_ == era::data_2018 || era_ == era::data_2018UL)) {
      Met *old_met = new Met(*mets);
      // up uncertainty is uncorrected met so first copy this and add to the event
      event->Add("pfMET_up", old_met);

      double new_px = func->Eval(mets->vector().Px(), tau_neutrinos.Px());
      double new_py = func->Eval(mets->vector().Py(), tau_neutrinos.Py());
      TVector3 vec3(new_px,new_py,0.);
      double new_pt = vec3.Pt();
      double new_phi = vec3.Phi();
      ROOT::Math::PtEtaPhiEVector new_met_vec(new_pt,0., new_phi, new_pt);
      mets->set_vector(new_met_vec);

      //std::cout << "----------" << std::endl;
      //std::cout << tau_neutrinos.Px() << "    " << tau_neutrinos.Py() << std::endl;
      //std::cout << old_met->pt() << "    " << old_met->phi() << "    " << old_met->vector().Px() << "    " << old_met->vector().Py() << std::endl;
      //std::cout << mets->pt() << "    " << mets->phi() << "    " << mets->vector().Px() << "    " << mets->vector().Py() << std::endl;

      // now repeat procedure again for down uncertainty variation (scale applied twice)
      double new_px_down = func->Eval(mets->vector().Px(), tau_neutrinos.Px());
      double new_py_down = func->Eval(mets->vector().Py(), tau_neutrinos.Py());
      TVector3 vec3_down(new_px_down,new_py_down, 0.);
      double new_pt_down = vec3_down.Pt();
      double new_phi_down = vec3_down.Phi();
      ROOT::Math::PtEtaPhiEVector new_met_vec_down(new_pt_down,0., new_phi_down, new_pt_down);
      Met *met_down = new Met(*mets);
      met_down->set_vector(new_met_vec_down);
      event->Add("pfMET_down", met_down);

    }

    if(!is_embedded_ && false) {
      // my version of a recopil corrections like corrections for MC events derived from Z->mumu - seems to have similar problems as recoil corrections so don't apply it but keep code for reference
      double pt = 50.; //for samples with no Z boson set pt equal to average value
      if(bosonpT>=0.) pt = bosonpT;  
      double new_px = func->Eval(mets->vector().Px(), neutrinos.Px(),pt);
      double new_py = func2->Eval(mets->vector().Py(), neutrinos.Py(),pt);
      TVector3 vec3(new_px,new_py,0.);
      double new_pt = vec3.Pt();
      double new_phi = vec3.Phi();
      ROOT::Math::PtEtaPhiEVector new_met_vec(new_pt,0., new_phi, new_pt);
      mets->set_vector(new_met_vec);
    }

    if(undecayed_taus.size()>0){
      gen_match_undecayed_1_pt = undecayed_taus[0]->pt();
      gen_match_undecayed_1_eta = undecayed_taus[0]->eta();
      gen_match_undecayed_1_=undecayed_taus[0];
      event->Add("gen_match_undecayed_1", gen_match_undecayed_1_);
    }
    if(undecayed_taus.size()>1){
      gen_match_undecayed_2_pt = undecayed_taus[1]->pt();
      gen_match_undecayed_2_eta = undecayed_taus[1]->eta();
      gen_match_undecayed_2_=undecayed_taus[1];
      event->Add("gen_match_undecayed_2", gen_match_undecayed_2_);
    }

    int tauFlag1 = 0;
    int tauFlag2 = 0;
    
    std::vector<GenJet> gen_taus = BuildTauJets(particles, false,true);
    std::vector<GenJet *> gen_taus_ptr;
    for (auto & x : gen_taus) gen_taus_ptr.push_back(&x);
    ic::erase_if(gen_taus_ptr, !boost::bind(MinPtMaxEta, _1, 15.0, 999.));
 
 
    std::vector<std::pair<Candidate*, GenParticle*> > leading_lepton_match = MatchByDR(leading_lepton, sel_particles, 0.2, true, true);
    std::vector<std::pair<Candidate*, GenJet*> > leading_tau_match  = MatchByDR(leading_lepton, gen_taus_ptr, 0.2, true, true);

    std::vector<std::pair<Candidate*, GenParticle*> > subleading_lepton_match = MatchByDR(subleading_lepton, sel_particles, 0.2, true, true);
    std::vector<std::pair<Candidate*, GenJet*> > subleading_tau_match  = MatchByDR(subleading_lepton, gen_taus_ptr, 0.2, true, true);

   double gen_nu_p_1=-9999.;
   double gen_nu_p_2=-9999.;
   double gen_nu_px_1=-9999.;
   double gen_nu_px_2=-9999.;
   double gen_nu_py_1=-9999.;
   double gen_nu_py_2=-9999.;
   double gen_nu_pz_1=-9999.;
   double gen_nu_pz_2=-9999.;
   double gen_nu_phi_1=-9999.;
   double gen_nu_phi_2=-9999.;
   double gen_nu_eta_1=-9999.;
   double gen_nu_eta_2=-9999.;

   double gen_vis_p_1=-9999.;
   double gen_vis_p_2=-9999.;
   double gen_vis_px_1=-9999.;
   double gen_vis_px_2=-9999.;
   double gen_vis_py_1=-9999.;
   double gen_vis_py_2=-9999.;
   double gen_vis_pz_1=-9999.;
   double gen_vis_pz_2=-9999.;

   double gen_vis_phi_1=-9999.;
   double gen_vis_phi_2=-9999.;
   double gen_vis_eta_1=-9999.;
   double gen_vis_eta_2=-9999.;
   double gen_vis_E_1=-9999.;
   double gen_vis_E_2=-9999.;

   double gen_neutral_p_1=-9999.;
   double gen_neutral_p_2=-9999.;
   double gen_neutral_px_1=-9999.;
   double gen_neutral_px_2=-9999.;
   double gen_neutral_py_1=-9999.;
   double gen_neutral_py_2=-9999.;
   double gen_neutral_pz_1=-9999.;
   double gen_neutral_pz_2=-9999.;
   double gen_neutral_phi_1=-9999.;
   double gen_neutral_phi_2=-9999.;
   double gen_neutral_eta_1=-9999.;
   double gen_neutral_eta_2=-9999.;
   double gen_neutral_E_1=-9999.;
   double gen_neutral_E_2=-9999.;

   mcorigin gen_match_1 = mcorigin::fake;
   mcorigin gen_match_2 = mcorigin::fake;
   int leptonsize = leading_lepton_match.size();
   int tausize = leading_tau_match.size();
   double gen_match_1_pt = -1;
   double gen_match_2_pt = -1;

   if(leptonsize!=0&&tausize!=0){
     DR(leading_lepton_match.at(0).first,leading_lepton_match.at(0).second) < DR(leading_tau_match.at(0).first,leading_tau_match.at(0).second) ? tausize=0 : leptonsize = 0;
   }

   // get gen IP for candidate 1
   if ((leptonsize==0&&tausize==0) && !foundboson) {
     gen_ip_1.SetXYZ(0.,0.,0.);
   } else if(leptonsize!=0  && foundboson) {
     GenParticle *t = leading_lepton_match.at(0).second;
     gen_ip_1 = GenIP(h,t);
   } else if (tausize!=0 && foundboson) {
     std::vector<std::size_t> daughter_ids = leading_tau_match.at(0).second->constituents();
     bool found=false;
     for (auto id: daughter_ids) {
       for(auto p: particles) {
         if(p->id() == id && p->charge()!=0 && !found) {
           gen_ip_1 = GenIP(h,p);
           found=true;
           break; 
         }
         if(p->id() == id && p->charge()==0 && fabs(p->pdgid())==22) {
           gen_neutral_1 += p->vector();
           break;
         }
       }
     } 
   }


   if(leptonsize==0&&tausize==0) gen_match_1 = mcorigin::fake;
   if(leptonsize!=0) {
      std::vector<bool> status_flags = leading_lepton_match.at(0).second->statusFlags();
      gen_match_1_pt = leading_lepton_match.at(0).second->pt();
//      leading_lepton_match_DR = DR(leading_lepton_match.at(0).first,leading_lepton_match.at(0).second);
      if(fs_ && write_plots_){ 
       hists_[0]->Fill("relpt_vs_drlep_lead",fabs(leading_lepton_match.at(0).second->pt()-leading_lepton_match.at(0).first->pt())/leading_lepton_match.at(0).first->pt(),DR(leading_lepton_match.at(0).first,leading_lepton_match.at(0).second));
      }
      if(status_flags[IsPrompt]){
        if(abs(leading_lepton_match.at(0).second->pdgid())==11){
          gen_match_1 = mcorigin::promptE;
         } else gen_match_1 = mcorigin::promptMu;
       }
       if(status_flags[IsDirectPromptTauDecayProduct]){
        if(abs(leading_lepton_match.at(0).second->pdgid())==11){
          gen_match_1 = mcorigin::tauE;
         } else gen_match_1 = mcorigin::tauMu;
        }
       if(status_flags[IsDirectHadronDecayProduct]){
         if(abs(leading_lepton_match.at(0).second->pdgid())==11){
           gen_match_1 = mcorigin::hadE;
          } else gen_match_1 = mcorigin::hadMu;
        } 
      }
     if(tausize!=0){
       if(fs_ && write_plots_){
       hists_[0]->Fill("relpt_vs_drtau_lead",fabs(leading_tau_match.at(0).second->pt()-leading_tau_match.at(0).first->pt())/leading_tau_match.at(0).first->pt(),DR(leading_tau_match.at(0).first,leading_tau_match.at(0).second));
       }
       gen_match_1_pt = leading_tau_match.at(0).second->pt();
  //     leading_lepton_match_DR = DR(leading_tau_match.at(0).first,leading_tau_match.at(0).second);
       gen_match_1 = mcorigin::tauHad;
       tauFlag1 = leading_tau_match.at(0).second->flavour();

       gen_nu_p_1=leading_tau_match.at(0).second->nu_vector().P();
       gen_nu_px_1=leading_tau_match.at(0).second->nu_vector().Px();
       gen_nu_py_1=leading_tau_match.at(0).second->nu_vector().Py();
       gen_nu_pz_1=leading_tau_match.at(0).second->nu_vector().Pz();
       gen_nu_phi_1=leading_tau_match.at(0).second->nu_vector().Phi();
       gen_nu_eta_1=leading_tau_match.at(0).second->nu_vector().Rapidity();

       gen_vis_p_1=leading_tau_match.at(0).second->vector().P();
       gen_vis_px_1=leading_tau_match.at(0).second->vector().Px();
       gen_vis_py_1=leading_tau_match.at(0).second->vector().Py();
       gen_vis_pz_1=leading_tau_match.at(0).second->vector().Pz();
       gen_vis_E_1=leading_tau_match.at(0).second->vector().E();
       gen_vis_phi_1=leading_tau_match.at(0).second->vector().Phi();
       gen_vis_eta_1=leading_tau_match.at(0).second->vector().Rapidity();

       gen_neutral_p_1=gen_neutral_1.P();
       gen_neutral_px_1=gen_neutral_1.Px();
       gen_neutral_py_1=gen_neutral_1.Py();
       gen_neutral_pz_1=gen_neutral_1.Pz();
       gen_neutral_E_1=gen_neutral_1.E();
       gen_neutral_phi_1=gen_neutral_1.Phi();
       gen_neutral_eta_1=gen_neutral_1.Rapidity();

      }
   
//Now for subleading lepton:
   leptonsize = subleading_lepton_match.size();
   tausize = subleading_tau_match.size();

   if(leptonsize!=0&&tausize!=0){
     DR(subleading_lepton_match.at(0).first,subleading_lepton_match.at(0).second) < DR(subleading_tau_match.at(0).first,subleading_tau_match.at(0).second) ? tausize=0 : leptonsize = 0;
   }

   // get gen IP for candidate 2
   if ((leptonsize==0&&tausize==0) && !foundboson) {
     gen_ip_2.SetXYZ(0.,0.,0.);
   } else if(leptonsize!=0 && foundboson) {
     GenParticle *t = subleading_lepton_match.at(0).second;
     gen_ip_2 = GenIP(h,t);
   } else if (tausize!=0 && foundboson) {
     std::vector<std::size_t> daughter_ids = subleading_tau_match.at(0).second->constituents();
     bool found=false;
     for (auto id: daughter_ids) {
       for(auto p: particles) {
         if(p->id() == id && p->charge()!=0 && !found) {
           gen_ip_2 = GenIP(h,p);
           found=true;
           break;
         }
         if(p->id() == id && p->charge()==0 && fabs(p->pdgid())==22) {
           gen_neutral_2 += p->vector();
           break;
         }
       }
     }
   }


   if(leptonsize==0&&tausize==0) gen_match_2 = mcorigin::fake;
   if(leptonsize!=0) {
      std::vector<bool> status_flags = subleading_lepton_match.at(0).second->statusFlags();
      gen_match_2_pt = subleading_lepton_match.at(0).second->pt();
    //  subleading_lepton_match_DR = DR(subleading_lepton_match.at(0).first,subleading_lepton_match.at(0).second);
      if(fs_ && write_plots_){
       hists_[0]->Fill("relpt_vs_drlep_sublead",fabs(subleading_lepton_match.at(0).second->pt()-subleading_lepton_match.at(0).first->pt())/subleading_lepton_match.at(0).first->pt(),DR(subleading_lepton_match.at(0).first,subleading_lepton_match.at(0).second));
      }
      if(status_flags[IsPrompt]){
        if(abs(subleading_lepton_match.at(0).second->pdgid())==11){
          gen_match_2 = mcorigin::promptE;
         } else gen_match_2 = mcorigin::promptMu;
       }
       if(status_flags[IsDirectPromptTauDecayProduct]){
        if(abs(subleading_lepton_match.at(0).second->pdgid())==11){
          gen_match_2 = mcorigin::tauE;
         } else gen_match_2 = mcorigin::tauMu;
        }
       if(status_flags[IsDirectHadronDecayProduct]){
         if(abs(subleading_lepton_match.at(0).second->pdgid())==11){
           gen_match_2 = mcorigin::hadE;
          } else gen_match_2 = mcorigin::hadMu;
        } 
      }
     if(tausize!=0){
       if(fs_&&write_plots_){
       hists_[0]->Fill("relpt_vs_drtau_sublead",fabs(subleading_tau_match.at(0).second->pt()-subleading_tau_match.at(0).first->pt())/subleading_tau_match.at(0).first->pt(),DR(subleading_tau_match.at(0).first,subleading_tau_match.at(0).second));
       }
       gen_match_2_pt = subleading_tau_match.at(0).second->pt();
      // subleading_lepton_match_DR = DR(subleading_tau_match.at(0).first,subleading_tau_match.at(0).second);
       gen_match_2 = mcorigin::tauHad;
       tauFlag2 = subleading_tau_match.at(0).second->flavour(); 
       gen_nu_p_2=subleading_tau_match.at(0).second->nu_vector().P();
       gen_nu_px_2=subleading_tau_match.at(0).second->nu_vector().Px();
       gen_nu_py_2=subleading_tau_match.at(0).second->nu_vector().Py();
       gen_nu_pz_2=subleading_tau_match.at(0).second->nu_vector().Pz();
       gen_nu_phi_2=subleading_tau_match.at(0).second->nu_vector().Phi();
       gen_nu_eta_2=subleading_tau_match.at(0).second->nu_vector().Rapidity();

       gen_vis_p_2=subleading_tau_match.at(0).second->vector().P();
       gen_vis_px_2=subleading_tau_match.at(0).second->vector().Px();
       gen_vis_py_2=subleading_tau_match.at(0).second->vector().Py();
       gen_vis_pz_2=subleading_tau_match.at(0).second->vector().Pz();
       gen_vis_E_2=subleading_tau_match.at(0).second->vector().E();
       gen_vis_phi_2=subleading_tau_match.at(0).second->vector().Phi();
       gen_vis_eta_2=subleading_tau_match.at(0).second->vector().Rapidity();

       gen_neutral_p_2=gen_neutral_2.P();
       gen_neutral_px_2=gen_neutral_2.Px();
       gen_neutral_py_2=gen_neutral_2.Py();
       gen_neutral_pz_2=gen_neutral_2.Pz();
       gen_neutral_E_2=gen_neutral_2.E();
       gen_neutral_phi_2=gen_neutral_2.Phi();
       gen_neutral_eta_2=gen_neutral_2.Rapidity();

      }

   if(gen_match_1 == mcorigin::tauHad) event->Add("leading_gen_tau", new ic::GenJet(*(leading_tau_match.at(0).second)));
   if(gen_match_2 == mcorigin::tauHad) event->Add("subleading_gen_tau", new ic::GenJet(*(subleading_tau_match.at(0).second)));

   int gen_match_1_mother=-9999;
   int gen_match_2_mother=-9999;

   if(leading_lepton_match.size()>0 && leading_lepton_match.at(0).second->mothers().size()>0) gen_match_1_mother = particles[leading_lepton_match.at(0).second->mothers()[0]]->pdgid();
   if(subleading_lepton_match.size()>0 && subleading_lepton_match.at(0).second->mothers().size()>0) gen_match_2_mother = particles[subleading_lepton_match.at(0).second->mothers()[0]]->pdgid();

   int gen_match_1_q=-9999;
   int gen_match_2_q=-9999;

   if(leading_lepton_match.size()>0) gen_match_1_q = leading_lepton_match.at(0).second->charge();
   if(subleading_lepton_match.size()>0) gen_match_2_q = subleading_lepton_match.at(0).second->charge();

   event->Add("gen_match_1_q",gen_match_1_q);
   event->Add("gen_match_2_q",gen_match_2_q);

   event->Add("gen_match_1_mother",gen_match_1_mother);
   event->Add("gen_match_2_mother",gen_match_2_mother);

   event->Add("gen_match_1",gen_match_1);
   event->Add("gen_match_2",gen_match_2);
   event->Add("gen_match_1_pt", gen_match_1_pt);
   event->Add("gen_match_2_pt", gen_match_2_pt);
   event->Add("gen_match_undecayed_1_pt", gen_match_undecayed_1_pt);
   event->Add("gen_match_undecayed_2_pt", gen_match_undecayed_2_pt);
   event->Add("gen_match_undecayed_1_eta", gen_match_undecayed_1_eta);
   event->Add("gen_match_undecayed_2_eta", gen_match_undecayed_2_eta);
   event->Add("tauFlag1", tauFlag1);
   event->Add("tauFlag2", tauFlag2);

   event->Add("gen_ip_1", gen_ip_1);

   event->Add("gen_ip_2", gen_ip_2);

   event->Add("gen_nu_p_1",gen_nu_p_1);
   event->Add("gen_nu_px_1",gen_nu_px_1);
   event->Add("gen_nu_py_1",gen_nu_py_1);
   event->Add("gen_nu_pz_1",gen_nu_pz_1);
   event->Add("gen_nu_phi_1",gen_nu_phi_1);
   event->Add("gen_nu_eta_1",gen_nu_eta_1);

   event->Add("gen_nu_p_2",gen_nu_p_2);
   event->Add("gen_nu_px_2",gen_nu_px_2);
   event->Add("gen_nu_py_2",gen_nu_py_2);
   event->Add("gen_nu_pz_2",gen_nu_pz_2);
   event->Add("gen_nu_phi_2",gen_nu_phi_2);
   event->Add("gen_nu_eta_2",gen_nu_eta_2);

   event->Add("gen_vis_p_1",  gen_vis_p_1);
   event->Add("gen_vis_px_1",  gen_vis_px_1);
   event->Add("gen_vis_py_1",  gen_vis_py_1);
   event->Add("gen_vis_pz_1",  gen_vis_pz_1);
   event->Add("gen_vis_E_1",  gen_vis_E_1);
   event->Add("gen_vis_phi_1",gen_vis_phi_1);
   event->Add("gen_vis_eta_1",gen_vis_eta_1);

   event->Add("gen_vis_p_2",  gen_vis_p_2);
   event->Add("gen_vis_px_2",  gen_vis_px_2);
   event->Add("gen_vis_py_2",  gen_vis_py_2);
   event->Add("gen_vis_pz_2",  gen_vis_pz_2);
   event->Add("gen_vis_E_2",  gen_vis_E_2);
   event->Add("gen_vis_phi_2",gen_vis_phi_2);
   event->Add("gen_vis_eta_2",gen_vis_eta_2);

   event->Add("gen_neutral_p_1",  gen_neutral_p_1);
   event->Add("gen_neutral_px_1",  gen_neutral_px_1);
   event->Add("gen_neutral_py_1",  gen_neutral_py_1);
   event->Add("gen_neutral_pz_1",  gen_neutral_pz_1);
   event->Add("gen_neutral_E_1",  gen_neutral_E_1);
   event->Add("gen_neutral_phi_1",gen_neutral_phi_1);
   event->Add("gen_neutral_eta_1",gen_neutral_eta_1);

   event->Add("gen_neutral_p_2",  gen_neutral_p_2);
   event->Add("gen_neutral_px_2",  gen_neutral_px_2);
   event->Add("gen_neutral_py_2",  gen_neutral_py_2);
   event->Add("gen_neutral_pz_2",  gen_neutral_pz_2);
   event->Add("gen_neutral_E_2",  gen_neutral_E_2);
   event->Add("gen_neutral_phi_2",gen_neutral_phi_2);
   event->Add("gen_neutral_eta_2",gen_neutral_eta_2);



    if(ngenjets_){
      //Get gen-jets collection, filter Higgs decay products and add Njets variable to event
      std::vector<ic::GenJet*> gen_jets = event->GetPtrVec<ic::GenJet>("genJets");

      // erase jets matched to our taus or leptons
      for(unsigned i=0; i<gen_jets.size(); ++i){
        ic::GenJet *genjet = gen_jets[i];
        bool MatchedToPrompt = false;
        for(unsigned j=0; j<sel_particles.size(); ++j){
         if(DRLessThan(std::make_pair(genjet, sel_particles[j]),0.5)) MatchedToPrompt = true;
        }
        for(unsigned j=0; j<gen_taus_ptr.size(); ++j){
         if(DRLessThan(std::make_pair(genjet, gen_taus_ptr[j]),0.5)) MatchedToPrompt = true;
        }
        //remove jets that are matched to Higgs decay products
        if(MatchedToPrompt) gen_jets.erase (gen_jets.begin()+i);
      }

      double genjetpt_1=-1., genjetpt_2=-1.;
      double genjeteta_1=-1., genjeteta_2=-1.;
      if(gen_jets.size()>0) {
        genjetpt_1 = gen_jets[0]->pt();
        genjeteta_1 = fabs(gen_jets[0]->eta());
      }
      if(gen_jets.size()>1) {
        genjetpt_2 = gen_jets[1]->pt();
        genjeteta_2 = fabs(gen_jets[1]->eta());
      }
      event->Add("genjetpt_1", genjetpt_1);
      event->Add("genjetpt_2", genjetpt_2);
      event->Add("genjeteta_1", genjeteta_1);
      event->Add("genjeteta_2", genjeteta_2);

      ic::erase_if(gen_jets,!boost::bind(MinPtMaxEta, _1, 20.0, 5.0));
      unsigned ngenjets20 = gen_jets.size();
      event->Add("ngenjets20", ngenjets20);   
      double gen_sjdphi_ = -999; 
      double gen_mjj_=-9999;
      if (gen_jets.size() >= 2) {
        gen_mjj_ = (gen_jets[0]->vector()+gen_jets[1]->vector()).M();
        if(gen_jets[0]->eta() > gen_jets[1]->eta())
          gen_sjdphi_ =  ROOT::Math::VectorUtil::DeltaPhi(gen_jets[0]->vector(), gen_jets[1]->vector());
        
        else
          gen_sjdphi_ =  ROOT::Math::VectorUtil::DeltaPhi(gen_jets[1]->vector(), gen_jets[0]->vector());
        
      }
      ic::erase_if(gen_jets,!boost::bind(MinPtMaxEta, _1, 30.0, 4.7));
      unsigned ngenjets = gen_jets.size();
      event->Add("ngenjets", ngenjets);
      event->Add("gen_sjdphi", gen_sjdphi_);
      event->Add("gen_mjj", gen_mjj_);
    }

    return 0;
  }
  int HTTPairGenInfo::PostAnalysis() {
    return 0;
  }

  void HTTPairGenInfo::PrintInfo() {
    ;
  }

}
