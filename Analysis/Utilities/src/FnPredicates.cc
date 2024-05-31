#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPredicates.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/FnPairs.h"
#include "UserCode/ICHiggsTauTau//interface/city.h"
#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/ElectronEffectiveArea.h"

#include <boost/lexical_cast.hpp>
#include "TRandom2.h"
#include "TMath.h"
#include <utility>
#include <set>

namespace ic {

  bool GoodVertex(Vertex const* vertex) {
    if (     vertex->ndof() > 4 
          && fabs(vertex->vz()) < 24.0 
          && vertex->point().rho() < 2.0) {
      return true;
    } else {
      return false; 
    }
  }

  bool GenParticleInMassBand(GenParticle const* part, int pdgid, double m_low, double m_high) {
    return (abs(part->pdgid()) == pdgid && part->vector().M() > m_low && part->vector().M() < m_high);
  }

  bool IsFilterMatched(Candidate const* cand, std::vector<TriggerObject *> const& objs, std::string const& filter, double const& max_dr) {
    std::size_t hash = CityHash64(filter);
    for (unsigned i = 0; i < objs.size(); ++i) {
      std::vector<std::size_t> const& labels = objs[i]->filters();
      if (std::find(labels.begin(),labels.end(), hash) == labels.end()) continue;
      if (DR(cand, objs[i]) < max_dr) return true;
    }
    return false;
  }

 std::pair<bool,unsigned> IsFilterMatchedWithIndex(Candidate const* cand, std::vector<TriggerObject*> const& objs, std::string const& filter, double const& max_dr){
    std::size_t hash = CityHash64(filter);
    for (unsigned i = 0; i < objs.size(); ++i) {
      std::vector<std::size_t> const& labels = objs[i]->filters();
      if (std::find(labels.begin(),labels.end(), hash) == labels.end()) continue;
      if (DR(cand, objs[i]) < max_dr){
         const unsigned index_val = i;
         return std::make_pair(true,index_val);
      }
    }
    return std::make_pair(false,0);
  }
  
  std::pair<bool,std::vector<unsigned>> IsFilterMatchedWithMultipleIndexs(Candidate const* cand, std::vector<TriggerObject*> const& objs, std::string const& filter, double const& max_dr){
    std::size_t hash = CityHash64(filter);
    std::vector<unsigned> index_vals;
    bool matched = false;
    for (unsigned i = 0; i < objs.size(); ++i) {
      std::vector<std::size_t> const& labels = objs[i]->filters();
      if (std::find(labels.begin(),labels.end(), hash) == labels.end()) continue;
      if (DR(cand, objs[i]) < max_dr){
         index_vals.push_back(i);
         matched = true;
      }
    }
    return std::make_pair(matched,index_vals);
  }
  
  std::set<std::size_t> ListTriggerFilters(std::vector<TriggerObject *> const& objs){
    std::set<std::size_t> filters;  
    for (unsigned i = 0; i < objs.size(); ++i) {
      std::vector<std::size_t> const& labels = objs[i]->filters();
      for(unsigned j=0; j<labels.size(); j++){
        filters.insert(labels[j]);    
      }
    }  
    return filters;
  }
  
  bool VertexDz(Tau const* cand, double const& vertexZ) {
    return ( fabs(cand->vz() - vertexZ)==0) ; 
  }

  bool MinPtMaxEta(Candidate const* cand, double const& minPt, double const& maxEta) {
    return ( (cand->pt() > minPt) && (std::fabs(cand->eta()) < maxEta) ); 
  }

  bool MinPtMaxSCEta(Electron const* cand, double const& minPt, double const& maxEta) {
    return ( (cand->pt() > minPt) && (std::fabs(cand->sc_eta()) < maxEta) ); 
  }

  bool EtaOutsideRange(Candidate const* cand, double const& minEtaCut, double const& maxEtaCut) {
    return ( (std::fabs(cand->eta()) < minEtaCut) || (std::fabs(cand->eta()) > maxEtaCut) ); 
  }


  bool CSVMediumWP(Jet const* cand, double const& minDiscri) {
    return (cand->GetBDiscriminator("combinedSecondaryVertexBJetTags") > minDiscri);
  }



  bool InEcalGap(Electron const* electron) {
    double scEta = fabs(electron->sc_eta());
    return (scEta >= 1.4442 && scEta <= 1.566);
  }

  bool MinPtMaxEtaSC(SuperCluster const* cand, double const& minPt, double const& maxEta) {
    return ( (cand->pt() > minPt) && (std::fabs(cand->eta()) < maxEta) ); 
  }

  bool InEcalGapSC(SuperCluster const* electron) {
    double scEta = fabs(electron->eta());
    return (scEta >= 1.4442 && scEta <= 1.566);
  }

  bool InEcalGapCand(Candidate const* cand) {
    double scEta = fabs(cand->eta());
    return (scEta >= 1.4442 && scEta <= 1.566);
  }

  bool MaxEtaSC(Electron const* electron, double const& maxEta) {
    double scEta = fabs(electron->sc_eta());
    return (scEta < maxEta);
  }

  bool PFJetID(PFJet const* jet) {
    int n_pu = jet->charged_multiplicity() + jet->neutral_multiplicity() + jet->HF_had_multiplicity() + jet->HF_em_multiplicity();
    double eta = fabs(jet->eta());
    bool result = false;

    double neutralFrac = ( jet->neutral_had_energy() + jet->HF_had_energy() ) / jet->uncorrected_energy();

    if (eta < 2.4) {
      result = neutralFrac   < 0.99
      && jet->neutral_em_energy_frac()    < 0.99
            && n_pu                             > 0
            && jet->charged_had_energy_frac()   > 0.0
            && jet->charged_multiplicity()      > 0
            && jet->charged_em_energy_frac()    < 0.99;
    } else {
      result = neutralFrac   < 0.99
            && jet->neutral_em_energy_frac()    < 0.99
            && n_pu                             > 0;
    }
    return result;
  }

  bool PFJetIDNoHFCut(PFJet const* jet) {
    int n_pu = jet->charged_multiplicity() + jet->neutral_multiplicity() + jet->HF_had_multiplicity() + jet->HF_em_multiplicity();
    double eta = fabs(jet->eta());
    bool result = false;

    double neutralFrac = jet->neutral_had_energy() / jet->uncorrected_energy();

    if (eta < 2.4) {
      result = neutralFrac   < 0.99
      && jet->neutral_em_energy_frac()    < 0.99
            && n_pu                             > 0
            && jet->charged_had_energy_frac()   > 0.0
            && jet->charged_multiplicity()      > 0
            && jet->charged_em_energy_frac()    < 0.99;
    } else {
      result = neutralFrac   < 0.99
            && jet->neutral_em_energy_frac()    < 0.99
            && n_pu                             > 0;
    }
    return result;
  }

  bool PFJetID2015(PFJet const* jet) {
    double eta = fabs(jet->eta());
    bool result = false;

    double neutralFrac = jet->neutral_had_energy() / jet->uncorrected_energy();
//    int n_pf = jet->charged_multiplicity() + jet->neutral_multiplicity() + jet->HF_had_multiplicity() + jet->HF_em_multiplicity();

    if (eta <= 2.4) {
      result = neutralFrac   < 0.99
      && jet->neutral_em_energy_frac()    < 0.99
            && jet->charged_multiplicity()+jet->neutral_multiplicity() > 1
            && jet->charged_had_energy_frac()   > 0.0
            && jet->charged_multiplicity()      > 0
            && jet->charged_em_energy_frac()    < 0.99;
    } else if(eta<=3.0){
      result = neutralFrac   < 0.99
            && jet->neutral_em_energy_frac()    < 0.99
            && jet->charged_multiplicity()+jet->neutral_multiplicity() > 1;
    }
    else{
      result = jet->neutral_em_energy_frac()    < 0.90
	    && jet->neutral_multiplicity()>10;
    }
    return result;
  }

  bool PFJetID2016(PFJet const* jet) {
    double eta = fabs(jet->eta());
    bool result = false;
    
    double neutralFrac = jet->neutral_had_energy() / jet->uncorrected_energy();
    //    int n_pf = jet->charged_multiplicity() + jet->neutral_multiplicity() + jet->HF_had_multiplicity() + jet->HF_em_multiplicity();
    
    if (eta <= 2.4) {
      result = neutralFrac   < 0.99
      && jet->neutral_em_energy_frac()    < 0.99
      && jet->charged_multiplicity()+jet->neutral_multiplicity() > 1
      && jet->charged_had_energy_frac()   > 0.0
      && jet->charged_multiplicity()      > 0
      && jet->charged_em_energy_frac()    < 0.99;
    } else if (eta <= 2.7){
      result = neutralFrac < 0.90 // 0.99 (loose)
      && jet->neutral_em_energy_frac()   < 0.90 // 0.99 (loose)
      && jet->charged_multiplicity()+jet->neutral_multiplicity() > 1;
    } else if(eta<=3.0){
      result = jet->neutral_em_energy_frac() > 0.01
      && neutralFrac                         < 0.98
      && jet->neutral_multiplicity()         > 2;
    }
    else{
      result = jet->neutral_em_energy_frac() < 0.90
      && jet->neutral_multiplicity()         > 10;
    }
    return result;
  } 
  
  bool PFJetID2017(PFJet const* jet) {
    double eta = fabs(jet->eta());
    bool result = false;
    
    double neutralFrac = jet->neutral_had_energy() / jet->uncorrected_energy();
    
    if (eta <= 2.4) {
      result = neutralFrac   < 0.90
      && jet->neutral_em_energy_frac()    < 0.90
      && jet->charged_multiplicity()+jet->neutral_multiplicity() > 1
      && jet->charged_had_energy_frac()   > 0.0
      && jet->charged_multiplicity()      > 0;
    } else if (eta <= 2.7){
      result = neutralFrac < 0.90
      && jet->neutral_em_energy_frac()   < 0.90
      && jet->charged_multiplicity()+jet->neutral_multiplicity() > 1;
    } else if(eta<=3.0){
      result = jet->neutral_em_energy_frac()    < 0.99
      &&jet->neutral_em_energy_frac()    > 0.02
      && jet->neutral_multiplicity() > 2;
    }
    else{
      result = jet->neutral_em_energy_frac()    < 0.90
      && neutralFrac>0.02
      && jet->neutral_multiplicity()>10;
    }
    return result;
  } 

  bool PFJetID2018(PFJet const* jet) {
    // Applying JetID criteria (not JetIDLepVeto) from 
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
    double eta = fabs(jet->eta());
    bool result = false;
    
    double neutralFrac = jet->neutral_had_energy() / jet->uncorrected_energy();
    
    if (eta <= 2.6) {
      result = neutralFrac   < 0.90
      && jet->neutral_em_energy_frac()    < 0.90
      && jet->charged_multiplicity()+jet->neutral_multiplicity() > 1.
      && jet->charged_had_energy_frac()   > 0.
      && jet->charged_multiplicity()      > 0.;
    } else if (eta <= 2.7){
      result = neutralFrac < 0.90
      && jet->neutral_em_energy_frac()   < 0.99
      && jet->charged_multiplicity()     > 0.;
    } else if(eta<=3.0){
      result = jet->neutral_em_energy_frac() < 0.99
      && jet->neutral_em_energy_frac()       > 0.02
      && jet->neutral_multiplicity()         > 2;
    }
    else{
      result = jet->neutral_em_energy_frac() < 0.90
      && neutralFrac                         > 0.02
      && jet->neutral_multiplicity()         >10;
    }
    return result;
  } 

  bool PUJetID(PFJet const* jet, bool is_2012) {
    // Pt2030_Loose   = cms.vdouble(-0.80,-0.85,-0.84,-0.85),
    // Pt3050_Loose   = cms.vdouble(-0.80,-0.74,-0.68,-0.77)
    // #4 Eta Categories  0-2.5 2.5-2.75 2.75-3.0 3.0-5.0
    double abs_eta = fabs(jet->eta());
    double pt = jet->pt();
    double pu_id_mva_value = jet->pu_id_mva_value();
    if (is_2012) {
      if (pt > 20. && pt <= 30) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.80);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.85);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.84);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.85);
            } else return true;
      } else if (pt > 30.) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.80);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.74);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.68);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.77);
            } else return true;
      } else return true;
    } else {
      if (pt > 20. && pt <= 30) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.40);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.85);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.7);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.6);
            } else return true;
      } else if (pt > 30.) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.40);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.85);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.7);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.6);
            } else return true;
      } else return true;
    }
  }


  bool PileupJetID(PFJet const* jet, unsigned training, bool doTight, bool doLoose) {
    double abs_eta = fabs(jet->eta());
    double pt = jet->pt();
    double pu_id_mva_value = jet->pu_id_mva_value();
    if (training == 0) {
      if (pt > 20. && pt <= 30.) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.80);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.85);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.84);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.85);
            } else return true;
      } else if (pt > 30. && pt <= 50.) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.80);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.74);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.68);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.77);
            } else return true;
      } else return true;
    } else if (training == 1) {
      if (pt > 20. && pt <= 30.) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.40);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.85);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.7);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.6);
            } else return true;
      } else if (pt > 30. && pt <= 50.) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.40);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.85);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.7);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.6);
            } else return true;
      } else return true;
    } else if (training == 2) {
      if (pt > 20. && pt <= 30.) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.63);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.60);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.55);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.45);
            } else return true;
      } else if (pt > 30. && pt <= 50.) {
            if (abs_eta < 2.5) {
                  return (pu_id_mva_value > -0.63);
            } else if (abs_eta < 2.75) {
                  return (pu_id_mva_value > -0.60);
            } else if (abs_eta < 3.0) {
                  return (pu_id_mva_value > -0.55);
            } else if (abs_eta < 5.0) {
                  return (pu_id_mva_value > -0.45);
            } else return true;
      } else return true;
    } else if (training == 4) {
      //reference: page 16 of https://indico.cern.ch/event/559594/contributions/2257924/attachments/1317046/1973307/PUID_JMAR_2016_07_26_v1.pdf
      if (pt <= 30.) {
        if (abs_eta < 2.5) {
          if (doTight) return (pu_id_mva_value > 0.69);
          else if (doLoose) return (pu_id_mva_value > -0.97);
          else return (pu_id_mva_value > 0.18);
        } 
        else if (abs_eta < 2.75) {
          if (doTight) return (pu_id_mva_value > -0.35);
          else if (doLoose) return (pu_id_mva_value > -0.68);
          else return (pu_id_mva_value > -0.55);
        }
        else if (abs_eta < 3.0) {
          if (doTight) return (pu_id_mva_value > -0.26);
          else if (doLoose) return (pu_id_mva_value > -0.53);
          else return (pu_id_mva_value > -0.42);
        }
        else if (abs_eta < 5.0) {
          if (doTight) return (pu_id_mva_value > -0.21);
          else if (doLoose) return (pu_id_mva_value > -0.47);
          else return (pu_id_mva_value > -0.36);
        }
        else return true;
      }
      else if (pt > 30. && pt <= 50.) {
        if (abs_eta < 2.5) {
          if (doTight) return (pu_id_mva_value > 0.86);
          else if (doLoose) return (pu_id_mva_value > -0.89);
          else return (pu_id_mva_value > 0.61);
        } 
        else if (abs_eta < 2.75) {
          if (doTight) return (pu_id_mva_value > -0.10);
          else if (doLoose) return (pu_id_mva_value > -0.52);
          else return (pu_id_mva_value > -0.35);
        }
        else if (abs_eta < 3.0) {
          if (doTight) return (pu_id_mva_value > -0.05);
          else if (doLoose) return (pu_id_mva_value > -0.38);
          else return (pu_id_mva_value > -0.23);
        }
        else if (abs_eta < 5.0) {
          if (doTight) return (pu_id_mva_value > -0.01);
          else if (doLoose) return (pu_id_mva_value > -0.30);
          else return (pu_id_mva_value > -0.17);
        }
        else return true;
      }
      else return true;
    }
    else return true;
  }

  bool IsReBTagged(Jet const* jet, std::map<std::size_t, bool> const& tag_map) {
    std::map<std::size_t,bool>::const_iterator it = tag_map.find(jet->id());
    if (it != tag_map.end()) {
      return it->second;
    } else {
      std::cout << "Warning in IsReBTagged: Jet with id " << jet->id() << " not found in result map" << std::endl;
      return true;
    }
  }


  double PZeta(CompositeCandidate const* cand, Candidate const* met, double const& alpha) {
    if (cand->AsVector().size() < 2) return 0.0;
    Candidate const* leg1 = cand->At(0);
    Candidate const* leg2 = cand->At(1);
    double leg1x = cos(leg1->phi());
    double leg1y = sin(leg1->phi());
    double leg2x = cos(leg2->phi());
    double leg2y = sin(leg2->phi());
    double zetaX = leg1x + leg2x;
    double zetaY = leg1y + leg2y;
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) {
     zetaX /= zetaR;
     zetaY /= zetaR;
   }

   double visPx = leg1->vector().px() + leg2->vector().px();
   double visPy = leg1->vector().py() + leg2->vector().py();
   double pZetaVis = visPx*zetaX + visPy*zetaY;

   double px = /*visPx + */met->vector().px();
   double py = /*visPy + */met->vector().py();
   double pZetaMiss = px*zetaX + py*zetaY;

   return (pZetaMiss - alpha * pZetaVis);
  }

  double PZetaVis(CompositeCandidate const* cand) {
     if (cand->AsVector().size() < 2) return 0.0;
     Candidate const* leg1 = cand->At(0);
     Candidate const* leg2 = cand->At(1); //MUG
     double leg1x = cos(leg1->phi());
     double leg1y = sin(leg1->phi());
     double leg2x = cos(leg2->phi());
     double leg2y = sin(leg2->phi());
     double zetaX = leg1x + leg2x;
     double zetaY = leg1y + leg2y;
     double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
     if ( zetaR > 0. ) {
      zetaX /= zetaR;
      zetaY /= zetaR;
    }

    double visPx = leg1->vector().px() + leg2->vector().px();
    double visPy = leg1->vector().py() + leg2->vector().py();
    double pZetaVis = visPx*zetaX + visPy*zetaY;
    return pZetaVis;
  }



  double MT(Candidate const* cand1, Candidate const* cand2) {
    double mt = 2. * cand1->pt() * cand2->pt() * (1. - cos(ROOT::Math::VectorUtil::DeltaPhi(cand1->vector(), cand2->vector())));
    if (mt > 0) {
      return std::sqrt(mt);
    } else {
      std::cerr << "Transverse mass would be negative! Returning 0.0" << std::endl;
    }
    return 0.0;
  }

  bool passAntiEMVA(Tau const * tau, int WP){
    
    double z_2 = tau->vz() + (130. / tan(tau->vector().theta()));
    if (z_2 > -1.5 && z_2 < 0.5) return false;
   
    int iCat=(int)(tau->GetTauID("againstElectronMVA3category")+0.5);
    float raw=tau->GetTauID("againstElectronMVA3raw");

    //Code from Ivo
    if(iCat<0) return false;
    if(iCat>15) return true;

    static float cutsMedium[16]={0.933,0.921,0.944,0.945,0.918,0.941,0.981,0.943,0.956,0.947,0.951,0.95,0.897,0.958,0.955,0.942};
    static float cutsTight[16]={ 0.96,0.968,0.971,0.972,0.969,0.959,0.981,0.965,0.975,0.972,0.974,0.971,0.897,0.971,0.961,0.97};
    static float cutsVeryTight[16]={0.978,0.98,0.982,0.985,0.977,0.974,0.989,0.977,0.986,0.983,0.984,0.983,0.971,0.987,0.977,0.981};
    float cut=0;
    if(WP==0) cut = cutsMedium[iCat];
    if(WP==1) cut = cutsTight[iCat];
    if(WP==2) cut = cutsVeryTight[iCat];
    return (raw>cut);
  }

  bool TauEoverP(Tau const* tau, double cut) {
    if (tau->decay_mode() != 0) {
      return true;
    } else {
      double e_over_p = (tau->lead_ecal_energy() + tau->lead_hcal_energy()) / tau->lead_p();
      if (e_over_p != e_over_p) return true;
      return e_over_p > cut;
    }
  }

  bool ElectronZbbID(Electron const * elec) {
    if (InEcalGap(elec)) return false;
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 1) return false;
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;
    double OneOverEOneOverP = fabs( (1./elec->ecal_energy()) - (1./(elec->ecal_energy()/elec->sc_e_over_p())) );
    return ( (in_barrel 
              && elec->sigma_IetaIeta()           < 0.01
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.06
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.004
              && OneOverEOneOverP                 < 0.05
              && elec->hadronic_over_em()         < 0.12)
          || (!in_barrel 
              && elec->sigma_IetaIeta()           < 0.03
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.03
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.007
              && OneOverEOneOverP                 < 0.05) );
  }
  
  bool MVAMETElectron(Electron const * elec) {
    if (fabs(elec->eta()) > 1.4442 && fabs(elec->eta()) < 1.566) return false;
    bool pass_kin = elec->pt() > 9.5 && fabs(elec->eta()) < 2.5;
    bool in_barrel = true;
    if (fabs(elec->eta()) > 1.4442) in_barrel = false;
    double OneOverEOneOverP = fabs( (1./elec->ecal_energy()) - (1./(elec->ecal_energy()/elec->sc_e_over_p())) );
    double pass_iso = (elec->dr03_tk_sum_pt()/elec->pt()) < 0.2;
    return pass_kin && pass_iso && 
          ( (in_barrel 
              && elec->gsf_tk_nhits()            == 0
              && elec->sigma_IetaIeta()           < 0.01
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.8
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.007
              && OneOverEOneOverP                 < 0.05
              && elec->hadronic_over_em()         < 0.15)
          || (!in_barrel 
              && elec->gsf_tk_nhits()            == 0
              && elec->sigma_IetaIeta()           < 0.03
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.10
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.009
              && OneOverEOneOverP                 < 0.05
              && elec->hadronic_over_em()         < 0.10) );
  }

  bool ElectronZbbIso(Electron const* elec, bool is_data, double const& rho, double const& cut) {
    double eta = fabs(elec->sc_eta());
    double photon_area, neutral_area;
    if (!is_data) {
      if (eta < 1.0)                  { neutral_area = 0.022; photon_area = 0.084; }    
      if (eta >= 1.0 && eta < 1.479)  { neutral_area = 0.039; photon_area = 0.090; }    
      if (eta >= 1.479 && eta < 2.0)  { neutral_area = 0.040; photon_area = 0.049; }    
      if (eta >= 2.0 && eta < 2.2)    { neutral_area = 0.028; photon_area = 0.099; }    
      if (eta >= 2.2 && eta < 2.3)    { neutral_area = 0.027; photon_area = 0.122; }    
      if (eta >= 2.3 && eta < 2.4)    { neutral_area = 0.024; photon_area = 0.132; }
      if (eta >= 2.4)                 { neutral_area = 0.030; photon_area = 0.155; }    
    } else {
      if (eta < 1.0)                  { neutral_area = 0.024; photon_area = 0.081; }    
      if (eta >= 1.0 && eta < 1.479)  { neutral_area = 0.037; photon_area = 0.084; }    
      if (eta >= 1.479 && eta < 2.0)  { neutral_area = 0.037; photon_area = 0.048; }    
      if (eta >= 2.0 && eta < 2.2)    { neutral_area = 0.023; photon_area = 0.089; }    
      if (eta >= 2.2 && eta < 2.3)    { neutral_area = 0.023; photon_area = 0.092; }    
      if (eta >= 2.3 && eta < 2.4)    { neutral_area = 0.021; photon_area = 0.097; }
      if (eta >= 2.4)                 { neutral_area = 0.021; photon_area = 0.110; }    
    }

    double iso = elec->dr03_pfiso_charged()
      + std::max(elec->dr03_pfiso_gamma() - rho * photon_area, 0.)
      + std::max(elec->dr03_pfiso_neutral() - rho * neutral_area, 0.);
      iso = iso / elec->pt();

    return (iso < cut);
  }



  bool Electron2011WP85Iso(Electron const* elec, double const& rho) {
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;
    // Lazy way
    //double hcal_sum = elec->dr03_hcal_tower_sum_et() + elec->hadronic_over_em() * elec->sc_energy() * fabs(sin(elec->sc_theta()));
    // Correct way
    double hcal_sum = elec->hcal_sum();
    double iso = 0;
    if (in_barrel) {
      iso = elec->dr03_tk_sum_pt() 
            + std::max(0., elec->dr03_ecal_rechit_sum_et() - 1.)
            + hcal_sum
            - rho * TMath::Pi() * 0.3 * 0.3;
    } else {
      iso = elec->dr03_tk_sum_pt()
            + elec->dr03_ecal_rechit_sum_et()
            + hcal_sum
            - rho * TMath::Pi() * 0.3 * 0.3;
    }

    iso = iso / elec->pt();
    return ( (in_barrel && iso < 0.053) || (!in_barrel && iso < 0.042) );
  }

  bool Electron2011WP85IsoNoRhoSub(Electron const* elec) {
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;
    // Lazy way
    //double hcal_sum = elec->dr03_hcal_tower_sum_et() + elec->hadronic_over_em() * elec->sc_energy() * fabs(sin(elec->sc_theta()));
    // Correct way
    double hcal_sum = elec->hcal_sum();
    double iso = 0;
    if (in_barrel) {
      iso = elec->dr03_tk_sum_pt() 
            + std::max(0., elec->dr03_ecal_rechit_sum_et() - 1.)
            + hcal_sum;
    } else {
      iso = elec->dr03_tk_sum_pt()
            + elec->dr03_ecal_rechit_sum_et()
            + hcal_sum;
    }

    iso = iso / elec->pt();
    return ( (in_barrel && iso < 0.053) || (!in_barrel && iso < 0.042) );
  }

  bool ElectronSimpleWP85Iso(Electron const* elec) {
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;
    double trk_iso = elec->dr03_tk_sum_pt() / elec->pt();
    double ecal_iso = elec->dr03_ecal_rechit_sum_et() / elec->pt();
    double hcal_iso = elec->dr03_hcal_tower_sum_et() / elec->pt();
    return ( (in_barrel 
              && trk_iso  < 0.09 
              && ecal_iso < 0.08 
              && hcal_iso < 0.10) 
          || (!in_barrel 
              && trk_iso  < 0.05 
              && ecal_iso < 0.05 
              && hcal_iso < 0.025) );
  }

  bool Electron2011WP85ID(Electron const* elec) {
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;

    return ( (in_barrel 
              && elec->sigma_IetaIeta()           < 0.01
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.039
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.005
              && elec->gsf_tk_nhits()            == 0
              && (fabs(elec->conv_dist()) >= 0.02 || fabs(elec->conv_dcot()) >= 0.02))
          || (!in_barrel 
              && elec->sigma_IetaIeta()           < 0.031
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.028
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.007
              && elec->gsf_tk_nhits()            == 0
              && (fabs(elec->conv_dist()) >= 0.02 || fabs(elec->conv_dcot()) >= 0.02)));
  }

  bool Electron2011WP70ID(Electron const* elec) {
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;

    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));

    return ( elec->gsf_tk_nhits()            <= 0
       && !elec->has_matched_conversion()
       && ( (in_barrel 
       && elec->sigma_IetaIeta()           < 0.01
       && fabs(elec->dphi_sc_tk_at_vtx())  < 0.03
       && fabs(elec->deta_sc_tk_at_vtx())  < 0.004
       && elec->hadronic_over_em() < 0.120
       && ooemoop < 0.050
       ) || 
      (!in_barrel 
       && elec->sigma_IetaIeta()           < 0.03
       && fabs(elec->dphi_sc_tk_at_vtx())  < 0.02
       && fabs(elec->deta_sc_tk_at_vtx())  < 0.005
       && elec->hadronic_over_em() < 0.100
       && ooemoop < 0.050
       )
      )
       );
  }

  bool TightElectronIDSpring15(Electron const* elec) {//function for spring15 ID 25ns !! needs updating
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.479) in_barrel = false;
    
    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));
//    double dbiso = elec->dr03_pfiso_charged() + std::max(0., elec->dr03_pfiso_neutral()+elec->dr03_pfiso_gamma() - 0.5*elec->dr03_pfiso_pu());
    
    return(
       !elec->has_matched_conversion()
        && ( (in_barrel       
	&& elec->full5x5_sigma_IetaIeta()   <0.0101
        && fabs(elec->deta_sc_tk_at_vtx())  <0.00926
        && fabs(elec->dphi_sc_tk_at_vtx())  <0.0336
        && elec->hadronic_over_em()         <0.0597
        && ooemoop                          <0.012
        && elec->gsf_tk_nhits()             <=2
	      //   && fabs(elec->dxy_vertex())         <0.060279
	      //  && fabs(elec->dz_vertex())          <0.800538
	      //  && dbiso/elec->pt()                 <0.164369
        ) ||
        (!in_barrel       
        && elec->full5x5_sigma_IetaIeta()   <0.0279
        && fabs(elec->deta_sc_tk_at_vtx())  <0.00724
        && fabs(elec->dphi_sc_tk_at_vtx())  <0.0918
        && elec->hadronic_over_em()         <0.0615
        && ooemoop                          <0.00999
        && elec->gsf_tk_nhits()             <=1
    //    && fabs(elec->dxy_vertex())         <0.273097
     //   && fabs(elec->dz_vertex())          <0.885860
   //     && dbiso/elec->pt()                 <0.212604
        ) )
        );

  }

  bool TightElectronFullIDSpring15(Electron const* elec, double const& rho) {//function for spring15 ID 25ns !! needs updating
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.479) in_barrel = false;
    
    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));
    //double dbiso = elec->dr03_pfiso_charged() + std::max(0., elec->dr03_pfiso_neutral()+elec->dr03_pfiso_gamma() - 0.5*elec->dr03_pfiso_pu());
    double lEA = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , elec->sc_eta() , ElectronEffectiveArea::kEleEAData2012);

    double relisoWithEA = (elec->dr03_pfiso_charged() + std::max(0., elec->dr03_pfiso_neutral()+elec->dr03_pfiso_gamma() - rho*lEA))/elec->pt();

    return(
	   !elec->has_matched_conversion()
	   && ( (in_barrel       
		 && elec->full5x5_sigma_IetaIeta()   <0.0101
		 && fabs(elec->deta_sc_tk_at_vtx())  <0.00926
		 && fabs(elec->dphi_sc_tk_at_vtx())  <0.0336
		 && elec->hadronic_over_em()         <0.0597
		 && ooemoop                          <0.012
		 && elec->gsf_tk_nhits()             <=2
		 && fabs(elec->dxy_vertex())         <0.0111
		 && fabs(elec->dz_vertex())          <0.0466
		 && relisoWithEA                     <0.0354
		 ) ||
		(!in_barrel       
		 && elec->full5x5_sigma_IetaIeta()   <0.0279
		 && fabs(elec->deta_sc_tk_at_vtx())  <0.00724
		 && fabs(elec->dphi_sc_tk_at_vtx())  <0.0918
		 && elec->hadronic_over_em()         <0.0615
		 && ooemoop                          <0.00999
		 && elec->gsf_tk_nhits()             <=1
		 && fabs(elec->dxy_vertex())         <0.0351
		 && fabs(elec->dz_vertex())          <0.417
		 && relisoWithEA                     <0.0646
		 ) )
	   );

  }

  bool VetoElectronFullIDSpring15(Electron const* elec, double const& rho) {//function for spring15 ID 25ns !! needs updating
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.479) in_barrel = false;
    
    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));
    //double dbiso = elec->dr03_pfiso_charged() + std::max(0., elec->dr03_pfiso_neutral()+elec->dr03_pfiso_gamma() - 0.5*elec->dr03_pfiso_pu());
    double lEA = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , elec->sc_eta() , ElectronEffectiveArea::kEleEAData2012);

    double relisoWithEA = (elec->dr03_pfiso_charged() + std::max(0., elec->dr03_pfiso_neutral()+elec->dr03_pfiso_gamma() - rho*lEA))/elec->pt();

    return(
	   !elec->has_matched_conversion()
	   && ( (in_barrel       
		 && elec->full5x5_sigma_IetaIeta()   <0.0114
		 && fabs(elec->deta_sc_tk_at_vtx())  <0.0152
		 && fabs(elec->dphi_sc_tk_at_vtx())  <0.216
		 && elec->hadronic_over_em()         <0.181
		 && ooemoop                          <0.207
		 && elec->gsf_tk_nhits()             <=2
		 && fabs(elec->dxy_vertex())         <0.0564
		 && fabs(elec->dz_vertex())          <0.472
		 && relisoWithEA                     <0.126
		 ) ||
		(!in_barrel       
		 && elec->full5x5_sigma_IetaIeta()   <0.0352
		 && fabs(elec->deta_sc_tk_at_vtx())  <0.0113
		 && fabs(elec->dphi_sc_tk_at_vtx())  <0.237
		 && elec->hadronic_over_em()         <0.116
		 && ooemoop                          <0.174
		 && elec->gsf_tk_nhits()             <=3
		 && fabs(elec->dxy_vertex())         <0.222
		 && fabs(elec->dz_vertex())          <0.921
		 && relisoWithEA                     <0.144
		 ) )
	   );

  }

  bool VetoElectronIDPhys14(Electron const* elec){
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.479) in_barrel = false;
    
    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));
    //double dbiso = elec->dr03_pfiso_charged() + std::max(0., elec->dr03_pfiso_neutral()+elec->dr03_pfiso_gamma() - 0.5*elec->dr03_pfiso_pu());
    
    return(
       !elec->has_matched_conversion()
        && ( (in_barrel       
        && elec->full5x5_sigma_IetaIeta()   <0.011586
        && fabs(elec->deta_sc_tk_at_vtx())  <0.013625
        && fabs(elec->dphi_sc_tk_at_vtx())  <0.230374
        && elec->hadronic_over_em()         <0.181130
        && ooemoop                          <0.295751
        && elec->gsf_tk_nhits()             <=2
     //   && fabs(elec->dxy_vertex())         <0.060279
      //  && fabs(elec->dz_vertex())          <0.800538
  //      && dbiso/elec->pt()                 <0.164369
        ) ||
        (!in_barrel       
        && elec->full5x5_sigma_IetaIeta()   <0.031849
        && fabs(elec->deta_sc_tk_at_vtx())  <0.011932
        && fabs(elec->dphi_sc_tk_at_vtx())  <0.255450
        && elec->hadronic_over_em()         <0.223870
        && ooemoop                          <0.155501
        && elec->gsf_tk_nhits()             <=3
    //    && fabs(elec->dxy_vertex())         <0.273097
     //   && fabs(elec->dz_vertex())          <0.885860
   //     && dbiso/elec->pt()                 <0.212604
        ) )
        );

   } 

  
  bool VetoElectronIDSpring15(Electron const* elec){//function for spring15 veto 25ns
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.479) in_barrel = false;
    
    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));
//    double dbiso = elec->dr03_pfiso_charged() + std::max(0., elec->dr03_pfiso_neutral()+elec->dr03_pfiso_gamma() - 0.5*elec->dr03_pfiso_pu());
    
    return(
       !elec->has_matched_conversion()
        && ( (in_barrel       
	&& elec->full5x5_sigma_IetaIeta()   <0.0114
        && fabs(elec->deta_sc_tk_at_vtx())  <0.0152
        && fabs(elec->dphi_sc_tk_at_vtx())  <0.216
        && elec->hadronic_over_em()         <0.181
        && ooemoop                          <0.207
        && elec->gsf_tk_nhits()             <=2
     //   && fabs(elec->dxy_vertex())         <0.060279
	      //  && fabs(elec->dz_vertex())          <0.800538
	      //  && dbiso/elec->pt()                 <0.164369
        ) ||
        (!in_barrel       
        && elec->full5x5_sigma_IetaIeta()   <0.0352
        && fabs(elec->deta_sc_tk_at_vtx())  <0.0113
        && fabs(elec->dphi_sc_tk_at_vtx())  <0.237
        && elec->hadronic_over_em()         <0.116
        && ooemoop                          <0.174
        && elec->gsf_tk_nhits()             <=3
    //    && fabs(elec->dxy_vertex())         <0.273097
     //   && fabs(elec->dz_vertex())          <0.885860
   //     && dbiso/elec->pt()                 <0.212604
        ) )
        );

   } 


  bool VetoElectronID(Electron const* elec) {
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;

    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));

    return ( elec->gsf_tk_nhits()            <= 999
       //&& !elec->has_matched_conversion()
       && ( (in_barrel 
       && elec->sigma_IetaIeta()           < 0.01
       && fabs(elec->dphi_sc_tk_at_vtx())  < 0.800
       && fabs(elec->deta_sc_tk_at_vtx())  < 0.007
       && elec->hadronic_over_em() < 0.150
       && ooemoop < 999.9
       ) || 
      (!in_barrel 
       && elec->sigma_IetaIeta()           < 0.03
       && fabs(elec->dphi_sc_tk_at_vtx())  < 0.700
       && fabs(elec->deta_sc_tk_at_vtx())  < 0.010
       && elec->hadronic_over_em() < 999.9
       && ooemoop < 999.9
       )
      )
       );
  }

  bool Electron2011WP95ID(Electron const* elec) {
    //if (elec->has_matched_conversion()) return false;
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.479) in_barrel = false;
//    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;

    return ( (in_barrel 
              && elec->sigma_IetaIeta()           < 0.01
              && elec->hadronic_over_em()         < 0.15
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.8
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.007
              /*&& elec->gsf_tk_nhits()            == 0*/)
          || (!in_barrel 
              && elec->sigma_IetaIeta()           < 0.03
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.7
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.01
              /*&& elec->gsf_tk_nhits()            == 0*/));
  }


  bool ElectronSimpleWP85ID(Electron const* elec) {

    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.4442) in_barrel = false;

    return ( (in_barrel 
              && elec->sigma_IetaIeta()           < 0.01
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.06
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.006
              && elec->hadronic_over_em()         < 0.04
              && elec->gsf_tk_nhits()            <= 1
              && (fabs(elec->conv_dist()) >= 0.02 || fabs(elec->conv_dcot()) >= 0.02))
          || (!in_barrel 
              && elec->sigma_IetaIeta()           < 0.03
              && fabs(elec->dphi_sc_tk_at_vtx())  < 0.04
              && fabs(elec->deta_sc_tk_at_vtx())  < 0.007
              && elec->hadronic_over_em()         < 0.025
              && elec->gsf_tk_nhits()            <= 1
              && (fabs(elec->conv_dist()) >= 0.02 || fabs(elec->conv_dcot()) >= 0.02)));
  }


  bool ElectronHTTIdIso(Electron const* elec, unsigned const& mode) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 0) return false;
    bool pass_mva = false;
    double pt = elec->pt();
    double eta = fabs(elec->sc_eta());
    double idmva = elec->GetIdIso("mvaNonTrigV0");
    if (pt <= 20.0 && eta <= 0.8                  && idmva > 0.925) pass_mva = true;
    if (pt <= 20.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.915) pass_mva = true;
    if (pt <= 20.0 && eta >  1.479                && idmva > 0.965) pass_mva = true;
    if (mode != 2) {
      if (pt >  20.0 && eta <= 0.8                  && idmva > 0.925) pass_mva = true;
      if (pt >  20.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.975) pass_mva = true;
      if (pt >  20.0 && eta >  1.479                && idmva > 0.985) pass_mva = true;
    } else {
      if (pt >  20.0 && eta <= 0.8                  && idmva > 0.905) pass_mva = true;
      if (pt >  20.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.955) pass_mva = true;
      if (pt >  20.0 && eta >  1.479                && idmva > 0.975) pass_mva = true;
    }
    return pass_mva;
  }

  bool ElectronHTTVHID(Electron const* elec) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 0) return false;
    bool pass_mva = false;
    double pt = elec->pt();
    double eta = fabs(elec->sc_eta());
    double idmva = elec->GetIdIso("mvaNonTrigV0");
    if (pt >= 10.0 && eta <= 0.8                  && idmva > 0.5) pass_mva = true;
    if (pt >= 10.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.12) pass_mva = true;
    if (pt >= 10.0 && eta >  1.479                && idmva > 0.6) pass_mva = true;
    return pass_mva;
  }

  bool ElectronHTTId(Electron const* elec, bool loose_wp) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 0) return false;
    bool pass_mva = false;
    double pt = elec->pt();
    double eta = fabs(elec->sc_eta());
    double idmva = elec->GetIdIso("mvaNonTrigV0");
    if (pt <= 20.0 && eta <= 0.8                  && idmva > 0.925) pass_mva = true;
    if (pt <= 20.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.915) pass_mva = true;
    if (pt <= 20.0 && eta >  1.479                && idmva > 0.965) pass_mva = true;
    if (!loose_wp) {
      if (pt >  20.0 && eta <= 0.8                  && idmva > 0.925) pass_mva = true;
      if (pt >  20.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.975) pass_mva = true;
      if (pt >  20.0 && eta >  1.479                && idmva > 0.985) pass_mva = true;
    } else {
      if (pt >  20.0 && eta <= 0.8                  && idmva > 0.905) pass_mva = true;
      if (pt >  20.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.955) pass_mva = true;
      if (pt >  20.0 && eta >  1.479                && idmva > 0.975) pass_mva = true;
    }
    return pass_mva;
  }
  
  bool ElectronHTTIdPhys14(Electron const* elec, bool loose_wp) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 0) return false;
    bool pass_mva = false;
    double eta = fabs(elec->sc_eta());
    double idmva = elec->GetIdIso("mvaNonTrigV025nsPHYS14");
    if (!loose_wp) {
      if (eta <= 0.8                  && idmva > 0.73) pass_mva = true;
      if (eta >  0.8 && eta <= 1.479  && idmva > 0.57) pass_mva = true;
      if (eta >  1.479                && idmva > 0.05) pass_mva = true;
    } else {
      if (eta <= 0.8                  && idmva > 0.35) pass_mva = true;
      if (eta >  0.8 && eta <= 1.479  && idmva > 0.20) pass_mva = true;
      if (eta >  1.479                && idmva > -0.52) pass_mva = true;
    }
    return pass_mva;
  }

  bool ElectronHTTIdSpring15(Electron const* elec, bool loose_wp) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 1) return false;
    bool pass_mva = false;
    double eta = fabs(elec->sc_eta());
    double pt = fabs(elec->pt());
    double idmva = elec->GetIdIso("mvaNonTrigSpring15");
    if (!loose_wp) {
      if (eta <= 0.8 && pt <= 10                  && idmva > 0.287435) pass_mva = true;
      if (eta >  0.8 && eta <= 1.479 && pt <=10   && idmva > 0.221846) pass_mva = true;
      if (eta >  1.479 && pt <= 10                && idmva > -0.303263) pass_mva = true;
      if (eta <= 0.8 && pt > 10                   && idmva > 0.967083) pass_mva = true;
      if (eta >  0.8 && eta <= 1.479 && pt > 10   && idmva > 0.929117) pass_mva = true;
      if (eta >  1.479 && pt > 10                 && idmva > 0.726311) pass_mva = true;

    } else {
      if (eta <= 0.8 && pt <= 10                  && idmva > -0.083313) pass_mva = true;
      if (eta >  0.8 && eta <= 1.479 && pt <=10   && idmva > -0.235222) pass_mva = true;
      if (eta >  1.479 && pt <= 10                && idmva > -0.67099) pass_mva = true;
      if (eta <= 0.8 && pt > 10                   && idmva > 0.913286) pass_mva = true;
      if (eta >  0.8 && eta <= 1.479 && pt > 10   && idmva > 0.805013) pass_mva = true;
      if (eta >  1.479 && pt > 10                 && idmva > 0.358969) pass_mva = true;
    }
    return pass_mva;
  }
  
  bool ElectronHTTIdSpring16(Electron const* elec, bool loose_wp) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;      
    if (elec->gsf_tk_nhits() > 1) return false;              
    bool pass_mva = false;
    double eta = fabs(elec->sc_eta());
    double pt = fabs(elec->pt());
    double idmva = elec->GetIdIso("generalPurposeMVASpring16");
    if (!loose_wp) {
      if (eta <= 0.8 && pt > 10                   && idmva > 0.940962684155 ) pass_mva = true;
      if (eta > 0.8 && eta <= 1.479 && pt > 10    && idmva > 0.899208843708 ) pass_mva = true;
      if (eta > 1.479 && pt > 10                  && idmva > 0.758484721184 ) pass_mva = true;
    } else {
      if (eta <= 0.8 && pt > 10                   && idmva > 0.836695742607 ) pass_mva = true;
      if (eta > 0.8 && eta <= 1.479 && pt > 10    && idmva > 0.715337944031 ) pass_mva = true;
      if (eta > 1.479 && pt > 10                  && idmva > 0.356799721718 ) pass_mva = true;
    }
    return pass_mva;
  }
  
   bool VetoElectronIDSpring16(Electron const* elec) {
    
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.479) in_barrel = false;

    double dEtaInSeed = elec->deta_sc_tk_at_vtx() - elec->sc_eta() + elec->sc_seed_eta(); 
    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));

    return (!elec->has_matched_conversion()
       && ( (in_barrel 
       && elec->sigma_IetaIeta()           < 0.0115
       && fabs(dEtaInSeed)                 < 0.00749
       && fabs(elec->dphi_sc_tk_at_vtx())  < 0.228
       && elec->hadronic_over_em()         < 0.356
       && ooemoop                          < 0.299
       && elec->gsf_tk_nhits()             <=2
       ) || 
      (!in_barrel 
       && elec->sigma_IetaIeta()           < 0.037
       && fabs(dEtaInSeed)                 < 0.00895
       && fabs(elec->dphi_sc_tk_at_vtx())  < 0.213
       && elec->hadronic_over_em()         < 0.211
       && ooemoop                          < 0.15
       && elec->gsf_tk_nhits()             <=3
       )
      )
    );
  }

  bool VetoElectronIDFall17(Electron const* elec, double const& rho) {
    bool in_barrel = true;
    if (fabs(elec->sc_eta()) > 1.479) in_barrel = false;
    double dEtaInSeed = elec->deta_sc_tk_at_vtx() - elec->sc_eta() + elec->sc_seed_eta();
    double ooemoop = fabs((1.0/elec->ecal_energy() - elec->sc_e_over_p()/elec->ecal_energy()));
    double sc_energy = elec->sc_energy();

    return (!elec->has_matched_conversion()
       && ( (in_barrel
       && elec->sigma_IetaIeta()           < 0.0128
       && fabs(dEtaInSeed)                 < 0.00523
       && fabs(elec->dphi_sc_tk_at_vtx())  < 0.159
       && elec->hadronic_over_em()         < 0.05 + 1.12/sc_energy + 0.0368*rho/sc_energy 
       && ooemoop                          < 0.193
       && elec->gsf_tk_nhits()             <=2
       ) ||
      (!in_barrel
       && elec->sigma_IetaIeta()           < 0.0445
       && fabs(dEtaInSeed)                 < 0.00984
       && fabs(elec->dphi_sc_tk_at_vtx())  < 0.157
       && elec->hadronic_over_em()         < 0.05 + 0.5/sc_energy + 0.201*rho/sc_energy
       && ooemoop                          < 0.0962
       && elec->gsf_tk_nhits()             <=3
       )
      )
    );
  }

 bool ElectronHTTIdFall17(Electron const* elec, bool loose_wp) {
   //Do some cut-based pre-selection
   if (elec->has_matched_conversion()) return false;
     if (elec->gsf_tk_nhits() > 1) return false;
     double eta = fabs(elec->sc_eta());
     double pt = fabs(elec->pt());
     double idmva = elec->GetIdIso("mvaRun2Fall17");
     double c=0;
     double tau=0;
     double A=0;
     if (eta <= 0.8) {
        if (pt > 10.){
          c = loose_wp ? 0.9616542816132922 : 0.9825268564943458;
          tau = loose_wp ? 8.757943837889817 : 8.702601455860762;
          A = loose_wp ? 3.1390200321591206 : 1.1974861596609097;
        } else {
          c = loose_wp ? 0.9165112826974601 : 0.9530240956555949;
          tau = loose_wp ? 2.7381703555094217 : 2.7591425841003647;
          A = loose_wp ? 1.03549199648109 : 0.4669644718545271; 
        }
      } 
      if (eta > 0.8 && eta <= 1.479) {
        if (pt > 10.){
          c = loose_wp ? 0.9319258011430132 : 0.9727509457929913;
          tau = loose_wp ? 8.846057432565809 : 8.179525631018565;
          A = loose_wp ? 3.5985063793347787 : 1.7111755094657688;
        } else {
          c = loose_wp ? 0.8655738322220173 : 0.9336564763961019;
          tau = loose_wp ? 2.4027944652597073 : 2.709276284272272;
          A = loose_wp ? 0.7975615613282494 : 0.33512286599215946;
        }
      }
      if (eta > 1.479){
        if (pt > 10.){
          c = loose_wp ? 0.8899260780999244 : 0.9562619539540145;
          tau = loose_wp ? 10.124234115859881 : 8.109845366281608;
          A = loose_wp ? 4.352791250718547 : 3.013927699126942;
        } else {
          c = loose_wp ? -3016.035055227131 : 0.9313133688365339;
          tau = loose_wp ? -52140.61856333602 : 1.5821934800715558;
          A = loose_wp ? -3016.3029387236506 : 3.8889462619659265;
        }
      } 
    double wp = c - A*exp(-pt/tau);
    bool pass_mva = idmva > wp;
    return pass_mva;
  }
  
 bool ElectronHTTIsoIdFall17(Electron const* elec, bool loose_wp) {
   //Do some cut-based pre-selection
   if (elec->has_matched_conversion()) return false;
     if (elec->gsf_tk_nhits() > 1) return false;
     double eta = fabs(elec->sc_eta());
     double pt = fabs(elec->pt());
     double idmva = elec->GetIdIso("mvaRun2IsoFall17");
     double c=0;
     double tau=0;
     double A=0;
     if (eta <= 0.8) {
          // c = C0, tau = C1,  A = C2 
        if (pt > 10.){
          c = loose_wp ?  0.9717674837607253 : 0.9896562087723659;
          tau = loose_wp ? 8.912850985100356 : 10.342490511998674;
          A = loose_wp ? 1.9712414940437244 : 0.40204156417414094;
        } else {
          c = loose_wp ? 0.9387070396095831 :  0.9725509559754997;
          tau = loose_wp ? 2.6525585228167636 : 2.976593261509491;
          A = loose_wp ? 0.8222647164151365 : 0.2653858736397496; 
        }
      } 
      if (eta > 0.8 && eta <= 1.479) {
        if (pt > 10.){
          c = loose_wp ? 0.9458745023265976 : 0.9819232656533827;
          tau = loose_wp ? 8.83104420392795 :  9.05548836482051;
          A = loose_wp ? 2.40849932040698 : 0.772674931169389;
        } else {
          c = loose_wp ? 0.8948802925677235 : 0.9508038141601247;
          tau = loose_wp ? 2.7645670358783523 : 2.6633500558725713;
          A = loose_wp ? 0.4123381218697539 : 0.2355820499260076;
        }
      }
      if (eta > 1.479){
        if (pt > 10.){
          c = loose_wp ? 0.8979112012086751 : 0.9625098201744635;
          tau = loose_wp ? 9.814082144168015 : 8.42589315557279;
          A = loose_wp ? 4.171581694893849 : 2.2916152615134173;
        } else {
          c = loose_wp ?  -1830.8583661119892 : 0.9365037167596238;
          tau = loose_wp ?  -36578.11055382301 : 1.5765442323949856;
          A = loose_wp ?  -1831.2083578116517 : 3.067015289215309;
        }
      } 
    double wp = c - A*exp(-pt/tau);
    bool pass_mva = idmva > wp;
    return pass_mva;
  }


  // Fall v2 electron ID noIso
  bool ElectronHTTIdFall17V2(Electron const* elec, bool loose_wp) {
    // loose_wp will correspond to 90% wp
    // there is a looser one (98% wp) but don't need it
    // 
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
      if (elec->gsf_tk_nhits() > 1) return false;
      double eta = fabs(elec->sc_eta());
      double pt = fabs(elec->pt());
      double idmva = elec->GetIdIso("mvaRun2Fall17V2");
      double c=0;
      double tau=0;
      double A=0;
      if (eta <= 0.8) {
        if (pt > 10.){ //EB1_10
          c = loose_wp ? 5.9175992258 : 7.1336238874;
          tau = loose_wp ? 13.4807294538 : 16.5605268797;
          A = loose_wp ? 9.31966232685 : 8.22531222391;
        } else { //EB1_5
          c = loose_wp ? 2.77072387339 : 3.26449620468;
          tau = loose_wp ? 3.81500912145 : 3.32657149223;
          A = loose_wp ? 8.16304860178 : 8.84669783568; 
        }
      } 
      if (eta > 0.8 && eta <= 1.479) {
        if (pt > 10.){ //EB2_10
          c = loose_wp ? 5.01598837255 : 6.18638275782;
          tau = loose_wp ? 13.1280451502 : 15.2694634284;
          A = loose_wp ? 8.79418193765 : 7.49764565324;
        } else { //EB2_5
          c = loose_wp ? 1.85602317813 : 2.83557838497;
          tau = loose_wp ? 2.18697654938 : 2.15150487651;
          A = loose_wp ? 11.8568936824 : 11.0978016567;
        }
      }
      if (eta > 1.479){
        if (pt > 10.){ //EE_10
          c = loose_wp ? 4.16921343208 : 5.43175865738;
          tau = loose_wp ? 13.2017224621 : 15.4290075949;
          A = loose_wp ? 9.00720913211 : 7.56899692285;
        } else { //EE_5
          c = loose_wp ? 1.73489307814 : 2.91994945177;
          tau = loose_wp ? 2.0163211971 : 1.69875477522;
          A = loose_wp ? 17.013880078 : 24.024807824;
        }
      } 
      double wp = c - A*exp(-pt/tau);
      bool pass_mva = idmva > wp;
      return pass_mva;
    }

  // Fall v2 electron IDIso
  bool ElectronHTTIsoIdFall17V2(Electron const* elec, bool loose_wp) {
    // loose_wp will correspond to 90% wp
    // there is a looser one (98% wp) but don't need it
    // 
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
      if (elec->gsf_tk_nhits() > 1) return false;
      double eta = fabs(elec->sc_eta());
      double pt = fabs(elec->pt());
      double idmva = elec->GetIdIso("mvaRun2IsoFall17V2");
      double c=0;
      double tau=0;
      double A=0;
      if (eta <= 0.8) {
        if (pt > 10.){ //EB1_10
          c = loose_wp ? 6.12931925263 : 7.35752275071;
          tau = loose_wp ? 13.281753835 : 15.87907864;
          A = loose_wp ? 8.71138432196 : 7.61288809226;
        } else { //EB1_5
          c = loose_wp ? 2.84704783417 : 3.53495358797;
          tau = loose_wp ? 3.32529515837 : 3.07272325141;
          A = loose_wp ? 9.38050947827 : 9.94262764352; 
        }
      } 
      if (eta > 0.8 && eta <= 1.479) {
        if (pt > 10.){ //EB2_10
          c = loose_wp ? 5.26289004857 : 6.41811074032;
          tau = loose_wp ? 13.2154971491 : 14.730562874;
          A = loose_wp ? 8.0997882835 : 6.96387331587;
        } else { //EB2_5
          c = loose_wp ? 2.03833922005 : 3.06015605623;
          tau = loose_wp ? 1.93288758682 : 1.95572234114;
          A = loose_wp ? 15.364588247 : 14.3091184421;
        }
      }
      if (eta > 1.479){
        if (pt > 10.){ //EE_10
          c = loose_wp ? 4.37338792902 : 5.64936312428;
          tau = loose_wp ? 14.0776094696 : 16.3664949747;
          A = loose_wp ? 8.48513324496 : 7.19607610311;
        } else { //EE_5
          c = loose_wp ? 1.82704158461 : 3.02052519639;
          tau = loose_wp ? 1.89796754399 : 1.59784164742;
          A = loose_wp ? 19.1236071158 : 28.719380105;
        }
      } 
      double wp = c - A*exp(-pt/tau);
      bool pass_mva = idmva > wp;
      return pass_mva;
   }
  
 bool PF03EAElecIsolation(Electron const* elec, double const rho, double const& cut) {
    double charged_iso = elec->dr03_pfiso_charged();
    double eff_area = GetEffectiveArea2017(elec);
    double iso =  charged_iso 
                  + std::max(elec->dr03_pfiso_neutral() + elec->dr03_pfiso_gamma() - rho*eff_area, 0.0);
    iso = iso / elec->pt();
    return (iso < cut);
  } 


  
 double PUW03IsolationVal(Muon const* muon){
   double charged_iso = muon->dr03_pfiso_charged();
   double neutral_weighted = muon->GetIdIso("neutral_pfweighted_iso_03");
   double gamma_weighted = muon->GetIdIso("gamma_pfweighted_iso_03");
   double iso = charged_iso + neutral_weighted + gamma_weighted;
    iso = iso / muon->pt();
    return iso;
 }
 
 double PUW04IsolationVal(Muon const* muon){
   double charged_iso = muon->dr04_pfiso_charged();
   double neutral_weighted = muon->GetIdIso("neutral_pfweighted_iso_04");
   double gamma_weighted = muon->GetIdIso("gamma_pfweighted_iso_04");
   double iso = charged_iso + neutral_weighted + gamma_weighted;
    iso = iso / muon->pt();
    return iso;
 }


  bool ElectronHTTIdTrigSpring15(Electron const* elec, bool loose_wp) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 1) return false;
    bool pass_mva = false;
    double eta = fabs(elec->sc_eta());
    double pt = fabs(elec->pt());
    bool pass_preselection = false;
/*    if(eta <= 1.479 && pt > 15 && elec->full5x5_sigma_IetaIeta()<0.012&&elec->hadronic_over_em()<0.09 && elec->dr03_tk_sum_pt()/pt <0.18&&elec->deta_sc_tk_at_vtx()<0.0095&&elec->dphi_sc_tk_at_vtx()<0.065) pass_preselection = true;
    if(eta > 1.479 && pt > 15 && elec->full5x5_sigma_IetaIeta()<0.033&&elec->hadronic_over_em()<0.09  && elec->dr03_tk_sum_pt()/pt <0.18) pass_preselection = true;
*/
     if(eta <= 1.479 && pt > 15 && elec->full5x5_sigma_IetaIeta()<0.012&&elec->hadronic_over_em()<0.09 && elec->ecal_pf_cluster_iso()/pt <0.37 && elec->hcal_pf_cluster_iso()/pt < 0.25 && elec->dr03_tk_sum_pt()/pt <0.18&&abs(elec->deta_sc_tk_at_vtx()<0.0095)&&abs(elec->dphi_sc_tk_at_vtx()<0.065)) pass_preselection = true;
    if(eta > 1.479 && pt > 15 && elec->full5x5_sigma_IetaIeta()<0.033&&elec->hadronic_over_em()<0.09 && elec->ecal_pf_cluster_iso()/pt <0.45 && elec->hcal_pf_cluster_iso()/pt < 0.28 && elec->dr03_tk_sum_pt()/pt <0.18) pass_preselection = true;
    
    if(!pass_preselection) return false;
    double idmva = elec->GetIdIso("mvaTrigSpring15");
    if (!loose_wp) {
      if (eta <= 0.8                    && idmva > 0.988153) pass_mva = true;
      if (eta >  0.8 && eta <= 1.479   && idmva > 0.967910) pass_mva = true;
      if (eta >  1.479                  && idmva > 0.841729) pass_mva = true;
    } else {
      if (eta <= 0.8                    && idmva > 0.972153) pass_mva = true;
      if (eta >  0.8 && eta <= 1.479    && idmva > 0.922126) pass_mva = true;
      if (eta >  1.479                 && idmva > 0.610764) pass_mva = true;
    }
    return pass_mva;
  }



  bool ElectronHTTTrigNoIPId(Electron const* elec, bool loose_wp) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 0) return false;
    bool pass_mva = false;
    double pt = elec->pt();
    double eta = fabs(elec->sc_eta());
    double idmva = elec->GetIdIso("mvaTrigNoIPV0");
    if (pt <= 20.0 && eta <= 0.8                  && idmva > -0.3500) pass_mva = true;
    if (pt <= 20.0 && eta >  0.8 && eta <= 1.479  && idmva >  0.0000) pass_mva = true;
    if (pt <= 20.0 && eta >  1.479                && idmva >  0.0250) pass_mva = true;
    if (!loose_wp) {
      if (pt >  20.0 && eta <= 0.8                  && idmva > 0.7875) pass_mva = true;
      if (pt >  20.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.9500) pass_mva = true;
      if (pt >  20.0 && eta >  1.479                && idmva > 0.9250) pass_mva = true;
    } else {
      if (pt >  20.0 && eta <= 0.8                  && idmva > 0.7000) pass_mva = true;
      if (pt >  20.0 && eta >  0.8 && eta <= 1.479  && idmva > 0.9000) pass_mva = true;
      if (pt >  20.0 && eta >  1.479                && idmva > 0.8375) pass_mva = true;
    }
    return pass_mva;
  }

  bool HttEMuFakeElectron(Electron const* elec) {
    if (elec->has_matched_conversion()) return false; // ConversionTools::hasMatchedConversion(...)
    bool in_barrel = fabs(elec->sc_eta()) < 1.479;
    if (elec->gsf_tk_nhits() > 0) return false;
    bool pass_iso = false;
    if (in_barrel) {
      pass_iso =        ((elec->dr03_tk_sum_pt() / elec->pt())                                < 0.2) // dr03TkSumPt()
                    &&  ((TMath::Max(elec->dr03_ecal_rechit_sum_et()-1.0, 0.0) / elec->pt())  < 0.2) // dr03EcalRecHitSumEt()
                    &&  ((elec->dr03_hcal_tower_sum_et() / elec->pt())                        < 0.2); // dr03HcalTowerSumEt()
    } else {
      pass_iso =        ((elec->dr03_tk_sum_pt() / elec->pt())          < 0.2) 
                    &&  ((elec->dr03_ecal_rechit_sum_et() / elec->pt()) < 0.2)
                    &&  ((elec->dr03_hcal_tower_sum_et() / elec->pt())  < 0.2);
    }
   if (!pass_iso) return false;     
   bool pass_id = ( 
            fabs(elec->dxy_vertex()) < 0.02 &&
            fabs(elec->dz_vertex()) < 0.1 &&
            ((in_barrel 
             && elec->sigma_IetaIeta()           < 0.01   // sigmaIetaIeta()
             && fabs(elec->dphi_sc_tk_at_vtx())  < 0.15   // deltaPhiSuperClusterTrackAtVtx()
             && fabs(elec->deta_sc_tk_at_vtx())  < 0.007
             && elec->hadronic_over_em()         < 0.12) // deltaEtaSuperClusterTrackAtVtx()
         || (!in_barrel 
             && elec->sigma_IetaIeta()           < 0.03
             && fabs(elec->dphi_sc_tk_at_vtx())  < 0.1
             && fabs(elec->deta_sc_tk_at_vtx())  < 0.009
             && elec->hadronic_over_em()         < 0.10)));
   return pass_id && !(ElectronHTTId(elec, true) && PF04IsolationEB(elec, 0.5, 0.15, 0.1) && fabs(elec->dxy_vertex()) < 0.02 && fabs(elec->dz_vertex()) < 0.1 ); // !(usual MVA ID && pf isolation < 0.15(0.1) barrel(endcap))
  }

  

  bool PF04IsolationEBElec(Electron const* cand, double const& dbeta, double const& cut_barrel, double const& cut_endcap) {
    double iso =  PF04IsolationVal(cand, dbeta);
    if (fabs(cand->sc_eta()) < 1.479) {
      return (iso < cut_barrel);
    } else {
      return (iso < cut_endcap);
    }
  }


  bool ElectronHTT2011IdIso(Electron const* elec) {
    //Do some cut-based pre-selection
    if (elec->has_matched_conversion()) return false;
    if (elec->gsf_tk_nhits() > 0) return false;
    bool pass_preid = false;
    if(fabs(elec->sc_eta()) < 1.479) {
      if(elec->sigma_IetaIeta()                       < 0.01 
        && elec->hadronic_over_em()                   < 0.12 
        && fabs(elec->deta_sc_tk_at_vtx())            < 0.007 
        && fabs(elec->dphi_sc_tk_at_vtx())            < 0.15 
        && elec->gsf_tk_nhits() == 0)
        pass_preid = true;
    }
    else {
      if(elec->sigma_IetaIeta()                       < 0.03 
        && elec->hadronic_over_em()                   < 0.10 
        && fabs(elec->deta_sc_tk_at_vtx())            < 0.009 
        && fabs(elec->dphi_sc_tk_at_vtx())            < 0.10
        && elec->gsf_tk_nhits() == 0)
        pass_preid = true;
    }
    if (!pass_preid) return false;
    bool pass_mva = false;
    double pt = elec->pt();
    double eta = fabs(elec->sc_eta());
    double idmva = elec->GetIdIso("id2011mva");
    if (pt <= 20.0 && eta <= 1.0                  && idmva > 0.133) pass_mva = true;
    if (pt <= 20.0 && eta >  1.0 && eta <= 1.5    && idmva > 0.465) pass_mva = true;
    if (pt <= 20.0 && eta >  1.5                  && idmva > 0.518) pass_mva = true;
    if (pt >  20.0 && eta <= 1.0                  && idmva > 0.942) pass_mva = true;
    if (pt >  20.0 && eta >  1.0 && eta <= 1.5    && idmva > 0.947) pass_mva = true;
    if (pt >  20.0 && eta >  1.5                  && idmva > 0.878) pass_mva = true;
    return pass_mva;
  }

  //Photon IDs taken from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
  bool TightPhotonIDSpring15(Photon const* photon,double const& rho) {//function for spring15 ID
    double eta = fabs(photon->eta());

    bool in_barrel = true;
    if (eta > 1.479) in_barrel = false;
    bool in_endcap = true;
    if (eta < 1.479) in_endcap = false;

    double photon_area, neutral_area, charged_area;
    if (eta < 1.0)                  { charged_area=0.0158; neutral_area=0.0143; photon_area=0.0725; }
    if (eta >= 1.0 && eta < 1.479)  { charged_area=0.0143; neutral_area=0.0210; photon_area=0.0604; }
    if (eta >= 1.479 && eta < 2.0)  { charged_area=0.0115; neutral_area=0.0147; photon_area=0.0320; }
    if (eta >= 2.0 && eta < 2.2)    { charged_area=0.0094; neutral_area=0.0082; photon_area=0.0512; }
    if (eta >= 2.2 && eta < 2.3)    { charged_area=0.0095; neutral_area=0.0124; photon_area=0.0766; }
    if (eta >= 2.3 && eta < 2.4)    { charged_area=0.0068; neutral_area=0.0186; photon_area=0.0949; }
    if (eta >= 2.4)                 { charged_area=0.0053; neutral_area=0.0320; photon_area=0.1160; }

    return(
	( (in_barrel       
	   && photon->had_tower_over_em()<0.05
	   && photon->sigma_IetaIeta()<0.01
	   && std::max(photon->dr03_pfiso_charged()-rho*charged_area,0.)<0.91
	   && std::max(photon->dr03_pfiso_neutral()-rho*neutral_area,0.)<(0.33+exp(0.0044*photon->pt()+0.5809))
           && std::max(photon->dr03_pfiso_gamma()-rho*photon_area,0.)<(0.61+0.0043*photon->pt())
	   ) ||
	  (in_endcap       
	   && photon->had_tower_over_em()<0.05
	   && photon->sigma_IetaIeta()<0.0267
	   && std::max(photon->dr03_pfiso_charged()-rho*charged_area,0.)<0.65
	   && std::max(photon->dr03_pfiso_neutral()-rho*neutral_area,0.)<(0.93+exp(0.004*photon->pt()+0.9402))
	   && std::max(photon->dr03_pfiso_gamma()-rho*photon_area,0.)<(0.54+0.0041*photon->pt())
	   ) )
        );

  }

  bool MediumPhotonIDSpring15(Photon const* photon,double const& rho) {//function for spring15 ID
    double eta = fabs(photon->eta());
    bool in_barrel = true;
    if (eta > 1.479) in_barrel = false;
    bool in_endcap = true;
    if (eta < 1.479) in_endcap = false;

    double photon_area, neutral_area, charged_area;
    if (eta < 1.0)                  { charged_area=0.0158; neutral_area=0.0143; photon_area=0.0725; }
    if (eta >= 1.0 && eta < 1.479)  { charged_area=0.0143; neutral_area=0.0210; photon_area=0.0604; }
    if (eta >= 1.479 && eta < 2.0)  { charged_area=0.0115; neutral_area=0.0147; photon_area=0.0320; }
    if (eta >= 2.0 && eta < 2.2)    { charged_area=0.0094; neutral_area=0.0082; photon_area=0.0512; }
    if (eta >= 2.2 && eta < 2.3)    { charged_area=0.0095; neutral_area=0.0124; photon_area=0.0766; }
    if (eta >= 2.3 && eta < 2.4)    { charged_area=0.0068; neutral_area=0.0186; photon_area=0.0949; }
    if (eta >= 2.4)                 { charged_area=0.0053; neutral_area=0.0320; photon_area=0.1160; }
    
    return(
	( (in_barrel       
	   && photon->had_tower_over_em()<0.05
	   && photon->sigma_IetaIeta()<0.01
	   && std::max(photon->dr03_pfiso_charged()-rho*charged_area,0.)<1.31
	   && std::max(photon->dr03_pfiso_neutral()-rho*neutral_area,0.)<(0.6+exp(0.0044*photon->pt()+0.5809))
           && std::max(photon->dr03_pfiso_gamma()-rho*photon_area,0.)<(1.33+0.0043*photon->pt())
	   ) ||
	  (in_endcap       
	   && photon->had_tower_over_em()<0.05
	   && photon->sigma_IetaIeta()<0.0267
	   && std::max(photon->dr03_pfiso_charged()-rho*charged_area,0.)<1.25
	   && std::max(photon->dr03_pfiso_neutral()-rho*neutral_area,0.)<(1.65+exp(0.004*photon->pt()+0.9402))
	   && std::max(photon->dr03_pfiso_gamma()-rho*photon_area,0.)<(1.02+0.0041*photon->pt())
	   ) )
	   );
    
  }
  

  bool LoosePhotonIDSpring15(Photon const* photon,double const& rho) {//function for spring15 ID
    double eta = fabs(photon->eta());
    bool in_barrel = true;
    if (eta > 1.479) in_barrel = false;
    bool in_endcap = true;
    if (eta < 1.479) in_endcap = false;
    
    double photon_area, neutral_area, charged_area;
    if (eta < 1.0)                  { charged_area=0.0158; neutral_area=0.0143; photon_area=0.0725; }
    if (eta >= 1.0 && eta < 1.479)  { charged_area=0.0143; neutral_area=0.0210; photon_area=0.0604; }
    if (eta >= 1.479 && eta < 2.0)  { charged_area=0.0115; neutral_area=0.0147; photon_area=0.0320; }
    if (eta >= 2.0 && eta < 2.2)    { charged_area=0.0094; neutral_area=0.0082; photon_area=0.0512; }
    if (eta >= 2.2 && eta < 2.3)    { charged_area=0.0095; neutral_area=0.0124; photon_area=0.0766; }
    if (eta >= 2.3 && eta < 2.4)    { charged_area=0.0068; neutral_area=0.0186; photon_area=0.0949; }
    if (eta >= 2.4)                 { charged_area=0.0053; neutral_area=0.0320; photon_area=0.1160; }

    return(
	( (in_barrel       
	   && photon->had_tower_over_em()<0.05
	   && photon->sigma_IetaIeta()<0.0103
	   && std::max(photon->dr03_pfiso_charged()-rho*charged_area,0.)<2.44
	   && std::max(photon->dr03_pfiso_neutral()-rho*neutral_area,0.)<(2.57+exp(0.0044*photon->pt()+0.5809))
           && std::max(photon->dr03_pfiso_gamma()-rho*photon_area,0.)<(1.92+0.0043*photon->pt())
	      ) ||
        (in_endcap       
	 && photon->had_tower_over_em()<0.05
	 && photon->sigma_IetaIeta()<0.0277
	 && std::max(photon->dr03_pfiso_charged()-rho*charged_area,0.)<1.84
	 && std::max(photon->dr03_pfiso_neutral()-rho*neutral_area,0.)<(4+exp(0.004*photon->pt()+0.9402))
	 && std::max(photon->dr03_pfiso_gamma()-rho*photon_area,0.)<(2.15+0.0041*photon->pt())
        ) )
        );

  }

  bool PairFirstLeadingSecond(CompositeCandidate const* cand){
    if (cand->At(0)->pt() > cand->At(1)->pt()) return true;
    return false;
  }

  bool PairOneWithPt(CompositeCandidate const* cand, double const& ptMin) {
    std::vector<Candidate*> const& vec = cand->AsVector();
    for (unsigned i = 0; i < vec.size(); ++i) {
      if (vec[i]->pt() > ptMin) return true;
    }
    return false;
  }
  
  bool PairPtSelection(CompositeCandidate const* cand, double const& jetpt1, double const& jetpt2){
    if (cand->At(0)->pt()>jetpt1 && cand->At(1)->pt()>jetpt2) return true;
    return false;
  }

  bool OrderedPairPtSelection(CompositeCandidate const* cand, double const& jetpt1, double const& jetpt2, bool const& aboveorbelow){
    double at0pt=cand->At(0)->pt();
    double at1pt=cand->At(1)->pt();
    if(at0pt>at1pt){
      if(aboveorbelow){
  if (at0pt>jetpt1 && at1pt>jetpt2) return true;
      }
      else{
  if (at0pt<jetpt1 && at1pt<jetpt2) return true;
      }
      return false;
    }
    else{
      if(aboveorbelow){
  if (at1pt>jetpt1 && at0pt>jetpt2) return true;
      }
      else{
  if (at1pt<jetpt1 && at0pt<jetpt2) return true;
      }
      return false;
    }
  }

  bool PairMassInRange(CompositeCandidate const* cand, double const& mLow, double const& mHigh) {
    double mass = (cand->At(0)->vector() + cand->At(1)->vector()).M();
    return (mass > mLow && mass < mHigh);
  }
  
  double PairMass(CompositeCandidate const* cand) {
    double mass = (cand->At(0)->vector() + cand->At(1)->vector()).M();
    return (mass);
  }
  //Return abolute difference between mass of candidate pair and given mass (e.g. Higgs), for sorting by closest candidate pair
  double PairMassDiff(CompositeCandidate const* cand, double const& mass) {
    double pairmass = (cand->At(0)->vector() + cand->At(1)->vector()).M();
    return abs(pairmass-mass);
  }

  bool PairEtaProdLessThan(CompositeCandidate const* cand, double const& max){
    if ((cand->At(0)->eta() * cand->At(1)->eta()) < max) return true;
    return false;
  }

  bool PairDEtaLessThan(CompositeCandidate const* cand, double const& max){
    if (fabs(cand->At(0)->eta() - cand->At(1)->eta()) < max) return true;
    return false;
  }

  bool PairDPhiLessThan(CompositeCandidate const* cand, double const& max){
    if (ROOT::Math::VectorUtil::DeltaPhi(cand->At(0)->vector(),cand->At(1)->vector()) < max) return true;
    return false;
  }

  bool PairAbsDPhiLessThan(CompositeCandidate const* cand, double const& max){
    if (fabs(ROOT::Math::VectorUtil::DeltaPhi(cand->At(0)->vector(),cand->At(1)->vector())) < max) return true;
    return false;
  }

  bool PairOppCharge(CompositeCandidate const& cand) {
    int charge = (cand.At(0)->charge() * cand.At(1)->charge());
    return (charge == -1);
  }

  bool PairOppSign(CompositeCandidate const* cand) {
    int charge = (cand->At(0)->charge() * cand->At(1)->charge());
    return (charge == -1 && abs(cand->At(0)->charge()) == 1 && abs(cand->At(1)->charge()) == 1);
  }

  bool PairSameSign(CompositeCandidate const* cand) {
    int charge = (cand->At(0)->charge() * cand->At(1)->charge());
    return (charge == 1 && abs(cand->At(0)->charge()) == 1 && abs(cand->At(1)->charge()) == 1);
  }

  bool MuonTight(Muon const* muon) {
    bool tightCut = ( 
        muon->is_global() &&
	muon->is_pf() &&
        //muon->is_tracker() &&
        //muon->it_tracker_hits() > 10 &&
        muon->it_pixel_hits() > 0 &&
        muon->gt_valid_muon_hits() > 0 &&
        muon->it_layers_with_measurement() > 5 &&
        muon->gt_normalized_chi2() < 10 &&
        muon->matched_stations() > 1
        );
    return tightCut;
  }

  //Defined in https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonId2015
  bool MuonMedium(Muon const* muon) {
    bool goodGlob = muon->is_global() && 
      muon->gt_normalized_chi2() < 3 &&
      muon->cq_chi2_localposition()<12 &&
      muon->cq_trk_kink()<20;

    bool isMedium =
      muon->is_pf() &&
      (muon->is_global()||muon->is_tracker()) && //Require loose muon except pf isolation should be done in individual analyses
      muon->it_valid_fraction() > 0.8 && 
      muon->segment_compatibility() > (goodGlob ? 0.303 : 0.451); 
    return isMedium;
  }

  bool MuonMediumHIPsafe(Muon const* muon) {
    bool goodGlob = muon->is_global() && 
      muon->gt_normalized_chi2() < 3 &&
      muon->cq_chi2_localposition()<12 &&
      muon->cq_trk_kink()<20;

    bool isMedium =
      muon->is_pf() &&
      (muon->is_global()||muon->is_tracker()) && //Require loose muon except pf isolation should be done in individual analyses
      muon->it_valid_fraction() > 0.49 && 
      muon->segment_compatibility() > (goodGlob ? 0.303 : 0.451); 
    return isMedium;
  }


  bool MuonLoose(Muon const* muon) {
    bool isLoose = 
      muon->is_pf() &&
      (muon->is_global()||muon->is_tracker());
    return isLoose;
  }

  bool HttEMuFakeMuon(Muon const* muon) {
    bool pass_iso = false;
    if ( muon->pt() > 20. &&
         (muon->dr03_tk_sum_pt()/muon->pt())          < 0.4 &&
         (muon->dr03_ecal_rechit_sum_et()/muon->pt()) < 0.4 &&
         (muon->dr03_hcal_tower_sum_et()/muon->pt())  < 0.4 ) pass_iso = true;
    if (  muon->pt() <= 20. &&
          muon->dr03_tk_sum_pt()          < 8. &&
          muon->dr03_ecal_rechit_sum_et() < 8. &&
          muon->dr03_hcal_tower_sum_et()  < 8.) pass_iso = true;
    bool pass = ( 
        muon->is_global() &&  // isGlobalMuon()
        // muon->gt_valid_muon_hits() > 0 &&  // globalTrack()->hitPattern().numberOfValidMuonHits()
        // muon->gt_normalized_chi2() < 10 && // globalTrack()->normalizedChi2()
        // muon->matched_stations() > 1 && // numberOfMatchedStations()
        fabs(muon->dxy_vertex()) < 0.2 &&
        pass_iso
        );
    return (pass && !( MuonTight(muon) && PF04IsolationEB(muon, 0.5, 0.15, 0.1) && fabs(muon->dxy_vertex()) < 0.02 && fabs(muon->dz_vertex()) < 0.1 )); // !(usual tight muon && pf iso < 0.15(0.1) barrel(endcap))
  }

  bool HttEMuFakeMuon2011(Muon const* muon) {
    bool pass_iso = false;
    if ( (muon->dr03_tk_sum_pt()/muon->pt())          < 0.2 &&
         (muon->dr03_ecal_rechit_sum_et()/muon->pt()) < 0.2 &&
         (muon->dr03_hcal_tower_sum_et()/muon->pt())  < 0.2 ) pass_iso = true;
    
    bool pass = ( 
        muon->is_global() &&  // isGlobalMuon()
        muon->gt_valid_muon_hits() > 0 &&  // globalTrack()->hitPattern().numberOfValidMuonHits()
        muon->gt_normalized_chi2() < 10 && // globalTrack()->normalizedChi2()
        muon->matched_stations() > 1 && // numberOfMatchedStations()
        fabs(muon->dxy_vertex()) < 0.2 &&
        fabs(muon->dz_vertex()) < 0.1 &&
        pass_iso
        );
    return (pass && !( MuonTight(muon) && PF04IsolationEB(muon, 0.5, 0.15, 0.1) && fabs(muon->dxy_vertex()) < 0.02 && fabs(muon->dz_vertex()) < 0.1 )); // !(usual tight muon && pf iso < 0.15(0.1) barrel(endcap))
  }


  bool MuonIso(Muon const* muon) {
    bool isoCut = (
        ((muon->dr03_tk_sum_pt() + 
          muon->dr03_ecal_rechit_sum_et() + 
          muon->dr03_hcal_tower_sum_et()) 
	 / muon->pt()) < 0.15);
    return isoCut;
  }

  bool MuonTkIso(Muon const* muon) {
    bool isoCut = (((muon->dr03_tk_sum_pt())/muon->pt())<0.1); 
    return isoCut;
  }

  double MuonTkIsoVal(Muon const* muon) {
    double isoCut = ((muon->dr03_tk_sum_pt())/muon->pt()); 
    return isoCut;
  }

  
  std::vector<Track *> GetTracksAtVertex(std::vector<Track *> const& trks, std::vector<Vertex *> const& vtxs, unsigned idx, double const& dz) {
    std::vector<Track *> result;
    std::vector<std::size_t> vtx_tracks = ExtractFirst(vtxs.at(idx)->tracks());
    for (unsigned i = 0; i < trks.size(); ++i) {
      //Is track in vtx fit -> keep
      if (std::find(vtx_tracks.begin(),vtx_tracks.end(),trks[i]->id()) != vtx_tracks.end()) {
        result.push_back(trks[i]);
        continue;
      }
      //If not, is track within dz and is it the shortest distance?
      double dist = fabs(trks[i]->vz() - vtxs.at(idx)->vz());
      if (dist < dz) {
        std::vector<double> dists(vtxs.size(),0.0);
        for (unsigned j = 0; j < vtxs.size(); ++j) {
          dists[j] = fabs(trks[i]->vz() - vtxs.at(j)->vz());
        }
        unsigned min_ele = unsigned(std::min_element(dists.begin(),dists.end()) - dists.begin());
        if (min_ele == idx) result.push_back(trks[i]);
      }
    }
    return result;
  }

  bool MassDiffCompare(Candidate const* p1, Candidate const* p2, 
      double const& mass) {
      double m1 = p1->vector().M();
      double m2 = p2->vector().M();
      double md1 = fabs(m1 - mass);
      double md2 = fabs (m2 - mass);
      return (md1 < md2);
  }

  std::vector<GenParticle *> ExtractStableDaughters(GenParticle * part, std::vector<GenParticle *> const& input) {
    std::vector<GenParticle *> tmp = ExtractDaughters(part, input);
    std::set<GenParticle *> result_set;
    bool all_found = false;
    while (!all_found) {
      std::vector<GenParticle *> tmp_other;
      unsigned n_st1 = 0;
      for (unsigned i = 0; i < tmp.size(); ++i) {
        if (tmp[i]->status() == 1) {
          ++n_st1;
          result_set.insert(tmp[i]);
        } else {
          std::vector<GenParticle *> tmp_daughters = ExtractDaughters(tmp[i], input);
          for(unsigned j = 0; j < tmp_daughters.size(); ++j) tmp_other.push_back(tmp_daughters[j]);
        }
      }
      tmp = tmp_other;
      if (tmp.size() == 0) all_found = true;
    }
    std::vector<GenParticle *> result;
    for (auto it : result_set) result.push_back(it);
    return result;
  }

  std::vector<GenParticle*> ExtractDaughters(
      GenParticle* part, std::vector<GenParticle*> const& input) {
    std::vector<GenParticle*> result;
    std::vector<int> daughters = part->daughters();
    for (unsigned i = 0; i < input.size(); ++i) {
      if (std::find(daughters.begin(), daughters.end(), input[i]->index()) !=
          daughters.end()) {
        result.push_back(input[i]);
      }
    }
    return result;
  }

  std::vector<GenParticle*> ExtractMothers(
      GenParticle* part, std::vector<GenParticle*> const& input) {
    std::vector<GenParticle*> result;
    std::vector<int> mothers = part->mothers();
    for (unsigned i = 0; i < input.size(); ++i) {
      if (std::find(mothers.begin(), mothers.end(), input[i]->index()) !=
          mothers.end()) {
        result.push_back(input[i]);
      }
    }
    return result;
  }


  std::vector<GenParticle*> ExtractDaughtersRecursive(
      GenParticle* part, std::vector<GenParticle*> const& input) {
    std::vector<GenParticle*> tmp = ExtractDaughters(part, input);
    std::set<GenParticle*> result_set;
    bool all_found = false;
    while (!all_found) {
      std::vector<GenParticle*> tmp_other;
      for (unsigned i = 0; i < tmp.size(); ++i) {
        result_set.insert(tmp[i]);
        std::vector<GenParticle*> tmp_daughters =
            ExtractDaughters(tmp[i], input);
        for (unsigned j = 0; j < tmp_daughters.size(); ++j)
          tmp_other.push_back(tmp_daughters[j]);
      }
      tmp = tmp_other;
      if (tmp.size() == 0) all_found = true;
    }
    std::vector<GenParticle*> result;
    for (auto it : result_set) result.push_back(it);
    std::sort(result.begin(), result.end(),
              [](GenParticle const* p1,
                 GenParticle const* p2) { return p1->index() < p2->index();
              });
    return result;
  }


  std::vector<GenJet> BuildTauJets(std::vector<GenParticle *> const& parts, bool include_leptonic, bool use_prompt) {
    // modified this to set the jet flavour as a flag indicating the tau decay mode, flag is integer defined as:
    // Warning: returned tau type works for taus decayed by Pythia8 but probably won't work for other generators e.g tauola!
    std::vector<GenJet> taus;
    for (unsigned i = 0; i < parts.size(); ++i) {
        std::vector<bool> status_flags;
        bool is_prompt=true; 
        if(use_prompt){
          status_flags = parts[i]->statusFlags();
          is_prompt=status_flags[IsPrompt];
        }
      if (abs(parts[i]->pdgid()) == 15 && is_prompt) {
        std::size_t tau_id = parts[i]->id();
        std::vector<GenParticle *> daughters = ExtractDaughters(parts[i], parts);
        bool has_tau_daughter = false;
        bool has_lepton_daughter = false;
        unsigned count_pi0 =  0, count_pi = 0, count_tot = 0, count_k = 0, count_rho = 0;
        ROOT::Math::PtEtaPhiEVector nuvec;
        for (unsigned j = 0; j < daughters.size(); ++j) {
          if (abs(daughters[j]->pdgid()) == 15) has_tau_daughter = true;
          if (abs(daughters[j]->pdgid()) == 11 || abs(daughters[j]->pdgid()) == 13) has_lepton_daughter = true;
          unsigned pdgId = abs(daughters[j]->pdgid());
          if(pdgId == 111) count_pi0++;
          if(pdgId == 211) count_pi++;
          if(pdgId == 213) count_rho++;
          if(pdgId == 321) count_k++;
          if(pdgId!=12&&pdgId!=14&&pdgId!=16) {
            count_tot++;
          } else {
            nuvec+=daughters[j]->vector();
          }
        }
        if (has_tau_daughter) continue;
        if (has_lepton_daughter && !include_leptonic) continue;
        int tauFlag = -1;
        if(count_tot==1 && count_pi==1 && count_pi0==0) tauFlag=0;
        if(count_tot==1 && count_pi==0 && count_k==1 && count_pi0==0) tauFlag=20;
        if(count_tot==2 && count_pi==1 && count_pi0==1) tauFlag=1; 
        if(count_tot==2 && count_rho==1) tauFlag=21; 
        if(count_tot==3 && count_pi==1 && count_pi0==2) tauFlag=2;
        if(count_tot==3 && count_pi==3 && count_pi0==0) tauFlag=10;
        if(count_tot==4 && count_pi==3 && count_pi0==1) tauFlag=11;
        std::vector<GenParticle *> jet_parts = ExtractStableDaughters(parts[i], parts);
        taus.push_back(GenJet());
        ROOT::Math::PtEtaPhiEVector vec;
        std::vector<std::size_t> id_vec;
        for (unsigned k = 0; k < jet_parts.size(); ++k) {
          if (  abs(jet_parts[k]->pdgid()) == 12 || 
                abs(jet_parts[k]->pdgid()) == 14 ||
                abs(jet_parts[k]->pdgid()) == 16 ) continue;
          vec += jet_parts[k]->vector();
          taus.back().set_charge(taus.back().charge() + jet_parts[k]->charge());
          id_vec.push_back(jet_parts[k]->id());
        }
        taus.back().set_vector(vec);
        taus.back().set_constituents(id_vec);
        taus.back().set_flavour(tauFlag);
        taus.back().set_id(tau_id);
        taus.back().set_nu_vector(nuvec);
   
      }
    }
    return taus;
  }

  std::vector<ic::PFCandidate*> GetTauGammas(ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands) {
    return GetTauGammas(tau, pfcands, 0.5, 0);
  }

  std::vector<ic::PFCandidate*> GetTauGammas(ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands, double pt_cut, int gammas_shift=0) {
    // shift gammas energy down: gammas_shift = 1
    // shift gammas energy up: gammas_shift = 2
    // No shift of gammas energy: gammas_shift = 0
    std::vector<ic::PFCandidate*> gammas = {}; 
    std::vector<std::size_t> gammas_ids = tau->sig_gamma_cands();
    for(auto id : gammas_ids){
      for(auto p : pfcands) {
        std::size_t pfid = p->id();
        if(pfid == id) {
          ic::PFCandidate *gam = new ic::PFCandidate(*p);
          if(gammas_shift!=0) {
            double uncert = sqrt(0.0009/gam->energy()+0.000001);
            double shift=1.;
            if(gammas_shift==1)      shift-=uncert;
            else if(gammas_shift==2) shift+=uncert;
            gam->set_pt(gam->pt()*shift);
            gam->set_energy(gam->energy()*shift);
          } 
          if(p->pt()>pt_cut) gammas.push_back(gam);
          break;
        }
      }
    }
    std::sort(gammas.begin(), gammas.end(), bind(&PFCandidate::pt, _1) > bind(&PFCandidate::pt, _2));
    return gammas;
  }


  std::vector<ic::PFCandidate*> GetTauIsoGammas(ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands, double pt_cut, int gammas_shift) {
    std::vector<ic::PFCandidate*> gammas = {};
    std::vector<std::size_t> gammas_ids = tau->iso_gamma_cands();
    for(auto id : gammas_ids){
      for(auto p : pfcands) {
        std::size_t pfid = p->id();
        if(pfid == id) {
          ic::PFCandidate *gam = new ic::PFCandidate(*p);
          if(gammas_shift!=0) {
            double uncert = sqrt(0.0009/gam->energy()+0.000001);
            double shift=1.;
            if(gammas_shift==1)      shift-=uncert;
            else if(gammas_shift==2) shift+=uncert;
            gam->set_pt(gam->pt()*shift);
            gam->set_energy(gam->energy()*shift);
          }

          if(gam->pt()>pt_cut) gammas.push_back(gam);
          break;
        }
      }
    }
    std::sort(gammas.begin(), gammas.end(), bind(&PFCandidate::pt, _1) > bind(&PFCandidate::pt, _2));
    return gammas;
  }

  std::vector<ic::PFCandidate*> GetTauIsoGammas(ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands) {
    return GetTauIsoGammas(tau, pfcands, 0.5, 0);
  }


  
  std::vector<ic::PFCandidate*> GetTauHads(ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands, int pi_shift) {
    std::vector<ic::PFCandidate*> hads = {};
    std::vector<std::size_t> hads_ids = tau->sig_charged_cands();
    for(auto id : hads_ids){
      for(auto p : pfcands) {
        std::size_t pfid = p->id();
        if(pfid == id) {
          ic::PFCandidate *pi = new ic::PFCandidate(*p);
          if(pi_shift!=0) {
            double uncert = sqrt(pow(0.00009*pi->pt(),2)+pow(0.0085/sqrt(sin(2*atan(exp(-pi->eta())))),2));
            double shift=1.;
            if(pi_shift==1)      shift-=uncert;
            else if(pi_shift==2) shift+=uncert;
            p->set_pt(pi->pt()*shift);
            p->set_energy(pi->energy()*shift);
          }

          hads.push_back(pi);
          break;
        }
      }
    }

    std::sort(hads.begin(), hads.end(), bind(&PFCandidate::pt, _1) > bind(&PFCandidate::pt, _2));
    return hads;
  }

  std::pair<ic::Candidate*,ic::Candidate*> GetRho (ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands, double pt_cut, int uncert_mode) {
    // uncert_mode = 0: no shifts to gammas or pi
    // uncert_mode = 1: shift gammas down
    // uncert mode = 2: shift gammas up
    // uncert_mode = 3: shift pi down
    // uncert_mode = 4: shift pi up
    ic::Candidate *pi = new ic::Candidate();
    ic::Candidate *pi0 = new ic::Candidate();
    int pi_uncert = 0;
    if(uncert_mode==3) pi_uncert = 1;
    if(uncert_mode==4) pi_uncert = 2;
    int gammas_uncert = 0;
    if(uncert_mode==1) gammas_uncert = 1;
    if(uncert_mode==2) gammas_uncert = 2;
    std::vector<ic::PFCandidate*> hads = GetTauHads(tau, pfcands, pi_uncert);
    if(hads.size()>0) pi = (ic::Candidate*)hads[0];
    std::vector<ic::PFCandidate*> gammas = GetTauGammas(tau, pfcands, pt_cut, gammas_uncert);
    // reconstruct strips from "signal" gammas
    double mass = 0.1349;     
    double cone_size = std::max(std::min(0.1, 3./tau->pt()),0.05); 
    std::vector<std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>>> strip_pairs = HPSGammas(gammas, 0, 0, 0, mass, 1);
    std::vector<std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>>> strips_incone;
    for(auto s : strip_pairs) {
      if(std::fabs(ROOT::Math::VectorUtil::DeltaR(s.first->vector(),tau->vector()))<cone_size) strips_incone.push_back(s);
    }
    if(tau->decay_mode()==0) {
      std::vector<std::size_t> signal_gammas = {}; 
      if(strips_incone.size()>0) {
        pi0 = (ic::Candidate*)GetPi0(strips_incone[0].second, true); 
        for (auto g : strips_incone[0].second) signal_gammas.push_back(g->id());
      } else if(strip_pairs.size()>0) {
        pi0 = (ic::Candidate*)GetPi0(strip_pairs[0].second, true);
        for (auto g : strip_pairs[0].second) signal_gammas.push_back(g->id());
      }
      Tau * t = const_cast<Tau*>(tau);
      t->set_sig_gamma_cands(signal_gammas);
    } else if(tau->decay_mode()==1 && strip_pairs.size()>0) {
      pi0 = (ic::Candidate*)GetPi0(strip_pairs[0].second, true);
      //pi0 = (ic::Candidate*)GetPi0(gammas, true); 
    } //else if(tau->decay_mode()==1) {
      //// energy of pi0 from sum of gammas, direction from lead gamma
      //pi0 = (ic::Candidate*)GetPi0(gammas, true);
    //}

    return std::make_pair(pi,pi0);
  }

  std::pair<ic::Candidate*,ic::Candidate*> GetRho (ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands) {
    return  GetRho (tau, pfcands, 0.5,0);
  }
 

  bool sortStrips (std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>> i,std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>> j) {
    return (i.first->pt() > j.first->pt());
  }

  std::vector<std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>>> HPSGammas (std::vector<ic::PFCandidate*> cands, double stripPtThreshold, double etaAssociationDistance, double phiAssociationDistance, double mass, unsigned mode) {

    // mode = 0 uses fixed strip sizes set by etaAssociationDistance and phiAssociationDistance
    // mode = 1 uses dynamic strip size
    // mode = 2 uses dynamic strip size but sets maximum sizes of strips by etaAssociationDistance and phiAssociationDistance
    std::vector<std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>>> strips;  
    while(!cands.empty()) {
      //save the non associated candidates to a different collection
      std::vector<ic::PFCandidate*> Associated = {};
      std::vector<ic::PFCandidate*> notAssociated = {};

      //Create a cluster from the Seed Photon
      ROOT::Math::PtEtaPhiEVector stripVector(0,0,0,0);
      stripVector=cands.at(0)->vector();
      Associated.push_back(cands.at(0));
      //Loop and associate
      bool repeat = true;
      unsigned int mini=1;
      while (repeat) {
        repeat = false;
        for(unsigned int i=mini;i<cands.size();++i) {
          if(mode==1) {
            etaAssociationDistance = 0.20*pow(cands[i]->pt(),-0.66) + 0.20*pow(stripVector.Pt(),-0.66);
            phiAssociationDistance = 0.35*pow(cands[i]->pt(),-0.71) + 0.35*pow(stripVector.Pt(),-0.71);
            etaAssociationDistance = std::min(etaAssociationDistance, 0.15);
            etaAssociationDistance = std::max(etaAssociationDistance, 0.05);
            phiAssociationDistance = std::min(phiAssociationDistance, 0.30);
            phiAssociationDistance = std::max(phiAssociationDistance, 0.05);
          } else if (mode==2) {
            etaAssociationDistance = std::min(0.20*pow(cands[i]->pt(),-0.66) + 0.20*pow(stripVector.Pt(),-0.66), etaAssociationDistance);
            phiAssociationDistance = std::min(0.35*pow(cands[i]->pt(),-0.71) + 0.35*pow(stripVector.Pt(),-0.71), phiAssociationDistance);
            etaAssociationDistance = std::min(etaAssociationDistance, 0.15);
            phiAssociationDistance = std::min(phiAssociationDistance, 0.30);
          }
          if(fabs(cands.at(i)->eta()-stripVector.eta())<etaAssociationDistance &&
            fabs(ROOT::Math::VectorUtil::DeltaPhi(cands.at(i)->vector(),stripVector))<phiAssociationDistance) {
            stripVector+=cands.at(i)->vector();
            Associated.push_back(cands.at(i));
            repeat = true;
          }
          else {
            notAssociated.push_back(cands.at(i));
          }
        }
        //Swap the candidate vector with the non associated vector
        cands.swap(notAssociated);
        //Clear
        notAssociated.clear(); 
        mini=0;
      }
      //Save the strip 
      ic::PFCandidate strip;
      strip.set_vector(stripVector);

      strip.set_vector(GetPi0(Associated, false)->vector());
      if(mass) {;}
      if(strip.pt()>=stripPtThreshold) strips.push_back(std::make_pair(new ic::PFCandidate(strip), Associated));
    }
    std::sort(strips.begin(), strips.end(), sortStrips);
    return strips;
  }


  std::vector<ic::PFCandidate*> HPS (std::vector<ic::PFCandidate*> cands, double stripPtThreshold, double etaAssociationDistance, double phiAssociationDistance, double mass, unsigned mode) {
    std::vector<std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>>> strip_pairs = HPSGammas(cands, stripPtThreshold, etaAssociationDistance, phiAssociationDistance, mass, mode);
    std::vector<ic::PFCandidate*> strips;
    for(auto s : strip_pairs) strips.push_back(s.first);
    return strips;
  }

  ic::Candidate*  GetPi0 (std::vector<ic::PFCandidate*> gammas, bool leadEtaPhi) {
    ic::Candidate *pi0 = new ic::Candidate;
    if(gammas.size()>0) {
      double E = 0.;
      double phi = 0.;
      double eta = 0.;
      for(auto g: gammas) {
        E+=g->energy();
        phi+=g->energy()*g->phi();
        eta+=g->energy()*g->eta();
      }
      eta/=E;
      phi/=E;

      if(leadEtaPhi){
        // if true sets the eta and phi of the pi0 to that of the leading gamma - improves resolution
        eta = gammas[0]->eta();
        phi = gammas[0]->phi(); 
      }

      double mass = 0.1349;
      double p = sqrt(E*E-mass*mass);
      double theta = atan(exp(-eta))*2;
      double pt = p*sin(theta);
      ROOT::Math::PtEtaPhiEVector pi0_vector(pt,eta,phi,E);
      pi0->set_vector(pi0_vector);
    }
    return pi0;
  }

  std::pair<std::vector<ic::PFCandidate*>, ic::Candidate*> GetA1 (ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands, double pt_cut, int uncert_mode) {
    // uncert_mode = 0: no shifts to gammas or pi
    // uncert_mode = 1: shift gammas down
    // uncert mode = 2: shift gammas up
    // uncert_mode = 3: shift pi down
    // uncert_mode = 4: shift pi up
    int pi_uncert = 0;
    if(uncert_mode==3) pi_uncert = 1;
    if(uncert_mode==4) pi_uncert = 2;
    int gammas_uncert = 0;
    if(uncert_mode==1) gammas_uncert = 1;
    if(uncert_mode==2) gammas_uncert = 2;

    std::vector<ic::PFCandidate*> hads = GetTauHads(tau, pfcands, pi_uncert);
    if(hads.size()==3) {
      // arrange hadrons so the oppositly charged hadron is contained in the first element
      if(hads[1]->charge()!=hads[0]->charge()&&hads[1]->charge()!=hads[2]->charge()){
        auto temp = hads[1];
        hads[1] = hads[0];
        hads[0] = temp;
      }
      else if(hads[2]->charge()!=hads[0]->charge()&&hads[2]->charge()!=hads[1]->charge()){
        auto temp = hads[2];
        hads[2] = hads[0];
        hads[0] = temp;
      } 
      // from the two same sign hadrons place the one that gives the mass most similar to the rho meson as the second element
      double rho_mass = 0.7755;
      double dM1 = std::fabs((hads[0]->vector()+hads[1]->vector()).M()-rho_mass);
      double dM2 = std::fabs((hads[0]->vector()+hads[2]->vector()).M()-rho_mass);
      if(dM2<dM1){
        auto temp = hads[2];
        hads[2] = hads[1];
        hads[1] = temp;
      }
    }
    std::vector<ic::PFCandidate*> gammas;
    if(tau->decay_mode()>10) gammas = GetTauGammas(tau, pfcands, pt_cut, gammas_uncert);
    else  {
      // need to modify isolation collection to exclude these gammas!
      gammas = GetTauGammas(tau, pfcands, pt_cut,gammas_uncert);//signal gammas
      std::vector<ic::PFCandidate*> iso_gammas = GetTauIsoGammas(tau, pfcands, pt_cut, gammas_uncert);//iso gammas 
      gammas.insert(gammas.end(), iso_gammas.begin(), iso_gammas.end()  ); 
      std::sort(gammas.begin(), gammas.end(), bind(&PFCandidate::pt, _1) > bind(&PFCandidate::pt, _2));
    }
    double cone_size = std::max(std::min(0.1, 3./tau->pt()),0.05);
    double mass = 0.1349;
    ic::Candidate *pi0 = new ic::Candidate();
    ic::erase_if(gammas,!boost::bind(MinPtMaxEta, _1, pt_cut, 9999.));
    std::vector<std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>>> strip_pairs = HPSGammas(gammas, 0, 0, 0, mass, 1); 
    std::vector<std::pair<ic::PFCandidate*,std::vector<ic::PFCandidate*>>> strips_incone; 
    for(auto s : strip_pairs) {
      if(std::fabs(ROOT::Math::VectorUtil::DeltaR(s.first->vector(),tau->vector()))<cone_size) strips_incone.push_back(s);
    }
    std::vector<std::size_t> signal_gammas = {};
    if(strips_incone.size()>0) {
      pi0 = (ic::Candidate*)GetPi0(strips_incone[0].second, true);
      for (auto g : strips_incone[0].second) signal_gammas.push_back(g->id());
    } else if(strip_pairs.size()>0) {
      pi0 = (ic::Candidate*)GetPi0(strip_pairs[0].second, true);
      for (auto g : strip_pairs[0].second) signal_gammas.push_back(g->id());
    } 

    Tau * t = const_cast<Tau*>(tau);
    t->set_sig_gamma_cands(signal_gammas);
    std::vector<std::size_t> iso_gammas = tau->iso_gamma_cands();
    std::vector<std::size_t> new_iso_gammas = {};
    for(auto g : iso_gammas) {
      if(std::find(signal_gammas.begin(), signal_gammas.end(), g) == signal_gammas.end()) new_iso_gammas.push_back(g);  
    }
    t->set_iso_gamma_cands(new_iso_gammas);
    return std::make_pair(hads, pi0);
  }

  std::pair<std::vector<ic::PFCandidate*>, ic::Candidate*> GetA1 (ic::Tau const* tau, std::vector<ic::PFCandidate*> pfcands) { 
    return  GetA1 (tau, pfcands, 0.5, 0); 
  }


  std::vector<GenParticle*> GetTauDaughters(std::vector<GenParticle *> const& parts, std::vector<std::size_t> id) {
    std::vector<GenParticle*> tau_daughters;
    for (unsigned i = 0; i < parts.size(); ++i) {
      for (unsigned j = 0; j < id.size(); ++j) {
        if(parts[i]->id() == id[j]){
          tau_daughters.push_back(parts[i]);
          continue;
        }
      }
    }
    return tau_daughters;
  }
  
    std::pair<bool, GenParticle*> GetTauPiDaughter(std::vector<GenParticle *> const& parts, std::vector<std::size_t> id) {
    std::vector<GenParticle*> tau_daughters;
    for (unsigned i = 0; i < parts.size(); ++i) {
      for (unsigned j = 0; j < id.size(); ++j) {
        if(parts[i]->id() == id[j]){
          tau_daughters.push_back(parts[i]);
          continue;
        }
      }
    }
    GenParticle* pi = new GenParticle();
    pi->set_pdgid(0);
    GenParticle* tau_rho_daughter = new GenParticle();
    int countpi = 0;
    if(tau_daughters.size()==1){
      for(unsigned i=0; i<tau_daughters.size(); ++i){
        if(fabs(tau_daughters[i]->pdgid()) == 211) { pi = tau_daughters[i]; countpi++;}
      }
    }
    bool isPi = countpi == 1;
    if (isPi) tau_rho_daughter = pi;

    return std::make_pair(isPi, tau_rho_daughter);
  }
  
    std::pair<bool, std::vector<GenParticle*>> GetTauRhoDaughter(std::vector<GenParticle *> const& parts, std::vector<std::size_t> id) {
    std::vector<GenParticle*> tau_daughters;
    for (unsigned i = 0; i < parts.size(); ++i) {
      for (unsigned j = 0; j < id.size(); ++j) {
        if(parts[i]->id() == id[j]){
          tau_daughters.push_back(parts[i]);
          continue;
        }
      }
    }
    GenParticle* pi0 = new GenParticle();
    GenParticle* pi = new GenParticle();
    pi->set_pdgid(0);
    pi0->set_pdgid(0);
    std::vector<GenParticle*> tau_rho_daughter;
    int countpi = 0;
    int countgamma = 0;
    if(tau_daughters.size()==3){
      for(unsigned i=0; i<tau_daughters.size(); ++i){
        if(fabs(tau_daughters[i]->pdgid()) == 22){ pi0->set_vector(pi0->vector()+tau_daughters[i]->vector()); pi0->set_pdgid(111); countgamma++;}
        if(fabs(tau_daughters[i]->pdgid()) == 211) { pi = tau_daughters[i]; countpi++;}
      }
    }
    bool isRho = countgamma == 2 && countpi == 1;
    if (isRho) tau_rho_daughter = {pi,pi0};
    return std::make_pair(isRho, tau_rho_daughter);
    }

    // new function for rho daughters
    /* std::vector<std::pair<ic::GenParticle*>> GetTauRhoDaughterFromGenParticles */
    /*     (std::vector<ic::GenParticle *> genParticles, std::vector<std::size_t> id) { */
    /*   unsigned count_taus=0; */
    /*   for (unsigned i = 0; i < genParticles.size(); ++i) { */
    /*     if(std::fabs(genParticles[i]->pdgid()) == 15 && genParticles[i]->statusFlags()[IsLastCopy]){ */
    /*       ic::GenParticle* tau = gen_particles[i]; */
    /*       std::cout << tau->daughters() << std::endl; */
    /*       break; */
    /*       std::cout << "" << std::endl; */
    /*       ++count_taus; */
    /*       bool isRho = false; */
    /*       std::pair<ic::GenParticle*,ic::GenParticle*> rho = */ 
    /*           std::make_pair(new ic::GenParticle(), new ic::GenParticle()); */
        
    /*       /1* for (unsigned d : tau->daughters()){ *1/ */
    /*       /1*   unsigned daughter_id = std::fabs(gen_particles[d]->pdgid()); *1/ */
    /*       /1*   if (daughter_id == 211) { *1/ */
    /*       /1*     rho.first = gen_particles[d]; *1/ */
    /*       /1*     ++count_pi; *1/ */
    /*       /1*     continue; *1/ */
    /*       /1*   } *1/ */
    /*       /1*   if (daughter_id == 111) { *1/ */
    /*       /1*     gammasP4 += gen_particles[d]->p4(); *1/ */
    /*       /1*     ++count_gamma; *1/ */
    /*       /1*     continue; *1/ */
    /*       /1*   } *1/ */
    /*       /1* } *1/ */
    /*       /1* isRho = (count_gamma==2 && count_pi==1); *1/ */
    /*       /1* rho.second = *1/ */ 

    /*     } */
    /*   } */
    /* } */

  
    std::pair<bool, std::vector<GenParticle*>> GetTauA1Daughter(std::vector<GenParticle *> const& parts, std::vector<std::size_t> id) {
    std::vector<GenParticle*> tau_daughters;
    for (unsigned i = 0; i < parts.size(); ++i) {
      for (unsigned j = 0; j < id.size(); ++j) {
        if(parts[i]->id() == id[j]){
          tau_daughters.push_back(parts[i]);
          continue;
        }
      }
    }
    GenParticle* pi1 = new GenParticle();
    GenParticle* pi2 = new GenParticle();
    GenParticle* pi3 = new GenParticle();
    pi1->set_pdgid(0);
    pi2->set_pdgid(0);
    pi3->set_pdgid(0);
    std::vector<GenParticle*> tau_a1_daughters;
    int countpiplus = 0;
    int countpiminus = 0;
    std::vector<GenParticle*> charged_pis;
    if(tau_daughters.size()==3){
      for(unsigned i=0; i<tau_daughters.size(); ++i){
        if(fabs(tau_daughters[i]->pdgid()) == 211) charged_pis.push_back(tau_daughters[i]);
        if(tau_daughters[i]->pdgid() == -211) countpiplus++;
        if(tau_daughters[i]->pdgid() == 211)  countpiminus++;
      }
    }
    bool isA1 = ( countpiplus==2 && countpiminus==1 ) || ( countpiplus==1 && countpiminus==2 );
    if (isA1) tau_a1_daughters = charged_pis;
    return std::make_pair(isA1, tau_a1_daughters);
  }
  
  TVector3 GetGenImpactParam (ic::Vertex primary_vtx, ic::Vertex secondary_vtx, ROOT::Math::PtEtaPhiEVector part_vec){
    TVector3 x(secondary_vtx.vx()-primary_vtx.vx(),secondary_vtx.vy()-primary_vtx.vy(),secondary_vtx.vz()-primary_vtx.vz());
    TVector3 unit_vec = ConvertToTVector3(part_vec).Unit();
    TVector3 u = -(x - x.Dot(unit_vec)*unit_vec);
    return u;    
  }
  
  TLorentzVector ConvertToLorentz(ROOT::Math::PtEtaPhiEVector input_vec){
      TLorentzVector out_vec;
      out_vec.SetXYZM(input_vec.Px(),input_vec.Py(),input_vec.Pz(),input_vec.M());
      return out_vec;
  }

  TVector3 ConvertToTVector3(ROOT::Math::PtEtaPhiEVector input_vec){
      TVector3 out_vec;
      out_vec.SetXYZ(input_vec.Px(),input_vec.Py(),input_vec.Pz());
      return out_vec;
  }

  ROOT::Math::PtEtaPhiEVector ConvertToPtEtaPhiEVector(TLorentzVector input_vec){
      ROOT::Math::PtEtaPhiEVector out_vec;
      out_vec.SetPt(input_vec.Pt());
      out_vec.SetEta(input_vec.Eta());
      out_vec.SetPhi(input_vec.Phi());
      out_vec.SetE(input_vec.E());
      return out_vec;
  }
  
  TLorentzVector SmearVectorVar(TLorentzVector input_vec, double smearVal, unsigned indexToSmear, double smearVal_Eta, double smearVal_Phi, double smearVal_E) {
    TLorentzVector out_vec = input_vec;
    if (indexToSmear == 1) 
      out_vec.SetPtEtaPhiE(input_vec.Pt()+(input_vec.P()*((1/std::cosh(input_vec.Eta() * smearVal))-(1/std::cosh(input_vec.Eta())))), input_vec.Eta() * smearVal, input_vec.Phi(), input_vec.E());
    if (indexToSmear == 2) 
      out_vec.SetPtEtaPhiE(input_vec.Pt(), input_vec.Eta(), input_vec.Phi() * smearVal, input_vec.E());
    if (indexToSmear == 3) 
      out_vec.SetPtEtaPhiE(input_vec.Pt() * smearVal, input_vec.Eta(), input_vec.Phi(), input_vec.E() * smearVal);

    if (indexToSmear == 0) {
      out_vec.SetPtEtaPhiE(input_vec.Pt()+(input_vec.P()*((1/std::cosh(input_vec.Eta() * smearVal_Eta))-(1/std::cosh(input_vec.Eta())))), input_vec.Eta() * smearVal_Eta, input_vec.Phi(), input_vec.E());
      out_vec.SetPtEtaPhiE(out_vec.Pt(), out_vec.Eta(), out_vec.Phi() * smearVal_Phi, out_vec.E());
      out_vec.SetPtEtaPhiE(out_vec.Pt() * smearVal_E, out_vec.Eta(), out_vec.Phi(), out_vec.E() * smearVal_E);
    }
    return out_vec;
  }

  double IPAcoSign(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool ZMF){
    //p1 = ip+, p2 = pi0-, p3 = pi+, p4 = pi-  

    TVector3 boost, n1, n2;
    if(ZMF) boost = (p1+p2+p3+p4).BoostVector();
    else boost = (p3+p4).BoostVector();
    p1.Boost(-boost);
    p2.Boost(-boost);
    p3.Boost(-boost);
    p4.Boost(-boost);

    n1 = p1.Vect() - p1.Vect().Dot(p3.Vect().Unit())*p3.Vect().Unit();
    n2 = p2.Vect() - p2.Vect().Dot(p4.Vect().Unit())*p4.Vect().Unit();

    n1 = n1.Unit();
    n2 = n2.Unit();

    double sign;
    sign = p4.Vect().Unit().Dot(n1.Cross(n2));

    return sign;
  }

  double IPAcoAngle(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool ZMF){

    TVector3 boost, n1, n2;
    if(ZMF) boost = (p1+p2+p3+p4).BoostVector();
    else boost = (p3+p4).BoostVector();

    p1.Boost(-boost);
    p2.Boost(-boost);
    p3.Boost(-boost);
    p4.Boost(-boost);
    
    n1 = p1.Vect() - p1.Vect().Dot(p3.Vect().Unit())*p3.Vect().Unit();    
    n2 = p2.Vect() - p2.Vect().Dot(p4.Vect().Unit())*p4.Vect().Unit();

    n1 = n1.Unit();
    n2 = n2.Unit();

    double angle = acos(n1.Dot(n2));
    double sign;

    sign = p4.Vect().Unit().Dot(n1.Cross(n2)); 
 
    if(sign<0) angle = 2*M_PI - angle;
    return angle;
  }

  double IPAcoAngle_Test(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool ZMF){

    TVector3 boost, n1, n2;

    if(ZMF) boost = (p1+p2+p3+p4).BoostVector();
    else boost = (p3+p4).BoostVector();
    p1.Boost(-boost);
    p2.Boost(-boost);
    p3.Boost(-boost);
    p4.Boost(-boost);


    n1 = p1.Vect() - p1.Vect().Dot(p3.Vect().Unit())*p3.Vect().Unit();
    n2 = p2.Vect() - p2.Vect().Dot(p4.Vect().Unit())*p4.Vect().Unit();

    n1 = n1.Unit();
    n2 = n2.Unit();

    double angle = acos(n1.Dot(n2));
    double sign;

    sign = p4.Vect().Unit().Dot(n1.Cross(n2));

    if(sign<0) angle = 2*M_PI - angle;
    return angle;
  }


  double IPAcoLinAngle(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool ZMF){
    TVector3 boost, n1, n2;
    if(ZMF) boost = (p1+p2+p3+p4).BoostVector();
    else boost = (p3+p4).BoostVector();
    p1.Boost(-boost);
    p2.Boost(-boost);
    p3.Boost(-boost);
    p4.Boost(-boost);

    n1 = p1.Vect().Dot(p3.Vect().Unit())*p3.Vect().Unit();
    n2 = p2.Vect().Dot(p4.Vect().Unit())*p4.Vect().Unit();

    n1 = n1.Unit();
    n2 = n2.Unit();

    double angle = acos(n1.Dot(n2));

    return angle;
  }

  double IPAngle(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool ZMF){
    TVector3 boost, n1, n2;
    if(ZMF) boost = (p1+p2+p3+p4).BoostVector();
    else boost = (p3+p4).BoostVector();
    p1.Boost(-boost);
    p2.Boost(-boost);
    p3.Boost(-boost);
    p4.Boost(-boost);

    n1 = p1.Vect();
    n2 = p2.Vect();

    n1 = n1.Unit();
    n2 = n2.Unit();

    double angle = acos(n1.Dot(n2));

    return angle;
  }

  std::pair<double,double> CTauAngle(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool ZMF){
    TVector3 boost, n1, n2 ,n3, n4;
    if(ZMF) boost = (p1+p2+p3+p4).BoostVector();
    else boost = (p3+p4).BoostVector();
    p1.Boost(-boost);
    p2.Boost(-boost);
    p3.Boost(-boost);
    p4.Boost(-boost);

    n1 = p1.Vect();
    n2 = p2.Vect();

    n1 = n1.Unit();
    n2 = n2.Unit();

    n3 = p3.Vect();
    n4 = p4.Vect();

    n3 = n3.Unit();
    n4 = n4.Unit();

    double c1 = n1.Dot(n3);
    double c2 = n2.Dot(n4);

    return std::make_pair(c1,c2);
  }

  std::pair<double,double> IPAcoAngleSep(TLorentzVector p1, TLorentzVector p2, TLorentzVector p3, TLorentzVector p4, bool ZMF){
    TVector3 boost, n1, n2;
    if(ZMF) boost = (p1+p2+p3+p4).BoostVector();
    else boost = (p3+p4).BoostVector();
    p1.Boost(-boost);
    p2.Boost(-boost);
    p3.Boost(-boost);
    p4.Boost(-boost);

    n1 = p1.Vect() - p1.Vect().Dot(p3.Vect().Unit())*p3.Vect().Unit();
    n2 = p2.Vect() - p2.Vect().Dot(p4.Vect().Unit())*p4.Vect().Unit();

    n1 = n1.Unit();
    n2 = n2.Unit();
    double angle = acos(n1.Dot(n2));
    double sign;
    sign = p4.Vect().Unit().Dot(n1.Cross(n2));

    if(sign<0) angle = 2*M_PI - angle;
    return std::make_pair(angle,angle);
  }





  double PolarimetricA1A1(TVector3 tau1, TVector3 tau2, TLorentzVector a1_1, TLorentzVector a1_2, std::vector<TLorentzVector> pis_1, std::vector<TLorentzVector> pis_2, std::vector<double> charges_1, std::vector<double> charges_2) {
    double angle = -9999.;
    TLorentzVector Tauminus, Tauplus;
    std::vector<TLorentzVector> solutions;
    TLorentzVector dummy1, dummy2;

    solutions=tauPairMomentumSolutions(tau1, a1_1, dummy1, true, tau2, a1_2, dummy2, true,true);

    Tauminus=solutions.at(3);
    Tauplus=solutions.at(7);

    SCalculator Scalc("a1");
    if(Scalc.isOk("a1", "a1", Tauminus, pis_1, charges_1, Tauplus, pis_2, charges_2)) //(security checks, if there is no problem with variables)
    {
       angle= Scalc.AcopAngle("a1", "a1", Tauminus, pis_1, charges_1, Tauplus, pis_2, charges_2);
    }

    return angle;
  }

  double AlphaAngle(TVector3 p, TVector3 ip) {
    TVector3 z(0,0,1);
    TVector3 p_unit = p.Unit();
    TVector3 ip_unit = ip.Unit();
    TVector3 cross1 = z.Cross(p_unit);
    TVector3 cross2 = ip_unit.Cross(p_unit);
    double angle = std::fabs(cross1.Dot(cross2)/ (cross1.Mag()*cross2.Mag()) );
    return angle;
  }

  double AlphaAngleRho(TVector3 rho, TVector3 pi) {
    TVector3 z(0,0,1);
    TVector3 rho_unit = rho.Unit();
    TVector3 pi_unit = pi.Unit();
    TVector3 cross1 = z.Cross(rho_unit);
    TVector3 cross2 = pi_unit.Cross(rho_unit);
    double angle = std::fabs(cross1.Dot(cross2)/ (cross1.Mag()*cross2.Mag()) );
    return angle;
  }


  std::vector<ic::PFCandidate*> GetTauGammaCands(ic::Tau const* tau, 
      std::map<std::size_t, ic::PFCandidate*> pfcands) {
    std::vector<ic::PFCandidate*> sig_gammas = {};
    auto const& sig_gammas_id = tau->sig_gamma_cands();
    for (auto id : sig_gammas_id) sig_gammas.push_back(pfcands[id]);
    std::sort(sig_gammas.begin(), sig_gammas.end(), 
        bind(&PFCandidate::pt, _1) > bind(&PFCandidate::pt, _2));
    return sig_gammas;
  }

  std::vector<ic::PFCandidate*> GetTauChargedHadrCands(ic::Tau const* tau, 
      std::map<std::size_t, ic::PFCandidate*> pfcands) {
    std::vector<ic::PFCandidate*> sig_charged = {};
    auto const& sig_charged_id = tau->sig_charged_cands();
    for (auto id : sig_charged_id) sig_charged.push_back(pfcands[id]);
    std::sort(sig_charged.begin(), sig_charged.end(), 
        bind(&PFCandidate::pt, _1) > bind(&PFCandidate::pt, _2));
    return sig_charged;
  }

  std::vector<ic::PFCandidate*> GetTauIsoGammaCands(ic::Tau const* tau, 
      std::map<std::size_t, ic::PFCandidate*> pfcands) {
    std::vector<ic::PFCandidate*> iso_gammas = {};
    auto const& iso_gammas_id = tau->iso_gamma_cands();
    for (auto id : iso_gammas_id) iso_gammas.push_back(pfcands[id]);
    std::sort(iso_gammas.begin(), iso_gammas.end(), 
        bind(&PFCandidate::pt, _1) > bind(&PFCandidate::pt, _2));
    return iso_gammas;
  }

  ic::Candidate* GetPiFromCands(ic::Tau const* tau, 
          std::map<std::size_t, ic::PFCandidate*> pfcands) {
    ic::Candidate* pi = new ic::Candidate();
    double pt = 0.; 
    double E = 0.; 
    double eta = 0.; 
    double phi = 0.; 
    std::vector<ic::PFCandidate*> hadrons = GetTauChargedHadrCands(tau, pfcands);
    for (auto hadron : hadrons){
      pt  += hadron->pt();
      E   += hadron->energy();
      eta += hadron->vector().Eta();
      phi += hadron->vector().Phi();
    }
    ROOT::Math::PtEtaPhiEVector vector(pt, eta, phi, E);
    pi->set_vector(vector);
    return pi;
  }

  ic::Candidate* GetPi0FromCands(ic::Tau const* tau, 
          std::map<std::size_t, ic::PFCandidate*> pfcands) {
    ic::Candidate* pi0 = new ic::Candidate();
    double sum_E = 0.; 
    double sum_EEta = 0.; 
    double sum_EPhi = 0.; 
    std::vector<ic::PFCandidate*> gammas = GetTauGammaCands(tau, pfcands);
    for (auto gamma : gammas){
      sum_E += gamma->energy();
      sum_EEta += gamma->energy()*gamma->vector().Eta();
      sum_EPhi += gamma->energy()*gamma->vector().Phi();
    }
    double pt = 0.;
    double eta = 0.;
    double phi = 0.;
    if (sum_E != 0) {
      eta = sum_EEta/sum_E;
      phi = sum_EPhi/sum_E;
      double p = sqrt(sum_E*sum_E-0.1349*0.1349);
      pt = p/cosh(eta); 
    }
    ROOT::Math::PtEtaPhiEVector vector(pt, eta, phi, sum_E);
    pi0->set_vector(vector);
    return pi0;
  }

  double boundPhi(double dphi) {
    if (dphi >= M_PI)
      dphi -= 2*M_PI;
    else if (dphi < -M_PI)
      dphi += 2*M_PI;

    return dphi;
  }

  ROOT::Math::PtEtaPhiEVector reconstructWboson(Candidate const*  lepton, Candidate const* met){


    const double metX = met->vector().px();
    const double metY = met->vector().py();
    ROOT::Math::PxPyPzEVector nuP4;
    nuP4.SetPx(metX);
    nuP4.SetPy(metY);
    double massW = 80.4;
    double bestW = -1000;
    double lowestMassDifference = -1;
    double step = 1;
    for (int i = -500; i<500;++i){

      double currentW = step*i;
      nuP4.SetPz(currentW);
      nuP4.SetE(std::sqrt(metX*metX+metY*metY+currentW*currentW ) );
  
      double M = (nuP4+lepton->vector()).M();
      double diff = std::abs(massW-M);
      //      std::cout << i << " " << diff << std::endl;
      if ( diff < lowestMassDifference  || lowestMassDifference <0 ) {
          lowestMassDifference = diff;
          bestW = currentW;
      }

    }

    if (bestW > -700) {
      nuP4.SetPz(bestW);
      nuP4.SetE(std::sqrt(metX*metX+metY*metY+bestW*bestW ) );
      //std::cout << " NUz=" << bestW << " M_W=" << (nuP4+lepton->vector()).M() <<std::endl;
    } else {
      std::cout << " Cannot reconstruct nu!\n";
      nuP4.SetPz(0);
      nuP4.SetE( std::sqrt(metX*metX+metY*metY) );
      //return;
    }

    ROOT::Math::PtEtaPhiEVector lWreco;
    lWreco = lepton->vector() + nuP4;
    return lWreco;
    
  }

  /*
  

  //---------------------------------------------------------------------------
  //Collection Generators
  //---------------------------------------------------------------------------
  std::vector<GenParticle*> MakeFinalBQuarksCollection(
      std::vector<GenParticle*> const& partVec){
    std::vector<GenParticle*> bquarks;
    std::vector<GenParticle*> bquarksFinal;
    std::vector<int> bquarkIdx;
    BOOST_FOREACH (GenParticle *particle, partVec) {
      if (unsigned(abs(particle->genVars().pdgId)) == 5) {
        bquarks.push_back(particle);
        bquarkIdx.push_back(particle->genVars().index);
      }
    }
    BOOST_FOREACH (GenParticle* bquark, bquarks) {
      bool has_bquark_daughter = false;
      std::vector<unsigned> daughterIdx = bquark->getChildren();
      BOOST_FOREACH (int Idx, daughterIdx) {
        if (std::count(bquarkIdx.begin(),bquarkIdx.end(), Idx) > 0){
          has_bquark_daughter = true;
          break;//No need to keep looping
        }//If this daughter index is in the list of bquark indices
      }//Loop through daughter indices, see if exists in bquarkIdx
      if (!has_bquark_daughter && (bquark->genVars().status != 3)) bquarksFinal.push_back(bquark);
    }//Loop through bquarks
    return bquarksFinal;
  }
  
  std::vector<GenParticle*> MakeInitialBQuarksCollection(
      std::vector<GenParticle*> const& partVec){
    std::vector<GenParticle*> bquarks;
    bool foundP5 = false;
    bool foundM5 = false;
    BOOST_FOREACH (GenParticle *particle, partVec) {
      if (foundP5 && foundM5) break;
      if ((particle->genVars().pdgId) == 5 && !foundP5
          && particle->genVars().status == 2) {
        bquarks.push_back(particle);
        foundP5 = true;
      }
      if ((particle->genVars().pdgId) == -5 && !foundM5
          && particle->genVars().status == 2) {
        bquarks.push_back(particle);
        foundM5 = true;
      }
    }
    return bquarks;
  }
  
  std::vector<GenParticle*> MakeFinalBHadronsCollection(
    std::vector<GenParticle*> const& partVec) {
    std::vector<GenParticle*> bhadrons;
    std::vector<GenParticle*> bhadronsFinal;
    std::vector<int> bhadronIdx;
    BOOST_FOREACH (GenParticle *particle, partVec) {
      if (IsBflavPdgId(particle->genVars().pdgId)) {
        bhadrons.push_back(particle);
        bhadronIdx.push_back(particle->genVars().index);
      }
    }
    BOOST_FOREACH (GenParticle* bhadron, bhadrons) {
      bool has_bhadron_daughter = false;
      std::vector<unsigned> daughterIdx = bhadron->getChildren();
      BOOST_FOREACH (int Idx, daughterIdx) {
        if (std::count(bhadronIdx.begin(),bhadronIdx.end(), Idx) > 0){
          has_bhadron_daughter = true;
          break;//No need to keep looping
        }//If this daughter index is in the list of bhadron indices
      }//Loop through daughter indices, see if exists in bhadronIdx
      if (!has_bhadron_daughter && bhadron->Pt() > 3) bhadronsFinal.push_back(bhadron);
    }//Loop through bhadrons
    return bhadronsFinal;
  }
  
  std::vector<GenParticle*> MakeDaughterCollection(
      std::vector<GenParticle*> const& partVec, unsigned const& pdgid) {
    std::vector<unsigned> dIdx;
    std::vector<GenParticle*> daughters;
    int maxDaughterIdx = 0;
    bool foundMother = false;
    BOOST_FOREACH (GenParticle *particle, partVec) {
      if (foundMother && particle->genVars().index > maxDaughterIdx) break;
      if (foundMother && std::count(dIdx.begin(), dIdx.end(), particle->genVars().index)) {
        daughters.push_back(particle);
      }
      if (!foundMother && unsigned(abs(particle->genVars().pdgId)) == pdgid) {
        dIdx = particle->getChildren();
        if (dIdx.size() == 0) {
          break;
        } else {
          foundMother = true;
          maxDaughterIdx = *std::max_element(dIdx.begin(),dIdx.end());
        }
      }
    }//Iterate through particles
    return daughters;
  }

  std::vector<GenParticle*> MakeDaughterCollectionMCATNLO(
      std::vector<GenParticle*> const& partVec, unsigned const& pdgid) {
    std::vector<GenParticle*> daughters;
    unsigned nFound = 0;
    BOOST_FOREACH (GenParticle *particle, partVec) {
      if (nFound == 2) break;
      if (unsigned(abs(particle->genVars().pdgId)) == pdgid &&
          particle->genVars().status == 3) {
        daughters.push_back(particle);
        ++nFound;
      }
    }//Iterate through particles
    return daughters;
  }

  std::vector<GenParticle*> GetStatus1Muons(
      std::vector<GenParticle*> const& partVec) {
    std::vector<GenParticle*> muons;
    BOOST_FOREACH (GenParticle *particle, partVec) {
      if (particle->genVars().status == 1 && abs(particle->genVars().pdgId) == 13) {
        muons.push_back(particle);
      }
    }//Iterate through particles
    std::vector<GenPGenP> muonPairs = MakePairs(muons);
    std::sort(muonPairs.begin(),muonPairs.end(), 
        boost::bind(MassDiffCompare<GenParticle*,GenParticle*>,_1, _2, 91.19));
    if (muonPairs.size() == 0) muons.clear();
    if (muonPairs.size() > 0) {
      muons.clear();
      muons.push_back(muonPairs.at(0).first);
      muons.push_back(muonPairs.at(0).second);
    }
    return muons;
  }

  bool IsBflavPdgId(int pdgid) {
    unsigned pdgidNoSign = unsigned(abs(pdgid));
    std::string asString = boost::lexical_cast<std::string>(pdgidNoSign);
    if (asString.size() > 1 && asString.at(0) == char('5')) return true;
    return false;
  }

  double ScaleReso(double pT, double eta, double jerUnc){
    unsigned int etaBin = 0;
    if (fabs(eta)> 1.1 && fabs(eta)<= 1.7) etaBin = 1;
    else if (fabs(eta)> 1.7 && fabs(eta)<= 2.3) etaBin = 2;
    else etaBin = 3;

    double N[4] = {-1.18,-1.50,-1.52,-0.34};
    //double Nerr[4] = {0.02,0.10,0.15,0.45};

    double S[4] = {0.41,0.63,0.77,0.57};
    //double Serr[4] = {0.00,0.01,0.03,0.03};
    double M[4] = {0.3521,0.2184,0.0229,0.1447};
    //double Merr[4] = {0.0021,0.0080,0.0181,0.0216};
    double C = 0;
    int signN = static_cast<int>(fabs(N[etaBin])/N[etaBin]);

    double sigPtoverPt_sq = (signN*pow(N[etaBin]/pT,2)+
        pow(S[etaBin],2)*pow(pT,M[etaBin]-1)+
        pow(C,2));
    double sigPtoverPt = sqrt(fabs(sigPtoverPt_sq));

    //std::cout << " -- pT = " << pT << ", sigPtoverPt = " << sigPtoverPt
    //  << " (" <<  sigPtoverPt*pT << ")"
    //  << std::endl;

    double fData = 0.10;

    if (fabs(eta) > 1.5) fData = 0.20;

    return sigPtoverPt*(sqrt(pow(fData,2)+pow(jerUnc,2)))*pT;
  }
  */
  
   std::set<int16_t> GetTriggerTypes(TriggerObject* obj){
    
     std::set<int16_t> types;
     ui64 packed_type;
     packed_type.one = obj->id();
     
     for(unsigned i=0; i<4; ++i){
       int16_t type = packed_type.four[i];
       if(type!=0) types.insert(type);
     }
     return types;
  }

  double quantile_mapping(double value, TH1D *input_cdf, TH1D *output_cdf) {
    double new_value=value;
    double p = input_cdf->Interpolate(value);
    double ylow = output_cdf->GetBinContent(1), yhigh = output_cdf->GetBinContent(2), xlow=output_cdf->GetBinLowEdge(1), xhigh=output_cdf->GetBinLowEdge(2);
    //double ydiff_hi=fabs(yhigh-p);
    double ydiff_lo=fabs(ylow-p);
    for (unsigned i=1; i<=(unsigned)output_cdf->GetNbinsX(); ++i) {
      double y = output_cdf->GetBinContent(i);
      if ((fabs(y-p)<ydiff_lo && y<p)){
         ydiff_lo = fabs(y-p);
         ylow = y;
         xlow = output_cdf->GetBinLowEdge(i);
      //}
      //if ((fabs(y-p)<ydiff_hi && y>p)){
         //ydiff_hi = fabs(y-p);
         //yhigh = y;
         xhigh = output_cdf->GetBinLowEdge(i+1);
         yhigh = output_cdf->GetBinContent(i+1);
      }
    }
    double x = (xhigh-xlow)/(yhigh-ylow)*(p-ylow);
    new_value=xlow+x;
    //std::cout << value << "    " << p << "    " << xlow << "    " << xhigh << "    " << ylow << "    " << yhigh << "    " << x << "    " << new_value << std::endl;
    return new_value;
  }

  std::pair<TVector3,double> IPAndSignificance(ic::Tau const *tau, ic::Vertex *vertex, std::vector<ic::PFCandidate*> pfcands) {
      std::vector<float> h_param = {};
      //ic::Tau  *tau = const_cast<Tau*>(intau);
      for(auto i :  tau->track_params()) h_param.push_back(i);
      double B = tau->bfield();
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > ref(tau->svx(),tau->svy(),tau->svz());
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > pv(vertex->vx(),vertex->vy(),vertex->vz());
      std::vector<ic::PFCandidate*> charged_cands = GetTauHads(tau, pfcands, 0);

      ImpactParameter IP;
      TVector3 ip = IP.CalculatePCA(B, h_param, ref, pv);

      ROOT::Math::SMatrix<double,5,5, ROOT::Math::MatRepSym<double,5>> helixCov = tau->track_params_covariance();
      SMatrixSym3D SigmaPrV = vertex->covariance();

      ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >> ip_cov = IP.CalculatePCACovariance(helixCov, SigmaPrV);

      double mag = ip.Mag();
      ROOT::Math::SVector<double, 3> ip_svec;
      ip_svec(0) = ip.X();
      ip_svec(1) = ip.Y();
      ip_svec(2) = ip.Z();

      ip_svec = ip_svec.Unit();

      double uncert = sqrt(ROOT::Math::Dot( ip_svec, ip_cov * ip_svec));
      double sig = mag/uncert;
      return std::make_pair(ip, sig);
  }

  std::pair<TVector3,ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> IPAndCov(ic::Tau const *tau, ic::Vertex *vertex, std::vector<ic::PFCandidate*> pfcands) {
      std::vector<float> h_param = {};
      //ic::Tau  *tau = const_cast<Tau*>(intau);
      for(auto i :  tau->track_params()) h_param.push_back(i);
      double B = tau->bfield();
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > ref(tau->svx(),tau->svy(),tau->svz());
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > pv(vertex->vx(),vertex->vy(),vertex->vz());
      std::vector<ic::PFCandidate*> charged_cands = GetTauHads(tau, pfcands, 0);

      ImpactParameter IP;
      TVector3 ip = IP.CalculatePCA(B, h_param, ref, pv);

      ROOT::Math::SMatrix<double,5,5, ROOT::Math::MatRepSym<double,5>> helixCov = tau->track_params_covariance();
      SMatrixSym3D SigmaPrV = vertex->covariance();

      ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >> ip_cov = IP.CalculatePCACovariance(helixCov, SigmaPrV);

      return std::make_pair(ip, ip_cov);
  }



  std::pair<TVector3,double> IPAndSignificanceMuon(ic::Muon const *muon, ic::Vertex *vertex) {
      std::vector<float> h_param = {};
      ////ic::Muon  *muon = const_cast<Muon*>(inmuon);
      for(auto i :  muon->track_params()) h_param.push_back(i);
      double B = muon->bfield();
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > ref(muon->vx(),muon->vy(),muon->vz());
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > pv(vertex->vx(),vertex->vy(),vertex->vz());

      ImpactParameter IP;
      TVector3 ip = IP.CalculatePCA(B, h_param, ref, pv);
      ROOT::Math::SMatrix<double,5,5, ROOT::Math::MatRepSym<double,5>> helixCov = muon->track_params_covariance();
      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SigmaPrV = vertex->covariance();

      ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >> ip_cov = IP.CalculatePCACovariance(helixCov, SigmaPrV);

      double mag = ip.Mag();
      ROOT::Math::SVector<double, 3> ip_svec;
      ip_svec(0) = ip.X();
      ip_svec(1) = ip.Y();
      ip_svec(2) = ip.Z();

      ip_svec = ip_svec.Unit();

      double uncert = sqrt(ROOT::Math::Dot( ip_svec, ip_cov * ip_svec));
      double sig = mag/uncert;
      return std::make_pair(ip, sig);
  }

  std::pair<TVector3,double> IPAndSignificanceElectron(ic::Electron const *elec, ic::Vertex *vertex) {
      std::vector<float> h_param = {};
      for(auto i :  elec->track_params()) h_param.push_back(i);
      double B = elec->bfield();
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > ref(elec->vx(),elec->vy(),elec->vz());
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > pv(vertex->vx(),vertex->vy(),vertex->vz());

      ImpactParameter IP;
      TVector3 ip = IP.CalculatePCA(B, h_param, ref, pv);
      ROOT::Math::SMatrix<double,5,5, ROOT::Math::MatRepSym<double,5>> helixCov = elec->track_params_covariance();
      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SigmaPrV = vertex->covariance();

      ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >> ip_cov = IP.CalculatePCACovariance(helixCov, SigmaPrV);

      double mag = ip.Mag();
      ROOT::Math::SVector<double, 3> ip_svec;
      ip_svec(0) = ip.X();
      ip_svec(1) = ip.Y();
      ip_svec(2) = ip.Z();

      ip_svec = ip_svec.Unit();

      double uncert = sqrt(ROOT::Math::Dot( ip_svec, ip_cov * ip_svec));
      double sig = mag/uncert;
      return std::make_pair(ip, sig);
  }

  std::pair<TVector3, ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> IPAndCovMuon(ic::Muon const *muon, ic::Vertex *vertex) {
      std::vector<float> h_param = {};
      ////ic::Muon  *muon = const_cast<Muon*>(inmuon);
      for(auto i :  muon->track_params()) h_param.push_back(i);
      double B = muon->bfield();
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > ref(muon->vx(),muon->vy(),muon->vz());
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > pv(vertex->vx(),vertex->vy(),vertex->vz());

      ImpactParameter IP;
      TVector3 ip = IP.CalculatePCA(B, h_param, ref, pv);
      ROOT::Math::SMatrix<double,5,5, ROOT::Math::MatRepSym<double,5>> helixCov = muon->track_params_covariance();
      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SigmaPrV = vertex->covariance();

      ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >> ip_cov = IP.CalculatePCACovariance(helixCov, SigmaPrV);

      return std::make_pair(ip, ip_cov);
  }

  std::pair<TVector3, ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >>> IPAndCovElec(ic::Electron const *elec, ic::Vertex *vertex) {
      std::vector<float> h_param = {};
      for(auto i :  elec->track_params()) h_param.push_back(i);
      double B = elec->bfield();
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > ref(elec->vx(),elec->vy(),elec->vz());
      ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float> > pv(vertex->vx(),vertex->vy(),vertex->vz());

      ImpactParameter IP;
      TVector3 ip = IP.CalculatePCA(B, h_param, ref, pv);
      ROOT::Math::SMatrix<double,5,5, ROOT::Math::MatRepSym<double,5>> helixCov = elec->track_params_covariance();
      ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SigmaPrV = vertex->covariance();

      ROOT::Math::SMatrix<double,3,3, ROOT::Math::MatRepStd< double, 3, 3 >> ip_cov = IP.CalculatePCACovariance(helixCov, SigmaPrV);

      return std::make_pair(ip, ip_cov);
  }

} //namespace
