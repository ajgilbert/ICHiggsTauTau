#include "UserCode/ICHiggsTauTau/Analysis/Utilities/interface/SCalculator.h"
#include <iostream>

SCalculator::SCalculator(){
}

SCalculator::SCalculator(string type):
  type_(type)
{
}

void 
SCalculator::Configure(vector<TLorentzVector> TauAndProd, TLorentzVector Frame, int charge){
  for(unsigned int i=0; i<TauAndProd.size(); i++){
    TauAndProd_HRF.push_back(Boost(TauAndProd.at(i), Frame));

  }
  charge_=charge;
}

bool
SCalculator::isConfigured(){
  if(TauAndProd_LF.size()!=2){ std::cout<<"Error:   SCalculator is not Configured! Check  the size of input vector!  Size =  "<< TauAndProd_LF.size() <<std::endl; return false;} return true;
}

SCalculator::~SCalculator(){
}

TLorentzVector 
SCalculator::Boost(TLorentzVector pB, TLorentzVector frame){
  TMatrixT<double> transform(4,4);
  TMatrixT<double> result(4,1);
  TVectorT<double> vec(4); 
  TVector3 b;
  
  if(frame.Vect().Mag()==0){ std::cout<<"RH Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
  if(frame.E()==0){ std::cout<<" Caution: Please check that you perform boost correctly!  " <<std::endl; return pB;} 
  else   b=frame.Vect()*(1/frame.E());
  
  
  vec(0)  = pB.E();    vec(1)  = pB.Px();
  vec(2)  = pB.Py();  vec(3)  = pB.Pz();

  double gamma  = 1/sqrt( 1 - b.Mag2());
 
  transform(0,0)=gamma; transform(0,1) =- gamma*b.X() ;  transform(0,2) =  - gamma*b.Y();  transform(0,3) = - gamma*b.Z(); 
 
  transform(1,0)=-gamma*b.X(); transform(1,1) =(1+ (gamma-1)*b.X()*b.X()/b.Mag2()) ;  transform(1,2) = ((gamma-1)*b.X()*b.Y()/b.Mag2());  transform(1,3) = ((gamma-1)*b.X()*b.Z()/b.Mag2());
  transform(2,0)=-gamma*b.Y(); transform(2,1) = ((gamma-1)*b.Y()*b.X()/b.Mag2());  transform(2,2) = (1 + (gamma-1)*b.Y()*b.Y()/b.Mag2());  transform(2,3) =  ((gamma-1)*b.Y()*b.Z()/b.Mag2()); 
  transform(3,0)=-gamma*b.Z(); transform(3,1) =((gamma-1)*b.Z()*b.X()/b.Mag2()) ;  transform(3,2) = ((gamma-1)*b.Z()*b.Y()/b.Mag2());  transform(3,3) = (1 + (gamma-1)*b.Z()*b.Z()/b.Mag2()); 
  result=transform*convertToMatrix(vec);
  return TLorentzVector(result(1,0), result(2,0) ,result(3,0), result(0,0));
}
TMatrixT<double> SCalculator::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}

TVector3
SCalculator::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}

TVector3
SCalculator::pv(){
  TVector3 out(0,0,0);
  if(type_=="pion") out = TauAndProd_HRF.at(1).Vect();
  if(type_=="rho"){
    TLorentzVector pi  = TauAndProd_HRF.at(1);
    TLorentzVector pi0 = TauAndProd_HRF.at(2);
    TLorentzVector Tau = TauAndProd_HRF.at(0);
    TLorentzVector q= pi  - pi0;
    TLorentzVector P= Tau;
    TLorentzVector N= Tau - pi - pi0;
    out = P.M()*(2*(q*N)*q.Vect() - q.Mag2()*N.Vect()) * (1/ (2*(q*N)*(q*P) - q.Mag2()*(N*P)));
  }
  if(type_=="a1"){
    PolarimetricA1  a1pol;
    a1pol.Configure(TauAndProd_HRF,charge_);
    out = -a1pol.PVC().Vect();
  }
  return out;
}

void SCalculator::SortPions(std::vector<TLorentzVector>& pionsvec, std::vector<double>& charges)
{
  
  int npim(0),npip(0), npin(0);
  int OSMCPionIndex(-99);
  int SSMCPion1Index(-99);
  int SSMCPion2Index(-99);
  int OSCharge(-99);
  int SS1Charge(-99);
  int SS2Charge(-99);
 
  TLorentzVector os;
  TLorentzVector ss1;
  TLorentzVector ss2;

  int MCNeutralPionIndex(-99);
  int MCChargedPionIndex(-99);
  int NeutralPionCharge(-99);
  int ChargedPionCharge(-99);

  TLorentzVector ChargedPion;
  TLorentzVector NeutralPion;
  //cout<<"test 0"<<endl;
  if(charges.size()==3) //A1
    {
      //cout<<"test 1"<<endl;
      for(unsigned int i=0; i<charges.size(); i++){
	if( charges.at(i)== 1) npip++;
	if( charges.at(i)==-1) npim++;
      }
      if(npip == 1 && npim == 2){
	//cout<<"test 2"<<endl;
	int nss=0;
	for(unsigned int i=0; i<charges.size(); i++){
	  if(charges.at(i)== 1){
	    OSCharge=1;
	    OSMCPionIndex=i;
	  }
	  if(charges.at(i)== -1 && nss == 0){
	    nss++;
	    SS1Charge=-1;
	    SSMCPion1Index=i;
	  }
	  if(charges.at(i)== -1 && nss == 1){
	    SS2Charge=-1;
	    SSMCPion2Index=i;
	  }
	}
      } 
      else if( npip== 2 && npim==1){
	//cout<<"test 3"<<endl;
	int nss=0;
	for(unsigned int i=0; i<charges.size(); i++){
	  if(charges.at(i) == -1){
	    OSCharge=-1;
	    OSMCPionIndex=i;
	  }
	  if(charges.at(i) == 1 && nss ==0){
	    nss++;
	    SS1Charge=1;
	    SSMCPion1Index=i;
	  }
	  if(charges.at(i) == 1 && nss == 1){
	    SS2Charge=1;
	    SSMCPion2Index=i;
	  }
	}
      }
      if(OSMCPionIndex!=-99 && SSMCPion1Index!=-99 && SSMCPion2Index!=-99){
	os=pionsvec.at(OSMCPionIndex);
	ss1=pionsvec.at(SSMCPion1Index);
	ss2=pionsvec.at(SSMCPion2Index);
      
	charges.clear();
	charges.push_back(OSCharge);
	charges.push_back(SS1Charge);
	charges.push_back(SS2Charge);
      
	pionsvec.clear();
	pionsvec.push_back(os);
	pionsvec.push_back(ss1);
	pionsvec.push_back(ss2);
      }
      else{
	charges.clear();
	charges.push_back(-99);
	charges.push_back(-99);
	charges.push_back(-99);
      }
    }
  
  if(charges.size()==2) //Rho
    {
      for(unsigned int i=0; i<charges.size(); i++){
	if( charges.at(i)== 1) npip++;
	if( charges.at(i)==-1) npim++;
	if( charges.at(i)==0) npin++;
      }
      if(npip == 1 && npin == 1){
	for(unsigned int i=0; i<charges.size(); i++){
	  if(charges.at(i)== 1){
	    ChargedPionCharge=1;
	    MCChargedPionIndex=i;
	  }
	  if(charges.at(i)== 0){
	    NeutralPionCharge=0;
	    MCNeutralPionIndex=i;
	  }
	}
      }
      if( npim== 1 && npin==1){
	for(unsigned int i=0; i<charges.size(); i++){
	  if(charges.at(i) == -1){
	    ChargedPionCharge=-1;
	    MCChargedPionIndex=i;
	  }
	  if(charges.at(i) == 0){
	    NeutralPionCharge=0;
	    MCNeutralPionIndex=i;
	  }
	}
      }
      
      ChargedPion=pionsvec.at(MCChargedPionIndex);
      NeutralPion=pionsvec.at(MCNeutralPionIndex);
      
      charges.clear();
      charges.push_back(ChargedPionCharge);
      charges.push_back(NeutralPionCharge);
      
      pionsvec.clear();
      pionsvec.push_back(ChargedPion);
      pionsvec.push_back(NeutralPion);
    }
  //return pionsvec;
}

bool SCalculator::isOk(TString type1, TString type2, TLorentzVector tauMinus, std::vector<TLorentzVector> sumPionsMinus, std::vector<double> sumPionsChargeMinus, TLorentzVector tauPlus, std::vector<TLorentzVector> sumPionsPlus, std::vector<double> sumPionsChargePlus)
{
  SCalculator Scalc1(type1.Data());
  SCalculator Scalc2(type2.Data());
  
  TLorentzVector zeroLV(0,0,0,0);
  TLorentzVector HadLVMinus(0,0,0,0);
  TLorentzVector HadLVPlus(0,0,0,0);

  if(type1!="pion") Scalc1.SortPions(sumPionsMinus, sumPionsChargeMinus);
  if(type2!="pion") Scalc2.SortPions(sumPionsPlus, sumPionsChargePlus);

  bool isgoodcharge=true;
  if(sumPionsChargeMinus.at(0)==-99 || sumPionsChargeMinus.at(1)==-99 || sumPionsChargeMinus.at(2)==-99 || sumPionsChargePlus.at(0)==-99 || sumPionsChargePlus.at(1)==-99 || sumPionsChargePlus.at(2)==-99)isgoodcharge=false;

  vector<TLorentzVector> tauandprodminus;
  vector<TLorentzVector> tauandprodplus;
  bool pionszero=false;
  tauandprodminus.push_back(tauMinus);
  for(unsigned int i=0; i<sumPionsMinus.size();i++) {tauandprodminus.push_back(sumPionsMinus.at(i))/*;sumPionsMinus.at(i).Print()*/;if (sumPionsMinus.at(i)==zeroLV)pionszero=true;}
  
  tauandprodplus.push_back(tauPlus); 
  for(unsigned int i=0; i<sumPionsPlus.size();i++) {tauandprodplus.push_back(sumPionsPlus.at(i))/*;sumPionsPlus.at(i).Print()*/;if (sumPionsPlus.at(i)==zeroLV)pionszero=true;}   

  Scalc1.Configure(tauandprodminus,tauandprodminus.at(0)+tauandprodplus.at(0), -1);

  TVector3 h1=Scalc1.pv();
  
  Scalc2.Configure(tauandprodplus,tauandprodminus.at(0)+tauandprodplus.at(0), +1);
  TVector3 h2=Scalc2.pv();
  
  if(std::isnan(h1.Mag())==true || std::isnan(h2.Mag())==true || tauMinus==zeroLV || tauPlus==zeroLV || tauMinus==tauPlus || pionszero ||sumPionsPlus==sumPionsMinus || !isgoodcharge) return false;
  else return true;
}
double SCalculator::AcopAngle(TString type1, TString type2, TLorentzVector tauMinus, std::vector<TLorentzVector> sumPionsMinus, std::vector<double> sumPionsChargeMinus, TLorentzVector tauPlus, std::vector<TLorentzVector> sumPionsPlus, std::vector<double> sumPionsChargePlus)
{  
  
  SCalculator Scalc1(type1.Data());
  SCalculator Scalc2(type2.Data());

  TLorentzVector zeroLV(0,0,0,0);
  TLorentzVector HadLVMinus(0,0,0,0);
  TLorentzVector HadLVPlus(0,0,0,0);

  if(type1!="pion") Scalc1.SortPions(sumPionsMinus, sumPionsChargeMinus);
  if(type2!="pion") Scalc2.SortPions(sumPionsPlus, sumPionsChargePlus);
  
  //if(type1=="pion" && type2!="a1")
  vector<TLorentzVector> tauandprodminus;
  vector<TLorentzVector> tauandprodplus;
  
  tauandprodminus.push_back(tauMinus);
  for(unsigned int i=0; i<sumPionsMinus.size();i++) {tauandprodminus.push_back(sumPionsMinus.at(i));}
  
  
  tauandprodplus.push_back(tauPlus); 
  for(unsigned int i=0; i<sumPionsPlus.size();i++) {tauandprodplus.push_back(sumPionsPlus.at(i));}  
 
  Scalc1.Configure(tauandprodminus,tauandprodminus.at(0)+tauandprodplus.at(0), -1);
  TVector3 h1=Scalc1.pv();

  Scalc2.Configure(tauandprodplus,tauandprodminus.at(0)+tauandprodplus.at(0), +1);
  TVector3 h2=Scalc2.pv();

  if(tauandprodminus.at(0)==zeroLV){ cout<<endl;tauandprodminus.at(0).Print();}
  if(tauandprodplus.at(0)==zeroLV) {cout<<endl;tauandprodplus.at(0).Print();}
  
  TLorentzVector tauminus_HRF = Scalc1.Boost(tauandprodminus.at(0),tauandprodminus.at(0)+tauandprodplus.at(0));

  TLorentzVector tauplus_HRF  = Scalc2.Boost(tauandprodplus.at(0),tauandprodminus.at(0)+tauandprodplus.at(0));
  
  // tauandprodminus.at(0).Print();
  // tauandprodplus.at(0).Print();
  // h1.Print();
  // h2.Print();
  
  double h1Norm=1./h1.Mag();
  double h2Norm=1./h2.Mag();
  h1=h1*h1Norm;
  h2=h2*h2Norm;
  double k1Norm=1./((h1.Cross(tauminus_HRF.Vect().Unit())).Mag());
  double k2Norm=1./((h2.Cross(tauplus_HRF.Vect().Unit())).Mag());
  TVector3 k1 = (h1.Cross(tauminus_HRF.Vect().Unit()))*k1Norm;
  TVector3 k2 = (h2.Cross(tauplus_HRF.Vect().Unit()))*k2Norm;
  
  for (unsigned int i=0;i<sumPionsMinus.size();i++){HadLVMinus+=sumPionsMinus.at(i);}
  for (unsigned int i=0;i<sumPionsPlus.size();i++){HadLVPlus+=sumPionsPlus.at(i);}
 
  // tauminus_HRF.Vect().Unit().Print();
  // tauplus_HRF.Vect().Unit().Print();
  // cout<<"k1norm: "<<k1Norm<<endl;
  // cout<<"k2norm: "<<k2Norm<<endl;
  // k1.Print();
  // k2.Print();
  // cout<<"atan2: "<<TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)<<endl;
  //if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0){cout<<"acos: "<<acos(k1*k2)<<endl; return acos(k1*k2);}
  //else{cout<<"acos: "<<2.*TMath::Pi()-acos(k1*k2)<<endl; return (2.*TMath::Pi()-acos(k1*k2));}
  if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0){// if( TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)>2.1) {tauandprodminus.at(0).Print();HadLVMinus.Print();
    // tauandprodplus.at(0).Print();HadLVPlus.Print();
    // h1.Print();
    // h2.Print();tauminus_HRF.Vect().Unit().Print();
    // tauplus_HRF.Vect().Unit().Print();
    // cout<<"k1norm: "<<k1Norm<<endl;
    // cout<<"k2norm: "<<k2Norm<<endl;
    // k1.Print();
    // k2.Print();cout<<"ATan2: "<<TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)<<endl;} 
    return TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2);}
  else{// if((2.*TMath::Pi()- TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2))<4.1){tauandprodminus.at(0).Print();HadLVMinus.Print();
    // tauandprodplus.at(0).Print();HadLVPlus.Print();
    // h1.Print();
    // h2.Print();tauminus_HRF.Vect().Unit().Print();
    // tauplus_HRF.Vect().Unit().Print();
    // cout<<"k1norm: "<<k1Norm<<endl;
    // cout<<"k2norm: "<<k2Norm<<endl;
    // k1.Print();
    // k2.Print();cout<<"ATan2: "<<2.*TMath::Pi()-TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)<<endl;} 
    return (2.*TMath::Pi()-TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2));}
}

double SCalculator::M(TLorentzVector LV) 
{
  Double_t mm = LV.T()*LV.T()-LV.X()*LV.X()-LV.Y()*LV.Y()-LV.Z()*LV.Z();
  return mm < 0.0 ? -TMath::Sqrt(-mm) : TMath::Sqrt(mm);
}
