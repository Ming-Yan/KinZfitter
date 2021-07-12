
// -*- C++ -*-
//
// Package:     Subsystem/Package
// Class  :     HelperFunction
// 
// Implementation:
//     [Notes on implementation]
//
// Original Author:  Tongguang Cheng
//         Created:  Mon, 21 Dec 2015 12:47:33 GMT
//

#ifndef HelperFunction_CC
#define HelperFunction_CC
// system include files

// user include files
#include "/afs/cern.ch/work/m/milee/MYcode/KinZfitter/HelperFunction/interface/HelperFunction.h"

// fileinPath
#include "FWCore/ParameterSet/interface/FileInPath.h"


HelperFunction::HelperFunction(int period, bool isData)
{

  //declarations
  debug_ = 0;

  // MORIOND 17
  /* TString s_corr_e_1 = TString("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_1.root");
  TString s_corr_e_2 = TString("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_2.root");
  TString s_corr_e_3 = TString("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_3.root");
  TString s_corr_mu = TString("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2muLUT_m2mu.root");

  f_corr_e_1 = boost::shared_ptr<TFile>( new TFile(s_corr_e_1)); 
  f_corr_e_2 = boost::shared_ptr<TFile>( new TFile(s_corr_e_2)); 
  f_corr_e_3 = boost::shared_ptr<TFile>( new TFile(s_corr_e_3)); 
  f_corr_mu = boost::shared_ptr<TFile>( new TFile(s_corr_mu));
        
  el_corr_1 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_1->Get("2e")->Clone() )) );
  el_corr_2 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_2->Get("2e")->Clone() )) );
  el_corr_3 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_3->Get("2e")->Clone() )) );

  mu_corr = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_mu->Get("2mu")->Clone() )) );
        
  x_elpTaxis_1 = el_corr_1->GetXaxis(); y_eletaaxis_1 = el_corr_1->GetYaxis();
  maxPtEl_1 = x_elpTaxis_1->GetXmax(); minPtEl_1 = x_elpTaxis_1->GetXmin();

  x_elpTaxis_2 = el_corr_2->GetXaxis(); y_eletaaxis_2 = el_corr_2->GetYaxis();
  maxPtEl_2 = x_elpTaxis_2->GetXmax(); minPtEl_2 = x_elpTaxis_2->GetXmin();

  x_elpTaxis_3 = el_corr_3->GetXaxis(); y_eletaaxis_3 = el_corr_3->GetYaxis();
  maxPtEl_3 = x_elpTaxis_3->GetXmax(); minPtEl_3 = x_elpTaxis_3->GetXmin();

  x_mupTaxis = mu_corr->GetXaxis(); y_muetaaxis = mu_corr->GetYaxis();
  maxPtMu = x_mupTaxis->GetXmax(); minPtMu = x_mupTaxis->GetXmin(); */
  mu_mc.clear(); trke_mc.clear();ecalp03_mc.clear();ecalp07_mc.clear();
  mu.clear(); trke.clear();ecalp03.clear();ecalp07.clear();
  if (period == 16)
  {
	mu_mc.push_back(1.28876);//0-0.9
	mu_mc.push_back(1.28944);//0.9to1.8
	mu_mc.push_back(1.1305);//1.8-2.4
	mu.push_back(1.22276);
	mu.push_back(1.07634);
	mu.push_back(1.19583);
	trke_mc.push_back(2.18089);//0-1.44
	trke_mc.push_back(5.66962);//1.44-1.6
	trke_mc.push_back(3.30158);//1.6-2
	trke_mc.push_back(2.82185);//2-2.5
	trke.push_back(2.33334);
	trke.push_back(4.44167);
	trke.push_back(3.36206);
	trke.push_back(2.80757);
	ecalp03_mc.push_back(1.23272);//0-0.8_<0.03
	ecalp03_mc.push_back(1.21213);//0.8-1_<0.03
	ecalp03_mc.push_back(1.22961);//0-1_>0.03
	ecalp03.push_back(1.15062);//0-0.8_<0.03
	ecalp03.push_back(1.25599);//0.8-1_<0.03
	ecalp03.push_back(1.18747);//0-1_>0.03
	ecalp07_mc.push_back(1.11581);//1-1.2_<0.07
	ecalp07_mc.push_back(1.26125);//1.2-1.44_<0.07
	ecalp07_mc.push_back(1.18524);//1.44-1.57_<0.07
	ecalp07_mc.push_back(1.14954);//1.57-2._<0.07
	ecalp07_mc.push_back(1.19635);//2-2.5_<0.07
	ecalp07_mc.push_back(0.864778);//1-2.5_>0.07
	ecalp07.push_back(1.09127);//1-1.2_<0.07
	ecalp07.push_back(1.18978);//1.2-1.44_<0.07
	ecalp07.push_back(1.33164);//1.44-1.57_<0.07
	ecalp07.push_back(1.12591);//1.57-2._<0.07
	ecalp07.push_back(1.1569);//2-2.5_<0.07
	ecalp07.push_back(0.812196);//1-2.5_>0.07
  }
if (period == 17||period==18)
  {
	mu_mc.push_back(1.2691);//0-0.9
	mu_mc.push_back(1.29502);//0.9-1.8
	mu_mc.push_back(1.16235);//1.8-2.4
	mu.push_back(1.21283);
	mu.push_back(1.055181);
	mu.push_back(1.2006);
	trke_mc.push_back(2.40465);//0-1.44
	trke_mc.push_back(6.54636);//1.44-1.6
	trke_mc.push_back(3.1319);//1.6-2
	trke_mc.push_back(4.84664);//2-2.5
	trke.push_back(2.29313);
	trke.push_back(4.456603);
	trke.push_back(3.21738);
	trke.push_back(4.16995);
	ecalp03_mc.push_back(1.15208);//0-0.8_<0.03
	ecalp03_mc.push_back(1.24718);//0.8-1_<0.03
	ecalp03_mc.push_back(1.20714);//0-1_>0.03
	ecalp03.push_back(1.14794);//0-0.8_<0.03
	ecalp03.push_back(1.15732);//0.8-1_<0.03
	ecalp03.push_back(1.18867);//0-1_>0.03
	ecalp07_mc.push_back(1.07811);//1-1.2_<0.07
	ecalp07_mc.push_back(1.18304);//1.2-1.44_<0.07
	ecalp07_mc.push_back(1.21457);//1.44-1.57_<0.07
	ecalp07_mc.push_back(1.17568);//1.57-2._<0.07
	ecalp07_mc.push_back(1.14786);//2-2.5_<0.07
	ecalp07_mc.push_back(0.95666);//1-2.5_>0.07
	ecalp07.push_back(1.06918);//1-1.2_<0.07
	ecalp07.push_back(1.23277);//1.2-1.44_<0.07
	ecalp07.push_back(1.400);//1.44-1.57_<0.07
	ecalp07.push_back(1.12866);//1.57-2._<0.07
	ecalp07.push_back(1.15218);//2-2.5_<0.07
	ecalp07.push_back(0.931127);//1-2.5_>0.07
  }

}


HelperFunction::~HelperFunction(){}

double HelperFunction:: masserrorFullCov(std::vector<TLorentzVector> p4s, TMatrixDSym covMatrix)
{
  int ndim = 3*p4s.size();
  //  if(debug_) cout<<""<<endl;

  TMatrixD jacobian(1,ndim);

  double e = 0; double mass = 0;
  double px = 0; double py = 0; double pz = 0;
  for (unsigned int ip = 0; ip < p4s.size(); ip++) {
         
    e = e + p4s[ip].E();
    px = px + p4s[ip].Px();
    py = py + p4s[ip].Py();
    pz = pz + p4s[ip].Pz();
  }

  mass = TMath::Sqrt(e*e-px*px-py*py-pz*pz);

  for (unsigned int i = 0, o = 0; i < p4s.size(); i++, o += 3) {

    double pxi = p4s[i].Px();
    double pyi = p4s[i].Py();
    double pzi = p4s[i].Pz();
    double ei = p4s[i].E();

    jacobian(0, o+0) = (e*(pxi/ei) - px)/mass;
    jacobian(0, o+1) = (e*(pyi/ei) - py)/mass;
    jacobian(0, o+2) = (e*(pzi/ei) - pz)/mass;
  }

  TMatrixDSym massCov = covMatrix.Similarity(jacobian);

  double dm2 = massCov(0,0);
  return (dm2 > 0 ? std::sqrt(dm2) : 0.0);
}


double HelperFunction::masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr){
  
  TLorentzVector compositeParticle ;
  for(unsigned int i=0; i<Lep.size(); i++){
    compositeParticle+=Lep[i];
  }
  double mass  =  compositeParticle.M();

  double masserr = 0;

  for(unsigned int i=0; i<Lep.size(); i++){
    TLorentzVector variedLep; // = Lep[i];

    variedLep.SetPtEtaPhiM(Lep[i].Pt()+ pterr[i], Lep[i].Eta(), Lep[i].Phi(), Lep[i].M());
    TLorentzVector compositeParticleVariation ;
    for(unsigned int j=0; j<Lep.size(); j++){
      if(i!=j)compositeParticleVariation+=Lep[j];
      else compositeParticleVariation+=variedLep;
    }

    masserr += (compositeParticleVariation.M()-mass)*(compositeParticleVariation.M()-mass);
  }

  return sqrt(masserr);
}


double HelperFunction::pterr(TreeReader &data, Int_t lepID, Int_t fsrID, Int_t channel,  bool isData,bool corr){
  // reco::GsfElectron *gsf; reco::Muon *mu;
  // reco::PFCandidate *pf;
  // pat::Muon *patmu;
  double pterrLep = 0.0;
  float *elePt = data.GetPtrFloat("elePt");
  float *muPt = data.GetPtrFloat("muPt");
  float *eleEta = data.GetPtrFloat("eleEta");
  float *muEta = data.GetPtrFloat("muEta");
  Int_t *eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
  //if(fsrID!=-99)data.GetPtrFloat("pfPhoEtErr");
  
  if (channel == 0) {
  double pT_e = elePt[lepID];
  double eta_e = eleEta[lepID];
  pterrLep = correlepterr(data,lepID);
  if(corr){
  if(eleEcalDrivenSeed[lepID]==1){		//ecal driven
     if(isData)
	 {
		if(pterrLep/pT_e < 0.03&&fabs(eta_e)<0.8)pterrLep*=ecalp03[0];
		if(pterrLep/pT_e < 0.03&&fabs(eta_e)>=0.8&&fabs(eta_e)<1.)pterrLep*=ecalp03[1];
		if(pterrLep/pT_e >= 0.03&&fabs(eta_e)<1.)pterrLep*=ecalp03[2];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=1.&&fabs(eta_e)<1.2)pterrLep*=ecalp07[0];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=1.2&&fabs(eta_e)<1.44)pterrLep*=ecalp07[1];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=1.44&&fabs(eta_e)<1.57)pterrLep*=ecalp07[2];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=1.57&&fabs(eta_e)<2.0)pterrLep*=ecalp07[3];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=2.0&&fabs(eta_e)<2.5)pterrLep*=ecalp07[4];
		if(pterrLep/pT_e >= 0.07&&fabs(eta_e)>=1.0&&fabs(eta_e)<2.5)pterrLep*=ecalp07[5];
	 }
	 else{
		if(pterrLep/pT_e < 0.03&&fabs(eta_e)<0.8)pterrLep*=ecalp03_mc[0];
		if(pterrLep/pT_e < 0.03&&fabs(eta_e)>=0.8&&fabs(eta_e)<1.)pterrLep*=ecalp03_mc[1];
		if(pterrLep/pT_e >= 0.03&&fabs(eta_e)<1.)pterrLep*=ecalp03_mc[2]; 
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=1.&&fabs(eta_e)<1.2)pterrLep*=ecalp07_mc[0];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=1.2&&fabs(eta_e)<1.44)pterrLep*=ecalp07_mc[1];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=1.44&&fabs(eta_e)<1.57)pterrLep*=ecalp07_mc[2];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=1.57&&fabs(eta_e)<2.0)pterrLep*=ecalp07_mc[3];
		if(pterrLep/pT_e < 0.07&&fabs(eta_e)>=2.0&&fabs(eta_e)<2.5)pterrLep*=ecalp07_mc[4];
		if(pterrLep/pT_e >= 0.07&&fabs(eta_e)>=1.0&&fabs(eta_e)<2.5)pterrLep*=ecalp07_mc[5];
	 }
    }
    else //tracker driven
      {
	if(isData)
	{
		if(fabs(eta_e)<1.44)pterrLep*=trke[0];
		else if(fabs(eta_e)<1.6)pterrLep*=trke[1];
		else if(fabs(eta_e)<2.)pterrLep*=trke[2];
		else if(fabs(eta_e)<2.5)pterrLep*=trke[3];
	}
	else 
	{
		if(fabs(eta_e)<1.44)pterrLep*=trke_mc[0];
		else if(fabs(eta_e)<1.6)pterrLep*=trke_mc[1];
		else if(fabs(eta_e)<2.)pterrLep*=trke_mc[2];
		else if(fabs(eta_e)<2.5)pterrLep*=trke_mc[3];
		}
  }} else {
	  pterrLep*=1.0;
	}
      }
	  
  
  else 
    {
      float* muBestTrkPtError = data.GetPtrFloat("muBestTrkPtError");
      pterrLep= muBestTrkPtError[lepID];
      if(debug_)cout<<"reco pt err is "<<pterrLep<<endl;
      if(corr){
		if(isData)
		{
			if(fabs(muEta[lepID])<0.9)pterrLep*=mu[0];
			else if(fabs(muEta[lepID])<1.8)pterrLep*=mu[1];
			else if(fabs(muEta[lepID])<2.4)pterrLep*=mu[2];
		}
		else
		{
			if(fabs(muEta[lepID])<0.9)pterrLep*=mu_mc[0];
			if(fabs(muEta[lepID])<1.8&&fabs(muEta[lepID])>0.9)pterrLep*=mu_mc[1];
			if(fabs(muEta[lepID])<2.4&&fabs(muEta[lepID])>1.8)pterrLep*=mu_mc[2];	
		}

      } else {
			pterrLep*=1.0;
      }
    }
  if(fsrID!=-99)//FSR constribution
	    {
	      float *pfPhoEtErr = data.GetPtrFloat("pfPhoEtErr");
	      pterrLep = pfPhoEtErr[fsrID];
	  }
  return pterrLep;
}
double HelperFunction::correlepterr(TreeReader &data, Int_t lepID)
{
	float *ecalEnergy = data.GetPtrFloat("eleCalibEn");
	float *eleSCEta = data.GetPtrFloat("eleSCEta");
	float *elePt = data.GetPtrFloat("elePt");
	float *eleEta = data.GetPtrFloat("eleEta");
	float *elePhi = data.GetPtrFloat("elePhi");
	Int_t *eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
	float *elePtError =data.GetPtrFloat("elePtError");
	double pterr = 0.;
	if(eleEcalDrivenSeed[lepID]==1&&elePtError[lepID]!=-99)
	  {
	    pterr =elePtError[lepID];
	  }
	else {
	TLorentzVector ele;
	ele.SetPtEtaPhiM(elePt[lepID],eleEta[lepID],elePhi[lepID],0.511*0.001);
	float err2 = 0.0;
	if(fabs(eleSCEta[lepID])< 1.479)
	{
		 err2 += (5.24e-02*5.24e-02)/ecalEnergy[lepID];
		 err2 += (2.01e-01*2.01e-01)/(ecalEnergy[lepID]*ecalEnergy[lepID]);
		 err2 += 1.00e-02*1.00e-02;
	}
	else 
	{
		err2 += (1.46e-01*1.46e-01)/ecalEnergy[lepID];
		err2 += (9.21e-01*9.21e-01)/(ecalEnergy[lepID]*ecalEnergy[lepID]);
		err2 += 1.94e-03*1.94e-03;	
	}
	float perr = ecalEnergy[lepID]*sqrt(err2);
	pterr = perr*elePt[lepID]/ele.P();
	}
	return pterr;
}
#endif
