Using Z mass kinematic constraint(s) to refit lepton momenta in H-to-ZZ2l2q

Code structure
KinZfitter : the class read the inputs from jets that form the hadronic Z decay,, do the refitting, and get the refitted results HelperFunction : the class that read jet pT errors by accessing the JME::JetResolution in CMSSW

To include the refit in your analyzer:

0.Check out package

cd $CMSSW_BASE/src

git clone https://github.com/tocheng/KinZfitter.git

cd KinZfitter

git checkout -b from-Zhadv1.0 Zhadv1.0

cd ../

scram b -j 8

0.Add the package and JER related package into your BuildFile.xml

In your main analyzer analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

1.include the head file

#include "KinZfitter/KinZfitter/interface/KinZfitter.h"

2.Declare and then initialize the KinZfitter class when initializing your analyzer i.e., declare and initialize in beginJob() in EDAnalyzer framework KinZfitter *kinZfitter; kinZfitter = new KinZfitter(isData); //(In data, (isData=true). In mc (isData=false))

3.Prepare inputs after Higgs decay to ZZ2l2q candidate is formed:

leptons: prepare the TLorentzVectors lep+fsr for leptontic Z decay (similar to Z1/Z2 lepton TLorentzVector for HZZ4L matrix element calculation) say TLorentzVector lep1, lep2;

jets: prepare the TLorentzVectors (resolved) jets from hadronic Z decay
say TLorentzVector jet1, jet2;

jet resolution: prepare the jet resolution objects for pt and phi resolution say JME::JetResolution resolution_pt, resolution_phi;

resolution_pt = JME::JetResolution::get(iSetup, "AK4PFchs_pt"); resolution_phi = JME::JetResolution::get(iSetup, "AK4PFchs_phi");

Rho : get rho parameter for isolation calculation (double precision), say rho

4.Setup, refit and get the refitted results:

In your analyzer, do

  kinZfitter->Setup2L2Q(Zlep,Zhad,resolution_pt,resolution_phi,rho);
  kinZfitter->KinRefitZlepZhad();

  // To get refit mZZ
  double massZZREFIT = kinZfitter->GetRefitMZZ2L2Q();
  // To get refit hadronic mZ (mjj)
  double mass4lErrREFIT = kinZfitter->GetRefitMZhad();
To access JER in local sqlite file in cms config file, add the following lines to you config file

 process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
 dBJERFile = os.environ.get('CMSSW_BASE')+"/src/KinZfitter/HelperFunction/hists/Summer15_25nsV6_MC_JER.db"
 # for data
 #dBJERFile = os.environ.get('CMSSW_BASE')+"/src/KinZfitter/HelperFunction/hists/Summer15_25nsV6_DATA_JER.db"
 process.jer = cms.ESSource("PoolDBESSource",
    CondDBSetup,
    connect = cms.string("sqlite_file:"+dBJERFile),
    toGet = cms.VPSet(

        cms.PSet(
            record = cms.string('JetResolutionRcd'),
            tag    = cms.string('JR_Summer15_25nsV6_MC_PtResolution_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs_pt') ),

        cms.PSet(
            record = cms.string('JetResolutionRcd'),
            tag    = cms.string('JR_Summer15_25nsV6_MC_PhiResolution_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs_phi') )
      )
 )

 process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

