# KinZfitter

Using Z mass kinematic constraint(s) to refit lepton momenta in H-to-Zgamma
-> code based on H->ZZ*

- Code structure

KinZfitter : the class read the inputs from leptons and fsr photons that form the Higgs Candidate,, do the refitting, and get the refitted results
HelperFunction : the class that read lepton/photon pT errors by accessing the electron/muon with [ggAnalysis][https://github.com/cmkuo/ggAnalysis/tree/master/ggNtuplizer] structure
To include the refit in your analyzer:

0.Check out package

  git clone <https://github.com/Ming-Yan/KinZfitter>

In your main analyzer:

1. Include the macro in the `xAna.C`, include after `untuplizer.h`, and load macro in `runAna.C` by

```
gROOT->LoadMacro("KinZfitter/KinZfitter/src/KinZfitter_v1.cpp");
gROOT->LoadMacro("KinZfitter/HelperFunction/src/HelperFunction.cc");
```

2. Initialize tool by `KinZfitter* kinZfitter = new KinZfitter(isData,period);` before event loop, periods can be 2016, 2017, 2018

3. After selecting the dilepton ID and FSR photon ID put into `vector<int>`, put empty `vector<float> pTerror`, channel is 0/1 to be electron/muon. 

```
kinZfitter->Setup(data,lepID_kin, fsrID_kin, pTerr, channel, isData, selectedLeptons, selectedFsrMap,corr)
```
The selectedLeptons, selectedFsrMap are the map `std::map<unsigned int, TLorentzVector> selectedLeptons;` by putting the TLorentzVecor `selectedLeptons[i] = kinlep[i]`

4. Do refit by `kinZfitter->KinRefitZ();` 

5. extract the dilepton TLorentzVector by `kinZfitter->GetRefitP4s()`, get Refit Z mass by `kinZfitter->GetRefitMZ1()`


