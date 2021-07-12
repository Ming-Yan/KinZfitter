#include <iostream>
using namespace RooFit;
void GENZfit(int emu)
{
  //set file name
  TString type;
  if(fabs(emu)==11)  type = "ele";
  if(fabs(emu)==13)  type = "mu";
  char fname[500];
  sprintf(fname,"HZg_ggF_125GeV_ext1_M125_13TeV_powheg2_pythia8_%s.txt",type.Data());
  ofstream outfile(fname);
  //import the data
  TFile *f = TFile::Open(Form("../../outfile/GEN/MC_wcut_%s_wFSR_HZg_ggF_125GeV.root",type.Data()));
  TTree *t = (TTree*)f->Get("tMC");
	
  //build GenZ model
  RooRealVar mcmll("mcmll" ,"mcmll", 40, 120, "GeV");
  RooRealVar eormu("eormu" ,"eormu", 10, 14, "GeV");
  RooDataSet Zmass("Zmass", "Zmass", RooArgSet(mcmll), Import(*t));
  /*RooRealVar meanCB("meanCB","meanCB", 8.64761e+01, 80,100);
  RooRealVar sigmaCB("sigmaCB","sigmaCB", 4.99069e+00, 4., 10.);
  RooRealVar alphaCB("alphaCB","alphaCB",5.51008e-01, 0.1,1.2);
  RooRealVar nCB("nCB","nCB", 4.99983e+01, 0, 100);
  RooRealVar meanGauss1("meanGauss1","meanGauss1",70., 40, 110);
  RooRealVar sigmaGauss1("sigmaGauss1","sigmaGauss1",8.08747e+00, 0, 20);
  RooRealVar f1("f1","f1", 0.00001,1.0);*/
  /*RooRealVar meanGauss2("meanGauss2","meanGauss2",9.10697e+01, 80, 100);
  RooRealVar sigmaGauss2("sigmaGauss2","sigmaGauss2",8.03251e-01, 0.1, 5.1);
  RooRealVar f2("f2","f2",7.40547e-01, 0.0001, 1.0);
  RooRealVar meanGauss3("meanGauss3","meanGauss3",9.05542e+01, 80, 100);
  RooRealVar sigmaGauss3("sigmaGauss3","sigmaGauss3", 1.86069e+00, 0, 10);
  RooRealVar f3("f3","f3",7.38652e-01,0.0001,1.0);*/
  RooRealVar bwMean("bwMean", "m_{Z^{0}}", 91.187);
  RooRealVar bwGamma("bwGamma", "#Gamma", 2.5,1.0,3.0);
  RooRealVar sg("sg", "sg", 4.99069e+00, 0.1, 10.);
  RooRealVar a("a", "a", 4.99069e+00, 0.1, 10.);
  RooRealVar n("n", "n", -0.5,2.0);
  //RooRealVar mean("mean","mean",0.,10);
  RooRealVar mean("mean","mean",-1,10);
  RooCBShape CB("CB","CB",mcmll,mean,sg,a,n);
  RooRealVar f1("f1","f1",0.01,1.0);

  RooRealVar mean2("mean2","mean2",0.,20);
  RooRealVar sigma("sigma","sigma",0.1,10.);
  RooRealVar f2("f2","f2",0.000001,1.0000);

  RooGenericPdf RelBW("RelBW","1/( pow(mcmll*mcmll-bwMean*bwMean,2)+pow(mcmll,4)*pow(bwGamma/bwMean,2) )", RooArgSet(mcmll,bwMean,bwGamma) );

  //RooAddPdf *RelBWxCB = new RooAddPdf("RelBWxCB","RelBWxCB", RelBW, CB, f1);
  RooFFTConvPdf *RelBWxCB = new RooFFTConvPdf("RelBWxCB","RelBWxCB", mcmll, RelBW,CB);
  RooGaussian gauss("gauss","gauss",mcmll,mean2,sigma);
  //RooAddPdf *RelBWxCBxgauss = new RooAddPdf("RelBWxCBxgauss","RelBWxCBxgauss", *RelBWxCB, gauss, f2);
  //RooFFTConvPdf *RelBWxCBxgauss = new RooAddPdf("RelBWxCBxgauss","RelBWxCBxgauss", mcmll,*RelBWxCB, gauss);
  //RooAddPdf *RelBWxCBxgauss = new RooAddPdf("RelBWxCBxgauss","RelBWxCBxgauss", RelBWxCB, gauss, f2);
  
  /*RooCBShape* singleCB = new RooCBShape("singleCB", "", mcmll, meanCB, sigmaCB, alphaCB, nCB);
  RooGaussian* gaussShape1 = new RooGaussian("gaussShape1", "", mcmll, meanGauss1, sigmaGauss1);
  RooAddPdf* CBplusGauss = new RooAddPdf("CBplusGauss", "", *singleCB, *gaussShape1, f1);*/
  /*RooGaussian* gaussShape2 = new RooGaussian("gaussShape2", "", mcmll, meanGauss2, sigmaGauss2);
  RooAddPdf* CBplusGaussplusGauss = new RooAddPdf("CBplusGaussplusGauss", "", *CBplusGauss, *gaussShape2, f2);
  RooGaussian* gaussShape3 = new RooGaussian("gaussShape3", "", mcmll, meanGauss3, sigmaGauss3);
  RooAddPdf* CBplusGaussplusGaussplusGauss = new RooAddPdf("CBplusGaussplusGaussplusGauss", "", *CBplusGaussplusGauss, *gaussShape3, f3);*/
  //RooFitResult* GENZ = CBplusGaussplusGaussplusGauss->fitTo(Zmass, Save(kTRUE));
  //RooFitResult* GENZ = CBplusGauss->fitTo(Zmass, Save(kTRUE));
  RooFitResult* GENZ = RelBWxCB->fitTo(Zmass, Save(kTRUE));
  //RooFitResult* GENZ = RelBWxCBxgauss->fitTo(Zmass, Save(kTRUE));
    RooPlot* xframe4 = mcmll.frame(40,120) ;
  Zmass.plotOn(xframe4,Binning(80), RooFit::Name("Zmass")) ;
  //CBplusGauss->plotOn(xframe4,RooFit::Name("CBplusGauss"),LineColor(TColor::GetColor("#0e95e9")));
  /*CBplusGaussplusGaussplusGauss->plotOn(xframe4,RooFit::Name("CBplusGaussplusGaussplusGauss"),LineColor(TColor::GetColor("#d9594c")));
  CBplusGaussplusGaussplusGauss->plotOn(xframe4,RooFit::Components("CBplusGaussplusGauss"),LineStyle(kDashed),LineColor(TColor::GetColor("#ce8d66")));
  CBplusGaussplusGaussplusGauss->plotOn(xframe4,RooFit::Components("gaussShape3"),LineStyle(kDashed),LineColor(TColor::GetColor("#b7b868")));
  CBplusGaussplusGaussplusGauss->paramOn(xframe4,Layout(0.6,0.96,0.95));*/
  //RelBWxCBxgauss->plotOn(xframe4,RooFit::Name("RelBWxCBxgauss"),LineColor(TColor::GetColor("#df9f1f")));
  RelBWxCB->plotOn(xframe4,RooFit::Name("RelBWxCB"),LineColor(TColor::GetColor("#df9f1f")));
  RelBWxCB->plotOn(xframe4,RooFit::Components("RelBW"),LineStyle(kDashed),LineColor(TColor::GetColor("#0e95e9")));
  RelBWxCB->plotOn(xframe4,RooFit::Components("CB"),LineStyle(kDashed),LineColor(TColor::GetColor("#f35c5d")));
  //RelBWxCB->plotOn(xframe4,RooFit::Components("gauss"),LineStyle(kDashed),LineColor(TColor::GetColor("#1aa450")));
  RelBWxCB->paramOn(xframe4,Layout(0.1,0.4,0.9));
 
  //CBplusGauss->paramOn(xframe4,Layout(0.1,0.4,0.9));
  xframe4->Draw();
  outfile<<"alpha\t"<<a.getValV()<<endl;
  outfile<<"mean\t"<<mean.getValV()<<endl;
  outfile<<"n\t"<<n.getValV()<<endl;
  outfile<<"sigma\t"<<sg.getValV()<<endl;
  outfile<<"bwgamma\t"<<bwGamma.getValV()<<endl;
  /*RooWorkspace *w  = new RooWorkspace("w","");
  w->import(Zmass);
  w->import(*CBplusGaussplusGaussplusGauss);
  w->import(*GENZ);
  w->Print() ;	*/
  // w->WritetoFile("ele_root");
}
