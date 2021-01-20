#include <iostream>
#include <TTree.h>
#include <vector>
#include <TLorentzVector.h>
#include <fstream>
#include <string>
#include <cstring>
#include <THStack.h>
#include <stdio.h>
#include <TLatex.h>
#include <TF1.h>
#include <TH1.h>
#include <TGaxis.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooAbsArg.h"
#include "RooStats/ModelConfig.h"
#include "RooExponential.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;


void Analise(){


// Declaration of leaf types

Double_t event;   
TLorentzVector *dimuon_p4;
TLorentzVector *muonP_p4;
TLorentzVector *muonN_p4;

// List of branches
   
TBranch *b_event;
TBranch *b_dimuon_p4;
TBranch *b_muonP_p4;
TBranch *b_muonN_p4;
	  
//List of histos	  
	  
TH1F *hUpiMassCut = new TH1F("UpiMass", "m_{#Upislon}", 100, 8.5, 9.86);
hUpiMassCut->GetXaxis()->SetTitle("m_{#Upislon} (GeV)");
hUpiMassCut->GetYaxis()->SetTitle("Events");
   
TH1F *hdimuonMass = new TH1F("dimuonMass", "m_{#mu^{+}#mu^{-}}", 10000,0.2,300);
hdimuonMass->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
hdimuonMass->GetYaxis()->SetTitle("Events");
   
//pointer to the analyzed TTree or TChain

TChain* fChain = new TChain("oniaTree");
fChain->Add("Skim4.root/oniaTree;5");
fChain->Add("Skim4.root/oniaTree;4");
fChain->SetMakeClass(1);

// Set object pointer

dimuon_p4 = 0;
muonP_p4 = 0;
muonN_p4 = 0;

// Set branch addresses and branch pointers

fChain->SetBranchAddress("event", &event, &b_event);
fChain->SetBranchAddress("dimuon_p4", &dimuon_p4, &b_dimuon_p4);
fChain->SetBranchAddress("muonP_p4", &muonP_p4, &b_muonP_p4);
fChain->SetBranchAddress("muonN_p4", &muonN_p4, &b_muonN_p4);

int nentries = fChain->GetEntries();	
	 
cout << "Numero de entradas " << nentries << endl;
 
for(int i=0; i<nentries; i++){

	fChain->GetEntry(i);
		
	double etacut = abs(dimuon_p4->Eta());
	double mass = dimuon_p4->M();
	double pt1 = muonN_p4->Pt();
	double pt2 = muonP_p4->Pt();

	if(pt1 > 3.5 && pt2 > 3.5 && etacut < 2.1 && mass >= 8.5 && mass <= 9.86)hUpiMassCut->Fill(mass);

	hdimuonMass->Fill(mass);

}  //for i

fChain->GetFile()->Close();	

TCanvas *c = new TCanvas("c","c",10,10,700,900);

hUpiMassCut->Draw("E1");
hUpiMassCut->SetMarkerStyle(8);
hUpiMassCut->SetMarkerColor(1);

// Declare observable x

RooRealVar x("x","m_{#Upsilon}(GeV)",8.5,9.86); 

// Fit a p.d.f to the data

RooRealVar mean("mean","m_{#Upsilon}",9.4,9.2,9.6) ;
RooRealVar sigma("sigma","#sigma_{m_{#Upsilon}}",0.01,0.5) ;
RooRealVar alpha("alpha","#alpha",-5,-0.000001);
RooRealVar constant("constant","Constant",0,10);
RooCBShape crystalBall("crystalBall", "crystalBall", x, mean, sigma, alpha, constant);

RooGaussian sig("sig","Signal component",x,mean,sigma);

//RooRealVar eff_upi("eff_upi","The upi efficiency",0.98,0.00001,1.);
//RooRealVar lumi_upi("lumi_psi","The CMS luminosity",900,0,1500,"pb-1");
//RooRealVar cross_upi("cross_upi","The upi xsec",0.013,0.,40.,"pb");


//RooRealVar a("a","expo",-1,1);
RooRealVar a0("a0","a0",-1,1) ;
RooRealVar a1("a1","a1",-1,1) ;
//RooExponential expo("expo", "expo", x, a);
RooChebychev bkg("bkg","Background",x,RooArgList(a0,a1));
RooRealVar sigfrac("sigfrac","fraction of signal component",0.5,0.,1.);

RooRealVar b("b", "N_{sig}",0,19000);
RooRealVar s("s", "N_{bg}",0,15000);


RooExtendPdf esig("esig","esig",sig,s) ;
RooExtendPdf ebkg("ebkg","ebkg",bkg,b) ;

RooAddPdf model("model", "model",RooArgList(esig,ebkg));

// Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'

RooDataHist dh("dh","dh",x,Import(*hUpiMassCut));
model.fitTo(dh);

RooFitResult* fitresult = model.fitTo(dh, Extended(true), Save() );

// Make plot of binned dataset showing Poisson error bars (RooFit default)

RooPlot* frame = x.frame();
dh.plotOn(frame);
model.plotOn(frame);
//model.plotOn(frame,Components(bkg),LineStyle(kDashed));
double n_param, reduced_chi_square;

n_param = fitresult->floatParsFinal().getSize();

reduced_chi_square = frame->chiSquare(n_param);
model.plotOn(frame,Components(bkg),LineStyle(kDashed));
frame->Draw();

TPaveText* paveText = new TPaveText(0.5,0.5,0.9,0.9,"brNDC");
//paveText->SetBorderSize(1);
paveText->SetFillColor(kWhite);
paveText->SetFillStyle(1001);
paveText->SetTextSize(0.025);
paveText->AddText(Form("#chi^{2}/ndof = %f ", reduced_chi_square));
paveText->AddText(Form("Mean_{m_{#Upsilon}} = %.5f #pm %.5f GeV", mean.getVal(), mean.getError()));
paveText->AddText(Form("#sigma_{m_{#Upsilon}} = %.5f #pm %.5f GeV", sigma.getVal(), sigma.getError()));
//paveText->AddText(Form("#alpha = %.5f #pm %.5f GeV", alpha.getVal(), alpha.getError()));
//paveText->AddText(Form("Constant = %.5f #pm %.5f GeV", constant.getVal(), constant.getError()));
//paveText->AddText(Form("chi square = %f ", reduced_chi_square));
//paveText->AddText(Form("sigfrac = %f #pm %f", sigfrac.getVal(), sigfrac.getError()));
paveText->AddText(Form("a0 = %.5f #pm %.5f GeV", a0.getVal(), a0.getError()));
paveText->AddText(Form("a1 = %f #pm %f", a1.getVal(), a1.getError()));
//paveText->AddText(Form("a = %f #pm %f", a.getVal(), a.getError()));
paveText->AddText(Form("s = %.5f #pm %.5f GeV", s.getVal(), s.getError()));
paveText->AddText(Form("b = %f #pm %f", b.getVal(), b.getError()));
paveText->Draw();
//frame->addObject(t1);



}
