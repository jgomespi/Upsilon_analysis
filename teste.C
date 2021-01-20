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


void teste(){


// Declaration of leaf types

Double_t event;   
TLorentzVector *dimuon_p4;
TLorentzVector *muonP_p4;
TLorentzVector *muonN_p4;
double mass;
TFile f("massupi.root", "NEW");
 TTree *tree = new TTree("T","T");

// List of branches
   
TBranch *b_event;
TBranch *b_dimuon_p4;
TBranch *b_muonP_p4;
TBranch *b_muonN_p4;

tree->Branch("mass",&mass,"mass/D");
	  
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
	 mass = dimuon_p4->M();
	double pt1 = muonN_p4->Pt();
	double pt2 = muonP_p4->Pt();

	if(pt1 > 3.5 && pt2 > 3.5 && etacut < 2.1 && mass >= 8.5 && mass <= 9.86){
	
		hUpiMassCut->Fill(mass);
		tree->Fill();
	}

	hdimuonMass->Fill(mass);

}  //for i

tree->Print();

f.Write();
 
    // Close the file. Note that this is automatically done when you leave
    // the application upon file destruction.
    f.Close();
fChain->GetFile()->Close();

}
