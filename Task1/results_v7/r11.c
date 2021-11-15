#include "iostream"
#include "fstream"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "Th2.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TApplication.h"
#include "stdio.h"
#include "cstring"
#include "TMath.h"
#include "TGaxis.h"
#include "TLatex.h"

using namespace std;


void r11() {
 
   
   TString inputfile1 = "./z+jets_mc_co_matrix.root";

   TFile file1(inputfile1);

   TH2D *h1;
   TH2D *h2;
   TH2D *h3;
  
   h1 = (TH2D*) file1.Get("z+jets reco_gen_z_energy_co");
   h2 = (TH2D*) file1.Get("z+jets gen_dressed_z_energy_co");
   h3 = (TH2D*) file1.Get("z+jets reco_dressed_z_energy_co");
  
   auto c1 = new TCanvas("c1", "A ratio example",800,800);
   c1->SetTitle("Z Mass Correlation");
   c1->Divide(2,2);
   
   
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 0.5, 1);
   pad1->SetTopMargin(0.15);
   pad1->SetBottomMargin(0.15); // Upper and lower plot are joined
   pad1->SetLeftMargin(0.15);
   pad1->SetRightMargin(0.15);
   pad1->SetLogz(1);
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();
   h1->Draw("colz");
   
   c1->cd();
   
   TPad *pad2 = new TPad("pad2", "pad2", 0.5, 0.5, 1, 1);
   pad2->SetTopMargin(0.15);
   pad2->SetBottomMargin(0.15); // Upper and lower plot are joined
   pad2->SetLeftMargin(0.15);
   pad2->SetRightMargin(0.15);
   pad2->SetLogz(1);
   pad2->Draw();             // Draw the upper pad: pad1
   pad2->cd();
   h2->Draw("colz");
   
   c1->cd();

   TPad *pad3 = new TPad("pad3", "pad3", 0, 0, 0.5, 0.5);
   pad3->SetTopMargin(0.15);
   pad3->SetBottomMargin(0.15); // Upper and lower plot are joined
   pad3->SetLeftMargin(0.15);
   pad3->SetRightMargin(0.15);
   pad3->SetLogz(1);
   pad3->Draw();             // Draw the upper pad: pad1
   pad3->cd();
   h3->Draw("colz");
   
   c1->cd();
   c1->SetTitle("z energy");
   
   c1->SaveAs("zenergy.png");
}
