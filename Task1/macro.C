#include "macro.h"


TGraphErrors* calculate(TProfile *prof_2pc, TCanvas *canvas, Color_t color)

TGraphErrors* calculate(TProfile *prof_2pc, TCanvas *canvas, Color_t color){
    
    Double_t x[prof_2pc->GetNcells()];
    Double_t y[prof_2pc->GetNcells()];
    
    Double_t ey[prof_2pc->GetNcells()];

    for (Int_t i=0;i<=prof_2pc->GetNcells()-4;++i) {
        x[i]= prof_2pc->GetBinCenter(i+1);
        y[i]= sqrt(prof_2pc->GetBinContent(i+1));
        ey[i]=prof_2pc->GetBinError(i+1)/(2*sqrt(prof_2pc->GetBinContent(i+1)));
    }

    auto *gr = new TGraphErrors(prof_2pc->GetNcells()-3,x,y,0,ey);
    
    gr->SetTitle("v2 as a function of centrality;" "centrality; ""v2");
    gr->SetMarkerStyle(kCircle);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(21);
    gr->DrawClone("PESame");

    return gr;
}

TGraphErrors* ratio_plot(TProfile *prof_2pc, Double_t y_pub[],  Double_t ey_pub[], TCanvas *canvas, Color_t color){
    
    Double_t x[prof_2pc->GetNcells()];
    Double_t y;
    Double_t ey;
    Double_t ry[prof_2pc->GetNcells()];
    Double_t rey[prof_2pc->GetNcells()];

    for (Int_t i=0;i<=prof_2pc->GetNcells()-4;++i) {
        x[i] = prof_2pc->GetBinCenter(i+1);
        y = sqrt(prof_2pc->GetBinContent(i+1));
        ey = prof_2pc->GetBinError(i+1)/(2*sqrt(prof_2pc->GetBinContent(i+1)));
        
        ry[i] = y/y_pub[i];
        rey[i] = sqrt(ey/y_pub[i]*ey/y_pub[i]+ (y/(y_pub[i]*y_pub[i])*ey_pub[i])*ey/y_pub[i]);
        cout<<xCross[i]<<"\t"<<x[i]<<"\t"<<y<<"\t"<<ey<<"\t"<<y_pub[i]<<"\t"<<ey_pub[i]<<"\n";
    }

    auto *gr = new TGraphErrors(prof_2pc->GetNcells()-3,xCross,ry,0,rey);
    
    gr->SetTitle("ratio_plot;" "centrality; ""v2");
    gr->SetMarkerStyle(kCircle);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(21);
    gr->DrawClone("PESame");
    gr->SetTitle(";Centrality ; obtained_results/real_results ");
    
    return gr;
}


TProfile* load_results(const char * fname){

    TFile* inputFile = TFile::Open(fname, "READ");

    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    TProfile *profile = (TProfile*)inputList->FindObject("results_v2");
    
    return profile;
}

void macro(){
    
    TProfile *prof_2pc= load_results("results_v2/AnalysisResults.root");
    TProfile *prof_2pc_Gap= load_results("results_v2/AnalysisResultsGap.root");
    TProfile *prof_2pc_Gap_Weight= load_results("results_v2/AnalysisResultsGapWeights.root");
    
    TProfile *real_prof_2pc= load_results("results_v2/AnalysisResults.root");
    TProfile *real_prof_2pc_Gap= load_results("results_v2/AnalysisResultsGap.root");
    TProfile *real_prof_2pc_Gap_Weight= load_results("results_v2/AnalysisResultsGapWeights.root");
    
    TCanvas *c1 = new TCanvas("c1","multigraph",200, 10,700,500);
    
    c1->Divide(1,2,0,0);
    c1->cd(1);
    c1->GetPad(1)->SetRightMargin(.01);
    

    TMultiGraph *mg = new TMultiGraph();


    TGraphErrors *grapherrors= calculate(prof_2pc, c1, kBlue);
    mg->Add(grapherrors);
    
    TGraphErrors *grapherrors_Gap= calculate(prof_2pc_Gap, c1, kRed);
    mg->Add(grapherrors_Gap);
    
    TGraphErrors *grapherrors_Gap_Weight= calculate(prof_2pc_Gap_Weight, c1, kBlack);
    mg->Add(grapherrors_Gap_Weight);
    
    mg->GetXaxis()->SetTitle("Centrality");
    mg->GetYaxis()->SetTitle("v_{n}");
    
    mg->Draw("ap");
    
    
    mg->SetMinimum(0. );

    
    TLegend leg(.1,.7,.3,.9,"Cases");
    leg.SetFillColor(0);
    leg.AddEntry(grapherrors,"No gap/weights");
    leg.AddEntry(grapherrors_Gap,"Gap, No weights ");
    leg.AddEntry(grapherrors_Gap_Weight,"Gap, Weights");
    leg.DrawClone("Same");
    
    c1->cd(2);
    
    c1->GetPad(2)->SetRightMargin(.01);

    TMultiGraph *mg2 = new TMultiGraph();

   /* TGraphErrors *ratio= ratio_plot(prof_2pc, real_prof_2pc, c1, kBlue);
    mg2->Add(ratio);
    
    TGraphErrors *ratio_Gap= ratio_plot(prof_2pc_Gap, real_prof_2pc_Gap, c1, kRed);
    mg2->Add(ratio_Gap);
    */
    
    TGraphErrors *ratio_Gap_Weight= ratio_plot(prof_2pc_Gap_Weight, v22Gap10Run2, v22Gap10Run2CombErr, c1, kBlack);
    mg2->Add(ratio_Gap_Weight);
    
    mg2->Draw("ap");

    mg2->SetMinimum(0. );
    
    mg2->GetXaxis()->SetTitle("Centrality");
    mg2->GetYaxis()->SetTitle("ratio");
    
    TLegend leg2(.1,.7,.3,.9,"Cases");
    leg2.SetFillColor(0);
    //!leg2.AddEntry(ratio,"No gap/weights");
    //!leg2.AddEntry(ratio_Gap,"Gap, No weights ");
    leg2.AddEntry(ratio_Gap_Weight,"Gap, Weights");
    leg2.DrawClone("Same");
    
    
    TList* export_v2 = new TList();
    export_v2->SetOwner(kTRUE);
    export_v2->Add(grapherrors_Gap_Weight);

    /*
    TH1D* weightsPlaceHolder = (TH1D*)fHistWeightsPhi->Clone("weightsPlaceHolder");
    weightsPlaceHolder->Scale(3.);
    weights->Add(weightsPlaceHolder);
*/
    TFile* output = new TFile("v2.root", "RECREATE");
    export_v2->Write("weightsList",TObject::kSingleKey+TObject::kOverwrite);

    
    
    return;
}

