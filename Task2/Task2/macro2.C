TGraphErrors* calculate(TProfile *prof_2pc, TCanvas *canvas, Color_t color){
    
    Double_t x[prof_2pc->GetNcells()];
    Double_t y[prof_2pc->GetNcells()];
    
    Double_t ey[prof_2pc->GetNcells()];

    for (Int_t i=0;i<=prof_2pc->GetNcells()-1;++i) {
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


TGraphErrors* ratio_plot(TProfile *prof_2pc, TProfile *prof_2pc_real, TCanvas *canvas, Color_t color){
    
    Double_t x[prof_2pc->GetNcells()];
    Double_t y[prof_2pc->GetNcells()];
    Double_t ey[prof_2pc->GetNcells()];

    for (Int_t i=0;i<=prof_2pc->GetNcells()-1;++i) {
        x[i]= prof_2pc->GetBinCenter(i+1);
        y[i]= sqrt(prof_2pc->GetBinContent(i+1));
        ey[i]=prof_2pc->GetBinError(i+1)/(2*sqrt(prof_2pc->GetBinContent(i+1)));
    }
    
    Double_t x_real[prof_2pc_real->GetNcells()];
    Double_t y_real[prof_2pc_real->GetNcells()];
    Double_t ey_real[prof_2pc_real->GetNcells()];

    for (Int_t i=0;i<=prof_2pc_real->GetNcells()-1;++i) {
        y_real[i]= sqrt(prof_2pc_real->GetBinContent(i+1));
        ey_real[i]=prof_2pc_real->GetBinError(i+1)/(2*sqrt(prof_2pc->GetBinContent(i+1)));
    }
    
    Double_t ry[prof_2pc->GetNcells()];
    Double_t rey[prof_2pc->GetNcells()];
    
    
    for (Int_t i=0;i<=prof_2pc_real->GetNcells()-1;++i) {
        ry[i]= y[i]/y_real[i];
        rey[i]= sqrt(ey[i]/y_real[i]*ey[i]/y_real[i]+ (y[i]/(y_real[i]*y_real[i])*ey_real[i])*ey[i]/y_real[i]);
    }

    auto *gr = new TGraphErrors(prof_2pc->GetNcells()-3,x,ry,0,rey);
    
    gr->SetTitle("ratio_plot;" "centrality; ""v2");
    gr->SetMarkerStyle(kCircle);
    gr->SetMarkerColor(color);
    gr->SetMarkerStyle(21);
    gr->DrawClone("PESame");
    gr->SetTitle(";Centrality ; obtained_results/real_results ");
    
    return gr;
}


TProfile2D* load_results(const char * fname){

    TFile* inputFile = TFile::Open(fname, "READ");

    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    TProfile2D *profile = (TProfile2D*)inputList->FindObject("results_BA");
    
    return profile;
}

TGraphErrors* load_v2(const char * fname ){

    TFile* inputFile = TFile::Open(fname, "READ");

    //!inputFile->cd("weightsList;1");
    TList* inputList = (TList*) gDirectory->Get("weightsList;1");
    TGraphErrors *profile = (TGraphErrors*)inputList->FindObject("Graph");
    
    return profile;
}

void macro2(){
    
    TProfile2D *prof_2pc= load_results("AnalysisResults.root");
    
    prof_2pc->ProfileY()->Draw();
    
    
    /*
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

    TGraphErrors *ratio= ratio_plot(prof_2pc, real_prof_2pc, c1, kBlue);
    mg2->Add(ratio);
    
    TGraphErrors *ratio_Gap= ratio_plot(prof_2pc_Gap, real_prof_2pc_Gap, c1, kRed);
    mg2->Add(ratio_Gap);
    
    TGraphErrors *ratio_Gap_Weight= ratio_plot(prof_2pc_Gap_Weight, real_prof_2pc_Gap_Weight, c1, kBlack);
    mg2->Add(ratio_Gap_Weight);
    
    mg2->Draw("ap");

    mg2->SetMinimum(0. );
    
    mg2->GetXaxis()->SetTitle("Centrality");
    mg2->GetYaxis()->SetTitle("obtained_results/real_results");
    
    TLegend leg2(.1,.7,.3,.9,"Cases");
    leg2.SetFillColor(0);
    leg2.AddEntry(ratio,"No gap/weights");
    leg2.AddEntry(ratio_Gap,"Gap, No weights ");
    leg2.AddEntry(ratio_Gap_Weight,"Gap, Weights");
    leg2.DrawClone("Same");
    */
    
    
    return;
}

