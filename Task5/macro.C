#include "macro.h"

void calculate_SC(TProfile *prof_results_SC, TProfile *prof_results_NSC, TProfile *prof_4pc_m_n, TProfile *prof_2pc_m,  TProfile *prof_2pc_n, TProfile *prof_2pc_m_gap,  TProfile *prof_2pc_n_gap){
    
    Double_t x;
    Double_t y;
    Double_t ey;
    
    for (Int_t i=1;i<=prof_4pc_m_n->GetNcells()-1;++i) {
        x= prof_4pc_m_n->GetBinCenter(i);
        
        y= prof_4pc_m_n->GetBinContent(i)-prof_2pc_m->GetBinContent(i)*prof_2pc_n->GetBinContent(i);
        ey= sqrt(prof_4pc_m_n->GetBinError(i)*prof_4pc_m_n->GetBinError(i) + prof_2pc_m->GetBinError(i)*prof_2pc_n->GetBinContent(i)*prof_2pc_m->GetBinError(i)*prof_2pc_n->GetBinContent(i)+ prof_2pc_m->GetBinContent(i)*prof_2pc_n->GetBinError(i)*prof_2pc_m->GetBinContent(i)*prof_2pc_n->GetBinError(i));
        
        
        //!cout<<x<<"\t"<<y<<"\t"<<ey<<"\n";
        if(ey>0.){
            //cout<<c2<<"\t"<<y<<"\t"<<ey<<"\n";
            prof_results_SC->Fill(x, y, 1/ey);
        }
        
        y/=(prof_2pc_m_gap->GetBinContent(i)*prof_2pc_n_gap->GetBinContent(i));
        ey = sqrt((ey/(prof_2pc_m_gap->GetBinContent(i)*prof_2pc_n_gap->GetBinContent(i)))*(ey/(prof_2pc_m_gap->GetBinContent(i)*prof_2pc_n_gap->GetBinContent(i))) +  (y*prof_2pc_m_gap->GetBinError(i)/(pow(prof_2pc_m_gap->GetBinContent(i),2)*prof_2pc_n_gap->GetBinContent(i)))*(y*prof_2pc_m_gap->GetBinError(i)/(pow(prof_2pc_m_gap->GetBinContent(i),2)*prof_2pc_n_gap->GetBinContent(i)))+ (y*prof_2pc_n_gap->GetBinError(i)/(pow(prof_2pc_n_gap->GetBinContent(i),2)*prof_2pc_m_gap->GetBinContent(i)))*(y*prof_2pc_n_gap->GetBinError(i)/(pow(prof_2pc_n_gap->GetBinContent(i),2)*prof_2pc_m_gap->GetBinContent(i))));
        
        if(ey>0.){
            //cout<<c2<<"\t"<<y<<"\t"<<ey<<"\n";
            prof_results_NSC->Fill(x, y, 1/ey);
        }
    }
    return;
}






TMultiGraph * represent(TList* fOutputList, TList* fErrorList){
    


    
    TMultiGraph *mg1 = new TMultiGraph();
    
    TLegend leg(.1,.7,.3,.9,"Alice Pb-Pb #sqrt{s_{NN}}=5.02 TeV, No gap ");
    leg.SetFillColor(0);
    
    //!Double_t mult_factor[]={1.,1.,100.,1.,1.};
    Double_t mult_factor[]={1.,1.,1.,1.,1.};  //!               For the NSC there is not rescaling
    
    Color_t colors[]={kRed,kBlue,kBlack, kGreen, kViolet};
    
    for(Int_t j=0; j<fOutputList->GetEntries();j++){
        TProfile* Mean_profile = (TProfile*)fOutputList->At(j);
        TProfile* Error_profile = (TProfile*)fErrorList->At(j);
        
        Double_t x[Mean_profile->GetNcells()];
        Double_t y[Mean_profile->GetNcells()];
        Double_t ey[Mean_profile->GetNcells()];
        
        for(Int_t i=0; i< Mean_profile->GetNcells(); i++){
            x[i]= Mean_profile->GetBinCenter(i+1);
            y[i]= Mean_profile->GetBinContent(i+1)*mult_factor[j];
            ey[i] = Error_profile->GetBinError(i+1)*mult_factor[j];
            
            
        }
        auto gr= new TGraphErrors(Mean_profile->GetNcells()-1,x,y,0,ey);
        
        gr->SetMarkerStyle(kCircle);
        gr->SetMarkerColor(colors[j]);
        gr->SetMarkerStyle(21);
        
        
        mg1->Add(gr);
        leg.AddEntry(gr,names[j]);
    }

    mg1->SetTitle("Symetric cumulants; V0M; SC");
    
    mg1->Draw("ap");
    leg.DrawClone("Same");
     
    

    return mg1;
}


TProfile* load_results(const char * fname, const char * particles){

    TFile* inputFile = TFile::Open(fname, "READ");

    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    TProfile *profile = (TProfile*)inputList->FindObject(particles);
    
    return profile;
}


void macro(){
    

    const Double_t bin_edges[] = {0.,5.,10.,20.,30.,40.,50.};
    
     Mean_profile_SC_5_2 = new TProfile(" Mean 5_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{2}",6, bin_edges);
     Mean_profile_SC_5_3 = new TProfile(" Mean 5_3 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",6, bin_edges);
     Mean_profile_SC_4_3 = new TProfile(" Mean 4_3 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{6}",6, bin_edges);
     Mean_profile_SC_4_2 = new TProfile(" Mean 4_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",6, bin_edges);
     Mean_profile_SC_3_2 = new TProfile(" Mean 3_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",6, bin_edges);
    
     Error_profile_SC_5_2 = new TProfile(" Error 5_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{2}",6, bin_edges);
     Error_profile_SC_5_3 = new TProfile(" Error 5_3 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",6, bin_edges);
     Error_profile_SC_4_3 = new TProfile(" Error 4_3 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{6}",6, bin_edges);
     Error_profile_SC_4_2 = new TProfile(" Error 4_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",6, bin_edges);
     Error_profile_SC_3_2 = new TProfile(" Error 3_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",6, bin_edges);
    
     Mean_profile_NSC_5_2 = new TProfile(" Mean 5_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{2}",6, bin_edges);
     Mean_profile_NSC_5_3 = new TProfile(" Mean 5_3 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",6, bin_edges);
     Mean_profile_NSC_4_3 = new TProfile(" Mean 4_3 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{6}",6, bin_edges);
     Mean_profile_NSC_4_2 = new TProfile(" Mean 4_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",6, bin_edges);
     Mean_profile_NSC_3_2 = new TProfile(" Mean 3_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",6, bin_edges);
    
     Error_profile_NSC_5_2 = new TProfile(" Error 5_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{2}",6, bin_edges);
     Error_profile_NSC_5_3 = new TProfile(" Error 5_3 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",6, bin_edges);
     Error_profile_NSC_4_3 = new TProfile(" Error 4_3 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{6}",6, bin_edges);
     Error_profile_NSC_4_2 = new TProfile(" Error 4_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",6, bin_edges);
     Error_profile_NSC_3_2 = new TProfile(" Error 3_2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",6, bin_edges);

    

    TFile* inputFile = TFile::Open("AnalysisResults_grid.root", "READ");

    TFile* inputFile3 = TFile::Open("~/Desktop/Thesis-training-/Task1/results_v2/AnalysisResultsGapWeights.root", "READ");
    TFile* inputFile4 = TFile::Open("~/Desktop/Thesis-training-/Task1/results_v3/AnalysisResultsGapWeights.root", "READ");
    TFile* inputFile5 = TFile::Open("~/Desktop/Thesis-training-/Task1/results_v4/AnalysisResultsGapWeights.root", "READ");
    TFile* inputFile6 = TFile::Open("~/Desktop/Thesis-training-/Task1/results_v5/AnalysisResultsGapWeights.root", "READ");
    
    
    
    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    inputFile3->cd("FlowExampleTask");
    TList* inputList3 = (TList*) gDirectory->Get("Output");
    inputFile4->cd("FlowExampleTask");
    TList* inputList4 = (TList*) gDirectory->Get("Output");
    inputFile5->cd("FlowExampleTask");
    TList* inputList5 = (TList*) gDirectory->Get("Output");
    inputFile6->cd("FlowExampleTask");
    TList* inputList6 = (TList*) gDirectory->Get("Output");
    
    for(Int_t i(0); i<10; i++){
        
        TProfile *prof_5_2 = (TProfile*)inputList->FindObject(name_SC_5_2[i]);
        TProfile *prof_5_3 = (TProfile*)inputList->FindObject(name_SC_5_3[i]);
        TProfile *prof_4_3 = (TProfile*)inputList->FindObject(name_SC_4_3[i]);
        TProfile *prof_3_2 = (TProfile*)inputList->FindObject(name_SC_3_2[i]);
        TProfile *prof_4_2 = (TProfile*)inputList->FindObject(name_SC_4_2[i]);
        
        TProfile *prof_v2 = (TProfile*)inputList->FindObject(name_2pc_2[i]);
        TProfile *prof_v3 = (TProfile*)inputList->FindObject(name_2pc_3[i]);
        TProfile *prof_v4 = (TProfile*)inputList->FindObject(name_2pc_4[i]);
        TProfile *prof_v5 = (TProfile*)inputList->FindObject(name_2pc_5[i]);
        
        TProfile *prof_v2_gap = (TProfile*)inputList3->FindObject(name_vn[i]);
        TProfile *prof_v3_gap = (TProfile*)inputList4->FindObject(name_vn[i]);
        TProfile *prof_v4_gap = (TProfile*)inputList5->FindObject(name_vn[i]);
        TProfile *prof_v5_gap = (TProfile*)inputList6->FindObject(name_vn[i]);
 
        calculate_SC(Error_profile_SC_3_2, Error_profile_NSC_3_2, prof_3_2, prof_v3, prof_v2, prof_v2_gap, prof_v3_gap);
        calculate_SC(Error_profile_SC_4_2, Error_profile_NSC_4_2, prof_4_2, prof_v2, prof_v4, prof_v2_gap, prof_v4_gap);
        calculate_SC(Error_profile_SC_5_2, Error_profile_NSC_5_2, prof_5_2, prof_v2, prof_v5, prof_v2_gap, prof_v5_gap);
        calculate_SC(Error_profile_SC_5_3, Error_profile_NSC_5_3, prof_5_3, prof_v3, prof_v5, prof_v3_gap, prof_v5_gap);
        calculate_SC(Error_profile_SC_4_3, Error_profile_NSC_4_3, prof_4_3, prof_v4, prof_v3, prof_v4_gap, prof_v3_gap);
        
    }
    
    
    TProfile *prof_5_2 = (TProfile*)inputList->FindObject("4_pc_5_2");
    TProfile *prof_5_3 = (TProfile*)inputList->FindObject("4_pc_5_3");
    TProfile *prof_4_3 = (TProfile*)inputList->FindObject("4_pc_4_3 ");
    TProfile *prof_4_2 = (TProfile*)inputList->FindObject("4_pc_4_2 ");
    TProfile *prof_3_2 = (TProfile*)inputList->FindObject("4_pc_3_2 ");
    
    TProfile *prof_v2 = (TProfile*)inputList->FindObject("2_pc_2 ");
    TProfile *prof_v3 = (TProfile*)inputList->FindObject("2_pc_3 ");
    TProfile *prof_v4 = (TProfile*)inputList->FindObject("2_pc_4 ");
    TProfile *prof_v5 = (TProfile*)inputList->FindObject("2_pc_5 ");
    
    TProfile *prof_v2_gap = (TProfile*)inputList3->FindObject("results_value_2");
    TProfile *prof_v3_gap = (TProfile*)inputList4->FindObject("results_value_2");
    TProfile *prof_v4_gap = (TProfile*)inputList5->FindObject("results_value_2");
    TProfile *prof_v5_gap = (TProfile*)inputList6->FindObject("results_value_2");
    
    
    calculate_SC(Mean_profile_SC_3_2, Mean_profile_NSC_3_2, prof_3_2, prof_v3, prof_v2, prof_v2_gap, prof_v4_gap);
    calculate_SC(Mean_profile_SC_4_2, Mean_profile_NSC_4_2, prof_4_2, prof_v2, prof_v4, prof_v2_gap, prof_v4_gap);
    calculate_SC(Mean_profile_SC_5_2, Mean_profile_NSC_5_2, prof_5_2, prof_v2, prof_v5, prof_v2_gap, prof_v5_gap);
    calculate_SC(Mean_profile_SC_5_3, Mean_profile_NSC_5_3, prof_5_3, prof_v3, prof_v5, prof_v3_gap, prof_v5_gap);
    calculate_SC(Mean_profile_SC_4_3, Error_profile_NSC_4_3, prof_4_3, prof_v4, prof_v3, prof_v4_gap, prof_v3_gap);
    
    
    
    fOutputListSC = new TList();
    fOutputListSC->Add(Mean_profile_SC_5_2);
    fOutputListSC->Add(Mean_profile_SC_5_3);
    fOutputListSC->Add(Mean_profile_SC_4_3);
    fOutputListSC->Add(Mean_profile_SC_4_2);
    fOutputListSC->Add(Mean_profile_SC_3_2);
    
    
    fErrorListSC = new TList();
    fErrorListSC->Add(Error_profile_SC_5_2);
    fErrorListSC->Add(Error_profile_SC_5_3);
    fErrorListSC->Add(Error_profile_SC_4_3);
    fErrorListSC->Add(Error_profile_SC_4_2);
    fErrorListSC->Add(Error_profile_SC_3_2);

    fOutputListNSC = new TList();
    fOutputListNSC->Add(Mean_profile_NSC_5_2);
    fOutputListNSC->Add(Mean_profile_NSC_5_3);
    fOutputListNSC->Add(Mean_profile_NSC_4_3);
    fOutputListNSC->Add(Mean_profile_NSC_4_2);
    fOutputListNSC->Add(Mean_profile_NSC_3_2);
    
    
    fErrorListNSC = new TList();
    fErrorListNSC->Add(Error_profile_NSC_5_2);
    fErrorListNSC->Add(Error_profile_NSC_5_3);
    fErrorListNSC->Add(Error_profile_NSC_4_3);
    fErrorListNSC->Add(Error_profile_NSC_4_2);
    fErrorListNSC->Add(Error_profile_NSC_3_2);

    TCanvas *c1 = new TCanvas("c1","multigraph",200, 10,700,500);
    
    /*
    c1->Divide(1,2,0,0);
    c1->cd(1);
    c1->GetPad(1)->SetRightMargin(.01);
     */

    auto *c_SC = represent(fOutputListSC, fErrorListSC);
    
    c1->Update();
    c1->Range( 0., -2e-07, 60., 2e-07 );
    
    auto *l = new TLine(0.,0.,60.,0.);
    
    l->SetLineWidth(2);
    l->Draw();
    
    TCanvas *c2 = new TCanvas("c2","multigraph",200, 10,700,500);
    
    /*
    c1->Divide(1,2,0,0);
    c1->cd(1);
    c1->GetPad(1)->SetRightMargin(.01);
     */

    auto *c_NSC = represent(fOutputListNSC, fErrorListNSC);
    
    c2->Update();
    c2->Range( 0., -2e-07, 60., 2e-07 );
    
    l->Draw();
    
    
    

    
    return;
}

