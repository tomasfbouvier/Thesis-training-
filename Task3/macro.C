#include "published.h"


Double_t c_4(Double_t pc_4_ref, Double_t pc_2_ref){
    return pc_4_ref - 2*pc_2_ref*pc_2_ref;
}

Double_t c_4_err(Double_t pc_4_ref, Double_t pc_2_ref, Double_t pc_4_ref_err, Double_t pc_2_ref_err){
    
    Double_t aux= pc_4_ref_err*pc_4_ref_err + (2*pc_2_ref_err*pc_2_ref)*(2*pc_2_ref*pc_2_ref_err)*(2*pc_2_ref_err*pc_2_ref)*(2*pc_2_ref*pc_2_ref_err);
    
    return sqrt(aux);
}

Double_t c_6(Double_t pc_6_ref, Double_t pc_4_ref, Double_t pc_2_ref){
    return pc_6_ref - 9*pc_2_ref*pc_4_ref + 12*pow(pc_2_ref, 3);
}

Double_t c_6_err(Double_t pc_6_ref, Double_t pc_4_ref, Double_t pc_2_ref, Double_t pc_6_ref_err, Double_t pc_4_ref_err, Double_t pc_2_ref_err){
    
    Double_t aux = pc_6_ref_err*pc_6_ref_err + (-9*pc_2_ref_err*pc_4_ref + 3*12*pow(pc_2_ref, 2)*pc_2_ref_err)*(-9*pc_2_ref_err*pc_4_ref + 3*12*pow(pc_2_ref, 2)*pc_2_ref_err) + 9*pc_2_ref*pc_4_ref_err*9*pc_2_ref*pc_4_ref_err ;
    
    return sqrt(aux);
}

Double_t c_8(Double_t pc_8_ref, Double_t pc_6_ref, Double_t pc_4_ref, Double_t pc_2_ref){
    return  - pc_8_ref +16*pc_6_ref*pc_2_ref +18*pc_4_ref*pc_4_ref-144*pc_4_ref*pc_2_ref*pc_2_ref+144*pow(pc_2_ref, 4.);
}

Double_t c_8_err(Double_t pc_8_ref, Double_t pc_6_ref, Double_t pc_4_ref, Double_t pc_2_ref, Double_t pc_8_ref_err, Double_t pc_6_ref_err, Double_t pc_4_ref_err, Double_t pc_2_ref_err){
    
    Double_t aux = pc_8_ref_err*pc_8_ref_err + pow(16*pc_6_ref_err*pc_2_ref,2.) + pow((18*2*pc_4_ref-144*pc_2_ref*pc_2_ref)*pc_4_ref_err,2.) + pow(( 16*pc_6_ref-144*pc_4_ref*pc_2_ref + 4*144*pow(pc_2_ref,3))*pc_2_ref_err ,2.);
    return sqrt(aux);
}

void calculate_v_2_2_int(TProfile *prof_results, TProfile *prof_2pc){
    
    Double_t x;
    Double_t y;
    Double_t ey;
    
    for (Int_t i=1;i<=prof_2pc->GetNcells()-3;++i) {
        x= prof_2pc->GetBinCenter(i);
        
        Double_t c2= prof_2pc->GetBinContent(i);
        Double_t c2_err=  prof_2pc->GetBinError(i);
        
        y= pow(c2,1./2.)  ;
        ey= 1./2.*c2_err/pow(c2, 1./2.) ;
        //cout<<c2<<"\t"<<y<<"\t"<<ey<<"\n";
        if(c2>0. && ey>0.){
            //cout<<c2<<"\t"<<y<<"\t"<<ey<<"\n";
            prof_results->Fill(x, y, 1/ey);
        }
    }
    
    return;
}

void calculate_v_2_4_int(TProfile *prof_results, TProfile *prof_2pc, TProfile *prof_4pc){
    
    Double_t x;
    Double_t y;
    Double_t ey;
    
    for (Int_t i=1;i<=prof_4pc->GetNcells()-4;++i) {
        x= prof_4pc->GetBinCenter(i);
        
        Double_t c4=c_4(prof_4pc->GetBinContent(i), prof_2pc->GetBinContent(i));
        Double_t c4_err= c_4_err(prof_4pc->GetBinContent(i), prof_2pc->GetBinContent(i), prof_4pc->GetBinError(i), prof_2pc->GetBinError(i));
        
        y= pow(-c4,1./4.)  ;
        ey= 1./4.*c4_err/pow(-c4, 3./4.) ;
        //cout<<c4<<"\t"<<y<<"\t"<<ey<<"\n";
        if(c4<0. && ey>0.){
            
            prof_results->Fill(x, y, 1/ey);
        }
    }
    
    return;
}

void calculate_v_2_6_int(TProfile *prof_results, TProfile *prof_2pc, TProfile *prof_4pc, TProfile *prof_6pc ){
    
    Double_t x;
    Double_t y;
    Double_t ey;
    
    for (Int_t i=1;i<=prof_4pc->GetNcells()-4;++i) {
        x= prof_4pc->GetBinCenter(i);
        
        Double_t c6=c_6(prof_6pc->GetBinContent(i), prof_4pc->GetBinContent(i), prof_2pc->GetBinContent(i));
        Double_t c6_err= c_6_err(prof_6pc->GetBinContent(i), prof_4pc->GetBinContent(i), prof_2pc->GetBinContent(i), prof_6pc->GetBinError(i), prof_4pc->GetBinError(i), prof_2pc->GetBinError(i));
        
        y= pow(1./4.*c6,1./6.)  ;
        ey= 1./6.*1./4.*pow(1./4.*c6,-5./6.)*c6_err;
        //cout<<c4<<"\t"<<y<<"\t"<<ey<<"\n";
        if(c6>0. && ey>0.){
            
            prof_results->Fill(x, y, 1/ey);
        }
    }
    
    return;
}

void calculate_v_2_8_int(TProfile *prof_results, TProfile *prof_2pc, TProfile *prof_4pc, TProfile *prof_6pc, TProfile *prof_8pc ){
    
    Double_t x;
    Double_t y;
    Double_t ey;
    
    for (Int_t i=1;i<=prof_4pc->GetNcells()-4;++i) {
        x= prof_4pc->GetBinCenter(i);
        
        Double_t c8=c_8(prof_8pc->GetBinContent(i),prof_6pc->GetBinContent(i), prof_4pc->GetBinContent(i), prof_2pc->GetBinContent(i));
        Double_t c8_err= c_8_err(prof_8pc->GetBinContent(i), prof_6pc->GetBinContent(i), prof_4pc->GetBinContent(i), prof_2pc->GetBinContent(i), prof_8pc->GetBinError(i), prof_6pc->GetBinError(i), prof_4pc->GetBinError(i), prof_2pc->GetBinError(i));
        
        y= pow(1./33.*c8,1./8.)  ;
        ey= 1./8.*1./33.*pow(1./33.*c8,-7./8.)*c8_err;
        cout<<c8<<"\t"<<c8_err<<"\t"<<y<<"\t"<<ey<<"\n";
        if(c8>0. && c8_err>0.){
            
            prof_results->Fill(x, y, 1/ey);
        }
    }
    
    return;
}


TCanvas* represent(TProfile* Mean_profile_int_2, TProfile* Mean_profile_int_4, TProfile* Mean_profile_int_6, TProfile* Mean_profile_int_8, TProfile* Error_profile_int_2, TProfile* Error_profile_int_4,  TProfile* Error_profile_int_6,  TProfile* Error_profile_int_8){
    
    TCanvas *c1 = new TCanvas("c1","multigraph",200, 10,700,500);
    
    c1->Divide(1,2,0,0);
    c1->cd(1);
    c1->GetPad(1)->SetRightMargin(.01);

    
    TMultiGraph *mg1 = new TMultiGraph();
  
  
    
    //!mg->GetXaxis()->SetTitle("Centrality");
    //!mg1->GetYaxis()->SetTitle("v_{2}");
    
    Double_t x[Mean_profile_int_2->GetNcells()-2];
    Double_t y2[Mean_profile_int_2->GetNcells()-2];
    Double_t ey2[Mean_profile_int_2->GetNcells()-2];
    Double_t y4[Mean_profile_int_4->GetNcells()-2];
    Double_t ey4[Mean_profile_int_4->GetNcells()-2];
    Double_t ry4[Mean_profile_int_4->GetNcells()-2];
    Double_t rye4[Mean_profile_int_4->GetNcells()-2];
    Double_t y6[Mean_profile_int_6->GetNcells()-2];
    Double_t ey6[Mean_profile_int_6->GetNcells()-2];
    Double_t ry6[Mean_profile_int_6->GetNcells()-2];
    Double_t rye6[Mean_profile_int_6->GetNcells()-2];
    Double_t y8[Mean_profile_int_8->GetNcells()-2];
    Double_t ey8[Mean_profile_int_8->GetNcells()-2];
    Double_t ry8[Mean_profile_int_8->GetNcells()-2];
    Double_t rye8[Mean_profile_int_8->GetNcells()-2];

    for(Int_t i=0; i< Mean_profile_int_2->GetNcells()-2; i++){
        x[i]= Mean_profile_int_2->GetBinCenter(i+1);
        
        y2[i]= Mean_profile_int_2->GetBinContent(i+1);
        ey2[i] = Error_profile_int_2->GetBinError(i+1);
        
        y4[i]= Mean_profile_int_4->GetBinContent(i+1);
        ey4[i] = Error_profile_int_4->GetBinError(i+1);
        ry4[i] = y4[i]/v24Run2[i];
        rye4[i] = sqrt(pow(ey4[i]/v24Run2[i],2.)+ pow(y4[i]*v24Run2CombErr[i]/(v24Run2[i]*v24Run2[i]),2.));
        
        y6[i]= Mean_profile_int_6->GetBinContent(i+1);
        ey6[i] = Error_profile_int_6->GetBinError(i+1);
        ry6[i] = y6[i]/v26Run2[i];
        rye6[i] = sqrt(pow(ey6[i]/v26Run2[i],2.)+ pow(y6[i]*v26Run2CombErr[i]/(v26Run2[i]*v26Run2[i]),2.));
        
        y8[i]= Mean_profile_int_8->GetBinContent(i+1);
        ey8[i] = Error_profile_int_8->GetBinError(i+1);
        ry8[i] = y8[i]/v28Run2[i];
        rye8[i] = sqrt(pow(ey8[i]/v28Run2[i],2.)+ pow(y8[i]*v28Run2CombErr[i]/(v28Run2[i]*v28Run2[i]),2.));
    }
    
    /*
    auto *gr2 = new TGraphErrors(Mean_profile_int_2->GetNcells()-2,x,y2,0,ey2);
    
    //!gr2->SetTitle("2" "centrality; pT [GeV]; v_2 ");
    gr2->SetMarkerStyle(kCircle);
    gr2->SetMarkerColor(kRed);
    gr2->SetMarkerStyle(21);
    
    */
    auto *gr4 = new TGraphErrors(Mean_profile_int_4->GetNcells()-2,x,y4,0,ey4);
    
    //!gr4->SetTitle("4" "centrality; pT [GeV]; v2");
    gr4->SetMarkerStyle(kCircle);
    gr4->SetMarkerColor(kBlue);
    gr4->SetMarkerStyle(21);
    
    auto *gr6 = new TGraphErrors(Mean_profile_int_6->GetNcells()-2,x,y6,0,ey6);
    
    //!gr6->SetTitle("6" "centrality; pT [GeV]; v2");
    gr6->SetMarkerStyle(kCircle);
    gr6->SetMarkerColor(kBlack);
    gr6->SetMarkerStyle(21);
    
    auto *gr8 = new TGraphErrors(Mean_profile_int_8->GetNcells()-2,x,y8,0,ey8);
    
    //!gr8->SetTitle("8" "centrality; pT [GeV]; v2");
    gr8->SetMarkerStyle(kCircle);
    gr8->SetMarkerColor(kRed);
    gr8->SetMarkerStyle(21);
    
    //mg1->Add(gr2);
    mg1->Add(gr4);
    mg1->Add(gr6);
    mg1->Add(gr8);
    
    mg1->SetTitle("Values; V0M; v_{2}");
  //!  mg1->GetYaxis()->SetTitle("v_{2}");
    mg1->GetYaxis()->SetTitleSize(0.03);
    //! mg1->SetTitle("30-40%");
    mg1->Draw("ap");
    
    TLegend leg(.1,.7,.3,.9,"Alice Pb-Pb #sqrt{s_{NN}}=5.02 TeV, No gap ");
    leg.SetFillColor(0);
    //!leg.AddEntry(gr2,"2 pc");
    leg.AddEntry(gr4,"4 pc");
    leg.AddEntry(gr6,"6 pc");
    leg.AddEntry(gr8,"8 pc");
    leg.DrawClone("Same");
    
    
    c1->cd(2);
    
    c1->GetPad(2)->SetRightMargin(.01);
    
    TMultiGraph *mg2 = new TMultiGraph();
    
    auto *gr4_r = new TGraphErrors(Mean_profile_int_2->GetNcells()-2,x,ry4,0,rye4);
    
    //!gr4->SetTitle("4" "centrality; pT [GeV]; v2");
    gr4_r->SetMarkerStyle(kCircle);
    gr4_r->SetMarkerColor(kBlue);
    gr4_r->SetMarkerStyle(21);

    auto *gr6_r = new TGraphErrors(Mean_profile_int_2->GetNcells()-2,x,ry6,0,rye6);
    
    //!gr4_r->SetTitle("4" "centrality; pT [GeV]; v2");
    gr6_r->SetMarkerStyle(kCircle);
    gr6_r->SetMarkerColor(kBlack);
    gr6_r->SetMarkerStyle(21);
    
    
    auto *gr8_r = new TGraphErrors(Mean_profile_int_2->GetNcells()-2,x,ry8,0,rye8);
    
    //!gr4_r->SetTitle("4" "centrality; pT [GeV]; v2");
    gr8_r->SetMarkerStyle(kCircle);
    gr8_r->SetMarkerColor(kRed);
    gr8_r->SetMarkerStyle(21);
    
    
    mg2->Add(gr4_r);
    mg2->Add(gr6_r);
    mg2->Add(gr8_r);
    
    mg2->SetTitle(" Ratio plot; V0M; y/y_pub");
  //!  mg2->GetYaxis()->SetTitle("v_{2}");
    mg2->GetYaxis()->SetTitleSize(0.03);
    //! mg2->SetTitle("30-40%");
    mg2->Draw("ap");
    
 
    c1->Update();
    c1->Range( 0.2, 0., 5., 2. );
    
    auto *l = new TLine(0.,1.,90.,1.);
    
    l->SetLineWidth(2);
    l->Draw();
     
    
    return c1;
}



TProfile2D* load_results_2D(const char * fname, const char * particles){

    TFile* inputFile = TFile::Open(fname, "READ");

    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    TProfile2D *profile = (TProfile2D*)inputList->FindObject(particles);
    
    return profile;
}

TProfile* load_results(const char * fname, const char * particles){

    TFile* inputFile = TFile::Open(fname, "READ");

    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    TProfile *profile = (TProfile*)inputList->FindObject(particles);
    
    return profile;
}

void macro(){
    
    
    Double_t bin_edges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
    
    const char *name_2[10]= {"2_1", "2_2", "2_3", "2_4","2_5", "2_6", "2_7", "2_8", "2_9" ,"2_10"}; // TODO: FIX THIS FUCKIN SHIT
    const char *name_4[10]= {"4_1", "4_2", "4_3", "4_4","4_5", "4_6", "4_7", "4_8", "4_9" ,"4_10"}; // TODO: FIX THIS FUCKIN SHIT
    const char *name_6[10]= {"6_1", "6_2", "6_3", "6_4","6_5", "6_6", "6_7", "6_8", "6_9" ,"6_10"}; // TODO: FIX THIS FUCKIN SHIT
    const char *name_8[10]= {"8_1", "8_2", "8_3", "8_4", "8_5", "8_6", "8_7", "8_8", "8_9" ,"8_10"}; // TODO: FIX THIS FUCKIN SHIT
    
    TProfile* Mean_profile_int_2 = new TProfile(" Mean 2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{2}",10, bin_edges);
    TProfile* Mean_profile_int_4 = new TProfile(" Mean 4 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",10, bin_edges);
    TProfile* Mean_profile_int_6 = new TProfile(" Mean 6 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{6}",10, bin_edges);
    TProfile* Mean_profile_int_8 = new TProfile(" Mean 8 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",10, bin_edges);
    
    
    TProfile* Error_profile_int_2 = new TProfile(" Mean 2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{2}",10, bin_edges);
    TProfile* Error_profile_int_4 = new TProfile(" Error 4 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",10, bin_edges);
    TProfile* Error_profile_int_6 = new TProfile(" Error 6 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{6}",10, bin_edges);
    TProfile* Error_profile_int_8 = new TProfile(" Error 8 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{8}",10, bin_edges);
    
    
   //! TCanvas* stuff2= calculate_v_2_4_diff(Error_Profile, projected_prof2,prof_2pc, projected_prof4, prof_4pc, centr, kRed, v2_4_pub, v2_4_err);
    
    TFile* inputFile = TFile::Open("AnalysisResults.root", "READ");

    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    
    for(Int_t i(0); i<10; i++){
        
        TProfile *prof_2pc_int = (TProfile*)inputList->FindObject(name_2[i]);
        TProfile *prof_4pc_int = (TProfile*)inputList->FindObject(name_4[i]);
        TProfile *prof_6pc_int = (TProfile*)inputList->FindObject(name_6[i]);
        TProfile *prof_8pc_int = (TProfile*)inputList->FindObject(name_8[i]);
        
        calculate_v_2_2_int(Error_profile_int_2, prof_2pc_int);
        calculate_v_2_4_int(Error_profile_int_4, prof_2pc_int, prof_4pc_int);
        calculate_v_2_6_int(Error_profile_int_6, prof_2pc_int, prof_4pc_int, prof_6pc_int);
        calculate_v_2_8_int(Error_profile_int_8, prof_2pc_int, prof_4pc_int, prof_6pc_int, prof_8pc_int);
    }
    
    TProfile *prof_2pc_int = (TProfile*)inputList->FindObject("results_value_2");
    TProfile *prof_4pc_int = (TProfile*)inputList->FindObject("results_value_4");
    TProfile *prof_6pc_int = (TProfile*)inputList->FindObject("results_value_6");
    TProfile *prof_8pc_int = (TProfile*)inputList->FindObject("results_value_8");
    
    //prof_4pc_int->Draw();
    
    calculate_v_2_2_int(Mean_profile_int_2, prof_2pc_int);
    calculate_v_2_4_int(Mean_profile_int_4, prof_2pc_int, prof_4pc_int);
    calculate_v_2_6_int(Mean_profile_int_6, prof_2pc_int, prof_4pc_int, prof_6pc_int);
    calculate_v_2_8_int(Mean_profile_int_8, prof_2pc_int, prof_4pc_int, prof_6pc_int, prof_8pc_int);
    
    //!Mean_profile_int_8->Draw();
    
    
    //!auto c1= represent_diff(Mean_profile_diff_2, Mean_profile_diff_4,Error_profile_diff_2, Error_profile_diff_4, v2_4_pub, v2_4_err);
    
    auto c2= represent(Mean_profile_int_2, Mean_profile_int_4, Mean_profile_int_6, Mean_profile_int_8, Error_profile_int_2, Error_profile_int_4, Error_profile_int_6, Error_profile_int_8);
    
    
    
    
    
    //!stuff2->Draw();
    
   //! Mean_profile_diff_2->Draw();
    
    
    
    

    return;
}

