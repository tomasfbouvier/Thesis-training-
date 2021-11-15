#include "published.h"

void calculate_v_2_2_diff(TProfile *prof_results, TProfile *prof_2pc_diff, TProfile *prof_2pc, Int_t centr, Color_t color){
    
    Double_t x;
    Double_t y;
    Double_t ey;

    for (Int_t i=0;i<=prof_2pc_diff->GetNcells()-3;++i) {
        x= prof_2pc_diff->GetBinCenter(i+1);
        y= prof_2pc_diff->GetBinContent(i+1)/sqrt(prof_2pc->GetBinContent(centr));
        ey= sqrt((prof_2pc_diff->GetBinError(i+1)/sqrt(prof_2pc->GetBinContent(centr)))*(prof_2pc_diff->GetBinError(i+1)/sqrt(prof_2pc->GetBinContent(centr)))+ (prof_2pc_diff->GetBinContent(i+1)*prof_2pc->GetBinError(centr)/(2*pow(prof_2pc->GetBinContent(centr),3./2.)))*(prof_2pc_diff->GetBinContent(i+1)*prof_2pc->GetBinError(centr)/(2*pow(prof_2pc->GetBinContent(centr),3./2.)))) ;
        
        
        //!cout<<y<<"\t"<<ey<<"\n";
        if(ey>0.){
            
            prof_results->Fill(x, y, 1/ey);
        }
    }

    return;
}

Double_t d_4(Double_t pc_4_diff,Double_t pc_2_diff, Double_t pc_2_ref){
    return pc_4_diff - 2*pc_2_diff*pc_2_ref;
}

Double_t c_4(Double_t pc_4_ref, Double_t pc_2_ref){
    return pc_4_ref - 2*pc_2_ref*pc_2_ref;
}

Double_t d_4_err(Double_t pc_4_diff, Double_t pc_2_diff, Double_t pc_2_ref, Double_t pc_4_diff_err, Double_t pc_2_diff_err, Double_t pc_2_ref_err){

    Double_t aux= pc_4_diff_err*pc_4_diff_err + (2*pc_2_diff_err*pc_2_ref)*(2*pc_2_diff*pc_2_ref_err)*(2*pc_2_diff_err*pc_2_ref)*(2*pc_2_diff*pc_2_ref_err);
    
    return sqrt(aux);
}

Double_t c_4_err(Double_t pc_4_ref, Double_t pc_2_ref, Double_t pc_4_ref_err, Double_t pc_2_ref_err){
    
    Double_t aux= pc_4_ref_err*pc_4_ref_err + (2*pc_2_ref_err*pc_2_ref)*(2*pc_2_ref*pc_2_ref_err)*(2*pc_2_ref_err*pc_2_ref)*(2*pc_2_ref*pc_2_ref_err);
    
    return sqrt(aux);
}



void calculate_v_2_4_diff(TProfile *prof_results, TProfile *prof_2pc_diff, TProfile *prof_2pc, TProfile *prof_4pc_diff, TProfile *prof_4pc, Int_t centr, Color_t color, Double_t published[], Double_t published_err[]){
    
    Double_t x;
    Double_t y;
    Double_t ey;
    
    for (Int_t i=0;i<=prof_4pc_diff->GetNcells()-3;++i) {
        x= prof_4pc_diff->GetBinCenter(i+1);
        
        Double_t c4=c_4(prof_4pc->GetBinContent(centr), prof_2pc->GetBinContent(centr));
        Double_t d4= d_4(prof_4pc_diff->GetBinContent(i+1), prof_2pc_diff->GetBinContent(i+1),  prof_2pc->GetBinContent(centr));
        Double_t d4_err= d_4_err(prof_4pc_diff->GetBinContent(i+1), prof_2pc_diff->GetBinContent(i+1), prof_2pc->GetBinContent(centr) , prof_4pc_diff->GetBinError(i+1), prof_2pc_diff->GetBinError(i+1),  prof_2pc->GetBinError(centr));
        Double_t c4_err= c_4_err(prof_4pc->GetBinContent(centr), prof_2pc->GetBinContent(centr), prof_4pc->GetBinError(centr), prof_2pc->GetBinError(centr));
        
        y= - d4/pow(-c4,3./4.)  ;
        
        

        
        ey= sqrt((d4_err/pow(-c4,3./4.))*(d4_err/pow(-c4,3./4.))+ (d4*c4_err/(3./4.*pow(-c4,7./4.)))*(d4*c4_err/(3./4.*pow(-c4,7./4.)))) ;
        
        if(c4<0. && ey>0.){
            //!cout<<c4<<"\t"<<d4<<"\t"<<y<<"\t"<<ey<<"\n";
            prof_results->Fill(x, y, 1/ey);
        }
    }
    
    return;
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
        cout<<x<<"\n";
        //cout<<c4<<"\t"<<y<<"\t"<<ey<<"\n";
        if(c4<0. && ey>0.){
            
            prof_results->Fill(x, y, 1/ey);
        }
    }
    
    return;
}

TCanvas* represent_diff(TProfile* Mean_profile_diff_2, TProfile* Mean_profile_diff_4, TProfile* Error_profile_diff_2, TProfile* Error_profile_diff_4, Double_t v2_4_pub[], Double_t v2_4_err[]){
    
    TCanvas *c1 = new TCanvas("c1","multigraph",200, 10,700,500);
    
    c1->Divide(1,2,0,0);
    c1->cd(1);
    c1->GetPad(1)->SetRightMargin(.01);

    
    TMultiGraph *mg2 = new TMultiGraph();
    
    //!mg->GetXaxis()->SetTitle("Centrality");
    //!mg2->GetYaxis()->SetTitle("v_{2}");
    
    Double_t x[Mean_profile_diff_2->GetNcells()-2];
    Double_t y2[Mean_profile_diff_2->GetNcells()-2];
    Double_t ey2[Mean_profile_diff_2->GetNcells()-2];
    Double_t y4[Mean_profile_diff_4->GetNcells()-2];
    Double_t ey4[Mean_profile_diff_4->GetNcells()-2];
    Double_t ry4[Mean_profile_diff_4->GetNcells()-2];
    Double_t rye4[Mean_profile_diff_4->GetNcells()-2];

    for(Int_t i=0; i< Mean_profile_diff_2->GetNcells()-2; i++){
        x[i]= Mean_profile_diff_2->GetBinCenter(i+1);
        y2[i]= Mean_profile_diff_2->GetBinContent(i+1);
        ey2[i] = Error_profile_diff_2->GetBinError(i+1);
        y4[i]= Mean_profile_diff_4->GetBinContent(i+1);
        ey4[i] = Error_profile_diff_4->GetBinError(i+1);
        
        ry4[i] = y4[i]/ v2_4_pub[i];
        rye4[i] = sqrt(ey4[i]/v2_4_pub[i]*ey4[i]/v2_4_pub[i]+ y4[i]*v2_4_err[i]/(v2_4_pub[i]*v2_4_pub[i])*y4[i]*v2_4_err[i]/(v2_4_pub[i]*v2_4_pub[i]));
        //cout<<x[i]<<"\t"<<y2[i]<<"\t"<<ey2[i]<<"\n";
    }
    
    auto *gr2 = new TGraphErrors(Mean_profile_diff_2->GetNcells()-2,x,y2,0,ey2);
    
    //!gr2->SetTitle("2" "centrality; pT [GeV]; v_2 ");
    gr2->SetMarkerStyle(kCircle);
    gr2->SetMarkerColor(kRed);
    gr2->SetMarkerStyle(21);
    
    auto *gr4 = new TGraphErrors(Mean_profile_diff_2->GetNcells()-2,x,y4,0,ey4);
    
    //!gr4->SetTitle("4" "centrality; pT [GeV]; v2");
    gr4->SetMarkerStyle(kCircle);
    gr4->SetMarkerColor(kBlue);
    gr4->SetMarkerStyle(21);
    
    mg2->Add(gr2);
    mg2->Add(gr4);
    
    mg2->SetTitle("30-40%; pT (GeV/c); v_{2}");
  //!  mg2->GetYaxis()->SetTitle("v_{2}");
    mg2->GetYaxis()->SetTitleSize(0.06);
    //! mg2->SetTitle("30-40%");
    mg2->Draw("ap");
    
    TLegend leg(.1,.7,.3,.9,"Alice Pb-Pb #sqrt{s_{NN}}=5.02 TeV, ");
    leg.SetFillColor(0);
    leg.AddEntry(gr2,"2 pc, |#Delta #eta| = 0");
    leg.AddEntry(gr4,"4 pc, |#Delta #eta| = 0");
    leg.DrawClone("Same");
    
    c1->cd(2);
    
    c1->GetPad(2)->SetRightMargin(.01);
    
    auto *gr_ratio = new TGraphErrors(Mean_profile_diff_2->GetNcells()-2,x,ry4,0,rye4);
    
    gr_ratio->SetTitle("ratio plot; pT (GeV/c); v_{2}{4}/pub.");
    gr_ratio->SetMarkerStyle(kCircle);
    gr_ratio->SetMarkerColor(kBlue);
    gr_ratio->SetMarkerStyle(21);
    
    gr_ratio->Draw("ap");

    c1->Update();
    c1->Range( 0.2, 0., 5., 2. );
    
    auto *l = new TLine(0.,1.,4.9,1.);
    
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
    
    Int_t centr= 5;
    
    Double_t bin_edges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
    Double_t bin_edges_y[] =  {0.2,0.4,0.6,0.8,1.,1.25,1.5,1.75,2.,2.5,3.,3.5,4.,5.};
    
    Double_t v2_4_pub[]= {0.0402, 0.0661, 0.0895, 0.1103, 0.1285, 0.1457, 0.1664, 0.1796, 0.1936, 0.1911, 0.2022, 0.1819, 0.1652};
    Double_t v2_4_err[]= {0.0038, 0.0032, 0.0056, 0.0062, 0.0078, 0.0088, 0.0142, 0.0119, 0.0123, 0.0143, 0.23, 0.0439, 0.0555};

    const char *name_2[10]= {"2_1", "2_2", "2_3", "2_4","2_5", "2_6", "2_7", "2_8", "2_9" ,"2_10"}; // TODO: FIX THIS FUCKIN SHIT
    const char *name_4[10]= {"4_1", "4_2", "4_3", "4_4","4_5", "4_6", "4_7", "4_8", "4_9" ,"4_10"}; // TODO: FIX THIS FUCKIN SHIT

    
    TProfile* Mean_profile_int_2 = new TProfile(" Error 2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",10, bin_edges);
    TProfile* Mean_profile_int_4 = new TProfile(" Error 4 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",10, bin_edges);
    TProfile* Error_profile_int_2 = new TProfile(" Mean 2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",10, bin_edges);
    TProfile* Error_profile_int_4 = new TProfile(" Mean 4 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",10, bin_edges);
    
    TProfile* Mean_profile_diff_2 = new TProfile(" Error 2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",13, bin_edges_y);
    TProfile* Mean_profile_diff_4 = new TProfile(" Error 4 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",13, bin_edges_y);
    TProfile* Error_profile_diff_2 = new TProfile(" Mean 2 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",13, bin_edges_y);
    TProfile* Error_profile_diff_4 = new TProfile(" Mean 4 Pb-Pb 5.02GeV","Profile eta vs value; pT; v_2{4}",13, bin_edges_y);
    
    
   //! TCanvas* stuff2= calculate_v_2_4_diff(Error_Profile, projected_prof2,prof_2pc, projected_prof4, prof_4pc, centr, kRed, v2_4_pub, v2_4_err);
    
    TFile* inputFile = TFile::Open("AnalysisResults.root", "READ");

    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    
    
    
    for(Int_t i(0); i<10; i++){
        
        TProfile2D *prof_2pc = (TProfile2D*)inputList->FindObject(name_2[i]);
        TProfile *prof_2pc_int= (TProfile*) prof_2pc->ProfileX();
        TProfile *prof_2pc_diff= (TProfile*) prof_2pc->ProfileY(" ",centr,centr+1);
        
        TProfile2D *prof_4pc = (TProfile2D*)inputList->FindObject(name_4[i]);
        TProfile *prof_4pc_int= (TProfile*) prof_4pc->ProfileX();
        TProfile *prof_4pc_diff= (TProfile*) prof_4pc->ProfileY(" ",centr,centr+1);
        
        calculate_v_2_2_int(Error_profile_int_2, prof_2pc_int);
        calculate_v_2_4_int(Error_profile_int_4, prof_2pc_int, prof_4pc_int);
    
        calculate_v_2_2_diff(Error_profile_diff_2, prof_2pc_diff , prof_2pc_int, centr, kRed);
        calculate_v_2_4_diff(Error_profile_diff_4, prof_2pc_diff , prof_2pc_int, prof_4pc_diff, prof_4pc_int, centr, kRed, v2_4_pub, v2_4_err);
    }
    
    TProfile2D *prof_2pc = (TProfile2D*)inputList->FindObject("Big_2");
    TProfile *prof_2pc_int= (TProfile*) prof_2pc->ProfileX();
    cout<<"Effective entries 2"<<prof_2pc_int->GetEffectiveEntries()<<"\n";
    TProfile *prof_2pc_diff= (TProfile*) prof_2pc->ProfileY(" ",centr,centr+1);
    
    
    
    TProfile2D *prof_4pc = (TProfile2D*)inputList->FindObject("Big_4");
    TProfile *prof_4pc_int= (TProfile*) prof_4pc->ProfileX();
    cout<<"Effective entries 4"<<prof_4pc_int->GetEffectiveEntries()<<"\n";
    TProfile *prof_4pc_diff= (TProfile*) prof_4pc->ProfileY(" ",centr,centr+1);
    
    //prof_4pc_int->Draw();
    
    calculate_v_2_2_int(Mean_profile_int_2, prof_2pc_int);
    calculate_v_2_4_int(Mean_profile_int_4, prof_2pc_int, prof_4pc_int);
    
    calculate_v_2_2_diff(Mean_profile_diff_2, prof_2pc_diff , prof_2pc_int, centr, kRed);
    calculate_v_2_4_diff(Mean_profile_diff_4, prof_2pc_diff , prof_2pc_int, prof_4pc_diff, prof_4pc_int, centr, kRed, v2_4_pub, v2_4_err);
    
    
    //!auto c1= represent_diff(Mean_profile_diff_2, Mean_profile_diff_4,Error_profile_diff_2, Error_profile_diff_4, v2_4_pub, v2_4_err);
    
    auto c2= represent_diff(Mean_profile_int_2, Mean_profile_int_4,Error_profile_int_2, Error_profile_int_4, v2_4_pub, v2_4_err);
    
    
    //!stuff2->Draw();
    
   //! Mean_profile_diff_2->Draw();
    
    
    
    

    return;
}

