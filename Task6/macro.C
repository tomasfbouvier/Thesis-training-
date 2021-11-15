#include "macro.h"

void calculate_SC(TProfile *prof_results_rho, TProfile *prof_results_chi, TProfile *prof_Nom,  TProfile *prof_Denom, TProfile *prof_vn_gap){
    
    Double_t x;
    Double_t y;
    Double_t ey;
    
    Double_t y_rho;
    Double_t ey_rho;
    
    Double_t y_chi;
    Double_t ey_chi;
    
    
    for (Int_t i=0;i<=prof_Nom->GetNcells()-3;++i) {
        
        x= prof_Nom->GetBinCenter(i+1);
        if(prof_Denom->GetBinContent(i+1)>0.){
            y= prof_Nom->GetBinContent(i+1)/sqrt(prof_Denom->GetBinContent(i+1));
            ey= sqrt(pow(prof_Nom->GetBinError(i+1)/sqrt(prof_Denom->GetBinContent(i+1)),2) + pow( 1./2.*prof_Nom->GetBinContent(i+1)/pow(prof_Denom->GetBinContent(i+1), 3./2.)*prof_Denom->GetBinError(i+1),2));
            
            y_rho= y/sqrt(prof_vn_gap->GetBinContent(i+1));
            ey_rho= sqrt( pow(ey/sqrt(prof_vn_gap->GetBinContent(i+1)),2)+  pow(y/pow(prof_vn_gap->GetBinContent(i+1), 3./2.)*prof_vn_gap->GetBinError(i+1) ,2)) ;
            
            if(ey_rho>0.){
                prof_results_rho->Fill(x, y_rho, 1.);
            }
            y_chi= y/sqrt(prof_Denom->GetBinContent(i+1));
            ey_chi= sqrt(pow(ey/sqrt(prof_Denom->GetBinContent(i+1)),2) + pow( 1./2.*y/pow(prof_Denom->GetBinContent(i+1), 3./2.)*prof_Denom->GetBinError(i+1),2));;
            
            
            if(ey_chi>0.){
                prof_results_chi->Fill(x, y_chi, 1.);
            }
        }
    }
        /*
        y/=(prof_2pc_m_gap->GetBinContent(i)*prof_2pc_n_gap->GetBinContent(i));
        ey = sqrt((ey/(prof_2pc_m_gap->GetBinContent(i)*prof_2pc_n_gap->GetBinContent(i)))*(ey/(prof_2pc_m_gap->GetBinContent(i)*prof_2pc_n_gap->GetBinContent(i))) +  (y*prof_2pc_m_gap->GetBinError(i)/(pow(prof_2pc_m_gap->GetBinContent(i),2)*prof_2pc_n_gap->GetBinContent(i)))*(y*prof_2pc_m_gap->GetBinError(i)/(pow(prof_2pc_m_gap->GetBinContent(i),2)*prof_2pc_n_gap->GetBinContent(i)))+ (y*prof_2pc_n_gap->GetBinError(i)/(pow(prof_2pc_n_gap->GetBinContent(i),2)*prof_2pc_m_gap->GetBinContent(i)))*(y*prof_2pc_n_gap->GetBinError(i)/(pow(prof_2pc_n_gap->GetBinContent(i),2)*prof_2pc_m_gap->GetBinContent(i))));
        
        if(ey>0.){
            //cout<<c2<<"\t"<<y<<"\t"<<ey<<"\n";
            prof_results_NSC->Fill(x, y, 1/ey);
        }
    }
         */
    return;
}

/*

void get_published_values(const char fname){
    
    TFile* inputFile = TFile::Open(fname, "READ");

    //!inputFile->cd();
    TList* inputList = (TList*) gDirectory->Get("fname");
    TProfile *profile_y = (TProfile*)inputList->FindObject("Hist1D_y1");
    TProfile *profile_ey = (TProfile*)inputList->FindObject("Hist1D_y1_e1");
    
    Double_t    y[profile_y->GetNcells()];
    Double_t    ey[profile_y->GetNcells()];
    
    for(Int_t i=0; i<profile_y->GetNcells(); i++){
        y[i] = profile_y->GetBinContent(i+1);
        ey[i] = profile_y->GetBinError(i+1);
        
    }
    
}
*/


void represent(TList* fOutputList, TList* fErrorList){
    


    /*
    TMultiGraph *mg1 = new TMultiGraph();
    
    TLegend leg(.1,.7,.3,.9,"Alice Pb-Pb #sqrt{s_{NN}}=5.02 TeV, No gap ");
    leg.SetFillColor(0);
    //!leg.AddEntry(gr2,"2 pc");
    
    //!mg1->GetYaxis()->SetLimits(-0.0000002, 0.0000002);
   
    */
    
    //!TCanvas *plots[fOutputList->GetEntries()];
    
    const Double_t mult_factors[] = {1.,1.,1.,1.,1.,1.};
    
    TCanvas *c1 = new TCanvas("c1","multipads", 800, 700);
    
    c1->Divide(3,2,0,0);
    
    Color_t colors[]={kRed,kBlue,kBlack, kGreen, kViolet, kOrange};
    
    for(Int_t j=0; j<6;j++){
        
        TProfile* Mean_profile = (TProfile*)fOutputList->At(j);
        TProfile* Error_profile = (TProfile*)fErrorList->At(j);
        
        Double_t x[Mean_profile->GetNcells()-2];
        Double_t y[Mean_profile->GetNcells()-2];
        Double_t ey[Mean_profile->GetNcells()-2];
        
        for(Int_t i=0; i< Mean_profile->GetNcells()-2; i++){
            x[i]= Mean_profile->GetBinCenter(i+1);
            y[i]= Mean_profile->GetBinContent(i+1)*mult_factors[j];
            ey[i] = Error_profile->GetBinError(i+1)*mult_factors[j];
            //!cout<<x[i]<<"\t"<<y[i]<<"\t"<<ey[i]<<"\t"<<"\n";
            //!cout<<x[i]<<"\t"<<y[i]<<"\t"<<ey[i]<<"\n";
        }
        
        c1->cd(j+1);
        
        auto gr= new TGraphErrors(Mean_profile->GetNcells()-2,x,y,0,ey);
        gr->SetTitle(profile_nom_names2[j]);
        gr->SetMarkerStyle(kCircle);
        gr->SetMarkerColor(colors[j]);
        gr->SetMarkerStyle(21);

        gr->Draw("AP");
    }


     
    

    return ;
}

void macro(){
    

    const Double_t bin_edges[] = {0.,5.,10.,20.,30.,40.,50., 60.};

    std::string names[]={"_1","_2","_3","_4","_5", "_6","_7","_8","_9","_10", "_MEAN" };

    TFile* inputFile = TFile::Open("AnalysisResults.root", "READ");
    
    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    
    for(Int_t j(0); j<6; j++){
        
        TProfile* Error_profile_chi= new TProfile(("Error_profile_chi"+profile_nom_names[j]).c_str(),"Profile eta vs value; V0M; rho",7, bin_edges);
        TProfile* Error_profile_rho= new TProfile(("Error_profile_rho"+profile_nom_names[j]).c_str(),"Profile eta vs value; V0M; rho",7, bin_edges);
        TProfile* Mean_profile_chi= new TProfile(("Mean_profile_chi"+profile_nom_names[j]).c_str(),"Profile eta vs value; V0M; rho",7, bin_edges);
        TProfile* Mean_profile_rho= new TProfile(("Mean_profile_rho"+profile_nom_names[j]).c_str(),"Profile eta vs value; V0M; rho",7, bin_edges);
        
        TFile* inputFile3 = TFile::Open(inputFile_names[j],"READ");
    
        
        inputFile3->cd("FlowExampleTask");
        TList* inputList3 = (TList*) gDirectory->Get("Output");
        
        for(Int_t i(0); i<11; i++){
            

            TProfile *prof_vn_gap = (TProfile*)inputList3->FindObject(name_vn[i]);
            TProfile *prof_nom = (TProfile*)inputList->FindObject((profile_nom_names[j]+ names[i]).c_str());
            TProfile *prof_denom = (TProfile*)inputList->FindObject((profile_denom_names[j]+ names[i]).c_str());
            
            if(i==10){
                calculate_SC(Mean_profile_rho, Mean_profile_chi, prof_nom, prof_denom, prof_vn_gap);
            }
            else{
                calculate_SC(Error_profile_rho, Error_profile_chi, prof_nom, prof_denom, prof_vn_gap);
            }
        }
        
        fMeanList_rho->Add(Mean_profile_rho);
        fErrorList_rho->Add(Error_profile_rho);
        fMeanList_chi->Add(Mean_profile_chi);
        fErrorList_chi->Add(Error_profile_chi);
        
        Error_profile_rho->Draw();
    }
    
    represent(fMeanList_rho, fErrorList_rho);
    return;
}
