void get_weights(){

    // opening file
    TFile* inputFile = TFile::Open("/Users/tomasfernandezbouvier/Desktop/Thesis-training-/histos/AnalysisResults.root", "READ");
    if(!inputFile || !inputFile->IsOpen()) { printf("Problem with the input file \n"); return;}


    // find and open tprofile with <<2>> (this is NOT v2)
    inputFile->cd("FlowExampleTask");
    TList* inputList = (TList*) gDirectory->Get("Output");
    TH1* fHistPhi = (TH1*)inputList->FindObject("fHistPhi");
    
    new TCanvas();
    
    fHistPhi->Draw();
    
    TH1D* fHistWeightsPhi =(TH1D*) fHistPhi->Clone("hWeight");
    
    //!TH1D*  fHistWeightsPhi=  new TH1D("WeightsPhi", "WeightsPhi; phi", 20, -0.5, 7);
    
    
   
    
    Double_t Nmax= fHistPhi->GetBinContent(fHistPhi->GetMaximumBin());
    
    for(Int_t i=0; i< fHistWeightsPhi->GetNbinsX()+2; i++){
        Double_t content= fHistPhi->GetBinContent(i);
        if(content>0){
            fHistWeightsPhi->SetBinContent(i, Nmax/content);
        }
        else{
            fHistWeightsPhi->SetBinContent(i, 1);
        }
        
    }
    
    new TCanvas();
     
    fHistWeightsPhi->Draw();
    /*
    TFile outputFile("/Users/tomasfernandezbouvier/Desktop/alice/Thesis-training-/Weights.root","RECREATE");
    fHistWeightsPhi->Write();

   // Closing the ROOT file.
    outputFile.Close();
    */
    TList* weights = new TList();
    weights->SetOwner(kTRUE);
    weights->Add(fHistWeightsPhi);

    TH1D* weightsPlaceHolder = (TH1D*)fHistWeightsPhi->Clone("weightsPlaceHolder");
    weightsPlaceHolder->Scale(3.);
    weights->Add(weightsPlaceHolder);

    TFile* output = new TFile("weight.root", "RECREATE");
    weights->Write("weightsList",TObject::kSingleKey+TObject::kOverwrite);

    
    
    
    return;
}
